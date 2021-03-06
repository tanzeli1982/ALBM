/*
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: subfile.c 1468 2013-10-26 16:53:18Z wkliao $ */
#define _GNU_SOURCE
#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif 
#include "subfile.h"
#ifdef TAU_SSON
#include "TAU.h"
#endif
#include <stdlib.h>
#include <math.h>
#include "ncmpidtype.h"

static int DEBUG = 1;
enum {ONE, BALANCED};
int min_ndims = 1;

/* set default values for the following values */
int delegate_scheme = BALANCED; /* default: any proc can be delegate proc */
int is_partitioned = 0;

#define check_err(fn_name_)                                                   \
    do {                                                                      \
        if (err != MPI_SUCCESS) {                                             \
            errs++;                                                           \
            if (DEBUG) {                                                      \
                int len_;                                                     \
                char err_str_[MPI_MAX_ERROR_STRING];                          \
                MPI_Error_string(err, err_str_, &len_);                       \
                fprintf(stderr, #fn_name_ " failed at line %d, err=%d: %s\n", \
                        __LINE__, err, err_str_);                             \
            }                                                                 \
        }                                                                     \
    } while (0)

static int ncmpii_itoa(int val, char* buf)
{
    const unsigned int radix = 10;

    char* p;
    unsigned int a;        //every digit
    int len;
    char* b;            //start of the digit char
    char temp;
    unsigned int u;

    p = buf;

    if (val < 0) {
        *p++ = '-';
        val = 0 - val;
    }
    u = (unsigned int)val;

    b = p;

    do {
        a = u % radix;
        u /= radix;

        *p++ = a + '0';

    } while (u > 0);

    len = (int)(p - buf);

    *p-- = 0;

    //swap
    do {
        temp = *p;
        *p = *b;
        *b = temp;
        --p;
        ++b;
    } while (b < p);

    return len;
}

int ncmpii_subfile_create(NC *ncp, int *ncidp)
{
    int status=NC_NOERR;
    int myrank, nprocs, color;
    char path_sf[1024];
    MPI_Comm comm_sf;
    MPI_Info info;
    MPI_Info_create(&info);
    double ratio;

    MPI_Comm_rank (ncp->nciop->comm, &myrank);
    MPI_Comm_size (ncp->nciop->comm, &nprocs);
#ifdef SUBFILE_DEBUG
    if (myrank == 0)
      printf("%s: rank(%d): nprocs=%d, ncp->nc_num_subfiles=%d\n", 
           __func__, myrank, nprocs, ncp->nc_num_subfiles);
#endif

    /* split the orignial comm to subcomm */
    if (nprocs > ncp->nc_num_subfiles) {
        ratio = (double)nprocs/(double)(ncp->nc_num_subfiles);
	color = (int)((double)myrank/ratio);
    }
    else 
        color = myrank%ncp->nc_num_subfiles;

#ifdef SUBFILE_DEBUG
    printf("rank(%d): color=%d\n", myrank, color);
#endif
    //key = myrank/comm_size;
    
    /* TODO: fix error when using generated key value.
     * for now, just passing 0 value. */
    MPI_Comm_split(ncp->nciop->comm, color, myrank, &comm_sf);

    sprintf(path_sf, "%s.subfile_%i.%s", ncp->nciop->path, color, "nc");

    //MPI_Info_set(info, "romio_lustre_start_iodevice", offset);
    //MPI_Info_set(info, "striping_factor", "1");

    status = ncmpi_create(comm_sf, path_sf, NC_CLOBBER|NC_64BIT_DATA, info, ncidp);  
    if (status != NC_NOERR) {
      if (myrank == 0) fprintf(stderr, "%s: error in creating file(%s): %s\n", 
		__func__, path_sf, ncmpi_strerror(status));
        return status;
    }
    MPI_Info_free(&info);

    return status;
}

int
ncmpii_subfile_open(NC *ncp, int *ncidp)
{
    int status=NC_NOERR;
    int myrank, nprocs, color;
    char path_sf[1024];
    MPI_Comm comm_sf;
    double ratio;

    MPI_Comm_rank (ncp->nciop->comm, &myrank);
    MPI_Comm_size (ncp->nciop->comm, &nprocs);
#ifdef SUBFILE_DEBUG
    if (myrank == 0)
      printf("%s: rank(%d): nprocs=%d, ncp->nc_num_subfiles=%d\n", __func__,
           myrank, nprocs, ncp->nc_num_subfiles);
#endif

    /* split the original comm to subcomm */
    if (nprocs > ncp->nc_num_subfiles) {
        ratio = (double)nprocs/(double)ncp->nc_num_subfiles;
	color = (int)((double)myrank/ratio);
    }
    else 
        color = myrank%ncp->nc_num_subfiles;

#ifdef SUBFILE_DEBUG
    if (myrank == 0)
      printf("%s: rank(%d): color=%d\n", __func__, myrank, color);
#endif
    //key = myrank/comm_size;
    
    /* TODO: fix error when using generated key value.
     * for now, just passing 0 value. */
    MPI_Comm_split(ncp->nciop->comm, color, myrank, &comm_sf);

    //char path[1024], file[1024];
    //find_path_and_fname (ncp->nciop->path, path, file);
    sprintf(path_sf, "%s.subfile_%i.%s", ncp->nciop->path, color, "nc");
    //sprintf(path_sf, "%s%d/%s", path, color, file);
    
    status = ncmpi_open(comm_sf, path_sf, NC_CLOBBER|NC_64BIT_DATA, MPI_INFO_NULL, ncidp);
    if (status != NC_NOERR) {
        fprintf(stderr, "error in %s: %s\n", __func__, ncmpi_strerror(status));
        return status;
    }

    return status;
}

int ncmpii_subfile_close (NC *ncp)
{
    int status = NC_NOERR;
    NC *ncp_sf;
    
    status = ncmpii_NC_check_id(ncp->ncid_sf, &ncp_sf);
    if (status != NC_NOERR)
        return status;
	
    status = ncmpii_NC_close (ncp_sf);
    if (status != NC_NOERR)
        return status;
    
    /* reset values to 0 */
    is_partitioned = 0;

    return status;
}

/*----< ncmpii_subfile_partition() >---------------------------------------*/
int ncmpii_subfile_partition (NC *ncp, int *ncidp)
{
    int i, color, myrank, nprocs, status=NC_NOERR;
    nc_type type;
    MPI_Offset attlen;
    NC *ncp_sf;
    double ratio;

    MPI_Comm_rank (ncp->nciop->comm, &myrank);
    MPI_Comm_size (ncp->nciop->comm, &nprocs);
#ifdef SUBFILE_DEBUG
    if (myrank==0)  /* debug */
    {
        printf("rank(%d): is_partitioned=%d\n", myrank, is_partitioned);
    }
#endif
    if (is_partitioned == 1) return NC_NOERR;

    if (nprocs > ncp->nc_num_subfiles) {
        ratio = (double)nprocs/(double)ncp->nc_num_subfiles;
	color = (int)((double)myrank/ratio);
    }
    else 
	color = myrank%ncp->nc_num_subfiles;

    if (myrank == 0) /* debug */
        printf("%s: rank(%d): color=%d\n", __func__, myrank, color);
#ifdef SUBFILE_DEBUG     
    printf("rank(%d): color=%d\n", myrank, color);
#endif

    /* check whether file is already partitioned */
    /* this is to handle when app has multiple ncmpi_enddef() */
    status = ncmpi_inq_att(ncp->nciop->fd, NC_GLOBAL, "num_subfiles", &type, &attlen);
    /* TODO: should check specific error like such attr doesn't exist? */
    if (status != NC_NOERR) {
	status = ncmpii_subfile_create(ncp, ncidp);
	TEST_HANDLE_ERR ("ncmpii_subfile_create", status);
	
	status = ncmpi_put_att_int (ncp->nciop->fd, NC_GLOBAL, 
				    "num_subfiles", NC_INT, 1, 
				    &ncp->nc_num_subfiles); 
	TEST_HANDLE_ERR ("ncmpi_put_att_int", status);
    }
    
    /* TODO: ignore UNLIMITED dims */
    /* NOTE: the following "for loop" should be before NC_begins() */
    /* NOTE: ncp->vars.nalloc shouldn't be used for loop bound
       because it will incremented by 4, not by 1 */
    status = ncmpii_NC_check_id(ncp->ncid_sf, &ncp_sf);
    TEST_HANDLE_ERR ("ncmpii_NC_check_id", status);
    
    /* adjust the hints to be used by PnetCDF; use the same value in master */
    ncp_sf->nciop->hints.header_align_size = ncp->nciop->hints.header_align_size;
    ncp_sf->nciop->hints.var_align_size    = ncp->nciop->hints.var_align_size;
    
    for(i=0; i<ncp->vars.ndefined; i++) { /* travere all variables */
	NC_var **vpp = ncp->vars.value;
	NC_dim **dpp = ncp->dims.value;
	
	/* check attr for subfiles */
	int par_dim_id = (IS_RECVAR(vpp[i])?1:0); /* default is the most significant dim excluding NC_UNLIMITED */

	/* skip if the var is partitioned already */
	if (vpp[i]->dimids == NULL) continue; 

        /* par_dim_id already exists? */
	status = ncmpi_inq_att(ncp->nciop->fd, i, "par_dim_id",
			       &type, &attlen);
	if (status == NC_NOERR) {
	    int par_dim_id_temp, i_temp;
	    status = ncmpi_get_att_int (ncp->nciop->fd, i, 
					"par_dim_id",
					&par_dim_id_temp);
#ifdef SUBFILE_DEBUG 
	    printf("par_dim_id_temp=%d\n", par_dim_id_temp);
#endif
	    /* convert par_dim_id to dim index defined in a var */
	    for (i_temp = 0; i_temp <vpp[i]->ndims ; i_temp++)
		if (vpp[i]->dimids[i_temp] == par_dim_id_temp) {
#ifdef SUBFILE_DEBUG
		    if (myrank == 0)
			printf("dimids[%d]=%d\n", i_temp, 
			       vpp[i]->dimids[i_temp]);
#endif
		    par_dim_id = i_temp;
		    break;
		}
#ifdef SUBFILE_DEBUG
	    printf("par_dim_id=%d\n", par_dim_id);
#endif
	}
	if (par_dim_id < vpp[i]->ndims)
	    ncmpi_put_att_int (ncp->nciop->fd, i, "par_dim_index", 
			       NC_INT, 1, &par_dim_id);
	
#ifdef SUBFILE_DEBUG
	if (myrank == 0) 
	    printf ("%s: var(%s): size of partitioning dim (id=%d)=%d\n", 
		    __func__, (*vpp[i]).name->cp, par_dim_id, 
		    dpp[par_dim_id]->size);
#endif
	/* divide only when dim is partitionable */
	/* 1. skip sizeof(par_dim_id) is smaller than num_subfiles */
	/* 2. skip if ndims < min_ndims */
	if ( (dpp[(*vpp[i]).dimids[par_dim_id]]->size)/(ncp->nc_num_subfiles) != 0 && (vpp[i]->ndims >= par_dim_id+1) && (vpp[i]->ndims >= min_ndims)) { 
	    int varid, j, jj, k;
	    int var_ndims = vpp[i]->ndims; /* keep org ndims */
	    int dimids[var_ndims];
	    char *key[ncp->nc_num_subfiles][var_ndims];
            
	    for (jj=0; jj < ncp->nc_num_subfiles; jj++) 
		for (k=0; k<var_ndims; k++)
		    key[jj][k] = calloc(100, sizeof(char));
	    
	    /* save original value ndims to attribute before making to 0 */
	    status = ncmpi_put_att_int (ncp->nciop->fd, i, "ndims_org", NC_INT, 1, &vpp[i]->ndims);
	    TEST_HANDLE_ERR ("ncmpi_put_att_int", status);
	    
	    int sf_range[ncp->nc_num_subfiles][var_ndims][3];
	    
            /* j: each dimension */
	    /* subfile: create a var with partitioned dim sizes */
	    for(j=0; j<var_ndims; j++) {
		MPI_Offset dim_sz;
		char str[80];
		double x; /* = org dim size / num_subfiles */
                int min, max;
                
		dim_sz = dpp[vpp[i]->dimids[j]]->size; /* init both to org dim sz */
		/* determine the most significant dimension */
		if (j == par_dim_id) {
		    x = (double)(dpp[vpp[i]->dimids[j]]->size)/(double)(ncp->nc_num_subfiles);
		}
                
		/* don't partition dim if dim size is less than ratio x */
		if ((int)x < 1 && j == par_dim_id)
		    continue;

                /* don't need to check? */
                if (j == par_dim_id)
                {
                    double xx = x*(double)color;
		    double yy = x*(double)(color+1);
                    min = xx+(color==0||(xx-(int)xx==0.0)?0:1);
                    max = yy-(yy-(int)yy==0.0?1:0);
		    if (max >= dpp[vpp[i]->dimids[j]]->size) max = dpp[vpp[i]->dimids[j]]->size-1;
		    dim_sz = max-min+1;
		}

		snprintf(str, strlen (dpp[vpp[i]->dimids[j]]->name->cp) +
			 strlen (vpp[i]->name->cp) + 2, "%s.%s", 
			 dpp[vpp[i]->dimids[j]]->name->cp, 
			 vpp[i]->name->cp);
#ifdef SUBFILE_DEBUG
		printf("rank(%d): new dim name = %s, dim_sz=%d\n", myrank, str, dim_sz);
#endif		    
		status = ncmpi_def_dim(ncp->ncid_sf, str, 
				       dim_sz, &dimids[j]);
		TEST_HANDLE_ERR ("ncmpi_def_dim", status);
		
		//dpp_sf[color][j] = ncp_sf->dims.value[j];
                
		for (jj=0; jj < ncp->nc_num_subfiles; jj++) {
		    char ss[80];
                    double xx, yy;
		    snprintf(key[jj][j], strlen(dpp[vpp[i]->dimids[j]]->name->cp) + ncmpii_itoa(jj, ss) + 17, "range(%s).subfile.%d", dpp[vpp[i]->dimids[j]]->name->cp, jj); /* dim name*/ 
                    xx = x*(double)jj;
                    min = xx+(jj==0||(xx-(int)xx==0.0)?0:1);
		    yy = x*(double)(jj+1);
                    max = yy-(yy-(int)yy==0.0?1:0);
		    if (max >= dpp[vpp[i]->dimids[j]]->size) max = dpp[vpp[i]->dimids[j]]->size-1;
#ifdef SUBFILE_DEBUG
                    if (myrank == 0) printf("subfile(%d): min=%d, max=%d\n", jj, min, max);
#endif
		    if (j == par_dim_id) { /* partitioning dims? */ 
			sf_range[jj][j][0] = min;
                        sf_range[jj][j][1] = max;
			sf_range[jj][j][2] = (max-min+1);;
		    } 
		    else { 
			sf_range[jj][j][0] = 0;
			sf_range[jj][j][2] = dpp[vpp[i]->dimids[j]]->size;
			sf_range[jj][j][1] = (dpp[vpp[i]->dimids[j]]->size!=0?(dpp[vpp[i]->dimids[j]]->size)-1:0);
		    }
		}
	    } /* for each dim */
	    
	    /* master file: replace the original var with scalar var */
	    vpp[i]->dimids_org = calloc(vpp[i]->ndims, sizeof(int));
	    memcpy(vpp[i]->dimids_org, vpp[i]->dimids, vpp[i]->ndims*sizeof(int));
	    vpp[i]->ndims_org = vpp[i]->ndims;
	    vpp[i]->ndims = 0;
	    vpp[i]->dimids = (IS_RECVAR(vpp[i])?&vpp[i]->dimids[0]:NULL);
	    vpp[i]->len = vpp[i]->xsz; /* size of type  */
	    vpp[i]->shape = NULL;
	    vpp[i]->num_subfiles = ncp->nc_num_subfiles;
            
	    status = ncmpi_put_att_int (ncp->nciop->fd, i, "num_subfiles", 
					NC_INT, 1, &ncp->nc_num_subfiles);
	    TEST_HANDLE_ERR ("ncmpi_put_att_int", status);
	    
	    for (jj=0; jj < ncp->nc_num_subfiles; jj++)
		for (k=0; k < var_ndims; k++) {
		    status = ncmpi_put_att_int (ncp->nciop->fd, i, 
						key[jj][k], NC_INT, 
						2, sf_range[jj][k]);
		    TEST_HANDLE_ERR ("ncmpi_put_att_int", status);
		}
	    
	    /* define a var with new dim */
	    status = ncmpi_def_var(ncp->ncid_sf, (*vpp[i]).name->cp, 
				   vpp[i]->type, var_ndims, dimids, &varid);
	    TEST_HANDLE_ERR ("ncmpi_def_var", status);
	    
	    /* add an attribute about each dim's range in subfile */
	    /* varid: var id in subfile */
	    for (k=0; k<var_ndims; k++) {
		status = ncmpi_put_att_int (ncp->ncid_sf, varid, key[color][k], 
					    NC_INT, 2, sf_range[color][k]); 
		TEST_HANDLE_ERR ("ncmpi_put_att_int", status);
	    }
	    
	    /* deallocate buffers */
	    for (jj=0; jj < ncp->nc_num_subfiles; jj++) 
		for (k=0; k<var_ndims; k++)
		    free(key[jj][k]);
	    
	    status = ncmpi_put_att_int (ncp->ncid_sf, varid, "subfile_index",
					NC_INT, 1, &color);
	    TEST_HANDLE_ERR ("ncmpi_put_att_int", status);
	} /* end if() */
    } /* for each variable */

    is_partitioned = 1;

    return NC_NOERR;
}

/*----< ncmpii_subfile_getput_vars() >---------------------------------------*/
int
ncmpii_subfile_getput_vars(NC               *ncp,
			   NC_var           *varp,
			   const MPI_Offset  start[],
			   const MPI_Offset  count[],
			   const MPI_Offset  stride[],
			   void             *buf,
			   MPI_Offset        bufcount,
			   MPI_Datatype      buftype, /* data type of the buffer */
			   int               rw_flag,
			   int               io_method)
{
    int err, errs=0, status;
    NC_subfile_access *my_req, *others_req;
    int i, j, k, myrank, nprocs;
    NC *ncp_sf;
    int *count_my_req_per_proc, *count_others_req_per_proc;
    int varid, varid_sf;
    MPI_Offset *buf_count_my, *buf_count_others;
    void **xbuf=NULL, *cbuf=NULL;
    MPI_Request *requests=NULL;
    MPI_Status *statuses=NULL;
    int ndims_org = varp->ndims_org;
    int color;
    int nasyncios=0;

    /* "API error" will abort this API call, but not the entire program */
    err = status = NC_NOERR;

    MPI_Comm_rank(ncp->nciop->comm, &myrank);
    MPI_Comm_size(ncp->nciop->comm, &nprocs);

#ifdef SUBFILE_DEBUG
    for (i=0; i<ndims_org; i++)
        printf("rank(%d): %s: var(%s): start[%d]=%ld, count[%d]=%ld, stride[%d]=%ld, bufcount=%ld\n",
               myrank, __func__, varp->name->cp, i, start[i], i, count[i], i, 
               ((stride != NULL)?stride[i]:1), bufcount);
#endif

    /* get ncp info for the subfile */
    status = ncmpii_NC_check_id(ncp->ncid_sf, &ncp_sf);
    TEST_HANDLE_ERR("ncmpii_NC_check_id", status);

    /* check attr for subfiles */
    nc_type type;
    MPI_Offset attlen;
    int par_dim_id = 0; /* default is the most significant dim */

    status = ncmpi_inq_varid(ncp->nciop->fd, varp->name->cp, &varid);
    TEST_HANDLE_ERR("ncmpi_inq_varid", status);

    status = ncmpi_inq_att(ncp->nciop->fd, varid, "par_dim_index",
				   &type, &attlen);
    if (status == NC_NOERR) {
	status = ncmpi_get_att_int (ncp->nciop->fd, varid,
				    "par_dim_index",
				    &par_dim_id);
    }
#ifdef SUBFILE_DEBUG
    if (myrank == 0)
        printf("par_dim_index of %s = %d\n", varp->name->cp, par_dim_id);
#endif
    status = ncmpi_inq_varid(ncp_sf->nciop->fd, varp->name->cp, &varid_sf);
    TEST_HANDLE_ERR("ncmpi_inq_varid", status);

    status = ncmpi_get_att_int (ncp_sf->nciop->fd, varid_sf, "subfile_index",
                                &color);

    count_my_req_per_proc = (int *)calloc(nprocs, sizeof(int));
    count_others_req_per_proc = (int *)calloc(nprocs, sizeof(int));
    buf_count_my = (MPI_Offset *)calloc(nprocs, sizeof(MPI_Offset));
    buf_count_others = (MPI_Offset *)calloc(nprocs, sizeof(MPI_Offset));

    /* TODO: shouldn't it be initialized to 0? */
    for (i=0; i<nprocs; i++) {
        buf_count_my[i] = 1;
        buf_count_others[i] = 1;
    }

    /* allocate space for my_req */
    my_req = (NC_subfile_access *)calloc (nprocs, sizeof(NC_subfile_access));

    /* init allocated my_req */
    for (i=0; i<nprocs; i++) {
        my_req[i].start = (MPI_Offset *)calloc(ndims_org, sizeof(MPI_Offset));
        my_req[i].count = (MPI_Offset *)calloc(ndims_org, sizeof(MPI_Offset));
        my_req[i].start_org = (MPI_Offset *)calloc(ndims_org, sizeof(MPI_Offset));
        /* init start/count to -1 */
        for (j=0; j<ndims_org; j++) {
            my_req[i].start[j] = -1;
            my_req[i].count[j] = -1;
            my_req[i].start_org[j] = -1;
        }
    }

    /* i: for each subfile */
    for (i=0; i<ncp->nc_num_subfiles; i++) {
        int flag = 0; /* set to 1 if par_dim_id is partitioned, initially 0 */
        double ratio = (double)nprocs/(double)ncp->nc_num_subfiles;
        int aproc=-1; /* I/O delegate proc in subfile group */

        if (delegate_scheme == BALANCED) {
            double scaled, xx, yy;
            int min, max;
            
            xx = (ratio)*(double)i; /* i: each subfile */
            min = xx+(i==0||(xx-(int)xx==0.0)?0:1);
	    yy = (ratio)*(double)(i+1);
            max = yy-(yy-(int)yy==0.0?1:0);
	    if (max >= nprocs) max = nprocs-1;
            //scaled = (double)rand()/RAND_MAX;
            scaled = (double)myrank/ratio-(double)color;
            aproc = (i==color)?myrank:(min+(max-min+1)*scaled);
        }
        else if (delegate_scheme == ONE)
        {
            double xx;
            int min;
            xx = (ratio)*(double)i;
            min = xx+(i==0||(xx-(int)xx==0.0)?0:1);
            aproc = (i==color)?myrank:min;
        }

        /* check out of range? */
        if (aproc >= nprocs)
            aproc = nprocs-1;

#ifdef SUBFILE_DEBUG
        printf("rank(%d): color=%d, subfile=%d, aproc=%d\n", myrank, color, i, aproc);
#endif

	/* j: for each dim starting from par_dim_id in round-robin manner */
        for (j=par_dim_id; j<par_dim_id+ndims_org; j++) {
	    int jx = j%ndims_org;
            NC_dim *dimp = ncp_sf->dims.value[ncp_sf->vars.value[varid_sf]->dimids[jx]];
            int sf_range[2];
            int ii, jj, kk, stride_count;
	    char key[80], *org_dim_name;

            org_dim_name = strtok(dimp->name->cp, ".");
	    sprintf(key, "range(%s).subfile.%d", org_dim_name, i); /* dim name*/
#ifdef SUBFILE_DEBUG
            if (myrank == 0) {
                printf("rank(%d): org_dim_name=%s\n", myrank, org_dim_name);
                printf("rank(%d): key=%s\n", myrank, key);
            }
#endif
            /* should inquire master file */
            status = ncmpi_inq_att (ncp->nciop->fd, varid, key, &type, &attlen);
            TEST_HANDLE_ERR("ncmpi_inq_att", status);
            
            /* should get range info from the master file */
	    status = ncmpi_get_att_int (ncp->nciop->fd, varid, key, sf_range);
	    TEST_HANDLE_ERR("ncmpi_get_att_int", status);

#ifdef SUBFILE_DEBUG
            //if (myrank == 0) 
                printf("rank(%d): ndims_org=%d, rstart=%d, rend=%d, start[%d]=%d, count[%d]=%d\n", 
                       myrank, ndims_org, sf_range[0], sf_range[1], jx, start[jx], jx, count[jx]);
#endif
            /* ii: traverse user's request range, incremented by stride count
               jj: traverse subfile range, incrementing sequentially
               kk: count belong to my subfile range */
            ii=start[jx], jj=sf_range[0], kk=0;
            stride_count = (stride == NULL?1:stride[jx]);
            //printf("stride_count[%d]=%d\n", jx, stride_count);
	    
	    /* TODO: if sf_range is 1, count[] value is incorrect 
	       e.g., size of par_dim is 4 and the nproc is also 4 */
	    int is_unlimited_dim = sf_range[1]-sf_range[0];

            /* skip the remaining dim if par_dim is not 
               intersected with user's request */
            if (jx != par_dim_id && !flag) continue;

            while (ii < (start[jx]+(count[jx]*stride_count)) && 
		   jj <= (is_unlimited_dim==0?(count[jx]*stride_count-1):sf_range[1]))
            {
                if (ii < jj)
                    ii+=stride_count;
                else if (jj < ii)
                    jj++;
                else {
#ifdef SUBFILE_DEBUG
		    printf("rank(%d): var(%s): i=%d, j=%d, ii=%d, jj=%d, kk=%d, jx=%d\n", myrank, varp->name->cp, i, j, ii, jj, kk, jx);
#endif
                    if (kk == 0) {
			my_req[aproc].start[jx] = ii;
#ifdef SUBFILE_DEBUG
			printf("rank(%d): var(%s): my_req[%d].start[%d]=%d\n",
			       myrank, varp->name->cp, aproc, jx, ii);
#endif
		    }
                    if (jx == par_dim_id) flag = 1;
                    ii+=stride_count; jj++; kk++;
                }

            }
            if (kk > 0 && flag == 1) my_req[aproc].count[jx] = kk;
            /* adjust start offset based on subfile's range start.
               otherwise, there will be an out of bound error during I/O */
            if (my_req[aproc].start[jx] != -1) {
                my_req[aproc].start_org[jx] = my_req[aproc].start[jx];
                my_req[aproc].start[jx] -= sf_range[0];
            }
#ifdef SUBFILE_DEBUG
            //if (myrank == 0) {
	    {
                printf("rank(%d): my_req[%d].start[%d]=%d\n", myrank, aproc, 
                       jx, my_req[aproc].start[jx]);
                printf("rank(%d): my_req[%d].count[%d]=%d\n", myrank, aproc,
                       jx, my_req[aproc].count[jx]);
            }
#endif          
        } /* for each dim, j */
        if (my_req[aproc].start[0] == -1)
            count_my_req_per_proc[aproc] = 0;
        else
            count_my_req_per_proc[aproc] = 1;
    } /* for each subfile, i */
    
#ifdef SUBFILE_DEBUG
    for (i=0; i<ncp->nc_num_subfiles; i++) {
        char str_st[100], str_st_org[100], str_ct[100], str_t1[10];
        sprintf(str_st, ">> rank(%d): subfile(%d): var(%s): start(", myrank, i, varp->name->cp); 
        sprintf(str_ct, ">> rank(%d): subfile(%d): count(", myrank, i);
        sprintf(str_st_org, "%d: start_org(", i);
        for (j=0; j<ndims_org; j++) {
            sprintf(str_t1, "%d", my_req[i].start[j]);
            strcat(str_st, str_t1);
            sprintf(str_t1, "%d", my_req[i].count[j]);
            strcat(str_ct, str_t1);
            sprintf(str_t1, "%d", my_req[i].start_org[j]);
            strcat(str_st_org, str_t1);
            if(j != (ndims_org-1)) {
                strcat(str_st, ",");
                strcat(str_ct, ",");
                strcat(str_st_org, ",");
            }
            
        }
        strcat(str_st, ")");
        strcat(str_ct, ")");
        strcat(str_st_org, ")");
        printf("%s\n", str_st);
        printf("%s\n", str_ct);
        printf("%s\n", str_st_org);
    }
#endif

    /* calculate buf_count based on subfile */
    /* TODO: should do only when count_my_req_per_proc[myrank] != 0??? */
    for (i=0; i<ndims_org; i++)
        buf_count_my[myrank] *= my_req[myrank].count[i];

    /* build other proc's request */
    others_req = (NC_subfile_access *)calloc (nprocs, sizeof(NC_subfile_access));
    /* init allocated others_req */
    for (i=0; i<nprocs; i++) {
        others_req[i].start = (MPI_Offset *)calloc(ndims_org, sizeof(MPI_Offset));
        others_req[i].count = (MPI_Offset *)calloc(ndims_org, sizeof(MPI_Offset));
        others_req[i].start_org = (MPI_Offset *)calloc(ndims_org, sizeof(MPI_Offset));
        /* init start/count to -1 */
        for (j=0; j<ndims_org; j++) {
            others_req[i].start[j] = -1;
            others_req[i].count[j] = -1;
            others_req[i].start_org[j] = -1;
        }
    }
    
#ifdef SUBFILE_DEBUG
    for (i=0; i<nprocs && myrank==0; i++) {
        printf("%d: count_my_req_per_proc[%d]=%d\n", myrank, i, count_my_req_per_proc[i]);
        printf("%d: count_others_req_per_proc[%d]=%d\n", myrank, i, count_others_req_per_proc[i]);
    }
#endif

#ifdef TAU_SSON
    TAU_PHASE_CREATE_STATIC(t51, "SSON --- getput_vars MPI_Alltoall", "", TAU_USER);
    TAU_PHASE_START(t51);
#endif

    MPI_Alltoall (count_my_req_per_proc, 1, MPI_INT, count_others_req_per_proc, 1, MPI_INT, ncp->nciop->comm);

#ifdef TAU_SSON
    TAU_PHASE_STOP(t51);
#endif

#ifdef SUBFILE_DEBUG
    for (i=0; i<nprocs && myrank==0; i++) {
        printf("=> %d: count_others_req_per_proc[%d]=%d\n", myrank, i, count_others_req_per_proc[i]);
    }
#endif

    MPI_Offset buf_offset[nprocs];
    for (i=0; i<nprocs; i++)
        buf_offset[i] = 0;

#ifdef SUBFILE_DEBUG
    if (myrank == 0) printf("varname=%s: etype=0x%x, etype_size=%d\n", varp->name->cp, varp->type, varp->xsz);
#endif

    /* find the ptype (primitive MPI data type) from buftype
     * el_size is the element size of ptype
     * bnelems is the total number of ptype elements in the I/O buffer, buf
     * fnelems is the number of nc variable elements in nc_type
     * nbytes is the amount of read/write in bytes
     */
    MPI_Datatype ptype;
    int el_size;
    int isderived, buftype_is_contig;
    MPI_Offset bnelems;
    MPI_Aint lb, extent;

    status = ncmpii_dtype_decode(buftype, &ptype, &el_size, &bnelems,
				 &isderived, &buftype_is_contig);
    /* bnelems now is the number of ptype in a buftype */
    TEST_HANDLE_ERR("ncmpii_dtype_decode", status);

    MPI_Type_get_extent(buftype, &lb, &extent);

#ifdef SUBFILE_DEBUG
    printf("rank(%d): var(%s): ptype=0x%x, el_size=%d, bnelems=%d, isderived=%d, buftype_is_contig=%d, lb=%d, extent=%d\n",
	   myrank, varp->name->cp, ptype, el_size, bnelems, isderived, buftype_is_contig, lb, extent);
#endif

    /* if buftype is non-contiguous, pack to contiguous buffer*/
    /* NOTE: no conversion and byte swap are performed here
       as they are done underneath layer */
    if (!buftype_is_contig) {
	bnelems = bnelems*bufcount; /* update bnelems */
	cbuf = (void *)calloc(bnelems, el_size);
	if (cbuf == NULL)
	    printf("calloc of cbuf error!!!!\n");
	if (rw_flag == WRITE_REQ) {
	    status = ncmpii_data_repack((void*)buf, bufcount, buftype,
					cbuf, bnelems, ptype);
	    TEST_HANDLE_ERR("ncmpii_data_repack", status);
	}
    }
    else
	cbuf = (void *)buf;

    int diff[ndims_org];
    for (i=0; i < ndims_org; i++) {
        int stride_count;
        stride_count = (stride == NULL?1:stride[i]);
        /* making diff is necessary?? */
        diff[i] = ABS(my_req[myrank].start_org[i]-start[i])/stride_count;
#ifdef SUBFILE_DEBUG
        if (myrank == 0) printf("rank(%d): my_req[%d].start_org[%d]=%d, start[%d]=%d, diff[%d]=%d\n", myrank, 
               myrank, i, my_req[myrank].start_org[i], i, start[i], i, diff[i]);
#endif
    }

    for (i=0; i<ndims_org; i++) {
        int l, tmp=1;
        for (l=i+1; l < ndims_org ; l++) {
            tmp *= my_req[myrank].count[l];
        }
        buf_offset[myrank] += tmp*diff[i];
#ifdef SUBFILE_DEBUG
        if (myrank == 0) printf("local: rank(%d): buf_offset[%d]=%d\n", myrank, myrank, buf_offset[myrank]);
#endif
    }

    buf_offset[myrank] *= el_size;

#ifdef SUBFILE_DEBUG
    printf("rank(%d): buf_offset[%d]=%d, buf_count_my[%d]=%d\n", myrank, 
	   myrank, buf_offset[myrank], myrank, buf_count_my[myrank]);

    printf ("rank(%d): %s: var(%s)\n", myrank, __func__, varp->name->cp);
    for (i=0; i<ndims_org; i++)
        printf("my_req[%d].start[%d]=%ld, my_req[%d].count[%d]=%ld, my_req[%d].stride[%d]=%ld, bufcount=%ld\n",
               myrank, i, my_req[myrank].start[i], myrank, i, 
	       my_req[myrank].count[i], myrank, i, 
	       ((stride != NULL)?stride[i]:1), bufcount);
#endif

    int *array_of_requests;
    int *array_of_statuses;
    /* TODO: each proc can't get more than nprocs-1?? */
    array_of_requests = (int *)calloc(nprocs, sizeof(int));
    for (i=0; i<nprocs; i++) 
	array_of_requests[i] = NC_REQ_NULL; 

#ifdef SUBFILE_DEBUG
    printf("buf_count_my[%d]=%d\n", myrank, buf_count_my[myrank]);
#endif
    /* doing my portion of I/O */
    if (count_my_req_per_proc[myrank] != 0) {
        status = ncmpii_igetput_varm(ncp_sf, 
				     ncp_sf->vars.value[varid_sf], 
				     my_req[myrank].start, 
				     my_req[myrank].count, 
                                     stride, 
                                     NULL,
				     cbuf+buf_offset[myrank], 
				     buf_count_my[myrank],
				     (!buftype_is_contig?ptype:buftype),
				     &array_of_requests[nasyncios++],
				     rw_flag, 
				     0);
	TEST_HANDLE_ERR("ncmpii_igetput_varm1", status);
    }
#ifdef SUBFILE_DEBUG
    printf("rank(%d): var(%s): pushed local I/O to async calls...\n", myrank, varp->name->cp);
#endif

    ////////// doing other proc's request to my subfile
    /* TODO: each proc can't get more than nprocs?? */
    requests = (MPI_Request *)malloc(2*nprocs*sizeof(MPI_Request));
    
    j = 0;
#ifdef TAU_SSON
    TAU_PHASE_CREATE_STATIC(t52, "SSON --- getput_vars: exch count_my_req", "", TAU_USER);
    TAU_PHASE_START(t52);
#endif
    /* TODO: need to exchange count_my_req_per_proc? */
    for (i=0; i<nprocs; i++)
        if (count_my_req_per_proc[i] != 0 && i != myrank)
            MPI_Irecv(&count_others_req_per_proc[i], 1, MPI_INT, i, i+myrank, ncp->nciop->comm, &requests[j++]);

    for (i=0; i<nprocs; i++)
        if (count_others_req_per_proc[i] != 0 && i != myrank)
            MPI_Isend(&count_my_req_per_proc[i], 1, MPI_INT, i, i+myrank, ncp->nciop->comm, &requests[j++]);

    statuses = (MPI_Status *)malloc(j*sizeof(MPI_Status));
    err = MPI_Waitall(j, requests, statuses);
#ifdef TAU_SSON
    TAU_PHASE_STOP(t52);
#endif
    check_err(MPI_Waitall);
    free(statuses);
    free(requests);

    j = 0;
    requests = (MPI_Request *)malloc(2*nprocs*sizeof(MPI_Request));

#ifdef TAU_SSON
    TAU_PHASE_CREATE_STATIC(t53, "SSON --- getput_vars: exch start,count,", "", TAU_USER);
    TAU_PHASE_START(t53);
#endif
    for(i=0; i<nprocs; i++) {
        if(count_others_req_per_proc[i] != 0 && i != myrank) {
            /* MPI_Offset == MPI_LONG_LONG_INT */
            MPI_Irecv(others_req[i].start, ndims_org, MPI_LONG_LONG_INT, i, i+myrank, 
                      ncp->nciop->comm, &requests[j++]);
            MPI_Irecv(others_req[i].count, ndims_org, MPI_LONG_LONG_INT, i, i+myrank, 
                      ncp->nciop->comm, &requests[j++]);
            MPI_Irecv(others_req[i].start_org, ndims_org, MPI_LONG_LONG_INT, i, 
                      i+myrank, ncp->nciop->comm, &requests[j++]);
        }
    }

    for(i=0; i<nprocs; i++) {
        if(count_my_req_per_proc[i] != 0 && i != myrank) {
            MPI_Isend(my_req[i].start, ndims_org, MPI_LONG_LONG_INT, i, i+myrank, 
                      ncp->nciop->comm, &requests[j++]);
            MPI_Isend(my_req[i].count, ndims_org, MPI_LONG_LONG_INT, i, i+myrank, 
                      ncp->nciop->comm, &requests[j++]);
            MPI_Isend(my_req[i].start_org, ndims_org, MPI_LONG_LONG_INT, i, i+myrank, 
                      ncp->nciop->comm, &requests[j++]);
        }
    }
    
    statuses = (MPI_Status *)malloc(j*sizeof(MPI_Status));
    err = MPI_Waitall(j, requests, statuses);
    free(statuses);
    free(requests);

#ifdef TAU_SSON
    TAU_PHASE_STOP(t53);
#endif
    check_err(MPI_Waitall);

#ifdef SUBFILE_DEBUG
    /* DEBUG: print out others_req.{start,count} */
    for (i=0; i<nprocs && myrank == 0; i++) {
        char str_st[100], str_st_org[100], str_ct[100], str_t1[10];
        sprintf(str_st, "%d: others.start(", i); 
        sprintf(str_ct, "%d: others.count(", i);
        sprintf(str_st_org, "%d: others.start_org(", i);
        
        for (j=0; j<ndims_org; j++) {
            sprintf(str_t1, "%d", others_req[i].start[j]);
            strcat(str_st, str_t1);
            sprintf(str_t1, "%d", others_req[i].count[j]);
            strcat(str_ct, str_t1);
            sprintf(str_t1, "%d", others_req[i].start_org[j]);
            strcat(str_st_org, str_t1);
            if(j != (ndims_org-1)) {
                strcat(str_st, ",");
                strcat(str_ct, ",");
                strcat(str_st_org, ",");
            }
                
        }
        strcat(str_st, ")");
        strcat(str_ct, ")");
        strcat(str_st_org, ")");
        printf("%s\n", str_st);
        printf("%s\n", str_ct);
        printf("%s\n", str_st_org);
    }
#endif

    j = 0;
    requests = (MPI_Request *)malloc(2*nprocs*sizeof(MPI_Request));
    
    xbuf = (void *)calloc(nprocs, sizeof(void *));

#ifdef TAU_SSON
    TAU_PHASE_CREATE_STATIC(t54, "SSON --- getput_vars: exch buf", "", TAU_USER);
    TAU_PHASE_START(t54);
#endif

    for (i=0; i<nprocs; i++) {
	//xbuf[i] = NULL;
        buf_count_others[i] = 1;
        if (count_others_req_per_proc[i] != 0 && i != myrank) {
            for (k=0; k<ndims_org; k++) {
                buf_count_others[i] *= others_req[i].count[k];
            }
#ifdef SUBFILE_DEBUG
            printf("rank(%d): recv from rank %d: buf_count_others[%d]=%d\n", myrank, i, i, buf_count_others[i]);
#endif
	    xbuf[i] = (void*)calloc(buf_count_others[i], el_size);
            if (xbuf[i]==NULL) {
                printf("Error allocating memory!\n"); //print an error message
                return 1; //return with failure
            }
            MPI_Irecv(xbuf[i], buf_count_others[i], (!buftype_is_contig?ptype:buftype), i, i+myrank, ncp->nciop->comm, &requests[j++]);
        }
    }

    for (i=0; i<nprocs; i++) {
        buf_count_my[i] = 1;
        if (count_my_req_per_proc[i] != 0 && i != myrank) {
            int diff[ndims_org];
            for (k=0; k < ndims_org; k++) {
                int stride_count, tmp=1, l;

                stride_count = (stride == NULL?1:stride[k]);
                diff[k] = ABS(my_req[i].start_org[k]-start[k])/stride_count;

                for (l=k+1; l < ndims_org; l++) 
                    tmp *= my_req[i].count[l];

                buf_offset[i] += tmp*diff[k];
#ifdef SUBFILE_DEBUG
                if (myrank == 0) printf("rank(%d): subfile(%d): diff[%d]=%d\n", myrank, i, k, diff[k]);
#endif
            }

            buf_offset[i] *= el_size;

            for (k=0; k<ndims_org; k++) {
                buf_count_my[i] *= my_req[i].count[k];
            }
#ifdef SUBFILE_DEBUG
            printf("rank(%d): send to rank %d: buf_offset[%d]=%d, buf_count_my[%d]=%d\n", myrank, i, i, buf_offset[i], i, buf_count_my[i]);
#endif

	    MPI_Isend(cbuf+buf_offset[i], buf_count_my[i], (!buftype_is_contig?ptype:buftype), i, i+myrank, ncp->nciop->comm, &requests[j++]);
        } /* end if() */
    } /* end for() */     

    statuses = (MPI_Status *)calloc(j, sizeof(MPI_Status));
    err = MPI_Waitall(j, requests, statuses);
    free(statuses);
    free(requests);

#ifdef TAU_SSON
    TAU_PHASE_STOP(t54);
#endif
    check_err(MPI_Waitall);

#ifdef TAU_SSON
    TAU_PHASE_CREATE_STATIC(t55, "SSON --- getput_vars igetput", "", TAU_USER);
    TAU_PHASE_START(t55);
#endif

    for (i=0; i<nprocs; i++) {
        if(count_others_req_per_proc[i] != 0 && i != myrank) {
            status = ncmpii_igetput_varm(ncp_sf, 
                                         ncp_sf->vars.value[varid_sf],
                                         others_req[i].start, 
                                         others_req[i].count,
                                         stride, 
                                         NULL,
                                         xbuf[i], 
					 buf_count_others[i],
					 (!buftype_is_contig?ptype:buftype),
                                         &array_of_requests[nasyncios++], 
                                         rw_flag, 
                                         0);
            TEST_HANDLE_ERR ("ncmpii_igetput_varm2", status);
        }
    }
#ifdef TAU_SSON
    TAU_PHASE_STOP(t55);
#endif

#ifdef SUBFILE_DEBUG
    printf("rank(%d): nasyncios=%d\n", myrank, nasyncios);
#endif

#ifdef TAU_SSON
    TAU_PHASE_CREATE_STATIC(t56, "SSON --- getput_vars ncmpi_wait_all", "", TAU_USER);
    TAU_PHASE_START(t56);
#endif
    //double stime, wait_time;
    //stime = MPI_Wtime();
    array_of_statuses = (int *)calloc(nasyncios, sizeof(int));
    status = ncmpii_wait(ncp_sf, COLL_IO, nasyncios, array_of_requests, array_of_statuses);
    TEST_HANDLE_ERR("ncmpii_wait", status);
    free(array_of_statuses);
    free(array_of_requests);

    //int pending_nreqs;
    //status = ncmpi_inq_nreqs(ncp->nciop->fd, &pending_nreqs);
    //printf("myrank(%d): pending_nreqs=%d\n", myrank, pending_nreqs);
    //wait_time = MPI_Wtime() - stime;
    //printf("myrank(%d): ncmpii_wait time = %f\n", myrank, wait_time);
#ifdef TAU_SSON
    TAU_PHASE_STOP(t56);
#endif

#ifdef SUBFILE_DEBUG
    printf("rank(%d): var(%s): after ncmpi_wait_all()\n", myrank, varp->name->cp);
#endif

    //MPI_Barrier(ncp->nciop->comm);

    /* free all allocated memories */

    for(i=0; i<nprocs; i++)
	if (xbuf[i] != NULL) free(xbuf[i]);
    free(xbuf);

    if (cbuf != NULL && cbuf != buf) {
	free(cbuf);
    }

    free(count_my_req_per_proc);
    free(count_others_req_per_proc);
    free(buf_count_my);
    free(buf_count_others);

    for (i=0; i<nprocs; i++) {
        free(my_req[i].start);
        free(my_req[i].count);
        free(my_req[i].start_org);
        free(others_req[i].start);
        free(others_req[i].count);
        free(others_req[i].start_org);
    }
    free(my_req);
    free(others_req);

    return status;
}
