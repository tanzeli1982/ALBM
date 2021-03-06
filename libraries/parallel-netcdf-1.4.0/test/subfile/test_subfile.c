/*  
 *  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: test_subfile.c 1479 2013-10-27 06:36:18Z wkliao $ */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <pnetcdf.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <dirent.h>

/* Prototype for functions used only in this file */
static void handle_error(int status);

static void handle_error(int status) {
  fprintf(stderr, "%s\n", ncmpi_strerror(status));
}

#define MAXLINE 128

/* The file name is taken as a command-line argument. */

/* Measures the I/O bandwidth for writing/reading a 3D
   block-distributed array to a file corresponding to the global array
   in row-major (C) order.
   Note that the file access pattern is noncontiguous.
  
   Array size 128^3. For other array sizes, change array_of_gsizes below.*/
#define TEST_HANDLE_ERR(status)                         \
    {                                                   \
        if ((status) != NC_NOERR) {                     \
            printf("Error at line %d (%s)\n", __LINE__, \
                   ncmpi_strerror((status)) );          \
            nerr++;                                     \
        }                                               \
    }

int main(int argc, char **argv)
{
    int opt;
    extern char *optarg;
    extern int optind;
    int i, j, array_of_gsizes[3],array_of_distribs[3];
    int order, nprocs, len, **buf, mynod;
    MPI_Offset bufcount;
    int array_of_dargs[3], array_of_psizes[3];
    int status;
    MPI_Offset sizes[3], array_of_starts[3];
    char *basename = NULL, *basename1 = NULL, filename[100];
    char dimname[20], varname[20];
    int ncid, dimids0[3], rank_dim[3], *varid;
    MPI_Info info;
    MPI_Offset **starts_list, **count_list;
    MPI_Offset *bufcount_list;
    int ndims=3, nvars=1, ngatts, unlimdimid;
    MPI_Datatype *datatype_list;
    int length = 128; /* 8MB per proc */
    double stim, write_tim, new_write_tim, write_bw;
    double read_tim, new_read_tim, read_bw;
    double open_tim, new_open_tim;
    double close_tim, new_close_tim;
    int num_sf = 2;
    int par_dim_id = 0; /* default is 0 */
    int do_read = 0;
    int nerr=0, sum_nerr;
    MPI_Info info_used;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mynod);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    /* process 0 takes the file name as a command-line argument and
       broadcasts it to other processes */
    if (!mynod) {
	while ((opt = getopt(argc, argv, "f:s:rp:n:l:")) != EOF) {
	    switch (opt) {
	    case 'f': basename = optarg;
		break;
	    case 's': num_sf = atoi(optarg);
		break;
	    case 'r': do_read = 1;
		break;
            case 'p': par_dim_id = atoi(optarg);
                break;
            case 'n': nvars = atoi(optarg);
                break;
            case 'l': length = atoi(optarg);
                break;
	    default:
		break;
	    }
	}
	if (basename == NULL) {
	    fprintf(stderr, "\n*#  Usage: test_subfile -f pathname -s num_sf -p par_dim_id \n\n");
	    MPI_Abort(MPI_COMM_WORLD, 1);
	}
	
	basename1 = (char *) malloc (MAXLINE);
	sprintf(basename1, "%s", basename);
	len = strlen(basename1);
	MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(basename, len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(&num_sf, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&par_dim_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&nvars, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&do_read, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else {
	basename1 = (char *) malloc (MAXLINE);
	MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(basename1, len+1, MPI_CHAR, 0, MPI_COMM_WORLD);
        MPI_Bcast(&num_sf, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&par_dim_id, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&nvars, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&do_read, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&length, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    array_of_gsizes[0] = array_of_gsizes[1] = array_of_gsizes[2] = length;

    order = MPI_ORDER_C;
    buf = (int **)malloc(nvars*sizeof(int*));
    if (buf == NULL){
        printf("buf malloc error\n");
        return 0;               
    }   
    bufcount_list = (MPI_Offset *)malloc(nvars*sizeof(MPI_Offset));
    if (bufcount_list == NULL){
        printf("bufcount_list malloc error\n");
        return 0;               
    }   
    starts_list = (MPI_Offset **)malloc(nvars*sizeof(MPI_Offset *));
    if (starts_list== NULL){
        printf("starts_list malloc error\n");
        return 0;               
    }   
    count_list = (MPI_Offset **)malloc(nvars*sizeof(MPI_Offset *));
    if (count_list == NULL){
        printf("count_list malloc error\n");
        return 0;               
    }   
    datatype_list = (MPI_Datatype*)malloc(nvars*sizeof(MPI_Datatype));
    if (datatype_list == NULL){
        printf("count_list malloc error\n");
        return 0;               
    } 

    for (i=0; i<nvars; i++) {
	starts_list[i] = (MPI_Offset *)malloc(ndims*sizeof(MPI_Offset));
	if (starts_list[i] == NULL){
	    printf("starts_list[%d] malloc error\n", i);
	    return 0;		
	}	
	count_list[i] = (MPI_Offset *)malloc(ndims*sizeof(MPI_Offset));
	if (count_list[i] == NULL){
	    printf("count_list[%d] malloc error\n", i);
	    return 0;		
	}
    }
  
    array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
    array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
    array_of_distribs[2] = MPI_DISTRIBUTE_BLOCK;

    array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
    array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
    array_of_dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;
    
    bufcount = 1;
    for (i=0; i<ndims; i++) {
        array_of_psizes[i] = 0;
        sizes[i] = length;
        bufcount *= length;
    }
    MPI_Dims_create(nprocs, ndims, array_of_psizes);
    for(i=0; i<ndims&&mynod==0; i++)
	printf("array_of_psizes[%d]=%d\n", i, array_of_psizes[i]);

    /* subarray in each process is len x len x len */
    for (i=0; i<ndims; i++)
        array_of_gsizes[i] = length * array_of_psizes[i];

    /* mynd's process rank in each dimension (in MPI_ORDER_C) */
    rank_dim[2] =  mynod %  array_of_psizes[2];
    rank_dim[1] = (mynod /  array_of_psizes[2]) % array_of_psizes[1];
    rank_dim[0] =  mynod / (array_of_psizes[2]  * array_of_psizes[1]);

    /* starting coordinates of the subarray in each dimension */
    for (i=0; i<ndims; i++)
        array_of_starts[i] = length * rank_dim[i];

    for (i=0; i<nvars; i++) {
	for (j=0; j<ndims; j++) {
	    starts_list[i][j] = array_of_starts[j];
	    count_list[i][j]  = length;
	}
        bufcount_list[i] = bufcount;
        datatype_list[i] = MPI_INT;
    }

    for (i=0; i<nvars; i++) {
	buf[i] = (int *) malloc(bufcount * sizeof(int));
	if (buf[i] == NULL){
	    printf("buf[i]malloc error\n");
	    return 0;		
	}	

	for (j=0; j<bufcount; j++)
	    buf[i][j]=mynod+1;
    }

    MPI_Info_create(&info);
    //set all non-record variable to be subfiled
    char tmp[10];
    sprintf(tmp, "%d", num_sf);
    MPI_Info_set(info, "nc_num_subfiles", tmp);

    sprintf(filename, "%s.%d.%d.%d.nc", basename1, length, 1, 0);

    if (do_read == 1) goto read;

    stim = MPI_Wtime();
    status = ncmpi_create(MPI_COMM_WORLD, filename, NC_CLOBBER|NC_64BIT_DATA,
                          info, &ncid);
    TEST_HANDLE_ERR(status);

    open_tim = MPI_Wtime() - stim;
    
    MPI_Allreduce(&open_tim, &new_open_tim, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    if (mynod == 0) {
        fprintf(stderr, "create time = %f sec\n", new_open_tim);
    }

    /* define dimensions */
    for (i=0; i<ndims; i++){
        sprintf(dimname, "dim0_%d", i);
        status = ncmpi_def_dim(ncid, dimname, array_of_gsizes[i], &dimids0[i]);
        TEST_HANDLE_ERR(status);
    }

    /* define variables */
    varid = (int *)malloc(nvars*sizeof(int));
    for (i=0; i<nvars; i++) {
	sprintf(varname, "var0_%d", i);
	status = ncmpi_def_var(ncid, varname, NC_INT, ndims, dimids0, &varid[i]);
	TEST_HANDLE_ERR(status);
    }
     
    if (par_dim_id != 0) {
        for (i=0; i<nvars; i++) {
            status = ncmpi_put_att_int (ncid, varid[i], "par_dim_id",
                                        NC_INT, 1, &dimids0[par_dim_id]);
            if (status != NC_NOERR) handle_error(status);
        }
    }

    //set all non-record variable to be subfiled
    //MPI_Info_set(info, "nc_num_subfiles", "2");
    //status = ncmpi_set_var_info(ncid, varid, info);
    //TEST_HANDLE_ERR(status);
    
    status = ncmpi_enddef(ncid);
    TEST_HANDLE_ERR(status);

#if 0    
    if (mynod == 0)
        printf("*** Testing to write 1 non-record variable by using ncmpi_put_vara_all() ...");
#endif    
    stim = MPI_Wtime();
    for (i=0; i<nvars; i++) {
        status = ncmpi_put_vara_all(ncid, varid[i],
                                    starts_list[i], count_list[i],
                                    (const void *)&buf[i][0],
                                    bufcount_list[i], MPI_INT);
        TEST_HANDLE_ERR(status);
    }
    write_tim = MPI_Wtime() - stim;
            
    MPI_Allreduce(&write_tim, &new_write_tim, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    if (mynod == 0) {
        write_bw = ((double)array_of_gsizes[0]*(double)array_of_gsizes[1]*(double)array_of_gsizes[2]*(double)sizeof(int)*(double)nvars)/(new_write_tim*1024.0*1024.0);
        fprintf(stderr, "Global array size %d x %d x %d integers\n", array_of_gsizes[0], array_of_gsizes[1], array_of_gsizes[2]);
        fprintf(stderr, "Collective write time = %f sec, Collective write bandwidth = %f Mbytes/sec\n", new_write_tim, write_bw);
    }

    status = ncmpi_get_file_info(ncid, &info_used);
    TEST_HANDLE_ERR(status);

    stim = MPI_Wtime();
    ncmpi_close(ncid);
    close_tim = MPI_Wtime() - stim;
    
    MPI_Allreduce(&close_tim, &new_close_tim, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    if (mynod == 0) {
        fprintf(stderr, "close time = %f sec\n", new_close_tim);
    }

    goto end;

read:
    status = ncmpi_open(MPI_COMM_WORLD, filename, NC_NOWRITE, MPI_INFO_NULL, &ncid);
    TEST_HANDLE_ERR(status);
    
    stim = MPI_Wtime();
    
    /**
     * Inquire the dataset definitions of input dataset AND
     * Add dataset definitions for output dataset.
     */

    status = ncmpi_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
    if (status != NC_NOERR) handle_error(status);

    for (i=0; i<nvars; i++) {
        status = ncmpi_get_vara_all(ncid, i,
                                    starts_list[i], count_list[i],
                                    &(buf[i][0]), bufcount_list[i], MPI_INT);
        TEST_HANDLE_ERR(status);
    }
    read_tim = MPI_Wtime() - stim;
    
    MPI_Allreduce(&read_tim, &new_read_tim, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);
    
    if (mynod == 0) {
        read_bw = ((double)array_of_gsizes[0]*(double)array_of_gsizes[1]*(double)array_of_gsizes[2]*sizeof(int)*(double)nvars)/(new_read_tim*1024.0*1024.0);
        fprintf(stderr, "Collective read time = %f sec, Collective read bandwidth = %f Mbytes/sec\n", new_read_tim, read_bw);
    }

    status = ncmpi_get_file_info(ncid, &info_used);
    TEST_HANDLE_ERR(status);

end:    
    MPI_Info_free(&info);

    for (i=0; i<nvars; i++){
        free(buf[i]);
        free(starts_list[i]);
        free(count_list[i]);
    }
    free(buf);
    free(bufcount_list);
    free(datatype_list);
    if (!do_read) free(varid);
    free(starts_list);
    free(count_list);
    free(basename1);

    MPI_Reduce(&nerr, &sum_nerr, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mynod == 0) {
        char cmd_str[80];
        sprintf(cmd_str, "*** TESTING C   %s for subfiling", argv[0]);
        if (sum_nerr == 0)
            printf("%-66s ------ pass\n", cmd_str);
        else
            printf("%-66s ------ failed\n", cmd_str);
    }

    MPI_Finalize();

    return 0;
}
