/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 *
 *  $Id: test_varm.c 1468 2013-10-26 16:53:18Z wkliao $
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <pnetcdf.h>


#define ERRCODE 2
#define ERR(e) {printf("Error at line %d: err=%d %s\n", __LINE__, e, ncmpi_strerror(e)); exit(ERRCODE);}

/*----< main() >------------------------------------------------------------*/
int main(int argc, char **argv)
{
    int i, j, retval, err=0, rank, nprocs, verbose;
    int ncid, dimid[2], varid, req, status;

    MPI_Offset start[2];
    MPI_Offset count[2];
    MPI_Offset stride[2];
    MPI_Offset imap[2];
    int   var[6][4];
    float rh[4][6];
    signed char  varT[4][6];
    char *filename="testfile.nc";

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    if (argc > 2) {
        if (!rank) printf("Usage: %s [filename]\n",argv[0]);
        MPI_Finalize();
        return 0;
    }
    if (argc == 2) filename = argv[1];

    verbose = 0;
    if (nprocs > 1 && rank == 0 && verbose)
        printf("Warning: %s is designed to run on 1 process\n", argv[0]);

    if (NC_NOERR != (retval = ncmpi_create(MPI_COMM_WORLD, filename,
        NC_CLOBBER | NC_64BIT_DATA, MPI_INFO_NULL, &ncid)))
       ERR(retval);

    /* define a variable of a 6 x 4 integer array in the nc file */
    if (NC_NOERR != (retval = ncmpi_def_dim(ncid, "Y", 6, &dimid[0]))) ERR(retval);
    if (NC_NOERR != (retval = ncmpi_def_dim(ncid, "X", 4, &dimid[1]))) ERR(retval);
    if (NC_NOERR != (retval = ncmpi_def_var(ncid, "var", NC_INT, 2, dimid, &varid)))
        ERR(retval);
    if (NC_NOERR != (retval = ncmpi_enddef(ncid))) ERR(retval);

    /* create a 6 x 4 integer variable in the file with contents:
           0,  1,  2,  3,
           4,  5,  6,  7,
           8,  9, 10, 11,
          12, 13, 14, 15,
          16, 17, 18, 19,
          20, 21, 22, 23
     */
    for (j=0; j<6; j++) for (i=0; i<4; i++) var[j][i] = j*4+i;

    start[0] = 0; start[1] = 0;
    count[0] = 6; count[1] = 4;
    if (rank > 0) count[0] = count[1] = 0;
    if (NC_NOERR != (retval = ncmpi_put_vara_int_all(ncid, varid, start, count, &var[0][0])))
        ERR(retval);

    /* read the variable back in the matrix transposed way, rh is 4 x 6 */
    count[0] = 6; count[1] = 4;
    stride[0] = 1; stride[1] = 1;
    imap[0]   = 1; imap[1] = 6;   /* would be {4, 1} if not transposing */
#define TEST_NON_BLOCKING_API
#ifdef TEST_NON_BLOCKING_API
    if (NC_NOERR != (retval = ncmpi_iget_varm_float(ncid, varid, start, count, stride, imap, &rh[0][0], &req)))
        ERR(retval);

    if (NC_NOERR != (retval = ncmpi_wait_all(ncid, 1, &req, &status)))
        ERR(retval);

    if (status != NC_NOERR) ERR(status);
#else
    if (NC_NOERR != (retval = ncmpi_get_varm_float_all(ncid, varid, start, count, stride, imap, &rh[0][0])))
        ERR(retval);
#endif

    /* check the contents of read */
    float k = 0.0;
    for (i=0; i<6; i++) {
        for (j=0; j<4; j++) {
            if (rh[j][i] != k) {
#ifdef PRINT_ERR_ON_SCREEN
                printf("get_varm unexpected value at j=%d i=%d\n",j,i);
#endif
                err++;
                break;
            }
            k += 1.0;
        }
    }
#ifdef PRINT_ON_SCREEN
    /* print the contents of read */
    for (j=0; j<4; j++) {
        printf("[%2d]: ",j);
        for (i=0; i<6; i++) {
            printf("%5.1f",rh[j][i]);
        }
        printf("\n");
    }
#endif
    /* the stdout should be:
           [ 0]:   0.0  4.0  8.0 12.0 16.0 20.0
           [ 1]:   1.0  5.0  9.0 13.0 17.0 21.0
           [ 2]:   2.0  6.0 10.0 14.0 18.0 22.0
           [ 3]:   3.0  7.0 11.0 15.0 19.0 23.0
     */

    /* testing get_varm(), first zero-out the variable in the file */
    memset(&var[0][0], 0, 6*4*sizeof(int));
    start[0] = 0; start[1] = 0;
    count[0] = 6; count[1] = 4;
    if (rank > 0) count[0] = count[1] = 0;
    if (NC_NOERR != (retval = ncmpi_put_vara_int_all(ncid, varid, start, count, &var[0][0])))
        ERR(retval);

    /* set the contents of the write buffer varT, a 4 x 6 char array
          50, 51, 52, 53, 54, 55,
          56, 57, 58, 59, 60, 61,
          62, 63, 64, 65, 66, 67,
          68, 69, 70, 71, 72, 73
     */
    for (j=0; j<4; j++) for (i=0; i<6; i++) varT[j][i] = j*6+i + 50;

    /* write varT to the NC variable in the matrix transposed way */
    start[0]  = 0; start[1]  = 0;
    count[0]  = 6; count[1]  = 4;
    stride[0] = 1; stride[1] = 1;
    imap[0]   = 1; imap[1]   = 6;   /* would be {4, 1} if not transposing */
    if (rank > 0) count[0] = count[1] = 0;
#ifdef TEST_NON_BLOCKING_API
    if (NC_NOERR != (retval = ncmpi_iput_varm_schar(ncid, varid, start, count, stride, imap, &varT[0][0], &req)))
        ERR(retval);

    if (NC_NOERR != (retval = ncmpi_wait_all(ncid, 1, &req, &status)))
        ERR(retval);

    if (status != NC_NOERR) ERR(status);
#else
    if (NC_NOERR != (retval = ncmpi_put_varm_schar_all(ncid, varid, start, count, stride, imap, &varT[0][0])))
        ERR(retval);
#endif

    /* the output from command "ncmpidump -v var test.nc" should be:
           var =
            50, 56, 62, 68,
            51, 57, 63, 69,
            52, 58, 64, 70,
            53, 59, 65, 71,
            54, 60, 66, 72,
            55, 61, 67, 73 ;
     */

    /* check if the contents of write buffer have been altered */
    for (j=0; j<4; j++) {
        for (i=0; i<6; i++) {
            if (varT[j][i] != j*6+i + 50) {
#ifdef PRINT_ERR_ON_SCREEN
                /* this error is a pntecdf internal error, if occurs */
                printf("Error: get_varm write buffer has been altered at j=%d i=%d\n",j,i);
#endif
                err++;
                break;
            }
        }
    }
    if (NC_NOERR != (retval = ncmpi_close(ncid))) ERR(retval);

    if (rank == 0) {
        char cmd_str[80];
        sprintf(cmd_str, "*** TESTING C   %s for get/put varm ", argv[0]);

        if (err)
            printf("%-66s ------ failed\n", cmd_str);
        else
            printf("%-66s ------ pass\n", cmd_str);
    }

    MPI_Finalize();

    return err;
}

