/***********************************************************
 *
 * This test program writes a netCDF file using the parallel
 * netCDF library using MPI-IO. 
 *
 * The output file is: "testwrite.nc"
 *
 *  The CDL notation for the test dataset is shown below:
 *
 *    netcdf test {
 *
 *       dimensions:
 *
 *            x = 100, y = 100, z = 100, time = NC_UNLIMITED;
 *
 *
 *       variables:  // variable types, names, shapes, attributes
 *
 *            double square(x, y);
 *                     squre: description = "2-D integer array";
 *
 *            double cube(x,y,z);
 *
 *            double time(time);  // coordinate & record variable
 *
 *            double xytime(time, x, y);  // record variable
 *
 *
 *      // global attributes
 *
 *           :title = "example netCDF dataset";
 *
 *
 *      data:  // data written for variables
 *          square  = 0, 1, 2, 3,  ... , 9999;
 *          cube    = 0, 1, 2, 3,  ... , 999999;
 *	    time    = 0, 1, 2, 3,  ... , 99;    // 100 records
 *          xytime  = 0, 1, 2, 3,  ... , 9999;  // 100 records
 *   }
 *
 *
 *
 * This test uses collective APIs to write variable data and only
 * deals with double variables. 
 *
 * This test assume # of processors = 4
 *
 **********************************************************/



#include <mpi.h>
#include <stdio.h>
#include <pnetcdf.h>
#include <string.h>
#include "testutils.h"

/* Prototype for functions used only in this file */
static void handle_error(int status);

static void handle_error(int status) {
  fprintf(stderr, "%s\n", ncmpi_strerror(status));
}

int main(int argc, char **argv) {

  int i, j, k;
  int status;
  int ncid;
  int dimid1, dimid2, dimid3, udimid;
  int square_dim[2], cube_dim[3], xytime_dim[3], time_dim[1];
  MPI_Offset square_start[2], cube_start[3] = {0, 0, 0};
  MPI_Offset square_count[2] = {50, 50}, cube_count[3] = {100, 50, 50};
  MPI_Offset square_stride[2] = {2, 2};
  MPI_Offset xytime_start[3] = {0, 0, 0};
  MPI_Offset xytime_count[3] = {100, 50, 50};
  MPI_Offset time_start[1], time_count[1] = {25};
  int square_id, cube_id, xytime_id, time_id;
  static char title[] = "example netCDF dataset";
  static char description[] = "2-D integer array";
  double data[100][50][50], buffer[100];
  double stride_2d_data[50][50];
  int rank;
  int nprocs;
  MPI_Comm comm = MPI_COMM_WORLD;
  params opts;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0) 
	  fprintf(stderr, "Testing write ... ");
  parse_write_args(argc, argv, rank, &opts);

  /**********  START OF NETCDF ACCESS **************/

  /**
   * Create the dataset
   *   File name: "testwrite.nc"
   *   Dataset API: Collective
   */

  status = ncmpi_create(comm, opts.outfname, NC_CLOBBER, MPI_INFO_NULL, &ncid);
  if (status != NC_NOERR) handle_error(status);


  /**
   * Create a global attribute:
   *    :title = "example netCDF dataset";
   */

  status = ncmpi_put_att_text (ncid, NC_GLOBAL, "title",
                          strlen(title), title);
  if (status != NC_NOERR) handle_error(status);

  /**
   * Add 4 pre-defined dimensions:
   *   x = 100, y = 100, z = 100, time = NC_UNLIMITED
   */

  status = ncmpi_def_dim(ncid, "x", 100L, &dimid1);
  if (status != NC_NOERR) handle_error(status);
  status = ncmpi_def_dim(ncid, "y", 100L, &dimid2);
  if (status != NC_NOERR) handle_error(status);
  status = ncmpi_def_dim(ncid, "z", 100L, &dimid3);
  if (status != NC_NOERR) handle_error(status);
  status = ncmpi_def_dim(ncid, "time", NC_UNLIMITED, &udimid);
  if (status != NC_NOERR) handle_error(status);

  /**
   * Define the dimensionality and then add 4 variables:
   *    square(x, y), cube(x,y,z), time(time), xytime(time, x, y)  
   */

  square_dim[0] = cube_dim[0] = xytime_dim[1] = dimid1;
  square_dim[1] = cube_dim[1] = xytime_dim[2] = dimid2;
  cube_dim[2] = dimid3;
  xytime_dim[0] = udimid;
  time_dim[0] = udimid;
  status = ncmpi_def_var (ncid, "square", NC_DOUBLE, 2, square_dim, &square_id);
  if (status != NC_NOERR) handle_error(status);
  status = ncmpi_def_var (ncid, "cube", NC_DOUBLE, 3, cube_dim, &cube_id);
  if (status != NC_NOERR) handle_error(status);
  status = ncmpi_def_var (ncid, "time", NC_DOUBLE, 1, time_dim, &time_id);
  if (status != NC_NOERR) handle_error(status);
  status = ncmpi_def_var (ncid, "xytime", NC_DOUBLE, 3, xytime_dim, &xytime_id);
  if (status != NC_NOERR) handle_error(status);

  /**
   * Add an attribute for variable: 
   *    square: decsription = "2-D integer array"
   */

  status = ncmpi_put_att_text (ncid, square_id, "description",
                          strlen(description), description);
  if (status != NC_NOERR) handle_error(status);


  /**
   * End Define Mode (switch to data mode)
   *   Dataset API: Collective
   */

  status = ncmpi_enddef(ncid);
  if (status != NC_NOERR) handle_error(status);

  /**
   * Data Partition (Assume 4 processors):
   *   square: 2-D, (Cyclic, Cyclic), 50*50 from 100*100, strided access
   *   cube:   3-D, (*, Block, Block), 100*50*50 from 100*100*100
   *   xytime: 3-D, (*, Block, Block), 100*50*50 from 100*100*100
   *   time:   1-D, Block-wise, 25 from 100
   */

  /* square_start[0] = */
  cube_start[1] = xytime_start[1] = (rank/2) * 50;
  /* square_start[1] = */
  cube_start[2] = xytime_start[2] = (rank%2) * 50;
  time_start[0] = (rank%4) * 25;
  square_start[0] = rank/2;
  square_start[1] = rank%2;


  /**
   * Packing data in the buffer 
   */

  /* Data for variable: time */
  for ( i = time_start[0]; i < time_start[0] + time_count[0]; i++ )
    buffer[i - time_start[0]] = i;   

  /* Data for variable: cube and xytime */
  for ( i = 0; i < 100; i++ )
    for ( j = cube_start[1]; j < cube_start[1]+cube_count[1]; j++ )
      for ( k = cube_start[2]; k < cube_start[2]+cube_count[2]; k++ )
        data[i][j-cube_start[1]][k-cube_start[2]] = i*100*100 + j*100 + k;

  /* Data for variable: square */
  for ( i = 0; i < 50; i ++ )
    for ( j = 0; j < 50; j++ )
      stride_2d_data[i][j] = (2*i + rank/2)*100 + (2*j + rank%2);

  /**
   * Write data into variables: square, cube, time and xytime  
   *   Access Method: subarray
   *   Data Mode API: collective
   */ 
  
  status = ncmpi_put_vars_double_all(ncid, square_id,
                    square_start, square_count, square_stride,
                    &stride_2d_data[0][0]);
  if (status != NC_NOERR) handle_error(status);
  status = ncmpi_put_vara_double_all(ncid, cube_id,
                    cube_start, cube_count,
                    &data[0][0][0]);
  if (status != NC_NOERR) handle_error(status);
  status = ncmpi_put_vara_double_all(ncid, time_id,
                    time_start, time_count,
                    (void *)buffer);
  if (status != NC_NOERR) handle_error(status);
  status = ncmpi_put_vara_double_all(ncid, xytime_id,
                    xytime_start, xytime_count,
                    &data[0][0][0]);
  if (status != NC_NOERR) handle_error(status);

  /**
   * Close the dataset
   *   Dataset API:  collective
   */

  status = ncmpi_close(ncid);
  if (status != NC_NOERR) handle_error(status);

  /*******************  END OF NETCDF ACCESS  ****************/

if (rank == 0)
  fprintf(stderr, "OK\nFile written to: %s!\n", opts.outfname);

  MPI_Finalize();
  return 0;
}

