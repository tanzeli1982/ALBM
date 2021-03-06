/* -*- Mode: C; c-basic-offset:4 ; -*- */
/*  
 *  (C) 2003 by Argonne National Laboratory and Northwestern University.
 *      See COPYRIGHT in top-level directory.
 *
 * This file is automatically generated by ./buildiface -infile=../lib/pnetcdf.h -deffile=defs
 * DO NOT EDIT
 */
#include "mpinetcdf_impl.h"


#ifdef F77_NAME_UPPER
#define nfmpi_iget_var_int8_ NFMPI_IGET_VAR_INT8
#elif defined(F77_NAME_LOWER_2USCORE)
#define nfmpi_iget_var_int8_ nfmpi_iget_var_int8__
#elif !defined(F77_NAME_LOWER_USCORE)
#define nfmpi_iget_var_int8_ nfmpi_iget_var_int8
/* Else leave name alone */
#endif


/* Prototypes for the Fortran interfaces */
#include "mpifnetcdf.h"
FORTRAN_API int FORT_CALL nfmpi_iget_var_int8_ ( int *v1, int *v2, long long * v3, MPI_Fint *v4 ){
    int ierr;
    int l2 = *v2 - 1;
    ierr = ncmpi_iget_var_longlong( *v1, l2, v3, v4 );
    return ierr;
}
