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
#define nfmpi_get_var_ NFMPI_GET_VAR
#elif defined(F77_NAME_LOWER_2USCORE)
#define nfmpi_get_var_ nfmpi_get_var__
#elif !defined(F77_NAME_LOWER_USCORE)
#define nfmpi_get_var_ nfmpi_get_var
/* Else leave name alone */
#endif


/* Prototypes for the Fortran interfaces */
#include "mpifnetcdf.h"
FORTRAN_API int FORT_CALL nfmpi_get_var_ ( int *v1, int *v2, void *v3, MPI_Offset *v4, MPI_Fint *v5 ){
    int ierr;
    int l2 = *v2 - 1;
    ierr = ncmpi_get_var( *v1, l2, v3, *v4, MPI_Type_f2c(*v5) );
    return ierr;
}
