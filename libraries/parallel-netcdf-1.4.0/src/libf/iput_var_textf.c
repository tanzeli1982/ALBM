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
#define nfmpi_iput_var_text_ NFMPI_IPUT_VAR_TEXT
#elif defined(F77_NAME_LOWER_2USCORE)
#define nfmpi_iput_var_text_ nfmpi_iput_var_text__
#elif !defined(F77_NAME_LOWER_USCORE)
#define nfmpi_iput_var_text_ nfmpi_iput_var_text
/* Else leave name alone */
#endif


/* Prototypes for the Fortran interfaces */
#include "mpifnetcdf.h"
FORTRAN_API int FORT_CALL nfmpi_iput_var_text_ ( int *v1, int *v2, char *v3 FORT_MIXED_LEN(d3), MPI_Fint *v4 FORT_END_LEN(d3) ){
    int ierr;
    int l2 = *v2 - 1;
    ierr = ncmpi_iput_var_text( *v1, l2, v3, v4 );
    return ierr;
}