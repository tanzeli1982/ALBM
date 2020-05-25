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
#define nfmpi_issyserr_ NFMPI_ISSYSERR
#elif defined(F77_NAME_LOWER_2USCORE)
#define nfmpi_issyserr_ nfmpi_issyserr__
#elif !defined(F77_NAME_LOWER_USCORE)
#define nfmpi_issyserr_ nfmpi_issyserr
/* Else leave name alone */
#endif


/* Prototypes for the Fortran interfaces */
#include "mpifnetcdf.h"
FORTRAN_API int FORT_CALL nfmpi_issyserr_ ( MPI_Fint *v1 ){
    if (*v1 > 0)
      return 1;
    else
      return 0;
}