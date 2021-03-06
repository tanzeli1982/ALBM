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
#define nfmpi_rename_att_ NFMPI_RENAME_ATT
#elif defined(F77_NAME_LOWER_2USCORE)
#define nfmpi_rename_att_ nfmpi_rename_att__
#elif !defined(F77_NAME_LOWER_USCORE)
#define nfmpi_rename_att_ nfmpi_rename_att
/* Else leave name alone */
#endif


/* Prototypes for the Fortran interfaces */
#include "mpifnetcdf.h"
FORTRAN_API int FORT_CALL nfmpi_rename_att_ ( int *v1, int *v2, char *v3 FORT_MIXED_LEN(d3), char *v4 FORT_MIXED_LEN(d4) FORT_END_LEN(d3) FORT_END_LEN(d4) ){
    int ierr;
    int l2 = *v2 - 1;
    char *p3;
    char *p4;

    {char *p = v3 + d3 - 1;
     int  li;
        while (*p == ' ' && p > v3) p--;
        p++;
        p3 = (char *)malloc( p-v3 + 1 );
        for (li=0; li<(p-v3); li++) { p3[li] = v3[li]; }
        p3[li] = 0; 
    }

    {char *p = v4 + d4 - 1;
     int  li;
        while (*p == ' ' && p > v4) p--;
        p++;
        p4 = (char *)malloc( p-v4 + 1 );
        for (li=0; li<(p-v4); li++) { p4[li] = v4[li]; }
        p4[li] = 0; 
    }
    ierr = ncmpi_rename_att( *v1, l2, p3, p4 );
    free( p3 );
    free( p4 );
    return ierr;
}
