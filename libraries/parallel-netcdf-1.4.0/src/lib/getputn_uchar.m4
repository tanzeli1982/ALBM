dnl Process this m4 file to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
/* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: getputn_uchar.m4 1468 2013-10-26 16:53:18Z wkliao $ */

#if HAVE_CONFIG_H
# include <ncconfig.h>
#endif

#include <mpi.h>

#include "nc.h"
#include "ncx.h"

/* ftype is the variable's nc_type defined in file, eg. int64
 * btype is the I/O buffer's C data type, eg. long long
 * buftype is I/O bufer's MPI data type, eg. MPI_UNSIGNED_LONG_LONG
 * apitype is data type appeared in the API names, eg. ncmpi_get_vara_longlong
 */

/*---- uchar ----------------------------------------------------------------*/

/*----< ncmpix_getn_uchar_schar() >------------------------------------------*/
int
ncmpix_getn_uchar_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
    /* file type is uchar, buffer type is schar */
    int status = NC_NOERR;
    uchar *xp = (uchar *) *xpp;

    while (nelems-- > 0) {
        if (*xp > X_SCHAR_MAX) /* uchar might be too big for schar */
            status = NC_ERANGE;
        *tp++ = *xp++;
    }

    return status;
}

/*----< ncmpix_pad_getn_uchar_schar() >--------------------------------------*/
int
ncmpix_pad_getn_uchar_schar(const void **xpp, MPI_Offset nelems, schar *tp)
{
    int status = NC_NOERR;
    uchar *xp = (uchar *) *xpp;

    MPI_Offset rndup = nelems % X_ALIGN;
    if (rndup) rndup = X_ALIGN - rndup;

    while (nelems-- > 0) {
        if (*xp > X_SCHAR_MAX) /* uchar might be too big for schar */
            status = NC_ERANGE;
        *tp++ = *xp++;
    }
    *xpp = (void *)(xp + rndup);

    return status;
}

/*----< ncmpix_getn_uchar_uchar() >------------------------------------------*/
int
ncmpix_getn_uchar_uchar(const void **xpp, MPI_Offset nelems, uchar *tp)
{
    /* file type is uchar, buffer type is uchar */
    memcpy(tp, *xpp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems);

    return NC_NOERR;
}

dnl
dnl GETN_UCHAR(xpp, nelems, tp)
dnl
define(`GETN_UCHAR',dnl
`dnl
/*----< ncmpix_getn_uchar_$1() >----------------------------------------------*/
int
ncmpix_getn_uchar_$1(const void **xpp, MPI_Offset nelems, $1 *tp)
{
    /* there is no ENDIANness issue, as uchar is 1 byte */
    uchar *xp = (uchar *) *xpp;
    while (nelems-- != 0)
        *tp++ = *xp++;
    *xpp = (const void *)xp;

    return NC_NOERR;
}
')dnl

GETN_UCHAR(short)
GETN_UCHAR(int)
GETN_UCHAR(long)
GETN_UCHAR(float)
GETN_UCHAR(double)
GETN_UCHAR(ushort)
GETN_UCHAR(uint)
GETN_UCHAR(int64)
GETN_UCHAR(uint64)

/*----< ncmpix_pad_getn_uchar_uchar() >--------------------------------------*/
int
ncmpix_pad_getn_uchar_uchar(const void **xpp, MPI_Offset nelems, uchar *tp)
{
    MPI_Offset rndup = nelems % X_ALIGN;
    if (rndup) rndup = X_ALIGN - rndup;

    memcpy(tp, *xpp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems + rndup);

    return NC_NOERR;
}

dnl
dnl PAD_GETN_UCHAR(xpp, nelems, tp)
dnl
define(`PAD_GETN_UCHAR',dnl
`dnl
/*----< ncmpix_getn_uchar_$1() >----------------------------------------------*/
int
ncmpix_pad_getn_uchar_$1(const void **xpp, MPI_Offset nelems, $1 *tp)
{
    MPI_Offset rndup = nelems % X_ALIGN;
    if (rndup) rndup = X_ALIGN - rndup;

    /* there is no ENDIANness issue, as uchar is 1 byte */
    uchar *xp = (uchar *) *xpp;
    while (nelems-- != 0)
        *tp++ = *xp++;
    *xpp = (void *)((char *)(*xpp) + rndup);

    return NC_NOERR;
}
')dnl

PAD_GETN_UCHAR(short)
PAD_GETN_UCHAR(int)
PAD_GETN_UCHAR(long)
PAD_GETN_UCHAR(float)
PAD_GETN_UCHAR(double)
PAD_GETN_UCHAR(ushort)
PAD_GETN_UCHAR(uint)
PAD_GETN_UCHAR(int64)
PAD_GETN_UCHAR(uint64)


/*----< ncmpix_putn_uchar_uchar() >------------------------------------------*/
int
ncmpix_putn_uchar_uchar(void **xpp, MPI_Offset nelems, const uchar *tp)
{
    memcpy(*xpp, tp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems);

    return NC_NOERR;
}

/*----< ncmpix_pad_putn_uchar_uchar() >--------------------------------------*/
int
ncmpix_pad_putn_uchar_uchar(void **xpp, MPI_Offset nelems, const uchar *tp)
{
    MPI_Offset rndup = nelems % X_ALIGN;
    if (rndup) rndup = X_ALIGN - rndup;

    memcpy(*xpp, tp, nelems);
    *xpp = (void *)((char *)(*xpp) + nelems);

    if (rndup) {
        memcpy(*xpp, nada, rndup);
        *xpp = (void *)((char *)(*xpp) + rndup);
    }

    return NC_NOERR;
}

dnl
dnl PUTN_UCHAR(xpp, nelems, tp)
dnl
define(`PUTN_UCHAR',dnl
`dnl
/*----< ncmpix_putn_uchar_$1() >----------------------------------------------*/
int
ncmpix_putn_uchar_$1(void **xpp, MPI_Offset nelems, const $1 *tp)
{
    int status=NC_NOERR;
    uchar *xp = (uchar *) *xpp;

    while (nelems-- != 0) {
        $2           /* check if can fit into uchar */
        *xp++ = (uchar) *tp++;
    }
    *xpp = (void *)xp;

    return status;
}
')dnl

PUTN_UCHAR(schar,  if (*tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(short,  if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(int,    if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(long,   if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(float,  if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(double, if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(int64,  if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PUTN_UCHAR(ushort, if (*tp > X_UCHAR_MAX) status = NC_ERANGE;)
PUTN_UCHAR(uint,   if (*tp > X_UCHAR_MAX) status = NC_ERANGE;)
PUTN_UCHAR(uint64, if (*tp > X_UCHAR_MAX) status = NC_ERANGE;)

dnl
dnl PAD_PUTN_UCHAR(xpp, nelems, tp)
dnl
define(`PAD_PUTN_UCHAR',dnl
`dnl
/*----< ncmpix_pad_putn_uchar_$1() >------------------------------------------*/
int
ncmpix_pad_putn_uchar_$1(void **xpp, MPI_Offset nelems, const $1 *tp)
{
    int status=NC_NOERR;
    uchar *xp = (uchar *) *xpp;

    MPI_Offset rndup = nelems % X_ALIGN;
    if (rndup) rndup = X_ALIGN - rndup;

    while (nelems-- != 0) {
        $2         /* check if can fit into uchar */
        *xp++ = (uchar) *tp++;
    }

    if (rndup) {
        memcpy(xp, nada, rndup);
        xp += rndup;
    }
    *xpp = (void *)xp;

    return status;
}
')dnl

PAD_PUTN_UCHAR(schar,  if (*tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(short,  if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(int,    if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(long,   if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(float,  if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(double, if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(int64,  if (*tp > X_UCHAR_MAX || *tp < 0) status = NC_ERANGE;)
PAD_PUTN_UCHAR(ushort, if (*tp > X_UCHAR_MAX) status = NC_ERANGE;)
PAD_PUTN_UCHAR(uint,   if (*tp > X_UCHAR_MAX) status = NC_ERANGE;)
PAD_PUTN_UCHAR(uint64, if (*tp > X_UCHAR_MAX) status = NC_ERANGE;)


