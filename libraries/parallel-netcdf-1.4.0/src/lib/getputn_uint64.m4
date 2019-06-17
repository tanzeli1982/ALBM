dnl Process this m4 file to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
/* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: getputn_uint64.m4 1468 2013-10-26 16:53:18Z wkliao $ */

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

/*---- x_uint64 -------------------------------------------------------------*/

/*
static void
get_ix_uint64(const void *xp, uint64 *ip)
{
    *ip = *((uint64*)xp);
#ifndef WORDS_BIGENDIAN
    SWAP8B(ip);
#endif
}

static void
put_ix_uint64(void *xp, const uint64 *ip)
{
    *((uint64*) xp) = *ip;
#ifndef WORDS_BIGENDIAN
    SWAP8B(xp);
#endif
}
*/

static void
get_ix_uint64(const void *xp, uint64 *ip)
{
    /* are these bit shifting faster than byte swap? */
    const uchar *cp = (const uchar *) xp;

    *ip  = ((uint64)(*cp++) << 56);
    *ip |= ((uint64)(*cp++) << 48);
    *ip |= ((uint64)(*cp++) << 40);
    *ip |= ((uint64)(*cp++) << 32);
    *ip |= ((uint64)(*cp++) << 24);
    *ip |= ((uint64)(*cp++) << 16);
    *ip |= ((uint64)(*cp++) <<  8);
    *ip |=  (uint64)*cp;
}

static void
put_ix_uint64(void *xp, const uint64 *ip)
{
    uchar *cp = (uchar *) xp;

    *cp++ = (*ip) >> 56;
    *cp++ = ((*ip) & 0x00ff000000000000ULL) >> 48;
    *cp++ = ((*ip) & 0x0000ff0000000000ULL) >> 40;
    *cp++ = ((*ip) & 0x000000ff00000000ULL) >> 32;
    *cp++ = ((*ip) & 0x00000000ff000000ULL) >> 24;
    *cp++ = ((*ip) & 0x0000000000ff0000ULL) >> 16;
    *cp++ = ((*ip) & 0x000000000000ff00ULL) >>  8;
    *cp   = ((*ip) & 0x00000000000000ffULL);
}

dnl
dnl GET_UINT64(xp, ip)
dnl
define(`GET_UINT64',dnl
`dnl
/*----< ncmpix_get_uint64_$1() >----------------------------------------------*/
static int
ncmpix_get_uint64_$1(const void *xp, $1 *ip)
{
    uint64 xx;
    get_ix_uint64(xp, &xx);
    *ip = xx;
    $2       /* check if can fit into $1 */
    return NC_NOERR;
}
')dnl

/* for smaller-sized   signed types, check if the got uint64 is too big (schar, short, int, long, float)
 * for smaller-sized unsigned types, check if the got uint64 is too big (uchar, ushort, uint)
 * for equal-sized     signed types, check if the got uint64 is too big (double, int64)
 * for equal-sized   unsigned types, no check is needed (uint64)
 */

GET_UINT64(schar,  if (xx > SCHAR_MAX) return NC_ERANGE;)
GET_UINT64(uchar,  if (xx > UCHAR_MAX) return NC_ERANGE;)
GET_UINT64(short,  if (xx > SHRT_MAX)  return NC_ERANGE;)
GET_UINT64(int,    if (xx > INT_MAX)   return NC_ERANGE;)
#if SIZEOF_LONG == X_SIZEOF_INT
static int
ncmpix_get_uint64_long(const void *xp, long *ip)
{
    return ncmpix_get_uint64_int(xp, (int*)ip);
}
#else
GET_UINT64(long,   if (xx > LONG_MAX)  return NC_ERANGE;)
#endif
GET_UINT64(float,  if (xx > FLT_MAX)   return NC_ERANGE;)
GET_UINT64(double, if (xx > DBL_MAX)   return NC_ERANGE;)
GET_UINT64(ushort, if (xx > USHRT_MAX) return NC_ERANGE;)
GET_UINT64(uint,   if (xx > UINT_MAX)  return NC_ERANGE;)
GET_UINT64(int64,  if (xx > LLONG_MAX) return NC_ERANGE;)
dnl GET_UINT64(uint64)

dnl
dnl PUT_UINT64(xp, ip)
dnl
define(`PUT_UINT64',dnl
`dnl
/*----< ncmpix_put_uint64_$1() >----------------------------------------------*/
static int
ncmpix_put_uint64_$1(void *xp, const $1 *ip)
{
    uint64 xx = (uint64) *ip;
    put_ix_uint64(xp, &xx);
    $2       /* check if can fit into uint64 */
    return NC_NOERR;
}
')dnl

/* for smaller-sized   signed types, check if the put value is negative (schar, short, int, long, float)
 * for smaller-sized unsigned types, no check is needed (uchar, ushort, uint)
 * for equal-sized     signed types, check if the put value is negative (double, int64)
 * for equal-sized   unsigned types, no check is needed (uint64)
 */

PUT_UINT64(uchar)
PUT_UINT64(uint)
PUT_UINT64(ushort)
dnl PUT_UINT64(uint64)
PUT_UINT64(schar,  if (*ip < 0) return NC_ERANGE;)
PUT_UINT64(short,  if (*ip < 0) return NC_ERANGE;)
PUT_UINT64(int,    if (*ip < 0) return NC_ERANGE;)
PUT_UINT64(long,   if (*ip < 0) return NC_ERANGE;)
PUT_UINT64(float,  if (*ip < 0) return NC_ERANGE;)
PUT_UINT64(double, if (*ip < 0) return NC_ERANGE;)
PUT_UINT64(int64,  if (*ip < 0) return NC_ERANGE;)


dnl
dnl GETN_UINT64(xpp, nelems, tp)
dnl
define(`GETN_UINT64',dnl
`dnl
/*----< ncmpix_getn_uint64_$1() >---------------------------------------------*/
int
ncmpix_getn_uint64_$1(const void **xpp, MPI_Offset nelems, $1 *tp)
{
    const char *xp = (const char *) *xpp;
    int status = NC_NOERR;

    for ( ; nelems != 0; nelems--, xp += X_SIZEOF_UINT64, tp++) {
        const int lstatus = ncmpix_get_uint64_$1(xp, tp);
        if (lstatus != NC_NOERR) status = lstatus;
    }

    *xpp = (void *)xp;
    return status;
}
')dnl

GETN_UINT64(schar)
GETN_UINT64(uchar)
GETN_UINT64(short)
GETN_UINT64(int)
GETN_UINT64(long)
GETN_UINT64(float)
GETN_UINT64(double)
GETN_UINT64(ushort)
GETN_UINT64(uint)
GETN_UINT64(int64)

/*----< ncmpix_getn_uint64_uint64() >----------------------------------------*/
/* optimized version */
int
ncmpix_getn_uint64_uint64(const void **xpp, MPI_Offset nelems, uint64 *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(tp, *xpp, nelems * sizeof(uint64));
# else
    ncmpii_swapn8b(tp, *xpp, nelems);
# endif
    *xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_UINT64);
    return NC_NOERR;
}

dnl
dnl PUTN_UINT64(xpp, nelems, tp)
dnl
define(`PUTN_UINT64',dnl
`dnl
/*----< ncmpix_putn_uint64_$1() >---------------------------------------------*/
int
ncmpix_putn_uint64_$1(void **xpp, MPI_Offset nelems, const $1 *tp)
{
    char *xp = (char *) *xpp;
    int status = NC_NOERR;

    for ( ; nelems != 0; nelems--, xp += X_SIZEOF_UINT64, tp++) {
        int lstatus = ncmpix_put_uint64_$1(xp, tp);
        if (lstatus != NC_NOERR) status = lstatus;
    }

    *xpp = (void *)xp;
    return status;
}
')dnl

PUTN_UINT64(schar)
PUTN_UINT64(uchar)
PUTN_UINT64(short)
PUTN_UINT64(int)
PUTN_UINT64(long)
PUTN_UINT64(float)
PUTN_UINT64(double)
PUTN_UINT64(ushort)
PUTN_UINT64(uint)
PUTN_UINT64(int64)

/*----< ncmpix_putn_uint64_uint64() >----------------------------------------*/
/* optimized version */
int
ncmpix_putn_uint64_uint64(void **xpp, MPI_Offset nelems, const uint64 *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(*xpp, tp, nelems * X_SIZEOF_UINT64);
# else
    ncmpii_swapn8b(*xpp, tp, nelems);
# endif
    *xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_UINT64);
    return NC_NOERR;
}


