dnl Process this m4 file to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
/* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: getputn_float.m4 1468 2013-10-26 16:53:18Z wkliao $ */

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

/*---- x_float --------------------------------------------------------------*/

#if X_SIZEOF_FLOAT == SIZEOF_FLOAT && !defined(NO_IEEE_FLOAT)

static void
get_ix_float(const void *xp, float *ip)
{
#ifdef WORDS_BIGENDIAN
    memcpy(ip, xp, sizeof(float));
#else
    ncmpii_swap4b(ip, xp);
#endif
}

static void
put_ix_float(void *xp, const float *ip)
{
#ifdef WORDS_BIGENDIAN
    memcpy(xp, ip, X_SIZEOF_FLOAT);
#else
    ncmpii_swap4b(xp, ip);
#endif
}

#elif defined(vax)

/* What IEEE single precision floating point looks like on a Vax */
struct ieee_single {
    unsigned int exp_hi       : 7;
    unsigned int sign         : 1;
    unsigned int mant_hi      : 7;
    unsigned int exp_lo       : 1;
    unsigned int mant_lo_hi   : 8;
    unsigned int mant_lo_lo   : 8;
};

/* Vax single precision floating point */
struct vax_single {
    unsigned int mantissa1 : 7;
    unsigned int exp       : 8;
    unsigned int sign      : 1;
    unsigned int mantissa2 : 16;
};

#define VAX_SNG_BIAS    0x81
#define IEEE_SNG_BIAS   0x7f

static struct sgl_limits {
    struct vax_single s;
    struct ieee_single ieee;
} max = {
    { 0x7f, 0xff, 0x0, 0xffff },        /* Max Vax */
    { 0x7f, 0x0,  0x0, 0x1, 0x0, 0x0 }  /* Max IEEE */
};
static struct sgl_limits min = {
    { 0x0, 0x0, 0x0, 0x0 },             /* Min Vax */
    { 0x0, 0x0, 0x0, 0x0, 0x0, 0x0 }    /* Min IEEE */
};

static void
get_ix_float(const void *xp, float *ip)
{
    struct vax_single *const vsp = (struct vax_single *) ip;
    const struct ieee_single *const isp =
             (const struct ieee_single *) xp;
    unsigned exp = isp->exp_hi << 1 | isp->exp_lo;

    switch(exp) {
        case 0 :
            /* ieee subnormal */
            if(isp->mant_hi == min.ieee.mant_hi
                && isp->mant_lo_hi == min.ieee.mant_lo_hi
                && isp->mant_lo_lo == min.ieee.mant_lo_lo) {
                *vsp = min.s;
            } else {
                unsigned mantissa = (isp->mant_hi << 16)
                     | isp->mant_lo_hi << 8
                     | isp->mant_lo_lo;
                unsigned tmp = mantissa >> 20;
                if (tmp >= 4) {
                    vsp->exp = 2;
                } else if (tmp >= 2) {
                    vsp->exp = 1;
                } else {
                    *vsp = min.s;
                    break;
                } /* else */
                tmp = mantissa - (1 << (20 + vsp->exp ));
                tmp <<= 3 - vsp->exp;
                vsp->mantissa2 = tmp;
                vsp->mantissa1 = (tmp >> 16);
            }
            break;
        case 0xfe :
        case 0xff :
            *vsp = max.s;
            break;
        default :
            vsp->exp = exp - IEEE_SNG_BIAS + VAX_SNG_BIAS;
            vsp->mantissa2 = isp->mant_lo_hi << 8 | isp->mant_lo_lo;
            vsp->mantissa1 = isp->mant_hi;
    }
    vsp->sign = isp->sign;
}


static void
put_ix_float(void *xp, const float *ip)
{
    const struct vax_single *const vsp =
             (const struct vax_single *)ip;
    struct ieee_single *const isp = (struct ieee_single *) xp;

    switch(vsp->exp){
        case 0 :
            /* all vax float with zero exponent map to zero */
            *isp = min.ieee;
            break;
        case 2 :
        case 1 :
        {
            /* These will map to subnormals */
            unsigned mantissa = (vsp->mantissa1 << 16)
                     | vsp->mantissa2;
            mantissa >>= 3 - vsp->exp;
            mantissa += (1 << (20 + vsp->exp));
            isp->mant_lo_lo = mantissa;
            isp->mant_lo_hi = mantissa >> 8;
            isp->mant_hi = mantissa >> 16;
            isp->exp_lo = 0;
            isp->exp_hi = 0;
        }
            break;
        case 0xff : /* max.s.exp */
            if ( vsp->mantissa2 == max.s.mantissa2 &&
                 vsp->mantissa1 == max.s.mantissa1) {
                /* map largest vax float to ieee infinity */
                *isp = max.ieee;
                break;
            } /* else, fall thru */
        default :
        {
            unsigned exp = vsp->exp - VAX_SNG_BIAS + IEEE_SNG_BIAS;
            isp->exp_hi = exp >> 1;
            isp->exp_lo = exp;
            isp->mant_lo_lo = vsp->mantissa2;
            isp->mant_lo_hi = vsp->mantissa2 >> 8;
            isp->mant_hi = vsp->mantissa1;
        }
    }
    isp->sign = vsp->sign;

}

    /* vax */
#else
#error "ix_float implementation"
#endif


dnl
dnl GET_FLOAT(xp, ip)
dnl
define(`GET_FLOAT',dnl
`dnl
/*----< ncmpix_get_float_$1() >-----------------------------------------------*/
static int
ncmpix_get_float_$1(const void *xp, $1 *ip)
{
    float xx;
    get_ix_float(xp, &xx);
    *ip = xx;
    $2              /* check if can fit into $1 */
    return NC_NOERR;
}
')dnl

GET_FLOAT(uchar,  if (xx > UCHAR_MAX  || xx < 0) return NC_ERANGE;)
GET_FLOAT(ushort, if (xx > USHRT_MAX  || xx < 0) return NC_ERANGE;)
GET_FLOAT(uint,   if (xx > UINT_MAX   || xx < 0) return NC_ERANGE;)
GET_FLOAT(uint64, if (xx > ULLONG_MAX || xx < 0) return NC_ERANGE;)

GET_FLOAT(schar,  if (xx > SCHAR_MAX  || xx < SCHAR_MIN) return NC_ERANGE;)
GET_FLOAT(short,  if (xx > SHRT_MAX   || xx < SHRT_MIN)  return NC_ERANGE;)
GET_FLOAT(int,    if (xx > INT_MAX    || xx < INT_MIN)   return NC_ERANGE;)
GET_FLOAT(long,   if (xx > LONG_MAX   || xx < LONG_MIN)  return NC_ERANGE;)
GET_FLOAT(int64,  if (xx > LLONG_MAX  || xx < LLONG_MIN) return NC_ERANGE;)
GET_FLOAT(double)

#define CHECK_FLOAT_RANGE \
    if ((float)(*ip) > X_FLOAT_MAX || (float)(*ip) < X_FLOAT_MIN) \
        return NC_ERANGE;

dnl
dnl PUT_FLOAT(xp, ip)
dnl
define(`PUT_FLOAT',dnl
`dnl
/*----< ncmpix_put_float_$1() >-----------------------------------------------*/
static int
ncmpix_put_float_$1(void *xp, const $1 *ip)
{
    float xx = (float) *ip;
    put_ix_float(xp, &xx);
    $2               /* check if can fit into float */
    return NC_NOERR;
}
')dnl

PUT_FLOAT(schar,)
PUT_FLOAT(uchar,)
PUT_FLOAT(short,)
PUT_FLOAT(ushort,)
PUT_FLOAT(int,)
PUT_FLOAT(uint,)

PUT_FLOAT(long,   CHECK_FLOAT_RANGE)
PUT_FLOAT(double, CHECK_FLOAT_RANGE)
PUT_FLOAT(int64,  CHECK_FLOAT_RANGE)
PUT_FLOAT(uint64, CHECK_FLOAT_RANGE)

dnl
dnl GETN_FLOAT(xp, ip)
dnl
define(`GETN_FLOAT',dnl
`dnl
/*----< ncmpix_getn_float_$1() >----------------------------------------------*/
int
ncmpix_getn_float_$1(const void **xpp, MPI_Offset nelems, $1 *tp)
{
    const char *xp = (const char *) *xpp;
    int status = NC_NOERR;

    for ( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++) {
        const int lstatus = ncmpix_get_float_$1(xp, tp);
        if (lstatus != NC_NOERR) status = lstatus;
    }

    *xpp = (void *)xp;
    return status;
}
')dnl

GETN_FLOAT(schar)
GETN_FLOAT(uchar)
GETN_FLOAT(short)
GETN_FLOAT(int)
GETN_FLOAT(long)
GETN_FLOAT(double)
GETN_FLOAT(ushort)
GETN_FLOAT(uint)
GETN_FLOAT(int64)
GETN_FLOAT(uint64)

/*----< ncmpix_getn_float_float() >------------------------------------------*/
#if X_SIZEOF_FLOAT == SIZEOF_FLOAT && !defined(NO_IEEE_FLOAT)
/* optimized version */
int
ncmpix_getn_float_float(const void **xpp, MPI_Offset nelems, float *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(tp, *xpp, nelems * sizeof(float));
# else
    ncmpii_swapn4b(tp, *xpp, nelems);
# endif
    *xpp = (const void *)((const char *)(*xpp) + nelems * X_SIZEOF_FLOAT);
    return NC_NOERR;
}
#elif defined(vax)
int
ncmpix_getn_float_float(const void **xpp, MPI_Offset nfloats, float *ip)
{
    float *const end = ip + nfloats;

    while (ip < end) {
        struct vax_single *const vsp = (struct vax_single *) ip;
        const struct ieee_single *const isp =
             (const struct ieee_single *) (*xpp);
        unsigned exp = isp->exp_hi << 1 | isp->exp_lo;

        switch (exp) {
            case 0 :
                /* ieee subnormal */
                if (isp->mant_hi == min.ieee.mant_hi       &&
                    isp->mant_lo_hi == min.ieee.mant_lo_hi &&
                    isp->mant_lo_lo == min.ieee.mant_lo_lo) {
                    *vsp = min.s;
                }
                else {
                    unsigned mantissa = (isp->mant_hi << 16)
                         | isp->mant_lo_hi << 8
                         | isp->mant_lo_lo;
                    unsigned tmp = mantissa >> 20;
                    if (tmp >= 4) {
                        vsp->exp = 2;
                    } else if (tmp >= 2) {
                        vsp->exp = 1;
                    } else {
                        *vsp = min.s;
                        break;
                    } /* else */
                    tmp = mantissa - (1 << (20 + vsp->exp ));
                    tmp <<= 3 - vsp->exp;
                    vsp->mantissa2 = tmp;
                    vsp->mantissa1 = (tmp >> 16);
                }
                break;
            case 0xfe :
            case 0xff :
                *vsp = max.s;
                break;
            default :
                vsp->exp = exp - IEEE_SNG_BIAS + VAX_SNG_BIAS;
                vsp->mantissa2 = isp->mant_lo_hi << 8 | isp->mant_lo_lo;
                vsp->mantissa1 = isp->mant_hi;
        }

        vsp->sign = isp->sign;

        ip++;
        *xpp = (char *)(*xpp) + X_SIZEOF_FLOAT;
    }
    return NC_NOERR;
}
#else

/*----< ncmpix_get_float_float() >-------------------------------------------*/
static int
ncmpix_get_float_float(const void *xp, float *ip)
{
    /* TODO */
    get_ix_float(xp, ip);
    return NC_NOERR;
}

GETN_FLOAT(float)
#endif

dnl
dnl PUTN_FLOAT(xp, ip)
dnl
define(`PUTN_FLOAT',dnl
`dnl
/*----< ncmpix_putn_float_$1() >----------------------------------------------*/
int
ncmpix_putn_float_$1(void **xpp, MPI_Offset nelems, const $1 *tp)
{
    char *xp = (char *) *xpp;
    int status = NC_NOERR;

    for ( ; nelems != 0; nelems--, xp += X_SIZEOF_FLOAT, tp++) {
        int lstatus = ncmpix_put_float_$1(xp, tp);
        if (lstatus != NC_NOERR) status = lstatus;
    }

    *xpp = (void *)xp;
    return status;
}
')dnl

PUTN_FLOAT(schar)
PUTN_FLOAT(uchar)
PUTN_FLOAT(short)
PUTN_FLOAT(int)
PUTN_FLOAT(long)
PUTN_FLOAT(double)
PUTN_FLOAT(ushort)
PUTN_FLOAT(uint)
PUTN_FLOAT(int64)
PUTN_FLOAT(uint64)

/*----< ncmpix_putn_float_float() >------------------------------------------*/
#if X_SIZEOF_FLOAT == SIZEOF_FLOAT && !defined(NO_IEEE_FLOAT)
/* optimized version */
int
ncmpix_putn_float_float(void **xpp, MPI_Offset nelems, const float *tp)
{
# ifdef WORDS_BIGENDIAN
    memcpy(*xpp, tp, nelems * X_SIZEOF_FLOAT);
# else
    ncmpii_swapn4b(*xpp, tp, nelems);
# endif
    *xpp = (void *)((char *)(*xpp) + nelems * X_SIZEOF_FLOAT);
    return NC_NOERR;
}
#elif defined(vax)
int
ncmpix_putn_float_float(void **xpp, MPI_Offset nfloats, const float *ip)
{
    const float *const end = ip + nfloats;

    while(ip < end) {
        const struct vax_single *const vsp =
             (const struct vax_single *)ip;
        struct ieee_single *const isp = (struct ieee_single *) (*xpp);

        switch(vsp->exp){
            case 0 :
                /* all vax float with zero exponent map to zero */
                *isp = min.ieee;
                break;
            case 2 :
            case 1 :
            {
                /* These will map to subnormals */
                unsigned mantissa = (vsp->mantissa1 << 16)
                         | vsp->mantissa2;
                mantissa >>= 3 - vsp->exp;
                mantissa += (1 << (20 + vsp->exp));
                isp->mant_lo_lo = mantissa;
                isp->mant_lo_hi = mantissa >> 8;
                isp->mant_hi = mantissa >> 16;
                isp->exp_lo = 0;
                isp->exp_hi = 0;
            }
                break;
            case 0xff : /* max.s.exp */
                if ( vsp->mantissa2 == max.s.mantissa2 &&
                     vsp->mantissa1 == max.s.mantissa1) {
                    /* map largest vax float to ieee infinity */
                    *isp = max.ieee;
                    break;
                } /* else, fall thru */
            default :
            {
                unsigned exp = vsp->exp - VAX_SNG_BIAS + IEEE_SNG_BIAS;
                isp->exp_hi = exp >> 1;
                isp->exp_lo = exp;
                isp->mant_lo_lo = vsp->mantissa2;
                isp->mant_lo_hi = vsp->mantissa2 >> 8;
                isp->mant_hi = vsp->mantissa1;
            }
        }

        isp->sign = vsp->sign;
    
        ip++;
        *xpp = (char *)(*xpp) + X_SIZEOF_FLOAT;
    }
    return NC_NOERR;
}
#else

static int
ncmpix_put_float_float(void *xp, const float *ip)
{
    put_ix_float(xp, ip);
#ifdef NO_IEEE_FLOAT
    if (*ip > X_FLOAT_MAX || *ip < X_FLOAT_MIN)
        return NC_ERANGE;
#endif
    return NC_NOERR;
}

PUTN_FLOAT(float)
#endif


