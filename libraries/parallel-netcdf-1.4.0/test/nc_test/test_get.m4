dnl This is m4 source.
dnl Process using m4 to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
/* Do not edit this file. It is produced from the corresponding .m4 source */
dnl
/*
 *  Copyright (C) 2003, Northwestern University and Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */
/* $Id: test_get.m4 1484 2013-11-02 14:57:28Z wkliao $ */


undefine(`index')dnl
dnl dnl dnl
dnl
dnl Macros
dnl
dnl dnl dnl
dnl
dnl Upcase(str)
dnl
define(`Upcase',dnl
`dnl
translit($1, abcdefghijklmnopqrstuvwxyz, ABCDEFGHIJKLMNOPQRSTUVWXYZ)')dnl
dnl dnl dnl
dnl
dnl NCT_ITYPE(type)
dnl
define(`NCT_ITYPE', ``NCT_'Upcase($1)')dnl
dnl

#include "tests.h"

dnl TEST_NC_GET_VAR1(TYPE)
dnl
define(`TEST_NC_GET_VAR1',dnl
`dnl
int
test_ncmpi_get_var1_$1(void)
{
    int ncid;
    int i;
    int j;
    int err;
    int nok = 0;      /* count of valid comparisons */
    MPI_Offset index[MAX_RANK];
    int canConvert;     /* Both text or both numeric */
    double expect;
    $1 value;

    err = ncmpi_open(comm, testfile, NC_NOWRITE, info, &ncid);
    IF (err != NC_NOERR)
        error("ncmpi_open: %s", ncmpi_strerror(err));

    ncmpi_begin_indep_data(ncid);
    for (i = 0; i < numVars; i++) {
        canConvert = (var_type[i] == NC_CHAR) == (NCT_ITYPE($1) == NCT_TEXT);
        for (j = 0; j < var_rank[i]; j++)
            index[j] = 0;

        /* check if pnetcdf can detect a bad file ID */
        err = ncmpi_get_var1_$1(BAD_ID, i, index, &value);
        IF (err != NC_EBADID)
            error("bad ncid: err = %d", err);

        /* check if pnetcdf can detect a bad variable ID */
        err = ncmpi_get_var1_$1(ncid, BAD_VARID, index, &value);
        IF (err != NC_ENOTVAR)
            error("bad var id: err = %d", err);

        /* check if pnetcdf can detect out of boundary requests */
        for (j = 0; j < var_rank[i]; j++) {
            index[j] = var_shape[i][j]; /* make an out-of-boundary starts[] */
            err = ncmpi_get_var1_$1(ncid, i, index, &value);
            if (!canConvert) {
                IF (err != NC_ECHAR)
                    error("conversion: err = %d", err);
            } else IF (err != NC_EINVALCOORDS)
                error("bad index: err = %d", err);
            index[j] = 0;
        }

        /* check if the contents are supposed to be */
        for (j = 0; j < var_nels[i]; j++) {
            err = toMixedBase(j, var_rank[i], var_shape[i], index);
            IF (err != NC_NOERR)
                error("error in toMixedBase 1");

            /* when file is created the variable contents are generated by
             * hash functions */
            expect = hash4( var_type[i], var_rank[i], index, NCT_ITYPE($1) );

            if (var_rank[i] == 0 && i%2 )
                /* this var has no dim, a scalar variable */
                err = ncmpi_get_var1_$1(ncid, i, NULL, &value);
            else
                err = ncmpi_get_var1_$1(ncid, i, index, &value);

            if (canConvert) {
                if (inRange3(expect,var_type[i], NCT_ITYPE($1))) {
                    if (expect >= $1_min && expect <= $1_max) {
                        IF (err != NC_NOERR) {
                            error("%s", ncmpi_strerror(err));
                        } else {
                            IF (!equal(value,expect,var_type[i],NCT_ITYPE($1))) {
                                error("expected: %G, got: %G", expect, (double) value);
                            } else {
                                nok++;
                            }
                        }
                    } else {
                        IF (err != NC_ERANGE)
                            error("Range error: err = %d", err);
                    }
                } else {
                    IF (err != NC_NOERR && err != NC_ERANGE)
                        error("OK or Range error: err = %d", err);
                }
            } else {
                IF (err != NC_ECHAR)
                    error("wrong type: err = %d", err);
            }
        }
    }
    ncmpi_end_indep_data(ncid);
    err = ncmpi_close(ncid);
    IF (err != NC_NOERR)
        error("ncmpi_close: %s", ncmpi_strerror(err));
    return nok;
}
')dnl

TEST_NC_GET_VAR1(text)
TEST_NC_GET_VAR1(schar)
TEST_NC_GET_VAR1(uchar)
TEST_NC_GET_VAR1(short)
TEST_NC_GET_VAR1(int)
TEST_NC_GET_VAR1(long)
TEST_NC_GET_VAR1(float)
TEST_NC_GET_VAR1(double)
TEST_NC_GET_VAR1(ushort)
TEST_NC_GET_VAR1(uint)
TEST_NC_GET_VAR1(longlong)
TEST_NC_GET_VAR1(ulonglong)


dnl TEST_NC_GET_VAR(TYPE)
dnl
define(`TEST_NC_GET_VAR',dnl
`dnl
int
test_ncmpi_get_var_$1(void)
{
    int ncid;
    int i;
    int j;
    int err;
    int allInExtRange;  /* all values within range of external data type */
    int allInIntRange;  /* all values within range of internal data type */
    int nels;
    int nok = 0;      /* count of valid comparisons */
    MPI_Offset index[MAX_RANK];
    int canConvert;     /* Both text or both numeric */
    double expect[MAX_NELS];
    $1 value[MAX_NELS];

    err = ncmpi_open(comm, testfile, NC_NOWRITE, info, &ncid);
    IF (err != NC_NOERR)
        error("ncmpi_open: %s", ncmpi_strerror(err));

    for (i = 0; i < numVars; i++) {
        canConvert = (var_type[i] == NC_CHAR) == (NCT_ITYPE($1) == NCT_TEXT);
        assert(var_rank[i] <= MAX_RANK);
        assert(var_nels[i] <= MAX_NELS);

        /* check if pnetcdf can detect a bad file ID */
        err = ncmpi_get_var_$1_all(BAD_ID, i, value);
        IF (err != NC_EBADID)
            error("bad ncid: err = %d", err);

        /* check if pnetcdf can detect a bad variable ID */
        err = ncmpi_get_var_$1_all(ncid, BAD_VARID, value);
        IF (err != NC_ENOTVAR)
            error("bad var id: err = %d", err);

        /* check if the contents are supposed to be */
        nels = 1;
        for (j = 0; j < var_rank[i]; j++) {
            nels *= var_shape[i][j];
        }
        allInExtRange = allInIntRange = 1;
        for (j = 0; j < nels; j++) {
            err = toMixedBase(j, var_rank[i], var_shape[i], index);
            IF (err != NC_NOERR)
                error("error in toMixedBase 1");
            expect[j] = hash4(var_type[i], var_rank[i], index, NCT_ITYPE($1));
            if (inRange3(expect[j],var_type[i], NCT_ITYPE($1))) {
                allInIntRange = allInIntRange && expect[j] >= $1_min
                            && expect[j] <= $1_max;
            } else {
                allInExtRange = 0;
            }
        }
        err = ncmpi_get_var_$1_all(ncid, i, value);
        if (canConvert) {
            if (allInExtRange) {
                if (allInIntRange) {
                    IF (err != NC_NOERR)
                        error("%s", ncmpi_strerror(err));
                } else {
                    IF (err != NC_ERANGE)
                        error("Range error: err = %d", err);
                }
            } else {
                IF (err != NC_NOERR && err != NC_ERANGE)
                    error("OK or Range error: err = %d", err);
            }
            for (j = 0; j < nels; j++) {
                if (inRange3(expect[j],var_type[i],NCT_ITYPE($1))
                        && expect[j] >= $1_min && expect[j] <= $1_max) {
                    IF (!equal(value[j],expect[j],var_type[i],NCT_ITYPE($1))){
                        error("value read not that expected");
                        if (verbose) {
                            error("\n");
                            error("varid: %d, ", i);
                            error("var_name: %s, ", var_name[i]);
                            error("var_type: %s, ", s_nc_type(var_type[i]));
                            error("element number: %d, ", j);
                            error("expect: %g, ", expect[j]);
                            error("got: %g", (double) value[j]);
                        }
                    } else {
                        nok++;
                    }
                }
            }
        } else {
            IF (nels > 0 && err != NC_ECHAR)
                error("wrong type: err = %d", err);
        }
    }
    err = ncmpi_close(ncid);
    IF (err != NC_NOERR)
        error("ncmpi_close: %s", ncmpi_strerror(err));
    return nok;
}
')dnl

TEST_NC_GET_VAR(text)
TEST_NC_GET_VAR(uchar)
TEST_NC_GET_VAR(schar)
TEST_NC_GET_VAR(short)
TEST_NC_GET_VAR(int)
TEST_NC_GET_VAR(long)
TEST_NC_GET_VAR(float)
TEST_NC_GET_VAR(double)
TEST_NC_GET_VAR(ushort)
TEST_NC_GET_VAR(uint)
TEST_NC_GET_VAR(longlong)
TEST_NC_GET_VAR(ulonglong)


dnl TEST_NC_GET_VARA(TYPE)
dnl
define(`TEST_NC_GET_VARA',dnl
`dnl
int
test_ncmpi_get_vara_$1(void)
{
    int ncid;
    int d;
    int i;
    int j;
    int k;
    int err;
    int allInExtRange;        /* all values within external range? */
    int allInIntRange;        /* all values within internal range? */
    int nels;
    int nslabs;
    int nok = 0;      /* count of valid comparisons */
    MPI_Offset start[MAX_RANK];
    MPI_Offset edge[MAX_RANK];
    MPI_Offset index[MAX_RANK];
    MPI_Offset mid[MAX_RANK];
    int canConvert;     /* Both text or both numeric */
    double expect[MAX_NELS];
    $1 value[MAX_NELS];

    err = ncmpi_open(comm, testfile, NC_NOWRITE, info, &ncid);
    IF (err != NC_NOERR)
        error("ncmpi_open: %s", ncmpi_strerror(err));
    for (i = 0; i < numVars; i++) {
        canConvert = (var_type[i] == NC_CHAR) == (NCT_ITYPE($1) == NCT_TEXT);
        assert(var_rank[i] <= MAX_RANK);
        assert(var_nels[i] <= MAX_NELS);
        for (j = 0; j < var_rank[i]; j++) {
            start[j] = 0;
            edge[j] = 1;
        }
        err = ncmpi_get_vara_$1_all(BAD_ID, i, start, edge, value);
        IF (err != NC_EBADID)
            error("bad ncid: err = %d", err);
        err = ncmpi_get_vara_$1_all(ncid, BAD_VARID, start, edge, value);
        IF (err != NC_ENOTVAR)
            error("bad var id: err = %d", err);
        for (j = 0; j < var_rank[i]; j++) {
            start[j] = var_shape[i][j]; /* causes NC_EINVALCOORDS */
            err = ncmpi_get_vara_$1_all(ncid, i, start, edge, value);
            IF (canConvert && err != NC_EINVALCOORDS)
                error("bad index: err = %d", err);
            start[j] = 0;
            edge[j] = var_shape[i][j] + 1; /* causes NC_EEDGE */
            err = ncmpi_get_vara_$1_all(ncid, i, start, edge, value);
            IF (canConvert && err != NC_EEDGE)
                error("bad edge: err = %d", err);
            edge[j] = 1;
        }
        /* Check non-scalars for correct error returned even when */
        /* there is nothing to get (edge[j]==0) */
        if (var_rank[i] > 0) {
            for (j = 0; j < var_rank[i]; j++) {
                edge[j] = 0;
            }
            err = ncmpi_get_vara_$1_all(BAD_ID, i, start, edge, value);
            IF (err != NC_EBADID) 
                error("bad ncid: err = %d", err);
            err = ncmpi_get_vara_$1_all(ncid, BAD_VARID, start, edge, value);
            IF (err != NC_ENOTVAR) 
                error("bad var id: err = %d", err);
            for (j = 0; j < var_rank[i]; j++) {
                if (var_dimid[i][j] > 0) {                /* skip record dim */
                    start[j] = var_shape[i][j];
                    err = ncmpi_get_vara_$1_all(ncid, i, start, edge, value);
                    IF (canConvert && err != NC_EINVALCOORDS)
                        error("bad start: err = %d", err);
                    start[j] = 0;
                }
            }
            err = ncmpi_get_vara_$1_all(ncid, i, start, edge, value);
            if (canConvert) {
                IF (err != NC_NOERR) 
                    error("%s", ncmpi_strerror(err));
            } else {
                IF (err != NC_ECHAR)
                    error("wrong type: err = %d", err);
            }
            for (j = 0; j < var_rank[i]; j++) {
                edge[j] = 1;
            }
        }
        /* Choose a random point dividing each dim into 2 parts */
        /* get 2^rank (nslabs) slabs so defined */
        nslabs = 1;
        for (j = 0; j < var_rank[i]; j++) {
            mid[j] = roll( var_shape[i][j] );
            nslabs *= 2;
        }
        /* bits of k determine whether to get lower or upper part of dim */
        for (k = 0; k < nslabs; k++) {
            nels = 1;
            for (j = 0; j < var_rank[i]; j++) {
                if ((k >> j) & 1) {
                    start[j] = 0;
                    edge[j] = mid[j];
                } else {
                    start[j] = mid[j];
                    edge[j] = var_shape[i][j] - mid[j];
                }
                nels *= edge[j];
            }
            allInExtRange = allInIntRange = 1;
            for (j = 0; j < nels; j++) {
                err = toMixedBase(j, var_rank[i], edge, index);
                IF (err != NC_NOERR)
                    error("error in toMixedBase 1");
                for (d = 0; d < var_rank[i]; d++)
                    index[d] += start[d];
                expect[j] = hash4(var_type[i], var_rank[i], index, NCT_ITYPE($1));
                if (inRange3(expect[j],var_type[i], NCT_ITYPE($1))) {
                    allInIntRange = allInIntRange && expect[j] >= $1_min
                                && expect[j] <= $1_max;
                } else {
                    allInExtRange = 0;
                }
            }
            if (var_rank[i] == 0 && i%2)
                err = ncmpi_get_vara_$1_all(ncid, i, NULL, NULL, value);
            else
                err = ncmpi_get_vara_$1_all(ncid, i, start, edge, value);
            if (canConvert) {
                if (allInExtRange) {
                    if (allInIntRange) {
                        IF (err != NC_NOERR)
                            error("%s", ncmpi_strerror(err));
                    } else {
                        IF (err != NC_ERANGE)
                            error("Range error: err = %d", err);
                    }
                } else {
                    IF (err != NC_NOERR && err != NC_ERANGE)
                        error("OK or Range error: err = %d", err);
                }
                for (j = 0; j < nels; j++) {
                    if (inRange3(expect[j],var_type[i],NCT_ITYPE($1))
                            && expect[j] >= $1_min && expect[j] <= $1_max) {
                        IF (!equal(value[j],expect[j],var_type[i],NCT_ITYPE($1))){
                            error("value read not that expected");
                            if (verbose) {
                                error("\n");
                                error("varid: %d, ", i);
                                error("var_name: %s, ", var_name[i]);
                                error("var_type: %s, ", s_nc_type(var_type[i]));
                                error("element number: %d, ", j);
                                error("expect: %g, ", expect[j]);
                                error("got: %g", (double) value[j]);
                            }
                        } else {
                            nok++;
                        }
                    }
                }
            } else {
                IF (nels > 0 && err != NC_ECHAR)
                    error("wrong type: err = %d", err);
            }
        }
    }
    err = ncmpi_close(ncid);
    IF (err != NC_NOERR)
        error("ncmpi_close: %s", ncmpi_strerror(err));
    return nok;
}
')dnl

TEST_NC_GET_VARA(text)
TEST_NC_GET_VARA(uchar)
TEST_NC_GET_VARA(schar)
TEST_NC_GET_VARA(short)
TEST_NC_GET_VARA(int)
TEST_NC_GET_VARA(long)
TEST_NC_GET_VARA(float)
TEST_NC_GET_VARA(double)
TEST_NC_GET_VARA(ushort)
TEST_NC_GET_VARA(uint)
TEST_NC_GET_VARA(longlong)
TEST_NC_GET_VARA(ulonglong)


dnl TEST_NC_GET_VARS(TYPE)
dnl
define(`TEST_NC_GET_VARS',dnl
`dnl
int
test_ncmpi_get_vars_$1(void)
{
    int ncid;
    int d;
    int i;
    int j;
    int k;
    int m;
    int err;
    int allInExtRange;        /* all values within external range? */
    int allInIntRange;        /* all values within internal range? */
    int nels;
    int nslabs;
    int nstarts;        /* number of different starts */
    int nok = 0;      /* count of valid comparisons */
    MPI_Offset start[MAX_RANK];
    MPI_Offset edge[MAX_RANK];
    MPI_Offset index[MAX_RANK];
    MPI_Offset index2[MAX_RANK];
    MPI_Offset mid[MAX_RANK];
    MPI_Offset count[MAX_RANK];
    MPI_Offset sstride[MAX_RANK];
    MPI_Offset stride[MAX_RANK];
    int canConvert;     /* Both text or both numeric */
    double expect[MAX_NELS];
    $1 value[MAX_NELS];

    err = ncmpi_open(comm, testfile, NC_NOWRITE, info, &ncid);
    IF (err != NC_NOERR)
        error("ncmpi_open: %s", ncmpi_strerror(err));
    for (i = 0; i < numVars; i++) {
        canConvert = (var_type[i] == NC_CHAR) == (NCT_ITYPE($1) == NCT_TEXT);
        assert(var_rank[i] <= MAX_RANK);
        assert(var_nels[i] <= MAX_NELS);
        for (j = 0; j < var_rank[i]; j++) {
            start[j] = 0;
            edge[j] = 1;
            stride[j] = 1;
        }
        err = ncmpi_get_vars_$1_all(BAD_ID, i, start, edge, stride, value);
        IF (err != NC_EBADID)
            error("bad ncid: err = %d", err);
        err = ncmpi_get_vars_$1_all(ncid, BAD_VARID, start, edge, stride, value);
        IF (err != NC_ENOTVAR)
            error("bad var id: err = %d", err);
        for (j = 0; j < var_rank[i]; j++) {
            start[j] = var_shape[i][j];
            err = ncmpi_get_vars_$1_all(ncid, i, start, edge, stride, value);
            if (!canConvert) {
                IF (err != NC_ECHAR)
                    error("conversion: err = %d", err);
            } else {
                IF (err != NC_EINVALCOORDS)
                    error("bad index: err = %d", err);
                start[j] = 0;
                edge[j] = var_shape[i][j] + 1;
                err = ncmpi_get_vars_$1_all(ncid, i, start, edge, stride, value);
                IF (err != NC_EEDGE)
                    error("bad edge: err = %d", err);
                edge[j] = 1;
                stride[j] = 0;
                err = ncmpi_get_vars_$1_all(ncid, i, start, edge, stride, value);
                IF (err != NC_ESTRIDE)
                    error("bad stride: err = %d", err);
                stride[j] = 1;
            }
        }
        /* Choose a random point dividing each dim into 2 parts */
        /* get 2^rank (nslabs) slabs so defined */
        nslabs = 1;
        for (j = 0; j < var_rank[i]; j++) {
            mid[j] = roll( var_shape[i][j] );
            nslabs *= 2;
        }
        /* bits of k determine whether to get lower or upper part of dim */
        /* choose random stride from 1 to edge */
        for (k = 0; k < nslabs; k++) {
            nstarts = 1;
            for (j = 0; j < var_rank[i]; j++) {
                if ((k >> j) & 1) {
                    start[j] = 0;
                    edge[j] = mid[j];
                } else {
                    start[j] = mid[j];
                    edge[j] = var_shape[i][j] - mid[j];
                }
                sstride[j] = stride[j] = edge[j] > 0 ? 1+roll(edge[j]) : 1;
                nstarts *= stride[j];
            }
            for (m = 0; m < nstarts; m++) {
                err = toMixedBase(m, var_rank[i], sstride, index);
                IF (err != NC_NOERR)
                    error("error in toMixedBase");
                nels = 1;
                for (j = 0; j < var_rank[i]; j++) {
                    count[j] = 1 + (edge[j] - index[j] - 1) / stride[j];
                    nels *= count[j];
                    index[j] += start[j];
                }
                /* Random choice of forward or backward */
/* TODO
                if ( roll(2) ) {
                    for (j = 0; j < var_rank[i]; j++) {
                        index[j] += (count[j] - 1) * stride[j];
                        stride[j] = -stride[j];
                    }
                }
*/
                allInExtRange = allInIntRange = 1;
                for (j = 0; j < nels; j++) {
                    err = toMixedBase(j, var_rank[i], count, index2);
                    IF (err != NC_NOERR)
                        error("error in toMixedBase 1");
                    for (d = 0; d < var_rank[i]; d++)
                        index2[d] = index[d] + index2[d] * stride[d];
                    expect[j] = hash4(var_type[i], var_rank[i], index2, 
                        NCT_ITYPE($1));
                    if (inRange3(expect[j],var_type[i],NCT_ITYPE($1))) {
                        allInIntRange = allInIntRange && expect[j] >= $1_min
                            && expect[j] <= $1_max;
                    } else {
                        allInExtRange = 0;
                    }
                }
                if (var_rank[i] == 0 && i%2 )
                    err = ncmpi_get_vars_$1_all(ncid, i, NULL, NULL, NULL, value);
                else
                    err = ncmpi_get_vars_$1_all(ncid, i, index, count, stride, value);
                if (canConvert) {
                    if (allInExtRange) {
                        if (allInIntRange) {
                            IF (err != NC_NOERR)
                                error("%s", ncmpi_strerror(err));
                        } else {
                            IF (err != NC_ERANGE)
                                error("Range error: err = %d", err);
                        }
                    } else {
                        IF (err != NC_NOERR && err != NC_ERANGE)
                            error("OK or Range error: err = %d", err);
                    }
                    for (j = 0; j < nels; j++) {
                        if (inRange3(expect[j],var_type[i],NCT_ITYPE($1))
                                && expect[j] >= $1_min && expect[j] <= $1_max) {
                            IF (!equal(value[j],expect[j],var_type[i],
                                    NCT_ITYPE($1))){
                                error("value read not that expected");
                                if (verbose) {
                                    error("\n");
                                    error("varid: %d, ", i);
                                    error("var_name: %s, ", var_name[i]);
                                    error("var_type: %s, ", s_nc_type(var_type[i]));
                                    error("element number: %d, ", j);
                                    error("expect: %g, ", expect[j]);
                                    error("got: %g", (double) value[j]);
                                }
                            } else {
                                nok++;
                            }
                        }
                    }
                } else {
                    IF (nels > 0 && err != NC_ECHAR)
                        error("wrong type: err = %d", err);
                }
            }
        }

    }
    err = ncmpi_close(ncid);
    IF (err != NC_NOERR)
        error("ncmpi_close: %s", ncmpi_strerror(err));
    return nok;
}
')dnl

TEST_NC_GET_VARS(text)
TEST_NC_GET_VARS(uchar)
TEST_NC_GET_VARS(schar)
TEST_NC_GET_VARS(short)
TEST_NC_GET_VARS(int)
TEST_NC_GET_VARS(long)
TEST_NC_GET_VARS(float)
TEST_NC_GET_VARS(double)
TEST_NC_GET_VARS(ushort)
TEST_NC_GET_VARS(uint)
TEST_NC_GET_VARS(longlong)
TEST_NC_GET_VARS(ulonglong)


dnl TEST_NC_GET_VARM(TYPE)
dnl
define(`TEST_NC_GET_VARM',dnl
`dnl
int
test_ncmpi_get_varm_$1(void)
{
    int ncid;
    int d;
    int i;
    int j;
    int k;
    int m;
    int err;
    int allInExtRange;        /* all values within external range? */
    int allInIntRange;        /* all values within internal range? */
    int nels;
    int nslabs;
    int nstarts;        /* number of different starts */
    int nok = 0;      /* count of valid comparisons */
    MPI_Offset start[MAX_RANK];
    MPI_Offset edge[MAX_RANK];
    MPI_Offset index[MAX_RANK];
    MPI_Offset index2[MAX_RANK];
    MPI_Offset mid[MAX_RANK];
    MPI_Offset count[MAX_RANK];
    MPI_Offset sstride[MAX_RANK];
    MPI_Offset stride[MAX_RANK];
    MPI_Offset imap[MAX_RANK];
    int canConvert;     /* Both text or both numeric */
    double expect[MAX_NELS];
    $1 value[MAX_NELS];

    err = ncmpi_open(comm, testfile, NC_NOWRITE, info, &ncid);
    IF (err != NC_NOERR)
        error("ncmpi_open: %s", ncmpi_strerror(err));
    for (i = 0; i < numVars; i++) {
        canConvert = (var_type[i] == NC_CHAR) == (NCT_ITYPE($1) == NCT_TEXT);
        assert(var_rank[i] <= MAX_RANK);
        assert(var_nels[i] <= MAX_NELS);
        for (j = 0; j < var_rank[i]; j++) {
            start[j] = 0;
            edge[j] = 1;
            stride[j] = 1;
            imap[j] = 1;
        }
        err = ncmpi_get_varm_$1_all(BAD_ID, i, start, edge, stride, imap, value);
        IF (err != NC_EBADID)
            error("bad ncid: err = %d", err);
        err = ncmpi_get_varm_$1_all(ncid, BAD_VARID, start, edge, stride, imap, value);
        IF (err != NC_ENOTVAR)
            error("bad var id: err = %d", err);
        for (j = 0; j < var_rank[i]; j++) {
            start[j] = var_shape[i][j];
            err = ncmpi_get_varm_$1_all(ncid, i, start, edge, stride, imap, value);
            if (!canConvert) {
                IF (err != NC_ECHAR)
                    error("conversion: err = %d", err);
            } else {
                IF (err != NC_EINVALCOORDS)
                    error("bad index: err = %d", err);
                start[j] = 0;
                edge[j] = var_shape[i][j] + 1;
                err = ncmpi_get_varm_$1_all(ncid, i, start, edge, stride, imap, value);
                IF (err != NC_EEDGE)
                    error("bad edge: err = %d", err);
                edge[j] = 1;
                stride[j] = 0;
                err = ncmpi_get_varm_$1_all(ncid, i, start, edge, stride, imap, value);
                IF (err != NC_ESTRIDE)
                    error("bad stride: err = %d", err);
                stride[j] = 1;
            }
        }
        /* Choose a random point dividing each dim into 2 parts */
        /* get 2^rank (nslabs) slabs so defined */
        nslabs = 1;
        for (j = 0; j < var_rank[i]; j++) {
            mid[j] = roll( var_shape[i][j] );
            nslabs *= 2;
        }
        /* bits of k determine whether to get lower or upper part of dim */
        /* choose random stride from 1 to edge */
        for (k = 0; k < nslabs; k++) {
            nstarts = 1;
            for (j = 0; j < var_rank[i]; j++) {
                if ((k >> j) & 1) {
                    start[j] = 0;
                    edge[j] = mid[j];
                } else {
                    start[j] = mid[j];
                    edge[j] = var_shape[i][j] - mid[j];
                }
                sstride[j] = stride[j] = edge[j] > 0 ? 1+roll(edge[j]) : 1;
                nstarts *= stride[j];
            }
            for (m = 0; m < nstarts; m++) {
                err = toMixedBase(m, var_rank[i], sstride, index);
                IF (err != NC_NOERR)
                    error("error in toMixedBase");
                nels = 1;
                for (j = 0; j < var_rank[i]; j++) {
                    count[j] = 1 + (edge[j] - index[j] - 1) / stride[j];
                    nels *= count[j];
                    index[j] += start[j];
                }
                /* Random choice of forward or backward */
/* TODO
                if ( roll(2) ) {
                    for (j = 0; j < var_rank[i]; j++) {
                        index[j] += (count[j] - 1) * stride[j];
                        stride[j] = -stride[j];
                    }
                }
 */
                if (var_rank[i] > 0) {
                    j = var_rank[i] - 1;
                    imap[j] = 1;
                    for (; j > 0; j--)
                        imap[j-1] = imap[j] * count[j];
                }
                allInExtRange = allInIntRange = 1;
                for (j = 0; j < nels; j++) {
                    err = toMixedBase(j, var_rank[i], count, index2);
                    IF (err != NC_NOERR)
                        error("error in toMixedBase 1");
                    for (d = 0; d < var_rank[i]; d++)
                        index2[d] = index[d] + index2[d] * stride[d];
                    expect[j] = hash4(var_type[i], var_rank[i], index2,
                        NCT_ITYPE($1));
                    if (inRange3(expect[j],var_type[i],NCT_ITYPE($1))) {
                        allInIntRange = allInIntRange && expect[j] >= $1_min
                            && expect[j] <= $1_max;
                    } else {
                        allInExtRange = 0;
                    }
                }
                if (var_rank[i] == 0 && i%2 )
                    err = ncmpi_get_varm_$1_all(ncid,i,NULL,NULL,NULL,NULL,value);
                else
                    err = ncmpi_get_varm_$1_all(ncid,i,index,count,stride,imap,value);
                if (canConvert) {
                    if (allInExtRange) {
                        if (allInIntRange) {
                            IF (err != NC_NOERR)
                                error("%s", ncmpi_strerror(err));
                        } else {
                            IF (err != NC_ERANGE)
                                error("Range error: err = %d", err);
                        }
                    } else {
                        IF (err != NC_NOERR && err != NC_ERANGE)
                            error("OK or Range error: err = %d", err);
                    }
                    for (j = 0; j < nels; j++) {
                        if (inRange3(expect[j],var_type[i],NCT_ITYPE($1))
                                && expect[j] >= $1_min 
                                && expect[j] <= $1_max) {
                            IF (!equal(value[j],expect[j],var_type[i],
                                    NCT_ITYPE($1))){
                                error("value read not that expected");
                                if (verbose) {
                                    error("\n");
                                    error("varid: %d, ", i);
                                    error("var_name: %s, ", var_name[i]);
                                    error("var_type: %s, ", s_nc_type(var_type[i]));
                                    error("element number: %d, ", j);
                                    error("expect: %g, ", expect[j]);
                                    error("got: %g", (double) value[j]);
                                }
                            } else {
                                nok++;
                            }
                        }
                    }
                } else {
                    IF (nels > 0 && err != NC_ECHAR)
                        error("wrong type: err = %d", err);
                }
            }
        }
    }
    err = ncmpi_close(ncid);
    IF (err != NC_NOERR)
        error("ncmpi_close: %s", ncmpi_strerror(err));
    return nok;
}
')dnl

TEST_NC_GET_VARM(text)
TEST_NC_GET_VARM(uchar)
TEST_NC_GET_VARM(schar)
TEST_NC_GET_VARM(short)
TEST_NC_GET_VARM(int)
TEST_NC_GET_VARM(long)
TEST_NC_GET_VARM(float)
TEST_NC_GET_VARM(double)
TEST_NC_GET_VARM(ushort)
TEST_NC_GET_VARM(uint)
TEST_NC_GET_VARM(longlong)
TEST_NC_GET_VARM(ulonglong)


dnl TEST_NC_GET_ATT(TYPE)
dnl
define(`TEST_NC_GET_ATT',dnl
`dnl
int
test_ncmpi_get_att_$1(void)
{
    int ncid;
    int i;
    int j;
    MPI_Offset k;
    int err;
    int allInExtRange;
    int allInIntRange;
    int canConvert;     /* Both text or both numeric */
    int nok = 0;        /* count of valid comparisons */
    double expect[MAX_NELS];
    $1 value[MAX_NELS];

    err = ncmpi_open(comm, testfile, NC_NOWRITE, info, &ncid);
    IF (err != NC_NOERR) 
        error("ncmpi_open: %s", ncmpi_strerror(err));

    for (i = -1; i < numVars; i++) {
        for (j = 0; j < NATTS(i); j++) {
            canConvert = (ATT_TYPE(i,j) == NC_CHAR) == (NCT_ITYPE($1) == NCT_TEXT);
            err = ncmpi_get_att_$1(BAD_ID, i, ATT_NAME(i,j), value);
            IF (err != NC_EBADID) 
                error("bad ncid: err = %d", err);
            err = ncmpi_get_att_$1(ncid, BAD_VARID, ATT_NAME(i,j), value);
            IF (err != NC_ENOTVAR) 
                error("bad var id: err = %d", err);
            err = ncmpi_get_att_$1(ncid, i, "noSuch", value);
            IF (err != NC_ENOTATT) 
                error("Bad attribute name: err = %d", err);
            allInExtRange = allInIntRange = 1;
            for (k = 0; k < ATT_LEN(i,j); k++) {
                expect[k] = hash4(ATT_TYPE(i,j), -1, &k, NCT_ITYPE($1));
                if (inRange3(expect[k],ATT_TYPE(i,j),NCT_ITYPE($1))) {
                    allInIntRange = allInIntRange && expect[k] >= $1_min
                                && expect[k] <= $1_max;
                } else {
                    allInExtRange = 0;
                }
            }
            err = ncmpi_get_att_$1(ncid, i, ATT_NAME(i,j), value);
            if (canConvert || ATT_LEN(i,j) == 0) {
                if (allInExtRange) {
                    if (allInIntRange) {
                        IF (err != NC_NOERR)
                            error("%s", ncmpi_strerror(err));
                    } else {
                        IF (err != NC_ERANGE)
                            error("Range error: err = %d", err);
                    }
                } else {
                    IF (err != NC_NOERR && err != NC_ERANGE)
                        error("OK or Range error: err = %d", err);
                }
                for (k = 0; k < ATT_LEN(i,j); k++) {
                    if (inRange3(expect[k],ATT_TYPE(i,j),NCT_ITYPE($1))
                            && expect[k] >= $1_min && expect[k] <= $1_max) {
                        IF (!equal(value[k],expect[k],ATT_TYPE(i,j),
                                NCT_ITYPE($1))){
                            error("value read not that expected");
                            if (verbose) {
                                error("\n");
                                error("varid: %d, ", i);
                                error("att_name: %s, ", ATT_NAME(i,j));
                                error("var_type: %s, ", s_nc_type(var_type[i]));
                                error("element number: %d, ", k);
                                error("expect: %g, ", expect[k]);
                                error("got: %g", (double) value[k]);
                            }
                        } else {
                            nok++;
                        }
                    }
                }
            } else {
                IF (err != NC_ECHAR)
                    error("wrong type: err = %d", err);
            }
        }
    }

    err = ncmpi_close(ncid);
    IF (err != NC_NOERR)
        error("ncmpi_close: %s", ncmpi_strerror(err));
    return nok;
}
')dnl

TEST_NC_GET_ATT(text)
TEST_NC_GET_ATT(uchar)
TEST_NC_GET_ATT(schar)
TEST_NC_GET_ATT(short)
TEST_NC_GET_ATT(int)
TEST_NC_GET_ATT(long)
TEST_NC_GET_ATT(float)
TEST_NC_GET_ATT(double)
TEST_NC_GET_ATT(ushort)
TEST_NC_GET_ATT(uint)
TEST_NC_GET_ATT(longlong)
TEST_NC_GET_ATT(ulonglong)

