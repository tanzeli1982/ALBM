!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id: api.f90 1468 2013-10-26 16:53:18Z wkliao $
!

!
! miscellaneous subroutines:
!
INTERFACE
    character*80   FUNCTION nfmpi_inq_libvers()
    END FUNCTION   nfmpi_inq_libvers

    character*80   FUNCTION nfmpi_strerror(ncerr)
                            INTEGER,                       INTENT(IN)  :: ncerr
    END FUNCTION   nfmpi_strerror

    logical        FUNCTION nfmpi_issyserr(ncerr)
                            INTEGER,                       INTENT(IN)  :: ncerr
    END FUNCTION   nfmpi_issyserr
!
! control subroutines:
!
    INTEGER        FUNCTION nfmpi_create(mpi_comm, path, cmode, mpi_info, ncid)
                            INTEGER,                       INTENT(IN)  :: mpi_comm
                            CHARACTER(len=*),              INTENT(IN)  :: path
                            INTEGER,                       INTENT(IN)  :: cmode
                            INTEGER,                       INTENT(IN)  :: mpi_info
                            INTEGER,                       INTENT(OUT) :: ncid
    END FUNCTION   nfmpi_create

    INTEGER        FUNCTION nfmpi_open(mpi_comm, path, mode, mpi_info, ncid)
                            INTEGER,                       INTENT(IN)  :: mpi_comm
                            CHARACTER(len=*),              INTENT(IN)  :: path
                            INTEGER,                       INTENT(IN)  :: mode
                            INTEGER,                       INTENT(IN)  :: mpi_info
                            INTEGER,                       INTENT(OUT) :: ncid
    END FUNCTION   nfmpi_open

    INTEGER        FUNCTION nfmpi_inq_format(ncid, format)
                            INTEGER,                       INTENT(IN)   :: ncid
                            INTEGER,                       INTENT(OUT)  :: format
    END FUNCTION   nfmpi_inq_format

    INTEGER        FUNCTION nfmpi_inq_file_format(path, format)
                            CHARACTER(len=*),              INTENT(IN)  :: path
                            INTEGER,                       INTENT(OUT)  :: format
    END FUNCTION   nfmpi_inq_file_format

    INTEGER        FUNCTION nfmpi_get_file_info(ncid, mpi_info)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(OUT) :: mpi_info
    END FUNCTION   nfmpi_get_file_info

    INTEGER        FUNCTION nfmpi_delete(path, mpi_info)
                            CHARACTER(len=*),              INTENT(IN)  :: path
                            INTEGER,                       INTENT(IN)  :: mpi_info
    END FUNCTION   nfmpi_delete

    INTEGER        FUNCTION nfmpi_enddef(ncid)
                            INTEGER,                       INTENT(IN)  :: ncid
    END FUNCTION   nfmpi_enddef

    INTEGER        FUNCTION nfmpi_redef(ncid)
                            INTEGER,                       INTENT(IN)  :: ncid
    END FUNCTION   nfmpi_redef

    INTEGER        FUNCTION nfmpi_set_default_format(fillmode, old_mode)
                            INTEGER,                       INTENT(IN)  :: fillmode
                            INTEGER,                       INTENT(IN)  :: old_mode
    END FUNCTION   nfmpi_set_default_format

    INTEGER        FUNCTION nfmpi_sync(ncid)
                            INTEGER,                       INTENT(IN)  :: ncid
    END FUNCTION   nfmpi_sync

    INTEGER        FUNCTION nfmpi_abort(ncid)
                            INTEGER,                       INTENT(IN)  :: ncid
    END FUNCTION   nfmpi_abort

    INTEGER        FUNCTION nfmpi_close(ncid)
                            INTEGER,                       INTENT(IN)  :: ncid
    END FUNCTION   nfmpi_close

    INTEGER        FUNCTION nfmpi_set_fill(ncid, fillmode, old_mode)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: fillmode
                            INTEGER,                       INTENT(OUT) :: old_mode
    END FUNCTION   nfmpi_set_fill

!
! general inquiry subroutines:
!

    INTEGER        FUNCTION nfmpi_inq(ncid, ndims, nvars, ngatts, unlimdimid)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(OUT) :: ndims
                            INTEGER,                       INTENT(OUT) :: nvars
                            INTEGER,                       INTENT(OUT) :: ngatts
                            INTEGER,                       INTENT(OUT) :: unlimdimid
    END FUNCTION   nfmpi_inq

    INTEGER        FUNCTION nfmpi_inq_ndims(ncid, ndims)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(OUT) :: ndims
    END FUNCTION   nfmpi_inq_ndims

    INTEGER        FUNCTION nfmpi_inq_nvars(ncid, nvars)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(OUT) :: nvars
    END FUNCTION   nfmpi_inq_nvars

    INTEGER        FUNCTION nfmpi_inq_natts(ncid, ngatts)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(OUT) :: ngatts
    END FUNCTION   nfmpi_inq_natts

    INTEGER        FUNCTION nfmpi_inq_unlimdim(ncid, unlimdimid)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(OUT) :: unlimdimid
    END FUNCTION   nfmpi_inq_unlimdim

!
! dimension subroutines:
!

    INTEGER        FUNCTION nfmpi_def_dim(ncid, name, len, dimid)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: len
                            INTEGER,                       INTENT(OUT) :: dimid
    END FUNCTION   nfmpi_def_dim

    INTEGER        FUNCTION nfmpi_inq_dimid(ncid, name, dimid)
                            INTEGER,                       INTENT(IN)  :: ncid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(OUT) :: dimid
    END FUNCTION   nfmpi_inq_dimid

    INTEGER        FUNCTION nfmpi_inq_dim(ncid, dimid, name, len)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: dimid
                            CHARACTER(len=*),              INTENT(OUT) :: name
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(OUT) :: len
    END FUNCTION   nfmpi_inq_dim

    INTEGER        FUNCTION nfmpi_inq_dimname(ncid, dimid, name)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: dimid
                            CHARACTER(len=*),              INTENT(OUT) :: name
    END FUNCTION   nfmpi_inq_dimname

    INTEGER        FUNCTION nfmpi_inq_dimlen(ncid, dimid, len)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: dimid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(OUT) :: len
    END FUNCTION   nfmpi_inq_dimlen

    INTEGER        FUNCTION nfmpi_rename_dim(ncid, dimid, name)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: dimid
                            CHARACTER(len=*),              INTENT(IN)  :: name
    END FUNCTION   nfmpi_rename_dim

!
! general attribute subroutines:
!

    INTEGER        FUNCTION nfmpi_inq_att(ncid, varid, name, xtype, len)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(OUT) :: xtype
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(OUT) :: len
    END FUNCTION   nfmpi_inq_att

    INTEGER        FUNCTION nfmpi_inq_attid(ncid, varid, name, attid)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(OUT) :: attid
    END FUNCTION   nfmpi_inq_attid

    INTEGER        FUNCTION nfmpi_inq_atttype(ncid, varid, name, xtype)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(OUT) :: xtype
    END FUNCTION   nfmpi_inq_atttype

    INTEGER        FUNCTION nfmpi_inq_attlen(ncid, varid, name, len)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(OUT) :: len
    END FUNCTION   nfmpi_inq_attlen

    INTEGER        FUNCTION nfmpi_inq_attname(ncid, varid, attid, name)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER,                       INTENT(IN)  :: attid
                            CHARACTER(len=*),              INTENT(OUT) :: name
    END FUNCTION   nfmpi_inq_attname

    INTEGER        FUNCTION nfmpi_copy_att(ncid_in, varid_in, name, ncid_out, varid_out)
                            INTEGER,                       INTENT(IN)  :: ncid_in
                            INTEGER,                       INTENT(IN)  :: varid_in
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(IN)  :: ncid_out
                            INTEGER,                       INTENT(IN)  :: varid_out
    END FUNCTION   nfmpi_copy_att

    INTEGER        FUNCTION nfmpi_rename_att(ncid, varid, curname, newname)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: curname
                            CHARACTER(len=*),              INTENT(IN)  :: newname
    END FUNCTION   nfmpi_rename_att

    INTEGER        FUNCTION nfmpi_del_att(ncid, varid, name)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
    END FUNCTION   nfmpi_del_att

!
! attribute put/get subroutines:
!

    INTEGER        FUNCTION nfmpi_put_att_text(ncid, varid, name, len, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: len
                            CHARACTER(len=*),              INTENT(IN)  :: text
    END FUNCTION   nfmpi_put_att_text

    INTEGER        FUNCTION nfmpi_get_att_text(ncid, varid, name, text)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            CHARACTER(len=*),              INTENT(OUT) :: text
    END FUNCTION   nfmpi_get_att_text

    INTEGER        FUNCTION nfmpi_put_att_int1(ncid, varid, name, xtype, len, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(IN)  :: xtype
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: len
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
    END FUNCTION   nfmpi_put_att_int1

    INTEGER        FUNCTION nfmpi_get_att_int1(ncid, varid, name, i1vals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
    END FUNCTION   nfmpi_get_att_int1

    INTEGER        FUNCTION nfmpi_put_att_int2(ncid, varid, name, xtype, len, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(IN)  :: xtype
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: len
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
    END FUNCTION   nfmpi_put_att_int2

    INTEGER        FUNCTION nfmpi_get_att_int2(ncid, varid, name, i2vals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
    END FUNCTION   nfmpi_get_att_int2

    INTEGER        FUNCTION nfmpi_put_att_int(ncid, varid, name, xtype, len, IVALS)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(IN)  :: xtype
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: len
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
    END FUNCTION   nfmpi_put_att_int

    INTEGER        FUNCTION nfmpi_get_att_int(ncid, varid, name, IVALS)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(OUT) :: ivals(*)
    END FUNCTION   nfmpi_get_att_int

    INTEGER        FUNCTION nfmpi_put_att_real(ncid, varid, name, xtype, len, RVALS)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(IN)  :: xtype
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: len
                            REAL,                          INTENT(INOUT) :: rvals(*)
    END FUNCTION   nfmpi_put_att_real

    INTEGER        FUNCTION nfmpi_get_att_real(ncid, varid, name, rvals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            REAL,                          INTENT(OUT) :: rvals(*)
    END FUNCTION   nfmpi_get_att_real

    INTEGER        FUNCTION nfmpi_put_att_double(ncid, varid, name, xtype, len, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(IN)  :: xtype
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: len
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
    END FUNCTION   nfmpi_put_att_double

    INTEGER        FUNCTION nfmpi_get_att_double(ncid, varid, name, dvals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
    END FUNCTION   nfmpi_get_att_double

    INTEGER        FUNCTION nfmpi_put_att_int8(ncid, varid, name, xtype, len, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(IN)  :: xtype
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: len
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
    END FUNCTION   nfmpi_put_att_int8

    INTEGER        FUNCTION nfmpi_get_att_int8(ncid, varid, name, i8vals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
    END FUNCTION   nfmpi_get_att_int8

!
! independent data mode subroutines:
!
    INTEGER        FUNCTION nfmpi_begin_indep_data(ncid)
                            INTEGER,                       INTENT(IN)  :: ncid
    END FUNCTION   nfmpi_begin_indep_data

    INTEGER        FUNCTION nfmpi_end_indep_data(ncid)
                            INTEGER,                       INTENT(IN)  :: ncid
    END FUNCTION   nfmpi_end_indep_data
!
! general variable subroutines:
!

    INTEGER        FUNCTION nfmpi_def_var(ncid, name, datatype, ndims, dimids, varid)
                            INTEGER,                       INTENT(IN)  :: ncid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(IN)  :: datatype
                            INTEGER,                       INTENT(IN)  :: ndims
                            INTEGER,                       INTENT(IN)  :: dimids(ndims)
                            INTEGER,                       INTENT(OUT) :: varid
    END FUNCTION   nfmpi_def_var

    INTEGER        FUNCTION nfmpi_inq_var(ncid, varid, name, datatype, ndims, dimids, natts)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(OUT) :: name
                            INTEGER,                       INTENT(OUT) :: datatype
                            INTEGER,                       INTENT(OUT) :: ndims
                            INTEGER,                       INTENT(OUT) :: dimids(*)
                            INTEGER,                       INTENT(OUT) :: natts
    END FUNCTION   nfmpi_inq_var

    INTEGER        FUNCTION nfmpi_inq_varid(ncid, name, varid)
                            INTEGER,                       INTENT(IN)  :: ncid
                            CHARACTER(len=*),              INTENT(IN)  :: name
                            INTEGER,                       INTENT(OUT) :: varid
    END FUNCTION   nfmpi_inq_varid

    INTEGER        FUNCTION nfmpi_inq_varname(ncid, varid, name)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(OUT) :: name
    END FUNCTION   nfmpi_inq_varname

    INTEGER        FUNCTION nfmpi_inq_vartype(ncid, varid, xtype)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER,                       INTENT(OUT) :: xtype
    END FUNCTION   nfmpi_inq_vartype

    INTEGER        FUNCTION nfmpi_inq_varndims(ncid, varid, ndims)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER,                       INTENT(OUT) :: ndims
    END FUNCTION   nfmpi_inq_varndims

    INTEGER        FUNCTION nfmpi_inq_vardimid(ncid, varid, dimids)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER,                       INTENT(OUT) :: dimids(*)
    END FUNCTION   nfmpi_inq_vardimid

    INTEGER        FUNCTION nfmpi_inq_varnatts(ncid, varid, natts)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER,                       INTENT(OUT) :: natts
    END FUNCTION   nfmpi_inq_varnatts

    INTEGER        FUNCTION nfmpi_rename_var(ncid, varid, name)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: name
    END FUNCTION   nfmpi_rename_var

!
! entire variable put/get subroutines:
!

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_put_var(ncid, varid, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_put_var

!   INTEGER        FUNCTION nfmpi_get_var(ncid, varid, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_get_var

!   INTEGER        FUNCTION nfmpi_get_var_all (ncid, varid, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_get_var_all

    INTEGER        FUNCTION nfmpi_put_var_text(ncid, varid, text)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: text
    END FUNCTION   nfmpi_put_var_text

    INTEGER        FUNCTION nfmpi_get_var_text(ncid, varid, text)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(OUT) :: text
    END FUNCTION   nfmpi_get_var_text

    INTEGER        FUNCTION nfmpi_get_var_text_all(ncid, varid, text)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(OUT) :: text
    END FUNCTION   nfmpi_get_var_text_all

    INTEGER        FUNCTION nfmpi_put_var_int1(ncid, varid, i1vals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
    END FUNCTION   nfmpi_put_var_int1

    INTEGER        FUNCTION nfmpi_get_var_int1(ncid, varid, i1vals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
    END FUNCTION   nfmpi_get_var_int1

    INTEGER        FUNCTION nfmpi_get_var_int1_all(ncid, varid, i1vals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
    END FUNCTION   nfmpi_get_var_int1_all

    INTEGER        FUNCTION nfmpi_put_var_int2(ncid, varid, i2vals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
    END FUNCTION   nfmpi_put_var_int2

    INTEGER        FUNCTION nfmpi_get_var_int2(ncid, varid, i2vals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
    END FUNCTION   nfmpi_get_var_int2

    INTEGER        FUNCTION nfmpi_get_var_int2_all(ncid, varid, i2vals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
    END FUNCTION   nfmpi_get_var_int2_all

    INTEGER        FUNCTION nfmpi_put_var_int(ncid, varid, ivals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
    END FUNCTION   nfmpi_put_var_int

    INTEGER        FUNCTION nfmpi_get_var_int(ncid, varid, ivals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER,                       INTENT(OUT) :: ivals(*)
    END FUNCTION   nfmpi_get_var_int

    INTEGER        FUNCTION nfmpi_get_var_int_all(ncid, varid, ivals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER,                       INTENT(OUT) :: ivals(*)
    END FUNCTION   nfmpi_get_var_int_all

    INTEGER        FUNCTION nfmpi_put_var_real(ncid, varid, rvals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            REAL,                          INTENT(INOUT) :: rvals(*)
    END FUNCTION   nfmpi_put_var_real

    INTEGER        FUNCTION nfmpi_get_var_real(ncid, varid, rvals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            REAL,                          INTENT(OUT) :: rvals(*)
    END FUNCTION   nfmpi_get_var_real

    INTEGER        FUNCTION nfmpi_get_var_real_all(ncid, varid, rvals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            REAL,                          INTENT(OUT) :: rvals(*)
    END FUNCTION   nfmpi_get_var_real_all

    INTEGER        FUNCTION nfmpi_put_var_double(ncid, varid, dvals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
    END FUNCTION   nfmpi_put_var_double

    INTEGER        FUNCTION nfmpi_get_var_double(ncid, varid, dvals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
    END FUNCTION   nfmpi_get_var_double

    INTEGER        FUNCTION nfmpi_get_var_double_all(ncid, varid, dvals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
    END FUNCTION   nfmpi_get_var_double_all

    INTEGER        FUNCTION nfmpi_put_var_int8(ncid, varid, i8vals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
    END FUNCTION   nfmpi_put_var_int8

    INTEGER        FUNCTION nfmpi_get_var_int8(ncid, varid, i8vals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
    END FUNCTION   nfmpi_get_var_int8

    INTEGER        FUNCTION nfmpi_get_var_int8_all(ncid, varid, i8vals)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
    END FUNCTION   nfmpi_get_var_int8_all

!
! single variable put/get subroutines:
!

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_put_var1(ncid, varid, index, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_put_var1

!   INTEGER        FUNCTION nfmpi_get_var1(ncid, varid, index, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_get_var1

    INTEGER        FUNCTION nfmpi_put_var1_text(ncid, varid, index, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            CHARACTER,                     INTENT(IN)  :: text
    END FUNCTION   nfmpi_put_var1_text

    INTEGER        FUNCTION nfmpi_get_var1_text(ncid, varid, index, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            CHARACTER,                     INTENT(OUT) :: text
    END FUNCTION   nfmpi_get_var1_text

    INTEGER        FUNCTION nfmpi_put_var1_int1(ncid, varid, index, i1val)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*1,                     INTENT(IN)  :: i1val
    END FUNCTION   nfmpi_put_var1_int1

    INTEGER        FUNCTION nfmpi_get_var1_int1(ncid, varid, index, i1val)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*1,                     INTENT(OUT) :: i1val
    END FUNCTION   nfmpi_get_var1_int1

    INTEGER        FUNCTION nfmpi_put_var1_int2(ncid, varid, index, i2val)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2val
    END FUNCTION   nfmpi_put_var1_int2

    INTEGER        FUNCTION nfmpi_get_var1_int2(ncid, varid, index, i2val)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*2,                     INTENT(OUT) :: i2val
    END FUNCTION   nfmpi_get_var1_int2

    INTEGER        FUNCTION nfmpi_put_var1_int(ncid, varid, index, ival)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER,                       INTENT(INOUT) :: ival
    END FUNCTION   nfmpi_put_var1_int

    INTEGER        FUNCTION nfmpi_get_var1_int(ncid, varid, index, ival)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER,                       INTENT(OUT) :: ival
    END FUNCTION   nfmpi_get_var1_int

    INTEGER        FUNCTION nfmpi_put_var1_real(ncid, varid, index, rval)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            REAL,                          INTENT(INOUT) :: rval
    END FUNCTION   nfmpi_put_var1_real

    INTEGER        FUNCTION nfmpi_get_var1_real(ncid, varid, index, rval)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            REAL,                          INTENT(OUT) :: rval
    END FUNCTION   nfmpi_get_var1_real

    INTEGER        FUNCTION nfmpi_put_var1_double(ncid, varid, index, dval)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dval
    END FUNCTION   nfmpi_put_var1_double

    INTEGER        FUNCTION nfmpi_get_var1_double(ncid, varid, index, dval)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            DOUBLE PRECISION,              INTENT(OUT) :: dval
    END FUNCTION   nfmpi_get_var1_double

    INTEGER        FUNCTION nfmpi_put_var1_int8(ncid, varid, index, i8val)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8val
    END FUNCTION   nfmpi_put_var1_int8

    INTEGER        FUNCTION nfmpi_get_var1_int8(ncid, varid, index, i8val)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*8,                     INTENT(OUT) :: i8val
    END FUNCTION   nfmpi_get_var1_int8

!
! variable array put/get subroutines:
!

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_put_vara(ncid, varid, start, count, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_put_vara

!   INTEGER        FUNCTION nfmpi_put_vara_all(ncid, varid, start, count, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_put_vara_all

!   INTEGER        FUNCTION nfmpi_get_vara(ncid, varid, start, count, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_get_vara

!   INTEGER        FUNCTION nfmpi_get_vara_all (ncid, varid, start, count, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_get_vara_all

    INTEGER        FUNCTION nfmpi_put_vara_text(ncid, varid, start, count, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            CHARACTER(len=*),              INTENT(IN)  :: text
    END FUNCTION   nfmpi_put_vara_text

    INTEGER        FUNCTION nfmpi_put_vara_text_all(ncid, varid, start, count, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            CHARACTER(len=*),              INTENT(IN)  :: text
    END FUNCTION   nfmpi_put_vara_text_all

    INTEGER        FUNCTION nfmpi_get_vara_text(ncid, varid, start, count, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            CHARACTER(len=*),              INTENT(OUT) :: text
    END FUNCTION   nfmpi_get_vara_text

    INTEGER        FUNCTION nfmpi_get_vara_text_all(ncid, varid, start, count, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            CHARACTER(len=*),              INTENT(OUT) :: text
    END FUNCTION   nfmpi_get_vara_text_all

    INTEGER        FUNCTION nfmpi_put_vara_int1(ncid, varid, start, count, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
    END FUNCTION   nfmpi_put_vara_int1

    INTEGER        FUNCTION nfmpi_put_vara_int1_all(ncid, varid, start, count, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
    END FUNCTION   nfmpi_put_vara_int1_all

    INTEGER        FUNCTION nfmpi_get_vara_int1(ncid, varid, start, count, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
    END FUNCTION   nfmpi_get_vara_int1

    INTEGER        FUNCTION nfmpi_get_vara_int1_all(ncid, varid, start, count, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
    END FUNCTION   nfmpi_get_vara_int1_all

    INTEGER        FUNCTION nfmpi_put_vara_int2(ncid, varid, start, count, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
    END FUNCTION   nfmpi_put_vara_int2

    INTEGER        FUNCTION nfmpi_put_vara_int2_all(ncid, varid, start, count, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
    END FUNCTION   nfmpi_put_vara_int2_all

    INTEGER        FUNCTION nfmpi_get_vara_int2(ncid, varid, start, count, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
    END FUNCTION   nfmpi_get_vara_int2

    INTEGER        FUNCTION nfmpi_get_vara_int2_all(ncid, varid, start, count, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
    END FUNCTION   nfmpi_get_vara_int2_all

    INTEGER        FUNCTION nfmpi_put_vara_int(ncid, varid, start, count, ivals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
    END FUNCTION   nfmpi_put_vara_int

    INTEGER        FUNCTION nfmpi_put_vara_int_all(ncid, varid, start, count, ivals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
    END FUNCTION   nfmpi_put_vara_int_all

    INTEGER        FUNCTION nfmpi_get_vara_int(ncid, varid, start, count, ivals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER,                       INTENT(OUT) :: ivals(*)
    END FUNCTION   nfmpi_get_vara_int

    INTEGER        FUNCTION nfmpi_get_vara_int_all(ncid, varid, start, count, ivals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER,                       INTENT(OUT) :: ivals(*)
    END FUNCTION   nfmpi_get_vara_int_all

    INTEGER        FUNCTION nfmpi_put_vara_real(ncid, varid, start, count, rvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            REAL,                          INTENT(INOUT) :: rvals(*)
    END FUNCTION   nfmpi_put_vara_real

    INTEGER        FUNCTION nfmpi_put_vara_real_all(ncid, varid, start, count, rvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            REAL,                          INTENT(INOUT) :: rvals(*)
    END FUNCTION   nfmpi_put_vara_real_all

    INTEGER        FUNCTION nfmpi_get_vara_real(ncid, varid, start, count, rvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            REAL,                          INTENT(OUT) :: rvals(*)
    END FUNCTION   nfmpi_get_vara_real

    INTEGER        FUNCTION nfmpi_get_vara_real_all(ncid, varid, start, count, rvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            REAL,                          INTENT(OUT) :: rvals(*)
    END FUNCTION   nfmpi_get_vara_real_all

    INTEGER        FUNCTION nfmpi_put_vara_double(ncid, varid, start, count, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
    END FUNCTION   nfmpi_put_vara_double

    INTEGER        FUNCTION nfmpi_put_vara_double_all(ncid, varid, start, count, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
    END FUNCTION   nfmpi_put_vara_double_all

    INTEGER        FUNCTION nfmpi_get_vara_double(ncid, varid, start, count, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
    END FUNCTION   nfmpi_get_vara_double

    INTEGER        FUNCTION nfmpi_get_vara_double_all(ncid, varid, start, count, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
    END FUNCTION   nfmpi_get_vara_double_all

    INTEGER        FUNCTION nfmpi_put_vara_int8(ncid, varid, start, count, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
    END FUNCTION   nfmpi_put_vara_int8

    INTEGER        FUNCTION nfmpi_put_vara_int8_all(ncid, varid, start, count, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
    END FUNCTION   nfmpi_put_vara_int8_all

    INTEGER        FUNCTION nfmpi_get_vara_int8(ncid, varid, start, count, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
    END FUNCTION   nfmpi_get_vara_int8

    INTEGER        FUNCTION nfmpi_get_vara_int8_all(ncid, varid, start, count, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
    END FUNCTION   nfmpi_get_vara_int8_all

!
! strided variable put/get subroutines:
!

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_put_vars(ncid, varid, start, count, stride, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_put_vars

!   INTEGER        FUNCTION nfmpi_put_vars_all(ncid, varid, start, count, stride, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_put_vars_all

!   INTEGER        FUNCTION nfmpi_get_vars(ncid, varid, start, count, stride, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_get_vars

!   INTEGER        FUNCTION nfmpi_get_vars_all(ncid, varid, start, count, stride, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_get_vars_all

    INTEGER        FUNCTION nfmpi_put_vars_text(ncid, varid, start, count, stride, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            CHARACTER(len=*),              INTENT(IN)  :: text
    END FUNCTION   nfmpi_put_vars_text

    INTEGER        FUNCTION nfmpi_put_vars_text_all(ncid, varid, start, count, stride, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            CHARACTER(len=*),              INTENT(IN)  :: text
    END FUNCTION   nfmpi_put_vars_text_all

    INTEGER        FUNCTION nfmpi_get_vars_text(ncid, varid, start, count, stride, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            CHARACTER(len=*),              INTENT(OUT) :: text
    END FUNCTION   nfmpi_get_vars_text

    INTEGER        FUNCTION nfmpi_get_vars_text_all(ncid, varid, start, count, stride, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            CHARACTER(len=*),              INTENT(OUT) :: text
    END FUNCTION   nfmpi_get_vars_text_all

    INTEGER        FUNCTION nfmpi_put_vars_int1(ncid, varid, start, count, stride, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
    END FUNCTION   nfmpi_put_vars_int1

    INTEGER        FUNCTION nfmpi_put_vars_int1_all(ncid, varid, start, count, stride, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
    END FUNCTION   nfmpi_put_vars_int1_all

    INTEGER        FUNCTION nfmpi_get_vars_int1(ncid, varid, start, count, stride, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
    END FUNCTION   nfmpi_get_vars_int1

    INTEGER        FUNCTION nfmpi_get_vars_int1_all(ncid, varid, start, count, stride, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
    END FUNCTION   nfmpi_get_vars_int1_all

    INTEGER        FUNCTION nfmpi_put_vars_int2(ncid, varid, start, count, stride, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
    END FUNCTION   nfmpi_put_vars_int2

    INTEGER        FUNCTION nfmpi_put_vars_int2_all(ncid, varid, start, count, stride, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
    END FUNCTION   nfmpi_put_vars_int2_all

    INTEGER        FUNCTION nfmpi_get_vars_int2(ncid, varid, start, count, stride, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
    END FUNCTION   nfmpi_get_vars_int2

    INTEGER        FUNCTION nfmpi_get_vars_int2_all(ncid, varid, start, count, stride, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
    END FUNCTION   nfmpi_get_vars_int2_all

    INTEGER        FUNCTION nfmpi_put_vars_int(ncid, varid, start, count, stride, ivals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
    END FUNCTION   nfmpi_put_vars_int

    INTEGER        FUNCTION nfmpi_put_vars_int_all(ncid, varid, start, count, stride, ivals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
    END FUNCTION   nfmpi_put_vars_int_all

    INTEGER        FUNCTION nfmpi_get_vars_int(ncid, varid, start, count, stride, ivals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER,                       INTENT(OUT) :: ivals(*)
    END FUNCTION   nfmpi_get_vars_int

    INTEGER        FUNCTION nfmpi_get_vars_int_all(ncid, varid, start, count, stride, ivals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER,                       INTENT(OUT) :: ivals(*)
    END FUNCTION   nfmpi_get_vars_int_all

    INTEGER        FUNCTION nfmpi_put_vars_real(ncid, varid, start, count, stride, rvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            REAL,                          INTENT(INOUT) :: rvals(*)
    END FUNCTION   nfmpi_put_vars_real

    INTEGER        FUNCTION nfmpi_put_vars_real_all(ncid, varid, start, count, stride, rvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            REAL,                          INTENT(INOUT) :: rvals(*)
    END FUNCTION   nfmpi_put_vars_real_all

    INTEGER        FUNCTION nfmpi_get_vars_real(ncid, varid, start, count, stride, rvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            REAL,                          INTENT(OUT) :: rvals(*)
    END FUNCTION   nfmpi_get_vars_real

    INTEGER        FUNCTION nfmpi_get_vars_real_all(ncid, varid, start, count, stride, rvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            REAL,                          INTENT(OUT) :: rvals(*)
    END FUNCTION   nfmpi_get_vars_real_all

    INTEGER        FUNCTION nfmpi_put_vars_double(ncid, varid, start, count, stride, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
    END FUNCTION   nfmpi_put_vars_double

    INTEGER        FUNCTION nfmpi_put_vars_double_all(ncid, varid, start, count, stride, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
    END FUNCTION   nfmpi_put_vars_double_all

    INTEGER        FUNCTION nfmpi_get_vars_double(ncid, varid, start, count, stride, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
    END FUNCTION   nfmpi_get_vars_double

    INTEGER        FUNCTION nfmpi_get_vars_double_all(ncid, varid, start, count, stride, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
    END FUNCTION   nfmpi_get_vars_double_all

    INTEGER        FUNCTION nfmpi_put_vars_int8(ncid, varid, start, count, stride, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
    END FUNCTION   nfmpi_put_vars_int8

    INTEGER        FUNCTION nfmpi_put_vars_int8_all(ncid, varid, start, count, stride, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
    END FUNCTION   nfmpi_put_vars_int8_all

    INTEGER        FUNCTION nfmpi_get_vars_int8(ncid, varid, start, count, stride, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
    END FUNCTION   nfmpi_get_vars_int8

    INTEGER        FUNCTION nfmpi_get_vars_int8_all(ncid, varid, start, count, stride, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
    END FUNCTION   nfmpi_get_vars_int8_all

!
! mapped variable put/get subroutines:
!

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_put_varm(ncid, varid, start, count, stride, imap, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_put_varm

!   INTEGER        FUNCTION nfmpi_put_varm_all(ncid, varid, start, count, stride, imap, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_put_varm_all

!   INTEGER        FUNCTION nfmpi_get_varm(ncid, varid, start, count, stride, imap, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_get_varm

!   INTEGER        FUNCTION nfmpi_get_varm_all(ncid, varid, start, count, stride, imap, buf, bufcount, datatype)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!   END FUNCTION   nfmpi_get_varm_all

    INTEGER        FUNCTION nfmpi_put_varm_text(ncid, varid, start, count, stride, imap, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            CHARACTER(len=*),              INTENT(IN)  :: text
    END FUNCTION   nfmpi_put_varm_text

    INTEGER        FUNCTION nfmpi_put_varm_text_all(ncid, varid, start, count, stride, imap, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            CHARACTER(len=*),              INTENT(IN)  :: text
    END FUNCTION   nfmpi_put_varm_text_all

    INTEGER        FUNCTION nfmpi_get_varm_text(ncid, varid, start, count, stride, imap, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            CHARACTER(len=*),              INTENT(OUT) :: text
    END FUNCTION   nfmpi_get_varm_text

    INTEGER        FUNCTION nfmpi_get_varm_text_all(ncid, varid, start, count, stride, imap, text)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            CHARACTER(len=*),              INTENT(OUT) :: text
    END FUNCTION   nfmpi_get_varm_text_all

    INTEGER        FUNCTION nfmpi_put_varm_int1(ncid, varid, start, count, stride, imap, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
    END FUNCTION   nfmpi_put_varm_int1

    INTEGER        FUNCTION nfmpi_put_varm_int1_all(ncid, varid, start, count, stride, imap, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
    END FUNCTION   nfmpi_put_varm_int1_all

    INTEGER        FUNCTION nfmpi_get_varm_int1(ncid, varid, start, count, stride, imap, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
    END FUNCTION   nfmpi_get_varm_int1

    INTEGER        FUNCTION nfmpi_get_varm_int1_all(ncid, varid, start, count, stride, imap, i1vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
    END FUNCTION   nfmpi_get_varm_int1_all

    INTEGER        FUNCTION nfmpi_put_varm_int2(ncid, varid, start, count, stride, imap, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
    END FUNCTION   nfmpi_put_varm_int2

    INTEGER        FUNCTION nfmpi_put_varm_int2_all(ncid, varid, start, count, stride, imap, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
    END FUNCTION   nfmpi_put_varm_int2_all

    INTEGER        FUNCTION nfmpi_get_varm_int2(ncid, varid, start, count, stride, imap, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
    END FUNCTION   nfmpi_get_varm_int2

    INTEGER        FUNCTION nfmpi_get_varm_int2_all(ncid, varid, start, count, stride, imap, i2vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
    END FUNCTION   nfmpi_get_varm_int2_all

    INTEGER        FUNCTION nfmpi_put_varm_int (ncid, varid, start, count, stride, imap, ivals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
    END FUNCTION   nfmpi_put_varm_int

    INTEGER        FUNCTION nfmpi_put_varm_int_all(ncid, varid, start, count, stride, imap, ivals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
    END FUNCTION   nfmpi_put_varm_int_all

    INTEGER        FUNCTION nfmpi_get_varm_int (ncid, varid, start, count, stride, imap, ivals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER,                       INTENT(OUT) :: ivals(*)
    END FUNCTION   nfmpi_get_varm_int

    INTEGER        FUNCTION nfmpi_get_varm_int_all(ncid, varid, start, count, stride, imap, ivals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER,                       INTENT(OUT) :: ivals(*)
    END FUNCTION   nfmpi_get_varm_int_all

    INTEGER        FUNCTION nfmpi_put_varm_real(ncid, varid, start, count, stride, imap, rvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            REAL,                          INTENT(INOUT) :: rvals(*)
    END FUNCTION   nfmpi_put_varm_real

    INTEGER        FUNCTION nfmpi_put_varm_real_all(ncid, varid, start, count, stride, imap, rvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            REAL,                          INTENT(INOUT) :: rvals(*)
    END FUNCTION   nfmpi_put_varm_real_all

    INTEGER        FUNCTION nfmpi_get_varm_real(ncid, varid, start, count, stride, imap, rvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            REAL,                          INTENT(OUT) :: rvals(*)
    END FUNCTION   nfmpi_get_varm_real

    INTEGER        FUNCTION nfmpi_get_varm_real_all(ncid, varid, start, count, stride, imap, rvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            REAL,                          INTENT(OUT) :: rvals(*)
    END FUNCTION   nfmpi_get_varm_real_all

    INTEGER        FUNCTION nfmpi_put_varm_double(ncid, varid, start, count, stride, imap, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
    END FUNCTION   nfmpi_put_varm_double

    INTEGER        FUNCTION nfmpi_put_varm_double_all(ncid, varid, start, count, stride, imap, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
    END FUNCTION   nfmpi_put_varm_double_all

    INTEGER        FUNCTION nfmpi_get_varm_double(ncid, varid, start, count, stride, imap, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
    END FUNCTION   nfmpi_get_varm_double

    INTEGER        FUNCTION nfmpi_get_varm_double_all(ncid, varid, start, count, stride, imap, dvals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
    END FUNCTION   nfmpi_get_varm_double_all

    INTEGER        FUNCTION nfmpi_put_varm_int8(ncid, varid, start, count, stride, imap, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
    END FUNCTION   nfmpi_put_varm_int8

    INTEGER        FUNCTION nfmpi_put_varm_int8_all(ncid, varid, start, count, stride, imap, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
    END FUNCTION   nfmpi_put_varm_int8_all

    INTEGER        FUNCTION nfmpi_get_varm_int8(ncid, varid, start, count, stride, imap, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
    END FUNCTION   nfmpi_get_varm_int8

    INTEGER        FUNCTION nfmpi_get_varm_int8_all(ncid, varid, start, count, stride, imap, i8vals)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
    END FUNCTION   nfmpi_get_varm_int8_all

!
! non-blocking variable array iput/iget subroutines:
!

!
! entire variable iput/iget subroutines:
!

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_iput_var(ncid, varid, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_iput_var

!   INTEGER        FUNCTION nfmpi_iget_var(ncid, varid, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_iget_var


    INTEGER        FUNCTION nfmpi_iput_var_text(ncid, varid, text, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var_text

    INTEGER        FUNCTION nfmpi_iget_var_text(ncid, varid, text, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(OUT) :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var_text

    INTEGER        FUNCTION nfmpi_iput_var_int1(ncid, varid, i1vals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var_int1

    INTEGER        FUNCTION nfmpi_iget_var_int1(ncid, varid, i1vals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var_int1

    INTEGER        FUNCTION nfmpi_iput_var_int2(ncid, varid, i2vals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var_int2

    INTEGER        FUNCTION nfmpi_iget_var_int2(ncid, varid, i2vals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var_int2

    INTEGER        FUNCTION nfmpi_iput_var_int(ncid, varid, ivals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var_int

    INTEGER        FUNCTION nfmpi_iget_var_int(ncid, varid, ivals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER,                       INTENT(OUT) :: ivals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var_int

    INTEGER        FUNCTION nfmpi_iput_var_real(ncid, varid, rvals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            REAL,                          INTENT(INOUT) :: rvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var_real

    INTEGER        FUNCTION nfmpi_iget_var_real(ncid, varid, rvals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            REAL,                          INTENT(OUT) :: rvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var_real

    INTEGER        FUNCTION nfmpi_iput_var_double(ncid, varid, dvals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var_double

    INTEGER        FUNCTION nfmpi_iget_var_double(ncid, varid, dvals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var_double

    INTEGER        FUNCTION nfmpi_iput_var_int8(ncid, varid, i8vals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var_int8

    INTEGER        FUNCTION nfmpi_iget_var_int8(ncid, varid, i8vals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var_int8

!
! single variable iput/iget subroutines:
!

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_iput_var1(ncid, varid, index, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_iput_var1

!   INTEGER        FUNCTION nfmpi_iget_var1(ncid, varid, index, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_iget_var1

    INTEGER        FUNCTION nfmpi_iput_var1_text(ncid, varid, index, text, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            CHARACTER,                     INTENT(IN)  :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var1_text

    INTEGER        FUNCTION nfmpi_iget_var1_text(ncid, varid, index, text, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            CHARACTER,                     INTENT(OUT) :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var1_text

    INTEGER        FUNCTION nfmpi_iput_var1_int1(ncid, varid, index, i1val, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*1,                     INTENT(IN)  :: i1val
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var1_int1

    INTEGER        FUNCTION nfmpi_iget_var1_int1(ncid, varid, index, i1val, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*1,                     INTENT(OUT) :: i1val
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var1_int1

    INTEGER        FUNCTION nfmpi_iput_var1_int2(ncid, varid, index, i2val, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2val
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var1_int2

    INTEGER        FUNCTION nfmpi_iget_var1_int2(ncid, varid, index, i2val, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*2,                     INTENT(OUT) :: i2val
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var1_int2

    INTEGER        FUNCTION nfmpi_iput_var1_int(ncid, varid, index, ival, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER,                       INTENT(INOUT) :: ival
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var1_int

    INTEGER        FUNCTION nfmpi_iget_var1_int(ncid, varid, index, ival, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER,                       INTENT(OUT) :: ival
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var1_int

    INTEGER        FUNCTION nfmpi_iput_var1_real(ncid, varid, index, rval, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            REAL,                          INTENT(INOUT) :: rval
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var1_real

    INTEGER        FUNCTION nfmpi_iget_var1_real(ncid, varid, index, rval, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            REAL,                          INTENT(OUT) :: rval
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var1_real

    INTEGER        FUNCTION nfmpi_iput_var1_double(ncid, varid, index, dval, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dval
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var1_double

    INTEGER        FUNCTION nfmpi_iget_var1_double(ncid, varid, index, dval, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            DOUBLE PRECISION,              INTENT(OUT) :: dval
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var1_double

    INTEGER        FUNCTION nfmpi_iput_var1_int8(ncid, varid, index, i8val, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8val
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_var1_int8

    INTEGER        FUNCTION nfmpi_iget_var1_int8(ncid, varid, index, i8val, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*8,                     INTENT(OUT) :: i8val
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_var1_int8

!
! variable array iput/iget subroutines:
!

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_iput_vara(ncid, varid, start, count, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_iput_vara

!   INTEGER        FUNCTION nfmpi_iget_vara(ncid, varid, start, count, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_iget_vara

    INTEGER        FUNCTION nfmpi_iput_vara_text(ncid, varid, start, count, text, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            CHARACTER(len=*),              INTENT(IN)  :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vara_text

    INTEGER        FUNCTION nfmpi_iget_vara_text(ncid, varid, start, count, text, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            CHARACTER(len=*),              INTENT(OUT) :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vara_text

    INTEGER        FUNCTION nfmpi_iput_vara_int1(ncid, varid, start, count, i1vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vara_int1

    INTEGER        FUNCTION nfmpi_iget_vara_int1(ncid, varid, start, count, i1vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vara_int1

    INTEGER        FUNCTION nfmpi_iput_vara_int2(ncid, varid, start, count, i2vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vara_int2

    INTEGER        FUNCTION nfmpi_iget_vara_int2(ncid, varid, start, count, i2vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vara_int2

    INTEGER        FUNCTION nfmpi_iput_vara_int(ncid, varid, start, count, ivals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vara_int

    INTEGER        FUNCTION nfmpi_iget_vara_int(ncid, varid, start, count, ivals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER,                       INTENT(OUT) :: ivals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vara_int

    INTEGER        FUNCTION nfmpi_iput_vara_real(ncid, varid, start, count, rvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            REAL,                          INTENT(INOUT) :: rvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vara_real

    INTEGER        FUNCTION nfmpi_iget_vara_real(ncid, varid, start, count, rvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            REAL,                          INTENT(OUT) :: rvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vara_real

    INTEGER        FUNCTION nfmpi_iput_vara_double(ncid, varid, start, count, dvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vara_double

    INTEGER        FUNCTION nfmpi_iget_vara_double(ncid, varid, start, count, dvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vara_double

    INTEGER        FUNCTION nfmpi_iput_vara_int8(ncid, varid, start, count, i8vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vara_int8

    INTEGER        FUNCTION nfmpi_iget_vara_int8(ncid, varid, start, count, i8vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vara_int8

!
! strided variable iput/iget subroutines:
!

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_iput_vars(ncid, varid, start, count, stride, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_iput_vars

!   INTEGER        FUNCTION nfmpi_iget_vars(ncid, varid, start, count, stride, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_iget_vars

    INTEGER        FUNCTION nfmpi_iput_vars_text(ncid, varid, start, count, stride, text, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            CHARACTER(len=*),              INTENT(IN)  :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vars_text

    INTEGER        FUNCTION nfmpi_iget_vars_text(ncid, varid, start, count, stride, text, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            CHARACTER(len=*),              INTENT(OUT) :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vars_text

    INTEGER        FUNCTION nfmpi_iput_vars_int1(ncid, varid, start, count, stride, i1vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vars_int1

    INTEGER        FUNCTION nfmpi_iget_vars_int1(ncid, varid, start, count, stride, i1vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vars_int1

    INTEGER        FUNCTION nfmpi_iput_vars_int2(ncid, varid, start, count, stride, i2vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vars_int2

    INTEGER        FUNCTION nfmpi_iget_vars_int2(ncid, varid, start, count, stride, i2vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vars_int2

    INTEGER        FUNCTION nfmpi_iput_vars_int(ncid, varid, start, count, stride, ivals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vars_int

    INTEGER        FUNCTION nfmpi_iget_vars_int(ncid, varid, start, count, stride, ivals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER,                       INTENT(OUT) :: ivals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vars_int

    INTEGER        FUNCTION nfmpi_iput_vars_real(ncid, varid, start, count, stride, rvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            REAL,                          INTENT(INOUT) :: rvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vars_real

    INTEGER        FUNCTION nfmpi_iget_vars_real(ncid, varid, start, count, stride, rvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            REAL,                          INTENT(OUT) :: rvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vars_real

    INTEGER        FUNCTION nfmpi_iput_vars_double(ncid, varid, start, count, stride, dvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vars_double

    INTEGER        FUNCTION nfmpi_iget_vars_double(ncid, varid, start, count, stride, dvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vars_double

    INTEGER        FUNCTION nfmpi_iput_vars_int8(ncid, varid, start, count, stride, i8vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_vars_int8

    INTEGER        FUNCTION nfmpi_iget_vars_int8(ncid, varid, start, count, stride, i8vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_vars_int8

!
! mapped variable iput/iget subroutines:
!

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_iput_varm(ncid, varid, start, count, stride, imap, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_iput_varm

!   INTEGER        FUNCTION nfmpi_iget_varm(ncid, varid, start, count, stride, imap, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_iget_varm

    INTEGER        FUNCTION nfmpi_iput_varm_text(ncid, varid, start, count, stride, imap, text, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            CHARACTER(len=*),              INTENT(IN)  :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_varm_text

    INTEGER        FUNCTION nfmpi_iget_varm_text(ncid, varid, start, count, stride, imap, text, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            CHARACTER(len=*),              INTENT(OUT) :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_varm_text

    INTEGER        FUNCTION nfmpi_iput_varm_int1(ncid, varid, start, count, stride, imap, i1vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_varm_int1

    INTEGER        FUNCTION nfmpi_iget_varm_int1(ncid, varid, start, count, stride, imap, i1vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*1,                     INTENT(OUT) :: i1vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_varm_int1

    INTEGER        FUNCTION nfmpi_iput_varm_int2(ncid, varid, start, count, stride, imap, i2vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_varm_int2

    INTEGER        FUNCTION nfmpi_iget_varm_int2(ncid, varid, start, count, stride, imap, i2vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*2,                     INTENT(OUT) :: i2vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_varm_int2

    INTEGER        FUNCTION nfmpi_iput_varm_int (ncid, varid, start, count, stride, imap, ivals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_varm_int

    INTEGER        FUNCTION nfmpi_iget_varm_int (ncid, varid, start, count, stride, imap, ivals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER,                       INTENT(OUT) :: ivals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_varm_int

    INTEGER        FUNCTION nfmpi_iput_varm_real(ncid, varid, start, count, stride, imap, rvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            REAL,                          INTENT(INOUT) :: rvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_varm_real

    INTEGER        FUNCTION nfmpi_iget_varm_real(ncid, varid, start, count, stride, imap, rvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            REAL,                          INTENT(OUT) :: rvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_varm_real

    INTEGER        FUNCTION nfmpi_iput_varm_double(ncid, varid, start, count, stride, imap, dvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_varm_double

    INTEGER        FUNCTION nfmpi_iget_varm_double(ncid, varid, start, count, stride, imap, dvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_varm_double

    INTEGER        FUNCTION nfmpi_iput_varm_int8(ncid, varid, start, count, stride, imap, i8vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iput_varm_int8

    INTEGER        FUNCTION nfmpi_iget_varm_int8(ncid, varid, start, count, stride, imap, i8vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*8,                     INTENT(OUT) :: i8vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_iget_varm_int8

!
! Begin buffered put non-blocking subroutines:
!

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_bput_var(ncid, varid, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_bput_var

!   INTEGER        FUNCTION nfmpi_bget_var(ncid, varid, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_bget_var

    INTEGER        FUNCTION nfmpi_bput_var_text(ncid, varid, text, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            CHARACTER(len=*),              INTENT(IN)  :: text(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var_text

    INTEGER        FUNCTION nfmpi_bput_var_int1(ncid, varid, i1vals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var_int1

    INTEGER        FUNCTION nfmpi_bput_var_int2(ncid, varid, i2vals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var_int2

    INTEGER        FUNCTION nfmpi_bput_var_int(ncid, varid, ivals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var_int

    INTEGER        FUNCTION nfmpi_bput_var_real(ncid, varid, rvals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            REAL,                          INTENT(INOUT) :: rvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var_real

    INTEGER        FUNCTION nfmpi_bput_var_double(ncid, varid, dvals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var_double

    INTEGER        FUNCTION nfmpi_bput_var_int8(ncid, varid, i8vals, req)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var_int8

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_bput_var1(ncid, varid, start, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_bput_var1

!   INTEGER        FUNCTION nfmpi_bget_var1(ncid, varid, start, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_bget_var1

    INTEGER        FUNCTION nfmpi_bput_var1_text(ncid, varid, index, text, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            CHARACTER,                     INTENT(IN)  :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var1_text

    INTEGER        FUNCTION nfmpi_bput_var1_int1(ncid, varid, index, i1val, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*1,                     INTENT(IN)  :: i1val
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var1_int1

    INTEGER        FUNCTION nfmpi_bput_var1_int2(ncid, varid, index, i2val, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2val
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var1_int2

    INTEGER        FUNCTION nfmpi_bput_var1_int(ncid, varid, index, ival, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER,                       INTENT(INOUT) :: ival
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var1_int

    INTEGER        FUNCTION nfmpi_bput_var1_real(ncid, varid, index, rval, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            REAL,                          INTENT(INOUT) :: rval
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var1_real

    INTEGER        FUNCTION nfmpi_bput_var1_double(ncid, varid, index, dval, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dval
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var1_double

    INTEGER        FUNCTION nfmpi_bput_var1_int8(ncid, varid, index, i8val, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: index(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8val
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_var1_int8

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_bput_vara(ncid, varid, start, count, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_bput_vara

!   INTEGER        FUNCTION nfmpi_bget_vara(ncid, varid, start, count, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_bget_vara

    INTEGER        FUNCTION nfmpi_bput_vara_text(ncid, varid, start, count, text, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            CHARACTER(len=*),              INTENT(IN)  :: text(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vara_text

    INTEGER        FUNCTION nfmpi_bput_vara_int1(ncid, varid, start, count, i1vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vara_int1

    INTEGER        FUNCTION nfmpi_bput_vara_int2(ncid, varid, start, count, i2vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vara_int2

    INTEGER        FUNCTION nfmpi_bput_vara_int(ncid, varid, start, count, ivals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vara_int

    INTEGER        FUNCTION nfmpi_bput_vara_real(ncid, varid, start, count, rvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            REAL,                          INTENT(INOUT) :: rvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vara_real

    INTEGER        FUNCTION nfmpi_bput_vara_double(ncid, varid, start, count, dvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vara_double

    INTEGER        FUNCTION nfmpi_bput_vara_int8(ncid, varid, start, count, i8vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vara_int8

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_bput_vars(ncid, varid, start, count, stride, imap, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_bput_vars

!   INTEGER        FUNCTION nfmpi_bget_vars(ncid, varid, start, count, stride, imap, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_bget_vars

    INTEGER        FUNCTION nfmpi_bput_vars_text(ncid, varid, start, count, stride, text, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            CHARACTER(len=*),              INTENT(IN)  :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vars_text

    INTEGER        FUNCTION nfmpi_bput_vars_int1(ncid, varid, start, count, stride, i1vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vars_int1

    INTEGER        FUNCTION nfmpi_bput_vars_int2(ncid, varid, start, count, stride, i2vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vars_int2

    INTEGER        FUNCTION nfmpi_bput_vars_int(ncid, varid, start, count, stride, ivals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vars_int

    INTEGER        FUNCTION nfmpi_bput_vars_real(ncid, varid, start, count, stride, rvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            REAL,                          INTENT(INOUT) :: rvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vars_real

    INTEGER        FUNCTION nfmpi_bput_vars_double(ncid, varid, start, count, stride, dvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vars_double

    INTEGER        FUNCTION nfmpi_bput_vars_int8(ncid, varid, start, count, stride, i8vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_vars_int8

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER        FUNCTION nfmpi_bput_varm(ncid, varid, start, count, stride, imap, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
!                           <type>,                        INTENT(IN)  :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_bput_varm

!   INTEGER        FUNCTION nfmpi_bget_varm(ncid, varid, start, count, stride, imap, buf, bufcount, datatype, req)
!                           INTEGER,                       INTENT(IN)  :: ncid
!                           INTEGER,                       INTENT(IN)  :: varid
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
!                           <type>,                        INTENT(OUT) :: buf(*)
!                           INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                           INTEGER,                       INTENT(IN)  :: datatype
!                           INTEGER,                       INTENT(OUT) :: req
!   END FUNCTION   nfmpi_bget_varm

    INTEGER        FUNCTION nfmpi_bput_varm_text(ncid, varid, start, count, stride, imap, text, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            CHARACTER(len=*),              INTENT(IN)  :: text
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_varm_text

    INTEGER        FUNCTION nfmpi_bput_varm_int1(ncid, varid, start, count, stride, imap, i1vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*1,                     INTENT(IN)  :: i1vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_varm_int1

    INTEGER        FUNCTION nfmpi_bput_varm_int2(ncid, varid, start, count, stride, imap, i2vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_varm_int2

    INTEGER        FUNCTION nfmpi_bput_varm_int (ncid, varid, start, count, stride, imap, ivals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER,                       INTENT(INOUT) :: ivals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_varm_int

    INTEGER        FUNCTION nfmpi_bput_varm_real(ncid, varid, start, count, stride, imap, rvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            REAL,                          INTENT(INOUT) :: rvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_varm_real

    INTEGER        FUNCTION nfmpi_bput_varm_double(ncid, varid, start, count, stride, imap, dvals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_varm_double

    INTEGER        FUNCTION nfmpi_bput_varm_int8(ncid, varid, start, count, stride, imap, i8vals, req)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: start(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: count(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: stride(*)
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: imap(*)
                            INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
                            INTEGER,                       INTENT(OUT) :: req
    END FUNCTION   nfmpi_bput_varm_int8

    INTEGER        FUNCTION nfmpi_buffer_attach(ncid, bufsize)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufsize
    END FUNCTION   nfmpi_buffer_attach

    INTEGER        FUNCTION nfmpi_buffer_detach(ncid)
                            INTEGER,                       INTENT(IN)  :: ncid
    END FUNCTION   nfmpi_buffer_detach

    INTEGER        FUNCTION nfmpi_inq_buffer_usage(ncid, usage)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(OUT) :: usage
    END FUNCTION   nfmpi_inq_buffer_usage

!
! End buffered put non-blocking subroutines:
!

    INTEGER        FUNCTION nfmpi_wait(ncid, count, req, status)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: count
                            INTEGER,                       INTENT(IN)  :: req(count)
                            INTEGER,                       INTENT(OUT) :: status(count)
    END FUNCTION   nfmpi_wait

    INTEGER        FUNCTION nfmpi_wait_all(ncid, count, req, status)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: count
                            INTEGER,                       INTENT(IN)  :: req(count)
                            INTEGER,                       INTENT(OUT) :: status(count)
    END FUNCTION   nfmpi_wait_all

    INTEGER        FUNCTION nfmpi_cancel(ncid, count, req, status)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: count
                            INTEGER,                       INTENT(IN)  :: req(count)
                            INTEGER,                       INTENT(OUT) :: status(count)
    END FUNCTION   nfmpi_cancel

    INTEGER        FUNCTION nfmpi_inq_put_size(ncid, size)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(OUT) :: size
    END FUNCTION   nfmpi_inq_put_size

    INTEGER        FUNCTION nfmpi_inq_get_size(ncid, size)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(OUT) :: size
    END FUNCTION   nfmpi_inq_get_size

    INTEGER        FUNCTION nfmpi_inq_header_size(ncid, size)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(OUT) :: size
    END FUNCTION   nfmpi_inq_header_size

    INTEGER        FUNCTION nfmpi_inq_header_extent(ncid, extent)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(OUT) :: extent
    END FUNCTION   nfmpi_inq_header_extent

    INTEGER        FUNCTION nfmpi_inq_varoffset(ncid, varid, offset)
                   use mpi, only: MPI_OFFSET_KIND
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(IN)  :: varid
                            INTEGER(KIND=MPI_OFFSET_KIND), INTENT(OUT) :: offset
    END FUNCTION   nfmpi_inq_varoffset

    INTEGER        FUNCTION nfmpi_inq_nreqs(ncid, nreqs)
                            INTEGER,                       INTENT(IN)  :: ncid
                            INTEGER,                       INTENT(OUT) :: nreqs
    END FUNCTION   nfmpi_inq_nreqs

!
! Begin of varn subroutines:
!

!
! flexible APIs, not ready yet for Fortran 77
!
!   INTEGER FUNCTION nfmpi_get_varn(ncid, varid, num, starts, counts, buf, bufcount, buftype)
!                    INTEGER,                       INTENT(IN)  :: ncid
!                    INTEGER,                       INTENT(IN)  :: varid
!                    INTEGER,                       INTENT(IN)  :: num
!                    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
!                    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
!                    <type>,                        INTENT(OUT) :: buf(*)
!                    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                    INTEGER,                       INTENT(IN)  :: buftype
!   END FUNCTION nfmpi_get_varn

!   INTEGER FUNCTION nfmpi_get_varn_all(ncid, varid, num, starts, counts, buf, bufcount, buftype)
!                    INTEGER,                       INTENT(IN)  :: ncid
!                    INTEGER,                       INTENT(IN)  :: varid
!                    INTEGER,                       INTENT(IN)  :: num
!                    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
!                    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
!                    <type>,                        INTENT(OUT) :: buf(*)
!                    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                    INTEGER,                       INTENT(IN)  :: buftype
!   END FUNCTION nfmpi_get_varn_all

!   INTEGER FUNCTION nfmpi_put_varn(ncid, varid, num, starts, counts, buf, bufcount, buftype)
!                    INTEGER,                       INTENT(IN)  :: ncid
!                    INTEGER,                       INTENT(IN)  :: varid
!                    INTEGER,                       INTENT(IN)  :: num
!                    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
!                    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
!                    <type>,                        INTENT(IN)  :: buf(*)
!                    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                    INTEGER,                       INTENT(IN)  :: buftype
!   END FUNCTION nfmpi_put_varn

!   INTEGER FUNCTION nfmpi_put_varn_all(ncid, varid, num, starts, counts, buf, bufcount, buftype)
!                    INTEGER,                       INTENT(IN)  :: ncid
!                    INTEGER,                       INTENT(IN)  :: varid
!                    INTEGER,                       INTENT(IN)  :: num
!                    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
!                    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
!                    <type>,                        INTENT(IN)  :: buf(*)
!                    INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: bufcount
!                    INTEGER,                       INTENT(IN)  :: buftype
!   END FUNCTION nfmpi_put_varn_all


    INTEGER FUNCTION nfmpi_get_varn_text(ncid, varid, num, starts, counts, text)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     CHARACTER(len=*),              INTENT(OUT) :: text
    END FUNCTION nfmpi_get_varn_text

    INTEGER FUNCTION nfmpi_get_varn_int1(ncid, varid, num, starts, counts, i1vals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER*1,                     INTENT(OUT) :: i1vals(*)
    END FUNCTION nfmpi_get_varn_int1

    INTEGER FUNCTION nfmpi_get_varn_int2(ncid, varid, num, starts, counts, i2vals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER*2,                     INTENT(OUT) :: i2vals(*)
    END FUNCTION nfmpi_get_varn_int2

    INTEGER FUNCTION nfmpi_get_varn_int(ncid, varid, num, starts, counts, ivals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER,                       INTENT(OUT) :: ivals(*)
    END FUNCTION nfmpi_get_varn_int

    INTEGER FUNCTION nfmpi_get_varn_real(ncid, varid, num, starts, counts, rvals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     REAL,                          INTENT(OUT) :: rvals(*)
    END FUNCTION nfmpi_get_varn_real

    INTEGER FUNCTION nfmpi_get_varn_double(ncid, varid, num, starts, counts, dvals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
    END FUNCTION nfmpi_get_varn_double

    INTEGER FUNCTION nfmpi_get_varn_int8(ncid, varid, num, starts, counts, i8vals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER*8,                     INTENT(OUT) :: i8vals(*)
    END FUNCTION nfmpi_get_varn_int8

    INTEGER FUNCTION nfmpi_get_varn_text_all(ncid, varid, num, starts, counts, text)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     CHARACTER(len=*),              INTENT(OUT) :: text
    END FUNCTION nfmpi_get_varn_text_all

    INTEGER FUNCTION nfmpi_get_varn_int1_all(ncid, varid, num, starts, counts, i1vals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER*1,                     INTENT(OUT) :: i1vals(*)
    END FUNCTION nfmpi_get_varn_int1_all

    INTEGER FUNCTION nfmpi_get_varn_int2_all(ncid, varid, num, starts, counts, i2vals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER*2,                     INTENT(OUT) :: i2vals(*)
    END FUNCTION nfmpi_get_varn_int2_all

    INTEGER FUNCTION nfmpi_get_varn_int_all(ncid, varid, num, starts, counts, ivals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER,                       INTENT(OUT) :: ivals(*)
    END FUNCTION nfmpi_get_varn_int_all

    INTEGER FUNCTION nfmpi_get_varn_real_all(ncid, varid, num, starts, counts, rvals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     REAL,                          INTENT(OUT) :: rvals(*)
    END FUNCTION nfmpi_get_varn_real_all

    INTEGER FUNCTION nfmpi_get_varn_double_all(ncid, varid, num, starts, counts, dvals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     DOUBLE PRECISION,              INTENT(OUT) :: dvals(*)
    END FUNCTION nfmpi_get_varn_double_all

    INTEGER FUNCTION nfmpi_get_varn_int8_all(ncid, varid, num, starts, counts, i8vals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER*8,                     INTENT(OUT) :: i8vals(*)
    END FUNCTION nfmpi_get_varn_int8_all


    INTEGER FUNCTION nfmpi_put_varn_text(ncid, varid, num, starts, counts, text)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     CHARACTER(len=*),              INTENT(IN)  :: text
    END FUNCTION nfmpi_put_varn_text

    INTEGER FUNCTION nfmpi_put_varn_int1(ncid, varid, num, starts, counts, i1vals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER*1,                     INTENT(IN)  :: i1vals(*)
    END FUNCTION nfmpi_put_varn_int1

    INTEGER FUNCTION nfmpi_put_varn_int2(ncid, varid, num, starts, counts, i2vals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
    END FUNCTION nfmpi_put_varn_int2

    INTEGER FUNCTION nfmpi_put_varn_int(ncid, varid, num, starts, counts, ivals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER,                       INTENT(INOUT) :: ivals(*)
    END FUNCTION nfmpi_put_varn_int

    INTEGER FUNCTION nfmpi_put_varn_real(ncid, varid, num, starts, counts, rvals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     REAL,                          INTENT(INOUT) :: rvals(*)
    END FUNCTION nfmpi_put_varn_real

    INTEGER FUNCTION nfmpi_put_varn_double(ncid, varid, num, starts, counts, dvals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
    END FUNCTION nfmpi_put_varn_double

    INTEGER FUNCTION nfmpi_put_varn_int8(ncid, varid, num, starts, counts, i8vals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
    END FUNCTION nfmpi_put_varn_int8

    INTEGER FUNCTION nfmpi_put_varn_text_all(ncid, varid, num, starts, counts, text)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     CHARACTER(len=*),              INTENT(IN)  :: text
    END FUNCTION nfmpi_put_varn_text_all

    INTEGER FUNCTION nfmpi_put_varn_int1_all(ncid, varid, num, starts, counts, i1vals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER*1,                     INTENT(IN)  :: i1vals(*)
    END FUNCTION nfmpi_put_varn_int1_all

    INTEGER FUNCTION nfmpi_put_varn_int2_all(ncid, varid, num, starts, counts, i2vals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER*2,                     INTENT(INOUT) :: i2vals(*)
    END FUNCTION nfmpi_put_varn_int2_all

    INTEGER FUNCTION nfmpi_put_varn_int_all(ncid, varid, num, starts, counts, ivals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER,                       INTENT(INOUT) :: ivals(*)
    END FUNCTION nfmpi_put_varn_int_all

    INTEGER FUNCTION nfmpi_put_varn_real_all(ncid, varid, num, starts, counts, rvals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     REAL,                          INTENT(INOUT) :: rvals(*)
    END FUNCTION nfmpi_put_varn_real_all

    INTEGER FUNCTION nfmpi_put_varn_double_all(ncid, varid, num, starts, counts, dvals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     DOUBLE PRECISION,              INTENT(INOUT) :: dvals(*)
    END FUNCTION nfmpi_put_varn_double_all

    INTEGER FUNCTION nfmpi_put_varn_int8_all(ncid, varid, num, starts, counts, i8vals)
            use mpi, only: MPI_OFFSET_KIND
                     INTEGER,                       INTENT(IN)  :: ncid
                     INTEGER,                       INTENT(IN)  :: varid
                     INTEGER,                       INTENT(IN)  :: num
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: starts(*)
                     INTEGER(KIND=MPI_OFFSET_KIND), INTENT(IN)  :: counts(*)
                     INTEGER*8,                     INTENT(INOUT) :: i8vals(*)
    END FUNCTION nfmpi_put_varn_int8_all

!
! End of varn subroutines:
!

END INTERFACE

! PnetCDF flexible APIs
      integer, external :: &
          nfmpi_put_var, &
          nfmpi_put_var1, &
          nfmpi_put_vara, &
          nfmpi_put_vars, &
          nfmpi_put_varm, &
          nfmpi_get_var, &
          nfmpi_get_var1, &
          nfmpi_get_vara, &
          nfmpi_get_vars, &
          nfmpi_get_varm

!         nfmpi_put_var_all
! collective put an entire variable does not make sense
!

      integer, external :: &
          nfmpi_put_vara_all, &
          nfmpi_put_vars_all, &
          nfmpi_put_varm_all, &
          nfmpi_get_var_all, &
          nfmpi_get_vara_all, &
          nfmpi_get_vars_all, &
          nfmpi_get_varm_all, &
          nfmpi_iput_var, &
          nfmpi_iput_var1, &
          nfmpi_iput_vara, &
          nfmpi_iput_vars, &
          nfmpi_iput_varm, &
          nfmpi_iget_var, &
          nfmpi_iget_var1, &
          nfmpi_iget_vara, &
          nfmpi_iget_vars, &
          nfmpi_iget_varm, &
          nfmpi_bput_var, &
          nfmpi_bput_var1, &
          nfmpi_bput_vara, &
          nfmpi_bput_vars, &
          nfmpi_bput_varm, &
          nfmpi_get_varn, &
          nfmpi_put_varn, &
          nfmpi_get_varn_all, &
          nfmpi_put_varn_all

