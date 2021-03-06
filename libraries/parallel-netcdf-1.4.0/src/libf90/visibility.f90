!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id: visibility.f90 1468 2013-10-26 16:53:18Z wkliao $
!
! This file is taken from netcdf_visibility.f90 with changes for PnetCDF use
!
!

  ! Library version, error string
  public :: nf90mpi_inq_libvers, nf90mpi_strerror
  
  ! Control routines 
  public :: nf90mpi_create, nf90mpi_open, &
            nf90mpi_set_fill, nf90mpi_redef, nf90mpi_enddef,                     &
            nf90mpi_sync, nf90mpi_abort, nf90mpi_close, nf90mpi_delete
            
  ! File level inquiry
  public :: nf90mpi_inquire
  
  ! Dimension routines
  public :: nf90mpi_def_dim, nf90mpi_inq_dimid, nf90mpi_rename_dim, nf90mpi_inquire_dimension
  
  ! attribute routines
  public :: nf90mpi_copy_att, nf90mpi_rename_att, nf90mpi_del_att, nf90mpi_inq_attname, &
            nf90mpi_inquire_attribute
  ! overloaded functions
  public :: nf90mpi_put_att, nf90mpi_get_att
  
  ! Variable routines
  public :: nf90mpi_def_var, nf90mpi_inq_varid, nf90mpi_rename_var, nf90mpi_inquire_variable 

  ! overloaded functions
  ! independent APIs
  public :: nf90mpi_put_var, nf90mpi_get_var, &
            nf90mpi_put_varn, nf90mpi_get_varn

  ! collective APIs
  public :: nf90mpi_put_var_all, nf90mpi_get_var_all, &
            nf90mpi_put_varn_all, nf90mpi_get_varn_all

  ! nonblocking APIs
  public :: nf90mpi_iput_var, nf90mpi_iget_var, &
            nf90mpi_bput_var, nf90mpi_inq_nreqs, &
            nf90mpi_buffer_attach, nf90mpi_inq_buffer_usage, nf90mpi_buffer_detach, &
            nf90mpi_wait, nf90mpi_wait_all, nf90mpi_cancel

  public :: nf90mpi_begin_indep_data, nf90mpi_end_indep_data, &
            nf90mpi_inq_put_size, nf90mpi_inq_get_size, &
            nf90mpi_inq_header_size, nf90mpi_inq_header_extent, &
            nf90mpi_inq_varoffset, nf90mpi_get_file_info

!
! F77 APIs
!

    public :: &
        nfmpi_inq_libvers, &
        nfmpi_strerror, &
        nfmpi_issyserr
!
! control subroutines:
!
    public :: &
        nfmpi_create, &
        nfmpi_open, &
        nfmpi_inq_format, &
        nfmpi_inq_file_format, &
        nfmpi_get_file_info, &
        nfmpi_delete, &
        nfmpi_enddef, &
        nfmpi_redef, &
        nfmpi_set_default_format, &
        nfmpi_sync, &
        nfmpi_abort, &
        nfmpi_close, &
        nfmpi_set_fill
!
! general inquiry subroutines:
!
    public :: &
        nfmpi_inq, &
        nfmpi_inq_ndims, &
        nfmpi_inq_nvars, &
        nfmpi_inq_natts, &
        nfmpi_inq_unlimdim
!
! dimension subroutines:
!
    public :: &
        nfmpi_def_dim, &
        nfmpi_inq_dimid, &
        nfmpi_inq_dim, &
        nfmpi_inq_dimname, &
        nfmpi_inq_dimlen, &
        nfmpi_rename_dim
!
! general attribute subroutines:
!
    public :: &
        nfmpi_inq_att, &
        nfmpi_inq_attid, &
        nfmpi_inq_atttype, &
        nfmpi_inq_attlen, &
        nfmpi_inq_attname, &
        nfmpi_copy_att, &
        nfmpi_rename_att, &
        nfmpi_del_att
!
! attribute put/get subroutines:
!
    public :: &
        nfmpi_put_att_text, &
        nfmpi_get_att_text, &
        nfmpi_put_att_int1, &
        nfmpi_get_att_int1, &
        nfmpi_put_att_int2, &
        nfmpi_get_att_int2, &
        nfmpi_put_att_int, &
        nfmpi_get_att_int, &
        nfmpi_put_att_real, &
        nfmpi_get_att_real, &
        nfmpi_put_att_double, &
        nfmpi_get_att_double, &
        nfmpi_put_att_int8, &
        nfmpi_get_att_int8
!
! independent data mode subroutines:
!
    public :: &
        nfmpi_begin_indep_data, &
        nfmpi_end_indep_data
!
! general variable subroutines:
!
    public :: &
        nfmpi_def_var, &
        nfmpi_inq_var, &
        nfmpi_inq_varid, &
        nfmpi_inq_varname, &
        nfmpi_inq_vartype, &
        nfmpi_inq_varndims, &
        nfmpi_inq_vardimid, &
        nfmpi_inq_varnatts, &
        nfmpi_rename_var
!
! entire variable put/get subroutines:
!
    public :: &
        nfmpi_put_var_text, &
        nfmpi_get_var_text, &
        nfmpi_get_var_text_all, &
        nfmpi_put_var_int1, &
        nfmpi_get_var_int1, &
        nfmpi_get_var_int1_all, &
        nfmpi_put_var_int2, &
        nfmpi_get_var_int2, &
        nfmpi_get_var_int2_all, &
        nfmpi_put_var_int, &
        nfmpi_get_var_int, &
        nfmpi_get_var_int_all, &
        nfmpi_put_var_real, &
        nfmpi_get_var_real, &
        nfmpi_get_var_real_all, &
        nfmpi_put_var_double, &
        nfmpi_get_var_double, &
        nfmpi_get_var_double_all, &
        nfmpi_put_var_int8, &
        nfmpi_get_var_int8, &
        nfmpi_get_var_int8_all
!
! single variable put/get subroutines:
!
    public :: &
        nfmpi_put_var1_text, &
        nfmpi_get_var1_text, &
        nfmpi_put_var1_int1, &
        nfmpi_get_var1_int1, &
        nfmpi_put_var1_int2, &
        nfmpi_get_var1_int2, &
        nfmpi_put_var1_int, &
        nfmpi_get_var1_int, &
        nfmpi_put_var1_real, &
        nfmpi_get_var1_real, &
        nfmpi_put_var1_double, &
        nfmpi_get_var1_double, &
        nfmpi_put_var1_int8, &
        nfmpi_get_var1_int8

!
! variable array put/get subroutines:
!
    public :: &
        nfmpi_put_vara_text, &
        nfmpi_put_vara_text_all, &
        nfmpi_get_vara_text, &
        nfmpi_get_vara_text_all, &
        nfmpi_put_vara_int1, &
        nfmpi_put_vara_int1_all, &
        nfmpi_get_vara_int1, &
        nfmpi_get_vara_int1_all, &
        nfmpi_put_vara_int2, &
        nfmpi_put_vara_int2_all, &
        nfmpi_get_vara_int2, &
        nfmpi_get_vara_int2_all, &
        nfmpi_put_vara_int, &
        nfmpi_put_vara_int_all, &
        nfmpi_get_vara_int, &
        nfmpi_get_vara_int_all, &
        nfmpi_put_vara_real, &
        nfmpi_put_vara_real_all, &
        nfmpi_get_vara_real, &
        nfmpi_get_vara_real_all, &
        nfmpi_put_vara_double, &
        nfmpi_put_vara_double_all, &
        nfmpi_get_vara_double, &
        nfmpi_get_vara_double_all, &
        nfmpi_put_vara_int8, &
        nfmpi_put_vara_int8_all, &
        nfmpi_get_vara_int8, &
        nfmpi_get_vara_int8_all

!
! strided variable put/get subroutines:
!
    public :: &
        nfmpi_put_vars_text, &
        nfmpi_put_vars_text_all, &
        nfmpi_get_vars_text, &
        nfmpi_get_vars_text_all, &
        nfmpi_put_vars_int1, &
        nfmpi_put_vars_int1_all, &
        nfmpi_get_vars_int1, &
        nfmpi_get_vars_int1_all, &
        nfmpi_put_vars_int2, &
        nfmpi_put_vars_int2_all, &
        nfmpi_get_vars_int2, &
        nfmpi_get_vars_int2_all, &
        nfmpi_put_vars_int, &
        nfmpi_put_vars_int_all, &
        nfmpi_get_vars_int, &
        nfmpi_get_vars_int_all, &
        nfmpi_put_vars_real, &
        nfmpi_put_vars_real_all, &
        nfmpi_get_vars_real, &
        nfmpi_get_vars_real_all, &
        nfmpi_put_vars_double, &
        nfmpi_put_vars_double_all, &
        nfmpi_get_vars_double, &
        nfmpi_get_vars_double_all, &
        nfmpi_put_vars_int8, &
        nfmpi_put_vars_int8_all, &
        nfmpi_get_vars_int8, &
        nfmpi_get_vars_int8_all

!
! mapped variable put/get subroutines:
!
    public :: &
        nfmpi_put_varm_text, &
        nfmpi_put_varm_text_all, &
        nfmpi_get_varm_text, &
        nfmpi_get_varm_text_all, &
        nfmpi_put_varm_int1, &
        nfmpi_put_varm_int1_all, &
        nfmpi_get_varm_int1, &
        nfmpi_get_varm_int1_all, &
        nfmpi_put_varm_int2, &
        nfmpi_put_varm_int2_all, &
        nfmpi_get_varm_int2, &
        nfmpi_get_varm_int2_all, &
        nfmpi_put_varm_int , &
        nfmpi_put_varm_int_all, &
        nfmpi_get_varm_int , &
        nfmpi_get_varm_int_all, &
        nfmpi_put_varm_real, &
        nfmpi_put_varm_real_all, &
        nfmpi_get_varm_real, &
        nfmpi_get_varm_real_all, &
        nfmpi_put_varm_double, &
        nfmpi_put_varm_double_all, &
        nfmpi_get_varm_double, &
        nfmpi_get_varm_double_all, &
        nfmpi_put_varm_int8, &
        nfmpi_put_varm_int8_all, &
        nfmpi_get_varm_int8, &
        nfmpi_get_varm_int8_all

!
! non-blocking variable array iput/iget subroutines:
!

!
! entire variable iput/iget subroutines:
!
    public :: &
        nfmpi_iput_var_text, &
        nfmpi_iget_var_text, &
        nfmpi_iput_var_int1, &
        nfmpi_iget_var_int1, &
        nfmpi_iput_var_int2, &
        nfmpi_iget_var_int2, &
        nfmpi_iput_var_int, &
        nfmpi_iget_var_int, &
        nfmpi_iput_var_real, &
        nfmpi_iget_var_real, &
        nfmpi_iput_var_double, &
        nfmpi_iget_var_double, &
        nfmpi_iput_var_int8, &
        nfmpi_iget_var_int8

!
! single variable iput/iget subroutines:
!
    public :: &
        nfmpi_iput_var1_text, &
        nfmpi_iget_var1_text, &
        nfmpi_iput_var1_int1, &
        nfmpi_iget_var1_int1, &
        nfmpi_iput_var1_int2, &
        nfmpi_iget_var1_int2, &
        nfmpi_iput_var1_int, &
        nfmpi_iget_var1_int, &
        nfmpi_iput_var1_real, &
        nfmpi_iget_var1_real, &
        nfmpi_iput_var1_double, &
        nfmpi_iget_var1_double, &
        nfmpi_iput_var1_int8, &
        nfmpi_iget_var1_int8

!
! variable array iput/iget subroutines:
!
    public :: &
        nfmpi_iput_vara_text, &
        nfmpi_iget_vara_text, &
        nfmpi_iput_vara_int1, &
        nfmpi_iget_vara_int1, &
        nfmpi_iput_vara_int2, &
        nfmpi_iget_vara_int2, &
        nfmpi_iput_vara_int, &
        nfmpi_iget_vara_int, &
        nfmpi_iput_vara_real, &
        nfmpi_iget_vara_real, &
        nfmpi_iput_vara_double, &
        nfmpi_iget_vara_double, &
        nfmpi_iput_vara_int8, &
        nfmpi_iget_vara_int8

!
! strided variable iput/iget subroutines:
!
    public :: &
        nfmpi_iput_vars_text, &
        nfmpi_iget_vars_text, &
        nfmpi_iput_vars_int1, &
        nfmpi_iget_vars_int1, &
        nfmpi_iput_vars_int2, &
        nfmpi_iget_vars_int2, &
        nfmpi_iput_vars_int, &
        nfmpi_iget_vars_int, &
        nfmpi_iput_vars_real, &
        nfmpi_iget_vars_real, &
        nfmpi_iput_vars_double, &
        nfmpi_iget_vars_double, &
        nfmpi_iput_vars_int8, &
        nfmpi_iget_vars_int8

!
! mapped variable iput/iget subroutines:
!
    public :: &
        nfmpi_iput_varm_text, &
        nfmpi_iget_varm_text, &
        nfmpi_iput_varm_int1, &
        nfmpi_iget_varm_int1, &
        nfmpi_iput_varm_int2, &
        nfmpi_iget_varm_int2, &
        nfmpi_iput_varm_int , &
        nfmpi_iget_varm_int , &
        nfmpi_iput_varm_real, &
        nfmpi_iget_varm_real, &
        nfmpi_iput_varm_double, &
        nfmpi_iget_varm_double, &
        nfmpi_iput_varm_int8, &
        nfmpi_iget_varm_int8

!
! Begin buffered put non-blocking subroutines:
!
    public :: &
        nfmpi_bput_var_text, &
        nfmpi_bput_var_int1, &
        nfmpi_bput_var_int2, &
        nfmpi_bput_var_int, &
        nfmpi_bput_var_real, &
        nfmpi_bput_var_double, &
        nfmpi_bput_var_int8, &
        nfmpi_bput_var1_text, &
        nfmpi_bput_var1_int1, &
        nfmpi_bput_var1_int2, &
        nfmpi_bput_var1_int, &
        nfmpi_bput_var1_real, &
        nfmpi_bput_var1_double, &
        nfmpi_bput_var1_int8, &
        nfmpi_bput_vara_text, &
        nfmpi_bput_vara_int1, &
        nfmpi_bput_vara_int2, &
        nfmpi_bput_vara_int, &
        nfmpi_bput_vara_real, &
        nfmpi_bput_vara_double, &
        nfmpi_bput_vara_int8, &
        nfmpi_bput_vars_text, &
        nfmpi_bput_vars_int1, &
        nfmpi_bput_vars_int2, &
        nfmpi_bput_vars_int, &
        nfmpi_bput_vars_real, &
        nfmpi_bput_vars_double, &
        nfmpi_bput_vars_int8, &
        nfmpi_bput_varm_text, &
        nfmpi_bput_varm_int1, &
        nfmpi_bput_varm_int2, &
        nfmpi_bput_varm_int , &
        nfmpi_bput_varm_real, &
        nfmpi_bput_varm_double, &
        nfmpi_bput_varm_int8, &
        nfmpi_buffer_attach, &
        nfmpi_buffer_detach, &
        nfmpi_inq_buffer_usage

!
! End buffered put non-blocking subroutines:
!
    public :: &
        nfmpi_wait, &
        nfmpi_wait_all, &
        nfmpi_cancel, &
        nfmpi_inq_put_size, &
        nfmpi_inq_get_size, &
        nfmpi_inq_header_size, &
        nfmpi_inq_header_extent, &
        nfmpi_inq_varoffset, &
        nfmpi_inq_nreqs

!
! Begin of varn subroutines:
!
    public :: &
        nfmpi_get_varn_text, &
        nfmpi_get_varn_int1, &
        nfmpi_get_varn_int2, &
        nfmpi_get_varn_int, &
        nfmpi_get_varn_real, &
        nfmpi_get_varn_double, &
        nfmpi_get_varn_int8, &
        nfmpi_get_varn_text_all, &
        nfmpi_get_varn_int1_all, &
        nfmpi_get_varn_int2_all, &
        nfmpi_get_varn_int_all, &
        nfmpi_get_varn_real_all, &
        nfmpi_get_varn_double_all, &
        nfmpi_get_varn_int8_all, &
        nfmpi_put_varn_text, &
        nfmpi_put_varn_int1, &
        nfmpi_put_varn_int2, &
        nfmpi_put_varn_int, &
        nfmpi_put_varn_real, &
        nfmpi_put_varn_double, &
        nfmpi_put_varn_int8, &
        nfmpi_put_varn_text_all, &
        nfmpi_put_varn_int1_all, &
        nfmpi_put_varn_int2_all, &
        nfmpi_put_varn_int_all, &
        nfmpi_put_varn_real_all, &
        nfmpi_put_varn_double_all, &
        nfmpi_put_varn_int8_all

!
! End of varn subroutines:
!

! PnetCDF flexible APIs
      public :: &
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

      public :: &
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
