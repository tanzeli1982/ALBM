!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id: overloads.f90 1468 2013-10-26 16:53:18Z wkliao $
!
! This file is taken from netcdf_overloads.f90 with changes for PnetCDF use
!  
!

  ! Overloaded variable functions
  interface nf90mpi_def_var
    module procedure nf90mpi_def_var_Scalar, nf90mpi_def_var_oneDim, nf90mpi_def_var_ManyDims
  end interface ! nf90mpi_def_var
  
  ! Overloaded attribute functions
  interface nf90mpi_put_att
    module procedure nf90mpi_put_att_text,                                      &
                     nf90mpi_put_att_OneByteInt,     nf90mpi_put_att_TwoByteInt,   &
                     nf90mpi_put_att_FourByteInt,    nf90mpi_put_att_EightByteInt, &
                     nf90mpi_put_att_FourByteReal,   nf90mpi_put_att_EightByteReal
    module procedure nf90mpi_put_att_one_OneByteInt,   nf90mpi_put_att_one_TwoByteInt,   &
                     nf90mpi_put_att_one_FourByteInt,  nf90mpi_put_att_one_EightByteInt, &
                     nf90mpi_put_att_one_FourByteReal, nf90mpi_put_att_one_EightByteReal
  end interface !nf90mpi_put_att
  interface nf90mpi_get_att
    module procedure nf90mpi_get_att_text,                                      &
                     nf90mpi_get_att_OneByteInt,     nf90mpi_get_att_TwoByteInt,   &
                     nf90mpi_get_att_FourByteInt,    nf90mpi_get_att_EightByteInt, &
                     nf90mpi_get_att_FourByteReal,   nf90mpi_get_att_EightByteReal
    module procedure nf90mpi_get_att_one_OneByteInt,   nf90mpi_get_att_one_TwoByteInt,   &
                     nf90mpi_get_att_one_FourByteInt,  nf90mpi_get_att_one_EightByteInt, &
                     nf90mpi_get_att_one_FourByteReal, nf90mpi_get_att_one_EightByteReal
  end interface ! nf90mpi_get_att

  ! Overloaded variable functions
  interface nf90mpi_put_var
    module procedure nf90mpi_put_var_text,                                   &
                     nf90mpi_put_var_OneByteInt, nf90mpi_put_var_TwoByteInt,    &
                     nf90mpi_put_var_FourByteInt, nf90mpi_put_var_EightByteInt, &
                     nf90mpi_put_var_FourByteReal, nf90mpi_put_var_EightByteReal
    module procedure nf90mpi_put_var_1D_text,                                      &
                     nf90mpi_put_var_1D_OneByteInt, nf90mpi_put_var_1D_TwoByteInt,    &
                     nf90mpi_put_var_1D_FourByteInt, nf90mpi_put_var_1D_EightByteInt, &
                     nf90mpi_put_var_1D_FourByteReal, nf90mpi_put_var_1D_EightByteReal
    module procedure nf90mpi_put_var_2D_text,                                       &
                     nf90mpi_put_var_2D_OneByteInt, nf90mpi_put_var_2D_TwoByteInt,     &
                     nf90mpi_put_var_2D_FourByteInt, nf90mpi_put_var_2D_EightByteInt,  &
                     nf90mpi_put_var_2D_FourByteReal, nf90mpi_put_var_2D_EightByteReal
    module procedure nf90mpi_put_var_3D_text,                                       &
                     nf90mpi_put_var_3D_OneByteInt, nf90mpi_put_var_3D_TwoByteInt,     &
                     nf90mpi_put_var_3D_FourByteInt, nf90mpi_put_var_3D_EightByteInt,  &
                     nf90mpi_put_var_3D_FourByteReal, nf90mpi_put_var_3D_EightByteReal
    module procedure nf90mpi_put_var_4D_text,                                       &
                     nf90mpi_put_var_4D_OneByteInt, nf90mpi_put_var_4D_TwoByteInt,     &
                     nf90mpi_put_var_4D_FourByteInt, nf90mpi_put_var_4D_EightByteInt,  &
                     nf90mpi_put_var_4D_FourByteReal, nf90mpi_put_var_4D_EightByteReal
    module procedure nf90mpi_put_var_5D_text,                                       &
                     nf90mpi_put_var_5D_OneByteInt, nf90mpi_put_var_5D_TwoByteInt,     &
                     nf90mpi_put_var_5D_FourByteInt, nf90mpi_put_var_5D_EightByteInt,  &
                     nf90mpi_put_var_5D_FourByteReal, nf90mpi_put_var_5D_EightByteReal
    module procedure nf90mpi_put_var_6D_text,                                       &
                     nf90mpi_put_var_6D_OneByteInt, nf90mpi_put_var_6D_TwoByteInt,     &
                     nf90mpi_put_var_6D_FourByteInt, nf90mpi_put_var_6D_EightByteInt,  &
                     nf90mpi_put_var_6D_FourByteReal, nf90mpi_put_var_6D_EightByteReal
    module procedure nf90mpi_put_var_7D_text,                                       &
                     nf90mpi_put_var_7D_OneByteInt, nf90mpi_put_var_7D_TwoByteInt,     &
                     nf90mpi_put_var_7D_FourByteInt, nf90mpi_put_var_7D_EightByteInt,  &
                     nf90mpi_put_var_7D_FourByteReal, nf90mpi_put_var_7D_EightByteReal
  end interface ! nf90mpi_put_var

  interface nf90mpi_get_var
    module procedure nf90mpi_get_var_text,                                   &
                     nf90mpi_get_var_OneByteInt, nf90mpi_get_var_TwoByteInt,    &
                     nf90mpi_get_var_FourByteInt, nf90mpi_get_var_EightByteInt, &
                     nf90mpi_get_var_FourByteReal, nf90mpi_get_var_EightByteReal
    module procedure nf90mpi_get_var_1D_text,                                      &
                     nf90mpi_get_var_1D_OneByteInt, nf90mpi_get_var_1D_TwoByteInt,    &
                     nf90mpi_get_var_1D_FourByteInt, nf90mpi_get_var_1D_EightByteInt, &
                     nf90mpi_get_var_1D_FourByteReal, nf90mpi_get_var_1D_EightByteReal
    module procedure nf90mpi_get_var_2D_text,                                      &
                     nf90mpi_get_var_2D_OneByteInt, nf90mpi_get_var_2D_TwoByteInt,    &
                     nf90mpi_get_var_2D_FourByteInt, nf90mpi_get_var_2D_EightByteInt, &
                     nf90mpi_get_var_2D_FourByteReal, nf90mpi_get_var_2D_EightByteReal
    module procedure nf90mpi_get_var_3D_text,                                      &
                     nf90mpi_get_var_3D_OneByteInt, nf90mpi_get_var_3D_TwoByteInt,    &
                     nf90mpi_get_var_3D_FourByteInt, nf90mpi_get_var_3D_EightByteInt, &
                     nf90mpi_get_var_3D_FourByteReal, nf90mpi_get_var_3D_EightByteReal
    module procedure nf90mpi_get_var_4D_text,                                      &
                     nf90mpi_get_var_4D_OneByteInt, nf90mpi_get_var_4D_TwoByteInt,    &
                     nf90mpi_get_var_4D_FourByteInt, nf90mpi_get_var_4D_EightByteInt, &
                     nf90mpi_get_var_4D_FourByteReal, nf90mpi_get_var_4D_EightByteReal
    module procedure nf90mpi_get_var_5D_text,                                      &
                     nf90mpi_get_var_5D_OneByteInt, nf90mpi_get_var_5D_TwoByteInt,    &
                     nf90mpi_get_var_5D_FourByteInt, nf90mpi_get_var_5D_EightByteInt, &
                     nf90mpi_get_var_5D_FourByteReal, nf90mpi_get_var_5D_EightByteReal
    module procedure nf90mpi_get_var_6D_text,                                      &
                     nf90mpi_get_var_6D_OneByteInt, nf90mpi_get_var_6D_TwoByteInt,    &
                     nf90mpi_get_var_6D_FourByteInt, nf90mpi_get_var_6D_EightByteInt, &
                     nf90mpi_get_var_6D_FourByteReal, nf90mpi_get_var_6D_EightByteReal
    module procedure nf90mpi_get_var_7D_text,                                      &
                     nf90mpi_get_var_7D_OneByteInt, nf90mpi_get_var_7D_TwoByteInt,    &
                     nf90mpi_get_var_7D_FourByteInt, nf90mpi_get_var_7D_EightByteInt, &
                     nf90mpi_get_var_7D_FourByteReal, nf90mpi_get_var_7D_EightByteReal
  end interface ! nf90mpi_get_var

  ! Overloaded variable functions varn
  interface nf90mpi_put_varn
    module procedure nf90mpi_put_varn_text,                                   &
                     nf90mpi_put_varn_OneByteInt, nf90mpi_put_varn_TwoByteInt,    &
                     nf90mpi_put_varn_FourByteInt, nf90mpi_put_varn_EightByteInt, &
                     nf90mpi_put_varn_FourByteReal, nf90mpi_put_varn_EightByteReal
    module procedure nf90mpi_put_varn_1D_text,                                      &
                     nf90mpi_put_varn_1D_OneByteInt, nf90mpi_put_varn_1D_TwoByteInt,    &
                     nf90mpi_put_varn_1D_FourByteInt, nf90mpi_put_varn_1D_EightByteInt, &
                     nf90mpi_put_varn_1D_FourByteReal, nf90mpi_put_varn_1D_EightByteReal
    module procedure nf90mpi_put_varn_2D_text,                                       &
                     nf90mpi_put_varn_2D_OneByteInt, nf90mpi_put_varn_2D_TwoByteInt,     &
                     nf90mpi_put_varn_2D_FourByteInt, nf90mpi_put_varn_2D_EightByteInt,  &
                     nf90mpi_put_varn_2D_FourByteReal, nf90mpi_put_varn_2D_EightByteReal
    module procedure nf90mpi_put_varn_3D_text,                                       &
                     nf90mpi_put_varn_3D_OneByteInt, nf90mpi_put_varn_3D_TwoByteInt,     &
                     nf90mpi_put_varn_3D_FourByteInt, nf90mpi_put_varn_3D_EightByteInt,  &
                     nf90mpi_put_varn_3D_FourByteReal, nf90mpi_put_varn_3D_EightByteReal
    module procedure nf90mpi_put_varn_4D_text,                                       &
                     nf90mpi_put_varn_4D_OneByteInt, nf90mpi_put_varn_4D_TwoByteInt,     &
                     nf90mpi_put_varn_4D_FourByteInt, nf90mpi_put_varn_4D_EightByteInt,  &
                     nf90mpi_put_varn_4D_FourByteReal, nf90mpi_put_varn_4D_EightByteReal
    module procedure nf90mpi_put_varn_5D_text,                                       &
                     nf90mpi_put_varn_5D_OneByteInt, nf90mpi_put_varn_5D_TwoByteInt,     &
                     nf90mpi_put_varn_5D_FourByteInt, nf90mpi_put_varn_5D_EightByteInt,  &
                     nf90mpi_put_varn_5D_FourByteReal, nf90mpi_put_varn_5D_EightByteReal
    module procedure nf90mpi_put_varn_6D_text,                                       &
                     nf90mpi_put_varn_6D_OneByteInt, nf90mpi_put_varn_6D_TwoByteInt,     &
                     nf90mpi_put_varn_6D_FourByteInt, nf90mpi_put_varn_6D_EightByteInt,  &
                     nf90mpi_put_varn_6D_FourByteReal, nf90mpi_put_varn_6D_EightByteReal
    module procedure nf90mpi_put_varn_7D_text,                                       &
                     nf90mpi_put_varn_7D_OneByteInt, nf90mpi_put_varn_7D_TwoByteInt,     &
                     nf90mpi_put_varn_7D_FourByteInt, nf90mpi_put_varn_7D_EightByteInt,  &
                     nf90mpi_put_varn_7D_FourByteReal, nf90mpi_put_varn_7D_EightByteReal
  end interface ! nf90mpi_put_varn

  interface nf90mpi_get_varn
    module procedure nf90mpi_get_varn_text,                                   &
                     nf90mpi_get_varn_OneByteInt, nf90mpi_get_varn_TwoByteInt,    &
                     nf90mpi_get_varn_FourByteInt, nf90mpi_get_varn_EightByteInt, &
                     nf90mpi_get_varn_FourByteReal, nf90mpi_get_varn_EightByteReal
    module procedure nf90mpi_get_varn_1D_text,                                      &
                     nf90mpi_get_varn_1D_OneByteInt, nf90mpi_get_varn_1D_TwoByteInt,    &
                     nf90mpi_get_varn_1D_FourByteInt, nf90mpi_get_varn_1D_EightByteInt, &
                     nf90mpi_get_varn_1D_FourByteReal, nf90mpi_get_varn_1D_EightByteReal
    module procedure nf90mpi_get_varn_2D_text,                                      &
                     nf90mpi_get_varn_2D_OneByteInt, nf90mpi_get_varn_2D_TwoByteInt,    &
                     nf90mpi_get_varn_2D_FourByteInt, nf90mpi_get_varn_2D_EightByteInt, &
                     nf90mpi_get_varn_2D_FourByteReal, nf90mpi_get_varn_2D_EightByteReal
    module procedure nf90mpi_get_varn_3D_text,                                      &
                     nf90mpi_get_varn_3D_OneByteInt, nf90mpi_get_varn_3D_TwoByteInt,    &
                     nf90mpi_get_varn_3D_FourByteInt, nf90mpi_get_varn_3D_EightByteInt, &
                     nf90mpi_get_varn_3D_FourByteReal, nf90mpi_get_varn_3D_EightByteReal
    module procedure nf90mpi_get_varn_4D_text,                                      &
                     nf90mpi_get_varn_4D_OneByteInt, nf90mpi_get_varn_4D_TwoByteInt,    &
                     nf90mpi_get_varn_4D_FourByteInt, nf90mpi_get_varn_4D_EightByteInt, &
                     nf90mpi_get_varn_4D_FourByteReal, nf90mpi_get_varn_4D_EightByteReal
    module procedure nf90mpi_get_varn_5D_text,                                      &
                     nf90mpi_get_varn_5D_OneByteInt, nf90mpi_get_varn_5D_TwoByteInt,    &
                     nf90mpi_get_varn_5D_FourByteInt, nf90mpi_get_varn_5D_EightByteInt, &
                     nf90mpi_get_varn_5D_FourByteReal, nf90mpi_get_varn_5D_EightByteReal
    module procedure nf90mpi_get_varn_6D_text,                                      &
                     nf90mpi_get_varn_6D_OneByteInt, nf90mpi_get_varn_6D_TwoByteInt,    &
                     nf90mpi_get_varn_6D_FourByteInt, nf90mpi_get_varn_6D_EightByteInt, &
                     nf90mpi_get_varn_6D_FourByteReal, nf90mpi_get_varn_6D_EightByteReal
    module procedure nf90mpi_get_varn_7D_text,                                      &
                     nf90mpi_get_varn_7D_OneByteInt, nf90mpi_get_varn_7D_TwoByteInt,    &
                     nf90mpi_get_varn_7D_FourByteInt, nf90mpi_get_varn_7D_EightByteInt, &
                     nf90mpi_get_varn_7D_FourByteReal, nf90mpi_get_varn_7D_EightByteReal
  end interface ! nf90mpi_get_varn

  !
  ! Collective APIs
  !

  ! Overloaded variable functions
  interface nf90mpi_put_var_all
    module procedure nf90mpi_put_var_text_all
    module procedure nf90mpi_put_var_1D_text_all,                                      &
                     nf90mpi_put_var_1D_OneByteInt_all, nf90mpi_put_var_1D_TwoByteInt_all,    &
                     nf90mpi_put_var_1D_FourByteInt_all, nf90mpi_put_var_1D_EightByteInt_all, &
                     nf90mpi_put_var_1D_FourByteReal_all, nf90mpi_put_var_1D_EightByteReal_all
    module procedure nf90mpi_put_var_2D_text_all,                                       &
                     nf90mpi_put_var_2D_OneByteInt_all, nf90mpi_put_var_2D_TwoByteInt_all,     &
                     nf90mpi_put_var_2D_FourByteInt_all, nf90mpi_put_var_2D_EightByteInt_all,  &
                     nf90mpi_put_var_2D_FourByteReal_all, nf90mpi_put_var_2D_EightByteReal_all
    module procedure nf90mpi_put_var_3D_text_all,                                       &
                     nf90mpi_put_var_3D_OneByteInt_all, nf90mpi_put_var_3D_TwoByteInt_all,     &
                     nf90mpi_put_var_3D_FourByteInt_all, nf90mpi_put_var_3D_EightByteInt_all,  &
                     nf90mpi_put_var_3D_FourByteReal_all, nf90mpi_put_var_3D_EightByteReal_all
    module procedure nf90mpi_put_var_4D_text_all,                                       &
                     nf90mpi_put_var_4D_OneByteInt_all, nf90mpi_put_var_4D_TwoByteInt_all,     &
                     nf90mpi_put_var_4D_FourByteInt_all, nf90mpi_put_var_4D_EightByteInt_all,  &
                     nf90mpi_put_var_4D_FourByteReal_all, nf90mpi_put_var_4D_EightByteReal_all
    module procedure nf90mpi_put_var_5D_text_all,                                       &
                     nf90mpi_put_var_5D_OneByteInt_all, nf90mpi_put_var_5D_TwoByteInt_all,     &
                     nf90mpi_put_var_5D_FourByteInt_all, nf90mpi_put_var_5D_EightByteInt_all,  &
                     nf90mpi_put_var_5D_FourByteReal_all, nf90mpi_put_var_5D_EightByteReal_all
    module procedure nf90mpi_put_var_6D_text_all,                                       &
                     nf90mpi_put_var_6D_OneByteInt_all, nf90mpi_put_var_6D_TwoByteInt_all,     &
                     nf90mpi_put_var_6D_FourByteInt_all, nf90mpi_put_var_6D_EightByteInt_all,  &
                     nf90mpi_put_var_6D_FourByteReal_all, nf90mpi_put_var_6D_EightByteReal_all
    module procedure nf90mpi_put_var_7D_text_all,                                       &
                     nf90mpi_put_var_7D_OneByteInt_all, nf90mpi_put_var_7D_TwoByteInt_all,     &
                     nf90mpi_put_var_7D_FourByteInt_all, nf90mpi_put_var_7D_EightByteInt_all,  &
                     nf90mpi_put_var_7D_FourByteReal_all, nf90mpi_put_var_7D_EightByteReal_all
  end interface ! nf90mpi_put_var_all

  interface nf90mpi_get_var_all
    module procedure nf90mpi_get_var_text_all
    module procedure nf90mpi_get_var_1D_text_all,                                      &
                     nf90mpi_get_var_1D_OneByteInt_all, nf90mpi_get_var_1D_TwoByteInt_all,    &
                     nf90mpi_get_var_1D_FourByteInt_all, nf90mpi_get_var_1D_EightByteInt_all, &
                     nf90mpi_get_var_1D_FourByteReal_all, nf90mpi_get_var_1D_EightByteReal_all
    module procedure nf90mpi_get_var_2D_text_all,                                      &
                     nf90mpi_get_var_2D_OneByteInt_all, nf90mpi_get_var_2D_TwoByteInt_all,    &
                     nf90mpi_get_var_2D_FourByteInt_all, nf90mpi_get_var_2D_EightByteInt_all, &
                     nf90mpi_get_var_2D_FourByteReal_all, nf90mpi_get_var_2D_EightByteReal_all
    module procedure nf90mpi_get_var_3D_text_all,                                      &
                     nf90mpi_get_var_3D_OneByteInt_all, nf90mpi_get_var_3D_TwoByteInt_all,    &
                     nf90mpi_get_var_3D_FourByteInt_all, nf90mpi_get_var_3D_EightByteInt_all, &
                     nf90mpi_get_var_3D_FourByteReal_all, nf90mpi_get_var_3D_EightByteReal_all
    module procedure nf90mpi_get_var_4D_text_all,                                      &
                     nf90mpi_get_var_4D_OneByteInt_all, nf90mpi_get_var_4D_TwoByteInt_all,    &
                     nf90mpi_get_var_4D_FourByteInt_all, nf90mpi_get_var_4D_EightByteInt_all, &
                     nf90mpi_get_var_4D_FourByteReal_all, nf90mpi_get_var_4D_EightByteReal_all
    module procedure nf90mpi_get_var_5D_text_all,                                      &
                     nf90mpi_get_var_5D_OneByteInt_all, nf90mpi_get_var_5D_TwoByteInt_all,    &
                     nf90mpi_get_var_5D_FourByteInt_all, nf90mpi_get_var_5D_EightByteInt_all, &
                     nf90mpi_get_var_5D_FourByteReal_all, nf90mpi_get_var_5D_EightByteReal_all
    module procedure nf90mpi_get_var_6D_text_all,                                      &
                     nf90mpi_get_var_6D_OneByteInt_all, nf90mpi_get_var_6D_TwoByteInt_all,    &
                     nf90mpi_get_var_6D_FourByteInt_all, nf90mpi_get_var_6D_EightByteInt_all, &
                     nf90mpi_get_var_6D_FourByteReal_all, nf90mpi_get_var_6D_EightByteReal_all
    module procedure nf90mpi_get_var_7D_text_all,                                      &
                     nf90mpi_get_var_7D_OneByteInt_all, nf90mpi_get_var_7D_TwoByteInt_all,    &
                     nf90mpi_get_var_7D_FourByteInt_all, nf90mpi_get_var_7D_EightByteInt_all, &
                     nf90mpi_get_var_7D_FourByteReal_all, nf90mpi_get_var_7D_EightByteReal_all
  end interface ! nf90mpi_get_var_all

  ! Overloaded variable functions varn
  interface nf90mpi_put_varn_all
    module procedure nf90mpi_put_varn_text_all
    module procedure nf90mpi_put_varn_1D_text_all,                                      &
                     nf90mpi_put_varn_1D_OneByteInt_all, nf90mpi_put_varn_1D_TwoByteInt_all,    &
                     nf90mpi_put_varn_1D_FourByteInt_all, nf90mpi_put_varn_1D_EightByteInt_all, &
                     nf90mpi_put_varn_1D_FourByteReal_all, nf90mpi_put_varn_1D_EightByteReal_all
    module procedure nf90mpi_put_varn_2D_text_all,                                       &
                     nf90mpi_put_varn_2D_OneByteInt_all, nf90mpi_put_varn_2D_TwoByteInt_all,     &
                     nf90mpi_put_varn_2D_FourByteInt_all, nf90mpi_put_varn_2D_EightByteInt_all,  &
                     nf90mpi_put_varn_2D_FourByteReal_all, nf90mpi_put_varn_2D_EightByteReal_all
    module procedure nf90mpi_put_varn_3D_text_all,                                       &
                     nf90mpi_put_varn_3D_OneByteInt_all, nf90mpi_put_varn_3D_TwoByteInt_all,     &
                     nf90mpi_put_varn_3D_FourByteInt_all, nf90mpi_put_varn_3D_EightByteInt_all,  &
                     nf90mpi_put_varn_3D_FourByteReal_all, nf90mpi_put_varn_3D_EightByteReal_all
    module procedure nf90mpi_put_varn_4D_text_all,                                       &
                     nf90mpi_put_varn_4D_OneByteInt_all, nf90mpi_put_varn_4D_TwoByteInt_all,     &
                     nf90mpi_put_varn_4D_FourByteInt_all, nf90mpi_put_varn_4D_EightByteInt_all,  &
                     nf90mpi_put_varn_4D_FourByteReal_all, nf90mpi_put_varn_4D_EightByteReal_all
    module procedure nf90mpi_put_varn_5D_text_all,                                       &
                     nf90mpi_put_varn_5D_OneByteInt_all, nf90mpi_put_varn_5D_TwoByteInt_all,     &
                     nf90mpi_put_varn_5D_FourByteInt_all, nf90mpi_put_varn_5D_EightByteInt_all,  &
                     nf90mpi_put_varn_5D_FourByteReal_all, nf90mpi_put_varn_5D_EightByteReal_all
    module procedure nf90mpi_put_varn_6D_text_all,                                       &
                     nf90mpi_put_varn_6D_OneByteInt_all, nf90mpi_put_varn_6D_TwoByteInt_all,     &
                     nf90mpi_put_varn_6D_FourByteInt_all, nf90mpi_put_varn_6D_EightByteInt_all,  &
                     nf90mpi_put_varn_6D_FourByteReal_all, nf90mpi_put_varn_6D_EightByteReal_all
    module procedure nf90mpi_put_varn_7D_text_all,                                       &
                     nf90mpi_put_varn_7D_OneByteInt_all, nf90mpi_put_varn_7D_TwoByteInt_all,     &
                     nf90mpi_put_varn_7D_FourByteInt_all, nf90mpi_put_varn_7D_EightByteInt_all,  &
                     nf90mpi_put_varn_7D_FourByteReal_all, nf90mpi_put_varn_7D_EightByteReal_all
  end interface ! nf90mpi_put_varn_all

  interface nf90mpi_get_varn_all
    module procedure nf90mpi_get_varn_text_all
    module procedure nf90mpi_get_varn_1D_text_all,                                      &
                     nf90mpi_get_varn_1D_OneByteInt_all, nf90mpi_get_varn_1D_TwoByteInt_all,    &
                     nf90mpi_get_varn_1D_FourByteInt_all, nf90mpi_get_varn_1D_EightByteInt_all, &
                     nf90mpi_get_varn_1D_FourByteReal_all, nf90mpi_get_varn_1D_EightByteReal_all
    module procedure nf90mpi_get_varn_2D_text_all,                                      &
                     nf90mpi_get_varn_2D_OneByteInt_all, nf90mpi_get_varn_2D_TwoByteInt_all,    &
                     nf90mpi_get_varn_2D_FourByteInt_all, nf90mpi_get_varn_2D_EightByteInt_all, &
                     nf90mpi_get_varn_2D_FourByteReal_all, nf90mpi_get_varn_2D_EightByteReal_all
    module procedure nf90mpi_get_varn_3D_text_all,                                      &
                     nf90mpi_get_varn_3D_OneByteInt_all, nf90mpi_get_varn_3D_TwoByteInt_all,    &
                     nf90mpi_get_varn_3D_FourByteInt_all, nf90mpi_get_varn_3D_EightByteInt_all, &
                     nf90mpi_get_varn_3D_FourByteReal_all, nf90mpi_get_varn_3D_EightByteReal_all
    module procedure nf90mpi_get_varn_4D_text_all,                                      &
                     nf90mpi_get_varn_4D_OneByteInt_all, nf90mpi_get_varn_4D_TwoByteInt_all,    &
                     nf90mpi_get_varn_4D_FourByteInt_all, nf90mpi_get_varn_4D_EightByteInt_all, &
                     nf90mpi_get_varn_4D_FourByteReal_all, nf90mpi_get_varn_4D_EightByteReal_all
    module procedure nf90mpi_get_varn_5D_text_all,                                      &
                     nf90mpi_get_varn_5D_OneByteInt_all, nf90mpi_get_varn_5D_TwoByteInt_all,    &
                     nf90mpi_get_varn_5D_FourByteInt_all, nf90mpi_get_varn_5D_EightByteInt_all, &
                     nf90mpi_get_varn_5D_FourByteReal_all, nf90mpi_get_varn_5D_EightByteReal_all
    module procedure nf90mpi_get_varn_6D_text_all,                                      &
                     nf90mpi_get_varn_6D_OneByteInt_all, nf90mpi_get_varn_6D_TwoByteInt_all,    &
                     nf90mpi_get_varn_6D_FourByteInt_all, nf90mpi_get_varn_6D_EightByteInt_all, &
                     nf90mpi_get_varn_6D_FourByteReal_all, nf90mpi_get_varn_6D_EightByteReal_all
    module procedure nf90mpi_get_varn_7D_text_all,                                      &
                     nf90mpi_get_varn_7D_OneByteInt_all, nf90mpi_get_varn_7D_TwoByteInt_all,    &
                     nf90mpi_get_varn_7D_FourByteInt_all, nf90mpi_get_varn_7D_EightByteInt_all, &
                     nf90mpi_get_varn_7D_FourByteReal_all, nf90mpi_get_varn_7D_EightByteReal_all
  end interface ! nf90mpi_get_varn_all

  !
  ! Nonblocking APIs
  !
  interface nf90mpi_iput_var
    module procedure nf90mpi_iput_var_text,                                   &
                     nf90mpi_iput_var_OneByteInt, nf90mpi_iput_var_TwoByteInt,    &
                     nf90mpi_iput_var_FourByteInt, nf90mpi_iput_var_EightByteInt, &
                     nf90mpi_iput_var_FourByteReal, nf90mpi_iput_var_EightByteReal
    module procedure nf90mpi_iput_var_1D_text,                                      &
                     nf90mpi_iput_var_1D_OneByteInt, nf90mpi_iput_var_1D_TwoByteInt,    &
                     nf90mpi_iput_var_1D_FourByteInt, nf90mpi_iput_var_1D_EightByteInt, &
                     nf90mpi_iput_var_1D_FourByteReal, nf90mpi_iput_var_1D_EightByteReal
    module procedure nf90mpi_iput_var_2D_text,                                       &
                     nf90mpi_iput_var_2D_OneByteInt, nf90mpi_iput_var_2D_TwoByteInt,     &
                     nf90mpi_iput_var_2D_FourByteInt, nf90mpi_iput_var_2D_EightByteInt,  &
                     nf90mpi_iput_var_2D_FourByteReal, nf90mpi_iput_var_2D_EightByteReal
    module procedure nf90mpi_iput_var_3D_text,                                       &
                     nf90mpi_iput_var_3D_OneByteInt, nf90mpi_iput_var_3D_TwoByteInt,     &
                     nf90mpi_iput_var_3D_FourByteInt, nf90mpi_iput_var_3D_EightByteInt,  &
                     nf90mpi_iput_var_3D_FourByteReal, nf90mpi_iput_var_3D_EightByteReal
    module procedure nf90mpi_iput_var_4D_text,                                       &
                     nf90mpi_iput_var_4D_OneByteInt, nf90mpi_iput_var_4D_TwoByteInt,     &
                     nf90mpi_iput_var_4D_FourByteInt, nf90mpi_iput_var_4D_EightByteInt,  &
                     nf90mpi_iput_var_4D_FourByteReal, nf90mpi_iput_var_4D_EightByteReal
    module procedure nf90mpi_iput_var_5D_text,                                       &
                     nf90mpi_iput_var_5D_OneByteInt, nf90mpi_iput_var_5D_TwoByteInt,     &
                     nf90mpi_iput_var_5D_FourByteInt, nf90mpi_iput_var_5D_EightByteInt,  &
                     nf90mpi_iput_var_5D_FourByteReal, nf90mpi_iput_var_5D_EightByteReal
    module procedure nf90mpi_iput_var_6D_text,                                       &
                     nf90mpi_iput_var_6D_OneByteInt, nf90mpi_iput_var_6D_TwoByteInt,     &
                     nf90mpi_iput_var_6D_FourByteInt, nf90mpi_iput_var_6D_EightByteInt,  &
                     nf90mpi_iput_var_6D_FourByteReal, nf90mpi_iput_var_6D_EightByteReal
    module procedure nf90mpi_iput_var_7D_text,                                       &
                     nf90mpi_iput_var_7D_OneByteInt, nf90mpi_iput_var_7D_TwoByteInt,     &
                     nf90mpi_iput_var_7D_FourByteInt, nf90mpi_iput_var_7D_EightByteInt,  &
                     nf90mpi_iput_var_7D_FourByteReal, nf90mpi_iput_var_7D_EightByteReal
  end interface ! nf90mpi_iput_var

  interface nf90mpi_iget_var
    module procedure nf90mpi_iget_var_text,                                   &
                     nf90mpi_iget_var_OneByteInt, nf90mpi_iget_var_TwoByteInt,    &
                     nf90mpi_iget_var_FourByteInt, nf90mpi_iget_var_EightByteInt, &
                     nf90mpi_iget_var_FourByteReal, nf90mpi_iget_var_EightByteReal
    module procedure nf90mpi_iget_var_1D_text,                                      &
                     nf90mpi_iget_var_1D_OneByteInt, nf90mpi_iget_var_1D_TwoByteInt,    &
                     nf90mpi_iget_var_1D_FourByteInt, nf90mpi_iget_var_1D_EightByteInt, &
                     nf90mpi_iget_var_1D_FourByteReal, nf90mpi_iget_var_1D_EightByteReal
    module procedure nf90mpi_iget_var_2D_text,                                      &
                     nf90mpi_iget_var_2D_OneByteInt, nf90mpi_iget_var_2D_TwoByteInt,    &
                     nf90mpi_iget_var_2D_FourByteInt, nf90mpi_iget_var_2D_EightByteInt, &
                     nf90mpi_iget_var_2D_FourByteReal, nf90mpi_iget_var_2D_EightByteReal
    module procedure nf90mpi_iget_var_3D_text,                                      &
                     nf90mpi_iget_var_3D_OneByteInt, nf90mpi_iget_var_3D_TwoByteInt,    &
                     nf90mpi_iget_var_3D_FourByteInt, nf90mpi_iget_var_3D_EightByteInt, &
                     nf90mpi_iget_var_3D_FourByteReal, nf90mpi_iget_var_3D_EightByteReal
    module procedure nf90mpi_iget_var_4D_text,                                      &
                     nf90mpi_iget_var_4D_OneByteInt, nf90mpi_iget_var_4D_TwoByteInt,    &
                     nf90mpi_iget_var_4D_FourByteInt, nf90mpi_iget_var_4D_EightByteInt, &
                     nf90mpi_iget_var_4D_FourByteReal, nf90mpi_iget_var_4D_EightByteReal
    module procedure nf90mpi_iget_var_5D_text,                                      &
                     nf90mpi_iget_var_5D_OneByteInt, nf90mpi_iget_var_5D_TwoByteInt,    &
                     nf90mpi_iget_var_5D_FourByteInt, nf90mpi_iget_var_5D_EightByteInt, &
                     nf90mpi_iget_var_5D_FourByteReal, nf90mpi_iget_var_5D_EightByteReal
    module procedure nf90mpi_iget_var_6D_text,                                      &
                     nf90mpi_iget_var_6D_OneByteInt, nf90mpi_iget_var_6D_TwoByteInt,    &
                     nf90mpi_iget_var_6D_FourByteInt, nf90mpi_iget_var_6D_EightByteInt, &
                     nf90mpi_iget_var_6D_FourByteReal, nf90mpi_iget_var_6D_EightByteReal
    module procedure nf90mpi_iget_var_7D_text,                                      &
                     nf90mpi_iget_var_7D_OneByteInt, nf90mpi_iget_var_7D_TwoByteInt,    &
                     nf90mpi_iget_var_7D_FourByteInt, nf90mpi_iget_var_7D_EightByteInt, &
                     nf90mpi_iget_var_7D_FourByteReal, nf90mpi_iget_var_7D_EightByteReal
  end interface ! nf90mpi_iget_var

  interface nf90mpi_bput_var
    module procedure nf90mpi_bput_var_text,                                   &
                     nf90mpi_bput_var_OneByteInt, nf90mpi_bput_var_TwoByteInt,    &
                     nf90mpi_bput_var_FourByteInt, nf90mpi_bput_var_EightByteInt, &
                     nf90mpi_bput_var_FourByteReal, nf90mpi_bput_var_EightByteReal
    module procedure nf90mpi_bput_var_1D_text,                                      &
                     nf90mpi_bput_var_1D_OneByteInt, nf90mpi_bput_var_1D_TwoByteInt,    &
                     nf90mpi_bput_var_1D_FourByteInt, nf90mpi_bput_var_1D_EightByteInt, &
                     nf90mpi_bput_var_1D_FourByteReal, nf90mpi_bput_var_1D_EightByteReal
    module procedure nf90mpi_bput_var_2D_text,                                       &
                     nf90mpi_bput_var_2D_OneByteInt, nf90mpi_bput_var_2D_TwoByteInt,     &
                     nf90mpi_bput_var_2D_FourByteInt, nf90mpi_bput_var_2D_EightByteInt,  &
                     nf90mpi_bput_var_2D_FourByteReal, nf90mpi_bput_var_2D_EightByteReal
    module procedure nf90mpi_bput_var_3D_text,                                       &
                     nf90mpi_bput_var_3D_OneByteInt, nf90mpi_bput_var_3D_TwoByteInt,     &
                     nf90mpi_bput_var_3D_FourByteInt, nf90mpi_bput_var_3D_EightByteInt,  &
                     nf90mpi_bput_var_3D_FourByteReal, nf90mpi_bput_var_3D_EightByteReal
    module procedure nf90mpi_bput_var_4D_text,                                       &
                     nf90mpi_bput_var_4D_OneByteInt, nf90mpi_bput_var_4D_TwoByteInt,     &
                     nf90mpi_bput_var_4D_FourByteInt, nf90mpi_bput_var_4D_EightByteInt,  &
                     nf90mpi_bput_var_4D_FourByteReal, nf90mpi_bput_var_4D_EightByteReal
    module procedure nf90mpi_bput_var_5D_text,                                       &
                     nf90mpi_bput_var_5D_OneByteInt, nf90mpi_bput_var_5D_TwoByteInt,     &
                     nf90mpi_bput_var_5D_FourByteInt, nf90mpi_bput_var_5D_EightByteInt,  &
                     nf90mpi_bput_var_5D_FourByteReal, nf90mpi_bput_var_5D_EightByteReal
    module procedure nf90mpi_bput_var_6D_text,                                       &
                     nf90mpi_bput_var_6D_OneByteInt, nf90mpi_bput_var_6D_TwoByteInt,     &
                     nf90mpi_bput_var_6D_FourByteInt, nf90mpi_bput_var_6D_EightByteInt,  &
                     nf90mpi_bput_var_6D_FourByteReal, nf90mpi_bput_var_6D_EightByteReal
    module procedure nf90mpi_bput_var_7D_text,                                       &
                     nf90mpi_bput_var_7D_OneByteInt, nf90mpi_bput_var_7D_TwoByteInt,     &
                     nf90mpi_bput_var_7D_FourByteInt, nf90mpi_bput_var_7D_EightByteInt,  &
                     nf90mpi_bput_var_7D_FourByteReal, nf90mpi_bput_var_7D_EightByteReal
  end interface ! nf90mpi_bput_var

