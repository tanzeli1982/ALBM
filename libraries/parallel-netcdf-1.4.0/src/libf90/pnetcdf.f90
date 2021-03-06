!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id: pnetcdf.f90 1468 2013-10-26 16:53:18Z wkliao $
!
! This file is taken from netcdf.f90 with changes for PnetCDF use
!
!
 module pnetcdf
  use mpi, only: MPI_OFFSET_KIND

  implicit none
  private

  integer, parameter ::   OneByteInt = selected_int_kind(2), &
                          TwoByteInt = selected_int_kind(4), &
                         FourByteInt = selected_int_kind(9), &
                        EightByteInt = selected_int_kind(18)

  integer, parameter ::                                          &
                        FourByteReal = selected_real_kind(P =  6, R =  37), &
                       EightByteReal = selected_real_kind(P = 13, R = 307)

  include "nfmpi_constants.f90"
  include "nf90_constants.f90"
  include "api.f90"
  include "overloads.f90"
  include "visibility.f90"

contains
  include "file.f90"
  include "dims.f90"
  include "attributes.f90"
  include "variables.f90"
  include "getput_text.f90"
  include "getput_var.f90"
  include "getput_varn.f90"
end module pnetcdf
