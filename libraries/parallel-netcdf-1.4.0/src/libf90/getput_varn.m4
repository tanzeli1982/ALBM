dnl Process this m4 file to produce 'C' language file.
dnl
dnl If you see this line, you can ignore the next one.
! Do not edit this file. It is produced from the corresponding .m4 source
dnl
!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
! $Id: getput_varn.m4 1468 2013-10-26 16:53:18Z wkliao $
!

dnl
dnl PUTVARN1(ncid, varid, values, num, start, count)
dnl
define(`PUTVARN1',dnl
`dnl
   function nf90mpi_put_varn_$1$4(ncid, varid, values, num, start, count)
     integer,                                        intent( in) :: ncid, varid, num
     $2 (kind=$1),                                   intent( in) :: values
     integer (kind=MPI_OFFSET_KIND), dimension(:,:), intent( in) :: start, count
     integer                                                     :: nf90mpi_put_varn_$1$4
     $2 (kind=$1),                   dimension(1)                :: tempValue
     tempValue(1) = values
     nf90mpi_put_varn_$1$4 = nfmpi_put_varn_$3$4(ncid, varid, num, start, count, tempValue)
   end function nf90mpi_put_varn_$1$4
')dnl

dnl
dnl GETVARN1(ncid, varid, values, num, start, count)
dnl
define(`GETVARN1',dnl
`dnl
   function nf90mpi_get_varn_$1$4(ncid, varid, values, num, start, count)
     integer,                                        intent( in) :: ncid, varid, num
     $2 (kind=$1),                                   intent(out) :: values
     integer (kind=MPI_OFFSET_KIND), dimension(:,:), intent( in) :: start, count
     integer                                                     :: nf90mpi_get_varn_$1$4
     $2 (kind=$1),                   dimension(1)                :: tempValue
     nf90mpi_get_varn_$1$4 = nfmpi_get_varn_$3$4(ncid, varid, num, start, count, tempValue)
     values = tempValue(1)
   end function nf90mpi_get_varn_$1$4
')dnl

!
! Independent put APIs
!

PUTVARN1(OneByteInt,    integer, int1)
PUTVARN1(TwoByteInt,    integer, int2)
PUTVARN1(FourByteInt,   integer, int)
PUTVARN1(FourByteReal,  real,    real)
PUTVARN1(EightByteReal, real,    double)
PUTVARN1(EightByteInt,  integer, int8)

!
! Independent get APIs
!

GETVARN1(OneByteInt,    integer, int1)
GETVARN1(TwoByteInt,    integer, int2)
GETVARN1(FourByteInt,   integer, int)
GETVARN1(FourByteReal,  real,    real)
GETVARN1(EightByteReal, real,    double)
GETVARN1(EightByteInt,  integer, int8)

!
! Collective put APIs
!

PUTVARN1(OneByteInt,    integer, int1,   _all)
PUTVARN1(TwoByteInt,    integer, int2,   _all)
PUTVARN1(FourByteInt,   integer, int,    _all)
PUTVARN1(FourByteReal,  real,    real,   _all)
PUTVARN1(EightByteReal, real,    double, _all)
PUTVARN1(EightByteInt,  integer, int8,   _all)

!
! Collective get APIs
!

GETVARN1(OneByteInt,    integer, int1,   _all)
GETVARN1(TwoByteInt,    integer, int2,   _all)
GETVARN1(FourByteInt,   integer, int,    _all)
GETVARN1(FourByteReal,  real,    real,   _all)
GETVARN1(EightByteReal, real,    double, _all)
GETVARN1(EightByteInt,  integer, int8,   _all)


dnl
dnl VARN(ncid, varid, values, num, start, count)
dnl
define(`VARN',dnl
`dnl
   function nf90mpi_$1_varn_$2_$3$8(ncid, varid, values, num, start, count)
     integer,                                        intent(in) :: ncid, varid, num
     $4 (kind=$3),   dimension($6),                  intent($7) :: values
     integer (kind=MPI_OFFSET_KIND), dimension(:,:), intent(in) :: start, count
     integer                                                    :: nf90mpi_$1_varn_$2_$3$8
     nf90mpi_$1_varn_$2_$3$8 = nfmpi_$1_varn_$5$8(ncid, varid, num, start, count, values)
   end function nf90mpi_$1_varn_$2_$3$8
')dnl

!
! put APIs
!

VARN(put, 1D, OneByteInt, integer, int1,  :,              in)
VARN(put, 2D, OneByteInt, integer, int1, `:,:',           in)
VARN(put, 3D, OneByteInt, integer, int1, `:,:,:',         in)
VARN(put, 4D, OneByteInt, integer, int1, `:,:,:,:',       in)
VARN(put, 5D, OneByteInt, integer, int1, `:,:,:,:,:',     in)
VARN(put, 6D, OneByteInt, integer, int1, `:,:,:,:,:,:',   in)
VARN(put, 7D, OneByteInt, integer, int1, `:,:,:,:,:,:,:', in)

VARN(put, 1D, TwoByteInt, integer, int2,  :,              inout)
VARN(put, 2D, TwoByteInt, integer, int2, `:,:',           inout)
VARN(put, 3D, TwoByteInt, integer, int2, `:,:,:',         inout)
VARN(put, 4D, TwoByteInt, integer, int2, `:,:,:,:',       inout)
VARN(put, 5D, TwoByteInt, integer, int2, `:,:,:,:,:',     inout)
VARN(put, 6D, TwoByteInt, integer, int2, `:,:,:,:,:,:',   inout)
VARN(put, 7D, TwoByteInt, integer, int2, `:,:,:,:,:,:,:', inout)

VARN(put, 1D, FourByteInt, integer, int,  :,              inout)
VARN(put, 2D, FourByteInt, integer, int, `:,:',           inout)
VARN(put, 3D, FourByteInt, integer, int, `:,:,:',         inout)
VARN(put, 4D, FourByteInt, integer, int, `:,:,:,:',       inout)
VARN(put, 5D, FourByteInt, integer, int, `:,:,:,:,:',     inout)
VARN(put, 6D, FourByteInt, integer, int, `:,:,:,:,:,:',   inout)
VARN(put, 7D, FourByteInt, integer, int, `:,:,:,:,:,:,:', inout)

VARN(put, 1D, FourByteReal, real,   real,  :,              inout)
VARN(put, 2D, FourByteReal, real,   real, `:,:',           inout)
VARN(put, 3D, FourByteReal, real,   real, `:,:,:',         inout)
VARN(put, 4D, FourByteReal, real,   real, `:,:,:,:',       inout)
VARN(put, 5D, FourByteReal, real,   real, `:,:,:,:,:',     inout)
VARN(put, 6D, FourByteReal, real,   real, `:,:,:,:,:,:',   inout)
VARN(put, 7D, FourByteReal, real,   real, `:,:,:,:,:,:,:', inout)

VARN(put, 1D, EightByteReal, real, double,  :,              inout)
VARN(put, 2D, EightByteReal, real, double, `:,:',           inout)
VARN(put, 3D, EightByteReal, real, double, `:,:,:',         inout)
VARN(put, 4D, EightByteReal, real, double, `:,:,:,:',       inout)
VARN(put, 5D, EightByteReal, real, double, `:,:,:,:,:',     inout)
VARN(put, 6D, EightByteReal, real, double, `:,:,:,:,:,:',   inout)
VARN(put, 7D, EightByteReal, real, double, `:,:,:,:,:,:,:', inout)

VARN(put, 1D, EightByteInt, integer, int8,  :,              inout)
VARN(put, 2D, EightByteInt, integer, int8, `:,:',           inout)
VARN(put, 3D, EightByteInt, integer, int8, `:,:,:',         inout)
VARN(put, 4D, EightByteInt, integer, int8, `:,:,:,:',       inout)
VARN(put, 5D, EightByteInt, integer, int8, `:,:,:,:,:',     inout)
VARN(put, 6D, EightByteInt, integer, int8, `:,:,:,:,:,:',   inout)
VARN(put, 7D, EightByteInt, integer, int8, `:,:,:,:,:,:,:', inout)

!
! get APIs
!

VARN(get, 1D, OneByteInt, integer, int1,  :,              out)
VARN(get, 2D, OneByteInt, integer, int1, `:,:',           out)
VARN(get, 3D, OneByteInt, integer, int1, `:,:,:',         out)
VARN(get, 4D, OneByteInt, integer, int1, `:,:,:,:',       out)
VARN(get, 5D, OneByteInt, integer, int1, `:,:,:,:,:',     out)
VARN(get, 6D, OneByteInt, integer, int1, `:,:,:,:,:,:',   out)
VARN(get, 7D, OneByteInt, integer, int1, `:,:,:,:,:,:,:', out)

VARN(get, 1D, TwoByteInt, integer, int2,  :,              out)
VARN(get, 2D, TwoByteInt, integer, int2, `:,:',           out)
VARN(get, 3D, TwoByteInt, integer, int2, `:,:,:',         out)
VARN(get, 4D, TwoByteInt, integer, int2, `:,:,:,:',       out)
VARN(get, 5D, TwoByteInt, integer, int2, `:,:,:,:,:',     out)
VARN(get, 6D, TwoByteInt, integer, int2, `:,:,:,:,:,:',   out)
VARN(get, 7D, TwoByteInt, integer, int2, `:,:,:,:,:,:,:', out)

VARN(get, 1D, FourByteInt, integer, int,  :,              out)
VARN(get, 2D, FourByteInt, integer, int, `:,:',           out)
VARN(get, 3D, FourByteInt, integer, int, `:,:,:',         out)
VARN(get, 4D, FourByteInt, integer, int, `:,:,:,:',       out)
VARN(get, 5D, FourByteInt, integer, int, `:,:,:,:,:',     out)
VARN(get, 6D, FourByteInt, integer, int, `:,:,:,:,:,:',   out)
VARN(get, 7D, FourByteInt, integer, int, `:,:,:,:,:,:,:', out)

VARN(get, 1D, FourByteReal, real,   real,  :,              out)
VARN(get, 2D, FourByteReal, real,   real, `:,:',           out)
VARN(get, 3D, FourByteReal, real,   real, `:,:,:',         out)
VARN(get, 4D, FourByteReal, real,   real, `:,:,:,:',       out)
VARN(get, 5D, FourByteReal, real,   real, `:,:,:,:,:',     out)
VARN(get, 6D, FourByteReal, real,   real, `:,:,:,:,:,:',   out)
VARN(get, 7D, FourByteReal, real,   real, `:,:,:,:,:,:,:', out)

VARN(get, 1D, EightByteReal, real, double,  :,              out)
VARN(get, 2D, EightByteReal, real, double, `:,:',           out)
VARN(get, 3D, EightByteReal, real, double, `:,:,:',         out)
VARN(get, 4D, EightByteReal, real, double, `:,:,:,:',       out)
VARN(get, 5D, EightByteReal, real, double, `:,:,:,:,:',     out)
VARN(get, 6D, EightByteReal, real, double, `:,:,:,:,:,:',   out)
VARN(get, 7D, EightByteReal, real, double, `:,:,:,:,:,:,:', out)

VARN(get, 1D, EightByteInt, integer, int8,  :,              out)
VARN(get, 2D, EightByteInt, integer, int8, `:,:',           out)
VARN(get, 3D, EightByteInt, integer, int8, `:,:,:',         out)
VARN(get, 4D, EightByteInt, integer, int8, `:,:,:,:',       out)
VARN(get, 5D, EightByteInt, integer, int8, `:,:,:,:,:',     out)
VARN(get, 6D, EightByteInt, integer, int8, `:,:,:,:,:,:',   out)
VARN(get, 7D, EightByteInt, integer, int8, `:,:,:,:,:,:,:', out)

!
! collective put APIs
!

VARN(put, 1D, OneByteInt, integer, int1,  :,              in, _all)
VARN(put, 2D, OneByteInt, integer, int1, `:,:',           in, _all)
VARN(put, 3D, OneByteInt, integer, int1, `:,:,:',         in, _all)
VARN(put, 4D, OneByteInt, integer, int1, `:,:,:,:',       in, _all)
VARN(put, 5D, OneByteInt, integer, int1, `:,:,:,:,:',     in, _all)
VARN(put, 6D, OneByteInt, integer, int1, `:,:,:,:,:,:',   in, _all)
VARN(put, 7D, OneByteInt, integer, int1, `:,:,:,:,:,:,:', in, _all)

VARN(put, 1D, TwoByteInt, integer, int2,  :,              inout, _all)
VARN(put, 2D, TwoByteInt, integer, int2, `:,:',           inout, _all)
VARN(put, 3D, TwoByteInt, integer, int2, `:,:,:',         inout, _all)
VARN(put, 4D, TwoByteInt, integer, int2, `:,:,:,:',       inout, _all)
VARN(put, 5D, TwoByteInt, integer, int2, `:,:,:,:,:',     inout, _all)
VARN(put, 6D, TwoByteInt, integer, int2, `:,:,:,:,:,:',   inout, _all)
VARN(put, 7D, TwoByteInt, integer, int2, `:,:,:,:,:,:,:', inout, _all)

VARN(put, 1D, FourByteInt, integer, int,  :,              inout, _all)
VARN(put, 2D, FourByteInt, integer, int, `:,:',           inout, _all)
VARN(put, 3D, FourByteInt, integer, int, `:,:,:',         inout, _all)
VARN(put, 4D, FourByteInt, integer, int, `:,:,:,:',       inout, _all)
VARN(put, 5D, FourByteInt, integer, int, `:,:,:,:,:',     inout, _all)
VARN(put, 6D, FourByteInt, integer, int, `:,:,:,:,:,:',   inout, _all)
VARN(put, 7D, FourByteInt, integer, int, `:,:,:,:,:,:,:', inout, _all)

VARN(put, 1D, FourByteReal, real,   real,  :,              inout, _all)
VARN(put, 2D, FourByteReal, real,   real, `:,:',           inout, _all)
VARN(put, 3D, FourByteReal, real,   real, `:,:,:',         inout, _all)
VARN(put, 4D, FourByteReal, real,   real, `:,:,:,:',       inout, _all)
VARN(put, 5D, FourByteReal, real,   real, `:,:,:,:,:',     inout, _all)
VARN(put, 6D, FourByteReal, real,   real, `:,:,:,:,:,:',   inout, _all)
VARN(put, 7D, FourByteReal, real,   real, `:,:,:,:,:,:,:', inout, _all)

VARN(put, 1D, EightByteReal, real, double,  :,              inout, _all)
VARN(put, 2D, EightByteReal, real, double, `:,:',           inout, _all)
VARN(put, 3D, EightByteReal, real, double, `:,:,:',         inout, _all)
VARN(put, 4D, EightByteReal, real, double, `:,:,:,:',       inout, _all)
VARN(put, 5D, EightByteReal, real, double, `:,:,:,:,:',     inout, _all)
VARN(put, 6D, EightByteReal, real, double, `:,:,:,:,:,:',   inout, _all)
VARN(put, 7D, EightByteReal, real, double, `:,:,:,:,:,:,:', inout, _all)

VARN(put, 1D, EightByteInt, integer, int8,  :,              inout, _all)
VARN(put, 2D, EightByteInt, integer, int8, `:,:',           inout, _all)
VARN(put, 3D, EightByteInt, integer, int8, `:,:,:',         inout, _all)
VARN(put, 4D, EightByteInt, integer, int8, `:,:,:,:',       inout, _all)
VARN(put, 5D, EightByteInt, integer, int8, `:,:,:,:,:',     inout, _all)
VARN(put, 6D, EightByteInt, integer, int8, `:,:,:,:,:,:',   inout, _all)
VARN(put, 7D, EightByteInt, integer, int8, `:,:,:,:,:,:,:', inout, _all)

!
! collective get APIs
!

VARN(get, 1D, OneByteInt, integer, int1,  :,              out, _all)
VARN(get, 2D, OneByteInt, integer, int1, `:,:',           out, _all)
VARN(get, 3D, OneByteInt, integer, int1, `:,:,:',         out, _all)
VARN(get, 4D, OneByteInt, integer, int1, `:,:,:,:',       out, _all)
VARN(get, 5D, OneByteInt, integer, int1, `:,:,:,:,:',     out, _all)
VARN(get, 6D, OneByteInt, integer, int1, `:,:,:,:,:,:',   out, _all)
VARN(get, 7D, OneByteInt, integer, int1, `:,:,:,:,:,:,:', out, _all)

VARN(get, 1D, TwoByteInt, integer, int2,  :,              out, _all)
VARN(get, 2D, TwoByteInt, integer, int2, `:,:',           out, _all)
VARN(get, 3D, TwoByteInt, integer, int2, `:,:,:',         out, _all)
VARN(get, 4D, TwoByteInt, integer, int2, `:,:,:,:',       out, _all)
VARN(get, 5D, TwoByteInt, integer, int2, `:,:,:,:,:',     out, _all)
VARN(get, 6D, TwoByteInt, integer, int2, `:,:,:,:,:,:',   out, _all)
VARN(get, 7D, TwoByteInt, integer, int2, `:,:,:,:,:,:,:', out, _all)

VARN(get, 1D, FourByteInt, integer, int,  :,              out, _all)
VARN(get, 2D, FourByteInt, integer, int, `:,:',           out, _all)
VARN(get, 3D, FourByteInt, integer, int, `:,:,:',         out, _all)
VARN(get, 4D, FourByteInt, integer, int, `:,:,:,:',       out, _all)
VARN(get, 5D, FourByteInt, integer, int, `:,:,:,:,:',     out, _all)
VARN(get, 6D, FourByteInt, integer, int, `:,:,:,:,:,:',   out, _all)
VARN(get, 7D, FourByteInt, integer, int, `:,:,:,:,:,:,:', out, _all)

VARN(get, 1D, FourByteReal, real,   real,  :,              out, _all)
VARN(get, 2D, FourByteReal, real,   real, `:,:',           out, _all)
VARN(get, 3D, FourByteReal, real,   real, `:,:,:',         out, _all)
VARN(get, 4D, FourByteReal, real,   real, `:,:,:,:',       out, _all)
VARN(get, 5D, FourByteReal, real,   real, `:,:,:,:,:',     out, _all)
VARN(get, 6D, FourByteReal, real,   real, `:,:,:,:,:,:',   out, _all)
VARN(get, 7D, FourByteReal, real,   real, `:,:,:,:,:,:,:', out, _all)

VARN(get, 1D, EightByteReal, real, double,  :,              out, _all)
VARN(get, 2D, EightByteReal, real, double, `:,:',           out, _all)
VARN(get, 3D, EightByteReal, real, double, `:,:,:',         out, _all)
VARN(get, 4D, EightByteReal, real, double, `:,:,:,:',       out, _all)
VARN(get, 5D, EightByteReal, real, double, `:,:,:,:,:',     out, _all)
VARN(get, 6D, EightByteReal, real, double, `:,:,:,:,:,:',   out, _all)
VARN(get, 7D, EightByteReal, real, double, `:,:,:,:,:,:,:', out, _all)

VARN(get, 1D, EightByteInt, integer, int8,  :,              out, _all)
VARN(get, 2D, EightByteInt, integer, int8, `:,:',           out, _all)
VARN(get, 3D, EightByteInt, integer, int8, `:,:,:',         out, _all)
VARN(get, 4D, EightByteInt, integer, int8, `:,:,:,:',       out, _all)
VARN(get, 5D, EightByteInt, integer, int8, `:,:,:,:,:',     out, _all)
VARN(get, 6D, EightByteInt, integer, int8, `:,:,:,:,:,:',   out, _all)
VARN(get, 7D, EightByteInt, integer, int8, `:,:,:,:,:,:,:', out, _all)

!
! text variable
!

dnl
dnl TEXTVARN1(ncid, varid, values, num, start, count)
dnl
define(`TEXTVARN1',dnl
`dnl
   function nf90mpi_$1_varn_text$3(ncid, varid, values, num, start, count)
     integer,                                        intent( in) :: ncid, varid, num
     character (len = *),                             intent($2) :: values
     integer (kind=MPI_OFFSET_KIND), dimension(:,:), intent( in) :: start, count
     integer                                                     :: nf90mpi_$1_varn_text$3
     nf90mpi_$1_varn_text$3 = nfmpi_$1_varn_text$3(ncid, varid, num, start, count, values)
   end function nf90mpi_$1_varn_text$3
')dnl

TEXTVARN1(put, in)
TEXTVARN1(get, out)

TEXTVARN1(put, in, _all)
TEXTVARN1(get, out, _all)

dnl
dnl TEXTVARN(ncid, varid, values, num, start, count)
dnl
define(`TEXTVARN',dnl
`dnl
   function nf90mpi_$1_varn_$2_text$6(ncid, varid, values, num, start, count)
     integer,                                        intent( in) :: ncid, varid, num
     character (len = *), dimension($3),             intent( $5) :: values
     integer (kind=MPI_OFFSET_KIND), dimension(:,:), intent( in) :: start, count
     integer                                                     :: nf90mpi_$1_varn_$2_text$6
 
     nf90mpi_$1_varn_$2_text$6 = nfmpi_$1_varn_text$6(ncid, varid, num, start, count, values($4))
   end function nf90mpi_$1_varn_$2_text$6
')dnl

TEXTVARN(put, 1D,  :,               1,              in)
TEXTVARN(put, 2D, `:,:',           `1,1',           in)
TEXTVARN(put, 3D, `:,:,:',         `1,1,1',         in)
TEXTVARN(put, 4D, `:,:,:,:',       `1,1,1,1',       in)
TEXTVARN(put, 5D, `:,:,:,:,:',     `1,1,1,1,1',     in)
TEXTVARN(put, 6D, `:,:,:,:,:,:',   `1,1,1,1,1,1',   in)
TEXTVARN(put, 7D, `:,:,:,:,:,:,:', `1,1,1,1,1,1,1', in)

TEXTVARN(get, 1D,  :,               1,              out)
TEXTVARN(get, 2D, `:,:',           `1,1',           out)
TEXTVARN(get, 3D, `:,:,:',         `1,1,1',         out)
TEXTVARN(get, 4D, `:,:,:,:',       `1,1,1,1',       out)
TEXTVARN(get, 5D, `:,:,:,:,:',     `1,1,1,1,1',     out)
TEXTVARN(get, 6D, `:,:,:,:,:,:',   `1,1,1,1,1,1',   out)
TEXTVARN(get, 7D, `:,:,:,:,:,:,:', `1,1,1,1,1,1,1', out)

!
! Collective APIs
!

TEXTVARN(put, 1D,  :,               1,              in, _all)
TEXTVARN(put, 2D, `:,:',           `1,1',           in, _all)
TEXTVARN(put, 3D, `:,:,:',         `1,1,1',         in, _all)
TEXTVARN(put, 4D, `:,:,:,:',       `1,1,1,1',       in, _all)
TEXTVARN(put, 5D, `:,:,:,:,:',     `1,1,1,1,1',     in, _all)
TEXTVARN(put, 6D, `:,:,:,:,:,:',   `1,1,1,1,1,1',   in, _all)
TEXTVARN(put, 7D, `:,:,:,:,:,:,:', `1,1,1,1,1,1,1', in, _all)

TEXTVARN(get, 1D,  :,               1,              out, _all)
TEXTVARN(get, 2D, `:,:',           `1,1',           out, _all)
TEXTVARN(get, 3D, `:,:,:',         `1,1,1',         out, _all)
TEXTVARN(get, 4D, `:,:,:,:',       `1,1,1,1',       out, _all)
TEXTVARN(get, 5D, `:,:,:,:,:',     `1,1,1,1,1',     out, _all)
TEXTVARN(get, 6D, `:,:,:,:,:,:',   `1,1,1,1,1,1',   out, _all)
TEXTVARN(get, 7D, `:,:,:,:,:,:,:', `1,1,1,1,1,1,1', out, _all)

