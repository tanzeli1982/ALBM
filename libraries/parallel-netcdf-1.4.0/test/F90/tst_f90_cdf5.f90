!
!  Copyright (C) 2013, Northwestern University and Argonne National Laboratory
!  See COPYRIGHT notice in top-level directory.
!
!     This is part of the PnetCDF package.
!
!     $Id: tst_f90_cdf5.f90 1468 2013-10-26 16:53:18Z wkliao $


program tst_f90_nc4
  use mpi
  use pnetcdf
  implicit none
  integer :: fh, cmode, err, dimid, varid, ndim, nvar
  character (len = *), parameter :: FILE_NAME = "tst_f90_nc4.nc"
  integer(KIND=MPI_OFFSET_KIND) :: ten=10
  character(LEN=128) filename, cmd, msg
  integer argc, iargc, my_rank, p

  call MPI_Init(err)
  call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, err)
  call MPI_Comm_size(MPI_COMM_WORLD, p, err)

  ! take filename from command-line argument if there is any
  call getarg(0, cmd)
  argc = IARGC()
  if (argc .GT. 1) then
     if (my_rank .EQ. 0) print*,'Usage: ',trim(cmd),' [filename]'
     goto 999
  endif
  filename = FILE_NAME
  if (argc .EQ. 1) call getarg(1, filename)

  if (p .ne. 1 .AND. my_rank .eq. 0) then
     print *, 'Warning: ',trim(cmd),' is design to run on 1 process'
  endif

  cmode = IOR(NF90_CLOBBER, NF90_64BIT_DATA)
  call check(nf90mpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, fh))
  call check(nf90mpi_def_dim(fh, 'fred', ten, dimid))
  call check(nf90mpi_def_var(fh, 'john', NF90_INT, (/dimid/), varid))
  call check(nf90mpi_close(fh))
  
  ! Check the file.
  call check(nf90mpi_open(MPI_COMM_WORLD, filename, NF90_WRITE, MPI_INFO_NULL, fh))
  call check(nf90mpi_inquire(fh, nDimensions = ndim, nVariables = nvar))
  if (nvar .ne. 1 .or. ndim .ne. 1) stop 3
  call check(nf90mpi_close(fh))

  call check(nf90mpi_create(MPI_COMM_WORLD, filename, cmode, MPI_INFO_NULL, fh))
  call check(nf90mpi_def_dim(fh, 'fred', ten, dimid))
  call check(nf90mpi_def_var(fh, 'john', NF90_INT, (/dimid/), varid))
  call check(nf90mpi_close(fh))
  
  ! Check the file.
  call check(nf90mpi_open(MPI_COMM_WORLD, filename, NF90_WRITE, MPI_INFO_NULL, fh))
  call check(nf90mpi_inquire(fh, nDimensions = ndim, nVariables = nvar))
  if (nvar .ne. 1 .or. ndim .ne. 1) stop 3
  call check(nf90mpi_close(fh))

   msg = '*** TESTING F90 '//trim(cmd)
   if (my_rank .eq. 0)   write(*,"(A67,A)") msg,'------ pass'

 999 call MPI_Finalize(err)

contains
!     This subroutine handles errors by printing an error message and
!     exiting with a non-zero status.
  subroutine check(errcode)
    implicit none
    integer, intent(in) :: errcode
    
    if(errcode /= nf90_noerr) then
       print *, 'Error: ', trim(nf90mpi_strerror(errcode))
       stop 2
    endif
  end subroutine check
end program tst_f90_nc4

