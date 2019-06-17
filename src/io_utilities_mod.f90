module io_utilities_mod
!---------------------------------------------------------------------------------
! Purpose:
!
! This module provides io-related utilities
!
!---------------------------------------------------------------------------------
   use shr_kind_mod,       only : cs => SHR_KIND_CS, cx => SHR_KIND_CX, r8, r4
   use shr_typedef_mod,    only : SimTime
   use shr_ctrl_mod,       only : archive_tstep, archive_dir, masterproc, DEBUG
   use math_utilities_mod, only : Mean
   use pnetcdf
   use mpi

   interface Endrun
      module procedure Endrun1arg
      module procedure EndrunInt
      module procedure EndrunStr
   end interface

   interface MstrMsg
      module procedure MstrMsgInt
      module procedure MstrMsgStr
   end interface

   interface check
      module procedure check1arg
      module procedure checkstr
   end interface

   interface GetArchiveFullname
      module procedure GetStaticArchiveFullname
      module procedure GetDynamicArchiveFullname
   end interface

   interface WriteData
      module procedure WriteReal4Data1D
      module procedure WriteReal4Data2D
      module procedure WriteReal4Data3D
      module procedure WriteStaticReal8Data1D
   end interface

   interface DefNcVariable
      module procedure DefNcVariable4R
      module procedure DefNcVariable8R
   end interface
 
contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: Get full file name
   !
   !------------------------------------------------------------------------------
   subroutine GetFullFileName(filename, fullname)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), intent(out) :: fullname
      character(cs) :: prefix
      integer :: indx

      if (filename(1:2) == '~/') then
         call get_environment_variable('HOME',prefix)
         fullname = trim(prefix) // filename(2:len_trim(filename))
      else if (filename(1:2) == './') then
         call GETCWD(prefix)
         fullname = trim(prefix) // filename(2:len_trim(filename))
      else if (filename(1:1) == '$') then
         indx = INDEX(filename, '/')
         call get_environment_variable(filename(2:indx-1),prefix)
         fullname = trim(prefix) // filename(indx:len_trim(filename))
      else
         fullname = filename
      end if
   end subroutine

   subroutine AppendFilename(dirname, filename, fullname)
      implicit none
      character(len=*), intent(in) :: dirname
      character(len=*), intent(in) :: filename
      character(cx), intent(out) :: fullname
      character(cs) :: prefix
      integer :: indx

      if (dirname(1:2) == '~/') then
         call get_environment_variable('HOME',prefix)
         fullname = trim(prefix) // dirname(2:len_trim(dirname)) &
            // trim(filename)
      else if (dirname(1:2) == './') then
         call GETCWD(prefix)
         fullname = trim(prefix) // dirname(2:len_trim(dirname)) &
            // trim(filename)
      else if (dirname(1:1) == '$') then
         indx = INDEX(dirname, '/')
         call get_environment_variable(dirname(2:indx-1),prefix)
         fullname = trim(prefix) // dirname(indx:len_trim(dirname)) &
            // trim(filename)
      else
         fullname = trim(dirname) // trim(filename)
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: End mpi execution and print error message
   !
   !------------------------------------------------------------------------------
   subroutine Endrun1arg(msg)
      implicit none
      character(len=*), intent(in) :: msg
      integer :: err

      print *, "Error: ", trim(msg), "!!"
      call MPI_ABORT(MPI_COMM_WORLD, 1, err)
   end subroutine

   subroutine EndrunInt(id, msg)
      implicit none
      integer, intent(in) :: id
      character(len=*), intent(in) :: msg
      integer :: err

      print "(A, I0, A)", "Error in Object ", id, ": " // trim(msg) // "!!"
      call MPI_ABORT(MPI_COMM_WORLD, 1, err)
   end subroutine

   subroutine EndrunStr(fid, msg)
      implicit none
      character(len=*), intent(in) :: fid
      character(len=*), intent(in) :: msg
      integer :: err

      print *, "Error in File " // trim(fid) // ": " // trim(msg) // "!!"
      call MPI_ABORT(MPI_COMM_WORLD, 1, err)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: print warning message by master process
   !
   !------------------------------------------------------------------------------
   subroutine MstrMsgInt(id, msg)
      implicit none
      integer, intent(in) :: id
      character(len=*), intent(in) :: msg

      print "(A, I0, A)", "Warning in Object ", id, ": " // trim(msg) // "!!"
   end subroutine

   subroutine MstrMsgStr(fid, msg)
      implicit none
      character(len=*), intent(in) :: fid
      character(len=*), intent(in) :: msg

      print *, "Warning in File " // trim(fid) //  ": " // trim(msg) // "!!"
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Check the status of netcdf-related function 
   !
   !------------------------------------------------------------------------------
   subroutine check1arg(istatus)
      implicit none
      integer, intent(in) :: istatus

      if(istatus /= NF90_NOERR) then
         call Endrun(nf90mpi_strerror(istatus)) 
      end if
   end subroutine

   subroutine checkstr(fid, istatus)
      implicit none
      character(len=*), intent(in) :: fid
      integer, intent(in) :: istatus

      if(istatus /= NF90_NOERR) then
         call Endrun(fid, nf90mpi_strerror(istatus))
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Convert a string to lower case or upper case
   !
   !------------------------------------------------------------------------------
   subroutine To_upper(str_in, str_out)
      implicit none
      character(len=*), intent(in) :: str_in
      character(len=*), intent(out) :: str_out
      integer :: i
 
      do i = 1, len_trim(str_in)
         select case(str_in(i:i))
            case("a":"z")
               str_out(i:i) = achar(iachar(str_in(i:i))-32)
            case default
               str_out(i:i) = str_in(i:i)
         end select
      end do 
   end subroutine
 
   subroutine To_lower(str_in, str_out)
      implicit none
      character(len=*), intent(in) :: str_in
      character(len=*), intent(out) :: str_out
      integer :: i
 
      do i = 1, len_trim(str_in)
         select case(str_in(i:i))
            case("A":"Z")
               str_out(i:i) = achar(iachar(str_in(i:i))+32)
            case default
               str_out(i:i) = str_in(i:i)
         end select
      end do  
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: form full file name for archive and backup files
   !
   !------------------------------------------------------------------------------
   subroutine GetStaticArchiveFullname(varname, fullname)
      implicit none
      character(len=*), intent(in) :: varname
      character(len=*), intent(out) :: fullname
      character(cx) :: fulldir, command
      character(cs) :: tmpstr
      logical :: isexist

      call GetFullFileName(archive_dir, fulldir)
      inquire(directory=trim(fulldir), exist=isexist)
      if (.not. isexist) then
         command = 'mkdir ' // fulldir
         call system(trim(command))
      end if
      fullname = trim(fulldir) // 'bLakeOut.' // trim(varname) // '.nc'
   end subroutine

   subroutine GetDynamicArchiveFullname(time, varname, fullname)
      implicit none
      Type(SimTime), intent(in) :: time
      character(len=*), intent(in) :: varname
      character(len=*), intent(out) :: fullname
      character(cx) :: fulldir, command
      character(cs) :: tmpstr
      logical :: isexist

      write(tmpstr,"(I4, I2.2, I2.2, A, I4, I2.2, I2.2)") time%year0, &
            time%month0, time%day0, '_', time%year1, time%month1, time%day1
      call GetFullFileName(archive_dir, fulldir)
      inquire(directory=trim(fulldir), exist=isexist)
      if (.not. isexist) then
         command = 'mkdir ' // fulldir
         call system(trim(command))
      end if
      fullname = trim(fulldir) // 'bLakeOut.' // trim(varname) // '.' // &
         trim(tmpstr) // '.nc'
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: define a variable in nc file
   !
   !------------------------------------------------------------------------------
   subroutine DefNcVariable8R(ncid, dimids, varname, longname, units, &
                              fillvalue, missvalue, varid)
      implicit none
      integer, intent(in) :: ncid
      integer, dimension(:), intent(in) :: dimids
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units
      real(r8), intent(in) :: fillvalue
      real(r8), intent(in) :: missvalue
      integer, intent(out) :: varid
      real(r8) :: fval, mval

      fval = fillvalue
      mval = missvalue
      call check( nf90mpi_def_var(ncid, trim(varname), NF90_DOUBLE, &
                  dimids, varid) )
      call check( nf90mpi_put_att(ncid, varid, "long_name", trim(longname)) )
      call check( nf90mpi_put_att(ncid, varid, "units", trim(units)) )
      call check( nf90mpi_put_att(ncid, varid, "_FillValue", fval) )
      call check( nf90mpi_put_att(ncid, varid, "missing_value", mval) )
   end subroutine
   
   subroutine DefNcVariable4R(ncid, dimids, varname, longname, units, &
                              fillvalue, missvalue, varid)
      implicit none
      integer, intent(in) :: ncid
      integer, dimension(:), intent(in) :: dimids
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units
      real(r4), intent(in) :: fillvalue
      real(r4), intent(in) :: missvalue
      integer, intent(out) :: varid
      real(r4) :: fval, mval

      fval = fillvalue
      mval = missvalue
      call check( nf90mpi_def_var(ncid, trim(varname), NF90_FLOAT, &
                  dimids, varid) )
      call check( nf90mpi_put_att(ncid, varid, "long_name", trim(longname)) )
      call check( nf90mpi_put_att(ncid, varid, "units", trim(units)) )
      call check( nf90mpi_put_att(ncid, varid, "_FillValue", fval) )
      call check( nf90mpi_put_att(ncid, varid, "missing_value", mval) )
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: write data in collective mode
   !
   !------------------------------------------------------------------------------
   subroutine WriteStaticReal8Data1D(lakeid, varname, odata)
      implicit none
      integer, intent(in) :: lakeid
      character(len=*), intent(in) :: varname
      real(r8), intent(inout) :: odata(:)
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND) :: nstart(2), ncount(2)
      integer :: ncid, varid, nlayer

      call GetArchiveFullname(varname, fullname)
      call check( nf90mpi_open(MPI_COMM_WORLD, trim(fullname), NF90_WRITE, &
                  MPI_INFO_NULL, ncid) )
      nlayer = size(odata)
      nstart = (/1, lakeid/)
      ncount = (/nlayer, 1/)
      call check( nf90mpi_inq_varid(ncid, trim(varname), varid) )
      call check( nf90mpi_put_var_all(ncid, varid, odata, nstart, ncount) )
      call check( nf90mpi_close(ncid) )
      if (masterproc .and. DEBUG) then
         print *, "Write Real8 1-D variable " // trim(varname)
      end if
   end subroutine

   subroutine WriteReal4Data1D(lakeid, time, varname, odata)
      implicit none
      integer, intent(in) :: lakeid
      type(SimTime), intent(in) :: time
      character(len=*), intent(in) :: varname
      real(r4), intent(inout) :: odata(:)
      real(r4), allocatable :: odata1d(:)
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND) :: nstart(2), ncount(2)
      integer :: ncid, varid, nt, ii

      call GetArchiveFullname(time, varname, fullname)
      if (trim(archive_tstep)=='day') then
         nt = INT( size(odata)/24 )
         allocate(odata1d(nt))
         do ii = 1, nt, 1
            call Mean(odata(24*ii-23:24*ii), odata1d(ii))
         end do
         nstart = (/1, lakeid/)
         ncount = (/nt, 1/)
         call check( nf90mpi_open(MPI_COMM_WORLD, trim(fullname), & 
                     NF90_WRITE, MPI_INFO_NULL, ncid) )
         call check( nf90mpi_inq_varid(ncid, trim(varname), varid) )
         call check( nf90mpi_put_var_all(ncid, varid, odata1d, &
                     nstart, ncount) )
         call check( nf90mpi_close(ncid) )
         deallocate(odata1d)
      else if (trim(archive_tstep)=='hour') then
         nt = size(odata)
         nstart = (/1, lakeid/)
         ncount = (/nt, 1/)
         call check( nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                     NF90_WRITE, MPI_INFO_NULL, ncid) )
         call check( nf90mpi_inq_varid(ncid, trim(varname), varid) )
         call check( nf90mpi_put_var_all(ncid, varid, odata, &
                     nstart, ncount) )
         call check( nf90mpi_close(ncid) )
      end if 
      if (masterproc .and. DEBUG) then
         print *, "Write Real4 1-D variable " // trim(varname)
      end if
   end subroutine

   subroutine WriteReal4Data2D(lakeid, time, varname, odata)
      implicit none
      integer, intent(in) :: lakeid
      type(SimTime), intent(in) :: time
      character(len=*), intent(in) :: varname
      real(r4), intent(inout) :: odata(:,:)
      real(r4), allocatable :: odata2d(:,:)
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND) :: nstart(3), ncount(3)
      integer :: ncid, varid, nt, nlayer, ii

      nlayer = size(odata,1)
      call GetArchiveFullname(time, varname, fullname)
      if (trim(archive_tstep)=='day') then
         nt = INT( size(odata,2)/24 )
         allocate(odata2d(nlayer,nt))
         do ii = 1, nt, 1
            call Mean(odata(:,24*ii-23:24*ii), 2, odata2d(:,ii))
         end do
         nstart = (/1, 1, lakeid/)
         ncount = (/nlayer, nt, 1/)
         call check( nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                     NF90_WRITE, MPI_INFO_NULL, ncid) )
         call check( nf90mpi_inq_varid(ncid, trim(varname), varid) )
         call check( nf90mpi_put_var_all(ncid, varid, odata2d, &
                     nstart, ncount) )
         call check( nf90mpi_close(ncid) ) 
         deallocate(odata2d)
      else if (trim(archive_tstep)=='hour') then
         nt = size(odata,2)
         nstart = (/1, lakeid, 1/)
         ncount = (/nlayer, 1, nt/)
         call check( nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                     NF90_WRITE, MPI_INFO_NULL, ncid) )
         call check( nf90mpi_inq_varid(ncid, trim(varname), varid) )
         call check( nf90mpi_put_var_all(ncid, varid, odata, &
                     nstart, ncount) )
         call check( nf90mpi_close(ncid) )
      end if
      if (masterproc .and. DEBUG) then
         print *, "Write Real4 2-D variable " // trim(varname)
      end if
   end subroutine

   subroutine WriteReal4Data3D(lakeid, time, varname, odata)
      implicit none
      integer, intent(in) :: lakeid
      type(SimTime), intent(in) :: time
      character(len=*), intent(in) :: varname
      real(r4), intent(inout) :: odata(:,:,:)
      real(r4), allocatable :: odata3d(:,:,:)
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND) :: nstart(4), ncount(4)
      integer :: ncid, varid, nt, ng, nlayer, ii

      ng = size(odata,1)
      nlayer = size(odata,2)
      if (trim(archive_tstep)=='day') then
         nt = INT( size(odata,3)/24 )
         allocate(odata3d(ng,nlayer,nt))
         do ii = 1, nt, 1
            call Mean(odata(:,:,24*ii-23:24*ii), 3, odata3d(:,:,ii))
         end do
         nstart = (/1, 1, 1, lakeid/)
         ncount = (/ng, nlayer, nt, 1/)
         call check( nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                     NF90_WRITE, MPI_INFO_NULL, ncid) )
         call check( nf90mpi_inq_varid(ncid, trim(varname), varid) )
         call check( nf90mpi_put_var_all(ncid, varid, odata3d, &
                     nstart, ncount) )
         call check( nf90mpi_close(ncid) )
         deallocate(odata3d)
      else if (trim(archive_tstep)=='hour') then
         nt = size(odata,3)
         nstart = (/1, 1, 1, lakeid/)
         ncount = (/ng, nlayer, nt, 1/)
         call check( nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                     NF90_WRITE, MPI_INFO_NULL, ncid) )
         call check( nf90mpi_inq_varid(ncid, trim(varname), varid) )
         call check( nf90mpi_put_var_all(ncid, varid, odata, &
                     nstart, ncount) )
         call check( nf90mpi_close(ncid) )
      end if
      if (masterproc .and. DEBUG) then
         print *, "Write Real4 3-D variable " // trim(varname)
      end if
   end subroutine
end module io_utilities_mod
