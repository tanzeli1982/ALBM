module radiation_io_mod
!---------------------------------------------------------------------------------
! Purpose: I/O utilities for SMARTS model
!---------------------------------------------------------------------------------
   use shr_kind_mod,       only : cs => SHR_KIND_CS, cx => SHR_KIND_CX, r8, r4
   use shr_ctrl_mod,       only : solar_dir, gas_dir, albedo_dir 
   use io_utilities_mod
   use math_utilities_mod, only : Interp1d

   public
   private :: fid
   integer, parameter :: fid = 15   ! file identifier

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: read spectrum file
   !
   !------------------------------------------------------------------------------
   subroutine ReadSpctrmFile(ispctr, ESC, vWVLN, vH0)
      implicit none
      integer, intent(in) :: ispctr          ! spectrum file specifier
      real(r8), intent(out) :: ESC           ! solar constant
      real(r8), intent(out) :: vWVLN(:)      ! spectrum vector (nm)
      real(r8), intent(out) :: vH0(:)        ! spectrum coefficient vector 
      character(len=32) :: fname = "ReadSpctrmFile"
      character(cx) :: fullname, filename
      character(cs) :: tempstr
      integer :: nline, error, ii, IWVLN1
      
      if (ispctr>=0) then
         write(filename, "(A, I0, A)") 'Spctrm_', ispctr, '.dat'
         filename = trim(solar_dir) // trim(filename)
      else
         filename = trim(solar_dir) // 'Spctrm_U.dat'
      end if
      call GetFullFileName(filename, fullname)
      open(unit=fid, file=fullname, status="old", action="read", iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun(fname, 'Cannot open file ' // trim(fullname))
      end if
      read(fid, "(A80)", iostat=error) tempstr
      read(unit=fid, fmt=*, iostat=error) ESC
      if (error==0) then
         nline = size(vWVLN)
         do ii = 1, nline, 1
            read(unit=fid, fmt=*, iostat=error) IWVLN1, vH0(ii)
            if (error/=0) exit
            vWVLN(ii) = DBLE(IWVLN1) / 10.0
         end do
      end if
      close(unit=fid)
      if (error/=0) then
         call Endrun(fname, 'Error in reading ' // trim(fullname))
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: read gas absorption band file
   !
   !------------------------------------------------------------------------------
   subroutine ReadGasFile(gstr, coefs)
      implicit none
      character(len=*), intent(in) :: gstr
      real(r8), intent(out) :: coefs(:,:)
      character(len=32) :: fname = "ReadGasFile"
      character(cx) :: fullname, filename
      character(cx) :: tempstr
      integer :: nline, error, ii

      filename = trim(gas_dir) // trim(gstr) // '.dat'
      call GetFullFileName(filename, fullname)
      open(unit=fid, file=fullname, status="old", action="read", iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun(fname, 'Cannot open file ' // trim(fullname))
      end if
      read(fid, "(A512)", iostat=error) tempstr 
      nline = size(coefs, 1)
      do ii = 1, nline, 1
         read(unit=fid, fmt=*, iostat=error) coefs(ii,:)
         if (error/=0) exit
      end do
      close(unit=fid)
      if (error/=0) then
         call Endrun(fname, 'Error in reading ' // trim(fullname))
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: read albedo band file.
   !        There are albedo band files for other surface cover types, see
   !        /home/tan80/packages/SMARTS_295_Linux/Albedo for details. 
   !
   !------------------------------------------------------------------------------
   subroutine ReadAlbedoFile(snowstr, vwvln, albedos)
      implicit none
      character(len=*), intent(in) :: snowstr
      real(r8), intent(in) :: vwvln(:)          ! nm
      real(r8), intent(out) :: albedos(:)
      character(len=32) :: fname = "ReadAlbedoFile"
      character(cx) :: fullname, filename
      real(r8), allocatable :: rawdata(:,:)
      real(r8) :: pair(2)
      integer :: nline, nwvln, error, ii

      filename = trim(albedo_dir) // trim(snowstr) // '.dat'
      call GetFullFileName(filename, fullname)
      open(unit=fid, file=fullname, status="old", action="read", iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun(fname, 'Cannot open file ' // trim(fullname))
      end if
      select case(trim(snowstr))
         case("FineSnow")
            nline = 1388
         case("MeltSnow")
            nline = 2151
      end select
      allocate(rawdata(nline,2)) 
      do ii = 1, nline, 1
         read(unit=fid, fmt=*, iostat=error) rawdata(ii,:)
         if (error/=0) exit
      end do
      close(unit=fid)
      if (error/=0) then
         call Endrun(fname, 'Error in reading ' // trim(fullname))
      end if
      ! The units of rawdata(:,1) is mm
      call Interp1d(rawdata(:,1),rawdata(:,2),1e-3*vwvln,albedos)
      deallocate(rawdata)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: write irradiance at earth surface
   !
   !------------------------------------------------------------------------------
   subroutine WriteIrradianceFile(vWVLN, vFRAD, vDFRAC)
      implicit none
      real(r8), intent(in) :: vWVLN(:)    ! nm
      real(r8), intent(in) :: vFRAD(:)    ! photons m-2 s-1
      real(r8), intent(in) :: vDFRAC(:)   ! diffuse fraction
      character(len=32) :: fname = "WriteIrradianceFile"
      character(cx) :: filename, fullname
      integer :: ii, nwvln, error

      filename = "./smarts295.ext.txt"
      call GetFullFileName(filename, fullname)
      open(unit=fid, file=fullname, status="replace", action="write", iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun(fname, 'Cannot create/open file ' // trim(fullname))
      end if
      write(unit=fid, fmt=*, iostat=error) "Wvlgth,Beam_normal_photon_irrad,", &
                        "Difuse_horiz_photn_irrad"
      nwvln = size(vWVLN)
      do ii = 1, nwvln, 1
         write(unit=fid, fmt='(E12.4E3, E12.4E3, E12.4E3)', iostat=error) vWVLN(ii), &
                  vFRAD(ii)*(1.0-vDFRAC(ii)), vFRAD(ii)*vDFRAC(ii)
         if (error/=0) exit
      end do
      close(unit=fid)
      if (error/=0) then
         call Endrun(fname, 'Error in writing ' // trim(fullname))
      end if
   end subroutine

end module radiation_io_mod
