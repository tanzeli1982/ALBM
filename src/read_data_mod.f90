module read_data_mod
!---------------------------------------------------------------------------------
! Purpose:
!
! This module serves for the data preparation and storage
! ECMWF data reanalysis: http://stommel.tamu.edu/~baum/ecmwf.html
!
!---------------------------------------------------------------------------------
   use shr_kind_mod,       only : cs => SHR_KIND_CS, cx => SHR_KIND_CX, &
                                  r8, r4, i8, i4, i1
   use shr_typedef_mod,    only : LakeInfo, SimTime
   use shr_ctrl_mod
   use shr_param_mod
   use io_utilities_mod
   use math_utilities_mod
   use phy_utilities_mod
   use bg_const_mod

   public
   private :: fid
   integer, parameter :: fid = 10   ! file identifier
   integer, parameter :: Ncorigin = 000, Ncinterp = 001

   interface ReadStaticData
      module procedure ReadStaticRealData
      module procedure ReadStaticByteData
   end interface

   interface ReadLakeInfoData
      module procedure ReadLakeInfoRealData
      module procedure ReadLakeInfoIntData
   end interface

   interface CreateOutputFile
      module procedure CreateOutputFile0d
      module procedure CreateOutputFile1d
      module procedure CreateOutputFile2d
      module procedure CreateOutputFile3d
   end interface

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: Read lake information (lat, lon, depth, area) and calculate discrete
   !          water layers
   !
   !------------------------------------------------------------------------------
   subroutine ReadLakeInfo(lakeId, sampleId)
      implicit none
      integer, intent(in) :: lakeId
      integer, intent(in) :: sampleId
      real(r8) :: area, elev, eta

      ! read lake file
      call ReadLakeInfoData(lake_file, 'type', lakeId, lake_info%itype)
      call ReadLakeInfoData(lake_file, 'lon', lakeId, lake_info%longitude) ! deg
      call ReadLakeInfoData(lake_file, 'lat', lakeId, lake_info%latitude)  ! deg
      call ReadLakeInfoData(lake_file, 'depth', lakeId, lake_info%depth)   ! m
      call ReadLakeInfoData(lake_file, 'maxdepth', lakeId, lake_info%maxdepth) ! m
      call ReadLakeInfoData(lake_file, 'Asurf', lakeId, area)  ! km^2
      call ReadLakeInfoData(lake_file, 'alt', lakeId, elev)    ! m
      call ReadLakeInfoData(lake_file, 'wrt', lakeId, lake_info%wrt)    ! day
      call ReadLakeInfoData(lake_file, 'pH', lakeId, lake_info%pH)      !
      call ReadLakeInfoData(lake_file, 'sal', lakeId, lake_info%sal)    !
      call ReadLakeInfoData(lake_file, 'SRP', lakeId, lake_info%srp)    ! g/m^3
      call ReadLakeInfoData(lake_file, 'Fe', lakeId, lake_info%sedFe)   ! g/m^3
      call ReadLakeInfoData(lake_file, 'csed', lakeId, lake_info%csed)  ! kgC/m^3
      call ReadLakeInfoData(lake_file, 'cdep', lakeId, lake_info%cdep)  ! gC/m^2/yr
      call ReadLakeInfoData(lake_file, 'ice', lakeId, lake_info%excice) ! fraction
      call ReadLakeInfoData(lake_file, 'hsed', lakeId, lake_info%hsed)  ! m
      call ReadLakeInfoData(lake_file, 'bfdep', lakeId, lake_info%bfdep)  ! m
      call ReadLakeInfoData(lake_file, 'thrmkst', lakeId, lake_info%thrmkst)
      call ReadLakeInfoData(lake_file, 'eta', lakeId, eta) ! m^-1
      call ReadLakeInfoData(lake_file, 'refPOC', lakeId, lake_info%refPOC) ! gC/m3
      lake_info%id = lakeId
      lake_info%sampleId = sampleId
      lake_info%Asurf = 1.d6 * area  ! m^2
      lake_info%zalt = 1.d-3 * elev  ! km
      lake_info%kext = eta
      lake_info%bfdep = max(0.1_r8, lake_info%bfdep)
      if (lake_info%maxdepth<=0.1*NWLAYER) then
         WATER_LAYER = INT(10 * lake_info%maxdepth)
      else
         WATER_LAYER = NWLAYER
      end if
      if (NINT(lake_info%hsed)==10) then
         SED_LAYER = 20
      else
         SED_LAYER = NSLAYER
      end if
      if (lake_info%itype==yedoma .and. lake_info%thrmkst/=2) then
         call Endrun('lake_info%itype and lake_info%thrmkst must be both yedoma') 
      end if
      if (DEBUG) then
         print *, "Lake information: ", lake_info
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Read lake parameter samples for Monte Carlo sensitivity analysis
   !
   !------------------------------------------------------------------------------
   subroutine ReadParameterSamples(nmaxsample, samples)
      implicit none
      integer, intent(in) :: nmaxsample
      real(r8), intent(out) :: samples(:,:)
      character(len=32) :: fname = "ReadParameterSamples"
      character(cx) :: fullname, msg
      integer :: ii, id, error

      call GetFullFileName(mc_file, fullname)
      open(unit=fid, file=fullname, status="old", action="read", iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun('Cannot open file ' // trim(fullname))
      end if
      do ii = 1, nmaxsample, 1
         read(unit=fid, fmt=*, iostat=error) id, samples(ii,:)
         if (error/=0) then
            close(unit=fid)
            write(msg, "(A, I0)") "Reading stops at line ", ii
            call Endrun(fname, msg)
         end if
      end do
      close(unit=fid)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Read soil column active flags for ebullition
   !
   !------------------------------------------------------------------------------
   !subroutine ReadSoilColBubFlag(soilColBub)
   !   implicit none
   !   integer, intent(out) :: soilColBub(:)
   !   integer :: nflag
   !   integer :: pos2, pos1
   !
   !   soilColBub = 1 ! default values
   !   nflag = 0
   !   pos1 = 1
   !   do while (.True.)
   !      pos2 = index(Bubble_Cols(pos1:), ",")
   !      nflag = nflag + 1
   !      if (pos2==0) then
   !         exit
   !      end if
   !      pos1 = pos2 + pos1
   !   end do
   !   if (nflag>NSCOL) then
   !      call Endrun("The size of Bubble_Cols cannot exceed NSCOL")
   !   end if
   !   read(Bubble_Cols,fmt=*) soilColBub(1:nflag)
   !end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Read lake bathymetry 
   !
   !------------------------------------------------------------------------------
   subroutine ReadLakeBathymetry(info, Zw, dZw, Az, dAz, SAL, soilColInd)
      implicit none
      type(LakeInfo), intent(inout) :: info
      real(r8), intent(in) :: Zw(:)
      real(r8), intent(in) :: dZw(:)
      real(r8), intent(out) :: Az(:)
      real(r8), intent(out) :: dAz(:)
      real(r8), intent(out) :: SAL(:)
      integer, intent(out) :: soilColInd(:)
      character(len=32) :: fname = "ReadLakeBathymetry"
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND) :: nstart(2)
      integer(kind=MPI_OFFSET_KIND) :: ncount(2)
      integer(kind=MPI_OFFSET_KIND) :: nlevel
      integer :: ncid, dimid
      integer :: h_varid, A_varid, S_varid
      integer :: ii, nz, indx, nLvalid
      real(r8), allocatable :: tmph(:,:)
      real(r8), allocatable :: tmpA(:,:)
      real(r8), allocatable :: tmpS(:,:)
      real(r8), allocatable :: tmph_corr(:)
      real(r8), allocatable :: tmpA_corr(:)
      real(r8), allocatable :: tmpS_corr(:)
      real(r8), allocatable :: tmpZw(:)
      real(r8), allocatable :: tmpAz(:)
      real(r8), allocatable :: tmpSAL(:)

      call GetFullFileName(bathy_file, fullname) 
      call check( fname, nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                  NF90_NOWRITE, MPI_INFO_NULL, ncid) )
      call check( fname, nf90mpi_inq_dimid(ncid, 'level', dimid) )
      call check( nf90mpi_inquire_dimension(ncid, dimid, len = nlevel) )
      call check( nf90mpi_inq_varid(ncid, 'h', h_varid) )
      call check( nf90mpi_inq_varid(ncid, 'A', A_varid) )
      call check( nf90mpi_inq_varid(ncid, 'sal', S_varid) )
      nstart = (/1, info%id/)
      ncount = (/nlevel, 1_i8/)
      allocate(tmph(nlevel,1))
      allocate(tmpA(nlevel,1))
      allocate(tmpS(nlevel,1))
      call check( nf90mpi_get_var_all(ncid, h_varid, tmph, nstart, ncount) )
      call check( nf90mpi_get_var_all(ncid, A_varid, tmpA, nstart, ncount) )
      call check( nf90mpi_get_var_all(ncid, S_varid, tmpS, nstart, ncount) )
      call check( nf90mpi_close(ncid) )

      ! remove all zero area above the maximum depth
      allocate(tmph_corr(nlevel))
      allocate(tmpA_corr(nlevel))
      allocate(tmpS_corr(nlevel))
      nLvalid = 0
      do ii = 1, nlevel, 1
         if (tmpA(ii,1)>e8 .or. abs(tmph(ii,1)-info%maxdepth)<e8) then
            nLvalid = nLvalid + 1
            tmph_corr(nLvalid) = tmph(ii,1)
            tmpA_corr(nLvalid) = tmpA(ii,1)
            tmpS_corr(nLvalid) = tmpS(ii,1)
         end if
      end do

      ! interpolate and link soil column
      nz = size(Zw)
      allocate(tmpZw(nz+1))
      allocate(tmpAz(nz+1))
      allocate(tmpSAL(nz+1))
      do ii = 1, nz+1, 1
         if (ii==1) then
            tmpZw(ii) = Zw(ii)
         else if (ii==nz) then
            tmpZw(ii) = Zw(ii) - dZw(ii)
         else if (ii==nz+1) then
            tmpZw(ii) = Zw(ii-1)
         else
            tmpZw(ii) = Zw(ii) - 0.5*dZw(ii)
         end if
      end do

      call Interp1d(tmph_corr(1:nLvalid), sqrt(tmpA_corr(1:nLvalid)), &
               tmpZw, tmpAz)
      call Interp1d(tmph_corr(1:nLvalid), tmpS_corr(1:nLvalid), tmpZw, tmpSAL)
      tmpAz = 1.d6 * tmpAz**2 ! km^2 > m^2
      Az = tmpAz(1:nz)
      SAL = tmpSAL(1:nz)
      do ii = 1, nz, 1
         if (ii<nz) then
            dAz(ii) = max( Az(ii)-Az(ii+1), 0._r8 )
         else
            dAz(ii) = max( Az(ii), 0._r8 )
         end if
         indx = NSCOL - INT(NSCOL*Az(ii)/Az(1))
         indx = max(min(indx,NSCOL),1)
         soilColInd(ii) = indx
      end do

      deallocate(tmph, tmpA, tmpS)
      deallocate(tmph_corr, tmpA_corr, tmpS_corr)
      deallocate(tmpZw, tmpAz, tmpSAL)
      if (DEBUG .and. masterproc) then
         print *, "Read bathymetry from " // trim(fullname)
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Read all model simulation settings
   !
   !------------------------------------------------------------------------------
   subroutine ReadSimulationSettings(filename)
      implicit none
      character(len=*), intent(in) :: filename
      character(cx) :: fullname
      integer :: error
      namelist /general/ run_mode, lake_file, lake_range, bathy_file, opt_file
      namelist /simulation/ Thermal_Module, Bubble_Module, Diagenesis_Module, &
                            Carbon_Module, Start_Year, Start_Month,  &
                            Start_Day, End_Year, End_Month, End_Day, &
                            Spinup_Month, Spinup_Day, nSpinup, Use_Leap
      namelist /resolution/ NWLAYER, NSLAYER, NRLAYER, NSCOL 
      namelist /bayesian/ NMAXSAMPLE, sample_range, mc_file, SA_Start_Year, &
                          SA_Start_Month, SA_Start_Day, SA_End_Year, &
                          SA_End_Month, SA_End_Day
      namelist /radiation/ solar_dir, gas_dir, albedo_dir, co2_file, &
                           o3_file, aod_file
      namelist /rundata/ forcing_tstep, forcing_time, tas_file, tasmax_file, &
                         tasmin_file, hurs_file, ps_file, pr_file, prsn_file, &
                         rsds_file, rlds_file, wind_file, hydro_file, &
                         tref_file, soc_file, veg_file, wlnd_file
      namelist /archive/ archive_tstep, archive_dir
      namelist /dbg/ DEBUG
      
      call GetFullFileName(filename, fullname)
      open(unit=fid, file=fullname, status="old", action="read", &
           delim='APOSTROPHE', iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun('Cannot open file ' // trim(fullname))
      end if
      read(unit=fid, NML=general, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("reading general group")
      end if
      read(unit=fid, NML=simulation, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("reading simulation group")
      end if
      read(unit=fid, NML=resolution, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("reading resolution group")
      end if
      read(unit=fid, NML=bayesian, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("reading bayesian group")
      end if
      read(unit=fid, NML=radiation, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("reading radiation group")
      end if
      read(unit=fid, NML=rundata, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("reading rundata group")
      end if
      read(unit=fid, NML=archive, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("reading archive group")
      end if
      read(unit=fid, NML=dbg, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("reading debug group")
      end if
      close(unit=fid)

      call To_lower(run_mode, run_mode)
      call To_upper(forcing_time, forcing_time)
      call To_lower(archive_tstep, archive_tstep)
   end subroutine

   subroutine BcastSimulationSettings()
      implicit none
      integer :: err

      ! general group
      call MPI_BCAST(run_mode, len(run_mode), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(lake_file, len(lake_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(opt_file, len(opt_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(bathy_file, len(bathy_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(lake_range, size(lake_range), MPI_INTEGER, 0, &
                     MPI_COMM_WORLD, err)
      ! simulation group
      call MPI_BCAST(Thermal_Module, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Bubble_Module, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Diagenesis_Module, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Carbon_Module, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Start_Year, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Start_Month, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Start_Day, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(End_Year, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(End_Month, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(End_Day, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Spinup_Month, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Spinup_Day, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(nSpinup, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Use_Leap, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
      ! resolution group
      call MPI_BCAST(NWLAYER, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(NSLAYER, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(NRLAYER, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(NSCOL, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      ! bayesian group
      call MPI_BCAST(NMAXSAMPLE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(sample_range, size(sample_range), MPI_INTEGER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(mc_file, len(mc_file), MPI_CHARACTER, 0, MPI_COMM_WORLD, &
                     err)
      call MPI_BCAST(SA_Start_Year, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(SA_Start_Month, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(SA_Start_Day, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(SA_End_Year, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(SA_End_Month, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(SA_End_Day, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      ! radiation group
      call MPI_BCAST(solar_dir, len(solar_dir), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(gas_dir, len(gas_dir), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(albedo_dir, len(gas_dir), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(co2_file, len(co2_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(o3_file, len(o3_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(aod_file, len(aod_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      ! rundata group
      call MPI_BCAST(forcing_tstep, len(forcing_tstep), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(forcing_time, len(forcing_time), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(tas_file, len(tas_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(tasmax_file, len(tasmax_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(tasmin_file, len(tasmin_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(hurs_file, len(hurs_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(ps_file, len(ps_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(pr_file, len(pr_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(prsn_file, len(prsn_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(rsds_file, len(rsds_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(rlds_file, len(rlds_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(wind_file, len(wind_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(tref_file, len(tref_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(soc_file, len(soc_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(wlnd_file, len(wlnd_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(veg_file, len(veg_file), MPI_CHARACTER, &
                     0, MPI_COMM_WORLD, err)
      call MPI_BCAST(hydro_file, len(hydro_file), MPI_CHARACTER, &
                     0, MPI_COMM_WORLD, err)
      ! archive group
      call MPI_BCAST(archive_tstep, len(archive_tstep), MPI_CHARACTER, &
                     0, MPI_COMM_WORLD, err)
      call MPI_BCAST(archive_dir, len(archive_dir), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      ! debug group
      call MPI_BCAST(DEBUG, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Write all model simulation settings for re-run job
   !
   !------------------------------------------------------------------------------
   subroutine WriteSimulationSettings(filename)
      implicit none
      character(len=*), intent(in) :: filename
      character(cx) :: fullname
      integer :: error
      namelist /general/ run_mode, lake_file, lake_range, opt_file
      namelist /simulation/ Thermal_Module, Bubble_Module, Diagenesis_Module, &
                            Carbon_Module, Start_Year, Start_Month,  &
                            Start_Day, End_Year, End_Month, End_Day, &
                            Spinup_Month, Spinup_Day, nSpinup, Use_Leap
      namelist /resolution/ NWLAYER, NSLAYER, NRLAYER, NSCOL 
      namelist /bayesian/ NMAXSAMPLE, sample_range, mc_file, SA_Start_Year, &
                          SA_Start_Month, SA_Start_Day, SA_End_Year, &
                          SA_End_Month, SA_End_Day
      namelist /radiation/ solar_dir, gas_dir, albedo_dir, co2_file, &
                           o3_file, aod_file
      namelist /rundata/ forcing_tstep, forcing_time, tas_file, tasmax_file, &
                         tasmin_file, hurs_file, ps_file, pr_file, prsn_file, &
                         rsds_file, rlds_file, wind_file, hydro_file, &
                         tref_file, soc_file, veg_file, wlnd_file
      namelist /archive/ archive_tstep, archive_dir
      namelist /dbg/ DEBUG

      call GetFullFileName(filename,fullname)
      open(unit=fid, file=fullname, status="replace", action="write", &
           delim='APOSTROPHE', iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun('Cannot open file ' // trim(fullname))
      end if
      write(unit=fid, NML=general, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("writing general group")
      end if
      write(unit=fid, NML=simulation, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("writing simulation group")
      end if
      write(unit=fid, NML=resolution, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("writing resolution group")
      end if
      write(unit=fid, NML=bayesian, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("writing bayesian group")
      end if
      write(unit=fid, NML=radiation, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("writing radiation group")
      end if
      write(unit=fid, NML=rundata, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("writing rundata group")
      end if
      write(unit=fid, NML=archive, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("writing archive group")
      end if
      write(unit=fid, NML=dbg, iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun("writing debug group")
      end if
      close(unit=fid)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Read the total number of simulation lakes 
   !
   !------------------------------------------------------------------------------
   subroutine ReadLakeTotalNumber(nlake) 
      implicit none
      integer, intent(out) :: nlake
      character(cx) :: fullname
      character(len=32) :: fname = "ReadLakeTotalNumber" 
      integer(kind=MPI_OFFSET_KIND) :: ncohort
      integer :: ncid, dimid

      call GetFullFileName(lake_file, fullname)
      call check( nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                  NF90_NOWRITE, MPI_INFO_NULL, ncid) )
      call check( nf90mpi_inq_dimid(ncid, "lake", dimid) )
      call check( nf90mpi_inquire_dimension(ncid, dimid, len=ncohort) )
      call check( nf90mpi_close(ncid) )
      nlake = INT(ncohort)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Read and write the calibrated parameters
   !
   !------------------------------------------------------------------------------
   subroutine ReadOptimumParameters(OptParams)
      implicit none
      real(r8), intent(out) :: OptParams(NPARAM)
      character(len=32) :: fname = "ReadOptimumParameters"
      character(cx) :: fullname, tempstr
      character(len=15) :: delimstr
      integer :: error, indx, ii

      if (len_trim(opt_file)==0) then
         if (lake_info%itype==temperate) then
            opt_file = './Parameters/optpar_temperate.dat'
         else if (lake_info%itype==boreal) then
            opt_file = './Parameters/optpar_boreal.dat'
         else if (lake_info%itype==tundra) then
            opt_file = './Parameters/optpar_tundra.dat'
         else if (lake_info%itype==yedoma) then
            opt_file = './Parameters/optpar_yedoma.dat'
         else if (lake_info%itype==arctic) then
            opt_file = './Parameters/optpar_arctic.dat'
         else if (lake_info%itype==alpine) then
            opt_file = './Parameters/optpar_alpine.dat'
         else if (lake_info%itype==eutrophic) then
            opt_file = './Parameters/optpar_eutrophic.dat'
         else if (lake_info%itype==subtropic) then
            opt_file = './Parameters/optpar_subtropic.dat'
         else if (lake_info%itype==amazon) then
            opt_file = './Parameters/optpar_amazon.dat'
         else if (lake_info%itype==rift) then
            opt_file = './Parameters/optpar_rift.dat'
         else if (lake_info%itype==saline) then
            opt_file = './Parameters/optpar_saline.dat'
         else
            call Endrun('Invalid lake_info%itype')
         end if
         call GetFullFileName(opt_file, fullname)
      else
         call GetFullFileName(opt_file, fullname)
      end if
      open(unit=fid, file=fullname, status="old", action="read", &
         iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun(fname, 'Cannot open file ' // trim(fullname))
      end if
      read(fid, "(A512)", iostat=error) tempstr
      if (error/=0) then
         close(unit=fid)
         call Endrun(fname, 'reading ' // trim(fullname))
      end if
      do ii = 1, NPARAM, 1
         read(unit=fid, fmt='(I2, A2, F)', iostat=error) &
            indx, delimstr, OptParams(ii)
         if (error/=0) then
            close(unit=fid)
            call Endrun(fname, 'reading ' // trim(fullname))
         end if
      end do
      close(unit=fid)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: This subroutine is used to read lake information data. 
   !
   !------------------------------------------------------------------------------
   subroutine ReadLakeInfoRealData(filename, varname, lakeId, odata)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      integer, intent(in) :: lakeId
      real(r8), intent(out) :: odata
      character(len=32) :: fname = "ReadLakeInfoRealData"
      integer(kind=MPI_OFFSET_KIND) :: nstart(1), ncount(1)
      character(cx) :: fullname
      integer :: ncid, varid
      real(r4) :: val(1)

      call GetFullFileName(filename, fullname)
      call check( fname, nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                  NF90_NOWRITE, MPI_INFO_NULL, ncid) )
      call check( fname, nf90mpi_inq_varid(ncid, trim(varname), varid) )
      nstart = lakeId 
      ncount = 1
      call check( fname, nf90mpi_get_var_all(ncid, varid, val, &
                  nstart, ncount) )
      call check( fname, nf90mpi_close(ncid) )
      odata = val(1)
      if (DEBUG .and. masterproc) then
         print *, "Read " // trim(varname) // " from " // trim(fullname)
      end if
   end subroutine

   subroutine ReadLakeInfoIntData(filename, varname, lakeId, odata)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      integer, intent(in) :: lakeId 
      integer, intent(out) :: odata
      character(len=32) :: fname = "ReadLakeInfoIntData"
      integer(kind=MPI_OFFSET_KIND) :: nstart(1), ncount(1)
      character(cx) :: fullname
      integer :: ncid, varid
      integer(i1) :: val(1)

      call GetFullFileName(filename, fullname)
      call check( fname, nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                  NF90_NOWRITE, MPI_INFO_NULL, ncid) )
      call check( fname, nf90mpi_inq_varid(ncid, trim(varname), varid) )
      nstart = lakeId
      ncount = 1
      call check( fname, nf90mpi_get_var_all(ncid, varid, val, &
                  nstart, ncount) )
      call check( fname, nf90mpi_close(ncid) )
      odata = val(1)
      if (DEBUG .and. masterproc) then
         print *, "Read " // trim(varname) // " from " // trim(fullname)
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Read air forcing data for a specific lake during given time
   !        Grid: 0.75*0.75, time: 1984~2009 (two data points per day)
   !        LON: 0 : 0.75 : 359.25; LAT: 90 : -0.75 : -90
   !
   !------------------------------------------------------------------------------
   subroutine GetTargetCoordinate(lonindx, latindx, nreslon, nreslat, nrange, &
                                  target_lonid, target_latid, b1start, b1count, &
                                  b2start, b2count, b3start, b3count, b4start, &
                                  b4count)
      implicit none
      integer(kind=MPI_OFFSET_KIND), intent(in) :: lonindx, latindx
      integer(kind=MPI_OFFSET_KIND), intent(in) :: nreslon
      integer(kind=MPI_OFFSET_KIND), intent(in) :: nreslat
      integer(kind=MPI_OFFSET_KIND), intent(in) :: nrange
      integer, intent(out) :: target_lonid(2*nrange)
      integer, intent(out) :: target_latid(2*nrange)
      integer(kind=MPI_OFFSET_KIND), intent(out) :: b1start(2)
      integer(kind=MPI_OFFSET_KIND), intent(out) :: b1count(2)
      integer(kind=MPI_OFFSET_KIND), intent(out) :: b2start(2)
      integer(kind=MPI_OFFSET_KIND), intent(out) :: b2count(2)
      integer(kind=MPI_OFFSET_KIND), intent(out) :: b3start(2)
      integer(kind=MPI_OFFSET_KIND), intent(out) :: b3count(2)
      integer(kind=MPI_OFFSET_KIND), intent(out) :: b4start(2)
      integer(kind=MPI_OFFSET_KIND), intent(out) :: b4count(2)
      integer :: indx0, ii

      if (lonindx>=nrange .and. lonindx+nrange<=nreslon .and. &
            latindx>=nrange .and. latindx+nrange<=nreslat) then
         target_lonid = (/(ii, ii = lonindx-nrange+1, lonindx+nrange)/)
         target_latid = (/(ii, ii = latindx-nrange+1, latindx+nrange)/)
         b1start = (/lonindx - nrange + 1, latindx - nrange + 1/)
         b1count = (/nrange, nrange/)
         b2start = (/lonindx + 1, latindx - nrange + 1/)
         b2count = (/nrange, nrange/)
         b3start = (/lonindx - nrange + 1, latindx + 1/)
         b3count = (/nrange, nrange/)
         b4start = (/lonindx + 1, latindx + 1/)
         b4count = (/nrange, nrange/)
      else if (lonindx<nrange .and. latindx>=nrange .and. &
            latindx+nrange<=nreslat) then
         indx0 = nrange - lonindx + 1
         target_lonid(1:indx0-1) = &
            (/(ii, ii = nreslon-nrange+lonindx+1, nreslon)/)
         target_lonid(indx0:2*nrange) = (/(ii, ii = 1, nrange+lonindx)/)
         target_latid = (/(ii, ii = latindx-nrange+1, latindx+nrange)/)
         b1start = (/nreslon - nrange + lonindx + 1, latindx - nrange + 1/)
         b1count = (/nrange - lonindx, nrange/)
         b2start = (/1_i8, latindx - nrange + 1/)
         b2count = (/nrange + lonindx, nrange/)
         b3start = (/nreslon - nrange + lonindx + 1, latindx + 1/)
         b3count = (/nrange - lonindx, nrange/)
         b4start = (/1_i8, latindx + 1/)
         b4count = (/nrange + lonindx, nrange/)
      else if (lonindx+nrange>nreslon .and. latindx>=nrange .and. &
            latindx+nrange<=nreslat) then
         indx0 = nreslon - lonindx + nrange
         target_lonid(1:indx0) = (/(ii, ii = lonindx-nrange+1, nreslon)/)
         target_lonid(indx0+1:2*nrange) = &
            (/(ii, ii = 1, nrange-nreslon+lonindx)/)
         target_latid = (/(ii, ii = latindx-nrange+1, latindx+nrange)/)
         b1start = (/lonindx - nrange + 1, latindx - nrange + 1/)
         b1count = (/nreslon - lonindx + nrange, nrange/)
         b2start = (/1_i8, latindx - nrange + 1/)
         b2count = (/nrange - nreslon + lonindx, nrange/)
         b3start = (/lonindx - nrange + 1, latindx + 1/)
         b3count = (/nreslon - lonindx + nrange, nrange/)
         b4start = (/1_i8, latindx + 1/)
         b4count = (/nrange - nreslon + lonindx, nrange/)
      else if (lonindx>=nrange .and. lonindx+nrange<=nreslon .and. &
            latindx<nrange) then
         indx0 = nrange - latindx + 1
         target_lonid = (/(ii, ii = lonindx-nrange+1, lonindx+nrange)/)
         target_latid(1:indx0-1) = (/(ii, ii = 1, nrange-latindx)/)
         target_latid(indx0:2*nrange) = (/(ii, ii = 1, nrange+latindx)/)
         b1start = (/lonindx - nrange + 1, 1_i8/)
         b1count = (/nrange, nrange - latindx/)
         b2start = (/lonindx + 1, 1_i8/)
         b2count = (/nrange, nrange - latindx/)
         b3start = (/lonindx - nrange + 1, 1_i8/)
         b3count = (/nrange, nrange + latindx/)
         b4start = (/lonindx + 1, 1_i8/)
         b4count = (/nrange, nrange + latindx/)
      else if (lonindx>=nrange .and. lonindx+nrange<=nreslon .and. &
            latindx+nrange>nreslat) then
         indx0 = nreslat - latindx + nrange
         target_lonid = (/(ii, ii = lonindx-nrange+1, lonindx+nrange)/)
         target_latid(1:indx0) = (/(ii, ii = latindx-nrange+1, nreslat)/)
         target_latid(indx0+1:2*nrange) = &
            (/(ii, ii = 2*nreslat-nrange-latindx+1, nreslat)/)
         b1start = (/lonindx - nrange + 1, latindx - nrange + 1/)
         b1count = (/nrange, nreslat - latindx + nrange/)
         b2start = (/lonindx + 1, latindx - nrange + 1/)
         b2count = (/nrange, nreslat - latindx + nrange/)
         b3start = (/lonindx - nrange + 1, 2*nreslat - nrange - latindx + 1/)
         b3count = (/nrange, nrange + latindx - nreslat/)
         b4start = (/lonindx + 1, 2*nreslat - nrange - latindx + 1/)
         b4count = (/nrange, nrange + latindx - nreslat/)
      end if
   end subroutine

   subroutine ReadSiteTSData(filename, varname, info, time, odata)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      type(LakeInfo), intent(in) :: info 
      type(SimTime), intent(in) :: time
      real(r8), intent(out) :: odata(:)
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND) :: nstart(2)
      integer(kind=MPI_OFFSET_KIND) :: ncount(2)
      integer(kind=MPI_OFFSET_KIND) :: ntottime
      integer, parameter :: ndays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      integer :: ncid, varid, date_varid
      integer :: time_dimid, t0, t1, indx
      integer :: h0, h1, hindx0, hindx1
      integer :: JDN0, JDN1, JDNb, date(1)
      integer :: year, month, day, ii
      integer :: noutval
      real(r8), allocatable :: tmpArr(:,:)
      real(r8) :: filled_value

      ! Units: air temperature (K), relative humidity (%),
      ! wind (m/s), air pressure (Pa), solar radiation (W/m2)
      ! precipitation (kg/m2/s),
      ! atmospheric radiation (W/m2)

      ! inquire data info
      call GetFullFileName(filename, fullname)
      call check( nf90mpi_open(MPI_COMM_WORLD, trim(fullname), NF90_NOWRITE, &
                  MPI_INFO_NULL, ncid) )
      call check( nf90mpi_inq_dimid(ncid, "time", time_dimid) )
      call check( nf90mpi_inquire_dimension(ncid, time_dimid, len = ntottime) )
      call check( nf90mpi_inq_varid(ncid, "date", date_varid) )
      call check( nf90mpi_inq_varid(ncid, trim(varname), varid) )
      call check( nf90mpi_get_att(ncid, varid, "_FillValue", filled_value) )

      ! read data start and end date
      call check( nf90mpi_get_var_all(ncid, date_varid, date, &
                  (/1_i8/), (/1_i8/)) )
      call YYMMDD2Date(date(1), year, month, day)
      call Date2JDN(year, month, day, JDNb)
      if (Use_Leap) then
         call Date2JDN(time%year0, time%month0, time%day0, JDN0)
         call Date2JDN(time%year1, time%month1, time%day1, JDN1)
      else
         JDN0 = JDNb + 365*(time%year0-year-1) + sum(ndays(month:12)) + &
            sum(ndays(1:time%month0-1)) - day + time%day0
         JDN1 = JDNb + 365*(time%year1-year-1) + sum(ndays(month:12)) + &  
            sum(ndays(1:time%month1-1)) - day + time%day1
      end if

      ! read data set
      t0 = max(JDN0 - JDNb + 1, 1)
      t1 = min(max(JDN1 - JDNb, 1), INT4(ntottime))
      nstart = (/info%id, 24 / forcing_nhour * (t0-1) + 1/)
      ncount = (/1, 24 / forcing_nhour * (t1-t0+1)/)
      allocate(tmpArr(1, ncount(2)))
      call check( nf90mpi_get_var_all(ncid, varid, tmpArr, nstart, ncount) )
      call check( nf90mpi_close(ncid) )
      noutval = 0
      do ii = 1, JDN1-JDN0, 1
         indx = JDN0 + ii - JDNb
         h0 = 24 / forcing_nhour * (ii-1) + 1
         h1 = 24 / forcing_nhour * ii
         if (indx>=1 .and. indx<=ntottime) then
            if (JDN0>=JDNb) then
               hindx0 = 24 / forcing_nhour * (ii-1) + 1
               hindx1 = 24 / forcing_nhour * ii
            else
               hindx0 = 24 / forcing_nhour * (indx-1) + 1
               hindx1 = 24 / forcing_nhour * indx
            end if
            odata(h0:h1) = tmpArr(1,hindx0:hindx1)
         else
            odata(h0:h1) = filled_value
            noutval = noutval + 1
         end if
      end do
      deallocate(tmpArr)

      if (noutval>0) then
         call Endrun("Out of valid time period " // trim(varname))
      end if

      if (DEBUG .and. masterproc) then
         print *, "Read " // trim(varname) // " from " // trim(fullname)
      end if
   end subroutine

   subroutine ReadSiteHydroTSData(filename, varname, info, time, odata)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      type(LakeInfo), intent(in) :: info
      type(SimTime), intent(in) :: time
      real(r8), intent(out) :: odata(:)
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND) :: nstart(2)
      integer(kind=MPI_OFFSET_KIND) :: ncount(2)
      integer(kind=MPI_OFFSET_KIND) :: ntottime
      integer, parameter :: ndays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      integer :: ncid, varid, date_varid
      integer :: time_dimid, t0, t1, indx
      integer :: JDN0, JDN1, JDNb, date(1)
      integer :: year, month, day, ii
      real(r8), allocatable :: tmpArr(:,:)
      real(r8) :: filled_value

      ! inquire data info
      call GetFullFileName(filename, fullname)
      call check( nf90mpi_open(MPI_COMM_WORLD, trim(fullname), NF90_NOWRITE, &
                  MPI_INFO_NULL, ncid) )
      call check( nf90mpi_inq_dimid(ncid, "time", time_dimid) )
      call check( nf90mpi_inquire_dimension(ncid, time_dimid, len = ntottime) )
      call check( nf90mpi_inq_varid(ncid, "date", date_varid) )
      call check( nf90mpi_inq_varid(ncid, trim(varname), varid) )
      call check( nf90mpi_get_att(ncid, varid, "_FillValue", filled_value) )

      ! read data start and end date
      call check( nf90mpi_get_var_all(ncid, date_varid, date, &
                  (/1_i8/), (/1_i8/)) )
      call YYMMDD2Date(date(1), year, month, day)
      call Date2JDN(year, month, day, JDNb)
      if (Use_Leap) then
         call Date2JDN(time%year0, time%month0, time%day0, JDN0)
         call Date2JDN(time%year1, time%month1, time%day1, JDN1)
      else
         JDN0 = JDNb + 365*(time%year0-year-1) + sum(ndays(month:12)) + &
            sum(ndays(1:time%month0-1)) - day + time%day0
         JDN1 = JDNb + 365*(time%year1-year-1) + sum(ndays(month:12)) + &
            sum(ndays(1:time%month1-1)) - day + time%day1
      end if

      ! read data set
      t0 = max(JDN0 - JDNb + 1, 1)
      t1 = min(max(JDN1 - JDNb, 1), INT4(ntottime))
      nstart = (/info%id, t0/)
      ncount = (/1, t1-t0+1/)
      allocate(tmpArr(1, ncount(2)))
      call check( nf90mpi_get_var_all(ncid, varid, tmpArr, nstart, ncount) )
      call check( nf90mpi_close(ncid) )
      do ii = 1, JDN1-JDN0, 1
         indx = JDN0 + ii - JDNb 
         if (indx>=1 .and. indx<=ntottime) then
            if (JDN0>=JDNb) then
               odata(ii) = tmpArr(1,ii)
            else
               odata(ii) = tmpArr(1,indx)
            end if
         else
            odata(ii) = filled_value
         end if
      end do

      deallocate(tmpArr)
      if (DEBUG .and. masterproc) then
         print *, "Read " // trim(varname) // " from " // trim(fullname)
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: This subroutine is used to read time-variant 2D data that does not
   !          have scale factor or add offset.
   !
   !------------------------------------------------------------------------------
   subroutine Read2DTSData(filename, varname, lonname, latname, &
                           loc, timeIndx, odata)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: lonname
      character(len=*), intent(in) :: latname
      real(r8), intent(in) :: loc(2)
      integer, intent(in) :: timeIndx(2)
      real(r8), intent(out) :: odata(:)
      character(len=32) :: fname = "Read2DTSData"
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND) :: b1start(3), b1count(3)
      integer(kind=MPI_OFFSET_KIND) :: nreslon, nreslat
      integer :: ncid, varid, lon_dimid, lat_dimid
      integer :: lon_varid, lat_varid, lonindx, latindx
      integer :: ii, ntime
      real(r8), allocatable :: tmplons(:), tmplats(:)
      real(r8), allocatable :: tmpArr(:,:,:)
      real(r8) :: plon, plat, filled_val

      ! Units: AOD (no_units), O3 (DU)
      
      plon = loc(1) 
      plat = loc(2) 
      call GetFullFileName(filename, fullname)
      call check( fname, nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                  NF90_NOWRITE, MPI_INFO_NULL, ncid) )
      call check( fname, nf90mpi_inq_dimid(ncid, trim(lonname), &
                  lon_dimid) )
      call check( fname, nf90mpi_inq_dimid(ncid, trim(latname), &
                  lat_dimid) )
      call check( nf90mpi_inquire_dimension(ncid, lon_dimid, len = nreslon) )
      call check( nf90mpi_inquire_dimension(ncid, lat_dimid, len = nreslat) )
      allocate(tmplons(nreslon))
      allocate(tmplats(nreslat))
      call check( fname, nf90mpi_inq_varid(ncid, trim(lonname), lon_varid) )
      call check( fname, nf90mpi_inq_varid(ncid, trim(latname), lat_varid) )
      call check( nf90mpi_get_var_all(ncid, lon_varid, tmplons) )
      call check( nf90mpi_get_var_all(ncid, lat_varid, tmplats) )
      ! search for corresponding indice
      call DichotomySectionSearch(tmplons, plon, lonindx)
      call DichotomySectionSearch(tmplats, plat, latindx)
      call check( fname, nf90mpi_inq_varid(ncid, trim(varname), varid) )
      call check( fname, nf90mpi_get_att(ncid, varid, "_FillValue", &
                  filled_val) )

      ntime = timeIndx(2) + timeIndx(1) - 1
      if (timeIndx(1)>0) then
         allocate(tmpArr(1,1,timeIndx(2)))
         b1start = (/lonindx, latindx, timeIndx(1)/)
         b1count = (/1, 1, timeIndx(2)/)
         call check( fname, nf90mpi_get_var_all(ncid, varid, tmpArr, &
                     b1start, b1count) )
         odata = tmpArr(1,1,:)
      else if (ntime>0) then
         ntime = timeIndx(2) + timeIndx(1) - 1
         allocate(tmpArr(1,1,ntime))
         b1start = (/lonindx, latindx, 1/)
         b1count = (/1, 1, ntime/)
         call check( fname, nf90mpi_get_var_all(ncid, varid, tmpArr, &
                     b1start, b1count) )
         do ii = 1, timeIndx(2), 1
            if (ii+timeIndx(1)-1>0) then
               odata(ii) = tmpArr(1,1,ii+timeIndx(1)-1)
            else
               odata(ii) = tmpArr(1,1,mod(ii-1,ntime)+1)
            end if
         end do
      else
         allocate(tmpArr(1,1,1))
         b1start = (/lonindx, latindx, 1/)
         b1count = (/1, 1, 1/)
         call check( fname, nf90mpi_get_var_all(ncid, varid, tmpArr, &
                     b1start, b1count) )
         odata = tmpArr(1,1,1)
      end if
      call check( nf90mpi_close(ncid) )

      deallocate(tmpArr)
      deallocate(tmplons)
      deallocate(tmplats)
      if (DEBUG .and. masterproc) then
         print *, "Read " // trim(varname) // " from " // trim(fullname)
      end if
   end subroutine
   
   !------------------------------------------------------------------------------
   !
   ! Purpose: This subroutine is used to read time-invariant data, such as 
   !          climatology data and soil characteristics data.
   !
   !------------------------------------------------------------------------------
   subroutine ReadStaticRealData(filename, varname, lonname, latname, loc, &
                                 res, interp, odata, defval)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: lonname
      character(len=*), intent(in) :: latname
      real(r8), intent(in) :: loc(2)
      real(r8), intent(in) :: res(2)
      integer, intent(in) :: interp
      real(r8), intent(out) :: odata
      real(r8), optional, intent(in) :: defval
      character(len=32) :: fname = "ReadStaticRealData"
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND) :: bstart(2), bcount(2)
      integer(kind=MPI_OFFSET_KIND) :: nreslon, nreslat
      integer(kind=MPI_OFFSET_KIND) :: nrange1, nrange2
      integer :: ncid, varid, lon_varid, lat_varid
      integer :: lon_dimid, lat_dimid, ii, jj
      integer :: ncount, lonindx, latindx
      real(r8), allocatable :: tmplons(:), tmplats(:)
      real(r8) :: lon, lat, plon, plat
      real(r8) :: numer, denom
      real(r4), allocatable :: tmpArr(:,:)
      real(r4) :: filled_val

      plon = loc(1) 
      plat = loc(2)
      call GetFullFileName(filename, fullname)
      call check( fname, nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                  NF90_NOWRITE, MPI_INFO_NULL, ncid) )

      ! data set scope
      call check( fname, nf90mpi_inq_dimid(ncid, trim(lonname), lon_dimid) )
      call check( fname, nf90mpi_inq_dimid(ncid, trim(latname), lat_dimid) )
      call check( nf90mpi_inquire_dimension(ncid, lon_dimid, len = nreslon) )
      call check( nf90mpi_inquire_dimension(ncid, lat_dimid, len = nreslat) )
      allocate(tmplons(nreslon))
      allocate(tmplats(nreslat))
      call check( fname, nf90mpi_inq_varid(ncid, trim(lonname), lon_varid) )
      call check( fname, nf90mpi_inq_varid(ncid, trim(latname), lat_varid) )
      call check( nf90mpi_get_var_all(ncid, lon_varid, tmplons) )
      call check( nf90mpi_get_var_all(ncid, lat_varid, tmplats) )
      call DichotomySectionSearch(tmplons, plon, lonindx)
      call DichotomySectionSearch(tmplats, plat, latindx)

      call check( fname, nf90mpi_inq_varid(ncid, trim(varname), varid) )
      call check( fname, nf90mpi_get_att(ncid, varid, "_FillValue", &
                  filled_val) )

      if (interp==Ncinterp) then

         nrange1 = NINT(0.5 * res(1) * nreslon / 360.0)
         nrange2 = NINT(0.5 * res(2) * nreslat / 180.0)
         allocate(tmpArr(2*nrange1,2*nrange2))
         bstart = (/lonindx-nrange1+1, latindx-nrange2+1/)
         bcount = (/2*nrange1, 2*nrange2/)
         call check( fname, nf90mpi_get_var_all(ncid, varid, tmpArr, &
                     bstart, bcount) )

         numer = 0.0
         denom = 0.0
         do ii = 1, 2*nrange1, 1
            do jj = 1, 2*nrange2, 1
               if (tmpArr(ii,jj)/=filled_val) then
                  numer = numer + tmpArr(ii,jj)
                  denom = denom + 1.0 
               end if
            end do
         end do

         call check( nf90mpi_close(ncid) )

         if (denom>0) then
            odata = numer / denom
         else if (present(defval)) then
            odata = defval
         else
            call Endrun(trim(varname) // ": cannot find valid scope") 
         end if
      else
         allocate(tmpArr(1,1))
         bstart = (/lonindx, latindx/)
         bcount = (/1, 1/)
         call check( fname, nf90mpi_get_var_all(ncid, varid, tmpArr, &
                     bstart, bcount) )
         odata = tmpArr(1,1)
         call check( nf90mpi_close(ncid) )
      end if

      if (allocated(tmpArr))        deallocate(tmpArr)
      if (allocated(tmplons))       deallocate(tmplons)
      if (allocated(tmplats))       deallocate(tmplats)
      if (DEBUG .and. masterproc) then
         print *, "Read " // trim(varname) // " from " // trim(fullname)
      end if
   end subroutine

   subroutine ReadStaticByteData(filename, varname, lonname, latname, &
                                 loc, odata)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: lonname
      character(len=*), intent(in) :: latname
      real(r8), intent(in) :: loc(2)
      real(r8), intent(out) :: odata
      character(len=32) :: fname = "ReadStaticByteData"
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND) :: nreslon, nreslat
      integer(kind=MPI_OFFSET_KIND) :: nstart(2), ncount(2)
      integer :: ncid, varid, lon_varid, lat_varid
      integer :: lon_dimid, lat_dimid
      integer :: lonindx, latindx
      integer(i1) :: tmpArr(1,1) 
      real(r8), allocatable :: tmplons(:), tmplats(:)
      real(r8) :: plon, plat

      plon = loc(1)
      plat = loc(2)
      call GetFullFileName(filename, fullname)
      call check( fname, nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                  NF90_NOWRITE, MPI_INFO_NULL, ncid) )

      ! data set scope
      call check( fname, nf90mpi_inq_dimid(ncid, trim(lonname), lon_dimid) )
      call check( fname, nf90mpi_inq_dimid(ncid, trim(latname), lat_dimid) )
      call check( nf90mpi_inquire_dimension(ncid, lon_dimid, len=nreslon) )
      call check( nf90mpi_inquire_dimension(ncid, lat_dimid, len=nreslat) )
      allocate(tmplons(nreslon))
      allocate(tmplats(nreslat))
      call check( fname, nf90mpi_inq_varid(ncid, trim(lonname), lon_varid) )
      call check( fname, nf90mpi_inq_varid(ncid, trim(latname), lat_varid) )
      call check( fname, nf90mpi_inq_varid(ncid, trim(varname), varid) )
      call check( nf90mpi_get_var_all(ncid, lon_varid, tmplons) )
      call check( nf90mpi_get_var_all(ncid, lat_varid, tmplats) )
      call DichotomySectionSearch(tmplons, plon, lonindx)
      call DichotomySectionSearch(tmplats, plat, latindx)
      deallocate(tmplons)
      deallocate(tmplats)

      nstart = (/lonindx, latindx/)
      ncount = (/1, 1/)
      call check( fname, nf90mpi_get_var_all(ncid, varid, tmpArr, &
                  nstart, ncount) )
      call check( nf90mpi_close(ncid) )
      odata = DBLE(tmpArr(1,1))

      if (DEBUG .and. masterproc) then
         print *, "Read " // trim(varname) // " from " // trim(fullname)
      end if
   end subroutine

   subroutine RescaleTreeCoverFraction(ftree)
      implicit none
      real(r8), intent(inout) :: ftree

      if (ftree==127) then
         ftree = 5.0d-2
      else if (ftree>=10 .and. ftree<=80) then
         ftree = 1.0d-2 * ftree
      else
         ftree = 0.0_r8
      end if
   end subroutine

   subroutine RescaleWetlandFraction(fwlnd)
      implicit none
      real(r8), intent(inout) :: fwlnd

      if (fwlnd==4 .or. fwlnd==5 .or. fwlnd==6 .or. fwlnd==8) then
         fwlnd = 1.0
      else if (fwlnd==10) then
         fwlnd = 0.75
      else if (fwlnd==11) then
         fwlnd = 0.375
      else if (fwlnd==12) then
         fwlnd = 0.125
      else
         fwlnd = 0.0_r8
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: This subroutine is used to read global mean time series. 
   !
   !------------------------------------------------------------------------------
   subroutine ReadGlobalTSData(filename, varname, timeIndx, odata)
      implicit none
      character(len=*), intent(in) :: filename
      character(len=*), intent(in) :: varname
      integer, intent(in) :: timeIndx(2)
      real(r8), intent(out) :: odata(:)
      character(len=32) :: fname = "ReadGlobalTSData"
      integer(kind=MPI_OFFSET_KIND) :: nstart(1), ncount(1)
      character(cx) :: fullname
      integer :: ncid, varid

      call GetFullFileName(filename, fullname)
      call check( fname, nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                  NF90_NOWRITE, MPI_INFO_NULL, ncid) )
      call check( fname, nf90mpi_inq_varid(ncid, trim(varname), varid) )
      nstart = timeIndx(1)
      ncount = timeIndx(2)
      call check( fname, nf90mpi_get_var_all(ncid, varid, odata, &
                  nstart, ncount) )
      call check( fname, nf90mpi_close(ncid) )
      if (DEBUG .and. masterproc) then
         print *, "Read " // trim(varname) // " from " // trim(fullname)
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Create a netcdf file for model output storage. It includes the
   !        creation of file, dimensions and variables. The dimension lengths should
   !        the maximum of all lakes.    
   !
   !------------------------------------------------------------------------------
   subroutine CreateOutputFile0d(time, nz, varname, longname, units)
      implicit none
      Type(SimTime), intent(in)  :: time
      integer, intent(in) :: nz
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units
      character(cx) :: fullname, histstr
      character(cx) :: str1, str2
      integer :: ncid, varid, cmode
      integer :: z_dimid, lake_dimid
      integer :: gmtime(6), err

      call GetArchiveFullname(time, varname, archive_tstep, fullname) 
      if (masterproc) then
         call GetGMTime(Use_Leap, gmtime)
         write(str1,"(I4, '-', I2.2, '-', I2.2)") gmtime(1), gmtime(2), &
               gmtime(3)
         write(str2,"(' ', I2.2, ':', I2.2, ':', I2.2, A)") gmtime(4), &
               gmtime(5), gmtime(6), " GMT from ALBM v2.0 by Zeli Tan"
         histstr = trim(str1) // trim(str2)
      end if
      !call MPI_BCAST(histstr, len(histstr), MPI_CHARACTER, 0, &
      !               MPI_COMM_WORLD, err)

      ! Create the file.
      cmode = IOR(NF90_CLOBBER, NF90_64BIT_OFFSET)
      call check( nf90mpi_create(MPI_COMM_SELF, trim(fullname), &
                  cmode, MPI_INFO_NULL, ncid) )

      ! Define the dimensions
      call check( nf90mpi_def_dim(ncid, 'Z', INT8(nz), z_dimid) )
      call check( nf90mpi_def_dim(ncid, 'Lake', NFMPI_UNLIMITED, lake_dimid) )

      ! Define the coordinate variables
      call check( nf90mpi_put_att(ncid, NF90_GLOBAL, "history", &
                  trim(histstr)) )

      ! ! Define the netCDF variables and assign attributes for state variables. 
      call DefNcVariable(ncid, (/z_dimid, lake_dimid/), trim(varname), &
                        trim(longname), trim(units), -9999.0_r8, &
                        -9999.0_r8, varid)

      ! End define mode.
      call check( nf90mpi_enddef(ncid) )

      ! Write the coordinate variable data. This will put our data grid 
      ! into the netCDF file.

      ! Close the file. This causes netCDF to flush all buffers and make
      ! sure your data are really written to disk.
      call check( nf90mpi_close(ncid) )

      if (masterproc) then
         print *, 'Create ouput file ' // trim(fullname)
      end if
   end subroutine

   subroutine CreateOutputFile1d(time, varname, longname, units, defval)
      implicit none
      Type(SimTime), intent(in)  :: time
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units
      real(r4), intent(in) :: defval
      character(cx) :: fullname, histstr, timestr
      character(cx) :: str1, str2
      integer(kind=MPI_OFFSET_KIND) :: ntime
      integer :: ncid, varid, cmode
      integer :: time_dimid, lake_dimid
      integer :: gmtime(6), simday, err

      call GetArchiveFullname(time, varname, archive_tstep, fullname)
      if (masterproc) then
         call GetGMTime(Use_Leap, gmtime)
         write(str1,"(I4, '-', I2.2, '-', I2.2)") gmtime(1), gmtime(2), &
               gmtime(3)
         write(str2,"(' ', I2.2, ':', I2.2, ':', I2.2, A)") gmtime(4), &
               gmtime(5), gmtime(6), " GMT from ALBM v2.0 by Zeli Tan"         
         histstr = trim(str1) // trim(str2)
      end if
      !call MPI_BCAST(histstr, len(histstr), MPI_CHARACTER, 0, &
      !               MPI_COMM_WORLD, err)

      write(str1,"('Records start since ', I4, '-', I2.2, '-', I2.2)") &
         time%year0, time%month0, time%day0
      timestr = trim(str1) // " 00:00:00 GMT"
      simday = CalcRunningDays(time, Use_Leap)
      if (trim(archive_tstep)=='day') then
         ntime = simday
      else if (trim(archive_tstep)=='hour') then
         ntime = 24 * simday
      end if

      ! Create the file.
      cmode = IOR(NF90_CLOBBER, NF90_64BIT_OFFSET)
      call check( nf90mpi_create(MPI_COMM_SELF, trim(fullname), &
                  cmode, MPI_INFO_NULL, ncid) )

      ! Define the dimensions
      call check( nf90mpi_def_dim(ncid, 'Time', ntime, time_dimid) )
      call check( nf90mpi_def_dim(ncid, 'Lake', NFMPI_UNLIMITED, lake_dimid) )

      ! Define the coordinate variables
      call check( nf90mpi_put_att(ncid, NF90_GLOBAL, "history", &
                  trim(histstr)) )
      call check( nf90mpi_put_att(ncid, NF90_GLOBAL, "description", &
                  trim(timestr)) )

      ! ! Define the netCDF variables and assign attributes for state variables. 
      call DefNcVariable(ncid, (/time_dimid, lake_dimid/), trim(varname), &
                        trim(longname), trim(units), defval, defval, varid)

      ! End define mode.
      call check( nf90mpi_enddef(ncid) )

      ! Write the coordinate variable data. This will put our data grid 
      ! into the netCDF file.

      ! Close the file. This causes netCDF to flush all buffers and make
      ! sure your data are really written to disk.
      call check( nf90mpi_close(ncid) )

      if (masterproc) then
         print *, 'Create ouput file ' // trim(fullname)
      end if
   end subroutine

   subroutine CreateOutputFile2d(time, ndim, varname, dimname, longname, &
                  units, defval)
      implicit none
      Type(SimTime), intent(in)  :: time
      integer, intent(in) :: ndim
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: dimname
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units
      real(r4), intent(in) :: defval
      character(cx) :: fullname, histstr, timestr
      character(cx) :: str1, str2
      integer(kind=MPI_OFFSET_KIND) :: ntime
      integer :: ncid, varid, cmode
      integer :: dimid, time_dimid, lake_dimid
      integer :: gmtime(6), simday, err

      call GetArchiveFullname(time, varname, archive_tstep, fullname)
      if (masterproc) then
         call GetGMTime(Use_Leap, gmtime)
         write(str1,"(I4, '-', I2.2, '-', I2.2)") gmtime(1), gmtime(2), &
               gmtime(3)
         write(str2,"(' ', I2.2, ':', I2.2, ':', I2.2, A)") gmtime(4), &
               gmtime(5), gmtime(6), " GMT from ALBM v2.0 by Zeli Tan"
         histstr = trim(str1) // trim(str2)
      end if
      !call MPI_BCAST(histstr, len(histstr), MPI_CHARACTER, 0, &
      !               MPI_COMM_WORLD, err)

      write(str1,"('Records start since ', I4, '-', I2.2, '-', I2.2)") &
         time%year0, time%month0, time%day0
      timestr = trim(str1) // " 00:00:00 GMT"
      simday = CalcRunningDays(time, Use_Leap)
      if (trim(archive_tstep)=='day') then
         ntime = simday
      else if (trim(archive_tstep)=='hour') then
         ntime = 24 * simday
      end if

      ! Create the file.
      cmode = IOR(NF90_CLOBBER, NF90_64BIT_OFFSET)
      call check( nf90mpi_create(MPI_COMM_SELF, trim(fullname), &
                  cmode, MPI_INFO_NULL, ncid) )

      ! Define the dimensions
      call check( nf90mpi_def_dim(ncid, trim(dimname), INT8(ndim), dimid) )
      call check( nf90mpi_def_dim(ncid, 'Time', ntime, time_dimid) )
      call check( nf90mpi_def_dim(ncid, 'Lake', NFMPI_UNLIMITED, lake_dimid) )

      ! Define the coordinate variables
      call check( nf90mpi_put_att(ncid, NF90_GLOBAL, "history", &
                  trim(histstr)) )
      call check( nf90mpi_put_att(ncid, NF90_GLOBAL, "description", &
                  trim(timestr)) )

      ! ! Define the netCDF variables and assign attributes for state variables. 
      call DefNcVariable(ncid, (/dimid, time_dimid, lake_dimid/), &
                        trim(varname), trim(longname), trim(units), &
                        defval, defval, varid)

      ! End define mode.
      call check( nf90mpi_enddef(ncid) )

      ! Write the coordinate variable data. This will put our data grid 
      ! into the netCDF file.

      ! Close the file. This causes netCDF to flush all buffers and make
      ! sure your data are really written to disk.
      call check( nf90mpi_close(ncid) )

      if (masterproc) then
         print *, 'Create ouput file ' // trim(fullname)
      end if
   end subroutine

   subroutine CreateOutputFile3d(time, ndim1, ndim2, varname, &
                  dimname1, dimname2, longname, units, defval)
      implicit none
      Type(SimTime), intent(in)  :: time
      integer, intent(in) :: ndim1
      integer, intent(in) :: ndim2
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: dimname1
      character(len=*), intent(in) :: dimname2
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units
      real(r4), intent(in) :: defval
      character(cx) :: fullname, histstr, timestr
      character(cx) :: str1, str2
      integer(kind=MPI_OFFSET_KIND) :: ntime
      integer :: ncid, varid, cmode
      integer :: dimid1, dimid2, time_dimid, lake_dimid
      integer :: gmtime(6), simday, err

      call GetArchiveFullname(time, varname, archive_tstep, fullname)
      if (masterproc) then
         call GetGMTime(Use_Leap, gmtime)
         write(str1,"(I4, '-', I2.2, '-', I2.2)") gmtime(1), gmtime(2), &
               gmtime(3)
         write(str2,"(' ', I2.2, ':', I2.2, ':', I2.2, A)") gmtime(4), &
               gmtime(5), gmtime(6), " GMT from ALBM v2.0 by Zeli Tan"
         histstr = trim(str1) // trim(str2)
      end if
      !call MPI_BCAST(histstr, len(histstr), MPI_CHARACTER, 0, &
      !               MPI_COMM_WORLD, err)

      write(str1,"('Records start since ', I4, '-', I2.2, '-', I2.2)") &
         time%year0, time%month0, time%day0
      timestr = trim(str1) // " 00:00:00 GMT"
      simday = CalcRunningDays(time, Use_Leap)
      if (trim(archive_tstep)=='day') then
         ntime = simday
      else if (trim(archive_tstep)=='hour') then
         ntime = 24 * simday
      end if

      ! Create the file.
      cmode = IOR(NF90_CLOBBER, NF90_64BIT_OFFSET)
      call check( nf90mpi_create(MPI_COMM_SELF, trim(fullname), &
                  cmode, MPI_INFO_NULL, ncid) )

      ! Define the dimensions
      call check( nf90mpi_def_dim(ncid, trim(dimname1), INT8(ndim1), dimid1) )
      call check( nf90mpi_def_dim(ncid, trim(dimname2), INT8(ndim2), dimid2) )
      call check( nf90mpi_def_dim(ncid, 'Time', ntime, time_dimid) )
      call check( nf90mpi_def_dim(ncid, 'Lake', NFMPI_UNLIMITED, lake_dimid) )

      ! Define the coordinate variables
      call check( nf90mpi_put_att(ncid, NF90_GLOBAL, "history", &
                  trim(histstr)) )
      call check( nf90mpi_put_att(ncid, NF90_GLOBAL, "description", &
                  trim(timestr)) )

      ! ! Define the netCDF variables and assign attributes for state variables. 
      call DefNcVariable(ncid, (/dimid1, dimid2, time_dimid, lake_dimid/), &
                        trim(varname), trim(longname), trim(units), &
                        defval, defval, varid)

      ! End define mode.
      call check( nf90mpi_enddef(ncid) )

      ! Write the coordinate variable data. This will put our data grid 
      ! into the netCDF file.

      ! Close the file. This causes netCDF to flush all buffers and make
      ! sure your data are really written to disk.
      call check( nf90mpi_close(ncid) )

      if (masterproc) then
         print *, 'Create ouput file ' // trim(fullname)
      end if
   end subroutine

end module read_data_mod
