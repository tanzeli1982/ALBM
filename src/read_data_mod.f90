module read_data_mod
!---------------------------------------------------------------------------------
! Purpose:
!
! This module serves for the data preparation and storage
! ECMWF data reanalysis: http://stommel.tamu.edu/~baum/ecmwf.html
!
!---------------------------------------------------------------------------------
   use shr_kind_mod,       only : cs => SHR_KIND_CS, cx => SHR_KIND_CX, &
                                  r8, r4, i8, i1
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

   interface CalcCostfunc4Var
      module procedure CalcCostfunc4Var0D
      module procedure CalcCostfunc4Var1D
   end interface

   interface ReadLakeInfoData
      module procedure ReadLakeInfoRealData
      module procedure ReadLakeInfoIntData
   end interface

   interface CreateOutputFile
      module procedure CreateOutputFile0d
      module procedure CreateOutputFile1d
      module procedure CreateOutputFile2d
   end interface

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: Read lake information (lat, lon, depth, area) and calculate discrete
   !          water layers
   !
   !------------------------------------------------------------------------------
   subroutine ReadLakeInfo(lakeId)
      implicit none
      integer, intent(in) :: lakeId
      real(r8) :: lon, lat, depth, area
      real(r8) :: basin, zalt, excice, kext
      integer :: thrmkst, hydroconn 

      ! read lake file
      call ReadLakeInfoData(lakeId, 'lon', lon)     ! deg
      call ReadLakeInfoData(lakeId, 'lat', lat)     ! deg
      call ReadLakeInfoData(lakeId, 'depth', depth) ! m
      call ReadLakeInfoData(lakeId, 'basin', basin) ! m
      call ReadLakeInfoData(lakeId, 'Asurf', area)  ! km^2
      call ReadLakeInfoData(lakeId, 'zalt', zalt)   ! m
      call ReadLakeInfoData(lakeId, 'kext', kext)   ! m^-1
      call ReadLakeInfoData(lakeId, 'ice', excice)  ! fraction
      call ReadLakeInfoData(lakeId, 'thrmkst', thrmkst)
      call ReadLakeInfoData(lakeId, 'hydroconn', hydroconn)
      lake_info%id = lakeId
      lake_info%latitude = lat
      lake_info%longitude = lon
      lake_info%depth = depth
      lake_info%basin = basin
      lake_info%Asurf = 1d6 * area  ! m^2
      lake_info%zalt = 1d-3 * zalt  ! km
      lake_info%excice = excice 
      lake_info%hsed = SED_DEPTH
      lake_info%margin = 0
      lake_info%thrmkst = thrmkst
      lake_info%hydroconn = hydroconn
      lake_info%itype = temperate_lake
      if (kext>0) then
         lake_info%kext = sa_params(Param_Feta) * kext
      else
         lake_info%kext = sa_params(Param_Feta) * CalcLightExtinction(depth)
      end if
      if (lake_info%depth<=0.1*NWLAYER) then
         WATER_LAYER = INT(10 * lake_info%depth)
      else
         WATER_LAYER = NWLAYER
      end if
      if (DEBUG) then
         print *, "Lake information: ", lake_info
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Read lake bathymetry.
   !
   !------------------------------------------------------------------------------
   subroutine ReadLakeBathymetry(info, Zw, Az, dAz)
      implicit none
      type(LakeInfo), intent(inout) :: info
      real(r8), intent(in) :: Zw(:)
      real(r8), intent(out) :: Az(:)
      real(r8), intent(out) :: dAz(:)
      character(cx) :: filename, fullname
      character(cx) :: tmpstr, msg
      real(r8), allocatable :: tmpZw(:), tmpAz(:)
      integer :: nline, error, ii, nz

      if (len_trim(lakeid_file)==0) then
         write(filename, "(A,I0,A)") trim(bthmtry_dir) // 'bathymetry_', &
            lake_info%id, '.dat'
      else
         filename = trim(bthmtry_dir) // 'bathymetry_' // &
            trim(lake_info%name) // '.dat'
      end if
      call GetFullFileName(filename, fullname)
      open(unit=fid, file=trim(fullname), status="old", action="read", &
           iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun('Cannot open file ' // trim(fullname))
      end if
      read(unit=fid, fmt="(A9,I)", iostat=error) tmpstr, nline
      allocate(tmpZw(nline))
      allocate(tmpAz(nline))
      read(fid, "(A512)", iostat=error) tmpstr
      read(fid, "(A512)", iostat=error) tmpstr
      do ii = 1, nline, 1
         read(unit=fid, fmt=*, iostat=error) tmpZw(ii), tmpAz(ii)
         if (error/=0) then
            close(unit=fid)
            write(msg, "(A, I0)") "Reading stops at line ", ii
            call Endrun(trim(msg))
         end if
      end do
      close(unit=fid)
      info%Abasin = tmpAz(nline)
      do ii = 1, nline, 1
         if (tmpZw(ii)>info%basin) then
            info%Abasin = tmpAz(ii)
            exit
         end if
      end do
      ! interpolate
      call Interp1d(tmpZw, tmpAz, Zw, Az)
      nz = size(Az)
      do ii = 1, nz, 1
         if (ii==1) then
            dAz(ii) = 0.5 * (Az(ii) - Az(ii+1))
         else if (ii==nz) then
            dAz(ii) = 0.5 * (Az(ii) + Az(ii-1))
         else
            dAz(ii) = 0.5 * (Az(ii-1) - Az(ii+1))
         end if
      end do
      deallocate(tmpZw)
      deallocate(tmpAz)
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
      open(unit=fid, file=trim(fullname), status="old", action="read", &
           iostat=error)
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
   ! Purpose: Read all model simulation settings
   !
   !------------------------------------------------------------------------------
   subroutine ReadSimulationSettings(filename)
      implicit none
      character(len=*), intent(in) :: filename
      character(cx) :: fullname
      integer :: error
      namelist /general/ run_mode, lake_file, lakeid_file, lake_range, &
                         bthmtry_dir, param_dir, opt_file
      namelist /simulation/ Thermal_Module, Bubble_Module, Diagenesis_Module, &
                            Carbon_Module, Hydro_Module, Start_Year, &
                            Start_Month,  Start_Day, End_Year, End_Month, &
                            End_Day, Spinup_Month, Spinup_Day, nSpinup
      namelist /resolution/ NWLAYER, NSLAYER, NRLAYER 
      namelist /bayesian/ NMAXSAMPLE, sample_range, obs_dir, obs_var, &
                          obs_weight, mc_file, sa_file
      namelist /radiation/ solar_dir, gas_dir, albedo_dir, co2_file, &
                           o3_file, aod_file
      namelist /rundata/ forcing_tstep, forcing_dir, hydro_dir, tref_file, & 
                         soc_file, veg_file, wlnd_file, tas_file, &
                         tasmax_file, tasmin_file, hurs_file, ps_file, &
                         pr_file, prsn_file, rsds_file, rlds_file, &
                         wind_file
      namelist /archive/ archive_tstep, archive_dir
      namelist /dbg/ DEBUG, RESUBMIT
      
      call GetFullFileName(filename, fullname)
      open(unit=fid, file=trim(fullname), status="old", action="read", &
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

      ! check the consistency of namelist
      if (len_trim(opt_file)==0 .and. len_trim(param_dir)==0) then
         call Endrun("Must set either opt_file or param_dir")
      end if
      if (len_trim(forcing_dir)==0 .and. len_trim(tas_file)==0) then
         call Endrun("Must set either forcing_dir or tas_file")
      end if
      if (len_trim(forcing_dir)==0 .and. len_trim(tasmax_file)==0) then
         call Endrun("Must set either forcing_dir or tasmax_file")
      end if
      if (len_trim(forcing_dir)==0 .and. len_trim(tasmin_file)==0) then
         call Endrun("Must set either forcing_dir or tasmin_file")
      end if
      if (len_trim(forcing_dir)==0 .and. len_trim(hurs_file)==0) then
         call Endrun("Must set either forcing_dir or hurs_file")
      end if
      if (len_trim(forcing_dir)==0 .and. len_trim(ps_file)==0) then
         call Endrun("Must set either forcing_dir or ps_file")
      end if
      if (len_trim(forcing_dir)==0 .and. len_trim(pr_file)==0) then
         call Endrun("Must set either forcing_dir or pr_file")
      end if
      if (len_trim(forcing_dir)==0 .and. len_trim(prsn_file)==0) then
         call Endrun("Must set either forcing_dir or prsn_file")
      end if
      if (len_trim(forcing_dir)==0 .and. len_trim(rsds_file)==0) then
         call Endrun("Must set either forcing_dir or rsds_file")
      end if
      if (len_trim(forcing_dir)==0 .and. len_trim(rlds_file)==0) then
         call Endrun("Must set either forcing_dir or rlds_file")
      end if
      if (len_trim(forcing_dir)==0 .and. len_trim(wind_file)==0) then
         call Endrun("Must set either forcing_dir or wind_file")
      end if
   end subroutine

   subroutine BcastSimulationSettings()
      implicit none
      integer :: err

      ! general group
      call MPI_BCAST(run_mode, len(run_mode), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(lake_file, len(lake_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(lakeid_file, len(lakeid_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(lake_range, size(lake_range), MPI_INTEGER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(bthmtry_dir, len(bthmtry_dir), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(param_dir, len(param_dir), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(opt_file, len(opt_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      ! simulation group
      call MPI_BCAST(Thermal_Module, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Bubble_Module, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Diagenesis_Module, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Carbon_Module, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Hydro_Module, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Start_Year, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Start_Month, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Start_Day, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(End_Year, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(End_Month, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(End_Day, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Spinup_Month, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(Spinup_Day, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(nSpinup, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      ! resolution group
      call MPI_BCAST(NWLAYER, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(NSLAYER, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(NRLAYER, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      ! bayesian group
      call MPI_BCAST(NMAXSAMPLE, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(sample_range, size(sample_range), MPI_INTEGER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(obs_dir, len(obs_dir), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(obs_var, len(obs_var), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(obs_weight, len(obs_weight), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(mc_file, len(mc_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(sa_file, len(sa_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
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
      call MPI_BCAST(forcing_dir, len(forcing_dir), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(hydro_dir, len(hydro_dir), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(tref_file, len(tref_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(soc_file, len(soc_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(wlnd_file, len(wlnd_file), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(veg_file,len(veg_file), MPI_CHARACTER, 0, &
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
      ! archive group
      call MPI_BCAST(archive_tstep, len(archive_tstep), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(archive_dir, len(archive_dir), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)
      ! debug group
      call MPI_BCAST(DEBUG, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
      call MPI_BCAST(RESUBMIT, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, err)
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
      namelist /general/ run_mode, lake_file, lakeid_file, lake_range, &
                         bthmtry_dir, param_dir, opt_file
      namelist /simulation/ Thermal_Module, Bubble_Module, Diagenesis_Module, &
                            Carbon_Module, Hydro_Module, Start_Year, &
                            Start_Month, Start_Day, End_Year, End_Month, &
                            End_Day, Spinup_Month, Spinup_Day, nSpinup 
      namelist /resolution/ NWLAYER, NSLAYER, NRLAYER 
      namelist /bayesian/ NMAXSAMPLE, sample_range, obs_dir, obs_var, &
                          obs_weight, mc_file, sa_file
      namelist /radiation/ solar_dir, gas_dir, albedo_dir, co2_file, &
                           o3_file, aod_file
      namelist /rundata/ forcing_tstep, forcing_dir, hydro_dir,  tref_file, &
                         soc_file, veg_file, wlnd_file, tas_file, &
                         tasmax_file, tasmin_file, hurs_file, ps_file, &
                         pr_file, prsn_file, rsds_file, rlds_file, &
                         wind_file
      namelist /archive/ archive_tstep, archive_dir
      namelist /dbg/ DEBUG, RESUBMIT

      call GetFullFileName(filename,fullname)
      open(unit=fid, file=trim(fullname), status="replace", action="write", &
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
   ! Purpose: Save results for Monte Carlo sensitivity analysis
   !
   !------------------------------------------------------------------------------
   subroutine CreateSampleResultFile(nrec, fillval)
      implicit none
      integer, intent(in) :: nrec
      real(r8), intent(in) :: fillval
      character(len=32) :: fname = "CreateSampleResultFile"
      character(cx) :: fullname, histstr
      integer :: gmtime(6), ncid, varid
      integer :: sample_dimid, rec_dimid
      integer :: err

      call GetFullFileName(sa_file, fullname)
      ! create the sample output file
      if (masterproc) then
         call GetGMTime(gmtime)
         write(histstr,"(I4, A, I2.2, A, I2.2, ' ', I2.2, ':', " // &
            "I2.2, ':', I2.2, A)") gmtime(1), '-', gmtime(2), '-', &
            gmtime(3), gmtime(4), gmtime(5), gmtime(6), &
            " GMT from ALBM by zeli tan"
      end if
      call MPI_BCAST(histstr, len(histstr), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)

      call check( nf90mpi_create(MPI_COMM_WORLD, trim(fullname), &
                  NF90_CLOBBER, MPI_INFO_NULL, ncid) )
      call check( nf90mpi_def_dim(ncid, "rec", INT8(nrec), rec_dimid) )
      call check( nf90mpi_def_dim(ncid, "id", NFMPI_UNLIMITED, &
                  sample_dimid) )
      call check( nf90mpi_put_att(ncid, NF90_GLOBAL, "history", &
                  trim(histstr)) )
      call DefNcVariable(ncid, (/rec_dimid, sample_dimid/), "sa", &
                  "sample return value", "n/a", fillval, fillval, varid)
      call check( nf90mpi_enddef(ncid) )
      call check( nf90mpi_close(ncid) )

      if (masterproc) then
         print *, 'Create sample result file ' // trim(fullname)
      end if
   end subroutine

   subroutine WriteSampleResults(sp_range, results)
      implicit none
      integer, intent(in) :: sp_range(2)
      real(r8), intent(inout) :: results(:,:)
      character(len=32) :: fname = "WriteSampleResults"
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND) :: nstart(2), ncount(2)
      integer :: ncid, varid, minid, maxid
      integer :: nrec

      maxid = maxval(sp_range)
      minid = minval(sp_range)
      nrec = size(results,1)
      call GetFullFileName(sa_file, fullname)
      nstart = (/1, minid/)
      ncount = (/nrec, maxid-minid+1/)
      call check( fname, nf90mpi_open(MPI_COMM_SELF, trim(fullname), &
                  NF90_WRITE, MPI_INFO_NULL, ncid) )
      call check( fname, nf90mpi_inq_varid(ncid, "sa", varid) )
      call check( fname, nf90mpi_put_var_all(ncid, varid, results, &
                  nstart, ncount) )
      call check( fname, nf90mpi_close(ncid) )
      if (masterproc) then
         print *, 'Save sample results in ' // trim(fullname)
      end if
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
      character(cx) :: filename, fullname
      character(cx) :: tempstr
      character(len=15) :: delimstr
      integer :: error, indx, ii

      if (len_trim(opt_file)==0) then
         if (len_trim(lakeid_file)==0) then
            write(filename, "(A,I0,A)") trim(param_dir) // 'optpar_', &
               lake_info%id, '.dat'
         else
            filename = trim(param_dir) // 'optpar_' // &
               trim(lake_info%name) // '.dat'
         end if
         call GetFullFileName(filename, fullname)
      else
         call GetFullFileName(opt_file, fullname)
      end if
      open(unit=fid, file=trim(fullname), status="old", action="read", &
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
         read(unit=fid, fmt='(I2, A2, E14.6)', iostat=error) &
            indx, delimstr, OptParams(ii)
         if (error/=0) then
            close(unit=fid)
            call Endrun(fname, 'reading ' // trim(fullname))
         end if
      end do
      close(unit=fid)
   end subroutine

   subroutine ReadCalibVariables(vars, weights, nvar)
      implicit none
      character(len=8), intent(inout) :: vars(:)
      real(r8), intent(inout) :: weights(:)
      integer, intent(out) :: nvar
      integer :: pos2, pos1

      nvar = 0
      pos1 = 1
      do while (.True.)
         pos2 = index(obs_var(pos1:), ",")
         if (pos2==0) then
            nvar = nvar + 1
            vars(nvar) = obs_var(pos1:)
            exit
         end if
         nvar = nvar + 1
         vars(nvar) = obs_var(pos1:pos1+pos2-2)
         pos1 = pos2 + pos1
      end do
      read(obs_weight,fmt=*) weights(1:nvar)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: This subroutine is used to read lake information data. 
   !
   !------------------------------------------------------------------------------
   subroutine ReadLakeInfoRealData(lakeId, varname, odata)
      implicit none
      integer, intent(in) :: lakeId
      character(len=*), intent(in) :: varname
      real(r8), intent(out) :: odata
      character(len=32) :: fname = "ReadLakeInfoRealData"
      integer(kind=MPI_OFFSET_KIND) :: nstart(1), ncount(1)
      character(cx) :: fullname
      integer :: ncid, varid
      real(r4) :: val(1)

      call GetFullFileName(lake_file, fullname)
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

   subroutine ReadLakeInfoIntData(lakeId, varname, odata)
      implicit none
      integer, intent(in) :: lakeId
      character(len=*), intent(in) :: varname
      integer, intent(out) :: odata
      character(len=32) :: fname = "ReadLakeInfoIntData"
      integer(kind=MPI_OFFSET_KIND) :: nstart(1), ncount(1)
      character(cx) :: fullname
      integer :: ncid, varid
      integer(i1) :: val(1)

      call GetFullFileName(lake_file, fullname)
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
   ! Purpose: Read the lake_id to lake_name map.
   !
   !------------------------------------------------------------------------------
   subroutine ReadLakeName(lakeId)
      implicit none
      integer, intent(in) :: lakeId
      character(len=32) :: fname = "ReadLakeName"
      character(cx) :: fullname, tmpstr, msg
      character(len=32) :: lname
      integer :: nline, ii, id, error

      if (len_trim(lakeid_file)==0) then
         return
      end if
      call GetFullFileName(lakeid_file, fullname)
      open(unit=fid, file=trim(fullname), status="old", action="read", &
           iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun('Cannot open file ' // trim(fullname))
      end if
      read(unit=fid, fmt='(A7, I)', iostat=error) tmpstr, nline
      read(fid, "(A512)", iostat=error) tmpstr
      read(fid, "(A512)", iostat=error) tmpstr
      do ii = 1, nline, 1
         read(unit=fid, fmt=*, iostat=error) id, lname
         if (error/=0) then
            close(unit=fid)
            write(msg, "(A, I0)") "Reading stops at line ", ii
            call Endrun(fname, trim(msg))
         end if
         if (id==lakeId) then
            exit
         end if
      end do
      close(unit=fid)
      if (lakeId==id) then
         lake_info%name = lname
      else
         write(msg, "(A, I0)") "Cannot find name for ", lakeId
         call Endrun(fname, trim(msg))
      end if
   end subroutine 

   !------------------------------------------------------------------------------
   !
   ! Purpose: Read observation data for calibration.
   !
   !------------------------------------------------------------------------------
   subroutine CalcCostfunc4Var0D(filename, timeHist, ntday, varHist, &
                                 std, fcost)
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(in) :: timeHist(:)
      integer, intent(in) :: ntday
      real(r4), intent(in) :: varHist(:)
      real(r8), intent(in) :: std
      real(r8), intent(out) :: fcost
      character(len=32) :: fname = "CalcCostfunc4Var0D"
      character(cx) :: fullname, tmpstr, msg
      real(r8) :: fobs, fobs_std, ferr, fsim
      integer :: nline, nobs, error, ii, jj
      integer :: nt, date, idx, nrec

      nt = size(timeHist)
      call GetFullFileName(filename, fullname)
      open(unit=fid, file=trim(fullname), status="old", action="read", &
           iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun('Cannot open file ' // trim(filename))
      end if
      read(unit=fid, fmt='(A14, I)', iostat=error) tmpstr, nline
      read(fid, "(A512)", iostat=error) tmpstr
      read(fid, "(A512)", iostat=error) tmpstr
      fcost = 0.0_r8
      nobs = 0
      do ii = 1, nline, 1
         read(unit=fid, fmt=*, iostat=error) date, fobs, fobs_std 
         if (error/=0) then
            close(unit=fid)
            write(msg, "(A, I0)") "Reading stops at line ", ii
            call Endrun(fname, trim(msg))
         end if
         if (date<timeHist(1) .or. date>timeHist(nt)) then
            cycle 
         end if
         idx = BSEARCHQQ(LOC(date),LOC(timeHist),INT8(nt),SRT$INTEGER4)
         fsim = 0.0_r8
         nrec = 0
         do jj = max(1,idx-24*ntday), min(nt,idx+24*ntday), 1
            if (abs(timeHist(jj)-date)<ntday) then
               fsim = fsim + varHist(jj)
               nrec = nrec + 1
            end if
         end do
         ferr = max( 0.01*fobs_std*fobs, std )
         fcost = fcost + ((fobs-fsim/DBLE(nrec)) / ferr)**2
         nobs = nobs + 1
      end do
      close(unit=fid)
      fcost = fcost / DBLE(nobs)
   end subroutine

   subroutine CalcCostfunc4Var1D(filename, timeHist, ntday, depthHist, &
                                 varHist, std, fcost)
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(in) :: timeHist(:)
      integer, intent(in) :: ntday
      real(r8), intent(in) :: depthHist(:,:)
      real(r4), intent(in) :: varHist(:,:)
      real(r8), intent(in) :: std
      real(r8), intent(out) :: fcost
      character(len=32) :: fname = "CalcCostfunc4Var1D"
      character(cx) :: fullname, tmpstr, msg
      real(r8) :: fobs, fobs_std, ferr
      real(r8) :: fsim, val, par
      real(r8) :: depth, z1, z2
      integer :: nline, nobs, error, ii, jj
      integer :: nt, nz, date, idx, nrec
      integer :: iz1, iz2 

      nt = size(timeHist)
      nz = size(depthHist,1)
      call GetFullFileName(filename, fullname)
      open(unit=fid, file=trim(fullname), status="old", action="read", &
           iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun('Cannot open file ' // trim(filename))
      end if
      read(unit=fid, fmt='(A14, I)', iostat=error) tmpstr, nline
      read(fid, "(A512)", iostat=error) tmpstr
      read(fid, "(A512)", iostat=error) tmpstr
      fcost = 0.0_r8
      nobs = 0
      do ii = 1, nline, 1
         read(unit=fid, fmt=*, iostat=error) date, depth, fobs, fobs_std
         if (error/=0) then
            close(unit=fid)
            write(msg, "(A, I0)") "Reading stops at line ", ii
            call Endrun(fname, trim(msg))
         end if
         if (date<timeHist(1) .or. date>timeHist(nt)) then
            cycle
         end if
         idx = BSEARCHQQ(LOC(date),LOC(timeHist),INT8(nt),SRT$INTEGER4)
         if (depth<=depthHist(1,idx)) then
            iz1 = 1
            iz2 = 1
            par = 0.0_r8
         else if (depth>=depthHist(nz,idx)) then
            iz1 = nz
            iz2 = nz
            par = 0.0_r8
         else
            call DichotomySectionSearch(depthHist(:,idx), depth, iz1)
            iz2 = iz1 + 1
            z1 = depthHist(iz1,idx)
            z2 = depthHist(iz2,idx)
            par = (depth - z1) / (z2 - z1)
         end if
         fsim = 0.0_r8
         nrec = 0
         do jj = max(1,idx-24*ntday), min(nt,idx+24*ntday), 1
            if (abs(timeHist(jj)-date)<ntday) then
               val = varHist(iz1,jj)*(1.0-par) + varHist(iz1,jj)*par 
               fsim = fsim + val
               nrec = nrec + 1
            end if
         end do
         ferr = max( 0.01*fobs_std*fobs, std )
         fcost = fcost + ((fobs-fsim/DBLE(nrec)) / ferr)**2
         nobs = nobs + 1
      end do
      close(unit=fid)
      fcost = fcost / DBLE(nobs)
   end subroutine

   subroutine CalcCostfunc4AnnVar(filename, timeHist, varHist, &
                                  std, fcost)
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(in) :: timeHist(:)
      real(r8), intent(in) :: varHist(:)
      real(r8), intent(in) :: std
      real(r8), intent(out) :: fcost
      character(len=32) :: fname = "CalcCostfunc4AnnVar"
      character(cx) :: fullname, tmpstr, msg
      real(r8) :: fobs, fobs_std, ferr, fsim
      integer :: nline, nobs, error, ii
      integer :: nt, date, idx

      nt = size(timeHist)
      call GetFullFileName(filename, fullname)
      open(unit=fid, file=trim(fullname), status="old", action="read", &
           iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun('Cannot open file ' // trim(filename))
      end if
      read(unit=fid, fmt='(A14, I)', iostat=error) tmpstr, nline
      read(fid, "(A512)", iostat=error) tmpstr
      read(fid, "(A512)", iostat=error) tmpstr
      fcost = 0.0_r8
      nobs = 0
      do ii = 1, nline, 1
         read(unit=fid, fmt=*, iostat=error) date, fobs, fobs_std 
         if (error/=0) then
            close(unit=fid)
            write(msg, "(A, I0)") "Reading stops at line ", ii
            call Endrun(fname, trim(msg))
         end if
         if (date<timeHist(1) .or. date>timeHist(nt)) then
            cycle
         end if
         idx = BSEARCHQQ(LOC(date),LOC(timeHist),INT8(nt),SRT$INTEGER4)
         fsim = varHist(idx) 
         ferr = max( 0.01*fobs_std*fobs, std )
         fcost = fcost + ((fobs-fsim) / ferr)**2
         nobs = nobs + 1
      end do
      close(unit=fid)
      fcost = fcost / DBLE(nobs)
   end subroutine

   subroutine CalcCostfunc4AvgVar(filename, timeHist, ntday, depthHist, &
                                 varHist, std, fcost)
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(in) :: timeHist(:)
      integer, intent(in) :: ntday
      real(r8), intent(in) :: depthHist(:,:)
      real(r4), intent(in) :: varHist(:,:)
      real(r8), intent(in) :: std
      real(r8), intent(out) :: fcost
      character(len=32) :: fname = "CalcCostfunc4AvgVar"
      character(cx) :: fullname, tmpstr, msg
      real(r8) :: fobs, fobs_std, ferr
      real(r8) :: secchi, depth, fsim
      real(r8) :: val, wt, wt_sum
      integer :: nline, error, ii, jj, kk
      integer :: nt, nz, date, idx, nrec
      integer :: iz1, iz2, nobs

      nt = size(timeHist)
      nz = size(depthHist,1)
      call GetFullFileName(filename, fullname)
      open(unit=fid, file=trim(fullname), status="old", action="read", &
           iostat=error)
      if (error/=0) then
         close(unit=fid)
         call Endrun('Cannot open file ' // trim(filename))
      end if
      read(unit=fid, fmt='(A14, I)', iostat=error) tmpstr, nline
      read(fid, "(A512)", iostat=error) tmpstr
      read(fid, "(A512)", iostat=error) tmpstr
      fcost = 0.0_r8
      nobs = 0
      do ii = 1, nline, 1
         read(unit=fid, fmt=*, iostat=error) date, secchi, fobs, &
            fobs_std
         if (error/=0) then
            close(unit=fid)
            write(msg, "(A, I0)") "Reading stops at line ", ii
            call Endrun(fname, trim(msg))
         end if
         if (date<timeHist(1) .or. date>timeHist(nt)) then
            cycle
         end if
         depth = 2.0 * secchi
         idx = BSEARCHQQ(LOC(date),LOC(timeHist),INT8(nt),SRT$INTEGER4)
         if (depth>=depthHist(nz,idx)) then
            iz2 = nz
         else
            call DichotomySectionSearch(depthHist(:,idx), depth, iz2)
         end if
         do jj = 1, nz, 1
            if (depthHist(jj,idx)>=0) then
               iz1 = jj
               exit
            end if
         end do
         fsim = 0.0_r8
         nrec = 0
         do jj = max(1,idx-24*ntday), min(nt,idx+24*ntday), 1
            if (abs(timeHist(jj)-date)<ntday) then
               val = 0.0_r8
               wt_sum = 0.0_r8
               do kk = iz1, iz2, 1
                  if (kk==iz1) then
                     wt = 0.5*(depthHist(kk+1,jj)-depthHist(kk,jj)) 
                  else if (kk==iz2) then
                     wt = 0.5*(depthHist(kk,jj)-depthHist(kk-1,jj))
                  else
                     wt = 0.5*(depthHist(kk+1,jj)-depthHist(kk-1,jj))
                  end if
                  val = val + varHist(kk,jj) * wt
                  wt_sum = wt_sum + wt
               end do
               val = val / wt_sum 
               fsim = fsim + val
               nrec = nrec + 1
            end if
         end do
         ferr = max( 0.01*fobs_std*fobs, std )
         fcost = fcost + ((fobs-fsim/DBLE(nrec)) / ferr)**2
         nobs = nobs + 1
      end do
      close(unit=fid)
      fcost = fcost / DBLE(nobs)
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

   subroutine ReadSiteTSData(info, time, tstep, varname, odata)
      implicit none
      type(LakeInfo), intent(in) :: info
      type(SimTime), intent(in) :: time
      character(len=*), intent(in) :: tstep
      character(len=*), intent(in) :: varname
      real(r8), intent(out) :: odata(:)
       character(len=32) :: fname = "ReadSiteDTSData"
      character(cx) :: filename, fullname
      integer(kind=MPI_OFFSET_KIND) :: nstart(1)
      integer(kind=MPI_OFFSET_KIND) :: ncount(1)
      integer :: ncid, varid, date_varid
      integer :: JDN0, JDN1, JDNb, date(1)
      integer :: year, month, day
      real(r8) :: filled_value

      ! inquire data info
      if (len_trim(lakeid_file)==0) then
         write(filename, "(A,I0,A)") trim(forcing_dir), info%id, '.nc'
      else
         filename = trim(forcing_dir) // trim(info%name) // '.nc'
      end if
      call GetFullFileName(filename, fullname)
      if (trim(run_mode)=='regular') then
         call check( fname, nf90mpi_open(MPI_COMM_SELF, trim(fullname), &
                     NF90_NOWRITE, MPI_INFO_NULL, ncid) )
      else
         call check( fname, nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                     NF90_NOWRITE, MPI_INFO_NULL, ncid) )
      end if
      call check( nf90mpi_inq_varid(ncid, "date", date_varid) )
      call check( nf90mpi_inq_varid(ncid, trim(varname), varid) )
      call check( nf90mpi_get_att(ncid, varid, "_FillValue", filled_value) )

      ! read data start and end date
      call check( nf90mpi_get_var_all(ncid, date_varid, date, &
                  (/1_i8/), (/1_i8/)) )
      call YYMMDD2Date(date(1), year, month, day)
      call Date2JDN(year, month, day, JDNb)
      call Date2JDN(time%year0, time%month0, time%day0, JDN0)
      call Date2JDN(time%year1, time%month1, time%day1, JDN1)

      ! read data set
      if (trim(tstep)=='hour') then
         nstart = 24 * (JDN0 - JDNb) + 1
         ncount = 24 * (JDN1 - JDN0)
      else if (trim(tstep)=='day') then
         nstart = JDN0 - JDNb + 1
         ncount = JDN1 - JDN0
      end if
      call check( nf90mpi_get_var_all(ncid, varid, odata, nstart, ncount) )
      call check( nf90mpi_close(ncid) )

      !if (count(odata<=filled_value)>0) then
      !   call Endrun("Invalid forcing data " // trim(varname))
      !end if

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
   subroutine Read2DTSData(filename, time, loc, varname, odata)
      implicit none
      character(len=*), intent(in) :: filename
      type(SimTime), intent(in) :: time
      real(r8), intent(in) :: loc(2)
      character(len=*), intent(in) :: varname
      real(r8), intent(out) :: odata(:)
      character(len=32) :: fname = "Read2DTSData"
      character(cx) :: fullname
      character(len=32) :: date_units
      integer(kind=MPI_OFFSET_KIND) :: b1start(3), b1count(3)
      integer(kind=MPI_OFFSET_KIND) :: nreslon, nreslat
      integer :: ncid, varid, lon_dimid, lat_dimid
      integer :: lon_varid, lat_varid, date_varid
      integer :: lonindx, latindx, lonerr, laterr
      integer :: date(1), JDN0, JDN1, JDNb, timeIndx(2)
      integer :: year, month, day
      real(r8), allocatable :: tmplons(:), tmplats(:)
      real(r8), allocatable :: tmpArr(:,:,:)
      real(r8) :: plon, plat, filled_val

      plon = loc(1) 
      plat = loc(2) 
      call GetFullFileName(filename, fullname)
      call check( fname, nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                  NF90_NOWRITE, MPI_INFO_NULL, ncid) )
      call check( nf90mpi_inq_varid(ncid, "date", date_varid) )
      call check( nf90mpi_get_att(ncid, date_varid, "units", date_units) )
      lonerr = nf90mpi_inq_dimid(ncid, "longitude", lon_dimid)
      if (lonerr/=NF90_NOERR) then
         call check( nf90mpi_inq_dimid(ncid, "lon", lon_dimid) )
      end if
      laterr = nf90mpi_inq_dimid(ncid, "latitude", lat_dimid)
      if (laterr/=NF90_NOERR) then
         call check( nf90mpi_inq_dimid(ncid, "lat", lat_dimid) )
      end if
      call check( nf90mpi_inquire_dimension(ncid, lon_dimid, len = nreslon) )
      call check( nf90mpi_inquire_dimension(ncid, lat_dimid, len = nreslat) )
      allocate(tmplons(nreslon))
      allocate(tmplats(nreslat))
      if (lonerr==NF90_NOERR) then
         call check( nf90mpi_inq_varid(ncid, "longitude", lon_varid) )
      else
         call check( nf90mpi_inq_varid(ncid, "lon", lon_varid) )
      end if
      if (laterr==NF90_NOERR) then
         call check( nf90mpi_inq_varid(ncid, "latitude", lat_varid) )
      else
         call check( nf90mpi_inq_varid(ncid, "lat", lat_varid) )
      end if
      call check( nf90mpi_get_var_all(ncid, lon_varid, tmplons) )
      call check( nf90mpi_get_var_all(ncid, lat_varid, tmplats) )
      ! search for corresponding indice
      call DichotomySectionSearch(tmplons, plon, lonindx)
      call DichotomySectionSearch(tmplats, plat, latindx)
      call check( fname, nf90mpi_inq_varid(ncid, trim(varname), varid) )
      call check( fname, nf90mpi_get_att(ncid, varid, "_FillValue", &
                  filled_val) )

      call check( nf90mpi_get_var_all(ncid, date_varid, date, &
                  (/1_i8/), (/1_i8/)) ) 
      if (trim(date_units)=='YYYYMMDD') then
         call YYMMDD2Date(date(1), year, month, day)
         call Date2JDN(year, month, day, JDNb)
         call Date2JDN(time%year0, time%month0, time%day0, JDN0)
         call Date2JDN(time%year1, time%month1, time%day1, JDN1)
         timeIndx = (/JDN0-JDNb+1, JDN1-JDN0/)
      else if (trim(date_units)=='YYYYMM') then
         date = 100 * date + 1
         call YYMMDD2Date(date(1), year, month, day)
         if (time%day1==1) then 
            timeIndx = (/12*(time%year0-year)+time%month0, &
               12*(time%year1-time%year0)-time%month0+time%month1/)
         else
            timeIndx = (/12*(time%year0-year)+time%month0, &
               12*(time%year1-time%year0)-time%month0+time%month1+1/)
         end if
      else if (trim(date_units)=='YYYY') then
         date = 10000 * date + 101
         call YYMMDD2Date(date(1), year, month, day)
         if (time%month1==1 .and. time%day1==1) then
            timeIndx = (/time%year0-year+1, time%year1-time%year0/)
         else
            timeIndx = (/time%year0-year+1, time%year1-time%year0+1/)
         end if
      end if

      allocate(tmpArr(1,1,timeIndx(2)))
      b1start = (/lonindx, latindx, timeIndx(1)/)
      b1count = (/1, 1, timeIndx(2)/)
      call check( fname, nf90mpi_get_var_all(ncid, varid, tmpArr, &
                  b1start, b1count) )
      odata = tmpArr(1,1,:)
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
   subroutine ReadStaticRealData(filename, loc, varname, interp, odata)
      implicit none
      character(len=*), intent(in) :: filename
      real(r8), intent(in) :: loc(2)
      character(len=*), intent(in) :: varname
      integer, intent(in) :: interp
      real(r8), intent(out) :: odata
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND), parameter :: nrange2 = 10
      integer(kind=MPI_OFFSET_KIND) :: b1start(2), b1count(2)
      integer(kind=MPI_OFFSET_KIND) :: b2start(2), b2count(2)
      integer(kind=MPI_OFFSET_KIND) :: b3start(2), b3count(2)
      integer(kind=MPI_OFFSET_KIND) :: b4start(2), b4count(2)
      integer(kind=MPI_OFFSET_KIND) :: nreslon, nreslat
      integer(kind=MPI_OFFSET_KIND) :: idx, jdx, nrange
      integer, allocatable :: target_lonid(:), target_latid(:)
      integer :: ncid, varid, lon_varid, lat_varid
      integer :: lon_dimid, lat_dimid, ii, jj
      integer :: ncount, lonindx, latindx
      integer :: lonerr, laterr
      real(r8), allocatable :: tmplons(:), tmplats(:)
      real(r8), allocatable :: weights(:,:)
      real(r8) :: lon, lat, plon, plat
      real(r8) :: numer, denom, dist
      real(r4), allocatable :: tmpArr(:,:)
      real(r4) :: filled_val

      plon = loc(1)
      plat = loc(2)
      call GetFullFileName(filename, fullname)
      call check( nf90mpi_open(MPI_COMM_WORLD, trim(fullname), NF90_NOWRITE, &
                  MPI_INFO_NULL, ncid) )

      ! data set scope
      lonerr = nf90mpi_inq_dimid(ncid, "longitude", lon_dimid)
      if (lonerr/=NF90_NOERR) then
         call check( nf90mpi_inq_dimid(ncid, "lon", lon_dimid) )
      end if
      laterr = nf90mpi_inq_dimid(ncid, "latitude", lat_dimid)
      if (laterr/=NF90_NOERR) then
         call check( nf90mpi_inq_dimid(ncid, "lat", lat_dimid) )
      end if
      call check( nf90mpi_inquire_dimension(ncid, lon_dimid, len = nreslon) )
      call check( nf90mpi_inquire_dimension(ncid, lat_dimid, len = nreslat) )
      allocate(tmplons(nreslon))
      allocate(tmplats(nreslat))
      if (lonerr==NF90_NOERR) then
         call check( nf90mpi_inq_varid(ncid, "longitude", lon_varid) )
      else
         call check( nf90mpi_inq_varid(ncid, "lon", lon_varid) )
      end if
      if (laterr==NF90_NOERR) then
         call check( nf90mpi_inq_varid(ncid, "latitude", lat_varid) )
      else
         call check( nf90mpi_inq_varid(ncid, "lat", lat_varid) )
      end if
      call check( nf90mpi_get_var_all(ncid, lon_varid, tmplons) )
      call check( nf90mpi_get_var_all(ncid, lat_varid, tmplats) )
      call DichotomySectionSearch(tmplons, plon, lonindx)
      call DichotomySectionSearch(tmplats, plat, latindx)

      call check( nf90mpi_inq_varid(ncid, trim(varname), varid) )
      call check( nf90mpi_get_att(ncid, varid, "_FillValue", filled_val) )

      if (interp==Ncinterp) then

         ncount = 0
         nrange = nrange2
         do while (ncount==0 .and. nrange<=8*nrange2)
            if (allocated(tmpArr))  deallocate(tmpArr)
            if (allocated(target_lonid))  deallocate(target_lonid)
            if (allocated(target_latid))  deallocate(target_latid)
            if (allocated(weights))  deallocate(weights)
            allocate(tmpArr(nrange,nrange))
            allocate(weights(nrange,nrange))
            allocate(target_lonid(nrange))
            allocate(target_latid(nrange))

            ! construct data blocks
            call GetTargetCoordinate(INT8(lonindx), INT8(latindx), nreslon, &
                     nreslat, nrange/2, target_lonid, target_latid, b1start, &
                     b1count, b2start, b2count, b3start, b3count, b4start, &
                     b4count)

            ! read block1 data 
            idx = b1count(1)
            jdx = b1count(2)
            call check( nf90mpi_get_var_all(ncid, varid, tmpArr(1:idx,1:jdx), &
                        b1start, b1count) )
            ! read block2 data
            idx = nrange - b2count(1) + 1
            jdx = b2count(2)
            call check( nf90mpi_get_var_all(ncid, varid, &
                        tmpArr(idx:nrange,1:jdx), b2start, b2count) )
            ! read block3 data
            idx = b3count(1)
            jdx = nrange - b3count(2) + 1
            call check( nf90mpi_get_var_all(ncid, varid, &
                        tmpArr(1:idx,jdx:nrange), b3start, b3count) )
            ! read block4 data
            idx = nrange - b4count(1) + 1
            jdx = nrange - b4count(2) + 1
            call check( nf90mpi_get_var_all(ncid, varid, &
                        tmpArr(idx:nrange,jdx:nrange), b4start, b4count) )

            ! construct weights matrix
            do ii = 1, nrange, 1
               lon = tmplons(target_lonid(ii))
               do jj = 1, nrange, 1
                  lat = tmplats(target_latid(jj))
                  if (lat==plat .and. lon==plon) then
                     if (tmpArr(ii,jj)/=filled_val) then
                        weights = 0.0
                        weights(ii,jj) = 1.0
                        exit
                     else
                        weights(ii,jj) = 0.0
                     end if
                  else
                     dist = CalcGreatCircleDistance(plon, plat, lon, lat)
                     weights(ii,jj) = dist**(-2)
                  end if
               end do
            end do

            ! interpolate
            numer = 0.0
            denom = 0.0
            ncount = 0
            do ii = 1, nrange, 1
               do jj = 1, nrange, 1
                  if (tmpArr(ii,jj)/=filled_val) then
                     numer = numer + tmpArr(ii,jj) * weights(ii,jj)
                     denom = denom + weights(ii,jj)
                     ncount = ncount + 1
                  end if
               end do
            end do

            nrange = 2 * nrange
         end do
         call check( nf90mpi_close(ncid) )

         if (ncount>0) then
            odata = numer / denom
         else
            call Endrun(trim(varname) // ": cannot find valid scope")
         end if
      else
         allocate(tmpArr(1,1))
         b1start = (/lonindx, latindx/)
         b1count = (/1, 1/)
         call check( nf90mpi_get_var_all(ncid, varid, tmpArr, b1start, &
                     b1count) )
         odata = tmpArr(1,1)
         call check( nf90mpi_close(ncid) )
      end if

      if (allocated(tmpArr))        deallocate(tmpArr)
      if (allocated(target_lonid))  deallocate(target_lonid)
      if (allocated(target_latid))  deallocate(target_latid)
      if (allocated(weights))       deallocate(weights)
      if (allocated(tmplons))       deallocate(tmplons)
      if (allocated(tmplats))       deallocate(tmplats)
      if (DEBUG .and. masterproc) then
         print *, "Read " // trim(varname) // " from " // trim(fullname)
      end if
   end subroutine

   subroutine ReadStaticByteData(filename, loc, varname, odata)
      implicit none
      character(len=*), intent(in) :: filename
      real(r8), intent(in) :: loc(2)
      character(len=*), intent(in) :: varname
      real(r8), intent(out) :: odata
      character(cx) :: fullname
      integer(kind=MPI_OFFSET_KIND) :: nreslon, nreslat
      integer(kind=MPI_OFFSET_KIND) :: nstart(2), ncount(2)
      integer :: ncid, varid, lon_varid, lat_varid
      integer :: lon_dimid, lat_dimid
      integer :: lonindx, latindx
      integer :: lonerr, laterr
      integer(i1) :: tmpArr(1,1)
      real(r8), allocatable :: tmplons(:), tmplats(:)
      real(r8) :: plon, plat

      plon = loc(1)
      plat = loc(2)
      call GetFullFileName(filename, fullname)
      call check( nf90mpi_open(MPI_COMM_WORLD, trim(fullname), NF90_NOWRITE, &
                  MPI_INFO_NULL, ncid) )

      ! data set scope
      lonerr = nf90mpi_inq_dimid(ncid, "longitude", lon_dimid)
      if (lonerr/=NF90_NOERR) then
         call check( nf90mpi_inq_dimid(ncid, "lon", lon_dimid) )
      end if
      laterr = nf90mpi_inq_dimid(ncid, "latitude", lat_dimid)
      if (laterr/=NF90_NOERR) then
         call check( nf90mpi_inq_dimid(ncid, "lat", lat_dimid) )
      end if
      call check( nf90mpi_inquire_dimension(ncid, lon_dimid, len = nreslon) )
      call check( nf90mpi_inquire_dimension(ncid, lat_dimid, len = nreslat) )
      allocate(tmplons(nreslon))
      allocate(tmplats(nreslat))
      if (lonerr==NF90_NOERR) then
         call check( nf90mpi_inq_varid(ncid, "longitude", lon_varid) )
      else
         call check( nf90mpi_inq_varid(ncid, "lon", lon_varid) )
      end if
      if (laterr==NF90_NOERR) then
         call check( nf90mpi_inq_varid(ncid, "latitude", lat_varid) )
      else
         call check( nf90mpi_inq_varid(ncid, "lat", lat_varid) )
      end if
      call check( nf90mpi_inq_varid(ncid, trim(varname), varid) )
      call check( nf90mpi_get_var_all(ncid, lon_varid, tmplons) )
      call check( nf90mpi_get_var_all(ncid, lat_varid, tmplats) )
      call DichotomySectionSearch(tmplons, plon, lonindx)
      call DichotomySectionSearch(tmplats, plat, latindx)
      deallocate(tmplons)
      deallocate(tmplats)

      nstart = (/lonindx, latindx/)
      ncount = (/1, 1/)
      call check( nf90mpi_get_var_all(ncid, varid, tmpArr, nstart, &
                  ncount) )
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
   subroutine ReadGlobalTSData(filename, time, varname, odata)
      implicit none
      character(len=*), intent(in) :: filename
      type(SimTime), intent(in) :: time
      character(len=*), intent(in) :: varname
      real(r8), intent(out) :: odata(:)
      character(len=32) :: fname = "ReadGlobalTSData"
      integer(kind=MPI_OFFSET_KIND) :: nstart(1), ncount(1)
      character(cx) :: fullname
      character(len=32) :: date_units
      integer :: ncid, varid, date_varid
      integer :: date(1), JDN0, JDN1, JDNb, timeIndx(2)
      integer :: year, month, day

      call GetFullFileName(filename, fullname)
      call check( fname, nf90mpi_open(MPI_COMM_WORLD, trim(fullname), &
                  NF90_NOWRITE, MPI_INFO_NULL, ncid) )
      call check( nf90mpi_inq_varid(ncid, "date", date_varid) )
      call check( nf90mpi_get_att(ncid, date_varid, "units", date_units) )
      call check( fname, nf90mpi_inq_varid(ncid, trim(varname), varid) )
      call check( nf90mpi_get_var_all(ncid, date_varid, date, &
                  (/1_i8/), (/1_i8/)) )
      if (trim(date_units)=='YYYYMMDD') then
         call YYMMDD2Date(date(1), year, month, day)
         call Date2JDN(year, month, day, JDNb)
         call Date2JDN(time%year0, time%month0, time%day0, JDN0)
         call Date2JDN(time%year1, time%month1, time%day1, JDN1)
         timeIndx = (/JDN0-JDNb+1, JDN1-JDN0/)
      else if (trim(date_units)=='YYYYMM') then
         date = 100 * date + 1
         call YYMMDD2Date(date(1), year, month, day)
         if (time%day1==1) then
            timeIndx = (/12*(time%year0-year)+time%month0, &
               12*(time%year1-time%year0)-time%month0+time%month1/)
         else
            timeIndx = (/12*(time%year0-year)+time%month0, &
               12*(time%year1-time%year0)-time%month0+time%month1+1/)
         end if
      else if (trim(date_units)=='YYYY') then
         date = 10000 * date + 101
         call YYMMDD2Date(date(1), year, month, day)
         if (time%month1==1 .and. time%day1==1) then
            timeIndx = (/time%year0-year+1, time%year1-time%year0/)
         else
            timeIndx = (/time%year0-year+1, time%year1-time%year0+1/)
         end if
      end if
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
   subroutine CreateOutputFile0d(nz, varname, longname, units)
      implicit none
      integer, intent(in) :: nz
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units
      character(cx) :: fullname, histstr
      character(cx) :: str1, str2
      integer :: ncid, varid, cmode
      integer :: z_dimid, lake_dimid
      integer :: gmtime(6), err

      call GetArchiveFullname(varname, fullname)
      if (masterproc) then
         call GetGMTime(gmtime)
         write(str1,"(I4, '-', I2.2, '-', I2.2)") gmtime(1), gmtime(2), &
               gmtime(3)
         write(str2,"(' ', I2.2, ':', I2.2, ':', I2.2, A)") gmtime(4), &
               gmtime(5), gmtime(6), " GMT from ALBM v2.0 by Zeli Tan"
         histstr = trim(str1) // trim(str2)
      end if
      call MPI_BCAST(histstr, len(histstr), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)

      ! Create the file.
      cmode = IOR(NF90_CLOBBER, NF90_64BIT_OFFSET)
      call check( nf90mpi_create(MPI_COMM_WORLD, trim(fullname), &
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
      integer :: ncid, varid, cmode
      integer :: time_dimid, lake_dimid
      integer :: gmtime(6), simday, err

      call GetArchiveFullname(time, varname, fullname)
      if (masterproc) then
         call GetGMTime(gmtime)
         write(str1,"(I4, '-', I2.2, '-', I2.2)") gmtime(1), gmtime(2), &
               gmtime(3)
         write(str2,"(' ', I2.2, ':', I2.2, ':', I2.2, A)") gmtime(4), &
               gmtime(5), gmtime(6), " GMT from ALBM v2.0 by Zeli Tan"         
         histstr = trim(str1) // trim(str2)
      end if
      call MPI_BCAST(histstr, len(histstr), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)

      write(str1,"('Records start since ', I4, '-', I2.2, '-', I2.2)") &
         time%year0, time%month0, time%day0
      timestr = trim(str1) // " 00:00:00 GMT"
      simday = CalcRunningDays(time)

      ! Create the file.
      cmode = IOR(NF90_CLOBBER, NF90_64BIT_OFFSET)
      call check( nf90mpi_create(MPI_COMM_WORLD, trim(fullname), &
                  cmode, MPI_INFO_NULL, ncid) )

      ! Define the dimensions
      call check( nf90mpi_def_dim(ncid, 'Time', INT8(simday), time_dimid) )
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

   subroutine CreateOutputFile2d(time, nz, varname, longname, units, defval)
      implicit none
      Type(SimTime), intent(in)  :: time
      integer, intent(in) :: nz
      character(len=*), intent(in) :: varname
      character(len=*), intent(in) :: longname
      character(len=*), intent(in) :: units
      real(r4), intent(in) :: defval
      character(cx) :: fullname, histstr, timestr
      character(cx) :: str1, str2
      integer :: ncid, varid, cmode
      integer :: z_dimid, time_dimid, lake_dimid
      integer :: gmtime(6), simday, err

      call GetArchiveFullname(time, varname, fullname)
      if (masterproc) then
         call GetGMTime(gmtime)
         write(str1,"(I4, '-', I2.2, '-', I2.2)") gmtime(1), gmtime(2), &
               gmtime(3)
         write(str2,"(' ', I2.2, ':', I2.2, ':', I2.2, A)") gmtime(4), &
               gmtime(5), gmtime(6), " GMT from ALBM v2.0 by Zeli Tan"
         histstr = trim(str1) // trim(str2)
      end if
      call MPI_BCAST(histstr, len(histstr), MPI_CHARACTER, 0, &
                     MPI_COMM_WORLD, err)

      write(str1,"('Records start since ', I4, '-', I2.2, '-', I2.2)") &
         time%year0, time%month0, time%day0
      timestr = trim(str1) // " 00:00:00 GMT"
      simday = CalcRunningDays(time)

      ! Create the file.
      cmode = IOR(NF90_CLOBBER, NF90_64BIT_OFFSET)
      call check( nf90mpi_create(MPI_COMM_WORLD, trim(fullname), &
                  cmode, MPI_INFO_NULL, ncid) )

      ! Define the dimensions
      call check( nf90mpi_def_dim(ncid, 'Z', INT8(nz), z_dimid) )
      call check( nf90mpi_def_dim(ncid, 'Time', INT8(simday), time_dimid) )
      call check( nf90mpi_def_dim(ncid, 'Lake', NFMPI_UNLIMITED, lake_dimid) )

      ! Define the coordinate variables
      call check( nf90mpi_put_att(ncid, NF90_GLOBAL, "history", &
                  trim(histstr)) )
      call check( nf90mpi_put_att(ncid, NF90_GLOBAL, "description", &
                  trim(timestr)) )

      ! ! Define the netCDF variables and assign attributes for state variables. 
      call DefNcVariable(ncid, (/z_dimid, time_dimid, lake_dimid/), &
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
