module assimilation_mod
!---------------------------------------------------------------------------------
! Purpose: This module runs data assimilation on each specified lake.
!
!---------------------------------------------------------------------------------
   use shr_ctrl_mod
   use shr_param_mod
   use shr_typedef_mod
   use sim_coupler_mod
   use math_utilities_mod, only : init_cond_random_seed
   use read_data_mod
   use epfm_da_mod
   use epfm_bridge_mod
#ifdef USE_INTEL_COMPILER
   use ifport
#endif
   use mpi

   private
   public :: RunAssimilation
   ! current particle Id
   integer :: cur_particle

contains
   subroutine RunAssimilation(taskid, numprocs, arg)
      implicit none
      integer, intent(in) :: taskid
      integer, intent(in) :: numprocs
      character(len=*), intent(in) :: arg
      integer, allocatable :: partIds(:)
      integer :: partId, err
      integer :: ii, jj, itmp

      ! only support numprocs == NPART
      if (numprocs /= NPART) then
         print *, "Number of particles and processors: ", NPART, numprocs
         call Endrun("Number of processors must be set to number of particles")
      end if

      allocate(partIds(numprocs))

      if (masterproc) then
         print "(A, I0, A)", 'Run ', NPART, ' particles'
      end if

      call DoAssimilationWarmup()

      if (masterproc) then
         partIds = (/(ii, ii = 1, numprocs)/)
      end if
      call MPI_BCAST(partIds, numprocs, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
      partId = partIds(taskid+1)
      if (masterproc .and. partId/=1) then
         call Endrun("partId 1 must be on the root processor")  
      end if
      call DASimulation(partId, taskid, err)
      print "(A, I0, A, I0, A, I0)", "Particle ", partId, ", processor ", &
            taskid, ", error ", err

      deallocate(partIds)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: create simulation archive files 
   !
   !------------------------------------------------------------------------------
   subroutine DoAssimilationWarmup()
      implicit none
      type(SimTime) :: time
      
      time = SimTime(Start_Year, Start_Month, Start_Day, 0, &
                     End_Year, End_Month, End_Day, 0)
      if (masterproc) then
         ! create output files
         call CreateOutputFile(time, NWLAYER+1, 'zw', 'Z', 'water layer depth', &
                              'm', -9999.0_r4) 
         call CreateOutputFile(time, NSLAYER+1, 'zs', 'sediment layer depth', 'm')
         call CreateOutputFile(time, NWLAYER+1, 'Az', 'Z', 'water layer ' // & 
                              'cross-section area', 'm^2', -9999.0_r4);
         call CreateOutputFile(time, NWLAYER+1, 'colindx', 'water layer ' // &
                              'connected sediment column index', 'index')
         if (Thermal_Module) then
            call CreateOutputFile(time, 'snowthick', 'Snow Thickness', &
                                 'm', -9999.0_r4)
            call CreateOutputFile(time, 'icethick', 'Ice Thickness', &
                                 'm', -9999.0_r4)
            call CreateOutputFile(time, 'sensheatf', 'Sensible Heat Flux ' // &
                                 'at Lake-Atmosphere Interface', 'W m-2', &
                                 -9999.0_r4)
            call CreateOutputFile(time, 'latentheatf', 'Latent Heat Flux ' // &
                                 'at Lake-Atmosphere Interface', 'W m-2', &
                                 -9999.0_r4)
            call CreateOutputFile(time, 'momf', 'Momentum Flux at ' // &
                                 'Lake-Atmosphere Interface', 'kg m-1 s-2', &
                                 -9999.0_r4)
            call CreateOutputFile(time, 'lwup', 'Upward Long-Wave ' // &
                                 'Radiation Flux at Lake-Atmosphere Interface', &
                                 'W m-2', -9999.0_r4)
            call CreateOutputFile(time, 'lakeheatf', 'Downward Heat Flux ' // &
                                 'at Lake-Atmosphere Interface', 'W m-2', &
                                 -9999.0_r4)
            !call CreateOutputFile(time, 'swdw', 'downward shortwave radiation', &
            !                     'W/m2', -9999.0_r4)
            call CreateOutputFile(time, 'swup', 'Upward Short-Wave ' // &
                                 'Radiation Flux at Lake-Atmosphere Interface', &
                                 'W m-2', -9999.0_r4)
            call CreateOutputFile(time, 'sedheatf', 'Sediment Upward Heat ' // &
                                 'Flux at Lake-Sediment Interface', 'W m-2', &
                                 -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'watertemp', 'Z', &
                                 'Temperature of Lake Water', 'K', -9999.0_r4)
            !call CreateOutputFile(time, NWLAYER+1, 'turbdiffheat', 'Z', &
            !                     'Turbulent diffusivity of heat', 'm2/s', &
            !                     -9999.0_r4)
            call CreateOutputFile(time, NSCOL, NSLAYER+1, 'sedtemp', 'COL', 'Z', &
                                 'Temperature of Lake Sediment', 'K', -9999.0_r4)
         end if
         if (Hydro_Module) then
            call CreateOutputFile(time, 'Sw', 'Water storage', 'm^3', -9999.0_r4)
            call CreateOutputFile(time, 'Qwt', 'Outflow water temperature', &
                                 'm', -9999.0_r4)
            call CreateOutputFile(time, 'Qch4', 'Outflow CH4 concentration', &
                                 'mol m-3', -9999.0_r4)
            call CreateOutputFile(time, 'Qco2', 'Outflow CO2 concentration', &
                                 'mol m-3', -9999.0_r4)
            call CreateOutputFile(time, 'Qo2', 'Outflow O2 concentration', &
                                 'mol m-3', -9999.0_r4)
            call CreateOutputFile(time, 'Qsrp', 'Outflow SRP concentration', &
                                 'mol m-3', -9999.0_r4)
         end if
         if (Carbon_Module) then
            call CreateOutputFile(time, 'ch4df', 'surface methane diffusion flux', &
                                 'mol m-2 s-1', -9999.0_r4)
            call CreateOutputFile(time, 'gpp', 'total gross primary production', &
                                 'gC m-2 s-1', -9999.0_r4)
            call CreateOutputFile(time, 'npp', 'total net primary production', &
                                 'gC m-2 s-1', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'ch4oxid', 'Z', &
                                 'Methane Oxidation', 'mol m-3 s-1', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'och4prod', 'Z', 'oxic ' // &
                                 'methane production', 'mol m-3 s-1', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'ch4conc', 'Z', 'dissolved ' // &
                                 'CH4 concentration', 'mol m-3', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'o2conc', 'Z', 'dissolved ' // &
                                 'O2 concentration', 'mol m-3', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'n2conc', 'Z', 'dissolved ' // &
                                 'N2 concentration', 'mol m-3', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'co2conc', 'Z', 'dissolved ' // &
                                 'CO2 concentration', 'mol m-3', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'srp', 'Z', 'soluable ' // &
                                 'reactive P concentration', 'mol m-3', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'chla', 'Z', 'chlorophyll ' // &
                                 'concentration', 'g m-3', -9999.0_r4)
            call CreateOutputFile(time, NSCOL, 'cdep', 'COL', 'active organic ' // &
                                 'carbon deposition', 'gC m-2 s-1', -9999.0_r4)
            call CreateOutputFile(time, NPOC, NWLAYER+1, 'biomass', 'POC', 'Z', &
                                 'phytoplankton biomass', 'gC m-3', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'bveg', 'Z', 'submerged ' // &
                                 'macrophyte biomass', 'gC m-2', -9999.0_r4)
         end if
         if (Diagenesis_Module) then
            call CreateOutputFile(time, NSCOL, 'sedch4df', 'COL', 'sediment methane ' // &
                                 'diffusion flux', 'mol m-2 s-1', -9999.0_r4)
            call CreateOutputFile(time, NSCOL, 'sedch4eb', 'COL', 'sediment methane ' // &
                                 'ebullition flux', 'mol m-2 s-1', -9999.0_r4)
            call CreateOutputFile(time, NSCOL, 'ch4prod', 'COL', 'total sediment  ' // &
                                 'methane production', 'mol m-2 s-1', -9999.0_r4)
            call CreateOutputFile(time, NSCOL, NSLAYER+1, 'pascarb', 'COL', 'Z', &
                                 'sediment available passive carbon density', &
                                 'gC m-3', -9999.0_r4)
            call CreateOutputFile(time, NSCOL, NSLAYER+1, 'actcarb', 'COL', 'Z', &
                                 'sediment available active carbon density', &
                                 'gC m-3', -9999.0_r4)
            call CreateOutputFile(time, NSCOL, NSLAYER+1, 'oldcarb', 'COL', 'Z', &
                                 'sediment available permafrost carbon density', &
                                 'gC m-3', -9999.0_r4)
            call CreateOutputFile(time, NSCOL, NSLAYER+1, 'sedrdx', 'COL', 'Z', &
                                 'sediment redox potential', 'mV', -9999.0_r4)
         end if
         if (Bubble_Module) then
            call CreateOutputFile(time, 'ch4eb', 'surface methane ebullition flux', &
                                 'mol m-2 s-1', -9999.0_r4)
            call CreateOutputFile(time, 'icebch4', 'total CH4 in ice-trapped ' // &
                                 'bubbles', 'mole', -9999.0_r4)
         end if
      end if
   end subroutine

   subroutine DASimulation(partId, rank, error)
      implicit none
      integer, intent(in) :: partId
      integer, intent(in) :: rank
      integer, intent(out) :: error
      type(SimTime) :: time, otime
      real(r8) :: OptParams(NPARAM)
      integer :: i4ret, lakeId

#ifdef USE_INTEL_COMPILER
      i4ret = SIGNALQQ(SIG$FPE, hand_fpe)
#endif
      ! read lake information (i.e. depth, location ...)
      lakeId = lake_range(1)
      call ReadLakeInfo(lakeId, partId)
      call ReadOptimumParameters(OptParams)
      call LoadSensitiveParameters(OptParams)
      time = SimTime(Start_Year, Start_Month, Start_Day, 0, &
                     End_Year, End_Month, End_Day, 0)
      otime = time 
      call InitializeSimulation()
      call EPFM_Initialize(lakeId, time) 
      call ModelRun(lakeId, partId, rank, time, otime, error)
      call ArchiveAsssimilationOutput(partId, otime)
      call EPFM_Finalize()
      call FinalizeSimulation()
   end subroutine

   subroutine ModelRun(lakeId, partId, rank, time, otime, error)
      implicit none
      integer, intent(in) :: lakeId
      integer, intent(in) :: partId 
      integer, intent(in) :: rank
      type(SimTime), intent(in) :: time
      type(SimTime), intent(in) :: otime
      integer, intent(out) :: error
      type(SimTime) :: window
      type(EPFM_Obs) :: obs_c, obs_p
      integer :: ii, n_assim_times

      if (len_trim(restart_file)>0) then
         call ExtractRestartStates(lakeId)
      else
         call Endrun('restart_file cannot be empty in the DA mode')
      end if
      
      call init_cond_random_seed(12345, partId)
      call UpdateEPFMConfig(lake_info%kext)
      call PerturbTemperatureProfile(epfm_cfg_ini, partId) 

      n_assim_times = size(epfm_obs4da)
      ! Run simulation during the interested period
      do ii = 1, n_assim_times, 1
         if (ii==1) then
            obs_c = epfm_obs4da(ii)
            window = SimTime(time%year0, time%month0, time%day0, time%hour0, &
               obs_c%year, obs_c%month, obs_c%day, obs_c%hour+1)
            call EvolveParameters(epfm_cfg_ini, obs_c%month)
         else
            obs_c = epfm_obs4da(ii)
            obs_p = epfm_obs4da(ii-1)
            window = SimTime(obs_p%year, obs_p%month, obs_p%day, obs_p%hour, &
               obs_c%year, obs_c%month, obs_c%day, obs_c%hour+1)
            call EvolveParameters(epfm_cfg, obs_c%month)
         end if
         call ModuleCoupler(partId, time, window, otime, .False., error)
         if (error/=0) then
            exit
         end if
         call CopyLakeStateToEPFM() 
         ! assimilation
         if (obs_c%has_lwst .or. obs_c%has_secchi) then
            call GatherParticlesToRoot(rank)
            if (masterproc) then
               call ArraysToEPFMParticles() 
               call EPFM_Assimilate(obs_c, ObsOperator)
               call EPFMParticlesToArrays()
            end if
            call ScatterParticlesFromRoot(rank)
         end if
         call SaveEPFMToLakeState()
      end do

      if (error/=0) then
         call InitializeModelOutputs()
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: some utilities for exceptions
   !
   !------------------------------------------------------------------------------
#ifdef USE_INTEL_COMPILER
   function hand_fpe(sigid, except)
      !DEC$ ATTRIBUTES C :: hand_fpe
      use ifport
      INTEGER(4) :: hand_fpe
      INTEGER(2) :: sigid, except

      if (sigid/=SIG$FPE) then
         print "('The hand_fpe is not for signal ', I0)", sigid
         hand_fpe = 1
         return
      end if
      select case(except)
         case( FPE$INVALID )
            print *, ' Floating point exception: Invalid number'
         case( FPE$DENORMAL )
            print *, ' Floating point exception: Denormalized number'
         case( FPE$ZERODIVIDE )
            print *, ' Floating point exception: Zero divide'
         case( FPE$OVERFLOW )
            print *, ' Floating point exception: Overflow'
         case( FPE$UNDERFLOW )
            print *, ' Floating point exception: Underflow'
         case( FPE$INEXACT )
            print *, ' Floating point exception: Inexact precision'
         case default
            print *, ' Floating point exception: Non-IEEE type'
      end select
      !CALL TRACEBACKQQ(trim(header), USER_EXIT_CODE=-1)
      print *, 'failed particle ', cur_particle
      hand_fpe = 1
   end function
#endif

end module assimilation_mod
