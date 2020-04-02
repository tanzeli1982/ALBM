module simulation_mod
!---------------------------------------------------------------------------------
! Purpose: This module runs simulations locally or regionally.
!  There are two main components in this model: thermal module and bubbling module.
!  Thermal module includes water thermal sub module and soil thermal sub module.
!  Bubbling module includes gas production module and gas transportation module.
!
!---------------------------------------------------------------------------------
   use shr_ctrl_mod
   use shr_param_mod
   use shr_typedef_mod
   use sim_coupler_mod
   use read_data_mod

   private
   public :: RunRegular

contains
   subroutine RunRegular(arg)
      implicit none
      character(len=*), intent(in) :: arg
      integer :: err, istep, ii, itmp
      integer :: ntotlake, jj

      print "(A, I0, A, I0)", 'Run lake ', lake_id

      call DoSimulationWarmup()

      call RegularSimulation(lake_id, err)
      print "(A, I0, A, I0)", "Lake ", lake_id, ", error ", err

   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: create simulation archive files 
   !
   !------------------------------------------------------------------------------
   subroutine DoSimulationWarmup()
      implicit none
      type(SimTime) :: time
      
      time = SimTime(Start_Year, Start_Month, Start_Day, End_Year, &
                     End_Month, End_Day)
      call CreateOutputFile(lake_id, NWLAYER+1, 'zw', 'water layer depth', 'm') 
      call CreateOutputFile(lake_id, time, 'snowthick', 'snow cover thickness', &
                            'm', -9999.0_r4)
      call CreateOutputFile(lake_id, time, 'icethick', 'ice cover thickness', &
                            'm', -9999.0_r4)
      call CreateOutputFile(lake_id, time, 'sensheatf', 'upward sensible ' // &
                            'heat flux', 'W m-2', -9999.0_r4)
      call CreateOutputFile(lake_id, time, 'latentheatf', 'upward latent ' // &
                            'heat flux', 'W m-2', -9999.0_r4)
      call CreateOutputFile(lake_id, time, 'momf', 'momentum energy flux', &
                            'kg m-1 s-2', -9999.0_r4)
      call CreateOutputFile(lake_id, time, 'lwup', 'upward longwave radiation', &
                            'W m-2', -9999.0_r4)
      call CreateOutputFile(lake_id, time, 'lakeheatf', 'downward net heat flux', &
                            'W m-2', -9999.0_r4)
      call CreateOutputFile(lake_id, time, 'swdw', 'downward shortwave radiation', &
                            'W m-2', -9999.0_r4)
      call CreateOutputFile(lake_id, time, 'swup', 'upward shortwave radiation', &
                            'W m-2', -9999.0_r4)
      call CreateOutputFile(lake_id, time, 'sedheatf', 'upward sediment heat ' // &
                            'flux', 'W m-2', -9999.0_r4)
      call CreateOutputFile(lake_id, time, NWLAYER+1, 'watertemp', 'water ' // &
                            'temperature', 'K', -9999.0_r4)
      call CreateOutputFile(lake_id, time, NWLAYER+1, 'turbdiffheat', &
                            'Turbulent diffusivity of heat', 'm2 s-1', &
                            -9999.0_r4)
   end subroutine

   subroutine RegularSimulation(lakeId, error)
      implicit none
      integer, intent(in) :: lakeId
      integer, intent(out) :: error
      type(SimTime) :: time, spinup
      real(r8) :: OptParams(NPARAM)

      ! read lake information (i.e. depth, location ...)
      call ReadLakeName(lakeId)
      call ReadOptimumParameters(OptParams)
      call LoadSensitiveParameters(OptParams)
      call ReadLakeInfo(lakeId)
      time = SimTime(Start_Year, Start_Month, Start_Day, End_Year, &
                     End_Month, End_Day)
      spinup = SimTime(Start_Year-nSpinup, Spinup_Month, Spinup_Day, &
                     Start_Year, Start_Month, Start_Day)
      call InitializeSimulation()
      call ModelRun(lakeId, time, spinup, error)
      call ArchiveModelOutput(lakeId, time)
      call FinalizeSimulation()
   end subroutine

end module
