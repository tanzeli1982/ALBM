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
      integer :: lake_next_range(2)
      integer :: lakeId, err, minid, maxid
      integer :: nlake, istep, ii, itmp
      integer :: ntotlake, jj

      minid = minval(lake_range)
      maxid = maxval(lake_range)
      nlake = maxid - minid + 1

      print "(A, I0, A, I0)", 'Run lakes from ', minid, ' to ', maxid

      call DoSimulationWarmup(minid, ntotlake)

      do istep = 1, nlake, 1 
         lakeId = minid + istep - 1 
         call RegularSimulation(lakeId, err)
         print "(A, I0, A, I0)", "Lake ", lakeId, ", error ", err
      end do

   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: create simulation archive files 
   !
   !------------------------------------------------------------------------------
   subroutine DoSimulationWarmup(minid, nlake)
      implicit none
      integer, intent(in) :: minid
      integer, intent(out) :: nlake
      type(SimTime) :: time
      
      call ReadLakeTotalNumber(nlake)
      if (minid==1 .or. DEBUG) then
         time = SimTime(Start_Year, Start_Month, Start_Day, End_Year, &
                     End_Month, End_Day)
         call CreateOutputFile(NWLAYER+1, 'zw', 'water layer depth', 'm') 
         call CreateOutputFile(time, 'snowthick', 'snow cover thickness', &
                              'm', -9999.0_r4)
         call CreateOutputFile(time, 'icethick', 'ice cover thickness', &
                              'm', -9999.0_r4)
         call CreateOutputFile(time, 'sensheatf', 'upward sensible ' // &
                              'heat flux', 'W m-2', -9999.0_r4)
         call CreateOutputFile(time, 'latentheatf', 'upward latent ' // &
                              'heat flux', 'W m-2', -9999.0_r4)
         call CreateOutputFile(time, 'momf', 'momentum energy flux', &
                              'kg m-1 s-2', -9999.0_r4)
         call CreateOutputFile(time, 'lwup', 'upward longwave radiation', &
                              'W m-2', -9999.0_r4)
         call CreateOutputFile(time, 'lakeheatf', 'downward net heat flux', &
                              'W m-2', -9999.0_r4)
         call CreateOutputFile(time, 'swdw', 'downward shortwave radiation', &
                              'W m-2', -9999.0_r4)
         call CreateOutputFile(time, 'swup', 'upward shortwave radiation', &
                              'W m-2', -9999.0_r4)
         call CreateOutputFile(time, 'sedheatf', 'upward sediment heat ' // &
                              'flux', 'W m-2', -9999.0_r4)
         call CreateOutputFile(time, NWLAYER+1, 'watertemp', 'water ' // &
                              'temperature', 'K', -9999.0_r4)
         call CreateOutputFile(time, NWLAYER+1, 'turbdiffheat', &
                               'Turbulent diffusivity of heat', 'm2 s-1', &
                               -9999.0_r4)
      end if
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
