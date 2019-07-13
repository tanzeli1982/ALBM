module sim_coupler_mod
!---------------------------------------------------------------------------------
! Purpose: Govern the simulation process and manage different modules
!          Currently, the modules work in a single CPU. But in the future,
!          different modules will run with the least interactions by MPI technique.
!
!---------------------------------------------------------------------------------
   use shr_ctrl_mod, e8 => SHR_CTRL_E8
   use shr_typedef_mod,       only : SimTime
   use math_utilities_mod,    only : RungeKutta4, adaptive_mode, fixed_mode
   use read_data_mod
   use data_buffer_mod
   use boundary_mod
   use radiation_mod
   use thermal_mod
   use soil_thermal_mod
   use bubble_mod
   use diagenesis_mod
   use carbon_cycle_mod

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: manage the simulation of different modules
   !
   !------------------------------------------------------------------------------
   subroutine ModelRun(lakeId, time, spinup, error)
      implicit none
      integer, intent(in) :: lakeId
      type(SimTime), intent(in) :: time
      type(SimTime), intent(in) :: spinup
      integer, intent(out) :: error

      ! run spinup at first
      call ModuleCoupler(spinup, .True., error)

      ! Run simulation during the interested period
      if (error==0) then
         call InitializeModelOutputs()
         call ConstructActCarbonPool()
         call ModuleCoupler(time, .False., error)
      end if

      if (error==1) then
         call SetNullModelOutputs()
      end if
   end subroutine

   subroutine ModuleCoupler(time, isspinup, error)
      implicit none
      type(SimTime), intent(in) :: time
      logical, intent(in) :: isspinup
      integer, intent(out) :: error
      real(r8) :: curstep, nextstep
      real(r8) :: curstep1, nextstep1
      real(r8) :: curstep2, nextstep2
      real(r8) :: curstep3, nextstep3
      real(r8) :: t, tf
      integer(i8) :: simhour, hindx
      integer :: simday, ncount
      integer :: year, month, day
      logical :: isHourNode

      ! simulation time length
      simday = CalcRunningDays(time)
      simhour = 24 * simday
      tf = 3.6d+3 * DBLE(simhour)

      error = 0                     ! error flag (/=0, error)
      t = 0.0_r8                    ! the timer of simulation
      hindx = 0                     ! simulation output index
      ncount = 0                    ! the count of consecutive small steps
      curstep = 50.0_r8             ! time step of simulation (sec)
      nextstep = MAX_OF_STEP        ! time step in the next cycle
      isHourNode = .False.          ! hourly node flag

      do while(t<tf .and. error==0)
         if(t>=3.6d+3*DBLE(hindx) .and. hindx<simhour) then
            if (DEBUG) then
               if (IsSpinup) then
                  print "(A, I0, A, I0)", "Lake ", lake_info%id, &
                        ": start spinup step ", hindx
               else
                  print "(A, I0, A, I0)", "Lake ", lake_info%id, &
                        ": start formal step ", hindx
               end if
            end if
            isHourNode = .True.
            hindx = hindx + 1
         end if

         if (isHourNode) then
            call GetBoundaryConditions(time, hindx, isspinup)
            call UpdateCatchmentThermalRegime(time, hindx)
            call GetSolarConditions(hindx)
            call CacheHourlyResults(hindx, isspinup)
         end if

         if (Thermal_Module) then
            curstep1 = curstep
            call ThermalModuleSetup()
            call SedThermalModuleSetup()
            call RungeKutta4(HeatEquation, mem_tw, adaptive_mode, Ttol, &
                             curstep1, nextstep1, m_waterTemp, m_tmpWaterTemp)
            call RungeKutta4(SedHeatEquation, mem_ts, fixed_mode, Ttol, &
                             curstep1, curstep1, m_sedTemp, m_tmpSedTemp)
            curstep = min(curstep,curstep1)
            nextstep = min(nextstep,nextstep1)
         end if
         if (Diagenesis_Module) then
            curstep2 = curstep
            call DiagenesisModuleSetup()
            call RungeKutta4(DiagenesisEquation, mem_ch4, adaptive_mode, SStol, &
                             curstep2, nextstep2, m_sedSubCon, m_tmpSedSubCon)
            curstep = min(curstep,curstep2)
            nextstep = min(nextstep,nextstep2)
         end if
         if (Carbon_Module) then
            curstep3 = curstep
            call CarbonModuleSetup(isHourNode)
            call RungeKutta4(CarbonCycleEquation, mem_sub, adaptive_mode, WStol, &
                             curstep3, nextstep3, m_waterSubCon, m_tmpWaterSubCon)
            call RungeKutta4(ParticulateEquation, mem_poc, fixed_mode, WPtol, &
                             curstep3, curstep3, m_waterPOC, m_tmpWaterPOC)
            curstep = min(curstep,curstep3)
            nextstep = min(nextstep,nextstep3)
         end if
         if (Bubble_Module) then
            call BubbleModuleSetup(isHourNode)
            call BubbleDynamics(isHourNode)
         end if
         if (curstep1>curstep+e8 .and. Thermal_Module) then
            call RungeKutta4(HeatEquation, mem_tw, fixed_mode, Ttol, curstep, &
                             curstep, m_waterTemp, m_tmpWaterTemp)
            call RungeKutta4(SedHeatEquation, mem_ts, fixed_mode, Ttol, &
                             curstep, curstep, m_sedTemp, m_tmpSedTemp)
         end if
         if (curstep2>curstep+e8 .and. Diagenesis_Module) then
            call RungeKutta4(DiagenesisEquation, mem_ch4, fixed_mode, SStol, &
                             curstep, curstep, m_sedSubCon, m_tmpSedSubCon)
         end if
         if (curstep3>curstep+e8 .and. Carbon_Module) then
            call RungeKutta4(CarbonCycleEquation, mem_sub, fixed_mode, WStol, &
                             curstep, curstep, m_waterSubCon, m_tmpWaterSubCon)
            call RungeKutta4(ParticulateEquation, mem_poc, fixed_mode, WPtol, &
                             curstep, curstep, m_waterPOC, m_tmpWaterPOC)
         end if
         if (Thermal_Module) then
            m_waterTemp = m_tmpWaterTemp
            m_sedTemp = m_tmpSedTemp
            call SedThermalModuleCallback()
            call ThermalModuleCallback(curstep)
         end if
         if (Diagenesis_Module) then
            m_sedSubCon = m_tmpSedSubCon
            call DiagenesisModuleCallback(curstep)
         end if
         if (Carbon_Module) then
            m_waterSubCon = m_tmpWaterSubCon
            m_waterPOC = m_tmpWaterPOC
            call CarbonModuleCallback(isHourNode, hindx, curstep)
         end if
         if (Bubble_Module) then
            call BubbleModuleCallback(curstep)
         end if

         isHourNode = .False.
         if (curstep<0.1) then
            ncount = ncount + 1
            if (ncount>100) then
               error = 1
               print "(A, I0, A, I0)", 'Lake ', lake_info%id, &
                     ': run diverges at step ', hindx
            end if
            nextstep = 50.0_r8
         else 
            ncount = 0   ! restart the counter
         end if
         t = t + curstep
         curstep = nextstep
         nextstep = MAX_OF_STEP
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: initialization and load restart environment
   !
   !------------------------------------------------------------------------------
   subroutine InitializeSimulation()
      implicit none

      call InitializeDataBuffer()
      call InitializeSmartsModule()
      call InitializeThermalModule()
      call InitializeSedThermalModule()
      call InitializeDiagenesisModule()
      call InitializeCarbonModule()
      call InitializeBubbleModule()
   end subroutine

   subroutine InitializeModelOutputs()
      implicit none

      m_timeHist = 0
      m_tempwHist = 0.0_r4
      m_snowHist = 0.0_r4
      m_iceHist = 0.0_r4
      m_shHist = 0.0_r4
      m_lhHist = 0.0_r4
      m_fmmHist = 0.0_r4
      m_lwHist = 0.0_r4
      m_hnetHist = 0.0_r4
      m_fsedHist = 0.0_r4
      m_swdwHist = 0.0_r4
      m_swupHist = 0.0_r4
      m_fturbHist = 0.0_r4
   end subroutine

   subroutine SetNullModelOutputs()
      implicit none

      m_timeHist = -9999
      m_tempwHist = -9999.0_r4
      m_snowHist = -9999.0_r4
      m_iceHist = -9999.0_r4
      m_shHist = -9999.0_r4
      m_lhHist = -9999.0_r4
      m_fmmHist = -9999.0_r4
      m_lwHist = -9999.0_r4
      m_hnetHist = -9999.0_r4
      m_fsedHist = -9999.0_r4
      m_swdwHist = -9999.0_r4 
      m_swupHist = -9999.0_r4
      m_fturbHist = -9999.0_r4
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: cache results for storage
   !
   !------------------------------------------------------------------------------
   subroutine CacheHourlyResults(hindx, isspinup)
      implicit none
      integer(i8), intent(in)  :: hindx
      logical, intent(in) :: isspinup
      real(r8) :: sh, lh, fmm, lw, hnet
      real(r8) :: fsed, turbdiff(NWLAYER+1)

      if (isspinup) then
         return
      end if

      if (Thermal_Module) then
         m_tempwHist(:,hindx) = m_waterTemp
         m_snowHist(hindx) = m_Hsnow
         m_iceHist(hindx) = m_Hice + m_Hgrayice
         call GetBoundaryOutputs(sh, lh, fmm, lw, hnet, fsed, turbdiff)
         m_shHist(hindx) = sh
         m_lhHist(hindx) = lh
         m_fmmHist(hindx) = fmm
         m_lwHist(hindx) = lw
         m_hnetHist(hindx) = hnet
         m_fsedHist(hindx) = fsed
         m_swdwHist(hindx) = m_surfData%sw_sim
         m_swupHist(hindx) = m_surfData%sw_sim - m_surfData%srd
         m_fturbHist(:,hindx) = turbdiff
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Archive model outputs
   !
   !------------------------------------------------------------------------------
   subroutine ArchiveModelOutput(lakeId, time)
      implicit none
      integer, intent(in) :: lakeId
      type(SimTime), intent(in) :: time

      call WriteData(lakeId, 'zw', m_Zw)
      if (Thermal_Module) then
         call WriteData(lakeId, time, 'watertemp', m_tempwHist)
         call WriteData(lakeId, time, 'snowthick', m_snowHist)
         call WriteData(lakeId, time, 'icethick', m_iceHist)
         call WriteData(lakeId, time, 'sensheatf', m_shHist)
         call WriteData(lakeId, time, 'latentheatf', m_lhHist)
         call WriteData(lakeId, time, 'momf', m_fmmHist)
         call WriteData(lakeId, time, 'lwup', m_lwHist)
         call WriteData(lakeId, time, 'lakeheatf', m_hnetHist)
         call WriteData(lakeId, time, 'sedheatf', m_fsedHist)
         call WriteData(lakeId, time, 'swup', m_swupHist)
         call WriteData(lakeId, time, 'swdw', m_swdwHist)
         call WriteData(lakeId, time, 'turbdiffheat', m_fturbHist)
      end if
   end subroutine

   subroutine FinalizeSimulation()
      implicit none

      call DestructThermalModule()
      call DestructSedThermalModule()
      call DestructDiagenesisModule()
      call DestructCarbonModule()
      call DestructBubbleModule()
      call DestructSmartsModule()
      call DestructDataBuffer()
   end subroutine

end module sim_coupler_mod
