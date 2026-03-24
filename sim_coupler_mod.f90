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
   subroutine ModelRun(rid, time, spinup, otime, error)
      implicit none
      integer, intent(in) :: rid  ! lake or sample id
      type(SimTime), intent(in) :: time
      type(SimTime), intent(in) :: spinup
      type(SimTime), intent(in) :: otime
      integer, intent(out) :: error

      ! run spinup at first
      call ModuleCoupler(rid, spinup, otime, .True., error)

      ! Run simulation during the interested period
      if (error==0) then
         call ConstructOldCarbonPool()
         call ConstructActCarbonPool()
         call ModuleCoupler(rid, time, otime, .False., error)
      end if

      if (error==1) then
         call SetNullModelOutputs()
      end if
   end subroutine

   subroutine ModuleCoupler(rid, time, otime, isspinup, error)
      implicit none
      integer, intent(in) :: rid ! lake or sample id
      type(SimTime), intent(in) :: time
      type(SimTime), intent(in) :: otime
      logical, intent(in) :: isspinup
      integer, intent(out) :: error
      type(SimTime) :: otime0, otime1
      real(r8) :: STtol_col(NSCOL)
      real(r8) :: curstep, nextstep
      real(r8) :: curstep1, nextstep1
      real(r8) :: curstep2, nextstep2
      real(r8) :: curstep3, nextstep3
      real(r8) :: curstep4, nextstep4
      real(r8) :: t, t_old, tf
      integer(i8) :: simhour, hindx, hindx_utc
      integer(i8) :: ohindx0, ohindx1
      integer :: simday, ncount, ndayout
      integer :: year, month, day
      logical :: isHourNode, isSubHourNode

      ! simulation time length
      simday = CalcRunningDays(time, Use_Leap)
      simhour = 24 * simday
      tf = 3.6d+3 * DBLE(simhour)

      ! model output time window
      otime0 = SimTime(time%year0, time%month0, time%day0, &
                       otime%year0, otime%month0, otime%day0)
      otime1 = SimTime(time%year0, time%month0, time%day0, &
                       otime%year1, otime%month1, otime%day1)
      ndayout = CalcRunningDays(otime0, Use_Leap)
      ohindx0 = 24 * ndayout
      ndayout = CalcRunningDays(otime1, Use_Leap)
      ohindx1 = 24 * ndayout


      error = 0                     ! error flag (/=0, error)
      t = 0.0_r8                    ! the timer of simulation
      t_old = 0.0_r8                ! the old time point of sub-hourly update
      hindx = 0                     ! simulation output index
      ncount = 0                    ! the count of consecutive small steps
      curstep = 50.0_r8             ! time step of simulation (sec)
      curstep1 = 50.0_r8
      curstep2 = 50.0_r8
      curstep3 = 50.0_r8
      curstep4 = 50.0_r8
      nextstep = MAX_OF_STEP        ! time step in the next cycle
      isHourNode = .True.           ! hourly node flag
      isSubHourNode = .True.        ! sub-hourly node flag

      STtol_col = STtol
      hindx_utc = GetUTCHourIndex(hindx, lake_info%longitude)

      do while(t<tf .and. error==0)
         if(t>=3.6d+3*DBLE(hindx) .and. hindx<simhour) then
            if (DEBUG) then
               if (IsSpinup) then
                  print "(A, I0, A, I0)", "Run ", rid, &
                        ": start spinup hour step ", hindx
               else
                  print "(A, I0, A, I0)", "Run ", rid, &
                        ": start formal hour step ", hindx
               end if
            end if
            isHourNode = .True.
            hindx = hindx + 1
            hindx_utc = GetUTCHourIndex(hindx, lake_info%longitude)
         end if

         if (t-t_old>=TSTEP_SUB .or. isHourNode) then
            t_old = t
            isSubHourNode = .True.
         end if

         if (isHourNode) then
            if (trim(forcing_time)=='UTC') then
               call GetBoundaryConditions(time, hindx_utc, isspinup)
            else
               call GetBoundaryConditions(time, hindx, isspinup)
            end if
            call GetSolarConditions(hindx)
         end if

         if (Thermal_Module) then
            curstep1 = curstep
            curstep4 = curstep
            call ThermalModuleSetup()
            call SedThermalModuleSetup()
            call RungeKutta4(HeatEquation, mem_tw, adaptive_mode, Ttol, &
                             curstep1, nextstep1, m_waterTemp, m_tmpWaterTemp)
            call RungeKutta4(SedHeatEquation, mem_ts, adaptive_mode, STtol_col, &
                             curstep4, nextstep4, m_sedTemp, m_tmpSedTemp)
            curstep = min(curstep,min(curstep1,curstep4))
            nextstep = min(nextstep,min(nextstep1,nextstep4))
         end if
         if (Diagenesis_Module) then
            curstep2 = curstep
            call DiagenesisModuleSetup(isSubHourNode)
            call RungeKutta4(DiagenesisEquation, mem_ch4, adaptive_mode, SStol, &
                             curstep2, nextstep2, m_sedSubCon, m_tmpSedSubCon)
            curstep = min(curstep,curstep2)
            nextstep = min(nextstep,nextstep2)
         end if
         if (Carbon_Module) then
            curstep3 = curstep
            call CarbonModuleSetup(isHourNode, isSubHourNode, hindx)
            call RungeKutta4(CarbonCycleEquation, mem_sub, adaptive_mode, WStol, &
                             curstep3, nextstep3, m_waterSubCon, m_tmpWaterSubCon)
            curstep = min(curstep,curstep3)
            nextstep = min(nextstep,nextstep3)
         end if
         if (Bubble_Module) then
            call BubbleModuleSetup(isHourNode)
            call BubbleDynamics(isHourNode)
         end if
         if (Thermal_Module) then
            if (curstep1>curstep+e8) then
               call RungeKutta4(HeatEquation, mem_tw, fixed_mode, Ttol, &
                     curstep, curstep, m_waterTemp, m_tmpWaterTemp)
            end if
            if (curstep4>curstep+e8) then
               call RungeKutta4(SedHeatEquation, mem_ts, fixed_mode, &
                     STtol_col, curstep, curstep, m_sedTemp, m_tmpSedTemp)
            end if
            m_waterTemp = m_tmpWaterTemp
            m_sedTemp = m_tmpSedTemp
         end if
         if (Diagenesis_Module) then
            if (curstep2>curstep+e8) then
               call RungeKutta4(DiagenesisEquation, mem_ch4, fixed_mode, &
                     SStol, curstep, curstep, m_sedSubCon, m_tmpSedSubCon)
            end if
            m_sedSubCon = m_tmpSedSubCon
         end if
         if (Carbon_Module) then
            if (curstep3>curstep+e8) then
               call RungeKutta4(CarbonCycleEquation, mem_sub, fixed_mode, &
                     WStol, curstep, curstep, m_waterSubCon, m_tmpWaterSubCon)
            end if
            m_waterSubCon = m_tmpWaterSubCon
         end if
         if (Thermal_Module) then
            call SedThermalModuleCallback()
            call ThermalModuleCallback(curstep)
         end if
         if (Carbon_Module) then
            call PhytoplanktonDynamics(curstep)
            call BedVegetationDynamics(curstep)
            call CarbonModuleCallback(curstep)
         end if
         if (Diagenesis_Module) then
            call DiagenesisModuleCallback(curstep)
         end if
         if (Bubble_Module) then
            call BubbleModuleCallback(curstep)
         end if

         if (isHourNode .and. (.NOT. isspinup)) then
            if (hindx>ohindx0 .and. hindx<=ohindx1) then
               call CacheHourlyResults(hindx-ohindx0)
            end if
         end if

         isHourNode = .False.
         isSubHourNode = .False.
         if (curstep<0.1) then
            ncount = ncount + 1
            if (ncount>100) then
               error = 1
               print "(A, I0, A, I0)", 'Run ', rid, ': diverges at step ', hindx
            end if
            nextstep = 50.0_r8
         else if (curstep>=0.1) then
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
      call InitializeModelOutputs()
   end subroutine

   subroutine InitializeModelOutputs()
      implicit none

      m_tempwHist = 0.0_r4
      m_tempsHist = 0.0_r4
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
      m_dfch4airHist = 0.0_r4
      m_ebch4airHist = 0.0_r4
      m_dfch4sedHist = 0.0_r4
      m_ebch4sedHist = 0.0_r4
      m_totGPPHist = 0.0_r4
      m_totNPPHist = 0.0_r4
      m_totpch4sedHist = 0.0_r4
      m_pch4watHist = 0.0_r4
      m_och4watHist = 0.0_r4
      m_depAtCHist = 0.0_r4
      m_ch4concHist = 0.0_r4
      m_o2concHist = 0.0_r4
      m_n2concHist = 0.0_r4
      m_co2concHist = 0.0_r4
      m_srpconcHist = 0.0_r4
      m_icebch4Hist = 0.0_r4
      m_bvegHist = 0.0_r4
      m_biomassHist = 0.0_r4
      m_chlaHist = 0.0_r4
      m_pascarbHist = 0.0_r4
      m_actcarbHist = 0.0_r4
      m_oldcarbHist = 0.0_r4
      m_sedEHLHist = 0.0_r4
   end subroutine

   subroutine SetNullModelOutputs()
      implicit none

      m_tempwHist = -9999.0_r4
      m_tempsHist = -9999.0_r4
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
      m_dfch4airHist = -9999.0_r4
      m_ebch4airHist = -9999.0_r4
      m_dfch4sedHist = -9999.0_r4
      m_ebch4sedHist = -9999.0_r4
      m_totGPPHist = -9999.0_r4
      m_totNPPHist = -9999.0_r4
      m_totpch4sedHist = -9999.0_r4
      m_pch4watHist = -9999.0_r4
      m_och4watHist = -9999.0_r4
      m_depAtCHist = -9999.0_r4
      m_ch4concHist = -9999.0_r4
      m_o2concHist = -9999.0_r4
      m_n2concHist = -9999.0_r4
      m_co2concHist = -9999.0_r4
      m_srpconcHist = -9999.0_r4
      m_icebch4Hist = -9999.0_r4
      m_bvegHist = -9999.0_r4
      m_biomassHist = -9999.0_r4
      m_chlaHist = -9999.0_r4
      m_pascarbHist = -9999.0_r4
      m_actcarbHist = -9999.0_r4
      m_oldcarbHist = -9999.0_r4
      m_sedEHLHist = -9999.0_r4
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: cache results for storage
   !
   !------------------------------------------------------------------------------
   subroutine CacheHourlyResults(hindx)
      implicit none
      integer(i8), intent(in)  :: hindx
      real(r8) :: sh, lh, fmm, lw, hnet
      real(r8) :: fsed, turbdiff(NWLAYER+1)
      real(r8) :: gpp, npp, pch4(NWLAYER+1)
      real(r8) :: och4(NWLAYER+1), sedpch4(NSCOL)

      if (Thermal_Module) then
         call GetBoundaryOutputs(sh, lh, fmm, lw, hnet, fsed, turbdiff)
         m_tempwHist(:,hindx) = REAL(m_waterTemp)
         m_tempsHist(:,:,hindx) = REAL(m_sedTemp)
         m_snowHist(hindx) = REAL(m_Hsnow)
         m_iceHist(hindx) = REAL(m_Hice + m_Hgrayice)
         m_shHist(hindx) = REAL(sh)
         m_lhHist(hindx) = REAL(lh)
         m_fmmHist(hindx) = REAL(fmm)
         m_lwHist(hindx) = REAL(lw)
         m_hnetHist(hindx) = REAL(hnet)
         m_fsedHist(hindx) = REAL(fsed)
         m_swdwHist(hindx) = REAL(m_surfData%sw_sim)
         m_swupHist(hindx) = REAL(m_surfData%sw_sim - m_surfData%srd)
         m_fturbHist(:,hindx) = REAL(turbdiff)
      end if
      if (Carbon_Module) then
         call GetProductionRates(gpp, npp)
         call GetRespirationRates(pch4, och4)
         m_dfch4airHist(hindx) = REAL(m_topdflux(Wch4))
         m_totGPPHist(hindx) = REAL(gpp)
         m_totNPPHist(hindx) = REAL(npp)
         m_pch4watHist(:,hindx) = REAL(pch4)
         m_och4watHist(:,hindx) = REAL(och4)
         m_ch4concHist(:,hindx) = REAL(m_waterSubCon(Wch4,:))
         m_o2concHist(:,hindx) = REAL(m_waterSubCon(Wo2,:))
         m_n2concHist(:,hindx) = REAL(m_waterSubCon(Wn2,:))
         m_co2concHist(:,hindx) = REAL(m_waterSubCon(Wco2,:))
         m_srpconcHist(:,hindx) = REAL(m_waterSubCon(Wsrp,:))
         m_chlaHist(:,hindx) = REAL(m_chla)
         m_biomassHist(:,:,hindx) = REAL(m_waterPOC)
         m_bvegHist(:,hindx) = REAL(cCPerDW * m_bedVegDW)
         m_depAtCHist(:,hindx) = REAL(m_aqDepAtCarb + m_trDepAtCarb + &
            m_vegDepAtCarb)
      end if
      if (Diagenesis_Module) then
         call GetSedCLossRate(sedpch4) 
         m_dfch4sedHist(:,hindx) = REAL(m_btmdflux(Wch4,:))
         m_ebch4sedHist(:,hindx) = REAL(m_btmbflux(Wch4,:))
         m_pascarbHist(:,:,hindx) = REAL(m_unfrzCarbPool(pasC,:,:))
         m_actcarbHist(:,:,hindx) = REAL(m_unfrzCarbPool(actC,:,:))
         m_oldcarbHist(:,:,hindx) = REAL(m_unfrzCarbPool(oldC,:,:))
         m_totpch4sedHist(:,hindx) = REAL(sedpch4)
         M_sedEHLHist(:,:,hindx) = REAL(m_sedEHL)
      end if
      if (Bubble_Module) then
         m_ebch4airHist(hindx) = REAL(m_topbflux(Wch4))
         m_icebch4Hist(hindx) = REAL(m_iceBubblePool(Wch4))
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

      call WriteData(lakeId, time, 'zw', archive_tstep, m_Zw)
      call WriteData(lakeId, time, 'zs', archive_tstep, m_Zs)
      call WriteData(lakeId, time, 'Az', archive_tstep, m_Az)
      call WriteData(lakeId, time, 'colindx', archive_tstep, m_soilColInd)
      if (Thermal_Module) then
         call WriteData(lakeId, time, 'watertemp', archive_tstep, m_tempwHist)
         call WriteData(lakeId, time, 'sedtemp', archive_tstep, m_tempsHist)
         call WriteData(lakeId, time, 'snowthick', archive_tstep, m_snowHist)
         call WriteData(lakeId, time, 'icethick', archive_tstep, m_iceHist)
         call WriteData(lakeId, time, 'sensheatf', archive_tstep, m_shHist)
         call WriteData(lakeId, time, 'latentheatf', archive_tstep, m_lhHist)
         call WriteData(lakeId, time, 'momf', archive_tstep, m_fmmHist)
         call WriteData(lakeId, time, 'lwup', archive_tstep, m_lwHist)
         call WriteData(lakeId, time, 'lakeheatf', archive_tstep, m_hnetHist)
         call WriteData(lakeId, time, 'sedheatf', archive_tstep, m_fsedHist)
         call WriteData(lakeId, time, 'swup', archive_tstep, m_swupHist)
         !call WriteData(lakeId, time, 'swdw', archive_tstep, m_swdwHist)
         !call WriteData(lakeId, time, 'turbdiffheat', archive_tstep, m_fturbHist)
      end if
      if (Carbon_Module) then
         call WriteData(lakeId, time, 'ch4df', archive_tstep, m_dfch4airHist)
         call WriteData(lakeId, time, 'gpp', archive_tstep, m_totGPPHist)
         call WriteData(lakeId, time, 'npp', archive_tstep, m_totNPPHist)
         call WriteData(lakeId, time, 'och4prod', archive_tstep, m_pch4watHist)
         call WriteData(lakeId, time, 'ch4oxid', archive_tstep, m_och4watHist)
         call WriteData(lakeId, time, 'ch4conc', archive_tstep, m_ch4concHist)
         call WriteData(lakeId, time, 'o2conc', archive_tstep, m_o2concHist)
         call WriteData(lakeId, time, 'n2conc', archive_tstep, m_n2concHist)
         call WriteData(lakeId, time, 'co2conc', archive_tstep, m_co2concHist)
         call WriteData(lakeId, time, 'srp', archive_tstep, m_srpconcHist)
         call WriteData(lakeId, time, 'chla', archive_tstep, m_chlaHist)
         call WriteData(lakeId, time, 'cdep', archive_tstep, m_depAtCHist)
         call WriteData(lakeId, time, 'biomass', archive_tstep, m_biomassHist)
         call WriteData(lakeId, time, 'bveg', archive_tstep, m_bvegHist)
      end if
      if (Diagenesis_Module) then
         call WriteData(lakeId, time, 'sedch4df', archive_tstep, m_dfch4sedHist)
         call WriteData(lakeId, time, 'sedch4eb', archive_tstep, m_ebch4sedHist)
         call WriteData(lakeId, time, 'ch4prod', archive_tstep, m_totpch4sedHist)
         call WriteData(lakeId, time, 'pascarb', archive_tstep, m_pascarbHist)
         call WriteData(lakeId, time, 'actcarb', archive_tstep, m_actcarbHist)
         call WriteData(lakeId, time, 'oldcarb', archive_tstep, m_oldcarbHist)
         call WriteData(lakeId, time, 'sedrdx', archive_tstep, m_sedEHLHist) 
      end if
      if (Bubble_Module) then
         call WriteData(lakeId, time, 'ch4eb', archive_tstep, m_ebch4airHist)
         call WriteData(lakeId, time, 'icebch4', archive_tstep, m_icebch4Hist)
      end if
   end subroutine

   subroutine ArchiveSensitivityOutput(sampleId, time)
      implicit none
      integer, intent(in) :: sampleId
      type(SimTime), intent(in) :: time

      call WriteData(sampleId, time, 'zw', archive_tstep, m_Zw)
      call WriteData(sampleId, time, 'zs', archive_tstep, m_Zs)
      call WriteData(sampleId, time, 'Az', archive_tstep, m_Az)
      call WriteData(sampleId, time, 'colindx', archive_tstep, m_soilColInd)
      if (Thermal_Module) then
         call WriteData(sampleId, time, 'watertemp', archive_tstep, m_tempwHist)
         call WriteData(sampleId, time, 'snowthick', archive_tstep, m_snowHist)
         call WriteData(sampleId, time, 'icethick', archive_tstep, m_iceHist)
      end if
      if (Carbon_Module) then
         call WriteData(sampleId, time, 'ch4df', archive_tstep, m_dfch4airHist)
         call WriteData(sampleId, time, 'gpp', archive_tstep, m_totGPPHist)
         call WriteData(sampleId, time, 'npp', archive_tstep, m_totNPPHist)
         call WriteData(sampleId, time, 'och4prod', archive_tstep, m_pch4watHist)
         call WriteData(sampleId, time, 'ch4oxid', archive_tstep, m_och4watHist)
         call WriteData(sampleId, time, 'ch4conc', archive_tstep, m_ch4concHist)
         call WriteData(sampleId, time, 'o2conc', archive_tstep, m_o2concHist)
         call WriteData(sampleId, time, 'chla', archive_tstep, m_chlaHist)
         call WriteData(sampleId, time, 'cdep', archive_tstep, m_depAtCHist)
         call WriteData(sampleId, time, 'biomass', archive_tstep, m_biomassHist)
         call WriteData(sampleId, time, 'bveg', archive_tstep, m_bvegHist)
      end if
      if (Diagenesis_Module) then
         call WriteData(sampleId, time, 'sedch4df', archive_tstep, m_dfch4sedHist)
         call WriteData(sampleId, time, 'sedch4eb', archive_tstep, m_ebch4sedHist)
         call WriteData(sampleId, time, 'ch4prod', archive_tstep, m_totpch4sedHist)
         call WriteData(sampleId, time, 'actcarb', archive_tstep, m_actcarbHist)
         call WriteData(sampleId, time, 'oldcarb', archive_tstep, m_oldcarbHist)
      end if
      if (Bubble_Module) then
         call WriteData(sampleId, time, 'ch4eb', archive_tstep, m_ebch4airHist)
         call WriteData(sampleId, time, 'icebch4', archive_tstep, m_icebch4Hist)
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
