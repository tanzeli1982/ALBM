module data_buffer_mod
!---------------------------------------------------------------------------------
! Purpose:
!
! This module serves for the buffer of large temporary data in the model
!
!---------------------------------------------------------------------------------
   use shr_kind_mod,          only : r8, r4, i8
   use shr_ctrl_mod, e8 => SHR_CTRL_E8
   use shr_typedef_mod
   use read_data_mod
   use phy_utilities_mod

   implicit none
   public
   ! YYYYMMDD of each record
   integer, allocatable :: m_timeHist(:)
   ! simulated hourly water and sediment temperature profiles (K)
   real(r4), allocatable :: m_tempwHist(:,:)
   ! simulated hourly lake snow and ice thickness
   real(r4), allocatable :: m_snowHist(:), m_iceHist(:)
   ! hourly boundary conditions
   real(r4), allocatable :: m_shHist(:)   ! sensible heat
   real(r4), allocatable :: m_lhHist(:)   ! latent heat
   real(r4), allocatable :: m_fmmHist(:)  ! momentum energy flux
   real(r4), allocatable :: m_lwHist(:)   ! upward longwave radiation
   real(r4), allocatable :: m_hnetHist(:) ! downward net heat flux
   real(r4), allocatable :: m_swdwHist(:) ! downward shortwave radiation
   real(r4), allocatable :: m_swupHist(:) ! upward shortwave radiation
   real(r4), allocatable :: m_fsedHist(:) ! upward sediment heat flux
   real(r4), allocatable :: m_fturbHist(:,:) ! turbulent diffusivity of heat
   ! incumbent water and sediment temperature profile (K)
   real(r8), allocatable :: m_waterTemp(:), m_sedTemp(:)
   ! temporary water and sediment temperature profile (K) 
   real(r8), allocatable :: m_tmpWaterTemp(:), m_tmpSedTemp(:)
   ! water and sediment ice profile (m)
   real(r8), allocatable :: m_waterIce(:), m_sedIce(:)
   ! incumbent and temporary bubble gas amount (umol/m3/mm)
   real(r8), allocatable :: m_bubbleGasCon(:,:,:), m_tmpBubbleGasCon(:,:,:)
   ! incumbent and temporary solvent concentrations in water (umol/m3)
   real(r8), allocatable :: m_waterSubCon(:,:), m_tmpWaterSubCon(:,:)
   ! incumbent and temporary solvent concentrations in sediments (umol/m3)
   real(r8), allocatable :: m_sedSubCon(:,:), m_tmpSedSubCon(:,:)
   ! total heat diffusivity (m2/s)
   real(r8), allocatable :: m_Kt(:)
   real(r8), allocatable :: m_Ks(:)
   ! total momentum diffusivity (m2/s)
   real(r8), allocatable :: m_Kv(:)
   ! water density (kg/m3)
   real(r8), allocatable :: m_wrho(:)
   ! water dynamic viscosity (kg/s/m)
   real(r8), allocatable :: m_dVsc(:)
   ! Solar radiation in water column (W/m2)
   real(r8), allocatable :: m_Iab(:)
   ! incumbent and temporary POC and DOC concentrations in water (umol/m3)
   real(r8), allocatable :: m_waterPOC(:,:)
   real(r8), allocatable :: m_tmpWaterPOC(:,:)
   ! Chl : C ratio profile (mg Chl mmol C-1)
   real(r8), allocatable :: m_rChl2C(:,:)
   ! gas exchange amount from bubble to water (umol/m3/s)
   real(r8), allocatable :: m_gasExchange(:,:)
   ! gas pools for bubbles trapped in the ice (umol) 
   real(r8), allocatable :: m_iceBubblePool(:)
   ! carbon pools for settling phytoplankton (umol) 
   real(r8), allocatable :: m_sinkPOCPool(:)
   ! 14C-enriched and 14C-depleted decomposable carbon pools (umol/m3)
   real(r8), allocatable :: m_frzCarbPool(:,:)
   real(r8), allocatable :: m_unfrzCarbPool(:,:)
   ! downward scalar irradiance (mol photons m-2 s-1 nm-1)
   real(r8), allocatable :: m_fsphot(:,:)
   ! downward irradiance spectrum (nm)
   real(r8), allocatable :: m_wvln(:)
   ! daily air climate records above the lake
   real(r8), allocatable :: m_airTemp(:)
   real(r8), allocatable :: m_airTempMax(:)
   real(r8), allocatable :: m_airTempMin(:)
   real(r8), allocatable :: m_airRH(:)
   real(r8), allocatable :: m_airWind(:)
   real(r8), allocatable :: m_airPr(:)
   real(r8), allocatable :: m_airPrsn(:) 
   real(r8), allocatable :: m_airPs(:)
   real(r8), allocatable :: m_airSWRad(:)
   real(r8), allocatable :: m_airLWRad(:)
   real(r8), allocatable :: m_Qsi(:)
   real(r8), allocatable :: m_tQsi(:)
   real(r8), allocatable :: m_dQsi(:)
   real(r8), allocatable :: m_DOQsi(:)
   real(r8), allocatable :: m_DICQsi(:)
   real(r8), allocatable :: m_DOCQsi(:)
   real(r8), allocatable :: m_POCQsi(:)
   real(r8), allocatable :: m_SRPQsi(:)
   real(r8), allocatable :: m_Qso(:)
   real(r8), allocatable :: m_Qgw(:)
   ! radiation-related conditions
   real(r8), allocatable :: m_aCO2(:)
   real(r8), allocatable :: m_aO3(:)
   real(r8), allocatable :: m_aAOD(:)
   ! water and sediment depth vector (m)
   real(r8), allocatable :: m_Zw(:), m_Zs(:)
   real(r8), allocatable :: m_dZw(:), m_dZs(:)
   ! lake depth-dependent area vector (m2)
   real(r8), allocatable :: m_Az(:)
   real(r8), allocatable :: m_dAz(:)
   ! bubble radius vector (m)
   real(r8), allocatable :: m_Rb0(:)
   ! bubble water interface fluxes (umol/m2/s)
   real(r8), allocatable :: m_btmbflux(:)
   ! autochthonous and allochthonous sedimentation rate (umol/m2/s)
   real(r8) :: m_burialAtCarb, m_burialAlCarb
   ! lake ice, gray ice and snow thickness (m)
   real(r8) :: m_Hice, m_Hsnow, m_Hgrayice 
   ! areal forest fraction
   real(r8) :: m_ftree, m_fwlnd
   ! SOC density (kg/m2)
   real(r8) :: m_soc
   ! initial DOC (mol/m3), SRP (mmol/m3) and algae (mg/m3)
   real(r8) :: m_chla0, m_DOC0, m_DOC1, m_SRP0
   ! bottom molecular diffusivity (m2/s)
   real(r8) :: m_Kbm
   ! mixing layer thickness (m)
   real(r8) :: m_Hmix
   ! surface and bottom boundary layer thickness (m)
   real(r8) :: m_HbLayer(2)
   ! effective heat flux of the open-water mixing layer (W/m2)
   real(r8) :: m_Heff
   ! lake remaining total kinetic energy (J)
   real(r8) :: m_Ekin
   ! user-defined surface forcing data
   type(SurfaceData) :: m_surfData
   ! user-defined radiation parameter set
   type(RadParaData) :: m_radPars
   ! lake free water top index
   integer :: m_lakeWaterTopIndex              
   ! lake water mixing indices
   integer :: m_mixTopIndex
   ! sediment free water top and bottom index
   integer :: m_sedWaterTopIndex, m_sedWaterBtmIndex 
   ! winter flag
   integer :: winter_flag, prewinter_flag

contains
   !------------------------------------------------------------------------------
   !
   ! Data buffer must be destructed before changing its size
   !
   !------------------------------------------------------------------------------
   subroutine InitializeDataBuffer()
      implicit none
      type(SimTime) :: time
      integer :: simday, simmon, simyr
      integer :: ntin, ntout
      real(r8) :: loc180(2)

      loc180 = (/lake_info%longitude, lake_info%latitude/)

      time = SimTime(Start_Year,Start_Month,Start_Day,End_Year, &
                     End_Month,End_Day)
      simday = CalcRunningDays(time)
      simmon = CalcRunningMonths(time)
      simyr = CalcRunningYears(time)

      if (trim(forcing_tstep)=='hour') then
         ntin = 24 * simday
      else if (trim(forcing_tstep)=='day') then
         ntin = simday
      end if
      ntout = 24 * simday

      ! allocate memory for the archive variables
      allocate(m_timeHist(ntout))
      allocate(m_tempwHist(WATER_LAYER+1,ntout))
      allocate(m_snowHist(ntout))
      allocate(m_iceHist(ntout))
      allocate(m_shHist(ntout))
      allocate(m_lhHist(ntout))
      allocate(m_fmmHist(ntout))
      allocate(m_lwHist(ntout))
      allocate(m_hnetHist(ntout))
      allocate(m_swdwHist(ntout))
      allocate(m_swupHist(ntout))
      allocate(m_fsedHist(ntout))
      allocate(m_fturbHist(WATER_LAYER+1,ntout))
      ! allocate memory for the state variables
      allocate(m_waterTemp(WATER_LAYER+1))
      allocate(m_sedTemp(NSLAYER+1))
      allocate(m_tmpWaterTemp(WATER_LAYER+1))
      allocate(m_tmpSedTemp(NSLAYER+1))
      allocate(m_waterIce(WATER_LAYER+1))
      allocate(m_sedIce(NSLAYER+1))
      allocate(m_bubbleGasCon(NGAS,NRLAYER+1,WATER_LAYER+1))
      allocate(m_tmpBubbleGasCon(NGAS,NRLAYER+1,WATER_LAYER+1))
      allocate(m_iceBubblePool(NGAS))
      allocate(m_waterSubCon(NWSUB,WATER_LAYER+1))
      allocate(m_tmpWaterSubCon(NWSUB,WATER_LAYER+1))
      allocate(m_gasExchange(NGAS,WATER_LAYER+1))
      allocate(m_sedSubCon(NSSUB,NSLAYER+1))
      allocate(m_tmpSedSubCon(NSSUB,NSLAYER+1)) 
      allocate(m_Kt(WATER_LAYER+1))
      allocate(m_Kv(WATER_LAYER+1))
      allocate(m_Ks(NSLAYER+1))
      allocate(m_wrho(WATER_LAYER+1))
      allocate(m_dVsc(WATER_LAYER+1))
      allocate(m_Iab(WATER_LAYER+2))
      allocate(m_sinkPOCPool(NPOC))
      allocate(m_waterPOC(NPOC,WATER_LAYER+1))
      allocate(m_tmpWaterPOC(NPOC,WATER_LAYER+1))
      allocate(m_rChl2C(NPOC,WATER_LAYER+1))
      allocate(m_frzCarbPool(NPOOL,NSLAYER+1))
      allocate(m_unfrzCarbPool(NPOOL,NSLAYER+1))
      ! allocate memory for air forcing data
      allocate(m_airTemp(ntin))
      allocate(m_airTempMax(ntin))
      allocate(m_airTempMin(ntin))
      allocate(m_airRH(ntin))
      allocate(m_airWind(ntin))
      allocate(m_airPrsn(ntin))
      allocate(m_airPr(ntin))
      allocate(m_airPs(ntin))
      allocate(m_airSWRad(ntin))
      allocate(m_airLWRad(ntin))
      allocate(m_Qsi(simday))
      allocate(m_tQsi(simday))
      allocate(m_dQsi(simday))
      allocate(m_DOQsi(simday))
      allocate(m_DICQsi(simday))
      allocate(m_DOCQsi(simday))
      allocate(m_POCQsi(simday))
      allocate(m_SRPQsi(simday))
      allocate(m_Qso(simday))
      allocate(m_Qgw(simday))
      ! allocate memory for irradiation data
      allocate(m_aCO2(simyr))
      allocate(m_aO3(simmon))
      allocate(m_aAOD(simmon))
      ! allocate memory for radiation vector
      allocate(m_fsphot(NSPCTM,WATER_LAYER+1))
      allocate(m_wvln(NSPCTM))
      ! allocate depth vector
      allocate(m_Zs(NSLAYER+1))
      allocate(m_Zw(WATER_LAYER+1))
      allocate(m_dZw(WATER_LAYER+1))
      allocate(m_dZs(NSLAYER+1))
      allocate(m_Az(WATER_LAYER+1))
      allocate(m_dAz(WATER_LAYER+1))
      allocate(m_Rb0(NRLAYER+1))
      ! allocate intermediate vector
      allocate(m_btmbflux(NGAS))

      ! initialize depth vectors
      call ConstructDepthVector(lake_info, m_Zw, m_Zs, m_dZw, m_dZs)
      if (len_trim(bthmtry_dir)==0) then
         call BuildLakeBathymetry(lake_info, m_Zw, m_Az, m_dAz)
      else
         call ReadLakeBathymetry(lake_info, m_Zw, m_Az, m_dAz)
      end if
      call ReadStaticData(veg_file, loc180, 'ftree', m_ftree)
      call ReadStaticData(wlnd_file, loc180, 'glwd', m_fwlnd)
      call ReadStaticData(soc_file, loc180, 'soc', Ncinterp, m_soc)
      call ReadStaticData(tref_file, loc180, 't2m', Ncinterp, m_radPars%tref)
      ! initialize formal air forcing data
      if (len_trim(forcing_dir)==0) then
         call Read2DTSData(tas_file, time, loc180, 'tas', m_airTemp)
         call Read2DTSData(tasmax_file, time, loc180, 'tasmax', m_airTempMax)
         call Read2DTSData(tasmin_file, time, loc180, 'tasmin', m_airTempMin)
         call Read2DTSData(hurs_file, time, loc180, 'hurs', m_airRH)
         call Read2DTSData(pr_file, time, loc180, 'pr', m_airPr)
         call Read2DTSData(prsn_file, time, loc180, 'prsn', m_airPrsn)
         call Read2DTSData(ps_file, time, loc180, 'ps', m_airPs)
         call Read2DTSData(wind_file, time, loc180, 'sfcWind', m_airWind)
         call Read2DTSData(rsds_file, time, loc180, 'rsds', m_airSWRad)
         call Read2DTSData(rlds_file, time, loc180, 'rlds', m_airLWRad)
      else
         call ReadSiteTSData(lake_info, time, forcing_tstep, 'tas',  m_airTemp)
         call ReadSiteTSData(lake_info, time, forcing_tstep, 'tasmax', m_airTempMax)
         call ReadSiteTSData(lake_info, time, forcing_tstep, 'tasmin',  m_airTempMin)
         call ReadSiteTSData(lake_info, time, forcing_tstep, 'hurs',  m_airRH)
         call ReadSiteTSData(lake_info, time, forcing_tstep, 'pr',  m_airPr)
         call ReadSiteTSData(lake_info, time, forcing_tstep, 'prsn',  m_airPrsn)
         call ReadSiteTSData(lake_info, time, forcing_tstep, 'ps',  m_airPs)
         call ReadSiteTSData(lake_info, time, forcing_tstep, 'sfcWind',  m_airWind)
         call ReadSiteTSData(lake_info, time, forcing_tstep, 'rsds',  m_airSWRad)
         call ReadSiteTSData(lake_info, time, forcing_tstep, 'rlds',  m_airLWRad)
      end if
      ! initialize long-term forcing data
      call ReadGlobalTSData(co2_file, time, 'co2_rcp26', m_aCO2) 
      call Read2DTSData(o3_file, time, loc180, 'tro3', m_aO3)
      call Read2DTSData(aod_file, time, loc180, 'AOD_550', m_aAOD)
      ! no hydrology and chemistry input data
      m_Qsi = 0.0_r8
      m_tQsi = T0 + 4.0
      m_dQsi = 1d3
      m_DOQsi = 0.0_r8
      m_DICQsi = 0.0_r8
      m_DOCQsi = 0.0_r8
      m_POCQsi = 0.0_r8
      m_SRPQsi = 0.0_r8
      m_Qso = 0.0_r8
      m_Qgw = 0.0_r8
      ! rescale tree cover and wetland fraction
      call RescaleTreeCoverFraction(m_ftree)
      call RescaleWetlandFraction(m_fwlnd)
      m_fwlnd = 1.0  ! assume all wetlands
      m_ftree = 0.0  ! assume no trees
      ! units conversion
      m_airPr = 1.0d-3 * m_airPr             ! convert to m/s (water)
      m_airPrsn = 1.0d-3 * m_airPrsn         ! convert to m/s (water)
      m_DOQsi = 1d+3 * m_DOQsi               ! convert to umol/m3
      m_DICQsi = 1d+3 * m_DICQsi             ! convert to umol/m3
      m_DOCQsi = 1d+3 * m_DOCQsi             ! convert to umol/m3
      m_POCQsi = 1d+3 * m_POCQsi             ! convert to umol/m3
      m_SRPQsi = 1d+3 * m_SRPQsi             ! convert to umol/m3
      m_aO3 = 1.0d-3 * m_aO3                 ! convert to 1000 DU
   end subroutine

   subroutine DestructDataBuffer()
      implicit none

      ! destroy the memory for the archive variables
      deallocate(m_timeHist)
      deallocate(m_tempwHist)
      deallocate(m_snowHist)
      deallocate(m_iceHist)
      deallocate(m_shHist)
      deallocate(m_lhHist)
      deallocate(m_fmmHist)
      deallocate(m_lwHist)
      deallocate(m_hnetHist)
      deallocate(m_swdwHist)
      deallocate(m_swupHist)
      deallocate(m_fsedHist)
      deallocate(m_fturbHist)
      ! destroy the memory for the state variables
      deallocate(m_waterTemp)
      deallocate(m_sedTemp)
      deallocate(m_tmpWaterTemp)
      deallocate(m_tmpSedTemp)
      deallocate(m_bubbleGasCon)
      deallocate(m_tmpBubbleGasCon)
      deallocate(m_waterSubCon)
      deallocate(m_tmpWaterSubCon)
      deallocate(m_gasExchange)
      deallocate(m_sedSubCon)  
      deallocate(m_tmpSedSubCon)
      deallocate(m_waterPOC)
      deallocate(m_tmpWaterPOC)
      deallocate(m_rChl2C)
      deallocate(m_Kt)
      deallocate(m_Kv)
      deallocate(m_Ks)
      deallocate(m_wrho)
      deallocate(m_dVsc)
      deallocate(m_Iab)
      deallocate(m_waterIce)
      deallocate(m_sedIce)
      deallocate(m_iceBubblePool)
      deallocate(m_sinkPOCPool)
      deallocate(m_frzCarbPool)
      deallocate(m_unfrzCarbPool)
      ! destroy the memory for air forcing data
      deallocate(m_airTemp)
      deallocate(m_airTempMax)
      deallocate(m_airTempMin)
      deallocate(m_airWind)
      deallocate(m_airRH)
      deallocate(m_airPrsn)
      deallocate(m_airPr)
      deallocate(m_airPs)
      deallocate(m_airSWRad)
      deallocate(m_airLWRad)
      ! destroy memory for irradiation data
      deallocate(m_aCO2)
      deallocate(m_aO3)
      deallocate(m_aAOD)
      ! destroy the memory for radiation vectors
      deallocate(m_fsphot)
      deallocate(m_wvln)
      ! destroy memory for discharge data
      deallocate(m_Qsi)
      deallocate(m_tQsi)
      deallocate(m_dQsi)
      deallocate(m_DOQsi)
      deallocate(m_DICQsi)
      deallocate(m_DOCQsi)
      deallocate(m_POCQsi)
      deallocate(m_SRPQsi)
      deallocate(m_Qso)
      deallocate(m_Qgw)
      ! destroy the memory for depth vectors
      deallocate(m_Zs)
      deallocate(m_Zw)               
      deallocate(m_dZw)
      deallocate(m_dZs)
      deallocate(m_Az)
      deallocate(m_dAz)
      deallocate(m_Rb0)
      ! destroy memory for intermediate variables
      deallocate(m_btmbflux)
   end subroutine
    
end module data_buffer_mod
