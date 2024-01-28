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
   ! simulated hourly water and sediment temperature profiles (K)
   real(r4), allocatable :: m_tempwHist(:,:)
   real(r4), allocatable :: m_tempsHist(:,:,:)
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
   ! simulated hourly BGC fluxes
   real(r4), allocatable :: m_dfch4airHist(:), m_ebch4airHist(:)
   real(r4), allocatable :: m_dfch4sedHist(:,:), m_ebch4sedHist(:,:)
   real(r4), allocatable :: m_totGPPHist(:), m_totNPPHist(:)
   real(r4), allocatable :: m_totpch4sedHist(:,:)
   real(r4), allocatable :: m_pch4watHist(:,:), m_och4watHist(:,:)
   real(r4), allocatable :: m_depAtCHist(:,:)
   ! simulated hourly solute concentrations
   real(r4), allocatable :: m_ch4concHist(:,:)
   real(r4), allocatable :: m_o2concHist(:,:)
   real(r4), allocatable :: m_n2concHist(:,:)
   real(r4), allocatable :: m_co2concHist(:,:)
   real(r4), allocatable :: m_srpconcHist(:,:)
   real(r4), allocatable :: m_icebch4Hist(:)
   ! simulated hourly unfrozen C budgets
   real(r4), allocatable :: m_bvegHist(:,:)
   real(r4), allocatable :: m_biomassHist(:,:,:) 
   real(r4), allocatable :: m_chlaHist(:,:)
   real(r4), allocatable :: m_pascarbHist(:,:,:)
   real(r4), allocatable :: m_actcarbHist(:,:,:)
   real(r4), allocatable :: m_oldcarbHist(:,:,:)
   real(r4), allocatable :: m_sedEHLHist(:,:,:)
   ! incumbent water and sediment temperature profile (K)
   real(r8), allocatable :: m_waterTemp(:), m_sedTemp(:,:)
   ! temporary water and sediment temperature profile (K) 
   real(r8), allocatable :: m_tmpWaterTemp(:), m_tmpSedTemp(:,:)
   ! water and sediment ice profile (m)
   real(r8), allocatable :: m_waterIce(:), m_sedIce(:,:)
   ! incumbent and temporary solution concentrations in water (mol/m3)
   real(r8), allocatable :: m_waterSubCon(:,:), m_tmpWaterSubCon(:,:)
   ! incumbent and temporary solution concentrations in sediments (mol/m3)
   real(r8), allocatable :: m_sedSubCon(:,:,:), m_tmpSedSubCon(:,:,:)
   real(r8), allocatable :: m_sedSRPCon(:)
   ! redox potential of lake sediment (mV)
   real(r8), allocatable :: m_sedEHL(:,:)
   ! total heat diffusivity (m2/s)
   real(r8), allocatable :: m_Kt(:)
   real(r8), allocatable :: m_Ks(:,:)
   ! total momentum diffusivity (m2/s)
   real(r8), allocatable :: m_Kv(:)
   ! water density (kg/m3)
   real(r8), allocatable :: m_wrho(:)
   ! sediment porosity (fraction)
   real(r8), allocatable :: m_sedpor(:,:)
   ! Solar radiation in water column (W/m2)
   real(r8), allocatable :: m_Iab(:)
   real(r8), allocatable :: m_Idw(:)
   ! phytoplankton biomass (gC/m3)
   real(r8), allocatable :: m_waterPOC(:,:)
   ! Chl : C ratio profile (g Chl gC-1)
   real(r8), allocatable :: m_rChl2C(:,:)
   ! Chla concentration (gchl/m3)
   real(r8), allocatable :: m_chla(:)
   ! submerged macrophyte biomass (gDW/m2)
   real(r8), allocatable :: m_bedVegDW(:)
   ! bed vegetation cover fraction
   real(r8), allocatable :: m_fcovBedVeg(:)
   ! vegetation P uptake (gP/m2/s)
   real(r8), allocatable :: m_vegPUptake(:)
   ! gas exchange amount from bubble to water (mol/m3/s)
   real(r8), allocatable :: m_gasExchange(:,:)
   ! gas pools for bubbles trapped in the ice (mol) 
   real(r8), allocatable :: m_iceBubblePool(:)
   ! sediment carbon pools (gC/m3)
   real(r8), allocatable :: m_frzCarbPool(:,:,:)
   real(r8), allocatable :: m_unfrzCarbPool(:,:,:)
   ! autochthonous and allochthonous active OC deposition rate (gC/m2/s)
   real(r8), allocatable :: m_aqDepAtCarb(:)
   real(r8), allocatable :: m_trDepAtCarb(:)
   real(r8), allocatable :: m_vegDepAtCarb(:)
   ! biology-based downward scalar irradiance (mol photons m-2 s-1 nm-1)
   real(r8), allocatable :: m_fsphot(:,:)
   ! downward irradiance spectrum (nm)
   real(r8), allocatable :: m_wvln(:)
   ! PAR radiation (mol/m2/s)
   real(r8), allocatable :: m_Ipar(:)
   ! hourly air climate records above the lake
   real(r8), allocatable :: m_airTemp(:)
   real(r8), allocatable :: m_airRH(:)
   real(r8), allocatable :: m_airWind(:)
   real(r8), allocatable :: m_airPr(:)
   real(r8), allocatable :: m_airPs(:)
   real(r8), allocatable :: m_airSWRad(:)
   real(r8), allocatable :: m_airLWRad(:)
   ! hydrological conditions
   real(r8), allocatable :: m_dzsurf(:)
   real(r8), allocatable :: m_srp(:)
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
   ! lake depth-dependent salinity (g/kg)
   real(r8), allocatable :: m_SAL(:)
   ! bubble radius vector (m)
   real(r8), allocatable :: m_Rb0(:)
   ! bottom heat diffusivity (m2/s)
   real(r8), allocatable :: m_Ktb(:)
   ! sediment-water interface bubble and diffusion fluxes (mol/m2/s)
   real(r8), allocatable :: m_btmbflux(:,:)
   real(r8), allocatable :: m_btmdflux(:,:)
   ! air-water interface bubble and diffusion fluxes (mol/m2/s)
   real(r8), allocatable :: m_topbflux(:)
   real(r8), allocatable :: m_topdflux(:)
   ! upwelling bubble CH4 flux upder ice (mol/m2/s)
   real(r8) :: m_upbflux
   ! sediment free water top and bottom index
   integer, allocatable :: m_sedWaterTopIndex(:)
   integer, allocatable :: m_sedWaterBtmIndex(:)
   ! sediment column index
   integer, allocatable :: m_soilColInd(:)
   !integer, allocatable :: m_soilColBub(:)
   ! lake ice, gray ice and snow thickness (m)
   real(r8) :: m_Hice, m_Hsnow, m_Hgrayice 
   ! areal forest fraction
   real(r8) :: m_ftree, m_fwlnd
   ! SOC density (kg/m2)
   real(r8) :: m_soc
   ! initial DOC (mol/m3), SRP (mmol/m3) and algae (mg/m3)
   real(r8) :: m_chla0, m_DOC0, m_DOC1, m_SRP0
   ! mixing layer thickness (m)
   real(r8) :: m_Hmix(2)
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
   integer :: m_mixTopIndex, m_mixBotIndex

contains
   !------------------------------------------------------------------------------
   !
   ! Data buffer must be destructed before changing its size
   !
   !------------------------------------------------------------------------------
   subroutine InitializeDataBuffer()
      implicit none
      type(SimTime) :: time, otime
      integer :: simday, simmon, simyr
      integer :: ntin, ntout, ndayout, indx
      integer :: o3Indx(2), aodIndx(2)
      integer :: co2Indx(2), timeIndx(2)
      real(r8) :: loc180(2)

      loc180 = (/lake_info%longitude, lake_info%latitude/)

      time = SimTime(Start_Year, Start_Month, Start_Day, End_Year, &
                     End_Month,End_Day)
      simday = CalcRunningDays(time, Use_Leap)
      simmon = CalcRunningMonths(time)
      simyr = CalcRunningYears(time)

      indx = INDEX(forcing_tstep, 'hour')
      if (indx==1) then
         forcing_nhour = 1
      else
         read(forcing_tstep(1:indx-1),fmt=*) forcing_nhour
      end if
      ntin = 24 * simday / forcing_nhour 

      if (trim(run_mode)=='sensitivity') then
         otime = SimTime(SA_Start_Year, SA_Start_Month, SA_Start_Day, &
                         SA_End_Year, SA_End_Month, SA_End_Day)
      else
         otime = SimTime(Start_Year, Start_Month, Start_Day, End_Year, &
                         End_Month,End_Day)
      end if
      ndayout = CalcRunningDays(otime, Use_Leap)
      ntout = 24 * ndayout

      co2Indx = (/time%year0-1764, simyr/)
      o3Indx = (/12*(time%year0-1978)+time%month0, simmon/)
      aodIndx = (/12*(time%year0-1978)+time%month0, simmon/)

      ! allocate memory for the archive variables
      allocate(m_tempwHist(WATER_LAYER+1,ntout))
      allocate(m_tempsHist(NSCOL,SED_LAYER+1,ntout))
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
      allocate(m_dfch4airHist(ntout))
      allocate(m_ebch4airHist(ntout))
      allocate(m_dfch4sedHist(NSCOL,ntout))
      allocate(m_ebch4sedHist(NSCOL,ntout))
      allocate(m_totGPPHist(ntout))
      allocate(m_totNPPHist(ntout))
      allocate(m_totpch4sedHist(NSCOL,ntout))
      allocate(m_pch4watHist(WATER_LAYER+1,ntout))
      allocate(m_och4watHist(WATER_LAYER+1,ntout))
      allocate(m_depAtCHist(NSCOL,ntout))
      allocate(m_ch4concHist(WATER_LAYER+1,ntout))
      allocate(m_o2concHist(WATER_LAYER+1,ntout))
      allocate(m_n2concHist(WATER_LAYER+1,ntout))
      allocate(m_co2concHist(WATER_LAYER+1,ntout))
      allocate(m_srpconcHist(WATER_LAYER+1,ntout))
      allocate(m_icebch4Hist(ntout))
      allocate(m_bvegHist(WATER_LAYER+1,ntout))
      allocate(m_biomassHist(NPOC,WATER_LAYER+1,ntout))
      allocate(m_chlaHist(WATER_LAYER+1,ntout))
      allocate(m_pascarbHist(NSCOL,SED_LAYER+1,ntout))
      allocate(m_actcarbHist(NSCOL,SED_LAYER+1,ntout))
      allocate(m_oldcarbHist(NSCOL,SED_LAYER+1,ntout))
      allocate(m_sedEHLHist(NSCOL,SED_LAYER+1,ntout))
      ! allocate memory for the state variables
      allocate(m_waterTemp(WATER_LAYER+1))
      allocate(m_sedTemp(NSCOL,SED_LAYER+1))
      allocate(m_tmpWaterTemp(WATER_LAYER+1))
      allocate(m_tmpSedTemp(NSCOL,SED_LAYER+1))
      allocate(m_waterIce(WATER_LAYER+1))
      allocate(m_sedIce(NSCOL,SED_LAYER+1))
      allocate(m_iceBubblePool(NGAS))
      allocate(m_waterSubCon(NWSUB,WATER_LAYER+1))
      allocate(m_tmpWaterSubCon(NWSUB,WATER_LAYER+1))
      allocate(m_gasExchange(NGAS,WATER_LAYER+1))
      allocate(m_sedSubCon(NSSUB,NSCOL,SED_LAYER+1))
      allocate(m_tmpSedSubCon(NSSUB,NSCOL,SED_LAYER+1)) 
      allocate(m_sedSRPCon(NSCOL))
      allocate(m_sedEHL(NSCOL,SED_LAYER+1))
      allocate(m_Kt(WATER_LAYER+1))
      allocate(m_Kv(WATER_LAYER+1))
      allocate(m_Ks(NSCOL,SED_LAYER+1))
      allocate(m_wrho(WATER_LAYER+1))
      allocate(m_sedpor(NSCOL,SED_LAYER+1))
      allocate(m_Iab(WATER_LAYER+1))
      allocate(m_Idw(WATER_LAYER+2))
      allocate(m_waterPOC(NPOC,WATER_LAYER+1))
      allocate(m_rChl2C(NPOC,WATER_LAYER+1))
      allocate(m_chla(WATER_LAYER+1))
      allocate(m_bedVegDW(WATER_LAYER+1))
      allocate(m_fcovBedVeg(WATER_LAYER+1))
      allocate(m_vegPUptake(WATER_LAYER+1))
      allocate(m_frzCarbPool(NPOOL,NSCOL,SED_LAYER+1))
      allocate(m_unfrzCarbPool(NPOOL,NSCOL,SED_LAYER+1))
      allocate(m_aqDepAtCarb(NSCOL))
      allocate(m_trDepAtCarb(NSCOL))
      allocate(m_vegDepAtCarb(NSCOL))
      ! allocate memory for air forcing data
      allocate(m_airTemp(ntin))
      allocate(m_airRH(ntin))
      allocate(m_airWind(ntin))
      allocate(m_airPr(ntin))
      allocate(m_airPs(ntin))
      allocate(m_airSWRad(ntin))
      allocate(m_airLWRad(ntin))
      ! allocate memory for hydrology data
      allocate(m_dzsurf(simday))
      allocate(m_srp(simday))
      ! allocate memory for irradiation data
      allocate(m_aCO2(simyr))
      allocate(m_aO3(simmon))
      allocate(m_aAOD(simmon))
      ! allocate memory for radiation vector
      allocate(m_fsphot(NSPCTM,WATER_LAYER+1))
      allocate(m_wvln(NSPCTM))
      allocate(m_Ipar(WATER_LAYER+1))
      ! allocate depth vector
      allocate(m_Zs(SED_LAYER+1))
      allocate(m_Zw(WATER_LAYER+1))
      allocate(m_dZw(WATER_LAYER+1))
      allocate(m_dZs(SED_LAYER+1))
      allocate(m_Az(WATER_LAYER+1))
      allocate(m_dAz(WATER_LAYER+1))
      allocate(m_SAL(WATER_LAYER+1))
      allocate(m_soilColInd(WATER_LAYER+1))
      allocate(m_Rb0(NRLAYER+1))
      ! allocate intermediate vector
      allocate(m_Ktb(WATER_LAYER+1))
      allocate(m_btmbflux(NGAS,NSCOL))
      allocate(m_btmdflux(NSSUB,NSCOL))
      allocate(m_topbflux(NGAS))
      allocate(m_topdflux(NGAS))
      allocate(m_sedWaterTopIndex(NSCOL))
      allocate(m_sedWaterBtmIndex(NSCOL))

      ! initialize depth vectors
      call ConstructDepthVector(lake_info, m_Zw, m_dZw, m_Zs, m_dZs)
      call ReadLakeBathymetry(lake_info, m_Zw, m_dZw, m_Az, m_dAz, &
               m_SAL, m_soilColInd)
      !call ReadSoilColBubFlag(m_soilColBub)
      call SetSedPorosity(m_Zs, m_dZs, m_sedpor)
      call ReadStaticData(veg_file, 'ftree', 'longitude', 'latitude', &
               loc180, m_ftree)
      call ReadStaticData(wlnd_file, 'glwd', 'longitude', 'latitude', &
               loc180, m_fwlnd)
      call ReadStaticData(soc_file, 'soc', 'longitude', 'latitude', &
               loc180, (/5.d-1,5.d-1/), Ncinterp, m_soc, 0d0)
      call ReadStaticData(tref_file, "t2m", 'longitude', 'latitude', &
               loc180, (/5.d-1,5.d-1/), Ncorigin, m_radPars%tref)
      ! initialize formal air forcing data
      call ReadSiteTSData(tas_file, 'tas', lake_info, time, m_airTemp)
      call ReadSiteTSData(hurs_file, 'hurs', lake_info, time, m_airRH)
      call ReadSiteTSData(pr_file, 'pr', lake_info, time, m_airPr)
      call ReadSiteTSData(ps_file, 'ps', lake_info, time, m_airPs)
      call ReadSiteTSData(wind_file, 'sfcwind', lake_info, time, m_airWind)
      call ReadSiteTSData(rsds_file, 'rsds', lake_info, time, m_airSWRad)
      call ReadSiteTSData(rlds_file, 'rlds', lake_info, time, m_airLWRad)
     
      call ReadGlobalTSData(co2_file, 'co2_rcp26', co2Indx, m_aCO2) 
      call Read2DTSData(o3_file, 'tro3', 'lon', 'lat', loc180, o3Indx, m_aO3)
      call Read2DTSData(aod_file, 'AOD_550', 'lon', 'lat', loc180, aodIndx, &
               m_aAOD)
      ! hydrology and chemistry input data
      call ReadSiteHydroTSData(hydro_file, 'dzsurf', lake_info, time, m_dzsurf)
      call ReadSiteHydroTSData(hydro_file, 'SRP', lake_info, time, m_srp)
      ! rescale tree cover and wetland fraction
      call RescaleTreeCoverFraction(m_ftree)
      call RescaleWetlandFraction(m_fwlnd)
      ! units conversion
      m_airPr = 1.0d-3 * m_airPr             ! convert to m/s (water)
      m_aO3 = 1.0d-3 * m_aO3                 ! convert to 1000 DU
   end subroutine

   subroutine DestructDataBuffer()
      implicit none

      ! destroy the memory for the archive variables
      deallocate(m_tempwHist)
      deallocate(m_tempsHist)
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
      deallocate(m_dfch4airHist)
      deallocate(m_ebch4airHist)
      deallocate(m_dfch4sedHist)
      deallocate(m_ebch4sedHist)
      deallocate(m_totGPPHist)
      deallocate(m_totNPPHist)
      deallocate(m_totpch4sedHist)
      deallocate(m_pch4watHist)
      deallocate(m_och4watHist)
      deallocate(m_depAtCHist)
      deallocate(m_ch4concHist)
      deallocate(m_o2concHist)
      deallocate(m_n2concHist)
      deallocate(m_co2concHist)
      deallocate(m_srpconcHist)
      deallocate(m_icebch4Hist)
      deallocate(m_bvegHist)
      deallocate(m_biomassHist)
      deallocate(m_chlaHist)
      deallocate(m_pascarbHist)
      deallocate(m_actcarbHist)
      deallocate(m_oldcarbHist)
      deallocate(m_sedEHLHist)
      ! destroy the memory for the state variables
      deallocate(m_waterTemp)
      deallocate(m_sedTemp)
      deallocate(m_tmpWaterTemp)
      deallocate(m_tmpSedTemp)
      deallocate(m_waterSubCon)
      deallocate(m_tmpWaterSubCon)
      deallocate(m_gasExchange)
      deallocate(m_sedSubCon)  
      deallocate(m_tmpSedSubCon)
      deallocate(m_sedSRPCon)
      deallocate(m_sedEHL)
      deallocate(m_waterPOC)
      deallocate(m_rChl2C)
      deallocate(m_chla)
      deallocate(m_bedVegDW)
      deallocate(m_vegPUptake)
      deallocate(m_fcovBedVeg)
      deallocate(m_Kt)
      deallocate(m_Kv)
      deallocate(m_Ks)
      deallocate(m_wrho)
      deallocate(m_sedpor)
      deallocate(m_Iab)
      deallocate(m_Idw)
      deallocate(m_waterIce)
      deallocate(m_sedIce)
      deallocate(m_iceBubblePool)
      deallocate(m_frzCarbPool)
      deallocate(m_unfrzCarbPool)
      deallocate(m_aqDepAtCarb)
      deallocate(m_trDepAtCarb)
      deallocate(m_vegDepAtCarb)
      ! destroy the memory for air forcing data
      deallocate(m_airTemp)
      deallocate(m_airWind)
      deallocate(m_airRH)
      deallocate(m_airPr)
      deallocate(m_airPs)
      deallocate(m_airSWRad)
      deallocate(m_airLWRad)
      ! destroy memory for hydrology data
      deallocate(m_dzsurf)
      deallocate(m_srp)
      ! destroy memory for irradiation data
      deallocate(m_aCO2)
      deallocate(m_aO3)
      deallocate(m_aAOD)
      ! destroy the memory for radiation vectors
      deallocate(m_fsphot)
      deallocate(m_wvln)
      deallocate(m_Ipar)
      ! destroy memory for discharge data
      ! destroy the memory for depth vectors
      deallocate(m_Zs)
      deallocate(m_Zw)               
      deallocate(m_dZw)
      deallocate(m_dZs)
      deallocate(m_Az)
      deallocate(m_dAz)
      deallocate(m_SAL)
      deallocate(m_soilColInd)
      deallocate(m_Rb0)
      ! destroy memory for intermediate variables
      deallocate(m_Ktb)
      deallocate(m_btmbflux)
      deallocate(m_btmdflux)
      deallocate(m_topbflux)
      deallocate(m_topdflux)
      deallocate(m_sedWaterTopIndex)
      deallocate(m_sedWaterBtmIndex)
   end subroutine
    
end module data_buffer_mod
