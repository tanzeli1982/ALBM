module bg_const_mod
!----------------------------------------------------------------------------
!  Purpose: biological and biogeochemical constants
!           Notice: normally [DOC] is measured on the basis of carbon, e.g. 
!           mg C L-1 or umol C L-1, not the total mass of DOM if not clearly 
!           indicated!!
!           All photochemical parameters are from "Seasonality of 
!           photochemical dissolved organic carbon mineralization and its 
!           relative contribution to pelagic CO2 production in northern 
!           lakes" (Vachon et al., 2016).
!----------------------------------------------------------------------------
   use shr_kind_mod,    only : r8
   use shr_ctrl_mod,    only : NPOC, NPOOL, NLAKTYPE
   public
   ! lake type identifier
   integer, parameter :: temperate = 001, boreal = 002, tundra = 003, &
                         yedoma = 004, arctic = 005, alpine = 006, &
                         eutrophic = 007, subtropic = 008, amazon = 009, &
                         rift = 010, saline = 011
   integer, parameter :: small_ppk = 001, large_ppk = 002
   integer, parameter :: pasC = 001, actC = 002, oldC = 003
   integer, parameter :: dormant_season = 0, onset_season = 1, &
                         grow_season = 2, offset_season = 3
   ! CH4 oxidation Reference temperature (K)
   real(r8), parameter :: Tor = 267.65
   ! O2 amplification coefficient for CH4 oxidation (m3/mol)
   real(r8), parameter :: betaO2 = 30._r8
   ! critical CH4 concentration for O2-inhibitation (mol/m3)
   real(r8), parameter :: Dch4_cr = 1.7d-3
   ! critical O2 concentration for developing anoxia (mol/m3)
   real(r8), parameter :: Do2_cr = 0.05_r8
   ! CH4 and O2 threshold for CH4 oxidation (mol/m3)
   real(r8), parameter :: Dch4min = 1.0d-5
   real(r8), parameter :: Do2min = 6.25d-3
   ! CH4 production reference temperature (K)
   real(r8), parameter :: Tpr(NPOOL) = (/276.65, 273.15, 273.15/)
   ! the initial density of 14C-depleted organic matter (kgC/m3)
   real(r8), parameter :: oldCarb0 = 29.3
   ! minimum and maximum allowable pH for CH4 production
   real(r8), parameter :: pHmin = 2.2, pHmax = 10.0
   ! chla-specific light absorption coefficient (self-shading) ((ug chl L-1)-1 m-1)
   real(r8), parameter :: aCHL = 0.0249_r8
   ! The carbon-specific absorption coefficient by NAP at 440 nm (m2 gC-1)
   real(r8), parameter :: aNAP_440 = 0.1
   ! The diameter (m) and density (kg/m3) of detritus
   real(r8), parameter :: daDetrs = 8.0d-5
   real(r8), parameter :: dsDetrs = 1.04d+3
   ! The diameter (m) of small and large phytoplankton
   real(r8), parameter :: daPico = 3.0d-6
   real(r8), parameter :: daMicro = 1.0d-5
   ! The minimum and maximum carbon to chlorophill ratio (g chl gC^-1)
   real(r8), parameter :: Chl2Cmax(NPOC) = (/0.030, 0.045/)
   real(r8), parameter :: Chl2Cmin(NPOC) = (/0.005, 0.008/)
   ! The slope of C:Chl ratio vs. growth rate (gC (g chl)-1 d)
   real(r8), parameter :: Kpc2chl(NPOC) = (/95, 70/)
   ! Apparant quantum yield of photochemical degradation
   ! units: AQYem1 (mol C mol photons-1), AQYm2 (nm-1)
   real(r8), parameter :: AQYem1(NLAKTYPE) = (/2.6d-4/)
   real(r8), parameter :: AQYm2(NLAKTYPE) = (/0.018/)
   ! The ratio of partial photo-oxidation to photo-mineralization
   !real(r8), parameter :: PPO2PM(NLAKTYPE) = (/2.547/)
   ! The spectral slope of CDOM absorption (nm-1)
   real(r8), parameter :: SCDOM(NLAKTYPE) = (/0.018/)
   ! The reference CDOM-specific UV absorption at 305 nm (m3 g C-1 m-1)
   real(r8), parameter :: SUVA305(NLAKTYPE) = (/2.0/)
   ! chlorophyll specific initial slope of P‐E curve (gC (g chl)^-1 (umol photons m-2 s-1 d)^-1) 
   real(r8), parameter :: alphaChl(NPOC) = (/0.28, 0.28/) 
   ! initial slope of P-E curve ((mol photons m-2 s-1)-1 d-1)
   real(r8), parameter :: phAlpha(NPOC) = (/3.0d3, 6.1d3/)
   ! photoinhibition parameter ((mol photons m-2 s-1)-1 d-1)
   real(r8), parameter :: phBeta(NPOC) = (/8.0d2, 2.0d2/)
   ! temperature multiplier for phytoplankton growth
   real(r8), parameter :: ThetaG = 1.08
   ! temperature function parameters for phytoplankton growth
   real(r8), parameter :: fTmax_ppk(NPOC) = (/2.2207, 2.4508/)
   real(r8), parameter :: Topt_ppk(NPOC) = (/30.0, 22.0/)
   real(r8), parameter :: Tmax_ppk(NPOC) = (/40.0, 35.0/)
   ! default sinking velocity (m/s)
   real(r8), parameter :: Vs0(NPOC) = (/9.838d-8, 8.68d-6/)
   ! swim velocity of small phytoplankton (m/s) 
   real(r8), parameter :: Vswim = 2.89d-6
   ! threshold of nutrient limitation for chemotaxis
   real(r8), parameter :: fsrp_vmdown = 0.67_r8
   ! threshold of nutrient limition for phototaxis
   real(r8), parameter :: fsrp_vmup = 0.75_r8
   ! threshold of light for phototaxis (umol/m2/s)
   real(r8), parameter :: ipar_crit = 0.1_r8
   ! threshold of hypoxia mortality (mol/m3)
   real(r8), parameter :: Ko2AR = 0.0625_r8
   ! temperature multiplier for metabolic loss
   real(r8), parameter :: ThetaML = 1.073
   ! fraction of respiratioin relative to total metabolic loss
   real(r8), parameter :: Fres(NPOC) = (/0.8, 0.5/)
   ! fraction of excretion relative to non-respiration metabolic loss
   real(r8), parameter :: Fdom(NPOC) = (/0.1, 0.5/)
   ! temperature multiplier for DOC microbial mineralization
   real(r8), parameter :: ThetaCM = 1.073
   ! The O2 half-saturation constant for microbial degradation (mol/m3)
   real(r8), parameter :: Ko2CM = 4.6875d-2
   ! The baseline heterotrophic respiration rate (mol/m3/s)
   real(r8), parameter :: rCO2_BOD = 3.62d-9
   ! The O2 half-saturation constant for POM degradation (mol/m3)
   !real(r8), parameter :: Ko2PM = 1.5625d-2
   ! temperature multiplier for POC decomposition
   real(r8), parameter :: ThetaPM = 1.073
   ! The CO2 half-saturation constant for photosynthesis (mol/m3)
   real(r8), parameter :: Kco2 = 6.163d-3
   ! stoichemistry of C:P in POM and in allochthonous DOM
   real(r8), parameter :: YC2P_POM = 106.0
   real(r8), parameter :: YC2P_DOM = 199.0
   ! Redfield C : DOM mass ratio (C106H175O42N16P) (g/g)
   real(r8), parameter :: YC2DOM = 0.5358_r8 
   ! Mole mass of DOM (g/mol)
   real(r8), parameter :: MasDOM = 2374.0_r8 
   ! The rate coefficient of density change to irradiance (kg/m3/s)
   real(r8), parameter :: dsc1(NPOC) = (/7.155d-3, 7.816d-3/)
   ! The minimum rate of density increase (kg/m3/s)
   real(r8), parameter :: dsc3 = 3.83d-4
   ! POC resuspension rate (umol/m2/s), if for chla this rate ranges
   ! at 0.1 to 1 mg/m2/d [Saloranta and Anderson, 2007]
   real(r8), parameter :: Scsed = 5.787d-3 
   ! P adsorption affinity at oxidized conditions (m3/gP)
   real(r8), parameter :: cKPAdsOxS = 1.5_r8
   ! max. reduction factor of P adsorption affinity
   real(r8), parameter :: fRedMaxS = 0.96_r8
   ! maximum amount of sorbed P in sediment (gP/m3)
   real(r8), parameter :: cRelPAdsDWS = 11.88_r8
   ! max. P adsorption per g Fe (gP/gFe)
   real(r8), parameter :: cRelPAdsFeS = 0.06_r8
   ! temperature q10 of macrophyte growth
   real(r8), parameter :: cQ10ProdVeg = 1.2_r8
   ! half-saturation PAR of macrophyte growth at 20 celsius (W/m2)
   real(r8), parameter :: hLRefVeg = 17.0_r8
   ! carrying capacity (maximum macrophyte biomass) (g/m2)
   real(r8), parameter :: cDCarrVeg = 300._r8
   ! dark respiration rate of submerged macrophytes at 20 celsius (d-1)
   real(r8), parameter :: kDRespVeg = 0.02_r8
   ! temperature q10 of respiration
   real(r8), parameter :: cQ10RespVeg = 2.0_r8
   ! fraction surviving outside of growing season
   real(r8), parameter :: fWinVeg = 0.3_r8
   ! mortality rate during growing season (d-1)
   real(r8), parameter :: kMortVegSum = 0.005_r8
   ! critical growing-degree-day summation value (d)  
   real(r8), parameter :: GDDsum_crit = 30._r8
   ! duration of onset/offset phase (d)
   real(r8), parameter :: cLengTrans = 15._r8
   ! root fraction outside of growing season
   real(r8), parameter :: fRootVegWin = 0.6_r8
   ! root fraction during growing season
   real(r8), parameter :: fRootVegSum = 0.1_r8   
   ! specific cover of macrophytes (1/(gDW/m2))
   real(r8), parameter :: cCovSpVeg = 0.005_r8
   ! C content of organic matter DW
   real(r8), parameter :: cCPerDW = 0.4_r8
   ! fraction of allochthonous OC into active C pool
   real(r8), parameter :: fTrActC = 0.1_r8
   ! fraction of phytoplankton debris into active C pool
   real(r8), parameter :: fPhyActC = 0.9_r8
   ! fraction of macrophytes debris into active C pool
   real(r8), parameter :: fVegActC = 0.75_r8
end module bg_const_mod
