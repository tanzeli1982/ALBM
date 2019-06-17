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
   use shr_ctrl_mod,    only : NPOC, NLAKTYPE
   public
   ! lake type identifier
   integer, parameter :: temperate_lake = 001
   integer, parameter :: small_ppk = 001, large_ppk = 002
   integer, parameter :: pasC = 001, actC = 002
   ! freshwater lake pH
   real(r8), parameter :: LKpH(NLAKTYPE) = (/6.5/)
   ! CH4 oxidation Reference temperature (K)
   real(r8), parameter :: Tor = 267.65
   ! CH4 production reference temperature for passive matter (K)
   real(r8), parameter :: Tpr = 276.65
   ! CH4 production reference temperature for active matter
   real(r8), parameter :: Tpr_act = 273.15
   ! O2, CH4 and CO2 minimum dissolved concentrations (umol/m3)
   real(r8), parameter :: minDo2 = 1.0d+2
   real(r8), parameter :: minDch4 = 1.0d+1
   real(r8), parameter :: minDco2 = 1.0d+1
   ! the initial density of 14C-depleted organic matter (kg/m3)
   real(r8), parameter :: oldCarb0 = 29.3
   ! the suppression coefficient of O2 to methanogenesis (m3 water mol-1)
   real(r8), parameter :: etaO2 = 400.0
   ! minimum, maximum and optimum pH for CH4 production
   real(r8), parameter :: pHmin = 5.5, pHmax = 9.0, pHopt = 7.5
   ! The carbon-specific absorption coefficient by NAP at 440 nm (m2 gC-1)
   real(r8), parameter :: aNAP_440 = 0.1
   ! The chlorophill specific coefficient of phytoplankton (m2 mg-1)
   real(r8), parameter :: apico_676 = 0.023
   real(r8), parameter :: amicro_674 = 0.0086
   ! The diameter (m) and density (kg/m3) of detritus
   real(r8), parameter :: daDetrs = 8.0d-5
   real(r8), parameter :: dsDetrs = 1.04d+3
   ! The diameter (m) of small and large phytoplankton
   real(r8), parameter :: daPico = 3.0d-6
   real(r8), parameter :: daMicro = 1.0d-5
   ! The background backscattering coefficient (m-1)
   real(r8), parameter :: Bbbg = 1.7d-4
   ! The backscattering ratio of sea water
   real(r8), parameter :: Bbsw = 0.5_r8
   ! Proportion of biosynthate allocated to synthesis of biomass of the
   ! biosynthetic machinery (Geider et al., 1996)
   real(r8), parameter :: kE = 0.6_r8
   ! The minimum and maximum carbon to chlorophill ratio (mmol C mg chl-1)
   ! Wang et al. (2009)
   real(r8), parameter :: C2Chlmin(NPOC) = (/12.0, 12.0/)
   real(r8), parameter :: C2Chlmax(NPOC) = (/40.0, 30.0/)
   ! The slope of C:Chl ratio vs. growth rate (mg C mg chl-1 d)
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
   ! the maximum growth rate of phytoplankton at 0 celcius (d-1) 
   real(r8), parameter :: Vm0(NPOC) = (/0.4, 1.0/)
   ! temperature multiplier for phytoplankton growth
   real(r8), parameter :: ThetaG = 1.08
   ! temperature function parameters for phytoplankton growth
   real(r8), parameter :: kt_ppk(NPOC) = (/1.93841, 12.76830/)
   real(r8), parameter :: at_ppk(NPOC) = (/29.27777, 21.67022/)
   real(r8), parameter :: bt_ppk(NPOC) = (/0.28991, 0.21632/)
   ! default sinking velocity (m/s)
   real(r8), parameter :: Vs0(NPOC) = (/9.838d-8, 8.68d-6/)
   ! temperature multiplier for metabolic loss
   real(r8), parameter :: ThetaML = 1.073
   ! fraction of respiratioin relative to total metabolic loss
   real(r8), parameter :: Fres(NPOC) = (/0.8, 0.5/)
   ! fraction of excretion relative to non-respiration metabolic loss
   real(r8), parameter :: Fdom(NPOC) = (/0.1, 0.5/)
   ! temperature multiplier for DOC microbial mineralization
   real(r8), parameter :: ThetaCM = 1.073
   ! The O2 half-saturation constant for DOC microbial degradation (umol/m3)
   real(r8), parameter :: Ko2CM = 4.6875d+4
   ! The O2 half-saturation constant for POM degradation (umol/m3)
   !real(r8), parameter :: Ko2PM = 1.5625d+4
   ! temperature multiplier for POC decomposition
   real(r8), parameter :: ThetaPM = 1.073
   ! algae mortality rate due to hypoxia (day-1)
   real(r8), parameter :: Kdhyp = 0.8
   ! The CO2 half-saturation constant for photosynthesis (umol/m3)
   real(r8), parameter :: Kco2 = 6.163d+4
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
   ! The concentration of DOC in precipitation (gC/m3)
   real(r8), parameter :: DOCrf = 2.0
   ! aerial OC loading factor (gC/m/d)
   real(r8), parameter :: DOCae = 2.0
   ! ratio of DIC to DOC loading
   real(r8), parameter :: rDIC2DOC = 0.25
   ! The ratio of Chla:POC in streams (mg chla mmol C-1)
   real(r8), parameter :: rChla2POC = 0.036
   ! POC resuspension rate (umol/m2/s), if for chla this rate ranges
   ! at 0.1 to 1 mg/m2/d [Saloranta and Anderson, 2007]
   real(r8), parameter :: Scsed = 5.787d-3 
   ! flocculation rate (s-1)
   real(r8), parameter :: floc = 2.2d-8
end module bg_const_mod
