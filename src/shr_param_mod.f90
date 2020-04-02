module shr_param_mod
!----------------------------------------------------------------------------
! Sensitive parameters list below (Param Id should be identical with 
! optpar.dat)
!
!1 Param_Ks  : thermal conductivity of soil constituent (W/(m*K))
!2 Param_Cps : heat capacity of soil constituent (J/(kg*K))
!3 Param_Por : porosity of lake sediment (0~1)
!4 Param_Rous: density of soil solid particle (kg/m3)
!5 Param_OQ10: Methane oxidation Q10
!6 Param_Qch4: Oxidation potential when substrates are not limited (umol/m3/s)
!7 Param_Kch4: Michaelis-Menten constant (unit: umol/m3)
!8 Param_Ko2 : Michaelis-Menten constant (unit: umol/m3)
!9 Param_Re  : ebullition rate (unit: s-1)
!10 Param_PQ10n: Methane production Q10 of acetate fermentation
!11 Param_Rcn : the fraction of decomposed recalcitrant carbon (s-1)
!12 Param_DMP : recalcitrant organic matter dampening rate (m-1)
!13 Param_Rca : the fraction of aerobic decomposed carbon (s-1)
!14 Param_Vchs : Chla-specific light saturated growth rate (mg C mg Chl-1 d-1)
!15 Param_Vchl : Chla-specific light saturated growth rate (mg C mg Chl-1 d-1)
!16 Param_Klrs : metabolic loss rate coefficient (day-1)
!17 Param_Klrl : metabolic loss rate coefficient (day-1)
!18 Param_RDOMaq: aquatic DOM microbial degradation rate (d-1)
!19 Param_RDOCtr: terrestrail DOM microbial degradation rate (d-1)
!20 Param_DOCwt: groundwater DOC concentration (mol/m3) 
!21 Param_phAlphas: initial slope of P-E curve ((mol photons m-2 s-1)-1 d-1)
!22 Param_phAlphal: initial slope of P-E curve ((mol photons m-2 s-1)-1 d-1)
!23 Param_phBetas: photoinhibition parameter ((mol photons m-2 s-1)-1 d-1)
!24 Param_phBetal: photoinhibition parameter ((mol photons m-2 s-1)-1 d-1)
!25 Param_Ksrps: half-saturation for phosphorus limitation (umol/m3)
!26 Param_Ksrpl: half-saturation for phosphorus limitation (umol/m3)
!27 Param_Roun: snow density (kg/m3)
!28 Param_Feta: light attenuation correction factor for chla and CDOM
!29 Param_Wstr: wind shielding factor of mixing 
!30 Param_Ktscale: turbulence diffusivity scaling factor 
!31 Param_Hscale: sensible and latent heat transfer coefficent scaling factor
!----------------------------------------------------------------------------
   use shr_kind_mod,    only : r8, cx => SHR_KIND_CX 
   
   implicit none

   integer, parameter :: NPARAM = 31
   ! thermal related parameters
   integer, parameter :: Param_Ks = 1, Param_Cps = 2, Param_Por = 3, &
                         Param_Rous = 4
   ! methane oxidation related parameters
   integer, parameter :: Param_OQ10 = 5, Param_Qch4 = 6,  Param_Kch4 = 7, &
                         Param_Ko2 = 8
   ! 14C-enriched pool related parameters
   integer, parameter :: Param_Re = 9, Param_PQ10n = 10, Param_Rcn = 11, &
                         Param_DMP = 12
   ! aerobic decomposition parameters
   integer, parameter :: Param_Rca = 13
   ! phytoplankton related parameters
   integer, parameter :: Param_Vchs = 14, Param_Vchl = 15
   integer, parameter :: Param_Klrs = 16, Param_Klrl = 17
   ! allochthonous carbon parameters
   integer, parameter :: Param_RDOMaq = 18, Param_RDOMtr = 19
   integer, parameter :: Param_DOCwt = 20
   ! photo-response parameters
   integer, parameter :: Param_phAlphas = 21, Param_phAlphal = 22
   integer, parameter :: Param_phBetas = 23, Param_phBetal = 24
   ! phosphorus limitation
   integer, parameter :: Param_Ksrps = 25, Param_Ksrpl = 26
   ! ice phenology
   integer, parameter :: Param_Roun = 27
   ! thermodynamics
   integer, parameter :: Param_Feta = 28
   integer, parameter :: Param_Wstr = 29, Param_Ktscale = 30
   integer, parameter :: Param_Hscale = 31
   ! sensitive parameter collection
   real(r8) :: sa_params(NPARAM)

contains
   subroutine LoadSensitiveParameters(params)
      implicit none
      real(r8), intent(in) :: params(NPARAM)

      sa_params = params
   end subroutine

end module shr_param_mod
