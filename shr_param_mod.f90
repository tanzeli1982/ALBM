module shr_param_mod
!----------------------------------------------------------------------------
! Sensitive parameters list below (Param Id should be identical with 
! optpar.dat)
!
!1  Param_Roun: snow density (kg/m3)
!2  Param_Feta: light attenuation correction factor
!3  Param_Wstr: wind shielding factor of mixing 
!4  Param_Hscale: heat transfer coefficeint scaling factor
!5  Param_Dscale: turbulent diffusivity scaling factor
!6  Param_TinDiff: temperature difference between inflow and forebay
!7  Param_Vm0s: the maximum growth rate of phytoplankton at 0 celcius (d-1) 
!8  Param_Vm0l: the maximum growth rate of phytoplankton at 0 celcius (d-1) 
!9  Param_Ksrps: half-saturation for phosphorus limitation (mol/m3)
!10 Param_Ksrpl: half-saturation for phosphorus limitation (mol/m3)
!11 Param_Klrs: metabolic loss rate coefficient (day-1)
!12 Param_Klrl: metabolic loss rate coefficient (day-1)
!13 Param_fDepMorts: mortality fraction of settled phytoplankton
!14 Param_fDepMortl: mortality fraction of settled phytoplankton
!15 Param_Re: ebullition rate (unit: s-1)
!16 Param_Ae: relative concentration at which ebullition begins
!17 Param_icebflux: ice bubble flux rate (s-1)
!18 Param_icebloss: ice bubble dissolution rate (s-1)
!19 Param_csedDMP: recalcitrant organic matter dampening rate (m-1)
!20 Param_Rcapas: passive OC decomposition rate through aerobic respiration (s-1)
!21 Param_Rcaact: active OC decomposition rate through aerobic respiration (s-1)
!22 Param_Rcpas: passive OC decomposition rate through methanogenesis (s-1)
!23 Param_Rcact: active OC decomposition rate through methanogenesis (s-1)
!24 Param_Rcold: old OC decomposition rate through methanogenesis (s-1)
!25 Param_PQ10pas: passive OC decomposition Q10 through methanogenesis
!26 Param_PQ10act: active OC decomposition Q10 through methanogenesis 
!27 Param_PQ10old: old OC decomposition Q10 through methanogenesis
!28 Param_cRx: change rate of sediment redox potential (mV d-1) 
!29 Param_Qch4: CH4 oxidation potential (mol/m3/s)
!30 Param_OQ10: Methane oxidation Q10
!31 Param_betaCH4: CH4 concentration exponent coefficient of CH4 oxidation (unitless)
!32 Param_lamdaO2: O2 inhibition coefficient of CH4 oxidation (unitless)
!33 Param_RcOMP: oxic CH4 production rate to photosynthesis (mol CH4/m3/s (gC/m3/s)-1)
!34 Param_MuMaxVeg: maximum growth rate of submerged macrophytes at 20 celsius (d-1)
!----------------------------------------------------------------------------
   use shr_kind_mod,    only : r8, cx => SHR_KIND_CX 
   
   implicit none

   integer, parameter :: NPARAM = 34
   ! thermal related parameters
   integer, parameter :: Param_Roun = 1, Param_Feta = 2, Param_Wstr = 3, &
                         Param_Hscale = 4, Param_Dscale = 5, &
                         Param_TinDiff = 6                     
   ! phytoplankton related parameters
   integer, parameter :: Param_Vm0s = 7, Param_Vm0l = 8, Param_Ksrps = 9, &
                         Param_Ksrpl = 10, Param_Klrs = 11, Param_Klrl = 12, &
                         Param_fDepMorts = 13, Param_fDepMortl = 14
   ! ebullition related parameters
   integer, parameter :: Param_Re = 15, Param_Ae = 16, Param_icebflux = 17, &
                         Param_icebloss = 18 
   ! sediment C pool related parameters
   integer, parameter :: Param_csedDMP = 19, Param_Rcapas = 20, Param_Rcaact = 21, &
                         Param_Rcpas = 22, Param_Rcact = 23, Param_Rcold = 24
   ! other sediment methanogenesis parameters
   integer, parameter :: Param_PQ10pas = 25, Param_PQ10act = 26, Param_PQ10old = 27, &
                         Param_cRx = 28
   ! methane oxidation related parameters
   integer, parameter :: Param_Qch4 = 29, Param_OQ10 = 30, Param_betaCH4 = 31, &
                         Param_lamdaO2 = 32
   ! oxic methane production parameter
   integer, parameter :: Param_RcOMP = 33
   ! submerged macrophytes parameters
   integer, parameter :: Param_MuMaxVeg = 34
   ! sensitive parameter collection
   real(r8) :: sa_params(NPARAM)
   ! maximum sample size
   integer :: NMAXSAMPLE
   ! first and last sample id to run
   integer :: sample_range(2)
   ! Monte Carlo sample file
   character(cx) :: mc_file
   ! sensitivity output time window
   integer :: SA_Start_Year, SA_Start_Month, SA_Start_Day
   integer :: SA_End_Year, SA_End_Month, SA_End_Day

contains
   subroutine LoadSensitiveParameters(params)
      implicit none
      real(r8), intent(in) :: params(NPARAM)

      sa_params = params
   end subroutine

end module shr_param_mod
