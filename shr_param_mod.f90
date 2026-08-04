module shr_param_mod
!----------------------------------------------------------------------------
! Sensitive parameters list below (Param Id should be identical with 
! optpar.dat)
!
!1  Param_Roun: snow density (kg/m3)
!2  Param_Feta: light attenuation correction factor
!3  Param_Wstr: wind shielding factor of mixing 
!4  Param_Hscale: heat transfer coefficeint scaling factor
!5  Param_Vm0s: the maximum growth rate of phytoplankton at 0 celcius (d-1) 
!6  Param_Vm0l: the maximum growth rate of phytoplankton at 0 celcius (d-1) 
!7  Param_Ksrps: half-saturation for phosphorus limitation (mol/m3)
!8  Param_Ksrpl: half-saturation for phosphorus limitation (mol/m3)
!9  Param_Klrs: metabolic loss rate coefficient (day-1)
!10 Param_Klrl: metabolic loss rate coefficient (day-1)
!11 Param_fDepMorts: mortality fraction of settled phytoplankton
!12 Param_fDepMortl: mortality fraction of settled phytoplankton
!13 Param_Re: ebullition rate (unit: s-1)
!14 Param_Ae: relative concentration at which ebullition begins
!15 Param_icebflux: ice bubble flux rate (s-1)
!16 Param_icebloss: ice bubble dissolution rate (s-1)
!17 Param_csedDMP: recalcitrant organic matter dampening rate (m-1)
!18 Param_Rcapas: passive OC decomposition rate through aerobic respiration (s-1)
!19 Param_Rcaact: active OC decomposition rate through aerobic respiration (s-1)
!20 Param_Rcpas: passive OC decomposition rate through methanogenesis (s-1)
!21 Param_Rcact: active OC decomposition rate through methanogenesis (s-1)
!22 Param_Rcold: old OC decomposition rate through methanogenesis (s-1)
!23 Param_PQ10pas: passive OC decomposition Q10 through methanogenesis
!24 Param_PQ10act: active OC decomposition Q10 through methanogenesis 
!25 Param_PQ10old: old OC decomposition Q10 through methanogenesis
!26 Param_cRx: change rate of sediment redox potential (mV d-1) 
!27 Param_Qch4: CH4 oxidation potential (mol/m3/s)
!28 Param_OQ10: Methane oxidation Q10
!29 Param_betaCH4: CH4 concentration exponent coefficient of CH4 oxidation (unitless)
!30 Param_lamdaO2: O2 inhibition coefficient of CH4 oxidation (unitless)
!31 Param_RcOMP: oxic CH4 production rate to photosynthesis (mol CH4/m3/s (gC/m3/s)-1)
!32 Param_MuMaxVeg: maximum growth rate of submerged macrophytes at 20 celsius (d-1)
!----------------------------------------------------------------------------
   use shr_kind_mod,    only : r8, cx => SHR_KIND_CX 
   
   implicit none

   integer, parameter :: NPARAM = 33
   ! thermal related parameters
   integer, parameter :: Param_Roun = 1, Param_Feta = 2, Param_Wstr = 3, &
                         Param_Hscale = 4, Param_Dscale = 5                      
   ! phytoplankton related parameters
   integer, parameter :: Param_Vm0s = 6, Param_Vm0l = 7, Param_Ksrps = 8, &
                         Param_Ksrpl = 9, Param_Klrs = 10, Param_Klrl = 11, &
                         Param_fDepMorts = 12, Param_fDepMortl = 13
   ! ebullition related parameters
   integer, parameter :: Param_Re = 14, Param_Ae = 15, Param_icebflux = 16, &
                         Param_icebloss = 17 
   ! sediment C pool related parameters
   integer, parameter :: Param_csedDMP = 18, Param_Rcapas = 19, Param_Rcaact = 20, &
                         Param_Rcpas = 21, Param_Rcact = 22, Param_Rcold = 23
   ! other sediment methanogenesis parameters
   integer, parameter :: Param_PQ10pas = 24, Param_PQ10act = 25, Param_PQ10old = 26, &
                         Param_cRx = 27        
   ! methane oxidation related parameters
   integer, parameter :: Param_Qch4 = 28, Param_OQ10 = 29, Param_betaCH4 = 30, &
                         Param_lamdaO2 = 31
   ! oxic methane production parameter
   integer, parameter :: Param_RcOMP = 32
   ! submerged macrophytes parameters
   integer, parameter :: Param_MuMaxVeg = 33
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
