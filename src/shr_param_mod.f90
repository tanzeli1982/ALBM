module shr_param_mod
!----------------------------------------------------------------------------
! Sensitive parameters list below and defined in optpar.dat 
!
!----------------------------------------------------------------------------
   use shr_kind_mod,    only : r8, cx => SHR_KIND_CX
   use shr_ctrl_mod,    only : NPTK, NSEDC
   
   implicit none

   ! thermal related parameters
   integer, parameter :: NPARAM_temp = 5     ! number of thermal related parameters
   real(r8) :: PAR_Roun = 200.0_r8           ! snow density (kg/m3)
   real(r8) :: PAR_Feta = 1.0_r8             ! light attenuation correction factor
   real(r8) :: PAR_Wstr = 1.0_r8             ! wind shielding factor of mixing
   real(r8) :: PAR_Hscale = 1.0_r8           ! heat transfer coefficeint scaling factor
   real(r8) :: PAR_Dscale = 1.0_r8           ! turbulent diffusivity scaling factor
   
   ! phytoplankton related parameters
   integer, parameter :: NPARAM_phyto = 6 * NPTK               ! number of phytoplankton related parameters
   real(r8) :: PAR_Vmax(NPTK) = (/2.4, 1.0, 2.0/)              ! the maximum growth rate of phytoplankton at 20 celsius (d-1)
   real(r8) :: PAR_Ksrp(NPTK) = (/4.d-3, 1.d-3, 2.d-3/)        ! half-saturation for phosphorus limitation (mol/m3)
   real(r8) :: PAR_Klr(NPTK) = (/0.1, 0.03, 0.075/)            ! respiration rate coefficient at 20 celsius (day-1) 
   real(r8) :: PAR_fDepMort(NPTK) = (/0.2, 0.2, 0.8/)          ! settled fraction of deposited phytoplankton
   real(r8) :: PAR_SigT(NPTK) = (/15.0, 11.5, 16.0/)           ! temperature sensitivity parameter of phytoplankton growth (celsius)
   real(r8) :: PAR_Topt(NPTK) = (/22.0, 28.0, 18.0/)           ! optimum temperature parameter of phytoplankton growth (celsius)

   ! zooplankton related parameters
   integer, parameter :: NPARAM_zoo = 5
   real(r8) :: PAR_Gmax = 2._r8                                ! zooplankton max grazing rate at 20 celsius (m3/gC/d)
   real(r8) :: PAR_zSigT = 13.0_r8                             ! temperature sensitivity parameter of zooplankton grazing (celsius)
   real(r8) :: PAR_zTopt = 22.0_r8                             ! optimum temperature parameter of zooplankton grazing (celsius)
   real(r8) :: PAR_kMortZoo = 0.06_r8                          ! zooplankton mortality rate (d^-1)
   real(r8) :: PAR_kDRespZoo = 0.15_r8                         ! zooplankton maintenance respiration rate at 20*C (d^-1) 

   ! ebullition related parameters
   integer, parameter :: NPARAM_ebu = 4      ! number of ebullition related parameters
   real(r8) :: PAR_Re = 2.06d-04             ! ebullition rate (unit: s-1)
   real(r8) :: PAR_Ae = 0.4                  ! relative concentration at which ebullition begins
   real(r8) :: PAR_icebflux = 1.70d-06       ! ice-trapped bubble flux rate (s-1)
   real(r8) :: PAR_icebloss = 1.90d-07       ! ice-trapped bubble dissolution rate (s-1)

   ! sediment C pool related parameters
   integer, parameter :: NPARAM_cpool = 2 + 3 * NSEDC                         ! number of C pool related parameters
   real(r8) :: PAR_csedDMP = 1.0_r8                                           ! sediment OC e-folding dampening rate in depth (m-1)
   real(r8) :: PAR_Rc_ar(NSEDC) = (/6.34d-11, 2.38d-09, 2.38d-09/)            ! OC decomposition rate through aerobic respiration (s-1)
   real(r8) :: PAR_Rc_ch4(NSEDC) = (/1.39d-12, 3.64d-10, 3.64d-10/)           ! OC decomposition rate through methanogenesis (s-1)
   real(r8) :: PAR_PQ10(NSEDC) = (/3.7, 2.0, 1.2/)                            ! Q10 of sediment OC methanogenesis
   real(r8) :: PAR_cRx = 100._r8                                              ! change rate of sediment redox potential (mV d-1)

   ! methane oxidation and OMP related parameters
   integer, parameter :: NPARAM_ch4 = 5      ! number of methane oxidation and OMP related parameters
   real(r8) :: PAR_Qch4 = 6.36d-07           ! CH4 oxidation potential (mol/m3/s)
   real(r8) :: PAR_OQ10 = 2.0_r8             ! CH4 oxidation Q10
   real(r8) :: PAR_betaCH4 = 0.79_r8         ! CH4 concentration exponent coefficient of CH4 oxidation (unitless)
   real(r8) :: PAR_lamdaO2 = 10._r8          ! O2 inhibition coefficient of CH4 oxidation (unitless)
   real(r8) :: PAR_RcOMP = 0._r8             ! oxic CH4 production rate (mol CH4/m3/s (gC/m3/s)-1)

   ! submerged macrophytes parameters
   integer, parameter :: NPARAM_smv = 1      ! number of submerged macrophytes related parameters
   real(r8) :: PAR_MuMaxVeg = 0.08_r8        ! maximum growth rate of submerged macrophytes at 20 celsius (d-1)

   ! sensitive parameter collection
   integer, parameter :: NPARAM = NPARAM_temp + NPARAM_phyto + NPARAM_zoo + &
                                  NPARAM_ebu + NPARAM_cpool + NPARAM_ch4 + &
                                  NPARAM_smv 

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
   subroutine LoadSensitiveParameters(sa_params)
      implicit none
      real(r8), intent(in) :: sa_params(NPARAM)
      integer :: indx, ii

      indx = 0
      PAR_Roun = sa_params(indx+1)
      PAR_Feta = sa_params(indx+2)
      PAR_Wstr = sa_params(indx+3)
      PAR_Hscale = sa_params(indx+4)
      PAR_Dscale = sa_params(indx+5) 

      indx = indx + NPARAM_temp
      do ii = 1, NPTK, 1
         PAR_Vmax(ii) = sa_params(indx+ii)
         PAR_Ksrp(ii) = sa_params(indx+NPTK+ii)
         PAR_Klr(ii) = sa_params(indx+2*NPTK+ii)
         PAR_fDepMort(ii) = sa_params(indx+3*NPTK+ii)
         PAR_SigT(ii) = sa_params(indx+4*NPTK+ii)
         PAR_Topt(ii) = sa_params(indx+5*NPTK+ii)
      end do

      indx = indx + NPARAM_phyto
      PAR_Gmax = sa_params(indx+1)
      PAR_zSigT = sa_params(indx+2)
      PAR_zTopt = sa_params(indx+3)
      PAR_kMortZoo = sa_params(indx+4)
      PAR_kDRespZoo = sa_params(indx+5)

      indx = indx + NPARAM_zoo
      PAR_Re = sa_params(indx+1)
      PAR_Ae = sa_params(indx+2)
      PAR_icebflux = sa_params(indx+3)
      PAR_icebloss = sa_params(indx+4)

      indx = indx + NPARAM_ebu
      PAR_csedDMP = sa_params(indx+1)
      do ii = 1, NSEDC, 1
         PAR_Rc_ar(ii) = sa_params(indx+1+ii)
         PAR_Rc_ch4(ii) = sa_params(indx+NSEDC+1+ii)
         PAR_PQ10(ii) = sa_params(indx+2*NSEDC+1+ii)
      end do 
      PAR_cRx = sa_params(indx+3*NSEDC+2)

      indx = indx + NPARAM_cpool
      PAR_Qch4 = sa_params(indx+1)
      PAR_OQ10 = sa_params(indx+2)
      PAR_betaCH4 = sa_params(indx+3)
      PAR_lamdaO2 = sa_params(indx+4)
      PAR_RcOMP = sa_params(indx+5)

      indx = indx + NPARAM_ch4
      PAR_MuMaxVeg = sa_params(indx+1)

   end subroutine

end module shr_param_mod
