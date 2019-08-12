module radiation_mod
!---------------------------------------------------------------------------------
!                        SMARTS MODEL, Version 2.9.5
! Purpose: CALCULATES CLEAR SKY SPECTRAL SOLAR IRRADIANCES FROM 280 TO 4000 nm.
!          Research code written by C. GUEYMARD, Solar Consulting Services
!
!          This module calculates solar radiations just below water surface.
!           References include "Sunlight-induced carbon dioxide emissions from
!           inland waters" (Koehler et al., 2014), "Sunlight controls water
!           column processing of carbon in arctic fresh waters" (Cory et al.,
!           2014), "Light and Photosynthesis in Aquatic Ecosystems" (Kirk,
!           2011), and etc.
!
!           When ice layers form, only the broad band absorption coefficients of
!           water, ice and snow will be counted.
!---------------------------------------------------------------------------------
   use shr_kind_mod,       only : i8, r8
   use shr_ctrl_mod, e8 => SHR_CTRL_E8, e30 => INFINITESIMAL_E8
   use math_utilities_mod, only : BinarySearch
   use phy_const_mod
   use bg_utilities_mod,   only : CalcAcWater, CalcAcIce
   use bg_utilities_mod,   only : CalcAcSnow, CalcAcGrayIce
   use radiation_io_mod
   use data_buffer_mod,    only : m_radPars, m_surfData, m_wvln

   implicit none
   private
   public :: InitializeSmartsModule, DestructSmartsModule
   public :: CalcClearSkyIrradiance, CorrIrradianceByCloud
   public :: CorrIrradianceByReflection
   public :: CorrIrradianceForSnow
   public :: CorrIrradianceForGrayIce
   public :: abW, abI, abN, abE
   public :: abCDOM, abAP, abNAP, bsP
   public :: Abd_fsnow, Abd_msnow
   public :: fgphot, frdif, fbphot
   public :: mem_pico, mem_micro
   real(r8), parameter :: NLosch = 2.6867775d+19
   real(r8), parameter :: epsilm = 1.0d-6
   real(r8), allocatable :: Zref50(:)  ! reference altitude (km)
   real(r8), allocatable :: Zref41(:)  ! reference altitude (km)
   real(r8), allocatable :: O3a0(:)    ! coefficient a1
   real(r8), allocatable :: O3a1(:)    ! coefficient a0
   real(r8), allocatable :: O3Min(:)   ! 
   real(r8), allocatable :: O3Max(:)   !
   real(r8), allocatable :: TO3Min(:)  ! minimum O3 temperature (K)
   real(r8), allocatable :: TO3Max(:)  ! maximum O3 temperature (K)
   real(r8), allocatable :: Tref(:,:)  ! reference temperature (K)
   real(r8), allocatable :: Pref(:,:)  ! reference air pressure (hPa)
   real(r8), allocatable :: RHref(:,:) ! reference relative humidity (%)
   real(r8), allocatable :: O3ref(:,:) ! reference ozone
   real(r8), allocatable :: Wref(:,:)  !
   real(r8), allocatable :: TO3(:,:)   ! reference ozone temperature (K)
   real(r8), dimension(4) :: C1, C2, C3, C4, C5, C6, D1, D2, D3, D4, &
      D5, D6, BP00, BP01, BP02, BP10, BP11, BP12, BP20, BP21, BP22, &
      BP30, BP31, BP32, BQ00, BQ01, BQ02, BQ10, BQ11, BQ12, BQ20, &
      BQ21, BQ22, AG41, AG42, AG00, AG01, AG02, AG10, AG11, AG12, &
      AG20, AG21, AG22, AG30, AG31, AG32, AG40
   real(r8), allocatable :: fgphot(:)     ! surface global downward irradiance
   real(r8), allocatable :: frdif(:)      ! diffuse fraction
   real(r8), allocatable :: fbphot(:)     ! backscattering irradiance
   real(r8), allocatable :: abCDOM(:)     ! CDOM absorption (m-1)
   real(r8), allocatable :: abAP(:)       ! algae absorption (m-1)
   real(r8), allocatable :: abNAP(:)      ! detritus absorption (m-1)
   real(r8), allocatable :: bsP(:)        ! particle backscattering (m-1) 
   real(r8), allocatable :: vH0(:)        ! spectrum coefficient
   real(r8), allocatable :: abW(:)        ! water absorption coefficient
   real(r8), allocatable :: abI(:)        ! ice absorption coefficient
   real(r8), allocatable :: abN(:)        ! snow absorption coefficient
   real(r8), allocatable :: abE(:)        ! gray ice absorption coefficient
   real(r8), allocatable :: mem_pico(:)   ! normalized highlight Apico
   real(r8), allocatable :: mem_micro(:)  ! normalized highlight Amicro
   real(r8), allocatable :: Abd_fsnow(:)  ! albedo band of fresh snow
   real(r8), allocatable :: Abd_msnow(:)  ! albedo band of melting snow
   real(r8), allocatable :: Abs_H2O(:,:)  ! H2O absorption band
   real(r8), allocatable :: Abs_O2(:,:)   ! O2 absorption band
   real(r8), allocatable :: Abs_CH4(:,:)  ! CH4 absorption band
   real(r8), allocatable :: Abs_CO(:,:)   ! CO absorption band
   real(r8), allocatable :: Abs_N2O(:,:)  ! N2O absorption band
   real(r8), allocatable :: Abs_CO2(:,:)  ! CO2 absorption band
   real(r8), allocatable :: Abs_N2(:,:)   ! N2 absorption band
   real(r8), allocatable :: Abs_O4(:,:)   ! O2-O2 absorption band
   real(r8), allocatable :: Abs_HNO3(:,:) ! HNO3 absorption band
   real(r8), allocatable :: Abs_NO(:,:)   ! NO absorption band
   real(r8), allocatable :: Abs_NO2(:,:)  ! NO2 absorption band
   real(r8), allocatable :: Abs_NO3(:,:)  ! NO3 absorption band
   real(r8), allocatable :: Abs_SO2U(:,:) ! SO2 UV absorption band
   real(r8), allocatable :: Abs_SO2I(:,:) ! SO2 IR absorption band
   real(r8), allocatable :: Abs_O3UV(:,:) ! O3 UV absorption band
   real(r8), allocatable :: Abs_O3IR(:,:) ! O3 UV absorption band
   real(r8), allocatable :: Abs_NH3(:,:)  ! NH3 absorption band
   real(r8), allocatable :: Abs_BrO(:,:)  ! BrO absorption band
   real(r8), allocatable :: Abs_CH2O(:,:) ! CH2O absorption band
   real(r8), allocatable :: Abs_HNO2(:,:) ! HNO2 absorption band
   real(r8), allocatable :: Abs_ClNO(:,:) ! ClNO3 absorption band
   real(r8) :: ESC                        ! solar constant (W/m2)

contains
   subroutine InitializeSmartsModule()
      implicit none

      allocate(Zref50(50))
      allocate(Zref41(41))
      allocate(O3a0(41))
      allocate(O3a1(41))
      allocate(O3Min(41))
      allocate(O3Max(41))
      allocate(TO3Min(41))
      allocate(TO3Max(41))
      allocate(Tref(50,10))
      allocate(Pref(50,10))
      allocate(RHref(50,10))
      allocate(O3ref(50,10))
      allocate(Wref(50,10))
      allocate(TO3(50,10))
      allocate(fgphot(NSPCTM))
      allocate(frdif(NSPCTM))
      allocate(fbphot(NSPCTM))
      allocate(abCDOM(NSPCTM))
      allocate(abAP(NSPCTM))
      allocate(abNAP(NSPCTM))
      allocate(bsP(NSPCTM))
      allocate(vH0(NSPCTM))
      allocate(abW(NSPCTM))
      allocate(abI(NSPCTM))
      allocate(abN(NSPCTM))
      allocate(abE(NSPCTM))
      allocate(Abd_fsnow(NSPCTM))
      allocate(Abd_msnow(NSPCTM))
      allocate(Abs_H2O(1722,17))
      allocate(Abs_O2(955,2))
      allocate(Abs_CH4(545,2))
      allocate(Abs_CO(20,2))
      allocate(Abs_N2O(411,2))
      allocate(Abs_CO2(1126,2))
      allocate(Abs_N2(72,2))
      allocate(Abs_O4(1434,2))
      allocate(Abs_HNO3(141,3))
      allocate(Abs_NO(21,2))
      allocate(Abs_NO2(767,3))
      allocate(Abs_NO3(304,3))
      allocate(Abs_SO2U(261,3))
      allocate(Abs_SO2I(10,2))
      allocate(Abs_O3UV(932,4))
      allocate(Abs_O3IR(307,2))
      allocate(Abs_NH3(421,2))
      allocate(Abs_BrO(177,2))
      allocate(Abs_CH2O(241,3))
      allocate(Abs_HNO2(193,2))
      allocate(Abs_ClNO(273,4))
      allocate(mem_pico(151))
      allocate(mem_micro(151))

      call SetSmartsConstants()
      call ReadSpctrmFile(0, ESC, m_wvln, vH0)
      call CalcAcWater(m_wvln, abW)
      call CalcAcIce(m_wvln, abI)
      call CalcAcSnow(m_wvln, abN)
      call CalcAcGrayIce(m_wvln, abE)
      call ReadAlbedoFile("FineSnow", m_wvln, Abd_fsnow)
      call ReadAlbedoFile("MeltSnow", m_wvln, Abd_msnow)
      call ReadGasFile("Abs_H2O", Abs_H2O)
      call ReadGasFile("Abs_O2", Abs_O2)
      call ReadGasFile("Abs_CH4", Abs_CH4)
      call ReadGasFile("Abs_CO", Abs_CO)
      call ReadGasFile("Abs_N2O", Abs_N2O)
      call ReadGasFile("Abs_CO2", Abs_CO2)
      call ReadGasFile("Abs_N2", Abs_N2)
      call ReadGasFile("Abs_O4", Abs_O4)
      call ReadGasFile("Abs_HNO3", Abs_HNO3)
      call ReadGasFile("Abs_NO", Abs_NO)
      call ReadGasFile("Abs_NO2", Abs_NO2)
      call ReadGasFile("Abs_NO3", Abs_NO3)
      call ReadGasFile("Abs_SO2U", Abs_SO2U)
      call ReadGasFile("Abs_SO2I", Abs_SO2I)
      call ReadGasFile("Abs_O3UV", Abs_O3UV)
      call ReadGasFile("Abs_O3IR", Abs_O3IR)
      call ReadGasFile("Abs_NH3", Abs_NH3)
      call ReadGasFile("Abs_BrO", Abs_BrO)
      call ReadGasFile("Abs_CH2O", Abs_CH2O)
      call ReadGasFile("Abs_HNO2", Abs_HNO2)
      call ReadGasFile("Abs_ClNO", Abs_ClNO)
   end subroutine

   subroutine DestructSmartsModule()
      implicit none

      deallocate(Zref50)
      deallocate(Zref41)
      deallocate(O3a0)
      deallocate(O3a1)
      deallocate(O3Min)
      deallocate(O3Max)
      deallocate(TO3Min)
      deallocate(TO3Max)
      deallocate(Tref)
      deallocate(Pref)
      deallocate(RHref)
      deallocate(O3ref)
      deallocate(Wref)
      deallocate(TO3)
      deallocate(fgphot)
      deallocate(frdif)
      deallocate(fbphot)
      deallocate(abCDOM)
      deallocate(abAP)
      deallocate(abNAP)
      deallocate(bsP)
      deallocate(vH0)
      deallocate(abW)
      deallocate(abI)
      deallocate(abN)
      deallocate(abE)
      deallocate(Abd_fsnow)
      deallocate(Abd_msnow)
      deallocate(Abs_H2O)
      deallocate(Abs_O2)
      deallocate(Abs_CH4)
      deallocate(Abs_CO)
      deallocate(Abs_N2O)
      deallocate(Abs_CO2)
      deallocate(Abs_N2)
      deallocate(Abs_O4)
      deallocate(Abs_HNO3)
      deallocate(Abs_NO)
      deallocate(Abs_NO2)
      deallocate(Abs_NO3)
      deallocate(Abs_SO2U)
      deallocate(Abs_SO2I)
      deallocate(Abs_O3UV)
      deallocate(Abs_O3IR)
      deallocate(Abs_NH3)
      deallocate(Abs_BrO)
      deallocate(Abs_CH2O)
      deallocate(Abs_HNO2)
      deallocate(Abs_ClNO)
      deallocate(mem_pico)
      deallocate(mem_micro)
   end subroutine

   !------------------------------------------------------------------------------
   !
   !  Purpose: Make the correction of downward irradiance for cloud cover.
   !
   !------------------------------------------------------------------------------
   subroutine CorrIrradianceByCloud(zenith, tcc, ofgphot, ofrdif)
      implicit none
      real(r8), intent(in) :: zenith         ! surface zenith angle
      real(r8), intent(in) :: tcc            ! total cloud cover
      real(r8), intent(inout) :: ofgphot(:)  ! total radiation spectrum
      real(r8), intent(inout) :: ofrdif(:)   ! diffuse fraction
      real(r8), parameter :: tau_ovc = 0.37
      real(r8), parameter :: alpha = 2.1
      real(r8) :: oktas, tau
      real(r8) :: duva, duvb, dpar

      ! no correction after sunset and before sunrise
      if (zenith>90) then
         return
      end if
      ! convert cloud cover fraction to oktas (see Boers et al., 2010, JGR)
      if (tcc<=0.01) then
         oktas = 0.0
      else if (tcc>=0.99) then
         oktas = 8.0
      else
         oktas = 8.0 * tcc
      end if
      ! correct radiation on cloudy days
      tau = 1.0 - (1.0 - tau_ovc)*((oktas/8.0)**alpha)
      ofgphot = tau * ofgphot
      ! calculate diffuse fraction increment for UV (Grant & Gao, 2003; JGR)
      ! the diffuse fraction for PAR and NIR is the same as UVA
      duva = -0.031 + 1.3 * exp( -0.5*(((oktas-17)/6.3)**2 + &
         ((zenith-28)/30.0)**2) )
      duvb = -0.009 + 3.73 * exp( -0.5*(((oktas-28)/9.4)**2 + &
         ((zenith-26)/26.0)**2) )
      dpar = duva
      if (NINT(oktas)==8) then
         ofrdif = 0.95 
      else
         where (m_wvln<=320)
            ofrdif = ofrdif + duvb
         elsewhere (m_wvln<=400)
            ofrdif = ofrdif + duva
         elsewhere
            ofrdif = ofrdif + dpar
         end where
         ofrdif = max(min(1.0, ofrdif), 0.0)
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   !  Purpose: Correct belowwater downward irradiance for snow reflection.
   !
   !------------------------------------------------------------------------------
   subroutine CorrIrradianceForSnow(zenith, temp, fgphot, ofrdif, zcos)
      implicit none
      real(r8), intent(in) :: zenith         ! abovewater zenith angle
      real(r8), intent(in) :: temp           ! air temperature
      real(r8), intent(inout) :: fgphot(:)   ! downward irradiance
      real(r8), intent(in) :: ofrdif(:)      ! diffuse fraction
      real(r8), intent(out) :: zcos          ! cosine of underwater zenith angle

      zcos = cos(zenith*Pi/180.0)
      if (temp<=T0 .and. m_radPars%season==0) then
         fgphot = (1.0 - (0.939*ofrdif + (1.0-0.176*zcos)/0.94* &
               (1.0-ofrdif)) * Abd_fsnow) * fgphot
      else
         if (zcos>1.0d-6) then
            fgphot = (1.0 - (1.167*ofrdif + (1.0-zcos*log(1.0+1.0/zcos))/ &
                  0.35*(1.0-ofrdif)) * Abd_msnow) * fgphot
         else
            fgphot = (1.0 - (1.167*ofrdif + (1.0-ofrdif)/0.35) * &
                  Abd_msnow) * fgphot
         end if
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   !  Purpose: Correct belowwater downward irradiance for gray ice reflection.
   !
   !------------------------------------------------------------------------------
   subroutine CorrIrradianceForGrayIce(zenith, rfrIndx, fgphot, ofrdif, zcos)
      implicit none
      real(r8), intent(in) :: zenith         ! abovewater zenith angle
      real(r8), intent(in) :: rfrIndx        ! refractive index
      real(r8), intent(inout) :: fgphot(:)   ! downward irradiance
      real(r8), intent(in) :: ofrdif(:)      ! diffuse fraction
      real(r8), intent(out) :: zcos          ! cosine of underwater zenith angle
      real(r8), parameter :: txd = 0.934     ! see Koehler et al., 2014
      real(r8) :: tha, thw, txb

      tha = zenith * Pi / 180.0
      thw = asin(sin(tha)/rfrIndx)
      zcos = cos(thw)
      ! Fresnel's equation for beam radiation reflectance
      txb = 1.0 - 0.5 * ( (sin(tha-thw)/sin(tha+thw))**2 + &
            (tan(tha-thw)/tan(tha+thw))**2 )
      ! radiation after transmitting through the gray ice-ice interface
      fgphot = (1.0 - (txd*ofrdif + txb*(1.0-ofrdif))*Alphae) &
            * fgphot
   end subroutine

   !------------------------------------------------------------------------------
   !
   !  Purpose: Correct belowwater downward irradiance for reflection.
   !
   !------------------------------------------------------------------------------
   subroutine CorrIrradianceByReflection(zenith, rfrIndx, fgphot, ofrdif, zcos)
      implicit none
      real(r8), intent(in) :: zenith         ! abovewater zenith angle
      real(r8), intent(in) :: rfrIndx        ! refractive index 
      real(r8), intent(inout) :: fgphot(:)   ! downward irradiance
      real(r8), intent(in) :: ofrdif(:)      ! diffuse fraction
      real(r8), intent(out) :: zcos          ! cosine of underwater zenith angle
      real(r8), parameter :: txd = 0.934     ! see Koehler et al., 2014
      real(r8) :: tha, thw, txb

      tha = zenith * Pi / 180.0
      thw = asin(sin(tha)/rfrIndx)
      zcos = cos(thw)
      ! Fresnel's equation for beam radiation reflectance
      txb = 1.0 - 0.5 * ( (sin(tha-thw)/sin(tha+thw))**2 + &
            (tan(tha-thw)/tan(tha+thw))**2 )
      ! radiation after transmitting through the water-air interface
      fgphot = (txd * ofrdif + txb * (1.0 - ofrdif)) * fgphot
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: The main algorithm to calculate transfer radiation.
   !
   !------------------------------------------------------------------------------
   subroutine CalcClearSkyIrradiance(hour, fgphot, ofrdif, zenith)
      implicit none
      real(r8), intent(in) :: hour        ! local hour
      real(r8), intent(out) :: fgphot(:)  ! global radiation (mol/m2/s/nm)
      real(r8), intent(out) :: ofrdif(:)  ! diffuse fraction
      real(r8), intent(out) :: zenith     ! zenith angle (decimal)
      real(r8), parameter :: epsiln = 1d-3
      real(r8), parameter :: pinb = 3.14159265d+0
      real(r8), parameter :: PHOT = 5.03411762d+24
      real(r8), parameter :: evolt = 1.602176463d-19
      real(r8), parameter :: Avogad = 6.02214199d+23
      real(r8) :: RHC, RHC2, wref, TK, TT0
      real(r8) :: TK1, EVS, EV, W, ROV, TT, HV
      real(r8) :: Tempa, Tmin, Tmax, Ozmin, Ozmax
      real(r8) :: pp0, qp, TO3ini, Tempo, Tempn
      real(r8) :: UOC, ApO3, RH2, pi2, pix4, RPD
      real(r8) :: xrh, xrh2, xrh3, alpha1, alpha2, alpha
      real(r8) :: beta, alf1, alf2, t550, RANGE1, Visi
      real(r8) :: cosx, Tilt, Rhog, Wazim, Suncor, SolarC
      real(r8) :: BETAC, EPSIRP
      real(r8) :: dec, Zenit, Azim, Julian, Radius, EOT
      real(r8) :: HourUT, ZR, ZCOS, AmH2O, AmO3, AmNO2
      real(r8) :: AmNO, AmN2, AmCO2, AmCO, AmCH4, AmN2O
      real(r8) :: AmO2, AmNH3, AmSO2, AmHNO3, AmAER, AmCH2O
      real(r8) :: AmClNO, AmHNO2, AmNO3, AmBrO, AmPOL, Amdif
      real(r8) :: ELEV, ZSIN, EPSIR, SIGMAR, EXPR
      real(r8) :: FHTcy, FHTcx, FHTcz, FHTdx, AbCO2, AbO2
      real(r8) :: AbO4, AbBrO, AbNO3, AbCH2O, AbHNO2, ABClNO
      real(r8) :: AbO3, AbCH4, AbCO, AbN2O, AbN2, AbHNO3, AbNO2
      real(r8) :: AbSO2, AbNO, AbNH3, pp0x2, pln, pln2, pln3
      real(r8) :: tau550, t5, t52, t53, t54, t505, t515, t503
      real(r8) :: tau5, t5032, t5033, t5034, t5035, t56, t562
      real(r8) :: Ama1, Ama11, Ama12, Ama15, Ama05, Ama2, Ama3
      real(r8) :: Ama4, Ama5, tra0, tra1, tra2, trb0, trb1
      real(r8) :: trb2, tzb0, tzb1, tzb2, fdifa0, fdifa1
      real(r8) :: fdifb0, fdifb1, fdifb2, tzc1b0, tzc2b0
      real(r8) :: tzc3b0, tzc4b0, tzd1b1, tzd2b1, tzd3b1, tzd4b1
      real(r8) :: tze1b2, tze2b2, tze3b2, tze4b2, alba00
      real(r8) :: alba01, alba0, alba1, alba2, albb0, albb1
      real(r8) :: albb2, albc0, albc1, albc2, albc3, albc4
      real(r8) :: Amo32, Amo33, FHTa0, FHTa1, FHTb0, FHTb1
      real(r8) :: FHTc0, FHTc1, FHTd0, FHTd1, FHTe0, FHTe1
      real(r8) :: FHTf0, FHTf1
      real(r8) :: ra00, ra01, ra10, ra11, ra12, ra02, ra13
      real(r8) :: ra03, ra14, ra04, ra15, ra16, ra05, ra17
      real(r8) :: ra06, ra18, ESCC, Scor, Scor2, wvold 
      real(r8) :: wvl, wvl2, wvl3, wvl4, TAURL, TR 
      real(r8) :: TH2o, TH2op, Tauw, wvlw, AW, bwa0, bwa1
      real(r8) :: bwa2, bma0, bma1, bma2, bmwa0, bmwa1, bmwa2
      real(r8) :: bpa1, bpa2, Bw, w0, ww0, ww02, Bm, Bmp
      real(r8) :: xamw, xamw1, xamw11, xamw12, xamw15, xamw25
      real(r8) :: Bmx(2), Bmw, Bmwp, wfw0, yamw, yamw1, yamw12
      real(r8) :: Bmwx(2), Bp, wamw, qp1, pp01, pp02, qp2 
      real(r8) :: wvln, wamp, Tautrc, Taumix, TO2, TO2P
      real(r8) :: wvlo, AO2, TauO2, TCH4, TCH4P
      real(r8) :: ACH4, TCO, TCOP, ACO, TN2O, TN2OP, AN2O
      real(r8) :: TCO2, TCO2P, tauco2, ACO2, TN2, TN2P, AN2
      real(r8) :: TO4, TO4P, xsO4, AO4, THNO3, THNO3P, xsHNO3
      real(r8) :: athno3, TNO2, TNO2P, xsNO2, atno2, xsNO3
      real(r8) :: atno3, TNO3, TNO3P, TNO, TNOP, ANO
      real(r8) :: TSO2, TSO2P, xsSO2, atso2, ASO2, TxO3
      real(r8) :: xo3, xso3, a0o3, a1o3, TO3, TAUZ3, TrO3, AO3
      real(r8) :: TNH3, TNH3P, ANH3, TBrO, TBrOP, xsBrO
      real(r8) :: ABrO, TCH2O, TCH2OP, xsCH2O, atCH2O
      real(r8) :: THNO2, THNO2P, xsHNO2, TClNO, TClNOP
      real(r8) :: TCl, xsClNO, a1tCl, a2tCl, AClNO
      real(r8) :: Tmixd, TmixdP, Trace, TraceP, TAUA, Tauaa
      real(r8) :: TAUAS, TAA, TAAP, TAS, TAT, OMEGL, BQ
      real(r8) :: TSCAT, TABS0, TABS0P, TABS, TDIR, DIR
      real(r8) :: DIRH, H0, H0H, FHTO, FHT1, Fda00, Fdazt
      real(r8) :: ssaro, Taurf, ra0, ra1, FHT1x, HT, HTa
      real(r8) :: GG, ALG, AFS, BFS, FA1, FA1P, DRAY, FR, FRP
      real(r8) :: DAER, Fdifz, Fdiftz, DIF0, Glob0, TRP, TAUAP
      real(r8) :: TASP, TTp5, GAMOZ, Rhob0, Rhob, Rhod, Rhor
      real(r8) :: Rhoa, Rhos, Roro, Dgrnd, Rho, Fatau, Upward
      real(r8) :: DIF, GLOB, WPHT, Rhox, AmR, SolarH
      real(r8) :: PFG, PFB, PFD, cosZ, qCO2
      real(r8) :: BP0, BP1, BP2, BP3, BQ0, BQ1, BQ2
      real(r8) :: AG0, AG1, AG2, AG3, AG4
      integer :: Iwarn(9), IAER, IREF, DayUT, IWVL, iabs_indx
      integer :: iband, ifitm, ifitw, ifitmw, im, ITURB
      integer :: isptm, DayoYr

      ! check whether there is sunshine
      HourUT = hour - m_radPars%Longit/15.0
      DayUT = m_radPars%day
      if (HourUT>24.) then
         HourUT = HourUT - 24.
         DayUT = DayUT + 1
      end if
      if (HourUT<0.0) then
         HourUT = HourUT + 24.
         DayUT = DayUT - 1
      end if

      pi2 = pinb/2.
      pix4 = pinb*4.
      RPD = pinb/180.

      pp0 = m_radPars%spr/1013.25   ! mb -> bar
      qp = 1.0 - pp0
      qCO2 = m_radPars%qCO2      ! ppm

      wref = 1.4164
      TK = m_radPars%tref 
      TT0 = Tk/T0
      IREF = 11
      IAER = 1
      ITURB = 5
      Iwarn = 0

      CALL SunPSA(HourUT, m_radPars%Latit, m_radPars%Longit, dec, &
                  Zenit, Azim, Julian, Radius, EOT, m_radPars%spr, TK, &
                  m_radPars%year, m_radPars%month, DayUT)
      zenith = Zenit
      if (zenith>90) then
         fgphot = 0.0_r8
         ofrdif = 1.0_r8
         return
      end if
      AmR = AMZ(Zenit)

      ! Estimate ozone temperature at sea level
      if(m_radPars%season==0) then
         TO3ini = 219.25
      else
         TO3ini = 230.87
      end if

      ! Converts Temperature at given altitude to ozone temperature
      Call Ozon(lake_info%zalt,TK,Tempa,Tmin,Tmax,Ozmin,Ozmax)

      !if (Tempa>=Tmin) then
         if (Tempa>Tmax) then
            Iwarn(6) = 1
            Tempa = Tmax
         end if
         
         if (Tempa<Tmin) then
            Iwarn(5) = 1
            Tempa = Tmin
         end if
         
         ! SATURATION VAPOR PRESSURE FROM GUEYMARD (J. Appl. Met. 1993)
         TK1 = TK/100.
         EVS = EXP(22.329699-49.140396/TK1-10.921853/TK1/TK1-.39015156*TK1)
         EV = EVS * m_radPars%RH / 100. 

         ! W=f(T,RH) USING EMPIRICAL MODEL OF GUEYMARD (SOLAR ENERGY 1994)
         ROV = 216.7*EV/TK
         TT = m_radPars%tair/T0
         HV = 0.4976 + 1.5265*TT + EXP(13.6897*TT-14.9188*TT**3)
         HV = max(min(HV,8.0_r8),1.0_r8)
         W = 0.1*HV*ROV
         if (W>12) then
            call Endrun(lake_info%id, 'Precipitable water W is above ' // &
               'the allowed maximum value of 12 cm') 
         end if
         TEMPO = TEMPA
         TEMPN = TEMPA
         
         UOC = m_radPars%AbO3
         ApO3 = 0.0

         if(lake_info%zalt>6) then
            IAER = 4
            Iwarn(7) = 1
         end if
         if (m_radPars%RH>100.) then
            RH2 = 50.0
         else
            RH2 = m_radPars%RH
         end if

         ! Angstrom's alpha =f(RH) for Shettle & Fenn aerosols
         xrh = COS(0.9*RH2*RPD)
         xrh2 = xrh*xrh
         xrh3= xrh**3
         alpha1 = (C1(IAER)+C2(IAER)*xrh+C3(IAER)*xrh2+C4(IAER)*xrh3)/ &
                  (1.+C5(IAER)*xrh+C6(IAER)*xrh2)
         alpha2 = (D1(IAER)+D2(IAER)*xrh+D3(IAER)*xrh2+D4(IAER)*xrh3)/ &
                  (1.+D5(IAER)*xrh+D6(IAER)*xrh2)

         RHC = max(50., RH2)
         RHC2 = RHC**2

         ! COEFFICIENTS FOR OMEGL (SINGLE SCATTERING ALBEDO)
         BP0 = BP00(IAER)+BP01(IAER)*RHC+BP02(IAER)*RHC2
         BP1 = BP10(IAER)+BP11(IAER)*RHC+BP12(IAER)*RHC2
         BP2 = BP20(IAER)+BP21(IAER)*RHC+BP22(IAER)*RHC2
         BP3 = BP30(IAER)+BP31(IAER)*RHC+BP32(IAER)*RHC2
         BQ0 = BQ00(IAER)+BQ01(IAER)*RHC+BQ02(IAER)*RHC2
         BQ1 = exp(BQ10(IAER)+BQ11(IAER)*RHC+BQ12(IAER)*RHC2)
         BQ2 = BQ20(IAER)+BQ21(IAER)*RHC+BQ22(IAER)*RHC2

         ! COEFFICIENTS FOR GG (ASYMMETRY FACTOR)
         AG0 = AG00(IAER)+AG01(IAER)*RHC+AG02(IAER)*RHC2
         AG1 = AG10(IAER)+AG11(IAER)*RHC+AG12(IAER)*RHC2
         AG2 = AG20(IAER)+AG21(IAER)*RHC+AG22(IAER)*RHC2
         AG3 = AG30(IAER)+AG31(IAER)*RHC+AG32(IAER)*RHC2
         AG4 = AG40(IAER)+AG41(IAER)*RHC+AG42(IAER)*RHC2

         ! Aerosol Turbidity
         tau550 = m_radPars%tau550
         tau5 = tau550*(1.1**alpha2)
         if (tau5<1d-4) then
            beta = tau5*(0.5**alpha2)
         else
            Call ALFA(m_radPars%season,IAER,ITURB,IREF,alpha1,alpha2,tau5, &
                     beta,alf1,alf2,t550,0)
            alpha1 = alf1
            alpha2 = alf2 
         end if
         if (tau550>5.0) then
            call Endrun(lake_info%id, 'Input value for turbidity is too large')
         else if (tau550<0) then
            call Endrun(lake_info%id, 'Input value for turbidity is too small')
         end if

         ! FUNCTION RANGE=F(BETA,ALPHA) FROM NEW FIT
         RANGE1 = 999.0
         Visi = 764.9
         beta = (0.5**1.33669)*tau5
         if (tau5>=0.001) then
            Call VISTAU(m_radPars%season,Range1,tau550,0)
            Visi = RANGE1/1.306
         end if
      !else
      !   Iwarn(5) = 1
      !   Tempa = Tmin

      !   RANGE1 = Min(RANGE1, 999.0)
      !   Visi = RANGE1/1.306
      
      !   Call VISTAU(m_radPars%season,RANGE1,tau550,1)
      !   tau5 = tau550*(1.1**alpha2)
      !   beta = (0.55**alpha2)*tau550
      !   Call ALFA(m_radPars%season,IAER,ITURB,IREF,alpha1,alpha2,tau5, &
      !            beta,alf1,alf2,tau550,1)
      !   alpha1 = alf1
      !   alpha2 = alf2
      !end if

      alpha = (alpha1+alpha2)/2.

      if (Iwarn(9)==1) then
         print *, 'Receiver is at more than 6 km above sea level, ' // &
                  'hence the aerosol optical depth has been fixed to ' // &
                  'a default value, dependent only on altitude.'
      end if

      ! Albedo will be counted in boundary_mod.f90
      Rhox = 0.0
      
      ! no tilt
      Tilt = 0.0
      Rhog = 0.0
      Wazim = 0.0

      ! PRELIMINARY CALCULATIONS FOR SPECTRAL FUNCTIONS
      BETAC = beta*(2.**(alpha2-alpha1))
      EPSIRP = 0.1686

      ! Solar Time
      !SolarH = hour + EOT/60.

      ZR = Zenit*RPD
      ZCOS = cos(ZR)
      ZSIN = sin(ZR)

      AmH2O = 1./(ZCOS+.10648*(Zenit**0.11423)*(93.781-Zenit)**(-1.9203))
      AmO3 = 1./(ZCOS+1.0651*(Zenit**.6379)*(101.8-Zenit)**(-2.2694))
      AmNO2 = 1./(ZCOS+1.1212*(Zenit**1.6132)*(111.55-Zenit)**(-3.2629))
      AmNO = 1./(ZCOS+.77738*(Zenit**.11075)*(100.34-Zenit)**(-1.5794))
      AmN2 = 1./(ZCOS+.38155*(Zenit**8.871d-5)*(95.195-Zenit)**(-1.8053))
      AmCO2 = 1./(ZCOS+.65786*(Zenit**.064688)*(96.974-Zenit)**(-1.8083))
      AmCO = 1./(ZCOS+.505*(Zenit**.063191)*(95.899-Zenit)**(-1.917))
      AmCH4 = 1./(ZCOS+.49381*(Zenit**.35569)*(98.23-Zenit)**(-2.1616))
      AmN2O = 1./(ZCOS+.61696*(Zenit**.060787)*(96.632-Zenit)**(-1.8279))
      AmO2 = 1./(ZCOS+.65779*(Zenit**.064713)*(96.974-Zenit)**(-1.8084))
      AmNH3 = 1./(ZCOS+.32101*(Zenit**.010793)*(94.337-Zenit)**(-2.0548))
      AmSO2 = 1./(ZCOS+.63454*(Zenit**.0099198)*(95.804-Zenit)**(-2.0573))
      AmHNO3 = 1./(ZCOS+1.044*(Zenit**.78456)*(103.15-Zenit)**(-2.4794))
      AmAER = 1./(ZCOS+.16851*(Zenit**.18198)*(95.318-Zenit)**(-1.9542))
      AmCH2O = AmN2O
      AmClNO = AmNO2
      AmHNO2 = AmHNO3
      AmNO3 = AmNO2
      AmBrO = AmO3
      AmPOL = 1./(.0001569+.9998431*zcos*zcos)**.5
      Amdif = 1.732

      Zenit = min(90.,Zenit)
      ELEV = 90. - Zenit
      EPSIR = 0.17*(1.-EXP(-8.*ZCOS))
      SIGMAR = 3.65-2.3*EXP(-4.*ZCOS)
      EXPR = 0.72 + ZCOS

      FHTcy = (-11.012+12.392*AMO3)/(1.+.23644*AMO3)
      FHTcx = 3.2656*(1.-EXP(-.46464*AMO3**1.25))-.965936*FHTcy
      FHTcz = 2.*FHTcx+1.93187*FHTcy
      FHTdx = EXP(.31045+.001684*AMO3-.28549/AMO3**4)

      SUNCOR = 1./RADIUS**2

      ! Gas Abundances for total column in atm-cm
      AbCO2 = .802685 * qCO2 * pp0
      AbO2 = 1.67766d5 * pp0
      AbO4 = 1.8171d4*NLosch*NLosch*(pp0**1.7984)/(TT0**.344)
      AbBrO = 2.5d-06
      AbNO3 = .00005
      AbCH2O = .0003
      AbHNO2 = .0001
      AbClNO = .00012
      AbO3 = uoc
      pp0x2 = pp0*pp0
      pln = log(pp0)
      pln2 = pln*pln
      pln3 = pln2*pln
      AbCH4 = 1.31195*(pp0*1.1245)*(TT0**.047283)
      AbCO = .31491*(pp0**2.6105)*exp(.82546-2.9753*pp0+.88437*pp0x2)
      AbN2O = .24344*(pp0**1.1625)
      AbN2 = 3.8298*(pp0**1.8643)/(TT0**.50342)
      AbHNO3 = 3.6739d-4*(pp0**.13568)/(TT0**.0714)
      AbNO2 = 1d-4*Min(1.864+.20314*pp0,41.693*pp0)
      AbNO = 1d-4*Min(.74039+2.4154*pp0,57.314*pp0)
      AbSO2 = 1.114d-5*(pp0**.81319)*exp(.81356+3.0448*pp0x2-1.5652*pp0x2*pp0)
      AbNH3 = exp(-8.6472+2.239*pLN-2.3822*PLN2-1.4408*PLN3-.46322*pln3*pln)

      ! Turbidity dependent coefficients for diffuse algorithm
      t5 = tau550
      t52 = t5*t5
      t53 = t5*t52
      t54 = t5*t53
      t505 = t5**0.5
      t515 = t5**1.5
      t503 = max(t5-.03, 0.0)
      t5032 = t503*t503
      t5033 = t503*t5032
      t5034 = t503*t5033
      t5035 = t503**.5
      t56 = t5 - 0.6
      t562 = t56*t56

      Ama1 = min(Amaer,20.)
      Ama11 = Ama1 - 1.0
      Ama12 = Ama11*Ama11
      Ama15 = Ama11**1.5
      Ama05 = Ama11**0.5
      Ama2 = Ama1*Ama1
      Ama3 = Ama1*Ama2
      Ama4 = Ama1*Ama3
      Ama5 = Ama1**.5
      tra0 = (.79996+5.2399*t5-.45849*t52)/(1.+5.6621*t5+.67258*t52)
      tra1 = (-.25538+5.4972*t5-2.484*t52)/(1.+3.0475*t5+3.1615*t52)
      tra2 = (1.7845+8.5655*t5-2.8046*t52)/(1.+1.9911*t5+1.0593*t52)
      trb0 = (.98+4.2022*t5-.36341*t52)/(1.+4.1566*t5-.24198*t52)
      trb1 = (-2.927-3.5797*t5+19.036*t52)/(1.+28.266*t5+1.4545*t52)
      trb2 = (.017854+6.9189*t5+3.1496*t52)/(1.+7.1728*t5-.75729*t52)
      tzb0 = 1.0
      tzb1 = 0.0
      tzb2 = 0.0
      if(Zenit>=45.) then
         fdifa0 = (-81.197-42.369*Ama1+1.6385*Ama2+117.14*Ama5)/ &
                  (1.+.0048737*Ama2)
         fdifa1 = (287.15+148.6*Ama1-5.6994*Ama2-410.8*Ama5)/ &
                  (1.+.005139*Ama2)
         fdifb0 = (3.1461+4.3549*t5-.023703*t52-5.3845*t505)/ &
                  (1.+1.8418*t52)
         fdifb1 = (-.34571-.1714*t5-.0022115*t52+.55363*t505)/ &
                  (1.+.046034*t52)
         fdifb2 = (.017953+.087827*t5-.016289*t52-.065605*t505)/ &
                  (1.+3.231*t52)
      else
         fdifa0 = (-405.89+669.51*Ama1-421.72*Ama3+163.26*Ama4)/ &
                  (1.-1.7483*Ama1+.86666*Ama2)
         fdifa1 = (1375.2-2268.7*Ama1+1430.2*Ama3-553.82*Ama4)/ &
                  (1.-1.7459*Ama1+.86337*Ama2)
      end if
      if(zenit>0.0 .and. t5>0.03) then
         if (zenit<=80) then
            tzc1b0 = (1.0303-.29853*t5+3.0145*t52-5.1787*t53+.75487*t54)/ &
                     (1.+48.21*t5-11.627*t505)
            tzc2b0 = (.0087736-.023675*t5+.038433*t52-.02068*t53+.0025949*t54)/ &
                     (1.+.75681*t5-1.102*t505)
            tzc3b0 = (-1.6395+4.6978*t5-9.2494*t52+5.9966*t53-.78497*t54)/ &
                     (1.+37.007*t5-3.1734*t505)
            tzd1b1 = (.79647-.34098*t56+.23205*t562)/(1.-.6048*t56+.16386*t562)
            tzd2b1 = (-.059924+.084758*t56-.025479*t562)/(1.-.60223*t56+.12697*t562)
            tzd3b1 = (-.042678+.34131*t56-.28495*t562)/(1.+.14978*t562)
            tzd4b1 = (-.013436+.23576*t56+.016945*t562)/(1.-.61716*t56+.13362*t562)
            tze1b2 = -1.4472+.67069*t5-.18447*t52+.031889*t53
            tze2b2 = .10793-.054098*t5+.072644*t52-.015276*t53
            tze3b2 = .4919-1.0133*t5+.35268*t52-.044266*t53
            tze4b2 = .073577-.036235*t5+.066935*t52-.014448*t53
            if (tau550<=1.75) then
               tze1b2 = (-1.7064-3.3834*t5-.47597*t54+2.6748*t505)/ &
                        (1.+.92521*t52)
               tze2b2 = (.16239+.5601*t5+.13857*t54-.39887*t505)/(1.+1.1706*t52)
               tze3b2 = (-.14892-.084698*t5+.020901*t54+.19723*t505)/(1.+.074613*t52)
               tze4b2 = (.11187+.36411*t5+.19632*t54-.168*t505)/(1.+.42997*t52)
               if (tau550<=0.6) then
                  tzc1b0 = (.068343-1.4788*t5+23.173*t52-50.108*t53+42.158*t54)/ &
                           (1.+14.968*t5-7.0445*t505)
                  tzc2b0 = (.010684-.10861*t5-4.3737*t52+9.6237*t53-8.1427*t54)/ &
                           (1.+70.477*t5-19.077*t505)
                  tzc3b0 = (.0055635+.15348*t5-2.3326*t52+6.547*t53-5.2644*t54)/ &
                           (1.+26.745*t5-9.9257*t505)
                  tzc4b0 = (.030341-1.0782*t5+23.098*t52-49.809*t53+44.96*t54)/ &
                           (1.+18.655*t5-7.7781*t505)
                  tzd1b1 = (2.4+24.478*t503-3.8623*t5032-10.047*t5035) &
                           /(1.+14.505*t503)
                  tzd2b1 = (-.068+.68664*t503-.40106*t5032-.26178*t5035) &
                           /(1.-1.6208*t503)
                  tzd3b1 = (-1.2-15.514*t503+10.004*t5032+8.1483*t5035) &
                           /(1.+19.234*t503)
                  tzd4b1 = (.5-.12601*t503+7.0083*t5032-3.805*t5035) &
                           /(1.+34.177*t503)
               end if
            end if
            tzb0 = 1.+tzc1b0*Ama11+tzc2b0*Ama12+tzc3b0*Ama15
            tzb1 = (tzd1b1*Ama11+tzd2b1*Ama12+tzd3b1*Ama05)/(1.+tzd4b1*Ama11)
            tzb2 = (tze1b2*Ama11+tze2b2*Ama12+tze3b2*Ama05)/(1.+tze4b2*Ama12)
            if(tau550<=1.75) then
               tzb2 = (tze1b2*Ama11+tze2b2*Ama12+tze3b2*Ama05)/(1.+tze4b2*Ama11)
               if (tau550<=0.6) then
                  tzb0 = (1.+tzc1b0*Ama11+tzc2b0*Ama12+tzc3b0*Ama05)/(1.+tzc4b0*Ama11)
                  tzb1 = (tzd1b1*Ama11+tzd2b1*Ama12+tzd3b1*Ama05)/(1.+tzd4b1*Ama11)
               end if
            end if
         else
            tzc1b0 = (-.21837+.42655*t5-.25804*t52+.053587*t53-.0034933*t54)/ &
                     (1.-1.0045*t5)
            tzc2b0 = (.007183-.011421*t5+.0048386*t52-.00060429*t53 &
                     -1.1538d-5*t54)/(1.-1.0044*t5)
            tzc3b0 = (.41391-.94615*t5+.66482*t52-.14276*t53+.010187*t54)/ &
                     (1.-1.0012*t5)
            tzd1b1 = (-1.28-10.727*t5-.84523*t52+1.6956*t53-.26123*t54)/ &
                     (1.+6.48*t5)
            tzd2b1 = (-.014423+5.2179*t5+.10566*t52-.52195*t53+.08956*t54)/ &
                     (1.+51.341*t5)
            tzd3b1 = 3.3705+.91411*t5-1.1374*t52+.27451*t53-.020812*t54
            tze1b2 = (4.8107+9.1856*t5-2.4024*t52+.29673*t53-10.278*t505)/ &
                     (1.+.16684*t52)
            tze2b2 = (-.18328-.39604*t5+.098351*t52-.011915*t53+.43926*t505)/ &
                     (1.+.17865*t52)
            tze3b2 = (-9.9946-14.963*t5+4.0066*t52-.50372*t53+18.124*t505)/ &
                     (1.+.13084*t52)
            if(tau550<=0.6) then
               tzc1b0 = (-.076486-2.9836*t5-3.6388*t52+5.2983*t53+1.0036*t505)/ &
                        (1.+31.865*t5-11.358*t505)
               tzc2b0 = (.0075+.060214*t503+.20831*t5032-.40488*t5033+ &
                        .20894*t5034)/(1.+15.819*t503+2.2681*t5035)
               tzc3b0 = (.072803+7.6593*t5+3.2585*t52-11.177*t53-1.9115*t505)/ &
                        (1.+37.274*t5-12.44*t505)
               tze1b2 = (-.09282+3.6079*t5+72.469*t52-42.168*t53-.81503*t505)/ &
                        (1.+34.684*t5-11.93*t505)
               tze2b2 = (-.051779+.36296*t5-2.8752*t52)/(1.+52.058*t52)
               tze3b2 = (-.94044-51.232*t5-120.6*t52+88.994*t53+15.858*t505)/ &
                        (1.+35.032*t5-11.999*t505)
               if (tau550<=0.1) then
                  tzd1b1 = (-3.6546-15.903*t5+15.443*t505)/(1.-7.4089*t5)
                  tzd2b1 = (.2838+3.3062*t5-1.3909*t505)/(1.+11.041*t5)
                  tzd3b1 = (8.0948+25.437*t5-31.606*t505)/(1.-8.1072*t5)
               end if
            end if
            tzb0 = 1.+tzc1b0*Ama11+tzc2b0*Ama12+tzc3b0*Ama05
            tzb1 = tzd1b1*Ama11+tzd2b1*Ama12+tzd3b1*Ama05
            tzb2 = tze1b2*Ama11+tze2b2*Ama12+tze3b2*Ama05
         end if
      end if

      alba00 = 21.712-74.917*t5
      alba01 = -79.895+263.33*t5
      alba0 = (-13.126+331.3*t5+48.496*t52)/(1.+4.5503*t515)
      alba1 = (113.53-2152.3*t5-280.*t52)/(1.+4.7907*t515)
      alba2 = (-260.19+3516.8*t5+438.98*t52)/(1.+5.1594*t515)
      albb0 = (-169.74+366.27*t5-1497.4*t52+417.66*t53)/(1.+211.44*t52)
      albb1 = (695.55-2080.1*t5+5470.9*t52-1385.8*t53)/(1.+178.71*t52)
      albb2 = (-724.21+2368.7*t5-5331.3*t52+1335.4*t53)/(1.+158.89*t52)
      albc0 = (.030827+2.1215*t5+.23068*t52)/(1.+.84396*t5)
      albc1 = (-.03264-2.1233*t5-.29271*t52)/(1.+1.1121*t5)
      albc2 = (.015472+.7404*t5+.12473*t52)/(1.+1.1626*t5)
      albc3 = (-.002198-.084679*t5-.015984*t52)/(1.+1.1628*t5)
      albc4 = (1.019+1.6635*t5+.47646*t52)/(1.+80.333*t5)

      ! Ozone-dependent coefficients for diffuse algorithm
      Amo32 = Amo3*Amo3
      Amo33 = Amo32*Amo3
      FHTa0 = (-.47724+2.4591*Amo3-1.5553*Amo32+.36215*Amo33)/ &
               (1.-1.5041*Amo3+.58992*Amo32)
      FHTa1 = (.018787-.010397*Amo3-.017601*Amo32+.0052716*Amo33)/ &
               (1.-1.4903*Amo3+.56957*Amo32)
      FHTb0 = (15.376-4.3125*Amo3+.83276*Amo32-.040636*Amo33)/ &
               (1.-4.6253*Amo3+2.8311*Amo32)
      FHTb1 = (-.13145+.035853*Amo3-.0076174*Amo32+.00038038*Amo33)/ &
               (1.-1.9406*Amo3+.93498*Amo32)
      FHTc0 = (-62.912+49.558*Amo3)/(1.-.027447*Amo3)
      FHTc1 = (3.3372-2.4508*Amo3)/(1.-.035041*Amo3)
      FHTd0 = (30.875-4.5802*Amo3+.50428*Amo32)/ &
               (1.-2.5567*Amo3+1.8116*Amo32)
      FHTd1 = (-.23641+.020167*Amo3-.004264*Amo32)/ &
               (1.-1.3537*Amo3+.59902*Amo32)
      FHTe0 = (2.5689-4.824*Amo3+2.07*Amo32)/ &
               (1.-1.1905*Amo3+.3643*Amo32)
      FHTe1 = (-.066605+.21397*Amo3-.10567*Amo32)/ &
               (1.-1.2116*Amo3+.37579*Amo32)
      FHTf0 = (2.2348-.44161*Amo3+.3772*Amo32-.020273*Amo33)/ &
               (1.-.97556*Amo3+.28047*Amo32)
      FHTf1 = (-.012594-.012008*Amo3-.015456*Amo32+.00085387*Amo33)/ &
               (1.-.96938*Amo3+.26425*Amo32)
      ra00 = (1.8973-2.8609*Amo3+1.4498*Amo32-.18485*Amo33)/ &
               (1.-.95212*Amo3+.24444*Amo32)
      ra01 = (.35236-.2446*Amo3+.24659*Amo32-.013065*Amo33)/ &
               (1.-.88635*Amo3+.22055*Amo32)
      ra10 = (-.58215+.31643*Amo3-.023724*Amo32+.00068713*Amo33)/ &
               (1.-.1444*Amo3-.11746*Amo32)
      ra11 = (18.015-121.17*Amo3+81.105*Amo32-13.644*Amo33)/ &
               (1.-7.7047*Amo3)
      ra12 = (.092338+1.1519*Amo3-2.3328*Amo32+1.1325*Amo33)/ &
               (1.-1.4379*Amo3+.7014*Amo32)
      ra02 = (1.4738-.90914*Amo3+.14322*Amo32)/ &
               (1.-.30469*Amo3+.027331*Amo32)
      ra13 = (-.20733+.19451*Amo3-.029374*Amo32)/ &
               (1.-.51985*Amo3+.081935*Amo32)
      ra03 = (-8.1831+3.2169*Amo3-.18812*Amo32)/(1.+.32473*Amo3)
      ra14 = (2.1533-.57263*Amo3+.03416*Amo32)/(1.+.0027812*Amo3)
      ra04 = (1.9915-.58894*Amo3+.05611*Amo32 &
               -.0013311*Amo33)/(1.-.26987*Amo3+.02113*Amo32)
      ra15 = (1.7389-3.3761*Amo3+2.1498*Amo32-.40236*Amo33)/ &
               (1.-.71275*Amo3+.16324*Amo32)
      ra16 = (1.4282-.61201*Amo3+.10302*Amo32 &
               -.0055147*Amo33)/(1.-.6277*Amo3+.12961*Amo32)
      ra05 = (18.214+65.305*Amo3-8.6308*Amo32+1.1603*Amo33)/ &
               (1.0+69.727*Amo3-.34374*Amo32)
      ra17 = (-.12865-.18023*Amo3+.084979*Amo32-.0068614*Amo33)/ &
               (1.-.18029*Amo3+.020348*Amo32)
      ra06 = (1.0122-.18077*Amo3+.43678*Amo32)/ &
               (1.-.15562*Amo3+.47075*Amo32)
      ra18 = (1.9286-2.0616*Amo3+.24389*Amo32)/ &
               (1.+.2096*Amo3+.030102*Amo32)

      SolarC = 1366.1 
      ESCC = SUNCOR*SolarC
      Scor = ESCC/ESC
      Scor2 = SolarC/ESC

      ! START SPECTRAL CALCULATIONS
      wvold = m_wvln(1) - 0.5
      
      IWVL = 0
      do isptm = 1, NSPCTM, 1
         wvln = m_wvln(isptm)
         wvl = wvln/1000.
         wvl2 = wvl*wvl
         wvl3 = wvl*wvl2
         WVL4 = WVL*WVL3
         if (wvln>wvold) then
            H0 = vH0(isptm)*Scor
            wvold = wvln
            IWVL = IWVL+1

            ! RAYLEIGH SCATTERING FUNCTION
            TAURL = pp0/(117.3405*wvl4-1.5107*wvl2+.017535-8.7743d-4/wvl2)
            TR = exp(-AmR*TAURL)

            ! WATER VAPOR ABSORPTION
            TH2o = 1.0d00
            TH2oP = 1.d00
            Tauw = 0.d00
            
            if (wvln>=440.0 .and. W>0.0) then
               call SearchWVLM(Abs_H2O(:,1), wvln, iabs_indx) ! see line 2246 to 2249
               wvlw = Abs_H2O(iabs_indx,1) 
               AW = Abs_H2O(iabs_indx,2)
               iband = int(Abs_H2O(iabs_indx,3))
               ifitw = int(Abs_H2O(iabs_indx,4))
               bwa0 = Abs_H2O(iabs_indx,5)
               bwa1 = Abs_H2O(iabs_indx,6)
               bwa2 = Abs_H2O(iabs_indx,7)
               ifitm = int(Abs_H2O(iabs_indx,8))
               bma0 = Abs_H2O(iabs_indx,9)
               bma1 = Abs_H2O(iabs_indx,10)
               bma2 = Abs_H2O(iabs_indx,11)
               ifitmw = int(Abs_H2O(iabs_indx,12))
               bmwa0 = Abs_H2O(iabs_indx,13)
               bmwa1 = Abs_H2O(iabs_indx,14)
               bmwa2 = Abs_H2O(iabs_indx,15)
               bpa1 = Abs_H2O(iabs_indx,16)
               bpa2 = Abs_H2O(iabs_indx,17)
               if (Aw>0.0) then
                  Bw = 1.0
                  select case (iband)
                     case (2)
                        w0 = 2.92232
                     case (3)
                        w0 = 1.41642
                     case (4)
                        w0 = 0.41612
                     case (5)
                        w0 = 0.05663
                     case default
                        w0 = 4.11467
                  end select
                  ww0 = W - w0
                  ww02 = ww0*ww0 
                  if (ifitw==1) then
                     Bw = Bw/(1.+bwa2*ww0)
                  else if (ifitw==2) then
                     Bw = Bw/(1.+bwa2*ww02)
                  else if (ifitw==6) then
                     Bw = bwa0+bwa1*ww0
                  else
                     Bw = 1.+bwa0*ww0+bwa1*ww02
                  end if
                  Bm = 1.0
                  Bmp = 1.0
                  do im = 1, 2, 1
                     xamw = AmH2O
                     if(im==2) xamw = Amdif
                     xamw1 = xamw-1.
                     xamw11 = xamw1**.1
                     xamw12 = xamw1*xamw1
                     xamw15 = xamw1**.5
                     xamw25 = xamw1**.25
                     
                     if(ifitm==0) then
                        Bmx(im) = bma1*(xamw**bma2)
                     else if(ifitm==1) then
                        Bmx(im) = (1.+bma0*xamw1+bma1*xamw12)/ &
                                 (1.+bma2*xamw1)
                     else if(ifitm==2) then
                        Bmx(im) = (1.+bma0*xamw1+bma1*xamw12)/ &
                                 (1.+bma2*xamw12)
                     else if(ifitm==3) then
                        Bmx(im) = (1.+bma0*xamw1+bma1*xamw12)/ &
                                 (1.+bma2*xamw15)
                     else if(ifitm==5) then
                        Bmx(im) = (1.+bma0*xamw25)/(1.+bma2*xamw11)
                     end if
                  end do

                  Bm = Bmx(1)
                  Bmp = Bmx(2)
                  Bm = min(max(Bm,.05),7.0)
                  Bmp= min(max(Bmp,.1),6.0)
                  Bw = min(max(Bw,.05),7.0)
                  Bmw = Bm*Bw
                  Bmwp = Bmp*Bw
                  if(abs(Bw-1.)>=epsilm .and. (ifitm==0 .or. abs(Bm-1.)>=epsilm) &
                        .and. (ifitm/=0 .or. Bm<=0.968 .or. Bm>=1.0441) .and. &
                        ifitmw/=-1) then
                     Bmw = 1.0
                     wfw0 = w/w0
                     do im = 1, 2, 1
                        xamw = AmH2O
                        if(im==2) xamw = Amdif
                        yamw = xamw*wfw0
                        yamw1 = yamw - 1.0
                        yamw12 = yamw1*yamw1
                        if(ifitmw==0) then
                           Bmwx(im) = bmwa1*(yamw**bmwa2)
                        else if(ifitmw==1) then
                           Bmwx(im) = (1.+bmwa0*yamw1+bmwa1*yamw12)/(1.+bmwa2*yamw1)
                        else if(ifitmw==2) then
                           Bmwx(im) = (1.+bmwa0*yamw1+bmwa1*yamw12)/(1.+bmwa2*yamw12)
                        end if
                     end do
                     Bmw = Bmwx(1)
                     Bmwp = Bmwx(2)
                     Bmw = min(max(Bmw,.05),7.0)
                     Bmwp= min(max(Bmwp,.1),6.0)
                  end if

                  Bp = 1.0
                  wamw = W*AmH2O
                  if(abs(qp)>=1d-5) then
                     qp1 = min(.35,qp)
                     pp01 = max(.65,pp0)
                     pp02 = pp01*pp01
                     qp2 = qp1*qp1
                     select case (iband)
                        case (2)
                           Bp = 1. + .08721*qp1
                        case (3)
                           Bp = 1. - bpa1*qp1 - bpa2*qp2
                        case (4)
                           Bp = (1.-exp(-.63486+6.9149*pp01-13.853*pp02)*wamw) &
                              *(1.-bpa1*qp1-bpa2*qp2)
                        case (5)
                           Bp = (1.-wamw*exp(8.9243-18.197*pp01+2.4141*pp02)) &
                              *(1.-bpa1*qp1-bpa2*qp2)
                        case default
                           Bp = 1.0 + 0.1623*qp
                     end select
                     Bp=min(max(Bp,.3),1.7)
                  end if
                  wamw = (W*AmH2O)**.9426
                  wamp = (Amdif*W)**.9426
                  
                  Tauw = Bmw*Bp*Aw*wamw
                  TH2O = exp(-Tauw)
                  TH2OP = exp(-Bmwp*Bp*Aw*wamp)
                  Tauw = Tauw/AmH2O
               end if
            end if

            ! Absorption from all other gases. Those in variable quantities 
            ! due to pollution are treated in subroutines
            Tauz3 = 0.0
            Tautrc = 0.0
            Taumix = 0.0

            ! Uniformly mixed gases
            ! 1. Oxygen (O2)
            TO2 = 1.0
            TO2P = 1.0
            if(wvln>=627.0 .and. wvln<=1581.0) then
               call SearchWVLM(Abs_O2(:,1), wvln, iabs_indx) 
               wvlo = Abs_O2(iabs_indx,1)
               AO2 = Abs_O2(iabs_indx,2)
               if (AO2>0.0) then
                  Tauo2 = AO2*AbO2
                  TO2 = exp(-(Tauo2*AmO2))
                  TO2P = exp(-(Tauo2*Amdif))
               end if
            end if

            ! 2. Methane (CH4)
            TCH4 = 1.0
            TCH4P = 1.0
            if (AbCH4>0.0 .and. wvln>=1617.0) then
               call SearchWVLM(Abs_CH4(:,1), wvln, iabs_indx) 
               wvlo = Abs_CH4(iabs_indx,1)
               ACH4 = Abs_CH4(iabs_indx,2)
               if(ACH4>0.0) then
                  Call GSCH4(ACH4,AbCH4,AmCH4,TCH4,TCH4P,Amdif)
               end if
            end if

            ! 3. Carbon Monoxide (CO)
            TCO = 1.0
            TCOP = 1.0
            if (AbCO>0 .and. wvln>=2310 .and. wvln<=2405) then
               call SearchWVLM(Abs_CO(:,1), wvln, iabs_indx)
               wvlo = Abs_CO(iabs_indx,1)
               ACO = Abs_CO(iabs_indx,2)
               Call GSCO(ACO,AbCO,AmCO,TCO,TCOP,Amdif)
            end if

            ! 4. Nitrous oxide (N2O)
            TN2O = 1.0
            TN2OP = 1.0
            if (wvln>=1950.0) then
               call SearchWVLM(Abs_N2O(:,1), wvln, iabs_indx)
               wvlo = Abs_N2O(iabs_indx,1)
               AN2O = Abs_N2O(iabs_indx,2)
               if (AN2O>0.0) then
                  TN2O = exp(-(AN2O*AbN2O*AmN2O))
                  TN2OP = exp(-(AN2O*AbN2O*Amdif))
               end if
            end if

            ! 5. Carbon Dioxide (CO2)
            TCO2 = 1.0
            TCO2P = 1.0
            if(wvln>=1036.0) then
               call SearchWVLM(Abs_CO2(:,1), wvln, iabs_indx)
               wvlo = Abs_CO2(iabs_indx,1)
               ACO2 = Abs_CO2(iabs_indx,2)
               if (ACO2>0.0) then
                  tauco2 = ACO2*AbCO2
                  TCO2 = exp(-(tauco2*AmCO2))
                  TCO2P = exp(-(tauco2*Amdif))
               end if
            end if

            ! 6. Nitrogen (N2)
            TN2 = 1.0
            TN2P = 1.0
            if(wvln>=3645.0) then
               call SearchWVLM(Abs_N2(:,1), wvln, iabs_indx)
               wvlo = Abs_N2(iabs_indx,1)
               AN2 = Abs_N2(iabs_indx,2)
               TN2 = exp(-(AN2*AbN2*AmN2))
               TN2P = exp(-(AN2*AbN2*Amdif))
            end if

            ! 7. Oxygen-oxygen (O2-O2 or O4) 
            ! [Collision-induced absorption; also includes O2-N2]
            TO4 = 1.0
            TO4P = 1.0
            if(wvln<=1593.0) then
               call SearchWVLM(Abs_O4(:,1), wvln, iabs_indx)
               wvlo = Abs_O4(iabs_indx,1)
               xsO4 = Abs_O4(iabs_indx,2)
               if(xsO4>0.0) then
                  AO4 = xsO4*1d-46
                  TO4 = exp(-AO4*AbO4*AmO2)
                  TO4P = exp(-AO4*AbO4*Amdif)
               end if
            end if

            ! Misc. Trace gases
            ! 1. Nitric acid (HNO3)
            THNO3 = 1.0
            THNO3P = 1.0
            if(AbHNO3>0.0 .and. wvln<=350.0) then
               call SearchWVLM(Abs_HNO3(:,1), wvln, iabs_indx)
               wvlo = Abs_HNO3(iabs_indx,1)
               xsHNO3 = Abs_HNO3(iabs_indx,2)
               athno3 = Abs_HNO3(iabs_indx,3)
               Call GSHNO3(234.2d0,xsHNO3,athno3,AbHNO3,AmHNO3,THNO3, &
                           THNO3P,Amdif)
            end if

            ! 2. Nitrogen dioxide (NO2)
            TNO2 = 1.0
            TNO2P = 1.0
            if(AbNO2>0.0 .and. wvln<=926.0) then
               call SearchWVLM(Abs_NO2(:,1), wvln, iabs_indx)
               wvlo = Abs_NO2(iabs_indx,1)
               xsNO2 = Abs_NO2(iabs_indx,2)
               atno2 = Abs_NO2(iabs_indx,3)
               Call GSNO2(TempN,xsNO2,atno2,AbNO2,AmNO2,TNO2,TNO2P,Amdif)
            end if

            ! 3. Nitrogen trioxide (NO3)
            TNO3 = 1.0
            TNO3P = 1.0
            if(AbNO3>0.0 .and. wvln>=400.0 .and. wvln<=703.0) then
               call SearchWVLM(Abs_NO3(:,1), wvln, iabs_indx)
               wvlo = Abs_NO3(iabs_indx,1)
               xsNO3 = Abs_NO3(iabs_indx,2)
               atno3 = Abs_NO3(iabs_indx,3)
               Call GSNO3(TempN,xsNO3,atno3,AbNO3,AmNO3,TNO3,TNO3P,Amdif)
            end if

            ! 4. Nitric oxide (NO)
            TNO = 1.0
            TNOP = 1.0
            if(AbNO>0.0 .and. wvln>=2645 .and. wvln<=2745) then
               call SearchWVLM(Abs_NO(:,1), wvln, iabs_indx)
               wvlo = Abs_NO(iabs_indx,1)
               ANO = Abs_NO(iabs_indx,2)
               Call GSNO(ANO,AbNO,AmNO,TNO,TNOP,Amdif)
            end if

            ! 5. Sulfur Dioxide (SO2)
            TSO2 = 1.0
            TSO2P = 1.0
            if(AbSO2>0.0) then
               ! 5a. Sulfur Dioxide (SO2) [UV band]
               if (wvln<=420.0) then
                  call SearchWVLM(Abs_SO2U(:,1), wvln, iabs_indx)
                  wvlo = Abs_SO2U(iabs_indx,1)
                  xsSO2 = Abs_SO2U(iabs_indx,2)
                  atso2 = Abs_SO2U(iabs_indx,3)
                  Call GSSO2U(247.1d0,xsSO2,atso2,AbSO2,AmSO2,TSO2, &
                              TSO2P,Amdif)
               end if

               ! 5b. Sulfur Dioxide (SO2) [IR band]
               if(wvln>=3955.0) then
                  call SearchWVLM(Abs_SO2I(:,1), wvln, iabs_indx)
                  wvlo = Abs_SO2I(iabs_indx,1)
                  ASO2 = Abs_SO2I(iabs_indx,2)
                  Call GSSO2I(ASO2,AbSO2,AmSO2,TSO2,TSO2P,Amdif)
               end if
            end if

            ! 6. Ozone (O3)
            TO3 = 1.0
            TXO3 = TO3
            TAUZ3 = 0.0
            if(AbO3>0.0) then
               ! 6a. Ozone (O3) [UV and VIS bands]
               if(wvln<=1091) then
                  call SearchWVLM(Abs_O3UV(:,1), wvln, iabs_indx)
                  wvlo = Abs_O3UV(iabs_indx,1)
                  xso3 = Abs_O3UV(iabs_indx,2)
                  a0o3 = Abs_O3UV(iabs_indx,3)
                  a1o3 = Abs_O3UV(iabs_indx,4)
                  if (wvl<345) then
                     TrO3 = 228.0
                  else
                     TrO3 = 223.0
                  end if
                  Call GSO3U(TrO3,Tempo,xso3,a0o3,a1o3,Abo3,AmO3,TAUZ3, &
                             TXO3,xo3,AO3)
                  TAUZ3 = max(TAUZ3, 0.0d+0)
               end if

               ! 6b. Ozone (O3) [IR bands]
               if(wvln>=2470.) then
                  call SearchWVLM(Abs_O3IR(:,1), wvln, iabs_indx)
                  wvlo = Abs_O3IR(iabs_indx,1)
                  AO3 = Abs_O3IR(iabs_indx,2)
                  TXO3 = exp(-AO3*AbO3*AmO3)
               end if
            end if
            TO3 = TXO3

            ! 7. Ammonia (NH3)
            TNH3 = 1.0
            TNH3P = 1.0
            if(wvln>=1900.0) then
               call SearchWVLM(Abs_NH3(:,1), wvln, iabs_indx)
               wvlo = Abs_NH3(iabs_indx,1)
               ANH3 = Abs_NH3(iabs_indx,2)
               if(ANH3>0.0) then
                  TNH3 = exp(-ANH3*AbNH3*AmNH3)
                  TNH3P = exp(-ANH3*AbNH3*Amdif)
               end if
            end if

            ! 8. Bromine monoxide (BrO)
            TBrO = 1.0
            TBrOP = 1.0
            if(wvln>=296.5 .and. wvln<=384.5) then
               call SearchWVLM(Abs_BrO(:,1), wvln, iabs_indx)
               wvlo = Abs_BrO(iabs_indx,1)
               xsBrO = Abs_BrO(iabs_indx,2)
               ABrO = xsBrO*NLosch
               TBrO = exp(-ABrO*AbBrO*AmBrO)
               TBrOP = exp(-ABrO*AbBrO*Amdif)
            end if

            ! 9. Formaldehyde (CH2O)
            TCH2O = 1.0
            TCH2OP = 1.0
            if(AbCH2O>0.0 .and. wvln<=400.0) then
               call SearchWVLM(Abs_CH2O(:,1), wvln, iabs_indx)
               wvlo = Abs_CH2O(iabs_indx,1)
               xsCH2O = Abs_CH2O(iabs_indx,2)
               atCH2O = Abs_CH2O(iabs_indx,3)
               Call GSCH2O(TK-24.,xsCH2O,atCH2O,AbCH2O,AmCH2O, &
                           TCH2O,TCH2OP,Amdif)
            end if

            ! 10. Nitrous acid (HNO2)
            THNO2 = 1.0
            THNO2P = 1.0
            if(AbHNO2>0.0 .and. wvln>=300.5 .and. wvln<=396.5) then
               call SearchWVLM(Abs_HNO2(:,1), wvln, iabs_indx)
               wvlo = Abs_HNO2(iabs_indx,1)
               xsHNO2 = Abs_HNO2(iabs_indx,2)
               Call GSHNO2(xsHNO2,AbHNO2,AmHNO2,THNO2,THNO2P,Amdif)
            end if

            ! 11. Chlorine nitrate (ClNO3)
            TClNO = 1.0
            TClNOP = 1.0
            TCl = 230.0
            if(wvln<=432.0) then
               call SearchWVLM(Abs_ClNO(:,1), wvln, iabs_indx)
               wvlo = Abs_ClNO(iabs_indx,1)
               xsClNO = Abs_ClNO(iabs_indx,2)
               a1tCl = Abs_ClNO(iabs_indx,3)
               a2tCl = Abs_ClNO(iabs_indx,4)
               AClNO = xsClNO*(1.+a1tCl*(TCl-296.)+a2tCl*(TCl-296.) &
                     *(TCl-296.))*NLosch
               TClNO = exp(-AClNO*AbClNO*AmClNO)
               TClNOP = exp(-AClNO*AbClNO*Amdif)
            end if

            ! Total gaseous absorption excluding H2O and O3
            ! Mixed gases
            Tmixd = TO2*TO4*TN2*TN2O*TCO*TCO2*TCH4
            TmixdP = TO2P*TO4P*TN2P*TN2OP*TCOP*TCO2P*TCH4P
            Trace = TNO*TNO2*TNO3*THNO3*TSO2*TNH3*TBrO*TCH2O*THNO2*TClNO
            TraceP = TNOP*TNO2P*TNO3P*THNO3P*TSO2P*TNH3P*TBrOP*TCH2OP*THNO2P*TClNOP
            if (Tmixd<e30 .or. Trace<e30) then
               fgphot(isptm) = 0.0_r8
               ofrdif(isptm) = 1.0_r8
               cycle
            end if
            Taumix = -log(Tmixd)/AmR
            Tautrc = -log(Trace)/AmR

            ! AEROSOL EXTINCTION
            TAUA = 0.0
            TAUAS = 0.0
            TAA = 1.0
            TAAP = 1.0
            TAS = 1.0
            TAT = 1.0

            if (wvln>1999.) then
               BQ = exp(BQ1*(WVL-BQ2))
               OMEGL = 1.0 - (BQ0*BQ)/(1.+BQ)**2
            else
               OMEGL = BP0+BP1*WVL+BP2*WVL2+BP3*WVL3
            end if

            if (beta>0.0) then
               if (wvln<=499.) then
                  TAUA = BETAC/(WVL**alpha1)
               else if (wvln>=1001.0 .and. IAER==8) then
                  TAUA = beta/WVL**.825
               else
                  TAUA = beta/WVL**alpha2
               end if
               if(IAER==1 .and. IREF==1) then
                  if (wvl<=0.2999) then
                     alpha = 1.1696-4.6814*wvl+12.96*wvl2
                  else if (abs(wvl-0.3)<1d-4) then
                     alpha = 0.93178
                  else if (wvl>=0.3001 .and. wvl<=0.3369) then
                     alpha = 1.443-4.6051*wvl+9.6723*wvl2
                  else if (abs(wvl-0.337)<1d-4) then
                     alpha = 0.989
                  else if (wvl>=0.3371 .and. wvl<=0.5499) then
                     alpha = 0.98264+.032539*wvl-.040251*wvl2
                  else if (abs(wvl-0.55)<1d-4) then
                     alpha = 0.98839
                  else if (wvl>=0.5501 .and. wvl<=0.6939) then
                     alpha = -32.0108+151.02*wvl-229.75*wvl2+116.83*wvl3
                  else if (abs(wvl-0.694)<1d-4) then
                     alpha = 1.192
                  else if (wvl>=0.6941 .and. wvl<=1.0599) then
                     alpha = -1.9669+9.576*wvl-9.4345*wvl2+3.1621*wvl3
                  else if (abs(wvl-1.06)<1d-4) then
                     alpha = 1.3485
                  else if (wvl>=1.0601 .and. wvl<=1.536) then
                     alpha = -.25628+3.0677*wvl-1.9011*wvl2+.41005*wvl3
                  else if (wvl>1.536 .and. wvl<=2.0) then
                     alpha = -1.3018+3.7405*wvl-1.6633*wvl2+.25856*wvl3
                  else if (wvl>2.0 .and. wvl<=2.25) then
                     alpha = 2.1665-.40189*wvl+.057873*wvl2
                  else if (wvl>2.25 .and. wvl<=2.5) then
                     alpha = 2.1188-.35073*wvl+.044553*wvl2
                  else if (wvl>2.5 .and. wvl<=2.7) then
                     alpha = 4.3108-1.5493*wvl+.17324*wvl2
                  else if (wvl>2.7 .and. wvl<=3.0) then
                     alpha = 2.1947-.33892*wvl+.015213*wvl2
                  else if (wvl>3.0 .and. wvl<=3.39) then
                     alpha = -2.993+3.3795*wvl- .86713*wvl2+.073101*wvl3
                  else if (wvl>3.39 .and. wvl<=3.75) then
                     alpha = 1.6801-.12171*wvl+.0068994*wvl2
                  else if(wvl>3.75) then
                     alpha = 2.0473-.27977*wvl+.022939*wvl2
                  end if
                  TAUA = Tau5/(2.*wvl)**alpha
               end if
               TAUAS = OMEGL*TAUA
               TAUAA = TAUA - TAUAS
               TAS = exp(-TAUAS*AmAER)
               TAT = exp(-TAUA*AmAER)
               TAA = exp(-TAUAA*AmAER)
            else
               OMEGL = 1.0
            end if

            ! BEAM RADIATION
            TSCAT = TR*TAS
            TABS0 = TH2O*Tmixd*Trace*TAA
            TAAp = exp(-TAUAA*Amdif)
            TABS0P = TH2OP*TmixdP*TraceP*TAAP
            TABS = TABS0*TO3
            TDIR = TABS*TSCAT
            DIR = H0*TDIR
            DIRH = DIR*ZCOS
            H0H = H0*ZCOS
            FHTO = 1.0
            FHT1 = 1.0
            if (TAUZ3>5d-6) then
               FHTO = exp(-FHTcz-FHTdx*(TAUZ3-2.0))
            end if
            if (TAUZ3<=2.0) then
               FHTO = exp(-FHTcx*TAUZ3-FHTcy*(TAUZ3**.95))
            end if

            ! Improved multiple scattering algorithm
            Fda00 = 1.0
            Fdazt = 1.0
            if (t5>0.03) then
               ssaro = Taurl/(Taurl+Tauz3)
               Taurf = (ssaro**.5)*Taurl*Taurl
               Fda00 = (trb0+trb1*Taurf)/(1.+trb2*(Taurf**.5))
               if (wvln>294) then
                  Fdazt = tzb0+tzb1*Taurf+tzb2*Taurf**.5
               end if
            end if
            if (wvln<=400) then
               Fda00 = 1.0
               if(wvln>294.0 .and. t5>0.03) then
                  Fda00 = (tra0+tra1*Taurf)/(1.+tra2*Taurf)   
               end if
               FHT1 = 0.962 - 9.1*tauz3

               if (Tauz3>0.01) then
                  if (Tauz3<=22.5) then
                     if (Tauz3<=15.5) then
                        if (Tauz3<=10) then
                           if (Tauz3>6.0) then
                              if (Amo3>2.2) then 
                                 ra0 = ra01
                              else
                                 ra0 = ra00
                              end if
                              if (Amo3<=1.72) then
                                 ra1 = ra12
                              else if (Amo3<=2.6) then
                                 ra1 = ra11
                              else
                                 ra1 = ra10
                              end if
                              FHT1 = min(10.D00, ra0+(Tauz3-6.D00)*ra1)
                           else
                              if (Tauz3<=1.0) then
                                 if (Tauz3<=0.1) then
                                    FHT1 = min(1.0, ra06+Tauz3*ra18)
                                 else
                                    if (Amo3>2.0) then
                                       FHT1 = min(1.6, ra05+Tauz3*ra17)
                                    else
                                       FHT1 = min(1.0, ra05+Tauz3*ra17)
                                    end if
                                 end if
                              else
                                 if (Amo3>3.2) then
                                    ra0 = ra03
                                    ra1 = ra14
                                    FHT1 = min(2.0, ra0+Tauz3*ra1)
                                    FHT1x = min(2.0, ra0+2.505*ra1) 
                                 else
                                    ra0 = ra02
                                    ra1 = ra13
                                    FHT1 = min(2., ra0+Tauz3*ra1)
                                    FHT1x = min(2., ra0+2.505*ra1)
                                 end if
                                 if (Tauz3>2.505) then
                                    if (Amo3>3.5) then
                                       ra0 = ra04
                                    else
                                       ra0 = FHT1x
                                    end if
                                    if (Amo3>2.4) then
                                       ra1 = ra16
                                    else
                                       ra1 = ra15
                                    end if
                                    FHT1 = min(7.5, ra0+(Tauz3-2.505D00)*ra1)
                                 end if
                              end if
                           end if
                        else
                           if(Amo3>1.9) then
                              FHT1 = min(8.D00, FHTf0+Tauz3*FHTf1)
                           else
                              FHT1 = min(9.D00, FHTe0+Tauz3*FHTe1)
                           end if
                        end if
                     else
                        if (Amo3>1.6) then
                           FHT1 = min(6.D00, FHTd0+Tauz3*FHTd1)
                        else
                           FHT1 = min(7.D00, FHTc0+Tauz3*FHTc1)
                        end if
                     end if
                  else
                     if (Amo3>2.0) then
                        FHT1 = min(12.0, FHTb0+Tauz3*FHTb1)
                     else
                        FHT1 = min(12.0, FHTa0+Tauz3*FHTa1)
                     end if
                  end if
               end if
            end if
            Fdazt = min( max(Fdazt, 0.0), 1.0 )
            Fda00 = min( max(Fda00, 0.0), 1.0 )
            if (Tauz3<5.0D00) then
               FHT1 = max(FHT1, 5.D-01)
            else
               FHT1 = max(FHT1, 2.D-03)
            end if
            FHTO = FHTO/FHT1
            if (Zenit>=89.) then
               HT = H0/AmR
            else
               HT = H0H
            end if
            HTa = HT*TABS0p*FHTO

            ! Diffuse radiation
            ! Asymmetry and forward scatterance
            GG = AG0+AG1*WVL+AG2*WVL2+AG3*WVL3+AG4*WVL4
            GG = min(0.99, GG)
            ALG = log(1.0-GG)
            AFS = ALG*(1.459+ALG*(.1595+ALG*.4129))
            BFS = ALG*(.0783+ALG*(-.3824-ALG*.5874))
            FA1 = 1.0 - 0.5*exp((AFS+BFS*ZCOS)*ZCOS)
            FA1P = 1.0 - 0.5*exp((AFS+BFS*.6)*.6)

            ! 1. DIFFUSE RADIATION FROM RAYLEIGH SCATTERING
            DRAY = 0.0
            FR = 0.5
            FRP = 0.5
            if (TAURL>=EPSIR) then
               FR = 0.5*exp(-((TAURL-EPSIR)/SIGMAR)**EXPR)
            end if 
            if (TAURL>=EPSIRP) then
               FRP = 0.5*exp(-.1957*(TAURL-0.0648)**1.32)
            end if
            DRAY = HTa*FR*(1.-TR)*(Max(Tas,1.d-10)**.167)

            ! 2. DIFFUSE RADIATION FROM AEROSOL SCATTERING
            if (beta>0.0) then
               DAER = HTa*FA1*(TR**.167)*(1.-TAS)*Fda00*Fdazt
            else
               DAER = 0.D00
            end if
            ! Sky diffuse before backscattering
            Fdifz = 1.0
            Fdiftz = 1.0
            if (wvln<=294.) then
               Fdifz = fdifa0 + fdifa1*wvl
               if (zenit>=45.) then
                  Fdiftz = fdifb0+fdifb1*Ama1+fdifb2*Ama2
               end if
            end if
            DIF0 = Fdifz*Fdiftz*(DRAY+DAER)
            Glob0 = DIRH + Dif0
               
            ! 3. BACKSCATTERING - reflection from ground to space and back
            TRP = exp(-Amdif*TAURL)
            TAUAP = TAUA*Amdif
            TASP = exp(-OMEGL*TAUAP)
            TTp5 = Tabs0p**.5
            if (wvln<=379.5) then
               GAMOZ = exp(-1D+5*(4.8344+23.088*(AbO3+ApO3))*(.38-WVL)**5.8)
            else
               GAMOZ = 1.0
            end if
            Rhob0 = Rhox
            Rhob = Rhox
            Rhod = Rhox
            Rhor = 0.0
            Rhoa = 0.0
            Rhos = 0.0
            Roro = 0.0
            Dgrnd = 0.0
            Rho = max(Rhob,Rhod)
            if (TTP5>1.0d-12) then
               Rhor = TTP5*((1.-FRP)**.85)*(TASP**.05)*(1.-TRP)*GAMOZ
               Fatau = 0.0
               if (tau550>=0.03) then
                  if (wvl<=0.35) then
                     Fatau = exp(alba00+alba01*wvl)
                     if (t5>0.2) then
                        Fatau = exp(alba0+alba1*wvl+alba2*wvl2)
                     end if
                  else
                     if(wvl>0.5) then
                        Fatau = (albc0+albc1*wvl+albc2*wvl2+albc3*wvl3)/ &
                                 (1.+albc4*wvl)
                     else
                        Fatau = exp(albb0+albb1*wvl+albb2*wvl2)
                     end if
                  end if
                  Rhoa = TTP5*((1.-FA1P)**.85)*GAMOZ*Fatau
               end if
               Rhos = Rhor + Rhoa
               Roro = Rho*Rhos
               Upward = Rho*Glob0
               Dgrnd = Upward*Rhos/(1.-Roro)
            end if
            DIF = DIF0 + Dgrnd
            GLOB = DIRH + DIF

            ! PHOTON FLUX per wavelength (mole photons m-2 s-1 nm-1)
            ! save for radiation and diffuse fraction 
            WPHT = 1.0d-6*WVL*PHOT/NA
            fgphot(isptm) = GLOB*WPHT
            if (GLOB>e30) then
               ofrdif(isptm) = DIF / GLOB
            else
               ofrdif(isptm) = 1.0
            end if
         end if
      end do   ! End of SPECTRAL CALCULATIONS

   end subroutine

   !------------------------------------------------------------------------------
   !
   ! some utilities
   ! 
   !------------------------------------------------------------------------------
   function AMZ(X)
      implicit none
      real(r8), intent(in) :: X
      real(r8) :: AMZ, XR, XCOS

      XR = X*0.017453293d+0
      XCOS = cos(XR)
      AMZ = 1./(XCOS+.48353*(X**.095846)*(96.741-X)**(-1.754))
      return
   end function

   subroutine SearchWVLM(vwvlm, wvln, iabs_indx)
      implicit none
      real(r8), intent(in) :: vwvlm(:)
      real(r8), intent(in) :: wvln
      integer, intent(out) :: iabs_indx
      integer :: indx, nn

      nn = size(vwvlm)
      call BinarySearch(vwvlm, wvln, indx)
      if (indx<nn) then
         if (abs(wvln-vwvlm(indx))<=epsilm) then
            iabs_indx = indx
         else if (abs(wvln-vwvlm(indx+1))<=epsilm) then
            iabs_indx = indx + 1
         end if
      else
         if (abs(wvln-vwvlm(indx))<=epsilm) then
            iabs_indx = indx
         end if
      end if
   end subroutine

   Subroutine GSCH4(ACH4,AbCH4,AmCH4,TCH4,TCH4P,Amdif)
      implicit none
      real(r8), intent(in) :: ACH4, AbCH4, AmCH4, Amdif
      real(r8), intent(out) :: TCH4, TCH4P

      TCH4 = exp(-(ACH4*AbCH4*AmCH4))
      TCH4P = exp(-(ACH4*AbCH4*Amdif))
   end subroutine

   subroutine GSCO(ACO,AbCO,AmCO,TCO,TCOP,Amdif)
      implicit none
      real(r8), intent(in) :: ACO, AbCO, AmCO, Amdif
      real(r8), intent(out) :: TCO, TCOP

      TCO = exp(-ACO*AbCO*AmCO)
      TCOP = exp(-ACO*AbCO*Amdif)
   end subroutine

   subroutine GSHNO3(T,xsHNO3,athno3,AbHNO3,AmHNO3, &
                     THNO3,THNO3P,Amdif)
      implicit none
      real(r8), intent(in) :: T, xsHNO3, athno3
      real(r8), intent(in) :: AbHNO3, AmHNO3, Amdif
      real(r8), intent(out) :: THNO3, THNO3P
      real(r8) :: xnl, AHNO3

      xnl = NLosch*1d-19
      AHNO3 = xsHNO3*0.1*exp(athno3*1d-3*(T-298.))*xnl
      THNO3 = exp(-AHNO3*AbHNO3*AmHNO3)
      THNO3P = exp(-AHNO3*AbHNO3*Amdif)
   end subroutine

   subroutine GSNO(ANO,AbNO,AmNO,TNO,TNOP,Amdif)
      implicit none
      real(r8), intent(in) :: ANO, AbNO, AmNO, Amdif
      real(r8), intent(out) :: TNO, TNOP

      TNO = exp(-ANO*AbNO*AmNO)
      TNOP = exp(-ANO*AbNO*Amdif)
   end subroutine

   subroutine GSNO2(T,xsNO2,atno2,AbNO2,AmNO2,TNO2,TNO2P,Amdif)
      implicit none
      real(r8), intent(in) :: T, xsNO2, atno2
      real(r8), intent(in) :: AbNO2, AmNO2, Amdif
      real(r8), intent(out) :: TNO2, TNO2P 
      real(r8) :: ANO2

      ANO2 = (xsNO2+atno2*(T-220.))*NLosch
      TNO2 = exp(-ANO2*AbNO2*AmNO2)
      TNO2P = exp(-ANO2*AbNO2*Amdif)
   end subroutine

   subroutine GSNO3(T,xsNO3,atno3,AbNO3,AmNO3,TNO3,TNO3P,Amdif)
      implicit none
      real(r8), intent(in) :: T, xsNO3, atno3
      real(r8), intent(in) :: AbNO3, AmNO3, Amdif
      real(r8), intent(out) :: TNO3, TNO3P
      real(r8) :: ANO3

      ANO3 = (xsNO3+atno3*(T-230.))*NLosch
      TNO3 = exp(-ANO3*AbNO3*AmNO3)
      TNO3P = exp(-ANO3*AbNO3*Amdif)
   end subroutine

   subroutine GSSO2U(T,xsSO2,atso2,AbSO2,AmSO2,TSO2,TSO2P,Amdif)
      implicit none
      real(r8), intent(in) :: T, xsSO2, atso2
      real(r8), intent(in) :: AbSO2, AmSO2, Amdif
      real(r8), intent(out) :: TSO2, TSO2P
      real(r8) :: ASO2

      ASO2 = (xsSO2+atso2*(T-213.))*NLosch
      TSO2 = exp(-ASO2*AbSO2*AmSO2)
      TSO2P = exp(-ASO2*AbSO2*Amdif)
   end subroutine

   subroutine GSSO2I(ASO2,AbSO2,AmSO2,TSO2,TSO2P,Amdif)
      implicit none
      real(r8), intent(in) :: ASO2, AbSO2, AmSO2, Amdif
      real(r8), intent(out) :: TSO2, TSO2P
      
      TSO2 = exp(-(ASO2*AbSO2*AmSO2))
      TSO2P = exp(-(ASO2*AbSO2*Amdif))
   end subroutine

   subroutine GSO3U(Tref,T,xs,a0o3,a1o3,Abo3,Amo3,TAUZ3,TO3,xso3, &
                    AO3)
      implicit none
      real(r8), intent(in) :: Tref, T, xs
      real(r8), intent(in) :: a0o3, a1o3, Abo3, Amo3
      real(r8), intent(out) :: TAUZ3, TO3, xso3, AO3
      real(r8) :: fto3, XCHECK

      fto3 = 1.0 - T/Tref
      xso3 = xs*(1.+a0o3*fto3)/(1.+a1o3*fto3)
      AO3 = NLosch*xso3
      TAUZ3 = AbO3*AO3
      XCHECK = AMO3*TAUZ3
      if (XCHECK<=499.) then
         TO3 = exp(-XCHECK)
      else
         TO3 = 0.D+000
      end if
   end subroutine

   subroutine GSCH2O(T,xsCH2O,atCH2O,AbCH2O,AmCH2O,TCH2O,TCH2OP, &
                     Amdif)
      implicit none
      real(r8), intent(in) :: T, xsCH2O, atCH2O
      real(r8), intent(in) :: AbCH2O, AmCH2O, Amdif
      real(r8), intent(out) :: TCH2O, TCH2OP
      real(r8) :: ACH2O

      ACH2O = (xsCH2O+atCH2O*(T-293.))*NLosch
      TCH2O = exp(-ACH2O*AbCH2O*AmCH2O)
      TCH2OP = exp(-ACH2O*AbCH2O*Amdif)
   end subroutine

   subroutine GSHNO2(xsHNO2,AbHNO2,AmHNO2,THNO2,THNO2P,Amdif)
      implicit none
      real(r8), intent(in) :: xsHNO2, AbHNO2, AmHNO2, Amdif
      real(r8), intent(out) :: THNO2, THNO2P
      real(r8) :: AHNO2

      AHNO2 = xsHNO2*NLosch
      THNO2 = exp(-AHNO2*AbHNO2*AmHNO2)
      THNO2P = exp(-AHNO2*AbHNO2*Amdif)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate alpha parameters for aerosol model.
   !          Overall alpha1, alpha2; considering TROPO aerosols in the free
   !          atmosphere for Shettle & Fenn models
   !
   !------------------------------------------------------------------------------ 
   Subroutine ALFA(season, Iaer, Iturb, Iref, alpha1, alpha2, tau5, &
                   beta,al1,al2,tau550,indx)
      implicit none
      integer, intent(in) :: season
      integer, intent(in) :: Iaer
      integer, intent(in) :: Iturb
      integer, intent(in) :: Iref
      real(r8), intent(in) :: alpha1
      real(r8), intent(in) :: alpha2
      real(r8), intent(inout) :: tau5
      real(r8), intent(inout) :: beta
      real(r8), intent(inout) :: al1
      real(r8), intent(inout) :: al2
      real(r8), intent(inout) :: tau550
      integer, intent(in) :: indx
      real(r8) :: alfa1s, alfa1w, alfa2s, alfa2w, frctn
      real(r8) :: alfa1, alfa2, t55, frctns, frctnw
      real(r8) :: t552

      alfa1s = .999
      alfa1w = .999
      alfa2s = 1.565
      alfa2w = 1.584
      frctn = 0.0
      alfa1 = 0.0
      alfa2 = 0.0
      t55 = tau5/(1.1**alpha2)
      if(indx==1) T55=tau550
      frctns = 1.0
      frctnw = 1.0
      if (Iaer/=4) then
         t552 = t55*t55
         frctns = Min(1.0, exp((-1.6435+1.4011*t55-8.2491*t552+.065552/t55)/ &
                  (1.+1.936*t552)))
         frctnw = Min(1.0, exp((-1.8181+.84983*t55-8.641*t552+.0436/t55)/ &
                  (1.+1.8959*t552)))
      end if
      frctn = frctns
      alfa1 = alfa1s
      alfa2 = alfa2s
      if(Season==0) then
         frctn = frctnw
         alfa1 = alfa1w
         alfa2 = alfa2w
      end if
      al1 = alpha1*(1.-frctn)+alfa1*frctn
      al2 = alpha2*(1.-frctn)+alfa2*frctn

      if(indx/=1) then
         ! Recalculate Tau550 from Tau5 and the new value of alpha2
         if(iTurb==1) tau5 = beta/(0.5**al2)
         if(iTurb==0 .or. iTurb==2 .or. Iturb==5) beta = (0.5**al2)*tau5
         if(iTurb/=5) tau550 = tau5/(1.1**al2)
         if(Iaer/=1 .and. Iref==1) then
            if(iTurb==1) tau5 = beta/(0.5**1.336688)
            if(iTurb==0 .or. iTurb==2 .or. Iturb==5) beta = (0.5**1.33669)*tau5
            if(iTurb/=5) tau550 = tau5/(1.1**.988415)
         end if
      else
         ! Recalculate Tau5 from Tau550 and the new value of alpha2
         tau5 = tau550*(1.1**al2)
         beta = (0.55**al2)*tau550
         if(Iaer==1 .and. Iref==1) then
            tau5 = tau550*(1.1**.9883)
            beta = tau550*(.55**1.39223) 
         end if
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: AVERAGE STRATOSPHERIC TEMPERATURE AND REFERENCE ATMOSPHERIC 
   !          CONDITIONS.
   !
   !------------------------------------------------------------------------------ 
   subroutine RefAtm(Z, Pm, Tm, TO3m, O3m, RHm, Wm, TO3ini, k)
      implicit none
      real(r8), intent(in) :: Z
      real(r8), intent(out) :: Pm
      real(r8), intent(out) :: Tm
      real(r8), intent(out) :: TO3m
      real(r8), intent(out) :: O3m
      real(r8), intent(out) :: RHm
      real(r8), intent(out) :: Wm
      real(r8), intent(out) :: TO3ini
      integer, intent(in) :: k
      integer, parameter :: NN = 50
      integer :: INTRP, ii, J1
      real(r8) :: xp, xp1
      
      TO3ini = TO3(1,k)
      INTRP = 0
      do ii = 1, NN, 1
         if (INTRP/=0) then
            exit
         end if
         if (abs(Z-Zref50(ii))<=1.d-4) then
            INTRP = 2
            RHm = RHref(ii,k)
            Tm = Tref(ii,k)
            Pm = Pref(ii,k)
            O3m = O3ref(ii,k)
            TO3m = TO3(ii,k)
            Wm = Wref(ii,k)
         else if (Z<Zref50(ii)) then
            if (ii<=2) then
               INTRP = 1
               J1 = max(ii-1,1)
            else if (ii>=NN-1) then
               INTRP = 1
               J1 = min(ii-1,NN-1)
            else
               INTRP = 1
               J1 = ii - 1
            end if
         end if
      end do

      if (INTRP/=2) then
         xp = Zref50(J1)
         xp1 = Zref50(J1+1)
         Call Interp(1,xp,xp1,Z,RHref(J1,k),RHref(J1+1,k),RHm)
         Call Interp(1,xp,xp1,Z,Tref(J1,k),Tref(J1+1,k),Tm)
         Call Interp(1,xp,xp1,Z,Pref(J1,k),Pref(J1+1,k),Pm)
         Call Interp(1,xp,xp1,Z,O3ref(J1,k),O3ref(J1+1,k),O3m)
         Call Interp(1,xp,xp1,Z,Wref(J1,k),Wref(J1+1,k),Wm)
         Call Interp(1,xp,xp1,Z,TO3(J1,k),TO3(J1+1,k),TO3m)
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the zone temperature and ozone conditions
   !
   !------------------------------------------------------------------------------
   subroutine Ozon(Z,Tz,Tm,Tmin,Tmax,Ozmin,Ozmax)
      implicit none
      real(r8), intent(in) :: Z
      real(r8), intent(in) :: Tz
      real(r8), intent(out) :: Tm
      real(r8), intent(out) :: Tmin
      real(r8), intent(out) :: Tmax
      real(r8), intent(out) :: Ozmin
      real(r8), intent(out) :: Ozmax
      integer, parameter :: NN = 41
      integer :: INTRP, IZ, J1
      real(r8) :: xp, xp1, Tm0, Tm1, Tmin0, Tmin1
      real(r8) :: Tmax0, Tmax1, Ozmin0, Ozmin1
      real(r8) :: Ozmax0, Ozmax1

      INTRP=0
      do IZ=1, NN, 1
         if(INTRP/=0) then
            exit
         end if
         if (abs(Z-Zref41(IZ))<=1.d-4) then
            INTRP = 2
            Tm = (1.-O3a1(IZ))*Tz-O3a0(IZ)
            Tmin = TO3min(IZ)
            Tmax = TO3max(IZ)
            Ozmin = O3min(IZ)
            Ozmax = O3max(IZ) 
         else
            if (Z<=Zref41(IZ)) then
               if (IZ>=NN-1) then
                  INTRP = 1
                  J1 = Min(IZ-1,NN-1)
               else if (IZ>2) then
                  INTRP = 1
                  J1 = IZ-1
               else
                  INTRP = 1
                  J1 = Max(IZ-1,1)
               end if
            end if
         end if
      end do

      if(INTRP/=2) then
         xp = Zref41(J1)
         xp1 = Zref41(J1+1)
         Tm0 = (1.-O3a1(J1))*Tz-O3a0(J1)
         Tm1 = (1.-O3a1(J1+1))*Tz-O3a0(J1+1)
         Tmin0 = TO3min(J1)
         Tmin1 = TO3min(J1+1)
         Tmax0 = TO3max(J1)
         Tmax1 = TO3max(J1+1)
         Ozmin0 = O3min(J1)
         Ozmin1 = O3min(J1+1)
         Ozmax0 = O3max(J1)
         Ozmax1 = O3max(J1+1)

         Call Interp(1,xp,xp1,Z,Tm0,Tm1,Tm)
         Call Interp(1,xp,xp1,Z,Tmin0,Tmin1,Tmin)
         Call Interp(1,xp,xp1,Z,Tmax0,Tmax1,Tmax)
         Call Interp(1,xp,xp1,Z,Ozmin0,Ozmin1,Ozmin)
         Call Interp(1,xp,xp1,Z,Ozmax0,Ozmax1,Ozmax)
      end if
   end subroutine

   Subroutine Ozon2(Z, Ozmin, Ozmax)
      implicit none
      real(r8), intent(in) :: Z
      real(r8), intent(out) :: Ozmin
      real(r8), intent(out) :: Ozmax
      integer, parameter :: NN = 41
      integer :: INTRP, ii, J1
      real(r8) :: xp, xp1, Ozmin0, Ozmin1
      real(r8) :: Ozmax0, Ozmax1

      INTRP = 0
      do ii = 1, NN, 1
         if (INTRP/=0) then
            exit
         end if
         if (abs(Z-Zref41(ii))<=1.d-4) then
            INTRP = 2
            Ozmin = O3min(ii)
            Ozmax = O3max(ii)
         else if (Z<Zref41(ii)) then
            if (ii<=2) then
               INTRP = 1
               J1 = max(ii-1,1)
            else if (ii>=NN-1) then
               INTRP = 1
               J1 = min(ii-1,NN-1)
            else
               INTRP = 1
               J1 = ii - 1
            end if
         end if
      end do

      if (INTRP/=2) then
         xp = Zref41(J1)
         xp1 = Zref41(J1+1)
         Ozmin0 = O3min(J1)
         Ozmin1 = O3min(J1+1)
         Ozmax0 = O3max(J1)
         Ozmax1 = O3max(J1+1)

         Call Interp(1,xp,xp1,Z,Ozmin0,Ozmin1,Ozmin)
         Call Interp(1,xp,xp1,Z,Ozmax0,Ozmax1,Ozmax)
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   !  Calculate sun position information ('PSA Algorithm' by Blanco-Muriel et
   !  al. (Solar Energy 2001))
   !
   !------------------------------------------------------------------------------
   subroutine SunPSA(dHour, dLat, dLong, decli, Zenith, Azimu, Julian, R, &
                     EOT, P, T, Year, Month, Day)
      implicit none
      integer, intent(in) :: Year
      integer, intent(in) :: Month
      integer, intent(in) :: Day
      real(r8), intent(in) :: P        ! pressure (mb)
      real(r8), intent(in) :: T        ! temperature (K)
      real(r8), intent(in) :: dHour    ! decimal Hour (Universal Time, UT) 
      real(r8), intent(in) :: dLat     ! decimal Latitude (deg., positive North)
      real(r8), intent(in) :: dLong    ! decimal Longitude (deg., positive East)
      real(r8), intent(out) :: Zenith
      real(r8), intent(out) :: Azimu
      real(r8), intent(out) :: Julian
      real(r8), intent(out) :: decli
      real(r8), intent(out) :: R
      real(r8), intent(out) :: EOT
      real(r8) :: pinb, twopi, rad, LMST, GMST, dRLat
      real(r8) :: HrAngl, cosHA, Elapsd, dX, dY, EclipL
      real(r8) :: Paralx, AUnit, Radius, EclipO, RightA
      real(r8) :: Anomly, Omega, sinELg, SunLng
      real(r8) :: cosLat, sinLat, xsun, ELD, ELD2
      real(r8) :: REFR, PT
      integer(i8) :: liAux1, liAux2 

      ! Declaration of some constants
      pinb = 3.14159265358979323846d+0
      twopi = 2.*pinb
      rad = pinb/180.
      Radius = 6371.01d+0
      AUnit = 149597890.0d+0

      ! Calculate current Julian Day
      liAux1 = (Month-14)/12
      liAux2 = (1461*(Year + 4800 + liAux1))/4 &
       + (367*(Month- 2-12*liAux1))/12 &
       - (3*((Year + 4900+ liAux1)/100))/4 + Day - 32075
      Julian = DBLE(liAux2) -0.5 + dHour/24.0

      ! Calculate difference in days between the current Julian Day 
      ! and JD 2451545.0, which is noon 1 January 2000 Universal Time 
      Elapsd = Julian - 2451545.0d+0

      ! Calculate ecliptic coordinates (ecliptic longitude and obliquity of 
      ! the ecliptic in radians) but without limiting the angle to be less 
      ! than 2*Pi (i.e., the result may be greater than 2*Pi)
      Omega = 2.1429d+0 - 0.0010394594d+0 * Elapsd
      SunLng = 4.8950630d+0 + 0.017202791698d+0 * Elapsd
      Anomly = 6.2400600d+0 + 0.0172019699d+0 * Elapsd
      EclipL = SunLng + 0.03341607d+0*sin(Anomly) + 0.00034894d+0* &
            sin(2.*Anomly) - 0.0001134d+0 - 0.0000203d+0*sin(Omega)
      EclipO = 0.4090928d+0 - 6.2140d-9*Elapsd + 0.0000396d+0*cos(Omega)

      ! Calculate celestial coordinates (right ascension and declination)       
      ! in radians but without limiting the angle to be less than 2*Pi 
      ! (i.e., the result may be greater than 2*Pi)
      SinELg = sin(EclipL)
      dY = cos(EclipO) * SinELg
      dX = cos(EclipL)
      RightA = atan2(dY,dX)
      if(RightA<0.0) RightA = RightA + twopi
      Decli = asin(sin(EclipO)*SinELg)

      ! Calculate local coordinates (azimuth and zenith angle) in degrees
      GMST = 6.6974243242d+0 + 0.0657098283d+0*Elapsd + dHour
      LMST = (GMST*15. + dLong)*rad
      HrAngl = LMST - RightA
      dRLat = dLat*rad
      cosLat = cos(dRLat)
      sinLat = sin(dRLat)
      cosHA= cos(HrAngl)
      Zenith = acos(cosLat*cosHA*cos(Decli) + sin(Decli)*sinLat)
      dY = -sin(HrAngl)
      dX = tan(Decli)*cosLat - sinLat*cosHA
      Azimu = atan2(dY, dX)
      if (Azimu<0.0) Azimu = Azimu + twopi
      Azimu = Azimu/rad

      ! Parallax Correction
      Paralx = (Radius/AUnit)*sin(Zenith)
      Zenith = (Zenith + Paralx)/rad

      ! Sun-Earth actual distance in AU (from Michalsky's paper)
      R = 1.00014 - 0.01671*cos(Anomly) - 0.00014*cos(2.*Anomly)

      ! Equation of Time (in min, from Michalsky's paper)
      RightA = RightA/rad
      SunLng = SunLng/rad
      xsun = -aint(abs(SunLng)/360.)
      if(Sunlng<0.0) xsun = -xsun + 1.0
      SunLng = SunLng + xsun*360.
      EOT = (SunLng-RightA)*4.

      ! REFRACTION CORRECTION FOR ACTUAL ATMOSPHERIC CONDITIONS (P,T)
      ELD = 90.-Zenith
      ELD2 = ELD*ELD
      REFR = 0.0
      PT = P/T
      if(ELD<15.0 .and. ELD>=-2.5) then
         REFR = PT*(.1594+.0196*ELD+2d-5*ELD2)/(1.+.505*ELD+.0845*ELD2)
      else if (ELD>=15.0 .and. ELD<90.0) then
         REFR = .00452*PT/tan(ELD*rad)
      end if
      Zenith = 90.0 - (ELD+REFR)

      ! Declination in degrees
      Decli = Decli/rad
   end subroutine

   subroutine VISTAU(Season, Range1, Tau, indx)
      implicit none
      integer, intent(in) :: Season
      real(r8), intent(inout) :: Range1
      real(r8), intent(inout) :: Tau
      integer, intent(in) :: indx
      real(r8), parameter :: vs1(5) = (/-3.2998,-5.37,156.14,42.389,48.957/)
      real(r8), parameter :: vs2(3) = (/.026483,7.133,-6.6238/)
      real(r8), parameter :: vs3(2) = (/.039987,.43928/)
      real(r8), parameter :: vw1(5) = (/-3.6629,-6.5109,165.85,44.857,51.968/)
      real(r8), parameter :: vw2(3) = (/.010149,6.7705,-1.7703/)
      real(r8), parameter :: vw3(2) = (/.023339,.27928/)
      real(r8) :: tln, delta2, delta1, Yvis, Yvis2

      if (indx==1) then
         ! (Index=1) Calculate Tau from Range
         Yvis = 1.0/RANGE1 - 1.0d-3
         Yvis2 = Yvis*Yvis
         if(Range1>=100.) then
            if(Range1>320.) then
               if (Season==0) then
                  Tau = vw3(1)+vw3(2)*Yvis
               else
                  Tau = vs3(1)+vs3(2)*Yvis
               end if
            else
               if (Season==0) then
                  Tau = vw2(1)+vw2(2)*Yvis+vw2(3)*Yvis2
               else
                  Tau = vs2(1)+vs2(2)*Yvis+vs2(3)*Yvis2
               end if
            end if
         else
            if (season==0) then
               Tau = exp((vw1(1)+vw1(2)*Yvis+vw1(3)*Yvis2)/ &
                     (1.+vw1(4)*Yvis+vw1(5)*Yvis2))
            else
               Tau = exp((vs1(1)+vs1(2)*Yvis+vs1(3)*Yvis2)/ &
                     (1.+vs1(4)*Yvis+vs1(5)*Yvis2))
            end if
         end if
      else
         ! (Index=0) Calculate Range from Tau
         tln = log(Tau)
         Range1 = 999.0
         if (Season/=0) then
            ! Calculations for SPRING/SUMMER conditions
            if (Tau>=0.0402) then
               if (Tau>0.0416) then
                  if(Tau<0.0901) then
                     delta2 = vs2(2)*vs2(2)+4.*vs2(3)*(Tau-vs2(1))
                     Range1 = min(999., 1./(.001+.5*(vs2(2)-(delta2**.5))/(-vs2(3)))) 
                  else
                     delta1 = (vs1(4)*tln-vs1(2))*(vs1(4)*tln-vs1(2))- &
                              4.*(vs1(5)*tln-vs1(3))*(tln-vs1(1))
                     Range1 = 1./(.001+.5*(vs1(2)-vs1(4)*tln-(delta1**.5))/ &
                              (vs1(5)*tln-vs1(3)))
                  end if
               else
                  Range1 = 1./(.001+(Tau-vs3(1))/vs3(2))
               end if
            end if
         else
            ! Calculations for FALL/WINTER conditions
            if (Tau>=0.0235) then
               if (Tau>0.0245) then
                  if (Tau>=0.0709) then
                     delta1 = (vw1(4)*tln-vw1(2))*(vw1(4)*tln-vw1(2))- &
                              4.*(vw1(5)*tln-vw1(3))*(tln-vw1(1))
                     Range1 = 1./(.001+.5*(vw1(2)-vw1(4)*tln-(delta1**.5))/ &
                              (vw1(5)*tln-vw1(3)))
                  else
                     delta2 = vw2(2)*vw2(2)+4.*vw2(3)*(Tau-vw2(1))
                     Range1 = min(999., 1./(.001+.5*(vw2(2)-(delta2**.5))/(-vw2(3))))
                  end if
               else
                  Range1 = 1./(.001+(Tau-vw3(1))/vw3(2))
               end if
            end if
         end if
         Range1 = min(Range1, 999.0)
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !  
   ! SATURATION VAPOR PRESSURE FROM GUEYMARD (J. Appl. Met. 1993)
   ! W=f(T,RH) USING EMPIRICAL MODEL OF GUEYMARD (SOLAR ENERGY 1994)
   !
   !------------------------------------------------------------------------------
   subroutine CalcPrecipitableWater(TK, TAIR, RH, W)
      implicit none
      real(r8), intent(in) :: TK    ! reference atmospheric temperature
      real(r8), intent(in) :: TAIR  ! surface air temperature
      real(r8), intent(in) :: RH    ! surface air humidity
      real(r8), intent(out) :: W    ! cm or g/cm^2
      real(r8) :: TK1, EVS, EV, ROV, TT, HV
      
      TK1 = TK/100.
      EVS = exp(22.329699-49.140396/TK1-10.921853/TK1/TK1-.39015156*TK1)
      EV = EVS*RH/100.
      ROV = 216.7*EV/TK
      TT = 1.0 + (TAIR/T0)
      HV = 0.4976 + 1.5265*TT + exp(13.6897*TT-14.9188*TT**3)
      W = 0.1*HV*ROV
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Interpolation routine using either a linear interpolation (IntTyp=1) 
   !          or logarithmic interpolation (IntTyp=2). The log interpolation 
   !          returns to a linear interpolation if the gradient is very small.
   !
   !------------------------------------------------------------------------------
   Subroutine Interp(IntTyp, x1, x2, X, y1, y2, y)
      implicit none
      integer, intent(in) :: IntTyp
      real(r8), intent(in) :: x1, x2
      real(r8), intent(in) :: y1, y2
      real(r8), intent(in) :: x
      real(r8), intent(out) :: y
      real(r8) :: a, a1, a2, grad
      integer :: it

      it = IntTyp
      If (y1<=1.d-5 .or. y2<=1.d-5) then
         it = 1
      end if

      if (it==1) then
         y = y1 + (x-x1)*(y2-y1)/(x2-x1)
      else
         a1 = log(y1)
         a2 = log(y2)
         grad = (a2-a1)/(x2-x1)
         if (abs(grad)<=1.d-2) then
            y = y1 + (x-x1)*(y2-y1)/(x2-x1)
         else
            a = a1 + (x-x1)*grad
            y = exp(a)
         end if
      end if
   end subroutine



   subroutine SetSmartsConstants()
      implicit none

      C1 = (/.4998,.27999,.049331,.57973/)
      C2 = (/45.236,55.642,7.9767,65.559/)
      C3 = (/96.233,1382.3,17.726,206.15/)
      C4 = (/-13.067,-132.47,-14.555,-26.911/)
      C5 = (/55.506,108.73,41.369,78.478/)
      C6 = (/83.115,1500.9,-18.384,166.38/)
      D1 = (/.86887,.69983,.039871,1.1194/)
      D2 = (/43.547,39.689,12.397,76.251/)
      D3 = (/-29.719,-26.736,98.641,129.32/)
      D4 = (/-1.8192,4.0596,-60.939,-17.537/)
      D5 = (/33.783,31.674,128.0,55.211/)
      D6 = (/-24.849,-16.936,-34.736,66.192/)
      BP00 = (/1.0151,.84946,.94016,.99926/)
      BP01 = (/-6.0574E-3,-9.7903E-3,-3.5957E-4,-5.0201E-3/)
      BP02 = (/5.5945E-5,1.0266E-4,9.8774E-6,4.8169E-5/)
      BP10 = (/-1.2901E-1,-2.0852E-1,1.2843E-1,-5.5311E-2/)
      BP11 = (/2.1565E-2,1.2935E-2,1.2117E-3,1.8072E-2/)
      BP12 = (/-1.95E-4,-9.4275E-5,-2.7557E-5,-1.693E-4/)
      BP20 = (/2.0622E-1,3.9371E-1,-1.4612E-1,9.0412E-2/)
      BP21 = (/-3.1109E-2,-2.3536E-2,-8.5631E-4,-2.3949E-2/)
      BP22 = (/2.8096E-4,1.8413E-4,2.7298E-5,2.2335E-4/)
      BP30 = (/-8.1528E-2,-1.3342E-1,3.9982E-2,-3.9868E-2/)
      BP31 = (/1.0582E-2,7.301E-3,3.7258E-4,7.5484E-3/)
      BP32 = (/-9.5007E-5,-5.7236E-5,-9.5415E-6,-6.9475E-5/)
      BQ00 = (/-3.0306,7.5308,-3.7748,-4.4981/)
      BQ01 = (/.12324,-.15526,.13631,.17798/)
      BQ02 = (/-6.408E-4,1.0762E-3,-7.6824E-4,-9.9386E-4/)
      BQ10 = (/1.0949,-.88621,1.5129,-5.0756/)
      BQ11 = (/5.4308E-3,-7.2508E-2,1.5867E-2,.13536/)
      BQ12 = (/1.7654E-5,9.8766E-4,-1.2999E-4,-6.7061E-4/)
      BQ20 = (/2.5572,2.2092,2.8725,6.6072/)
      BQ21 = (/7.2117E-3,2.9849E-2,2.6098E-3,-8.1503E-2/)
      BQ22 = (/-2.5712E-5,-2.2029E-4,-9.2133E-6,4.5423E-4/)
      AG00 = (/.75831,.65473,.77681,.77544/)
      AG01 = (/9.5376E-4,6.0975E-3,-2.7558E-3,-3.1632E-3/)
      AG02 = (/-2.3126E-6,-4.3907E-5,3.635E-5,3.577E-5/)
      AG10 = (/6.5007E-2,1.0582E-2,-3.07E-1,-2.3927E-3/)
      AG11 = (/-1.9238E-2,-2.0473E-2,5.5554E-3,-3.8837E-3/)
      AG12 = (/1.6785E-4,1.9499E-4,-4.014E-5,2.8519E-5/)
      AG20 = (/-2.5092E-2,7.2283E-2,1.1744E-1,-9.6464E-3/)
      AG21 = (/1.5397E-2,1.3209E-2,3.7471E-4,5.8684E-4/)
      AG22 = (/-1.3813E-4,-1.3393E-4,-1.5242E-6,-4.3942E-6/)
      AG30 = (/-4.7607E-4,-3.3056E-2,-7.4695E-3,0./)
      AG31 = (/-4.0963E-3,-3.0744E-3,-1.0596E-3,0./)
      AG32 = (/3.6814E-5,3.191E-5,6.5979E-6,0./)
      AG40 = (/7.4163E-4,3.6485E-3,-1.381E-3,0./)
      AG41 = (/3.5332E-4,2.4708E-4,1.7037E-4,0./)
      AG42 = (/-3.146E-6,-2.544E-6,-1.0431E-6,0./)
      Zref50 = (/0.,0.5,1.,1.5,2.,2.5,3.,3.5,4.,5.,6.,7.,8.,9.,10.,11.,12.,13., &
         14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,27.5,30.,32.5,35., &
         37.5,40.,42.5,45.,47.5,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100./)
      Zref41 = (/0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17., &
         18.,19.,20.,21.,22.,23.,24.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70., &
         75.,80.,85.,90.,95.,100./)
      O3a0 = (/-128.24,-118.47,-99.894,-99.881,-102.06,-106.79,-112.54,-115.88, &
         -117.69,-115.02,-91.984,-18.072,68.685,7.1018,-83.208,-145.22,-165.41, &
         -166.81,-145.86,-113.39,-74.301,-26.34,23.488,54.405,60.253,51.106, &
         7.4044,-16.679,-21.216,-43.058,-44.575,-55.374,-99.21,44.285,43.945, &
         32.81,-19.095,-16.587,-9.3195,-6.7324,-10.863/)
      O3a1 = (/0.65139,0.61388,0.54085,0.53352,0.5325,0.54076,0.55315,0.55586, &
         0.55201,0.52965,0.41801,0.078528,-0.32603,-0.058176,0.34991,0.6326, &
         0.72425,0.72921,0.63055,0.47874,0.2971,0.075867,-0.15234,-0.29358, &
         -0.31898,-0.27613,-0.073547,0.034212,0.060706,0.1574,0.18328,0.23388, &
         0.42174,-0.15026,-0.14293,-0.096462,0.13032,0.087694,0.029694,0.0187, &
         0.028993/)
      O3Min = (/0.2434,0.2416,0.2389,0.2364,0.2339,0.2319,0.2297,0.2276,0.2259, &
         0.2247,0.2245,0.2248,0.2243,0.2232,0.2236,0.2239,0.226,0.2279,0.2256, &
         0.211,0.1866,0.1607,0.1362,0.1145,0.0961,0.0807,0.0409,0.0201,7.99E-03, &
         2.63E-03,8.89E-04,3.36E-04,9.51E-05,3.31E-05,1.37E-05,7.50E-06,4.38E-06, &
         2.37E-06,8.94E-07,2.59E-07,8.00E-08/)
      O3Max = (/0.4444,0.4428,0.4406,0.4383,0.4359,0.4336,0.4308,0.4278,0.4235, &
         0.4182,0.4092,0.3979,0.3856,0.3727,0.3552,0.336,0.3132,0.289,0.2669, &
         0.2537,0.2464,0.239,0.2302,0.219,0.2048,0.1887,0.113,0.0482,1.69E-02, &
         5.43E-03,1.91E-03,8.24E-04,3.03E-04,9.99E-05,2.93E-05,1.04E-05,5.92E-06, &
         4.39E-06,1.95E-06,6.75E-07,1.50E-07/)
      TO3Min = (/211.,210.8,210.5,210.3,210.1,209.9,209.7,209.6,209.4,209.3, &
         209.2,209.1,208.9,208.8,208.7,208.6,208.5,208.4,208.4,208.4,208.6,208.9, &
         209.2,209.6,210.1,210.7,214.8,223.1,234.,244.9,250.3,246.5,237.7,218.8, &
         188.2,163.1,151.6,149.2,154.5,166.6,173./)
      TO3Max = (/241.4,241.1,240.5,240.,239.6,239.3,238.9,238.7,238.5,238.3, &
         238.3,238.4,238.6,238.8,239.1,239.5,240.,240.5,241.1,241.8,242.7,243.8, &
         244.9,246.3,247.8,249.3,257.2,267.5,277.3,282.3,277.3,268.1,253.1,241.5, &
         235.6,226.9,219.7,216.8,213.5,218.4,230./)
      Tref(:,1) = (/288.15,284.90,281.65,278.40,275.15,271.91,268.66,265.41, &
         262.17,255.68,249.19,242.70,236.22,229.73,223.25,216.77,216.65, &
         216.65,216.65,216.65,216.65,216.65,216.65,216.65,216.65,217.58, &
         218.57,219.57,220.56,221.55,224.03,226.51,229.58,236.51,243.44, &
         250.35,257.25,264.16,270.65,270.65,260.77,247.02,233.29,219.59, &
         208.40,198.64,188.89,186.87,188.42,195.08/)
      Tref(:,2) = (/294.15,291.90,289.65,287.40,285.15,282.15,279.15,276.15, &
         273.15,267.15,261.15,254.65,248.15,241.65,235.15,228.80,222.30, &
         215.80,215.70,215.70,215.70,215.70,216.80,217.90,219.20,220.40, &
         221.60,222.80,223.90,225.10,228.50,233.70,239.00,245.20,251.30, &
         257.50,263.70,269.90,275.20,275.70,269.30,257.10,240.10,218.10, &
         196.10,174.10,165.10,165.00,178.30,190.50/)
      Tref(:,3) = (/272.15,270.40,268.65,266.90,265.15,263.40,261.65,258.65, &
         255.65,249.65,243.65,237.65,231.65,225.65,219.65,219.20,218.70, &
         218.20,217.70,217.20,216.70,216.20,215.70,215.20,215.20,215.20, &
         215.20,215.20,215.20,215.20,215.50,217.40,220.40,227.90,235.50, &
         243.20,250.80,258.50,265.10,265.70,260.60,250.80,240.90,230.70, &
         220.40,210.10,199.80,199.50,208.30,218.60/)
      Tref(:,4) = (/287.15,284.45,281.75,279.05,276.35,273.65,270.95,268.25, &
         265.55,260.15,253.15,246.15,239.15,232.15,225.15,225.15,225.15, &
         225.15,225.15,225.15,225.15,225.15,225.15,225.15,225.15,225.15, &
         225.15,225.15,226.60,228.10,231.00,235.10,240.00,247.20,254.60, &
         262.10,269.50,273.60,276.20,277.20,274.00,262.70,239.70,216.60, &
         193.60,170.60,161.70,161.60,176.80,190.40/)
      Tref(:,5) = (/257.15,258.15,259.15,257.55,255.95,254.35,252.75,251.15, &
         247.75,240.95,234.15,227.30,220.55,217.15,217.15,217.15,217.15, &
         217.15,217.15,217.15,216.60,216.00,215.40,214.80,214.20,213.60, &
         213.00,212.40,211.80,211.20,213.60,216.00,218.50,222.30,228.50, &
         234.70,240.80,247.00,253.20,259.30,259.10,250.90,248.40,245.40, &
         234.70,223.90,213.10,202.30,211.00,218.50/)
      Tref(:,6) = (/299.65,296.65,293.65,290.65,287.65,286.95,283.65,280.32, &
         277.00,270.30,263.60,257.00,250.30,243.60,237.00,230.10,223.60, &
         217.00,210.30,203.70,197.00,194.80,198.80,202.70,206.70,210.70, &
         214.60,217.00,219.20,221.40,227.00,232.30,237.70,243.10,248.50, &
         254.00,259.40,264.80,269.60,270.20,263.40,253.10,236.00,218.90, &
         201.80,184.80,177.10,177.00,184.30,190.70/)
      Tref(:,7) = (/301.15,297.40,293.65,290.90,288.15,285.40,282.65,279.90, &
         277.15,271.65,266.15,259.15,252.15,245.15,238.15,231.15,224.15, &
         217.15,210.15,203.50,203.15,205.20,207.38,209.57,211.75,213.93, &
         215.94,217.92,219.90,221.88,226.85,231.80,237.70,243.10,248.50, &
         254.00,259.40,264.80,269.60,270.20,263.40,253.10,236.00,218.90, &
         201.80,184.80,177.10,177.00,184.30,190.70/)
      Tref(:,8) = (/287.15,285.65,284.15,282.65,281.15,277.90,274.65,271.40, &
         268.15,261.65,255.15,248.65,242.15,235.65,229.15,222.65,216.15, &
         213.55,210.95,208.35,205.75,203.15,203.15,205.65,208.15,210.65, &
         213.15,215.15,217.15,219.15,224.15,229.15,234.48,239.80,245.13, &
         250.56,255.88,261.21,265.94,266.54,259.83,249.67,232.80,215.93, &
         199.06,182.29,174.70,174.60,181.80,188.11/)
      Tref(:,9) = (/278.15,276.85,275.55,274.25,272.95,271.65,268.40,265.15, &
         261.45,255.51,248.90,242.40,235.90,229.36,226.66,227.66,228.65, &
         229.65,230.15,230.15,230.15,230.15,230.15,230.15,230.15,230.15, &
         230.15,230.15,230.71,231.90,234.88,237.86,240.00,247.20,254.60, &
         262.10,269.50,273.60,276.20,277.20,274.00,262.70,239.70,216.60, &
         193.60,170.60,161.70,161.60,176.80,190.40/)
      Tref(:,10) = (/249.15,250.65,252.15,253.65,250.90,248.15,245.40, &
         242.65,239.90,234.38,228.87,223.36,217.86,214.90,214.40,213.90, &
         213.25,212.45,211.65,210.85,210.05,209.26,208.46,207.66,207.65, &
         207.65,207.65,207.65,207.65,207.65,207.65,207.65,210.50,214.30, &
         220.50,226.70,232.80,239.00,245.20,251.30,251.10,242.90,240.40, &
         237.40,226.70,215.90,205.10,194.30,203.00,210.50/)
      Pref(:,1) = (/1013.25,954.61,898.76,845.59,795.01,746.91,701.21, &
         657.8,616.6,540.48,472.17,411.05,356.51,308.,264.99,226.99, &
         193.99,165.79,141.7,121.11,103.52,88.5,75.65,64.67,55.29,47.29, &
         40.47,34.67,29.72,25.49,17.43,11.97,8.01,5.746,4.15,2.871,2.06, &
         1.491,1.09,0.7978,0.425,0.219,0.109,0.0522,0.024,0.0105, &
         0.00446,0.00184,0.00076,0.00032/)
      Pref(:,2) = (/1013.5,956.5,902.2,850.6,801.6,754.9,710.5,668.2, &
         628.,553.6,486.6,426.4,372.4,324.,280.9,242.6,208.6,178.6, &
         152.5,130.3,111.3,95.0,81.2,69.5,59.5,51.0,43.7,37.6,32.2,27.7, &
         19.07,13.2,9.3,6.52,4.64,3.33,2.41,1.76,1.29,0.951,0.515, &
         0.272,0.139,0.067,0.03,0.012,0.00448,0.00164,0.000625,0.000258/)
      Pref(:,3) = (/1018.,956.,897.4,842.,789.7,740.4,693.8,649.8, &
         608.1,531.3,462.7,401.6,347.3,299.3,256.8,219.9,188.2,161.1, &
         137.8,117.8,100.7,86.1,73.6,62.8,53.7,45.8,39.1,33.4,28.6, &
         24.4,16.46,11.1,7.56,5.18,3.6,2.53,1.8,1.29,0.94,0.683,0.362, &
         0.188,0.095,0.047,0.0222,0.0103,0.00456,0.00198,0.000877, &
         0.0004074/)
      Pref(:,4) = (/1010.,951.6,896.,843.2,793.,745.3,700.,657.1,616.4, &
         541.4,474.,413.4,359.2,310.8,267.7,230.1,197.8,170.,146.1, &
         125.6,108.0,92.85,79.83,68.64,59.0,50.7,43.6,37.5,32.28,27.8, &
         19.23,13.4,9.4,6.61,4.72,3.4,2.48,1.82,1.34,0.987,0.537, &
         0.288,0.147,0.071,0.032,0.0125,0.00451,0.00161,0.000606, &
         0.000248/)
      Pref(:,5) = (/1013.5,948.4,887.8,831.,777.5,727.1,679.8,635.2, &
         593.2,515.8,446.7,385.3,330.8,282.9,241.8,206.7,176.6,151., &
         129.1,110.3,94.31,80.58,68.82,58.75,50.14,42.77,36.47,31.09, &
         26.49,22.56,15.13,10.2,6.91,4.701,3.23,2.243,1.57,1.113,0.79, &
         0.5719,0.299,0.155,0.079,0.04,0.02,0.00966,0.0045,0.002022, &
         0.000907,0.000423/)
      Pref(:,6) = (/1013.25,957.5,904.2,853.3,804.8,758.6,714.8,673., &
         633.2,559.2,492.4,432.1,377.9,329.3,285.9,247.2,212.8,182.4, &
         155.6,132.1,111.5,93.7,78.9,66.6,56.5,48.,40.9,35.,30.,25.7, &
         17.63,12.2,8.52,6.0,4.26,3.05,2.2,1.59,1.16,0.854,0.456,0.239, &
         0.121,0.058,0.026,0.011,0.0044,0.00172,0.000688,0.000289/)
      Pref(:,7) = (/1013.5,957.9,904.6,853.6,805.1,758.8,714.8,672.9, &
         633.1,559.3,492.9,433.1,379.1,330.7,287.3,248.6,214.1,183.6, &
         156.6,132.9,112.5,95.3,80.8,68.7,58.5,49.9,42.6,36.4,31.2, &
         26.8,18.45,14.,8.52,6.,4.26,3.05,2.2,1.59,1.16,0.854,0.456, &
         0.239,0.121,0.058,0.026,0.011,0.0044,0.00172,0.000688,0.000289/)
      Pref(:,8) = (/1021.,962.1,906.4,853.5,803.5,756.,710.7,667.7, &
         626.8,551.,482.8,421.6,366.8,317.9,274.4,235.9,201.9,172.2, &
         146.6,124.6,105.6,89.4,75.5,63.9,54.2,46.,39.2,33.4,28.5, &
         24.4,16.65,12.6,8.52,6.,4.26,3.05,2.2,1.59,1.16,0.854,0.456, &
         0.239,0.121,0.058,0.026,0.011,0.0044,0.00172,0.000688,0.000289/)
      Pref(:,9) = (/1012.5,952.1,895.,841.1,790.2,742.2,696.7,653.5, &
         612.5,536.7,468.6,407.8,353.5,305.2,262.6,226.,194.6,167.7, &
         144.6,124.7,107.5,92.7,80.,69.,59.5,51.3,44.3,38.2,33.,28.5, &
         19.85,13.8,9.4,6.61,4.72,3.4,2.48,1.82,1.34,0.987,0.537,0.288, &
         0.147,0.071,0.032,0.0125,0.00451,0.00161,0.000606,0.000248/)
      Pref(:,10) = (/1013.5,946.4,884.1,826.3,772.1,721.,672.7,627.2, &
         584.3,505.8,436.4,375.2,321.4,274.3,234.,199.5,170.1,144.9, &
         123.4,105.,89.3,75.9,64.49,54.75,46.48,39.45,33.49,28.43, &
         24.14,20.5,13.67,9.05,6.91,4.701,3.23,2.243,1.57,1.113,0.79, &
         0.5719,0.299,0.155,0.079,0.04,0.02,0.00966,0.0045,0.002022, &
         0.000907,0.000423/)
      RHref(:,1) = (/46.04,47.61,49.18,50.67,52.16,51.53,50.91,50.56, &
         50.21,48.49,49.28,48.17,50.50,36.98,28.87,27.54,12.61,6.134, &
         2.864,2.065,1.394,1.162,0.987,0.849,0.735,0.572,0.444,0.349, &
         0.272,0.214,1.38E-01,6.10E-02,3.58E-02,1.06E-02,6.04E-03, &
         1.47E-03,8.61E-04,2.51E-04,5.36E-04,8.20E-04,9.15E-05, &
         1.43E-04,2.90E-04,4.37E-04,6.43E-04,8.48E-04,6.38E-04, &
         4.28E-04,2.19E-04,8.68E-06/)
      RHref(:,2) = (/76.18,71.10,66.02,60.61,55.20,50.25,45.29,42.17, &
         39.05,31.42,29.98,30.31,29.63,30.15,29.44,19.48,10.69,5.423, &
         2.933,1.695,1.404,1.166,0.856,0.651,0.492,0.382,0.297,0.238, &
         0.186,0.147,8.90E-02,3.10E-02,1.81E-02,5.23E-03,3.08E-03, &
         9.37E-04,5.68E-04,1.99E-04,1.35E-04,7.13E-05,5.99E-05, &
         7.75E-05,3.93E-04,7.08E-04,3.66E-02,7.24E-02,5.43E-02, &
         3.62E-02,1.81E-02,1.46E-05/)
      RHref(:,3) = (/77.07,73.79,70.50,67.97,65.45,61.10,56.74,53.32, &
         49.89,47.17,44.02,31.03,23.01,19.66,17.93,5.505,3.001,2.274, &
         1.984,1.765,1.57,1.396,1.27,1.153,0.986,0.841,0.723,0.62, &
         0.537,0.463,3.13E-01,1.64E-01,9.36E-02,2.33E-02,1.29E-02, &
         2.45E-03,1.39E-03,3.28E-04,2.12E-04,9.67E-05,7.51E-05, &
         8.27E-05,9.44E-05,1.06E-04,1.34E-04,1.62E-04,1.22E-04, &
         8.12E-05,4.08E-05,4.38E-07/)
      RHref(:,4) = (/75.23,72.64,70.04,69.96,69.89,67.51,65.12,62.80, &
         60.47,53.46,50.43,49.11,41.19,23.59,14.17,3.819,1.481,0.944, &
         0.729,0.629,0.539,0.469,0.428,0.385,0.339,0.298,0.261,0.226, &
         0.167,0.123,7.60E-02,2.90E-02,1.67E-02,4.46E-03,2.55E-03, &
         6.45E-04,3.95E-04,1.44E-04,1.02E-04,5.99E-05,4.01E-05, &
         4.69E-05,4.25E-04,8.04E-04,7.49E-02,1.49E-01,1.12E-01, &
         7.45E-02,3.73E-02,1.43E-05/)
      RHref(:,5) = (/80.46,74.88,69.30,69.62,69.93,67.79,65.64,63.02, &
         60.40,54.09,50.71,55.97,23.73,26.84,15.42,6.589,3.378,2.142, &
         1.852,1.6,1.488,1.384,1.287,1.197,1.113,1.035,0.963,0.895, &
         0.832,0.766,4.77E-01,1.88E-01,1.14E-01,4.09E-02,2.30E-02, &
         5.05E-03,2.91E-03,7.65E-04,4.50E-04,1.35E-04,7.01E-05, &
         6.76E-05,4.43E-05,2.10E-05,2.45E-05,2.79E-05,2.10E-05, &
         1.42E-05,7.32E-06,4.61E-07/)
      RHref(:,6) = (/75.59,74.21,72.83,73.71,74.58,61.44,48.29,41.61, &
         34.94,37.74,34.82,32.01,29.49,25.37,19.52,13.16,9.26,5.886, &
         7.416,9.931,16.78,19.25,8.328,3.755,1.82,0.923,0.504,0.332, &
         0.24,0.161,9.46E-02,2.82E-02,1.68E-02,5.45E-03,3.31E-03, &
         1.17E-03,7.24E-04,2.78E-04,1.91E-04,1.04E-04,9.36E-05, &
         1.15E-04,3.96E-04,6.77E-04,4.68E-03,8.68E-03,6.51E-03, &
         4.35E-03,2.18E-03,1.58E-05/)
      RHref(:,7) = (/80.00,70.00,65.00,62.50,60.00,60.00,60.00,55.00, &
         50.00,45.00,40.00,40.00,40.00,35.00,30.00,15.69,9.83,5.701, &
         5.623,6.637,10.630,12.016,5.339,2.513,1.289,0.707,0.421, &
         0.294,0.218,0.155,9.24E-02,2.93E-02,1.73E-02,5.36E-03, &
         3.22E-03,1.08E-03,6.62E-04,2.46E-04,1.69E-04,9.09E-05, &
         8.01E-05,1.00E-04,3.95E-04,6.89E-04,1.74E-02,3.42E-02, &
         2.56E-02,1.71E-02,8.55E-03,1.53E-05/)
      RHref(:,8) = (/80.00,75.00,70.00,60.00,50.00,47.50,45.00,40.00, &
         35.00,32.50,30.00,30.00,30.00,30.00,30.00,13.16,9.26,5.886, &
         7.416,9.931,16.780,19.250,8.328,3.755,1.820,0.923,0.504, &
         0.332,0.240,0.161,9.46E-02,2.82E-02,1.68E-02,5.45E-03, &
         3.31E-03,1.17E-03,7.24E-04,2.78E-04,1.91E-04,1.04E-04, &
         9.36E-05,1.15E-04,3.96E-04,6.77E-04,4.68E-03,8.68E-03, &
         6.51E-03,4.35E-03,2.18E-03,1.58E-05/)
      RHref(:,9) = (/85.00,80.00,75.00,70.00,65.00,65.00,60.00,57.50, &
         55.00,50.00,45.00,40.00,35.00,32.00,20.00,9.22,4.73,2.999, &
         2.593,2.240,2.083,1.938,1.802,1.676,1.558,1.449,1.348,1.253, &
         1.165,1.072,6.68E-01,2.63E-01,1.60E-01,5.73E-02,3.22E-02, &
         7.07E-03,4.07E-03,1.07E-03,6.30E-04,1.89E-04,9.81E-05, &
         9.46E-05,6.20E-05,2.94E-05,3.42E-05,3.91E-05,2.95E-05, &
         1.99E-05,1.02E-05,6.45E-07/)
      RHref(:,10) = (/80.00,70.00,65.00,60.00,60.00,57.50,55.00, &
         52.50,50.00,47.50,45.00,42.50,40.00,21.23,12.75,3.44,1.33, &
         0.850,0.656,0.566,0.485,0.422,0.385,0.347,0.305,0.268,0.235, &
         0.203,0.150,0.111,6.84E-02,2.61E-02,1.51E-02,4.01E-03, &
         2.30E-03,5.81E-04,3.55E-04,1.30E-04,9.18E-05,5.39E-05, &
         3.61E-05,4.22E-05,3.83E-04,7.24E-04,6.74E-02,1.34E-01, &
         1.01E-01,6.71E-02,3.35E-02,1.29E-05/)
      O3ref(:,1) = (/0.3438,0.3425,0.3413,0.3400,0.3387,0.3375, &
         0.3363,0.3352,0.3341,0.3319,0.3298,0.3276,0.3253,0.3224, &
         0.3186,0.3135,0.3067,0.2990,0.2906,0.2819,0.2708,0.2586, &
         0.2446,0.2290,0.2120,0.1942,0.1762,0.1583,0.1410,0.1247, &
         8.981E-02,6.308E-02,4.291E-02,2.795E-02,1.743E-02,1.045E-02, &
         6.049E-03,3.430E-03,1.921E-03,1.082E-03,3.459E-04,1.059E-04, &
         3.238E-05,1.253E-05,5.770E-06,2.578E-06,1.101E-06,4.354E-07, &
         1.353E-07,3.028E-07/)
      O3ref(:,2) = (/0.3320,0.3306,0.3292,0.3278,0.3263,0.3249, &
         0.3235,0.3220,0.3206,0.3175,0.3144,0.3110,0.3074,0.3036, &
         0.2995,0.2948,0.2894,0.2832,0.2757,0.2674,0.2592,0.2509, &
         0.2417,0.2302,0.2164,0.2016,0.1864,0.1709,0.1555,0.1397, &
         1.037E-01,7.417E-02,5.066E-02,3.248E-02,1.959E-02,1.132E-02, &
         6.370E-03,3.659E-03,2.142E-03,1.271E-03,4.590E-04,1.554E-04, &
         4.782E-05,1.565E-05,5.579E-06,2.206E-06,9.322E-07,3.650E-07, &
         1.124E-07,4.303E-07/)
      O3ref(:,3) = (/0.3771,0.3757,0.3744,0.3732,0.3720,0.3709, &
         0.3697,0.3686,0.3674,0.3649,0.3621,0.3588,0.3549,0.3500, &
         0.3435,0.3348,0.3239,0.3108,0.2970,0.2836,0.2701,0.2557, &
         0.2399,0.2226,0.2038,0.1840,0.1644,0.1459,0.1284,0.1122, &
         7.872E-02,5.412E-02,3.623E-02,2.335E-02,1.441E-02,8.445E-03, &
         4.757E-03,2.672E-03,1.503E-03,8.540E-04,2.750E-04,8.590E-05, &
         2.758E-05,1.037E-05,4.523E-06,2.090E-06,9.869E-07,4.275E-07, &
         1.438E-07,2.633E-07/)
      O3ref(:,4) = (/0.3451,0.3439,0.3427,0.3414,0.3402,0.3388, &
         0.3375,0.3361,0.3347,0.3318,0.3287,0.3253,0.3217,0.3173, &
         0.3117,0.3045,0.2954,0.2855,0.2751,0.2646,0.2539,0.2428, &
         0.2310,0.2178,0.2035,0.1878,0.1710,0.1541,0.1378,0.1223, &
         8.979E-02,6.476E-02,4.518E-02,2.959E-02,1.821E-02,1.064E-02, &
         6.019E-03,3.454E-03,2.020E-03,1.211E-03,4.470E-04,1.540E-04, &
         4.918E-05,1.600E-05,5.447E-06,2.110E-06,8.941E-07,3.504E-07, &
         1.081E-07,4.893E-07/)
      O3ref(:,5) = (/0.3759,0.3750,0.3740,0.3730,0.3721,0.3711, &
         0.3701,0.3691,0.3681,0.3659,0.3637,0.3609,0.3571,0.3513, &
         0.3432,0.3342,0.3253,0.3148,0.3015,0.2861,0.2691,0.2507, &
         0.2306,0.2086,0.1855,0.1631,0.1426,0.1240,0.1074,0.0929, &
         6.423E-02,4.378E-02,2.916E-02,1.867E-02,1.146E-02,6.736E-03, &
         3.824E-03,2.156E-03,1.214E-03,6.885E-04,2.260E-04,7.800E-05, &
         2.770E-05,8.960E-06,2.905E-06,1.320E-05,7.201E-07,3.566E-07, &
         1.351E-07,2.399E-07/)
      O3ref(:,6) = (/0.2776,0.2763,0.2750,0.2737,0.2724,0.2712, &
         0.2700,0.2688,0.2677,0.2655,0.2635,0.2615,0.2597,0.2578, &
         0.2560,0.2542,0.2522,0.2501,0.2480,0.2459,0.2437,0.2410, &
         0.2367,0.2298,0.2204,0.2098,0.1981,0.1845,0.1692,0.1528, &
         1.131E-01,7.803E-02,5.034E-02,3.095E-02,1.819E-02,1.036E-02, &
         5.802E-03,3.305E-03,1.909E-03,1.110E-03,3.700E-04,1.151E-04, &
         3.629E-05,1.449E-05,6.783E-06,2.951E-06,1.182E-06,4.409E-07, &
         1.302E-07,3.135E-07/)
      O3ref(:,7) = (/0.3104,0.3090,0.3077,0.3064,0.3050,0.3037, &
         0.3024,0.3011,0.2999,0.2980,0.2955,0.2930,0.2905,0.2879, &
         0.2852,0.2824,0.2792,0.2757,0.2717,0.2669,0.2621,0.2572, &
         0.2518,0.2445,0.2341,0.2219,0.2091,0.1953,0.1800,0.1640, &
         1.346E-01,9.576E-02,6.165E-02,3.954E-02,2.394E-02,1.382E-02, &
         7.734E-03,4.362E-03,2.541E-03,1.534E-03,6.398E-04,2.198E-04, &
         6.960E-05,2.053E-05,8.388E-06,5.192E-06,2.962E-06,1.194E-06, &
         3.920E-07,1.001E-07/)
      O3ref(:,8) = (/0.2798,0.2784,0.2770,0.2757,0.2744,0.2731, &
         0.2718,0.2706,0.2694,0.2678,0.2656,0.2636,0.2616,0.2598, &
         0.2580,0.2562,0.2543,0.2523,0.2503,0.2484,0.2463,0.2443, &
         0.2414,0.2364,0.2284,0.2186,0.2080,0.1961,0.1818,0.1666, &
         1.383E-01,9.922E-02,6.476E-02,4.063E-02,2.428E-02,1.397E-02, &
         7.819E-03,4.400E-03,2.554E-03,1.540E-03,6.327E-04,2.069E-04, &
         6.500E-05,1.946E-05,8.597E-06,5.428E-06,2.744E-06,1.047E-06, &
         3.565E-07,1.015E-07/)
      O3ref(:,9) = (/0.3296,0.3287,0.3278,0.3268,0.3259,0.3249, &
         0.3239,0.3229,0.3218,0.3203,0.3180,0.3156,0.3127,0.3093, &
         0.3039,0.2974,0.2899,0.2819,0.2722,0.2614,0.2496,0.2371, &
         0.2237,0.2088,0.1924,0.1752,0.1581,0.1414,0.1257,0.1111, &
         8.823E-02,6.321E-02,4.452E-02,2.994E-02,1.911E-02,1.153E-02, &
         6.658E-03,3.854E-03,2.253E-03,1.380E-03,5.915E-04,2.161E-04, &
         7.877E-05,2.699E-05,9.977E-06,4.957E-06,3.631E-06,1.375E-06, &
         3.911E-07,7.327E-08/)
      O3ref(:,10) = (/0.4496,0.4484,0.4472,0.4459,0.4447,0.4435, &
         0.4423,0.4410,0.4397,0.4377,0.4349,0.4321,0.4280,0.4228, &
         0.4137,0.4027,0.3917,0.3810,0.3644,0.3456,0.3207,0.2926, &
         0.2623,0.2315,0.2006,0.1731,0.1493,0.1289,0.1113,0.0963, &
         7.345E-02,5.127E-02,3.600E-02,2.339E-02,1.466E-02,8.741E-03, &
         5.074E-03,2.938E-03,1.714E-03,1.041E-03,4.194E-04,1.305E-04, &
         4.881E-05,2.002E-05,8.666E-06,4.742E-06,3.958E-06,1.740E-06, &
         6.179E-07,1.362E-07/)
      Wref(:,1) = (/1.416E+00,1.145E+00,9.163E-01,7.247E-01,5.655E-01, &
         4.365E-01,3.349E-01,2.552E-01,1.928E-01,1.079E-01,5.802E-02, &
         2.937E-02,1.329E-02,5.580E-03,2.596E-03,1.350E-03,7.853E-04, &
         5.217E-04,3.958E-04,3.242E-04,2.741E-04,2.367E-04,2.053E-04, &
         1.784E-04,1.552E-04,1.351E-04,1.176E-04,1.023E-04,8.884E-05, &
         7.707E-05,5.387E-05,3.762E-05,2.631E-05,1.851E-05,1.310E-05, &
         9.338E-06,6.685E-06,4.788E-06,3.424E-06,2.438E-06,1.193E-06, &
         5.385E-07,2.252E-07,8.616E-08,2.959E-08,9.154E-09,2.513E-09, &
         6.605E-10,1.440E-10,2.737E-11/)
      Wref(:,2) = (/2.922E+00,2.289E+00,1.774E+00,1.358E+00,1.027E+00, &
         7.709E-01,5.796E-01,4.355E-01,3.262E-01,1.860E-01,1.071E-01, &
         5.916E-02,3.093E-02,1.485E-02,5.951E-03,2.020E-03,7.893E-04, &
         4.701E-04,3.632E-04,3.041E-04,2.636E-04,2.300E-04,2.020E-04, &
         1.782E-04,1.575E-04,1.392E-04,1.229E-04,1.081E-04,9.490E-05, &
         8.307E-05,5.926E-05,4.220E-05,3.008E-05,2.153E-05,1.550E-05, &
         1.123E-05,8.161E-06,5.905E-06,4.253E-06,3.046E-06,1.512E-06, &
         7.007E-07,3.008E-07,1.155E-07,3.861E-08,1.118E-08,2.734E-09, &
         6.447E-10,1.278E-10,2.429E-11/)
      Wref(:,3) = (/8.517E-01,6.907E-01,5.547E-01,4.394E-01,3.417E-01, &
         2.603E-01,1.938E-01,1.419E-01,1.035E-01,5.281E-02,2.416E-02, &
         1.034E-02,4.709E-03,2.283E-03,1.162E-03,7.322E-04,5.733E-04, &
         4.783E-04,4.056E-04,3.452E-04,2.946E-04,2.521E-04,2.162E-04, &
         1.854E-04,1.590E-04,1.365E-04,1.173E-04,1.007E-04,8.649E-05, &
         7.418E-05,5.051E-05,3.446E-05,2.360E-05,1.628E-05,1.131E-05, &
         7.919E-06,5.579E-06,3.950E-06,2.802E-06,1.983E-06,9.629E-07, &
         4.361E-07,1.853E-07,7.327E-08,2.650E-08,8.642E-09,2.505E-09, &
         6.931E-09,1.580E-10,3.002E-11/)
      Wref(:,4) = (/2.081E+00,1.671E+00,1.337E+00,1.063E+00,8.329E-01, &
         6.446E-01,4.936E-01,3.732E-01,2.776E-01,1.457E-01,7.110E-02, &
         3.092E-02,1.098E-02,3.529E-03,1.368E-03,7.591E-04,5.689E-04, &
         4.770E-04,4.129E-04,3.606E-04,3.157E-04,2.769E-04,2.423E-04, &
         2.110E-04,1.832E-04,1.587E-04,1.372E-04,1.184E-04,1.022E-04, &
         8.832E-05,6.128E-05,4.260E-05,2.975E-05,2.095E-05,1.487E-05, &
         1.064E-05,7.649E-06,5.510E-06,3.966E-06,2.842E-06,1.418E-06, &
         6.625E-07,2.880E-07,1.126E-07,3.849E-08,1.122E-08,2.710E-09, &
         6.314E-10,1.239E-10,2.354E-11/)
      Wref(:,5) = (/4.161E-01,3.562E-01,2.962E-01,2.398E-01,1.898E-01, &
         1.465E-01,1.096E-01,7.955E-02,5.624E-02,2.699E-02,1.270E-02, &
         5.324E-03,2.622E-03,1.658E-03,1.013E-03,6.885E-04,5.381E-04, &
         4.532E-04,3.908E-04,3.369E-04,2.902E-04,2.498E-04,2.148E-04, &
         1.845E-04,1.583E-04,1.357E-04,1.161E-04,9.920E-05,8.460E-05, &
         7.207E-05,4.831E-05,3.252E-05,2.198E-05,1.494E-05,1.021E-05, &
         7.031E-06,4.874E-06,3.393E-06,2.369E-06,1.656E-06,7.932E-07, &
         3.571E-07,1.517E-07,6.133E-08,2.316E-08,7.862E-09,2.371E-09, &
         6.804E-10,1.601E-10,3.042E-11/)
      Wref(:,6) = (/4.115E+00,3.250E+00,2.534E+00,1.936E+00,1.430E+00, &
         1.036E+00,7.557E-01,5.603E-01,4.266E-01,2.439E-01,1.295E-01, &
         6.546E-02,3.064E-02,1.293E-02,4.942E-03,1.884E-03,8.286E-04, &
         4.799E-04,3.438E-04,2.679E-04,2.222E-04,1.889E-04,1.620E-04, &
         1.411E-04,1.242E-04,1.100E-04,9.765E-05,8.682E-05,7.701E-05, &
         6.820E-05,5.048E-05,3.726E-05,2.746E-05,2.026E-05,1.495E-05, &
         1.101E-05,8.085E-06,5.922E-06,4.317E-06,3.122E-06,1.578E-06, &
         7.346E-07,3.082E-07,1.126E-07,3.520E-08,1.012E-08,2.633E-09, &
         6.586E-10,1.375E-10,2.613E-11/)
      Wref(:,7) = (/4.229E+00,3.302E+00,2.630E+00,2.106E+00,1.679E+00, &
         1.326E+00,1.028E+00,7.877E-01,6.045E-01,3.504E-01,1.955E-01, &
         1.034E-01,4.968E-02,2.106E-02,7.606E-03,2.394E-03,9.766E-04, &
         5.661E-04,4.431E-04,3.886E-04,3.412E-04,2.699E-04,2.043E-04, &
         1.648E-04,1.393E-04,1.216E-04,1.086E-04,9.785E-05,8.824E-05, &
         7.954E-05,5.833E-05,4.666E-05,3.006E-05,2.156E-05,1.559E-05, &
         1.134E-05,8.339E-06,6.128E-06,4.506E-06,3.334E-06,2.024E-06, &
         1.311E-06,8.065E-07,5.496E-07,2.235E-07,5.156E-08,2.538E-08, &
         1.649E-08,4.038E-09,4.038E-09/)
      Wref(:,8) = (/2.096E+00,1.649E+00,1.269E+00,9.573E-01,7.184E-01, &
         5.303E-01,3.851E-01,2.822E-01,2.098E-01,1.171E-01,6.459E-02, &
         3.496E-02,1.794E-02,8.485E-03,3.434E-03,1.371E-03,8.733E-04, &
         6.892E-04,5.747E-04,4.675E-04,3.500E-04,2.365E-04,1.650E-04, &
         1.295E-04,1.065E-04,9.062E-05,7.937E-05,7.040E-05,6.252E-05, &
         5.568E-05,4.034E-05,3.214E-05,2.301E-05,1.690E-05,1.243E-05, &
         9.125E-06,6.685E-06,4.847E-06,3.462E-06,2.451E-06,1.309E-06, &
         6.930E-07,3.030E-07,1.209E-07,3.859E-08,8.343E-09,4.117E-09, &
         2.706E-09,6.839E-10,6.839E-10/)
      Wref(:,9) = (/1.484E+00,1.215E+00,9.843E-01,7.863E-01,6.177E-01, &
         4.695E-01,3.471E-01,2.562E-01,1.888E-01,1.007E-01,5.207E-02, &
         2.719E-02,1.495E-02,9.138E-03,6.382E-03,5.053E-03,4.343E-03, &
         3.902E-03,3.551E-03,3.239E-03,2.960E-03,2.700E-03,2.459E-03, &
         2.235E-03,2.026E-03,1.831E-03,1.652E-03,1.484E-03,1.325E-03, &
         1.163E-03,7.841E-04,5.275E-04,3.574E-04,2.442E-04,1.481E-04, &
         8.482E-05,4.661E-05,2.287E-05,1.181E-05,6.077E-06,2.165E-06, &
         5.152E-07,4.169E-08,1.743E-09,9.137E-11,9.077E-12,7.235E-12, &
         6.812E-12,3.317E-12,3.317E-12/)
      Wref(:,10) = (/2.164E-01,1.858E-01,1.545E-01,1.216E-01,9.159E-02, &
         6.824E-02,5.057E-02,3.734E-02,2.751E-02,1.459E-02,7.450E-03, &
         3.626E-03,1.650E-03,7.449E-04,3.564E-04,1.847E-04,1.372E-04, &
         1.170E-04,1.044E-04,9.512E-05,8.790E-05,8.225E-05,7.772E-05, &
         7.402E-05,7.089E-05,6.814E-05,6.573E-05,6.363E-05,6.194E-05, &
         6.070E-05,5.862E-05,5.755E-05,5.717E-05,5.684E-05,5.660E-05, &
         5.641E-05,5.628E-05,5.616E-05,5.604E-05,5.593E-05,5.573E-05, &
         5.561E-05,5.512E-05,5.424E-05,3.369E-05,1.103E-05,2.561E-06, &
         1.068E-06,2.766E-07,2.766E-07/)
      TO3(:,1) = (/225.3,225.2,225.1,224.8,224.6,224.4,224.3,224.1, &
         224.0,223.7,223.5,223.3,223.2,223.1,223.0,223.0,223.2,223.3, &
         223.5,223.7,224.0,224.3,224.7,225.2,225.8,226.6,227.4,228.4, &
         229.4,230.5,232.9,236.4,240.8,246.2,251.8,257.3,262.2,265.7, &
         266.9,264.2,254.0,240.3,226.3,209.8,199.5,192.7,188.6,188.2, &
         190.6,198.3/)
      TO3(:,2) = (/231.9,231.8,231.7,231.4,231.2,230.9,230.7,230.5, &
         230.3,229.9,229.5,229.2,228.9,228.7,228.5,228.4,228.4,228.5, &
         228.8,229.3,229.7,230.1,230.6,231.2,231.9,232.8,233.7,234.8, &
         235.9,237.1,239.8,244.1,248.5,253.3,258.3,263.3,267.9,271.0, &
         271.8,269.6,261.5,248.3,229.7,204.7,180.7,169.5,167.1,170.1, &
         182.2,196.9/)
      TO3(:,3) = (/220.4,220.4,220.3,220.1,219.9,219.8,219.7,219.5, &
         219.4,219.2,218.9,218.7,218.6,218.4,218.3,218.3,218.2,218.2, &
         218.2,218.2,218.3,218.4,218.5,218.7,219.0,219.4,219.9,220.4, &
         221.1,222.0,223.9,227.6,232.3,238.6,244.8,250.9,256.6,261.2, &
         263.2,261.8,255.8,245.4,233.6,221.2,211.0,204.5,202.0,204.2, &
         212.4,225.1/)
      TO3(:,4) = (/233.4,233.4,233.3,233.1,232.9,232.7,232.6,232.4, &
         232.3,232.0,231.7,231.5,231.4,231.3,231.3,231.4,231.6,231.8, &
         232.0,232.3,232.6,232.9,233.2,233.7,234.2,234.9,235.8,236.9, &
         238.2,239.6,242.4,246.7,251.3,256.7,262.3,267.6,271.9,273.6, &
         273.6,271.9,265.3,250.9,228.6,202.3,177.8,165.8,163.8,166.8, &
         180.7,197.5/)
      TO3(:,5) = (/217.3,217.2,217.2,217.1,217.0,216.9,216.8,216.7, &
         216.6,216.4,216.3,216.2,216.0,216.0,216.0,215.9,215.9,215.9, &
         215.8,215.8,215.7,215.6,215.6,215.6,215.7,215.9,216.2,216.7, &
         217.3,218.2,220.3,223.3,226.8,231.2,236.4,241.8,247.1,252.0, &
         255.9,257.8,255.4,248.0,243.1,235.5,222.5,212.4,210.1,206.3, &
         213.5,222.1/)
      TO3(:,6) = (/229.6,229.4,229.2,228.9,228.6,228.3,228.1,227.8, &
         227.6,227.2,226.8,226.6,226.3,226.2,226.1,226.0,226.0,226.0, &
         226.1,226.2,226.4,226.6,227.1,227.7,228.6,229.6,230.6,231.6, &
         232.8,234.1,236.8,240.9,245.4,250.0,254.6,259.1,263.1,265.9, &
         266.8,264.9,257.3,244.9,226.9,205.7,189.0,181.6,178.4,180.5, &
         187.4,195.2/)
      TO3(:,7) = (/229.6,229.3,229.0,228.7,228.4,228.2,227.9,227.7, &
         227.5,227.2,226.8,226.5,226.2,225.9,225.8,225.6,225.6,225.6, &
         225.7,226.0,226.4,226.8,227.3,227.9,228.7,229.6,230.6,231.6, &
         232.8,234.0,236.7,240.7,245.6,250.1,254.6,259.1,263.1,265.9, &
         266.7,264.7,257.1,244.9,227.4,206.7,189.0,181.1,178.4,180.3, &
         187.1,195.0/)
      TO3(:,8) = (/227.7,227.4,227.1,226.8,226.6,226.3,226.1,225.8, &
         225.6,225.4,225.1,224.9,224.7,224.6,224.5,224.4,224.5,224.5, &
         224.6,224.7,224.9,225.0,225.3,225.7,226.4,227.3,228.1,229.0, &
         230.1,231.3,233.8,237.6,242.1,246.6,251.1,255.6,259.5,262.3, &
         263.2,261.3,253.9,241.6,223.9,202.9,186.5,179.1,176.0,178.1, &
         184.9,192.8/)
      TO3(:,9) = (/235.4,235.3,235.2,235.1,234.9,234.8,234.7,234.6, &
         234.5,234.4,234.2,234.1,234.0,234.0,234.1,234.3,234.4,234.6, &
         234.8,235.0,235.2,235.5,235.8,236.2,236.7,237.3,238.1,239.0, &
         240.2,241.4,243.9,247.4,251.4,257.0,262.5,267.8,271.9,273.6, &
         273.6,272.0,265.1,249.7,227.1,202.9,179.7,165.6,163.7,167.0, &
         180.7,188.5/)
      TO3(:,10) = (/211.6,211.5,211.4,211.3,211.2,211.0,210.9,210.8, &
         210.7,210.6,210.5,210.3,210.2,210.1,210.0,209.9,209.8,209.7, &
         209.6,209.5,209.3,209.3,209.3,209.4,209.7,210.0,210.3,210.8, &
         211.3,211.8,213.1,215.5,218.8,223.3,228.6,234.1,239.4,244.3, &
         248.0,249.8,247.6,240.0,235.1,227.5,214.5,204.4,202.1,198.3, &
         205.5,213.3/)
      mem_pico = (/1.7439,1.8264,1.9128,1.9992,2.0895,2.1799,2.2702,2.3684, &
         2.4666,2.5687,2.6669,2.7612,2.8437,2.9183,2.9890,3.0479,3.1029, &
         3.1500,3.1854,3.2089,3.2247,3.2325,3.2286,3.2168,3.1932,3.1540, &
         3.1029,3.0361,2.9576,2.8712,2.7848,2.6944,2.5137,2.4273,2.3488, &
         2.2781,2.2486,2.2192,2.1720,2.1328,2.1013,2.0660,2.0267,1.9835, &
         1.9285,1.8657,1.7989,1.7203,1.6339,1.5357,1.4336,1.3276,1.2176, &
         1.1076,1.0016,0.8994,0.8013,0.7109,0.6284,0.5538,0.4870,0.4320, &
         0.3782,0.3307,0.2875,0.2486,0.2137,0.1842,0.1599,0.1402,0.1233, &
         0.1080,0.0935,0.0789,0.0656,0.0530,0.0424,0.0344,0.0290,0.0260, &
         0.0258,0.0268,0.0304,0.0320,0.0331,0.0347,0.0355,0.0363,0.0382, &
         0.0401,0.0416,0.0428,0.0432,0.0432,0.0432,0.0424,0.0416,0.0408, &
         0.0408,0.0424,0.0452,0.0503,0.0562,0.0628,0.0695,0.0758,0.0821, &
         0.0880,0.0939,0.1002,0.1060,0.1123,0.1178,0.1229,0.1261,0.1280, &
         0.1288,0.1296,0.1308,0.1331,0.1371,0.1422,0.1493,0.1591,0.1728, &
         0.1909,0.2137,0.2416,0.2757,0.3178,0.3692,0.4281,0.5499,0.6009, &
         0.6324,0.6402,0.6324,0.6245,0.5892,0.5342,0.4674,0.3967,0.3276, &
         0.2635,0.2078,0.1618,0.1249,0.0958,0.0746,0.0601,0.0503/)
      mem_micro = (/1.574,1.584,1.600,1.617,1.633,1.654,1.669,1.674,1.684, &
         1.697,1.708,1.710,1.716,1.737,1.763,1.793,1.812,1.827,1.830, &
         1.834,1.824,1.800,1.771,1.741,1.712,1.685,1.667,1.650,1.641, &
         1.631,1.631,1.623,1.616,1.606,1.592,1.568,1.542,1.509,1.481, &
         1.459,1.437,1.415,1.399,1.387,1.377,1.367,1.349,1.338,1.319, &
         1.301,1.271,1.242,1.222,1.196,1.169,1.141,1.118,1.096,1.075, &
         1.057,1.035,1.013,0.992,0.977,0.959,0.944,0.927,0.909,0.888, &
         0.868,0.847,0.826,0.806,0.785,0.764,0.737,0.711,0.682,0.653, &
         0.626,0.604,0.580,0.555,0.535,0.514,0.501,0.487,0.478,0.475, &
         0.468,0.464,0.459,0.452,0.452,0.449,0.443,0.433,0.424,0.416, &
         0.406,0.401,0.400,0.403,0.408,0.416,0.429,0.443,0.458,0.473, &
         0.487,0.495,0.499,0.504,0.514,0.521,0.525,0.532,0.535,0.534, &
         0.535,0.532,0.528,0.526,0.528,0.538,0.549,0.574,0.605,0.655, &
         0.720,0.798,0.889,0.979,1.068,1.147,1.207,1.243,1.249,1.227, &
         1.174,1.096,1.004,0.893,0.767,0.635,0.516,0.409,0.323,0.253, &
         0.200,0.158/)
   end subroutine

end module radiation_mod 
