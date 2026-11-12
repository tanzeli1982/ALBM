module bg_utilities_mod
!---------------------------------------------------------------------------------
!
! Purpose: this module contains utilities related to aquatic biogeosciences. 
!
!---------------------------------------------------------------------------------
   use shr_ctrl_mod,    only: e8 => SHR_CTRL_E8
   use bg_const_mod
   use phy_const_mod
   use shr_param_mod
   use math_utilities_mod, only: SolveQuadraticEquation

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the absorption coefficient of CDOM (Cory et al., 2014
   !           Science; Cory et al., 2013 PNAS; and Koehler et al., 2014 GBC)
   !           
   !  The spectral slopes of yedoma and tundra lakes are fitted according to
   !  Table S2 of Cory et al. (2014). The spectral slope of boreal lake is set
   !  to a value of Fujii et al. (2007).
   !
   !------------------------------------------------------------------------------
   subroutine CalcAcCDOM(ltype, DOC, vwvln, aCDOM)
      implicit none
      integer, intent(in) :: ltype        ! lake type identifier 
      real(r8), intent(in) :: DOC         ! DOC concentration (mol C/m3)
      real(r8), intent(in) :: vwvln(:)    ! wavelength (nm)
      real(r8), intent(out) :: aCDOM(:)   ! absorption coefficient (m-1)
      real(r8) :: suva, sac

      suva = SUVA305(ltype)
      sac = SCDOM(ltype)
      aCDOM = suva * DOC * MasC * exp(-sac*(vwvln-305.0))
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the absorption coefficient of Non-algal particles, e.g. 
   !           detritus (Fujii et al., 2007 Biogeosciences).
   !
   !------------------------------------------------------------------------------
   subroutine CalcAcNAP(DPOC, vwvln, aNAP)
      implicit none
      real(r8), intent(in) :: DPOC        ! dead algae biomass (mol C/m3)
      real(r8), intent(in) :: vwvln(:)    ! wavelength (nm)
      real(r8), intent(out) :: aNAP(:)    ! absorption coefficient (m-1)

      aNAP = aNAP_440 * DPOC * MasC * exp(-0.011*(vwvln-440.0))
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the total backscattering coefficient of water and POC
   !           (Morel, 1974; Fujii et al., 2007 Biogeosciences; Twardowski et
   !           al., 2007 Biogeosciences).
   !
   !------------------------------------------------------------------------------
   subroutine CalcBbWaterPOC(mPOC, vwvln, Bp)
      implicit none
      real(r8), parameter :: Bbbg = 1.7d-4   ! background backscattering coefficient (m-1)
      real(r8), parameter :: Bbsw = 0.5_r8   ! The backscattering ratio of sea water
      real(r8), intent(in) :: mPOC(:)     ! phytoplankton biomass (gC/m3)
      real(r8), intent(in) :: vwvln(:)    ! wavelength (nm)
      real(r8), intent(out) :: Bp(:)      ! backscattering coefficient (m-1)
      real(r8) :: aa, bb, tPPOC, tMPOC

      tPPOC = 1.0d3 * mPOC(small_ppk) ! mg C/m3
      tMPOC = 1.0d3 * mPOC(large_ppk) ! mg C/m3 
      aa = (tPPOC/476935.8)**(1.0/1.277)
      bb = (tMPOC/17069.0)**(1.0/0.859)
      Bp = Bbbg + aa * (vwvln/510.0)**(-0.5) + bb + Bbsw * 3.5d-3 * &
            (vwvln/450.0)**(-4.32)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: DOC photo degradation
   !
   !------------------------------------------------------------------------------
   subroutine PhotoDegradation(ltype, DOC, vwvln, irradiance, rate)
      implicit none
      integer, intent(in) :: ltype              ! lake type identifier
      real(r8), intent(in) :: DOC               ! mol/m3
      real(r8), intent(in) :: vwvln(:)          ! nm
      real(r8), intent(in) :: irradiance(:)     ! mol/m2/s/nm
      real(r8), intent(out) :: rate             ! mol C m-3 s-1
      real(r8) :: AQYCDOM0, AQYCDOM1, avgAQY
      real(r8) :: aCDOM0, aCDOM1, wvln, suva, sac
      integer :: nwvln, ii

      suva = SUVA305(ltype)
      sac = SCDOM(ltype)
      nwvln = size(vwvln)
      rate = 0.0_r8
      do ii = 1, nwvln-1, 1
         ! apparent quantum yield (mol C mol photons-1)
         if (ii==1) then
            wvln = vwvln(ii)
            AQYCDOM0 = AQYem1(ltype) * exp(-AQYm2(ltype)*(wvln-290.0))
            aCDOM0 = suva * DOC * MasC * exp(-sac*(wvln-305.0)) 
         end if
         if (vwvln(ii+1)>700) then
            exit
         end if
         wvln = vwvln(ii+1)
         AQYCDOM1 = AQYem1(ltype) * exp(-AQYm2(ltype)*(wvln-290.0))
         aCDOM1 = suva * DOC * MasC * exp(-sac*(wvln-305.0))
         avgAQY = 0.5 * (AQYCDOM0*aCDOM0*irradiance(ii) + AQYCDOM1*aCDOM1* &
                  irradiance(ii+1))
         rate = rate + 1.0d+6 * avgAQY * (vwvln(ii+1)-vwvln(ii))
         AQYCDOM0 = AQYCDOM1
         aCDOM0 = aCDOM1
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: There are many marine ecosystem models including two algae 
   !          function types (diatom and non-diatom) [Li et al., 2010].
   !          non-diatom algae dominates as Arctic LTER data shows.
   !
   !------------------------------------------------------------------------------
   subroutine Photosynthesis(temp, Dco2, Dsrp, Ipar, PCO2)
      implicit none
      real(r8), intent(in) :: temp              ! water temperature (K)
      real(r8), intent(in) :: Dco2              ! CO2 (mol/m3)
      real(r8), intent(in) :: Dsrp              ! soluble reactive P (mol m-3)
      real(r8), intent(in) :: Ipar              ! mol m-2 s-1
      real(r8), intent(out) :: PCO2(NPTK)       ! CO2 fixation rate (s-1)
      real(r8) :: fpar(NPTK), ftemp(NPTK)
      real(r8) :: fsrp(NPTK), fco2(NPTK)
      
      ! temperature factor
      call PhytoGrowthTempFactor(temp, ftemp)
      ! CO2 factor
      fco2 = Dco2 / (Dco2 + Kco2)
      ! nutrient limitation factor
      fsrp = Dsrp / (Dsrp + PAR_Ksrp)
      ! light factor
      fpar =  (1.0 - exp(-phAlpha*Ipar/VmaxConst)) * exp(-phBeta*Ipar/VmaxConst)

      PCO2 = PAR_Vmax * fpar * ftemp * min(fsrp,fco2) / SECOND_OF_DAY
   end subroutine

   subroutine PhytoGrowthTempFactor(temp, ftemp)
      implicit none
      real(r8), intent(in) :: temp              ! water temperature (K)
      real(r8), intent(out) :: ftemp(NPTK)      ! unitless
      real(r8), parameter :: Tref = 20._r8      ! reference temperature
      real(r8) :: Tm

      Tm = temp - T0    ! celsius
      ftemp = exp( -0.5_r8/PAR_SigT**2 * ( (Tm-PAR_Topt)**2 - &
               (TRef-PAR_Topt)**2) )
   end subroutine

   !------------------------------------------------------------------------------
   !
   !  Purpose: calculate the ratio of Chlorophill-a to carbon biomass 
   !           1) Chla:C vs depth "Regulation of phytoplankton carbon to 
   !           chlorophyll ratio by light, nutrients and temperature in the
   !           equatorial Pacific Ocean: a basin-scale model" (Wang et al., 2009)
   !           2) Dynamic model of phytoplankton growth and acclimation: Response
   !           of the balanced growth rate and the chlorophyll a: Carbon ratio to
   !           light, nutrient‐limitation and temperature (Geider et al., 1996)
   !
   !------------------------------------------------------------------------------
   subroutine CalcChl2CRatio(Ipar, temp, Dsrp, Dco2, topIndx, Chl2C)
      implicit none
      real(r8), intent(in) :: Ipar(:)        ! PAR radiation (mol/m2/s)
      real(r8), intent(in) :: temp(:)        ! water temperature (K)
      real(r8), intent(in) :: Dsrp(:)        ! SRP (mol/m3)
      real(r8), intent(in) :: Dco2(:)        ! CO2 (mol/m3)
      integer, intent(in) :: topIndx         ! top water layer index
      real(r8), intent(out) :: Chl2C(:,:)    ! g Chl gC-1
      real(r8) :: fsrp(NPTK), ftemp(NPTK)
      real(r8) :: Pmc(NPTK), fco2(NPTK)
      real(r8) :: fnutri(NPTK)
      real(r8) :: omega, Ipar0, Tm
      integer :: ii, kk, nz

      nz = size(Ipar)
      Ipar0 = Ipar(topIndx)
      if (Ipar0>e8) then
         Pmc = PAR_Vmax * (phAlpha/(phAlpha+phBeta)) * &
            (phBeta/(phAlpha+phBeta))**(phBeta/phAlpha)
         do ii = 1, nz, 1
            if (ii>=topIndx) then
               fsrp = Dsrp(ii) / (Dsrp(ii) + PAR_Ksrp)
               fco2 = Dco2(ii) / (Dco2(ii) + Kco2) 
               fnutri = min(fsrp,fco2)
               call PhytoGrowthTempFactor(temp(ii), ftemp)
               do kk = 1, NPTK, 1
                  if (ftemp(kk)<e8 .or. fnutri(kk)<e8) then
                     Chl2C(kk,ii) = Chl2Cmin(kk) 
                  else
                     omega = Chl2Cmax(kk) * alphaChl(kk) * (1.d6*Ipar(ii)) / &
                        (Pmc(kk)*ftemp(kk)*fnutri(kk)) 
                     Chl2C(kk,ii) = max( 1./(1.+0.5*omega)*Chl2Cmax(kk), &
                        Chl2Cmin(kk) )
                  end if
               end do
            else
               Chl2C(:,ii) = Chl2Cmin
            end if
         end do
      else
         do ii = 1, nz, 1
            Chl2C(:,ii) = Chl2Cmax
         end do
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate algae metabolic loss, including autotrophic
   !          respiration (DIC), mortality (detritus) and excretion (DOC).
   !
   !------------------------------------------------------------------------------
   subroutine AutotrophicR(temp, rresp, rmort)
      implicit none
      real(r8), intent(in) :: temp           ! water temperature (K)
      real(r8), intent(out) :: rresp(NPTK)   ! repsiration (s-1) 
      real(r8), intent(out) :: rmort(NPTK)   ! mortality (s-1)
      real(r8) :: ftemp(NPTK)

      call PhytoGrowthTempFactor(temp, ftemp)
      rresp = PAR_Klr * ftemp / SECOND_OF_DAY
      rmort = KMort / SECOND_OF_DAY
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Deposition-induced phytoplankton mortality rate.
   !
   !------------------------------------------------------------------------------
   subroutine DepositionMortality(Vs, lensink, fsed, rmort)
      implicit none
      real(r8), intent(in) :: Vs(NPTK)       ! sink velocity (m/s)
      real(r8), intent(in) :: lensink        ! sinking length to sediment bed (m)
      real(r8), intent(in) :: fsed           ! ratio of water-sediment interface to cross section
      real(r8), intent(out) :: rmort(NPTK)   ! mortality rate (s-1)
      integer :: ii

      do ii = 1, NPTK, 1
         rmort(ii) = fsed * max(Vs(ii),0._r8) / lensink * PAR_fDepMort(ii)
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Zooplankton grazing induced phytoplankton mortality rate based on
   !          Lotka–Volterra predator–prey model.
   !
   !------------------------------------------------------------------------------
   subroutine GrazingMortality(temp, Do2, zooC, phytoC, rmort)
      implicit none
      real(r8), intent(in) :: temp           ! water temperature (K)
      real(r8), intent(in) :: Do2            ! dissolved O2 (mol/m3)
      real(r8), intent(in) :: zooC           ! zooplankton C (gC/m3)
      real(r8), intent(in) :: phytoC(NPTK)   ! phytoplankton C (gC/m3)
      real(r8), intent(out) :: rmort(NPTK)   ! mortality (s-1)
      real(r8) :: ftemp, fo2, ffood

      call ZooGrowthTempFactor(temp, ftemp)
      ! grazing reduced at hypoxia conditions
      fo2 = (Do2/(Do2+Ko2AR))
      ! food saturation factor
      ffood = 1._r8 !hFoodZoo / (hFoodZoo + sum(phytoC))
      rmort = PrefZoo * PAR_Gmax * ftemp * fo2 * ffood * max(zooC,0._r8) &
         / SECOND_OF_DAY
   end subroutine

   subroutine ZooGrowthTempFactor(temp, ftemp)
      implicit none
      real(r8), intent(in) :: temp              ! water temperature (K)
      real(r8), intent(out) :: ftemp            ! unitless
      real(r8), parameter :: Tref = 20._r8      ! reference temperature
      real(r8) :: Tm

      Tm = temp - T0    ! celsius
      ftemp = exp( -0.5_r8/PAR_zSigT**2 * ( (Tm-PAR_zTopt)**2 - &
               (TRef-PAR_zTopt)**2) )
   end subroutine

   subroutine ZooplanktonLoss(temp, Do2, rresp, rmort)
      implicit none
      real(r8), intent(in) :: temp           ! water temperature (K)
      real(r8), intent(in) :: Do2            ! dissolved O2 (mol/m3)
      real(r8), intent(out) :: rresp         ! respiration rate (s-1)
      real(r8), intent(out) :: rmort         ! mortality rate (s-1)
      real(r8) :: ftemp, fo2, fo2r

      call ZooGrowthTempFactor(temp, ftemp)
      fo2 = (Do2/(Do2+Ko2AR))
      rresp = PAR_kDRespZoo * ftemp * fo2 / SECOND_OF_DAY
      if (Do2 < Ko2AR) then
         fo2r = 1._r8 + Ko2AR / (Do2+0.1_r8*Ko2AR)
         rmort = PAR_kMortZoo * fo2r / SECOND_OF_DAY 
      else
         rmort = PAR_kMortZoo / SECOND_OF_DAY
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate CH4 oxidation rate. 
   !
   !------------------------------------------------------------------------------
   subroutine Methanotrophy(Dch4, Do2, temp, oCH4)
      implicit none
      real(r8), intent(in) :: Dch4     ! CH4 concentration (mol/m3)
      real(r8), intent(in) :: Do2      ! O2 concentration (mol/m3)
      real(r8), intent(in) :: temp     ! water temperature (K)
      real(r8), intent(out) :: oCH4    ! mol CH4 m-3 s-1
      real(r8) :: lamdaO2

      lamdaO2 = PAR_lamdaO2 / (1._r8 + Dch4/Dch4_cr)
      oCH4 = PAR_Qch4 * (PAR_OQ10**(0.1*(temp-Tor))) * &
         (Dch4**PAR_betaCH4) * exp(-lamdaO2*Do2) * (1.0 - exp(-betaO2*Do2))
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate microbial mineralization rate in water column.
   !
   !------------------------------------------------------------------------------
   subroutine HeterotrophicR(temp, Do2, trophic, rCO2)
      implicit none
      real(r8), intent(in) :: temp        ! water temperature (K)
      real(r8), intent(in) :: Do2         ! dissolved O2 (mol/m3)
      integer, intent(in) :: trophic      ! trophic class
      real(r8), intent(out) :: rCO2       ! mol m-3 s-1
      real(r8) :: rRDOM(2), ftemp, fo2

      ftemp = ThetaCM**(temp-T0-20.0)
      fo2 = Do2 / (Do2 + Ko2CM)
      rCO2 = rCO2_BOD(trophic) * ftemp * fo2
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate carbon mineralization rate in sediments due to 
   !          aerobic and anaerobic reactions. 
   !          CH4:CO2 ratio is based on Zhu et al., 2023 L&O Letters
   !
   !------------------------------------------------------------------------------
   subroutine Methanogenesis(carb, temp, sal, eHL, rCH4, rCO2)
      implicit none
      real(r8), intent(in) :: carb(NSEDC)    ! unfrozen carbon pool (gC/m3)
      real(r8), intent(in) :: temp           ! sediment temperature (K)
      real(r8), intent(in) :: sal            ! porewater salinity (g/kg)
      real(r8), intent(in) :: eHL            ! redox potential (mV)
      real(r8), intent(out) :: rCH4(NSEDC)   ! mol CH4 m-3 s-1
      real(r8), intent(out) :: rCO2(NSEDC)   ! mol CO2 m-3 s-1
      real(r8), parameter :: frcCH4(NSEDC) = (/0.25, 0.5, 0.5/)
      real(r8), parameter :: frcCO2(NSEDC) = (/1.0, 4.0, 4.0/)
      real(r8) :: fRx, ftemp, fph, fsal
      integer :: ii

      ! pH factor (borrowed from ELM wetland CH4 model)
      !if (pH>pHmin .and. pH<pHmax) then
      !   fph = 10._r8**(-0.2235_r8*pH**2._r8 + 2.7727_r8*pH - 8.6_r8)
      !   fph = min(1._r8, max(0._r8, fph))
      !else
      !   fph = 0._r8
      !end if
      fph = 1.0_r8

      ! salinity factor (Poffenbarger et al., 2011)
      fsal = exp(-0.056*sal) 

      ! redox potential effects (Zhuang et al., 2004, GBC) 
      if (eHL<=-200.0_r8) then
         fRx = 1.0_r8
      else if (eHL>=-100.0_r8) then
         fRx = 0.0_r8
      else
         fRx = -0.01 * eHL - 1.0 
      end if

      ! passive C mainly through hydrogenotrophic methanogenesis
      ! active and old C mainly through acetoclastic methanogenesis
      do ii = 1, NSEDC, 1
         if (carb(ii)>e8) then
            ftemp = PAR_PQ10(ii)**(0.1*(temp-Tpr(ii)))
            rCH4(ii) = frcCH4(ii) * PAR_Rc_ch4(ii) * carb(ii) * &
               ftemp * fRx * fsal / MasC
            rCO2(ii) = frcCO2(ii) * rCH4(ii)
         else
            rCH4(ii) = 0._r8
            rCO2(ii) = 0._r8
         end if
      end do

   end subroutine

   subroutine OxicMethanogenesis(gpp, rCH4)
      implicit none
      real(r8), intent(in) :: gpp(:)   ! GPP (gC/m3/s)
      real(r8), intent(out) :: rCH4    ! oxic methane production (mol/m3/s)

      rCH4 = PAR_RcOMP * sum(gpp)
   end subroutine

   subroutine AerobicRespiration(carb, temp, Do2, rCO2)
      implicit none
      real(r8), intent(in) :: carb(NSEDC)    ! unfrozen carbon pool (gC/m3)
      real(r8), intent(in) :: temp           ! K
      real(r8), intent(in) :: Do2            ! mol/m3
      real(r8), intent(out) :: rCO2(NSEDC)   ! mol/m3/s
      real(r8) :: Rca, ftemp, fo2
      integer :: ii

      ftemp = ThetaCM**(temp-T0-20.0)
      fo2 = Do2 / (Do2 + Ko2CM)
      do ii = 1, NSEDC, 1
         if (carb(ii)>e8) then
            rCO2(ii) = PAR_Rc_ar(ii) * ftemp * fo2 * carb(ii) / MasC 
         else
            rCO2(ii) = 0._r8 
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Partition sediment mineral P into SRP based on dissolved O2
   !           (PCLake; Osaka et al., 2021, Limnol Oceanogr Methods).
   !
   !------------------------------------------------------------------------------
   function CalcSedSRPConc(minP, minFe, Do2)
      implicit none
      real(r8), intent(in) :: minP  ! mineral P (mol/m3)
      real(r8), intent(in) :: minFe ! mineral Fe (g/m3) 
      real(r8), intent(in) :: Do2   ! mol/m3
      real(r8) :: CalcSedSRPConc    ! mol/m3
      real(r8) :: afOxySed, aKPAdsS
      real(r8) :: aPAdsMaxS
      real(r8) :: A(3), roots(2)

      afOxySed = min(Do2*MasO2/13.0_r8, 1.0_r8)  ! aerobic fraction
      aKPAdsS = cKPAdsOxS * (1.0 - fRedMaxS*(1.0-afOxySed)) 
      aPAdsMaxS = cRelPAdsDWS + afOxySed*cRelPAdsFeS*minFe
      A(1) = aKPAdsS
      A(2) = 1.0 + aPAdsMaxS*aKPAdsS - aKPAdsS*(minP*MasP)
      A(3) = -minP*MasP
      call SolveQuadraticEquation(A, roots)
      CalcSedSRPConc = maxval(roots)/MasP
      return
   end function

   function CalcSedMinPConc(SRP, minFe)
      implicit none
      real(r8), intent(in) :: SRP   ! soluable reactive P (gP/m3)
      real(r8), intent(in) :: minFe ! mineral Fe (g/m3)
      real(r8) :: CalcSedMinPConc   ! mol/m3
      real(r8) :: afOxySed, aKPAdsS
      real(r8) :: aPAdsMaxS

      afOxySed = 1.0_r8  ! aerobic fraction
      aKPAdsS = cKPAdsOxS * (1.0 - fRedMaxS*(1.0-afOxySed))
      aPAdsMaxS = cRelPAdsDWS + afOxySed*cRelPAdsFeS*minFe
      CalcSedMinPConc = SRP * (aKPAdsS*SRP + (1._r8+aPAdsMaxS*aKPAdsS)) / &
         (1._r8 + aKPAdsS*SRP) / MasP
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate submerged macrophyte production, including 
   !          photorespiration, maintenance respiration, and mortality 
   !          following PCLake. 
   !
   !------------------------------------------------------------------------------   
   subroutine MacrophyteMetabolism(biomass, Ipar, temp, tonfset, phenol, &
                                   shootfrac, rProd, rResp, rMort)
      implicit none
      real(r8), intent(in) :: biomass     ! gDW/m2
      real(r8), intent(in) :: Ipar        ! mol/m2/s
      real(r8), intent(in) :: temp        ! K
      real(r8), intent(in) :: tonfset     ! day 
      integer, intent(in) :: phenol       ! phenology flag
      real(r8), intent(out) :: shootfrac  ! fraction
      real(r8), intent(out) :: rProd      ! gDW/m2/s
      real(r8), intent(out) :: rResp      ! gDW/m2/s
      real(r8), intent(out) :: rMort      ! gDW/m2/s
      real(r8) :: ftemp_p, ftemp_r, fpar, pPAR 
      real(r8) :: rootfrac, donfset, cMuMaxVeg
      real(r8) :: akDIncrVeg, aMuVeg, tDEnvVeg
      real(r8) :: tDEnvProdVeg, bkMortVeg

      pPAR = Ipar * 1.d6 / fconvPAR    ! W/m2
      ! temperature factor
      ftemp_p = cQ10ProdVeg ** (0.1_r8 * (temp - T0 - 20._r8))
      ftemp_r = cQ10RespVeg ** (0.1_r8 * (temp - T0 - 20._r8))
      ! light factor
      fpar = pPAR / (pPAR + hLRefVeg) 
      ! shoot fraction
      donfset = cLengTrans - tonfset
      if (phenol==dormant_season) then
         rootfrac = fRootVegWin
      else if (phenol==grow_season) then
         rootfrac = fRootVegSum
      else if (phenol==onset_season) then
         rootfrac = 0.5 * (fRootVegWin + fRootVegSum) + 0.5 * (fRootVegWin - &
            fRootVegSum) * cos(Pi*donfset/cLengTrans)
      else if (phenol==offset_season) then
         rootfrac = 0.5 * (fRootVegWin + fRootVegSum) - 0.5 * (fRootVegWin - &
            fRootVegSum) * cos(Pi*donfset/cLengTrans)
      end if
      shootfrac = 1.0_r8 - rootfrac

      if (phenol==offset_season) then
         bkMortVeg = -log(fWinVeg) / cLengTrans
      else
         bkMortVeg = kMortVegSum
      end if

      cMuMaxVeg = PAR_MuMaxVeg
      if (phenol==dormant_season) then
         aMuVeg = 0._r8
         tDEnvVeg = 0._r8
         tDEnvProdVeg = 0._r8
      else
         ! max. growth rate at current temp. AND light (d-1)
         aMuVeg = ftemp_p * fpar * shootfrac * cMuMaxVeg
         ! intrinsc growth rate (d-1)
         akDIncrVeg = aMuVeg - kDRespVeg * ftemp_r - bkMortVeg
         ! growth rate correction due to capacity control (g/m2/d)
         tDEnvVeg = max( akDIncrVeg / cDCarrVeg * biomass**2._r8, 0._r8 )
         ! production reduction proportion from growth rate correction (g/m2/d)
         tDEnvProdVeg = aMuVeg / cMuMaxVeg * tDEnvVeg
      end if

      rProd = max(0._r8, aMuVeg * biomass - tDEnvProdVeg) / SECOND_OF_DAY
      rResp = kDRespVeg * ftemp_r * biomass / SECOND_OF_DAY
      rMort = (bkMortVeg * biomass + tDEnvVeg - tDEnvProdVeg) / SECOND_OF_DAY
   end subroutine

end module bg_utilities_mod
