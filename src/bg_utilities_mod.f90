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
      real(r8), intent(out) :: PCO2(NPOC)       ! CO2 fixation rate (s-1)
      real(r8) :: fpar(NPOC), ftemp(NPOC)
      real(r8) :: fsrp(NPOC), fco2(NPOC)
      real(r8) :: Ksrp(NPOC), Vm0(NPOC)
      
      ! temperature factor
      call PhotosynTempFactor(temp-T0, ftemp)
      ! CO2 factor
      fco2 = Dco2 / (Dco2 + Kco2)
      ! nutrient limitation factor
      Ksrp = (/sa_params(Param_Ksrps), sa_params(Param_Ksrpl)/)
      fsrp = Dsrp / (Dsrp + Ksrp)
      ! light factor
      Vm0 = (/sa_params(Param_Vm0s), sa_params(Param_Vm0l)/)
      fpar =  (1.0 - exp(-phAlpha*Ipar/(Vm0+e8))) * exp(-phBeta*Ipar/(Vm0+e8))

      PCO2 = Vm0 * fpar * ftemp * min(fsrp,fco2) / SECOND_OF_DAY
   end subroutine

   subroutine PhotosynTempFactor(temp, ftemp)
      implicit none
      real(r8), intent(in) :: temp           ! celsius
      real(r8), intent(out) :: ftemp(NPOC)
      real(r8) :: thetaH, tfrac
      integer :: kk
      
      ! maximum temperature factor is 3.0
      do kk = 1, NPOC, 1
         if (temp<0._r8) then
            ftemp(kk) = 0._r8
         else if (temp<=Topt_ppk(kk)) then
            ftemp(kk) = fTmax_ppk(kk) * ( ThetaG**(temp-Topt_ppk(kk)) - &
               ThetaG**(-Topt_ppk(kk)) ) + 1._r8
         else
            thetaH = 2._r8 - ThetaG**(-Topt_ppk(kk))
            tfrac = (temp - Topt_ppk(kk)) / (Tmax_ppk(kk) - Topt_ppk(kk))
            ftemp(kk) = fTmax_ppk(kk) * ( thetaH - thetaH**tfrac ) + 1._r8
            ftemp(kk) = max(ftemp(kk), 0._r8)
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   !  Purpose: calculate the ratio of Chlorophill-a to carbon biomass 
   !           1) Chla:C vs depth "Regulation of phytoplankton carbon to 
   !           chlorophyll ratio by light, nutrients and temperature in the
   !           equatorial Pacific Ocean: a basin-scale model" (Wang et al., 2009)
   !           2) Dynamic model ofphytoplankton growth and acclimation: Response
   !           of the balanced growthrate and the chlorophylla: Carbon ratio to
   !           light, nutrient‐limitation andtemperature (Geider et al., 1996)
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
      real(r8) :: fsrp(NPOC), ftemp(NPOC)
      real(r8) :: Vm0(NPOC), Ksrp(NPOC)
      real(r8) :: Pmc(NPOC), fco2(NPOC)
      real(r8) :: fnutri(NPOC)
      real(r8) :: omega, Ipar0
      integer :: ii, kk, nz

      nz = size(Ipar)
      Ipar0 = Ipar(topIndx)
      if (Ipar0>e8) then
         Vm0 = (/sa_params(Param_Vm0s), sa_params(Param_Vm0l)/)
         Ksrp = (/sa_params(Param_Ksrps), sa_params(Param_Ksrpl)/)
         Pmc = Vm0 * (phAlpha/(phAlpha+phBeta)) * (phBeta/(phAlpha+phBeta))** &
            (phBeta/phAlpha)
         do ii = 1, nz, 1
            if (ii>=topIndx) then
               fsrp = Dsrp(ii) / (Dsrp(ii) + Ksrp)
               fco2 = Dco2(ii) / (Dco2(ii) + Kco2) 
               fnutri = min(fsrp,fco2)
               call PhotosynTempFactor(temp(ii)-T0, ftemp)
               do kk = 1, NPOC, 1
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
   ! Purpose: Calculate CH4 oxidation rate. 
   !
   !------------------------------------------------------------------------------
   subroutine Methanotrophy(Dch4, Do2, temp, oCH4)
      implicit none
      real(r8), intent(in) :: Dch4     ! CH4 concentration (mol/m3)
      real(r8), intent(in) :: Do2      ! O2 concentration (mol/m3)
      real(r8), intent(in) :: temp     ! water temperature (K)
      real(r8), intent(out) :: oCH4    ! mol CH4 m-3 s-1
      real(r8) :: Qch4, OQ10, betaCH4, lamdaO2

      Qch4 = sa_params(Param_Qch4)
      OQ10 = sa_params(Param_OQ10)
      betaCH4 = sa_params(Param_betaCH4)
      lamdaO2 = sa_params(Param_lamdaO2) / (1._r8 + Dch4/Dch4_cr)
      oCH4 = Qch4 * (OQ10**(0.1*(temp-Tor))) * (Dch4**betaCH4) * &
         exp(-lamdaO2*Do2) * (1.0 - exp(-betaO2*Do2))
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
      real(r8), intent(in) :: carb(:)     ! unfrozen carbon pool (gC/m3)
      real(r8), intent(in) :: temp        ! sediment temperature (K)
      real(r8), intent(in) :: sal         ! porewater salinity (g/kg)
      real(r8), intent(in) :: eHL         ! redox potential (mV)
      real(r8), intent(out) :: rCH4(:)    ! mol CH4 m-3 s-1
      real(r8), intent(out) :: rCO2(:)    ! mol CO2 m-3 s-1
      real(r8) :: fRx, ftemp, fph, fsal
      real(r8) :: Rc, PQ10

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
      if (carb(pasC)>e8) then
         Rc = sa_params(Param_Rcpas) 
         PQ10 = sa_params(Param_PQ10pas)
         ftemp = PQ10**(0.1*(temp-Tpr(pasC)))
         rCH4(pasC) = 0.25 * Rc * carb(pasC) * ftemp * fRx * fsal / MasC
         rCO2(pasC) = rCH4(pasC)
      else
         rCH4(pasC) = 0.0_r8
         rCO2(pasC) = 0.0_r8
      end if

      ! active C mainly through acetoclastic methanogenesis
      if (carb(actC)>e8) then
         Rc = sa_params(Param_Rcact)
         PQ10 = sa_params(Param_PQ10act)
         ftemp = PQ10**(0.1*(temp-Tpr(actC)))
         rCH4(actC) = 0.5 * Rc * carb(actC) * ftemp * fRx * fsal / MasC
         rCO2(actC) = 4.0 * rCH4(actC)
      else
         rCH4(actC) = 0.0_r8
         rCO2(actC) = 0.0_r8
      end if

      ! yedoma old C mainly through acetoclastic methanogenesis
      if (carb(oldC)>e8) then
         Rc = sa_params(Param_Rcold)
         PQ10 = sa_params(Param_PQ10old)
         ftemp = PQ10**(0.1*(temp-Tpr(oldC)))
         rCH4(oldC) = 0.5 * Rc * carb(oldC) * ftemp * fRx * fsal / MasC
         rCO2(oldC) = 4.0 * rCH4(oldC)
      else
         rCH4(oldC) = 0.0_r8
         rCO2(oldC) = 0.0_r8
      end if
   end subroutine

   subroutine OxicMethanogenesis(gpp, rCH4)
      implicit none
      real(r8), intent(in) :: gpp(:)   ! GPP (gC/m3/s)
      real(r8), intent(out) :: rCH4    ! oxic methane production (mol/m3/s)
      real(r8) :: Rc 

      Rc = sa_params(Param_RcOMP)
      rCH4 = Rc * sum(gpp)
   end subroutine

   subroutine AerobicRespiration(carb, temp, Do2, rCO2)
      implicit none
      real(r8), intent(in) :: carb(:)     ! unfrozen carbon pool (gC/m3)
      real(r8), intent(in) :: temp        ! K
      real(r8), intent(in) :: Do2         ! mol/m3
      real(r8), intent(out) :: rCO2(:)    ! mol/m3/s
      real(r8) :: Rca, ftemp, fo2
     
      if (carb(pasC)>e8) then
         Rca = sa_params(Param_Rcapas)
         ftemp = ThetaCM**(temp-T0-20.0)
         fo2 = Do2 / (Do2 + Ko2CM)
         rCO2(pasC) = Rca * ftemp * fo2 * carb(pasC) / MasC
      else
         rCO2(pasC) = 0.0_r8
      end if

      if (carb(actC)>e8) then
         Rca = sa_params(Param_Rcaact)
         ftemp = ThetaCM**(temp-T0-20.0)
         fo2 = Do2 / (Do2 + Ko2CM)
         rCO2(actC) = Rca * ftemp * fo2 * carb(actC) / MasC
      else
         rCO2(actC) = 0.0_r8
      end if

      rCO2(oldC) = 0.0_r8
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate algae metabolic loss, including autotrophic
   !          respiration (DIC), mortality (detritus) and excretion (DOC).
   !
   !------------------------------------------------------------------------------
   subroutine AutotrophicR(temp, Do2, rloss, rresp, rdoc, rmort)
      implicit none
      real(r8), intent(in) :: temp           ! water temperature (K)
      real(r8), intent(in) :: Do2            ! dissolved O2 (mol/m3)
      real(r8), intent(out) :: rloss(NPOC)   ! metabolic loss rate (s-1)
      real(r8), intent(out) :: rresp(NPOC)   ! repsiration (s-1) 
      real(r8), intent(out) :: rdoc(NPOC)    ! excretion (s-1)
      real(r8), intent(out) :: rmort(NPOC)   ! mortality (s-1)
      real(r8) :: Klr(NPOC)
      real(r8) :: fo2AR, fo2HR

      Klr = (/sa_params(Param_Klrs), sa_params(Param_Klrl)/)
      if (Do2<Ko2AR) then
         fo2AR = 0._r8
         fo2HR = 0._r8
      else
         fo2AR = Do2 / (Do2 + Ko2AR)
         fo2HR = Do2 / (Do2 + Ko2CM)
      end if
      rloss = Klr * (ThetaML**(temp-T0-20.)) / SECOND_OF_DAY
      rresp = Fres * fo2AR * rloss
      rdoc = Fdom * (1._r8 - Fres) * fo2HR * rloss
      rmort = rloss - rresp - rdoc
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Deposition-induced phytoplankton mortality rate.
   !
   !------------------------------------------------------------------------------
   subroutine DepositionMortality(dist, hhyp, rmort)
      implicit none
      real(r8), intent(in) :: dist           ! sinking path to bottom (m)
      real(r8), intent(in) :: hhyp           ! hypolimnion thickness (m)
      real(r8), intent(out) :: rmort(NPOC)   ! mortality rate (s-1)
      real(r8) :: fDepMort(NPOC)

      fDepMort = (/sa_params(Param_fDepMorts), sa_params(Param_fDepMortl)/)
      if (dist>hhyp .or. hhyp<e8) then
         rmort = 0._r8
      else
         rmort = Vs0 / hhyp * fDepMort 
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate microbial mineralization rate.
   !
   !------------------------------------------------------------------------------
   subroutine HeterotrophicR(temp, Do2, rCO2)
      implicit none
      real(r8), intent(in) :: temp        ! water temperature (K)
      real(r8), intent(in) :: Do2         ! dissolved O2 (mol/m3)
      real(r8), intent(out) :: rCO2       ! mol m-3 s-1
      real(r8) :: rRDOM(2), ftemp, fo2

      ftemp = ThetaCM**(temp-T0-20.0)
      fo2 = Do2 / (Do2 + Ko2CM)
      rCO2 = rCO2_BOD * ftemp * fo2
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

      cMuMaxVeg = sa_params(Param_MuMaxVeg)
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
