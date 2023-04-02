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

contains
   !------------------------------------------------------------------------------
   !  
   ! Purpose: the integral of PAR and total irradiance
   !
   !------------------------------------------------------------------------------
   subroutine GetIncidentPAR(vwvln, fphot, rPAR, pPAR)
      implicit none
      real(r8), intent(in) :: vwvln(:)          ! nm
      real(r8), intent(in) :: fphot(:)          ! mol m-2 s-1 nm-1
      real(r8), intent(out) :: rPAR             ! mol m-2 s-1
      real(r8), optional, intent(out) :: pPAR   ! W/m2
      integer :: ii, nwvln

      nwvln = size(fphot)
      rPAR = 0.0d+0
      do ii = 1, nwvln-1, 1
         if (vwvln(ii)>=400 .and. vwvln(ii+1)<=700) then
            rPAR = rPAR + 0.5 * (vwvln(ii+1)-vwvln(ii)) * &
                  (fphot(ii)+fphot(ii+1))
         end if
      end do
      if (present(pPAR)) then
         pPAR = rPAR * 1d6 / fconvPAR 
      end if
   end subroutine

   subroutine GetIncidentSRD(vwvln, fphot, rTOT)
      implicit none
      real(r8), intent(in) :: vwvln(:)       ! nm
      real(r8), intent(in) :: fphot(:)       ! mol m-2 s-1 nm-1
      real(r8), intent(out) :: rTOT          ! W/m2
      real(r8), parameter :: PHOT = 0.119626565801574d+0  ! units: J m mol-1
      real(r8) :: swv, dwv
      integer :: ii, nwvln

      nwvln = size(fphot)
      rTOT = 0.0d+0
      do ii = 1, nwvln-1, 1
         dwv = vwvln(ii+1) - vwvln(ii)
         swv = 0.5d-9 * (vwvln(ii) + vwvln(ii+1))  ! m
         rTOT = rTOT + 0.5 * dwv * (fphot(ii)+fphot(ii+1)) * PHOT / swv
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the absorption coefficient of pure water (Absorption
   !           spectrum (380-700 nm) of pure water. II. Integrating cavity
   !           measurements [Pope & Fry, 1997]).
   !
   !------------------------------------------------------------------------------
   subroutine CalcAcWater(vwvln, aw)
      implicit none
      real(r8), intent(in) :: vwvln(:)       ! wavelength (nm)
      real(r8), intent(out) :: aw(:)         ! absorption coefficient (m-1)
      real(r8), allocatable :: aw_data(:)
      real(r8) :: par, wvln, wvmin, wvmax, dwv
      integer :: ii, nwvln, indx

      allocate(aw_data(140))
      aw_data = (/0.01137,0.01044,0.00941,0.00917,0.00851,0.00829,0.00813, &
         0.00775,0.00663,0.00579,0.00530,0.00503,0.00473,0.00452,0.00444, &
         0.00442,0.00454,0.00474,0.00478,0.00482,0.00495,0.00504,0.00530, &
         0.00580,0.00635,0.00696,0.00751,0.00830,0.00922,0.00969,0.00962, &
         0.00957,0.00979,0.01005,0.01011,0.0102,0.0106,0.0109,0.0114,0.0121, &
         0.0127,0.0131,0.0136,0.0144,0.0150,0.0162,0.0173,0.0191,0.0204, &
         0.0228,0.0256,0.0280,0.0325,0.0372,0.0396,0.0399,0.0409,0.0416, &
         0.0417,0.0428,0.0434,0.0447,0.0452,0.0466,0.0474,0.0489,0.0511, &
         0.0537,0.0565,0.0593,0.0596,0.0606,0.0619,0.0640,0.0642,0.0672, &
         0.0695,0.0733,0.0772,0.0836,0.0896,0.0989,0.1100,0.1220,0.1351, &
         0.1516,0.1672,0.1925,0.2224,0.2470,0.2577,0.2629,0.2644,0.2665, &
         0.2678,0.2707,0.2755,0.2810,0.2834,0.2904,0.2916,0.2995,0.3012, &
         0.3077,0.3108,0.322,0.325,0.335,0.340,0.358,0.371,0.393,0.410, &
         0.424,0.429,0.436,0.439,0.448,0.448,0.461,0.465,0.478,0.486,0.502, &
         0.516,0.538,0.559,0.592,0.624,0.663,0.704,0.756,0.827,0.914,1.007, &
         1.119,1.231,1.356,1.489,1.678/)

      wvmin = 380.0
      wvmax = 727.5
      dwv = 2.5
      nwvln = size(vwvln)
      do ii = 1, nwvln, 1
         if (vwvln(ii)<wvmin) then
            aw(ii) = aw_data(1) - 1.0d-3*(vwvln(ii)-wvmin)   ! approximately
         else if (vwvln(ii)>wvmax) then
            aw(ii) = min(100.0, 10**(5.61*(log10(vwvln(ii))-2.861833)))
         else
            indx = NINT( (vwvln(ii) - wvmin) / dwv ) + 1
            indx = min(139, indx)
            wvln = wvmin + (indx-1) * dwv
            par = (vwvln(ii) - wvln) / dwv
            aw(ii) = (1.0-par)*aw_data(indx) + par*aw_data(indx+1)
         end if
      end do
      deallocate(aw_data)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the absorption coefficient of pure ice
   !          "Optical constants of ice from the ultraviolet to the microwave: A 
   !          revised compilation" [Warren & Brandt, 2008]
   !          "Visible and near-ultraviolet absorption spectrum of ice
   !          from transmission of solar radiation into snow" [Warren et al.,
   !          2006]
   !
   !          K_ice = 4*Pi*mim/wvln (m_im: imaginary index of refraction)
   !
   !------------------------------------------------------------------------------
   subroutine CalcAcIce(vwvln, aI)
      implicit none
      real(r8), intent(in) :: vwvln(:)    ! wavelength (nm)
      real(r8), intent(out) :: aI(:)      ! absorption coefficient (m-1)
      real(r8), allocatable :: aI_data(:)
      real(r8) :: wvln, wvmin, wvmax, par
      integer :: nwvln, ii, indx

      allocate(aI_data(181))
      aI_data = (/1.78d-3,1.45d-3,1.25d-3,1.08d-3,9.21d-4,8.11d-4,7.31d-4, &
         6.60d-4,6.38d-4,6.79d-4,7.43d-4,7.82d-4,8.18d-4,8.62d-4,9.38d-4, &
         1.04d-3,1.21d-3,1.47d-3,1.79d-3,2.15d-3,2.58d-3,3.05d-3,3.62d-3, &
         4.33d-3,5.23d-3,6.29d-3,7.49d-3,8.91d-3,1.07d-2,1.27d-2,1.48d-2, &
         1.72d-2,1.98d-2,2.28d-2,2.60d-2,2.96d-2,3.34d-2,3.76d-2,4.22d-2, &
         4.71d-2,5.23d-2,5.78d-2,6.37d-2,6.99d-2,7.63d-2,8.31d-2,9.01d-2, &
         9.73d-2,1.05d-1,1.12d-1,1.20d-1,1.42d-1,1.74d-1,2.07d-1,2.40d-1, &
         2.76d-1,3.16d-1,3.54d-1,3.86d-1,4.37d-1,5.21d-1,6.09d-1,7.03d-1, &
         7.40d-1,8.35d-1,9.84d-1,1.17d+0,1.40d+0,1.64d+0,1.88d+0,2.10d+0, &
         2.17d+0,2.19d+0,2.20d+0,2.26d+0,2.71d+0,3.14d+0,3.83d+0,4.78d+0, &
         5.53d+0,5.86d+0,6.13d+0,6.47d+0,6.90d+0,7.39d+0,7.96d+0,9.88d+0, &
         1.20d+1,1.44d+1,1.69d+1,2.04d+1,2.77d+1,2.82d+1,2.32d+1,2.02d+1, &
         1.94d+1,2.04d+1,2.48d+1,3.29d+1,5.08d+1,7.03d+1,1.05d+2,1.24d+2, &
         1.32d+2,1.31d+2,1.28d+2,1.26d+2,1.26d+2,1.31d+2,1.44d+2,1.78d+2, &
         5.27d+2,1.32d+3,2.53d+3,4.11d+3,4.55d+3,4.17d+3,3.70d+3,3.18d+3, &
         2.65d+3,2.27d+3,2.02d+3,1.87d+3,1.71d+3,1.53d+3,1.39d+3,1.26d+3, &
         1.16d+3,1.09d+3,1.04d+3,9.85d+2,9.24d+2,8.95d+2,9.67d+2,1.40d+3, &
         2.36d+3,4.29d+3,6.90d+3,9.13d+3,1.02d+4,1.03d+4,8.44d+3,4.90d+3, &
         2.35d+3,1.46d+3,1.14d+3,1.55d+3,2.45d+3,3.03d+3,3.40d+3,3.78d+3, &
         3.77d+3,3.65d+3,5.26d+3,1.13d+4,3.73d+4,1.08d+5,2.83d+5,7.46d+5, &
         1.33d+6,1.80d+6,2.41d+6,2.36d+6,1.80d+6,1.12d+6,6.56d+5,3.92d+5, &
         2.39d+5,1.46d+5,9.37d+4,6.13d+4,4.25d+4,3.25d+4,2.69d+4,2.42d+4, &
         2.28d+4,2.45d+4,2.71d+4,3.05d+4,3.39d+4,3.76d+4/)

      wvmin = 350.0
      wvmax = 4000.0
      nwvln = size(vwvln)
      do ii = 1, nwvln, 1
         wvln = vwvln(ii)
         if (wvln<=wvmin) then
            aI(ii) = aI_data(1)
         else if (wvln>=wvmax) then
            aI(ii) = aI_data(181)
         else
            if (wvln<=600) then
               indx = int((wvln-350)/5) + 1
               par = (wvln-350-5*(indx-1)) / 5.0
            else if (wvln<=1000) then
               indx = int((wvln-600)/10) + 51
               par = (wvln-600-10*(indx-51)) / 10.0
            else if (wvln<=2000) then
               indx = int((wvln-1000)/20) + 91
               par = (wvln-1000-20*(indx-91)) / 20.0
            else
               indx = int((wvln-2000)/50) + 141
               par = (wvln-2000-50*(indx-141)) / 50.0
            end if
            aI(ii) = par*aI_data(indx+1) + (1.0-par)*aI_data(indx)
         end if
      end do
      deallocate(aI_data)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the absorption coefficient of snow and gray ice
   !          (Rogers et al., 1995; Limnol. Oceanogr.)
   !
   !------------------------------------------------------------------------------
   subroutine CalcAcSnow(vwvln, aN)
      implicit none
      real(r8), intent(in) :: vwvln(:)
      real(r8), intent(out) :: aN(:)
      
      where (vwvln<=700)
         aN = EtanVIS 
      elsewhere
         aN = EtanIR
      end where
   end subroutine

   subroutine CalcAcGrayIce(vwvln, aN)
      implicit none
      real(r8), intent(in) :: vwvln(:)
      real(r8), intent(out) :: aN(:)

      where (vwvln<=700)
         aN = EtaeVIS
      elsewhere
         aN = EtaeIR
      end where
   end subroutine

   subroutine CalcCoverRadAbsfrc(hsnow, hgrayice, fs, fe)
      implicit none
      real(r8), intent(in) :: hsnow       ! units: m
      real(r8), intent(in) :: hgrayice    ! units: m
      real(r8), intent(out) :: fs         ! snow abs fraction
      real(r8), intent(out) :: fe         ! gray ice abs fraction
      real(r8) :: vis1, vis2, ir1, ir2

      ! assume VIS:IR radiation = 1 : 1
      if (hsnow<e8 .and. hgrayice<e8) then
         fs = 0.0_r8
         fe = 0.0_r8
      else
         vis1 = exp(-hsnow*EtanVIS)
         vis2 = exp(-hsnow*EtanVIS-hgrayice*EtaeVIS)
         ir1 = exp(-hsnow*EtanIR)
         ir2 = exp(-hsnow*EtanIR-hgrayice*EtaeIR)
         fe = 0.5*(vis1-vis2)/(1.0-vis2) + 0.5*(ir1-ir2)/(1.0-ir2)
         fs = 1.0 - fe 
      end if
   end subroutine

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
      real(r8), intent(in) :: DOC         ! DOC concentration (umol C/m3)
      real(r8), intent(in) :: vwvln(:)    ! wavelength (nm)
      real(r8), intent(out) :: aCDOM(:)   ! absorption coefficient (m-1)
      real(r8) :: suva, sac

      suva = SUVA305(ltype)
      sac = SCDOM(ltype)
      aCDOM = suva * DOC * MasC * 1.0d-6 * exp(-sac*(vwvln-305.0))
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the absorption coefficient of Non-algal particles, e.g. 
   !           detritus (Fujii et al., 2007 Biogeosciences).
   !
   !------------------------------------------------------------------------------
   subroutine CalcAcNAP(DPOC, vwvln, aNAP)
      implicit none
      real(r8), intent(in) :: DPOC        ! dead algae biomass (umol C/m3)
      real(r8), intent(in) :: vwvln(:)    ! wavelength (nm)
      real(r8), intent(out) :: aNAP(:)    ! absorption coefficient (m-1)

      aNAP = aNAP_440 * DPOC * MasC * 1.0d-6 * exp(-0.011*(vwvln-440.0))
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the absorption coefficient of algae photosynthesis
   !
   !          Chlorophill-a specific absorption coefficient is calculated
   !          according to Fujii et al. (2007, Biogeosciences); Ciotti et al.
   !          (2002, Limnol. Oceanogr.); Ciotti & Bricaud (2006, Limnol.
   !          Oceanogr.: Methods).
   !          
   !          maximum and minimum Chla:C (mg Chl mmol C-1) in Fujii et al. (2007)
   !          Chl2Cmin = 0.036, Chl2Cmax = 1.2
   !
   !------------------------------------------------------------------------------
   subroutine CalcAcAlgae(LPOC, chla, vwvln, rapico_h, ramicro_h, aAP)
      implicit none
      real(r8), intent(in) :: LPOC(:)        ! living algae biomass (umol C/m3)
      real(r8), intent(in) :: chla(:)        ! Chla (mg/m3) 
      real(r8), intent(in) :: vwvln(:)       ! wavelength (nm)
      real(r8), intent(in) :: rapico_h(:)    ! normalized highlight Apico
      real(r8), intent(in) :: ramicro_h(:)   ! normalized highlight Amicro
      real(r8), intent(out) :: aAP(:)        ! absorption coefficient (m-1)
      real(r8), parameter :: wvmin = 400.0
      real(r8), parameter :: wvmax = 700.0
      real(r8), parameter :: dwv = 2.0
      real(r8) :: wvln, par, rChl2C(NPOC)
      real(r8) :: apico_h, apico_l, amicro_h, amicro_l
      real(r8) :: apico, amicro, qChl2C1, qChl2C2
      integer :: nwvln, ii, indx

      if (sum(chla)<e8) then
         aAP = 0.0_r8
         return
      end if

      ! Chl:C ratio (mg Chl mmol C-1)
      rChl2C = 1.0d3 * chla / (LPOC + e8)
      nwvln = size(vwvln)
      qChl2C1 = (rChl2C(small_ppk) - 1.0/C2Chlmax(small_ppk)) / &
         (1.0/C2Chlmin(small_ppk) - 1.0/C2Chlmax(small_ppk))
      qChl2C2 = (rChl2C(large_ppk) - 1.0/C2Chlmax(large_ppk)) / &
         (1.0/C2Chlmin(large_ppk) - 1.0/C2Chlmax(large_ppk))
      qChl2C1 = min(1.0, max(0.0, qChl2C1))
      qChl2C2 = min(1.0, max(0.0, qChl2C2))

      do ii = 1, nwvln, 1
         if (vwvln(ii)>=wvmin .and. vwvln(ii)<=wvmax) then
            indx = NINT( (vwvln(ii) - wvmin)/dwv ) + 1
            indx = min(150, indx)
            wvln = wvmin + (indx-1) * dwv
            par = (vwvln(ii) - wvln) / dwv
            ! high-light and low-light chlorophill-a specific absorption
            ! coefficient (m2 mg chl-1)
            apico_h = (rapico_h(indx)*(1.0-par)+rapico_h(indx+1)*par) &
                     / rapico_h(139) * apico_676
            apico_l = apico_h / 1.5
            amicro_h = (ramicro_h(indx)*(1.0-par)+ramicro_h(indx+1)*par) &
                     / ramicro_h(138) * amicro_674
            amicro_l = amicro_h / 1.5
            apico = apico_h * (1.0 - qChl2C1) + apico_l * qChl2C1
            amicro = amicro_h * (1.0 - qChl2C2) + amicro_l * qChl2C2
            aAP(ii) = apico * rChl2C(small_ppk) * 1.0d-3 * LPOC(small_ppk) + &
                     amicro * rChl2C(large_ppk) * 1.0d-3 * LPOC(large_ppk)
         else
            aAP(ii) = 0.0
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the total backscattering coefficient of water and POC
   !           (Morel, 1974; Fujii et al., 2007 Biogeosciences; Twardowski et
   !           al., 2007 Biogeosciences).
   !
   !------------------------------------------------------------------------------
   subroutine CalcBbWaterPOC(PPOC, MPOC, vwvln, Bp)
      implicit none
      real(r8), intent(in) :: PPOC        ! small phytoplankton (umol C/m3)
      real(r8), intent(in) :: MPOC        ! large phytoplankton (umol C/m3)
      real(r8), intent(in) :: vwvln(:)    ! wavelength (nm)
      real(r8), intent(out) :: Bp(:)      ! backscattering coefficient (m-1)
      real(r8) :: aa, bb, tPPOC, tMPOC

      tPPOC = PPOC * 1.0d-3 * MasC   ! mg C/m3
      tMPOC = MPOC * 1.0d-3 * MasC   ! mg C/m3 
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
      real(r8), intent(in) :: DOC               ! umol/m3
      real(r8), intent(in) :: vwvln(:)          ! nm
      real(r8), intent(in) :: irradiance(:)     ! mol/m2/s/nm
      real(r8), intent(out) :: rate             ! umol C m-3 s-1
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
            aCDOM0 = suva * DOC * MasC * 1.0d-6 * exp(-sac*(wvln-305.0)) 
         end if
         if (vwvln(ii+1)>700) then
            exit
         end if
         wvln = vwvln(ii+1)
         AQYCDOM1 = AQYem1(ltype) * exp(-AQYm2(ltype)*(wvln-290.0))
         aCDOM1 = suva * DOC * MasC * 1.0d-6 * exp(-sac*(wvln-305.0))
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
   subroutine Photosynthesis(Chla, Chl2C, Dco2, Temp, SRP, Ipar, PCO2)
      implicit none
      real(r8), intent(in) :: Chla(NPOC)  ! mg chla m-3
      real(r8), intent(in) :: Chl2C(NPOC) ! mg chla (mg C)-1
      real(r8), intent(in) :: Dco2        ! umol/m3
      real(r8), intent(in) :: Temp        ! water temperature (K)
      real(r8), intent(in) :: SRP         ! soluble reactive P (umol m-3)
      real(r8), intent(in) :: Ipar        ! mol m-2 s-1
      real(r8), intent(out) :: PCO2(NPOC) ! CO2 fixation rate (umol C m-3 s-1)
      real(r8) :: fpar(NPOC), ftemp(NPOC)
      real(r8) :: fsrp(NPOC), fco2(NPOC)
      real(r8) :: phAlpha(NPOC), phBeta(NPOC) 
      real(r8) :: Ksrp(NPOC), Vch(NPOC), Vm0(NPOC)
      real(r8) :: Tw
      
      Tw = Temp - T0
      Vch = (/sa_params(Param_Vchs), sa_params(Param_Vchl)/)
      Vch = Vch * 1.0d+3 / MasC / SECOND_OF_DAY
      Vm0 = (/sa_params(Param_Vchs), sa_params(Param_Vchl)/)
      Vm0 = Vm0 * Chl2C 
      phAlpha = (/sa_params(Param_phAlphas), sa_params(Param_phAlphal)/)
      phBeta = (/sa_params(Param_phBetas), sa_params(Param_phBetal)/)
      fpar =  (1.0 - exp(-phAlpha*Ipar/Vm0)) * exp(-phBeta*Ipar/Vm0)
      ftemp = ThetaG**(Tw-20) - ThetaG**(kt_ppk*(Tw-at_ppk)) + bt_ppk 
      ftemp = max(0.0, ftemp)
      !ftemp = ThetaG**Tw
      fco2 = Dco2 / (Dco2 + Kco2)
      Ksrp = (/sa_params(Param_Ksrps), sa_params(Param_Ksrpl)/)
      fsrp = SRP / (SRP + Ksrp) 
      PCO2 = Vch * fpar * ftemp * fco2 * fsrp * Chla
   end subroutine

   !------------------------------------------------------------------------------
   !
   !  Purpose: calculate the ratio of Chlorophill-a to carbon biomass 
   !           1) Chla:C vs depth "Regulation of phytoplankton carbon to 
   !           chlorophyll ratio by light, nutrients and temperature in the
   !           equatorial Pacific Ocean: a basin-scale model" (Wang et al., 2009)
   !           2) surface Chla:C "A dynamic model of photoadaptation in
   !           phytoplankton" (Geider et al., 1996)
   !
   !           previously, from Wang et al. (2009), surface Chla:C is governed
   !           by growth rate by: C2Chlmax - Kpc2chl * rGrowth   
   !
   !------------------------------------------------------------------------------
   subroutine CalcChl2CRatio(Ipar, wIce, Temp, SRP, Chl2C)
      implicit none
      real(r8), intent(in) :: Ipar(:)        ! PAR radiation (mol/m2/s)
      real(r8), intent(in) :: wIce(:)        ! water ice fraction
      real(r8), intent(in) :: Temp           ! water temperature (K)
      real(r8), intent(in) :: SRP            ! soluble reactive P (umol/m3)
      real(r8), intent(out) :: Chl2C(:,:)    ! mg Chl mmol C-1
      real(r8) :: C2Chl0(NPOC), C2Chl(NPOC)
      real(r8) :: fsrp(NPOC), ftemp(NPOC), Ksrp(NPOC)
      real(r8) :: Ipar0, Tw
      integer :: ii, nn

      nn = size(Ipar)
      Ipar0 = Ipar(1)
      if (Ipar0>e8) then
         Ksrp = (/sa_params(Param_Ksrps), sa_params(Param_Ksrpl)/)
         fsrp = SRP / (SRP + Ksrp)
         Tw = Temp - T0
         ftemp = ThetaG**(Tw-20) - ThetaG**(kt_ppk*(Tw-at_ppk)) + bt_ppk
         ftemp = max(0.0, ftemp)
         C2Chl0 = C2Chlmax - Kpc2chl * mu0 * ftemp * fsrp 
         C2Chl0 = max(C2Chl0, C2Chlmin)
         do ii = 1, nn, 1
            if (wIce(ii)<1.0 .and. Ipar(ii)>0.01*Ipar0) then
               C2Chl = C2Chl0 - (C2Chl0 - C2Chlmin) * log(Ipar0/Ipar(ii)) &
                  / 4.605
               C2Chl = min(max(C2Chl, C2Chlmin), C2Chlmax)
               Chl2C(:,ii) = 1.0 / C2Chl
            else
               Chl2C(:,ii) = 1.0 / C2Chlmax
            end if
         end do
      else
         do ii = 1, nn, 1
            Chl2C(:,ii) = 1.0 / C2Chlmax
         end do 
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate CH4 oxidation rate by Michaelis-Menten kinetics. 
   !
   !------------------------------------------------------------------------------
   subroutine Methanotrophy(Dch4, Do2, temp, rate)
      implicit none
      real(r8), intent(in) :: Dch4     ! CH4 concentration (umol/m3)
      real(r8), intent(in) :: Do2      ! O2 concentration (umol/m3)
      real(r8), intent(in) :: temp     ! water temperature (K)
      real(r8), intent(out) :: rate    ! umol CH4 m-3 s-1
      real(r8) :: Kch4, Ko2, Qch4, OQ10

      Qch4 = sa_params(Param_Qch4)
      OQ10 = sa_params(Param_OQ10)
      Kch4 = sa_params(Param_Kch4)
      Ko2 = sa_params(Param_Ko2)
      rate = Qch4 * (OQ10**(0.1*(temp-Tor))) * (Dch4/(Kch4+Dch4)) * &
            (Do2/(Ko2+Do2))
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate carbon mineralization rate in sediments due to 
   !           aerobic and anaerobic reactions. 
   !
   !------------------------------------------------------------------------------
   subroutine Methanogenesis(carb, temp, Do2, rCH4, rCO2)
      implicit none
      real(r8), intent(in) :: carb(:)     ! unfrozen labile carbon (umol/m3)
      real(r8), intent(in) :: temp        ! sediment temperature (K)
      real(r8), intent(in) :: Do2         ! O2 concentration (umol/m3)
      real(r8), intent(out) :: rCH4(:)    ! umol CH4 m-3 s-1
      real(r8), intent(out) :: rCO2(:)    ! umol CO2 m-3 s-1
      real(r8) :: fo2, ftemp, Rcn, Rco
      real(r8) :: PQ10_act, PQ10_pas

      Rcn = sa_params(Param_Rcn)
      Rco = 2.149145d-10 
      PQ10_pas = sa_params(Param_PQ10n)
      PQ10_act = 1.006950d+00 
      ! O2 suppression (Tang et al., 2010; Biogeosciences)
      fo2 = 1.0 / (1.0 + etaO2*1.0d-6*Do2)
      ! 14C-enriched pool
      ftemp = PQ10_pas**(0.1*(temp-Tpr))
      rCH4(pasC) = 0.25 * Rcn * carb(pasC) * ftemp * fo2
      rCO2(pasC) = 3.0 * rCH4(pasC)
      ! 14C-depleted pool
      ftemp = PQ10_act**(0.1*(temp-Tpr_act))
      rCH4(actC) = 0.5 * Rco * carb(actC) * ftemp * fo2
      rCO2(actC) = rCH4(actC)
   end subroutine

   subroutine AerobicDegradation(pool, temp, Do2, rCO2)
      implicit none
      real(r8), intent(in) :: pool(:)     ! umol/m3
      real(r8), intent(in) :: temp        ! K
      real(r8), intent(in) :: Do2         ! umol/m3
      real(r8), intent(out) :: rCO2(:)    ! umol/m3/s
      real(r8) :: Rca, ftemp, fo2
      
      Rca = sa_params(Param_Rca)
      ftemp = ThetaCM**(temp-T0-20)
      fo2 = Do2 / (Do2 + Ko2CM)
      rCO2 = Rca * ftemp * fo2 * pool
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate pH factor for CH4 production. (A2, Zhuang et al., 2004)
   !
   !------------------------------------------------------------------------------
   subroutine GetCH4ProdpHFactor(pH, factor)
      implicit none
      real(r8), intent(in) :: pH
      real(r8), intent(out) :: factor
      real(r8) :: f1, f2

      f1 = (pH - pHmin) * (pH - pHmax)
      f2 = (pH - pHopt)**2
      factor = f1 / (f1 - f2)
      factor = min(1.0, max(0.0, factor))
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate algae metabolic loss, including autotrophic
   !          respiration (DIC), mortality (detritus) and excretion (DOC).
   !
   !------------------------------------------------------------------------------
   subroutine MetabolicLoss(biomass, temp, rloss, rrloss)
      implicit none
      real(r8), intent(in) :: biomass(NPOC)  ! umol C m-3
      real(r8), intent(in) :: temp           ! water temperature (K)
      real(r8), intent(out) :: rloss(NPOC)   ! loss rate (umol C m-3 s-1)
      real(r8), intent(out) :: rrloss(NPOC)  ! respiration loss
      real(r8) :: Klr(NPOC), Tw

      Tw = temp - T0
      Klr = (/sa_params(Param_Klrs), sa_params(Param_Klrl)/)
      rloss = Klr / SECOND_OF_DAY * (ThetaML**(Tw-20)) * biomass 
      rrloss = Fres * rloss
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate algae mortality rate due to hypoxia.
   !
   !------------------------------------------------------------------------------
   subroutine HypoxiaDeath(biomass, temp, rdeath)
      implicit none
      real(r8), intent(in) :: biomass(:)     ! umol C m-3
      real(r8), intent(in) :: temp           ! water temperature (K)
      real(r8), intent(out) :: rdeath(:)     ! loss rate (umol C m-3 s-1)
      real(r8) :: Tw

      Tw = temp - T0
      rdeath = Kdhyp / SECOND_OF_DAY * (ThetaML**(Tw-20)) * biomass
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate DOC microbial mineralization rate.
   !
   !------------------------------------------------------------------------------
   subroutine MicrobialRespiration(DOC, temp, Do2, rDOC)
      implicit none
      real(r8), intent(in) :: DOC(2)      ! umol C m-3
      real(r8), intent(in) :: temp        ! water temperature (K)
      real(r8), intent(in) :: Do2         ! dissolved O2 (umol/m3)
      real(r8), intent(out) :: rDOC(2)    ! umol C m-3 s-1
      real(r8) :: rRDOM(2), ftemp, fo2

      rRDOM(1) = sa_params(Param_RDOMaq)
      rRDOM(2) = sa_params(Param_RDOMtr)
      ftemp = ThetaCM**(temp-T0-20.0)
      fo2 = Do2 / (Do2 + Ko2CM)
      rDOC = rRDOM / SECOND_OF_DAY * ftemp * fo2 * DOC
   end subroutine

   subroutine UpdateAlgaeDensity(Chl2C, Ipar, dt, rho)
      implicit none
      real(r8), intent(in) :: Chl2C(NPOC)    ! mg chla (mg C)-1
      real(r8), intent(in) :: Ipar           ! PAR (mol/m2/s)
      real(r8), intent(in) :: dt             ! time interval (s)
      real(r8), intent(inout) :: rho(NPOC)   ! density (kg/m3)
      real(r8) :: Vm0(NPOC), drho(NPOC), fpar(NPOC)
      real(r8) :: phAlpha(NPOC), phBeta(NPOC)

      phAlpha = (/sa_params(Param_phAlphas), sa_params(Param_phAlphal)/)
      phBeta = (/sa_params(Param_phBetas), sa_params(Param_phBetal)/)
      Vm0 = (/sa_params(Param_Vchs), sa_params(Param_Vchl)/)
      Vm0 = Vm0 * Chl2C
      fpar =  (1.0 - exp(-phAlpha*Ipar/Vm0)) * exp(-phBeta*Ipar/Vm0)
      drho = (dsc1 * fpar - dsc3) * dt
      rho = rho + drho
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the degradation rate from POM to DOM by hydrolysis
   !           (Hanson et al., 2011).
   !
   !------------------------------------------------------------------------------
   !subroutine POCDecomposition(POC, temp, Do2, rPOC)
   !   implicit none
   !   real(r8), intent(in) :: POC         ! umol/m3
   !   real(r8), intent(in) :: temp        ! water temperature (K)
   !   real(r8), intent(in) :: Do2         ! dissolved O2 (umol/m3)
   !   real(r8), intent(out) :: rPOC       ! umol C m-3 s-1
   !   real(r8) :: rDPOM, ftemp, fo2
   !
   !   rDPOM = sa_params(Param_DPOM) 
   !   ftemp = ThetaPM**(temp-T0-20.0)
   !   fo2 = Do2 / (Do2 + Ko2CM)
   !   rPOC = rDPOM / SECOND_OF_DAY * ftemp * fo2 * POC
   !end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the net carbon fixation rate from Stefan & Fang (1994)
   !
   !------------------------------------------------------------------------------
   !subroutine NetCarbonFixation(Temp, Chla, SRP, Ipar, PCO2)
   !   implicit none
   !   real(r8), intent(in) :: Temp     ! water temperature (K)
   !   real(r8), intent(in) :: Chla     ! mg chla m-3
   !   real(r8), intent(in) :: SRP      ! soluble reactive P (umol m-3)
   !   real(r8), intent(in) :: Ipar     ! mol m-2 s-1
   !   real(r8), intent(out) :: PCO2    ! CO2 fixation rate (umol C m-3 s-1)
   !   real(r8) :: K1, K2               ! mol m-2 s-1
   !   real(r8) :: Pmax                 ! umol C (mg Chla)-1 s-1
   !   real(r8) :: Tw, minL
   !
   !   ! photosynthesis
   !   Tw = Temp - T0
   !   K1 = 1.908333333333d-4 * (1.086**(Tw-20.0))
   !   if (Tw<=10) then
   !      K2 = 1.388888888889d-3
   !   else
   !      K2 = 4.166666666667d-3
   !   end if
   !   minL = Ipar*(1.0+2.0*sqrt(K1/K2))/(Ipar+K1+Ipar*Ipar/K2)
   !   Pmax = 8.333d-2 * (1.036**(Tw-20.0)) * SRP / (SRP+Ksrp)
   !   PCO2 = Pmax * minL * Chla
   !   ! autotrophic respiration
   !   PCO2 = PCO2 - 3765.06 * 1.1574d-6 * (1.047**(Tw-20.0)) * Chla
   !end subroutine

end module bg_utilities_mod
