module boundary_mod
!---------------------------------------------------------------------------------
! Purpose: setup the real-time boundary and heat conditions.
!---------------------------------------------------------------------------------
   use data_buffer_mod
   use phy_utilities_mod
   use bg_utilities_mod
   use radiation_mod

   private
   public :: GetBoundaryConditions, GetSolarConditions 
   public :: UpdateCatchmentThermalRegime

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: get lake model boundary conditions
   !
   !------------------------------------------------------------------------------
   subroutine GetBoundaryConditions(time, hindx, isspinup)
      implicit none
      type(SimTime), intent(in) :: time
      integer(i8), intent(in) :: hindx
      logical, intent(in) :: isspinup
      integer :: nyr, nmon, nday, nt
      integer :: year, month, day, dofy
      integer :: JDN0, JDN1
      real(r8) :: Roun, hour, lat, elev
      real(r8) :: Tmean, Tdiff, Tmax, Tmin

      ! Get the model running date
      if (.NOT. isspinup) then
         call GetSolarDate(time, 3.6d3*(hindx-1), year, month, day)
         nday = GetDay(3.6d3*(hindx-1)) + 1
         nyr = year - time%year0 + 1
         nmon = 12 * (year-time%year0) - time%month0 + month + 1
         m_timeHist(hindx) = 10000*year + 100*month + day
      else
         call Date2JDN(time%year1, time%month1, time%day1, JDN0)
         call Date2JDN(time%year1, time%month0, time%day0, JDN1)
         nday = JDN1 - JDN0 + GetDay(3.6d3*(hindx-1)) + 1
         nday = mod(nday-1, 365) + 1
         call JDN2Date(JDN0+nday-1, year, month, day)
         nyr = 1
         nmon = mod(month-time%month1+12, 12) + 1
      end if

      if (trim(forcing_tstep)=='hour') then
         nt = 24 * (nday-1) + mod(hindx-1,24) + 1
      else if (trim(forcing_tstep)=='day') then
         nt = nday
      end if
      hour = mod(hindx-1,24) + 0.5 
      lat = lake_info%latitude
      dofy = GetSolarDay(year, month, day)

      ! air conditions
      Roun = sa_params(Param_Roun)
      Tmean = m_airTemp(nt)
      Tmax = m_airTempMax(nt)
      Tmin = m_airTempMin(nt)
      Tdiff = min(Tmax-Tmean, Tmean-Tmin)
      elev = 1d3 * lake_info%zalt
      call GetHourlyTair(dofy, hour, lat, elev, Tmean, Tdiff, &
         m_surfData%temp)
      m_surfData%RH = m_airRH(nt)
      m_surfData%wind = max(m_airWind(nt), 0.1)
      m_surfData%rainfall = max(0.0, m_airPr(nt)-m_airPrsn(nt))
      m_surfData%snowfall = max(0.0, Roul/Roun*m_airPrsn(nt))
      m_surfData%pressure = m_airPs(nt) 
      m_surfData%sw = m_airSWRad(nt)
      m_surfData%lw = m_airLWRad(nt) 
      m_surfData%Qsi = m_Qsi(nday)
      m_surfData%Qso = m_Qso(nday)
      m_surfData%Qgw = m_Qgw(nday)
      m_surfData%tQsi = m_tQsi(nday)
      m_surfData%dQsi = m_dQsi(nday) 
      m_surfData%DICQsi = m_DICQsi(nday)
      m_surfData%DOCQsi = m_DOCQsi(nday)
      m_surfData%POCQsi = m_POCQsi(nday)
      m_surfData%SRPQsi = m_SRPQsi(nday)
      m_surfData%DOQsi = m_DOQsi(nday) 
      m_radPars%qCO2 = m_aCO2(nyr)
      m_radPars%AbO3 = m_aO3(nmon)
      m_radPars%tau550 = m_aAOD(nmon)

      ! radiation parameters
      m_radPars%spr = 1.0d-2 * m_surfData%pressure
      m_radPars%tair = m_surfData%temp
      m_radPars%RH = m_surfData%RH
      m_radPars%Latit = lake_info%latitude
      m_radPars%Longit = lake_info%longitude
      m_radPars%year = year
      m_radPars%month = month
      m_radPars%day = day
      if (m_radPars%Latit>=0) then
         if (month>=3 .and. month<=8) then
            m_radPars%season = 1
         else
            m_radPars%season = 0
         end if
      else
         if (month>=3 .and. month<=8) then
            m_radPars%season = 0
         else
            m_radPars%season = 1
         end if
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: get water column radiation
   !
   !          Important: for the calculation of irradiance attenuation (vertical 
   !          attenuate coeffcient) and energy balance, the downward irradiance 
   !          should be used!!
   !
   !------------------------------------------------------------------------------
   subroutine GetSolarConditions(hindx) 
      implicit none
      integer(i8), intent(in) :: hindx
      real(r8) :: LPOC(NPOC), Chla(NPOC)
      real(r8) :: PPOC, MPOC, DPOC, trDOC, rhour
      real(r8) :: hour, dzi, dzw, zcos, rBsc
      real(r8) :: Iab0, Iab1, zenith, srd_daily
      real(r8) :: tcc, tair
      real(r8), save :: rTot = 0.0_r8
      integer :: ii

      if (m_surfData%sw<e8) then
         m_fsphot = 0.0_r8
         m_Iab = 0.0_r8
         m_surfData%sw_sim = 0.0_r8
         m_surfData%srd = 0.0_r8
         return
      end if

      ! correct solar radiation by measured solar radiation
      !tcc = m_surfData%cloud
      if (trim(forcing_tstep)=='hour') then
         hour = DBLE( mod(hindx-1,24) ) + 0.5
         call CalcClearSkyIrradiance(hour, fgphot, frdif, zenith)
         !call CorrIrradianceByCloud(zenith, tcc, fgphot, frdif)
         call GetIncidentSRD(m_wvln, fgphot, rTot)
      else if (trim(forcing_tstep)=='day') then
         if (mod(hindx-1,24)==0) then
            rTot = 0.0_r8
            do ii = 0, 23, 1
               hour = DBLE(ii) + 0.5
               call CalcClearSkyIrradiance(hour, fgphot, frdif, zenith)
               !call CorrIrradianceByCloud(zenith, tcc, fgphot, frdif)
               call GetIncidentSRD(m_wvln, fgphot, rhour)
               rTot = rTot + rhour
            end do
            rTot = rTot / 24.0
         end if
         hour = DBLE( mod(hindx-1,24) ) + 0.5
         call CalcClearSkyIrradiance(hour, fgphot, frdif, zenith)
         !call CorrIrradianceByCloud(zenith, tcc, fgphot, frdif)
      end if
      if (rTot>e8) then
         fgphot = (m_surfData%sw/rTot) * fgphot
         call GetIncidentSRD(m_wvln, fgphot, m_surfData%sw_sim)
      else
         m_fsphot = 0.0_r8
         m_Iab = 0.0_r8
         m_surfData%sw_sim = 0.0_r8
         m_surfData%srd = 0.0_r8
         return
      end if

      if (m_Hsnow>e8) then
         ! correct for snow reflection
         tair = m_surfData%temp
         call CorrIrradianceForSnow(zenith, tair, fgphot, frdif, zcos)
      else if (m_Hgrayice>e8) then
         ! correct for gray ice reflection
         call CorrIrradianceForGrayIce(zenith, Rfre, fgphot, frdif, zcos)
      else if (m_Hice>e8) then
         ! correct for ice reflection
         call CorrIrradianceByReflection(zenith, Rfri, fgphot, frdif, zcos)
      else
         ! correct for water reflection
         call CorrIrradianceByReflection(zenith, Rfrw, fgphot, frdif, zcos)
      end if
      call GetIncidentSRD(m_wvln, fgphot, m_surfData%srd)

      ! correct for snow and gray ice absorption
      if (m_Hsnow>e8 .or. m_Hgrayice>e8) then
         fgphot = fgphot * exp(-abN*m_Hsnow-abE*m_Hgrayice)   
      end if
      call GetIncidentSRD(m_wvln, fgphot, Iab0)

      do ii = 1, WATER_LAYER+1, 1
         ! scalar irradiance (Fichot & Miler, 2010)
         m_fsphot(:,ii) = ((1.0 - frdif)/zcos + frdif/0.859) * fgphot

         dzi = m_waterIce(ii) * m_dZw(ii)
         dzw = m_dZw(ii) - dzi
         if (m_waterIce(ii)<1.0) then
            ! absorption and scattering coefficients
            LPOC = m_waterPOC(:,ii)
            DPOC = 0.0_r8  ! dead phytoplankton biomass
            ! if DPOC included, 0.5 should be replace by the fraction of
            ! small and large phytoplankton biomass fraction
            PPOC = LPOC(small_ppk) + DPOC * 0.5
            MPOC = LPOC(large_ppk) + DPOC * 0.5
            trDOC = m_waterSubCon(Wtrdoc,ii)
            Chla = 1.0d-3 * m_rChl2C(:,ii) * LPOC 
            !call CalcAcCDOM(lake_info%itype, trDOC, m_wvln, abCDOM) 
            !call CalcAcAlgae(LPOC, Chla, m_wvln, mem_pico, mem_micro, abAP)
            ! irradiance attenuation
            fgphot = fgphot * exp(-lake_info%kext*dzw-abI*dzi)
         else
            fgphot = fgphot * exp(-abI*dzi)
         end if
         call GetIncidentSRD(m_wvln, fgphot, Iab1)
         m_Iab(ii) = Iab0 - Iab1
         Iab0 = Iab1
         if (m_Hice>0.2 .and. m_waterIce(ii)<e8) then
            m_Iab(ii) = 0.0_r8
         end if
         ! absorbed by wet dark sediments
         if (ii==WATER_LAYER+1) then
            m_Iab(ii+1) = 0.92 * Iab0
         end if
         if (Iab0<e8) then
            m_fsphot(:,ii+1:WATER_LAYER+1) = 0.0_r8
            m_Iab(ii+1:WATER_LAYER+2) = 0.0_r8
            exit
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   !  Purpose: update the seasonal flags
   !           
   !------------------------------------------------------------------------------
   subroutine UpdateSeasonalFlags(time, hindx, year, month, day)
      implicit none
      type(SimTime), intent(in) :: time
      integer(i8), intent(in) :: hindx
      integer, intent(out) :: year
      integer, intent(out) :: month
      integer, intent(out) :: day
      real(r8) :: dayl, prev_dayl
      integer :: dofy

      call GetSolarDate(time, 3.6d3*(hindx-1), year, month, day) 
      dofy = GetSolarDay(year, month, day)
      dayl = CalcDaylength(dofy, lake_info%latitude)
      prev_dayl = CalcDaylength(max(1,dofy-1), lake_info%latitude)
      ! set flag for solstice period (winter->summer = 1, summer->winter = 0)
      if (dayl>28800.0_r8 .and. m_Hice<e8) then
         if (dayl>=prev_dayl) then
            winter_flag = 0
         else if (lake_info%latitude>=0) then
            if (month==7 .or. month==8) then
               winter_flag = 0
            else
               winter_flag = prewinter_flag
            end if
         else
            if (month==12 .or. month==1) then
               winter_flag = 0
            else
               winter_flag = prewinter_flag
            end if
         end if
      else
         if (dayl<prev_dayl) then
            winter_flag = 1
         else
            winter_flag = prewinter_flag
         end if
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   !  Purpose: update the active layer thickness of lake catchement soils
   !           
   !------------------------------------------------------------------------------
   subroutine UpdateCatchmentThermalRegime(time, hindx)
      implicit none
      type(SimTime), intent(in) :: time
      integer(i8), intent(in) :: hindx
      integer :: year, month, day
      integer :: ii
      real(r8) :: rsdl 

      ! update daily
      if (mod(hindx,24)/=1) then
         return
      end if

      ! update seasonal flags
      call UpdateSeasonalFlags(time, hindx, year, month, day)
      if (winter_flag==0 .and. prewinter_flag==1) then
         prewinter_flag = winter_flag
         !m_rChl2C = 0.24
         !m_waterPOC = 5d2 * m_chla0 / 0.24 
      else if (winter_flag==1 .and. prewinter_flag==0) then
         prewinter_flag = winter_flag
         do ii = 1, NPOC, 1
            rsdl = max(m_sinkPOCPool(ii)-2.5d3, 0.0)
            m_burialAtCarb = m_burialAtCarb + rsdl
            m_sinkPOCPool(ii) = 2.5d3
         end do
      end if
   end subroutine

end module boundary_mod
