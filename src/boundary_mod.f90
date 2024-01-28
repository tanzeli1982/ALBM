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
      integer :: nyr, nmon, nday
      integer :: year, month, day
      integer :: nt, ntot, nt_day
      integer :: nytot, nmtot, ndtot
      integer :: JDN0, JDN1, doy
      real(r8) :: Roun

      ! Get the model running date
      if (.NOT. isspinup) then
         call GetSolarDate(time, 3.6d3*(hindx-1), Use_Leap, year, month, day)
         nday = GetDay(3.6d3*(hindx-1)) + 1
         nyr = year - time%year0 + 1
         nmon = 12 * (year-time%year0) - time%month0 + month + 1
      else
         call Date2JDN(time%year1, time%month1, time%day1, JDN0)
         call Date2JDN(time%year1, time%month0, time%day0, JDN1)
         nday = JDN1 - JDN0 + GetDay(3.6d3*(hindx-1)) + 1
         nday = mod(nday-1, 365) + 1
         call JDN2Date(JDN0+nday-1, year, month, day)
         nyr = 1
         nmon = mod(month-time%month1+12, 12) + 1
      end if

      ntot = size(m_airTemp)
      if (hindx<=0 .and. nday<=1) then
         nt = int(mod(hindx+23,24)/forcing_nhour) + 1
      else if (hindx>ntot) then
         nt = 24/forcing_nhour*(nday-2) + int(mod(hindx-1,24)/forcing_nhour) + 1
      else
         nt = 24/forcing_nhour*(nday-1) + int(mod(hindx-1,24)/forcing_nhour) + 1
      end if

      nytot = size(m_aCO2)
      nmtot = size(m_aO3)
      nyr = max(min(nyr, nytot), 1)
      nmon = max(min(nmon, nmtot), 1)
      
      ! air conditions
      Roun = sa_params(Param_Roun)
      m_surfData%temp = m_airTemp(nt)
      m_surfData%RH = m_airRH(nt)
      m_surfData%wind = max(m_airWind(nt), 0.1_r8)
      if (m_surfData%temp>=T0) then
         m_surfData%rainfall = max(0.0_r8, m_airPr(nt))
         m_surfData%snowfall = 0.0_r8
      else
         m_surfData%rainfall = 0.0_r8
         m_surfData%snowfall = max(0.0_r8, Roul/Roun*m_airPr(nt))
      end if
      m_surfData%pressure = m_airPs(nt) 
      m_surfData%sw = m_airSWRad(nt)
      m_surfData%lw = m_airLWRad(nt) 
      m_radPars%qCO2 = m_aCO2(nyr)
      m_radPars%AbO3 = m_aO3(nmon)
      m_radPars%tau550 = m_aAOD(nmon)

      ! hydrology conditions
      ndtot = size(m_dzsurf)
      nday = max(min(nday, ndtot), 1)
      m_surfData%dzsurf = m_dzsurf(nday)
      if (m_srp(nday)>e8) then
         m_surfData%srp = m_srp(nday) 
      else
         m_surfData%srp = lake_info%srp
      end if

      call Date2JDN(2001, 1, 1, JDN0)
      call Date2JDN(2001, month, day, JDN1)
      doy = JDN1 - JDN0 + 1

      ! radiation parameters
      m_radPars%spr = 1.0d-2 * m_surfData%pressure
      m_radPars%tair = m_surfData%temp
      m_radPars%RH = m_surfData%RH
      m_radPars%Latit = lake_info%latitude
      m_radPars%Longit = lake_info%longitude
      m_radPars%dayl = CalcDaylength(doy, lake_info%latitude)
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
      real(r8) :: hour, dzi, dzw, zcos, rhour 
      real(r8) :: Iab0, Iab1, zenith, srd_daily
      real(r8) :: tcc, tair, chla, kext
      real(r8), save :: rTot = 0.0_r8
      integer :: ii, hour0, top

      ! correct solar radiation by measured solar radiation
      !tcc = m_surfData%cloud
      hour0 = mod(hindx-1,24)
      if (mod(hindx-1,forcing_nhour)==0) then
         rTot = 0.0_r8
         do ii = hour0, hour0+forcing_nhour-1, 1
            hour = DBLE(ii) + 0.5
            call CalcClearSkyIrradiance(hour, fgphot, frdif, zenith)
            !call CorrIrradianceByCloud(zenith, tcc, fgphot, frdif)
            call GetIncidentSRD(m_wvln, fgphot, rhour)
            rTot = rTot + rhour
         end do
         rTot = rTot / float(forcing_nhour)
      end if
      if (rTot>e8 .and. m_surfData%sw>e8) then
         hour = DBLE( hour0 ) + 0.5
         call CalcClearSkyIrradiance(hour, fgphot, frdif, zenith)
         !call CorrIrradianceByCloud(zenith, tcc, fgphot, frdif)
         fgphot = (m_surfData%sw/rTot) * fgphot
         call GetIncidentSRD(m_wvln, fgphot, m_surfData%sw_sim)
      else
         m_fsphot = 0.0_r8
         m_Ipar = 0.0_r8
         m_Iab = 0.0_r8
         m_Idw = 0.0_r8
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
         !call CorrIrradianceByReflection(zenith, Rfri, fgphot, frdif, zcos)
         call CorrIrradianceForIce(zenith, Rfri, m_Hice, fgphot, frdif, zcos)
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
      fgphot_tmp = fgphot

      m_Idw(1) = Iab0
      top = m_lakeWaterTopIndex
      kext = sa_params(Param_Feta) * lake_info%kext
      do ii = 1, WATER_LAYER+1, 1
         dzi = m_waterIce(ii) * m_dZw(ii)
         dzw = m_dZw(ii) - dzi
         if (m_waterIce(ii)<1.0) then
            ! irradiance attenuation
            chla = 1.d3 * m_chla(ii)  ! ug/L
            fgphot_tmp = fgphot_tmp * exp(-kext*dzw - abI*dzi)
            if (ii>top .and. m_Idw(ii)<0.01*m_Idw(top)) then
               fgphot = fgphot_tmp   
            else
               fgphot = fgphot * exp(-abW*dzw - abI*dzi - aCHL*chla*dzw)
            end if
         else
            fgphot_tmp = fgphot_tmp * exp(-abI*dzi)
            fgphot = fgphot * exp(-abI*dzi)
         end if
         ! scalar irradiance (Fichot & Miler, 2010)
         m_fsphot(:,ii) = ((1.0 - frdif)/zcos + frdif/0.859) * fgphot
         if (m_Hsnow>e8 .or. m_Hgrayice>e8) then
            m_Ipar(ii) = 0.0_r8
         else
            call GetIncidentPAR(m_wvln, m_fsphot(:,ii), m_Ipar(ii))
         end if
         ! incident solar radiation
         call GetIncidentSRD(m_wvln, fgphot_tmp, Iab1)
         m_Iab(ii) = Iab0 - Iab1
         m_Idw(ii+1) = Iab1
         Iab0 = Iab1 
         if (Iab0<e8) then
            m_fsphot(:,ii+1:WATER_LAYER+1) = 0.0_r8
            m_Ipar(ii+1:WATER_LAYER+1) = 0.0_r8
            m_Iab(ii+1:WATER_LAYER+1) = 0.0_r8
            m_Idw(ii+2:WATER_LAYER+2) = 0.0_r8
            exit
         end if
      end do
   end subroutine

end module boundary_mod
