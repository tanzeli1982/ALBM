module phy_utilities_mod
!---------------------------------------------------------------------------------
! Purpose:
!
! This module serves for the implementation of some basic physical utilities
!
!---------------------------------------------------------------------------------
   use shr_kind_mod,          only: r8, i8
   use shr_ctrl_mod,          only: e8 => SHR_CTRL_E8, inf => INFINITE_E8, &
                                    inft => INFINITESIMAL_E8
   use phy_const_mod
   use shr_typedef_mod,       only: SimTime, LakeInfo

   interface CalcPistonVelocity
      module procedure CalcPistonVelocitySR
      module procedure CalcPistonVelocityCC
      module procedure CalcPistonVelocityUW
   end interface

   interface CalcLatentHeatWater
      module procedure CalcLatentHeatWaterAero
      module procedure CalcLatentHeatWaterPM
   end interface

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: Get which day t is in a year and which hour t is in a day
   !
   !------------------------------------------------------------------------------
   subroutine GetSolarDate(base_time, t, useleap, year, month, day)
      implicit none
      type(SimTime), intent(in) :: base_time
      real(r8), intent(in) :: t
      logical, intent(in) :: useleap
      integer, intent(out) :: year
      integer, intent(out) :: month
      integer, intent(out) :: day
      integer :: JDN, nyear, nday
      integer :: year0, month0, day0

      if (useleap) then
         call Date2JDN(base_time%year0, base_time%month0, &
               base_time%day0, JDN)
         JDN = JDN + GetDay(t)
         call JDN2Date(JDN, year, month, day)
      else
         nyear = INT(GetDay(t)/365.)
         year0 = base_time%year0 + nyear
         month0 = base_time%month0
         day0 = base_time%day0
         nday = GetDay(t) - 365*nyear
         call Date2JDN(year0, month0, day0, JDN)
         JDN = JDN + nday
         call JDN2Date(JDN, year, month, day) 
      end if
   end subroutine

   function GetSolarDay(year, month, day, useleap)
      implicit none
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: day
      logical, intent(in) :: useleap
      integer :: GetSolarDay
      logical :: leap
      real(r8) :: Xdm

      if (useleap) then
         leap = IsLeapYear(year)
      else
         leap = .False.
      end if
      if (leap .and. month>2) then
         Xdm = 31.8
      else if (month<=2) then
         Xdm = 30.6
      else
         Xdm = 32.8
      end if
      GetSolarDay = INT(30.6*month + day - Xdm + 0.5)
      return
   end function

   function GetDay(t)
      implicit none
      real(r8), intent(in) :: t     ! units: second
      integer :: GetDay

      GetDay = floor(t/86400.0_r8)
      return
   end function

   function GetHour(t)
      implicit none
      real(r8), intent(in) :: t     ! units: second
      real(r8) :: GetHour
      integer :: day

      day = GetDay(t)
      GetHour = (t-day*86400.0_r8)/3600.0_r8
      return
   end function

   subroutine YYMMDD2Date(YYYYMMDD, year, month, day)
      implicit none
      integer, intent(in) :: YYYYMMDD
      integer, intent(out) :: year
      integer, intent(out) :: month
      integer, intent(out) :: day

      year = INT( YYYYMMDD / 1.0d4 )
      month = INT( ( YYYYMMDD - 1.0d4*year ) / 1.0d2 )
      day = INT( YYYYMMDD - 1.0d4*year - 1.0d2*month )
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Check whether the input is a leap year
   !
   !------------------------------------------------------------------------------
   function IsLeapYear(year)
      implicit none
      integer, intent(in) :: year
      logical :: IsLeapYear

      if (mod(year,4)/=0) then
         IsLeapYear = .false.
      else if (mod(year,100)/=0) then
         IsLeapYear = .true.
      else if (mod(year,400)==0) then
         IsLeapYear = .true.
      else
         IsLeapYear = .false.
      end if
      return
   end function

   function CalcRunningDays(time, useleap)
      implicit none
      type(SimTime), intent(in) :: time
      logical, intent(in) :: useleap
      integer :: CalcRunningDays
      integer, parameter :: ndays(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
      integer :: JDN0, JDN1

      if (useleap) then
         call Date2JDN(time%year0, time%month0, time%day0, JDN0)
         call Date2JDN(time%year1, time%month1, time%day1, JDN1)
         CalcRunningDays = JDN1 - JDN0
      else
         if (time%year1>time%year0) then
            CalcRunningDays = 365*(time%year1-time%year0-1) + &
               sum(ndays(time%month0:12)) + sum(ndays(1:time%month1-1)) - &
               time%day0 + time%day1
         else if (time%month1>time%month0) then
            CalcRunningDays = sum(ndays(time%month0:time%month1-1)) - &
               time%day0 + time%day1
         else
            CalcRunningDays = time%day1 - time%day0
         end if
      end if
      return
   end function

   function CalcRunningMonths(time)
      implicit none
      type(SimTime), intent(in) :: time
      integer :: CalcRunningMonths

      if (time%day1==1) then
         CalcRunningMonths = 12*(time%year1-time%year0) - &
            time%month0 + time%month1
      else
         CalcRunningMonths = 12*(time%year1-time%year0) - &
            time%month0 + time%month1 + 1
      end if
      return
   end function

   function CalcRunningYears(time)
      implicit none
      type(SimTime), intent(in) :: time
      integer :: CalcRunningYears

      if (time%month1==1 .and. time%day1==1) then
         CalcRunningYears = time%year1 - time%year0
      else
         CalcRunningYears = time%year1 - time%year0 + 1
      end if
      return
   end function

   function GetUTCHourIndex(hindx, longitude)
      implicit none
      integer(i8), intent(in) :: hindx    ! local time hour index
      real(r8), intent(in) :: longitude   ! degree
      integer(i8) :: GetUTCHourIndex      ! UTC hour index
      
      GetUTCHourIndex = hindx - NINT(longitude/15.0)
      return
   end function
   
   !------------------------------------------------------------------------------
   !
   ! Purpose: Convert between Julian day number and Gregorian calendar
   !          see http://en.wikipedia.org/wiki/Julian_day
   !
   !------------------------------------------------------------------------------
   subroutine Date2JDN(year, month, day, JDN)
      implicit none
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: day
      integer, intent(out) :: JDN
      integer :: a, y, m
      real(r8) :: date

      a = (14-month) / 12
      y = year + 4800 - a
      m = month + 12*a - 3
      date = 1.d4 * year + 1.d2 * month + day
      if (date>=15821015) then
         JDN = day + (153*m+2)/5 + 365*y + y/4 - y/100 + y/400 - 32045
      else
         JDN = day + (153*m+2)/5 + 365*y + y/4 - 32083
      end if
   end subroutine 

   subroutine JDN2Date(JDN, year, month, day)
      implicit none
      integer, intent(in) :: JDN
      integer, intent(out) :: year
      integer, intent(out) :: month
      integer, intent(out) :: day
      integer :: nf, ne, ng, nh

      nf = JDN + 1401 + (((4*JDN+274277)/146097)*3)/4 - 38
      ne = 4*nf + 3
      ng = mod(ne,1461)/4
      nh = 5*ng + 2
      day = mod(nh,153)/5 + 1
      month = mod(nh/153+2,12) + 1
      year = ne/1461 - 4716 + (14-month)/12
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Convert a year-based floating time to Gregorian calendar
   !
   !------------------------------------------------------------------------------
   subroutine CalendarConversion(t, useleap, year, month, day)
      implicit none
      real(r8), intent(in) :: t
      logical, intent(in) :: useleap
      integer, intent(out) :: year
      integer, intent(out) :: month
      integer, intent(out) :: day
      integer :: nday, JDN
      logical :: isleap

      year = INT(t)
      isleap = IsLeapYear(year)
      if (isleap) then
         if (useleap) then
            nday = NINT(3.66d2*(t-year))
            call Date2JDN(year,1,1,JDN)
            JDN = JDN + nday
            call JDN2Date(JDN, year, month, day)
         else
            nday = NINT(3.65d2*(t-year))
            call Date2JDN(year+1,1,1,JDN)
            JDN = JDN + nday
            call JDN2Date(JDN, year, month, day)
            year = year - 1
         end if
      else
         nday = NINT(3.65d2*(t-year))
         call Date2JDN(year,1,1,JDN)
         JDN = JDN + nday
         call JDN2Date(JDN, year, month, day)
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Get a GMT date and time
   !
   !------------------------------------------------------------------------------
   subroutine GetGMTime(useleap, time)
      implicit none
      logical, intent(in) :: useleap
      integer, intent(out) :: time(6)
      integer :: date_time(8), day
      logical :: leap

      call date_and_time(values=date_time)
      if (useleap) then
         leap = IsLeapYear(date_time(1))
      else
         leap = .False.
      end if
      date_time(5) = date_time(5) - date_time(4)/60   ! adjust hour to GMT zone
      ! The adjust algorithm is actually the addition/subtraction about date
      if (date_time(5)>=24) then
         date_time(5) = date_time(5) - 24
         date_time(3) = date_time(3) + 1
         if (leap .and. date_time(2)==2) then
            day = DAY_IN_MONTH(date_time(2)) + 1
         else
            day = DAY_IN_MONTH(date_time(2))
         end if
         if (date_time(3)>day) then
            date_time(3) = date_time(3) - day
            date_time(2) = date_time(2) + 1
            if (date_time(2)>12) then
               date_time(2) = date_time(2) - 12
               date_time(1) = date_time(1) + 1
            end if
         end if
      else if (date_time(5)<0) then
         date_time(5) = date_time(5) + 24
         date_time(3) = date_time(3) - 1
         if (date_time(3)<=0) then
            date_time(2) = date_time(2) - 1
            if (date_time(2)<=0) then
               date_time(2) = date_time(2) + 12
               date_time(1) = date_time(1) - 1
            end if
            if (leap .and. date_time(2)==2) then
               day = DAY_IN_MONTH(date_time(2)) + 1
            else
               day = DAY_IN_MONTH(date_time(2))
            end if
            date_time(3) = date_time(3) + day
         end if
      end if
      time(1) = date_time(1)  ! year
      time(2) = date_time(2)  ! month
      time(3) = date_time(3)  ! day
      time(4) = date_time(5)  ! hour
      time(5) = date_time(6)  ! minute
      time(6) = date_time(7)  ! second
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate thermal diffusivity of water (liquid or ice)
   !          "Thermal diffusivity of thermokarst lake ice in the Beiluhe basin 
   !           of the Qinghai–Tibetan Plateau"
   !
   !------------------------------------------------------------------------------
   function CalcWaterHeatDiffusivity(temp)
      implicit none
      real(r8), intent(in) :: temp           ! temperature (K)
      real(r8) :: CalcWaterHeatDiffusivity   ! units: m^2/s
      real(r8) :: Tw, param

      Tw = temp - T0
      if (Tw>e8) then
         CalcWaterHeatDiffusivity = (0.558 + 2.223d-3*Tw - 1.797d-5*(Tw**2)) / 4.2d6
      else if (Tw<e8 .and. Tw>-5._r8) then
         param = Tw / -5._r8 
         CalcWaterHeatDiffusivity = Ktw_ref*(1.0 - param) + Kti_ref*param 
      else
         CalcWaterHeatDiffusivity = 1.16 * (1.91 - 8.66d-3*Tw + &
            2.97d-5*(Tw**2)) / 1.9551d6
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate molecular heat conductivity of some common materials
   !          1) the equations for frozen and unfrozen saturated soils are from
   !             Omar T. Farouki (1981);
   !          2) the equation for phase-changing soils is from Klas Hansson (2004).
   !
   !------------------------------------------------------------------------------
   function CalcLWHeatConductivity(temp, cita, icita)
      implicit none
      real(r8), intent(in) :: temp           ! temperature (K)
      real(r8), intent(in) :: cita           ! water content
      real(r8), intent(in) :: icita          ! ice content
      real(r8) :: CalcLWHeatConductivity     ! units: W/(m*K)
      real(r8) :: Tw

      Tw = temp - T0
      if (Tw < e8 .and. Tw > -e8) then
         CalcLWHeatConductivity = 1.0/(cita/Kw0+icita/Ki0)
      else if (Tw > e8) then
         CalcLWHeatConductivity = 0.558 + 2.223d-3*Tw - 1.797d-5*(Tw**2)
      else
         CalcLWHeatConductivity = 1.16*(1.91 - 8.66d-3*Tw + 2.97d-5*(Tw**2))
      end if
      return
   end function

   function CalcSedHeatConductivity(poro, satLW, Kss)
      implicit none
      real(r8), intent(in) :: poro           ! porosity
      real(r8), intent(in) :: satLW          ! liquid water saturation
      real(r8), intent(in) :: Kss            ! solid thermal conductivity
      real(r8) :: CalcSedHeatConductivity    ! units: W/m/K
      real(r8) :: satIce

      satIce = 1.0 - satLW
      CalcSedHeatConductivity = (Kw0**(poro*satLW)) * (Ki0**(poro*satIce)) &
            * (Kss**(1.0-poro))
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: liquid water saturation to subfreezing soil temperature.
   !
   !------------------------------------------------------------------------------
   function CalcSoilWaterSaturation(temp)
      implicit none
      real(r8), intent(in) :: temp  ! units: K
      real(r8) :: CalcSoilWaterSaturation
      real(r8), parameter :: aa = 0.15
      real(r8), parameter :: bb = -0.4

      if (temp<T0) then
         CalcSoilWaterSaturation = min( 1.0, aa*((T0 - temp)**bb) )
      else
         CalcSoilWaterSaturation = 1.0_r8
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate saturation vapor pressure (units: mb) and specific humidity, 
   !          which are from "The computation of equivalent potential temperature", 
   !          Monthly Weather Review, [Bolton, D., 1980]
   !          
   !          0.3% bias within -35C to 35C
   !------------------------------------------------------------------------------
   function CalcSatVP(T)
      implicit none
      real(r8), intent(in) :: T     ! units: K
      real(r8) :: CalcSatVP         ! units: mb
      real(r8) :: Tc

      Tc = T - T0
      Tc = max(min(Tc, 1.0d2), -3.5d1)
      CalcSatVP = 6.112*exp(17.67*Tc/(Tc+243.5))
      return
   end function

   function CalcSatVPSlope(T)
      implicit none
      real(r8), intent(in) :: T     ! units: K
      real(r8) :: CalcSatVPSlope    ! units: mb/K
      real(r8) :: Tc, es

      Tc = T - T0
      es = CalcSatVP(T)
      CalcSatVPSlope = 4.3026d3*es/(Tc+243.5)**2
      return
   end function

   function CalcSpecificHumidity(vp)
      implicit none
      real(r8), intent(in) :: vp                   ! units: mb
      real(r8) :: CalcSpecificHumidity             ! units: g/g
      real(r8), parameter :: epsilon = 0.622       ! H2O/Air molecular weight

      CalcSpecificHumidity = epsilon*vp/(0.01*P0-(1-epsilon)*vp)
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: calculate the wet-bub temperature from "Moored Observations of 
   !          Precipitation Temperature" [Anderson et al., 1998].
   !
   !------------------------------------------------------------------------------
   function CalcWetBubTemp(temp, RH, pressure)
      implicit none
      real(r8), intent(in) :: temp        ! air temperature (K)
      real(r8), intent(in) :: RH          ! relative humidity (%)
      real(r8), intent(in) :: pressure    ! air pressure (Pa)
      real(r8) :: CalcWetBubTemp          ! K
      real(r8) :: Tair, Pair, ev, lambda
      real(r8) :: Twet, Tnew, GTwet, dGTwet

      Tair = temp - T0
      Pair = 0.01 * pressure
      ev = 0.01 * RH * CalcSatVP(temp)

      ! Newton method (maybe better Newton downhill method)
      Twet = Tair
      lambda = 1.0
      do while (.True.)
         GTwet = 6.6d-4*Pair*(1.+1.15d-3*Twet)*(Tair-Twet) + ev - &
                  6.112*exp(17.67*Twet/(Twet+243.5))
         dGTWet = 6.6d-4*Pair*(-2.3d-3*Twet+1.15d-3*Tair-1.) - 2.6298d+4* &
                  exp(17.67*Twet/(Twet+243.5))/((Twet+243.5)**2)
         Tnew = Twet - lambda * GTwet / dGTWet
         if (abs(Tnew-Twet)<1.d-3) then
            Twet = Tnew
            exit
         end if
         Twet = Tnew
      end do

      CalcWetBubTemp = Twet + T0
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate solar declination, zenith angle and downward solar 
   !          radiation from "A Large-Scale Numerical Model of Sea Ice" 
   !          [Parkinson et al., 1979] and downward long-wave radiation from 
   !          "Simulation of lake evaporation with application to modeling lake 
   !          level variations of Harney-Malheur Lake, Oregon" [Hostetler, 1990]. 
   !          Dwelling Longwave radiation (DLR) can also be evaluated by Eq. (1) 
   !          of "The Global Character of the Flux of Downward Longwave Radiation"
   !          [Stephens et al., 2012]
   !
   !------------------------------------------------------------------------------
   function GetSolarDeclination(day)
      implicit none
      integer, intent(in) :: day       ! DOY
      real(r8) :: GetSolarDeclination  ! units: degree 

      GetSolarDeclination = -23.44*cos(2*(day+10)*Pi/365.0)
      return
   end function

   function GetZenithAngle(day, hour, lat)
      implicit none
      integer, intent(in) :: day       ! solar day
      real(r8), intent(in) :: hour     ! local hour
      real(r8), intent(in) :: lat      ! units: decimal
      real(r8) :: GetZenithAngle
      real(r8) :: hourAngle, decl
      real(r8) :: rlat, rdec

      rlat = lat * Pi / 180.0
      hourAngle = (12-hour)*Pi/12
      decl = GetSolarDeclination(day)
      rdec = decl * Pi / 180.0
      GetZenithAngle = acos(sin(rlat)*sin(rdec)+cos(rlat)*cos(rdec)* &
                        cos(hourAngle))
      return
   end function

   function CalcSunriseHour(day, lat, elev)
      implicit none
      integer, intent(in) :: day    ! DOY
      real(r8), intent(in) :: lat   ! units: decimal 
      real(r8), intent(in) :: elev  ! units: m
      real(r8) :: CalcSunriseHour
      real(r8) :: hourAngle, decl
      real(r8) :: rdecl, rlat, rcorr

      decl = GetSolarDeclination(day)
      if (lat>-90+abs(decl) .and. lat<90-abs(decl)) then
         rlat = lat / 180.0 * Pi
         rdecl = decl * Pi / 180.0
         !rcorr = -0.83 - 2.076 * sqrt(elev) / 60.0 
         !rcorr = rcorr * Pi / 180.0
         rcorr = 0.0_r8
         hourAngle = acos( (sin(rcorr)-sin(rlat)*sin(rdecl)) / &
            (cos(rlat)*cos(rdecl)) )
         CalcSunriseHour = 12 * (1-hourAngle/Pi)
      else
         if (lat*decl>0) then
            CalcSunriseHour = 0.0_r8
         else
            CalcSunriseHour = 12.0_r8
         end if
      end if
      return
   end function

   function CalcDaylength(day, lat)
      implicit none
      integer, intent(in) :: day    ! day of year
      real(r8), intent(in) :: lat   ! units: decimal
      real(r8) :: CalcDaylength     ! units: seconds
      real(r8) :: decl, dayl, rlat
      real(r8) :: rdecl

      rlat = lat * Pi / 180.0
      decl = GetSolarDeclination(day)
      rdecl = decl * Pi / 180.0
      dayl = -(sin(rlat)*sin(rdecl))/(cos(rlat) * cos(rdecl))
      dayl = min(1._r8, max(-1._r8,dayl))
      CalcDaylength = 2.0_r8 * 13750.9871_r8 * acos(dayl)
      return
   end function

   function CalcSolarRadiation(day, hour, lat, vp, cloud)
      implicit none
      integer, intent(in) :: day       ! solar day
      real(r8), intent(in) :: hour     ! local hour
      real(r8), intent(in) :: lat      ! units: decimal
      real(r8), intent(in) :: vp       ! units: mb
      real(r8), intent(in) :: cloud    ! units: fraction
      real(r8) :: CalcSolarRadiation
      real(r8) :: zenith, cosz, rad
      real(r8) :: time, rlat

      zenith = GetZenithAngle(day,hour,lat)
      cosZ = cos(zenith)
      if (cosZ<e8) then
         CalcSolarRadiation = 0.0     ! No sunshine
      else
         rad = Solar*(cosZ**2)/((cosZ+2.7)*vp*0.001+1.085*cosZ+0.1)
         CalcSolarRadiation = rad*(1-0.6*cloud**3)
      end if
      return
   end function

   function CalcThermalRadiation(surfTemp, vp, cloud)
      implicit none
      real(r8), intent(in) :: surfTemp ! units: K
      real(r8), intent(in) :: vp       ! units: mb
      real(r8), intent(in) :: cloud    ! units: fraction
      real(r8) :: CalcThermalRadiation
      real(r8) :: factor

      if (cloud<=0.6) then
         factor = 0.84-(0.1-9.973d-7*vp)*(1.0-cloud)+3.491d-6*vp
      else
         factor = 0.87-(0.175-29.92d-7*vp)*(1.0-cloud)+2.693d-6*vp
      end if
      CalcThermalRadiation = factor*Stefan*(surfTemp**4)
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate latent heat, sensible heat, which are from
   !          "A Large-Scale Numerical Model of Sea Ice" [Parkinson et al., 1979]
   !
   !------------------------------------------------------------------------------
   function GetSpecificLatentHeat4Evap(temp)
      implicit none
      real(r8), intent(in) :: temp              ! units: K
      real(r8) :: GetSpecificLatentHeat4Evap    ! for evaporation
      real(r8) :: SLH, Tw

      ! in the water temperature range from -25 deg to 40 deg
      Tw = temp - T0
      SLH = 2500.8 - 2.36*Tw + 0.0016*Tw**2 - 0.00006*Tw**3
      GetSpecificLatentHeat4Evap = 1.d3 * SLH   ! J/g --> J/kg
      return
   end function

   function CalcLatentHeatWaterAero(waterTemp, airTemp, RH, wind)
      implicit none
      real(r8), intent(in) :: waterTemp      ! units: K
      real(r8), intent(in) :: airTemp        ! units: K
      real(r8), intent(in) :: RH             ! units: %
      real(r8), intent(in) :: wind           ! units: m/s
      real(r8) :: CalcLatentHeatWaterAero
      real(r8) :: vps, vap, qs, q, Lv

      vps = CalcSatVP(waterTemp)
      vap = 0.01 * RH * CalcSatVP(airTemp)
      qs = CalcSpecificHumidity(vps)
      q = CalcSpecificHumidity(vap)
      Lv = GetSpecificLatentHeat4Evap(waterTemp)
      CalcLatentHeatWaterAero = max( Roua*Lv*Ce*wind*(qs-q), 0d0 )
      return
   end function

   function CalcLatentHeatWaterPM(waterTemp, airTemp, RH, wind, ps, Rn)
      implicit none
      real(r8), intent(in) :: waterTemp      ! units: K
      real(r8), intent(in) :: airTemp        ! units: K
      real(r8), intent(in) :: RH             ! units: %
      real(r8), intent(in) :: wind           ! units: m/s
      real(r8), intent(in) :: ps             ! units: pascal
      real(r8), intent(in) :: Rn             ! units: W/m2
      real(r8) :: CalcLatentHeatWaterPM      ! units: W/m2
      real(r8) :: delta, vps, vap, de
      real(r8) :: Ea, gama, Lv

      vps = CalcSatVP(waterTemp)
      vap = 0.01 * RH * CalcSatVP(airTemp) 
      de = 0.1 * (vps - vap)           ! vapor pressure deficit (kPa)
      delta = 0.1 * CalcSatVPSlope(waterTemp)  ! units: kPa/K
      Ea = 6.43*(1+0.536*wind)*de      ! bulk aerodynamic expression
      Lv = GetSpecificLatentHeat4Evap(waterTemp)
      gama = Cpa*1.d-3*ps/(0.622*Lv)   ! pyschrometric constant (kPa/K)
      CalcLatentHeatWaterPM = (delta*Rn+gama*11.574*Ea)/(delta+gama)
      return
   end function

   function CalcLatentHeatIce(waterTemp, airTemp, RH, wind)
      implicit none
      real(r8), intent(in) :: waterTemp      ! units: K
      real(r8), intent(in) :: airTemp        ! units: K 
      real(r8), intent(in) :: RH             ! units: %
      real(r8), intent(in) :: wind           ! units: m/s
      real(r8) :: CalcLatentHeatIce
      real(r8) :: vap, vps, qs, q

      vps = CalcSatVP(waterTemp)
      vap = 0.01 * RH * CalcSatVP(airTemp)
      qs = CalcSpecificHumidity(vps)
      q = CalcSpecificHumidity(vap)
      CalcLatentHeatIce = max( Roua*Ls*Ce*wind*(qs-q), 0d0 )
      return
   end function

   function CalcSensibleHeat(waterTemp, airTemp, wind)
      implicit none
      real(r8), intent(in) :: waterTemp      ! units: K
      real(r8), intent(in) :: airTemp        ! units: K
      real(r8), intent(in) :: wind           ! units: m/s
      real(r8) :: CalcSensibleHeat

      CalcSensibleHeat = Roua*Cpa*Ch*wind*(waterTemp-airTemp)
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate eddy diffusivity, gradient Richardson number, Brunt-
   !          Vaisala frequency, and water density profile, which are from
   !          "Simulation of lake evaporation with application to modeling 
   !          lake level variations of Harney-Malheur Lake, Oregon" [Hostetler 
   !          et al., 1990].
   !          The minimum value of N2 is 7.5d-5 in Fang and Stefan (1996) and
   !          1.0d-9 in Wuest and Lorke (2003). 
   !
   !------------------------------------------------------------------------------
   subroutine CalcEddyDiffusivity(depth, freq, wind, lat, keddy)
      implicit none
      real(r8), intent(in) :: depth(:)       ! units: m
      real(r8), intent(in) :: freq(:)        ! units: s-2
      real(r8), intent(in) :: wind           ! units: m/s
      real(r8), intent(in) :: lat            ! units: decimal
      real(r8), intent(out) :: keddy(:)      ! units: m2/s
      real(r8) :: uts, ubs, ekman, rlat
      real(r8) :: tmp, tmp1, Ri, w10
      real(r8) :: zhs, Nsqrt, maxdepth 
      integer :: ii, nn

      if (wind<0.1_r8) then
         ! on no-wind days
         keddy = 0.0_r8
         return
      end if
      nn = size(depth)
      maxdepth = depth(nn)
      ! surface and bottom friction velocity
      w10 = ConvertWindSpeed10(wind, 2.0d0)
      uts = 1.2d-3 * wind
      ubs = sqrt(Cb1m) * 1.0d-2 * w10
      rlat = lat * Pi / 180.0
      ! a latitudinally dependent parameter
      ekman = 6.6*sqrt(sin(abs(rlat)))*(wind**(-1.84))
      do ii = 1, nn, 1
         Nsqrt = max(7.5d-5, freq(ii))
         tmp = uts * exp(-ekman*depth(ii))
         tmp1 = (40*Nsqrt*(Karman*depth(ii))**2) / (tmp**2+inft)
         Ri = 0.05 * (-1 + sqrt(1.0+tmp1))
         keddy(ii) = (Karman*depth(ii)*tmp/Prandtl)/(1+37*Ri**2)
         if (maxdepth>5._r8) then
            zhs = maxdepth - depth(ii) + DeltaD
            tmp = ubs * exp(-ekman*zhs)
            tmp1 = (40*Nsqrt*(Karman*zhs)**2) / (tmp**2+inft)
            Ri = 0.05 * (-1 + sqrt(1.0+tmp1))
            keddy(ii) = keddy(ii) + (Karman*zhs*tmp/Prandtl)/(1+37*Ri**2)
         end if
      end do
   end subroutine

   ! The square of Brunt-Vaisala frequency
   subroutine CalcBruntVaisalaFreq(deltaH, density, freq)
      implicit none
      real(r8), intent(in) :: deltaH(:)
      real(r8), intent(in) :: density(:)
      real(r8), intent(out) :: freq(:)
      real(r8) :: drho
      integer :: ii, nn

      nn = size(deltaH)
      do ii = 1, nn, 1
         if (ii==1) then
            drho = 0.5 * (density(ii+1) - density(ii))
         else if (ii==nn) then
            drho = 0.5 * (density(ii) - density(ii-1))
         else
            drho = 0.5 * (density(ii+1) - density(ii-1))
         end if
         freq(ii) = G * drho / deltaH(ii) / density(ii)
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: KPP-based eddy diffusivity scheme.
   !     Zhang, Q., et al.: Improving lake mixing process simulations in the 
   !     Community Land Model by using K profile parameterization, Hydrol. Earth 
   !     Syst. Sci., 23, 4969–4982, 2019.
   !
   !------------------------------------------------------------------------------
   subroutine CalcKPPDiffusivity(depth, temp, freq, wrho, w10, lat, Heff, Keddy)
      implicit none
      real(r8), intent(in) :: depth(:)       ! units: m
      real(r8), intent(in) :: temp(:)        ! units: K
      real(r8), intent(in) :: freq(:)        ! units: s-2
      real(r8), intent(in) :: wrho(:)        ! units: kg/m3
      real(r8), intent(in) :: w10            ! units: m/s
      real(r8), intent(in) :: lat            ! units: decimal
      real(r8), intent(in) :: Heff           ! units: W/m2
      real(r8), intent(out) :: keddy(:)      ! units: m2/s
      real(r8), dimension(size(depth)) :: Vh_z     ! horizontal velocity of water
      real(r8), dimension(size(depth)) :: Buoy_z   ! buoyancy
      real(r8), dimension(size(depth)) :: nv_z     ! total diffusivity
      real(r8) :: uts, eps, wind, kml
      real(r8) :: Ri0, Ri, Rig, dRih
      real(r8) :: zbnd, fk, Zmax, Vh_r
      real(r8) :: Vt2, ws, phi, xi, dVh
      real(r8) :: Cd10, Lscale, Bf
      real(r8) :: sigma, Gsigma, Gsigma1
      real(r8) :: dGsigma1, ws1, dnv1, dws1
      real(r8) :: phi0, phi1
      real(r8) :: a0, a1, a2, a3
      integer :: ii, nn, ibnd

      nn = size(depth)
      Zmax = maxval(depth)
      wind = ConvertWindSpeed(w10, 2.0_r8)
      Cd10 = 1.0d-3*(2.7/w10+0.142+0.0764*w10)
      uts = sqrt(Cd10*Roua/Roul)*w10
      ! calculate bulk Richardson number
      eps = 0.1_r8
      do ii = 1, nn, 1
         Buoy_z(ii) = G * (1._r8 - wrho(1)/wrho(ii))
         Vh_z(ii) = 0.028_r8 * wind * (3.0*(depth(ii)/Zmax)**2 - &
            4.0*(depth(ii)/Zmax) + 1.0)
      end do
      kml = 1.4d-7
      do ii = 1, nn, 1
         if (ii==1) then
            dVh = (Vh_z(ii+1) - Vh_z(ii)) / (depth(ii+1) - depth(ii))
         else if (ii==nn) then
            dVh = (Vh_z(ii) - Vh_z(ii-1)) / (depth(ii) - depth(ii-1))
         else
            dVh = (Vh_z(ii+1) - Vh_z(ii-1)) / (depth(ii+1) - depth(ii-1))
         end if
         Rig = freq(ii) / dVh**2
         if (Rig<0._r8) then
            nv_z(ii) = 1.0d-5
         else if (Rig<0.7_r8) then
            nv_z(ii) = 1.0d-5 * (1.0 - (Rig/0.7)**2)**3
         else
            nv_z(ii) = 0._r8
         end if
         nv_z(ii) = nv_z(ii) + 1.0d-7 + kml
      end do
      Vh_r = 0.028_r8 * wind * (3.0*(0.1/Zmax)**2 - 4.0*(0.1/Zmax) + 1.0)
      if (Heff<e8) then
         Bf = G**2/(temp(1)*wrho(1)*Cpl)
      else
         Bf = CalcBuoyantFlux(temp(1), wrho(1), Heff)
      end if
      Lscale = uts**3._r8 / (Karman*Bf)
      fk = 2.0 * 7.2921d-5 * sin(lat/180.0*Pi)
      ! calculate the boundary layer depth
      zbnd = depth(2)
      Ri0 = 0._r8
      do ii = 3, nn, 1
         xi = depth(ii) / Lscale
         phi = CalcPhiFunction(xi)
         ws = Karman * uts / phi
         Vt2 = 1.6*depth(ii)*sqrt(abs(freq(ii)))*ws*sqrt(0.2/98.96/eps) / &
            0.25 / Karman**2
         Ri = (Buoy_z(ii) - Buoy_z(2))*depth(ii)/((Vh_r-Vh_z(ii))**2 + Vt2)
         if (Ri>=0.25_r8) then
            zbnd = depth(ii-1)*(Ri-0.25)/(Ri-Ri0) + depth(ii)*(0.25-Ri0)/(Ri-Ri0)
            exit
         else
            Ri0 = Ri
         end if
      end do
      if (Lscale>0._r8) then
         zbnd = min(min(zbnd,Lscale),0.7*uts/fk)
      end if
      zbnd = max(min(zbnd,Zmax),depth(2))   ! boundary layer depth
      ibnd = COUNT(depth<=zbnd)
      ! calculate shape function parameters
      a0 = 0._r8
      a1 = 1._r8
      if (ibnd<nn) then
         xi = zbnd / Lscale
         phi = CalcPhiFunction(xi)
         phi0 = CalcPhiFunction(0.999_r8*xi)
         phi1 = CalcPhiFunction(1.001_r8*xi)
         ws1 = Karman * uts / phi
         dws1 = Karman * uts * (1.0/phi1 - 1.0/phi0) / 2.0d-3
         dnv1 = (nv_z(ibnd+1) - nv_z(ibnd-1)) / (depth(ibnd+1) - depth(ibnd-1))
      else
         xi = zbnd / Lscale
         phi = CalcPhiFunction(xi)
         phi0 = CalcPhiFunction(0.999_r8*xi)
         ws1 = Karman * uts / phi
         dws1 = Karman * uts * (1.0/phi - 1.0/phi0) / 1.0d-3
         dnv1 = (nv_z(ibnd) - nv_z(ibnd-1)) / (depth(ibnd) - depth(ibnd-1))
      end if
      Gsigma1 = nv_z(ibnd)/zbnd/ws1
      dGsigma1 = -dnv1/ws1 - nv_z(ibnd)*dws1/zbnd/ws1**2
      a2 = -2.0 + 3.0*Gsigma1 - dGsigma1
      a3 = 1.0 - 2.0*Gsigma1 + dGsigma1

      ! calculate eddy diffusivity
      do ii = 1, nn, 1
         if (depth(ii)<=zbnd) then
            sigma = depth(ii)/zbnd
            if (sigma>eps .and. sigma<1.0 .and. Lscale<0) then
               xi = eps * zbnd / Lscale
            else
               xi = sigma * zbnd / Lscale
            end if
            phi = CalcPhiFunction(xi)
            ws = Karman * uts / phi
            Gsigma = a0 + a1*sigma + a2*sigma**2 + a3*sigma**3
            keddy(ii) = zbnd * ws * Gsigma
         else
            keddy(ii) = nv_z(ii) - kml
         end if
      end do
   end subroutine

   function CalcPhiFunction(xi)
      implicit none
      real(r8), intent(in) :: xi
      real(r8) :: CalcPhiFunction

      if (xi>=0._r8) then
         CalcPhiFunction = 1.0 + 5.0*xi
      else if (xi>=-1._r8) then
         CalcPhiFunction = (1.0 - 16.0*xi)**(-0.5)
      else
         CalcPhiFunction = (-28.86 - 98.96*xi)**(-1./3.)
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate enhanced diffusivity from "Long-term lake water temperature 
   !          and ice cover simulations/measurements" [Fang et al., 1996].
   !          Deprecate this diffusivity because it overestimates the heat transfer. 
   !
   !------------------------------------------------------------------------------
   subroutine CalcEnhancedDiffusivity(freq, ice, Kds)
      implicit none
      real(r8), intent(in) :: freq(:)        ! units: s-2
      real(r8), intent(in) :: ice(:)         ! fraction
      real(r8), intent(out) :: Kds(:)        ! units: m2/s
      real(r8) :: Nsqrt
      integer :: ii, nn

      nn = size(freq)
      do ii = 1, nn, 1
         if (ice(ii)>e8) then
            Kds(ii) = 0.0_r8
         else
            Nsqrt = max(7.5d-5, freq(ii))
            Kds(ii) = 1.04d-8 * (Nsqrt**(-0.43))
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Build a relationship between 10-meter wind speed and wind speed at
   !          any other height, which are from "Gas transfer velocities measured at
   !          low wind speed over a lake" [Crusius et al., 2003]
   !
   !------------------------------------------------------------------------------
   function ConvertWindSpeed(w10, height)
      implicit none
      real(r8), intent(in) :: w10      ! 10-meter wind speed (m/s)
      real(r8), intent(in) :: height   ! aimed height (m)
      real(r8) :: ConvertWindSpeed
      real(r8) :: coef

      coef = 1.0 + sqrt(Cd10)/Karman*log(10.0/height)
      ConvertWindSpeed = w10/coef
      return
   end function

   function ConvertWindSpeed10(wind, height)
      implicit none
      real(r8), intent(in) :: wind     ! wind speed at given height (m/s)
      real(r8), intent(in) :: height   ! height of known wind speed (m)
      real(r8) :: ConvertWindSpeed10
      real(r8) :: coef

      coef = 1.0 + sqrt(Cd10)/Karman*log(10.0/height)
      ConvertWindSpeed10 = wind*coef
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate water extinction coefficient (m^-1), which are from 
   !          "An improved lake model for climate simulations: Model structure, 
   !          evaluation, and sensitivity analyses in CESM1" [Subin et al., 2012]
   !
   !------------------------------------------------------------------------------
   function CalcLightExtinction(depth)
      implicit none
      real(r8), intent(in) :: depth    ! units: m
      real(r8) :: CalcLightExtinction

      CalcLightExtinction = 1.1925_r8*(max(depth,1.0)**(-0.424_r8))
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate gas solubility (henry law constant), which are from
   !          "Compilation of Henry law constants for inorganic and organic
   !          species of potential importance in environmental chemistry (version 3)" 
   !          [Sander, 1999]
   !   For Bunsen coefficient, C[gas] = Bunsen * (P[gas] / (R * T)) * water_content
   !         (if it is pure water, water_content = 1)
   !
   !------------------------------------------------------------------------------
   function CalcHenrySolubility(gas, temp, pH)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: temp        ! units: K
      real(r8), intent(in) :: pH          ! units: n/a
      real(r8) :: CalcHenrySolubility     ! units: mol/(m3*Pa) (M = mole/L)
      real(r8) :: hi, kc1, kc2, par
      integer :: indx

      if (gas==Wn2) then
         CalcHenrySolubility = 6.1d-6*exp(-1300*(1/temp-1/298.0))
      else if (gas==Wo2) then
         if (temp>=T0 .and. temp<=T0+50) then
            indx = int((temp-T0)/5) + 1
            indx = min(indx, 10)
            par = (temp - T0 - 5*indx + 5) / 5.0
            CalcHenrySolubility = (SOLO2(indx+1)*par + SOLO2(indx)*(1-par)) &
                                    / MasO2 / P0 / Xo2 
         else
            CalcHenrySolubility = 1.3d-5*exp(-1500*(1/temp-1/298.0))
         end if
      else if (gas==Wco2) then
         CalcHenrySolubility = 3.4d-4*exp(-2400*(1/temp-1/298.0))
         hi = 10**(-pH)    ! Concentration of hydrogen ion
         ! rate constant of dissolved CO2 for first and second dissolution
         kc1 = 4.3d-7*exp(-921.4*(1/temp-1/298.0))
         kc2 = 4.7d-11*exp(-1787.4*(1/temp-1/298.0))
         CalcHenrySolubility = CalcHenrySolubility*(1.0+kc1/hi+kc1*kc2/hi**2)
      else if (gas==Wch4) then
         CalcHenrySolubility = 1.3d-5*exp(-1700*(1/temp-1/298.0))
      end if
      return
   end function

   function CalcBunsenSolubility(gas, temp, pH)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: temp        ! units: K
      real(r8), intent(in) :: pH          ! units: n/a
      real(r8) :: CalcBunsenSolubility    ! units: m3 / m3
      real(r8) :: henry                   ! units: mol / (m3 * Pa)

      henry = CalcHenrySolubility(gas, temp, pH)
      CalcBunsenSolubility = henry * temp * R
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate water surface tension from Wikipedia "Surface_tension"
   !  Surface tension is expressed as a surface energy per unit area or the work
   !  required to increase unit surface area at constant temperature.
   !
   !------------------------------------------------------------------------------
   function CalcSurfaceTension(temp)
      implicit none
      real(r8), intent(in) :: temp
      real(r8) :: CalcSurfaceTension               ! units: N/m
      real(r8), parameter :: r0 = 75.64d-3, r25 = 71.97d-3, r50 = 67.91d-3
      real(r8) :: par

      if (temp<T0) then
         CalcSurfaceTension = 0.0_r8
      else if (temp>=T0 .and. temp<T0+25) then
         par = (temp-T0)/25.0
         CalcSurfaceTension = r0*(1-par) + r25*par
      else if (temp>=T0+25 .and. temp<T0+50) then
         par = (temp-T0-25.0)/25.0
         CalcSurfaceTension = r25*(1-par) + r50*par
      else
         CalcSurfaceTension = r50
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate water dynamic Viscosity from wikipedia "Viscosity". 
   !          kinematic viscosity = dynamic viscosity / density
   !
   !------------------------------------------------------------------------------
   function CalcDynamicViscosity(temp)
      implicit none
      real(r8), intent(in) :: temp        ! units: K
      real(r8) :: CalcDynamicViscosity    ! units: kg/s/m

      if (temp<T0-e8) then
         CalcDynamicViscosity = inf 
      else
         CalcDynamicViscosity = 2.414d-5 * (10.0**(247.8/(temp-140.0)))
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate buoyant rising speed of bubble from "Bubbles and the 
   !          air-sea exchange of gases in near-saturation conditions" [Woolf, 
   !          D. and Thorpe, S., 1991].
   !
   !------------------------------------------------------------------------------
   function CalcBuoyantVelocity(radius, vsc)
      implicit none
      real(r8), intent(in) :: radius      ! units: m
      real(r8), intent(in) :: vsc         ! kinematic viscosity (m2/s)
      real(r8) :: CalcBuoyantVelocity     ! units: m/s               
      real(r8) :: rr, xx, yy

      if (vsc>1.0d10) then
         CalcBuoyantVelocity = 1.e-30_r8
      else 
         xx = G * (radius**3.0_r8) / (vsc**2.0_r8)
         yy = 10.82_r8 / xx
         CalcBuoyantVelocity = (2.0_r8*(radius**2.0_r8)*G/9.0_r8/vsc) * &
            ((yy**2.0_r8+2.0_r8*yy)**0.5_r8-yy)
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate Schmidt number for N2, O2, CO2 and CH4
   !          Reference: "Relationship between wind-speed and gas-exchange over 
   !          the ocean", R. Wanninkhof, 1992 (N2 and O2) and "Measurement of the 
   !          diffusion coefficients of sparingly soluble gases in water", B. 
   !          Jahne, et al., 1987 (CO2 and CH4).
   !
   !------------------------------------------------------------------------------
   function CalcSchmidtNumber(gas, temp)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: temp     ! units: K
      real(r8) :: CalcSchmidtNumber
      real(r8) :: T

      T = min(temp-T0,30.0)     ! to celcius
      if (T<0) then
         CalcSchmidtNumber = inf
      else
         if (gas==Wn2) then
            CalcSchmidtNumber = 1970.7-131.45*T+4.139*T**2-0.052106*T**3
         else if (gas==Wo2) then
            CalcSchmidtNumber = 1800.6-120.1*T+3.7818*T**2-0.047608*T**3
         else if (gas==Wco2) then
            CalcSchmidtNumber = 1911-113.7*T+2.967*T**2-0.02943*T**3
         else if (gas==Wch4) then
            CalcSchmidtNumber = 1898-110.1*T+2.834*T**2-0.02791*T**3
         end if
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate gas diffusivity in water and air, from Appendix A of J. 
   !          Tang's thesis. In D. Woolf's "Bubbles and the air-sea exchange of 
   !          gases in near saturation conditions", the diffusivity have different 
   !          values for those gases (Table 1):
   !          Dn2 = 1.8d-9
   !          Do2 = 1.7d-9
   !          Dco2 = 1.3d-9
   !          Dch4 = 0.5d-9
   !          For ice: "Change in CO2 concentration and O2/N2 ratio in ice
   !          cores due to molecular diffusion" (Bereiter et al., 2009)
   !
   !------------------------------------------------------------------------------
   function CalcGasDiffusivityInWater(gas, temp)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: temp                 ! units: K
      real(r8) :: CalcGasDiffusivityInWater        ! units: m2/s

      if (gas==Wn2) then
         CalcGasDiffusivityInWater = 2.57d-7*(temp/273.0)
      else if (gas==Wo2) then
         CalcGasDiffusivityInWater = 2.4d-7*(temp/298.0)
      else if (gas==Wco2) then
         CalcGasDiffusivityInWater = 1.81d-3*exp(-2032.6/temp)
      else if (gas==Wch4) then
         CalcGasDiffusivityInWater = 1.5d-7*(temp/298.0)
      end if
      return
   end function

   function CalcGasDiffusivityInAir(gas, temp)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: temp                 ! units: K   
      real(r8) :: CalcGasDiffusivityInAir          ! units: m2/s

      if (gas==Wn2) then
         CalcGasDiffusivityInAir = 1.93d-5*(temp/273.0)**1.82
      else if (gas==Wo2) then
         CalcGasDiffusivityInAir = 1.8d-5*(temp/273.0)**1.82
      else if (gas==Wco2) then
         CalcGasDiffusivityInAir = 1.47d-5*(temp/273.15)**1.792
      else if (gas==Wch4) then
         CalcGasDiffusivityInAir = 1.9d-5*(temp/298.0)**1.82
      end if
      return
   end function

   function CalcGasDiffusivityInIce(gas, temp)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: temp                 ! units: K 
      real(r8) :: CalcGasDiffusivityInIce          ! units: m2/s

      if (gas==Wn2) then
         CalcGasDiffusivityInIce = 2.0d-10*exp(-5.1d+3/R/temp)
      else if (gas==Wo2) then
         CalcGasDiffusivityInIce = 3.5d-9*exp(-9.7d+3/R/temp)
      else if (gas==Wco2) then
         CalcGasDiffusivityInIce = 9.1d-8*exp(-14.7d+3/R/temp)
      else if (gas==Wch4) then
         CalcGasDiffusivityInIce = 2.2d-9*exp(-9.7d+3/R/temp)
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate Nusselt number defined as the ratio between the total gas 
   !          flux and the molecular diffusivity flux across the bubble surface, 
   !          which is from "Bubbles and the air-sea exchange of gases in near-
   !          saturation conditions" [Woolf D. and Thorpe S., 1991].
   !
   !------------------------------------------------------------------------------
   function CalcPecletNumber(gas, radius, temp, vsc)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: radius         ! units: m
      real(r8), intent(in) :: temp           ! units: K
      real(r8), intent(in) :: vsc            ! units: m2/s
      real(r8) :: CalcPecletNumber
      real(r8) :: wb, diffusivity

      wb = CalcBuoyantVelocity(radius,vsc)
      diffusivity = CalcGasDiffusivityInWater(gas,temp)
      CalcPecletNumber = radius*wb/(diffusivity+inft)
      return
   end function

   function CalcReynoldNumber(gas, radius, temp, vsc)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: radius         ! units: m
      real(r8), intent(in) :: temp           ! units: K
      real(r8), intent(in) :: vsc            ! units: m2/s
      real(r8) :: CalcReynoldNumber
      real(r8) :: Pe, Schmidt                ! Pe = Re*Sc

      Pe = CalcPecletNumber(gas,radius,temp,vsc)
      Schmidt = CalcSchmidtNumber(gas,temp)
      CalcReynoldNumber = Pe/Schmidt
      return
   end function

   function CalcNusseltNumber(gas, radius, temp, vsc)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: radius         ! units: m
      real(r8), intent(in) :: temp           ! units: K
      real(r8), intent(in) :: vsc            ! units: ms/s
      real(r8) :: CalcNusseltNumber
      real(r8) :: PeN, ReN

      ReN = CalcReynoldNumber(gas,radius,temp,vsc)
      PeN = CalcPecletNumber(gas,radius,temp,vsc)
      if (ReN<=1.0) then
         CalcNusseltNumber = sqrt(2*Pi*PeN/3)
      else
         CalcNusseltNumber = 0.45*(ReN**(1.0/6.0))*(PeN**(1.0/3.0))
      end if
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate gas transfer velocity from the wind-based model 
   !          [Cole and Caraco, 1998 Limnol. Oceanogr.] and the surface renewal 
   !          model [MacIntyre et al., 2010 GRL; Heiskanen et al., 2014 Tellus B].
   !
   !------------------------------------------------------------------------------
   function CalcPistonVelocityCC(wind)
      implicit none
      real(r8), intent(in) :: wind           ! 10-meter wind (m/s)
      real(r8) :: CalcPistonVelocityCC       ! units: m/s
      real(r8) :: k600

      k600 = 2.778d-6 * (2.07 + 0.215*wind**1.7)   ! cm/hr to m/s
      CalcPistonVelocityCC = k600
      return
   end function

   function CalcPistonVelocitySR(wind, temp, rho0, vv, Heff, zaml) 
      implicit none
      real(r8), intent(in) :: wind           ! 10-m wind (m/s)
      real(r8), intent(in) :: temp           ! units: K
      real(r8), intent(in) :: rho0           ! units: kg/m3
      real(r8), intent(in) :: vv             ! kinematic viscosity
      real(r8), intent(in) :: Heff           ! W/m2
      real(r8), intent(in) :: zaml           ! m
      real(r8) :: CalcPistonVelocitySR       ! m/s
      real(r8) :: k600, uts, beta, epslon

      uts = CalcWaterFrictionVelocity(wind) 
      beta = CalcBuoyantFlux(temp, rho0, Heff) 
      beta = min(beta, 0.0)
      epslon = -0.77*beta + 0.3*uts**3/zaml/Karman
      k600 = 0.5*(vv*epslon)**0.25
      CalcPistonVelocitySR = k600
      return
   end function

   function CalcPistonVelocityUW(info, wind, temp, rho0, Heff, zaml)
      implicit none
      type(LakeInfo), intent(in) :: info
      real(r8), intent(in) :: wind           ! 10-m wind (m/s)
      real(r8), intent(in) :: temp           ! units: K
      real(r8), intent(in) :: rho0           ! units: kg/m3
      real(r8), intent(in) :: Heff           ! W/m2
      real(r8), intent(in) :: zaml           ! m
      real(r8) :: CalcPistonVelocityUW       ! m/s
      real(r8) :: beta, wts, Wstr, k600

      beta = CalcBuoyantFlux(temp, rho0, Heff)
      beta = min(beta, 0.0)
      ! penetrative convective velocity
      wts = (-beta*zaml)**(1.0/3.0)
      Wstr = 1.0 - exp(-0.3*1.0d-6*info%Asurf)
      k600 = sqrt( Wstr*(0.00015*wind)**2 + (0.07*wts)**2 )
      CalcPistonVelocityUW = k600
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate equilibrium gas concentration
   !
   !------------------------------------------------------------------------------
   function CalcEQConc(gas, temp, pH, pressure)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: temp           ! units: Kelvin
      real(r8), intent(in) :: pH
      real(r8), intent(in) :: pressure       ! partial pressure (Pa)
      real(r8) :: CalcEQConc                 ! units: mol/m3
      real(r8) :: solubility

      solubility = CalcHenrySolubility(gas,temp,pH)
      CalcEQConc = solubility*pressure
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate airflow rate 
   !
   !------------------------------------------------------------------------------
   function CalcBubbleAirflow(mflow, area, depth, temp, pr0)
      implicit none
      real(r8), intent(in) :: mflow    ! units: mol/m2/s
      real(r8), intent(in) :: area     ! units: m2
      real(r8), intent(in) :: depth    ! units: m
      real(r8), intent(in) :: temp     ! units: K
      real(r8), intent(in) :: pr0      ! units: Pascal
      real(r8) :: CalcBubbleAirflow    ! units: m3/s
      real(r8) :: pr

      pr = pr0 + Roul * G * depth
      CalcBubbleAirflow = mflow * area * R * temp / pr
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Construct depth and area grids
   !
   !------------------------------------------------------------------------------
   subroutine ConstructDepthVector(info, Zw, dZw, Zs, dZs)
      implicit none
      type(LakeInfo), intent(in) :: info     ! lake information object
      real(r8), intent(out) :: Zw(:)         ! water layer depth vector
      real(r8), intent(out) :: dZw(:)        ! water layer thickness
      real(r8), intent(out) :: Zs(:)         ! sediment layer depth vector
      real(r8), intent(out) :: dZs(:)        ! sediment layer thickness
      real(r8), parameter :: KArr(20) = (/0.0733, 0.0914, 0.1017, &
            0.1088, 0.1143, 0.1188, 0.1225, 0.1257, 0.1285, &
            0.1310, 0.1333, 0.1354, 0.1373, 0.1390, 0.1407, &
            0.1422, 0.1436, 0.1449, 0.1462, 0.1474/)
      real(r8) :: K1, K2, denom, pp
      integer :: nw, ns, nz, ii, idx

      ! for water column
      nw = size(Zw) - 1
      if (info%maxdepth<=5.0) then
         denom = info%maxdepth / dble(nw)
         do ii = 1, nw+1, 1
            Zw(ii) = (ii-1) * denom
         end do
      else
         if (info%maxdepth<=20) then
            K1 = -1.572d-4*info%maxdepth**2 + 6.861d-3*info%maxdepth - 2.686d-2
         else if (info%maxdepth<=50) then
            K1 = -1.372d-5*info%maxdepth**2 + 1.794d-3*info%maxdepth + 1.757d-2
         else if (info%maxdepth<1000) then
            idx = INT(info%maxdepth/50.0)
            pp = (info%maxdepth - 50*idx) / 50.0
            K1 = (1.0-pp)*KArr(idx) + pp*KArr(idx+1)
         else
            K1 = KArr(20) + (info%maxdepth - 1000) / 50.0 * 1.0d-3
         end if
         denom = info%maxdepth / (exp(K1*nw) - 1)
         do ii = 1, nw+1, 1
            Zw(ii) = (exp(K1*(ii-1)) - 1) * denom
         end do
      end if
      do ii = 1, nw+1, 1
         if (ii==1) then
            dZw(ii) = 0.5 * (Zw(ii+1) - Zw(ii))
         else if (ii==nw+1) then
            dZw(ii) = 0.5 * (Zw(ii) - Zw(ii-1))
         else
            dZw(ii) = 0.5 * (Zw(ii+1) - Zw(ii-1))
         end if
      end do
      ! for sediment column
      ns = size(Zs) - 1
      K2 = 0.05_r8
      denom = info%hsed / (exp(K2*ns) - 1)
      do ii = 1, ns+1, 1
         Zs(ii) = (exp(K2*(ii-1)) - 1) * denom
      end do
      do ii = 1, ns+1, 1
         if (ii==1) then
            dZs(ii) = 0.5 * (Zs(ii+1) - Zs(ii))
         else if (ii==ns+1) then
            dZs(ii) = 0.5 * (Zs(ii) - Zs(ii-1))
         else
            dZs(ii) = 0.5 * (Zs(ii+1) - Zs(ii-1))
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Set sediment porosity vector 
   !
   !------------------------------------------------------------------------------
   subroutine SetSedPorosity(Zs, dZs, sedpor)
      implicit none
      real(r8), intent(in) :: Zs(:)          ! sediment layer depth 
      real(r8), intent(in) :: dZs(:)         ! sediment layer thickness
      real(r8), intent(out) :: sedpor(:,:)   ! sediment layer porosity
      integer :: ns, ii
      real(r8) :: tdepth, bdepth, par

      ns = size(Zs)
      do ii = 1, ns, 1
         if (ii==1) then
            tdepth = Zs(ii)
            bdepth = Zs(ii) + dZs(ii)
         else if (ii==ns) then
            tdepth = Zs(ii) - dZs(ii)
            bdepth = Zs(ii)
         else
            tdepth = Zs(ii) - 0.5 * dZs(ii)
            bdepth = Zs(ii) + 0.5 * dZs(ii)
         end if
         if (tdepth<=0.15 .and. bdepth<=0.15) then
            sedpor(:,ii) = 0.8_r8
         else if (tdepth<=0.15) then
            par = (0.15 - tdepth) / (bdepth - tdepth)
            sedpor(:,ii) = 0.4_r8 + 0.4_r8 * par 
         else
            sedpor(:,ii) = 0.4_r8
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate lake sediment bottom temperature 
   !          Uas air temperature climatology data as a reference of sediment 
   !          bottom temperature (Fang et al., 1998)
   !
   !------------------------------------------------------------------------------
   subroutine CalcSedBottomTemp(info, t2mref, margin, tsb)
      implicit none
      type(LakeInfo), intent(in) :: info
      real(r8), intent(in) :: t2mref      ! units: K
      logical, intent(in) :: margin       ! margin flag 
      real(r8), intent(out) :: tsb        ! units: K
      real(r8) :: tgw, zl, geom

      if (info%thrmkst>0 .and. margin) then
         tgw = t2mref - T0 - 0.5_r8       ! ground water temperature
      else
         tgw = t2mref - T0 + 2.0_r8       ! ground water temperature
      end if
      geom = (info%Asurf)**0.25 / info%depth    ! lake geometry ratio 
      geom = max( min(geom,4.0), 0.8 )    ! Fig. 9
      tsb = T0 + 0.5 * (tgw - 2) + 0.5 * (tgw + 6) * log(geom)
   end subroutine

   function CalcGreatCircleDistance(lon0, lat0, lon1, lat1)
      implicit none
      real(r8), intent(in) :: lon0, lat0  ! units: decimal
      real(r8), intent(in) :: lon1, lat1  ! units: decimal
      real(r8) :: CalcGreatCircleDistance
      real(r8), parameter :: radius =  6371.009_r8   ! units: km
      real(r8) :: dlon, dlat, rlon0, rlon1
      real(r8) :: rlat0, rlat1, dcigma

      rlon0 = lon0 / 180.0 * PI
      rlon1 = lon1 / 180.0 * PI
      rlat0 = lat0 / 180.0 * PI
      rlat1 = lat1 / 180.0 * PI
      dlon = abs(rlon0 - rlon1)
      dlat = abs(rlat0 - rlat1)
      dcigma = 2.0 * asin( sqrt( sin(0.5*dlat)**2 + cos(rlat0) * &
               cos(rlat1) * sin(0.5*dlon)**2 ) )
      CalcGreatCircleDistance = radius * dcigma
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate settling velocity of detritus and algae (Stoke's Law).
   !           (Hanson et al., 2011)
   !
   !------------------------------------------------------------------------------
   function CalcSettlingVelocity(rhow, rhos, da, dVsc)
      implicit none
      real(r8), intent(in) :: rhow           ! water density (kg/m3)
      real(r8), intent(in) :: rhos           ! particulate density (kg/m3)
      real(r8), intent(in) :: da             ! particulate diameter (m)
      real(r8), intent(in) :: dVsc           ! viscosity (kg/s/m)
      real(r8) :: CalcSettlingVelocity       ! m/s

      CalcSettlingVelocity = G * da * da * (rhos - rhow) / (18.0*dVsc)
      return
   end function

   !------------------------------------------------------------------------------
   ! 
   ! Purpose: Calculate turbulent kinetic energy dissipation.
   !          "Small-scale hydrodynamics in lakes" (Wuest et al., 2003)
   !
   !------------------------------------------------------------------------------
   function CalcSignificantWaveHeight(w10, area)
      implicit none
      real(r8), intent(in) :: w10      ! units: m/s
      real(r8), intent(in) :: area     ! units: m2
      real(r8) :: CalcSignificantWaveHeight
      real(r8) :: fetch, wstar, C10

      C10 = CalcWindDragCoefficient(w10)
      wstar = sqrt(C10) * w10
      fetch = 2.0 * sqrt(area/Pi)
      CalcSignificantWaveHeight = 0.051*wstar*sqrt(fetch/G)
      return
   end function

   function CalcWindDragCoefficient(w10)
      implicit none
      real(r8), intent(in) :: w10
      real(r8) :: CalcWindDragCoefficient
      real(r8) :: C10

      if (w10<=3.0) then
         C10 = 0.0044 * w10**(-1.15) 
      else if (w10>=5.0) then
         C10 = 5.7616d-4 * w10**0.4018
      else
         C10 = 1.0d-3 
      end if
      CalcWindDragCoefficient = C10
      return
   end function

   function CalcWaterFrictionVelocity(w10)
      implicit none
      real(r8), intent(in) :: w10
      real(r8) :: CalcWaterFrictionVelocity
      real(r8) :: ua

      ua = sqrt(Cd10) * w10
      CalcWaterFrictionVelocity = ua * sqrt(Roua/Roul)
      return
   end function

   !------------------------------------------------------------------------------
   ! 
   ! Purpose: calculate pure water density (Schmid and Wuest, 2012).
   !
   !------------------------------------------------------------------------------
   function CalcPureWaterDensity(temp)
      implicit none
      real(r8), intent(in) :: temp        ! units: celcius
      real(r8) :: CalcPureWaterDensity    ! units: kg/m3

      if (temp>=0) then
         CalcPureWaterDensity = 9.99868d+2 + 1.0d-3*(65.185*temp - &
            8.4878*temp**2. + 0.05607*temp**3.)
      else
         CalcPureWaterDensity = Roui 
      end if
      return
   end function

   !------------------------------------------------------------------------------
   ! 
   ! Purpose: Calculate water thermal expansion coefficent that is the relative
   !          gradient of water density and buoyant flux (Heff is positive when
   !          gaining heat and negative when losing heat) [Heiskanen et al.,
   !          2014].
   !
   !------------------------------------------------------------------------------
   function CalcThermalExpansionCoef(temp)
      implicit none
      real(r8), intent(in) :: temp           ! untis: K
      real(r8) :: CalcThermalExpansionCoef   ! units: fraction/K
      real(r8) :: Tw, denom, numer

      Tw = temp - T0
      denom = 9.99868d+2 + 1.0d-3*(65.185*Tw - 8.4878*Tw**2. + 0.05607*Tw**3.)
      numer = 1.0d-3 * (65.185 - 16.9756*Tw + 0.16821*Tw**2.)
      CalcThermalExpansionCoef = -numer / denom 
      return
   end function

   function CalcBuoyantFlux(temp, rho0, Heff)
      implicit none
      real(r8), intent(in) :: temp     ! mixing layer temperature (K)
      real(r8), intent(in) :: rho0     ! mixing layer density (kg/m3) 
      real(r8), intent(in) :: Heff     ! mixing layer heat flux (W/m2)
      real(r8) :: CalcBuoyantFlux      ! units: m^2 s^-3
      real(r8) :: alpha

      alpha = CalcThermalExpansionCoef(temp)
      CalcBuoyantFlux = G*alpha*Heff/rho0/Cpl
      return
   end function

   !------------------------------------------------------------------------------
   ! 
   ! Purpose: Calculate kinetic energy and potential energy for mixing.
   !          "MyLake-A multi-year lake simulation model code suitable for 
   !          uncertainty and sensitivity analysis simulations" (Saloranta &
   !          Andersen, 2007)
   !
   !------------------------------------------------------------------------------
   subroutine CalcTotalKineticPower(info, w10, Pkin)
      implicit none
      type(LakeInfo), intent(in) :: info
      real(r8), intent(in) :: w10         ! 10-m wind (m/s)
      real(r8), intent(out) :: Pkin       ! units: J/s
      real(r8) :: stress, As

      stress = Roua * Cd10 * (w10**2)  ! wind stress (N/m2)
      As = 1.0d-6 * info%Asurf
      Pkin =  As * sqrt(stress**3/Roul)
   end subroutine

   !------------------------------------------------------------------------------
   ! 
   ! Purpose: Calculate the snow burden capacity of ice layers according to the
   !          Archimedes' theorem (Rogers et al., 1995; Vavrus et al., 1996) 
   !
   !------------------------------------------------------------------------------
   subroutine CalcSnowBurdenLimit(Hice, Hgrayice, roun, snowMax)
      implicit none
      real(r8), intent(in) :: Hice           ! ice thickness (m)
      real(r8), intent(in) :: Hgrayice       ! gray ice thickness (m)
      real(r8), intent(in) :: roun           ! snow density (kg/m3)
      real(r8), intent(out) :: snowMax       ! units: m
      
      snowMax = ((Roul-Roui)*Hice+(Roul-Roue)*Hgrayice)/roun
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Construct hourly temperature profile
   !          "Estimating average daytime and daily temperature profiles within 
   !           Europe" (Huld et al., 2006)
   !
   !------------------------------------------------------------------------------
   subroutine GetHourlyTair(day, hour, lat, elev, Tmean, Tdiff, Thr)
      implicit none
      integer, intent(in) :: day       ! DOY
      real(r8), intent(in) :: hour     ! local hour
      real(r8), intent(in) :: lat      ! units: decimal
      real(r8), intent(in) :: elev     ! units: m
      real(r8), intent(in) :: Tmean    ! mean air temperature
      real(r8), intent(in) :: Tdiff    ! half of temperature variation
      real(r8), intent(out) :: Thr     ! hourly air temperature
      real(r8), parameter :: tpeak = 15.0_r8
      real(r8) :: tdawn

      tdawn = CalcSunriseHour(day, lat, elev)
      if (hour<=tdawn) then
         Thr = Tmean - Tdiff * cos((tdawn-hour)/(24+tdawn-tpeak)*Pi)
      else if (hour<=tpeak) then
         Thr = Tmean + Tdiff * cos((tpeak-hour)/(tpeak-tdawn)*Pi)
      else
         Thr = Tmean - Tdiff * cos((24+tdawn-hour)/(24+tdawn-tpeak)*Pi)
      end if
   end subroutine

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
         pPAR = rPAR * 1.d6 / fconvPAR
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

end module phy_utilities_mod
