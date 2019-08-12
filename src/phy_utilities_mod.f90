module phy_utilities_mod
!---------------------------------------------------------------------------------
! Purpose:
!
! This module serves for the implementation of some basic physical utilities
!
!---------------------------------------------------------------------------------
   use shr_kind_mod,          only: r8
   use shr_ctrl_mod,          only: e8 => SHR_CTRL_E8, inf => INFINITE_E8, &
                                    inft => INFINITESIMAL_E8
   use phy_const_mod
   use shr_typedef_mod,       only: SimTime, LakeInfo
   use shr_param_mod

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
   subroutine GetSolarDate(base_time, t, year, month, day)
      implicit none
      type(SimTime), intent(in) :: base_time
      real(r8), intent(in) :: t
      integer, intent(out) :: year
      integer, intent(out) :: month
      integer, intent(out) :: day
      integer :: JDN

      call Date2JDN(base_time%year0, base_time%month0, &
                    base_time%day0, JDN)
      JDN = JDN + GetDay(t)
      call JDN2Date(JDN, year, month, day)
   end subroutine

   function GetSolarDay(year, month, day)
      implicit none
      integer, intent(in) :: year
      integer, intent(in) :: month
      integer, intent(in) :: day
      integer :: GetSolarDay
      real(r8) :: Xdm

      if (IsLeapYear(year) .and. month>2) then
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

   function CalcRunningDays(time)
      implicit none
      type(SimTime), intent(in) :: time
      integer :: CalcRunningDays
      integer :: JDN0, JDN1

      call Date2JDN(time%year0, time%month0, time%day0, JDN0)
      call Date2JDN(time%year1, time%month1, time%day1, JDN1)
      CalcRunningDays = JDN1 - JDN0
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
      date = 1d4 * year + 1d2 * month + day
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
   subroutine CalendarConversion(t, year, month, day)
      implicit none
      real(r8), intent(in) :: t
      integer, intent(out) :: year
      integer, intent(out) :: month
      integer, intent(out) :: day
      integer :: nday, JDN
      logical :: isleap

      year = INT(t)
      isleap = IsLeapYear(year)
      if (isleap) then
         nday = NINT(3.66d2*(t-year))
      else
         nday = NINT(3.65d2*(t-year))
      end if
      call Date2JDN(year,1,1,JDN)
      JDN = JDN + nday
      call JDN2Date(JDN, year, month, day)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Get a GMT date and time
   !
   !------------------------------------------------------------------------------
   subroutine GetGMTime(time)
      implicit none
      integer, intent(out) :: time(6)
      integer :: date_time(8), day
      logical :: leap

      call date_and_time(values=date_time)
      leap = IsLeapYear(date_time(1))
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
   ! Purpose: Calculate molecular heat conductivity of some common materials
   !          1) the equations for frozen and unfrozen saturated soils are from
   !             Omar T. Farouki (1981);
   !          2) the equation for phase-changing soils is from Klas Hansson (2004).
   !
   !------------------------------------------------------------------------------
   function CalcLWHeatConductivity(temp, cita, icita)
      implicit none
      real(r8), intent(in) :: temp           ! temperature
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
      real(r8), intent(in) :: surfTemp
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
      GetSpecificLatentHeat4Evap = 1d3 * SLH    ! J/g --> J/kg
      return
   end function

   function CalcLatentHeatWaterAero(waterTemp, RH, wind)
      implicit none
      real(r8), intent(in) :: waterTemp      ! units: K
      real(r8), intent(in) :: RH             ! units: %
      real(r8), intent(in) :: wind           ! units: m/s
      real(r8) :: CalcLatentHeatWaterAero
      real(r8) :: vps, vap, qs, q, Lv
      real(r8) :: adjCe

      adjCe = sa_params(Param_Hscale) * Ce
      vps = CalcSatVP(waterTemp)
      vap = 0.01 * RH * vps
      qs = CalcSpecificHumidity(vps)
      q = CalcSpecificHumidity(vap)
      Lv = GetSpecificLatentHeat4Evap(waterTemp)
      CalcLatentHeatWaterAero = max( Roua*Lv*adjCe*wind*(qs-q), 0d0 )
      return
   end function

   function CalcLatentHeatWaterPM(waterTemp, RH, wind, ps, Rn)
      implicit none
      real(r8), intent(in) :: waterTemp      ! units: K
      real(r8), intent(in) :: RH             ! units: %
      real(r8), intent(in) :: wind           ! units: m/s
      real(r8), intent(in) :: ps             ! units: pascal
      real(r8), intent(in) :: Rn             ! units: W/m2
      real(r8) :: CalcLatentHeatWaterPM      ! units: W/m2
      real(r8) :: delta, vps, vpa, de
      real(r8) :: Ea, gama, Lv

      vps = CalcSatVP(waterTemp)
      vpa = 0.01 * RH * vps
      de = 0.1 * (vps - vpa)           ! vapor pressure deficit (kPa)
      delta = 0.1 * CalcSatVPSlope(waterTemp)  ! units: kPa/K
      Ea = 6.43*(1+0.536*wind)*de      ! bulk aerodynamic expression
      Lv = GetSpecificLatentHeat4Evap(waterTemp)
      gama = Cpa*1d-3*ps/(0.622*Lv)    ! pyschrometric constant (kPa/K)
      CalcLatentHeatWaterPM = (delta*Rn+gama*11.574*Ea)/(delta+gama)
      return
   end function

   function CalcLatentHeatIce(waterTemp, RH, wind)
      implicit none
      real(r8), intent(in) :: waterTemp      ! units: K
      real(r8), intent(in) :: RH             ! units: %
      real(r8), intent(in) :: wind           ! units: m/s
      real(r8) :: CalcLatentHeatIce
      real(r8) :: vap, vps, qs, q
      real(r8) :: adjCe

      adjCe = sa_params(Param_Hscale) * Ce
      vps = CalcSatVP(waterTemp)
      vap = 0.01 * RH * vps
      qs = CalcSpecificHumidity(vps)
      q = CalcSpecificHumidity(vap)
      CalcLatentHeatIce = max( Roua*Ls*adjCe*wind*(qs-q), 0d0 )
      return
   end function

   function CalcSensibleHeat(waterTemp, surfTemp, wind)
      implicit none
      real(r8), intent(in) :: waterTemp      ! units: K
      real(r8), intent(in) :: surfTemp       ! units: K
      real(r8), intent(in) :: wind           ! units: m/s
      real(r8) :: adjCh
      real(r8) :: CalcSensibleHeat

      adjCh = sa_params(Param_Hscale) * Ch
      CalcSensibleHeat = Roua*Cpa*adjCh*wind*(waterTemp-surfTemp)
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
      real(r8) :: zhs, Nsqrt 
      integer :: ii, nn

      if (wind<0.1_r8) then
         ! on no-wind days
         keddy = 0.0_r8
         return
      end if
      nn = size(depth)
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
         zhs = depth(nn) - depth(ii) + DeltaD
         tmp = ubs * exp(-ekman*zhs)
         tmp1 = (40*Nsqrt*(Karman*zhs)**2) / (tmp**2+inft)
         Ri = 0.05 * (-1 + sqrt(1.0+tmp1))
         keddy(ii) = keddy(ii) + (Karman*zhs*tmp/Prandtl)/(1+37*Ri**2)
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
      real(r8) :: CalcHenrySolubility     ! units: umol/(m3*Pa) (M = mole/L)
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
      CalcHenrySolubility = 1.0d+6*CalcHenrySolubility
      return
   end function

   function CalcBunsenSolubility(gas, temp, pH)
      implicit none
      integer, intent(in) :: gas
      real(r8), intent(in) :: temp        ! units: K
      real(r8), intent(in) :: pH          ! units: n/a
      real(r8) :: CalcBunsenSolubility    ! units: m3 / m3
      real(r8) :: henry                   ! units: umol / (m3 * Pa)

      henry = CalcHenrySolubility(gas, temp, pH)
      CalcBunsenSolubility = henry * temp * R * 1.0d-6
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
   subroutine CalcDynamicViscosity(temp, dVsc)
      implicit none
      real(r8), intent(in) :: temp(:)     ! units: K
      real(r8), intent(out) :: dVsc(:)    ! units: kg/s/m

      where (temp<T0-e8)
         dVsc = inf 
      elsewhere
         dVsc = 2.414d-5 * (10.0**(247.8/(temp-140.0)))
      end where
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate buoyant rising speed of bubble from "Bubbles and the 
   !          air-sea exchange of gases in near-saturation conditions" [Woolf, 
   !          D. and Thorpe, S., 1991].
   !
   !------------------------------------------------------------------------------
   function CalcBuoyantVelocity(radius, vv)
      implicit none
      real(r8), intent(in) :: radius      ! units: m
      real(r8), intent(in) :: vv          ! kinematic viscosity (m2/s)
      real(r8) :: CalcBuoyantVelocity     ! units: m/s               
      real(r8) :: rr, xx, par
      real(r8) :: tmp1, tmp2

      if (vv>1.0d+10) then
         CalcBuoyantVelocity = 0.0_r8
         return
      end if
      rr = radius*1.0d+6      ! change meter to micron
      if (rr<80) then
         CalcBuoyantVelocity = G*(radius**2)/(3*vv)
      else if (rr>150) then
         xx = G*(radius**3)/(vv**2)
         CalcBuoyantVelocity = G*(radius**2)/ &
            (vv*18*(1-2/(1+sqrt(1+0.091*xx))))
      else
         xx = G*(radius**3)/(vv**2)
         tmp1 = G*(radius**2)/(3*vv)
         tmp2 = G*(radius**2)/(vv*18*(1-2/(1+sqrt(1+0.091*xx))))
         par = (rr-80.0)/70.0
         CalcBuoyantVelocity = tmp1*(1-par) + tmp2*par
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

      T = temp - T0     ! to celcius
      if (T<0) then
         CalcSchmidtNumber = inf
         return
      else if (T>30) then
         T = 30
      end if
      if (gas==Wn2) then
         CalcSchmidtNumber = 1970.7-131.45*T+4.139*T**2-0.052106*T**3
      else if (gas==Wo2) then
         CalcSchmidtNumber = 1800.6-120.1*T+3.7818*T**2-0.047608*T**3
      else if (gas==Wco2) then
         CalcSchmidtNumber = 1911-113.7*T+2.967*T**2-0.02943*T**3
      else if (gas==Wch4) then
         CalcSchmidtNumber = 1898-110.1*T+2.834*T**2-0.02791*T**3
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
         CalcGasDiffusivityInWater = 2.57d-9*(temp/273.0)
      else if (gas==Wo2) then
         CalcGasDiffusivityInWater = 2.4d-9*(temp/298.0)
      else if (gas==Wco2) then
         CalcGasDiffusivityInWater = 1.81d-6*exp(-2032.6/temp)
      else if (gas==Wch4) then
         CalcGasDiffusivityInWater = 1.5d-9*(temp/298.0)
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
      real(r8), intent(in) :: pH             ! units: n/a
      real(r8), intent(in) :: pressure       ! partial pressure (Pa)
      real(r8) :: CalcEQConc                 ! units: umol/m3
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
      real(r8), intent(in) :: mflow    ! units: umol/m2/s
      real(r8), intent(in) :: area     ! units: m2
      real(r8), intent(in) :: depth    ! units: m
      real(r8), intent(in) :: temp     ! units: K
      real(r8), intent(in) :: pr0      ! units: Pascal
      real(r8) :: CalcBubbleAirflow    ! units: m3/s
      real(r8) :: pr

      pr = pr0 + Roul * G * depth
      CalcBubbleAirflow = 1d-6 * mflow * area * R * temp / pr
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: Construct the lake depth grid and bathymetry
   !
   !------------------------------------------------------------------------------
   subroutine ConstructDepthVector(info, Zw, Zs, dZw, dZs)
      implicit none
      type(LakeInfo), intent(in) :: info     ! lake information object
      real(r8), intent(out) :: Zw(:)         ! water depth vector
      real(r8), intent(out) :: Zs(:)         ! sediment depth vector
      real(r8), intent(out) :: dZw(:)        ! water grid thickness
      real(r8), intent(out) :: dZs(:)        ! sediment grid thickness
      real(r8), parameter :: KArr(20) = (/0.0733, 0.0914, 0.1017, &
            0.1088, 0.1143, 0.1188, 0.1225, 0.1257, 0.1285, &
            0.1310, 0.1333, 0.1354, 0.1373, 0.1390, 0.1407, &
            0.1422, 0.1436, 0.1449, 0.1462, 0.1474/)
      real(r8) :: K1, K2, denom, pp
      integer :: nw, ns, ii, idx

      ! for water column
      nw = size(Zw) - 1
      if (info%depth<=5.0) then
         denom = info%depth / dble(nw)
         do ii = 1, nw+1, 1
            Zw(ii) = (ii-1) * denom
         end do
      else
         if (info%depth<=20) then
            K1 = -1.572d-4*info%depth**2 + 6.861d-3*info%depth - 2.686d-2
         else if (info%depth<=50) then
            K1 = -1.372d-5*info%depth**2 + 1.794d-3*info%depth + 1.757d-2
         else if (info%depth<1000) then
            idx = INT(info%depth/50.0)
            pp = (info%depth - 50*idx) / 50.0
            K1 = (1.0-pp)*KArr(idx) + pp*KArr(idx+1)
         else
            K1 = KArr(20) + (info%depth - 1000) / 50.0 * 1.0d-3
         end if
         denom = info%depth / (exp(K1*nw) - 1)
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
      K2 = 0.03_r8
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

   subroutine BuildLakeBathymetry(info, Zw, Az, dAz)
      implicit none
      type(LakeInfo), intent(in) :: info     ! lake information object
      real(r8), intent(in) :: Zw(:)          ! water depth vector
      real(r8), intent(out) :: Az(:)         ! water depth cross-section
      real(r8), intent(out) :: dAz(:)        ! cross-section difference
      integer :: ii, nz

      ! for cross section
      Az = info%Asurf
      nz = size(Az)
      do ii = 1, nz, 1
         if (ii==1) then
            dAz(ii) = 0.5 * (Az(ii) - Az(ii+1))
         else if (ii==nz) then
            dAz(ii) = 0.5 * (Az(ii) + Az(ii-1))
         else
            dAz(ii) = 0.5 * (Az(ii-1) - Az(ii+1))
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
   subroutine CalcSedBottomTemp(info, t2mref, tsb)
      implicit none
      type(LakeInfo), intent(in) :: info
      real(r8), intent(in) :: t2mref      ! units: K
      real(r8), intent(out) :: tsb        ! units: K
      real(r8) :: tgw, zl, geom

      if (info%margin==1) then
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
   ! Purpose: Calculate active layer thickness based on Nixon and McRoberts
   !          Model (two-layer soils: upper layer 1 and lower layer 2)
   !  "An Analytical Model of the Ground Surface Temperature Under Snowcover
   !  with Soil Freezing" [RISEBOROUGH, 2001]
   !
   !------------------------------------------------------------------------------
   subroutine CalcActiveLayerThickness(winter, temp, snow, k1, k2, L1, L2, &
                                       t, dlt)
      implicit none
      integer, intent(in)  :: winter   ! winter flag
      real(r8), intent(in) :: temp     ! air temperature (K)
      real(r8), intent(in) :: snow     ! sonw thickness (m)
      real(r8), intent(in) :: k1       ! heat conductivity of layer 1 (W/m/K)
      real(r8), intent(in) :: k2       ! heat conductivity of layer 2 (W/m/K)
      real(r8), intent(in) :: L1       ! volumetric latent heat of layer 1 (J/m3)
      real(r8), intent(in) :: L2       ! volumetric latent heat of layer 2 (J/m3)
      real(r8), intent(in) :: t        ! time since freezing/thawing starts (s)
      real(r8), intent(inout) :: dlt   ! dynamic layer thickness (m)
      real(r8), parameter :: H1 = 0.2_r8  ! organic layer thickness (m)
      real(r8) :: Ts, rs, ra, XT

      if (winter==1) then
         if (temp>=T0) then
            return
         end if
         ! calculate temperature at snow-soil interface
         if (snow<0.01 .or. dlt<0.01) then
            Ts = temp
         else
            rs = Kn0 / snow
            ra = k1 / dlt
            Ts = T0 + (temp-T0) * rs / (rs + ra)
         end if
         XT = sqrt( 2 * k1 * (T0-Ts) * t / L1 )
         if (XT>H1) then
            XT = sqrt( (H1*(k2/k1))**2 + 2*(k2/k1)*(L2/L1)*(XT**2-H1**2) ) &
                  - (k2/k1-1)*H1
         end if
      else
         if (temp<=T0) then
            return
         end if
         XT = sqrt( 2 * k1 * (temp-T0) * t / L1 )
         if (XT>H1) then
            XT = sqrt( (H1*(k2/k1))**2 + 2*(k2/k1)*(L2/L1)*(XT**2-H1**2) ) &
                  - (k2/k1-1)*H1
         end if
      end if
      dlt = XT
   end subroutine

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
            8.4878*temp**2 + 0.05607*temp**3)
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
      denom = 9.99868d+2 + 1.0d-3*(65.185*Tw - 8.4878*Tw**2 + 0.05607*Tw**3)
      numer = 1.0d-3 * (65.185 - 16.9756*Tw + 0.16821*Tw**2)
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

end module phy_utilities_mod
