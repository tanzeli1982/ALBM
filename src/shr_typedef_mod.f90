module shr_typedef_mod
!----------------------------------------------------------------------------
! user self-defined types
!----------------------------------------------------------------------------
   use shr_kind_mod,          only: r8, r4
   public
   type SurfaceData
      real(r8) :: temp                 ! 2-m air temperature (K)
      real(r8) :: RH                   ! 2-m relative humidity (%)
      real(r8) :: pressure             ! surface air pressure (Pa)
      real(r8) :: wind                 ! 10-m air wind speed (m/s)
      real(r8) :: rainfall             ! rainfall (m/s)
      real(r8) :: snowfall             ! snowfall (m/s)
      real(r8) :: sw                   ! observed incoming short-wave radiation (W/m2)
      real(r8) :: sw_sim               ! simulated incoming short-wave radiation (W/m2)
      real(r8) :: srd                  ! simulated net SW irradiance (W/m2)
      real(r8) :: lw                   ! incoming long-wave radiation (W/m2)
      real(r8) :: Qsi                  ! inlet discharge (m3/s)
      real(r8) :: tQsi                 ! inlet water temperature (K)
      real(r8) :: dQsi                 ! inlet water density (kg/m3)
      real(r8) :: DICQsi               ! inlet water DIC (umol/m3)
      real(r8) :: DOCQsi               ! inlet water DOC (umol/m3)
      real(r8) :: POCQsi               ! inlet water POC (umol/m3)
      real(r8) :: SRPQsi               ! inlet water SRP (umol/m3)
      real(r8) :: DOQsi                ! inlet water O2 (umol/m3)
      real(r8) :: Qso                  ! water outflow (m3/s)
      real(r8) :: Qgw                  ! groundwater flow (m3/s)
   end type
   type LakeInfo
      integer  :: id                   ! lake id
      integer  :: itype                ! lake type identifier
      real(r8) :: latitude             ! lake latitude
      real(r8) :: longitude            ! lake longitude
      real(r8) :: depth                ! lake depth (m)
      real(r8) :: Asurf                ! lake surface area (m2)
      real(r8) :: zalt                 ! lake altitude (km)
      real(r8) :: kext                 ! chla and CDOM light attenuation (m-1)
      real(r8) :: excice               ! land excessive ice fraction
      real(r8) :: hsed                 ! sediment column thickness (m)
      integer  :: thrmkst              ! 1=thermokarst, 0=vice
      integer  :: margin               ! 1=margin, 2=slope, 0=center
      integer  :: hydroconn            ! 1=hydro connect, 0=vice
   end type
   type SimTime
      integer :: year0                 ! the starting year
      integer :: month0                ! the starting month
      integer :: day0                  ! the starting day
      integer :: year1                 ! the ending year
      integer :: month1                ! the ending month
      integer :: day1                  ! the ending day
   end type
   type RadParaData
      real(r8) :: spr                  ! surface air pressure (mb)
      real(r8) :: tair                 ! air temperature (K)
      real(r8) :: tref                 ! average air temperature (K)
      real(r8) :: RH                   ! relative humidity (%)
      real(r8) :: tau550               ! aerosol optical depth at 550 nm
      real(r8) :: AbO3                 ! ozone total-column density (1000 DU)
      real(r8) :: qCO2                 ! CO2 concentration (ppm)
      real(r8) :: Latit                ! site latitude (+/-)
      real(r8) :: Longit               ! site longitude (+/-)
      integer  :: year, month, day     ! local time date
      integer  :: season               ! 0: winter or fall; 1: summer or spring
   end type
   !-------------------------------------------------------------------------
   !  memory caches for Runge-Kutta methods
   !     K1-K6: Runge-Kutta parameters for 4th Runge-Kutta Method
   !     nxt4th: the guess value of 4th iteration
   !     nxt5th: the guess value of 5th iteration
   !     interim: interimediate guess
   !     rerr: relative error between nxt4th and nxt5th
   !-------------------------------------------------------------------------
   type RungeKuttaCache1D
      real(r8), pointer :: K1(:)
      real(r8), pointer :: K2(:)
      real(r8), pointer :: K3(:)
      real(r8), pointer :: K4(:)
      real(r8), pointer :: K5(:)
      real(r8), pointer :: K6(:)
      real(r8), pointer :: nxt4th(:)
      real(r8), pointer :: nxt5th(:)
      real(r8), pointer :: interim(:)
      real(r8), pointer :: rerr(:)
   end type
   type RungeKuttaCache2D
      real(r8), pointer :: K1(:,:)
      real(r8), pointer :: K2(:,:)
      real(r8), pointer :: K3(:,:)
      real(r8), pointer :: K4(:,:)
      real(r8), pointer :: K5(:,:)
      real(r8), pointer :: K6(:,:)
      real(r8), pointer :: nxt4th(:,:)
      real(r8), pointer :: nxt5th(:,:)
      real(r8), pointer :: interim(:,:)
      real(r8), pointer :: rerr(:,:)
   end type
end module shr_typedef_mod
