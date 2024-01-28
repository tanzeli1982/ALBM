module shr_ctrl_mod
!----------------------------------------------------------------------------
! shared precision control constants
!----------------------------------------------------------------------------
   use shr_kind_mod,       only : r8, cx => SHR_KIND_CX, cs => SHR_KIND_CS
   use shr_typedef_mod,    only : LakeInfo

   public
   ! Real infinitestimal valu
   real(r8), parameter :: SHR_CTRL_E8 = 1.d-8
   ! Infinite numble
   real(r8), parameter :: INFINITE_E8 = 1.d+30
   ! Infinitesimal number
   real(r8), parameter :: INFINITESIMAL_E8 = 1.d-30
   ! Relative tolerance of Runge-Kutta method
   real(r8), parameter :: TOL_E8 = 1.d-6
   ! sub-hourly update time step (s)
   real(r8), parameter :: TSTEP_SUB = 300.0_r8
   ! the maximum step length (s)
   integer, parameter :: MAX_OF_STEP = 1800
   ! N2, O2, CO2 and CH4
   integer, parameter :: NGAS = 4
   ! N2, O2, CO2, CH4, TP
   integer, parameter :: NWSUB = 5
   ! N2, O2, CO2, CH4, TP
   integer, parameter :: NSSUB = 5
   ! # of phytoplankton function groups
   integer, parameter :: NPOC = 2
   ! # of dissolved organic carbon groups
   integer, parameter :: NDOC = 2
   ! # of carbon pools
   integer, parameter :: NPOOL = 3
   ! # of lake types
   integer, parameter :: NLAKTYPE = 1
   ! Wavelength number
   integer, parameter :: NSPCTM = 2002
   ! Tolerance of temperature (units: K)
   real(r8), parameter :: Ttol = 1.0d-4
   ! Tolerance of soil temperature (units: K)
   real(r8), parameter :: STtol = 1.0d-4 
   ! Tolerance of sediment substance concentration (units: mol/m3)
   real(r8), parameter :: SStol(NSSUB) = (/1.d-8, 1.d-8, 1.d-8, &
                                    1.d-8, 1.d-8/)
   ! Tolerance of water substance concentration (units: mol/m3)
   real(r8), parameter :: WStol(NWSUB) = (/1.d-7, 1.d-7, 1.d-8, &
                                    1.d-8, 1.d-8/)
   ! general group
   integer :: lake_range(2)
   character(cx) :: lake_file
   character(cx) :: opt_file
   character(cx) :: bathy_file
   character(len=32) :: run_mode    ! "regular" or "sensitivity"
   ! debug group
   logical :: DEBUG
   ! simulation group
   logical :: Thermal_Module, Bubble_Module, Diagenesis_Module
   logical :: Carbon_Module
   logical :: Use_Leap
   integer :: Start_Year, Start_Month, Start_Day, End_Year
   integer :: End_Month, End_Day
   integer :: Spinup_Month, Spinup_Day, nSpinup
   ! resolution group
   integer :: NWLAYER, NSLAYER, NRLAYER, NSCOL 
   integer :: WATER_LAYER, SED_LAYER
   ! run data group
   character(len=32) :: forcing_tstep
   character(len=32) :: forcing_time   ! "UTC" or "LT"
   integer :: forcing_nhour
   character(cx) :: tas_file, tasmax_file, tasmin_file
   character(cx) :: hurs_file
   character(cx) :: ps_file
   character(cx) :: pr_file, prsn_file
   character(cx) :: rsds_file, rlds_file
   character(cx) :: wind_file
   character(cx) :: tref_file
   character(cx) :: soc_file
   character(cx) :: wlnd_file
   character(cx) :: veg_file
   character(cx) :: hydro_file
   ! archive group
   character(len=32) :: archive_tstep  ! "day" or "hour"
   character(cx) :: archive_dir
   ! radiation group
   character(cx) :: solar_dir, gas_dir, albedo_dir
   character(cx) :: co2_file, o3_file, aod_file
   ! lake information object
   type(LakeInfo) :: lake_info
   ! master processor flag
   integer :: masterproc
end module shr_ctrl_mod
