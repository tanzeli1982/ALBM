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
   ! thickness of sediments (m)
   real(r8), parameter :: SED_DEPTH = 25.0
   ! the maximum step length (s)
   integer, parameter :: MAX_OF_STEP = 1800
   ! N2, O2, CO2 and CH4
   integer, parameter :: NGAS = 4
   ! N2, O2, CO2, CH4, TP, aqDOC, trDOC
   integer, parameter :: NWSUB = 7
   ! N2, O2, CO2, CH4, and TP
   integer, parameter :: NSSUB = 5
   ! # of phytoplankton function groups
   integer, parameter :: NPOC = 2
   ! # of dissolved organic carbon groups
   integer, parameter :: NDOC = 2
   ! # of carbon pools
   integer, parameter :: NPOOL = 2
   ! # of lake types
   integer, parameter :: NLAKTYPE = 1
   ! Wavelength number
   integer, parameter :: NSPCTM = 2002
   ! Tolerance of temperature (units: K)
   real(r8), parameter :: Ttol = 1.0d-4
   ! Tolerance of sediment substance concentration (units: umol/m3)
   real(r8), parameter :: SStol(NSSUB) = (/1.d-2, 1.d-2, 1.d-2, &
                                    1.d-4, 1.d-2/)
   ! Tolerance of water substance concentration (units: umol/m3)
   real(r8), parameter :: WStol(NWSUB) = (/1.d-1, 1.d-1, 1.d-2, &
                                    1.d-4, 1.d-2, 1.0d-1, 1.0d-1/)
   ! Tolerance of water particulate concentration (units: umol/m3)
   real(r8), parameter :: WPtol(NPOC) = (/1.d-1, 1.d-1/)
   ! Tolerance of bubble gas concentration (units: umol/m3/mm)
   real(r8), parameter :: Bubtol = 1.0d-6
   ! general group
   integer :: lake_id 
   character(cx) :: lake_file = ""
   character(cx) :: lakeid_file = ""
   character(cx) :: bthmtry_dir = ""
   character(cx) :: param_dir = ""
   character(cx) :: param_file = ""
   character(len=32) :: run_mode = ""
   ! debug group
   logical :: DEBUG
   ! simulation group
   integer :: Start_Year, Start_Month, Start_Day, End_Year
   integer :: End_Month, End_Day
   integer :: Spinup_Month, Spinup_Day, nSpinup
   ! resolution group
   integer :: NWLAYER, NSLAYER, NRLAYER 
   integer :: WATER_LAYER
   ! run data group
   character(len=32) :: forcing_tstep = ""
   character(cx) :: forcing_dir = ""
   character(cx) :: hydro_dir = ""
   character(cx) :: tref_file = ""
   character(cx) :: tas_file = ""
   character(cx) :: tasmax_file = ""
   character(cx) :: tasmin_file = ""
   character(cx) :: hurs_file = ""
   character(cx) :: ps_file = ""
   character(cx) :: pr_file = ""
   character(cx) :: prsn_file = ""
   character(cx) :: rsds_file = ""
   character(cx) :: rlds_file = ""
   character(cx) :: wind_file = ""
   ! archive group
   character(len=32) :: archive_tstep = ""
   character(cx) :: archive_dir = ""
   ! radiation group
   character(cx) :: solar_dir = ""
   character(cx) :: gas_dir = ""
   character(cx) :: albedo_dir = ""
   character(cx) :: co2_file = ""
   character(cx) :: o3_file = ""
   character(cx) :: aod_file = ""
   ! lake information object
   type(LakeInfo) :: lake_info
end module shr_ctrl_mod
