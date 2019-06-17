module phy_const_mod
!----------------------------------------------------------------------------
!  Purpose: physical constants
!  Hydrodynamics parameters are from Wuest and Lorke (2003),
!  Goudsmit et al. (2002) and etc.
!----------------------------------------------------------------------------
   use shr_kind_mod,          only: r8
   public
   ! substance identifier
   integer, parameter :: Wn2 = 1, Wo2 = 2, Wco2 = 3, Wch4 = 4, Wsrp = 5, &
                         Waqdoc = 6, Wtrdoc = 7
   ! number of days in months (non-leap year)
   integer, parameter :: DAY_IN_MONTH(12) = (/31, 28, 31, 30, 31, 30, 31, &
                        31, 30, 31, 30, 31/)
   ! seconds in a day and year
   real(r8), parameter :: SECOND_OF_YEAR = 3.1536d+7
   real(r8), parameter :: SECOND_OF_DAY = 8.64d+4
   ! Avogadro constant
   real(r8), parameter :: NA = 6.0221412927d+023
   ! gravitational acceleration (m/s2)
   real(r8), parameter :: G = 9.8
   ! density of water, ice and snow (kg/m3)
   real(r8), parameter :: Roul = 1000.0, Roui = 931.0
   ! density of gray ice (kg/m3)
   real(r8), parameter :: Roue = 890.0
   ! density of air (kg/m3)
   real(r8), parameter :: Roua = 1.225
   ! Heat capacity of water, ice and snow at constant pressure (J/(K*kg))
   real(r8), parameter :: Cpl = 4.2d+3, Cpi = 2.1d+3, Cpn = 2.1d+3
   ! Heat capacity of dry air at constant pressure (J/(K*kg))
   real(r8), parameter :: Cpa = 1.006d+3
   ! latent heat of freezing and sublimation (J/kg)
   real(r8), parameter :: Lf = 3.337d+5, Ls = 2.834d+6
   ! latent heat of peat and mineral soil fusion (J/m3)
   real(r8), parameter :: Lpet = 2.672d+8, Lmnr = 1.002d+8
   ! Kelvin degree of 0 Celsius (K)
   real(r8), parameter :: T0 = 273.15
   ! Standard atmosphere pressure (Pa)
   real(r8), parameter :: P0 = 98150.0 
   ! Thermal conductivity of ice, water and snow (W/(m*K))
   real(r8), parameter :: Ki0 = 2.2156, Kw0 = 0.558, Kn0 = 0.27
   real(r8), parameter :: Ke0 = 2.0
   ! Thermal conductivity of air (W/(m*K))
   real(r8), parameter :: Ka0 = 0.0234
   ! Thermal conductivity of peat and mineral soils (W/m/K)
   real(r8), parameter :: Kpet = 0.25, Kmnr = 1.3
   ! Stefan-Boltzmann constant (W/(m2*K4))
   real(r8), parameter :: Stefan = 5.6704d-8
   ! Ideal-gas constant (J/(mol*K))
   real(r8), parameter :: R = 8.314462175d+0
   ! Solar constant (W/m2)
   real(r8), parameter :: Solar = 1366.1
   ! Pi constant
   real(r8), parameter :: Pi = 3.14159265d+0
   ! Euler's constant
   real(r8), parameter :: e = 2.71828183d+0
   ! Transfer coefficeint for sensible heat and latent heat
   real(r8), parameter :: Ch = 1.75d-3, Ce = 1.75d-3
   ! Emissivity of water, ice and snow for long-wave radiation
   real(r8), parameter :: Epsw = 0.97, Epsi = 0.97, Epsn = 0.94
   real(r8), parameter :: Epse = 0.94
   ! Absorption extinction coefficient (m-1) of snow and gray ice
   real(r8), parameter :: EtanVIS = 6.0, EtanIR = 40.0
   real(r8), parameter :: EtaeVIS = 3.75, EtaeIR = 40.0
   ! Refractive index of water, ice and gray ice
   real(r8), parameter :: Rfrw = 1.333, Rfri = 1.309, Rfre = 1.309
   ! Surface albedo of gray ice
   real(r8), parameter :: Alphae = 0.40
   ! N2, O2, CO2 and CH4 mixing ratio in the atmosphere
   real(r8), parameter :: Xn2 = 0.78, Xo2 = 0.21
   real(r8), parameter :: Xco2 = 380.0d-6, Xch4 = 1800.0d-9
   ! Mole mass of N2, O2, CO2 and CH4 (g/mol)
   real(r8), parameter :: MasN2 = 28, MasO2 = 32, MasCO2 = 44, MasCH4 = 16
   ! Mole mass of carbon and phosphorus
   real(r8), parameter :: MasC = 12, MasP = 31
   ! O2 solubility from Engineering toolbox (mg/L)
   real(r8), parameter :: SOLO2(11) = (/14.6, 12.8, 11.3, 10.1, 9.1, 8.3, &
                           7.6, 7.0, 6.5, 6.0, 5.6/)
   ! relative concentration at which ebullition begins
   real(r8), parameter :: Ae = 0.4
   ! bubble gas release rate after thawing or ice breakage (units: s-1)
   real(r8), parameter :: Blr = 7.7d-7
   ! Hydrodynamics parameters are from Wuest and Lorke (2003),
   ! Goudsmit et al. (2002) and etc
   ! diffusive boundary layer thickness (m) 
   real(r8), parameter :: DeltaD = 1.0d-3
   ! bottom drag coefficient at 1-m height
   real(r8), parameter :: Cb1m = 2.0d-3
   ! top drag coefficient at 10-m height
   real(r8), parameter :: Cd10 = 1.3d-3
   ! von Karman's constant
   real(r8), parameter :: Karman = 0.41
   ! turbulent Prandtl number in neutral conditions
   real(r8), parameter :: Prandtl = 1.0
   ! PAR radiation conversion factor from W/m2 to umol/m2/s
   real(r8), parameter :: fconvPAR = 4.6
end module phy_const_mod
