module soil_thermal_mod

!---------------------------------------------------------------------------------
! Purpose: this module calculates the temperature profile in the sediments
!          from "Temperature variability in lake sediments" [Fang et al., 1998].
!          The top boundary is driven by the bottom water temperature.
!          The bottom boundary is assumed no heat transfer.
!          
!---------------------------------------------------------------------------------
   use shr_kind_mod,          only : r8
   use shr_ctrl_mod,          only : WATER_LAYER, NSLAYER, DEBUG, lake_info
   use shr_param_mod
   use shr_typedef_mod,       only : RungeKuttaCache1D
   use phy_utilities_mod
   use data_buffer_mod

   implicit none
   private
   public :: InitializeSedThermalModule, DestructSedThermalModule
   public :: SedThermalModuleSetup, SedThermalModuleCallback
   public :: mem_ts, SedHeatEquation 
   ! memory cache for Runge-Kutta
   type(RungeKuttaCache1D) :: mem_ts
   ! volumetric heat capacity (J/K/m3)
   real(r8), allocatable :: Cvt(:)

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: Private data initialization, get/set and destruction
   !
   !------------------------------------------------------------------------------
   subroutine InitializeSedThermalModule()
      implicit none

      allocate(Cvt(NSLAYER+1))
      allocate(mem_ts%K1(NSLAYER+1))
      allocate(mem_ts%K2(NSLAYER+1))
      allocate(mem_ts%K3(NSLAYER+1))
      allocate(mem_ts%K4(NSLAYER+1))
      allocate(mem_ts%K5(NSLAYER+1))
      allocate(mem_ts%K6(NSLAYER+1))
      allocate(mem_ts%nxt4th(NSLAYER+1))
      allocate(mem_ts%nxt5th(NSLAYER+1))
      allocate(mem_ts%interim(NSLAYER+1))
      allocate(mem_ts%rerr(NSLAYER+1))

      call InitializeSoilStateVariables() 
      Cvt = 0.0_r8
   end subroutine

   subroutine DestructSedThermalModule()
      implicit none

      deallocate(Cvt)
      deallocate(mem_ts%K1)
      deallocate(mem_ts%K2)
      deallocate(mem_ts%K3)
      deallocate(mem_ts%K4)
      deallocate(mem_ts%K5)
      deallocate(mem_ts%K6)
      deallocate(mem_ts%nxt4th)
      deallocate(mem_ts%nxt5th)
      deallocate(mem_ts%interim)
      deallocate(mem_ts%rerr)
   end subroutine

   subroutine SedThermalModuleSetup()
      implicit none
      real(r8) :: Rous, Cps, Porosity
      real(r8) :: Kss, satLW
      integer :: ii

      Rous = sa_params(Param_Rous)
      Cps = sa_params(Param_Cps)
      Kss = sa_params(Param_Ks)
      Porosity = sa_params(Param_Por)
      do ii = 1, NSLAYER+1, 1
         satLW = CalcSoilWaterSaturation(m_sedTemp(ii))
         m_Ks(ii) = CalcSedHeatConductivity(Porosity, satLW, Kss)
      end do
      Cvt = Cps*Rous*(1.0-Porosity) + Cpl*Roul*(Porosity-m_sedIce) &
            + Cpi*Roui*m_sedIce
      m_Ks = m_Ks / Cvt
   end subroutine

   subroutine SedThermalModuleCallback()
      implicit none

      call AdjustIceTempForSediment()
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Discretize Heat differential equation by Finite Element Method
   !          This sub-routine returns time derivatives of lake sediment temperature
   ! 1. The first point is controlled by bottom lake water and it is excluded from
   !    sediment heat calculation but just used to driven top heat flux
   !
   !------------------------------------------------------------------------------
   subroutine SedHeatEquation(temp, dtemp)
      implicit none
      real(r8), intent(in) :: temp(NSLAYER+1)
      real(r8), intent(out) :: dtemp(NSLAYER+1)
      real(r8) :: qt, qb, dT1, dT2
      real(r8) :: dT, aa, bb, rad
      integer  :: ii

      dT = temp(1) - m_waterTemp(WATER_LAYER+1)
      qt = m_Kbm * dT / DeltaD   ! top boundary condition
      qb = 0.0_r8       ! bottom boundary condition
      do ii = 1, NSLAYER+1, 1
         if(ii==1) then
            ! suppress the rad item because it causes unexpected warming
            ! rad = m_Iab(WATER_LAYER+2) / m_dZs(ii) / Cvt(ii)
            rad = 0.0_r8
            aa = 0.5*(m_Ks(ii) + m_Ks(ii+1))
            dT1 = (temp(ii+1) - temp(ii)) / (m_Zs(ii+1) - m_Zs(ii))
            dtemp(ii) = (aa*dT1 - qt) / m_dZs(ii) + rad 
         else if(ii==NSLAYER+1) then
            bb = 0.5*(m_Ks(ii-1) + m_Ks(ii))
            dT2 = (temp(ii) - temp(ii-1)) / (m_Zs(ii) - m_Zs(ii-1))
            dtemp(ii) = (qb - bb*dT2) / m_dZs(ii)
         else
            aa = 0.5*(m_Ks(ii+1) + m_Ks(ii))
            bb = 0.5*(m_Ks(ii-1) + m_Ks(ii))
            dT1 = (temp(ii+1) - temp(ii)) / (m_Zs(ii+1) - m_Zs(ii))
            dT2 = (temp(ii) - temp(ii-1)) / (m_Zs(ii) - m_Zs(ii-1))
            dtemp(ii) = (aa*dT1 - bb*dT2) / m_dZs(ii)
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Adjust temperature and ice percentage of every layer
   !
   !------------------------------------------------------------------------------
   subroutine AdjustIceTempForSediment()
      implicit none
      real(r8) :: tmp1, tmp2, Porosity, Cps
      integer :: ii

      Porosity = sa_params(Param_Por)
      Cps = sa_params(Param_Cps)
      do ii = 1, NSLAYER+1, 1             ! Sediment freezing
         if (m_sedIce(ii)<Porosity-e8 .and. m_sedTemp(ii)<T0-e8) then
            tmp1 = Cpi*(T0-m_sedTemp(ii))*m_sedIce(ii) +       &
            Cpl*(T0-m_sedTemp(ii))*(Porosity-m_sedIce(ii)) +   &
            Cps*(T0-m_sedTemp(ii))*(1.0-Porosity)
            tmp2 = Lf*(Porosity-m_sedIce(ii))
            if (tmp2>tmp1) then
               m_sedIce(ii) = m_sedIce(ii) + tmp1/Lf
               m_sedTemp(ii) = T0
            else
               m_sedIce(ii) = Porosity
               m_sedTemp(ii) = T0 - (tmp1-tmp2)/(Cpi*Porosity+Cps*(1-Porosity))
            end if
         end if
      end do                           
      do ii = 1, NSLAYER+1, 1             ! Sediment thawing
         if (m_sedIce(ii)>e8 .and. m_sedTemp(ii)>T0+e8) then
            tmp1 = Cpi*(m_sedTemp(ii)-T0)*m_sedIce(ii) +       &
            Cpl*(m_sedTemp(ii)-T0)*(Porosity-m_sedIce(ii)) +   &
            Cps*(m_sedTemp(ii)-T0)*(1.0-Porosity)
            tmp2 = Lf*m_sedIce(ii)
            if (tmp2>tmp1) then
               m_sedIce(ii) = m_sedIce(ii) - tmp1/Lf
               m_sedTemp(ii) = T0
            else
               m_sedIce(ii) = 0.0_r8
               m_sedTemp(ii) = T0 + (tmp1-tmp2)/(Cpl*Porosity+Cps*(1-Porosity))
            end if
         end if
      end do
      call UpdateSedWaterBoundIndex()
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update the top and bottom index of sediment free water layer
   !
   !------------------------------------------------------------------------------
   subroutine UpdateSedWaterBoundIndex()
      implicit none
      logical, save :: sfirstime = .True.
      integer :: ii

      m_sedWaterTopIndex = NSLAYER + 2
      do ii = 1, NSLAYER+1, 1
         if (m_sedIce(ii)<e8) then
            m_sedWaterTopIndex = ii
            exit
         end if
      end do
      m_sedWaterBtmIndex = 0
      do ii = NSLAYER+1, 1, -1
         if (m_sedIce(ii)<e8) then
            m_sedWaterBtmIndex = ii
            exit
         end if
      end do
      if (m_sedWaterTopIndex>m_sedWaterBtmIndex) then
         if (sfirstime) then
            print "(A, I0, A)", "Lake ", lake_info%id, &
                  ": lake sediment is totally frozen"
            sfirstime = .False.
         end if
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Initialize soil temperature and ice profiles
   !
   !------------------------------------------------------------------------------
   subroutine InitializeSoilStateVariables()
      implicit none
      real(r8) :: MAGT, excess_ice
      real(r8) :: talik, dist, pMAGT, Porosity, Kss
      real(r8) :: Tbot, Ttop, Ttk, Tuf, Tlf
      real(r8) :: Ks_tk, Ks_uf, Ks_lf, dTtk, dTuf
      real(r8) :: satLW_tk, satLW_uf, satLW_lf
      integer  :: ii

      call CalcSedBottomTemp(lake_info, m_radPars%tref, MAGT)

      ! calculate talik depth
      if (lake_info%thrmkst==2 .and. lake_info%excice>e8) then
         talik = (lake_info%depth-2.0) / lake_info%excice 
      else if (lake_info%thrmkst==1 .and. lake_info%excice>e8) then
         talik = (lake_info%depth-1.0) / lake_info%excice
      else
         talik = lake_info%hsed 
      end if
      talik = max(0._r8, min(talik, lake_info%hsed))

      ! construct thermal profile
      Porosity = sa_params(Param_Por)
      Kss = sa_params(Param_Ks)
      if (talik>lake_info%hsed-e8) then
         Tbot = max(MAGT, T0)
      else
         if (MAGT<T0) then
            Tbot = MAGT
         else
            Tbot = T0 - 1.0
         end if
      end if
      if (DEBUG) then
         print "(A, I0, A, F10.4)", "Lake ", lake_info%id, &
               ": mean ground temperature ", Tbot-T0
      end if
      if (talik>e8) then
         Ttop = T0 + 4.6 
         Ttk = 0.5 * (Ttop + T0)
         Tuf = 0.5 * (T0 + Tbot)
         Tlf = Tbot
         satLW_tk = CalcSoilWaterSaturation(Ttk)
         satLW_uf = CalcSoilWaterSaturation(Tuf)
         satLW_lf = CalcSoilWaterSaturation(Tlf)
         Ks_tk = CalcSedHeatConductivity(Porosity, satLW_tk, Kss)
         Ks_uf = CalcSedHeatConductivity(Porosity, satLW_uf, Kss)
         Ks_lf = CalcSedHeatConductivity(Porosity, satLW_lf, Kss)
         dist = Ks_uf / Ks_tk * talik * (T0 - Tbot) / (Ttop - T0)
         if (dist>0) then
            pMAGT = talik + dist
            dTtk = (T0 - Ttop) / talik
            dTuf = (Tbot - T0) / dist
            do ii = 1, NSLAYER+1, 1
               if (m_Zs(ii)<=talik+m_Zs(1)) then
                  m_sedTemp(ii) = Ttop + dTtk * (m_Zs(ii)-m_Zs(1))
               else if (m_Zs(ii)<=pMAGT+m_Zs(1)) then
                  m_sedTemp(ii) = T0  + dTuf * (m_Zs(ii)-talik)
               else
                  m_sedTemp(ii) = Tbot
               end if
            end do
         else
            dTtk = (Tbot - Ttop) / talik
            do ii = 1, NSLAYER+1, 1
               m_sedTemp(ii) = Ttop + dTtk * (m_Zs(ii)-m_Zs(1))
            end do
         end if
      else
         Ttop = T0
         dTtk = (Tbot - Ttop) / lake_info%hsed
         do ii = 1, NSLAYER+1, 1
            m_sedTemp(ii) = Ttop + dTtk * (m_Zs(ii)-m_Zs(1))
         end do
      end if

      ! construct ice profile
      where (m_sedTemp<T0) 
         m_sedIce = Porosity
      elsewhere
         m_sedIce = 0.0_r8
      end where
      m_Ks = 0.0_r8
      call UpdateSedWaterBoundIndex()
   end subroutine

end module soil_thermal_mod
