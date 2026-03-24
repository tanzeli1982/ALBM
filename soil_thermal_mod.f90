module soil_thermal_mod

!---------------------------------------------------------------------------------
! Purpose: this module calculates the temperature profile in the sediments
!          from "Temperature variability in lake sediments" [Fang et al., 1998].
!          The top boundary is driven by the bottom water temperature.
!          The bottom boundary is assumed no heat transfer.
!          
!---------------------------------------------------------------------------------
   use shr_kind_mod,          only : r8
   use shr_ctrl_mod,          only : WATER_LAYER, SED_LAYER, DEBUG, lake_info
   use shr_param_mod
   use shr_typedef_mod,       only : RungeKuttaCache2D
   use phy_utilities_mod
   use data_buffer_mod

   implicit none
   private
   public :: InitializeSedThermalModule, DestructSedThermalModule
   public :: SedThermalModuleSetup, SedThermalModuleCallback
   public :: mem_ts, SedHeatEquation 
   ! memory cache for Runge-Kutta
   type(RungeKuttaCache2D) :: mem_ts
   ! volumetric heat capacity (J/K/m3)
   real(r8), allocatable :: Cvt(:,:)
   ! radiation-induced temperature change (K/m3/s)
   real(r8), allocatable :: rad(:)

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: Private data initialization, get/set and destruction
   !
   !------------------------------------------------------------------------------
   subroutine InitializeSedThermalModule()
      implicit none

      allocate(Cvt(NSCOL,SED_LAYER+1))
      allocate(rad(NSCOL))
      allocate(mem_ts%K1(NSCOL,SED_LAYER+1))
      allocate(mem_ts%K2(NSCOL,SED_LAYER+1))
      allocate(mem_ts%K3(NSCOL,SED_LAYER+1))
      allocate(mem_ts%K4(NSCOL,SED_LAYER+1))
      allocate(mem_ts%K5(NSCOL,SED_LAYER+1))
      allocate(mem_ts%K6(NSCOL,SED_LAYER+1))
      allocate(mem_ts%nxt4th(NSCOL,SED_LAYER+1))
      allocate(mem_ts%nxt5th(NSCOL,SED_LAYER+1))
      allocate(mem_ts%interim(NSCOL,SED_LAYER+1))
      allocate(mem_ts%rerr(NSCOL,SED_LAYER+1))

      call InitializeSoilStateVariables() 
      Cvt = 0.0_r8
      rad = 0.0_r8
   end subroutine

   subroutine DestructSedThermalModule()
      implicit none

      deallocate(Cvt)
      deallocate(rad)
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
      real(r8) :: Porosity, satLW, sumAz
      integer :: ii, jj, icol
      integer :: indx, indx0

      do icol = 1, NSCOL, 1
         indx0 = COUNT(m_soilColInd<=icol-1) + 1
         indx = COUNT(m_soilColInd<=icol)
         ! inactive sediment column
         if (indx0>indx) then 
            m_Ks(icol,:) = Ks0 / Cps / Rous 
            rad(icol) = 0._r8
            cycle
         end if
         sumAz = sum(m_dAz(indx0:indx))
         do ii = 1, SED_LAYER+1, 1
            Porosity = m_sedpor(icol,ii)
            satLW = CalcSoilWaterSaturation(m_sedTemp(icol,ii))
            m_Ks(icol,ii) = CalcSedHeatConductivity(Porosity, satLW, Ks0)
            Cvt(icol,ii) = Cps*Rous*(1.0-Porosity) + Cpl*Roul* &
                  (Porosity-m_sedIce(icol,ii)) + Cpi*Roui*m_sedIce(icol,ii)
            m_Ks(icol,ii) = m_Ks(icol,ii) / Cvt(icol,ii)
         end do
         if (icol<NSCOL) then
            rad(icol) = 0._r8
            do jj = indx0, indx, 1
               rad(icol) = rad(icol) + (1.0 - Alphas) * (1.0 - &
                  m_fcovBedVeg(jj)) * m_Idw(jj+1) * m_dAz(jj) / sumAz
            end do
            rad(icol) = rad(icol) / m_dZs(1) / Cvt(icol,1)
         else
            rad(icol) = 0._r8
         end if
      end do
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
      real(r8), intent(in) :: temp(NSCOL,SED_LAYER+1)
      real(r8), intent(out) :: dtemp(NSCOL,SED_LAYER+1)
      real(r8) :: qt, qb, dT1, dT2
      real(r8) :: dT, aa, bb, sumAz
      integer  :: ii, jj, icol, top
      integer  :: indx, indx0

      top = m_lakeWaterTopIndex
      qb = 0.0_r8       ! bottom boundary condition
      do icol = 1, NSCOL, 1
         indx0 = COUNT(m_soilColInd<=icol-1) + 1
         indx = COUNT(m_soilColInd<=icol)
         if (indx0>indx) then
            dtemp(icol,:) = 0._r8
            cycle
         end if
         sumAz = sum(m_dAz(indx0:indx))

         ! top boundary condition
         if (icol<NSCOL) then
            qt = 0._r8
            do jj = indx0, indx, 1
               dT = temp(icol,1) - m_waterTemp(jj)
               qt = qt + m_Ktb(jj) * dT * m_dAz(jj) / sumAz / DeltaD
            end do
         else
            dT = temp(icol,1) - m_waterTemp(indx)
            qt = m_Ktb(indx) * dT / DeltaD
         end if
         do ii = 1, SED_LAYER+1, 1
            if(ii==1) then
               aa = 0.5*(m_Ks(icol,ii) + m_Ks(icol,ii+1))
               dT1 = (temp(icol,ii+1) - temp(icol,ii)) / (m_Zs(ii+1) - m_Zs(ii))
               dtemp(icol,ii) = (aa*dT1 - qt) / m_dZs(ii) + rad(icol) 
            else if(ii==SED_LAYER+1) then
               bb = 0.5*(m_Ks(icol,ii-1) + m_Ks(icol,ii))
               dT2 = (temp(icol,ii) - temp(icol,ii-1)) / (m_Zs(ii) - m_Zs(ii-1))
               dtemp(icol,ii) = (qb - bb*dT2) / m_dZs(ii)
            else
               aa = 0.5*(m_Ks(icol,ii+1) + m_Ks(icol,ii))
               bb = 0.5*(m_Ks(icol,ii-1) + m_Ks(icol,ii))
               dT1 = (temp(icol,ii+1) - temp(icol,ii)) / (m_Zs(ii+1) - m_Zs(ii))
               dT2 = (temp(icol,ii) - temp(icol,ii-1)) / (m_Zs(ii) - m_Zs(ii-1))
               dtemp(icol,ii) = (aa*dT1 - bb*dT2) / m_dZs(ii)
            end if
         end do
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Adjust temperature and ice percentage of every layer
   !
   !------------------------------------------------------------------------------
   subroutine AdjustIceTempForSediment()
      implicit none
      real(r8) :: tmp1, tmp2, Porosity
      integer :: ii, icol

      do icol = 1, NSCOL, 1
         ! inactive sediment column
         if (COUNT(m_soilColInd==icol)==0) then
            cycle
         end if

         do ii = 1, SED_LAYER+1, 1             ! Sediment freezing
            Porosity = m_sedpor(icol,ii)
            if (m_sedIce(icol,ii)<Porosity-e8 .and. m_sedTemp(icol,ii)<T0-e8) then
               tmp1 = Cpi*(T0-m_sedTemp(icol,ii))*m_sedIce(icol,ii) + &
                     Cpl*(T0-m_sedTemp(icol,ii))*(Porosity-m_sedIce(icol,ii)) + &
                     Cps*(T0-m_sedTemp(icol,ii))*(1.0-Porosity)
               tmp2 = Lf*(Porosity-m_sedIce(icol,ii))
               if (tmp2>tmp1) then
                  m_sedIce(icol,ii) = m_sedIce(icol,ii) + tmp1/Lf
                  m_sedTemp(icol,ii) = T0
               else
                  m_sedIce(icol,ii) = Porosity
                  m_sedTemp(icol,ii) = T0 - (tmp1-tmp2)/(Cpi*Porosity+ &
                        Cps*(1-Porosity))
               end if
            end if
         end do                           
         do ii = 1, SED_LAYER+1, 1             ! Sediment thawing
            if (m_sedIce(icol,ii)>e8 .and. m_sedTemp(icol,ii)>T0+e8) then
               tmp1 = Cpi*(m_sedTemp(icol,ii)-T0)*m_sedIce(icol,ii) + &
                     Cpl*(m_sedTemp(icol,ii)-T0)*(Porosity-m_sedIce(icol,ii)) + &
                     Cps*(m_sedTemp(icol,ii)-T0)*(1.0-Porosity)
               tmp2 = Lf*m_sedIce(icol,ii)
               if (tmp2>tmp1) then
                  m_sedIce(icol,ii) = m_sedIce(icol,ii) - tmp1/Lf
                  m_sedTemp(icol,ii) = T0
               else
                  m_sedIce(icol,ii) = 0.0_r8
                  m_sedTemp(icol,ii) = T0 + (tmp1-tmp2)/(Cpl*Porosity+ &
                        Cps*(1-Porosity))
               end if
            end if
         end do
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
      integer :: ii, icol

      do icol = 1, NSCOL, 1 
         if (COUNT(m_soilColInd==icol)==0) then
            m_sedWaterTopIndex(icol) = 1
            m_sedWaterBtmIndex(icol) = SED_LAYER + 1
            cycle
         end if

         m_sedWaterTopIndex(icol) = SED_LAYER + 2
         do ii = 1, SED_LAYER+1, 1
            if (m_sedIce(icol,ii)<e8) then
               m_sedWaterTopIndex(icol) = ii
               exit
            end if
         end do
         m_sedWaterBtmIndex(icol) = 0
         do ii = SED_LAYER+1, 1, -1
            if (m_sedIce(icol,ii)<e8) then
               m_sedWaterBtmIndex(icol) = ii
               exit
            end if
         end do
      end do
      if (ANY(m_sedWaterTopIndex>m_sedWaterBtmIndex)) then
         if (sfirstime) then
            print "(A, I0, A, I0, A)", "Lake ", lake_info%id, & 
                  ": at least one sediment column is totally frozen"
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
      real(r8) :: MAGT, talik, talik_theory
      real(r8) :: dist, hMAGT
      real(r8) :: Tbot, Ttop, Ttk, Tfz
      real(r8) :: Ks_tk, Ks_fz, dTtk, dTfz
      real(r8) :: satLW_tk, satLW_fz
      integer  :: ii, icol, indx
      logical  :: margin

      do icol = 1, NSCOL, 1
         ! inactive sediment column
         if (COUNT(m_soilColInd==icol)==0) then
            m_sedTemp(icol,:) = T0 + 4.0_r8
            m_sedIce(icol,:) = 0._r8
            cycle
         end if

         indx = COUNT(m_soilColInd<=icol)
         ! calculate sediment bottom temperature
         margin = (icol==1)
         call CalcSedBottomTemp(lake_info, m_radPars%tref, margin, MAGT)
         ! calculate talik depth
         if (lake_info%thrmkst>0) then
            talik = m_Zw(indx) / lake_info%excice
         else
            talik = lake_info%hsed
         end if
         talik = min(talik, lake_info%hsed)

         ! construct temperature profile
         if (talik>lake_info%hsed-e8) then
            Tbot = max(MAGT,T0)
            Ttop = T0 + 4.0
         else
            Tbot = min(MAGT,T0-1.0)
            if (m_Zw(indx)>=2.0) then
               Ttop = T0 + 4.0
            else
               Ttop = min(max(m_radPars%tref,T0+0.5),T0+4.0)
            end if
         end if
         if (DEBUG) then
            print "(A, I0, A, I0, A, F10.4)", "Lake ", lake_info%id, &
                  ", soil column ", icol, ": bottom temperature ", Tbot-T0
         end if

         Ttk = 0.5 * (Ttop + T0)
         Tfz = 0.5 * (T0 + Tbot)
         satLW_tk = CalcSoilWaterSaturation(Ttk)
         satLW_fz = CalcSoilWaterSaturation(Tfz)
         Ks_tk = CalcSedHeatConductivity(0.4_r8, satLW_tk, Ks0)
         Ks_fz = CalcSedHeatConductivity(0.4_r8, satLW_fz, Ks0)
         dist = Ks_fz / Ks_tk * talik * (T0 - Tbot) / (Ttop - T0)
         if (dist>e8) then
            hMAGT = talik + dist
            dTtk = (T0 - Ttop) / talik
            dTfz = (Tbot - T0) / dist
            do ii = 1, SED_LAYER+1, 1
               if (m_Zs(ii)<=talik+m_Zs(1)) then
                  m_sedTemp(icol,ii) = Ttop + dTtk * (m_Zs(ii)-m_Zs(1))
               else if (m_Zs(ii)<=hMAGT+m_Zs(1)) then
                  m_sedTemp(icol,ii) = T0  + dTfz * (m_Zs(ii)-talik)
               else
                  m_sedTemp(icol,ii) = Tbot
               end if
            end do
         else
            dTtk = (Tbot - Ttop) / talik
            do ii = 1, SED_LAYER+1, 1
               m_sedTemp(icol,ii) = Ttop + dTtk * (m_Zs(ii)-m_Zs(1))
            end do
         end if

         ! construct ice profile
         do ii = 1, SED_LAYER+1, 1
            if (m_sedTemp(icol,ii)<T0) then
               m_sedIce(icol,ii) = m_sedpor(icol,ii)
            else
               m_sedIce(icol,ii) = 0.0_r8
            end if
         end do
      end do

      m_Ks = Ks0 / Cps / Rous 

      call UpdateSedWaterBoundIndex()
   end subroutine

end module soil_thermal_mod
