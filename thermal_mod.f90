module thermal_mod
!---------------------------------------------------------------------------------
! Purpose: Govern the heat diffusion of lake water
!
! "Simulation of lake evaporation with application to modeling lake level 
! variations of Harney-Malheur Lake, Oregon" [Hostetler et al., 1990]
!
!---------------------------------------------------------------------------------
   use shr_kind_mod,          only : r8
   use shr_param_mod
   use shr_ctrl_mod,          only : WATER_LAYER, lake_info, run_mode
   use shr_typedef_mod,       only : RungeKuttaCache1D
   use phy_utilities_mod
   use data_buffer_mod

   implicit none
   private
   public :: mem_tw
   public :: InitializeThermalModule, DestructThermalModule 
   public :: ThermalModuleSetup, ThermalModuleCallback
   public :: UpdateLakeWaterTopIndex, UpdateLakeIceThickness
   public :: GetBoundaryOutputs, HeatEquation
   ! module cache for Runge-Kutta
   type(RungeKuttaCache1D) :: mem_tw
   ! volumetric heat capacity (J/K/m3)
   real(r8), allocatable :: Cvt(:)
   ! buoyancy frequency square (s-2)
   real(r8), allocatable :: freq(:)
   ! turbulent diffusivity (m2/s)
   real(r8), allocatable :: df_turb(:)
   ! radiation-induced temperature change (K/m3/s)
   real(r8), allocatable :: rad(:)
   ! boundary conditions
   real(r8) :: lw_hr, sh_hr, lh_hr
   real(r8) :: hnet_hr, fmm_hr
   real(r8) :: hf_top, hf_btm
   real(r8) :: hf_bubble

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: Private data initialization, get/set and destruction
   !
   !------------------------------------------------------------------------------
   subroutine InitializeThermalModule()
      implicit none

      allocate(Cvt(WATER_LAYER+1))
      allocate(freq(WATER_LAYER+1))
      allocate(df_turb(WATER_LAYER+1))
      allocate(rad(WATER_LAYER+1))
      allocate(mem_tw%K1(WATER_LAYER+1))
      allocate(mem_tw%K2(WATER_LAYER+1))
      allocate(mem_tw%K3(WATER_LAYER+1))
      allocate(mem_tw%K4(WATER_LAYER+1))
      allocate(mem_tw%K5(WATER_LAYER+1))
      allocate(mem_tw%K6(WATER_LAYER+1))
      allocate(mem_tw%nxt4th(WATER_LAYER+1))
      allocate(mem_tw%nxt5th(WATER_LAYER+1))
      allocate(mem_tw%interim(WATER_LAYER+1))
      allocate(mem_tw%rerr(WATER_LAYER+1))

      call InitializeWaterStateVariables()
      Cvt = 0.0_r8
      freq = 0.0_r8
      df_turb = 0.0_r8
      rad = 0.0_r8
      hf_top = 0.0_r8
      hf_btm = 0.0_r8
      hf_bubble = 0.0_r8
      lw_hr = 0.0_r8
      sh_hr = 0.0_r8
      lh_hr = 0.0_r8
      hnet_hr = 0.0_r8
      fmm_hr = 0.0_r8
   end subroutine

   subroutine DestructThermalModule()
      implicit none

      deallocate(Cvt)
      deallocate(freq)
      deallocate(df_turb)
      deallocate(rad)
      deallocate(mem_tw%K1)
      deallocate(mem_tw%K2)
      deallocate(mem_tw%K3)
      deallocate(mem_tw%K4)
      deallocate(mem_tw%K5)
      deallocate(mem_tw%K6)
      deallocate(mem_tw%nxt4th)
      deallocate(mem_tw%nxt5th)
      deallocate(mem_tw%interim)
      deallocate(mem_tw%rerr)
   end subroutine

   subroutine ThermalModuleSetup()
      implicit none
      real(r8) :: cita, icita, temp
      real(r8) :: Kml, wind, w10, lat
      integer :: ii, top, icol

      top = m_lakeWaterTopIndex
      call UpdateWaterDensity()
      w10 = m_surfData%wind
      wind = ConvertWindSpeed(w10, 2.0d0) 
      lat = lake_info%latitude
      fmm_hr = Roua * Cd10 * (w10**3)
      ! wind-driven eddy diffusivity
      call CalcBruntVaisalaFreq(m_dZw, m_wrho, freq)
      if (m_Hice<e8) then
         call CalcEddyDiffusivity(m_Zw, freq, wind, lat, m_Kt)
         m_Kt(1) = max(m_Kt(1), 1.0d-5)
      else
         call CalcEnhancedDiffusivity(freq, m_waterIce, m_Kt)
         if (top<=WATER_LAYER+1) then
            m_Kt(top) = max(m_Kt(top), 2.0d-6)
         end if
      end if
      ! bubble-induced heat flux
      if (m_waterTemp(top)<T0+4._r8) then
         hf_bubble = min(1.d9*m_upbflux, Lf*Roul*m_dZw(1)/3.6d3)
      else
         hf_bubble = 0._r8
      end if
      df_turb = m_Kt 
      do ii = 1, WATER_LAYER+1, 1
         cita = 1.0 - m_waterIce(ii)
         icita = m_waterIce(ii)
         temp = m_waterTemp(ii)
         Cvt(ii) = Cpl*Roul*cita + Cpi*Roui*icita
         ! molecular heat conductivity
         Kml = CalcWaterHeatDiffusivity(temp)
         m_Kv(ii) = 1.25 * m_Kt(ii)
         m_Kt(ii) = m_Kt(ii) + Kml
         if (ii<top) then
            m_Ktb(ii) = Kml * DeltaD / 0.1_r8
         else
            m_Ktb(ii) = Kml
         end if
         ! radiation
         icol = m_soilColInd(ii)
         if ((m_Hice>e8 .and. ii==top-1) .or. (m_Hice<e8 .and. ii==top)) then
            rad(ii) = (m_Iab(ii) + hf_bubble) / m_dZw(ii) / Cvt(ii)
         else
            rad(ii) = m_Iab(ii) / m_dZw(ii) / Cvt(ii)
         end if
      end do
      call CalcBoundaryHeatFlux()
   end subroutine

   subroutine ThermalModuleCallback(dt)
      implicit none
      real(r8), intent(in) :: dt
      real(r8) :: rateS, rateE, fs, fe 
      real(r8) :: snowMax, mS2E, Roun

      call AdjustIceTemp()
      call UpdateWaterDensity()
      call ConvectiveMixing(dt)
      if (m_Hice>e8) then
         Roun = sa_params(Param_Roun)
         call CalcSnowBurdenLimit(m_Hice, m_Hgrayice, Roun, snowMax)
         call CalcSnowLayerChangeRate(rateS)
         call CalcGrayIceLayerMeltRate(rateE)
         m_Hsnow = max(m_Hsnow+rateS*dt, 0.0) 
         if (m_Hsnow>=0.1) then
            mS2E = max(0d0, (m_Hsnow-snowMax)*Roun/Roue)
            m_Hgrayice = max(m_Hgrayice+mS2E+rateE*dt, 0d0)
            m_Hsnow = max(0.0, min(m_Hsnow, snowMax))
         else
            m_Hgrayice = max(m_Hgrayice+rateE*dt, 0d0)
         end if
      else
         m_Hsnow = 0.0_r8
         m_Hgrayice = 0.0_r8
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Discretize Heat differential equation by Finite Difference Method.
   !          This sub-routine returns time derivatives of lake water temperature
   !          Cv*dT/dt = d(K*dT/dz)/dz - dQ/dz  (Q: incident solar energy) 
   !          Mean annual lake bottom temperature is applied at the periphery
   !          of a lake (West and Plug, 2008).
   !
   !------------------------------------------------------------------------------
   subroutine HeatEquation(temp, dtemp)
      implicit none
      real(r8), intent(in)  :: temp(WATER_LAYER+1)     ! temperature column vector
      real(r8), intent(out) :: dtemp(WATER_LAYER+1)    ! temperature time derivative
      real(r8) :: qb, qt, aa, bb, Az1, Az2 
      real(r8) :: dT, dT1, dT2
      integer  :: ii, icol, top

      ! Top and Bottom boundary conditions
      top = m_lakeWaterTopIndex
      qt = hf_top
      do ii = 1, WATER_LAYER+1, 1
         icol = m_soilColInd(ii)
         if (m_Hice<e8) then
            if (icol<NSCOL) then
               qb = 0._r8
            else
               dT = m_sedTemp(icol,1) - temp(WATER_LAYER+1)
               qb = m_Ktb(WATER_LAYER+1) * dT / DeltaD
            end if
         else
            qb = m_Ktb(ii) * (m_sedTemp(icol,1) - temp(ii)) / DeltaD
         end if
         if(ii==1) then
            aa = 0.5 * (m_Kt(ii+1) + m_Kt(ii))
            dT1 = (temp(ii+1) - temp(ii)) / (m_Zw(ii+1) - m_Zw(ii))
            Az1 = m_Az(ii)
            dtemp(ii) = (aa*Az1*dT1 + qb*m_dAz(ii) - qt*m_Az(ii)) / &
                  m_dZw(ii) / m_Az(ii) + rad(ii)
         else if(ii==WATER_LAYER+1) then
            bb = 0.5 * (m_Kt(ii-1) + m_Kt(ii))
            dT2 = (temp(ii) - temp(ii-1)) / (m_Zw(ii) - m_Zw(ii-1))
            Az2 = m_Az(ii) 
            dtemp(ii) = (qb*m_dAz(ii) - bb*Az2*dT2) / m_dZw(ii) / &
                  m_Az(ii) + rad(ii)
         else
            aa = 0.5 * (m_Kt(ii+1) + m_Kt(ii))
            bb = 0.5 * (m_Kt(ii-1) + m_Kt(ii))
            dT1 = (temp(ii+1) - temp(ii)) / (m_Zw(ii+1) - m_Zw(ii))
            dT2 = (temp(ii) - temp(ii-1)) / (m_Zw(ii) - m_Zw(ii-1))
            Az1 = m_Az(ii)
            Az2 = m_Az(ii) 
            dtemp(ii) = (aa*Az1*dT1 + qb*m_dAz(ii) - bb*Az2*dT2) / &
                  m_dZw(ii) / m_Az(ii) + rad(ii)
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate top and bottom boundary heat flux (W/m2)
   !          For snow scenario, see Appendix B of "Temperature and dissolved
   !          oxygen simulations in a lake with ice cover" (Fang, 1994)
   !          Roughness length of mud is 0.005 m.
   !
   !          Hsa: bulk heat transfer coefficient at air/snow interface
   !          (Fertuck et al., 1971)
   !
   !------------------------------------------------------------------------------
   subroutine CalcBoundaryHeatFlux()
      implicit none
      real(r8) :: lw_net, lh_net, sh_net, hgray
      real(r8) :: runoff_net, Hsa, hice, hsnow
      real(r8) :: Ttop, Twet, w10, wind, RH, Tair
      real(r8) :: icita, dtemp, temp, Rn, ps
      real(r8) :: sw_top, sw_btm, sw_remain
      real(r8) :: Roun, Hscaler
      integer :: ii, icol, top

      Ttop = m_waterTemp(1)
      w10 = m_surfData%wind
      wind = ConvertWindSpeed(w10, 2.0d0)
      Tair = m_surfData%temp
      RH = m_surfData%RH
      ps = m_surfData%pressure
      Roun = sa_params(Param_Roun)
      Hscaler = sa_params(Param_Hscale)
      !Hsa = 4.29 * max(1.d-1, wind)  ! units: W/m2/K 
      Hsa = 3.90 * max(1.d-1, wind)  ! units: W/m2/K
      if (m_Hsnow>e8 .or. m_Hgrayice>e8) then
         if (Tair<T0) then
            lw_net = 0.0_r8
         else
            if (m_Hsnow>e8) then
               lw_net = Epsn*Stefan*(T0**4) - Epsn*m_surfData%lw
            else
               lw_net = Epse*Stefan*(T0**4) - Epse*m_surfData%lw
            end if
            lw_net = min(lw_net, 0.0_r8)
         end if
         lh_net = 0.0_r8
         hsnow = m_Hsnow
         hice = m_Hice
         hgray = m_Hgrayice
         sh_net = (T0 - Tair) / (hsnow/Kn0 + hice/Ki0 + &
            hgray/ke0 + 1.0/Hsa)
         runoff_net = 0.0_r8

         lh_hr = lh_net 
         sh_hr = sh_net 
         if (m_Hsnow>e8) then
            lw_hr = min(Epsn*Stefan*(T0**4), Epsn*m_surfData%lw)
         else
            lw_hr = min(Epse*Stefan*(T0**4), Epse*m_surfData%lw)
         end if
         hnet_hr = m_surfData%srd - lw_net - sh_hr - lh_hr
      else if (m_Hice>e8) then
         lw_net = Epsi*Stefan*(Ttop**4) - Epsi*m_surfData%lw
         lh_net = CalcLatentHeatIce(Ttop, Tair, RH, wind)
         sh_net = CalcSensibleHeat(Ttop, Tair, wind)
         runoff_net = 0.0_r8

         lh_hr = lh_net
         sh_hr = sh_net
         lw_hr = Epsi*Stefan*(Ttop**4)
         hnet_hr = m_surfData%srd - lw_net - sh_hr - lh_hr 
      else
         lw_net = Epsw*Stefan*(Ttop**4) - Epsw*m_surfData%lw
         lh_net = Hscaler * CalcLatentHeatWater(Ttop, Tair, RH, wind)
         !Rn = m_surfData%srd - lw_net
         !lh_net = CalcLatentHeatWater(Ttop, Tair, RH, wind, ps, Rn)
         sh_net = Hscaler * CalcSensibleHeat(Ttop, Tair, wind)
         Twet = CalcWetBubTemp(Tair, RH, m_surfData%pressure)
         if (m_surfData%rainfall>e8) then
            runoff_net = Cpl*Roul*(Twet-Ttop)*m_surfData%rainfall
         else if (m_surfData%snowfall>e8) then
            runoff_net = (Cpn*(Twet-T0) + Cpl*(T0-Ttop) - Lf) * &
                        m_surfData%snowfall * Roun
         else
            runoff_net = 0.0_r8
         end if
         
         lh_hr = lh_net
         sh_hr = sh_net
         lw_hr = Epsw*Stefan*(Ttop**4)
         hnet_hr = m_surfData%srd - lw_net - sh_hr - lh_hr
      end if
      hf_top = (lw_net + sh_net + lh_net - runoff_net) / Cvt(1)
      ! bottom heat flux
      hf_btm = 0._r8
      do ii = 1, WATER_LAYER+1, 1
         icol = m_soilColInd(ii) 
         temp = m_waterTemp(ii)
         hf_btm = hf_btm + m_Ktb(ii) * (m_sedTemp(icol,1) - temp) / DeltaD * &
            m_dAz(ii) / m_Az(1)
      end do
      if (m_Hice<e8) then
         ! see Imberger [1985] for the calculation of trapped shortwave
         ! radiation but it should be minor when buoyant flux is negative
         m_Heff = -sh_net - lh_net - lw_net
      else
         m_Heff = 0.0_r8
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Adjust temperature and ice percentage of every layer
   !          In "Measuring the sensitivity of southern Wisconsin lake ice to
   !          climate variations and lake depth using a numerical model" (Vavrus
   !          et al., 1996), the first or last 1 cm ice is formed or melted
   !          instantaneously during ice onset or break.
   !
   !------------------------------------------------------------------------------
   subroutine AdjustIceTemp()
      implicit none
      real(r8) :: tmp1, tmp2, dH1, dH2
      real(r8) :: dv1, dv2
      integer :: ii, jj

      ! Water Freezing
      do ii = 1, WATER_LAYER+1, 1
         if (m_Hice+e8>lake_info%maxdepth) then
            exit  ! no liquid
         end if
         if (ii<m_lakeWaterTopIndex) then
            if (m_waterIce(ii)<1.0-e8 .and. m_waterTemp(ii)<T0-e8) then
               tmp1 = Cpi*(T0-m_waterTemp(ii))*m_waterIce(ii) + &
                     Cpl*(T0-m_waterTemp(ii))*(1.0-m_waterIce(ii))
               tmp2 = Lf*(1.0-m_waterIce(ii))
               if (tmp2>tmp1) then
                  m_waterIce(ii) = m_waterIce(ii) + tmp1/Lf
                  m_waterTemp(ii) = T0
               else
                  m_waterIce(ii) = 1.0
                  m_waterTemp(ii) = T0 - (tmp1-tmp2)/Cpi
               end if
            end if
         else
            jj = max(1,m_lakeWaterTopIndex-1)
            if (m_waterIce(jj)>1.0-e8) then
               jj = jj + 1
            end if
            do while (m_waterIce(jj)<1.0-e8 .and. m_waterTemp(ii)<T0-e8)
               tmp1 = Cpi*(T0-m_waterTemp(ii))*m_waterIce(ii) + &
                     Cpl*(T0-m_waterTemp(ii))*(1.0-m_waterIce(ii))
               tmp2 = Lf*(1.0-m_waterIce(jj))
               dv1 = m_dZw(ii) * m_Az(ii)
               dv2 = m_dZw(jj) * m_Az(jj)
               dH1 = tmp1 * dv1
               dH2 = tmp2 * dv2
               if (dH2>dH1) then
                  m_waterIce(jj) = m_waterIce(jj) + dH1/dv2/Lf
                  m_waterTemp(ii) = T0
                  exit
               else
                  m_waterIce(jj) = 1.0
                  if (ii==jj) then
                     m_waterTemp(ii) = T0 - (dH1-dH2)/dv1/Cpi
                  else
                     m_waterTemp(ii) = T0 - (dH1-dH2)/dv1/Cpl
                  end if
                  jj = jj + 1
                  if (jj>ii) then
                     exit
                  end if
               end if
            end do
            call UpdateLakeWaterTopIndex()
         end if
      end do
      ! Ice Thawing
      do ii = 1, WATER_LAYER+1, 1
         if (m_Hice<e8) then
            exit  ! no ice or freezing cycle
         end if
         jj = m_lakeWaterTopIndex - 1
         if (jj<ii) then
            exit
         end if
         do while (m_waterTemp(ii)>T0+e8 .and. m_waterIce(jj)>e8)
            tmp1 = Cpi*(m_waterTemp(ii)-T0)*m_waterIce(ii) + &
                  Cpl*(m_waterTemp(ii)-T0)*(1.0-m_waterIce(ii))
            tmp2 = Lf*m_waterIce(jj)
            dv1 = m_dZw(ii) * m_Az(ii)
            dv2 = m_dZw(jj) * m_Az(jj)
            dH1 = tmp1 * dv1
            dH2 = tmp2 * dv2
            if (dH2>dH1) then
               m_waterIce(jj) = m_waterIce(jj) - dH1/dv2/Lf
               m_waterTemp(ii) = T0
               exit
            else
               m_waterIce(jj) = 0.0
               if (ii==jj) then
                  m_waterTemp(ii) = T0 + (dH1-dH2)/dv1/Cpl
               else
                  m_waterTemp(ii) = T0 + (dH1-dH2)/dv1/Cpi
               end if
               jj = jj - 1
               if (jj<ii) then
                  exit
               end if
            end if
         end do
         call UpdateLakeWaterTopIndex()
      end do
      call UpdateLakeIceThickness()
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Adjust temperature triggered by convective mixing. The convective
   !          mixing is judged by the comparison of TKE dissipation and
   !          stability [Wuest et al., 2000]. And TKE dissipation is calculated
   !          according to Henderson-Sellers (1985).
   !
   !------------------------------------------------------------------------------
   subroutine ConvectiveMixing(dt)
      implicit none
      real(r8), intent(in) :: dt
      real(r8) :: tavg, tmv, Vz, Mv, mavg
      real(r8) :: Wstr, w10, As, Pkin
      real(r8) :: Epot, drho, maxrho, zmz, mz
      real(r8) :: zepi, zmepi, mepi, Vepi
      real(r8) :: zhpi, zmhpi, mhpi, Vhpi
      real(r8) :: tbot, tcur
      integer :: ii, jj, top, bottom 
      integer :: indx, btmIdx(2)

      top = m_lakeWaterTopIndex
      bottom = WATER_LAYER + 1

      if (m_Hice<e8) then
         w10 = m_surfData%wind 
         As = 1.d-6 * lake_info%Asurf
         Wstr = sa_params(Param_Wstr) * (1.0 - exp(-0.3*As))
         call CalcTotalKineticPower(lake_info, w10, Pkin)
         m_Ekin = m_Ekin + Wstr * Pkin * dt
      else
         m_Ekin = 0.0_r8
      end if

      ! surface boundary layer
      m_mixTopIndex = top
      ii = top + 1 
      do while (ii<=bottom)
         zepi = sum( m_dZw(top:ii-1) )
         zmepi = 0.0_r8
         tmv = 0.0_r8
         mepi = 0.0_r8
         do jj = top, ii-1, 1
            Mv = m_wrho(jj) * m_dZw(jj) * m_Az(jj)
            mepi = mepi + Mv
            tmv = tmv + m_dZw(jj) * m_Az(jj)
            zmepi = zmepi + (sum(m_dZw(top:jj))-0.5*m_dZw(jj)) * Mv
         end do
         zmepi = zmepi / mepi
         zmz = m_dZw(ii)
         mepi = mepi / tmv 
         Vepi = 1.d-6 * tmv 
         Vz = 1.d-6 * m_Az(ii) * zmz
         drho = m_wrho(ii) - mepi
         Epot = G*drho*Vepi*Vz/(Vepi+Vz)*(zepi - zmepi + 0.5*zmz)
         if (m_Ekin<Epot) then
            exit
         end if
         if (Epot>0) then
            m_Ekin = m_Ekin - Epot
         end if
         m_mixTopIndex = ii
         ii = ii + 1
      end do
      ! TKE totally dissipation if fully mixing
      if (m_mixTopIndex==bottom) then
         m_Ekin = 0.0_r8
      end if

      m_mixTopIndex = min(m_mixTopIndex, bottom)
      ! convective mixing
      if (m_mixTopIndex>top) then
         tavg = 0.0_r8 
         tmv = 0.0_r8
         mavg = 0.0_r8
         do ii = top, m_mixTopIndex, 1
            Mv = m_dZw(ii) * m_Az(ii)
            tmv = tmv + Mv
            tavg = tavg + m_waterTemp(ii) * Mv
            mavg = mavg + m_wrho(ii) * Mv 
         end do
         tavg = tavg / tmv
         mavg = mavg / tmv
         m_waterTemp(top:m_mixTopIndex) = tavg
         m_wrho(top:m_mixTopIndex) = mavg
         m_Hmix(1) = sum( m_dZw(top:m_mixTopIndex) )
      else
         m_Hmix(1) = 0.0_r8
      end if

      ! density-induced bottom mixing
      do while (m_mixTopIndex<bottom)
         maxrho = m_wrho(bottom)
         tbot = m_waterTemp(bottom) - T0 - 4._r8
         indx = bottom
         do ii = top, bottom-1, 1
            tcur = m_waterTemp(ii) - T0 - 4._r8
            if (maxrho+1.d-6<m_wrho(ii)) then
               maxrho = m_wrho(ii)
               indx = ii
            else if (((tbot<-e8 .and. tcur>e8) .or. (tbot>e8 .and. tcur<-e8)) &
               .and. abs(tbot-tcur)>0.1) then
               maxrho = maxval(m_wrho) + 1.d-3 
               indx = ii
            end if
         end do
         ! no inverse density gradient
         if (indx==bottom) then
            exit
         end if
         btmIdx(1) = indx
         btmIdx(2) = bottom
         ! convective mixing
         tavg = 0.0_r8
         tmv = 0.0_r8
         mavg = 0.0_r8
         do ii = btmIdx(1), btmIdx(2), 1
            Mv = m_dZw(ii) * m_Az(ii)
            tmv = tmv + Mv
            tavg = tavg + m_waterTemp(ii) * Mv
            mavg = mavg + m_wrho(ii) * Mv
         end do
         tavg = tavg / tmv
         mavg = mavg / tmv
         m_waterTemp(btmIdx(1):btmIdx(2)) = tavg
         m_wrho(btmIdx(1):btmIdx(2)) = mavg
      end do

      ! search the range of bottom mixing layer
      m_mixBotIndex = bottom
      do ii = bottom-1, m_mixTopIndex+1, -1
         if (abs(m_wrho(ii)-m_wrho(bottom))>=0.02_r8) then
            exit
         end if
         m_mixBotIndex = ii
      end do
      if (m_mixBotIndex<bottom) then
         m_Hmix(2) = sum( m_dZw(m_mixBotIndex:bottom) )
      else
         m_Hmix(2) = 0.0_r8
      end if

      m_HbLayer(1) = sum( m_dZw(1:m_mixTopIndex) ) 
      m_HbLayer(2) = lake_info%maxdepth - sum( m_dZw(m_mixBotIndex:bottom) )
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate the depth of snow cover depth from "Temperature and 
   !          Dissolved Oxygen Simulations for a Lake with Ice Cover" [Fang et 
   !          al., 1994]
   !          According to Vavrus et al. (1996), the maximum snow depth is
   !          decided by Archimedes' principle and if this maximum snow depth is
   !          exceeded, then the weight of the snow is assumed to depress the
   !          ice below the surface water level, resulting a conversion of the
   !          submerged snow layer into gray ice.
   !          Boike et al. (2015) observed that in spring, the frozen surfaces
   !          of the lakes were normally kept snow free by wind action.
   !
   !------------------------------------------------------------------------------
   subroutine CalcSnowLayerChangeRate(rate)
      implicit none
      real(r8), intent(out) :: rate       ! unit: m/s
      real(r8) :: rad, lh, sh, runoff
      real(r8) :: lw, hnet, Roun, wind

      if (m_surfData%temp>T0+e8 .or. m_radPars%season==1) then
         ! absorbed solar energy, latent heat, sensible heat and runoff heat
         wind = ConvertWindSpeed(m_surfData%wind, 2.0d0)
         rad = m_surfData%srd - m_Idw(1)
         lw = max( Epsn*m_surfData%lw-Epsn*Stefan*(T0**4), 0.0_r8 )
         sh = CalcSensibleHeat(T0, m_surfData%temp, wind)                               
         lh = 0.0_r8 
         runoff = Cpl*Roul*(m_surfData%temp-T0)*m_surfData%rainfall                         
         runoff = max( runoff, 0.0_r8 )
      else
         rad = 0.0_r8 
         lw = 0.0_r8
         lh = 0.0_r8
         sh = 0.0_r8
         runoff = 0.0_r8
      end if
      Roun = sa_params(Param_Roun)
      hnet = max(0d0, rad+lw-sh-lh+runoff)
      rate = m_surfData%snowfall - hnet/(Roun*Lf)
   end subroutine

   subroutine CalcGrayIceLayerMeltRate(rate)
      implicit none
      real(r8), intent(out) :: rate       ! units: m/s
      real(r8) :: rad, lh, sh, runoff
      real(r8) :: lw, hnet, wind

      if ((m_surfData%temp>T0+e8 .or. m_radPars%season==1) &
            .and. m_Hsnow<e8) then
         wind = ConvertWindSpeed(m_surfData%wind, 2.0d0)
         rad = m_surfData%srd - m_Idw(1)
         lw = max( Epse*m_surfData%lw-Epse*Stefan*(T0**4), 0.0_r8 )
         sh = CalcSensibleHeat(T0, m_surfData%temp, wind)
         lh = 0.0_r8 
         runoff = Cpl*Roul*(m_surfData%temp-T0)*m_surfData%rainfall
         runoff = max( runoff, 0.0_r8 )
      else
         rad = 0.0_r8
         lw = 0.0_r8
         sh = 0.0_r8
         lh = 0.0_r8
         runoff = 0.0_r8
      end if
      hnet = max(0d0, rad+lw-sh-lh+runoff)
      rate = -hnet/(Roue*Lf)
   end subroutine

   subroutine GetBoundaryOutputs(sh, lh, fmm, lw, hnet, fsed, turbdiff)
      implicit none
      real(r8), intent(out) :: sh      ! upward sensible heat (W/m2)
      real(r8), intent(out) :: lh      ! upward latent heat (W/m2)
      real(r8), intent(out) :: fmm     ! momentum energy flux (W/m2)
      real(r8), intent(out) :: lw      ! upward longwave radiation (W/m2)
      real(r8), intent(out) :: hnet    ! net downward heat flux (W/m2)
      real(r8), intent(out) :: fsed    ! net upward sed heat flux (W/m2)
      real(r8), intent(out) :: turbdiff(:) ! turbulent diffusivity (m2/s)

      sh = sh_hr
      lh = lh_hr
      fmm = fmm_hr
      lw = lw_hr
      hnet = hnet_hr
      fsed = hf_btm
      turbdiff = df_turb
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: initialize thermal module related state variables.
   !
   !------------------------------------------------------------------------------
   subroutine InitializeWaterStateVariables()
      implicit none
      real(r8) :: tair_avg

      if (m_radPars%tref>T0+4) then 
         m_waterTemp = m_radPars%tref
      else
         m_waterTemp = T0 + 4
      end if
      m_waterIce = 0.0_r8
      m_Hsnow = 0.0_r8
      m_Hgrayice = 0.0_r8
      m_Kt = 1.5d-7
      m_Kv = 1.5d-6
      m_Ktb = 1.5d-7
      m_Iab = 0.0_r8
      m_Idw = 0.0_r8
      m_wrho = Roul
      m_Hmix = 0.0_r8
      m_HbLayer = 0.0_r8
      m_Heff = 0.0_r8
      m_Ekin = 0.0_r8
      m_mixTopIndex = 1
      m_mixBotIndex = WATER_LAYER+1
      call UpdateLakeWaterTopIndex()
      call UpdateLakeIceThickness()
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update lake water density depending on temperature and salinity.
   !          "Stratification, mixing and transport processes in Lake Kivu" 
   !          (Schmid and Wuest, 2012).
   !
   !------------------------------------------------------------------------------
   subroutine UpdateWaterDensity()
      implicit none
      real(r8) :: DCO2(WATER_LAYER+1)  ! units: g/kg
      real(r8) :: DCH4(WATER_LAYER+1)  ! units: g/kg

      if (count(m_waterTemp>T0+100)>0) then
         if (trim(run_mode)=='sensitivity') then
            call Endrun(lake_info%sampleId, "Water temperature out of range")
         else
            call Endrun(lake_info%id, "Water temperature out of range")
         end if
      end if

      ! 2.75 is a compensation for non-CO2 substances
      DCO2 = 2.75 * m_waterSubCon(Wco2,:) * MasCO2 * 1.d-3
      DCH4 = m_waterSubCon(Wch4,:) * MasCH4 * 1.d-3
      !salinity = 0._r8 ! 15.0d-3 for freshwater (not used due to DCO2)
      ! update water density
      where (m_waterTemp<T0)
         m_wrho = Roui
      elsewhere
         m_wrho = ( 9.99868d+2 + 1.0d-3 * (65.185*(m_waterTemp-T0)- &
            8.4878*(m_waterTemp-T0)**2+0.05607*(m_waterTemp-T0)**3) ) * &
            ( 1.0 + 0.75d-3*m_SAL + 0.284d-3*DCO2 - 1.25d-3*DCH4 )
      end where
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update lake ice thickness and top index of lake free water layer.
   !
   !------------------------------------------------------------------------------
   subroutine UpdateLakeIceThickness()
      implicit none
      integer :: ii

      m_Hice = 0.0_r8
      do ii = 1, WATER_LAYER, 1
         m_Hice = m_Hice + m_waterIce(ii)*m_dZw(ii)
         if (m_waterIce(ii)>e8 .and. m_lakeWaterTopIndex<=ii) then
            if (trim(run_mode)=='sensitivity') then
               call Endrun(lake_info%sampleId, 'incorrect ice layers ' // &
                  'under the ice-water interface')
            else
               call Endrun(lake_info%id, 'incorrect ice layers under ' // &
                  'the ice-water interface')
            end if
         end if
      end do
   end subroutine

   subroutine UpdateLakeWaterTopIndex()
      implicit none
      logical, save :: wfirstime = .True.
      integer :: ii

      do ii = 1, WATER_LAYER+1, 1
         if (m_waterIce(ii)<e8) then
            m_lakeWaterTopIndex = ii
            return
         end if
      end do
      m_lakeWaterTopIndex = WATER_LAYER + 2
      if (wfirstime) then
         print "(I0, A, A)", lake_info%id, ': lake water is ", &
                               "totally frozen!!'
         wfirstime = .False.
      end if
   end subroutine

end module thermal_mod
