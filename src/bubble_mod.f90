module bubble_mod
!---------------------------------------------------------------------------------
! Purpose: this module governs the transportation of bubble (ebullition process) 
!          in the water, which is mainly from "Modeling bubbles and dissolved gases 
!          in the ocean" [Liang et al., 2011]. And some parameters are calculated
!          according to "Bubbles and the air-sea exchange of gases in near-saturation 
!          conditions" [Woolf, D. and Thorpe, S., 1991].
!          The bubble dynamics run in single bubble mode.
!
!---------------------------------------------------------------------------------
   use shr_kind_mod,          only : r8
   use shr_ctrl_mod, e8 => SHR_CTRL_E8, inft => INFINITESIMAL_E8
   use phy_utilities_mod
   use data_buffer_mod
   use shr_param_mod

   implicit none
   private
   public :: bubble_fluxes
   public :: InitializeBubbleModule, DestructBubbleModule 
   public :: BubbleModuleSetup, BubbleModuleCallback
   public :: BubbleDynamics
   ! bubble gas exchange cache vector (umol/(m3*mm*s))
   real(r8), allocatable :: Vex(:,:,:)
   ! bubble radius cache vector (m)
   real(r8), allocatable :: Vr(:,:)
   ! bubble gas concentration vector (umol/(m3*mm))
   real(r8), allocatable :: Vbg(:,:)
   ! buoyant velocity auxiliary vector
   real(r8), allocatable :: Vab(:,:)
   ! bubble radius gradient auxiliary vector
   real(r8), allocatable :: Var1(:,:)
   ! bubble radius gradient auxiliary vector
   real(r8), allocatable :: Var2(:,:)
   ! gas exchange auxiliary vector
   real(r8), allocatable :: Vag1(:,:,:)
   real(r8), allocatable :: Vag2(:,:,:)
   ! top surface bubble fluxes (umol/m2/s)
   real(r8), allocatable :: bubble_fluxes(:)
   ! time step (s)
   real(r8) :: tdelta

contains
   subroutine InitializeBubbleModule()
      implicit none     

      allocate(Vex(NGAS,NRLAYER+1,WATER_LAYER+1))
      allocate(Vr(NRLAYER+1,WATER_LAYER+1))
      allocate(Vbg(NGAS,NRLAYER+1))
      allocate(Vab(NRLAYER+1,WATER_LAYER+1))
      allocate(Var1(NRLAYER+1,WATER_LAYER+1))
      allocate(Var2(NRLAYER+1,WATER_LAYER+1))
      allocate(Vag1(NGAS,NRLAYER+1,WATER_LAYER+1))
      allocate(Vag2(NGAS,NRLAYER+1,WATER_LAYER+1))
      allocate(bubble_fluxes(NGAS))
      
      call InitializeBubbleStateVariables()
      call InitBubbleRadiusField()
      Vex = 0.0_r8
      Vbg = 0.0_r8
      Vab = 0.0_r8
      Var1 = 0.0_r8
      Var2 = 0.0_r8
      Vag1 = 0.0_r8
      Vag2 = 0.0_r8
      bubble_fluxes = 0.0_r8
   end subroutine

   subroutine DestructBubbleModule()
      implicit none

      deallocate(Vex)
      deallocate(Vr)
      deallocate(Vbg)
      deallocate(Vab)
      deallocate(Var1)
      deallocate(Var2)
      deallocate(Vag1)
      deallocate(Vag2)
      deallocate(bubble_fluxes)
   end subroutine

   subroutine BubbleModuleSetup(isHourNode)
      implicit none
      logical, intent(in) :: isHourNode
      real(r8) :: temp, pH, rr, wb, tmp
      real(r8) :: depth, vsc, gama, solubility
      real(r8) :: diffusivity, nusselt, pressure
      integer :: ii, jj, kk

      if (.NOT. isHourNode) then
         return
      end if

      ! Producing bubbles based on sediment fluxes
      call GetBubbleReleaseMagnitude(Vbg)
      ! Updating transient variables of bubble dynamics
      do ii = 1, WATER_LAYER+1, 1
         temp = max(m_waterTemp(ii), T0)
         vsc = m_dVsc(ii) / Roul
         depth = m_Zw(ii) 
         gama = CalcSurfaceTension(temp)
         do jj = 1, NRLAYER+1, 1
            rr = Vr(jj,ii) + inft
            Vab(jj,ii) = CalcBuoyantVelocity(rr, vsc)
            tmp = 3.0*m_surfData%pressure + 3.0*Roul*G*depth + 4.0*gama/rr
            Var1(jj,ii) = 0.75*R*temp/(Pi*rr**2)/tmp
            Var2(jj,ii) = rr*Roul*G*Vab(jj,ii)/tmp
         end do
      end do
      pH = LKpH(lake_info%itype)
      do kk = Wn2, Wch4, 1
         ! no gas exchange at water-sediment interface for code stability
         do ii = 1, WATER_LAYER, 1
            temp = max(m_waterTemp(ii), T0)
            vsc = m_dVsc(ii) / Roul
            depth = m_Zw(ii)
            gama = CalcSurfaceTension(temp)
            solubility = CalcHenrySolubility(kk,temp,pH)
            diffusivity = CalcGasDiffusivityInWater(kk,temp)
            do jj = 1, NRLAYER+1, 1
               rr = Vr(jj,ii) + inft
               nusselt = CalcNusseltNumber(kk,rr,temp,vsc)
               pressure = m_surfData%pressure + Roul*G*depth + 2.0*gama/rr
               Vag1(kk,jj,ii) = -4.0*Pi*rr*diffusivity*nusselt* &
                                 solubility*pressure
               Vag2(kk,jj,ii) = 4.0*Pi*rr*diffusivity*nusselt* &
                                 m_waterSubCon(kk,ii)
            end do
         end do
      end do
   end subroutine

   subroutine BubbleModuleCallback(dt)
      implicit none
      real(r8), intent(in) :: dt
      integer :: ii, jj, top

      top = m_lakeWaterTopIndex
      do ii = top, WATER_LAYER+1, 1
         ! gas transfer from bubble to water
         m_gasExchange(:,ii) = 0._r8
         do jj = 1, NRLAYER+1, 1
            m_gasExchange(:,ii) = m_gasExchange(:,ii) - Vex(:,jj,ii)
         end do
      end do
      ! assume no transfer in ice layers
      m_gasExchange(:,1:top-1) = 0.0_r8
      call UpdateBubbleFlux(dt)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Discretize single bubble evolution by Finite Difference Method.
   !          This sub-routine returns time derivatives of bubble gas concentration.
   !          Four gases are included: N2, O2, CO2 and CH4.
   !          
   !
   !------------------------------------------------------------------------------
   subroutine BubbleDynamics(isHourNode)
      implicit none
      logical, intent(in) :: isHourNode
      real(r8) :: Ndm(NGAS)               ! gas exchange rate
      real(r8) :: ratio(NGAS)             ! gas ratio
      real(r8) :: Vbg_old(NGAS)
      real(r8) :: pressure, gama, pos
      real(r8) :: wb, rr, rr_old, temp
      real(r8) :: Cb, tmp1, tmp2, yn
      real(r8) :: ice_water_interface
      real(r8) :: dt, rratio              ! time step and radius rate
      integer :: rIndx, locIndx, top, ii

      if (.NOT. isHourNode) then
         return
      end if
      
      if (m_lakeWaterTopIndex>WATER_LAYER+1) then
         m_bubbleGasCon = 0.0_r8
         return
      end if
      top = m_lakeWaterTopIndex
      ice_water_interface = -m_Zw(top) 
      do rIndx = 1, NRLAYER+1, 1
         ! bubble release from sediment-water interface
         locIndx = WATER_LAYER + 1
         pos = -lake_info%depth
         rr = m_Rb0(rIndx) 
         rr_old = rr
         do while (pos+e8<ice_water_interface .and. rr>e8)
            if (pos>-m_Zw(locIndx)+e8) then
               Vr(rIndx,locIndx) = rr_old
               Vex(:,rIndx,locIndx) = Cb*Ndm
               m_bubbleGasCon(:,rIndx,locIndx) = Vbg_old
               locIndx = locIndx - 1
            end if
            call GetWaterTemperature(locIndx, -pos, temp)
            call GetBuoyantVelocity(locIndx, -pos, rIndx, wb)
            gama = CalcSurfaceTension(temp)
            pressure = m_surfData%pressure - Roul*G*pos + 2*gama/rr
            call CalcGasRatioInSingleBubble(Vbg(:,rIndx), ratio)
            Cb = GetBubbleAmount(Vbg(:,rIndx), rr, temp, pressure)
            if (Cb<e8) then
               exit
            end if
            call CalcGasExchangeRate(locIndx, -pos, rIndx, ratio, Ndm) 
            call GetRadiusChangeRates(locIndx, -pos, rIndx, tmp1, tmp2)
            Vbg_old = Vbg(:,rIndx)
            Vbg(:,rIndx) = Vbg_old + Cb*Ndm*tdelta
            yn = minval(Vbg(:,rIndx))
            if (yn<-100*Bubtol) then
               rr_old = rr
               exit
            else
               rr_old = rr
               rr = rr + (1.0e-6*tmp1*sum(Ndm)+tmp2)*tdelta
            end if
            pos = pos + wb*tdelta
         end do
         if (rr>e8 .and. Cb>e8 .and. locIndx==top) then
            Vr(rIndx,top) = rr_old
            Vex(:,rIndx,top) = Cb*Ndm
            m_bubbleGasCon(:,rIndx,top) = Vbg_old 
         else
            Vr(rIndx,top:locIndx) = rr_old
            Vex(:,rIndx,top:locIndx) = 0.0_r8
            m_bubbleGasCon(:,rIndx,top:locIndx) = 0.0_r8
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: some utility functions
   !
   !------------------------------------------------------------------------------ 
   subroutine InitializeBubbleStateVariables()
      implicit none
      real(r8) :: minR, maxR, rratio, dr
      integer :: ii

      minR = 2.5_r8     ! mm
      maxR = 10.0_r8    ! mm
      dr = (maxR - minR) / dble(NRLAYER)
      m_Rb0 = (/(minR+dr*(ii-1), ii = 1, NRLAYER+1)/)
      m_Rb0 = 1.0e-3_r8 * m_Rb0  ! mm => m
      rratio = minR / 0.1_r8
      tdelta = 0.4_r8 / (0.1418*rratio**2 + 0.05579*rratio + 0.7794)
      m_bubbleGasCon = 0.0_r8
      m_iceBubblePool = 0.0_r8
   end subroutine

   ! bubbles of different radius have the equal surface tension (surface areas)
   subroutine GetBubbleReleaseMagnitude(con)
      implicit none
      real(r8), intent(out) :: con(NGAS,NRLAYER+1)    ! unit: umol/(m3*mm)
      real(r8) :: temp, gama, pressure, wb
      real(r8) :: rr, Atmp 
      real(r8) :: vsc
      integer :: rIndx

      if (m_lakeWaterTopIndex>WATER_LAYER+1) then
         con = 0.0_r8
         return
      end if

      temp = m_waterTemp(WATER_LAYER+1)
      vsc = m_dVsc(WATER_LAYER+1) / Roul
      gama = CalcSurfaceTension(temp)
      pressure = m_surfData%pressure + Roul*G*lake_info%depth
      Atmp = 0._r8
      do rIndx = 1, NRLAYER+1, 1
         rr = m_Rb0(rIndx)
         wb = CalcBuoyantVelocity(rr, vsc)
         Atmp = Atmp + (pressure*rr + 2*gama) * wb
         con(:,rindx) = m_btmbflux * (pressure*rr + 2*gama)
      end do
      con(:,rindx) = con(:,rindx) / Atmp
   end subroutine

   ! Initialize the bubble radius field in water column without considering dissolution
   subroutine InitBubbleRadiusField()
      implicit none
      real(r8) :: temp, vsc, wb, rr, rr_old
      real(r8) :: gama, pos, denom, numer
      real(r8) :: ice_water_interface
      integer :: rIndx, locIndx, top

      do locIndx = 1, WATER_LAYER+1, 1
         Vr(:,locIndx) = m_Rb0
      end do
      if (m_lakeWaterTopIndex>WATER_LAYER+1) then
         return
      end if
      top = m_lakeWaterTopIndex
      ice_water_interface = -m_Zw(top)
      do rIndx = 1, NRLAYER+1, 1
         ! bubble release from sediment-water interface
         locIndx = WATER_LAYER + 1
         pos = -lake_info%depth
         rr = m_Rb0(rIndx)
         rr_old = rr
         do while (pos+e8<ice_water_interface .and. rr>e8)
            if (pos>-m_Zw(locIndx)) then
               Vr(rIndx,locIndx) = rr_old
               locIndx = locIndx - 1
            end if
            call GetWaterTemperature(locIndx, -pos, temp)
            call GetWaterViscosity(locIndx, -pos, vsc)
            wb = CalcBuoyantVelocity(rr, vsc/Roul)
            gama = CalcSurfaceTension(temp)
            denom = 3*m_surfData%pressure - 3*Roul*G*pos + 4*gama/rr
            numer = rr*Roul*G*wb
            rr_old = rr
            rr = rr + numer/denom*tdelta
            pos = pos + wb*tdelta
         end do
         if (rr>e8 .and. locIndx==top) then
            Vr(rIndx,locIndx) = rr_old
         else
            Vr(rIndx,top:locIndx) = rr_old
         end if
      end do
   end subroutine

   subroutine GetWaterTemperature(idx, pos, temp)
      implicit none
      integer, intent(in) :: idx
      real(r8), intent(in) :: pos
      real(r8), intent(out) :: temp
      real(r8) :: tt

      if (idx==WATER_LAYER+1) then
         temp = m_waterTemp(idx)
      else
         tt = (pos - m_Zw(idx)) / (m_Zw(idx+1) - m_Zw(idx))
         temp = (1.0-tt)*m_waterTemp(idx) + tt*m_waterTemp(idx+1)
      end if
      temp = max(temp, T0)
   end subroutine

   subroutine GetWaterViscosity(idx, pos, vsc)
      implicit none
      integer, intent(in) :: idx
      real(r8), intent(in) :: pos
      real(r8), intent(out) :: vsc
      real(r8) :: tt

      if (idx==WATER_LAYER+1) then
         vsc = m_dVsc(idx)
      else
         tt = (pos - m_Zw(idx)) / (m_Zw(idx+1) - m_Zw(idx))
         vsc = (1.0-tt)*m_dVsc(idx) + tt*m_dVsc(idx+1)
      end if
   end subroutine

   subroutine GetWaterGasConcentrations(idx, pos, con)
      implicit none
      integer, intent(in) :: idx
      real(r8), intent(in) :: pos
      real(r8), intent(out) :: con(NGAS)
      real(r8) :: tt

      if (idx==WATER_LAYER+1) then
         con = m_waterSubCon(:,idx)
      else
         tt = (pos - m_Zw(idx)) / (m_Zw(idx+1) - m_Zw(idx))
         con = (1.0-tt)*m_waterSubCon(:,idx) + tt*m_waterSubCon(:,idx+1)
      end if
   end subroutine

   subroutine GetBuoyantVelocity(idx, pos, rIndx, wb)
      implicit none
      integer, intent(in) :: idx
      real(r8), intent(in) :: pos
      integer, intent(in) :: rIndx
      real(r8), intent(out) :: wb
      real(r8) :: tt

      if (idx==WATER_LAYER+1) then
         wb = Vab(rIndx,idx)
      else
         tt = (pos - m_Zw(idx)) / (m_Zw(idx+1) - m_Zw(idx))
         wb = (1.0-tt)*Vab(rIndx,idx) + tt*Vab(rIndx,idx+1)
      end if
   end subroutine

   subroutine CalcGasExchangeRate(idx, pos, rIndx, ratio, exchange)
      implicit none
      integer, intent(in) :: idx
      real(r8), intent(in) :: pos
      integer, intent(in) :: rIndx
      real(r8), intent(in) :: ratio(NGAS)
      real(r8), intent(out) :: exchange(NGAS)
      real(r8) :: tmp1(NGAS), tmp2(NGAS), tt

      if (idx==WATER_LAYER+1) then
         exchange = Vag1(:,rIndx,idx)*ratio + Vag2(:,rIndx,idx)
      else
         tt = (pos - m_Zw(idx)) / (m_Zw(idx+1) - m_Zw(idx))
         tmp1 = (1.0-tt)*Vag1(:,rIndx,idx) + tt*Vag1(:,rIndx,idx+1)
         tmp2 = (1.0-tt)*Vag2(:,rIndx,idx) + tt*Vag2(:,rIndx,idx+1)
         exchange = tmp1*ratio + tmp2
      end if
      exchange(Wn2:Wo2) = min(0.0_r8, exchange(Wn2:Wo2))
   end subroutine

   subroutine GetRadiusChangeRates(idx, pos, rIndx, rate1, rate2)
      implicit none
      integer, intent(in) :: idx
      real(r8), intent(in) :: pos
      integer, intent(in) :: rIndx
      real(r8), intent(out) :: rate1
      real(r8), intent(out) :: rate2
      real(r8) :: tt

      if (idx==WATER_LAYER+1) then
         rate1 = Var1(rIndx,idx)
         rate2 = Var2(rIndx,idx)
      else
         tt = (pos - m_Zw(idx)) / (m_Zw(idx+1) - m_Zw(idx))
         rate1 = (1.0-tt)*Var1(rIndx,idx) + tt*Var1(rIndx,idx+1)
         rate2 = (1.0-tt)*Var2(rIndx,idx) + tt*Var2(rIndx,idx+1)
      end if
   end subroutine

   subroutine CalcGasRatioInSingleBubble(bubbleGasCon, ratio)
      implicit none
      real(r8), intent(in) :: bubbleGasCon(NGAS)
      real(r8), intent(out) :: ratio(NGAS)
      real(r8) :: gases

      gases = sum(bubbleGasCon)
      if (abs(gases)<e8) then
         ratio = (/0.0, 0.0, 0.0, 0.0/)
      else
         ratio = bubbleGasCon / abs(gases)
      end if
   end subroutine

   function GetBubbleAmount(bubbleGasCon, radius, temp, pressure)
      implicit none
      real(r8), intent(in) :: bubbleGasCon(NGAS)
      real(r8), intent(in) :: radius
      real(r8), intent(in) :: temp
      real(r8), intent(in) :: pressure
      real(r8) :: GetBubbleAmount   ! unit: bubble/(m3*mm)
      real(r8) :: gases

      gases = sum(bubbleGasCon) 
      GetBubbleAmount = 1.0e-6_r8*gases*0.75*R*temp/ &
                  (Pi*pressure*(radius**3))
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: update gas flux via ebullition and bubble gas pool in winter
   !
   !------------------------------------------------------------------------------
   subroutine UpdateBubbleFlux(dt)
      implicit none
      real(r8), intent(in) :: dt
      real(r8) :: flux(NGAS), cb(NGAS)
      real(r8) :: rr, vsc, wb
      integer :: ii, top

      bubble_fluxes = 0.0_r8
      if (m_lakeWaterTopIndex>WATER_LAYER+1) then
         return
      end if

      top = m_lakeWaterTopIndex
      vsc = m_dVsc(top) / Roul
      flux = 0.0_r8
      do ii = 1, NRLAYER+1, 1
         cb = m_bubbleGasCon(:,ii,top)
         rr = Vr(ii,top)
         wb = CalcBuoyantVelocity(rr, vsc)
         flux = flux + cb * wb 
      end do
      if ( (lake_info%thrmkst==2 .and. lake_info%margin==1) .or. &
            (m_Hice<e8) ) then
         bubble_fluxes = flux + Blr * m_iceBubblePool
         m_iceBubblePool = m_iceBubblePool * (1.0 - Blr*dt)
         do ii = 1, NGAS, 1
            m_iceBubblePool(ii) = max(0d0, m_iceBubblePool(ii))
         end do
      else if (m_Hice>e8) then
         ! only 10% of bubble remains in gaseous state after winter
         ! Greene et al. (2014)
         bubble_fluxes = 0.0_r8
         m_iceBubblePool = m_iceBubblePool + 0.1 * flux * dt
      end if
   end subroutine
   
end module bubble_mod
