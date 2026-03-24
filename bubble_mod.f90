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
   public :: InitializeBubbleModule, DestructBubbleModule 
   public :: BubbleModuleSetup, BubbleModuleCallback
   public :: BubbleDynamics
   ! bubble gas concentration vector (mol/m3)
   real(r8), allocatable :: Vbg0(:,:,:)
   real(r8), allocatable :: Cbg(:,:,:) 
   ! buoyant velocity auxiliary vector
   real(r8), allocatable :: Vab(:,:)
   ! kinematic viscosity (m^2/s)
   real(r8), allocatable :: dVsc(:)
   ! bubble surface tension (N/m)
   real(r8), allocatable :: bGama(:)
   ! bubble gas solubility (mol/m3/Pa)
   real(r8), allocatable :: bSolu(:,:)
   ! bubble gas diffusivity (m^2/s)
   real(r8), allocatable :: bDiff(:,:)
   ! schmidt number (m^2/s)
   real(r8), allocatable :: schmidt(:,:)
   ! bubble amount
   real(r8), allocatable :: bNumb(:,:)
   ! water layer boundaries (m)
   real(r8), allocatable :: topZw(:)
   real(r8), allocatable :: btmZw(:)
   ! time step (s)
   real(r8) :: tdelta

contains
   subroutine InitializeBubbleModule()
      implicit none     

      allocate(Vbg0(NGAS,NRLAYER+1,NSCOL))
      allocate(Cbg(NGAS,NRLAYER+1,WATER_LAYER+2))
      allocate(Vab(NRLAYER+1,NSCOL))
      allocate(dVsc(WATER_LAYER+1))
      allocate(bGama(WATER_LAYER+1))
      allocate(bSolu(NGAS,WATER_LAYER+1))
      allocate(bDiff(NGAS,WATER_LAYER+1))
      allocate(schmidt(NGAS,WATER_LAYER+1))
      allocate(bNumb(NRLAYER+1,NSCOL))
      allocate(topZw(WATER_LAYER+1))
      allocate(btmZw(WATER_LAYER+1))
      
      call InitializeBubbleStateVariables()
      Vbg0 = 0.0_r8
      Cbg = 0.0_r8
      Vab = 0.0_r8
   end subroutine

   subroutine DestructBubbleModule()
      implicit none

      deallocate(Vbg0)
      deallocate(Cbg)
      deallocate(Vab)
      deallocate(dVsc)
      deallocate(bGama)
      deallocate(bSolu)
      deallocate(bDiff)
      deallocate(schmidt)
      deallocate(topZw)
      deallocate(btmZw)
      deallocate(bNumb)
   end subroutine

   subroutine BubbleModuleSetup(isUpdate)
      implicit none
      logical, intent(in) :: isUpdate
      real(r8) :: temp, pH, vsc
      integer :: ii, kk

      if (.NOT. isUpdate) then
         return
      end if

      pH = lake_info%pH
      do ii = 1, WATER_LAYER+1, 1
         temp = max(m_waterTemp(ii), T0)
         vsc = CalcDynamicViscosity(m_waterTemp(ii))
         dVsc(ii) = vsc / Roul
         bGama(ii) = CalcSurfaceTension(temp)
         do kk = Wn2, Wch4, 1
            bSolu(kk,ii) = CalcHenrySolubility(kk,temp,pH)
            bDiff(kk,ii) = CalcGasDiffusivityInWater(kk,temp)
            schmidt(kk,ii) = CalcSchmidtNumber(kk,temp)
         end do
      end do

      call UpdateBubbleReleaseMagnitude()
   end subroutine

   subroutine BubbleModuleCallback(dt)
      implicit none
      real(r8), intent(in) :: dt
      real(r8) :: icebflux, icebloss
      real(r8) :: sumAz, topAz
      integer :: ii, kk, icol, top
      integer :: indx, indx0

      top = m_lakeWaterTopIndex
      topAz = 0._r8
      do icol = 1, NSCOL, 1
         indx = COUNT(m_soilColInd<=icol)
         indx0 = COUNT(m_soilColInd<=icol-1) + 1
         !if (indx0<=indx .and. m_soilColBub(icol)>0) then
         if (indx0<=indx) then
            sumAz = sum(m_dAz(indx0:indx))
            topAz = topAz + sumAz
            if (indx<top) then
               ! effective emitting area
               m_iceBubblePool = m_iceBubblePool + m_btmbflux(:,icol)* &
                  sumAz*dt
            end if
         end if
      end do

      ! bubble trapping and release
      if (m_Hice>e8) then
         m_iceBubblePool = m_iceBubblePool + m_topbflux*topAz*dt
         m_upbflux = m_topbflux(Wch4)
         m_topbflux = 0.0_r8
      else
         icebflux = sa_params(Param_icebflux)
         m_upbflux = m_topbflux(Wch4)
         m_topbflux = m_topbflux + m_iceBubblePool*icebflux/topAz
         do kk = Wn2, Wch4, 1
            m_iceBubblePool(kk) = m_iceBubblePool(kk) * max(0.0_r8, &
               (1.0-icebflux*dt))
         end do
      end if

      ! CH4 leaking from ice bubbles
      icebloss = sa_params(Param_icebloss)
      if (m_iceBubblePool(Wch4)>e8 .and. top<=WATER_LAYER+1) then
         m_iceBubblePool(Wch4) = m_iceBubblePool(Wch4) * max(0.0_r8, &
            (1.0-icebloss*dt))
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Discretize single bubble evolution by Finite Difference Method.
   !          This sub-routine returns time derivatives of bubble gas concentration.
   !          Four gases are included: N2, O2, CO2 and CH4.
   !          
   !
   !------------------------------------------------------------------------------
   subroutine BubbleDynamics(isUpdate)
      implicit none
      logical, intent(in) :: isUpdate
      real(r8) :: ratio(NGAS), Vbg(NGAS), exlayer(NGAS)
      real(r8) :: k_ndm(NGAS), ndm(NGAS)
      real(r8) :: peclet(NGAS), reynold(NGAS), nusselt(NGAS)
      real(r8) :: pressure, gama, pos, pos_frz
      real(r8) :: wb, rr, temp, icebloss, ex_iceb 
      real(r8) :: pr_tmp, dr_tmp1, dr_tmp2
      real(r8) :: sumAz, topAz
      integer :: rIndx, locIndx, indx, indx0
      integer :: top, ii, jj, icol

      if (.NOT. isUpdate) then
         return
      end if

      top = m_lakeWaterTopIndex
      topAz = 0._r8
      do icol = 1, NSCOL, 1
         indx = COUNT(m_soilColInd<=icol)
         indx0 = COUNT(m_soilColInd<=icol-1) + 1
         !if (indx0<=indx .and. m_soilColBub(icol)>0) then
         if (indx0<=indx) then
            sumAz = sum(m_dAz(indx0:indx))
            topAz = topAz + sumAz
         end if
      end do

      m_gasExchange = 0.0_r8
      m_topbflux = 0.0_r8
      do icol = 1, NSCOL, 1
         indx = COUNT(m_soilColInd<=icol)
         indx0 = COUNT(m_soilColInd<=icol-1) + 1
         ! frozen water column or inactive sediment column
         !if (indx<top .or. indx0>indx .or. m_soilColBub(icol)==0) then
         if (indx<top .or. indx0>indx) then
            cycle
         end if
         ! effective emitting area
         sumAz = sum(m_dAz(indx0:indx))

         ! too small sediment flux
         if (m_btmbflux(Wch4,icol)<1.d-12) then
            do jj = indx0, indx, 1
               m_gasExchange(Wch4,jj) = m_gasExchange(Wch4,jj) + &
                  m_btmbflux(Wch4,icol)*m_dAz(jj)/m_Az(jj)/m_dZw(jj)
            end do
            cycle
         end if

         m_topbflux = m_topbflux + m_btmbflux(:,icol)*sumAz/topAz

         Cbg(:,:,indx+2:WATER_LAYER+2) = 0.0_r8 
         Cbg(:,:,indx+1) = Vbg0(:,:,icol)
         Cbg(:,:,1:indx) = 0.0_r8

         pos_frz = -topZw(top)
         do rIndx = 1, NRLAYER+1, 1
            ! bubble release from sediment-water interface
            pos = -btmZw(indx)
            rr = m_Rb0(rIndx) 
            Vbg = Vbg0(:,rIndx,icol)
            locIndx = indx
            do while (pos<pos_frz .and. rr>e8)
               if (pos>-topZw(locIndx)) then
                  Cbg(:,rIndx,locIndx) = Vbg
                  locIndx = locIndx - 1
               end if
               temp = m_waterTemp(locIndx) 
               wb = CalcBuoyantVelocity(rr, dVsc(locIndx))
               gama = bGama(locIndx) 
               pressure = m_surfData%pressure + 2.0*gama/rr - &
                  Roul*G*(pos+m_surfData%dzsurf)
               ! gas exchange rate
               peclet = rr * wb / bDiff(:,locindx)
               reynold = peclet / schmidt(:,locindx)
               where (reynold<=1.0)
                  nusselt = sqrt(2.0*Pi*peclet/3.0)
               elsewhere
                  nusselt = 0.45*(reynold**(1.0/6.0))*(peclet**(1.0/3.0))
               end where
               k_ndm = -4.0 * Pi * rr * bDiff(:,locindx) * nusselt
               call CalcGasRatioInSingleBubble(Vbg, ratio)
               ndm = k_ndm * (bSolu(:,locindx)*pressure*ratio - &
                  m_waterSubCon(1:NGAS,locindx))   ! gas exchange rate
               ! new bubble gas
               Vbg = Vbg + bNumb(rIndx,icol) * ndm * tdelta
               where (Vbg<0._r8) Vbg = 0._r8    ! remove small negatives
               ! new radius
               pr_tmp = 3.0*m_surfData%pressure - 3.0*Roul*G*(pos+ &
                  m_surfData%dzsurf) + 4.0*gama/rr
               dr_tmp1 = 0.75d-3*R*temp/(Pi*rr**2.0)/pr_tmp
               dr_tmp2 = rr*Roul*G*wb/pr_tmp
               rr = rr + (dr_tmp1*sum(ndm)+dr_tmp2)*tdelta
               ! new position
               pos = pos + wb*tdelta
            end do
            ! including ice layers
            do jj = 1, locIndx, 1
               Cbg(:,rIndx,jj) = Vbg
            end do
         end do

         ! update water-bubble gas exchange 
         do ii = indx, top, -1
            ! size-integrated bubble gas and gas exchange
            do rIndx = 1, NRLAYER+1, 1
               exlayer = (Cbg(:,rIndx,ii+1) - Cbg(:,rIndx,ii))*Vab(rIndx,icol)
               if (ii>=indx0) then
                  m_gasExchange(:,ii) = m_gasExchange(:,ii) + exlayer/m_dZw(ii)
               else
                  m_gasExchange(:,ii) = m_gasExchange(:,ii) + exlayer*sumAz/ &
                     m_Az(ii)/m_dZw(ii)
               end if
            end do
         end do
      end do

      ! update bubble flux
      do ii = WATER_LAYER+1, top, -1
         m_topbflux = m_topbflux - m_gasExchange(:,ii)*m_dZw(ii)* &
            m_Az(ii)/topAz
      end do

      ! CH4 leaking from ice bubbles
      icebloss = sa_params(Param_icebloss)
      if (m_iceBubblePool(Wch4)>e8 .and. top<=WATER_LAYER+1) then
         ex_iceb = m_iceBubblePool(Wch4) * icebloss / m_dZw(top) / m_Az(top)
         m_gasExchange(Wch4,top) = m_gasExchange(Wch4,top) + ex_iceb
      end if

   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: some utility functions
   !
   !------------------------------------------------------------------------------ 
   subroutine InitializeBubbleStateVariables()
      implicit none
      real(r8) :: minR, maxR, rratio, dr
      integer :: ii, kk, icol

      minR = 2.5_r8     ! mm
      maxR = 10.0_r8    ! mm
      dr = (maxR - minR) / dble(NRLAYER)
      m_Rb0 = (/(minR+dr*(ii-1), ii = 1, NRLAYER+1)/)
      m_Rb0 = 1.0e-3_r8 * m_Rb0  ! mm => m
      rratio = minR / 0.1_r8
      tdelta = 2.0_r8 / (0.1418*rratio**2 + 0.05579*rratio + 0.7794)
      m_iceBubblePool = 0.0_r8
      m_topbflux = 0.0_r8
      m_upbflux = 0.0_r8

      ! initialize intermediate variables
      do ii = 1, WATER_LAYER+1, 1
         if (ii==1) then
            topZw(ii) = m_Zw(ii)
            btmZw(ii) = m_Zw(ii) + m_dZw(ii)
         else if (ii==WATER_LAYER+1) then
            topZw(ii) = m_Zw(ii) - m_dZw(ii)
            btmZw(ii) = m_Zw(ii)
         else
            topZw(ii) = m_Zw(ii) - 0.5*m_dZw(ii)
            btmZw(ii) = m_Zw(ii) + 0.5*m_dZw(ii)
         end if
      end do
   end subroutine

   ! bubbles of different radius have the equal surface tension (surface areas)
   subroutine UpdateBubbleReleaseMagnitude()
      implicit none
      real(r8) :: temp, pressure, wb
      real(r8) :: rr, tmp1, tmp2
      integer :: rIndx, indx, indx0
      integer :: top, icol

      top = m_lakeWaterTopIndex
      do icol = 1, NSCOL, 1
         indx = COUNT(m_soilColInd<=icol)
         indx0 = COUNT(m_soilColInd<=icol-1) + 1

         ! frozen water column or inactive sediment column
         !if (indx<top .or. indx0>indx .or. m_soilColBub(icol)==0) then
         if (indx<top .or. indx0>indx) then
            Vbg0(:,:,icol) = 0.0_r8
            bNumb(:,icol) = 0.0_r8
            cycle
         end if

         pressure = m_surfData%pressure + Roul*G*(btmZw(indx)- &
            m_surfData%dzsurf)
         tmp1 = 0.0_r8
         do rIndx = 1, NRLAYER+1, 1
            rr = m_Rb0(rIndx)
            wb = CalcBuoyantVelocity(rr, dVsc(indx))
            Vab(rIndx,icol) = wb
            tmp1 = tmp1 + (pressure*rr+2.0*bGama(indx))*wb
            Vbg0(:,rIndx,icol) = m_btmbflux(:,icol)*(pressure*rr+ &
               2.0*bGama(indx))
         end do
         Vbg0(:,:,icol) = Vbg0(:,:,icol) / (tmp1+inft)

         temp = m_waterTemp(indx)
         do rIndx = 1, NRLAYER+1, 1
            rr = m_Rb0(rIndx)
            pressure = m_surfData%pressure + Roul*G*(btmZw(indx)- &
               m_surfData%dzsurf) + 2.0*bGama(indx)/rr
            bNumb(rIndx,icol) = 0.75d-3*R*temp/(Pi*pressure*rr**3.0)* &
               sum(Vbg0(:,rIndx,icol))
         end do
      end do
   end subroutine
   
   subroutine CalcGasRatioInSingleBubble(bubbleGasCon, ratio)
      implicit none
      real(r8), intent(in) :: bubbleGasCon(NGAS)
      real(r8), intent(out) :: ratio(NGAS)
      real(r8) :: mgas_tot

      mgas_tot = abs(sum(bubbleGasCon))
      if (mgas_tot<1.d-12) then
         ratio = (/0.0, 0.0, 0.0, 1.0/)
      else
         ratio = bubbleGasCon / mgas_tot
      end if
   end subroutine

end module bubble_mod
