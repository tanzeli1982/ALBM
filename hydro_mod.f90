module hydro_mod
!---------------------------------------------------------------------------------
! Purpose: Govern the change of lake states due to hydrologic fluxes.
!
!--------------------------------------------------------------------------------- 
   use shr_kind_mod,          only : r8
   use shr_ctrl_mod,          only : WATER_LAYER, SED_LAYER, NSCOL, NZWD
   use shr_ctrl_mod,          only : lake_info
   use shr_ctrl_mod,          only : Thermal_Module, Carbon_Module
   use phy_utilities_mod
   use data_buffer_mod

   implicit none
   private
   public :: InitializeHydroModule, DestructHydroModule
   public :: UpdateLakeStatesForHydroFluxes
   public :: GetOutflowFluxes
   ! withdrawal weights, supply, and real rate
   real(r8), allocatable :: weightsWD(:,:)
   real(r8), allocatable :: supplyWD(:)      ! m3
   real(r8), allocatable :: realWD(:)        ! m3
   ! layer volume
   real(r8), allocatable :: zVol(:)          ! m3
   real(r8), allocatable :: zVol0(:)         ! m3
   ! temporal state variables
   real(r8), allocatable :: tmpwaterTemp(:)
   real(r8), allocatable :: tmpwaterSubCon(:,:)

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: Private data initialization, get/set and destruction
   !
   !------------------------------------------------------------------------------
   subroutine InitializeHydroModule()
      implicit none

      allocate(weightsWD(NZWD,WATER_LAYER+1))
      allocate(supplyWD(WATER_LAYER+1))
      allocate(realWD(WATER_LAYER+1))
      allocate(zVol(WATER_LAYER+1))
      allocate(zVol0(WATER_LAYER+1))
      allocate(tmpwaterTemp(WATER_LAYER+1))
      allocate(tmpwaterSubCon(NWSUB,WATER_LAYER+1))
   end subroutine

   subroutine DestructHydroModule()
      implicit none

      deallocate(weightsWD)
      deallocate(supplyWD)
      deallocate(realWD)
      deallocate(zVol)
      deallocate(zVol0)
      deallocate(tmpwaterTemp)
      deallocate(tmpwaterSubCon)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Update lake structure and thermal & chemical states.
   !
   !------------------------------------------------------------------------------
   subroutine UpdateLakeStatesForHydroFluxes(dt)
      implicit none
      real(r8), intent(in) :: dt
      real(r8) :: dV, Vtot, Vtot0 
      real(r8) :: over_request
      real(r8) :: zw, za, zb, rw
      real(r8) :: overlap, sumw
      real(r8) :: rdiff, maxZ
      real(r8) :: Vres, Vnew
      real(r8), parameter :: rw_min = 0.2_r8
      real(r8), parameter :: dZw_min = 0.1_r8
      real(r8), parameter :: err_max = 1.0_r8
      integer :: ii, jj, indx, top
      logical :: to_regrid

      ! skip ice-on period
      if (m_Hice>e8) then
         return
      end if

      top = m_lakeWaterTopIndex

      ! increase surface layer volume by inflow
      zVol = m_dZw * m_Az        ! current volume
      zVol0 = zVol
      dV = m_hydroData%Qwi * dt
      zVol(top) = zVol(top) + dV
      
      ! update state variables 
      !if (Thermal_Module) then
      !   m_waterTemp(top) = (m_waterTemp(top)*zVol0(top) + &
      !      m_hydroData%WTi*dV) / zVol(top)
      !end if
      !if (Carbon_Module) then
      !   m_waterSubCon(Wo2,top) = (m_waterSubCon(Wo2,top)*zVol0(top) + &
      !      m_hydroData%o2i*dV) / zVol(top)
      !   m_waterSubCon(Wco2,top) = (m_waterSubCon(Wco2,top)*zVol0(top) + &
      !      m_hydroData%co2i*dV) / zVol(top)
      !   m_waterSubCon(Wch4,top) = (m_waterSubCon(Wch4,top)*zVol0(top) + &
      !      m_hydroData%ch4i*dV) / zVol(top)
      !   m_waterSubCon(Wsrp,top) = (m_waterSubCon(Wsrp,top)*zVol0(top) + &
      !      m_hydroData%srpi*dV) / zVol(top)
      !end if

      do ii = 1, WATER_LAYER+1, 1
         supplyWD(ii) = max(zVol(ii)-dZw_min*m_Az(ii), 0._r8)
      end do

      ! calculate withdrawal weights
      do ii = 1, NZWD, 1
         zw = m_Zw(WATER_LAYER+1) - m_hydroData%zWD(ii)
         indx = COUNT(m_Zw<=zw)
         if (indx==0) then
            weightsWD(ii,:) = 0._r8
            cycle
         end if
         if (indx<WATER_LAYER+1) then
            rw = max( 0.5*(m_dZw(indx)+m_dZw(indx+1)), rw_min )
         else
            rw = max( m_dZw(indx), rw_min )
         end if
         za = zw - rw
         zb = zw + rw
         do jj = 1, WATER_LAYER+1, 1
            if (jj==1) then
               overlap = min(m_Zw(jj)+m_dZw(jj), zb) - max(m_Zw(jj), za)
            else if (jj==WATER_LAYER+1) then
               overlap = min(m_Zw(jj), zb) - max(m_Zw(jj)-m_dZw(jj), za)
            else
               overlap = min(m_Zw(jj)+0.5*m_dZw(jj), zb) - &
                  max(m_Zw(jj)-0.5*m_dZw(jj), za)
            end if
            if (overlap>e8) then
               weightsWD(ii,jj) = overlap
            else
               weightsWD(ii,jj) = 0._r8
            end if
         end do
         sumw = sum(weightsWD(ii,:))
         if (sumw > e8) then
            weightsWD(ii,:) = weightsWD(ii,:) / sumw
         else
            call Endrun('Zero sum of withdrawal weights')
         end if
      end do 

      ! calculate requested outflow
      realWD = 0._r8
      do ii = 1, NZWD, 1
         realWD = realWD + m_hydroData%Qwo(ii) * dt * weightsWD(ii,:)
      end do

      ! cap and redistribute over-requested outflow 
      over_request = 0._r8
      do ii = 1, WATER_LAYER+1, 1
         if (supplyWD(ii)>=realWD(ii) .and. over_request>e8) then
            if (supplyWD(ii)>=realWD(ii)+over_request) then
               realWD(ii) = realWD(ii) + over_request
               over_request = 0._r8
            else
               over_request = over_request + realWD(ii) - supplyWD(ii)
               realWD(ii) = supplyWD(ii)
            end if
         else if (supplyWD(ii)<realWD(ii)) then
            over_request = over_request + realWD(ii) - supplyWD(ii)
            realWD(ii) = supplyWD(ii)
         end if
      end do
      do ii = WATER_LAYER+1, 1, -1
         if (supplyWD(ii)>=realWD(ii) .and. over_request>e8) then
            if (supplyWD(ii)>=realWD(ii)+over_request) then
               realWD(ii) = realWD(ii) + over_request
               over_request = 0._r8
            else
               over_request = over_request + realWD(ii) - supplyWD(ii)
               realWD(ii) = supplyWD(ii)
            end if
         else if (supplyWD(ii)<realWD(ii)) then
            over_request = over_request + realWD(ii) - supplyWD(ii)
            realWD(ii) = supplyWD(ii)
         end if
      end do

      if (over_request>err_max) then
         call Endrun('Requested withdrawal exceeds available water')
      end if

      ! decrease layer volume by outflow
      zVol = zVol - realWD

      ! update layer thickness
      m_dZw = zVol / m_Az
      do ii = WATER_LAYER, 1, -1
         if (ii==WATER_LAYER) then
            m_Zw(ii) = m_Zw(ii+1) - m_dZw(ii+1) - 0.5*m_dZw(ii)
         else if (ii==1) then
            m_Zw(ii) = m_Zw(ii+1) - 0.5*m_dZw(ii+1) - m_dZw(ii)
         else
            m_Zw(ii) = m_Zw(ii+1) - 0.5*m_dZw(ii+1) - 0.5*m_dZw(ii)
         end if
      end do

      ! check whether regridding is needed
      to_regrid = .False.
      do ii = 1, WATER_LAYER+1, 1
         rdiff = zVol(ii)/zVol0(ii) - 1._r8
         if (abs(rdiff)>0.01_r8) then
            to_regrid = .True.
            exit
         end if
      end do 

      ! regridding is not needed
      if (.NOT. to_regrid) then
         return
      end if

      ! adjust lake grid 
      Vtot0 = sum(m_dZw*m_Az)
      maxZ = sum(m_dZw)
      !print *, maxZ, lake_info%maxdepth, m_Zw(1), m_Zw(WATER_LAYER+1)

      call BuildDepthVector(maxZ, lake_info%hsed, m_Zw, m_dZw, m_Zs, m_dZs)
      m_Zw = m_Zw + (lake_info%maxdepth - maxZ)
      call BuildLakeStructure(lake_info, m_hypsocurve, Vtot0, m_Zw, m_dZw, &
               m_Az, m_dAz, m_SAL, m_soilColZ, m_soilColInd)
      Vtot = sum(m_dZw*m_Az)

      if (abs(Vtot - Vtot0) > 1.d-3) then
         print *, "New vs. old lake volume: ", Vtot, Vtot0
         call Endrun('Lake volume does not conserve')
      end if

      ! remap state variables onto the new grids
      indx = 1
      zVol0 = zVol
      do ii = 1, WATER_LAYER+1, 1
         Vnew = m_dZw(ii) * m_Az(ii)
         Vres = Vnew
         tmpwaterTemp(ii) = 0._r8
         tmpwaterSubCon(:,ii) = 0._r8
         do jj = indx, WATER_LAYER+1, 1
            if (zVol0(jj)>=Vres) then
               tmpwaterTemp(ii) = tmpwaterTemp(ii) + Vres*m_waterTemp(jj)
               tmpwaterSubCon(:,ii) = tmpwaterSubCon(:,ii) + &
                  Vres*m_waterSubCon(:,jj)
               zVol0(jj) = zVol0(jj) - Vres
               indx = jj
               exit
            else
               tmpwaterTemp(ii) = tmpwaterTemp(ii) + zVol0(jj)*m_waterTemp(jj)
               tmpwaterSubCon(:,ii) = tmpwaterSubCon(:,ii) + &
                  zVol0(jj)*m_waterSubCon(:,jj)
               Vres = Vres - zVol0(jj)
               zVol0(jj) = 0._r8
            end if
         end do
         ! average
         tmpwaterTemp(ii) = tmpwaterTemp(ii) / Vnew
         tmpwaterSubCon(:,ii) = tmpwaterSubCon(:,ii) / Vnew
      end do

      ! update state variables
      m_waterTemp = tmpwaterTemp
      m_waterSubCon = tmpwaterSubCon

   end subroutine 

   !------------------------------------------------------------------------------
   !
   ! Purpose: Get outflow thermal and chemical fluxes.
   !
   !------------------------------------------------------------------------------
   subroutine GetOutflowFluxes(Qwt, Qo2, Qco2, Qch4, Qsrp)
      implicit none
      real(r8), intent(out) :: Qwt     ! K
      real(r8), intent(out) :: Qo2     ! mol/m3
      real(r8), intent(out) :: Qco2    ! mol/m3
      real(r8), intent(out) :: Qch4    ! mol/m3
      real(r8), intent(out) :: Qsrp    ! mol/m3
      real(r8) :: Qwo

      Qwo = sum(realWD)    ! m3 
      if (Qwo > e8) then
         Qwt = sum( realWD * m_waterTemp ) / Qwo 
         Qo2 = sum( realWD * m_waterSubCon(Wo2,:) ) / Qwo
         Qco2 = sum( realWD * m_waterSubCon(Wco2,:) ) / Qwo
         Qch4 = sum( realWD * m_waterSubCon(Wch4,:) ) / Qwo
         Qsrp = sum( realWD * m_waterSubCon(Wsrp,:) ) / Qwo
      else
         Qwt = -9999._r8
         Qo2 = -9999._r8 
         Qco2 = -9999._r8
         Qch4 = -9999._r8
         Qsrp = -9999._r8
      end if
   end subroutine

end module hydro_mod
