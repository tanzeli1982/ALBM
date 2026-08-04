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
   use thermal_mod,           only : EnforceThermalConsistency
   use carbon_cycle_mod,      only : EnforceBGCConsistency 

   implicit none
   private
   public :: InitializeHydroModule, DestructHydroModule
   public :: UpdateLakeStatesForHydroFluxes
   public :: GetWithdrawFlow
   ! withdrawal weights, supply, and real rate
   real(r8), allocatable :: weightsWD(:,:)
   real(r8), allocatable :: supplyWD(:)      ! m3
   real(r8), allocatable :: realWD(:)        ! m3
   real(r8), allocatable :: ratioWD(:)
   ! layer volume
   real(r8), allocatable :: zVol(:)          ! m3
   real(r8), allocatable :: zVol0(:)         ! m3
   real(r8), allocatable :: accVol(:)        ! m3
   ! temporal state variables
   real(r8), allocatable :: tmpwaterHeat(:)
   real(r8), allocatable :: tmpWaterTemp(:)
   real(r8), allocatable :: tmpwaterIce(:)
   real(r8), allocatable :: tmpwaterSubCon(:,:)

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: Private data initialization, get/set and destruction
   !
   !------------------------------------------------------------------------------
   subroutine InitializeHydroModule()
      implicit none

      allocate(weightsWD(NZWD+1,WATER_LAYER+1))
      allocate(supplyWD(WATER_LAYER+1))
      allocate(realWD(WATER_LAYER+1))
      allocate(ratioWD(WATER_LAYER+1))
      allocate(zVol(WATER_LAYER+1))
      allocate(zVol0(WATER_LAYER+1))
      allocate(accVol(WATER_LAYER+1))
      allocate(tmpwaterHeat(WATER_LAYER+1))
      allocate(tmpWaterTemp(WATER_LAYER+1))
      allocate(tmpwaterIce(WATER_LAYER+1))
      allocate(tmpwaterSubCon(NWSUB,WATER_LAYER+1))

      zVol = m_dZw * m_Az        ! current volume
   end subroutine

   subroutine DestructHydroModule()
      implicit none

      deallocate(weightsWD)
      deallocate(supplyWD)
      deallocate(realWD)
      deallocate(ratioWD)
      deallocate(zVol)
      deallocate(zVol0)
      deallocate(accVol)
      deallocate(tmpwaterHeat)
      deallocate(tmpWaterTemp)
      deallocate(tmpwaterIce)
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
      real(r8) :: vw, va, vb
      real(r8) :: overlap, sumw
      real(r8) :: rw, rdiff, maxZ
      real(r8) :: Vres, Vnew, Ebg, Eed
      real(r8) :: Vsurf, Vhuman
      real(r8) :: tmpEres1, tmpEres2
      real(r8), parameter :: rw_min = 0.2_r8
      real(r8), parameter :: dZw_min = 0.1_r8
      real(r8), parameter :: err_max = 1.0_r8
      integer :: ii, jj, indx, top
      logical :: to_regrid

      top = m_lakeWaterTopIndex

      ! increase surface layer volume by inflow
      zVol0 = zVol
      Vtot0 = sum(zVol0)
      dV = m_hydroData%Qwi * dt
      zVol(top) = zVol(top) + dV
      
      ! update state variables 
      if (Thermal_Module) then
         if (m_hydroData%WTi>-e8) then
            m_waterTemp(top) = (m_waterTemp(top)*zVol0(top) + &
               m_hydroData%WTi*dV) / zVol(top)
         end if
      end if
      if (Carbon_Module) then
         m_waterSubCon(Wo2,top) = (m_waterSubCon(Wo2,top)*zVol0(top) + &
            m_hydroData%o2i*dV) / zVol(top)
         m_waterSubCon(Wco2,top) = (m_waterSubCon(Wco2,top)*zVol0(top) + &
            m_hydroData%co2i*dV) / zVol(top)
         m_waterSubCon(Wch4,top) = (m_waterSubCon(Wch4,top)*zVol0(top) + &
            m_hydroData%ch4i*dV) / zVol(top)
         m_waterSubCon(Wsrp,top) = (m_waterSubCon(Wsrp,top)*zVol0(top) + &
            m_hydroData%srpi*dV) / zVol(top)
      end if

      do ii = 1, WATER_LAYER+1, 1
         supplyWD(ii) = max(zVol(ii)-dZw_min*m_Az(ii), 0._r8)
      end do

      ! update cumulative volume 
      do jj = WATER_LAYER+1, 1, -1
         if (jj==WATER_LAYER+1) then
            accVol(jj) = zVol(jj)
         else
            accVol(jj) = accVol(jj+1) + zVol(jj)
         end if
      end do

      ! calculate withdrawal weights
      weightsWD = 0._r8
      do ii = 1, NZWD, 1
         vw = m_hydroData%vWD(ii)
         indx = COUNT(accVol>=vw)
         ! water surface is below the gate
         if (indx==0) then
            weightsWD(ii,1) = 1._r8
            cycle
         end if
         ! water surface is above the gate
         if (indx<WATER_LAYER+1) then
            rw = max( 0.5*(zVol(indx)+zVol(indx+1)), rw_min*m_Az(indx) )
         else
            rw = max( zVol(indx), rw_min*m_Az(indx) )
         end if
         va = vw + rw
         vb = vw - rw
         do jj = 1, WATER_LAYER+1, 1
            if (jj<WATER_LAYER+1) then
               overlap = min(accVol(jj), va) - max(accVol(jj+1), vb)
            else if (jj==WATER_LAYER+1) then
               overlap = min(accVol(jj), va) - vb
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
      weightsWD(NZWD+1,1) = 1._r8

      ! calculate requested outflow
      realWD = 0._r8
      do ii = 1, NZWD, 1
         realWD = realWD + m_hydroData%Qwo(ii) * dt * weightsWD(ii,:)
      end do
      realWD = realWD + m_hydroData%Qws * dt * weightsWD(NZWD+1,:) 

      do ii = 1, WATER_LAYER+1, 1
         Vsurf = m_hydroData%Qws * weightsWD(NZWD+1,ii)
         Vhuman = sum( m_hydroData%Qwo * weightsWD(:,ii) )
         if (Vsurf>e8) then
            ratioWD(ii) = Vhuman / (Vhuman + Vsurf) 
         else
            ratioWD(ii) = 1._r8
         end if
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
      Vtot = sum(zVol)

      dV = Vtot - Vtot0 - m_hydroData%Qwi*dt + m_hydroData%Qws*dt + &
         sum(m_hydroData%Qwo)*dt
      if (abs(dV) > err_max) then
         print *, " Lake storage imbalance: ", dV
         call Endrun('Lake storage does not conserve after inflow/outflow')
      end if

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
         ! only update layer thickness 
         m_dZw = zVol / m_Az
         do ii = WATER_LAYER+1, 1, -1
            if (ii==WATER_LAYER+1) then
               m_Zw(ii) = 0._r8
            else if (ii==WATER_LAYER) then
               m_Zw(ii) = m_Zw(ii+1) + m_dZw(ii+1) + 0.5*m_dZw(ii)
            else if (ii==1) then
               m_Zw(ii) = m_Zw(ii+1) + 0.5*m_dZw(ii+1) + m_dZw(ii)
            else
               m_Zw(ii) = m_Zw(ii+1) + 0.5*m_dZw(ii+1) + 0.5*m_dZw(ii)
            end if
         end do

         return
      end if

      ! adjust lake grid 
      Vtot0 = sum(zVol)
      call GetHeightFromVolume(Vtot0, m_hypsocurve, maxZ)
      call BuildDepthVector(maxZ, lake_info%hsed, m_hypsocurve, m_Zw, m_dZw, &
               m_Zs, m_dZs)
      call BuildLakeStructure(lake_info, m_hypsocurve, Vtot0, m_Zw, m_dZw, &
               m_Az, m_dAz, m_SAL, m_soilColZ, m_soilColInd)
      Vtot = sum(m_dZw*m_Az)

      if (abs(Vtot - Vtot0) > 1.d-3) then
         call Endrun('Lake storage does not conserve after regridding')
      end if

      ! remap state variables onto the new grids
      zVol0 = zVol
      indx = 1
      !Ebg = sum(zVol0*m_waterTemp) / sum(zVol0)
      do ii = 1, WATER_LAYER+1, 1
         Vres = m_dZw(ii) * m_Az(ii)
         Vnew = 0._r8
         tmpwaterHeat(ii) = 0._r8
         tmpwaterIce(ii) = 0._r8
         tmpwaterSubCon(:,ii) = 0._r8
         do while (Vres>e8 .and. indx<=WATER_LAYER+1)
            if (Vres<=zVol0(indx)) then
               tmpwaterHeat(ii) = tmpwaterHeat(ii) + Vres* &
                  Cpi*(m_waterTemp(indx)-T0)*m_waterIce(indx) + Vres* &
                  Cpl*(m_waterTemp(indx)-T0)*(1.0-m_waterIce(indx))
               tmpwaterIce(ii) = tmpwaterIce(ii) + Vres*m_waterIce(indx)
               tmpwaterSubCon(:,ii) = tmpwaterSubCon(:,ii) + Vres* &
                  m_waterSubCon(:,indx)
               Vnew = Vnew + Vres
               zVol0(indx) = zVol0(indx) - Vres
               Vres = 0._r8
            else
               tmpwaterHeat(ii) = tmpwaterHeat(ii) + zVol0(indx)* &
                  Cpi*(m_waterTemp(indx)-T0)*m_waterIce(indx) + zVol0(indx)* &
                  Cpl*(m_waterTemp(indx)-T0)*(1.0-m_waterIce(indx))
               tmpwaterIce(ii) = tmpwaterIce(ii) + zVol0(indx)* &
                  m_waterIce(indx)
               tmpwaterSubCon(:,ii) = tmpwaterSubCon(:,ii) + zVol0(indx)* &
                  m_waterSubCon(:,indx)
               Vnew = Vnew + zVol0(indx)
               Vres = Vres - zVol0(indx)
               zVol0(indx) = 0._r8
               indx = indx + 1
            end if
         end do

         ! average
         tmpwaterIce(ii) = tmpwaterIce(ii) / Vnew
         tmpEres1 = Lf * (1.0-tmpwaterIce(ii)) * Vnew
         tmpEres2 = Lf * tmpwaterIce(ii) * Vnew
         if (tmpwaterIce(ii)<e8) then
            tmpwaterTemp(ii) = T0 + tmpwaterHeat(ii) / Vnew / Cpl
         else if (tmpwaterHeat(ii)>-e8) then
            if (tmpEres2>tmpwaterHeat(ii)) then
               tmpWaterIce(ii) = tmpWaterIce(ii) - tmpwaterHeat(ii) / Lf / Vnew 
               tmpwaterTemp(ii) = T0
            else
               tmpWaterIce(ii) = 0._r8
               tmpwaterTemp(ii) = T0 + (tmpwaterHeat(ii)-tmpEres2) / Cpl / Vnew
            end if
         else
            if (tmpEres1>-tmpwaterHeat(ii)) then
               tmpWaterIce(ii) = tmpWaterIce(ii) - tmpwaterHeat(ii) / Lf / Vnew
               tmpwaterTemp(ii) = T0
            else
               tmpWaterIce(ii) = 1._r8
               tmpwaterTemp(ii) = T0 + (tmpwaterHeat(ii)+tmpEres1) / Cpi / Vnew
            end if
         end if
         tmpwaterSubCon(:,ii) = tmpwaterSubCon(:,ii) / Vnew
      end do

      ! update state variables
      zVol = m_dZw * m_Az        ! update current volume

      if (Thermal_Module) then
         m_waterTemp = tmpwaterTemp
         m_waterIce = tmpwaterIce
         call EnforceThermalConsistency() 

         ! check energy balance
         !Eed = sum(zVol*m_waterTemp) / sum(zVol)
         !if (abs(Ebg - Eed) > 1.d-6) then
         !   print *, "Ebg and Eed: ", Ebg, Eed, Ebg-Eed
         !   call Endrun('Heat budget does not conserve after remapping')
         !end if
      end if

      if (Carbon_Module) then
         m_waterSubCon = tmpwaterSubCon
         call EnforceBGCConsistency() 
      end if

   end subroutine 

   !------------------------------------------------------------------------------
   !
   ! Purpose: Get outflow thermal and chemical fluxes.
   !
   !------------------------------------------------------------------------------
   subroutine GetWithdrawFlow(Sw, Qwt, Qo2, Qco2, Qch4, Qsrp)
      implicit none
      real(r8), intent(out) :: Sw      ! m^3
      real(r8), intent(out) :: Qwt     ! K
      real(r8), intent(out) :: Qo2     ! mol/m3
      real(r8), intent(out) :: Qco2    ! mol/m3
      real(r8), intent(out) :: Qch4    ! mol/m3
      real(r8), intent(out) :: Qsrp    ! mol/m3
      real(r8) :: Qwo

      Sw = sum(zVol)
      Qwo = sum(ratioWD*realWD)    ! m3
      if (Qwo > e8 .and. Thermal_Module) then
         Qwt = sum( ratioWD * realWD * m_waterTemp ) / Qwo
      else
         Qwt = -9999._r8
      end if
      if (Qwo > e8 .and. Carbon_Module) then
         Qo2 = sum( ratioWD * realWD * m_waterSubCon(Wo2,:) ) / Qwo
         Qco2 = sum( ratioWD * realWD * m_waterSubCon(Wco2,:) ) / Qwo
         Qch4 = sum( ratioWD * realWD * m_waterSubCon(Wch4,:) ) / Qwo
         Qsrp = sum( ratioWD * realWD * m_waterSubCon(Wsrp,:) ) / Qwo
      else
         Qo2 = -9999._r8
         Qco2 = -9999._r8
         Qch4 = -9999._r8
         Qsrp = -9999._r8
      end if
   end subroutine

end module hydro_mod
