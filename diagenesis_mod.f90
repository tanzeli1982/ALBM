module diagenesis_mod
!---------------------------------------------------------------------------------
! Purpose: this module governs the diagenetic reactions of organic carbon in the 
!          sediments, which converts organic carbon to either CO2 and CH4.
!
!          [Stepanenko et al., 2011; Izvestiya, Atmospheric and 
!           Oceanic Physics], [Kessler et al., 2012; JGR] and [Hanson et al.,
!           2011; PLus One].
!
!---------------------------------------------------------------------------------
   use shr_ctrl_mod, inft => INFINITESIMAL_E8
   use shr_typedef_mod,       only : RungeKuttaCache3D
   use shr_param_mod
   use phy_utilities_mod
   use bg_utilities_mod
   use data_buffer_mod

   implicit none
   private
   public :: InitializeDiagenesisModule, DestructDiagenesisModule
   public :: DiagenesisModuleSetup, DiagenesisModuleCallback 
   public :: DiagenesisEquation, GetSedCLossRate 
   public :: ConstructOldCarbonPool, ConstructActCarbonPool
   public :: mem_ch4
   ! memory cache for Runge-Kutta
   type(RungeKuttaCache3D) :: mem_ch4
   ! CH4 anaerobic production rates (mol/m3/s)
   real(r8), allocatable :: rPCH4(:,:,:)
   ! CO2 anaerobic production rates (mol/m3/s) 
   real(r8), allocatable :: rnPCO2(:,:,:)
   ! CO2 aerobic production rates (mol/m3/s) 
   real(r8), allocatable :: rPCO2(:,:,:)
   ! pore-water substance dynamcis (mol/m3/s)
   real(r8), allocatable :: rDYN(:,:,:)
   ! gas ebullition rates in the sediments (mol/m3/s)
   real(r8), allocatable :: Ebch4(:,:,:)
   ! lagged ice percentage
   real(r8), allocatable :: lagged_ice(:,:)
   ! diffusivity in the sediments (m2/s)
   real(r8), allocatable :: Kch4(:,:)
   ! water-sediment diffusion (mol/m2/s)
   real(r8), allocatable :: qt(:,:)
   ! active carbon layer
   real(r8) :: ztot_act
   integer :: indx_act

contains
   subroutine InitializeDiagenesisModule()
      implicit none
      integer :: ii

      allocate(rDYN(NSSUB,NSCOL,SED_LAYER+1))
      allocate(rPCH4(NPOOL,NSCOL,SED_LAYER+1))
      allocate(rPCO2(NPOOL,NSCOL,SED_LAYER+1))
      allocate(rnPCO2(NPOOL,NSCOL,SED_LAYER+1))
      allocate(Ebch4(NGAS,NSCOL,SED_LAYER+1))
      allocate(lagged_ice(NSCOL,SED_LAYER+1))
      allocate(Kch4(NSCOL,SED_LAYER+1))
      allocate(qt(NSSUB,NSCOL))
      allocate(mem_ch4%K1(NSSUB,NSCOL,SED_LAYER+1))
      allocate(mem_ch4%K2(NSSUB,NSCOL,SED_LAYER+1))
      allocate(mem_ch4%K3(NSSUB,NSCOL,SED_LAYER+1))
      allocate(mem_ch4%K4(NSSUB,NSCOL,SED_LAYER+1))
      allocate(mem_ch4%K5(NSSUB,NSCOL,SED_LAYER+1))
      allocate(mem_ch4%K6(NSSUB,NSCOL,SED_LAYER+1))
      allocate(mem_ch4%nxt4th(NSSUB,NSCOL,SED_LAYER+1))
      allocate(mem_ch4%nxt5th(NSSUB,NSCOL,SED_LAYER+1))
      allocate(mem_ch4%interim(NSSUB,NSCOL,SED_LAYER+1))
      allocate(mem_ch4%rerr(NSSUB,NSCOL,SED_LAYER+1))

      call InitializeDiagenesisStateVariables()
      Ebch4 = 0.0_r8
      rDYN = 0.0_r8
      rPCH4 = 0.0_r8
      rPCO2 = 0.0_r8
      rnPCO2 = 0.0_r8
      Kch4 = 0.0_r8
      qt = 0.0_r8
      indx_act = 2
      ztot_act = sum(m_dZs(1:indx_act))
   end subroutine

   subroutine DestructDiagenesisModule()
      implicit none

      deallocate(rDYN)
      deallocate(rPCH4)
      deallocate(rPCO2)
      deallocate(rnPCO2)
      deallocate(Ebch4)
      deallocate(lagged_ice)
      deallocate(Kch4)
      deallocate(qt)
      deallocate(mem_ch4%K1)
      deallocate(mem_ch4%K2)
      deallocate(mem_ch4%K3)
      deallocate(mem_ch4%K4)
      deallocate(mem_ch4%K5)
      deallocate(mem_ch4%K6)
      deallocate(mem_ch4%nxt4th)
      deallocate(mem_ch4%nxt5th)
      deallocate(mem_ch4%interim)
      deallocate(mem_ch4%rerr)
   end subroutine
   
   subroutine DiagenesisModuleSetup(isUpdate)
      implicit none
      logical, intent(in) :: isUpdate
      real(r8) :: minP, minFe, Do2
      integer :: icol, top
      integer :: indx0, indx, newindx0

      top = m_lakeWaterTopIndex
      if (isUpdate) then
         ! update sediment SRP
         do icol = 1, NSCOL, 1
            indx0 = COUNT(m_soilColInd<=icol-1) + 1
            indx = COUNT(m_soilColInd<=icol)
            ! invalid sediment column
            if (indx0>indx .or. indx<top) then
               cycle
            end if
            newindx0 = max(indx0, top)
            minP = m_sedSubCon(Wsrp,icol,1)
            minFe = lake_info%sedFe
            call Mean(m_waterSubCon(Wo2,newindx0:indx), Do2)
            Do2 = max(Do2, 0.0625_r8)
            m_sedSRPCon(icol) = CalcSedSRPConc(minP, minFe, Do2)
         end do
         call CalcDiagenesisRates()
         call UpdateSubDiffusivity()
      end if

      call CalcSurfaceExchange()
      call UpdateBubbleFlux()
   end subroutine

   subroutine DiagenesisModuleCallback(dt)
      implicit none
      real(r8), intent(in) :: dt

      lagged_ice = m_sedIce / m_sedpor
      call AdjustNegativeConcentration()
      call UpdateActCarbonPool(dt)
      call UpdatePasCarbonPool()
      call UpdateOldCarbonPool(dt)
      call UpdateSurfaceminP()
      call UpdateSediRedox(dt)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Govern the production of methane bubbles
   !
   !------------------------------------------------------------------------------
   subroutine DiagenesisEquation(con, dcon)
      implicit none
      real(r8), intent(in) :: con(NSSUB,NSCOL,SED_LAYER+1)
      real(r8), intent(out) :: dcon(NSSUB,NSCOL,SED_LAYER+1)
      real(r8), dimension(NSSUB) :: qb, dC1, dC2
      real(r8) :: a, b
      integer :: Top_Index, Bottom_Index
      integer :: ii, icol

      do icol = 1, NSCOL, 1
         ! inactive sediment column
         if (COUNT(m_soilColInd==icol)==0) then
            dcon(:,icol,:) = 0.0_r8
            cycle
         end if
         ! frozen sediment column
         Top_Index = m_sedWaterTopIndex(icol)
         Bottom_Index = m_sedWaterBtmIndex(icol)
         if (Top_Index>Bottom_Index) then
            dcon(:,icol,:) = 0.0_r8
            cycle 
         end if
         qb = 0.0_r8
         do ii = 1, SED_LAYER+1, 1
            if (ii==1) then
               a = 0.5*(Kch4(icol,ii) + Kch4(icol,ii+1))
               dC1 = (con(:,icol,ii+1) - con(:,icol,ii)) / (m_Zs(ii+1) - m_Zs(ii))
               dcon(:,icol,ii) = (a*dC1 - qt(:,icol)) / m_dZs(ii) + rDYN(:,icol,ii)
            else if (ii==SED_LAYER+1) then
               b = 0.5*(Kch4(icol,ii-1) + Kch4(icol,ii))
               dC2 = (con(:,icol,ii) - con(:,icol,ii-1)) / (m_Zs(ii) - m_Zs(ii-1))
               dcon(:,icol,ii) = (qb - b*dC2) / m_dZs(ii) + rDYN(:,icol,ii)
            else
               a = 0.5*(Kch4(icol,ii) + Kch4(icol,ii+1))
               b = 0.5*(Kch4(icol,ii-1) + Kch4(icol,ii))
               dC1 = (con(:,icol,ii+1) - con(:,icol,ii)) / (m_Zs(ii+1) - m_Zs(ii))
               dC2 = (con(:,icol,ii) - con(:,icol,ii-1)) / (m_Zs(ii) - m_Zs(ii-1))
               dcon(:,icol,ii) = (a*dC1 - b*dC2) / m_dZs(ii) + rDYN(:,icol,ii)
            end if
         end do
      end do
      where(con<=0 .and. dcon<0) dcon = 0.0_r8
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Adjust negative sediment gas concentration (project negative values 
   !          to zero) - "Non-negative solutions of ODEs" (Shampine, L., 2005).
   !
   !------------------------------------------------------------------------------
   subroutine AdjustNegativeConcentration()
      implicit none

      where (m_sedSubCon<0) m_sedSubCon = 0.0_r8
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate substance prod/loss rates
   !
   !------------------------------------------------------------------------------
   subroutine CalcDiagenesisRates()
      implicit none
      real(r8) :: Prtot, Prnet, Ps, Psat
      real(r8) :: Hco2, Ho2, Hch4, Hn2 
      real(r8) :: Porosity, Re, Ae, pH, sal
      real(r8) :: pch4, pn2, temp, depth
      real(r8) :: Do2, eHL, carb(NPOOL)
      integer :: ii, icol, indx
      logical :: isCH4sat

      Re = sa_params(Param_Re)
      Ae = sa_params(Param_Ae)
      Ps = m_surfData%pressure - Roul*G*m_surfData%dzsurf
      pH = lake_info%pH
      sal = lake_info%sal
      isCH4sat = .False.
      do icol = 1, NSCOL, 1
         ! inactive sediment column
         if (COUNT(m_soilColInd==icol)==0) then
            rPCH4(:,icol,:) = 0.0_r8
            rPCO2(:,icol,:) = 0.0_r8
            rnPCO2(:,icol,:) = 0.0_r8
            Ebch4(:,icol,:) = 0.0_r8
            cycle
         end if
         ! active sediment column
         indx = COUNT(m_soilColInd<=icol)
         do ii = 1, SED_LAYER+1, 1
            Porosity = m_sedpor(icol,ii)
            if (m_sedIce(icol,ii)>=Porosity) then
               rPCH4(:,icol,ii) = 0.0_r8
               rPCO2(:,icol,ii) = 0.0_r8
               rnPCO2(:,icol,ii) = 0.0_r8
               Ebch4(:,icol,ii) = 0.0_r8
            else
               ! carbon mineralization rates
               temp = m_sedTemp(icol,ii)
               Do2 = m_sedSubCon(Wo2,icol,ii)
               eHL = m_sedEHL(icol,ii)
               carb = m_unfrzCarbPool(:,icol,ii)
               call Methanogenesis(carb, temp, sal, eHL, rPCH4(:,icol,ii), &
                     rnPCO2(:,icol,ii)) 
               call AerobicRespiration(carb, temp, Do2, rPCO2(:,icol,ii))
               ! ebullition rates (Eq. (A.1), Stepanenko et al.)
               depth = m_Zs(ii) + m_Zw(indx) 
               Hch4 = CalcHenrySolubility(Wch4, temp, pH)
               Hn2 = CalcHenrySolubility(Wn2, temp, pH)
               pch4 = m_sedSubCon(Wch4,icol,ii)/Hch4
               pn2 = m_sedSubCon(Wn2,icol,ii)/Hn2
               Prtot = pn2 + pch4
               Psat = Ae * Porosity * (Ps + Roul*G*depth)
               Prnet = Prtot - Psat
               if (icol==NSCOL .and. ii==1 .and. Prnet>e8) then
                  isCH4sat = .True.
               end if
               if (Prnet>e8) then
                  Ebch4(Wn2,icol,ii) = Re * (pn2/Prtot) * Prnet * Hn2
                  Ebch4(Wo2,icol,ii) = 0.0_r8
                  Ebch4(Wch4,icol,ii) = Re * (pch4/Prtot) * Prnet * Hch4
                  Ebch4(Wco2,icol,ii) = 0.0_r8
               else
                  Ebch4(:,icol,ii) = 0.0_r8 
               end if
            end if
         end do
         ! update rDYN
         rDYN(Wn2,icol,:) = -Ebch4(Wn2,icol,:)
         rDYN(Wo2,icol,:) = -sum(rPCO2(:,icol,:),1) - Ebch4(Wo2,icol,:)
         rDYN(Wco2,icol,:) = sum(rPCO2(:,icol,:),1) + sum(rnPCO2(:,icol,:),1) - &
               Ebch4(Wco2,icol,:)
         rDYN(Wch4,icol,:) = sum(rPCH4(:,icol,:),1) - Ebch4(Wch4,icol,:)
         rDYN(Wsrp,icol,:) = 0._r8  ! fix mineral P
         !rDYN(Wsrp,icol,:) = ( sum(rPCH4(:,icol,:),1) + sum(rPCO2(:,icol,:),1) + &
         !      sum(rnPCO2(:,icol,:),1) ) / YC2P_POM
         
         if (icol==NSCOL .and. (.NOT. isCH4sat)) then
            rDYN(Wch4,icol,1) = rDYN(Wch4,icol,1) + sum(Ebch4(Wch4,icol,:)* &
               m_dZs)/m_dZs(1)
            rDYN(Wn2,icol,1) = rDYN(Wn2,icol,1) + sum(Ebch4(Wn2,icol,:)* &
               m_dZs)/m_dZs(1)
            Ebch4(Wch4,icol,:) = 0._r8
            Ebch4(Wn2,icol,:) = 0._r8
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: get sediment diffusive flux and P uptake.
   !
   !------------------------------------------------------------------------------
   subroutine CalcSurfaceExchange()
      implicit none
      real(r8) :: dflux_sed, sumAz, sumV
      real(r8) :: conc_sed, conc_wat, conc_tmp
      real(r8) :: Porosity
      integer :: kk, ii, icol, top
      integer :: indx, indx0

      qt = 0._r8
      top = m_lakeWaterTopIndex
      do icol = 1, NSCOL, 1
         indx0 = COUNT(m_soilColInd<=icol-1) + 1
         indx = COUNT(m_soilColInd<=icol)
         ! inactive sediment column 
         if (indx0>indx .or. indx<top) then
            cycle
         end if
         sumAz = sum(m_dAz(indx0:indx))
         sumV = sum(m_dZw(indx0:indx)*m_Az(indx0:indx))
         Porosity = m_sedpor(icol,1)

         do kk = 1, NSSUB, 1
            if (icol<NSCOL) then
               conc_wat = sum(m_waterSubCon(kk,indx0:indx)* &
                  m_dZw(indx0:indx)*m_Az(indx0:indx)) / sumV 
            else
               conc_wat = m_waterSubCon(kk,indx)
            end if
            if (kk==Wsrp) then
               conc_sed = m_sedSRPCon(icol)
               !dflux_sed = m_Ktb(indx) * (conc_sed - conc_wat) / DeltaD + &
               !  m_vegPUptake(ii)
               dflux_sed = 0._r8   ! fix mineral P
            else
               conc_sed = m_sedSubCon(kk,icol,1) / Porosity
               dflux_sed = m_Ktb(indx) * (conc_sed - conc_wat) * &
                  Porosity / DeltaD
            end if
            if (kk==Wn2 .or. kk==Wo2) then
               dflux_sed = min(0._r8, dflux_sed)
            else if (kk==Wco2 .or. kk==Wch4) then
               dflux_sed = max(0._r8, dflux_sed)
            end if
            do ii = indx0, indx, 1
               conc_tmp = m_waterSubCon(kk,ii)
               if (ii>=top) then
                  if (dflux_sed>0._r8 .and. conc_sed>conc_tmp) then
                     qt(kk,icol) = qt(kk,icol) + dflux_sed*m_dAz(ii)/sumAz
                  else if (dflux_sed<0._r8 .and. conc_sed<conc_tmp) then
                     qt(kk,icol) = qt(kk,icol) + dflux_sed*m_dAz(ii)/sumAz
                  end if
               end if
            end do
         end do
      end do
      m_btmdflux = qt
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update bubble flux rates from sediment top.
   !
   !------------------------------------------------------------------------------
   subroutine UpdateBubbleFlux()
      implicit none
      integer :: gas, icol

      do icol = 1, NSCOL, 1
         if (COUNT(m_soilColInd==icol)==0) then
            m_btmbflux(:,icol) = 0._r8
         else
            do gas = 1, NGAS, 1
               m_btmbflux(gas,icol) = sum(Ebch4(gas,icol,:)*m_dZs)
            end do
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate dissolved substance diffusivity in sediments
   !
   !------------------------------------------------------------------------------
   subroutine UpdateSubDiffusivity()
      implicit none
      real(r8) :: Porosity, diffst
      integer :: ii, icol

      do icol = 1, NSCOL, 1
         if (COUNT(m_soilColInd==icol)==0) then
            Kch4(icol,:) = 0.528d-7
         else
            do ii = 1, SED_LAYER+1, 1
               Porosity = m_sedpor(icol,ii)
               if (m_sedIce(icol,ii)>=Porosity) then
                  diffst = CalcGasDiffusivityInIce(Wch4, m_sedTemp(icol,ii))
                  Kch4(icol,ii) = diffst*0.66*Porosity 
               else
                  ! Eq. (8) in Walter et al. (2000) gives Kch4 at the magnitude
                  ! of 1e-9 but Goto et al. (2017; Marine Geophysical Research) shows
                  ! this value at the magnitude of 1e-7
                  Kch4(icol,ii) = 2.0d-7*0.66*Porosity   ! [Eq. (8) in Walter 2000]
               end if
            end do
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Construct active and passive carbon pool profile in sediments.
   !          Active carbon pool: [West and plug, 2008 JGR]
   !          Passive carbon pool: [Stepanenko et al., 2011]
   !  
   !------------------------------------------------------------------------------
   subroutine ConstructOldCarbonPool()
      implicit none
      real(r8) :: Ht, Zs, time, Ctotal
      real(r8) :: Ct, Rco, Porosity, depth
      integer :: ii, itk, icol, indx
      logical :: margin 

      if (lake_info%thrmkst==0) then
         m_frzCarbPool(oldC,:,:) = 0.0_r8
         m_unfrzCarbPool(oldC,:,:) = 0.0_r8
         return
      end if
      
      Rco = sa_params(Param_Rcold)
      do icol = 1, NSCOL, 1
         ! inactive sediment column
         if (COUNT(m_soilColInd==icol)==0) then
            m_frzCarbPool(oldC,icol,:) = 0.0_r8
            m_unfrzCarbPool(oldC,icol,:) = 0.0_r8
            cycle
         end if

         indx = COUNT(m_soilColInd<=icol)
         margin = (icol==1)
         ! permafrost thawing rate (m/(yr^0.5))
         if (margin) then
            Ct = 0.67_r8
         else
            Ct = 0.77_r8
         end if
         Ht = 0.0_r8
         do ii = 1, SED_LAYER+1, 1
            Porosity = m_sedpor(icol,ii)
            if (m_sedIce(icol,ii)>=Porosity) then
               exit
            end if
            Ht = Ht + m_dZs(ii)*(1.0 - m_sedIce(icol,ii)/Porosity)
            itk = ii
         end do
         if (abs(lake_info%hsed-Ht)<e8) then
            m_frzCarbPool(oldC,icol,:) = 0.0_r8
            m_unfrzCarbPool(oldC,icol,:) = 0.0_r8
         else
            if (lake_info%thrmkst==1) then
               depth = 5.0
            else if (lake_info%thrmkst==2) then
               depth = 25.0
            end if
            do ii = 1, SED_LAYER+1, 1
               Zs = m_Zs(ii)
               if (Zs+m_Zw(indx)>depth) then
                  m_frzCarbPool(oldC,icol,ii) = 0.0_r8
                  m_unfrzCarbPool(oldC,icol,ii) = 0.0_r8
               else
                  ! mean labile carbon density (gC/m3)
                  if (ii<=itk) then
                     ! time that talik has developed
                     time = SECOND_OF_YEAR * (Ht**2 - Zs**2) / Ct**2
                     Ctotal = max(0.0_r8, oldcarb0*(1.0-Rco*time))
                  else
                     Ctotal = oldcarb0
                  end if
                  Ctotal = Ctotal * 1.0d3 / 3.0
                  m_frzCarbPool(oldC,icol,ii) = Ctotal * lagged_ice(icol,ii)
                  m_unfrzCarbPool(oldC,icol,ii) = Ctotal - m_frzCarbPool(oldC,icol,ii)
               end if
            end do
         end if
      end do
   end subroutine

   subroutine ConstructPasCarbonPool()
      implicit none
      real(r8) :: csed, carbon, dampen 
      integer :: icol

      csed = 1.0d3 * lake_info%csed ! gC/m3
      dampen = sa_params(Param_csedDMP)
      do icol = 1, NSCOL, 1
         ! inactive sediment column
         if (COUNT(m_soilColInd==icol)==0) then
            m_frzCarbPool(pasC,icol,:) = 0.0_r8
            m_unfrzCarbPool(pasC,icol,:) = 0.0_r8
         else 
            carbon = csed * dampen / (1._r8 - exp(-dampen))
            m_frzCarbPool(pasC,icol,:) = carbon * lagged_ice(icol,:) * &
               exp(-dampen*m_Zs)
            m_unfrzCarbPool(pasC,icol,:) = carbon * (1.0 - &
               lagged_ice(icol,:)) * exp(-dampen*m_Zs)
         end if
      end do
   end subroutine

   subroutine ConstructActCarbonPool()
      implicit none
      real(r8) :: csed, Porosity
      integer :: icol, ii

      csed = 0.02_r8 * 1.0d3 * lake_info%csed ! gC/m3
      do icol = 1, NSCOL, 1
         ! inactive sediment column
         if (COUNT(m_soilColInd==icol)==0) then
            m_frzCarbPool(actC,icol,:) = 0._r8
            m_unfrzCarbPool(actC,icol,:) = 0._r8
            cycle
         end if
         do ii = 1, SED_LAYER+1, 1
            if (ii<=indx_act) then
               Porosity = m_sedpor(icol,ii)
               m_frzCarbPool(actC,icol,ii) = csed * m_sedIce(icol,ii) / &
                     Porosity
               m_unfrzCarbPool(actC,icol,ii) = csed * (1.0 - &
                     m_sedIce(icol,ii) / Porosity)
            else
               m_frzCarbPool(actC,icol,ii) = 0._r8
               m_unfrzCarbPool(actC,icol,ii) = 0._r8
            end if
         end do
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update the active and passive carbon profile in sediments.
   !          Active carbon pool: [West and plug, 2008 JGR]
   !          Passive carbon pool: [Stepanenko et al., 2011]
   !  
   !------------------------------------------------------------------------------
   subroutine UpdateOldCarbonPool(dt)
      implicit none
      real(r8), intent(in) :: dt       ! time interval (sec)
      real(r8) :: frice_lag, frice, ratio
      real(r8) :: Porosity, decay, Ctotal
      integer  :: ii, icol

      if (lake_info%thrmkst==0) then
         m_frzCarbPool(oldC,:,:) = 0.0_r8
         m_unfrzCarbPool(oldC,:,:) = 0.0_r8
         return
      end if

      do icol = 1, NSCOL, 1
         ! inactive sediment column
         if (COUNT(m_soilColInd==icol)==0) then
            m_frzCarbPool(oldC,icol,:) = 0.0_r8
            m_unfrzCarbPool(oldC,icol,:) = 0.0_r8
            cycle
         end if
         do ii = 1, SED_LAYER+1, 1
            ! Shrink C pool due to degradation
            Porosity = m_sedpor(icol,ii)
            if (m_unfrzCarbPool(oldC,icol,ii)>e8) then
               decay = dt * (rPCH4(oldC,icol,ii) + rnPCO2(oldC,icol,ii) + &
                  rPCO2(oldC,icol,ii)) * MasC
               m_unfrzCarbPool(oldC,icol,ii) = max(0.0_r8, &
                  m_unfrzCarbPool(oldC,icol,ii) - decay)
            end if
            ! redistribute frozen and unfrozen carbon pools
            Ctotal = m_frzCarbPool(oldC,icol,ii) + m_unfrzCarbPool(oldC,icol,ii)
            if (Ctotal>e8) then
               frice = m_sedIce(icol,ii) / Porosity 
               frice_lag = lagged_ice(icol,ii)
               if (frice>frice_lag) then
                  ! ice expand and freezing carbon
                  ratio = (frice - frice_lag) / (1.0 - frice_lag)
                  m_frzCarbPool(oldC,icol,ii) = m_frzCarbPool(oldC,icol,ii) + &
                     ratio * m_unfrzCarbPool(oldC,icol,ii)
                  m_unfrzCarbPool(oldC,icol,ii) = m_unfrzCarbPool(oldC,icol,ii) * &
                     (1.0 - ratio)
               else if (frice<frice_lag) then
                  ! ice shrink and thawing carbon
                  ratio = (frice_lag - frice) / frice_lag
                  m_unfrzCarbPool(oldC,icol,ii) = m_unfrzCarbPool(oldC,icol,ii) + &
                     ratio * m_frzCarbPool(oldC,icol,ii)
                  m_frzCarbPool(oldC,icol,ii) = m_frzCarbPool(oldC,icol,ii) * &
                     (1.0 - ratio)
               end if
            else
               m_frzCarbPool(oldC,icol,ii) = 0.0_r8
               m_unfrzCarbPool(oldC,icol,ii) = 0.0_r8
            end if
         end do
      end do
   end subroutine

   subroutine UpdatePasCarbonPool()
      implicit none
      real(r8) :: Porosity, ratio
      real(r8) :: frice_lag, frice, Ctotal
      integer :: ii, icol

      do icol = 1, NSCOL, 1
         ! inactive sediment column
         if (COUNT(m_soilColInd==icol)==0) then
            cycle
         end if
         ! freeze/thaw
         do ii = 1, SED_LAYER+1, 1
            Porosity = m_sedpor(icol,ii)
            Ctotal = m_frzCarbPool(pasC,icol,ii) + m_unfrzCarbPool(pasC,icol,ii)
            if (Ctotal>e8) then
               ! redistribute frozen and unfrozen carbon pools
               frice = m_sedIce(icol,ii) / Porosity
               frice_lag = lagged_ice(icol,ii)
               if (frice>frice_lag) then
                  ! ice expand and freezing carbon
                  ratio = (frice - frice_lag) / (1.0 - frice_lag)
                  m_frzCarbPool(pasC,icol,ii) = m_frzCarbPool(pasC,icol,ii) + &
                     ratio * m_unfrzCarbPool(pasC,icol,ii)
                  m_unfrzCarbPool(pasC,icol,ii) = m_unfrzCarbPool(pasC,icol,ii) * &
                     (1.0 - ratio)
               else if (frice<frice_lag) then     
                  ! ice shrink and thawing carbon
                  ratio = (frice_lag - frice) / frice_lag
                  m_unfrzCarbPool(pasC,icol,ii) = m_unfrzCarbPool(pasC,icol,ii) + &
                     ratio * m_frzCarbPool(pasC,icol,ii)
                  m_frzCarbPool(pasC,icol,ii) = m_frzCarbPool(pasC,icol,ii) * &
                     (1.0 - ratio)
               end if
            else
               m_frzCarbPool(pasC,icol,ii) = 0.0_r8
               m_unfrzCarbPool(pasC,icol,ii) = 0.0_r8
            end if
         end do
      end do
   end subroutine

   subroutine UpdateActCarbonPool(dt)
      implicit none
      real(r8), intent(in) :: dt ! time step (s)
      real(r8) :: ctot_dep, ctot_decay
      real(r8) :: ctot_act, cavg
      integer :: icol, ii

      do icol = 1, NSCOL, 1
         ! inactive sediment column
         if (COUNT(m_soilColInd==icol)==0) then
            cycle
         end if
         ! active C deposition
         ctot_dep = (m_trDepAtCarb(icol) + m_aqDepAtCarb(icol) + &
            m_vegDepAtCarb(icol)) * dt
         ! carbon decay
         ctot_decay = 0.0_r8
         ctot_act = 0.0_r8
         do ii = 1, SED_LAYER+1, 1
            if (m_unfrzCarbPool(actC,icol,ii)>e8) then
               ctot_decay = ctot_decay + dt * m_dZs(ii) * (rPCH4(actC,icol,ii) + &
                  rnPCO2(actC,icol,ii) + rPCO2(actC,icol,ii)) * MasC
               ctot_act = ctot_act + m_unfrzCarbPool(actC,icol,ii)*m_dZs(ii)
            end if
            ctot_act = ctot_act + m_frzCarbPool(actC,icol,ii)*m_dZs(ii)
         end do
         ! redistribution
         cavg = max( (ctot_act + ctot_dep - ctot_decay)/ztot_act, 0.0_r8 )
         do ii = 1, indx_act, 1
            m_unfrzCarbPool(actC,icol,ii) = cavg * (1._r8 - m_sedIce(icol,ii) / &
               m_sedpor(icol,ii))
            m_frzCarbPool(actC,icol,ii) = cavg * m_sedIce(icol,ii) / &
               m_sedpor(icol,ii)
         end do
         m_unfrzCarbPool(actC,icol,indx_act+1:SED_LAYER+1) = 0.0_r8
         m_frzCarbPool(actC,icol,indx_act+1:SED_LAYER+1) = 0.0_r8
      end do
   end subroutine

   subroutine UpdateSurfaceminP()
      implicit none
      real(r8) :: minP_avg
      integer :: icol

      do icol = 1, NSCOL, 1
         ! inactive sediment column
         if (COUNT(m_soilColInd==icol)==0) then
            cycle
         end if
         minP_avg = sum( m_sedSubCon(Wsrp,icol,1:indx_act) * &
            m_dZs(1:indx_act) ) / ztot_act
         m_sedSubCon(Wsrp,icol,1:indx_act) = minP_avg
      end do
   end subroutine

   subroutine UpdateSediRedox(dt)
      implicit none
      real(r8), intent(in) :: dt
      real(r8) :: cRx, Do2, Do2_crit 
      integer :: icol, top, ii 
      integer :: indx0, indx

      top = m_lakeWaterTopIndex
      cRx = sa_params(Param_cRx)
      ! update sediment redox potential
      do icol = 1, NSCOL, 1
         indx0 = COUNT(m_soilColInd<=icol-1) + 1
         indx = COUNT(m_soilColInd<=icol)
         ! invalid sediment column
         if (indx0>indx .or. indx<top) then
            cycle
         end if
         do ii = 1, SED_LAYER+1, 1
            Do2 = m_sedSubCon(Wo2,icol,ii)
            if (Do2<=Do2_cr) then
               m_sedEHL(icol,ii) = m_sedEHL(icol,ii) - cRx / 8.64d4 * dt
            else
               m_sedEHL(icol,ii) = m_sedEHL(icol,ii) + cRx / 8.64d4 * &
                  min(1._r8, (Do2-Do2_cr)/0.45625) * dt
            end if
            ! Mobilian, C., & Craft, C. B (2022). Wetland Soils: Physical and 
            ! Chemical Properties and Biogeochemical Processes. 157-168.
            m_sedEHL(icol,ii) = max(min(m_sedEHL(icol,ii),1.d2),-3.d2)
         end do
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate CH4/CO2 production rate.
   !
   !------------------------------------------------------------------------------
   subroutine GetSedCLossRate(pch4)
      implicit none
      real(r8), intent(out) :: pch4(:)   ! units: mol/m2/s
      integer :: icol, kk 

      pch4 = 0._r8
      do icol = 1, NSCOL, 1
         if (COUNT(m_soilColInd==icol)>0) then
            do kk = 1, NPOC, 1
               pch4(icol) = pch4(icol) + sum( rPCH4(kk,icol,:)*m_dZs )
            end do
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Initialize sediment gas profiles
   !
   !------------------------------------------------------------------------------
   subroutine InitializeDiagenesisStateVariables()
      implicit none
      real(r8) :: ps, temp, Porosity, pH
      integer :: ii, icol

      pH = lake_info%pH
      do icol = 1, NSCOL, 1
         ! invalid sediment column
         if (COUNT(m_soilColInd==icol)==0) then
            m_sedSubCon(:,icol,:) = 0._r8
            m_sedSRPCon(icol) = 0._r8
            m_sedEHL(icol,:) = 0._r8
            cycle
         end if
         do ii = 1, SED_LAYER+1, 1
            temp = m_sedTemp(icol,ii)
            Porosity = m_sedpor(icol,ii)
            ps = Xch4 * P0    ! CH4 
            m_sedSubCon(Wch4,icol,ii) = Porosity*CalcEQConc(Wch4,temp,pH,ps)
            ps = Xco2 * P0   ! CO2
            m_sedSubCon(Wco2,icol,ii) = Porosity*CalcEQConc(Wco2,temp,pH,ps)
            m_sedSubCon(Wn2,icol,ii) = 0.0_r8   ! initially no n2
            m_sedSubCon(Wo2,icol,ii) = 0.0_r8   ! initially no o2
            m_sedEHL(icol,ii) = -200.0_r8       ! fully anoxic 
         end do
         m_sedSRPCon(icol) = lake_info%srp / MasP
         m_sedSubCon(Wsrp,icol,:) = CalcSedMinPConc(lake_info%srp, &
            lake_info%sedFe)
      end do
      m_btmbflux = 0.0_r8
      m_btmdflux = 0.0_r8

      ! initialize three carbon pools
      lagged_ice = m_sedIce / m_sedpor
      call ConstructPasCarbonPool()
      call ConstructActCarbonPool()  
      call ConstructOldCarbonPool()
   end subroutine

end module diagenesis_mod
