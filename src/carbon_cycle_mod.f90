module carbon_cycle_mod
!---------------------------------------------------------------------------------
! Purpose: this module governs the concentration of water dissolved substances.
!          
!          References include Tang's thesis; Stepanenko et al., 2011;
!           and Hanson et al., 2004, 2011 & 2014.
!
!---------------------------------------------------------------------------------
   use shr_kind_mod,          only : r8
   use shr_ctrl_mod,          only : WATER_LAYER, NGAS, NWSUB, lake_info, &
                                     e8 => SHR_CTRL_E8
   use shr_typedef_mod,       only : RungeKuttaCache2D
   use shr_param_mod
   use phy_utilities_mod 
   use bg_utilities_mod
   use data_buffer_mod

   implicit none
   private
   public :: mem_sub
   public :: InitializeCarbonModule, DestructCarbonModule
   public :: CarbonModuleSetup, CarbonModuleCallback 
   public :: CarbonCycleEquation, PhytoplanktonDynamics 
   public :: GetProductionRates, GetRespirationRates
   public :: BedVegetationDynamics
   ! memory cache for Runge-Kutta
   type(RungeKuttaCache2D) :: mem_sub
   ! gas transfer velocity (m/s)
   real(r8), allocatable :: Kg(:)
   ! particulate settling velocity (m/s)
   real(r8), allocatable :: Vsettl(:,:)
   ! heterotrophic respiration (mol/m3/s)
   real(r8), allocatable :: rRDOC(:)
   real(r8), allocatable :: rRDOCnew(:)
   ! carbon fixation rate (s-1)
   real(r8), allocatable :: rGPP(:,:)
   ! carbon fixation rate (gC/m3/s)
   real(r8), allocatable :: rsumGPP(:)
   ! phytoplankton metabolic loss (s-1)
   real(r8), allocatable :: rLPOC(:,:)
   real(r8), allocatable :: rRLPOC(:,:)
   real(r8), allocatable :: rRLDOC(:,:)
   real(r8), allocatable :: rRLMort(:,:)
   ! phytoplankton metabolic loss (gC/m3/s)
   real(r8), allocatable :: rsumResp(:)
   real(r8), allocatable :: rsumLPOC(:)
   ! sinking-induced mortality (gC/m3/s)
   real(r8), allocatable :: rPOCMort(:,:)
   ! CH4 oxidation (mol/m3/s)
   real(r8), allocatable :: rOCH4(:)
   real(r8), allocatable :: rOCH4new(:)
   ! oxic CH4 production (mol/m3/s)
   real(r8), allocatable :: rOMP(:)
   ! dissolved elements dynamcis (mol/m3/s) 
   real(r8), allocatable :: rDYN(:,:)
   ! phytoplankton biomass dynamics (s-1)
   real(r8), allocatable :: rPDYN(:,:)
   ! submerged macrophyte biomass dynamics (gDW/m2/s)
   real(r8), allocatable :: rVegProd(:)
   real(r8), allocatable :: rVegResp(:)
   real(r8), allocatable :: rVegMort(:)
   ! phytoplankton movement (gC/m2/s)
   real(r8), allocatable :: bmotion(:,:)
   ! equilibrium gas concentration (mol/m3)
   real(r8), allocatable :: EQConc(:)
   ! bottom fluxes (mol/m2/s)
   real(r8), allocatable :: qb(:,:)
   real(r8), allocatable :: qt(:)
   ! growing-degree-day summation (d)
   real(r8), allocatable :: GDDsum(:)
   ! onset/offset timer (d)
   real(r8), allocatable :: tonfset(:)
   ! growing-degree-day summation start
   integer, allocatable :: flagGDDsum(:)
   ! submerged macrophyte phenology state
   integer, allocatable :: vegPhenol(:)
   ! phytoplankton phenology state
   integer :: phytoPhenol
   logical :: checkPhenol

contains
   subroutine InitializeCarbonModule()
      implicit none

      allocate(Kg(NGAS))
      allocate(Vsettl(NPOC,WATER_LAYER+1))
      allocate(rDYN(NWSUB,WATER_LAYER+1))
      allocate(rPDYN(NPOC,WATER_LAYER+1))
      allocate(rVegProd(WATER_LAYER+1))
      allocate(rVegResp(WATER_LAYER+1))
      allocate(rVegMort(WATER_LAYER+1))
      allocate(rGPP(NPOC,WATER_LAYER+1))
      allocate(rsumGPP(WATER_LAYER+1))
      allocate(rRDOC(WATER_LAYER+1))
      allocate(rRDOCnew(WATER_LAYER+1))
      allocate(rLPOC(NPOC,WATER_LAYER+1))
      allocate(rRLPOC(NPOC,WATER_LAYER+1))   
      allocate(rRLDOC(NPOC,WATER_LAYER+1))
      allocate(rRLMort(NPOC,WATER_LAYER+1))
      allocate(rsumResp(WATER_LAYER+1))
      allocate(rsumLPOC(WATER_LAYER+1))
      allocate(rPOCMort(NPOC,WATER_LAYER+1))
      allocate(rOCH4(WATER_LAYER+1))
      allocate(rOCH4new(WATER_LAYER+1))
      allocate(rOMP(WATER_LAYER+1))
      allocate(bmotion(NPOC,WATER_LAYER+1))
      allocate(EQConc(NGAS))
      allocate(qb(NWSUB,WATER_LAYER+1))
      allocate(qt(NWSUB))
      allocate(GDDsum(WATER_LAYER+1))
      allocate(tonfset(WATER_LAYER+1))
      allocate(flagGDDsum(WATER_LAYER+1))
      allocate(vegPhenol(WATER_LAYER+1))
      allocate(mem_sub%K1(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%K2(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%K3(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%K4(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%K5(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%K6(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%nxt4th(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%nxt5th(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%interim(NWSUB,WATER_LAYER+1))
      allocate(mem_sub%rerr(NWSUB,WATER_LAYER+1))

      call InitializeCarbonStateVariables()
      rVegProd = 0.0_r8
      rVegResp = 0.0_r8
      rVegMort = 0.0_r8
      rGPP = 0.0_r8
      rsumGPP = 0.0_r8
      rRDOC = 0.0_r8
      rRDOCnew = 0.0_r8
      rLPOC = 0.0_r8
      rRLPOC = 0.0_r8
      rRLDOC = 0.0_r8
      rRLMort = 0.0_r8
      rsumResp = 0.0_r8
      rsumLPOC = 0.0_r8
      rPOCMort = 0.0_r8
      rOCH4 = 0.0_r8
      rOCH4new = 0.0_r8
      rOMP = 0.0_r8
      rDYN = 0.0_r8
      rPDYN = 0.0_r8
      Kg = 0.0_r8
      Vsettl = 0.0_r8
      bmotion = 0.0_r8
      qb = 0.0_r8
      qt = 0.0_r8
      GDDsum = 0.0_r8
      tonfset = 0.0_r8
      flagGDDsum = 0
      vegPhenol = grow_season
      phytoPhenol = grow_season
      checkPhenol = .False.
   end subroutine

   subroutine DestructCarbonModule()
      implicit none

      deallocate(Kg)
      deallocate(Vsettl)
      deallocate(rDYN)
      deallocate(rPDYN)
      deallocate(rVegProd)
      deallocate(rVegResp)
      deallocate(rVegMort)
      deallocate(rGPP)
      deallocate(rsumGPP)
      deallocate(rRDOC)
      deallocate(rRDOCnew)
      deallocate(rLPOC)
      deallocate(rRLPOC)
      deallocate(rRLDOC)
      deallocate(rRLMort)
      deallocate(rsumResp)
      deallocate(rsumLPOC)
      deallocate(rPOCMort)
      deallocate(rOCH4)
      deallocate(rOCH4new)
      deallocate(rOMP)
      deallocate(bmotion)
      deallocate(EQConc)
      deallocate(qb)
      deallocate(qt)
      deallocate(GDDsum)
      deallocate(tonfset)
      deallocate(flagGDDsum)
      deallocate(vegPhenol)
      deallocate(mem_sub%K1)
      deallocate(mem_sub%K2)
      deallocate(mem_sub%K3)
      deallocate(mem_sub%K4)
      deallocate(mem_sub%K5)
      deallocate(mem_sub%K6)
      deallocate(mem_sub%nxt4th)
      deallocate(mem_sub%nxt5th)
      deallocate(mem_sub%interim)
      deallocate(mem_sub%rerr)
   end subroutine

   subroutine CarbonModuleSetup(isHourNode, isSubHourNode, hindx)
      implicit none
      logical, intent(in) :: isHourNode
      logical, intent(in) :: isSubHourNode
      integer(i8), intent(in) :: hindx
      real(r8) :: pH, temp, pressure
      real(r8) :: mixing_ratio(NGAS)
      integer :: top, gas

      top = m_lakeWaterTopIndex

      ! update EQConc
      if (isSubHourNode) then
         pH = lake_info%pH
         if (top<=WATER_LAYER+1) then
            temp = m_waterTemp(top)
         else
            temp = T0
         end if
         ! maximum saturation of O2 is 50% [Holgerson et al., 2016]
         mixing_ratio = (/Xn2, Xo2, 1.0d-6*m_radPars%qCO2, Xch4/)
         do gas = Wn2, Wch4, 1
            pressure = mixing_ratio(gas) * m_surfData%pressure
            EQConc(gas) = CalcEQConc(gas, temp, pH, pressure)
         end do          
      end if
      call UpdateGasTransferVelocity(isSubHourNode)

      ! supply O2 from surface water and groundwater
      if (m_Hice>e8 .and. m_Hice<=lake_info%bfdep) then
         m_waterSubCon(Wo2,top) = EQConc(Wo2)
      end if
      ! assume dissolved N2 always replete
      !if (m_Hice<e8) then
      !   m_waterSubCon(Wn2,top:bottom) = EQConc(Wn2) 
      !end if
      ! fix surface SRP
      if (m_Hice<e8) then
         m_waterSubCon(Wsrp,top) = m_surfData%srp / MasP 
      end if
      ! No dissolved chemicals in ice
      m_waterSubCon(:,1:top-1) = 0._r8

      ! update diffusivity for bottom mixing layer
      if (m_Hmix(2)>5._r8 .and. m_Hice<e8) then
         m_Kv(m_mixBotIndex:WATER_LAYER+1) = 1.0d-5
      end if

      call CalcSurfaceExchange()
      call CalcBottomExchange()
      call UpdateChla2POCRatio(isHourNode, hindx)
      call UpdatePOMSettlingVelocity()
      call UpdateCarbonCycleRates(isSubHourNode)
   end subroutine

   subroutine CarbonModuleCallback(dt)
      implicit none
      real(r8), intent(in) :: dt
      real(r8) :: redrPOC(NPOC), sumAz
      real(r8) :: deadPOC(NPOC)
      integer :: top, bottom, ii, kk
      integer :: icol, indx, indx0
 
      top = m_lakeWaterTopIndex
      bottom = WATER_LAYER + 1

      call ConvectiveMixing()
      call AdjustNegativeConcentrations()

      ! Count deposited phytoplankton for each sediment column
      do icol = 1, NSCOL, 1
         indx = COUNT(m_soilColInd<=icol) 
         indx0 = COUNT(m_soilColInd<=icol-1) + 1
         if (indx0<=indx) then
            sumAz = sum(m_dAz(indx0:indx))
            ! dead phytoplankton
            m_aqDepAtCarb(icol) = 0.0_r8
            do ii = top, indx, 1
               deadPOC = (rRLMort(:,ii) + rPOCMort(:,ii)) * m_waterPOC(:,ii)
               if (ii<indx0) then
                  m_aqDepAtCarb(icol) = m_aqDepAtCarb(icol) + &
                     fPhyActC * sum(deadPOC) * m_dZw(ii)
               else
                  m_aqDepAtCarb(icol) = m_aqDepAtCarb(icol) + &
                     fPhyActC * sum(deadPOC) * m_dZw(ii) * &
                     sum(m_dAz(ii:indx)) / sumAz
               end if
            end do
            ! submerged macrophyte mortality
            m_vegDepAtCarb(icol) = 0._r8
            do ii = indx0, indx, 1
               m_vegDepAtCarb(icol) = m_vegDepAtCarb(icol) + fVegActC * &
                  cCPerDW * rVegMort(ii) * m_dAz(ii) / sumAz      
            end do
            ! terrigenous source
            m_trDepAtCarb(icol) = fTrActC*lake_info%cdep/SECOND_OF_YEAR 
         end if
      end do    
      ! redistribute phytoplankton in winter
      if (top<=bottom) then
         redrPOC = 0.0_r8
         do ii = 1, top-1, 1
            redrPOC = redrPOC + m_waterPOC(:,ii)*m_dZw(ii)*m_Az(ii)
         end do
         m_waterPOC(:,top) = m_waterPOC(:,top) + redrPOC/m_dZw(top)/m_Az(top)
         m_waterPOC(:,1:top-1) = 0.0_r8
      end if

      call UpdatePhenologyState(dt)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Discretize substance transport equation by Finite Difference Method.
   !          This sub-routine returns time derivatives of dissolved substance
   !          concentration.
   !
   !------------------------------------------------------------------------------
   subroutine CarbonCycleEquation(con, dcon)
      implicit none
      real(r8), intent(in) :: con(NWSUB,WATER_LAYER+1)      ! Unit: mol/m3
      real(r8), intent(out) :: dcon(NWSUB,WATER_LAYER+1)    ! Unit: mol/(m3*s)
      real(r8), dimension(NWSUB) :: dCa, dCb
      real(r8) :: aa, bb, Az1, Az2
      integer :: ii, top

      top = m_lakeWaterTopIndex
      do ii = 1, WATER_LAYER+1, 1
         ! frozen layers
         if (ii<top) then
            dcon(:,ii) = 0._r8
            cycle
         end if
         ! unfrozen layers
         if(ii==1) then
            aa = 0.5 * (m_Kv(ii) + m_Kv(ii+1))
            Az1 = m_Az(ii+1)
            dca = (con(:,ii+1) - con(:,ii)) / (m_Zw(ii+1) - m_Zw(ii))
            dcon(:,ii) = (aa*Az1*dca + qb(:,ii)*m_dAz(ii) - qt*m_Az(ii)) / &
                  m_dZw(ii) / m_Az(ii) + rDYN(:,ii)
         else if(ii==WATER_LAYER+1) then
            if (ii-1>=top) then
               bb = 0.5 * (m_Kv(ii-1) + m_Kv(ii))
            else
               bb = 0._r8
            end if
            Az2 = m_Az(ii) 
            dcb = (con(:,ii) - con(:,ii-1)) / (m_Zw(ii) - m_Zw(ii-1))
            dcon(:,ii) = (qb(:,ii)*m_dAz(ii) - bb*Az2*dcb) / m_dZw(ii) / &
                  m_Az(ii) + rDYN(:,ii)
         else
            aa = 0.5 * (m_Kv(ii) + m_Kv(ii+1))
            if (ii-1>=top) then
               bb = 0.5 * (m_Kv(ii) + m_Kv(ii-1))
            else
               bb = 0._r8
            end if
            Az1 = m_Az(ii+1) 
            Az2 = m_Az(ii) 
            dca = (con(:,ii+1) - con(:,ii)) / (m_Zw(ii+1) - m_Zw(ii))
            dcb = (con(:,ii) - con(:,ii-1)) / (m_Zw(ii) - m_Zw(ii-1))
            dcon(:,ii) = (aa*Az1*dca + qb(:,ii)*m_dAz(ii) - bb*Az2*dcb) / &
                  m_dZw(ii) / m_Az(ii) + rDYN(:,ii)
         end if
      end do
      ! fixed boundary conditions
      if (m_Hice>e8 .and. m_Hice<=lake_info%bfdep) then
         dcon(Wo2,top) = 0._r8
      end if
      if (m_Hice<e8) then
         dcon(Wsrp,top) = 0._r8
      end if
      where(con<=0 .and. dcon<0) dcon = 0.0_r8
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate the dynamics of phytoplankton.
   !
   !------------------------------------------------------------------------------
   subroutine PhytoplanktonDynamics(dt)
      implicit none
      real(r8), intent(in) :: dt
      real(r8) :: blamda, aa, bb, cc
      integer :: kk, ii, top, mindx, bottom

      top = m_lakeWaterTopIndex
      mindx = m_mixBotIndex 
      bottom = WATER_LAYER+1
      bmotion = 0.0_r8
      do kk = 1, NPOC, 1
         do ii = top, bottom, 1
            ! biomass change due to metabolism
            blamda = -rPDYN(kk,ii)
            m_waterPOC(kk,ii) = m_waterPOC(kk,ii) * exp(-blamda*dt)
            ! phytoplankton movement
            if (Vsettl(kk,ii)>=0.0) then
               bmotion(kk,ii) = min(Vsettl(kk,ii),m_dZw(ii)/dt) * &
                  m_waterPOC(kk,ii)
            else
               bmotion(kk,ii) = max(Vsettl(kk,ii),-m_dZw(ii)/dt) * &
                  m_waterPOC(kk,ii)
            end if
         end do
         
         ! vertical movement including actively swim
         do ii = top, bottom, 1
            if (ii==top) then
               aa = 0.0_r8
               bb = -abs(bmotion(kk,ii))*m_Az(ii+1)
               if (ii<bottom) then
                  cc = max(-bmotion(kk,ii+1),0.0_r8)*m_Az(ii+1)
               else
                  cc = 0.0_r8
               end if
            else if (ii==bottom) then
               aa = max(bmotion(kk,ii-1),0.0_r8)*m_Az(ii)
               bb = -abs(bmotion(kk,ii))*m_Az(ii)
               cc = 0.0_r8
            else
               aa = max(bmotion(kk,ii-1),0.0_r8)*m_Az(ii)
               if (bmotion(kk,ii)>=0._r8) then
                  bb = -abs(bmotion(kk,ii))*m_Az(ii+1)
               else
                  bb = -abs(bmotion(kk,ii))*m_Az(ii)
               end if
               cc = max(-bmotion(kk,ii+1),0.0_r8)*m_Az(ii+1)
            end if
            m_waterPOC(kk,ii) = max( m_waterPOC(kk,ii) + (aa+bb+cc)/ &
               m_dZw(ii)/m_Az(ii)*dt, 0._r8 )
         end do
      end do

   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate the dynamics of submerged macrophytes.
   !
   !------------------------------------------------------------------------------
   subroutine BedVegetationDynamics(dt)
      implicit none
      real(r8), intent(in) :: dt
      real(r8) :: shootfrac
      integer :: ii

      do ii = 1, WATER_LAYER+1, 1
         ! no macrophyte survived or no bed area
         if (m_bedVegDW(ii)<e8 .or. m_dAz(ii)<e8) then
            rVegProd(ii) = 0._r8
            rVegResp(ii) = 0._r8
            rVegMort(ii) = 0._r8
            m_bedVegDW(ii) = 0._r8
            m_fcovBedVeg(ii) = 0._r8
            m_vegPUptake(ii) = 0._r8
            cycle
         end if
         ! update metabolic rates
         call MacrophyteMetabolism(m_bedVegDW(ii), m_Ipar(ii), m_waterTemp(ii), &
                  tonfset(ii), vegPhenol(ii), shootfrac, rVegProd(ii), &
                  rVegResp(ii), rVegMort(ii))
         m_bedVegDW(ii) = m_bedVegDW(ii) + (rVegProd(ii) - rVegResp(ii) - &
               rVegMort(ii)) * dt
         m_bedVegDW(ii) = max(m_bedVegDW(ii), 0._r8)
         ! update P uptake
         m_vegPUptake(ii) = fVegActC * cCPerDW * rVegProd(ii) / YC2P_POM 
         ! update bed vegeation cover
         m_fcovBedVeg(ii) = min(cCovSpVeg*shootfrac*m_bedVegDW(ii), 1._r8)
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Adjust negative dissolved gas concentrations (project negative 
   !          values to zero) ("Non-negative solutions of ODEs", Shampine, L., 
   !          2005).
   !
   !------------------------------------------------------------------------------
   subroutine AdjustNegativeConcentrations()
      implicit none
      real(r8) :: deficitCum, deficit
      integer :: ii, kk, top, btm

      top = m_lakeWaterTopIndex
      btm = WATER_LAYER + 1

      ! correct negative solute concentrations
      do kk = 1, NWSUB, 1
         deficitCum = 0.0_r8
         do ii = 1, WATER_LAYER+1, 1
            if (m_waterSubCon(kk,ii)<0._r8) then
               deficitCum = deficitCum - m_waterSubCon(kk,ii)*m_dZw(ii)*m_Az(ii)
               m_waterSubCon(kk,ii) = 0.0_r8
            end if
         end do
         ! add negative CH4 and CO2 to bottom and N2 and O2 to surface
         if (deficitCum>e8 .and. top<=btm) then
            if (kk==Wn2 .or. kk==Wo2) then
               m_waterSubCon(kk,top) = max(0._r8, m_waterSubCon(kk,top) - &
                  deficitCum/m_dZw(top)/m_Az(top)) 
            else if (kk==Wco2 .or. kk==Wch4) then
               m_waterSubCon(kk,btm) = max(0._r8, m_waterSubCon(kk,btm) - &
                  deficitCum/m_dZw(btm)/m_Az(btm))
            end if
         end if
      end do

      ! correct negative POC
      where (m_waterPOC<0._r8) m_waterPOC = 0.0_r8

      !where (m_waterSubCon<0) m_waterSubCon = 0.0_r8
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update carbon recycle reaction rates.
   !
   !          Important: for photosynthesis and CDOM photodegradation, the
   !          scalar irradiance should be used to calculate the rates.
   !
   !------------------------------------------------------------------------------
   subroutine UpdateCarbonCycleRates(isUpdate)
      implicit none
      logical, intent(in) :: isUpdate
      real(r8) :: LPOC(NPOC)
      real(r8) :: temp, dist, hhyp
      real(r8) :: Dsrp, Dco2, Do2, Dch4
      real(r8) :: aGPP, newGPP, minIpar
      real(r8) :: och4dist
      integer :: ii, jj, top, btm

      top = m_lakeWaterTopIndex
      btm = WATER_LAYER + 1 
      hhyp = sum(m_dZw(m_mixBotIndex:btm))
      minIpar = 1.d-3 * m_Ipar(1)
      och4dist = 1.0 - m_surfData%dzsurf/lake_info%maxdepth

      rGPP(:,1:top-1) = 0.0_r8
      rLPOC(:,1:top-1) = 0.0_r8
      rRLPOC(:,1:top-1) = 0.0_r8
      rRLDOC(:,1:top-1) = 0.0_r8
      rRLMort(:,1:top-1) = 0.0_r8
      rPOCMort(:,1:top-1) = 0.0_r8
      rRDOC(1:top-1) = 0.0_r8
      rRDOCnew(1:top-1) = 0.0_r8
      rOCH4(1:top-1) = 0.0_r8
      rOCH4new(1:top-1) = 0.0_r8
      rOMP(1:top-1) = 0.0_r8
      rPDYN(:,1:top-1) = 0.0_r8
      rsumGPP(1:top-1) = 0.0_r8
      rsumResp(1:top-1) = 0.0_r8
      rsumLPOC(1:top-1) = 0.0_r8
      do ii = top, WATER_LAYER+1, 1
         temp = max(m_waterTemp(ii), T0)
         dist = sum(m_dZw(ii:btm))
         Dch4 = m_waterSubCon(Wch4,ii)
         Do2 = m_waterSubCon(Wo2,ii)
         Dco2 = m_waterSubCon(Wco2,ii)
         Dsrp = m_waterSubCon(Wsrp,ii)
         LPOC = m_waterPOC(:,ii)

         if (isUpdate) then
            call Photosynthesis(temp, Dco2, Dsrp, m_Ipar(ii), rGPP(:,ii))
            call AutotrophicR(temp, Do2, rLPOC(:,ii), rRLPOC(:,ii), &
                     rRLDOC(:,ii), rRLMort(:,ii))
            call HeterotrophicR(temp, Do2, rRDOCnew(ii))
            call Methanotrophy(Dch4, Do2, temp, rOCH4new(ii))
            call OxicMethanogenesis(rGPP(:,ii)*LPOC, rOMP(ii))
            call DepositionMortality(dist, hhyp, rPOCMort(:,ii)) 
         end if
         
         ! avoid the unrealistic all death of phytoplankton
         if (m_Hice<e8 .and. m_Ipar(ii)>=minIpar) then
            do jj = 1, NPOC, 1
               if (LPOC(jj)<lake_info%refPOC) then
                  rLPOC(jj,ii) = 0.0_r8
                  rRLPOC(jj,ii) = 0.0_r8
                  rRLDOC(jj,ii) = 0.0_r8
                  rRLMort(jj,ii) = 0.0_r8
               end if
            end do
         end if

         ! carbon-based metabolism rates
         rsumGPP(ii) = sum( rGPP(:,ii)*LPOC )
         rsumResp(ii) = sum( (rRLPOC(:,ii)+rRLDOC(:,ii))*LPOC )
         rsumLPOC(ii) = sum( rLPOC(:,ii)*LPOC )

         rRDOC(ii) = min(rRDOCnew(ii), Do2)
         rOCH4(ii) = min(och4dist*rOCH4new(ii), 0.5*min(Do2,Dch4))
         newGPP = min(rsumGPP(ii), min(Dsrp*YC2P_POM,Dco2)*MasC)
         if (rsumGPP(ii)<1.d-15) then
            aGPP = 1.0_r8
         else
            aGPP = newGPP / rsumGPP(ii)
         end if
         rsumResp(ii) = min(rsumResp(ii), Do2*MasC)

         ! dynamics rates of phytoplankton
         rPDYN(:,ii) = aGPP*rGPP(:,ii) - rLPOC(:,ii) - rPOCMort(:,ii)
      end do

      ! dynamics rates of dissolved substances
      rDYN(Wn2,:) = m_gasExchange(Wn2,:)
      rDYN(Wo2,:) = (rsumGPP - rsumResp)/MasC - rRDOC - 2.0*rOCH4 + &
                     m_gasExchange(Wo2,:)
      rDYN(Wco2,:) = (rsumResp - rsumGPP)/MasC + rRDOC + rOCH4 + &
                     m_gasExchange(Wco2,:)
      rDYN(Wch4,:) = -rOCH4 + rOMP + m_gasExchange(Wch4,:)
      rDYN(Wsrp,:) = (rsumResp - rsumGPP)/MasC/YC2P_POM
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update POM settling velocity (Stokes law)
   !          Phytoplankton motion in SBL follows MacIntyre (1998).
   !
   !------------------------------------------------------------------------------
   subroutine UpdatePOMSettlingVelocity()
      implicit none
      real(r8) :: Ksrps, Ksrpl
      real(r8) :: fsrp, Dsrp, ipar_z 
      integer :: ii, top, bottom

      Ksrps = sa_params(Param_Ksrps)
      Ksrpl = sa_params(Param_Ksrpl)
      top = m_lakeWaterTopIndex
      bottom = WATER_LAYER + 1
      Vsettl(:,1:top-1) = 0.0_r8
      do ii = top, bottom, 1
         Dsrp = m_waterSubCon(Wsrp,ii)
         ipar_z = 1.d6 * m_Ipar(max(ii-1,1))
         if (ii<=m_mixTopIndex) then
            Vsettl(:,ii) = 0.0_r8 
         else if (ii<m_mixBotIndex) then
            fsrp = Dsrp / (Ksrps + Dsrp)
            if (fsrp<fsrp_vmdown .and. ipar_z>ipar_crit) then
               Vsettl(small_ppk,ii) = Vswim
            else if (fsrp>fsrp_vmup .and. ipar_z>ipar_crit) then
               Vsettl(small_ppk,ii) = -Vswim
            else
               Vsettl(small_ppk,ii) = 0.0_r8
            end if
            fsrp = Dsrp / (Ksrpl + Dsrp)
            if (fsrp<fsrp_vmdown .and. ipar_z>ipar_crit) then
               Vsettl(large_ppk,ii) = Vs0(large_ppk)
            else
               Vsettl(large_ppk,ii) = 0.0_r8
            end if
         else
            Vsettl(:,ii) = 0.0_r8 ! replaced by a mortality factor
         end if
      end do 
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update surface gas transfer velocity. 
   !
   !------------------------------------------------------------------------------
   subroutine UpdateGasTransferVelocity(isUpdate)
      implicit none
      logical, intent(in) :: isUpdate
      real(r8) :: temp, rho0, k600UW, Schmidt
      integer :: top, gas

      top = 1
      temp = m_waterTemp(top)
      rho0 = m_wrho(top)
      if (m_Hice<e8) then
         if (isUpdate) then
            k600UW = CalcPistonVelocity(lake_info, m_surfData%wind, temp, &
               rho0, m_Heff, m_Hmix(1))
            do gas = 1, NGAS, 1
               Schmidt = CalcSchmidtNumber(gas, temp)
               Kg(gas) = k600UW / sqrt(Schmidt)
            end do
         end if
      else
         Kg = 0.0_r8
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate gas flux via diffusive transfer
   !      The surface gas concentration is calculated from the top 20 cm layers
   !      (Schubert et al., 2012).
   !      When the maximum snow depth is exceeded, the accumulated gas can be
   !      emitted to the atmosphere.
   !
   !------------------------------------------------------------------------------
   subroutine CalcSurfaceExchange()
      implicit none
      real(r8) :: Ceq, pH, temp, pressure
      real(r8) :: mixing_ratio(NGAS)
      real(r8) :: Kg_srp, Ceq_srp
      integer :: gas

      if (m_Hice>e8) then
         qt = 0.0_r8
         m_topdflux = 0.0_r8
      else
         do gas = 1, NGAS, 1
            qt(gas) = Kg(gas) * (m_waterSubCon(gas,1) - EQConc(gas))
         end do
         m_topdflux = qt(1:NGAS)
         Ceq_srp = m_surfData%srp / MasP
         Kg_srp = lake_info%depth / lake_info%wrt / 8.64d4    ! m/s
         qt(Wsrp) = Kg_srp * (m_waterSubCon(Wsrp,1) - Ceq_srp) 
      end if
   end subroutine

   subroutine CalcBottomExchange()
      implicit none
      real(r8) :: flux_sed, sumV, Porosity
      real(r8) :: conc_sed, conc_wat, conc_tmp
      integer :: ii, kk, icol
      integer :: top, indx, indx0

      top = m_lakeWaterTopIndex
      qb = 0._r8
      do icol = 1, NSCOL, 1
         indx = COUNT(m_soilColInd<=icol)
         indx0 = COUNT(m_soilColInd<=icol-1) + 1
         if (indx0>indx) then
            cycle
         end if
         sumV = sum(m_dZw(indx0:indx)*m_Az(indx0:indx))
         Porosity = m_sedpor(icol,1)
         do kk = 1, NWSUB, 1
            if (kk==Wsrp) then
               conc_sed = m_sedSRPCon(icol)
            else
               conc_sed = m_sedSubCon(kk,icol,1) / Porosity
            end if
            if (icol<NSCOL) then
               conc_wat = sum(m_waterSubCon(kk,indx0:indx)* &
                  m_dZw(indx0:indx)*m_Az(indx0:indx)) / sumV
            else
               conc_wat = m_waterSubCon(kk,indx)
            end if
            flux_sed = m_Ktb(indx) * (conc_sed - conc_wat) * Porosity / DeltaD
            if (kk==Wn2 .or. kk==Wo2) then
               flux_sed = min(0._r8, flux_sed)
            else if (kk==Wco2 .or. kk==Wch4) then
               flux_sed = max(0._r8, flux_sed)
            end if
            do ii = indx0, indx, 1
               conc_tmp = m_waterSubCon(kk,ii)
               if (ii<top) then
                  qb(kk,ii) = 0._r8
               else if (flux_sed<0._r8 .and. conc_sed>conc_tmp) then
                  qb(kk,ii) = 0._r8
               else if (flux_sed>0._r8 .and. conc_sed<conc_tmp) then
                  qb(kk,ii) = 0._r8
               else
                  qb(kk,ii) = flux_sed
               end if
            end do
         end do
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Update phytoplankton Chla : C ratio 
   !
   !------------------------------------------------------------------------------
   subroutine UpdateChla2POCRatio(isHourNode, hindx)
      implicit none
      logical, intent(in) :: isHourNode
      integer(i8), intent(in) :: hindx
      integer :: ii, top
      logical :: isUpdate

      top = m_lakeWaterTopIndex
      if (isHourNode .and. mod(hindx-1,24)==12 .and. top<=WATER_LAYER+1) then
         isUpdate = .True.
      else
         isUpdate = .False.
      end if

      if (isUpdate) then
         call CalcChl2CRatio(m_Ipar, m_waterTemp, m_waterSubCon(Wsrp,:), &
               m_waterSubCon(Wco2,:), top, m_rChl2C)
      end if
      m_chla = sum(m_waterPOC*m_rChl2C,1) 
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Get GPP and NPP and other BGC rates.
   !
   !------------------------------------------------------------------------------
   subroutine GetProductionRates(gpp, npp)
      implicit none
      real(r8), intent(out) :: gpp   ! units: gC/m2/s
      real(r8), intent(out) :: npp   ! units: gC/m2/s

      gpp = sum(rsumGPP*m_dZw*m_Az) / m_Az(1)
      npp = gpp - sum(rsumLPOC*m_dZw*m_Az) / m_Az(1)
   end subroutine

   subroutine GetRespirationRates(rch4p, rch4o)
      implicit none
      real(r8), intent(out) :: rch4p(:)   ! units: mol/m3/s
      real(r8), intent(out) :: rch4o(:)   ! units: mol/m3/s

      rch4p = rOMP
      rch4o = rOCH4
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: convective mixing of dissolved substances.
   !
   !------------------------------------------------------------------------------
   subroutine ConvectiveMixing()
      implicit none
      real(r8) :: cavg, biomas_avg
      integer :: top, btm, kk, mindx

      ! convective mixing of solutes
      top = m_lakeWaterTopIndex
      btm = WATER_LAYER + 1
      mindx = m_mixTopIndex
      if (mindx>top) then
         do kk = 1, NWSUB, 1
            cavg = sum(m_waterSubCon(kk,top:mindx)*m_dZw(top:mindx)* &
               m_Az(top:mindx)) / sum(m_dZw(top:mindx)*m_Az(top:mindx))
            m_waterSubCon(kk,top:mindx) = cavg
         end do
      end if
      !mindx = m_mixBotIndex
      !if (m_Hmix(2)>e8 .and. m_Hice<e8) then
      !   do kk = 1, NWSUB, 1
      !      cavg = sum(m_waterSubCon(kk,mindx:btm)*m_dZw(mindx:btm)* &
      !         m_Az(mindx:btm)) / sum(m_dZw(mindx:btm)*m_Az(mindx:btm))
      !      m_waterSubCon(kk,mindx:btm) = cavg
      !   end do
      !end if

      ! convective mixing of phytoplankton
      mindx = m_mixTopIndex
      if (mindx>top) then
         do kk = 1, NPOC, 1
            biomas_avg = sum(m_waterPOC(kk,top:mindx)*m_dZw(top:mindx)* &
               m_Az(top:mindx)) / sum(m_dZw(top:mindx)*m_Az(top:mindx))
            m_waterPOC(kk,top:mindx) = biomas_avg
         end do
      end if
      !mindx = m_mixBotIndex
      !if (m_Hmix(2)>e8 .and. m_Hice<e8) then
      !   do kk = 1, NPOC, 1
      !      biomas_avg = sum(m_waterPOC(kk,mindx:btm)*m_dZw(mindx:btm)* &
      !         m_Az(mindx:btm)) / sum(m_dZw(mindx:btm)*m_Az(mindx:btm))
      !      m_waterPOC(kk,mindx:btm) = biomas_avg
      !   end do
      !end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update phenology state of submerged macrophytes.
   !
   !------------------------------------------------------------------------------
   subroutine UpdatePhenologyState(dt)
      implicit none
      real(r8), intent(in) :: dt
      integer :: ii, kk

      ! macrophyte phenology
      do ii = 1, WATER_LAYER+1, 1
         ! skip layers where submerged vegetation doesn't survive
         if (m_bedVegDW(ii)<e8) then
            vegPhenol(ii) = dormant_season
            cycle
         end if
         if (vegPhenol(ii)==dormant_season) then
            if (m_radPars%Latit>=0) then
               if (m_radPars%month==12 .and. m_radPars%day==21) then
                  GDDsum(ii) = 0._r8
                  flagGDDsum(ii) = 1
               else if (m_radPars%month==6 .and. m_radPars%day==21) then
                  GDDsum(ii) = 0._r8
                  flagGDDsum(ii) = 0
               end if
            else
               if (m_radPars%month==6 .and. m_radPars%day==21) then
                  GDDsum(ii) = 0._r8
                  flagGDDsum(ii) = 1
               else if (m_radPars%month==12 .and. m_radPars%day==21) then
                  GDDsum(ii) = 0._r8
                  flagGDDsum(ii) = 0
               end if
            end if
            if (flagGDDsum(ii)==1) then
               if (m_radPars%dayl>=39300._r8 .and. m_Hsnow<e8 .and. &
                     m_Hgrayice<e8 .and. m_waterTemp(ii)>T0) then
                  GDDsum(ii) = GDDsum(ii) + max(m_waterTemp(ii)-T0,3.0)* &
                     dt/8.64d4
                  if (GDDsum(ii)>GDDsum_crit) then
                     vegPhenol(ii) = onset_season
                     tonfset(ii) = cLengTrans
                     GDDsum(ii) = 0._r8
                     flagGDDsum(ii) = 0
                  end if
               end if
            end if
         else if (vegPhenol(ii)==onset_season) then
            tonfset(ii) = tonfset(ii) - dt/SECOND_OF_DAY
            if (tonfset(ii)<e8) then
               vegPhenol(ii) = grow_season
               tonfset(ii) = 0._r8
            end if
         else if (vegPhenol(ii)==grow_season) then
            if (m_radPars%Latit>=0 .and. m_radPars%month>=7 .and. &
                  m_radPars%dayl<39300._r8) then
               vegPhenol(ii) = offset_season
               tonfset(ii) = cLengTrans
            else if (m_radPars%Latit<0 .and. m_radPars%month<=6 .and. &
                  m_radPars%dayl<39300._r8) then
               vegPhenol(ii) = offset_season
               tonfset(ii) = cLengTrans
            end if
         else if (vegPhenol(ii)==offset_season) then
            tonfset(ii) = tonfset(ii) - dt/SECOND_OF_DAY
            if (tonfset(ii)<e8) then
               vegPhenol(ii) = dormant_season
               tonfset(ii) = 0._r8
            end if
         end if
      end do

      ! phytoplankton phenology
      if (phytoPhenol==dormant_season) then
         if (m_radPars%Latit>=0) then
            if (m_radPars%month==3 .and. m_radPars%day==20) then
               checkPhenol = .True.   
            end if
         else
            if (m_radPars%month==9 .and. m_radPars%day==23) then
               checkPhenol = .True.
            end if
         end if
         if (checkPhenol .and. m_Hice<e8) then
            phytoPhenol = grow_season
            checkPhenol = .False.
            ! prevent all dead
            do kk = 1, NPOC, 1
               m_waterPOC(kk,:) = m_waterPOC(kk,:) + 0.001
            end do
         end if
      else if (phytoPhenol==grow_season) then
         if (m_radPars%Latit>=0) then
            if (m_radPars%month==9 .and. m_radPars%day==23) then
               phytoPhenol = dormant_season
            end if
         else
            if (m_radPars%month==3 .and. m_radPars%day==20) then
               phytoPhenol = dormant_season
            end if
         end if 
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Initialize water gas profiles.
   !
   !------------------------------------------------------------------------------
   subroutine InitializeCarbonStateVariables()
      implicit none
      real(r8) :: pGas, temp, pH
      real(r8) :: mixing_ratio(NGAS)
      integer :: gas, ii

      pH = lake_info%pH
      mixing_ratio = (/Xn2, Xo2, Xco2, Xch4/)
      do ii = 1, WATER_LAYER+1, 1
         temp = m_waterTemp(ii)
         do gas = Wn2, Wch4, 1
            pGas = mixing_ratio(gas) * P0   ! gas partial pressure
            m_waterSubCon(gas,ii) = CalcEQConc(gas, temp, pH, pGas)
            if (ii==1) then
               EQConc(gas) = m_waterSubCon(gas,ii) 
            end if
         end do
      end do
      m_waterSubCon(Wsrp,:) = lake_info%srp / MasP
      
      m_rChl2C = 0.02_r8
      m_chla = 0.001_r8
      m_waterPOC = 0.05_r8
      m_bedVegDW = 10.0_r8 
      m_fcovBedVeg = 0.0_r8
      m_vegPUptake = 0.0_r8

      m_gasExchange = 0.0_r8
      m_aqDepAtCarb = 0.0_r8
      m_trDepAtCarb = 0.0_r8
      m_vegDepAtCarb = 0.0_r8
      m_topdflux = 0.0_r8
   end subroutine

end module carbon_cycle_mod 
