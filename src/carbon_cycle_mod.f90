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
   public :: CarbonCycleEquation
   public :: GetProductionRates, GetRespirationRates
   public :: PhytoplanktonDynamics
   public :: ZooplanktonDynamics
   public :: BedVegetationDynamics
   ! memory cache for Runge-Kutta
   type(RungeKuttaCache2D) :: mem_sub
   ! gas transfer velocity (m/s)
   real(r8), allocatable :: Kg(:)
   ! particulate settling velocity (m/s)
   real(r8), allocatable :: Vsettl(:,:)
   ! heterotrophic respiration (mol/m3/s)
   real(r8), allocatable :: rHR(:)
   real(r8), allocatable :: rHRbase(:)
   ! carbon fixation rate (s-1)
   real(r8), allocatable :: rGPP(:,:)
   real(r8), allocatable :: rGPPbase(:,:)
   ! phytoplankton metabolic loss (s-1)
   real(r8), allocatable :: rPhytoResp(:,:)
   real(r8), allocatable :: rPhytoRespbase(:,:)
   real(r8), allocatable :: rPhytoMort(:,:)
   ! sinking-induced mortality (s-1)
   real(r8), allocatable :: rPhytoDep(:,:)
   ! grazing-induced mortality (s-1)
   real(r8), allocatable :: rPhytoGraze(:,:)
   ! summary of phytoplanton C fluxes (gC/m3/s)
   real(r8), allocatable :: rsumGPP(:)
   real(r8), allocatable :: rsumResp(:)
   ! zooplankton metabolism (s-1)
   real(r8), allocatable :: rZooMort(:)
   real(r8), allocatable :: rZooResp(:)
   ! CH4 oxidation (mol/m3/s)
   real(r8), allocatable :: rOCH4(:)
   real(r8), allocatable :: rOCH4base(:)
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
   ! temporal zooplankton biomass (gC/m3)
   real(r8), allocatable :: zooC_new(:)
   ! equilibrium gas concentration (mol/m3)
   real(r8), allocatable :: EQConc(:)
   ! Chl : C ratio profile (g Chl gC-1)
   real(r8), allocatable :: rChl2C(:,:)
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
      allocate(Vsettl(NPTK,WATER_LAYER+1))
      allocate(rDYN(NWSUB,WATER_LAYER+1))
      allocate(rPDYN(NPTK,WATER_LAYER+1))
      allocate(rVegProd(WATER_LAYER+1))
      allocate(rVegResp(WATER_LAYER+1))
      allocate(rVegMort(WATER_LAYER+1))
      allocate(rGPP(NPTK,WATER_LAYER+1))
      allocate(rGPPbase(NPTK,WATER_LAYER+1))
      allocate(rHR(WATER_LAYER+1))
      allocate(rHRbase(WATER_LAYER+1))
      allocate(rPhytoResp(NPTK,WATER_LAYER+1))   
      allocate(rPhytoRespbase(NPTK,WATER_LAYER+1))
      allocate(rPhytoMort(NPTK,WATER_LAYER+1))
      allocate(rPhytoDep(NPTK,WATER_LAYER+1))
      allocate(rPhytoGraze(NPTK,WATER_LAYER+1))
      allocate(rsumGPP(WATER_LAYER+1))
      allocate(rsumResp(WATER_LAYER+1))
      allocate(rZooMort(WATER_LAYER+1))
      allocate(rZooResp(WATER_LAYER+1))
      allocate(rOCH4(WATER_LAYER+1))
      allocate(rOCH4base(WATER_LAYER+1))
      allocate(rOMP(WATER_LAYER+1))
      allocate(bmotion(NPTK,WATER_LAYER+1))
      allocate(zooC_new(WATER_LAYER+1))
      allocate(EQConc(NGAS))
      allocate(rChl2C(NPTK,WATER_LAYER+1))
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
      rGPPbase = 0.0_r8
      rHR = 0.0_r8
      rHRbase = 0.0_r8
      rPhytoResp = 0.0_r8
      rPhytoRespbase = 0.0_r8
      rPhytoMort = 0.0_r8
      rPhytoDep = 0.0_r8
      rPhytoGraze = 0.0_r8
      rsumGPP = 0.0_r8
      rsumResp = 0.0_r8
      rZooMort = 0.0_r8
      rZooResp = 0.0_r8
      rOCH4 = 0.0_r8
      rOCH4base = 0.0_r8
      rOMP = 0.0_r8
      rDYN = 0.0_r8
      rPDYN = 0.0_r8
      Kg = 0.0_r8
      Vsettl = 0.0_r8
      bmotion = 0.0_r8
      zooC_new = 0.0_r8
      rChl2C = 0.02_r8 
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
      deallocate(rGPPbase)
      deallocate(rHR)
      deallocate(rHRbase)
      deallocate(rPhytoResp)
      deallocate(rPhytoRespbase)
      deallocate(rPhytoMort)
      deallocate(rPhytoDep)
      deallocate(rPhytoGraze)
      deallocate(rsumGPP)
      deallocate(rsumResp)
      deallocate(rZooMort)
      deallocate(rZooResp)
      deallocate(rOCH4)
      deallocate(rOCH4base)
      deallocate(rOMP)
      deallocate(bmotion)
      deallocate(zooC_new)
      deallocate(EQConc)
      deallocate(rChl2C)
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
      integer :: top, gas, ii

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
         m_waterSubCon(Wsrp,top) = m_surfData%srp
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
      real(r8) :: redrPOC(NPTK), redrZOC
      integer :: top, bottom, ii
      integer :: icol, indx, indx0
 
      top = m_lakeWaterTopIndex
      bottom = WATER_LAYER + 1

      call ConvectiveMixing()
      call AdjustNegativeConcentrations()

      ! Count deposited terrigenous carbon
      do icol = 1, NSCOL, 1
         indx = COUNT(m_soilColInd<=icol) 
         indx0 = COUNT(m_soilColInd<=icol-1) + 1
         if (indx0<=indx) then
            m_trDepAtCarb(icol) = fTrActC*lake_info%cdep/SECOND_OF_YEAR 
         else
            m_trDepAtCarb(icol) = 0._r8
         end if
      end do    

      ! redistribute phytoplankton in winter
      if (top<=bottom) then
         redrPOC = 0._r8
         redrZOC = 0._r8
         do ii = 1, top-1, 1
            redrPOC = redrPOC + m_phytoC(:,ii)*m_dZw(ii)*m_Az(ii)
            redrZOC = redrZOC + m_zooC(ii)*m_dZw(ii)*m_Az(ii)
         end do
         
         m_phytoC(:,top) = m_phytoC(:,top) + redrPOC/m_dZw(top)/m_Az(top)
         m_phytoC(:,1:top-1) = 0.0_r8

         m_zooC(top) = m_zooC(top) + redrZOC/m_dZw(top)/m_Az(top)
         m_zooC(1:top-1) = 0.0_r8
      end if
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
      real(r8) :: deadPOC
      integer :: kk, ii, jj, icol
      integer :: top, bottom

      top = m_lakeWaterTopIndex
      bottom = WATER_LAYER+1

      m_aqDepAtCarb = 0.0_r8  ! initialize
      bmotion = 0.0_r8
      do kk = 1, NPTK, 1
         ! Count phytoplankton deposition
         do ii = top, bottom, 1
            icol = m_soilColInd(ii)
            if (m_phytoC(kk,ii)>phytoCmin) then
               deadPOC = rPhytoDep(kk,ii) * m_phytoC(kk,ii) * m_dZw(ii)
            else
               deadPOC = 0._r8
            end if
            do jj = top, ii, 1
               if (m_phytoC(kk,jj)>phytoCmin) then
                  deadPOC = deadPOC + (rPhytoMort(kk,jj) + (1._r8 - fDAssZoo)* &
                     rPhytoGraze(kk,jj)) * m_phytoC(kk,jj) * m_dZw(jj)
               end if
            end do
            m_aqDepAtCarb(icol) = m_aqDepAtCarb(icol) + fPhyActC * deadPOC * &
               (m_dAz(ii)/m_colAz(icol))
         end do

         ! biomass change due to metabolism
         do ii = top, bottom, 1
            blamda = -rPDYN(kk,ii)
            m_phytoC(kk,ii) = m_phytoC(kk,ii) * exp(-blamda*dt)
            ! phytoplankton movement
            if (Vsettl(kk,ii)>=0.0) then
               bmotion(kk,ii) = min(Vsettl(kk,ii),m_dZw(ii)/dt) * &
                  m_phytoC(kk,ii)
            else
               bmotion(kk,ii) = max(Vsettl(kk,ii),-m_dZw(ii)/dt) * &
                  m_phytoC(kk,ii)
            end if
         end do
         
         ! vertical movement including actively swim
         do ii = top, bottom, 1
            if (ii==top) then
               aa = 0.0_r8
               if (ii<bottom) then
                  bb = -abs(bmotion(kk,ii))*m_Az(ii+1)
                  cc = max(-bmotion(kk,ii+1),0.0_r8)*m_Az(ii+1)
               else
                  bb = 0.0_r8
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
            m_phytoC(kk,ii) = max( m_phytoC(kk,ii) + (aa+bb+cc)/ &
               m_dZw(ii)/m_Az(ii)*dt, phytoCmin )
         end do
      end do

      ! update phenology state
      call UpdatePhytoPhenology()

   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate the dynamics of zooplankton.
   !
   !------------------------------------------------------------------------------
   subroutine ZooplanktonDynamics(dt)
      real(r8), intent(in) :: dt
      real(r8) :: deadZOC
      integer :: ii, jj, icol
      integer :: top, bottom

      top = m_lakeWaterTopIndex
      bottom = WATER_LAYER+1

      ! zooplankton metabolism
      zooC_new = m_zooC
      do ii = top, bottom, 1
         zooC_new(ii) = m_zooC(ii) + ( -(rZooMort(ii)+rZooResp(ii))*m_zooC(ii) & 
            + fDAssZoo*sum(rPhytoGraze(:,ii)*m_phytoC(:,ii)) ) * dt
         zooC_new(ii) = max(zooC_new(ii), zooCmin)
      end do

      ! Count deposited zooplankton due to mortality and overgrowth 
      do ii = top, bottom, 1
         icol = m_soilColInd(ii)
         deadZOC = 0._r8
         do jj = top, ii, 1
            if (m_zooC(jj)>zooCmin) then
               deadZOC = deadZOC + min(rZooMort(jj),1._r8/dt) * m_zooC(jj) * &
                  m_dZw(jj)
            end if
            deadZOC = deadZOC + max(zooC_new(jj)-zooC_cap, 0._r8) / dt * &
               m_dZw(jj)
         end do
         m_aqDepAtCarb(icol) = m_aqDepAtCarb(icol) + fZooActC * deadZOC * &
            (m_dAz(ii)/m_colAz(icol))
      end do

      ! restrict by carrying capacity 
      do ii = top, bottom, 1
         m_zooC(ii) = min(zooC_new(ii), zooC_cap)
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
      real(r8) :: shootfrac, bedVegDW_next
      integer :: ii, icol

      ! initialize
      m_vegDepAtCarb = 0._r8

      do ii = 1, WATER_LAYER+1, 1
         if (m_dAz(ii)>e8) then
            ! update metabolic rates
            call MacrophyteMetabolism(m_bedVegDW(ii), m_Ipar(ii), &
                     m_waterTemp(ii), tonfset(ii), vegPhenol(ii), shootfrac, &
                     rVegProd(ii), rVegResp(ii), rVegMort(ii))

            ! count macrophyte mortality induced deposition
            icol = m_soilColInd(ii)
            m_vegDepAtCarb(icol) = m_vegDepAtCarb(icol) + fVegActC * &
               cCPerDW * rVegMort(ii) * m_dAz(ii) / m_colAz(icol)

            ! update biomass
            bedVegDW_next = m_bedVegDW(ii) + (rVegProd(ii) - rVegResp(ii) - &
                     rVegMort(ii)) * dt
            if (m_bedVegDW(ii)<=vegDWmin .and. bedVegDW_next<vegDWmin) then
               m_bedVegDW(ii) = vegDWmin
               rVegProd(ii) = 0._r8
               rVegResp(ii) = 0._r8
               rVegMort(ii) = 0._r8
            else
               m_bedVegDW(ii) = bedVegDW_next
            end if

            ! update P uptake
            m_vegPUptake(ii) = fVegActC * cCPerDW * rVegProd(ii) / YC2P_POM 

            ! update bed vegeation cover
            m_fcovBedVeg(ii) = min(cCovSpVeg*shootfrac*m_bedVegDW(ii), 1._r8)
         end if
      end do

      ! update phenology state
      call UpdateVegPhenology(dt)
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

      ! correct invalid phytoplankton biomass
      where (m_phytoC<phytoCmin) m_phytoC = phytoCmin

      ! correct invalid zooplankton biomass
      where (m_zooC<zooCmin) m_zooC = zooCmin

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
      real(r8) :: temp, lensink, fsed
      real(r8) :: depth, real_maxdepth, och4r
      real(r8) :: Dsrp, Dco2, Do2, Dch4
      real(r8) :: tmpGPP, tmpResp, dzfactor 
      real(r8) :: phytoC(NPTK), zooC
      integer :: ii, jj, top, btm
      integer :: trophic

      top = m_lakeWaterTopIndex
      btm = WATER_LAYER + 1 
      trophic = lake_info%trophic
      real_maxdepth = lake_info%maxdepth - m_surfData%dzsurf

      rGPP(:,1:top-1) = 0.0_r8
      rGPPbase(:,1:top-1) = 0.0_r8
      rPhytoResp(:,1:top-1) = 0.0_r8
      rPhytoRespbase(:,1:top-1) = 0.0_r8
      rPhytoMort(:,1:top-1) = 0.0_r8
      rPhytoDep(:,1:top-1) = 0.0_r8
      rPhytoGraze(:,1:top-1) = 0.0_r8
      rZooMort(1:top-1) = 0.0_r8
      rZooResp(1:top-1) = 0.0_r8
      rHR(1:top-1) = 0.0_r8
      rHRbase(1:top-1) = 0.0_r8
      rOCH4(1:top-1) = 0.0_r8
      rOCH4base(1:top-1) = 0.0_r8
      rOMP(1:top-1) = 0.0_r8
      rPDYN(:,1:top-1) = 0.0_r8
      rsumGPP(1:top-1) = 0.0_r8
      rsumResp(1:top-1) = 0.0_r8
      do ii = top, WATER_LAYER+1, 1
         temp = max(m_waterTemp(ii), T0)
         depth = sum(m_dZw(1:ii-1))
         lensink = 0.5_r8 * m_dZw(ii)
         fsed = m_dAz(ii) / m_Az(ii)
         Dch4 = m_waterSubCon(Wch4,ii)
         Do2 = m_waterSubCon(Wo2,ii)
         Dco2 = m_waterSubCon(Wco2,ii)
         Dsrp = m_waterSubCon(Wsrp,ii)
         phytoC = m_phytoC(:,ii)
         zooC = m_zooC(ii)

         if (isUpdate) then
            ! food web metabolism
            call Photosynthesis(temp, Dco2, Dsrp, m_Ipar(ii), rGPPbase(:,ii))
            call AutotrophicR(temp, rPhytoRespbase(:,ii), rPhytoMort(:,ii))
            call GrazingMortality(temp, Do2, zooC, phytoC, rPhytoGraze(:,ii))
            call DepositionMortality(Vsettl(:,ii), lensink, fsed, rPhytoDep(:,ii)) 
            call ZooplanktonLoss(temp, Do2, rZooResp(ii), rZooMort(ii)) 
            ! microbial mineralization
            call HeterotrophicR(temp, Do2, trophic, rHRbase(ii))
            call Methanotrophy(Dch4, Do2, temp, rOCH4base(ii))
            call OxicMethanogenesis(rGPPbase(:,ii)*phytoC, rOMP(ii))
         end if

         rHR(ii) = min(rHRbase(ii), Do2)

         if (Hydro_Module) then
            if (m_surfData%dzsurf>e8) then
               if (depth>real_maxdepth) then
                  och4r = 0._r8
               else if (real_maxdepth < depth + m_dZw(ii)) then
                  och4r = (real_maxdepth - depth) / m_dZw(ii)
               else
                  och4r = 1._r8
               end if
            else
               if (ii==WATER_LAYER+1) then
                  och4r = 1._r8 - m_surfData%dzsurf/m_dZw(ii)
               else
                  och4r = 1._r8
               end if
            end if 
            rOCH4(ii) = min(och4r*rOCH4base(ii), 0.5*min(Do2,Dch4))
         else
            rOCH4(ii) = min(rOCH4base(ii), 0.5*min(Do2,Dch4))
         end if

         ! carbon-based metabolism rates
         rsumGPP(ii) = sum( rGPPbase(:,ii)*phytoC )
         tmpGPP = min(rsumGPP(ii), min(Dsrp*YC2P_POM,Dco2)*MasC)
         if (rsumGPP(ii)>1.d-15) then
            rGPP(:,ii) = rGPPbase(:,ii) * tmpGPP / rsumGPP(ii)
            rsumGPP(ii) = tmpGPP
         else
            rGPP(:,ii) = rGPPbase(:,ii)
         end if

         rsumResp(ii) = sum( rPhytoRespbase(:,ii)*phytoC )
         tmpResp = min(rsumResp(ii), Do2*MasC)
         if (rsumResp(ii)>1.d-15) then
            rPhytoResp(:,ii) = rPhytoRespbase(:,ii) * tmpResp / rsumResp(ii)
            rsumResp(ii) = tmpResp
         else
            rPhytoResp(:,ii) = rPhytoRespbase(:,ii)
         end if

         ! dynamics rates of phytoplankton
         rPDYN(:,ii) = rGPP(:,ii) - rPhytoResp(:,ii) - rPhytoMort(:,ii) - &
            rPhytoDep(:,ii) - rPhytoGraze(:,ii)
      end do

      ! dynamics rates of dissolved substances
      rDYN(Wn2,:) = m_gasExch(Wn2,:)
      rDYN(Wo2,:) = (rsumGPP - rsumResp)/MasC - rZooResp*m_zooC/MasC - rHR - &
                     2.0*rOCH4 + m_gasExch(Wo2,:)
      rDYN(Wco2,:) = (rsumResp - rsumGPP)/MasC + rZooResp*m_zooC/MasC + rHR + &
                     rOCH4 - rOMP + m_gasExch(Wco2,:)
      rDYN(Wch4,:) = -rOCH4 + rOMP + m_gasExch(Wch4,:)
      rDYN(Wsrp,:) = (rsumResp - rsumGPP + rZooResp*m_zooC)/MasC/YC2P_POM
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update POM settling velocity (Stokes law)
   !          Phytoplankton motion in SBL follows MacIntyre (1998).
   !
   !------------------------------------------------------------------------------
   subroutine UpdatePOMSettlingVelocity()
      implicit none
      real(r8) :: fsrp(NPTK), Dsrp, ipar_z 
      integer :: ii, top, bottom

      top = m_lakeWaterTopIndex
      bottom = WATER_LAYER + 1
      Vsettl(:,1:top-1) = 0.0_r8
      do ii = top, bottom, 1
         Dsrp = m_waterSubCon(Wsrp,ii)
         ipar_z = 1.d6 * m_Ipar(max(ii-1,1))
         if (ii<=m_mixTopIndex) then
            Vsettl(:,ii) = 0.0_r8
         else
            ! nutrient- or light-driven movement
            fsrp = Dsrp / (PAR_Ksrp + Dsrp)
            if (ipar_z>ipar_crit) then
               where (fsrp<fsrp_vmdown)
                  Vsettl(:,ii) = Vswim
               elsewhere (fsrp>fsrp_vmup)
                  Vsettl(:,ii) = -Vswim
               elsewhere
                  Vsettl(:,ii) = 0.0_r8
               end where
               Vsettl(diatom_ppk,ii) = Vs0(diatom_ppk)
            else
               where (DsPhyto>m_wrho(ii))
                  Vsettl(:,ii) = Vs0
               elsewhere   ! within high-density waters
                  Vsettl(:,ii) = 0.1_r8 * Vs0
               end where
            end if
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
         Ceq_srp = m_surfData%srp
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
         call CalcChl2CRatio(m_Ipar(1:WATER_LAYER+1), m_waterTemp, &
               m_waterSubCon(Wsrp,:), m_waterSubCon(Wco2,:), top, rChl2C)
      end if
      m_chla = sum(m_phytoC*rChl2C,1) 
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
      npp = gpp - sum(rsumResp*m_dZw*m_Az) / m_Az(1)
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

      ! convective mixing of phytoplankton and zooplankton
      mindx = m_mixTopIndex
      if (mindx>top) then
         do kk = 1, NPTK, 1
            biomas_avg = sum(m_phytoC(kk,top:mindx)*m_dZw(top:mindx)* &
               m_Az(top:mindx)) / sum(m_dZw(top:mindx)*m_Az(top:mindx))
            m_phytoC(kk,top:mindx) = biomas_avg
         end do
         biomas_avg = sum(m_zooC(top:mindx)*m_dZw(top:mindx)* &
               m_Az(top:mindx)) / sum(m_dZw(top:mindx)*m_Az(top:mindx))
         m_zooC(top:mindx) = biomas_avg
      end if
      !mindx = m_mixBotIndex
      !if (m_Hmix(2)>e8 .and. m_Hice<e8) then
      !   do kk = 1, NPTK, 1
      !      biomas_avg = sum(m_phytoC(kk,mindx:btm)*m_dZw(mindx:btm)* &
      !         m_Az(mindx:btm)) / sum(m_dZw(mindx:btm)*m_Az(mindx:btm))
      !      m_phytoC(kk,mindx:btm) = biomas_avg
      !   end do
      !   biomas_avg = sum(m_zooC(mindx:btm)*m_dZw(mindx:btm)* &
      !      m_Az(mindx:btm)) / sum(m_dZw(mindx:btm)*m_Az(mindx:btm))
      !   m_zooC(mindx:btm) = biomas_avg
      !end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: update phenology state of submerged macrophytes.
   !
   !------------------------------------------------------------------------------
   subroutine UpdateVegPhenology(dt)
      implicit none
      real(r8), intent(in) :: dt
      integer :: ii, kk

      ! macrophyte phenology
      do ii = 1, WATER_LAYER+1, 1
         if (m_dAz(ii)<e8) then
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
   end subroutine

   subroutine UpdatePhytoPhenology()
      implicit none

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
      
      m_chla = 0.005_r8
      m_phytoC = 0.5_r8
      m_zooC = 0.25_r8
      where (m_dAz>e8)
         m_bedVegDW = 10.0_r8
      elsewhere
         m_bedVegDW = 0._r8
      end where
      m_fcovBedVeg = 0.0_r8
      m_vegPUptake = 0.0_r8

      m_gasExch = 0.0_r8
      m_aqDepAtCarb = 0.0_r8
      m_trDepAtCarb = 0.0_r8
      m_vegDepAtCarb = 0.0_r8
      m_topdflux = 0.0_r8
   end subroutine

end module carbon_cycle_mod 
