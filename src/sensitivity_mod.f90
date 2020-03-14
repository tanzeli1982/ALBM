module sensitivity_mod
!---------------------------------------------------------------------------------
! Purpose: this module governs the sensitivity analysis of the bLake model
!
!---------------------------------------------------------------------------------
   use shr_ctrl_mod
   use shr_typedef_mod
   use shr_param_mod
   use sim_coupler_mod
   use read_data_mod
   use io_utilities_mod

   implicit none
   private
   public :: RunSensitivity 
   ! # of output variables
   integer, parameter :: NOUT = 4
   ! current sample Id 
   integer :: cur_sample

contains
   subroutine RunSensitivity(arg)
      implicit none
      character(len=*), intent(in) :: arg
      real(r8), allocatable :: samples(:,:)     ! generated samples
      real(r8), allocatable :: results(:,:)
      integer :: ii, jj, cnt, minid, maxid, err
      integer :: sampleId, idx, nsample, itmp
      integer :: sample_next_range(2)
      real(r8) :: sir(NOUT)

      minid = minval(sample_range)
      maxid = maxval(sample_range)
      nsample = maxid - minid + 1

      print "(A, I0, A, I0)", 'Run samples from ', minid, ' to ', maxid

      allocate(samples(NMAXSAMPLE,NPARAM))
      allocate(results(NOUT,nsample))

      call ReadParameterSamples(NMAXSAMPLE, samples)
      
      if (minid==1) then
         call CreateSampleResultFile(NOUT, -9999.0_r8)
      end if

      do cnt = 1, nsample, 1 
         sampleId = minid + cnt - 1 
         cur_sample = sampleId
         call MonteCarloSimulation( sampleId, samples(sampleId,:), sir )
         results(:,cnt) = sir
      end do

      call WriteSampleResults(sample_range, results)

      deallocate(samples)
      deallocate(results)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: For slave ndoes, run sensitivity simulations for each sample
   !
   !------------------------------------------------------------------------------
   subroutine MonteCarloSimulation(sampleId, sample, odata)
      implicit none
      integer, intent(in) :: sampleId
      real(r8), intent(in) :: sample(NPARAM)
      real(r8), intent(out) :: odata(NOUT)
      type(SimTime) :: time, spinup
      integer :: lakeId, error
      real(r8) :: sir(NOUT)

      ! read lake information (i.e. depth, location ...)
      lakeId = lake_range(1)
      call ReadLakeName(lakeId)
      call LoadSensitiveParameters(sample)
      call ReadLakeInfo(lakeId)
      time = SimTime(Start_Year,Start_Month,Start_Day,End_Year, &
                     End_Month,End_Day)
      spinup = SimTime(Start_Year-nSpinup, Spinup_Month, Spinup_Day, &
                     Start_Year, Start_Month, Start_Day)
      call InitializeSimulation()
      call ModelRun(lakeId, time, spinup, error)
      if (error==0) then
         call GetSensitivityReturn(sir)  
      else
         sir = -9999.0_r8
      end if
      odata = sir
      print "(A, I0, A, I0, A, E14.6)", "Sample ", sampleId, &
            ": Error ", error
      call FinalizeSimulation()
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Get simulation returns for sensitivity test
   !     1. surface Twater (above 5 m) during 5/15-11/15
   !     2. bottom Twater (below 10 m) during 5/15-11/15
   !     3. ice-on DOY
   !     4. ice-off DOY
   !
   !------------------------------------------------------------------------------
   subroutine GetSensitivityReturn(odata)
      implicit none
      real(r8), intent(out) :: odata(NOUT)
      real(r8), allocatable :: tmp_zs(:)
      real(r8), allocatable :: tmp_zb(:)
      real(r8), allocatable :: tmp_zc(:)
      real(r8) :: avg_zs, avg_zb, avg_zc 
      real(r8) :: avg_co2, co2_hr(1), zco2(1)
      integer :: JDN0, JDN1, JDNb, JDNe
      integer :: ii, jj, izs, izb, izc
      integer :: idx0, idx1, nt1, nt2
      integer :: iceon, iceoff

      odata = 0.0_r8
      ! get the z indices
      izs = count(m_Zw<=5)
      izb = count(m_Zw>=10)
      izc = count(m_Zw<=8.42)
      allocate(tmp_zs(izs))
      allocate(tmp_zb(izb))
      allocate(tmp_zc(izc))
      izb = WATER_LAYER + 2 - izb
      ! get the mean annual values
      call Date2JDN(Start_Year, Start_Month, Start_Day, JDN0)
      call Date2JDN(End_Year, End_Month, End_Day, JDN1)
      nt1 = 0
      nt2 = 0
      do ii = Start_Year, End_Year, 1
         call Date2JDN(ii, 5, 15, JDNb)          
         call Date2JDN(ii, 11, 15, JDNe)
         if (JDNb>=JDN0 .and. JDNe<=JDN1) then
            nt1 = nt1 + 1
            idx0 = 24 * (JDNb - JDN0) + 1
            idx1 = 24 * (JDNe - JDN0)
            ! temperature
            call Mean(DBLE(m_tempwHist(1:izs,idx0:idx1)), 2, tmp_zs)
            call Mean(DBLE(m_tempwHist(izb:WATER_LAYER+1,idx0:idx1)), &
               2, tmp_zb)
            call WeightMean(tmp_zs, m_dZw(1:izs), avg_zs)
            call WeightMean(tmp_zb, m_dZw(izb:WATER_LAYER+1), avg_zb)
            odata(1) = odata(1) + avg_zs 
            odata(2) = odata(2) + avg_zb
         end if
         call Date2JDN(ii, 1, 1, JDNb)
         call Date2JDN(ii, 12, 31, JDNe)
         if (JDNb>=JDN0 .and. JDNe<=JDN1) then
            nt2 = nt2 + 1
            idx0 = 24 * (JDNb - JDN0) + 1
            idx1 = 24 * (JDNe - JDN0) + 24
            ! ice-on and ice-off DOY 
            iceoff = 1 
            do jj = idx0, idx1, 1
               if (m_iceHist(jj)<1d-6) then
                  iceoff = INT((jj-idx0)/24.0) + 1
                  exit 
               end if
            end do
            iceon = JDNe - JDNb + 1 
            do jj = idx1, idx0, -1
               if (m_iceHist(jj)<1d-6) then
                  iceon = INT((jj-idx0)/24.0) + 1 
                  exit
               end if
            end do
            odata(3) = odata(3) + DBLE(iceon)
            odata(4) = odata(4) + DBLE(iceoff)
         end if
      end do
      if (nt1>0) then
         odata(1:2) = odata(1:2) / DBLE(nt1)
      end if
      if (nt2>0) then
         odata(3:4) = odata(3:4) / DBLE(nt2)
      end if
      deallocate(tmp_zs)
      deallocate(tmp_zb)
      deallocate(tmp_zc)
   end subroutine

end module sensitivity_mod
