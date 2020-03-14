module bayesian_mod
!---------------------------------------------------------------------------------
! Purpose: this module governs the Monte Carlo analysis.
!
!---------------------------------------------------------------------------------
   use shr_ctrl_mod
   use shr_typedef_mod
   use shr_param_mod
   use sim_coupler_mod
   use costfunc_mod
   use read_data_mod
   use io_utilities_mod

   implicit none
   private
   public :: RunMonteCarlo
   ! current sample Id
   integer :: cur_sample

contains
   subroutine RunMonteCarlo(arg)
      implicit none
      character(len=*), intent(in) :: arg
      real(r8), allocatable :: samples(:,:)     ! generated samples
      real(r8), allocatable :: results(:,:)
      real(r8), allocatable :: sir(:)
      real(r8) :: weights(20)
      character(len=8) :: varnames(20)
      integer :: ii, cnt, minid, maxid, err
      integer :: sampleId, idx, nsample, itmp
      integer :: nvar, sample_next_range(2)

      minid = minval(sample_range)
      maxid = maxval(sample_range)
      nsample = maxid - minid + 1

      print "(A, I0, A, I0)", 'Run samples from ', minid, ' to ', maxid

      allocate(samples(NMAXSAMPLE,NPARAM))

      call ReadParameterSamples(NMAXSAMPLE, samples)
      call ReadCalibVariables(varnames, weights, nvar)

      allocate(results(nvar,nsample))
      allocate(sir(nvar))
      
      if (minid==1) then
         call CreateSampleResultFile(nvar, 9999.0_r8)
      end if

      do cnt = 1, nsample, 1
         sampleId = minid + cnt - 1
         cur_sample = sampleId
         call MonteCarloSimulation( sampleId, samples(sampleId,:), &
            varnames(1:nvar), weights(1:nvar), sir )
         results(:,cnt) = sir
      end do

      call WriteSampleResults(sample_range, results)

      deallocate(samples)
      deallocate(results)
      deallocate(sir)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: For slave ndoes, run sensitivity simulations for each sample
   !
   !------------------------------------------------------------------------------
   subroutine MonteCarloSimulation(sampleId, sample, vars, weights, &
                                   odata)
      implicit none
      integer, intent(in) :: sampleId
      real(r8), intent(in) :: sample(NPARAM)
      character(len=8), intent(in) :: vars(:)
      real(r8), intent(in) :: weights(:)
      real(r8), intent(out) :: odata(:)
      type(SimTime) :: time, spinup
      integer :: lakeId, error
      real(r8) :: sir

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
         call EvaluateSampleLikelihood(time, vars, odata)
      else
         odata = 9999.0_r8
      end if
      sir = sum(weights*odata) / sum(weights)
      print "(A, I0, A, I0, A, E14.6)", "Sample ", sampleId, &
            ": Error ", error, ", Likelihood ", sir
      call FinalizeSimulation()
   end subroutine

end module bayesian_mod
