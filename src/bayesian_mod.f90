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
   use ifport 
   use mpi

   implicit none
   private
   public :: RunMonteCarlo
   ! current sample Id
   integer :: cur_sample

contains
   subroutine RunMonteCarlo(taskid, numprocs, arg)
      implicit none
      integer, intent(in) :: taskid
      integer, intent(in) :: numprocs
      character(len=*), intent(in) :: arg
      real(r8), allocatable :: samples(:,:)     ! generated samples
      real(r8), allocatable :: results(:,:)
      real(r8), allocatable :: sirs(:)
      real(r8), allocatable :: sir(:)
      real(r8) :: weights(20)
      character(len=8) :: varnames(20)
      integer, allocatable :: sampleIds(:)
      integer, parameter :: ndid = 500
      integer :: ii, cnt, minid, maxid, err
      integer :: sampleId, idx, nsample, itmp
      integer :: nvar, sample_next_range(2)
      character(cx) :: script

      allocate(sampleIds(numprocs))
      minid = minval(sample_range)
      maxid = maxval(sample_range)
      nsample = maxid - minid + 1

      if (masterproc) then
         print "(A, I0, A, I0)", 'Run samples from ', minid, ' to ', maxid
      end if

      ! check the end sample id, if not end, refresh the sample range
      ! and restart with a new job
      if (masterproc .and. (.NOT. DEBUG)) then
         if (maxid<NMAXSAMPLE) then
            sample_next_range = sample_range
            sample_range = (/minid+ndid, min(maxid+ndid,NMAXSAMPLE)/)
            call WriteSimulationSettings(arg)
            call GetFullFileName('bLakeJob.sub', script)
            script = "qsub " // trim(script)
            err = system(trim(script))
            print "(A, I0)", "A new job is submitted. Return = ", err
            sample_range = sample_next_range
         end if
      end if

      allocate(samples(NMAXSAMPLE,NPARAM))

      if (masterproc) then
         call ReadParameterSamples(NMAXSAMPLE, samples)
         call ReadCalibVariables(varnames, weights, nvar)
      end if
      call MPI_BCAST(nvar, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, err) 
      call MPI_BCAST(weights, size(weights), MPI_REAL8, 0, &
                     MPI_COMM_WORLD, err)
      call MPI_BCAST(varnames, size(varnames)*len(varnames), &
                     MPI_CHARACTER, 0, MPI_COMM_WORLD, err)
      do ii = 1, NPARAM, 1
         call MPI_BCAST(samples(:,ii), NMAXSAMPLE, MPI_REAL8, 0, &
                        MPI_COMM_WORLD, err)
      end do

      if (masterproc) then
         allocate(results(nvar,nsample))
         allocate(sirs(nvar*numprocs))
      end if
      allocate(sir(nvar))
      
      if (minid==1) then
         call CreateSampleResultFile(nvar, 9999.0_r8)
      end if

      do cnt = 1, nsample, numprocs
         if (masterproc) then
            if (nsample-cnt+1>=numprocs) then
               sampleIds = (/(ii, ii = minid+cnt-1, minid+cnt+numprocs-2)/)
            else
               itmp = nsample - cnt + 1
               sampleIds(1:itmp) = (/(ii, ii = minid+cnt-1, maxid)/)
               sampleIds(itmp+1:numprocs) = (/(ii, ii = minid, &
                  minid+numprocs-itmp-1)/)
            end if
         end if
         call MPI_BCAST(sampleIds, numprocs, MPI_INTEGER, 0, &
                        MPI_COMM_WORLD, err)
         sampleId = sampleIds(taskid+1)
         cur_sample = sampleId
         call MonteCarloSimulation( sampleId, samples(sampleId,:), &
            varnames(1:nvar), weights(1:nvar), sir )
         call MPI_GATHER(sir, nvar, MPI_REAL8, sirs, nvar, MPI_REAL8, &
                         0, MPI_COMM_WORLD, err)
         if (masterproc) then
            do ii = 1, numprocs, 1
               idx = sampleIds(ii) - minid + 1 
               results(:,idx) = sirs(nvar*(ii-1)+1:nvar*ii)
            end do
         end if
      end do

      if (masterproc) then
         call WriteSampleResults(sample_range, results)
      end if

      deallocate(samples)
      deallocate(sampleIds)
      if (masterproc) then
         deallocate(results)
         deallocate(sirs)
      end if
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
      integer :: i4ret, lakeId, error
      real(r8) :: sir

      i4ret = SIGNALQQ(SIG$FPE, hand_fpe)
      ! read lake information (i.e. depth, location ...)
      lakeId = lake_range(1)
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

   !------------------------------------------------------------------------------
   !
   ! Purpose: some utilities for exceptions: SIG$FPE, SIG$ABORT, SIG$SEGV
   !
   !------------------------------------------------------------------------------
   function hand_fpe(sigid, except)
      !DEC$ ATTRIBUTES C :: hand_fpe
      use ifport
      !use ifcore
      INTEGER(4) :: hand_fpe
      INTEGER(2) :: sigid, except

      if (sigid/=SIG$FPE) then
         hand_fpe = 1
         return
      end if
      select case(except)
         case( FPE$INVALID )
            print *, ' Floating point exception: Invalid number'
         case( FPE$DENORMAL )
            print *, ' Floating point exception: Denormalized number'
         case( FPE$ZERODIVIDE )
            print *, ' Floating point exception: Zero divide'
         case( FPE$OVERFLOW )
            print *, ' Floating point exception: Overflow'
         case( FPE$UNDERFLOW )
            print *, ' Floating point exception: Underflow'
         case( FPE$INEXACT )
            print *, ' Floating point exception: Inexact precision'
         case default
            print *, ' Floating point exception: Non-IEEE type'
      end select
      !CALL TRACEBACKQQ(trim(header), USER_EXIT_CODE=-1)
      print *, 'failed sample ', cur_sample, sa_params 
      hand_fpe = 1
   end function

end module bayesian_mod
