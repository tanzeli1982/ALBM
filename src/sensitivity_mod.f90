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
   use ifport 
   use mpi

   implicit none
   private
   public :: RunSensitivity 
   ! current sample Id 
   integer :: cur_sample

contains
   subroutine RunSensitivity(taskid, numprocs, arg)
      implicit none
      integer, intent(in) :: taskid
      integer, intent(in) :: numprocs
      character(len=*), intent(in) :: arg
      real(r8), allocatable :: samples(:,:)     ! generated samples
      integer, allocatable :: sampleIds(:)
      integer, parameter :: ndid = 400
      integer :: ii, jj, cnt, minid, maxid, err
      integer :: sampleId, idx, nsample, itmp
      integer :: sample_next_range(2)
      character(cx) :: script

      allocate(sampleIds(numprocs))
      minid = minval(sample_range)
      maxid = maxval(sample_range)
      nsample = maxid - minid + 1

      if (masterproc) then
         print "(A, I0, A, I0)", 'Run samples from ', minid, ' to ', maxid
      end if

      call DoSensitivityWarmup(minid)

      ! check the end sample id, if not end, refresh the sample range
      ! and restart with a new job
      !if (masterproc .and. (.NOT. DEBUG)) then
      !   if (maxid<NMAXSAMPLE) then
      !      sample_next_range = sample_range
      !      sample_range = (/minid+ndid, min(maxid+ndid,NMAXSAMPLE)/)
      !      call WriteSimulationSettings(arg)
      !      call GetFullFileName('bLakeJob.sub', script)
      !      script = "sbatch " // trim(script)
      !      err = system(trim(script))
      !      print "(A, I0)", "A new job is submitted. Return = ", err
      !      sample_range = sample_next_range
      !   end if
      !end if

      allocate(samples(NMAXSAMPLE,NPARAM))
      if (masterproc) then
         call ReadParameterSamples(NMAXSAMPLE, samples)
      end if
      do ii = 1, NPARAM, 1
         call MPI_BCAST(samples(:,ii), NMAXSAMPLE, MPI_REAL8, 0, &
                        MPI_COMM_WORLD, err)
      end do
      
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
         call MPI_BCAST(sampleIds, numprocs, MPI_INTEGER, 0, MPI_COMM_WORLD, err)
         sampleId = sampleIds(taskid+1)
         cur_sample = sampleId
         call MonteCarloSimulation( sampleId, samples(sampleId,:), err )
         print "(A, I0, A, I0, A, I0)", "Sample ", sampleId, ", processor ", &
               taskid, ", error ", err
      end do

      deallocate(samples)
      deallocate(sampleIds)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: For slave ndoes, run sensitivity simulations for each sample
   !
   !------------------------------------------------------------------------------
   subroutine MonteCarloSimulation(sampleId, sample, error)
      implicit none
      integer, intent(in) :: sampleId
      real(r8), intent(in) :: sample(NPARAM)
      integer, intent(out) :: error
      type(SimTime) :: time, spinup, otime
      integer :: i4ret, lakeId

      i4ret = SIGNALQQ(SIG$FPE, hand_fpe)
      ! read lake information (i.e. depth, location ...)
      lakeId = lake_range(1)
      call ReadLakeInfo(lakeId, sampleId)
      call LoadSensitiveParameters(sample)
      time = SimTime(Start_Year,Start_Month,Start_Day,End_Year, &
                     End_Month,End_Day)
      spinup = SimTime(Start_Year-nSpinup, Spinup_Month, Spinup_Day, &
                       Start_Year, Start_Month, Start_Day)
      otime = SimTime(SA_Start_Year, SA_Start_Month, SA_Start_Day, &
                      SA_End_Year, SA_End_Month, SA_End_Day)  
      call InitializeSimulation()
      call ModelRun(sampleId, time, spinup, otime, error)
      call ArchiveSensitivityOutput(sampleId, otime)
      call FinalizeSimulation()
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: create sensitivity run archive files 
   !
   !------------------------------------------------------------------------------
   subroutine DoSensitivityWarmup(minid)
      implicit none
      integer, intent(in) :: minid
      type(SimTime) :: time 

      time = SimTime(SA_Start_Year, SA_Start_Month, SA_Start_Day, &
                     SA_End_Year, SA_End_Month, SA_End_Day)
      if (minid==1 .and. masterproc) then
         call CreateOutputFile(time, NWLAYER+1, 'zw', 'water layer depth', 'm')
         call CreateOutputFile(time, NSLAYER+1, 'zs', 'sediment layer depth', 'm')
         call CreateOutputFile(time, NWLAYER+1, 'Az', 'water layer ' // &
                              'cross-section area', 'm^2');
         call CreateOutputFile(time, NWLAYER+1, 'colindx', 'water layer ' // &
                              'connected sediment column index', 'index')
         if (Thermal_Module) then
            call CreateOutputFile(time, 'snowthick', 'snow cover thickness', &
                                 'm', -9999.0_r4)
            call CreateOutputFile(time, 'icethick', 'ice cover thickness', &
                                 'm', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'watertemp', 'Z', 'water ' // &
                                 'temperature', 'K', -9999.0_r4)
         end if
         if (Carbon_Module) then
            call CreateOutputFile(time, 'ch4df', 'surface methane diffusion flux', &
                                 'mol/m2/s', -9999.0_r4)
            call CreateOutputFile(time, 'gpp', 'total gross primary production', &
                                 'gC/m2/s', -9999.0_r4)
            call CreateOutputFile(time, 'npp', 'total net primary production', &
                                 'gC/m2/s', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'ch4oxid', 'Z', 'methane oxidation', &
                                 'mol/m3/s', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'och4prod', 'Z', 'oxic methane ' // &
                                 'production', 'mol/m3/s', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'ch4conc', 'Z', 'dissolved ' // &
                                 'CH4 concentration', 'mol/m3', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'o2conc', 'Z', 'dissolved ' // &
                                 'O2 concentration', 'mol/m3', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'chla', 'Z', 'chlorophyll ' // &
                                 'concentration', 'g/m3', -9999.0_r4)
            call CreateOutputFile(time, NSCOL, 'cdep', 'COL', 'active organic ' // &
                                 'carbon deposition', 'gC/m2/s', -9999.0_r4)
            call CreateOutputFile(time, NPOC, NWLAYER+1, 'biomass', 'POC', 'Z', &
                                 'phytoplankton biomass', 'gC/m3', -9999.0_r4)
            call CreateOutputFile(time, NWLAYER+1, 'bveg', 'Z', 'submerged ' // &
                                 'macrophyte biomass', 'gC m-2', -9999.0_r4)
         end if
         if (Diagenesis_Module) then
            call CreateOutputFile(time, NSCOL, 'sedch4df', 'COL', 'sediment methane ' // &
                                 'diffusion flux', 'mol/m2/s', -9999.0_r4)
            call CreateOutputFile(time, NSCOL, 'sedch4eb', 'COL', 'sediment methane ' // &
                                 'ebullition flux', 'mol/m2/s', -9999.0_r4)
            call CreateOutputFile(time, NSCOL, 'ch4prod', 'COL', 'total sediment  ' // &
                                 'methane production', 'mol/m2/s', -9999.0_r4)
            call CreateOutputFile(time, NSCOL, NSLAYER+1, 'actcarb', 'COL', 'Z', &
                                 'sediment available active carbon density', 'gC/m3', &
                                 -9999.0_r4)
            call CreateOutputFile(time, NSCOL, NSLAYER+1, 'oldcarb', 'COL', 'Z', &
                                 'sediment available permafrost carbon density', &
                                 'gC/m3', -9999.0_r4)
         end if
         if (Bubble_Module) then
            call CreateOutputFile(time, 'ch4eb', 'surface methane ebullition flux', &
                                 'mol/m2/s', -9999.0_r4)
            call CreateOutputFile(time, 'icebch4', 'total CH4 in ice-trapped ' // &
                                 'bubbles', 'mole', -9999.0_r4)
         end if
      end if
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

end module sensitivity_mod
