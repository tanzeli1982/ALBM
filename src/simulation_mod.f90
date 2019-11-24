module simulation_mod
!---------------------------------------------------------------------------------
! Purpose: This module runs simulations locally or regionally.
!  There are two main components in this model: thermal module and bubbling module.
!  Thermal module includes water thermal sub module and soil thermal sub module.
!  Bubbling module includes gas production module and gas transportation module.
!
!---------------------------------------------------------------------------------
   use shr_ctrl_mod
   use shr_param_mod
   use shr_typedef_mod
   use sim_coupler_mod
   use read_data_mod
   use ifport
   use mpi

   private
   public :: RunRegular

contains
   subroutine RunRegular(taskid, numprocs, arg)
      implicit none
      integer, intent(in) :: taskid
      integer, intent(in) :: numprocs
      character(len=*), intent(in) :: arg
      integer, allocatable :: lakeIds(:)
      integer, parameter :: ndid = 200
      integer :: lake_next_range(2)
      integer :: lakeId, err, minid, maxid
      integer :: nlake, istep, ii, itmp
      integer :: ntotlake, jj
      character(cx) :: script

      allocate(lakeIds(numprocs))
      minid = minval(lake_range)
      maxid = maxval(lake_range)
      nlake = maxid - minid + 1

      if (masterproc) then
         print "(A, I0, A, I0)", 'Run lakes from ', minid, ' to ', maxid
      end if

      call DoSimulationWarmup(minid, ntotlake)

      ! check the end sample id, if not end, refresh the sample range
      ! and restart with a new job
      if (masterproc .and. RESUBMIT) then
         if (maxid<ntotlake) then
            lake_next_range = lake_range
            lake_range = (/minid+ndid, min(maxid+ndid,ntotlake)/)
            call WriteSimulationSettings(arg)
            call GetFullFileName('bLakeJob.sub', script)
            script = "qsub " // trim(script)
            err = system(trim(script))
            print "(A, I0)", "A new job is submitted. Return = ", err
            lake_range = lake_next_range
         end if
      end if
      
      do istep = 1, nlake, numprocs
         if (masterproc) then
            if (nlake-istep+1>=numprocs) then
               lakeIds = (/(ii, ii = minid+istep-1, minid+istep+numprocs-2)/)
            else
               itmp = nlake - istep + 1
               lakeIds(1:itmp) = (/(ii, ii = minid+istep-1, maxid)/)
               do jj = itmp+1, numprocs, 1
                  lakeIds(jj) = minid + mod(jj-itmp-1, itmp)
               end do
            end if
         end if
         call MPI_BCAST(lakeIds, size(lakeIds), MPI_INTEGER, 0, &
                        MPI_COMM_WORLD, err)
         lakeId = lakeIds(taskid+1)
         call RegularSimulation(lakeId, err)
         print "(A, I0, A, I0, A, I0)", "Lake ", lakeId, ", processor ", &
               taskid, ", error ", err
      end do

      deallocate(lakeIds)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: create simulation archive files 
   !
   !------------------------------------------------------------------------------
   subroutine DoSimulationWarmup(minid, nlake)
      implicit none
      integer, intent(in) :: minid
      integer, intent(out) :: nlake
      type(SimTime) :: time
      
      call ReadLakeTotalNumber(nlake)
      if (minid==1 .or. DEBUG) then
         time = SimTime(Start_Year, Start_Month, Start_Day, End_Year, &
                     End_Month, End_Day)
         call CreateOutputFile(NWLAYER+1, 'zw', 'water layer depth', 'm') 
         call CreateOutputFile(time, 'snowthick', 'snow cover thickness', &
                              'm', -9999.0_r4)
         call CreateOutputFile(time, 'icethick', 'ice cover thickness', &
                              'm', -9999.0_r4)
         call CreateOutputFile(time, 'sensheatf', 'upward sensible ' // &
                              'heat flux', 'W m-2', -9999.0_r4)
         call CreateOutputFile(time, 'latentheatf', 'upward latent ' // &
                              'heat flux', 'W m-2', -9999.0_r4)
         call CreateOutputFile(time, 'momf', 'momentum energy flux', &
                              'kg m-1 s-2', -9999.0_r4)
         call CreateOutputFile(time, 'lwup', 'upward longwave radiation', &
                              'W m-2', -9999.0_r4)
         call CreateOutputFile(time, 'lakeheatf', 'downward net heat flux', &
                              'W m-2', -9999.0_r4)
         call CreateOutputFile(time, 'swdw', 'downward shortwave radiation', &
                              'W m-2', -9999.0_r4)
         call CreateOutputFile(time, 'swup', 'upward shortwave radiation', &
                              'W m-2', -9999.0_r4)
         call CreateOutputFile(time, 'sedheatf', 'upward sediment heat ' // &
                              'flux', 'W m-2', -9999.0_r4)
         call CreateOutputFile(time, NWLAYER+1, 'watertemp', 'water ' // &
                              'temperature', 'K', -9999.0_r4)
         call CreateOutputFile(time, NWLAYER+1, 'turbdiffheat', &
                               'Turbulent diffusivity of heat', 'm2 s-1', &
                               -9999.0_r4)
      end if
   end subroutine

   subroutine RegularSimulation(lakeId, error)
      implicit none
      integer, intent(in) :: lakeId
      integer, intent(out) :: error
      type(SimTime) :: time, spinup
      real(r8) :: OptParams(NPARAM)
      integer :: i4ret

      i4ret = SIGNALQQ(SIG$FPE, hand_fpe)
      ! read lake information (i.e. depth, location ...)
      call ReadLakeName(lakeId)
      call ReadOptimumParameters(OptParams)
      call LoadSensitiveParameters(OptParams)
      call ReadLakeInfo(lakeId)
      time = SimTime(Start_Year, Start_Month, Start_Day, End_Year, &
                     End_Month, End_Day)
      spinup = SimTime(Start_Year-nSpinup, Spinup_Month, Spinup_Day, &
                     Start_Year, Start_Month, Start_Day)
      call InitializeSimulation()
      call ModelRun(lakeId, time, spinup, error)
      call ArchiveModelOutput(lakeId, time)
      call FinalizeSimulation()
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: some utilities for exceptions
   !
   !------------------------------------------------------------------------------
   function hand_fpe(sigid, except)
      !DEC$ ATTRIBUTES C :: hand_fpe
      use ifport
      INTEGER(4) :: hand_fpe
      INTEGER(2) :: sigid, except

      if (sigid/=SIG$FPE) then
         print "('The hand_fpe is not for signal ', I0)", sigid
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
      print *, 'lake failed: ', lake_info 
      hand_fpe = 1
   end function

end module
