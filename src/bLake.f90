!****************************************************************************
!
!  PROGRAM: Advanced Lake Biogeochemistry Model (ALBM)
!
!  PURPOSE:  Entry point for the ALBM model application.
!            It can be run in two modes: regular and sensitivity simulation.
!
!****************************************************************************

program bLake
   use read_data_mod,   only : ReadSimulationSettings, BcastSimulationSettings 
   use simulation_mod,  only : RunSimulation
   use sensitivity_mod, only : RunSensitivity
   use shr_ctrl_mod
   use mpi

   implicit none
   integer :: err, numprocs, taskid
   character(len=32) :: arg

   call MPI_INIT(err)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,err)
   call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,err)

   if (taskid==0) then
      masterproc = .True.
   else
      masterproc = .False.
   end if

   call get_command_argument(1, arg)
   if (len_trim(arg)==0) then
      arg = 'namelist.bLake'
   end if

   if (masterproc) then
      call ReadSimulationSettings(arg)
   end if
   call BcastSimulationSettings()

   if (trim(run_mode)=='sensitivity') then
      call RunSensitivity(taskid, numprocs, arg)
   else if (trim(run_mode)=='regular') then
      call RunSimulation(taskid, numprocs, arg)
   else
      print *, "Error: unrecognized run_mode!! Must be " // &
         "'regular' or 'sensitivity'."
      call MPI_ABORT(MPI_COMM_WORLD, 1, err)
   end if

   call MPI_FINALIZE(err)
   stop
end program bLake

