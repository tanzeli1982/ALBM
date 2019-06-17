!****************************************************************************
!
!  PROGRAM: Generalized Lake Biogeochemical Model for temperate lakes
!
!  PURPOSE:  Entry point for the GLBM model application.
!            It can be run in three modes: regular simulation, model calibration,
!            and model sensitivity analysis.
!
!****************************************************************************

program bLake
   use read_data_mod,   only : ReadSimulationSettings, BcastSimulationSettings 
   use simulation_mod,  only : RunRegular
   use bayesian_mod,    only : RunMonteCarlo
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

   if (trim(run_mode)=='bayesian') then
      call RunMonteCarlo(taskid, numprocs, arg) 
   else if (trim(run_mode)=='sensitivity') then
      call RunSensitivity(taskid, numprocs, arg)
   else
      call RunRegular(taskid, numprocs, arg)
   end if

   call MPI_FINALIZE(err)
   stop
end program bLake

