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
   use read_data_mod,   only : ReadSimulationSettings
   use simulation_mod,  only : RunRegular
   use bayesian_mod,    only : RunMonteCarlo
   use sensitivity_mod, only : RunSensitivity
   use shr_ctrl_mod

   implicit none
   character(len=32) :: arg

   call get_command_argument(1, arg)
   if (len_trim(arg)==0) then
      arg = 'namelist.bLake'
   end if

   call ReadSimulationSettings(arg)

   if (trim(run_mode)=='bayesian') then
      call RunMonteCarlo(arg) 
   else if (trim(run_mode)=='sensitivity') then
      call RunSensitivity(arg)
   else
      call RunRegular(arg)
   end if

   stop
end program bLake

