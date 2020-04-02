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
   use shr_ctrl_mod

   implicit none
   character(len=32) :: arg

   call get_command_argument(1, arg)
   if (len_trim(arg)==0) then
      arg = 'namelist.bLake'
   end if

   call ReadSimulationSettings(arg)
   call RunRegular(arg)

   stop
end program bLake

