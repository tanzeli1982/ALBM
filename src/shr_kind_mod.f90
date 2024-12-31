module shr_kind_mod
!----------------------------------------------------------------------------
! precision/kind constants
!----------------------------------------------------------------------------
   use iso_fortran_env, only: r4 => REAL32, r8 => REAL64, i1 => INT8, &
                              i2 => INT16, i4 => INT32, i8 => int64
   public
   integer, parameter :: SHR_KIND_CS = 80          ! short char
   integer, parameter :: SHR_KIND_CL = 256         ! long char
   integer, parameter :: SHR_KIND_CX = 512         ! extra-long char
#ifdef USE_INTEL_COMPILER
   real(r8), parameter :: NaN = O'0777610000000000000000'
#else
   real(r8), parameter :: NaN = transfer(Z'7FF8000000000000', 0.0_r8)
#endif
end module shr_kind_mod
