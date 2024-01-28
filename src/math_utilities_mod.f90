module math_utilities_mod
!---------------------------------------------------------------------------------
! Purpose:
!
! This module serves for the implementation of some basic physical utilities
!
!---------------------------------------------------------------------------------
   use shr_kind_mod,       only : r8, r4, i8, NaN
   use shr_ctrl_mod,       only : inft => INFINITESIMAL_E8, inf => INFINITE_E8, &
                                  TOL_E8
   use shr_typedef_mod,    only : RungeKuttaCache1D, RungeKuttaCache2D, &
                                  RungeKuttaCache3D
   use ifport

   implicit none
   integer, parameter :: adaptive_mode = 101, fixed_mode = 102
   integer, parameter :: MAXITER = 100

   interface RungeKutta4
      module procedure RungeKutta4_1D
      module procedure RungeKutta4_2D
      module procedure RungeKutta4_3D
   end interface

   interface Norm
      module procedure Norm1d
      module procedure Norm2d
      module procedure Norm3d
   end interface

   interface Divide
      module procedure Divide1d
      module procedure Divide2d
      module procedure Divide3d
   end interface

   interface CalcFiniteDifference
      module procedure CalcFiniteDifference1D
      module procedure CalcFiniteDifference2D
      module procedure CalcFiniteDifference3D
   end interface

   interface Mean
      module procedure Mean1d
      module procedure Mean2d
      module procedure Mean3d
      module procedure Mean1d_r4
      module procedure Mean2d_r4
      module procedure Mean3d_r4
   end interface

   interface WeightMean
      module procedure Mean1d_weight
      module procedure Mean2d_weight
      module procedure Mean3d_weight
   end interface

   interface GetLeastValue
      module procedure GetLeastValue1d
      module procedure GetLeastValue2d
      module procedure GetLeastValue3d
   end interface

contains
   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate the mean value of a vector
   !
   !------------------------------------------------------------------------------
   subroutine Mean1d(vec, mean)
      implicit none
      real(r8), intent(in) :: vec(:)
      real(r8), intent(out) :: mean

      mean = sum(vec) / DBLE(size(vec))
   end subroutine

   subroutine Mean2d(vec, dir, mean)
      implicit none
      real(r8), intent(in) :: vec(:,:)
      integer, intent(in) :: dir
      real(r8), intent(out) :: mean(:)

      ! dir indicates the dimension to be averaged 
      mean = sum(vec,dir) / DBLE(size(vec,dir))
   end subroutine

   subroutine Mean3d(vec, dir, mean)
      implicit none
      real(r8), intent(in) :: vec(:,:,:)
      integer, intent(in) :: dir
      real(r8), intent(out) :: mean(:,:)

      ! dir indicates the dimension to be averaged
      mean = sum(vec,dir) / DBLE(size(vec,dir))
   end subroutine

   subroutine Mean1d_r4(vec, mean)
      implicit none
      real(r4), intent(in) :: vec(:)
      real(r4), intent(out) :: mean

      mean = sum(vec) / REAL(size(vec),4)
   end subroutine

   subroutine Mean2d_r4(vec, dir, mean)
      implicit none
      real(r4), intent(in) :: vec(:,:)
      integer, intent(in) :: dir
      real(r4), intent(out) :: mean(:)

      ! dir indicates the dimension to be averaged 
      mean = sum(vec,dir) / REAL(size(vec,dir),4)
   end subroutine

   subroutine Mean3d_r4(vec, dir, mean)
      implicit none
      real(r4), intent(in) :: vec(:,:,:)
      integer, intent(in) :: dir
      real(r4), intent(out) :: mean(:,:)

      ! dir indicates the dimension to be averaged
      mean = sum(vec,dir) / REAL(size(vec,dir),4)
   end subroutine

   subroutine Mean1d_weight(vec, weight, mean)
      implicit none
      real(r8), intent(in) :: vec(:)
      real(r8), intent(in) :: weight(:)
      real(r8), intent(out) :: mean

      mean = sum(vec*weight) / sum(weight)
   end subroutine

   subroutine Mean2d_weight(vec, weight, dir, mean)
      implicit none
      real(r8), intent(in) :: vec(:,:)
      real(r8), intent(in) :: weight(:)
      integer, intent(in) :: dir
      real(r8), intent(out) :: mean(:)
      real(r8) :: wsum
      integer :: ii, nn

      ! dir indicates the dimension to be averaged 
      wsum = sum(weight)
      if (dir==1) then
         nn = size(vec,2)
         do ii = 1, nn, 1
            mean(ii) = sum(vec(:,ii)*weight) / wsum
         end do
      else if (dir==2) then
         nn = size(vec,1)
         do ii = 1, nn, 1
            mean(ii) = sum(vec(ii,:)*weight) / wsum
         end do
      end if
   end subroutine

   subroutine Mean3d_weight(vec, weight, dir, mean)
      implicit none
      real(r8), intent(in) :: vec(:,:,:)
      real(r8), intent(in) :: weight(:)
      integer, intent(in) :: dir
      real(r8), intent(out) :: mean(:,:)
      integer :: ii, jj, nn, mm
      real(r8) :: wsum

      ! dir indicates the dimension to be averaged.
      wsum = sum(weight)
      if (dir==1) then
         nn = size(vec,2)
         mm = size(vec,3)
         do ii = 1, nn, 1
            do jj = 1, mm, 1
               mean(ii,jj) = sum(vec(:,ii,jj)*weight) / wsum
            end do
         end do
      else if (dir==2) then
         nn = size(vec,1)
         mm = size(vec,3)
         do ii = 1, nn, 1
            do jj = 1, mm, 1
               mean(ii,jj) = sum(vec(ii,:,jj)*weight) / wsum
            end do
         end do
      else if (dir==3) then
         nn = size(vec,1)
         mm = size(vec,2)
         do ii = 1, nn, 1
            do jj = 1, mm, 1
               mean(ii,jj) = sum(vec(ii,jj,:)*weight) / wsum
            end do
         end do
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Calculate the infinite norm of a vector
   !
   !------------------------------------------------------------------------------
   subroutine Norm1d(vec, norm)
      implicit none
      real(r8), intent(in) :: vec(:)
      real(r8), intent(out) :: norm

      norm = max(abs(minval(vec)), abs(maxval(vec)))
   end subroutine

   subroutine Norm2d(matrix, dir, norm)
      implicit none
      real(r8), intent(in) :: matrix(:,:)
      integer, intent(in) :: dir
      real(r8), intent(out) :: norm(:)
      integer :: ii, nn

      ! dir indicates the dimension to be kept.
      nn = size(matrix,dir)
      if (dir==1) then
         do ii = 1, nn, 1
            norm(ii) = max( abs(minval(matrix(ii,:))), &
               abs(maxval(matrix(ii,:))) )
         end do
      else if (dir==2) then
         do ii = 1, nn, 1
            norm(ii) = max( abs(minval(matrix(:,ii))), &
               abs(maxval(matrix(:,ii))) )
         end do
      end if
   end subroutine

   subroutine Norm3d(matrix, norm)
      implicit none
      real(r8), intent(in) :: matrix(:,:,:)
      real(r8), intent(out) :: norm 
      real(r8) :: tmp1, tmp2

      tmp1 = minval(matrix)
      tmp2 = maxval(matrix)
      norm = max(abs(tmp1), abs(tmp2))
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Solve tridiagonal linear equations by Chase After Method
   !
   !------------------------------------------------------------------------------
   subroutine ChaseAfter(mA, vf, vx)
      implicit none
      real(r8), intent(in) :: mA(:,:)
      real(r8), intent(in) :: vf(:)
      real(r8), intent(out) :: vx(:)
      real(r8) :: v1(size(vf)), v2(size(vf))
      integer :: ii, nn

      nn = size(vf)
      do ii = 1, nn, 1
         if (ii == 1) then
            v1(ii) = mA(3,ii)/mA(2,ii)
            v2(ii) = vf(ii)/mA(2,ii)
         else
            v1(ii) = mA(3,ii)/(mA(2,ii)-mA(1,ii)*v1(ii-1))
            v2(ii) = (vf(ii)-mA(1,ii)*v2(ii-1))/(mA(2,ii)-mA(1,ii)*v1(ii-1))
         end if
      end do
      do ii = nn, 1, -1
         if (ii == nn) then
            vx(ii) = v2(ii)
         else
            vx(ii) = v2(ii) - v1(ii)*vx(ii+1)
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Solve tridiagonal linear equations by Jacobi Iterative Method
   !
   !------------------------------------------------------------------------------
   subroutine Jacobi(mA, vf, vx)
      implicit none
      real(r8), intent(in) :: mA(:,:)
      real(r8), intent(in) :: vf(:)
      real(r8), intent(out) :: vx(:)
      real(r8), parameter :: tol = 1.0d-4
      real(r8) :: vxl(size(vf))
      real(r8) :: dx, tmp
      integer :: ii, nn

      nn = size(vf)
      do ii = 1, nn, 1
         vx(ii) = 0.0
      end do
      dx = 1.0
      do while (dx>tol)
         do ii = 1, nn, 1
            if (ii == 1) then
               tmp = mA(3,ii)*vx(ii+1)
            else if (ii == nn) then
               tmp = mA(1,ii)*vx(ii-1)
            else
               tmp = mA(1,ii)*vx(ii-1)+mA(3,ii)*vx(ii+1)
            end if
            vxl(ii) = (vf(ii)-tmp)/mA(2,ii)
         end do
         call Norm(vxl-vx,dx)
         vx = vxl
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: one matrix divided by another matrix
   !
   !------------------------------------------------------------------------------
   subroutine Divide1d(matrix1, matrix2, matrix)
      implicit none
      real(r8), intent(in) :: matrix1(:)
      real(r8), intent(in) :: matrix2(:)
      real(r8), intent(out) :: matrix(:)

      matrix = matrix1 / (inft + matrix2)
   end subroutine

   subroutine Divide2d(matrix1, matrix2, matrix)
      implicit none
      real(r8), intent(in) :: matrix1(:,:)
      real(r8), intent(in) :: matrix2(:,:)
      real(r8), intent(out) :: matrix(:,:)

      matrix = matrix1 / (inft + matrix2)
   end subroutine

   subroutine Divide3d(matrix1, matrix2, matrix)
      implicit none
      real(r8), intent(in) :: matrix1(:,:,:)
      real(r8), intent(in) :: matrix2(:,:,:)
      real(r8), intent(out) :: matrix(:,:,:)

      matrix = matrix1 / (inft + matrix2)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Finite difference
   !
   !------------------------------------------------------------------------------
   subroutine CalcFiniteDifference1D(X0, X1, dt, dX)
      implicit none
      real(r8), intent(in) :: X0(:)
      real(r8), intent(in) :: X1(:)
      real(r8), intent(in) :: dt
      real(r8), intent(out) :: dX(:)

      dX = (X1 - X0) / dt
   end subroutine

   subroutine CalcFiniteDifference2D(X0, X1, dt, dX)
      implicit none
      real(r8), intent(in) :: X0(:,:)
      real(r8), intent(in) :: X1(:,:)
      real(r8), intent(in) :: dt
      real(r8), intent(out) :: dX(:,:)

      dX = (X1 - X0) / dt 
   end subroutine

   subroutine CalcFiniteDifference3D(X0, X1, dt, dX)
      implicit none
      real(r8), intent(in) :: X0(:,:,:)
      real(r8), intent(in) :: X1(:,:,:)
      real(r8), intent(in) :: dt
      real(r8), intent(out) :: dX(:,:,:)

      dX = (X1 - X0) / dt 
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Find the most negative values
   !
   !------------------------------------------------------------------------------
   subroutine GetLeastValue1d(array, val)
      implicit none
      real(r8), intent(in) :: array(:)
      real(r8), intent(out) :: val

      val = minval(array)
   end subroutine

   subroutine GetLeastValue2d(matrix, dir, values)
      implicit none
      real(r8), intent(in) :: matrix(:,:)
      integer, intent(in) :: dir
      real(r8), intent(out) :: values(:)
      integer :: ii, nn

      ! dir indicates the dimension to be kept.
      if (dir==1) then
         nn = size(matrix,dir)
         do ii = 1, nn, 1
            values(ii) = minval(matrix(ii,:))
         end do
      else if (dir==2) then
         nn = size(matrix,dir)
         do ii = 1, nn, 1
            values(ii) = minval(matrix(:,ii))
         end do
      end if
   end subroutine

   subroutine GetLeastValue3d(matrix, val)
      implicit none
      real(r8), intent(in) :: matrix(:,:,:)
      real(r8), intent(out) :: val

      val = minval(matrix)
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: 1-dimensional Adaptive 4th Runge-Kutta-Fehlberg Method,
   !          <Numerical Analysis> (Richard L. Burden), pp254.
   !
   !  1-dimensional Adaptive 3th Bogacki-Shampine Method (wikipedia)
   !
   !------------------------------------------------------------------------------
   subroutine RungeKutta4_1D(odeFunc, mem, mode, tol, curstep, nextstep, &
                             invars, outvars)
      implicit none
      external :: odeFunc                          ! ODEs 
      type(RungeKuttaCache1D) :: mem               ! memory caches
      integer, intent(in) :: mode                  ! adaptive or fixed
      real(r8), intent(in) :: tol                  ! absolute tolerance
      real(r8), intent(inout) :: curstep           ! current time step
      real(r8), intent(out) :: nextstep            ! recommended next time step
      real(r8), intent(in) :: invars(:)            ! last state variables
      real(r8), intent(out) :: outvars(:)          ! current state variables
      real(r8) :: step, dy, dyn, rdy, delta, rate
      logical :: isLargeErr, isConstraintViolated
      integer :: iter

      isLargeErr = .True.
      isConstraintViolated = .False.
      step = curstep
      iter = 1
      call odeFunc(invars, mem%K1)
      do while (isLargeErr .or. isConstraintViolated)
         if (iter>MAXITER) then
            print *, "Runge-Kutta 1D iteration number is more than 100!!"
            exit
         end if
         curstep = step
         mem%interim = invars + step*0.25*mem%K1
         call odeFunc(mem%interim, mem%K2)
         mem%interim = invars + step*(0.09375*mem%K1 + 0.28125*mem%K2)
         call odeFunc(mem%interim, mem%K3)
         mem%interim = invars + step*(0.87938*mem%K1 - 3.27720*mem%K2 + &
                     3.32089*mem%K3)
         call odeFunc(mem%interim, mem%K4)
         mem%interim = invars + step*(2.03241*mem%K1 - 8.0*mem%K2 + &
                     7.17349*mem%K3 - 0.20590*mem%K4)
         call odeFunc(mem%interim, mem%K5)
         mem%nxt4th = invars + step*(0.11574*mem%K1 + 0.54893*mem%K3 + &
                     0.53533*mem%K4 - 0.2*mem%K5)
         if (mode==fixed_mode) then       ! fixed step
            nextstep = step
            outvars = mem%nxt4th
            return                                                         
         end if
         mem%interim = invars + step*(-0.29630*mem%K1 + 2.0*mem%K2 - &
                     1.38168*mem%K3 + 0.45297*mem%K4 - 0.275*mem%K5)
         call odeFunc(mem%interim, mem%K6)
         mem%nxt5th = invars + step*(0.11852*mem%K1 + 0.51899*mem%K3 + &
                     0.50613*mem%K4 - 0.18*mem%K5 + 0.03636*mem%K6)
         call Divide(mem%nxt4th-mem%nxt5th, mem%nxt4th, mem%rerr)
         call Norm(mem%rerr, rdy)
         call Norm(mem%nxt4th-mem%nxt5th, dy)
         call GetLeastValue(mem%nxt4th, dyn)
         isLargeErr = (dy>tol .and. rdy>TOL_E8)
         isConstraintViolated = (dyn<-100*tol)
         if (isConstraintViolated) then
            step = 0.5*step
         else
            rate = max(tol/(inft+dy), TOL_E8/(inft+rdy))
            delta = 0.84*rate**(0.25)
            if (delta<=0.1) then
               step = 0.1*step
            else if (delta>=4.0) then
               step = 4.0*step
            else
               step = delta*step
            end if
         end if
         iter = iter + 1
      end do
      nextstep = step
      outvars = mem%nxt4th
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: 2-dimensional Adaptive 4th Runge-Kutta-Fehlberg Method
   !
   !------------------------------------------------------------------------------
   subroutine RungeKutta4_2D(odeFunc, mem, mode, tol, curstep, nextstep, &
                             invars, outvars)
      implicit none
      external :: odeFunc                       ! ODEs 
      type(RungeKuttaCache2D) :: mem            ! memory caches
      integer, intent(in) :: mode               ! adaptive or fixed
      real(r8), intent(in) :: tol(:)            ! absolute tolerances
      real(r8), intent(inout) :: curstep        ! current time step
      real(r8), intent(out) :: nextstep         ! next time step
      real(r8), intent(in) :: invars(:,:)       ! last state variables
      real(r8), intent(out) :: outvars(:,:)     ! current state variables
      real(r8), dimension(size(invars,1)) :: dy
      real(r8), dimension(size(invars,1)) :: rdy
      real(r8), dimension(size(invars,1)) :: dyn
      real(r8), dimension(size(invars,1)) :: rlt_tol
      real(r8), dimension(size(invars,1)) :: abs_rate
      real(r8), dimension(size(invars,1)) :: rlt_rate
      real(r8) :: step, rate, delta
      logical :: isLargeErr, isConstraintViolated
      integer :: iter

      isLargeErr = .True.
      isConstraintViolated = .False.
      step = curstep
      iter = 1
      rlt_tol = TOL_E8
      call odeFunc(invars, mem%K1)
      do while (isLargeErr .or. isConstraintViolated)
         if (iter>MAXITER) then
            print *, "Runge-Kutta 2D iteration number is more than 100!!"
            exit
         end if 
         curstep = step
         mem%interim = invars + step*0.25*mem%K1
         call odeFunc(mem%interim, mem%K2)
         mem%interim = invars + step*(0.09375*mem%K1 + 0.28125*mem%K2)
         call odeFunc(mem%interim, mem%K3)
         mem%interim = invars + step*(0.87938*mem%K1 - 3.27720*mem%K2 + &
                     3.32089*mem%K3)
         call odeFunc(mem%interim, mem%K4)
         mem%interim = invars + step*(2.03241*mem%K1 - 8.0*mem%K2 + &
                     7.17349*mem%K3 - 0.20590*mem%K4)
         call odeFunc(mem%interim, mem%K5)
         mem%nxt4th = invars + step*(0.11574*mem%K1 + 0.54893*mem%K3 + &
                     0.53533*mem%K4 - 0.2*mem%K5)
         if (mode==fixed_mode) then       ! fixed step
            nextstep = step
            outvars = mem%nxt4th
            return                                                         
         end if
         mem%interim = invars + step*(-0.29630*mem%K1 + 2.0*mem%K2 - &
                     1.38168*mem%K3 + 0.45297*mem%K4 - 0.275*mem%K5)
         call odeFunc(mem%interim, mem%K6)
         mem%nxt5th = invars + step*(0.11852*mem%K1 + 0.51899*mem%K3 + &
                     0.50613*mem%K4 - 0.18*mem%K5 + 0.03636*mem%K6)
         call Divide(mem%nxt4th-mem%nxt5th, mem%nxt4th, mem%rerr)
         call Norm(mem%rerr, 1, rdy)
         call Norm(mem%nxt4th-mem%nxt5th, 1, dy)
         call GetLeastValue(mem%nxt4th, 1, dyn)
         isLargeErr = IsErrorStillLarge(dy,tol,rdy,rlt_tol)
         isConstraintViolated = IsTooLargeNegative(dyn,-100*tol)
         if (isConstraintViolated) then
            step = 0.5*step
         else
            call Divide(tol,dy,abs_rate)
            call Divide(rlt_tol,rdy,rlt_rate)
            rate = max(minval(abs_rate), minval(rlt_rate))
            delta = 0.84*rate**0.25
            if (delta<=0.1) then
               step = 0.1*step
            else if (delta>=4.0) then
               step = 4.0*step
            else
               step = delta*step
            end if
         end if
         iter = iter + 1
      end do
      nextstep = step
      outvars = mem%nxt4th
   end subroutine

   function IsErrorStillLarge(abs_diff, abs_tol, rlt_diff, rlt_tol)
      implicit none
      real(r8), intent(in) :: abs_diff(:)
      real(r8), intent(in) :: abs_tol(:)
      real(r8), intent(in) :: rlt_diff(:)
      real(r8), intent(in) :: rlt_tol(:)
      logical :: IsErrorStillLarge
      integer :: ii, nn

      IsErrorStillLarge = .False.
      nn = size(abs_diff)
      do ii = 1, nn, 1
         if (abs_diff(ii)>abs_tol(ii) .and. &
               rlt_diff(ii)>rlt_tol(ii)) then
            IsErrorStillLarge = .True.
            return
         end if
      end do
      return
   end function

   function IsTooLargeNegative(abs_diff, abs_tol)
      implicit none
      real(r8), intent(in) :: abs_diff(:)
      real(r8), intent(in) :: abs_tol(:)
      logical :: IsTooLargeNegative
      integer :: ii, nn

      IsTooLargeNegative = .False.
      nn = size(abs_diff)
      do ii = 1, nn, 1
         if (abs_diff(ii)<abs_tol(ii)) then
            IsTooLargeNegative = .True.
            return
         end if
      end do
      return
   end function

   !------------------------------------------------------------------------------
   !
   ! Purpose: 3-dimensional Adaptive 4th Runge-Kutta-Fehlberg Method
   !
   !------------------------------------------------------------------------------
   subroutine RungeKutta4_3D(odeFunc, mem, mode, tol, curstep, nextstep, &
                             invars, outvars)
      external :: odeFunc                       ! ODEs 
      type(RungeKuttaCache3D) :: mem            ! memory caches
      integer, intent(in) :: mode               ! adaptive or fixed
      real(r8), intent(in) :: tol(:)            ! absolute tolerances
      real(r8), intent(inout) :: curstep        ! current time step
      real(r8), intent(out) :: nextstep         ! next time step
      real(r8), intent(in) :: invars(:,:,:)     ! last state variables
      real(r8), intent(out) :: outvars(:,:,:)   ! current state variables
      real(r8), dimension(size(invars,1)) :: dy, dytmp
      real(r8), dimension(size(invars,1)) :: rdy, rdytmp
      real(r8), dimension(size(invars,1)) :: dyn, dyntmp
      real(r8), dimension(size(invars,1)) :: rlt_tol
      real(r8), dimension(size(invars,1)) :: abs_rate
      real(r8), dimension(size(invars,1)) :: rlt_rate
      real(r8) :: step, rate, delta
      logical :: isLargeErr, isConstraintViolated
      integer :: iter, kk, ii, ndim1, ndim2

      ndim1 = size(invars,1)
      ndim2 = size(invars,2)
      isLargeErr = .True.
      isConstraintViolated = .False.
      step = curstep
      iter = 1
      rlt_tol = TOL_E8
      call odeFunc(invars, mem%K1)
      do while (isLargeErr .or. isConstraintViolated)
         if (iter>MAXITER) then
            print *, "Runge-Kutta 3D iteration number is more than 100!!"
            exit
         end if
         curstep = step
         mem%interim = invars + step*0.25*mem%K1
         call odeFunc(mem%interim, mem%K2)
         mem%interim = invars + step*(0.09375*mem%K1 + 0.28125*mem%K2)
         call odeFunc(mem%interim, mem%K3)
         mem%interim = invars + step*(0.87938*mem%K1 - 3.27720*mem%K2 + &
                     3.32089*mem%K3)
         call odeFunc(mem%interim, mem%K4)
         mem%interim = invars + step*(2.03241*mem%K1 - 8.0*mem%K2 + &
                     7.17349*mem%K3 - 0.20590*mem%K4)
         call odeFunc(mem%interim, mem%K5)
         mem%nxt4th = invars + step*(0.11574*mem%K1 + 0.54893*mem%K3 + &
                     0.53533*mem%K4 - 0.2*mem%K5)
         if (mode==fixed_mode) then       ! fixed step
            nextstep = step
            outvars = mem%nxt4th
            return
         end if
         mem%interim = invars + step*(-0.29630*mem%K1 + 2.0*mem%K2 - &
                     1.38168*mem%K3 + 0.45297*mem%K4 - 0.275*mem%K5)
         call odeFunc(mem%interim, mem%K6)
         mem%nxt5th = invars + step*(0.11852*mem%K1 + 0.51899*mem%K3 + &
                     0.50613*mem%K4 - 0.18*mem%K5 + 0.03636*mem%K6)
         call Divide(mem%nxt4th-mem%nxt5th, mem%nxt4th, mem%rerr)
         dy = 0.0_r8
         rdy = 0.0_r8
         dyn = 1.0d36
         do kk = 1, ndim2, 1
            call Norm(mem%rerr(:,kk,:), 1, rdytmp)
            call Norm(mem%nxt4th(:,kk,:)-mem%nxt5th(:,kk,:), 1, dytmp)
            call GetLeastValue(mem%nxt4th(:,kk,:), 1, dyntmp)
            do ii = 1, ndim1, 1
               rdy(ii) = max(rdytmp(ii),rdy(ii))
               dy(ii) = max(dytmp(ii),dy(ii))
               dyn(ii) = min(dyntmp(ii),dyn(ii))
            end do
         end do
         isLargeErr = IsErrorStillLarge(dy,tol,rdy,rlt_tol)
         isConstraintViolated = IsTooLargeNegative(dyn,-100*tol)
         if (isConstraintViolated) then
            step = 0.5*step
         else
            call Divide(tol,dy,abs_rate)
            call Divide(rlt_tol,rdy,rlt_rate)
            rate = max(minval(abs_rate), minval(rlt_rate))
            delta = 0.84*rate**0.25
            if (delta<=0.1) then
               step = 0.1*step
            else if (delta>=4.0) then
               step = 4.0*step
            else
               step = delta*step
            end if
         end if
         iter = iter + 1
      end do
      nextstep = step
      outvars = mem%nxt4th
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Sorting and search method
   !  example: result = BSEARCHQQ(LOC(target),LOC(array),length,SRT$INTEGER4)
   !  example: Call SORTQQ (LOC(array), length, SRT$INTEGER2)
   !
   !------------------------------------------------------------------------------
   subroutine DichotomySectionSearch(array, obj, idx)
      implicit none
      real(r8), intent(in) :: array(:)
      real(r8), intent(in) :: obj
      integer, intent(out) :: idx
      integer :: middle, first, last
      logical :: ascend

      first = 1
      last = size(array)
      ascend = (array(first)<array(last))
      if (ascend) then
         if (obj<=array(first)) then
            idx = 1
            return
         else if (obj>=array(last)) then
            idx = last
            return
         end if
      else
         if (obj>=array(first)) then
            idx = 1
            return
         else if (obj<=array(last)) then
            idx = last
            return
         end if
      end if
      do while (last>first)
         middle = (first+last)/2
         if (array(middle)==obj) then
            last = middle
            exit
         else if (array(middle)<obj) then
            if (ascend) then
               first = middle + 1
            else
               last = middle
            end if
         else
            if (ascend) then
               last = middle
            else
               first = middle + 1
            end if
         end if
      end do
      idx = last - 1
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: solve nonlinear equation by Newton method.
   !
   !------------------------------------------------------------------------------
   subroutine SolveNonlinearEquation(ofunc, odfunc, x0, root)
      implicit none
      external :: ofunc                ! function object
      external :: odfunc               ! function derivative object 
      real(r8), intent(in) :: x0       ! initial value
      real(r8), intent(out) :: root    ! real root
      real(r8) :: fVal, fdVal, xcur, dx
      integer :: iter

      xcur = x0
      dx = 1.0_r8
      iter = 1
      do while (dx>1.0d-6)
         if (iter>MAXITER) then
            print *, "Newton iteration number is more than 100!!" 
            root = NaN 
            exit
         end if
         call ofunc(xcur, fVal)
         call odfunc(xcur, fdVal)
         root = xcur - fVal/fdVal
         if (abs(root)<1) then
            dx = abs( root - xcur )
         else
            dx = abs( 1.0d0 - xcur/root ) 
         end if
         xcur = root
         iter = iter + 1
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: solve quadratic equation a*x^2 + b*x + c = 0 for real roots.
   !
   !------------------------------------------------------------------------------
   subroutine SolveQuadraticEquation(A, root)
      implicit none
      real(r8), intent(in) :: A(3)
      real(r8), intent(out) :: root(2)
      real(r8) :: vD

      vD = A(2)**2.0 - 4.0*A(1)*A(3)
      ! If vD>0, quadratic equation has 2 real roots; vD==0, it has 1 real root
      if (vD>=0) then
         root(1) = 0.5 * (-A(2) + sqrt(vD)) / A(1) 
         root(2) = 0.5 * (-A(2) - sqrt(vD)) / A(1)
      else
         root(1) = inf
         root(2) = inf
      end if
   end subroutine
   
   !------------------------------------------------------------------------------
   !
   ! Purpose: solve cubic equation x^3 + a*(x^2) + b*x + c = 0 for real roots.
   !
   !------------------------------------------------------------------------------
   subroutine SolveCubicEquation(A, root)
      implicit none
      real(r8), intent(in) :: A(3)
      real(r8), intent(out) :: root(3)
      real(r8) :: vQ, vR, vD, rsum, rdif
      complex(r8) :: vS, vT

      vQ = (3.0*A(2) - A(1)**2) / 9.0
      vR = (9.0*A(1)*A(2) - 27.0*A(3) - 2.0*A(1)**3) / 54.0
      vD = vR**2 + vQ**3
      ! If vD<=0, cubic equation has 3 real roots; vD>0, it has 1 real root
      if (vD<=0) then
         vS = cmplx(vR, sqrt(-vD))
         vT = cmplx(vR, -sqrt(-vD))
         vS = vS**(1.0/3.0)
         vT = vT**(1.0/3.0)
         rsum = real(vS+vT)
         rdif = aimag(vS-vT)
         root(1) = rsum - A(1)/3.0
         root(2) = -0.5*rsum - A(1)/3.0 - 0.5*sqrt(3.0)*rdif
         root(3) = -0.5*rsum - A(1)/3.0 + 0.5*sqrt(3.0)*rdif
      else
         root(1) = (vR+sqrt(vD))**(1.0/3.0) + (vR-sqrt(vD))**(1.0/3.0) &
                     - A(1)/3.0
         root(2) = inf
         root(3) = inf
      end if
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: 1-D and 2-D interpolation methods 
   !          Interp1d and Interp2d are linear interpolation
   !          Interp1dcubic is Monotone cubic interpolation
   !
   ! https://en.wikipedia.org/wiki/Monotone_cubic_interpolation
   ! https://en.wikipedia.org/wiki/Cubic_Hermite_spline
   !
   !------------------------------------------------------------------------------
   subroutine Interp1d(xx, yy, xi, yi, vmin, vmax)
      implicit none
      real(r8), intent(in) :: xx(:)    ! defined point coordinate vector
      real(r8), intent(in) :: yy(:)    ! defined point value vector
      real(r8), intent(in) :: xi(:)    ! target point coordinate vector
      real(r8), intent(out) :: yi(:)   ! target point value vector
      real(r8), optional, intent(in) :: vmin
      real(r8), optional, intent(in) :: vmax
      integer :: ii, ntgt, ndef, idx
      real(r8) :: par

      ntgt = size(xi)
      ndef = size(xx)
      do ii = 1, ntgt, 1
         if (xi(ii)<=xx(1)) then
            if (present(vmin)) then
               yi(ii) = vmin
            else
               yi(ii) = yy(1)
            end if
         else if (xi(ii)>=xx(ndef)) then
            if (present(vmax)) then
               yi(ii) = vmax
            else
               yi(ii) = yy(ndef)
            end if
         else
            call DichotomySectionSearch(xx, xi(ii), idx)
            par = (xi(ii) - xx(idx)) / (xx(idx+1) - xx(idx))
            yi(ii) = yy(idx) * (1.0 - par) + yy(idx+1) * par
         end if
      end do
   end subroutine

   subroutine Interp2d(xx, yy, xi, yi, vmin, vmax)
      implicit none
      real(r8), intent(in) :: xx(:)    ! defined point coordinate vector
      real(r8), intent(in) :: yy(:,:)  ! defined point value time series
      real(r8), intent(in) :: xi(:)    ! target point coordinates
      real(r8), intent(out) :: yi(:,:) ! target point value time series
      real(r8), optional, intent(in) :: vmin
      real(r8), optional, intent(in) :: vmax
      integer :: ii, ntgt, ndef, idx
      real(r8) :: par

      ntgt = size(xi)
      ndef = size(xx)
      do ii = 1, ntgt, 1
         if (xi(ii)<=xx(1)) then
            if (present(vmin)) then
               yi(ii,:) = vmin
            else
               yi(ii,:) = yy(1,:)
            end if
         else if (xi(ii)>=xx(ndef)) then
            if (present(vmax)) then
               yi(ii,:) = vmax
            else
               yi(ii,:) = yy(ndef,:)
            end if
         else
            call DichotomySectionSearch(xx, xi(ii), idx)
            par = (xi(ii) - xx(idx)) / (xx(idx+1) - xx(idx))
            yi(ii,:) = yy(idx,:) * (1.0 - par) + yy(idx+1,:) * par
         end if
      end do
   end subroutine

   subroutine Interp1dcubic(xx, yy, xi, yi)
      implicit none
      real(r8), intent(in) :: xx(:)    ! defined point coordinate vector
      real(r8), intent(in) :: yy(:)    ! defined point value vector
      real(r8), intent(in) :: xi(:)    ! target point coordinate vector
      real(r8), intent(out) :: yi(:)   ! target point value vector
      integer :: ii, ntgt, ndef, idx
      real(r8) :: par, dkl, dkm, dkh
      real(r8) :: mkl, mkh, xd
      real(r8) :: h00, h01, h10, h11

      ntgt = size(xi)
      ndef = size(xx)
      do ii = 1, ntgt, 1
         if (xi(ii)<=xx(1)) then
            yi(ii) = yy(1)
         else if (xi(ii)>=xx(ndef)) then
            yi(ii) = yy(ndef)
         else
            call DichotomySectionSearch(xx, xi(ii), idx)
            if (idx==1) then
               dkm = (yy(idx+1) - yy(idx)) / (xx(idx+1) - xx(idx))
               dkh = (yy(idx+2) - yy(idx+1)) / (xx(idx+2) - xx(idx+1))
               mkl = dkm
               if (dkm*dkh<=0) then
                  mkh = 0.0
               else
                  mkh = min( min(0.5*(dkm+dkh), 3.0*dkh), 3.0*dkm )
               end if
            else if (idx==ndef-1) then
               dkl = (yy(idx) - yy(idx-1)) / (xx(idx) - xx(idx-1))
               dkm = (yy(idx+1) - yy(idx)) / (xx(idx+1) - xx(idx))
               mkh = dkm
               if (dkl*dkm<=0) then
                  mkl = 0.0
               else
                  mkl = min( min(0.5*(dkl+dkm), 3.0*dkm), 3.0*dkl )
               end if
            else
               dkl = (yy(idx) - yy(idx-1)) / (xx(idx) - xx(idx-1))
               dkm = (yy(idx+1) - yy(idx)) / (xx(idx+1) - xx(idx))
               dkh = (yy(idx+2) - yy(idx+1)) / (xx(idx+2) - xx(idx+1))
               if (dkl*dkm<=0) then
                  mkl = 0.0
               else
                  mkl = min( min(0.5*(dkl+dkm), 3.0*dkm), 3.0*dkl )
               end if
               if (dkm*dkh<=0) then
                  mkh = 0.0
               else
                  mkh = min( min(0.5*(dkm+dkh), 3.0*dkh), 3.0*dkm )
               end if
            end if
            xd = xx(idx+1) - xx(idx)
            par = (xi(ii) - xx(idx)) / (xx(idx+1) - xx(idx))
            h00 = 2*par**3 - 3*par**2 + 1
            h10 = par**3 - 2*par**2 + par
            h01 = -2*par**3 + 3*par**2
            h11 = par**3 - par**2
            yi(ii) = yy(idx)*h00 + xd*mkl*h10 + yy(idx+1)*h01 + xd*mkh*h11
         end if
      end do
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Generate sampling space using Latin Hypercube method
   !  randperm: random permutation
   !
   !------------------------------------------------------------------------------
   subroutine LHSampling(priors, nparam, nsample, samples)
      implicit none
      real(r8), intent(in) :: priors(:,:)    ! uniform prior bounds
      integer, intent(in) :: nparam          ! number of parameters
      integer, intent(in) :: nsample         ! number of samples
      real(r8), intent(out) :: samples(:,:)  ! generated samples
      real(r8), allocatable :: rands(:,:)    ! random numbers
      real(r8), allocatable :: param(:)
      integer, allocatable :: perm(:)
      integer :: ii

      ! Generate uniform distributed 
      allocate(rands(nsample,nparam))
      call random_number(rands)
      
      allocate(perm(nsample))
      allocate(param(nsample))
      do ii = 1, nparam, 1
         call randperm(nsample,perm)
         param = (perm - rands(:,ii))/nsample
         samples(:,ii) = priors(ii,1) + (priors(ii,2)-priors(ii,1))*param
      end do
      deallocate(param)
      deallocate(perm)
      deallocate(rands)
   end subroutine

   subroutine randperm(num, perm)
      implicit none
      integer, intent(in) :: num
      integer, intent(out) :: perm(:)
      integer :: ii, jj, kk, tmp
      real(r8) :: rrand

      perm = (/(ii, ii=1,num)/)
      do jj = num, 2, -1
         call random_number(rrand)
         kk = floor(jj*rrand) + 1
         ! exchange perm(kk) and perm(jj)
         tmp = perm(kk)
         perm(kk) = perm(jj)
         perm(jj) = tmp
      end do
   end subroutine

   subroutine init_random_seed_old()
      integer :: i, n, clock
      integer, allocatable :: seed(:)

      call RANDOM_SEED(size = n)
      allocate(seed(n))

      call SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call RANDOM_SEED(PUT = seed)

      deallocate(seed)
   end subroutine

   subroutine init_random_seed()
      implicit none
      integer, allocatable :: seed(:)
      integer :: i, n, un, istat, dt(8), pid, t(2), s
      integer(8) :: ncount, tms

      call random_seed(size = n)
      allocate(seed(n))
      ! First try if the OS provides a random number generator
      open(newunit=un, file="/dev/urandom", access="stream", &
           form="unformatted", action="read", status="old", iostat=istat)
      if (istat == 0) then
         read(un) seed
         close(un)
      else
         ! Fallback to XOR:ing the current time and pid. The PID is
         ! useful in case one launches multiple instances of the same
         ! program in parallel.
         call system_clock(ncount)
         if (ncount /= 0) then
            t = transfer(ncount, t)
         else
            call date_and_time(values=dt)
            tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
                 + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
                 + dt(3) * 24 * 60 * 60 * 60 * 1000 &
                 + dt(5) * 60 * 60 * 1000 &
                 + dt(6) * 60 * 1000 + dt(7) * 1000 &
                 + dt(8)
            t = transfer(tms, t)
         end if
         s = ieor(t(1), t(2))
         pid = GETPID() + 1099279 ! Add a prime
         s = ieor(s, pid)
         if (n >= 3) then
            seed(1) = t(1) + 36269
            seed(2) = t(2) + 72551
            seed(3) = pid
            if (n > 3) then
               seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
            end if
         else
            seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
         end if
      end if
      call random_seed(put=seed)
      deallocate(seed)
   end subroutine init_random_seed

   subroutine trapezoid_random_generator(ub, randn)
      implicit none
      integer, intent(in) :: ub
      integer, intent(out) :: randn
      real(r8) :: rrand
      
      call random_number(rrand)
      randn = 1 + int(ub + 0.5 - sqrt( (ub+0.5)**2 - ub*(ub+1)*rrand ))
   end subroutine

   subroutine gaussian_random_generator(avg, std, rand)
      implicit none
      real(r8), intent(in) :: avg
      real(r8), intent(in) :: std
      real(r8), intent(out) :: rand
      real(r8) :: rrand1, rrand2, fpi

      ! algorithm from http://pastebin.com/4bjBvZAD
      fpi = 4.0d+0 * atan(1.0d+0)   ! definition of pi
      call random_number(rrand1)
      call random_number(rrand2)
      rand = std * sqrt(-2.0_r8*log(rrand1)) * cos(2.0_r8*fpi*rrand2) + avg
   end subroutine

   !------------------------------------------------------------------------------
   !
   ! Purpose: Check population convergence
   !
   !------------------------------------------------------------------------------
   function IsPopulationConverged(samples, priors)
      implicit none
      real(r8), intent(in) :: samples(:,:)
      real(r8), intent(in) :: priors(:,:)
      logical :: IsPopulationConverged
      real(r8), parameter :: criteria = 1.0d-3
      real(r8), parameter :: delta = 1.0d-20
      real(r8), dimension(size(priors,1)) :: xmax, xmin, bound
      real(r8), dimension(size(priors,1)) :: xmean, xstd
      integer :: nsample, nparam, ii, jj
      real(r8) :: gsum, gnrng, xsum1, xsum2

      nsample = size(samples, 1)
      nparam = size(samples, 2)
      do ii = 1, nparam, 1
         bound(ii) = priors(ii,2) - priors(ii,1)
      end do
      gsum = 0.0_r8
      do ii = 1, nparam, 1
         xmax(ii) = maxval(samples(:,ii))
         xmin(ii) = minval(samples(:,ii))
         xsum1 = sum(samples(:,ii))
         xsum2 = 0.0_r8
         do jj = 1, nsample, 1
            xsum2 = xsum2 + samples(jj,ii)*samples(jj,ii)
         end do
         xmean(ii) = xsum1 / dble(nsample)
         xstd(ii) = xsum2 / dble(nsample) - xmean(ii)*xmean(ii)
         if (xstd(ii)<delta) then
            xstd(ii) = delta
         end if
         xstd(ii) = sqrt(xstd(ii))
         xstd(ii) = xstd(ii) / bound(ii)
         gsum = gsum + log( delta + (xmax(ii)-xmin(ii))/bound(ii) )
      end do
      gnrng = exp(gsum/dble(nparam))
      if (gnrng<=criteria) then
         IsPopulationConverged = .True.
      else
         IsPopulationConverged = .False.
      end if
      return
   end function

   subroutine CalcPopulationSTD(samples, priors, std)
      implicit none
      real(r8), intent(in) :: samples(:,:)
      real(r8), intent(in) :: priors(:,:)
      real(r8), intent(out) :: std(:)
      real(r8), parameter :: delta = 1.0d-20
      integer :: nsample, nparam, ii, jj
      real(r8) :: xmean, xmax, xmin, bound
      real(r8) :: xsum1, xsum2

      nsample = size(samples, 1)
      nparam = size(samples, 2)
      do ii = 1, nparam, 1
         bound = priors(ii,2) - priors(ii,1)
         xmax = maxval(samples(:,ii))
         xmin = minval(samples(:,ii))
         xsum1 = sum(samples(:,ii))
         xsum2 = 0.0_r8
         do jj = 1, nsample, 1
            xsum2 = xsum2 + samples(jj,ii)*samples(jj,ii)
         end do
         xmean = xsum1 / dble(nsample)
         std(ii) = xsum2 / dble(nsample) - xmean*xmean
         if (std(ii)<delta) then
            std(ii) = delta
         end if
         std(ii) = sqrt(std(ii))
         std(ii) = std(ii) / bound
      end do
   end subroutine

end module math_utilities_mod
