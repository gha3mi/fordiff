!===============================================================================
module mod_func5

   use fordiff, only: rk

   implicit none

contains

   function func5(x) result(f)
      real(rk), dimension(:), intent(in) :: x
      real(rk)                           :: f

      f = x(1)**2 + 0.5_rk*x(2)**2

   end function func5

end module mod_func5
!===============================================================================



!===============================================================================
program test5

   use fordiff, only: rk, derivative
   use mod_func5, only: func5
   use forunittest, only: unit_test

   implicit none

   real(rk), dimension(:), allocatable :: dfdx_f, dfdx_b, dfdx_c, expected_dfdx
   integer                             :: i
   type(unit_test) :: ut

   dfdx_f = derivative(f=func5, x=[1.0_rk, -1.0_rk], h=1e-5_rk, method='forward')

   dfdx_b = derivative(f=func5, x=[1.0_rk, -1.0_rk], h=1e-5_rk, method='backward')

   dfdx_c = derivative(f=func5, x=[1.0_rk, -1.0_rk], h=1e-5_rk, method='central')

   ! compute reference derivative
   allocate(expected_dfdx(size(dfdx_f)))
   expected_dfdx(1) = 2.0_rk*(1.0_rk)
   expected_dfdx(2) = 1.0_rk*(-1.0_rk)

   ! check if derivative is correct
   call ut%check(dfdx_f, expected_dfdx, 1.0e-2_rk, 'test5.1' )
   call ut%check(dfdx_b, expected_dfdx, 1.0e-2_rk, 'test5.1' )
   call ut%check(dfdx_c, expected_dfdx, 1.0e-2_rk, 'test5.1' )

   deallocate(dfdx_f, dfdx_b, dfdx_c, expected_dfdx)

end program test5
!===============================================================================