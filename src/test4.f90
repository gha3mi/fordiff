!===============================================================================
module mod_func4

   use kinds, only: rk
   implicit none

contains

   function func4(x) result(f)
      real(rk), intent(in) :: x
      real(rk)             :: f

      f = x**2 + 2.0_rk*x

   end function func4

end module mod_func4
!===============================================================================



!===============================================================================
program test4

   use kinds, only: rk
   use fordiff, only: derivative
   use mod_func4, only: func4
   use forunittest, only: unit_test

   implicit none

   real(rk) :: dfdx_f, dfdx_b, dfdx_c, expected_dfdx
   type(unit_test) :: ut

   ! compute derivative using forward mode
   dfdx_f = derivative(f=func4, x=1.0_rk, h=1e-5_rk, method='forward')
   
   ! compute derivative using backward mode
   dfdx_b = derivative(f=func4, x=1.0_rk, h=1e-5_rk, method='backward')

   ! compute derivative using central mode
   dfdx_c = derivative(f=func4, x=1.0_rk, h=1e-5_rk, method='central')

   ! compute derivative analytically
   expected_dfdx = 2.0_rk*(1.0_rk) + 2.0_rk

   ! check if derivative is correct
   call ut%check(dfdx_f, expected_dfdx, 1.0e-2_rk, 'test4.1' )
   call ut%check(dfdx_b, expected_dfdx, 1.0e-2_rk, 'test4.2' )
   call ut%check(dfdx_c, expected_dfdx, 1.0e-2_rk, 'test4.3' )

end program test4
!===============================================================================