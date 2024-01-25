!===============================================================================
module mod_func6

   use kinds
   implicit none

contains

   function func6(x) result(f)
      real(rk), dimension(:), intent(in)  :: x
      real(rk), dimension(:), allocatable :: f

      allocate(f(3))

      f(1) = x(1)**2 + 0.5*x(2)**2
      f(2) = x(1)**3 + 0.5*x(2)**3
      f(3) = x(1)**4 + 0.5*x(2)**4

   end function func6

end module mod_func6
!===============================================================================



!===============================================================================
program test6

   use fordiff
   use mod_func6
   use forunittest

   implicit none

   real(rk), dimension(:,:), allocatable :: dfdx_f, dfdx_b, dfdx_c, expected_dfdx
   integer                               :: i
   type(unit_test) :: ut

   ! compute derivative using forward, backward and central difference
   dfdx_f = derivative(f=func6, x=[1.0_rk, -1.0_rk], h=1e-5_rk, method='forward')
   dfdx_b = derivative(f=func6, x=[1.0_rk, -1.0_rk], h=1e-5_rk, method='backward')
   dfdx_c = derivative(f=func6, x=[1.0_rk, -1.0_rk], h=1e-5_rk, method='central')

   ! compute reference derivative
   allocate(expected_dfdx(size(dfdx_f,1), size(dfdx_f,2)))
   expected_dfdx(1,1) = 2.0_rk*(1.0_rk)
   expected_dfdx(1,2) = (-1.0_rk)
   expected_dfdx(2,1) = 3.0_rk*(1.0_rk)**2
   expected_dfdx(2,2) = 1.5_rk*(-1.0_rk)**2
   expected_dfdx(3,1) = 4.0_rk*(1.0_rk)**3
   expected_dfdx(3,2) = 2.0_rk*(-1.0_rk)**3   

   ! check if derivative is correct
   call ut%check(dfdx_f, expected_dfdx, 1.0e-2_rk, 'test6.1' )
   call ut%check(dfdx_b, expected_dfdx, 1.0e-2_rk, 'test6.1' )
   call ut%check(dfdx_c, expected_dfdx, 1.0e-2_rk, 'test6.1' )

   deallocate(dfdx_f, dfdx_b, dfdx_c, expected_dfdx)

end program test6
!===============================================================================