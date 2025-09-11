!===============================================================================
module mod_func3

   use kinds, only: rk
   implicit none

contains

   function func3(x) result(f)
      complex(rk), dimension(:), intent(in)  :: x
      complex(rk), dimension(:), allocatable :: f

      allocate(f(3))

      f(1) = x(1)**2 + 0.5_rk*x(2)**2
      f(2) = x(1)**3 + 0.5_rk*x(2)**3
      f(3) = x(1)**4 + 0.5_rk*x(2)**4

   end function func3

end module mod_func3
!===============================================================================



!===============================================================================
program test3

   use kinds, only: rk
   use mod_func3, only: func3
   use fordiff, only: derivative
   use forunittest, only: unit_test

   implicit none

   real(rk), dimension(:,:), allocatable :: dfdx, expected_dfdx
   integer                               :: i
   type(unit_test) :: ut

   ! compute derivative
   dfdx = derivative(f=func3, x=[1.0_rk, -1.0_rk], h=tiny(0.0_rk))

   ! reference derivative
   allocate(expected_dfdx(size(dfdx,1), size(dfdx,2)))
   expected_dfdx(1,1) = 2.0_rk*(1.0_rk)
   expected_dfdx(1,2) = (-1.0_rk)
   expected_dfdx(2,1) = 3.0_rk*(1.0_rk)**2
   expected_dfdx(2,2) = 1.5_rk*(-1.0_rk)**2
   expected_dfdx(3,1) = 4.0_rk*(1.0_rk)**3
   expected_dfdx(3,2) = 2.0_rk*(-1.0_rk)**3

   ! check if derivative is correct
   call ut%check(dfdx, expected_dfdx, 1.0e-6_rk, 'test3' )

   deallocate(dfdx, expected_dfdx)

end program test3
!===============================================================================