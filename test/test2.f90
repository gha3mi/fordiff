!===============================================================================
module mod_func2

   use fordiff, only: rk

   implicit none

contains

   function func2(x) result(f)
      complex(rk), dimension(:), intent(in) :: x
      complex(rk)                           :: f

      f = x(1)**2 + 0.5_rk*x(2)**2

   end function func2

end module mod_func2
!===============================================================================



!===============================================================================
program test2

   use fordiff, only: rk, derivative
   use mod_func2, only: func2
   use forunittest, only: unit_test

   implicit none

   real(rk), dimension(:), allocatable :: dfdx, expected_dfdx
   integer                             :: i
   type(unit_test) :: ut

   ! compute derivative
   dfdx = derivative(f=func2, x=[1.0_rk, -1.0_rk], h=tiny(0.0_rk))

   ! compute reference derivative
   allocate(expected_dfdx(size(dfdx)))
   expected_dfdx(1) = 2.0_rk*(1.0_rk)
   expected_dfdx(2) = 1.0_rk*(-1.0_rk)

   ! check if derivative is correct
   call ut%check(dfdx, expected_dfdx, 1.0e-6_rk, 'test2' )

   deallocate(dfdx, expected_dfdx)

end program test2
!===============================================================================