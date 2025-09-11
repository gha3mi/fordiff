!===============================================================================
module mod_func1

   use kinds, only: rk
   implicit none

contains

   function func1(x) result(f)
      complex(rk), intent(in)  :: x
      complex(rk)              :: f

      f = x**2 + 2.0_rk*x

   end function func1

end module mod_func1
!===============================================================================



!===============================================================================
program test1

   use kinds, only: rk
   use fordiff, only: derivative
   use mod_func1, only: func1
   use forunittest, only: unit_test

   implicit none

   real(rk) :: dfdx, expected_dfdx
   type(unit_test) :: ut

   ! compute derivative
   dfdx = derivative(f=func1, x=1.0_rk, h=tiny(0.0_rk))

   ! reference value
   expected_dfdx = 2.0_rk*(1.0_rk) + 2.0_rk

   ! check if derivative is correct
   call ut%check(dfdx, expected_dfdx, 1.0e-6_rk, 'test1' )

end program test1
!===============================================================================