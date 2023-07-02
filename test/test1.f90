!===============================================================================
module mod_func1

   use kinds
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

   use kinds
   use mod_func1
   use fordiff

   implicit none

   real(rk) :: dfdx

   dfdx = derivative(f=func1, x=1.0_rk, h=1e-100_rk)

   print*,dfdx

end program test1
!===============================================================================