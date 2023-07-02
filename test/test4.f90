!===============================================================================
module mod_func4

   use kinds
   implicit none

contains

   function func4(x) result(f)
      real(rk), intent(in) :: x
      real(rk)             :: f

      f = x**2

   end function func4

end module mod_func4
!===============================================================================



!===============================================================================
program test4

   use kinds
   use mod_func4
   use fordiff

   implicit none

   real(rk) :: dfdx

   dfdx = derivative(f=func4, x=1.0_rk, h=1e-10_rk, method='forward')
   print*, dfdx
   
   dfdx = derivative(f=func4, x=1.0_rk, h=1e-10_rk, method='backward')
   print*, dfdx

   dfdx = derivative(f=func4, x=1.0_rk, h=1e-10_rk, method='central')
   print*, dfdx

end program test4
!===============================================================================