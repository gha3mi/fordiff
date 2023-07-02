!===============================================================================
module mod_func5

   use kinds
   implicit none

contains

   function func5(x) result(f)
      real(rk), dimension(:), intent(in) :: x
      real(rk)                           :: f

      f = x(1)**2 + 0.5*x(2)**2

   end function func5

end module mod_func5
!===============================================================================



!===============================================================================
program test2

   use kinds
   use mod_func5
   use fordiff

   implicit none

   real(rk), dimension(:), allocatable :: dfdx
   integer                             :: i

   dfdx = derivative(f=func5, x=[1.0_rk, -1.0_rk], h=1e-10_rk, method='forward')
   do i = 1, size(dfdx)
      print*, (dfdx(i))
   end do

   dfdx = derivative(f=func5, x=[1.0_rk, -1.0_rk], h=1e-10_rk, method='backward')
   do i = 1, size(dfdx)
      print*, (dfdx(i))
   end do

   dfdx = derivative(f=func5, x=[1.0_rk, -1.0_rk], h=1e-10_rk, method='central')
   do i = 1, size(dfdx)
      print*, (dfdx(i))
   end do

   deallocate(dfdx)

end program test2
!===============================================================================