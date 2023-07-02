!===============================================================================
module mod_func2

   use kinds
   implicit none

contains

   function func2(x) result(f)
      complex(rk), dimension(:), intent(in) :: x
      complex(rk)                           :: f

      f = x(1)**2 + 0.5*x(2)**2

   end function func2

end module mod_func2
!===============================================================================



!===============================================================================
program test2

   use kinds
   use mod_func2
   use fordiff

   implicit none

   real(rk), dimension(:), allocatable :: dfdx
   integer                             :: i

   dfdx = derivative(f=func2, x=[1.0_rk, -1.0_rk], h=1e-100_rk)

   do i = 1, size(dfdx)
      print*, (dfdx(i))
   end do

   deallocate(dfdx)

end program test2
!===============================================================================