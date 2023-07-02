!===============================================================================
module mod_func3

   use kinds
   implicit none

contains

   function func3(x) result(f)
      complex(rk), dimension(:), intent(in)  :: x
      complex(rk), dimension(:), allocatable :: f

      allocate(f(3))

      f(1) = x(1)**2 + 0.5*x(2)**2
      f(2) = x(1)**3 + 0.5*x(2)**3
      f(3) = x(1)**4 + 0.5*x(2)**4

   end function func3

end module mod_func3
!===============================================================================



!===============================================================================
program test3

   use kinds
   use mod_func3
   use fordiff

   implicit none

   real(rk), dimension(:,:), allocatable :: dfdx
   integer                               :: i

   dfdx = derivative(f=func3, x=[1.0_rk, -1.0_rk], h=1e-100_rk)

   do i = 1, size(dfdx,1)
      print*, (dfdx(i, :))
   end do
   
   deallocate(dfdx)

end program test3
!===============================================================================