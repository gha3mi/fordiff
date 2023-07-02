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

   use kinds
   use mod_func6
   use fordiff

   implicit none

   real(rk), dimension(:,:), allocatable :: dfdx
   integer                               :: i

   dfdx = derivative(f=func6, x=[1.0_rk, -1.0_rk], h=1e-10_rk, method='forward')

   do i = 1, size(dfdx,1)
      print*, (dfdx(i, :))
   end do

   dfdx = derivative(f=func6, x=[1.0_rk, -1.0_rk], h=1e-10_rk, method='backward')

   do i = 1, size(dfdx,1)
      print*, (dfdx(i, :))
   end do

   dfdx = derivative(f=func6, x=[1.0_rk, -1.0_rk], h=1e-10_rk, method='central')

   do i = 1, size(dfdx,1)
      print*, (dfdx(i, :))
   end do

   deallocate(dfdx)

end program test6
!===============================================================================