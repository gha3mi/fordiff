program demo2

   ! This program demonstrates how to compute the derivative of a scalar-valued function w.r.t to a vector variable
   ! using complex-step and finite-difference methods.

   use fordiff, only: rk, derivative

   implicit none

   real(rk), dimension(:), allocatable :: dfdx

   print'(a)', 'f(x_1,x_2) = x_1**2 + 0.5*x_2**2'
   print'(a)', 'compute derivative of f(x_1,x_2) w.r.t. x_1 and x_2 at x_1 = 1, x_2 = -1'

   ! Compute derivative of a function using complex-step
   dfdx = derivative(f=f1, x=[1.0_rk, -1.0_rk], h=tiny(0.0_rk))
   print'(a,g0,", ",g0,a)', 'dfdx = ', dfdx, ' (complex-step)'

   ! Compute derivative of a function using forward finite-difference
   dfdx = derivative(f=f2, x=[1.0_rk, -1.0_rk], h=1e-5_rk, method='forward')
   print'(a,g0,", ",g0,a)', 'dfdx = ', dfdx, ' (forward finite-difference)'

   ! Compute derivative of a function using backward finite-difference
   dfdx = derivative(f=f2, x=[1.0_rk, -1.0_rk], h=1e-5_rk, method='backward')
   print'(a,g0,", ",g0,a)', 'dfdx = ', dfdx, ' (backward finite-difference)'

   ! Compute derivative of a function using central finite-difference
   dfdx = derivative(f=f2, x=[1.0_rk, -1.0_rk], h=1e-5_rk, method='central')
   print'(a,g0,", ",g0,a)', 'dfdx = ', dfdx, ' (central finite-difference)'

contains

   ! Define a scalar function of a vector variable (complex)
   function f1(x) result(f)
      complex(rk), dimension(:), intent(in) :: x
      complex(rk)                           :: f
      f = x(1)**2 + 0.5_rk*x(2)**2
   end function f1

   ! Define a scalar function of a vector variable (real)
   function f2(x) result(f)
      real(rk), dimension(:), intent(in) :: x
      real(rk)                           :: f
      f = x(1)**2 + 0.5_rk*x(2)**2
   end function f2

end program demo2
