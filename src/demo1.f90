program demo1

   ! This program demonstrates how to compute the derivative of a scalar-valued function w.r.t to a scalar variable
   ! using complex-step and finite-difference methods.

   use fordiff, only: rk, derivative

   implicit none

   real(rk) :: dfdx

   print'(a)', 'f(x) = x**2 + 2.0*x'
   print'(a)', 'compute derivative of f(x) w.r.t to x at x = 1.0'

   ! Compute derivative of a function using complex-step
   dfdx = derivative(f=f1, x=1.0_rk, h=tiny(0.0_rk))
   print'(a,g0,a)', 'dfdx = ', dfdx, ' (complex-step)'

   ! Compute derivative of a function using forward finite-difference
   dfdx = derivative(f=f2, x=1.0_rk, h=1e-5_rk, method='forward')
   print'(a,g0,a)', 'dfdx = ', dfdx, ' (forward finite-difference)'

   ! Compute derivative of a function using backward finite-difference
   dfdx = derivative(f=f2, x=1.0_rk, h=1e-5_rk, method='backward')
   print'(a,g0,a)', 'dfdx = ', dfdx, ' (backward finite-difference)'

   ! Compute derivative of a function using central finite-difference
   dfdx = derivative(f=f2, x=1.0_rk, h=1e-5_rk, method='central')
   print'(a,g0,a)', 'dfdx = ', dfdx, ' (central finite-difference)'

contains

   ! Define a scalar function of a scalar variable (complex)
   function f1(x) result(f)
      complex(rk), intent(in) :: x
      complex(rk)             :: f
      f = x**2 + 2.0_rk*x
   end function f1

   ! Define a scalar function of a scalar variable (real)
   function f2(x) result(f)
      real(rk), intent(in) :: x
      real(rk)             :: f
      f = x**2 + 2.0_rk*x
   end function f2

end program demo1
