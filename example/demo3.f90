program demo3

   ! This program demonstrates how to compute the derivative of a vector-valued function w.r.t. a vector of variables
   ! using complex-step and finite-difference methods.

   use kinds, only: rk
   use fordiff, only: derivative

   implicit none

   real(rk), dimension(:,:), allocatable :: dfdx

   print'(a)', 'f_1(x_1,x_2) = x_1**2 + 0.5*x_2**2'
   print'(a)', 'f_2(x_1,x_2) = x_1**3 + 0.5*x_2**3'
   print'(a)', 'f_3(x_1,x_2) = x_1**4 + 0.5*x_2**4'
   print'(a)', 'compute derivative of f w.r.t. x at x_1 = 1, x_2 = -1'

   ! Compute derivative of a function using complex-step
   dfdx = derivative(f=f1, x=[1.0_rk, -1.0_rk], h=tiny(0.0_rk))
   print'(a,*(g0,", "),a)', 'dfdx = ', dfdx, ' (complex-step)'

   ! Compute derivative of a function using forward finite-difference
   dfdx = derivative(f=f2, x=[1.0_rk, -1.0_rk], h=1e-5_rk, method='forward')
   print'(a,*(g0,", "),a)', 'dfdx = ', dfdx, ' (forward finite-difference)'

   ! Compute derivative of a function using backward finite-difference
   dfdx = derivative(f=f2, x=[1.0_rk, -1.0_rk], h=1e-5_rk, method='backward')
   print'(a,*(g0,", "),a)', 'dfdx = ', dfdx, ' (backward finite-difference)'

   ! Compute derivative of a function using central finite-difference
   dfdx = derivative(f=f2, x=[1.0_rk, -1.0_rk], h=1e-5_rk, method='central')
   print'(a,*(g0,", "),a)', 'dfdx = ', dfdx, ' (central finite-difference)'

contains

   ! Define a vector-valued function of a vector variable (complex)
   function f1(x) result(f)
      complex(rk), dimension(:), intent(in)  :: x
      complex(rk), dimension(:), allocatable :: f
      allocate(f(3))
      f(1) = x(1)**2 + 0.5_rk*x(2)**2
      f(2) = x(1)**3 + 0.5_rk*x(2)**3
      f(3) = x(1)**4 + 0.5_rk*x(2)**4
   end function f1

   ! Define a vector-valued function of a vector variable (real)
   function f2(x) result(f)
      real(rk), dimension(:), intent(in)  :: x
      real(rk), dimension(:), allocatable :: f
      allocate(f(3))
      f(1) = x(1)**2 + 0.5*x(2)**2
      f(2) = x(1)**3 + 0.5*x(2)**3
      f(3) = x(1)**4 + 0.5*x(2)**4
   end function f2

end program demo3
