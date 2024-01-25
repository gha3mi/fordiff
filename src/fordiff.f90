module fordiff
!! author: Seyed Ali Ghasemi
!! license: BSD 3-Clause License
!! Module for numerical differentiation using complex step differentiation and finite difference methods
!!

   use kinds  !! for real(kind=rk) and complex(kind=rk). Use -DREAL32 for real kind 4 and -DREAL64 for real kind 8. Default is real kind 8.

   implicit none

   private
   public :: derivative

   !===============================================================================
   interface derivative
      !! Derivative interface
      procedure :: complex_step_derivative_T0_T0
      procedure :: complex_step_derivative_T0_T1
      procedure :: complex_step_derivative_T1_T1
      procedure :: finite_difference_T0_T0
      procedure :: finite_difference_T0_T1
      procedure :: finite_difference_T1_T1
   end interface
   !===============================================================================

contains

   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a scalar-valued function f w.r.t. a scalar-valued variable x
   !> using complex step differentiation.
   impure function complex_step_derivative_T0_T0(f, x, h) result(dfdx)
      real(rk), intent(in) :: x    !! scalar variable
      real(rk), intent(in) :: h    !! perturbation for complex step differentiation
      real(rk)             :: dfdx !! derivative of f w.r.t. x

      interface
         !! scalar-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            complex(rk), intent(in) :: z  !! scalar complex variable
            complex(rk)             :: fz !! scalar complex function
         end function f
      end interface

      ! Check for potential division by zero
      if (h == 0.0_rk) error stop 'Division by zero. Please provide a non-zero value for h.'

      dfdx = aimag(f(cmplx(x, h, rk))) / h

   end function complex_step_derivative_T0_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a scalar-valued function f
   !> w.r.t. a vector-valued variable x using complex step differentiation.
   impure function complex_step_derivative_T0_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in) :: x      !! vector variable
      real(rk),               intent(in) :: h      !! perturbation for complex step differentiation
      real(rk), dimension(size(x))       :: dfdx   !! derivative of f w.r.t. x
      real(rk), dimension(size(x))       :: temp_x !! temporary vector variable
      integer                            :: i      !! loop index

      interface
         !! scalar-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            complex(rk), dimension(:), intent(in) :: z  !! vector complex variable
            complex(rk)                           :: fz !! scalar complex function
         end function f
      end interface

      if (h == 0.0_rk) error stop 'Division by zero. Please provide a non-zero value for h.'

      do i = 1, size(x)
         temp_x    = 0.0_rk
         temp_x(i) = x(i)
         dfdx(i) = aimag(f(cmplx(temp_x, h, rk))) / h
      end do

   end function complex_step_derivative_T0_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a vector-valued function f
   !> w.r.t. a vector-valued variable x using complex step differentiation.
   impure function complex_step_derivative_T1_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in)    :: x      !! vector variable
      real(rk),               intent(in)    :: h      !! perturbation for complex step differentiation
      real(rk), dimension(:,:), allocatable :: dfdx   !! derivative of f w.r.t. x
      real(rk), dimension(size(x))          :: temp_x !! temporary vector variable
      integer                               :: i      !! loop index

      interface
         !! vector-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            complex(rk), dimension(:), intent(in)  :: z  !! vector complex variable
            complex(rk), dimension(:), allocatable :: fz !! vector complex function
         end function f
      end interface

      allocate(dfdx(size(f(cmplx(x,kind=rk))),size(x)))

      if (h == 0.0_rk) error stop 'Division by zero. Please provide a non-zero value for h.'

      do i = 1, size(x)
         temp_x    = 0.0_rk
         temp_x(i) = x(i)
         dfdx(:,i) = aimag(f(cmplx(temp_x, h, rk))) / h
      end do

   end function complex_step_derivative_T1_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a scalar-valued function f w.r.t. a scalar-valued variable x
   !> using finite difference methods (forward, backward, central).
   impure function finite_difference_T0_T0(f,x,h,method) result(dfdx)
      real(rk),     intent(in) :: x      !! scalar variable
      real(rk),     intent(in) :: h      !! perturbation for finite difference methods
      character(*), intent(in) :: method !! finite difference method (forward, backward, central)
      real(rk)                 :: dfdx   !! derivative of f w.r.t. x

      interface
         !! scalar-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            real(rk), intent(in) :: z  !! scalar variable
            real(rk)             :: fz !! scalar function
         end function f
      end interface

      if (h == 0.0_rk) error stop 'Division by zero. Please provide a non-zero value for h.'

      select case (method)
       case('forward')
         dfdx = finite_difference_forward_T0_T0(f,x,h)
       case('backward')
         dfdx = finite_difference_backward_T0_T0(f,x,h)
       case('central')
         dfdx = finite_difference_central_T0_T0(f,x,h)
      end select

   end function finite_difference_T0_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a scalar-valued function f w.r.t. a scalar-valued variable x
   !> using the central finite difference method.
   impure function finite_difference_central_T0_T0(f, x, h) result(dfdx)
      real(rk), intent(in) :: x    !! scalar variable
      real(rk), intent(in) :: h    !! perturbation for finite difference methods
      real(rk)             :: dfdx !! derivative of f w.r.t. x

      interface
         !! scalar-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            real(rk), intent(in) :: z  !! scalar variable
            real(rk)             :: fz !! scalar function
         end function f
      end interface

      dfdx = (f(x+h) - f(x-h)) / (2.0_rk * h)

   end function finite_difference_central_T0_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a scalar-valued function f w.r.t. a scalar-valued variable x
   !> using the forward finite difference method.
   impure function finite_difference_forward_T0_T0(f, x, h) result(dfdx)
      real(rk), intent(in) :: x    !! scalar variable
      real(rk), intent(in) :: h    !! perturbation for finite difference methods
      real(rk)             :: dfdx !! derivative of f w.r.t. x

      interface
         !! scalar-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            real(rk), intent(in) :: z  !! scalar variable
            real(rk)             :: fz !! scalar function
         end function f
      end interface

      dfdx = (f(x+h) - f(x)) / h

   end function finite_difference_forward_T0_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a scalar-valued function f w.r.t. a scalar-valued variable x
   !> using the backward finite difference method.
   impure function finite_difference_backward_T0_T0(f, x, h) result(dfdx)
      real(rk), intent(in) :: x    !! scalar variable
      real(rk), intent(in) :: h    !! perturbation for finite difference methods
      real(rk)             :: dfdx !! derivative of f w.r.t. x

      interface
         !! scalar-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            real(rk), intent(in) :: z  !! scalar variable
            real(rk)             :: fz !! scalar function
         end function f
      end interface

      dfdx = (f(x) - f(x-h)) / h

   end function finite_difference_backward_T0_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a scalar-valued function f w.r.t. a vector-valued variable x
   !> using finite difference methods (forward, backward, central).
   impure function finite_difference_T0_T1(f,x,h,method) result(dfdx)
      real(rk), dimension(:), intent(in) :: x      !! vector variable
      real(rk),               intent(in) :: h      !! perturbation for finite difference methods
      character(*),           intent(in) :: method !! finite difference method (forward, backward, central)
      real(rk), dimension(size(x))       :: dfdx   !! derivative of f w.r.t. x

      interface
         !! scalar-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in) :: z  !! vector variable
            real(rk)                           :: fz !! scalar function        
         end function f
      end interface

      if (h == 0.0_rk) error stop 'Division by zero. Please provide a non-zero value for h.'

      select case (method)
       case('forward')
         dfdx = finite_difference_forward_T0_T1(f,x,h)
       case('backward')
         dfdx = finite_difference_backward_T0_T1(f,x,h)
       case('central')
         dfdx = finite_difference_central_T0_T1(f,x,h)
      end select

   end function finite_difference_T0_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a scalar-valued function f w.r.t. a vector-valued variable x
   !> using the central finite difference method.
   impure function finite_difference_central_T0_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in) :: x       !! vector variable
      real(rk),               intent(in) :: h       !! perturbation for finite difference methods
      real(rk), dimension(size(x))       :: dfdx    !! derivative of f w.r.t. x
      real(rk), dimension(size(x))       :: temp_x  !! temporary vector variable
      integer                            :: i       !! loop index

      interface
         !! scalar-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in) :: z  !! vector variable
            real(rk)                           :: fz !! scalar function
         end function f
      end interface

      do i = 1, size(x)
         temp_x    = 0.0_rk
         temp_x(i) = x(i)
         dfdx(i) = (f(temp_x+h) - f(temp_x-h)) / (2.0_rk * h)
      end do

   end function finite_difference_central_T0_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a scalar-valued function f w.r.t. a vector-valued variable x
   !> using the forward finite difference method.
   impure function finite_difference_forward_T0_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in) :: x       !! vector variable
      real(rk),               intent(in) :: h       !! perturbation for finite difference methods
      real(rk), dimension(size(x))       :: dfdx    !! derivative of f w.r.t. x
      real(rk), dimension(size(x))       :: temp_x  !! temporary vector variable
      integer                            :: i       !! loop index

      interface
         !! scalar-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in) :: z  !! vector variable
            real(rk)                           :: fz !! scalar function
         end function f
      end interface

      do i = 1, size(x)
         temp_x    = 0.0_rk
         temp_x(i) = x(i)
         dfdx(i) = (f(temp_x+h) - f(temp_x)) / h
      end do

   end function finite_difference_forward_T0_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a scalar-valued function f w.r.t. a vector-valued variable x
   !> using the backward finite difference method.
   impure function finite_difference_backward_T0_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in) :: x      !! vector variable
      real(rk),               intent(in) :: h      !! perturbation for finite difference methods
      real(rk), dimension(size(x))       :: dfdx   !! derivative of f w.r.t. x
      real(rk), dimension(size(x))       :: temp_x !! temporary vector variable
      integer                            :: i      !! loop index

      interface
         !! scalar-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in) :: z  !! vector variable
            real(rk)                           :: fz !! scalar function
         end function f
      end interface

      do i = 1, size(x)
         temp_x    = 0.0_rk
         temp_x(i) = x(i)
         dfdx(i) = (f(temp_x) - f(temp_x-h)) / h
      end do

   end function finite_difference_backward_T0_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a vector-valued function f w.r.t. a vector-valued variable x
   !> using finite difference methods (forward, backward, central).
   impure function finite_difference_T1_T1(f,x,h,method) result(dfdx)
      real(rk), dimension(:), intent(in)    :: x       !! vector variable
      real(rk),               intent(in)    :: h       !! perturbation for finite difference methods
      character(*),           intent(in)    :: method  !! finite difference method (forward, backward, central)
      real(rk), dimension(:,:), allocatable :: dfdx    !! derivative of f w.r.t. x

      interface
         !! vector-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in)  :: z  !! vector variable
            real(rk), dimension(:), allocatable :: fz !! vector function
         end function f
      end interface

      if (h == 0.0_rk) error stop 'Division by zero. Please provide a non-zero value for h.'

      select case (method)
       case('forward')
         dfdx = finite_difference_forward_T1_T1(f,x,h)
       case('backward')
         dfdx = finite_difference_backward_T1_T1(f,x,h)
       case('central')
         dfdx = finite_difference_central_T1_T1(f,x,h)
      end select

   end function finite_difference_T1_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a vector-valued function f w.r.t. a vector-valued variable x
   !> using the central finite difference method.
   impure function finite_difference_central_T1_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in)    :: x      !! vector variable
      real(rk),               intent(in)    :: h      !! perturbation for finite difference methods
      real(rk), dimension(:,:), allocatable :: dfdx   !! derivative of f w.r.t. x
      real(rk), dimension(size(x))          :: temp_x !! temporary vector variable
      integer                               :: i      !! loop index

      interface
         !! vector-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in)  :: z  !! vector variable
            real(rk), dimension(:), allocatable :: fz !! vector function
         end function f
      end interface

      allocate(dfdx(size(f(x)),size(x)))

      do i = 1, size(x)
         temp_x    = 0.0_rk
         temp_x(i) = x(i)
         dfdx(:,i) = (f(temp_x+h) - f(temp_x-h)) / (2.0_rk*h)
      end do

   end function finite_difference_central_T1_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a vector-valued function f w.r.t. a vector-valued variable x
   !> using the forward finite difference method.
   impure function finite_difference_forward_T1_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in)    :: x      !! vector variable
      real(rk),               intent(in)    :: h      !! perturbation for finite difference methods
      real(rk), dimension(:,:), allocatable :: dfdx   !! derivative of f w.r.t. x
      real(rk), dimension(size(x))          :: temp_x !! temporary vector variable
      integer                               :: i      !! loop index

      interface
         !! vector-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in)  :: z  !! vector variable
            real(rk), dimension(:), allocatable :: fz !! vector function
         end function f
      end interface

      allocate(dfdx(size(f(x)),size(x)))

      do i = 1, size(x)
         temp_x    = 0.0_rk
         temp_x(i) = x(i)
         dfdx(:,i) = (f(temp_x+h) - f(temp_x)) / h
      end do

   end function finite_difference_forward_T1_T1
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a vector-valued function f w.r.t. a vector-valued variable x
   !> using the backward finite difference method.
   impure function finite_difference_backward_T1_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in)    :: x       !! vector variable
      real(rk),               intent(in)    :: h       !! perturbation for finite difference methods
      real(rk), dimension(:,:), allocatable :: dfdx    !! derivative of f w.r.t. x
      real(rk), dimension(size(x))          :: temp_x  !! temporary vector variable
      integer                               :: i       !! loop index

      interface
         !! vector-valued function to differentiate
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in)  :: z  !! vector variable
            real(rk), dimension(:), allocatable :: fz !! vector function
         end function f
      end interface

      allocate(dfdx(size(f(x)),size(x)))

      do i = 1, size(x)
         temp_x    = 0.0_rk
         temp_x(i) = x(i)
         dfdx(:,i) = (f(temp_x) - f(temp_x-h)) / h
      end do

   end function finite_difference_backward_T1_T1
   !===============================================================================

end module fordiff