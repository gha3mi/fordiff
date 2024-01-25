module fordiff

   use kinds

   implicit none

   private
   public :: derivative

   !===============================================================================
   interface derivative
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
   !> Calculates the derivative of a scalar-valued function f
   !> w.r.t. a scalar-valued variable x using complex step differentiation.
   impure function complex_step_derivative_T0_T0(f, x, h) result(dfdx)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: h
      real(rk)             :: dfdx

      interface
         impure function f(z) result(fz)
            use kinds
            complex(rk), intent(in) :: z
            complex(rk)             :: fz
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
      real(rk), dimension(:), intent(in) :: x
      real(rk),               intent(in) :: h
      real(rk), dimension(size(x))       :: dfdx
      real(rk), dimension(size(x))       :: temp_x
      integer                            :: i

      interface
         impure function f(z) result(fz)
            use kinds
            complex(rk), dimension(:), intent(in) :: z
            complex(rk)                           :: fz
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
      real(rk), dimension(:), intent(in)    :: x
      real(rk),               intent(in)    :: h
      real(rk), dimension(:,:), allocatable :: dfdx
      real(rk), dimension(size(x))          :: temp_x
      integer                               :: i

      interface
         impure function f(z) result(fz)
            use kinds
            complex(rk), dimension(:), intent(in)  :: z
            complex(rk), dimension(:), allocatable :: fz
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
   impure function finite_difference_T0_T0(f,x,h,method) result(dfdx)
      real(rk),     intent(in) :: x
      real(rk),     intent(in) :: h
      character(*), intent(in) :: method
      real(rk)                 :: dfdx

      interface
         impure function f(z) result(fz)
            use kinds
            real(rk), intent(in) :: z
            real(rk)             :: fz
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
   impure function finite_difference_central_T0_T0(f, x, h) result(dfdx)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: h
      real(rk)             :: dfdx

      interface
         impure function f(z) result(fz)
            use kinds
            real(rk), intent(in) :: z
            real(rk)             :: fz
         end function f
      end interface

      dfdx = (f(x+h) - f(x-h)) / (2.0_rk * h)

   end function finite_difference_central_T0_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure function finite_difference_forward_T0_T0(f, x, h) result(dfdx)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: h
      real(rk)             :: dfdx

      interface
         impure function f(z) result(fz)
            use kinds
            real(rk), intent(in) :: z
            real(rk)             :: fz
         end function f
      end interface

      dfdx = (f(x+h) - f(x)) / h

   end function finite_difference_forward_T0_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure function finite_difference_backward_T0_T0(f, x, h) result(dfdx)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: h
      real(rk)             :: dfdx

      interface
         impure function f(z) result(fz)
            use kinds
            real(rk), intent(in) :: z
            real(rk)             :: fz
         end function f
      end interface

      dfdx = (f(x) - f(x-h)) / h

   end function finite_difference_backward_T0_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   impure function finite_difference_T0_T1(f,x,h,method) result(dfdx)
      real(rk), dimension(:), intent(in) :: x
      real(rk),               intent(in) :: h
      character(*),           intent(in) :: method
      real(rk), dimension(size(x))       :: dfdx

      interface
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in) :: z
            real(rk)             :: fz
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
   impure function finite_difference_central_T0_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in) :: x
      real(rk),               intent(in) :: h
      real(rk), dimension(size(x))       :: dfdx
      real(rk), dimension(size(x))       :: temp_x
      integer                            :: i

      interface
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in) :: z
            real(rk)                           :: fz
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
   impure function finite_difference_forward_T0_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in) :: x
      real(rk),               intent(in) :: h
      real(rk), dimension(size(x))       :: dfdx
      real(rk), dimension(size(x))       :: temp_x
      integer                            :: i

      interface
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in) :: z
            real(rk)                           :: fz
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
   impure function finite_difference_backward_T0_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in) :: x
      real(rk),               intent(in) :: h
      real(rk), dimension(size(x))       :: dfdx
      real(rk), dimension(size(x))       :: temp_x
      integer                            :: i

      interface
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in) :: z
            real(rk)                           :: fz
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
   impure function finite_difference_T1_T1(f,x,h,method) result(dfdx)
      real(rk), dimension(:), intent(in)    :: x
      real(rk),               intent(in)    :: h
      character(*),           intent(in)    :: method
      real(rk), dimension(:,:), allocatable :: dfdx

      interface
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in)  :: z
            real(rk), dimension(:), allocatable :: fz
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
   impure function finite_difference_central_T1_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in)    :: x
      real(rk),               intent(in)    :: h
      real(rk), dimension(:,:), allocatable :: dfdx
      real(rk), dimension(size(x))          :: temp_x
      integer                               :: i

      interface
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in)  :: z
            real(rk), dimension(:), allocatable :: fz
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
   impure function finite_difference_forward_T1_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in)    :: x
      real(rk),               intent(in)    :: h
      real(rk), dimension(:,:), allocatable :: dfdx
      real(rk), dimension(size(x))          :: temp_x
      integer                               :: i

      interface
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in)  :: z
            real(rk), dimension(:), allocatable :: fz
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
   impure function finite_difference_backward_T1_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in)    :: x
      real(rk),               intent(in)    :: h
      real(rk), dimension(:,:), allocatable :: dfdx
      real(rk), dimension(size(x))          :: temp_x
      integer                               :: i

      interface
         impure function f(z) result(fz)
            use kinds
            real(rk), dimension(:), intent(in)  :: z
            real(rk), dimension(:), allocatable :: fz
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