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
   end interface
   !===============================================================================

contains

   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a scalar-valued function f
   !> w.r.t. a scalar-valued variable x using complex step differentiation.
   function complex_step_derivative_T0_T0(f, x, h) result(dfdx)
      real(rk), intent(in) :: x
      real(rk), intent(in) :: h
      real(rk)             :: dfdx

      interface
         function f(z) result(fz)
            use kinds
            complex(rk), intent(in) :: z
            complex(rk)             :: fz
         end function f
      end interface

      dfdx = aimag(f(cmplx(x, h, rk))) / h

   end function complex_step_derivative_T0_T0
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the derivative of a scalar-valued function f
   !> w.r.t. a vector-valued variable x using complex step differentiation.
   function complex_step_derivative_T0_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in) :: x
      real(rk), intent(in)               :: h
      real(rk), dimension(size(x))       :: dfdx
      real(rk), dimension(size(x)) :: temp_x
      integer :: i
      interface
         function f(z) result(fz)
            use kinds
            complex(rk), dimension(:), intent(in) :: z
            complex(rk)                           :: fz
         end function f
      end interface

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
   function complex_step_derivative_T1_T1(f, x, h) result(dfdx)
      real(rk), dimension(:), intent(in)    :: x
      real(rk), intent(in)                  :: h
      real(rk), dimension(:,:), allocatable :: dfdx
      real(rk), dimension(size(x))          :: temp_x
      integer                               :: i

      interface
         function f(z) result(fz)
            use kinds
            complex(rk), dimension(:), intent(in)  :: z
            complex(rk), dimension(:), allocatable :: fz
         end function f
      end interface

      allocate(dfdx(size(f(cmplx(x,kind=rk))),size(x)))

      do i = 1, size(x)
         temp_x    = 0.0_rk
         temp_x(i) = x(i)
         dfdx(:,i) = aimag(f(cmplx(temp_x, h, rk))) / h
      end do

   end function complex_step_derivative_T1_T1
   !===============================================================================

end module fordiff