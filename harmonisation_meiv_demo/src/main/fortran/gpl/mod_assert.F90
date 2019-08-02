!> @file mod_assert.F90
!! Assertion of data validity.
!! Copyright (C) 2019 FastOpt GmbH, Hamburg, Germany (info@fastopt.de)
!!
!! This code was developed for the EC project "Fidelity and Uncertainty in
!! Climate Data Records from Earth Observations (FIDUCEO)".
!! Grant Agreement: 638822
!!
!! This program is free software; you can redistribute it and/or modify it
!! under the terms of the GNU General Public License as published by the Free
!! Software Foundation; either version 3 of the License, or (at your option)
!! any later version.
!! This program is distributed in the hope that it will be useful, but WITHOUT
!! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
!! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
!! more details.
!!
!! A copy of the GNU General Public License should have been supplied along
!! with this program; if not, see http://www.gnu.org/licenses/

!> @brief Module to assert the validity of data.
!! @author Ralf Quast
!! @copyright GNU Public License.
module mod_assert
  use, intrinsic :: ieee_arithmetic
  use mod_base, only: dp, sp, long
  use mod_logger

  implicit none

  private

  public assert_equals
  public assert_true
  public assert_not_huge
  public assert_not_infinite
  public assert_not_nan
  public assert_not_zero
  public assert_positive

  interface assert_equals
    module procedure assert_equals_0__integer
  end interface

  interface assert_true
    module procedure assert_true_0
  end interface

  interface assert_not_infinite
    module procedure assert_not_infinite_1__dp
    module procedure assert_not_infinite_1__sp
    module procedure assert_not_infinite_2__dp
    module procedure assert_not_infinite_2__sp
  end interface

  interface assert_not_nan
    module procedure assert_not_nan_1__dp
    module procedure assert_not_nan_1__sp
    module procedure assert_not_nan_2__dp
    module procedure assert_not_nan_2__sp
  end interface

  interface assert_not_huge
    module procedure assert_not_huge_1__dp
    module procedure assert_not_huge_1__sp
    module procedure assert_not_huge_2__dp
    module procedure assert_not_huge_2__sp
  end interface

  interface assert_not_zero
    module procedure assert_not_zero_1__dp
    module procedure assert_not_zero_1__sp
    module procedure assert_not_zero_2__dp
    module procedure assert_not_zero_2__sp
  end interface

  interface assert_positive
    module procedure assert_positive_1__dp
    module procedure assert_positive_1__sp
    module procedure assert_positive_2__dp
    module procedure assert_positive_2__sp
  end interface

  integer, parameter :: VALUE_IS_HUGE = 10
  integer, parameter :: VALUE_IS_INIFINITE = 11
  integer, parameter :: VALUE_IS_NAN = 12
  integer, parameter :: VALUE_IS_INCORRECT = 13
  integer, parameter :: VALUE_IS_ZERO = 14
  integer, parameter :: VALUE_IS_NEGATIVE_OR_ZERO = 14

contains

  logical function all_1( array )
    logical, intent(in) :: array(:)
    integer(kind=long) :: i

    if (get_log_level() < LOG_LEVEL_DEBUG) then
      all_1 = all( array )
    else
      all_1 = .true.
      do i = 1, size( array, kind=long )
        if (.not. array(i)) then
          all_1 = .false.
          call log_debug( 'assert', 'array(i) element i = ', i )
        end if
      end do
    end if
  end function all_1

  logical function all_2( array )
    logical, intent(in) :: array(:,:)
    integer :: i, j

    if (get_log_level() < LOG_LEVEL_DEBUG) then
      all_2 = all( array )
    else
      all_2 = .true.
      do i = 1, size( array, 1 )
        do j = 1, size( array, 2 )
          if (.not. array(i,j)) then
            all_2 = .false.
            call log_debug( 'assert', 'array(i,j) element i = ', i )
            call log_debug( 'assert', 'array(i,j) element j = ', j )
          end if
        end do
      end do
    end if
  end function all_2

  logical function any_1( array )
    logical, intent(in) :: array(:)
    integer(kind=long) :: i

    if (get_log_level() < LOG_LEVEL_DEBUG) then
      any_1 = any( array )
    else
      any_1 = .false.
      do i = 1, size( array, kind=long )
        if (array(i)) then
          any_1 = .true.
          call log_debug( 'assert', "array(i) element i = ", i )
        end if
      end do
    end if
  end function any_1

  logical function any_2( array )
    logical, intent(in) :: array(:,:)
    integer :: i, j

    if (get_log_level() < LOG_LEVEL_DEBUG) then
      any_2 = any( array )
    else
      any_2 = .false.
      do i = 1, size( array, 1 )
        do j = 1, size( array, 2 )
          if (array(i,j)) then
            any_2 = .true.
            call log_debug( "assert", "array(i,j) element i = ", i )
            call log_debug( "assert", "array(i,j) element j = ", j )
          end if
        end do
      end do
    end if
  end function any_2

  subroutine assert_equals_0__integer( actual, expected, name )
    integer, intent(in) :: actual
    integer, intent(in) :: expected
    character(len=*), intent(in) :: name

    if (.not. (actual == expected)) then
      call log_error( "assert", name//" != ", expected )
      error stop VALUE_IS_INCORRECT
    end if
  end subroutine assert_equals_0__integer

  subroutine assert_true_0( actual, name )
    logical, intent(in) :: actual
    character(len=*), intent(in) :: name

    if (.not. actual) then
      call log_error( "assert", name//" != .true." )
      error stop VALUE_IS_INCORRECT
    end if
  end subroutine assert_true_0

  subroutine assert_not_infinite_1__dp( array, array_name )
    real(kind=dp), intent(in) :: array(:)
    character(len=*), intent(in) :: array_name

    if (.not. all_1( ieee_is_finite( array ) )) then
      call log_error( "assert", "some value in array '"//array_name//"' is infinite" )
      error stop VALUE_IS_INIFINITE
    end if
  end subroutine assert_not_infinite_1__dp

  subroutine assert_not_infinite_1__sp( array, array_name )
    real(kind=sp), intent(in) :: array(:)
    character(len=*), intent(in) :: array_name

    if (.not. all_1( ieee_is_finite( array ) )) then
      call log_error( "assert", "some value in array '"//array_name//"' is infinite" )
      error stop VALUE_IS_INIFINITE
    end if
  end subroutine assert_not_infinite_1__sp

  subroutine assert_not_infinite_2__dp( array, array_name )
    real(kind=dp), intent(in) :: array(:,:)
    character(len=*), intent(in) :: array_name

    if (.not. all_2( ieee_is_finite( array ) )) then
      call log_error( "assert", "some value in array '"//array_name//"' is infinite" )
      error stop VALUE_IS_INIFINITE
    end if
  end subroutine assert_not_infinite_2__dp

  subroutine assert_not_infinite_2__sp( array, array_name )
    real(kind=sp), intent(in) :: array(:,:)
    character(len=*), intent(in) :: array_name

    if (.not. all_2( ieee_is_finite( array ) )) then
      call log_error( "assert", "some value in array '"//array_name//"' is infinite" )
      error stop VALUE_IS_INIFINITE
    end if
  end subroutine assert_not_infinite_2__sp

  subroutine assert_not_nan_1__dp( array, array_name )
    real(kind=dp), intent(in) :: array(:)
    character(len=*), intent(in) :: array_name

    if (any_1( ieee_is_nan( array ) )) then
      call log_error( "assert", "some value in array '"//array_name//"' is not a number (NaN)" )
      error stop VALUE_IS_NAN
    end if
  end subroutine assert_not_nan_1__dp

  subroutine assert_not_nan_1__sp( array, array_name )
    real(kind=sp), intent(in) :: array(:)
    character(len=*), intent(in) :: array_name

    if (any_1( ieee_is_nan( array ) )) then
      call log_error( "assert", "some value in array '"//array_name//"' is not a number (NaN)" )
      error stop VALUE_IS_NAN
    end if
  end subroutine assert_not_nan_1__sp

  subroutine assert_not_nan_2__dp( array, array_name )
    real(kind=dp), intent(in) :: array(:,:)
    character(len=*), intent(in) :: array_name

    if (any_2( ieee_is_nan( array ) )) then
      call log_error( "assert", "some value in array '"//array_name//"' is not a number (NaN)" )
      error stop VALUE_IS_NAN
    end if
  end subroutine assert_not_nan_2__dp

  subroutine assert_not_nan_2__sp( array, array_name )
    real(kind=sp), intent(in) :: array(:,:)
    character(len=*), intent(in) :: array_name

    if (any_2( ieee_is_nan( array ) )) then
      call log_error( "assert", "some value in array '"//array_name//"' is not a number (NaN)" )
      error stop VALUE_IS_NAN
    end if
  end subroutine assert_not_nan_2__sp

  subroutine assert_not_huge_1__dp( array, array_name )
    real(kind=dp), intent(in) :: array(:)
    character(len=*), intent(in) :: array_name

    if (any_1( abs( array ) > sqrt( huge( array ) ) )) then
      call log_error( "assert", "some value in array '"//array_name//"' is huge" )
      error stop VALUE_IS_HUGE
    end if
  end subroutine assert_not_huge_1__dp

  subroutine assert_not_huge_1__sp( array, array_name )
    real(kind=sp), intent(in) :: array(:)
    character(len=*), intent(in) :: array_name

    if (any_1( abs( array ) > sqrt( huge( array ) ) )) then
      call log_error( "assert", "some value in array '"//array_name//"' is huge" )
      error stop VALUE_IS_HUGE
    end if
  end subroutine assert_not_huge_1__sp

  subroutine assert_not_huge_2__dp( array, array_name )
    real(kind=dp), intent(in) :: array(:,:)
    character(len=*), intent(in) :: array_name

    if (any_2( abs( array ) > sqrt( huge( array ) ) )) then
      call log_error( "assert", "some value in array '"//array_name//"' is huge" )
      error stop VALUE_IS_HUGE
    end if
  end subroutine assert_not_huge_2__dp

  subroutine assert_not_huge_2__sp( array, array_name )
    real(kind=sp), intent(in) :: array(:,:)
    character(len=*), intent(in) :: array_name

    if (any_2( abs( array ) > sqrt( huge( array ) ) )) then
      call log_error( "assert", "some value in array '"//array_name//"' is huge" )
      error stop VALUE_IS_HUGE
    end if
  end subroutine assert_not_huge_2__sp

  subroutine assert_not_zero_1__dp( array, array_name )
    real(kind=dp), intent(in) :: array(:)
    character(len=*), intent(in) :: array_name

    if (any_1( array == 0.0_dp )) then
      call log_error( "assert", "some value in array '"//array_name//"' is zero" )
      error stop VALUE_IS_ZERO
    end if
  end subroutine assert_not_zero_1__dp

  subroutine assert_not_zero_1__sp( array, array_name )
    real(kind=sp), intent(in) :: array(:)
    character(len=*), intent(in) :: array_name

    if (any_1( array == 0.0_sp )) then
      call log_error( "assert", "some value in array '"//array_name//"' is zero" )
      error stop VALUE_IS_ZERO
    end if
  end subroutine assert_not_zero_1__sp

  subroutine assert_not_zero_2__dp( array, array_name )
    real(kind=dp), intent(in) :: array(:,:)
    character(len=*), intent(in) :: array_name

    if (any_2( array == 0.0_dp )) then
      call log_error( "assert", "some value in array '"//array_name//"' is zero" )
      error stop VALUE_IS_ZERO
    end if
  end subroutine assert_not_zero_2__dp

  subroutine assert_not_zero_2__sp( array, array_name )
    real(kind=sp), intent(in) :: array(:,:)
    character(len=*), intent(in) :: array_name

    if (any_2( array == 0.0_sp )) then
      call log_error( "assert", "some value in array '"//array_name//"' is zero" )
      error stop VALUE_IS_ZERO
    end if
  end subroutine assert_not_zero_2__sp

  subroutine assert_positive_1__dp( array, array_name )
    real(kind=dp), intent(in) :: array(:)
    character(len=*), intent(in) :: array_name

    if (any_1( array <= 0.0_dp )) then
      call log_error( "assert", "some value in array '"//array_name//"' is negative or zero" )
      error stop VALUE_IS_NEGATIVE_OR_ZERO
    end if
  end subroutine assert_positive_1__dp

  subroutine assert_positive_1__sp( array, array_name )
    real(kind=sp), intent(in) :: array(:)
    character(len=*), intent(in) :: array_name

    if (any_1( array <= 0.0_sp )) then
      call log_error( "assert", "some value in array '"//array_name//"' is negative or zero" )
      error stop VALUE_IS_NEGATIVE_OR_ZERO
    end if
  end subroutine assert_positive_1__sp

  subroutine assert_positive_2__dp( array, array_name )
    real(kind=dp), intent(in) :: array(:,:)
    character(len=*), intent(in) :: array_name

    if (any_2( array <= 0.0_dp )) then
      call log_error( "assert", "some value in array '"//array_name//"' is negative or zero" )
      error stop VALUE_IS_NEGATIVE_OR_ZERO
    end if
  end subroutine assert_positive_2__dp

  subroutine assert_positive_2__sp( array, array_name )
    real(kind=sp), intent(in) :: array(:,:)
    character(len=*), intent(in) :: array_name

    if (any_2( array <= 0.0_sp )) then
      call log_error( "assert", "some value in array '"//array_name//"' is negative or zero" )
      error stop VALUE_IS_NEGATIVE_OR_ZERO
    end if
  end subroutine assert_positive_2__sp

end module mod_assert
