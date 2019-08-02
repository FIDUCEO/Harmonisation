!> @file mod_logger.F90
!! Handling of log messages.
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

!> @brief Module to handle log messages.
!! @author Ralf Quast
!! @copyright GNU Public License.
module mod_logger
  use mod_base, only: dp, sp, long
  use mod_time, only: get_iso_date_and_time

  implicit none

  private

  public :: logger_initialize
  public :: get_log_level
  public :: log_error
  public :: log_warning
  public :: log_info
  public :: log_fine
  public :: log_debug

  interface log_debug
    module procedure log_debug_int
    module procedure log_debug_long
    module procedure log_debug_real__dp
    module procedure log_debug_real__sp
    module procedure log_debug_text
  end interface

  interface log_error
    module procedure log_error_int
    module procedure log_error_real__dp
    module procedure log_error_real__sp
    module procedure log_error_text
  end interface

  interface log_fine
    module procedure log_fine_int
    module procedure log_fine_real__dp
    module procedure log_fine_real__sp
    module procedure log_fine_text
  end interface

  interface log_info
    module procedure log_info_int
    module procedure log_info_real__dp
    module procedure log_info_real__sp
    module procedure log_info_text
  end interface

  interface log_warning
    module procedure log_warning_text
  end interface

  integer, public, parameter :: LOG_LEVEL_ERROR    = 1
  integer, public, parameter :: LOG_LEVEL_WARNING  = 2
  integer, public, parameter :: LOG_LEVEL_PROGRESS = 3
  integer, public, parameter :: LOG_LEVEL_INFO     = 4
  integer, public, parameter :: LOG_LEVEL_FINE     = 6
  integer, public, parameter :: LOG_LEVEL_DEBUG    = 7

  integer :: log_level = LOG_LEVEL_INFO
contains

  subroutine logger_initialize( log_level_in )
    integer, intent(in) :: log_level_in

    log_level = log_level_in
  end subroutine logger_initialize

  integer function get_log_level()
    get_log_level = log_level
  end function get_log_level

  subroutine log_debug_int( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    integer, intent(in) :: value

    if (log_level >= LOG_LEVEL_DEBUG) then
      write( *, "(a,a,a,a,a,i0)" ) get_iso_date_and_time(), " [debug] ", name, "::", message, value
    end if
  end subroutine log_debug_int

  subroutine log_debug_long( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    integer(kind=long), intent(in) :: value

    if (log_level >= LOG_LEVEL_DEBUG) then
      write( *, "(a,a,a,a,a,i0)" ) get_iso_date_and_time(), " [debug] ", name, "::", message, value
    end if
  end subroutine log_debug_long

  subroutine log_debug_real__dp( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    real(kind=dp), intent(in) :: value

    if (log_level >= LOG_LEVEL_DEBUG) then
      write( *, "(a,a,a,a,a,f0.4)" ) get_iso_date_and_time(), " [debug] ", name, "::", message, value
    end if
  end subroutine log_debug_real__dp

  subroutine log_debug_real__sp( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    real(kind=sp), intent(in) :: value

    if (log_level >= LOG_LEVEL_DEBUG) then
      write( *, "(a,a,a,a,a,f0.4)" ) get_iso_date_and_time(), " [debug] ", name, "::", message, value
    end if
  end subroutine log_debug_real__sp

  subroutine log_debug_text( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: value

    if (log_level >= LOG_LEVEL_DEBUG) then
      if (present( value )) then
        write( *, "(a,a,a,a,a,a)" ) get_iso_date_and_time(), " [debug] ", name, "::", message, trim( value )
      else
        write( *, "(a,a,a,a,a)" )   get_iso_date_and_time(), " [debug] ", name, "::", message
      end if
    end if
  end subroutine log_debug_text

  subroutine log_error_int( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    integer, intent(in) :: value

    if (log_level >= LOG_LEVEL_ERROR) then
      write( *, "(a,a,a,a,a,i0)" ) get_iso_date_and_time(), " [error] ", name, "::", message, value
    end if
  end subroutine log_error_int

  subroutine log_error_real__dp( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    real(kind=dp), intent(in) :: value

    if (log_level >= LOG_LEVEL_ERROR) then
      write( *, "(a,a,a,a,a,f0.4)" ) get_iso_date_and_time(), " [error] ", name, " :: ", message, value
    end if
  end subroutine log_error_real__dp

  subroutine log_error_real__sp( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    real(kind=sp), intent(in) :: value

    if (log_level >= LOG_LEVEL_ERROR) then
      write( *, "(a,a,a,a,a,f0.4)" ) get_iso_date_and_time(), " [error] ", name, "::", message, value
    end if
  end subroutine log_error_real__sp

  subroutine log_error_text( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: value

    if (log_level >= LOG_LEVEL_ERROR) then
      if (present( value )) then
        write( *, "(a,a,a,a,a,a)" ) get_iso_date_and_time(), " [error] ", name, "::", message, trim( value )
      else
        write( *, "(a,a,a,a,a)" )   get_iso_date_and_time(), " [error] ", name, "::", message
      end if
    end if
  end subroutine log_error_text

  subroutine log_fine_int( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    integer, intent(in) :: value

    if (log_level >= LOG_LEVEL_FINE) then
      write( *, "(a,a,a,a,a,i0)" ) get_iso_date_and_time(), " [fine] ", name, "::", message, value
    end if
  end subroutine log_fine_int

  subroutine log_fine_real__dp( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    real(kind=dp), intent(in) :: value

    if (log_level >= LOG_LEVEL_FINE) then
      write( *, "(a,a,a,a,a,f0.4)" ) get_iso_date_and_time(), " [fine] ", name, "::", message, value
    end if
  end subroutine log_fine_real__dp

  subroutine log_fine_real__sp( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    real(kind=sp), intent(in) :: value

    if (log_level >= LOG_LEVEL_FINE) then
      write( *, "(a,a,a,a,a,f0.4)" ) get_iso_date_and_time(), " [fine] ", name, "::", message, value
    end if
  end subroutine log_fine_real__sp

  subroutine log_fine_text( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: value

    if (log_level >= LOG_LEVEL_FINE) then
      if (present( value )) then
        write( *, "(a,a,a,a,a,a)" ) get_iso_date_and_time(), " [fine] ", name, "::", message, trim( value )
      else
        write( *, "(a,a,a,a,a)" )   get_iso_date_and_time(), " [fine] ", name, "::", message
      end if
    end if
  end subroutine log_fine_text

  subroutine log_info_int( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    integer, intent(in) :: value

    if (log_level >= LOG_LEVEL_INFO) then
      write( *, "(a,a,a,a,a,i0)" ) get_iso_date_and_time(), " [info] ", name, "::", message, value
    end if
  end subroutine log_info_int

  subroutine log_info_real__dp( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    real(kind=dp), intent(in) :: value

    if (log_level >= LOG_LEVEL_INFO) then
      write( *, "(a,a,a,a,a,f0.4)" ) get_iso_date_and_time(), " [info] ", name, "::", message, value
    end if
  end subroutine log_info_real__dp

  subroutine log_info_real__sp( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    real(kind=sp), intent(in) :: value

    if (log_level >= LOG_LEVEL_INFO) then
      write( *, "(a,a,a,a,a,f0.4)" ) get_iso_date_and_time(), " [info] ", name, "::", message, value
    end if
  end subroutine log_info_real__sp

  subroutine log_info_text( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: value

    if (log_level >= LOG_LEVEL_INFO) then
      if (present( value )) then
        write( *, "(a,a,a,a,a,a)" ) get_iso_date_and_time(), " [info] ", name, "::", message, trim( value )
      else
        write( *, "(a,a,a,a,a)" )   get_iso_date_and_time(), " [info] ", name, "::", message
      end if
    end if
  end subroutine log_info_text

  subroutine log_warning_text( name, message, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: message
    character(len=*), intent(in), optional :: value

    if (log_level >= LOG_LEVEL_WARNING) then
      if (present( value )) then
        write( *, "(a,a,a,a,a,a)" ) get_iso_date_and_time(), " [warning] ", name, "::", message, trim( value )
      else
        write( *, "(a,a,a,a,a)" )   get_iso_date_and_time(), " [warning] ", name, "::", message
      end if
    end if
  end subroutine log_warning_text

end module mod_logger
