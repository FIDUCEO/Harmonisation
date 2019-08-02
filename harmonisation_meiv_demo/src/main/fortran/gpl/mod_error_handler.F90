!> @file mod_error_handler.F90
!! Runtime error handling.
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

!> @brief Module to handle runtime errors.
!! @author Ralf Quast
!! @copyright GNU Public License.
module mod_error_handler
  use, intrinsic :: iso_fortran_env, only: output_unit, error_unit

  implicit none

  private

  public handle_error

  interface handle_error
    module procedure handle_error_impl
  end interface

  integer :: out = error_unit

contains

  subroutine handle_error_impl( status, message, what )
    integer, intent(in), optional :: status
    character(len=*), intent(in), optional :: message
    character(len=*), intent(in), optional :: what

    if (present( status )) then
      call error_message_name_value_int( "status", status )
    end if

    if (present( message )) then
      call error_message_text( message )
    end if

    if (present( what )) then
      call error_message_name_value_text( "what", what )
    end if

    error stop "An error has occurred."
  end subroutine handle_error_impl

  subroutine error_message_name_value_int( name, value )
    character(len=*), intent(in) :: name
    integer, intent(in) :: value

    write (out, '("Error: ",a," = ",i0)') name, value
  end subroutine error_message_name_value_int

  subroutine error_message_name_value_text( name, value )
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: value

    write (out, '("Error: ",a," = ",a)') name, value
  end subroutine error_message_name_value_text

  subroutine error_message_text( text )
    character(len=*), intent(in) :: text

    write (out, '("Error: ",a)') text
  end subroutine error_message_text

end module mod_error_handler
