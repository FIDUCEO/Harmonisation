!> @file mod_info.F90
!! Software build information.
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

!> @brief Module to define software build information.
!! @author Ralf Quast
!! @copyright GNU Public License.
module mod_info
  implicit none

  private

  !> @brief The project name.
  character(len=*), public, parameter :: PROJECT_NAME = "@PROJECT_NAME@"

  !> @brief The project version.
  character(len=*), public, parameter :: PROJECT_VERSION = "@PROJECT_VERSION@"

  !> @brief The project tag.
  character(len=*), public, parameter :: PROJECT_TAG = "@PROJECT_TAG@"

  !> @brief The project URL.
  character(len=*), public, parameter :: PROJECT_URL = "@PROJECT_URL@"

  !> @brief The project name and version identifier.
  character(len=*), public, parameter :: PROJECT_LONG_NAME = "@PROJECT_NAME@-@PROJECT_VERSION@ @PROJECT_TAG@"

  !> @brief The vendor of the Fortran compiler used to compile the project.
  character(len=*), public, parameter :: FORTRAN_COMPILER = "@CMAKE_Fortran_COMPILER_ID@"

  !> @brief The version of the Fortran compiler used to compile the project.
  character(len=*), public, parameter :: FORTRAN_COMPILER_VERSION = "@CMAKE_Fortran_COMPILER_VERSION@"

  !> @brief The composite name of the operating system the project is compiled for.
  character(len=*), public, parameter :: SYSTEM = "@CMAKE_SYSTEM@"

  !> The name of the CPU the project is build for.
  character(len=*), public, parameter :: SYSTEM_PROCESSOR = "@CMAKE_SYSTEM_PROCESSOR@"
    
end module mod_info
