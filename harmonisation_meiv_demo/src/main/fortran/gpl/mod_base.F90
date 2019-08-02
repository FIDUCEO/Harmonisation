!> @file mod_base.F90
!! Base kind parameters.
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

!> @brief Module to define base kind parameters.
!! @authors Ralf Giering, Ralf Quast
!! @copyright GNU Public License.
module mod_base
  implicit none

  private

  !> @brief The byte integer kind (one byte).
  integer, public, parameter :: byte = selected_int_kind(2)

  !> @brief The long integer kind (eight bytes).
  integer, public, parameter :: long = selected_int_kind(18)

  !> @brief At least IEEE S_floating ("single precision").
  integer, public, parameter :: sp = selected_real_kind(p=6,r=37)

  !> @brief At least IEEE T_floating ("double precision").
  integer, public, parameter :: dp = selected_real_kind(p=15,r=307)

  !> @brief The working precision -- use this to set/switch precision globally.
  integer, public, parameter :: wp = dp

  !> @brief The precision of coordinates or measurements.
  integer, public, parameter :: cp = sp

  !> @brief The precision of time data.
  integer, public, parameter :: tp = dp

end module mod_base
