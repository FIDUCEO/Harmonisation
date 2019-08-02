!> @file mod_time.F90
!! Time and date functions.
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

!> @brief Module providing time and date functions.
!! @author Ralf Quast
!! @copyright GNU Public License.
module mod_time
  use mod_base, only: dp, sp

  implicit none

  private

  public :: get_iso_date_and_time
  public :: iso_date_and_time

contains

  !> @brief Returns the system date and time in ISO format.
  !!
  !! @return the system date and time in ISO format.
  function get_iso_date_and_time() result (time)
    character(len=8) :: d
    character(len=10) :: t
    character(len=5) :: z
    character(len=29) :: time

    call date_and_time( d, t, z )
    time = d(1:4)//'-'//d(5:6)//'-'//d(7:8)//'T'//t(1:2)//':'//t(3:4)//':'//t(5:10)//z(1:3)//':'//z(4:5)
  end function get_iso_date_and_time

  !> @brief Expresses a Unix time stamp as date and time in ISO 8601 compliant format.
  !!
  !! @param[in] uts The unix time stamp
  !! @return the date and time in ISO 8601 compliant format.
  !!
  !! @remark Further reading:
  !! Richards, E. G. (2013). Calendars. Explanatory Supplement to the Astronomical Almanac, pp. 585-624
  function iso_date_and_time( uts ) result (iso)
    real(kind=dp), intent(in) :: uts
    character(len=16) :: iso
    real(kind=dp) :: jd
    integer :: day, month, year, date, time, hour, minute, secs
    integer :: jdn
    integer :: e
    integer :: f
    integer :: g
    integer :: h
    !! The Unix epoch expressed as Julian date
    real(kind=dp), parameter  :: EPOCH_19700101T000000Z = 2440587.5_dp
    !! Constants from Richards (2013)
    integer, parameter :: B = 274277
    integer, parameter :: J = 1401
    integer, parameter :: C = -38
    integer, parameter :: M = 2
    integer, parameter :: N = 12
    integer, parameter :: P = 1461
    integer, parameter :: R = 4
    integer, parameter :: S = 153
    integer, parameter :: U = 5
    integer, parameter :: V = 3
    integer, parameter :: W = 2
    integer, parameter :: Y = 4716

    !! Compute Julian date and Julian day number
    jd = EPOCH_19700101T000000Z + uts / 86400.0_dp
    jdn = int( floor( jd + 0.5_dp ) )

    !! Compute years, months, days (Richards 2013)
    f = jdn + J + (((4 * jdn + B) / 146097) * 3) / 4 + C
    e = R * f + V
    g = mod( e, P ) / R
    h = U * g + W
    day = mod( h, S ) / U + 1
    month = mod( h / S + M, N ) + 1
    year = (e / P) - Y + (N + M - month) / N

    !! Compute hours, minutes, and seconds
    secs = mod( int( floor( uts ) ), 86400 )
    if (secs < 0) then
      secs = secs + 86400
    end if
    hour = secs / 3600
    minute = (secs - hour * 3600) / 60
    secs = secs - hour * 3600 - minute * 60

    !! Compute date and time
    date = year * 10000 + month * 100 + day
    time = hour * 10000 + minute * 100 + secs

    !! Write ISO 8601 compliant time format
    write( iso, '(i8.8,a1,i6.6,a1)' ) date, 'T', time, 'Z'
  end function iso_date_and_time

end module mod_time
