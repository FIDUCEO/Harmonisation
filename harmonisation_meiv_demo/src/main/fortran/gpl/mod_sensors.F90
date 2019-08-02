!> @file mod_sensors.F90
!! Sensor data structure.
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

!> @brief Module providing the sensor data structure.
!! @author Ralf Quast
!! @copyright GNU Public License.
module mod_sensors
  use mod_type

  implicit none

  private

  public :: sensors_initialize
  public :: sensors_finalize
  public :: add_sensor
  public :: get_sensor_capacity
  public :: get_sensor_count
  public :: get_sensor_equation_id
  public :: get_sensor_name
  public :: get_sensor_n
  public :: get_sensor_r
  public :: get_sensor_prior_c
  public :: set_sensor_prior_c
  public :: get_sensor_prior_x
  public :: set_sensor_prior_x
  public :: get_start_index_x
  public :: get_final_index_x

  !> @brief The maximum number of sensors.
  integer, public, protected :: capacity = 0
  !> @brief The actual number of sensors.
  integer, public, protected :: nsensors = 0
  !> @brief The sensors.
  type(sensor_t), allocatable, public, protected :: sensors(:)
  !> @brief The start indices into the array of sensor calibration parameters.
  !! The ith entry points to the first calibration parameter of the ith sensor.
  integer, allocatable, public, protected :: istart(:)
  !> @brief The final indices into the array of sensor calibration parameters.
  !! The ith entry points to the final calibration parameter of the ith sensor.
  integer, allocatable, public, protected :: ifinal(:)

contains

  !> @brief Initializes the module.
  !! @param[in] capacity_in The maximum number of sensors that can be added.
  subroutine sensors_initialize( capacity_in )
    integer, intent(in) :: capacity_in

    call sensors_finalize
    capacity = capacity_in
    nsensors = 0
    allocate (sensors(capacity))
    allocate (istart(capacity))
    allocate (ifinal(capacity))
  end subroutine sensors_initialize

  !> @brief Finalizes the module and deallocates all memory used.
  subroutine sensors_finalize
    integer :: i

    do i = nsensors, 1, -1
      call delete_sensor( sensors(i) )
    end do
    if (allocated( sensors )) then
      deallocate (sensors)
    end if
    if (allocated( istart )) then
      deallocate (istart)
    end if
    if (allocated( ifinal )) then
      deallocate (ifinal)
    end if
    capacity = 0
    nsensors = 0
  end subroutine sensors_finalize

  subroutine delete_sensor( sensor )
    type(sensor_t), intent(inout) :: sensor

    sensor%name = ""
    sensor%id = 0
    sensor%n = 0
    sensor%r = 0
    sensor%o = 0
    if (allocated( sensor%prior_x )) then
      sensor%prior_x = 0.0_wp
      deallocate (sensor%prior_x)
    end if
    if (allocated( sensor%prior_c )) then
      sensor%prior_c = 0.0_wp
      deallocate (sensor%prior_c)
    end if
    if (allocated( sensor%prior_c_inverse )) then
      sensor%prior_c_inverse = 0.0_wp
      deallocate (sensor%prior_c_inverse)
    end if
    if (allocated( sensor%y )) then
      sensor%y = 0.0_wp
      deallocate (sensor%y)
    end if
  end subroutine delete_sensor

  !> @brief Adds a new sensor to the module.
  !! @param[in] name The sensor name.
  !! @param[in] id The sensor name.
  !! @param[in] n The number of sensor calibration parameters.
  !! @param[in] r The number of sensor telemetry data.
  !! @param[in] o The number of sensor configuration parameters.
  !! @param[in] prior_x The a priori expected values of the calibration parameters.
  !! @param[in] prior_u The a priori expected uncertainty of the calibration parameters.
  !! @param[out] sensor_index The sensor index.
  subroutine add_sensor( name, id, n, r, o, prior_x, prior_u, sensor_index )
    use mod_config, only: config_job_sensor_configurations

    character(len=*), intent(in) :: name
    integer, intent(in) :: id
    integer, intent(in) :: n
    integer, intent(in) :: r
    integer, intent(in) :: o
    real(kind=wp), intent(in), optional :: prior_x(n)
    real(kind=wp), intent(in), optional :: prior_u(n)
    integer, intent(out) :: sensor_index

    associate (i => nsensors)
      if (i < capacity) then
        i = i + 1
        sensors(i)%name = name
        sensors(i)%id = id
        sensors(i)%n = n
        sensors(i)%r = r

        if (i > 1) then
          istart(i) = ifinal(i - 1) + 1
          ifinal(i) = ifinal(i - 1) + n
        else
          istart(i) = 1
          ifinal(i) = n
        end if

        if (o == 0) then !! the sensor equation does not have configuration parameters, set add offsets and scale factors to neutral values
          sensors(i)%o = n + n
          allocate (sensors(i)%y(n + n))
          sensors(i)%y(1:n + n:2) = 0.0_wp !! the add offsets
          sensors(i)%y(2:n + n:2) = 1.0_wp !! the scale factors
        else
          sensors(i)%o = o
          allocate (sensors(i)%y(o))
          sensors(i)%y = config_job_sensor_configurations(i,1:o)
        end if

        allocate (sensors(i)%prior_x(n))
        allocate (sensors(i)%prior_c(n,n))
        allocate (sensors(i)%prior_c_inverse(n,n))

        sensors(i)%prior_x = 0.0_wp
        sensors(i)%prior_c = 0.0_wp
        sensors(i)%prior_c_inverse = 0.0_wp

        if (present( prior_x )) then
          call set_sensor_prior_x( i, n, prior_x )
        end if
        if (present( prior_u )) then
          call set_sensor_prior_c( i, n, prior_u )
        end if

        sensor_index = i
      else
        sensor_index = 0
      end if
    end associate
  end subroutine add_sensor

  pure integer function get_sensor_capacity()
    get_sensor_capacity = capacity
  end function get_sensor_capacity

  pure integer function get_sensor_count()
    get_sensor_count = nsensors
  end function get_sensor_count

  integer function get_sensor_equation_id( i )
    integer, intent(in) :: i

    call assert_sensor( i )
    get_sensor_equation_id = sensors(i)%id
  end function get_sensor_equation_id

  character(len=TYPE_MAX_LEN_NAME) function get_sensor_name( i )
    integer, intent(in) :: i

    call assert_sensor( i )
    get_sensor_name = sensors(i)%name
  end function get_sensor_name

  integer function get_sensor_n( i )
    integer, intent(in) :: i

    call assert_sensor( i )
    get_sensor_n = sensors(i)%n
  end function get_sensor_n

  !> Returns the number of sensor telemetry data columns.
  !! @param[in] i The sensor index.
  !! @return the number of sensor telemetry data columns for sensor @c i.
  integer function get_sensor_r( i )
    integer, intent(in) :: i

    call assert_sensor( i )
    get_sensor_r = sensors(i)%r
  end function get_sensor_r

  subroutine get_sensor_prior_c( i, n, c )
    integer, intent(in) :: i
    integer, intent(in) :: n
    real(kind=wp), intent(out) :: c(n,n)

    call assert_sensor( i )
    c = sensors(i)%prior_c
  end subroutine get_sensor_prior_c

  subroutine set_sensor_prior_c( i, n, u )
    integer,       intent(in) :: i
    integer,       intent(in) :: n
    real(kind=wp), intent(in) :: u(n)
    integer :: j

    call assert_sensor( i )

    if (n > 0) then
      sensors(i)%prior_c = 0.0_wp
      sensors(i)%prior_c_inverse = 0.0_wp

      do concurrent (j = 1:n)
        if (u(j) > 0.0_wp) then
          sensors(i)%prior_c(j,j) = u(j)**2
          sensors(i)%prior_c_inverse(j,j) = 1.0_wp / u(j)**2
        end if
      end do
    end if
  end subroutine set_sensor_prior_c

  subroutine get_sensor_prior_x( i, n, x )
    integer, intent(in) :: i
    integer, intent(in) :: n
    real(kind=wp), intent(out) :: x(n)

    call assert_sensor( i )
    x = sensors(i)%prior_x
  end subroutine get_sensor_prior_x

  subroutine set_sensor_prior_x( i, n, x )
    integer,       intent(in) :: i
    integer,       intent(in) :: n
    real(kind=wp), intent(in) :: x(n)

    call assert_sensor( i )
    sensors(i)%prior_x = x
  end subroutine set_sensor_prior_x

  pure integer function get_start_index_x( i )
    integer, intent(in) :: i

    get_start_index_x = istart(i)
  end function get_start_index_x

  pure integer function get_final_index_x( i )
    integer, intent(in) :: i

    get_final_index_x = ifinal(i)
  end function get_final_index_x

  subroutine assert_sensor( i )
    integer, intent(in) :: i

    if (i < 1 .or. i > nsensors) then
      !! @todo handle error
      error stop "Invalid sensor index"
    end if
  end subroutine assert_sensor

end module mod_sensors
