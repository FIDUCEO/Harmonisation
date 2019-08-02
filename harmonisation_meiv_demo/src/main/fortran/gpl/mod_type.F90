!> @file mod_type.F90
!! Constants and derived types.
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

!> @brief Module providing constants and derived types.
!! @author Ralf Quast
!! @copyright GNU Public License.
module mod_type
  use mod_base

  implicit none

  !> @brief The maximum length of a description text.
  integer, public, parameter :: TYPE_MAX_LEN_DESC = 512

  !> @brief The maximum length of a name.
  integer, public, parameter :: TYPE_MAX_LEN_NAME = 80

  !> @brief The maximum length of a path.
  integer, public, parameter :: TYPE_MAX_LEN_PATH = 512

  !> @brief The maximum length of a tag.
  integer, public, parameter :: TYPE_MAX_LEN_TAGS = 20

  !> @brief The type to represent a matchup dataset.
  type, public :: matchup_t
    !> @brief The name of the matchup dataset.
    character(len=TYPE_MAX_LEN_NAME) :: name = ""
    !> @brief The UUID of the matchup dataset.
    character(len=36) :: uuid = ""
    !> @brief The index number of sensor i.
    integer :: i = 0
    !> @brief The index number of sensor j.
    integer :: j = 0
    !> @brief The number of sensor telemetry data records.
    integer :: m = 0
    !> @brief The number of sensor telemetry data columns for sensor i.
    integer :: ri = 0
    !> @brief The number of sensor telemetry data columns for sensor j.
    integer :: rj = 0
    !> @brief The state observable measurements for sensor i.
    real(kind=cp), allocatable :: qi(:,:)
    !> @brief The state observable measurements for sensor j.
    real(kind=cp), allocatable :: qj(:,:)
    !> @brief The measurand adjustments for sensor i.
    real(kind=cp), allocatable :: hi(:)
    !> @brief The measurand adjustments for sensor j.
    real(kind=cp), allocatable :: hj(:)
    !> @brief The measurand adjustment error variance for sensor i.
    real(kind=cp), allocatable :: vhi(:)
    !> @brief The measurand adjustment error variance for sensor j.
    real(kind=cp), allocatable :: vhj(:)
    !> @brief The times of measurement for sensor i.
    real(kind=tp), allocatable :: ti(:)
    !> @brief The times of measurement for sensor j.
    real(kind=tp), allocatable :: tj(:)
    !> @brief The pixel positions for sensor i.
    integer, allocatable :: pi(:)
    !> @brief The pixel positions for sensor j.
    integer, allocatable :: pj(:)
    !> @brief The uncertainty types for sensor i.
    integer, allocatable :: ui(:)
    !> @brief The uncertainty types for sensor j.
    integer, allocatable :: uj(:)
    !> @brief The C matrix indices for sensor i.
    integer, allocatable :: ic(:)
    !> @brief The C matrix indices for sensor j.
    integer, allocatable :: jc(:)
    !> @brief The V matrix indices for sensor i.
    integer, allocatable :: iv(:)
    !> @brief The V matrix indices for sensor j.
    integer, allocatable :: jv(:)
    !> @brief The W matrix indices for sensor i.
    integer, allocatable :: iw(:)
    !> @brief The W matrix indices for sensor j.
    integer, allocatable :: jw(:)
    !> @brief The Z matrix index for sensor i.
    integer :: iz
    !> @brief The Z matrix index for sensor j.
    integer :: jz
  end type matchup_t

  !> @brief The type to represent matchup dataset statistics.
  type, public :: matchup_statistics_t
    !> @brief The number of matches.
    integer :: m
    !> @brief The cost.
    real(kind=wp) :: cost = 0.0_wp
    !> @brief The average cost per match.
    real(kind=wp) :: cost_reduced = 0.0_wp
    !> @brief The p-value of the cost.
    real(kind=wp) :: cost_p_value = 0.0_wp
    !> @brief The K-residual term of the cost function.
    real(kind=wp) :: cost_matchup = 0.0_wp
    !> @brief The sensor calibration term of the cost function.
    real(kind=wp) :: cost_sensors = 0.0_wp
    !> @brief The residual mean.
    real(kind=wp) :: residual_mean = 0.0_wp
    !> @brief The residual standard deviation.
    real(kind=wp) :: residual_sdev = 0.0_wp
  end type matchup_statistics_t

  !> @brief The type to represent a sensor.
  type, public :: sensor_t
    !> @brief The name of the sensor.
    character(len=TYPE_MAX_LEN_NAME) :: name = ""
    !> @brief The ID of the measurement equation
    integer :: id = 0
    !> @brief The number of calibration parameters used by the measurement equation.
    integer :: n = 0
    !> @brief The number of state observables used by the measurement equation.
    integer :: r = 0
    !> @brief The prior values of the calibration coefficients.
    real(kind=wp), allocatable :: prior_x(:)
    !> @brief The prior covariance matrix of the calibration coefficients.
    real(kind=wp), allocatable :: prior_c(:,:)
    !> @brief The inverse prior covariance matrix of the calibration coefficients.
    real(kind=wp), allocatable :: prior_c_inverse(:,:)
    !> @brief The number of configuration constants used by the measurement equation.
    integer :: o = 0
    !> @brief The configuration constants.
    real(kind=wp), allocatable :: y(:)
  end type sensor_t

  !> @brief The type to represent the calibration parameter dataset.
  type, public :: parameter_dataset_t
    !> @brief The number of calibration parameters.
    integer :: n = 0
    !> @brief The number of matchup datasets.
    integer :: nmatchup = 0
    !> @brief The number of sensors.
    integer :: nsensor = 0
    !> @brief The number of measurement equation configuration constants.
    integer :: nconfig = 0
    !> @brief The calibration parameters.
    real(kind=wp), allocatable :: x(:)
    !> @brief The calibration parameter uncertainties.
    real(kind=wp), allocatable :: xu(:)
    !> @brief The calibration parameter error covariance matrix.
    real(kind=wp), allocatable :: xc(:,:)
    !> @brief The calibration parameter error correlation matrix.
    real(kind=wp), allocatable :: xr(:,:)
    !> @brief The calibration parameter add-offsets used internally.
    real(kind=wp), allocatable :: x_add_offsets(:)
    !> @brief The calibration parameter scale factors used internally.
    real(kind=wp), allocatable :: x_scale_factors(:)
    !> @brief The names of the calibration parameters.
    character(len=TYPE_MAX_LEN_NAME), allocatable :: x_names(:)
    !> @brief The measurement equation identifiers.
    integer, allocatable :: sensor_eqn_ids(:)
    !> @brief The number of measurement equation parameters for each sensor.
    integer, allocatable :: sensor_eqn_n(:)
    !> @brief The number of measurement equation configuration constants for each sensor.
    integer, allocatable :: sensor_eqn_o(:)
    !> @brief The measurement equation configuration constants.
    real(kind=wp), allocatable :: sensor_eqn_y(:)
    !> @brief The sensor names.
    character(len=TYPE_MAX_LEN_NAME), allocatable :: sensor_names(:)
    !> @brief The matchup statistics.
    type(matchup_statistics_t), allocatable :: matchup_statistics(:)
    !> @brief The names of the sensors included with a matchup dataset.
    character(len=TYPE_MAX_LEN_NAME), allocatable :: matchup_sensor_names(:)
    !> @brief The begin date.
    character(len=8) :: begin_date = ""
    !> @brief The end date.
    character(len=8) :: end_date = ""
    !> @brief The total matchup statistics.
    type(matchup_statistics_t) :: statistics
    !> @brief The UUIDs of the source matchup datasets.
    character(len=36), allocatable :: matchup_dataset_uuids(:)
  end type parameter_dataset_t

  !> @brief The type to represent the calibration parameter dataset.
  type, public :: residuals_dataset_t
    !> @brief The number of matches.
    integer :: m = 0
    !> @brief The times of measurement for the first sensor.
    real(kind=tp), allocatable :: t(:)
    !> @brief The K residuals.
    real(kind=wp), allocatable :: k_res(:)
    !> @brief The uncertainty of K the residuals.
    real(kind=wp), allocatable :: k_res_unc(:)
    !> @brief The weights of the K residuals.
    real(kind=wp), allocatable :: k_res_weights(:)
    !> @brief The measurand for the 1st sensor.
    real(kind=wp), allocatable :: meas_i(:)
    !> @brief The measurand for the 2nd sensor.
    real(kind=wp), allocatable :: meas_j(:)
    !> @brief The measurand uncertainty for the 1st sensor (due uncertainty of measurement equation).
    real(kind=wp), allocatable :: meas_i_unc(:)
    !> @brief The measurand uncertainty for the 2nd sensor (due uncertainty of measurement equation).
    real(kind=wp), allocatable :: meas_j_unc(:)
    !> @brief The measurand uncertainty for the 1st sensor (due uncertainty of sensor telemetry).
    real(kind=wp), allocatable :: meas_i_unc_q(:)
    !> @brief The measurand uncertainty for the 2nd sensor (due uncertainty of sensor telemetry).
    real(kind=wp), allocatable :: meas_j_unc_q(:)
    !> @brief The measurand uncertainty for the 1st sensor (due to uncertainty of sensor calibration).
    real(kind=wp), allocatable :: meas_i_unc_x(:)
    !> @brief The measurand uncertainty for the 2nd sensor (due to uncertainty of sensor calibration).
    real(kind=wp), allocatable :: meas_j_unc_x(:)
    !> @brief The begin date.
    character(len=8) :: begin_date = ""
    !> @brief The end date.
    character(len=8) :: end_date = ""
    !> @brief The name of the 1st sensor.
    character(len=TYPE_MAX_LEN_NAME) :: sensor_i_name = ""
    !> @brief The name of the 2nd sensor.
    character(len=TYPE_MAX_LEN_NAME) :: sensor_j_name = ""
    !> @brief The UUID of the source matchup dataset.
    character(len=36) :: matchup_dataset_uuid = ""
  end type residuals_dataset_t

end module mod_type
