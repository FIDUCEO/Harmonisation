!> @file demo.F90
!! Reader and writer demonstration.
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

!> @brief Read demonstration.
!! @param[out] n The number of model parameters.
!! @author Ralf Quast
!! @copyright GNU Public License.
subroutine read( n )
  use mod_base
  use mod_config
  use mod_matrices
  use mod_matchup
  use mod_sensors
  use mod_reader

  implicit none
  integer, intent(out) :: n
  integer :: sensor_index
  integer :: i, k

  !! The full software adds sensors as specified in the job configuration automatically.
  !! Here, we add sensors manually.
  call add_sensor( config_job_sensor_names(1),   0, 0, 1, 0, sensor_index=sensor_index )
  call add_sensor( config_job_sensor_names(2), 101, 3, 4, 7, sensor_index=sensor_index )
  call add_sensor( config_job_sensor_names(3), 101, 3, 4, 7, sensor_index=sensor_index )
  call add_sensor( config_job_sensor_names(4), 101, 3, 4, 7, sensor_index=sensor_index )

  do k = 1, config_job_matchup_count
    call read_matchup_dataset( trim( config_job_matchup_dataset_path )//trim( config_job_matchup_datasets(k) ), &
      config_job_matchup_names(k), &
      config_job_matchup(k, 1), &
      config_job_matchup(k, 2), &
      config_job_matchup_samplings(k) )

    !! The full software configures preconditioning matrices ... call matchup_ops_configure_z_matrices( i )
    !! Here, we do not do that because no optimization is performed.

    if (.not. config_job_enable_covariance_matrix) then !! remove unused matrices to save computer memory
      call remove_b_matrices
      call remove_d_matrices
      call remove_w_matrices
      call remove_v_matrices
    end if
  end do

  n = 0
  do i = 1, nsensors
    n = n + get_sensor_n( i )
  end do
end subroutine read

!> @brief Write demonstration.
!! @param[in] n The number of model parameters.
!! @param[in] x The optimised model parameters.
!! @param[in] xc The inverse of the hessian matrix of the cost function at @c x.
!! @author Ralf Quast
!! @copyright GNU Public License.
subroutine write( n, x, xc )
  use mod_base
  use mod_type
  use mod_writer

  implicit none
  integer, intent(in) :: n
  real(kind=wp), intent(in) :: x(n), xc(n,n)
  type(parameter_dataset_t) :: parameter_dataset
  type(matchup_statistics_t) :: statistics
  character(len=36), parameter :: uuid = "02b650e8-d476-48fc-881b-71ad0ea8c025"
  character(len=36), parameter :: other_uuid = "9cab83ea-4f89-45d1-91ad-d9a17f225583"
  integer :: i, j, k

  !! The full software configures a parameter dataset from the optimisation results.
  !! Here, we configure a dummy parameter dataset manually.
  parameter_dataset%n = n
  parameter_dataset%nmatchup = 3
  parameter_dataset%nsensor = 4
  parameter_dataset%x = x
  parameter_dataset%xu = [(sqrt( xc(i,i) ), i = 1, n)]
  parameter_dataset%xc = xc
  parameter_dataset%xr = reshape( [((xc(i,j) / sqrt( xc(i,i) * xc(j,j) ), i = 1, n), j = 1, n)], [n, n] )
  parameter_dataset%x_add_offsets = [(0.0_wp, i = 1, n)]
  parameter_dataset%x_scale_factors = [(1.0_wp, i = 1, n)]
  parameter_dataset%x_names = [ "1/1", "1/2", "1/3", "2/1", "2/2", "2/3", "3/1", "3/2", "3/3" ]
  parameter_dataset%sensor_eqn_ids = [ 0, 101, 101, 101 ]
  parameter_dataset%sensor_eqn_n = [ 0, 3, 3, 3 ]
  parameter_dataset%sensor_names = [ "0", "1", "2", "3" ]
  parameter_dataset%matchup_sensor_names = ["0 1", "1 2", "2 3"]
  parameter_dataset%matchup_dataset_uuids = ["", "", ""]
  statistics%m = 30000
  statistics%cost = 15000.0_wp
  statistics%cost_reduced = 0.5_wp
  statistics%cost_p_value = 1.0_wp
  statistics%cost_matchup = 15000.0_wp
  statistics%residual_mean = 0.0_wp
  statistics%residual_sdev = 1.0_wp
  parameter_dataset%statistics = statistics
  statistics%m = 10000
  statistics%cost = 5000.0_wp
  statistics%cost_matchup = 5000.0_wp
  parameter_dataset%matchup_statistics = [ statistics, statistics, statistics ]
  parameter_dataset%matchup_dataset_uuids = [ other_uuid, other_uuid, other_uuid ]

  call write_parameter_dataset( parameter_dataset, uuid )

  !! The full software also writes a diagnostic residuals dataset for each input matchup dataset.
  !! Here, we do not write residuals datasets.
end subroutine write

program demo
  use mod_base
  use mod_config
  use mod_matrices
  use mod_matchup
  use mod_sensors

  implicit none
  integer :: i, n
  real(kind=wp), allocatable :: x(:), xc(:,:)

  call config_initialize
  call matrices_initialize( CONFIG_CAPACITY_MATRICES )
  call sensors_initialize( CONFIG_CAPACITY_SENSORS )
  call matchup_initialize( CONFIG_CAPACITY_MATCHUP )

  call read( n )

  !! The full software configures and initializes the calibration model, and then optimizes the model parameters.
  !! Here, we merely set the model parameters and the associated error covariance matrix to dummy values.
  allocate (x(n))
  allocate (xc(n,n))
  x = 0.0_wp
  xc = 0.1_wp
  do i = 1, n
    xc(i,i) = 0.5_wp
  end do

  call write( n, x, xc )
end program demo

