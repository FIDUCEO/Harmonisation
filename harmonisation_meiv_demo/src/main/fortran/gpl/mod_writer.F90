!> @file mod_writer.F90
!! Writing of harmonisation result datasets.
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

!> @brief Module to write harmonisation result datasets.
!! @author Ralf Quast
!! @copyright GNU Public License.
module mod_writer
  use mod_type
  use mod_config
  use mod_logger

  implicit none

  private

  public write_parameter_dataset
  public write_residuals_dataset

  interface write_parameter_dataset
    module procedure write_parameter_dataset_nc
  end interface

  interface write_residuals_dataset
    module procedure write_residuals_dataset_nc
  end interface

  interface add_var_vector
    module procedure add_var_vector_int
    module procedure add_var_vector_real__dp
    module procedure add_var_vector_real__sp
    module procedure add_var_vector_text
  end interface

contains

  subroutine write_parameter_dataset_nc( ds, uuid )
    use mod_config
    use mod_nc

    type(parameter_dataset_t), intent(in) :: ds
    character(len=36), intent(in) :: uuid

    integer :: ncid
    character(len=TYPE_MAX_LEN_PATH) :: dsname
    integer :: n_id
    integer :: npair_id
    integer :: lname_id
    integer :: luuid_id
    integer :: nsensor_id
    integer :: nconfig_id
    integer :: x_id
    integer :: xu_id
    integer :: xc_id
    integer :: xr_id
    integer :: x_add_offset_id
    integer :: x_scale_factor_id
    integer :: x_name_id
    integer :: sensor_eqn_id
    integer :: sensor_eqn_n_id
    integer :: sensor_eqn_o_id
    integer :: sensor_eqn_y_id
    integer :: sensor_name_id
    integer :: stat_m_id
    integer :: stat_cost_id
    integer :: stat_cost_reduced_id
    integer :: stat_mean_id
    integer :: stat_sdev_id
    integer :: matchup_sensor_pairs_id
    integer :: matchup_dataset_uuids_id
    logical :: must_write_config

    must_write_config = ds%nconfig > 0

    !! 1) Create the dataset
    dsname = get_parameter_dataset_name( ".nc", ds%begin_date, ds%end_date )
    call log_info( "mod_writer", "writing parameter dataset: ", dsname )
    call create_dataset( dsname, ncid )

    !! 2) Add global attributes
    call add_att_matchup_dataset( ncid, ds%begin_date, ds%end_date )
    call add_att_statistics( ncid, ds%statistics )
    call add_att_software_and_job( ncid )
    call add_att_uuid( ncid, uuid )

    !! 3) Add dimensions
    call add_dimension( ncid, config_name_dim_l_name, TYPE_MAX_LEN_NAME, lname_id )
    call add_dimension( ncid, config_name_dim_l_uuid, 36, luuid_id )
    call add_dimension( ncid, config_name_dim_n, ds%n, n_id )
    call add_dimension( ncid, config_name_dim_p, ds%nmatchup, npair_id )
    call add_dimension( ncid, config_name_dim_s, ds%nsensor, nsensor_id )
    if (must_write_config) then
      call add_dimension( ncid, config_name_dim_c, ds%nconfig, nconfig_id )
    end if

    !! 4) Add variables
    call add_var_vector( ncid, ds%x, n_id, config_name_var_x, config_desc_var_x, x_id )
    call add_var_vector( ncid, ds%xu, n_id, config_name_var_x_unc, config_desc_var_x_unc, xu_id )
    call add_var_matrix( ncid, ds%xc, n_id, config_name_var_x_cov, config_desc_var_x_cov, xc_id )
    call add_var_matrix( ncid, ds%xr, n_id, config_name_var_x_cor, config_desc_var_x_cor, xr_id )
    call add_var_vector( ncid, ds%x_add_offsets, n_id, config_name_var_x_add_offset, &
      config_desc_var_x_add_offset, x_add_offset_id )
    call add_var_vector( ncid, ds%x_scale_factors, n_id, config_name_var_x_scale_factor, &
      config_desc_var_x_scale_factor, x_scale_factor_id )
    call add_var_vector( ncid, ds%x_names, lname_id, n_id, config_name_var_x_names, &
      config_desc_var_x_names, x_name_id )
    call add_var_vector( ncid, ds%sensor_eqn_ids, nsensor_id, config_name_var_sensor_equation_id, &
      config_desc_var_sensor_equation_id, sensor_eqn_id )
    call add_var_vector( ncid, ds%sensor_eqn_n, nsensor_id, config_name_var_sensor_equation_x_count, &
      config_desc_var_sensor_equation_x_count, sensor_eqn_n_id )
    if (must_write_config) then
      call add_var_vector( ncid, ds%sensor_eqn_o, nsensor_id, config_name_var_sensor_equation_y_count, &
        config_desc_var_sensor_equation_y_count, sensor_eqn_o_id )
      call add_var_vector( ncid, ds%sensor_eqn_y, nconfig_id, config_name_var_sensor_equation_y, &
        config_desc_var_sensor_equation_y, sensor_eqn_y_id )
    end if
    call add_var_vector( ncid, ds%sensor_names, lname_id, nsensor_id, config_name_var_sensor_names, &
      config_desc_var_sensor_names, sensor_name_id )
    call add_statistics( ncid, ds%matchup_statistics, npair_id, &
      stat_m_id, stat_cost_id, stat_cost_reduced_id, stat_mean_id, stat_sdev_id )
    call add_var_vector( ncid, ds%matchup_sensor_names, lname_id, npair_id, config_name_var_matchup_sensor_names, &
      config_desc_var_matchup_sensor_names, matchup_sensor_pairs_id )
    call add_var_vector( ncid, ds%matchup_dataset_uuids, luuid_id, npair_id, config_name_var_matchup_dataset_uuids, &
      config_desc_var_matchup_dataset_uuids, matchup_dataset_uuids_id )

    !! 5) Finish definiton of the dataset
    call finish_dataset( ncid )

    !! 6) Write variables
    call write_variable( ncid, x_id, ds%x )
    call write_variable( ncid, xu_id, ds%xu )
    call write_variable( ncid, xc_id, ds%xc )
    call write_variable( ncid, xr_id, ds%xr )
    call write_variable( ncid, x_add_offset_id, ds%x_add_offsets )
    call write_variable( ncid, x_scale_factor_id, ds%x_scale_factors )
    call write_variable( ncid, x_name_id, ds%x_names )
    call write_variable( ncid, sensor_eqn_id, ds%sensor_eqn_ids )
    call write_variable( ncid, sensor_eqn_n_id, ds%sensor_eqn_n )
    if (must_write_config) then
      call write_variable( ncid, sensor_eqn_o_id, ds%sensor_eqn_o )
      call write_variable( ncid, sensor_eqn_y_id, ds%sensor_eqn_y )
    end if
    call write_variable( ncid, sensor_name_id, ds%sensor_names )
    call write_statistics( ncid, stat_m_id, stat_cost_id, stat_cost_reduced_id, stat_mean_id, stat_sdev_id, ds%matchup_statistics )
    call write_variable( ncid, matchup_sensor_pairs_id, ds%matchup_sensor_names )
    call write_variable( ncid, matchup_dataset_uuids_id, ds%matchup_dataset_uuids )

    !! 7) Close the dataset
    call close_dataset( ncid )
    call log_info( "mod_writer", "completed writing parameter dataset" )
  end subroutine write_parameter_dataset_nc

  subroutine write_residuals_dataset_nc( ds, uuid )
    use mod_config
    use mod_nc

    type(residuals_dataset_t), intent(in) :: ds
    character(len=36), intent(in) :: uuid
    integer :: ncid
    character(len=TYPE_MAX_LEN_PATH) :: dsname

    integer :: m_id
    integer :: t_id
    integer :: k_res_id
    integer :: k_res_unc_id
    integer :: k_res_weights_id
    integer :: meas_i_id
    integer :: meas_j_id
    integer :: meas_i_unc_q_id
    integer :: meas_j_unc_q_id
    integer :: meas_i_unc_x_id
    integer :: meas_j_unc_x_id

    !! 1) Create the dataset
    dsname = get_residuals_dataset_name( ".nc", ds%begin_date, ds%end_date, ds%sensor_i_name, ds%sensor_j_name )
    call log_info( "mod_writer", "writing residuals dataset: ", dsname )
    call create_dataset( dsname, ncid )

    !! 2) Add global attributes
    call add_att_matchup_dataset( ncid, ds%begin_date, ds%end_date, ds%matchup_dataset_uuid )
    call add_global_attribute( ncid, config_name_att_sensor_i, ds%sensor_i_name )
    call add_global_attribute( ncid, config_name_att_sensor_j, ds%sensor_j_name )
    call add_att_software_and_job( ncid )
    call add_att_uuid( ncid, uuid )

    !! 3) Add dimensions
    call add_dimension( ncid, config_name_dim_m, ds%m, m_id )

    !! 4) Add variables
    call add_var_vector( ncid, real( ds%t, cp ), m_id, config_name_var_time, &
      config_desc_var_time, t_id )
    call add_attribute( ncid, t_id, "units", config_unit_var_time )
    call add_var_vector( ncid, real( ds%k_res, cp ), m_id, config_name_var_k_res, &
      config_desc_var_k_res, k_res_id )
    call add_var_vector( ncid, real( ds%k_res_unc, cp ), m_id, config_name_var_k_res_unc, &
      config_desc_var_k_res_unc, k_res_unc_id )
    call add_var_vector( ncid, real( ds%k_res_weights, cp ), m_id, config_name_var_k_res_weights, &
      config_desc_var_k_res_weights, k_res_weights_id )
    call add_var_vector( ncid, real( ds%meas_i, cp ), m_id, "measurand1", &
      "The measurand (e.g. radiance) of the first sensor", meas_i_id )
    call add_var_vector( ncid, real( ds%meas_j, cp ), m_id, "measurand2", &
      "The measurand (e.g. radiance) of the other sensor", meas_j_id )
    call add_var_vector( ncid, real( ds%meas_i_unc_q, cp ), m_id, "measurand1_uncertainty_tel", &
      "The uncertainty of the measurand of the first sensor due to the uncertainty of telemetry", meas_i_unc_q_id )
    call add_var_vector( ncid, real( ds%meas_j_unc_q, cp ), m_id, "measurand2_uncertainty_tel", &
      "The uncertainty of the measurand of the other sensor due to the uncertainty of telemetry", meas_j_unc_q_id )
    call add_var_vector( ncid, real( ds%meas_i_unc_x, cp ), m_id, "measurand1_uncertainty_cal", &
      "The uncertainty of the measurand of the first sensor due to the uncertainty of calibration", meas_i_unc_x_id )
    call add_var_vector( ncid, real( ds%meas_j_unc_x, cp ), m_id, "measurand2_uncertainty_cal", &
      "The uncertainty of the measurand of the other sensor due to the uncertainty of calibration", meas_j_unc_x_id )
    if (config_job_enable_unit /= "") then
      call add_attribute( ncid, k_res_id, "units", config_job_enable_unit )
      call add_attribute( ncid, k_res_unc_id, "units", config_job_enable_unit )
      call add_attribute( ncid, k_res_unc_id, "units", "1 / "//config_job_enable_unit )
      call add_attribute( ncid, meas_i_id, "units", config_job_enable_unit )
      call add_attribute( ncid, meas_j_id, "units", config_job_enable_unit )
      call add_attribute( ncid, meas_i_unc_q_id, "units", config_job_enable_unit )
      call add_attribute( ncid, meas_j_unc_q_id, "units", config_job_enable_unit )
      call add_attribute( ncid, meas_i_unc_x_id, "units", config_job_enable_unit )
      call add_attribute( ncid, meas_j_unc_x_id, "units", config_job_enable_unit )
    end if

    !! 5) Finish definiton of the dataset
    call finish_dataset( ncid )

    !! 6) Write variables
    call write_variable( ncid, t_id, ds%t )
    call write_variable( ncid, k_res_id, real( ds%k_res, cp ) )
    call write_variable( ncid, k_res_unc_id, real( ds%k_res_unc, cp ) )
    call write_variable( ncid, k_res_weights_id, real( ds%k_res_weights, cp ) )
    call write_variable( ncid, meas_i_id, real( ds%meas_i, cp ) )
    call write_variable( ncid, meas_j_id, real( ds%meas_j, cp ) )
    call write_variable( ncid, meas_i_unc_q_id, real( ds%meas_i_unc_q, cp ) )
    call write_variable( ncid, meas_j_unc_q_id, real( ds%meas_j_unc_q, cp ) )
    call write_variable( ncid, meas_i_unc_x_id, real( ds%meas_i_unc_x, cp ) )
    call write_variable( ncid, meas_j_unc_x_id, real( ds%meas_j_unc_x, cp ) )

    !! 7) Close the dataset
    call close_dataset( ncid )
    call log_info( "mod_writer", "completed writing residuals dataset" )
  end subroutine write_residuals_dataset_nc

  subroutine add_att_matchup_dataset( ncid, begin_date, end_date, uuid )
    use mod_nc
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: begin_date
    character(len=*), intent(in) :: end_date
    character(len=36), intent(in), optional :: uuid

    call add_global_attribute( ncid, config_name_att_matchup_dataset_id, config_job_matchup_dataset_id )
    call add_global_attribute( ncid, config_name_att_matchup_dataset_begin, begin_date )
    call add_global_attribute( ncid, config_name_att_matchup_dataset_end, end_date )
    if (present( uuid )) then
      call add_global_attribute( ncid, config_name_att_matchup_dataset_uuid, uuid )
    end if
  end subroutine

  subroutine add_att_statistics( ncid, statistics )
    use mod_nc
    integer, intent(in) :: ncid
    type(matchup_statistics_t), intent(in) :: statistics

    call add_global_attribute( ncid, config_name_att_number_of_matchups, statistics%m )
    call add_global_attribute( ncid, config_name_att_cost, statistics%cost )
    call add_global_attribute( ncid, config_name_att_cost_p_value, statistics%cost_p_value )
    call add_global_attribute( ncid, config_name_att_cost_reduced, statistics%cost_reduced )
    call add_global_attribute( ncid, config_name_att_cost_matchup, statistics%cost_matchup )
    call add_global_attribute( ncid, config_name_att_cost_sensors, statistics%cost_sensors )
    call add_global_attribute( ncid, config_name_att_k_res_mean, statistics%residual_mean )
    call add_global_attribute( ncid, config_name_att_k_res_sdev, statistics%residual_sdev )
  end subroutine

  subroutine add_att_software_and_job( ncid )
    use mod_nc
    integer, intent(in) :: ncid

    call add_global_attribute( ncid, config_name_att_software, config_value_att_software )
    call add_global_attribute( ncid, config_name_att_software_version, config_value_att_software_version )
    call add_global_attribute( ncid, config_name_att_software_tag, config_value_att_software_tag )
    call add_global_attribute( ncid, config_name_att_job, config_job )
    call add_global_attribute( ncid, config_name_att_job_id, config_job_id )
  end subroutine

  subroutine add_att_uuid( ncid, uuid )
    use mod_nc
    integer, intent(in) :: ncid
    character(len=36), intent(in) :: uuid

    call add_global_attribute( ncid, config_name_att_uuid, uuid )
  end subroutine

  subroutine add_var_vector_int( ncid, prototype, dimid, varname, desc, varid )
    use mod_nc
    integer, intent(in) :: ncid
    integer, intent(in) :: prototype(:)
    integer, intent(in) :: dimid
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: desc
    integer, intent(out) :: varid

    call add_variable( ncid, varname, prototype, dimid, varid )
    call add_attribute( ncid, varid, "description", desc )
  end subroutine

  subroutine add_var_vector_real__dp( ncid, prototype, dimid, varname, desc, varid )
    use mod_nc
    integer, intent(in) :: ncid
    real(kind=dp), intent(in) :: prototype(:)
    integer, intent(in) :: dimid
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: desc
    integer, intent(out) :: varid

    call add_variable( ncid, varname, prototype, dimid, varid )
    call add_attribute( ncid, varid, "description", desc )
  end subroutine

  subroutine add_var_vector_real__sp( ncid, prototype, dimid, varname, desc, varid )
    use mod_nc
    integer, intent(in) :: ncid
    real(kind=sp), intent(in) :: prototype(:)
    integer, intent(in) :: dimid
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: desc
    integer, intent(out) :: varid

    call add_variable( ncid, varname, prototype, dimid, varid )
    call add_attribute( ncid, varid, "description", desc )
  end subroutine

  subroutine add_var_vector_text( ncid, prototype, lenid, dimid, varname, desc, varid )
    use mod_nc
    integer, intent(in) :: ncid
    character(len=*), intent(in) :: prototype(:)
    integer, intent(in) :: lenid
    integer, intent(in) :: dimid
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: desc
    integer, intent(out) :: varid

    call add_variable( ncid, varname, prototype, lenid, dimid, varid )
    call add_attribute( ncid, varid, "description", desc )
  end subroutine

  subroutine add_var_matrix( ncid, prototype, dimid, varname, desc, varid )
    use mod_nc
    integer, intent(in) :: ncid
    real(kind=wp), intent(in) :: prototype(:,:)
    integer, intent(in) :: dimid
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: desc
    integer, intent(out) :: varid

    call add_variable( ncid, varname, prototype, dimid, dimid, varid )
    call add_attribute( ncid, varid, "description", desc )
  end subroutine

  subroutine add_statistics( ncid, statistics, dimid, m_id, cost_id, cost_reduced_id, mean_id, sdev_id )
    use mod_config
    use mod_nc

    integer, intent(in) :: ncid
    type(matchup_statistics_t), intent(in) :: statistics(:)
    integer, intent(in) :: dimid
    integer, intent(out) :: m_id
    integer, intent(out) :: cost_id
    integer, intent(out) :: cost_reduced_id
    integer, intent(out) :: mean_id
    integer, intent(out) :: sdev_id

    integer :: counts(size( statistics ))
    real(kind=wp) :: others(size( statistics ))

    call add_var_vector( ncid, counts, dimid, config_name_var_number_of_matchups, config_desc_var_number_of_matchups, m_id )
    call add_var_vector( ncid, others, dimid, config_name_var_cost, config_desc_var_cost, cost_id )
    call add_var_vector( ncid, others, dimid, config_name_var_cost_reduced, config_desc_var_cost_reduced, cost_reduced_id )
    call add_var_vector( ncid, others, dimid, config_name_var_k_res_mean, config_desc_var_k_res_mean, mean_id )
    call add_var_vector( ncid, others, dimid, config_name_var_k_res_sdev, config_desc_var_k_res_sdev, sdev_id )
    if (config_job_enable_unit /= "") then
      call add_attribute( ncid, mean_id, "units", config_job_enable_unit )
      call add_attribute( ncid, sdev_id, "units", config_job_enable_unit )
    end if
  end subroutine

  subroutine write_statistics( ncid, m_id, cost_id, cost_reduced_id, mean_id, sdev_id, statistics )
    use mod_nc
    integer, intent(in) :: ncid
    integer, intent(in) :: m_id
    integer, intent(in) :: cost_id
    integer, intent(in) :: cost_reduced_id
    integer, intent(in) :: mean_id
    integer, intent(in) :: sdev_id
    type(matchup_statistics_t), intent(in) :: statistics(:)

    integer :: counts(size( statistics ))
    real(kind=wp) :: others(size( statistics ))
    integer :: k

    associate (nmatchup => size( statistics ))
      do k = 1, nmatchup
        counts(k) = statistics(k)%m
      end do
      call write_variable( ncid, m_id, counts )
      do k = 1, nmatchup
        others(k) = statistics(k)%cost
      end do
      call write_variable( ncid, cost_id, others )
      do k = 1, nmatchup
        others(k) = statistics(k)%cost_reduced
      end do
      call write_variable( ncid, cost_reduced_id, others )
      do k = 1, nmatchup
        others(k) = statistics(k)%residual_mean
      end do
      call write_variable( ncid, mean_id, others )
      do k = 1, nmatchup
        others(k) = statistics(k)%residual_sdev
      end do
      call write_variable( ncid, sdev_id, others )
    end associate
  end subroutine

end module mod_writer
