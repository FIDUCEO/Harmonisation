!> @file mod_config.F90
!! Run and job configurations.
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

!> @brief Module to hold run and job configurations.
!! @authors Ralf Quast
!! @copyright GNU Public License.
module mod_config
  use mod_type
  use mod_info
  use mod_logger

  implicit none

  private

  public config_initialize
  public get_parameter_dataset_name
  public get_residuals_dataset_name

  !> @brief The capacity of the sensor list.
  integer, public, parameter :: CONFIG_CAPACITY_SENSORS = 25

  !> @brief The capacity of the matchup list.
  integer, public, parameter :: CONFIG_CAPACITY_MATCHUP = (CONFIG_CAPACITY_SENSORS * (CONFIG_CAPACITY_SENSORS + 1)) / 2

  !> @brief The capacity of the matrix lists.
  integer, public, parameter :: CONFIG_CAPACITY_MATRICES = CONFIG_CAPACITY_MATCHUP * 20

  !> @brief The maximum number of sensor data columns
  integer, public, parameter :: CONFIG_MAX_SENSOR_COLUMN_COUNT = 32

  !> @brief The maximum length of a sensor configuration
  integer, public, parameter :: CONFIG_MAX_LEN_SENSOR_CONFIG = 512

  !> @brief The pathname of the gradient data file.
  character(len=*), public, parameter :: CONFIG_PATH_POSTADM_DATA = "postadm.dat"

  !> @brief The pathname of the optimum data file.
  character(len=*), public, parameter :: CONFIG_PATH_OPTIMUM_DATA = "optimum.dat"

  !> @brief The pathname of the prior data file.
  character(len=*), public, parameter :: CONFIG_PATH_PRIOR_DATA = "prior.dat"

  !> @brief The pathname of the Jacobian data file.
  character(len=*), public, parameter :: CONFIG_PATH_POSTJAC_DATA = "postjac.dat"

  !> @brief The pathname of the Hessian data file.
  character(len=*), public, parameter :: CONFIG_PATH_POSTHESS_DATA = "posthess.dat"

  !> @brief The pathname of the job configuration.
  character(len=TYPE_MAX_LEN_PATH), public, protected :: config_job = "job.nml"

  !> @brief The job identifier.
  character(len=TYPE_MAX_LEN_TAGS), public, protected :: config_job_id = "01"

  !> @brief To enable the use of covariance matrices.
  logical, public, protected :: config_job_enable_covariance_matrix = .false.

  !> @brief To enable the iterative update of variance or covariance information.
  logical, public, protected :: config_job_enable_covariance_update = .false.

  !> @brief To enable the use of covariance matrices which represent common random errors.
  logical, public, protected :: config_job_enable_covariance_common = .false.

  !> @brief To evaluate the sensor cost function.
  logical, public, protected :: config_job_enable_cost_sensors = .true.

  !> @brief The unit of the measurand.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_job_enable_unit = ""

  !> @brief To enable use of prior data on calibration parameters.
  logical, public, protected :: config_job_enable_prior_data = .false.

  !> @brief The pathname of the prior data file.
  character(len=TYPE_MAX_LEN_PATH), public, protected :: config_job_path_prior_data = CONFIG_PATH_PRIOR_DATA

  !> @brief To enable use of start data.
  logical, public, protected :: config_job_enable_start_data = .false.

  !> @brief The pathname of the start data file.
  character(len=TYPE_MAX_LEN_PATH), public, protected :: config_job_path_start_data = "start.dat"

  !> @brief The number of sensors.
  integer, public, protected :: config_job_sensor_count = 0

  !> @brief The sensor indices.
  integer, public, protected :: config_job_sensors(CONFIG_CAPACITY_SENSORS) = 0

  !> @brief The sensor names.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_job_sensor_names(CONFIG_CAPACITY_SENSORS) = ""

  !> @brief The sensor data columns selected.
  integer, public, protected :: config_job_sensor_columns(CONFIG_CAPACITY_SENSORS,CONFIG_MAX_SENSOR_COLUMN_COUNT) = 0

  !> @brief Select sensor data columns?
  logical, public, protected :: config_job_sensor_columns_select(CONFIG_CAPACITY_SENSORS) = .false.

  !> @brief The sensor configurations.
  real(kind=wp), public, protected :: config_job_sensor_configurations(CONFIG_CAPACITY_SENSORS, &
    CONFIG_MAX_LEN_SENSOR_CONFIG) = 0.0_wp

  !> @brief The number of matchup datasets.
  integer, public, protected :: config_job_matchup_count = 0

  !> @brief The matchup indices.
  integer, public, protected :: config_job_matchup(CONFIG_CAPACITY_MATCHUP, 2) = 0

  !> @brief The matchup datasets
  character(len=TYPE_MAX_LEN_PATH), public, protected :: config_job_matchup_datasets(CONFIG_CAPACITY_MATCHUP) = ""

  !> @brief The matchup names.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_job_matchup_names(CONFIG_CAPACITY_MATCHUP) = ""

  !> @brief The matchup dataset samplings.
  integer, public, protected :: config_job_matchup_samplings(CONFIG_CAPACITY_MATCHUP) = 1

  !> @brief The matchup dataset identifier.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_job_matchup_dataset_id = ""

  !> @brief The matchup dataset directory path.
  character(len=TYPE_MAX_LEN_PATH), public, protected :: config_job_matchup_dataset_path = "./"

  !> @brief The accuracy goal for the linear system solver.
  real(kind=wp), public, protected :: config_solver_accuracy_goal = 1.0E-4_wp

  !> @brief The maximum number of iterations for the conjugate gradient linear system solver.
  integer, public, protected :: config_solver_max_iteration = 100

  !> @brief To enable output of log information at a certain level.
  integer, public, protected :: config_enable_log_level = LOG_LEVEL_INFO

  !> @brief The name of the software identifier attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_software = "software"

  !> @brief The value of the software identifier attribute.
  character(len=TYPE_MAX_LEN_TAGS), public, protected :: config_value_att_software = "FO"

  !> @brief The name of the software version attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_software_version = "software_version"

  !> @brief The value of the software version attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_value_att_software_version = PROJECT_VERSION

  !> @brief The name of the software tag attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_software_tag = "software_tag"

  !> @brief The value of the software tag attribute.
  character(len=TYPE_MAX_LEN_TAGS), public, protected :: config_value_att_software_tag = PROJECT_TAG

  !> @brief The name of the job attribute.
  character(len=TYPE_MAX_LEN_PATH), public, protected :: config_name_att_job = "job"

  !> @brief The name of the job identifier attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_job_id = "job_id"

  !> @brief The name of the (semantic) matchup dataset identifier attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_matchup_dataset_id = "matchup_dataset"

  !> @brief The name of the matchup dataset UUID attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_matchup_dataset_uuid = "matchup_dataset_uuid"

  !> @brief The name of the matchup dataset begin date attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_matchup_dataset_begin = "matchup_dataset_begin"

  !> @brief The name of the matchup dataset end date attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_matchup_dataset_end = "matchup_dataset_end"

  !> @brief The name of the number of matchups attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_number_of_matchups = "number_of_matchups"

  !> @brief The name of the UUID attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_uuid = "uuid"

  !> @brief The name of the cost attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_cost = "cost"

  !> @brief The name of the reduced matchup cost attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_cost_reduced = "cost_reduced"

  !> @brief The name of the cost p-value attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_cost_p_value = "cost_p_value"

  !> @brief The name of the matchup cost attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_cost_matchup = "cost_matchup"

  !> @brief The name of the sensors cost attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_cost_sensors = "cost_sensors"

  !> @brief The name of the mean K residual attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_k_res_mean = "K_res_mean"

  !> @brief The name of the K residual standard deviation attribute.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_k_res_sdev = "K_res_sdev"

  !> @brief The name of the first sensor in a matchup.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_sensor_i = "sensor_1_name"

  !> @brief The name of the other sensor in a matchup.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_att_sensor_j = "sensor_2_name"

  !> @brief The name of the configuration constants dimension.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_dim_c = "n_config"

  !> @brief The name of the name length dimension.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_dim_l_name = "l_name"

  !> @brief The name of the UUID length dimension.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_dim_l_uuid = "l_uuid"

  !> @brief The name of the matchup dimension.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_dim_m = "M"

  !> @brief The name of the harmonisation parameter dimension.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_dim_n = "n"

  !> @brief The name of the sensor pair dimension.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_dim_p = "n_sensor_pair"

  !> @brief The name of the sensor dimension.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_dim_s = "n_sensor"

  !> @brief The name of the parameter variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_x = "parameter"

  !> @brief The description of the parameter variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_x = &
    "The harmonisation parameters"

  !> @brief The name of the parameter add-offset variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_x_add_offset = "parameter_add_offset"

  !> @brief The description of the parameter add-offset variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_x_add_offset = &
    "The harmonisation parameter add offsets (used internally only)"
  
  !> @brief The name of the parameter scaling variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_x_scale_factor = "parameter_scale_factor"

  !> @brief The description of the parameter scaling variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_x_scale_factor = &
    "The harmonisation parameter scale factors (used internally only)"

  !> @brief The name of the parameter uncertainty variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_x_unc = "parameter_uncertainty"

  !> @brief The description of the parameter uncertainty variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_x_unc = &
    "The harmonisation parameter uncertainties"

  !> @brief The name of the parameter covariance matrix variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_x_cov = "parameter_covariance_matrix"

  !> @brief The description of the parameter covariance matrix variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_x_cov = &
    "The harmonisation parameter error covariance matrix"

  !> @brief The name of the parameter correlation matrix variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_x_cor = "parameter_correlation_matrix"

  !> @brief The description of the parameter correlation matrix variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_x_cor = &
    "The harmonisation parameter error correlation matrix"

  !> @brief The name of the parameter names variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_x_names = "parameter_name"

  !> @brief The description of the parameter names variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_x_names = &
    "The names of the harmonisation parameters or the names of the sensors associated with them"

  !> @brief The name of the K residual variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_k_res = "K_res"

  !> @brief The description of the K residual variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_k_res = &
    "The harmonisation residual"

  !> @brief The name of the K residual uncertainty variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_k_res_unc = "K_res_uncertainty"

  !> @brief The description of the K residual uncertainty variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_k_res_unc = &
    "The uncertainty of the harmonisation residual"

  !> @brief The name of the K residual weights variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_k_res_weights = "K_res_weights"

  !> @brief The description of the K residual weights variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_k_res_weights = &
    "The weights of the harmonisation residual"

  !> @brief The name of the number of matchups variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_number_of_matchups = "number_of_matchups"

  !> @brief The description of the number of matchups variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_number_of_matchups = &
    "The number of matchups in a matchup dataset"

  !> @brief The name of the cost variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_cost = "cost"

  !> @brief The description of the cost variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_cost = &
    "The cost associated with a matchup dataset"

  !> @brief The name of the reduced matchup cost variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_cost_reduced = "cost_reduced"

  !> @brief The description of the reduced matchup cost variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_cost_reduced = &
    "The reduced cost associated with a matchup dataset"

  !> @brief The name of the mean K residual variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_k_res_mean = "K_res_mean"

  !> @brief The description of the mean K residual variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_k_res_mean = &
    "The mean harmonisation residual"

  !> @brief The name of the K residual standard deviation variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_k_res_sdev = "K_res_sdev"

  !> @brief The description of the K residual standard deviation variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_k_res_sdev = &
    "The standard deviation of the harmonisation residual"

  !> @brief The name of the matchup sensor names variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_matchup_sensor_names = "matchup_sensor_names"

  !> @brief The description of the matchup sensor names variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_matchup_sensor_names = &
    "The sensors associated with a matchup dataset"

  !> @brief The name of the matchup dataset UUIDs variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_matchup_dataset_uuids = "matchup_dataset_uuid"

  !> @brief The description of the matchup dataset UUIDs variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_matchup_dataset_uuids = &
    "The UUID associated with a matchup dataset"

  !> @brief The name of the time variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_time = "t"

  !> @brief The description of the time variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_time = &
    "The time of measurement"

  !> @brief The unit of the time variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_unit_var_time = &
    "seconds since 1970-1-1 0:0:0"

  !> @brief The name of the sensor equation identifier variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_sensor_equation_id = &
    "sensor_equation_id"

  !> @brief The description of the sensor equation identifier variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_sensor_equation_id = &
    "The sensor equation identifiers"

  !> @brief The name of the sensor equation parameter count variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_sensor_equation_x_count = &
    "sensor_equation_parameter_count"

  !> @brief The description of the sensor equation parameter count variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_sensor_equation_x_count = &
    "The number of sensor equation parameters"

  !> @brief The name of the sensor equation configuration constant count variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_sensor_equation_y_count = &
    "sensor_equation_config_count"

  !> @brief The description of the sensor equation configuration constant count variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_sensor_equation_y_count = &
    "The number of sensor equation configuration constants"

  !> @brief The name of the sensor equation configuration constants variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_sensor_equation_y = &
    "sensor_equation_config"

  !> @brief The description of the sensor equation configuration constants.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_sensor_equation_y = &
    "The sensor equation configuration constants"

  !> @brief The name of the sensor names variable.
  character(len=TYPE_MAX_LEN_NAME), public, protected :: config_name_var_sensor_names = &
    "sensor_name"

  !> @brief The description of the sensor names variable.
  character(len=TYPE_MAX_LEN_DESC), public, protected :: config_desc_var_sensor_names = &
    "The sensor names"

  !> @brief The name of the parameter dataset.
  character(len=TYPE_MAX_LEN_PATH) :: config_name_parameter_dataset = ""

  !> @brief The name of the residuals dataset.
  character(len=TYPE_MAX_LEN_PATH) :: config_name_residuals_dataset = ""

  !> @brief The namelist to read the job preferences.
  namelist /job/                             &
    config_job_id,                           &
    config_job_enable_cost_sensors,          &
    config_job_enable_covariance_matrix,     &
    config_job_enable_covariance_update,     &
    config_job_enable_covariance_common,     &
    config_job_enable_unit,                  &
    config_job_enable_prior_data,            &
    config_job_enable_start_data,            &
    config_job_path_prior_data,              &
    config_job_path_start_data,              &
    config_job_sensor_count,                 &
    config_job_sensors,                      &
    config_job_sensor_names,                 &
    config_job_sensor_columns_select,        &
    config_job_sensor_columns,               &
    config_job_sensor_configurations,        &
    config_job_matchup_count,                &
    config_job_matchup,                      &
    config_job_matchup_datasets,             &
    config_job_matchup_names,                &
    config_job_matchup_samplings,            &
    config_job_matchup_dataset_id,           &
    config_job_matchup_dataset_path

  !> @brief The namelist to read the run configuration.
  namelist /run/                             &
    config_job,                              &
    config_enable_log_level,                 &
    config_solver_accuracy_goal,             &
    config_solver_max_iteration,             &
    !! Attribute and variable names and descriptions
    config_name_att_cost,                    &
    config_name_att_cost_reduced,            &
    config_name_att_cost_p_value,            &
    config_name_att_cost_matchup,            &
    config_name_att_cost_sensors,            &
    config_name_att_job,                     &
    config_name_att_job_id,                  &
    config_name_att_uuid,                    &
    config_name_att_k_res_mean,              &
    config_name_att_k_res_sdev,              &
    config_name_att_matchup_dataset_id,      &
    config_name_att_matchup_dataset_begin,   &
    config_name_att_matchup_dataset_end,     &
    config_name_att_matchup_dataset_uuid,    &
    config_name_att_number_of_matchups,      &
    config_name_att_sensor_i,                &
    config_name_att_sensor_j,                &
    config_name_att_software,                &
    config_name_att_software_version,        &
    config_name_att_software_tag,            &
    config_name_dim_c,                       &
    config_name_dim_l_name,                  &
    config_name_dim_l_uuid,                  &
    config_name_dim_m,                       &
    config_name_dim_n,                       &
    config_name_var_cost,                    &
    config_desc_var_cost,                    &
    config_name_var_cost_reduced,            &
    config_desc_var_cost_reduced,            &
    config_name_var_time,                    &
    config_desc_var_time,                    &
    config_unit_var_time,                    &
    config_name_var_x,                       &
    config_desc_var_x,                       &
    config_name_var_x_add_offset,            &
    config_desc_var_x_add_offset,            &
    config_name_var_x_scale_factor,          &
    config_desc_var_x_scale_factor,          &
    config_name_var_x_cor,                   &
    config_desc_var_x_cor,                   &
    config_name_var_x_cov,                   &
    config_desc_var_x_cov,                   &
    config_name_var_x_names,                 &
    config_desc_var_x_names,                 &
    config_name_var_x_unc,                   &
    config_desc_var_x_unc,                   &
    config_name_var_k_res,                   &
    config_desc_var_k_res,                   &
    config_name_var_k_res_unc,               &
    config_desc_var_k_res_unc,               &
    config_name_var_k_res_weights,           &
    config_desc_var_k_res_weights,           &
    config_name_var_k_res_mean,              &
    config_desc_var_k_res_mean,              &
    config_name_var_k_res_sdev,              &
    config_desc_var_k_res_sdev,              &
    config_name_var_matchup_sensor_names,    &
    config_desc_var_matchup_sensor_names,    &
    config_name_var_matchup_dataset_uuids,   &
    config_desc_var_matchup_dataset_uuids,   &
    config_name_var_number_of_matchups,      &
    config_desc_var_number_of_matchups,      &
    config_name_var_sensor_equation_id,      &
    config_desc_var_sensor_equation_id,      &
    config_name_var_sensor_equation_x_count, &
    config_desc_var_sensor_equation_x_count, &
    config_name_var_sensor_equation_y,       &
    config_desc_var_sensor_equation_y,       &
    config_name_var_sensor_equation_y_count, &
    config_desc_var_sensor_equation_y_count, &
    config_name_var_sensor_names,            &
    config_desc_var_sensor_names,            &
    !! Parameter and residual dataset names
    config_name_parameter_dataset,           &
    config_name_residuals_dataset,           &
    !! Value of the software attribute
    config_value_att_software

contains

  !> @brief Initializes the module.
  subroutine config_initialize
    integer, parameter :: unit = 12

    call logger_initialize( config_enable_log_level )
    call read_preferences( unit )
    call read_run_configuration( unit )
    call logger_initialize( config_enable_log_level )
    call read_job_configuration( unit )
  end subroutine config_initialize

  !> Returns the name of the harmonisation parameter dataset.
  !!
  !! @param[in] ext The extension of the file name.
  !! @param[in] begin_date The date of the beginning of the harmonisation period.
  !! @param[in] end_date The date of the end of the harmonisation period.
  !! @return the name of the harmonisation parameter dataset.
  function get_parameter_dataset_name( ext, begin_date, end_date ) result (name)
    character(len=*), intent(in) :: ext
    character(len=*), intent(in) :: begin_date
    character(len=*), intent(in) :: end_date
    character(len=TYPE_MAX_LEN_PATH) :: name

    character(len=*), parameter  :: pre = "harm"
    character(len=1), parameter  :: sep = "_"

    if (len( config_name_parameter_dataset ) > 0) then
      name = pre
      name = trim( name )//sep//trim( config_value_att_software )
      name = trim( name )//sep//trim( config_value_att_software_version )
      name = trim( name )//sep//trim( config_value_att_software_tag )
      name = trim( name )//sep//trim( config_job_id )
      name = trim( name )//sep//trim( config_job_matchup_dataset_id )
      name = trim( name )//sep//trim( begin_date )
      name = trim( name )//sep//trim( end_date )
      name = trim( name )//trim( ext )
    else
      name = trim( config_name_parameter_dataset )//trim( ext )
    end if
  end function get_parameter_dataset_name

  !> Returns the file name of the harmonisation residuals dataset.
  !!
  !! @param[in] ext The extension of the file name.
  !! @param[in] begin_date The date of the beginning of the harmonisation period.
  !! @param[in] end_date The date of the end of the harmonisation period.
  !! @param[in] sensor_i_name The name of the first sensor.
  !! @param[in] sensor_j_name The name of the other sensor.
  !! @return the name of the harmonisation residuals dataset.
  function get_residuals_dataset_name( ext, begin_date, end_date, sensor_i_name, sensor_j_name ) result (name)
    character(len=*), intent(in) :: ext
    character(len=*), intent(in) :: begin_date
    character(len=*), intent(in) :: end_date
    character(len=*), intent(in) :: sensor_i_name
    character(len=*), intent(in) :: sensor_j_name
    character(len=TYPE_MAX_LEN_PATH)  :: name

    character(len=*), parameter  :: pre = "harm"
    character(len=1), parameter  :: sep = "_"

    if (config_name_residuals_dataset == "") then
      name = pre
      name = trim( name )//sep//trim( config_value_att_software )
      name = trim( name )//sep//trim( config_value_att_software_version )
      name = trim( name )//sep//trim( config_value_att_software_tag )
      name = trim( name )//sep//trim( config_job_id )
      name = trim( name )//sep//trim( config_job_matchup_dataset_id )
      name = trim( name )//sep//trim( begin_date )
      name = trim( name )//sep//trim( end_date )
      name = trim( name )//sep//"res"
      name = trim( name )//sep//trim( sensor_i_name )
      name = trim( name )//sep//trim( sensor_j_name )
      name = trim( name )//trim( ext )
    else
      name = config_name_residuals_dataset//sep//trim( sensor_i_name )//"_"//trim( sensor_j_name )//trim( ext )
    end if
  end function get_residuals_dataset_name

  subroutine read_preferences( unit )
    integer, intent(in) :: unit
    character(len=*), parameter :: path = "preferences.nml"
    logical :: answer

    inquire (file=path, exist=answer)
    if (answer) then
      call log_debug( "", "reading preferences ", path )

      open (unit, file=path)
      read (unit, nml=run)
      close (unit)

      call log_debug( "", "reading preferences completed" )
    end if
  end subroutine

  subroutine read_run_configuration( unit )
    integer, intent(in) :: unit
    character(len=*), parameter :: path = "run.nml"
    logical :: answer

    inquire (file=path, exist=answer)
    if (answer) then
      call log_debug( "", "reading run configuration ", path )

      open (unit, file=path)
      read (unit, nml=run)
      close (unit)

      call log_debug( "", "reading run configuration completed" )
    else
      call log_warning( "", "no run configuration specified, using default run configuration" )
    end if
  end subroutine

  subroutine read_job_configuration( unit )
    integer, intent(in) :: unit
    logical :: answer
    integer :: status

    inquire (file=trim( config_job ), exist=answer)
    if (answer) then
      call log_debug( "mod_config", "reading job configuration ", config_job )

      open (unit, file=trim( config_job ))
      read (unit, nml=job, iostat=status)
      if (status /= 0) then
        call log_error( "mod_config", "faulty job configuration ", config_job )
        error stop "Faulty job configuration"
      end if
      close (unit)

      call log_debug( "mod_config", "reading job configuration completed" )
    else
      call log_warning( "mod_config", "no job configuration specified or found ", config_job )
      call log_warning( "mod_config", "reading job configuration from standard input (pipe or keyboard)" )

      read (*, nml=job, iostat=status)
      if (status /= 0) then
        call log_error( "mod_config", "faulty job configuration" )
        error stop "Faulty job configuration"
      end if

      call log_debug( "mod_config", "reading job configuration completed" )
    end if
  end subroutine

end module mod_config
