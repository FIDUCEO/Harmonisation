# Harmonisation of satellite data series (for demonstration only)

The code included with this project is an excerpt from the full harmonisation software, which is described in detail by

Giering, R.; Quast, R.; Mittaz, J.P.D.; Hunt, S.E.; Harris, P.M.; Woolliams, E.R.; Merchant, C.J.
A Novel Framework to Harmonise Satellite Data Series for Climate Applications. Remote Sens. 2019, 11, 1002
<https://doi.org/10.3390/rs11091002>

The code included with this project is for demonstration only and does not include essential computational parts.
Please contact [FastOpt GmbH](<www.fastopt.de) for further information on the full harmonisation software.

## Building and installing
 
Building this software requires [CMake](https://cmake.org) and a compiler that implements the Fortran 2008 standard and [netCDF](https://www.unidata.ucar.edu/software/netcdf/). To use a specific Fortran compiler set the `FC` and `CC` environment variables before you build this software. For example:

    export FC=gfortran
    export CC=gcc

To build and test this software `cd` into the project root directory and type:

    mkdir cmake-build
    cd cmake-build
    cmake ..
    make clean all test

Typing `make test` shall yield console output like

    Running tests...
    Test project harmonisation_meiv_demo/cmake-build
        Start 1: demo_test
    1/1 Test #1: demo_test ........................   Passed    0.11 sec

    100% tests passed, 0 tests failed out of 1

    Total Test time (real) =   0.11 sec

and produce a file name like

    harm_FO_3.0_ggggggg_EIV_101_01_R_W__PH______TEST_01_RSW___.nc

which contains (dummy) calibration coefficients, calibration uncertainties, calibration error covariance and correlation matrices, and some result statistics.

### Adding a specific directory to search for [netCDF](https://www.unidata.ucar.edu/software/netcdf/)
       
To search for a specific netCDF installation add the directory to be searched to the `CMAKE_PREFIX_PATH`. For example:

    cmake -DCMAKE_PREFIX_PATH=/usr/local/netcdf-4.7-gfortran ..

### Editing the CMake configuration

To edit the CMake configuration `cd` into the project build directory and type:

    ccmake .
    
and follow the instruction displayed on the screen.

## Executables

The project produces a single executable program

    harmonisation-demo

The demonstration merely reads a certain run and job configuration, reads the matchup datasets associated with the job, and writes essentially empty output datasets. *No optimisation is performed.*

## Job configuration

The job configuration file is a Fortran namelist that consists entries that define the processing job. For example:

    config_job_id                        = "EIV_101_01_R_W_"
    
Defines the job identifier.

    config_job_enable_covariance_matrix  = T   
    
Defines whether the job uses the full covariance matrix (`T`) or its diagonal elements only (`F`).
 
    config_job_enable_covariance_update  = T
    
Defines whether the job updates the covariance matrix (or its diagonal elements) every iteration (`T`) or never (`F`).
    
    config_job_sensor_count              = 4
    
Defines the number of sensors considered (including the reference sensor).
    
    config_job_sensors(1)                = 0
    config_job_sensors(2)                = 101
    config_job_sensors(3)                = 101
    config_job_sensors(4)                = 101
    
Define the measurement equation to use for each sensor. Each measurement equation is defined by a numeric identifier (see measurement equation identifiers further below).

    config_job_sensor_configurations(2,:) = 0.0, 1.0, 0.0, 1.0E-02, 0.0, 1.0E-06
    config_job_sensor_configurations(3,:) = 0.0, 1.0, 0.0, 1.0E-02, 0.0, 1.0E-06
    config_job_sensor_configurations(4,:) = 0.0, 1.0, 0.0, 1.0E-02, 0.0, 1.0E-06
    
Define the configuration for each sensor. The sensor configuration consists of a sequence of pairs of add-offset and scale factor values for each calibration coefficient (three calibration coefficients in the example above). Any additional sensor configuration constants required by the measurement equation may be appended (no additional sensor configuration constants are required in the example above). 
     
    config_job_sensor_names(1)           = "0"
    config_job_sensor_names(2)           = "1"
    config_job_sensor_names(3)           = "2"
    config_job_sensor_names(4)           = "3"
    
Define the name of each sensor.

    config_job_matchup_count             = 3
    
Defines the number of matchup sensor pairs.

    config_job_matchup(1,:)              = 1, 2
    config_job_matchup(2,:)              = 2, 3
    config_job_matchup(3,:)              = 3, 4
    
Define the matchup sensor pairs. Each sensor is referenced by its index number.

    config_job_matchup_datasets(1)       = "0_1.nc"
    config_job_matchup_datasets(2)       = "1_2.nc"
    config_job_matchup_datasets(3)       = "2_3.nc"
    
Define the matchup data files.

    config_job_matchup_dataset_id        = "PH______TEST_01_RSW_"
    
Defines the matchup dataset identifier.

    config_job_matchup_dataset_path      = "@PROC@/source/test/ph_v1/"

Defines the path where the matchup dataset is stored. `@PROC@` is a variable that is replaced with the actual processing directory during configuration of the harmonisation software.

Job configurations for a few standard tests are predefined and included with this demonstration software. The full software allows you to select specific columns of telemetry data from your matchup input datasets. Thus, it is feasible to re-order columns and it is also feasible to treat the uncertainty (of common random errors) associated with a telemetry data column in a special way, which allows you to determine an otherwise unknown common random error as result of the optimisation. For example: 

    config_job_sensor_columns(2,:)       = 1, 2, 3, 4, 8
    config_job_sensor_columns(3,:)       = 1, 2, 3, 4, 8
    config_job_sensor_columns(4,:)       = 1, 2, 3, 4, 8
    config_job_sensor_columns_select(2)  = T
    config_job_sensor_columns_select(3)  = T
    config_job_sensor_columns_select(4)  = T

Here, the measurement equation uses five data columns for input, where columns `1-4` correspond to sensor telemetry data and column `5` in fact corresponds to the uncertainty (of common random errors) asociated with the telemetry in column `4`. Internally, data columns are enumerated like `1, ..., m, m + 1, ..., m + m` where the first `m` columns correspond to sensor telemetry and the second `m` columns correspond to the uncertainty (of common random errors) associated with the first `m` columns. In practice, however, this feature, has been applied to simulated test datasets only. 

## Measurement equations

The full version of the harmonisation software includes measurement equations for several sensors. In particular

| ID  | Sensor    | Channel | Number of calibration coefficients | Remark                                                                                                         |
|-----|-----------|---------|------------------------------------|----------------------------------------------------------------------------------------------------------------|
| 0   | reference | any     | 0                                  | no calibration coefficients                                                                                    |
| 101 | AVHRR     | 11, 12  | 3                                  | add-offset, emissivity correction, non-linearity                                                               |
| 102 | AVHRR     | 11, 12  | 4                                  | add-offset, emissivity correction, non-linearity, environment temperature coefficient                          |
| 105 | AVHRR     | 3.7     | 2                                  | add-offset, emissivity correction                                                                              |
| 106 | AVHRR     | 3.7     | 3                                  | add-offset, emissivity correction, environment temperature coefficient                                         |
| 111 | AVHRR     | 11, 12  | 4                                  | add-offset, emissivity correction, non-linearity, ICT radiance correction                                      |
| 112 | AVHRR     | 11, 12  | 5                                  | add-offset, emissivity correction, non-linearity, environment temperature coefficient, ICT radiance correction |
| 115 | AVHRR     | 3.7     | 2                                  | add-offset, emissivity correction, ICT radiance correction                                                     |
| 116 | AVHRR     | 3.7     | 3                                  | add-offset, emissivity correction, environment temperature coefficient, ICT radiance correction                |
| 201 | HIRS      | any     | 2                                  | add-offset, non-linearity                                                                                      |
| 202 | HIRS      | any     | 3                                  | add-offset, emissivity correction, non-linearity                                                               |
| 301 | MHS       | any     | 2                                  | non-linearity, warm temperature correction                                                                     |
| 302 | MHS       | any     | 3                                  | non-linearity, warm temperature correction, polarisation correction                                            |
| 303 | MHS       | any     | 4                                  | non-linearity, warm temperature correction, cold temperature correction, polarisation correction               |
| 401 | AMSU-B    | any     | 2                                  | non-linearity, warm temperature correction                                                                     |
| 402 | AMSU-B    | any     | 3                                  | non-linearity, warm temperature correction, polarisation correction                                            |
| 403 | AMSU-B    | any     | 4                                  | non-linearity, warm temperature correction, cold temperature correction, polarisation correction               |

It is straightforward to add more measurement equations to the Fortran code. Fortran code to evaluate partial derivatives of all measurement equations is generated automatically by means of [Transformation of Algorithms in Fortran](http://www.fastopt.de/products/taf/taf.shtml) which is integrated into the software build process. The demonstration software includes about 4,000 lines of Fortran code. The full software includes about 14,000 lines of hand-written Fortran code and about 36,000 lines of generated derivative Fortran code.
