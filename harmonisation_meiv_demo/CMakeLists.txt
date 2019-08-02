#
# Copyright (C) 2019 FastOpt GmbH, Hamburg, Germany (info@fastopt.de)
# All rights reserved.
#

cmake_minimum_required(VERSION 3.5.0)
include(src/main/cmake/project.cmake)
include(src/main/cmake/fortran.cmake)
include(src/main/cmake/unix.cmake)
include(src/main/cmake/netcdf_required.cmake)
include(src/main/cmake/openmp_optional.cmake)


##
# Project
#
project(harmonisation_meiv_demo VERSION 3.0 LANGUAGES Fortran)
project_tag(develop)
project_install_prefix($ENV{HOME}/apps)
project_url(http://http://www.fiduceo.eu)

set(CMAKE_CONFIGURATION_TYPES Release Debug)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set_property(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES mod_info.F90 processing resources)

netcdf_required()
openmp_optional()

##
# Directories
#
set(MAIN ${CMAKE_SOURCE_DIR}/src/main)
set(TEST ${CMAKE_SOURCE_DIR}/src/test)
set(PROC ${CMAKE_BINARY_DIR}/processing CACHE PATH "The path to the data processing directory")


##
# Configuration
#
configure_file(${MAIN}/fortran-filtered/gpl/mod_info.F90 mod_info.F90 @ONLY)

configure_file(${MAIN}/resources-filtered/jobs/job_ph____v1_EIV_101_01.nml resources/jobs/job_ph____v1_EIV_101_01.nml @ONLY)
configure_file(${MAIN}/resources-filtered/jobs/job_ph____v1_EIV_111_01.nml resources/jobs/job_ph____v1_EIV_111_01.nml @ONLY)
configure_file(${MAIN}/resources-filtered/jobs/job_ph____v1_ODR_101_01.nml resources/jobs/job_ph____v1_ODR_101_01.nml @ONLY)
configure_file(${MAIN}/resources-filtered/jobs/job_ph____v1_ODR_111_01.nml resources/jobs/job_ph____v1_ODR_111_01.nml @ONLY)

configure_file(${MAIN}/resources-filtered/templates/template_job.nml resources/templates/template_job.nml @ONLY)
configure_file(${MAIN}/resources-filtered/templates/template_run.nml resources/templates/template_run.nml @ONLY)

configure_file(${TEST}/resources/ph_v1/0_1.nc ${PROC}/source/test/ph_v1/0_1.nc COPYONLY)
configure_file(${TEST}/resources/ph_v1/1_2.nc ${PROC}/source/test/ph_v1/1_2.nc COPYONLY)
configure_file(${TEST}/resources/ph_v1/2_3.nc ${PROC}/source/test/ph_v1/2_3.nc COPYONLY)
configure_file(${TEST}/resources/run.nml ${CMAKE_BINARY_DIR}/run.nml COPYONLY)


##
# Main targets
#
add_executable(harmonisation-demo
        ${MAIN}/fortran/gpl/mod_base.F90
        ${MAIN}/fortran/gpl/mod_type.F90
        ${MAIN}/fortran/gpl/mod_time.F90
        ${MAIN}/fortran/gpl/mod_logger.F90
        ${CMAKE_BINARY_DIR}/mod_info.F90
        ${MAIN}/fortran/gpl/mod_config.F90
        ${MAIN}/fortran/gpl/mod_sensors.F90
        ${MAIN}/fortran/gpl/mod_assert.F90
        ${MAIN}/fortran/gpl/mod_error_handler.F90
        ${MAIN}/fortran/gpl/mod_matrices.F90
        ${MAIN}/fortran/gpl/mod_matchup.F90
        ${MAIN}/fortran/gpl/mod_nc.F90
        ${MAIN}/fortran/gpl/mod_reader.F90
        ${MAIN}/fortran/gpl/mod_writer.F90
        ${MAIN}/fortran/gpl/demo.F90)
target_link_libraries(harmonisation-demo ${NETCDF_LIBRARIES})


##
# Test targets
#
enable_testing()
add_test(NAME demo_test COMMAND harmonisation-demo)
set_tests_properties(demo_test PROPERTIES TIMEOUT 10)


##
# Custom targets
#
add_custom_target(tar
        COMMAND ${TAR} -f ${PROJECT_NAME}.tar.gz -C ${CMAKE_SOURCE_DIR} --exclude .git --exclude doc -X ${CMAKE_SOURCE_DIR}/.gitignore .
        BYPRODUCTS ${PROJECT_NAME}.tar.gz
        VERBATIM)
