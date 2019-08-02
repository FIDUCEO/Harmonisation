# Copyright (C) 2019 FastOpt GmbH, Hamburg, Germany (info@fastopt.de)

macro(netcdf_required)
    find_library(NETCDF_Fortran_LIBRARY netcdff
            PATHS $ENV{HOME}/lib $ENV{HOME}/usr/lib /usr/local/lib /usr/local/lib64 /usr/lib /usr/lib64
            DOC "The path to the netCDF Fortran library")
    if (NETCDF_Fortran_LIBRARY)
        message(STATUS "NetCDF Fortran library found: ${NETCDF_Fortran_LIBRARY}")
        get_filename_component(NETCDF_Fortran_LIBRARY_DIR ${NETCDF_Fortran_LIBRARY} DIRECTORY)
        get_filename_component(NETCDF_Fortran_ROOT ${NETCDF_Fortran_LIBRARY_DIR} DIRECTORY)
        find_file(NETCDF_Fortran_MODULE netcdf.mod
                PATHS ${NETCDF_Fortran_ROOT} ${NETCDF_Fortran_ROOT}/include
                DOC "The path to the netCDF Fortran module"
                NO_DEFAULT_PATH)
        if (NETCDF_Fortran_MODULE)
            message(STATUS "NetCDF Fortran module found: ${NETCDF_Fortran_MODULE}")
            get_filename_component(NETCDF_Fortran_INCLUDE ${NETCDF_Fortran_MODULE} DIRECTORY)
            include_directories(${NETCDF_Fortran_INCLUDE})
        else ()
            message(SEND_ERROR "NetCDF Fortran module netcdf.mod not found.")
        endif ()
        find_library(NETCDF_LIBRARY netcdf
                PATHS $ENV{HOME}/lib $ENV{HOME}/usr/lib /usr/local/lib /usr/local/lib64 /usr/lib /usr/lib64
                DOC "The path to the netCDF library")
        if (NETCDF_LIBRARY)
            message(STATUS "NetCDF library found: ${NETCDF_LIBRARY}")
            set(NETCDF_LIBRARIES ${NETCDF_Fortran_LIBRARY} ${NETCDF_LIBRARY} CACHE STRING
                    "The libraries for reading and writing netCDF")
        else ()
            message(SEND_ERROR "NetCDF library netcdf not found.")
        endif ()
    else ()
        message(SEND_ERROR "NetCDF Fortran library netcdff not found.")
    endif ()
endmacro()
