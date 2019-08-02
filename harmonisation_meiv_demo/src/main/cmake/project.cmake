# Copyright (C) 2019 FastOpt GmbH, Hamburg, Germany (info@fastopt.de)

function(project_doi DOI)
    set(PROJECT_DOI https://dx.doi.org/${DOI} PARENT_SCOPE)
endfunction()

function(project_url URL)
    set(PROJECT_URL ${URL} PARENT_SCOPE)
endfunction()

function(project_install_prefix PREFIX)
    if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
        set(CMAKE_INSTALL_PREFIX "${PREFIX}" CACHE PATH
                "This directory is prepended onto all install directories"
                FORCE)
    endif ()
endfunction()

function(project_tag TAG)
    set(PROJECT_TAG ${TAG} PARENT_SCOPE)

    if (${TAG} STREQUAL develop)
        include(${CMAKE_BINARY_DIR}/id.cmake OPTIONAL RESULT_VARIABLE VCS_ID)
        if (VCS_ID)       
            set(PROJECT_TAG ${ID} PARENT_SCOPE)
        else ()
            find_program(GIT git)
            if (GIT)
                execute_process(COMMAND ${GIT} rev-parse --short HEAD
                        OUTPUT_VARIABLE ID
                        OUTPUT_STRIP_TRAILING_WHITESPACE
                        RESULT_VARIABLE STATUS
                        ERROR_QUIET)
                if (${STATUS} EQUAL 0)
            	        if (DEFINED ID)
                        set(PROJECT_TAG ${ID} PARENT_SCOPE)
                    endif ()
                endif ()
            endif ()
        endif ()
    endif ()
endfunction()
