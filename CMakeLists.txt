cmake_minimum_required(VERSION 3.16)

project(agnsf
    VERSION 0.1.0
    LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# ----------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------

find_package(Threads REQUIRED)

find_package(Python REQUIRED COMPONENTS
    Interpreter
    Development
)

execute_process(
    COMMAND
        ${Python_EXECUTABLE}
        -m pybind11
        --cmakedir
    OUTPUT_VARIABLE
        pybind11_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

find_package(pybind11 CONFIG REQUIRED)

# cfitsio is used by the agnsf_io library to read FITS light curves.
#
# Detection order:
#   1. pkg-config (conda / homebrew / system)
#   2. manual find_path / find_library (conda prefix, homebrew, /usr)
#   3. download and build cfitsio (ExternalProject) so that a plain
#      `pip install .` works without a pre-installed cfitsio
find_package(PkgConfig QUIET)

set(CFITSIO_FOUND FALSE)

if(PkgConfig_FOUND)
    pkg_check_modules(CFITSIO QUIET IMPORTED_TARGET cfitsio)

    if(TARGET PkgConfig::CFITSIO)
        # Unify on cfitsio::cfitsio so that agnsf_io links the same
        # target regardless of which detection path found cfitsio.
        add_library(cfitsio::cfitsio INTERFACE IMPORTED)

        set_target_properties(cfitsio::cfitsio PROPERTIES
            INTERFACE_LINK_LIBRARIES PkgConfig::CFITSIO
        )

        set(CFITSIO_FOUND TRUE)
    endif()
endif()

if(NOT CFITSIO_FOUND)
    find_path(CFITSIO_INCLUDE_DIR fitsio.h
        HINTS
            ${CONDA_PREFIX}/include
            /opt/homebrew/include
            /usr/local/include
            /usr/include
    )

    find_library(CFITSIO_LIBRARY NAMES cfitsio
        HINTS
            ${CONDA_PREFIX}/lib
            /opt/homebrew/lib
            /usr/local/lib
            /usr/lib
    )

    if(CFITSIO_INCLUDE_DIR AND CFITSIO_LIBRARY)
        add_library(cfitsio::cfitsio UNKNOWN IMPORTED)

        set_target_properties(cfitsio::cfitsio PROPERTIES
            IMPORTED_LOCATION "${CFITSIO_LIBRARY}"
            INTERFACE_INCLUDE_DIRECTORIES "${CFITSIO_INCLUDE_DIR}"
        )

        set(CFITSIO_FOUND TRUE)
    endif()
endif()

if(NOT CFITSIO_FOUND)
    option(
        AGNSF_FETCH_CFITSIO
        "Download and build cfitsio when it is not found on the system"
        ON
    )

    if(AGNSF_FETCH_CFITSIO)
        message(STATUS "cfitsio not found; downloading and building it")

        include(ExternalProject)

        set(AGNSF_CFITSIO_INSTALL
            "${CMAKE_BINARY_DIR}/_deps/cfitsio/install"
        )

        ExternalProject_Add(agnsf_cfitsio
            GIT_REPOSITORY "https://github.com/heasarc/cfitsio.git"
            GIT_TAG "v4.6.3"
            PREFIX "${CMAKE_BINARY_DIR}/_deps/cfitsio"
            UPDATE_COMMAND ""
            CONFIGURE_COMMAND
                <SOURCE_DIR>/configure
                --prefix=<INSTALL_DIR>
                --disable-shared
                --enable-static
            BUILD_COMMAND ${CMAKE_MAKE_PROGRAM} -j
            INSTALL_COMMAND ${CMAKE_MAKE_PROGRAM} install
            INSTALL_DIR "${AGNSF_CFITSIO_INSTALL}"
        )

        add_library(cfitsio::cfitsio UNKNOWN IMPORTED)

        set_target_properties(cfitsio::cfitsio PROPERTIES
            IMPORTED_LOCATION
                "${AGNSF_CFITSIO_INSTALL}/lib/libcfitsio.a"
            INTERFACE_INCLUDE_DIRECTORIES
                "${AGNSF_CFITSIO_INSTALL}/include"
        )

        set(CFITSIO_FOUND TRUE)
    else()
        message(FATAL_ERROR
            "cfitsio is required to build agnsf. Install it "
            "(e.g. homebrew: brew install cfitsio) or re-run with "
            "-DAGNSF_FETCH_CFITSIO=ON."
        )
    endif()
endif()


# ----------------------------------------------------------------------
# Optimization dependency (NLopt)
#
# agnsf::optimization is a thin generic wrapper around NLopt. NLopt is 
# found via its CMake config when available; otherwise it is downloaded and built.
# ----------------------------------------------------------------------

find_package(nlopt CONFIG QUIET)

if(NOT nlopt_FOUND)
    message(STATUS "nlopt not found; downloading and building it")

    include(FetchContent)

    set(NLOPT_CXX ON CACHE BOOL "" FORCE)
    set(NLOPT_PYTHON OFF CACHE BOOL "" FORCE)
    set(NLOPT_TESTS OFF CACHE BOOL "" FORCE)
    set(NLOPT_EXAMPLES OFF CACHE BOOL "" FORCE)
    set(NLOPT_GUILE OFF CACHE BOOL "" FORCE)
    set(NLOPT_MATLAB OFF CACHE BOOL "" FORCE)
    set(NLOPT_OCTAVE OFF CACHE BOOL "" FORCE)
    set(NLOPT_SWIG OFF CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
        nlopt
        GIT_REPOSITORY "https://github.com/stevengj/nlopt.git"
        GIT_TAG "v2.9.1"
    )

    FetchContent_MakeAvailable(nlopt)
endif()

if(TARGET nlopt::nlopt)
    set(AGNSF_NLOPT_TARGET nlopt::nlopt)
elseif(TARGET nlopt)
    set(AGNSF_NLOPT_TARGET nlopt)
else()
    message(FATAL_ERROR "nlopt target not found")
endif()


# ----------------------------------------------------------------------
# Core Library (shared data model; no analysis-specific dependencies)
# ----------------------------------------------------------------------

add_library(agnsf_core
    src/core/light_curve.cpp
)

target_include_directories(agnsf_core
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)


# ----------------------------------------------------------------------
# Optimization Library (agnsf::optimization; generic, no domain deps)
# ----------------------------------------------------------------------

add_library(agnsf_optimization
    src/optimization/optimization.cpp
)

target_include_directories(agnsf_optimization
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(agnsf_optimization
    PRIVATE
        ${AGNSF_NLOPT_TARGET}
)


# ----------------------------------------------------------------------
# Time-delay Library (agnsf::timedelay; depends on agnsf_core only)
# ----------------------------------------------------------------------

add_library(agnsf_timedelay
    src/timedelay/cross_correlation.cpp
    src/timedelay/fr_rss.cpp
    src/timedelay/transfer_function.cpp
)

target_include_directories(agnsf_timedelay
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(agnsf_timedelay
    PRIVATE
        agnsf_core
        agnsf_optimization
)


# ----------------------------------------------------------------------
# IO Library (independent of esf; reusable by future modules)
# ----------------------------------------------------------------------

add_library(agnsf_io
    src/io/light_curve.cpp
    src/io/path_list.cpp
    src/io/table_writer.cpp
)

target_include_directories(agnsf_io
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(agnsf_io
    PRIVATE
        agnsf_core
        cfitsio::cfitsio
)

# When cfitsio is downloaded, make sure it is built before agnsf_io
# is compiled and linked against it.
if(TARGET agnsf_cfitsio)
    add_dependencies(agnsf_io agnsf_cfitsio)

    find_package(ZLIB QUIET)

    if(ZLIB_FOUND)
        target_link_libraries(agnsf_io
            PRIVATE
                ZLIB::ZLIB
        )
    endif()
endif()


# ----------------------------------------------------------------------
# ESF Library (agnsf::esf namespace)
# ----------------------------------------------------------------------

add_library(agnsf_esf
    src/esf/lag_bins.cpp
    src/esf/bin_accumulator.cpp
    src/esf/sf_result.cpp
    src/esf/sf_uncertainty_estimator.cpp
    src/esf/esf_uncertainty.cpp
    src/esf/sf_calculator.cpp
    src/esf/pooled_esf_calculator.cpp
    src/esf/sf_ensemble_calculator.cpp
    src/esf/pair_accumulator.cpp
    src/esf/sf_io.cpp
)

target_include_directories(agnsf_esf
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(agnsf_esf
    PRIVATE
        Threads::Threads
        agnsf_core
        agnsf_io
)


# ----------------------------------------------------------------------
# Python Module
# ----------------------------------------------------------------------

pybind11_add_module(_agnsf
    python/module.cpp
    python/bindings/core.cpp
    python/bindings/esf.cpp
    python/bindings/io.cpp
    python/bindings/sf.cpp
    python/bindings/timedelay.cpp
    python/bindings/uncertainty.cpp
)

target_link_libraries(_agnsf
    PRIVATE
        agnsf_esf
        agnsf_timedelay
)


# ----------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------

# Install the Python extension into the agnsf package.
install(
    TARGETS _agnsf
    LIBRARY DESTINATION agnsf
)

install(
    DIRECTORY python/agnsf/
    DESTINATION agnsf
    FILES_MATCHING
        PATTERN "*.py"
        PATTERN "*.pyi"
        PATTERN "py.typed"
)


# ----------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------

enable_testing()


add_executable(test_lag_bins
    tests/test_lag_bins.cpp
)

target_link_libraries(test_lag_bins
    PRIVATE
        agnsf_esf
)

add_test(
    NAME lag_bins
    COMMAND test_lag_bins
)


add_executable(test_light_curve
    tests/test_light_curve.cpp
)

target_link_libraries(test_light_curve
    PRIVATE
        agnsf_core
)

add_test(
    NAME light_curve
    COMMAND test_light_curve
)


add_executable(test_bin_accumulator
    tests/test_bin_accumulator.cpp
)

target_link_libraries(test_bin_accumulator
    PRIVATE
        agnsf_esf
)

add_test(
    NAME bin_accumulator
    COMMAND test_bin_accumulator
)


add_executable(test_sf_calculator
    tests/test_sf_calculator.cpp
)

target_link_libraries(test_sf_calculator
    PRIVATE
        agnsf_esf
)

add_test(
    NAME sf_calculator
    COMMAND test_sf_calculator
)


add_executable(test_pooled_esf_calculator
    tests/test_pooled_esf_calculator.cpp
)

target_link_libraries(test_pooled_esf_calculator
    PRIVATE
        agnsf_esf
)

add_test(
    NAME pooled_esf_calculator
    COMMAND test_pooled_esf_calculator
)


add_executable(test_sf_ensemble_calculator
    tests/test_sf_ensemble_calculator.cpp
)

target_link_libraries(test_sf_ensemble_calculator
    PRIVATE
        agnsf_esf
)

add_test(
    NAME sf_ensemble_calculator
    COMMAND test_sf_ensemble_calculator
)


add_executable(test_light_curve_io
    tests/test_light_curve_io.cpp
)

target_link_libraries(test_light_curve_io
    PRIVATE
        agnsf_io
        PkgConfig::CFITSIO
)

add_test(
    NAME light_curve_io
    COMMAND test_light_curve_io
)


add_executable(test_path_list
    tests/test_path_list.cpp
)

target_link_libraries(test_path_list
    PRIVATE
        agnsf_io
)

add_test(
    NAME path_list
    COMMAND test_path_list
)


add_executable(test_table_writer
    tests/test_table_writer.cpp
)

target_link_libraries(test_table_writer
    PRIVATE
        agnsf_io
)

add_test(
    NAME table_writer
    COMMAND test_table_writer
)


add_executable(test_sf_io
    tests/test_sf_io.cpp
)

target_link_libraries(test_sf_io
    PRIVATE
        agnsf_esf
)

add_test(
    NAME sf_io
    COMMAND test_sf_io
)


add_executable(test_sf_uncertainty
    tests/test_sf_uncertainty.cpp
)

target_link_libraries(test_sf_uncertainty
    PRIVATE
        agnsf_esf
)

add_test(
    NAME sf_uncertainty
    COMMAND test_sf_uncertainty
)


add_executable(test_esf_uncertainty
    tests/test_esf_uncertainty.cpp
)

target_link_libraries(test_esf_uncertainty
    PRIVATE
        agnsf_esf
)

add_test(
    NAME esf_uncertainty
    COMMAND test_esf_uncertainty
)


add_executable(test_redshift
    tests/test_redshift.cpp
)

target_link_libraries(test_redshift
    PRIVATE
        agnsf_esf
)

add_test(
    NAME redshift
    COMMAND test_redshift
)


add_executable(test_cross_correlation
    tests/test_cross_correlation.cpp
)

target_link_libraries(test_cross_correlation
    PRIVATE
        agnsf_timedelay
)

add_test(
    NAME cross_correlation
    COMMAND test_cross_correlation
)


add_executable(test_optimization
    tests/test_optimization.cpp
)

target_link_libraries(test_optimization
    PRIVATE
        agnsf_optimization
)

add_test(
    NAME optimization
    COMMAND test_optimization
)


add_executable(test_transfer_function
    tests/test_transfer_function.cpp
)

target_link_libraries(test_transfer_function
    PRIVATE
        agnsf_timedelay
)

add_test(
    NAME transfer_function
    COMMAND test_transfer_function
)




# ----------------------------------------------------------------------
# Command-line interface
# ----------------------------------------------------------------------

add_executable(agnsf
    src/cli/main.cpp
)

target_link_libraries(agnsf
    PRIVATE
        agnsf_esf
)

if(DEFINED SKBUILD)
    # Install the CLI into the Python environment's bin/Scripts directory.
    install(
        TARGETS agnsf
        RUNTIME DESTINATION "${SKBUILD_SCRIPTS_DIR}"
    )
else()
    install(
        TARGETS agnsf
        RUNTIME DESTINATION bin
    )
endif()