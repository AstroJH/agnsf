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


# ----------------------------------------------------------------------
# C++ Library
# ----------------------------------------------------------------------

add_library(esf
    src/lag_bins.cpp
    src/light_curve.cpp
    src/light_curve_view.cpp
    src/bin_accumulator.cpp
    src/sf_result.cpp
    src/sf_calculator.cpp
    src/pooled_esf_calculator.cpp
    src/sf_ensemble_calculator.cpp
    src/pair_accumulator.cpp
)

target_link_libraries(esf
    PRIVATE
        Threads::Threads
)

target_include_directories(esf
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)


# ----------------------------------------------------------------------
# Python Module
# ----------------------------------------------------------------------

pybind11_add_module(_agnsf
    python/module.cpp
)

target_link_libraries(_agnsf
    PRIVATE
        esf
)


# ----------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------

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
        esf
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
        esf
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
        esf
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
        esf
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
        esf
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
        esf
)

add_test(
    NAME sf_ensemble_calculator
    COMMAND test_sf_ensemble_calculator
)