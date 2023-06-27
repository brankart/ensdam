cmake_minimum_required(VERSION 3.15)

set(project_name "pyensdam")

if(USE_MPI)
  set(project_name "${project_name}_mpi")
endif()

project(${project_name} LANGUAGES C Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find target Python interpreter.
#set(Python_FIND_STRATEGY LOCATION)
#set(Python3_FIND_STRATEGY LOCATION)
find_package(Python3 COMPONENTS Interpreter Development NumPy)

add_library(censdam ensscores_wrapper.F90 obserror_wrapper.F90 stochtools_wrapper.F90 ensstat_wrapper.F90 ensanam_wrapper.F90 interp_wrapper.F90 ensaugm_wrapper.F90 ensupdate_wrapper.F90 transpho_wrapper.F90 )
target_link_libraries(censdam ensdam_static)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/scores.c
    COMMAND cython scores.pyx
    MAIN_DEPENDENCY scores.pyx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
Python3_add_library(scores MODULE scores.c)
target_link_libraries(scores PRIVATE censdam Python3::NumPy)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/obserror.c
    COMMAND cython obserror.pyx
    MAIN_DEPENDENCY obserror.pyx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
Python3_add_library(obserror MODULE obserror.c)
target_link_libraries(obserror PRIVATE censdam Python3::NumPy)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/random.c
    COMMAND cython random.pyx
    MAIN_DEPENDENCY random.pyx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
Python3_add_library(random MODULE random.c)
target_link_libraries(random PRIVATE censdam Python3::NumPy)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/probability.c
    COMMAND cython probability.pyx
    MAIN_DEPENDENCY probability.pyx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
Python3_add_library(probability MODULE probability.c)
target_link_libraries(probability PRIVATE censdam Python3::NumPy)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/statistics.c
    COMMAND cython statistics.pyx
    MAIN_DEPENDENCY statistics.pyx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
Python3_add_library(statistics MODULE statistics.c)
target_link_libraries(statistics PRIVATE censdam Python3::NumPy)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/anamorphosis.c
    COMMAND cython anamorphosis.pyx
    MAIN_DEPENDENCY anamorphosis.pyx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
Python3_add_library(anamorphosis MODULE anamorphosis.c)
target_link_libraries(anamorphosis PRIVATE censdam Python3::NumPy)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/interpolation.c
    COMMAND cython interpolation.pyx
    MAIN_DEPENDENCY interpolation.pyx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
Python3_add_library(interpolation MODULE interpolation.c)
target_link_libraries(interpolation PRIVATE censdam Python3::NumPy)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/augment.c
    COMMAND cython augment.pyx
    MAIN_DEPENDENCY augment.pyx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
Python3_add_library(augment MODULE augment.c)
target_link_libraries(augment PRIVATE censdam Python3::NumPy)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/update.c
    COMMAND cython update.pyx
    MAIN_DEPENDENCY update.pyx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
Python3_add_library(update MODULE update.c)
target_link_libraries(update PRIVATE censdam Python3::NumPy)

add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/transpho.c
    COMMAND cython transpho.pyx
    MAIN_DEPENDENCY transpho.pyx
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
Python3_add_library(transpho MODULE transpho.c)
target_link_libraries(transpho PRIVATE censdam Python3::NumPy)

