cmake_minimum_required(VERSION 3.15)
project(pyensdam LANGUAGES C Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find target Python interpreter.
#set(Python_FIND_STRATEGY LOCATION)
#set(Python3_FIND_STRATEGY LOCATION)
find_package(Python3 COMPONENTS Interpreter Development NumPy)

add_library(censdam ensscores_wrapper.F90 obserror_wrapper.F90 )
#add_library(censdam ensscore_wrapper.F90 pyensdam_wrapper.F90)
target_link_libraries(censdam ensdam_static)

#KBadd_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/anam.c
#KB  COMMAND cython anam.pyx
#KB  MAIN_DEPENDENCY anam.pyx
#KB  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#KB)
#KBPython3_add_library(anam MODULE anam.c)
#KBtarget_link_libraries(anam PRIVATE censdam Python3::NumPy)

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

#kbadd_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/stat.c
#KB    COMMAND cython stat.pyx
#KB    MAIN_DEPENDENCY stat.pyx
#KB    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#KB)
#KBPython3_add_library(stat MODULE stat.c)
#KBtarget_link_libraries(stat PRIVATE censdam Python3::NumPy)