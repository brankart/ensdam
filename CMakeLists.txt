cmake_minimum_required( VERSION 3.15 )

#list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

project( ensdam VERSION 1.0 LANGUAGES Fortran )

add_subdirectory(src EXCLUDE_FROM_ALL)
add_subdirectory(python EXCLUDE_FROM_ALL)

#add_custom_target(fortran DEPENDS ${PROJECT_NAME}_static ${PROJECT_NAME}_shared)
add_custom_target(fortran DEPENDS ${PROJECT_NAME}_static)
add_custom_target(python DEPENDS wheel )
add_dependencies(python ${PROJECT_NAME}_static)

#add_subdirectory(examples 2l_FROM_ALL)
