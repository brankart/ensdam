# Find target Python interpreter.
set(Python_FIND_STRATEGY LOCATION)
set(Python3_FIND_STRATEGY LOCATION)
find_package(Python3 COMPONENTS Interpreter Development NumPy)

add_subdirectory(src)

set(project_name "pyensdam")
set(setup_file "setup.py")

if(USE_MPI)
  set(project_name "${project_name}_mpi")
  set(setup_file "setup_mpi.py")
endif()

# Build wheel
add_custom_command(OUTPUT wheel.stamp POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/${project_name}" ${project_name}
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:scores>" ${project_name}/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:obserror>" ${project_name}/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:random>" ${project_name}/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:probability>" ${project_name}/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:statistics>" ${project_name}/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:anamorphosis>" ${project_name}/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:interpolation>" ${project_name}/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:augment>" ${project_name}/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:update>" ${project_name}/
                   COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/${setup_file}" bdist_wheel
                   VERBATIM)
add_custom_target(wheel ALL DEPENDS wheel.stamp)

# Let pip install the wheel
option(PYENSDAM_USER_INSTALL "Install with pip --user" ON)
if(PYENSDAM_USER_INSTALL)
  SET(EXTRA_PIP_ARGS "--user")
endif()
install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install --no-index --find-links=${CMAKE_CURRENT_BINARY_DIR}/dist ${EXTRA_PIP_ARGS} --force-reinstall --upgrade ${project_name} WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist)")

