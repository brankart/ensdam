# Find target Python interpreter.
set(Python_FIND_STRATEGY LOCATION)
set(Python3_FIND_STRATEGY LOCATION)
find_package(Python3 COMPONENTS Interpreter Development NumPy)

add_subdirectory(src)

# Build wheel
add_custom_command(OUTPUT wheel.stamp POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory "${CMAKE_CURRENT_SOURCE_DIR}/pyensdam" pyensdam
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:scores>" pyensdam/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:obserror>" pyensdam/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:random>" pyensdam/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:probability>" pyensdam/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:statistics>" pyensdam/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:anamorphosis>" pyensdam/
                   COMMAND ${CMAKE_COMMAND} -E copy "$<TARGET_FILE:interpolation>" pyensdam/
                   COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/setup.py" bdist_wheel
                   VERBATIM)
add_custom_target(wheel ALL DEPENDS wheel.stamp)

# Let pip install the wheel
option(PYENSDAM_USER_INSTALL "Install with pip --user" ON)
if(PYENSDAM_USER_INSTALL)
  SET(EXTRA_PIP_ARGS "--user")
endif()
install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install --no-index --find-links=${CMAKE_CURRENT_BINARY_DIR}/dist ${EXTRA_PIP_ARGS} --force-reinstall --upgrade pyensdam WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/dist)")

