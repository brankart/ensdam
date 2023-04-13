# compile object files and create static and shared libraries

set (EnsAnam
     EnsAnam/anaobs.f90
     EnsAnam/anaqua.f90
     EnsAnam/anatra.f90
     EnsAnam/anautil.f90
    )

set (EnsAugm
     EnsAugm/ensaugm.F90
     EnsAugm/schurprod.f90
    )

set (EnsScores
     EnsScores/score_crps.F90
     EnsScores/score_entropy.F90
     EnsScores/score_optimality.F90
     EnsScores/score_ranks.F90
     EnsScores/score_rcrv.F90
    )

set (EnsStat
     EnsStat/covariance.F90
     EnsStat/meanstd.F90
    )

set (EnsUpdate
     EnsUpdate/mcmc_update.F90
    )

set (InterpTools
     InterpTools/interp.f90
    )

set (ObsError
     ObsError/obserror.F90
    )

set (StochTools
     StochTools/constant.f90
     StochTools/qxgs.f90
     StochTools/stoanam.f90
     StochTools/stogprod.f90
     StochTools/storfg.F90
     StochTools/storng.f90
     StochTools/stotge.f90
     StochTools/stoutil.f90
    )

set (TranSpHO
     TranSpHO/spharea.f90
     TranSpHO/sphylm.F90
    )

set (SRC_ENSDAM
     ${EnsAnam}
     ${EnsAugm}
     ${EnsScores}
     ${EnsStat}
     ${EnsUpdate}
     ${InterpTools}
     ${ObsError}
     ${StochTools}
     ${TranSpHO}
    )

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

# Use a common object library for building shared and static targets
add_library( ${PROJECT_NAME}_objects OBJECT ${SRC_ENSDAM} )
target_compile_definitions( ${PROJECT_NAME}_objects PUBLIC ${PUBLIC_FLAGS}
                                                    PRIVATE ${PRIVATE_FLAGS} )
target_include_directories(${PROJECT_NAME}_objects
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
)

add_library( ${PROJECT_NAME}_static STATIC $<TARGET_OBJECTS:${PROJECT_NAME}_objects> )
target_include_directories(${PROJECT_NAME}_static
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
)
set_target_properties(${PROJECT_NAME}_static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})

#add_library( ${PROJECT_NAME}_shared SHARED $<TARGET_OBJECTS:${PROJECT_NAME}_objects> )
#target_include_directories(${PROJECT_NAME}_shared
#PUBLIC
#$<INSTALL_INTERFACE:include>
#$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
#)
#set_target_properties(${PROJECT_NAME}_shared PROPERTIES OUTPUT_NAME ${PROJECT_NAME})

#install(TARGETS ${PROJECT_NAME}_static ${PROJECT_NAME}_shared
install(TARGETS ${PROJECT_NAME}_static
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
)
install(DIRECTORY ${CMAKE_BINARY_DIR}/modules/
        DESTINATION include
)

