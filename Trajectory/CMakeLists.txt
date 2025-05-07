set(TARGET Trajectory)
add_library(${TARGET} SHARED src/${TARGET}.cpp include/${TARGET}.hpp
)
set_target_properties(${TARGET}  PROPERTIES VERSION ${VERSION} SOVERSION ${SOV})
target_link_libraries(${TARGET}
xdrfile
GSL::gsl
GSL::gslcblas
)
set_target_properties(${TARGET} PROPERTIES PUBLIC_HEADER include/${TARGET}.hpp)
install(TARGETS ${TARGET})
