set(TARGET ConvertTrajectoryType) 
add_library(${TARGET} SHARED src/${TARGET}.cpp include/${TARGET}.hpp)
add_executable( comput${TARGET}  src/${TARGET}_main.cpp include/${TARGET}_main.hpp )
set_target_properties(${TARGET}  PROPERTIES VERSION ${VERSION} SOVERSION ${SOV})
set_target_properties(${TARGET} PROPERTIES PUBLIC_HEADER include/${TARGET}.hpp)
include_directories(include)
target_link_libraries(${TARGET}
xdrfile
GSL::gsl
GSL::gslcblas
MeanSquaredDisplacement
Trajectory
)
target_link_libraries(comput${TARGET} ${TARGET} )
install(TARGETS ${TARGET})
install(TARGETS comput${TARGET} )

