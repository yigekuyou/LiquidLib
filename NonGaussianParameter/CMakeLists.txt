set(TARGET NonGaussianParameter) 
add_library(${TARGET} SHARED src/${TARGET}.cpp include/${TARGET}.hpp)
add_executable( comput${TARGET}  src/${TARGET}_main.cpp include/${TARGET}_main.hpp )
target_sources(${TARGET}  PUBLIC include/${TARGET}.hpp)
target_sources(comput${TARGET}  PRIVATE include/${TARGET}_main.hpp)
target_link_libraries(${TARGET}
xdrfile
GSL::gsl
GSL::gslcblas
MeanSquaredDisplacement
Trajectory
)
install(TARGETS ${TARGET})
install(TARGETS comput${TARGET} )

