set(TARGET Trajectory)
add_library(${TARGET} SHARED src/${TARGET}.cpp include/${TARGET}.hpp
)
set_target_properties(${TARGET}  PROPERTIES VERSION ${VERSION} SOVERSION ${SOV})

if(xdrfile_FOUND)
target_link_libraries(${TARGET}
xdrfile)
endif()

if(Boost_FOUND)
target_link_libraries(${TARGET}
Boost::math)
endif()
if(GSL_FOUND)
target_link_libraries(${TARGET}
GSL::gsl
GSL::gslcblas)
endif()
set_target_properties(${TARGET} PROPERTIES PUBLIC_HEADER include/${TARGET}.hpp)
install(TARGETS ${TARGET})
