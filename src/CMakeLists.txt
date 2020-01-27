find_library(OPENBLAS_LIB openblas HINT ${CGET_PREFIX}/lib)
find_library(GSL_LIB gsl HINT ${CGET_PREFIX}/lib)
find_package(Threads REQUIRED)

if(CGET_PREFIX)
   include_directories(${CGET_PREFIX}/include)
endif()

message(STATUS "OpenBLAS = ${OPENBLAS_LIB}")
message(STATUS "GSL = ${GSL_LIB}")

add_executable(laser laser.cpp)
target_link_libraries(laser ${OPENBLAS_LIB} ${GSL_LIB} ${CMAKE_THREAD_LIBS_INIT} gfortran)

add_executable(trace trace.cpp)
target_link_libraries(trace ${OPENBLAS_LIB} ${GSL_LIB} ${CMAKE_THREAD_LIBS_INIT} gfortran)

install(TARGETS laser trace DESTINATION bin)
