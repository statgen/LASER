execute_process(COMMAND bash -c "gfortran -print-search-dirs | sed -n -e 's/libraries: =//p' | sed -e 's/:/;/g'"  OUTPUT_VARIABLE GFORTRAN_PATHS)
string(STRIP "${GFORTRAN_PATHS}" GFORTRAN_PATHS)
foreach(PATH ${GFORTRAN_PATHS})
   find_library(GFORTRAN_LIB gfortran PATHS ${PATH})
   message(STATUS "${PATH} -- ${GFORTRAN_LIB}")
   if(GFORTRAN_LIB)
      break()
   endif()
endforeach()

find_library(OPENBLAS_LIB openblas HINT ${CGET_PREFIX}/lib)
find_library(GSL_LIB gsl HINT ${CGET_PREFIX}/lib)
find_package(Threads REQUIRED)
find_package(OpenMP REQUIRED)

if(CGET_PREFIX)
   include_directories(${CGET_PREFIX}/include)
endif()

# For CMake < 3.9, we need to make the target ourselves
if(NOT TARGET OpenMP::OpenMP_CXX)
    find_package(Threads REQUIRED)
    add_library(OpenMP::OpenMP_CXX IMPORTED INTERFACE)
    set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_COMPILE_OPTIONS ${OpenMP_CXX_FLAGS})
    # Only works if the same flag is passed to the linker; use CMake 3.9+ otherwise (Intel, AppleClang)
    set_property(TARGET OpenMP::OpenMP_CXX PROPERTY INTERFACE_LINK_LIBRARIES ${OpenMP_CXX_FLAGS} Threads::Threads)
endif()

message(STATUS "gfortran = ${GFORTRAN_LIB}")
message(STATUS "OpenBLAS = ${OPENBLAS_LIB}")
message(STATUS "GSL = ${GSL_LIB}")

add_executable(laser laser.cpp)
target_link_libraries(laser OpenMP::OpenMP_CXX ${GSL_LIB} ${OPENBLAS_LIB} ${GFORTRAN_LIB} Threads::Threads)

add_executable(trace trace.cpp)
target_link_libraries(trace OpenMP::OpenMP_CXX ${GSL_LIB} ${OPENBLAS_LIB} ${GFORTRAN_LIB} Threads::Threads)

install(TARGETS laser trace DESTINATION bin)
