project(Multi)

cmake_minimum_required (VERSION 2.6)

find_package(CGAL QUIET COMPONENTS Core )
  include( ${CGAL_USE_FILE} )

  include_directories (BEFORE "/include")
  include_directories (BEFORE "/dt")

find_package(MPI)
  include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  SET(CMAKE_C_COMPILER mpicc)
  SET(CMAKE_CXX_COMPILER mpicxx)
  SET(CXXFLAGS -lmetis)
  SET(MPIEXEC_MAX_NUMPROCS "4" CACHE STRING "Maximum number of processors available to run MPI applications.")

find_package(METIS)
  SET(CMAKE_C_COMPILER mpicc)
  SET(CMAKE_CXX_COMPILER mpicxx)
  SET(CXXFLAGS -lmetis) 

add_executable(example example.cpp)              
