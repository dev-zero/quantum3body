cmake_minimum_required(VERSION 2.8)

project(quantum3body)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wold-style-cast -Wredundant-decls -Wstrict-null-sentinel -Wmissing-noreturn -Woverloaded-virtual -Winit-self -Wunused -Wunused-function -Wshadow -Wwrite-strings -Wfloat-equal -Wconversion -pedantic -Wno-long-long")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2 -mtune=native")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb -pg")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_OPENMP -fopenmp") # enable OpenMP

find_package(Qt4 REQUIRED)
set(QT_USE_QTGUI TRUE)
include(${QT_USE_FILE})

# mark all 3rdparty headers as system headers to disable warnings generated
# for them by the compiler
include_directories(SYSTEM ${QT_INCLUDES})

# append the cmake build and source directories to include genereated headers
include_directories(${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

QT4_WRAP_UI(quantum3body_simulation_UIS_H
    quantum3body_simulation_window.ui)

QT4_WRAP_CPP(quantum3body_simulation_MOC_SRCS
    quantum3body_simulation_window.hh
    )

add_executable(quantum3body-simulation
    ${quantum3body_simulation_MOC_SRCS} ${quantum3body_simulation_UIS_H}
    two_dim_spo.cc
    quantum3body_simulation_window.cc
    color_maps.cc
    main.cc
    )

# remove -lfftw3_omp if you don't have OpenMP or replace
# it by -lfftw3_threads if you use FFTW < 3.3 and have FFTW compiled with OpenMP
target_link_libraries(quantum3body-simulation
    ${QT_LIBRARIES}
    -lfftw3_omp -lfftw3 -lm
    )

add_executable(two_dim_spo_TEST
    two_dim_spo.cc
    two_dim_spo_TEST.cc
    )
target_link_libraries(two_dim_spo_TEST
    -lfftw3_omp -lfftw3)

