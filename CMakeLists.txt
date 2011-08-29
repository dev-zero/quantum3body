cmake_minimum_required(VERSION 2.8)

project(quantum3body)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wold-style-cast -Wredundant-decls -Wstrict-null-sentinel -Wmissing-noreturn -Woverloaded-virtual -Winit-self -Wunused -Wunused-function -Wshadow -Wwrite-strings -Wfloat-equal -Wconversion -pedantic -Wno-long-long -ggdb -std=c++0x -O2")

find_package(Qt4 REQUIRED)
set(QT_USE_QTGUI TRUE)
set(QT_USE_QTOPENGL TRUE)
set(QT_USE_QTWEBKIT TRUE)
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
    quantum_pixel_plot.hh
    )

add_executable(quantum3body-simulation
    ${quantum3body_simulation_MOC_SRCS} ${quantum3body_simulation_UIS_H}
    quantum3body_simulation.cc
    quantum3body_simulation_window.cc
    quantum_pixel_plot.cc
    main.cc
    )
target_link_libraries(quantum3body-simulation
    ${QT_LIBRARIES}
    -lfftw3
    )

add_executable(quantum3body_simulation_TEST
    quantum3body_simulation.cc
    quantum3body_simulation_TEST.cc
    )
target_link_libraries(quantum3body_simulation_TEST
    -lfftw3)

add_executable(fftw_separated_computation_TEST
    fftw_separated_computation_TEST.cc
    )
target_link_libraries(fftw_separated_computation_TEST
    -lfftw3)
