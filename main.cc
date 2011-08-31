/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *                    Christian Reinhardt
 *
 *
 */

/** \mainpage
 * The purpose of this application is to simulate the quantum restricted
 * 3-body problem. This means that for the Hamiltonian in the time evolution
 * the Hamiltonian from the restricted 3-body problem is being used, transformed
 * to quantum mechanics. To be able to effectively time evolve the state,
 * the Split-Operator-Method has been used.
 * 
 * Please contact the authors for an in-depth theoretical description.
 * 
 * \section requirements Requirements
 * 
 * The application is written in C++, using constructs and syntax from the
 * new C++11 standard. You therefore need a recent compiler supporting the
 * new standard. The application has been tested using gcc-4.6.1 with
 * OpenMP support.
 * 
 * Furthermore you need:
 * - qt 4.7.1 or newer
 * - fftw 3.3.1_beta1 or newer, built with OpenMP support
 * - cmake 2.8 or newer
 * 
 * You could use an older fftw and/or one without OpenMP-support you have to
 * adjust the main.cc and the build flags.
 * 
 * \section build Build instructions
 * \code
 *   cd build/
 *   cmake ..
 *   make
 * \endcode
 * The executables are placed in the build/ directory:
 * - quantum3body-simulation .. this is the main application
 * - quantum3body_simulation_TEST.cc .. a simple test application for the initial SPO code
 * - fftw_separated_computation_TEST .. a test application to compare the 2x1D with the 2D FFT
 * 
 * \section running Running the application
 * 
 * Start the main application using:
 * \code
 *   ./quantum3body-simulation
 * \endcode
 * 
 * It should be mostly self-explanetary.
 * If you want to use multiple threads for computation, set the environment
 * variable OMP_NUM_THREADS to the number of desired threads before.
 *
 * For example using bash:
 * \code
 *   export OMP_NUM_THREADS=8
 *   ./quantum3body-simulation
 * \endcode
 *
 */

#include "quantum3body_simulation_window.hh"

#include <cassert>
#include <complex>
// make sure <complex> is included before fftw3.h
// to have FFTW use the datatype defined by gcc's complex
// rather than it's own
#include <fftw3.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <QtGui/QApplication>

int main(int argc, char* argv[])
{
    QApplication app(argc, argv);

#ifdef HAVE_OPENMP
    int ret(fftw_init_threads());
    assert(ret && "initializing fftw threads failed");
    fftw_plan_with_nthreads(omp_get_num_threads());
#endif

    Quantum3BodySimulationWindow window;
    window.show();

    return app.exec();
}

