/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *                    Christian Reinhardt
 *
 *
 */

#include "quantum3body_simulation_window.hh"

#include <complex>
// make sure <complex> is included before fftw3.h
// to have FFTW use the datatype defined by gcc's complex
// rather than it's own
#include <fftw3.h>
#include <omp.h>
#include <cassert>

#include <QtGui/QApplication>

int main(int argc, char* argv[])
{
    QApplication app(argc, argv);

    int ret(fftw_init_threads());
    assert(ret && "initializing fftw threads failed");
    fftw_plan_with_nthreads(omp_get_num_threads());

    Quantum3BodySimulationWindow window;
    window.show();

    return app.exec();
}

