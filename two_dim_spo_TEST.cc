/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *                    Christian Reinhardt
 *
 *
 *
 */

#include "two_dim_spo.hh"
#include "time_evolutions.hh"

#include <iostream>
#include <cmath>
#include <cassert>
#include <complex>
// make sure <complex> is included before fftw3.h
// to have FFTW use the datatype defined by gcc's complex
// rather than it's own
#include <fftw3.h>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

int main(int, char**)
{
#ifdef HAVE_OPENMP
    int ret(fftw_init_threads());
    assert(ret && "initializing fftw threads failed");
    fftw_plan_with_nthreads(omp_get_num_threads());
#endif

    auto phi0 = [&](const double& x, const double& y)->complex {
        double kx(0.0);
        double ky(2.0*M_PI);
        double dx(0.0);
        double dy(0.0);
        return exp(-0.5*((x-dx)*(x-dx)+(y-dy)*(y-dy)) - complex(0,1)*ky*y - complex(0,1)*kx*x);
    };
    const double dt(0.005);
    const size_t numberOfSteps(100);

    TwoDimSPO sim(1024, 1024);
    Quantum3BodyTimeEvolution te;
    sim.setTimeEvolution(&te);
    sim.initialize(phi0);

    std::cout << "Please wait..." << std::endl;
    for (size_t i(0); i < numberOfSteps; ++i)
    {
        std::cout << ".";
        std::flush(std::cout);
        sim.evolveStep(dt);
    }
    std::cout << " done." << std::endl;

    return 0;
}

