/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *
 *
 *
 */

#include "quantum3body_simulation.hh"

#include <limits>
#include <iostream>

int main(int, char**)
{
    auto potential = [](const double& x, const double& y) { return 0.5*(x*x+y*y); };
    auto phi0 = [](const double& x, const double& y)->complex { return exp(-0.5*(x*x+y*y)); };

    Quantum3BodySimulation sim(64, 64, potential);
    sim.setInitial(phi0);

    for (size_t i(0); i < 1000; ++i)
        sim.evolveStep(0.005);

    return 0;
}

