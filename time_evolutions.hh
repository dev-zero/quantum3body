/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *                    Christian Reinhardt
 *
 *
 */

#ifndef TIME_EVOLUTIONS_HH
#define TIME_EVOLUTIONS_HH

#include "two_dim_spo.hh"

#include <limits>

/**
 * This is the "default" time evolution which is being used for
 * to simulated a particle in a harmonic or flat potential.
 */
struct DefaultTimeEvolution :
    public TimeEvolution
{
    /**
     * Set this to false if you want a flat zero potential
     * instead of the default harmonic potential.
     */
    bool useHarmonicPotential;

    DefaultTimeEvolution() :
        useHarmonicPotential(true)
    {
    }

    /**
     * Returns the value of the potential at a given (x,y).
     * This will not be called by the SPO code, only by the
     * evolve functions defined in this class. It has been
     * written explicitly to ease the understanding
     */
    inline double V(const double& x, const double& y) const
    {
        if (useHarmonicPotential)
            return 0.5*(x*x+y*y);
        else
            return 0.0;
    }

    /**
     * The time evolution in x-y-coordinates (aka "spatial evolution").
     * The factor 0.5 in the implementation stems from the SPO method.
     */
    complex x_y_evolve(const double& x, const double& y, const double& dt) const
    {
        return exp(-0.5 * V(x, y) * dt * I);
    }

    /**
     * The time evolution in kx-y-coordinates (aka "momentum evolution for x").
     */
    complex kx_y_evolve(const double& kx, const double& /* y */, const double& dt) const
    {
        return exp(-0.5 * kx*kx * dt * I);
    }

    /**
     * The time evolution in x-ky-coordinates (aka "momentum evolution for y").
     */
    complex x_ky_evolve(const double& /* x */, const double& ky, const double& dt) const
    {
        return exp(-0.5 * ky*ky * dt * I);
    }
};

/**
 * This is the quantum restricted 3-body time evolution.
 * For a theoretical understanding, please consult the separate theoretical
 * description.
 */

struct Quantum3BodyTimeEvolution :
    public TimeEvolution
{
    /**
     * This parameter is the epsilon in the formula for the potential of the
     * restricted 3-body problem. It is initially set to 0.2.
     * Setting the initial conditions to (0,0,0,2*pi) for (x,y,px,py) should
     * give a bound orbit around the primary body at (x,y)=(0.2,0).
     */
    double e;

    Quantum3BodyTimeEvolution() :
        e(0.2)
    {
    }

    /**
     * This is the same potential as for the restricted 3-body problem.
     */
    inline double V(const double& x, const double& y) const
    {
        double a(sqrt((x+e)*(x+e) + y*y));
        double b(sqrt((x-1.0+e)*(x-1.0+e) + y*y));
        if (fabs(a) < std::numeric_limits<double>::epsilon())
            a = std::numeric_limits<double>::epsilon();
        if (fabs(b) < std::numeric_limits<double>::epsilon())
            b = std::numeric_limits<double>::epsilon();
        return -(1.0-e)/a - e/b;
    }

    /**
     * The time evolution in kx-y-coordinates (aka "momentum evolution for x").
     */
    complex x_y_evolve(const double& x, const double& y, const double& dt) const
    {
        return exp(-0.5 * V(x, y) * dt * I);
    }

    /**
     * The time evolution in kx-y-coordinates (aka "momentum evolution for x").
     */
    complex kx_y_evolve(const double& kx, const double& y, const double& dt) const
    {
        return exp(-I*(0.5*kx*kx - y*kx)*dt);
    }

    /**
     * The time evolution in x-ky-coordinates (aka "momentum evolution for y").
     */
    complex x_ky_evolve(const double& x, const double& ky, const double& dt) const
    {
        return exp(-I*(0.5*ky*ky + x*ky)*dt);
    }
};

#endif // TIME_EVOLUTIONS_HH
