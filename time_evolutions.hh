/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *
 *
 *
 */

#ifndef TIME_EVOLUTIONS_HH
#define TIME_EVOLUTIONS_HH

#include "two_dim_spo.hh"

struct DefaultTimeEvolution :
    public TimeEvolution
{
    bool useHarmonicPotential;

    DefaultTimeEvolution() :
        useHarmonicPotential(true)
    {
    }
    inline double V(const double& x, const double& y) const
    {
        if (useHarmonicPotential)
            return 0.5*(x*x+y*y);
        else
            return 0.0;
    }
    complex x_y_evolve(const double& x, const double& y, const double& dt) const
    {
        return exp(-0.5 * V(x, y) * dt * I);
    }
    complex kx_y_evolve(const double& kx, const double& /* y */, const double& dt) const
    {
        return exp(-0.5 * kx*kx * dt * I);
    }
    complex x_ky_evolve(const double& /* x */, const double& ky, const double& dt) const
    {
        return exp(-0.5 * ky*ky * dt * I);
    }
};

struct Quantum3BodyTimeEvolution :
    public TimeEvolution
{
    double e;

    Quantum3BodyTimeEvolution() :
        e(0.2)
    {
    }
    inline double V(const double& x, const double& y) const
    {
        return -(1.0-e)/sqrt((x+e)*(x+e) + y*y) - e/sqrt((x-1.0+e)*(x-1.0+e) + y*y);
    }
    complex x_y_evolve(const double& x, const double& y, const double& dt) const
    {
        return exp(-0.5 * V(x, y) * dt * I);
    }
    complex kx_y_evolve(const double& kx, const double& y, const double& dt) const
    {
        return exp(-I*(0.5*kx*kx - y*kx)*dt);
    }
    complex x_ky_evolve(const double& x, const double& ky, const double& dt) const
    {
        return exp(-I*(0.5*ky*ky + x*ky)*dt);
    }
};

#endif // TIME_EVOLUTIONS_HH
