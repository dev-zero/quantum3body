/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *
 *
 *
 */

#ifndef QUANTUM3BODY_SIMULATION_HH
#define QUANTUM3BODY_SIMULATION_HH

#include <cstring>
#include <vector>
#include <functional>
#include <complex>
// make sure <complex> is included before fftw3.h
// to have FFTW use the datatype defined by gcc's complex
// rather than it's own
#include <fftw3.h>

typedef std::complex<double> complex;

class Quantum3BodySimulation
{
public:
    typedef std::function<double (const double&, const double&)> PotentialFunction;

    Quantum3BodySimulation(size_t xSize, size_t ySize, PotentialFunction V);
    ~Quantum3BodySimulation();

    void setInitial(std::function<complex (const double&, const double&)> initialPhi);
    
    double binSizeX() const;
    double binSizeY() const;

    void evolveStep(const double& dt);

    const std::vector<complex>& f_spatial() const { return _spatial; }
    const std::vector<complex>& f_momentum() const { return _momentum; }

    const std::vector<double>& x() const { return _x; }
    const std::vector<double>& y() const { return _y; }
    const std::vector<double>& kx() const { return _kx; }
    const std::vector<double>& ky() const { return _ky; }

private:
    size_t _sizeX, _sizeY;

    PotentialFunction _V;
    std::vector<double> _x, _y, _kx, _ky;
    std::vector<complex> _spatial, _momentum;
    fftw_plan _fftPlanForwardX, _fftPlanBackwardX;
    fftw_plan _fftPlanForwardY, _fftPlanBackwardY;
 
};

#endif // QUANTUM3BODY_SIMULATION_HH
