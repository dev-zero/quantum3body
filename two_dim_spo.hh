/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *
 *
 *
 */

#ifndef TWO_DIM_SPO_HH
#define TWO_DIM_SPO_HH

#include <vector>
#include <functional>
#include <complex>
// make sure <complex> is included before fftw3.h
// to have FFTW use the datatype defined by gcc's complex
// rather than it's own
#include <fftw3.h>

typedef std::complex<double> complex;
const static complex I(complex(0,1));

struct TimeEvolution
{
    virtual complex x_y_evolve(const double& x, const double& y, const double& dt) const = 0;
    virtual complex kx_y_evolve(const double& kx, const double& y, const double& dt) const = 0;
    virtual complex x_ky_evolve(const double& x, const double& ky, const double& dt) const = 0;
};

class TwoDimSPO
{
public:
    TwoDimSPO(size_t sizeX, size_t sizeY);
    ~TwoDimSPO();
    void initialize(std::function<complex (const double&, const double&)> initialPhi);
    void evolveStep(const double& dt);
    void setTimeEvolution(const TimeEvolution* te);

    double binSizeX() const;
    double binSizeY() const;

    const complex* phi() const
    {
        return _phi;
    }

    size_t sizeX() const { return _sizeX; }
    size_t sizeY() const { return _sizeY; }

private:
    const size_t _sizeX, _sizeY;
    std::vector<double> _x, _y, _kx, _ky;
    complex* _phi;
    complex* _Phi;
    const TimeEvolution* _te;

    fftw_plan _fftPlanForwardX, _fftPlanBackwardX;
    fftw_plan _fftPlanForwardY, _fftPlanBackwardY;    
};

#endif // TWO_DIM_SPO_HH
