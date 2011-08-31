/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *                    Christian Reinhardt
 *
 *
 */

#include "two_dim_spo.hh"

#include <cassert>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

TwoDimSPO::TwoDimSPO(size_t sizeX, size_t sizeY) :
    _sizeX(sizeX),
    _sizeY(sizeY),
    _x(sizeX),
    _y(sizeY),
    _kx(sizeX),
    _ky(sizeY)
{
    assert(!(sizeX % 2) && "sizeX must be even");
    assert(!(sizeY % 2) && "sizeY must be even");

    // allocate the arrays using fftw_malloc to ensure the array are on memory boundaries suited for using SIMD
    _phi = static_cast<complex*>(fftw_malloc(sizeof(fftw_complex) * sizeX * sizeY));
    _Phi = static_cast<complex*>(fftw_malloc(sizeof(fftw_complex) * sizeX * sizeY));

    fftw_complex* phi = reinterpret_cast<fftw_complex*>(_phi);
    fftw_complex* Phi = reinterpret_cast<fftw_complex*>(_Phi);

    int fftw_method(FFTW_MEASURE); // eventually use FFTW_PATIENCE to get faster

    // Initialize the transformation in X
    {
        int rank = 1; // 1-Dimensional
        int n[] = { static_cast<int>(sizeX)}; // number of elements in each rank, 1 rank => one element in this array, value = sizeX
        int howmany = static_cast<int>(sizeY); // how many arrays with constant y are there? => sizeY
        int* inembed = NULL; // if the input array are embedded in larger rank arrays (???)
        int istride = static_cast<int>(sizeY); // the values of x are located at 1*i_y + sizeY*i_x, this is therefore sizeY
        int idist = 1; // the distance between the first element of the first array and the first element of the second array (1)

        // the output array is organized the same way
        int odist = idist;
        int ostride = istride;
        int* onembed = NULL;

        _fftPlanForwardX  = fftw_plan_many_dft(rank, n, howmany, phi, inembed, istride, idist, Phi, onembed, ostride, odist, FFTW_FORWARD, fftw_method);
        _fftPlanBackwardX = fftw_plan_many_dft(rank, n, howmany, Phi, inembed, istride, idist, phi, onembed, ostride, odist, FFTW_BACKWARD, fftw_method);
    }
    // initialize the transformation in Y
    {
        int rank = 1; // 1-Dimensional
        int n[] = {static_cast<int>(sizeY)}; // number of elements in each rank, 1 rank => one element in this array, value = sizeY
        int howmany = static_cast<int>(sizeX); // how many arrays with constant x are there? => sizeX
        int* inembed = NULL; // if the input array are embedded in larger rank arrays (???)
        int istride = 1; // the values of y are located at 1*i_y + sizeY*i_x, this is therefore 1
        int idist = static_cast<int>(sizeY); // the distance between the first element of the first array and the first element of the second array (sizeY)

        // the output array is organized the same way
        int odist = idist;
        int ostride = istride;
        int* onembed = NULL;

        _fftPlanForwardY  = fftw_plan_many_dft(rank, n, howmany, phi, inembed, istride, idist, Phi, onembed, ostride, odist, FFTW_FORWARD, fftw_method);
        _fftPlanBackwardY = fftw_plan_many_dft(rank, n, howmany, Phi, inembed, istride, idist, phi, onembed, ostride, odist, FFTW_BACKWARD, fftw_method);
    }


    const double delta_kx(2.0*M_PI*binSizeX());
    const double delta_ky(2.0*M_PI*binSizeY());
    const double delta_x(binSizeX());
    const double delta_y(binSizeY());

    // initialize x and kx
    for (size_t i(0); i < _sizeX; ++i)
    {
        _x[i] = (i - _sizeX*0.5)*delta_x;
        // FFTW stores the frequencies "in-order":
        if (i < _sizeX/2)
            _kx[i] = i*delta_kx;
        else
            _kx[i] = -static_cast<double>(_sizeX - i)*delta_kx;
    }

    // initialize y and ky
    for (size_t i(0); i < _sizeY; ++i)
    {
        _y[i] = (i - _sizeY*0.5)*delta_y;
        // FFTW stores the frequencies "in-order":
        if (i < _sizeY/2)
            _ky[i] = i*delta_ky;
        else
            _ky[i] = -static_cast<double>(_sizeY - i)*delta_ky;
    }

}

TwoDimSPO::~TwoDimSPO()
{
    fftw_free(_phi);
    fftw_free(_Phi);
}

void TwoDimSPO::initialize(std::function<complex (const double&, const double&)> initialPhi)
{
    for (size_t i(0); i < _sizeX; ++i)
    {
        for (size_t j(0); j < _sizeY; ++j)
        {
            _phi[j + _sizeY*i] = initialPhi(_x[i], _y[j]);
        }
    }
}

void TwoDimSPO::evolveStep(const double& dt)
{
    /**
     * The basic scheme here is:
     * - (half) time evolution in x-y
     * - forward FFT for x
     * - time evolution in kx-y
     * - backward FFT for x
     * - forward FFT for y
     * - time evolution in x-ky
     * - backward FFT for y
     * - (half) time evolution in x-y
     *
     * Please note that the FFTW library does not normalize
     * the values, requiring us to do it. The factor is 1/sqrt(sizeX)
     * after transformations in x, 1/sqrt(sizeY) resp. for transformations
     * in y. This normalization is mangled in the time evolution algorithm
     * (but quiet obvious) to reduce the number of loops.
     */

#pragma omp parallel for
    for (size_t i = 0; i < _sizeX; ++i)
    {
        for (size_t j(0); j < _sizeY; ++j)
        {
            _phi[j + _sizeY*i] *= _te->x_y_evolve(_x[i], _y[j], dt);
        }
    }

    fftw_execute(_fftPlanForwardX);

#pragma omp parallel for
    for (size_t i = 0; i < _sizeX; ++i)
    {
        for (size_t j(0); j < _sizeY; ++j)
        {
            _Phi[j + _sizeY*i] /= sqrt(static_cast<double>(_sizeX));
            _Phi[j + _sizeY*i] *= _te->kx_y_evolve(_kx[i], _y[j], dt);
        }
    }

    fftw_execute(_fftPlanBackwardX);

#pragma omp parallel for
    for (size_t i = 0; i < _sizeX*_sizeY; ++i) { _phi[i] /= sqrt(static_cast<double>(_sizeX)); }

    fftw_execute(_fftPlanForwardY);

    for (size_t i = 0; i < _sizeX; ++i)
    {
        for (size_t j(0); j < _sizeY; ++j)
        {
            _Phi[j + _sizeY*i] /= sqrt(static_cast<double>(_sizeY));
            _Phi[j + _sizeY*i] *= _te->x_ky_evolve(_x[i], _ky[j], dt);
        }
    }

    fftw_execute(_fftPlanBackwardY);

#pragma omp parallel for
    for (size_t i = 0; i < _sizeX; ++i)
    {
        for (size_t j(0); j < _sizeY; ++j)
        {
            _phi[j + _sizeY*i] /= sqrt(static_cast<double>(_sizeY));
            _phi[j + _sizeY*i] *= _te->x_y_evolve(_x[i], _y[j], dt);
        }
    }
}

double TwoDimSPO::binSizeX() const
{
    return 1.0/sqrt(static_cast<double>(_sizeX));
}

double TwoDimSPO::binSizeY() const
{
    return 1.0/sqrt(static_cast<double>(_sizeY));
}

void TwoDimSPO::setTimeEvolution(const TimeEvolution* te)
{
    _te = te;
}
