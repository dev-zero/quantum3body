/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *
 *
 *
 */

#include "quantum3body_simulation.hh"

#include <cassert>
#include <algorithm>

#include <iostream>
#include <iomanip>
#include <fstream>

Quantum3BodySimulation::Quantum3BodySimulation(size_t sizeX, size_t sizeY, PotentialFunction V) :
    _sizeX(sizeX),
    _sizeY(sizeY),
    _V(V),
    _x(sizeX),
    _y(sizeY),
    _kx(sizeX),
    _ky(sizeY),
    _spatial(sizeX*sizeY),
    _momentum(sizeX*sizeY)
{
    assert(!(sizeX % 2) && "sizeX must be even");
    assert(!(sizeY % 2) && "sizeY must be even");

    fftw_complex* spatial  = reinterpret_cast<fftw_complex*>(&_spatial.front());
    fftw_complex* momentum = reinterpret_cast<fftw_complex*>(&_momentum.front());

    // initialize the plan and determine the optimal method to calculate based on the size
    _fftPlanForward  = fftw_plan_dft_2d(static_cast<int>(_sizeX), static_cast<int>(_sizeY), spatial, momentum, FFTW_FORWARD, FFTW_ESTIMATE); // use FFTW_MEASURE once it gets real
    _fftPlanBackward = fftw_plan_dft_2d(static_cast<int>(_sizeX), static_cast<int>(_sizeY), momentum, spatial, FFTW_BACKWARD, FFTW_ESTIMATE); // use FFTW_MEASURE once it gets real


    // Trafo in X
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

        _fftPlanForwardX  = fftw_plan_many_dft(rank, n, howmany, spatial, inembed, istride, idist, momentum, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
        _fftPlanBackwardX = fftw_plan_many_dft(rank, n, howmany, momentum, inembed, istride, idist, spatial, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
    }
    // Trafo in Y
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

        _fftPlanForwardY  = fftw_plan_many_dft(rank, n, howmany, spatial, inembed, istride, idist, momentum, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
        _fftPlanBackwardY = fftw_plan_many_dft(rank, n, howmany, momentum, inembed, istride, idist, spatial, onembed, ostride, odist, FFTW_BACKWARD, FFTW_ESTIMATE);
    }
    
}

Quantum3BodySimulation::~Quantum3BodySimulation()
{
    fftw_destroy_plan(_fftPlanForward);
    fftw_destroy_plan(_fftPlanBackward);
    fftw_destroy_plan(_fftPlanForwardX);
    fftw_destroy_plan(_fftPlanBackwardX);
    fftw_destroy_plan(_fftPlanForwardY);
    fftw_destroy_plan(_fftPlanBackwardY);
}

void Quantum3BodySimulation::setInitial(std::function<complex (const double& x, const double& y)> initialPhi)
{
    const double delta_kx(2.0*M_PI*binSizeX());
    const double delta_ky(2.0*M_PI*binSizeY());
    const double delta_x(binSizeX());
    const double delta_y(binSizeY());

    // initialize x and k
    for (size_t i(0); i < _sizeX; ++i)
    {
        _x[i] = (i - _sizeX*0.5)*delta_x;
        // FFTW stores the frequencies "in-order":
        if (i < _sizeX/2)
           _kx[i] = i*delta_kx;
        else
           _kx[i] = -static_cast<double>(_sizeX - i)*delta_kx;
    }

    for (size_t i(0); i < _sizeY; ++i)
    {
        _y[i] = (i - _sizeY*0.5)*delta_y;
        // FFTW stores the frequencies "in-order":
        if (i < _sizeY/2)
           _ky[i] = i*delta_ky;
        else
           _ky[i] = -static_cast<double>(_sizeY - i)*delta_ky;
    }

    for (size_t i(0); i < _sizeX; ++i)
    {
        for (size_t j(0); j < _sizeY; ++j)
        {
            _spatial[j + (_sizeY*i)] = initialPhi(_x[i], _y[j]);
        }
    }
}

double Quantum3BodySimulation::binSizeX() const
{
    return 1.0/sqrt(static_cast<double>(_sizeX));
}

double Quantum3BodySimulation::binSizeY() const
{
    return 1.0/sqrt(static_cast<double>(_sizeY));
}

void Quantum3BodySimulation::evolveStep(const double& dt)
{
    auto spatialEvolve = [&](const complex& f, const double& x, const double& y)->complex {
        return f*exp(-0.5 * _V(x, y) * dt * complex(0, 1));
    };
    auto momentumEvolve = [&](const complex& f, const double& kx, const double& ky)->complex {
        return f*exp(-0.5 * (kx*kx+ky*ky) * dt * complex(0, 1));
    };
    auto momentumEvolveX = [&](const complex& f, const double& kx, const double&)->complex {
        return f*exp(-0.5 * kx*kx * dt * complex(0, 1));
    };
    auto momentumEvolveY = [&](const complex& f, const double&, const double& ky)->complex {
        return f*exp(-0.5 * ky*ky * dt * complex(0, 1));
    };
#if 1
    auto dump_spatial = [&](const std::string& identifier) {
        std::ofstream file("spatial.dump." + identifier);
//        std::cout << "Spatial Data:" << std::endl;
        for (size_t i(0); i < _sizeX; ++i)
        {
            for (size_t j(0); j < _sizeY; ++j)
            {
                file << _x[i] << " " << _y[j] << " " << abs(_spatial[j + _sizeY*i]) << std::endl;
            }
            file << std::endl;
        }
    };
#else
    auto dump_spatial = [&](const std::string&) {};
#endif

#if 0
    auto dump_momentum = [&]() {
        std::cout << "Momentum Data:" << std::endl;
        for (size_t i(0); i < _size/2; ++i)
        {
//            std::cout << std::fixed << std::setprecision(6) << "k[" << i << "]=" << _k[i] << ", F(k)=" << _momentum[i] << std::endl;
            std::cout << std::scientific
                << _k[_size/2+i] << ","
                << _momentum[_size/2+i].real() << ","
                << _momentum[_size/2+i].imag() << ","
                << std::endl;
        }
        for (size_t i(0); i < _size/2; ++i)
        {
//            std::cout << std::fixed << std::setprecision(6) << "k[" << i << "]=" << _k[i] << ", F(k)=" << _momentum[i] << std::endl;
            std::cout << std::scientific
                << _k[i] << ","
                << _momentum[i].real() << ","
                << _momentum[i].imag() << ","
                << std::endl;
        }
    };

#else
    auto dump_momentum = [&](const std::string&) {};
#endif

    dump_spatial("01");

//    std::cout << "spatial evolve" << std::endl;
    // f(x) *= exp(-0.5 * _V(x) * dt * i)
    for (size_t i(0); i < _sizeX; ++i)
    {
        for (size_t j(0); j < _sizeY; ++j)
        {
            spatialEvolve(_spatial[j + (_sizeY*i)], _x[i], _y[j]);
        }
    }

    dump_spatial("02");
    // draw spatial here (or below with the momentum)

//    std::cout << "forward fourier transformation" << std::endl;

#if 0
    fftw_execute(_fftPlanForward);
    for (auto& f: _momentum) { f /= sqrt(static_cast<double>(_sizeX*_sizeY)); } // renormalize, TODO: integrate in evolve function

    dump_momentum("01");
    // draw momentum (or both)

//    std::cout << "momentum evolve" << std::endl;
    // f(x) *= exp(-0.5 * k^2 * dt * i)
    for (size_t i(0); i < _sizeX; ++i)
    {
        for (size_t j(0); j < _sizeY; ++j)
        {
            momentumEvolve(_momentum[j + (_sizeY*i)], _kx[i], _ky[j]);
        }
    }

    dump_momentum("02");

//    std::cout << "backward fourier transformation" << std::endl;
    fftw_execute(_fftPlanBackward);
    for (auto& f: _spatial) { f /= sqrt(static_cast<double>(_sizeX*_sizeY)); } // renormalize, TODO: integrate in evolve function

#else
    fftw_execute(_fftPlanForwardX);
    for (auto& f: _momentum) { f /= sqrt(static_cast<double>(_sizeX)); }
    
    for (size_t i(0); i < _sizeX; ++i)
    {
        for (size_t j(0); j < _sizeY; ++j)
        {
            momentumEvolveX(_momentum[j + (_sizeY*i)], _kx[i], _y[j]);
        }
    }

    fftw_execute(_fftPlanBackwardX);
    for (auto& f: _spatial) { f /= sqrt(static_cast<double>(_sizeX)); }

    fftw_execute(_fftPlanForwardY);
    for (auto& f: _momentum) { f /= sqrt(static_cast<double>(_sizeY)); }

    for (size_t i(0); i < _sizeX; ++i)
    {
        for (size_t j(0); j < _sizeY; ++j)
        {
            momentumEvolveY(_momentum[j + (_sizeY*i)], _x[i], _ky[j]);
        }
    }
    fftw_execute(_fftPlanBackwardY);
    for (auto& f: _spatial) { f /= sqrt(static_cast<double>(_sizeY)); }

#endif

    dump_spatial("03");

//    std::cout << "spatial evolve" << std::endl;
    for (size_t i(0); i < _sizeX; ++i)
    {
        for (size_t j(0); j < _sizeY; ++j)
        {
            spatialEvolve(_spatial[j + (_sizeY*i)], _x[i], _y[j]);
        }
    }

    dump_spatial("04");
}
