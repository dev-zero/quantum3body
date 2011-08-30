/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *
 *
 *
 */

#include <iostream>
#include <limits>
#include <complex>
// make sure <complex> is included before fftw3.h
// to have FFTW use the datatype defined by gcc's complex
// rather than it's own
#include <fftw3.h>
#include <vector>

typedef std::complex<double> complex;

int main(int, char**)
{
    fftw_plan planForward;
    fftw_plan planForwardX, planForwardY;
    
    auto initialPhi = [](const double& x, const double& y)->complex { return exp(-0.5*(x*x+y*y)); };
    const size_t sizeX(256), sizeY(256);

    std::vector<double> x(sizeX), y(sizeY); 
    std::vector<complex> spatial(sizeX*sizeY), momentum(sizeX*sizeY);
    std::vector<complex> momentumSeparatedX(sizeX*sizeY), momentumSeparatedY(sizeX*sizeY);
    std::vector<complex> spatialAfterX(sizeX*sizeY);

    fftw_complex* fftw_spatial            = reinterpret_cast<fftw_complex*>(&spatial.front()); // common to all transformations

    fftw_complex* fftw_momentum           = reinterpret_cast<fftw_complex*>(&momentum.front());
    fftw_complex* fftw_momentumSeparatedX = reinterpret_cast<fftw_complex*>(&momentumSeparatedX.front());
    fftw_complex* fftw_momentumSeparatedY = reinterpret_cast<fftw_complex*>(&momentumSeparatedY.front());

    // initialize the plan and determine the optimal method to calculate based on the size
    planForward = fftw_plan_dft_2d(sizeX, sizeY, fftw_spatial, fftw_momentum, FFTW_FORWARD, FFTW_ESTIMATE);

    // Forward Trafo in X
    {
        int rank = 1; // 1-Dimensional
        int n[] = {sizeX}; // number of elements in each rank, 1 rank => one element in this array, value = sizeX
        int howmany = sizeY; // how many arrays with constant y are there? => sizeY
        int* inembed = NULL; // if the input array are embedded in larger rank arrays (???)
        int istride = sizeY; // the values of x are located at 1*i_y + sizeY*i_x, this is therefore sizeY
        int idist = 1; // the distance between the first element of the first array and the first element of the second array (1)

        // the output array is organized the same way
        int odist = idist;
        int ostride = istride;
        int* onembed = NULL;

        planForwardX = fftw_plan_many_dft(rank, n, howmany, fftw_spatial, inembed, istride, idist, fftw_momentumSeparatedX, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    // Forward Trafo in Y
    {
        int rank = 1; // 1-Dimensional
        int n[] = {sizeY}; // number of elements in each rank, 1 rank => one element in this array, value = sizeY
        int howmany = sizeX; // how many arrays with constant x are there? => sizeX
        int* inembed = NULL; // if the input array are embedded in larger rank arrays (???)
        int istride = 1; // the values of y are located at 1*i_y + sizeY*i_x, this is therefore 1
        int idist = sizeY; // the distance between the first element of the first array and the first element of the second array (sizeY)

        // the output array is organized the same way
        int odist = idist;
        int ostride = istride;
        int* onembed = NULL;

        planForwardY = fftw_plan_many_dft(rank, n, howmany, fftw_momentumSeparatedX, inembed, istride, idist, fftw_momentumSeparatedY, onembed, ostride, odist, FFTW_FORWARD, FFTW_ESTIMATE);
    }
     
    const double delta_x(1.0/sqrt(static_cast<double>(sizeX)));
    const double delta_y(1.0/sqrt(static_cast<double>(sizeY)));    

    // initialize x and k
    for (size_t i(0); i < sizeX; ++i)
        x[i] = (i - sizeX*0.5)*delta_x;

    for (size_t i(0); i < sizeY; ++i)
        y[i] = (i - sizeY*0.5)*delta_y;

    for (size_t i(0); i < sizeX; ++i)
    {
        for (size_t j(0); j < sizeY; ++j)
        {
            spatial[j + (sizeY*i)] = initialPhi(x[i], y[j]);
        }
    }
    
    fftw_execute(planForward);
    for (auto& f: momentum) { f /= sqrt(static_cast<double>(sizeX*sizeY)); }

    fftw_execute(planForwardX);
    for (auto& f: momentumSeparatedX) { f /= sqrt(static_cast<double>(sizeX)); }

    fftw_execute(planForwardY);
    for (auto& f: momentumSeparatedY) { f /= sqrt(static_cast<double>(sizeY)); }
    
    size_t differences(0);
    for (size_t i(0); i < sizeX*sizeY; ++i)
    {
        double abs_diff(fabs(momentumSeparatedY[i] - momentum[i]));
        if (abs_diff > std::numeric_limits<double>::epsilon())
        {
            std::cout << "difference found: " << abs_diff << " @" << i << std::endl;
            ++differences;
        }
    }
    std::cout << "number of differences: " << differences << std::endl;

    fftw_destroy_plan(planForward);
    fftw_destroy_plan(planForwardX);
    fftw_destroy_plan(planForwardY);
    return 0;
}
