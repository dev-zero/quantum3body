/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *                    Christian Reinhardt
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

/**
 * This is the base interface for a TimeEvolution.
 * Since the split of the operator can not be done at runtime, the
 * user of this class has to do it manually and provide the implementations.
 */
struct TimeEvolution
{
    /**
     * The time evolution in x-y-coordinates (aka "spatial evolution").
     */
    virtual complex x_y_evolve(const double& x, const double& y, const double& dt) const = 0;
    /**
     * The time evolution in kx-y-coordinates (aka "momentum evolution for x").
     */
    virtual complex kx_y_evolve(const double& kx, const double& y, const double& dt) const = 0;
    /**
     * The time evolution in x-ky-coordinates (aka "momentum evolution for y").
     */
    virtual complex x_ky_evolve(const double& x, const double& ky, const double& dt) const = 0;
};

/**
 * This is the driver for doing simulations using the Split-Operator (SPO) Method.
 * The dimensions have to be defined at the time of creation. The reason is that the libraries
 * used have to allocate and initialize storage.
 */

class TwoDimSPO
{
public:
    /**
     * Construct a new object.
     * \param sizeX the number of bins in X-direction
     * \param sizeY the number of bins in Y-direction
     */
    TwoDimSPO(size_t sizeX, size_t sizeY);

    ~TwoDimSPO();

    /**
     * Set phi0 by providing a function which is being called for each (x,y).
     * \param initialPhi This function should return the initial values for phi0,
     *                   phi0(x_i,x_j) to be precise.
     */
    void initialize(std::function<complex (const double&, const double&)> initialPhi);

    /**
     * Do one step in the simulation.
     * \param dt The time delta being used to calculate the next iteration
     */
    void evolveStep(const double& dt);

    /**
     * The time evolution functions can be exchanged at runtime, please
     * read the documentation for the TimeEvolution struct to know what
     * you have to provide
     * \param te Provide a TimeEvolution object. Ownership remains at your side.
     */
    void setTimeEvolution(const TimeEvolution* te);

    /**
     * Returns dx used in the calculations (this depends on the sizeX).
     */
    double binSizeX() const;

    /**
     * Returns dy used in the calculations (this depends on the sizeY).
     */
    double binSizeY() const;

    /**
     * Returns the current discretized phi.
     * \return This returns the current phi in row-major format. This means
     *         that the the value of phi at (x_i, y_j) can be found at
     *         array[j + sizeY*i].
     */
    const complex* phi() const
    {
        return _phi;
    }

    /**
     * Returns the same size as given at initialization.
     */
    size_t sizeX() const { return _sizeX; }

    /**
     * Returns the same size as given at initialization.
     */
    size_t sizeY() const { return _sizeY; }

private:
    const size_t _sizeX, _sizeY;
    std::vector<double> _x, _y, _kx, _ky;

    complex* _phi;
    complex* _Phi;
    complex* _spatialEvolutionValues;

    const TimeEvolution* _te;

    fftw_plan _fftPlanForwardX, _fftPlanBackwardX;
    fftw_plan _fftPlanForwardY, _fftPlanBackwardY;    

    void evolve_x_y_first(const double& dt);
    void evolve_x_y_second(const double& dt);
    void evolve_kx_y(const double& dt);
    void evolve_x_ky(const double& dt);

};

#endif // TWO_DIM_SPO_HH
