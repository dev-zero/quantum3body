/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *
 *
 *
 */

#include "quantum3body_simulation_window.hh"
#include "quantum3body_simulation.hh"
#include "quantum_pixel_plot.hh"

#include "ui_quantum3body_simulation_window.h"

#include <QtCore/QDebug>
#include <QtCore/QTimer>

#include <algorithm>
#include <cassert>

const static size_t gridSizeX(256);
const static size_t gridSizeY(256);

enum PotentialSelection
{
    HARMONIC_POTENTIAL     = 0,
    ZERO_POTENTIAL         = 1,
    QUANTUM3BODY_POTENTIAL = 2
};
enum TimeEvolutionSelection
{
    DEFAULT_EVOLUTION      = 0,
    QUANTUM3BODY_EVOLUTION = 1
};

Quantum3BodySimulationWindow::Quantum3BodySimulationWindow(QWidget* p) :
    QMainWindow(p),
    _timer(new QTimer(this)),
    _ui(new Ui::Quantum3BodyWindow),
    _currentIteration(0)
{
    _ui->setupUi(this);

    _simulation = new Quantum3BodySimulation(gridSizeX, gridSizeY);

    _spatialPlot = new QuantumPixelPlot(this);
    _ui->spatialPlot->setModel(_spatialPlot);
    PixelDelegate* delegate(new PixelDelegate(this));
    _ui->spatialPlot->setItemDelegate(delegate);

    connect(_timer, SIGNAL(timeout()), SLOT(evolve()));
    connect(_ui->run, SIGNAL(toggled(bool)), SLOT(runSimulation(bool)));
    connect(_ui->pixelSize, SIGNAL(valueChanged(int)), delegate, SLOT(setPixelSize(int)));
    connect(_ui->reset, SIGNAL(pressed()), SLOT(resetSimulation()));

    resetSimulation();
}

Quantum3BodySimulationWindow::~Quantum3BodySimulationWindow()
{
    _timer->stop();
    delete _simulation;
}

void Quantum3BodySimulationWindow::resetSimulation()
{
//    auto phi0 = [](const double& x, const double& y)->complex { return exp(-0.5*(x*x+y*y)); };

    auto phi0 = [&](const double& x, const double& y)->complex {
        double kx(_ui->initialPropagationX->value());
        double ky(_ui->initialPropagationY->value());
        double dx(_ui->initialPositionX->value());
        double dy(_ui->initialPositionY->value());
        return exp(-0.5*((x-dx)*(x-dx)+(y-dy)*(y-dy)) - complex(0,1)*ky*y - complex(0,1)*kx*x);
    };

    Quantum3BodySimulation::PotentialFunction initialPotential;

    switch (_ui->initialPotential->currentIndex())
    {
        case HARMONIC_POTENTIAL:
            initialPotential = [](const double& x, const double& y) { return 0.5*(x*x+y*y); }; 
            break;
        case ZERO_POTENTIAL:
            initialPotential = [](const double&, const double&) { return 0.0; };
            break;
        case QUANTUM3BODY_POTENTIAL:
            initialPotential = [](const double& x, const double& y)->double {
                const double e(0.2);
                return -(1.0-e)/sqrt((x+e)*(x+e) + y*y) - e/sqrt((x-1.0+e)*(x-1.0+e) + y*y);
            }; 
            break;
        default:
            assert(!"Selected potential not defined.");
    }

    switch (_ui->initialTimeEvolution->currentIndex())
    {
        case DEFAULT_EVOLUTION:
            _simulation->useDefaultEvolution();
            break;
        case QUANTUM3BODY_EVOLUTION:
            _simulation->useQuantum3BodyEvolution();
            break;
        default:
            assert(!"Selected time evolution not defined.");
    }
    _simulation->setInitial(phi0);
    _simulation->setPotential(initialPotential);

    _currentIteration = 0;
    _ui->currentIteration->setValue(_currentIteration);

    plot();
}

void Quantum3BodySimulationWindow::runSimulation(bool run)
{
    if (run)
        _timer->start();
    else
        _timer->stop();
}

void Quantum3BodySimulationWindow::evolve()
{
    ++_currentIteration;
    _ui->currentIteration->setValue(_currentIteration);

    _simulation->evolveStep(_ui->dtParameter->value());

    if (_currentIteration % _ui->updatePlotSteps->value() == 0)
        plot();

    if (_currentIteration == _ui->numberOfIterations->value())
        _ui->run->setChecked(false);
}

void Quantum3BodySimulationWindow::plot()
{
    _spatialPlot->setSpatialData(_simulation->f_spatial(), gridSizeX, gridSizeY);
    _ui->spatialPlot->resizeColumnsToContents();
    _ui->spatialPlot->resizeRowsToContents();
}

