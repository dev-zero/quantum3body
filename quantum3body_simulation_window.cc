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
#include <QtGui/QFileDialog>
#include <QtGui/QImageWriter>
#include <QtCore/QFile>
#include <QtCore/QDateTime>

#include <algorithm>
#include <cassert>

const static size_t gridSizeX(1024);
const static size_t gridSizeY(1024);

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
    _currentIteration(0),
    _lastResetTimestamp(0)
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
    connect(_ui->browsePictureFolder, SIGNAL(pressed()), SLOT(browsePictureFolder()));

    QByteArray format;
    foreach (format, QImageWriter::supportedImageFormats())
        _ui->pictureFormat->insertItem(0, format);
    _ui->pictureFormat->setCurrentIndex(_ui->pictureFormat->findText("png")); // set the selected item to png or nothing

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
            initialPotential = [](const double& x, const double& y)->double { return 0.5*(x*x+y*y); }; 
            break;
        case ZERO_POTENTIAL:
            initialPotential = [](const double&, const double&)->double { return 0.0; };
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

    _lastResetTimestamp = QDateTime::currentDateTime().toTime_t();

    plot();
}

void Quantum3BodySimulationWindow::browsePictureFolder()
{
    _ui->pictureFolder->setText(QFileDialog::getExistingDirectory(this,
            tr("Open base directory for saving pictures"),
            QString()));
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

    double totalProbability(0.0);
    const double binSizeX(_simulation->binSizeX()), binSizeY(_simulation->binSizeY());
    const std::vector<complex>& f(_simulation->f_spatial());
    for (size_t i(0); i < gridSizeX; ++i)
    {
        for (size_t j(0); j < gridSizeY; ++j)
        {
            totalProbability += binSizeX*binSizeY*fabs(f[j + gridSizeY*i]);
        }
    }
    _ui->totalProbability->setValue(totalProbability);

    if (_ui->savePictures->isChecked())
    {
        QDir pictureFolder(_ui->pictureFolder->text());
        if (pictureFolder.exists())
        {
            QString filepath(
                    pictureFolder.filePath(
                        QString("%1.%2_%3.%4")
                        .arg(_ui->pictureBasename->text())
                        .arg(_lastResetTimestamp)
                        .arg(_currentIteration, 4, 10, QLatin1Char('0'))
                        .arg(_ui->pictureFormat->currentText().toLower())));

            QImage image(gridSizeY, gridSizeX, QImage::Format_RGB32);
            for(size_t i(0); i < gridSizeX; ++i)
            {
                for(size_t j(0); j < gridSizeY; ++j)
                {
                    QModelIndex idx(_spatialPlot->index(j, i));
                    image.setPixel(i, j, _spatialPlot->data(idx).value<QRgb>());
                }
            }
            qDebug() << "dumping picture to " << filepath;
            image.save(filepath, _ui->pictureFormat->currentText().toAscii().constData());
        }
    }
}
