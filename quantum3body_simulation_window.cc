/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *                    Christian Reinhardt
 *
 *
 *
 */

#include "quantum3body_simulation_window.hh"
#include "two_dim_spo.hh"
#include "time_evolutions.hh"
#include "color_maps.hh"

#include "ui_quantum3body_simulation_window.h"

#include <QtCore/QDebug>
#include <QtCore/QTimer>
#include <QtGui/QFileDialog>
#include <QtGui/QImageWriter>
#include <QtCore/QFile>
#include <QtCore/QDateTime>

#include <algorithm>
#include <cassert>

/**
 * This should match the entries in the initialPotential QComboBox
 */
enum PotentialSelection
{
    HARMONIC_POTENTIAL     = 0,
    ZERO_POTENTIAL         = 1
};

/**
 * This should match the entries in the initialTimeEvolution QComboBox
 */
enum TimeEvolutionSelection
{
    DEFAULT_EVOLUTION      = 0,
    QUANTUM3BODY_EVOLUTION = 1
};

/**
 * This should match the entries in the resolution QComboBox
 */
enum Resolution
{
    RESOLUTION_64x64     = 0,
    RESOLUTION_128x128   = 1,
    RESOLUTION_256x256   = 2,
    RESOLUTION_512x512   = 3,
    RESOLUTION_1024x1024 = 4,
    RESOLUTION_2048x2048 = 5,
    RESOLUTION_4096x4096 = 6,
    RESOLUTION_8192x8192 = 7
};

Quantum3BodySimulationWindow::Quantum3BodySimulationWindow(QWidget* p) :
    QMainWindow(p),
    _simulation(nullptr),
    _timer(new QTimer(this)),
    _ui(new Ui::Quantum3BodyWindow),
    _currentIteration(0),
    _lastResetTimestamp(0),
    _image(nullptr)
{
    _ui->setupUi(this);

    _imageLabel = new QLabel;
    _imageLabel->setBackgroundRole(QPalette::Base);
    _imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    _imageLabel->setScaledContents(true);
    
    _ui->plotScrollArea->setBackgroundRole(QPalette::Dark);
    _ui->plotScrollArea->setWidget(_imageLabel);

    _ui->fitToWindowAction->setCheckable(true);
    _ui->fitToWindowAction->setChecked(true);
    updateViewActions();

     connect(_ui->fitToWindowAction, SIGNAL(triggered()), this, SLOT(viewFitToWindow()));
     connect(_ui->normalSizeAction, SIGNAL(triggered()), this, SLOT(viewNormalSize()));

    _defaultTimeEvolution = new DefaultTimeEvolution;
    _quantum3bodyTimeEvolution = new Quantum3BodyTimeEvolution;

    connect(_timer, SIGNAL(timeout()), SLOT(evolve()));
    connect(_ui->run, SIGNAL(toggled(bool)), SLOT(runSimulation(bool)));
    connect(_ui->reset, SIGNAL(pressed()), SLOT(resetSimulation()));
    connect(_ui->browsePictureFolder, SIGNAL(pressed()), SLOT(browsePictureFolder()));
    connect(_ui->initialTimeEvolution, SIGNAL(currentIndexChanged(int)), SLOT(setTimeEvolution(int)));

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
    delete _defaultTimeEvolution;
    delete _quantum3bodyTimeEvolution;
    delete _image;
}

void Quantum3BodySimulationWindow::resetSimulation()
{
    statusBar()->showMessage("resetting the simulation, this may take a while...");

    size_t gridSizeX(0), gridSizeY(0);

    switch (_ui->resolution->currentIndex())
    {
        case RESOLUTION_64x64:
            gridSizeX = gridSizeY = 64;
            break;
        case RESOLUTION_128x128:
            gridSizeX = gridSizeY = 128;
            break;
        case RESOLUTION_256x256:
            gridSizeX = gridSizeY = 256;
            break;
        case RESOLUTION_512x512:
            gridSizeX = gridSizeY = 512;
            break;
        case RESOLUTION_1024x1024:
            gridSizeX = gridSizeY = 1024;
            break;
        case RESOLUTION_2048x2048:
            gridSizeX = gridSizeY = 2048;
            break;
        case RESOLUTION_4096x4096:
            gridSizeX = gridSizeY = 4096;
            break;
        case RESOLUTION_8192x8192:
            gridSizeX = gridSizeY = 8192;
            break;
        default:
            assert(!"the selected grid size is not implemented");
    }

    // only reinitialize the simulation if the resolution changed
    if (_simulation == nullptr || gridSizeX != _simulation->sizeX() || gridSizeY != _simulation->sizeY())
    {
        delete _simulation;
        delete _image;

        _simulation = new TwoDimSPO(gridSizeX, gridSizeY);
        _image = new QImage(static_cast<int>(gridSizeY), static_cast<int>(gridSizeX), QImage::Format_RGB32);
    }

    switch (_ui->initialPotential->currentIndex())
    {
        case HARMONIC_POTENTIAL:
            _defaultTimeEvolution->useHarmonicPotential = true;
            break;
        case ZERO_POTENTIAL:
            _defaultTimeEvolution->useHarmonicPotential = false;
            break;
        default:
            assert(!"Selected potential not defined.");
    }

    _quantum3bodyTimeEvolution->e = _ui->initialEpsilon->value();

    switch (_ui->initialTimeEvolution->currentIndex())
    {
        case DEFAULT_EVOLUTION:
            _simulation->setTimeEvolution(_defaultTimeEvolution);
            break;
        case QUANTUM3BODY_EVOLUTION:
            _simulation->setTimeEvolution(_quantum3bodyTimeEvolution);
            break;
        default:
            assert(!"Selected time evolution not defined.");
    }

    auto phi0 = [&](const double& x, const double& y)->complex {
        double kx(_ui->initialPropagationX->value());
        double ky(_ui->initialPropagationY->value());
        double dx(_ui->initialPositionX->value());
        double dy(_ui->initialPositionY->value());
        return exp(-0.5*((x-dx)*(x-dx)+(y-dy)*(y-dy)) + complex(0,1)*ky*y + complex(0,1)*kx*x);
    };

    _simulation->initialize(phi0);

    _currentIteration = 0;
    _ui->currentIteration->setValue(static_cast<int>(_currentIteration));

    _lastResetTimestamp = QDateTime::currentDateTime().toTime_t();

    plot();
    statusBar()->showMessage("ready...");
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
    _ui->currentIteration->setValue(static_cast<int>(_currentIteration));

    _simulation->evolveStep(_ui->dtParameter->value());

    if (_currentIteration % _ui->updatePlotSteps->value() == 0)
        plot();

    if (_currentIteration == static_cast<size_t>(_ui->numberOfIterations->value()))
        _ui->run->setChecked(false);
}

void Quantum3BodySimulationWindow::plot()
{
    const size_t gridSizeX(_simulation->sizeX()), gridSizeY(_simulation->sizeY());

    double totalProbability(0.0);
    const double binSizeX(_simulation->binSizeX()), binSizeY(_simulation->binSizeY());
    const complex* f(_simulation->phi());
    for (size_t i(0); i < gridSizeX; ++i)
    {
        for (size_t j(0); j < gridSizeY; ++j)
        {
            totalProbability += binSizeX*binSizeY*fabs(f[j + gridSizeY*i])*fabs(f[j + gridSizeY*i]);
        }
    }
    _ui->totalProbability->setValue(totalProbability);

    for(size_t i(0); i < gridSizeX; ++i)
    {
        for(size_t j(0); j < gridSizeY; ++j)
        {
            _image->setPixel(static_cast<int>(i), static_cast<int>(gridSizeY-1-j), rainbowColorMap(fabs(f[j + gridSizeY*i])));
        }
    }
    _imageLabel->setPixmap(QPixmap::fromImage(*_image));
    if (!_ui->fitToWindowAction->isChecked())
        _imageLabel->adjustSize();

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

            qDebug() << "dumping picture to " << filepath;
            _image->save(filepath, _ui->pictureFormat->currentText().toAscii().constData());
        }
    }
}

void Quantum3BodySimulationWindow::setTimeEvolution(int selectedEvolutionAlgorithm)
{
    switch (selectedEvolutionAlgorithm)
    {
        case DEFAULT_EVOLUTION:
            _ui->initialPotential->setEnabled(true);
            _ui->initialEpsilon->setEnabled(false);
            break;
        case QUANTUM3BODY_EVOLUTION:
            _ui->initialPotential->setEnabled(false);
            _ui->initialEpsilon->setEnabled(true);
            break;
        default:
            assert(!"Selected time evolution not defined.");
    }
}

void Quantum3BodySimulationWindow::viewNormalSize()
{
    _imageLabel->adjustSize();
}

void Quantum3BodySimulationWindow::viewFitToWindow()
{
    bool fitToWindow = _ui->fitToWindowAction->isChecked();
    _ui->plotScrollArea->setWidgetResizable(fitToWindow);
    if (!fitToWindow) {
        viewNormalSize();
    }
    updateViewActions();
}

void Quantum3BodySimulationWindow::updateViewActions()
{
    _ui->normalSizeAction->setEnabled(!_ui->fitToWindowAction->isChecked());
}

