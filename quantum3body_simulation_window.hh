/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *                    Christian Reinhardt
 *
 *
 *
 */

#ifndef QUANTUM3BODY_SIMULATION_WINDOW_HH
#define QUANTUM3BODY_SIMULATION_WINDOW_HH

#include <QtGui/QMainWindow>

class QTimer;
class TwoDimSPO;
class DefaultTimeEvolution;
class Quantum3BodyTimeEvolution;
class QuantumPixelPlot;
class QLabel;
class QImage;

namespace Ui
{
    class Quantum3BodyWindow;
}
class Quantum3BodySimulationWindow :
    public QMainWindow
{
    Q_OBJECT
public:
    Quantum3BodySimulationWindow(QWidget* p = nullptr);
    ~Quantum3BodySimulationWindow();

public slots:
    void evolve();
    void plot();
    void runSimulation(bool);
    void resetSimulation();
    void browsePictureFolder();
    void setTimeEvolution(int);

    void viewNormalSize();
    void viewFitToWindow();
    void updateViewActions();
private:
    TwoDimSPO* _simulation;
    DefaultTimeEvolution* _defaultTimeEvolution;
    Quantum3BodyTimeEvolution* _quantum3bodyTimeEvolution;

    QTimer* _timer;
    Ui::Quantum3BodyWindow* _ui;
//    QuantumPixelPlot* _spatialPlot;
    size_t _currentIteration;
    unsigned int _lastResetTimestamp;
    QLabel* _imageLabel;
    QImage* _image;
};

#endif // QUANTUM3BODY_SIMULATION_WINDOW_HH
