/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *
 *
 *
 */

#ifndef QUANTUM3BODY_SIMULATION_WINDOW_HH
#define QUANTUM3BODY_SIMULATION_WINDOW_HH

#include <QtGui/QMainWindow>

class QTimer;
class Quantum3BodySimulation;
class QuantumPixelPlot;

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

private:
    Quantum3BodySimulation* _simulation;
    QTimer* _timer;
    Ui::Quantum3BodyWindow* _ui;
    QuantumPixelPlot* _spatialPlot;
    size_t _currentIteration;
};

#endif // QUANTUM3BODY_SIMULATION_WINDOW_HH
