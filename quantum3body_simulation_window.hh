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

// start forward declarations to reduce header dependencies
class QTimer;
class TwoDimSPO;
struct DefaultTimeEvolution;
struct Quantum3BodyTimeEvolution;
class QLabel;
class QImage;

namespace Ui
{
    class Quantum3BodyWindow;
}
// end forward declarations

/**
 * The main class for all the GUI stuff.
 * This is an implementation of a QMainWindow which sets itself up using the
 * generated code from the Qt Designer.
 */

class Quantum3BodySimulationWindow :
    public QMainWindow
{
    Q_OBJECT
public:

    /**
     * Construct a new MainWindow object.
     */
    Quantum3BodySimulationWindow(QWidget* p = NULL);
    ~Quantum3BodySimulationWindow();

public slots:
    /**
     * Run one step in the simulation using the dt specified in the interface
     */
    void evolve();

    /**
     * Update/Redraw the plots using the current data found in the simulation.
     */
    void plot();

    /**
     * Start/stop the continuous simulation.
     */
    void runSimulation(bool);

    /**
     * Reset the simulation using the values given in the interface
     */
    void resetSimulation();

    /**
     * Open a dialog to browse for a folder to save the rendered pictures.
     */
    void browsePictureFolder();

    /**
     * Select the time evolution algorithm to be used.
     */
    void setTimeEvolution(int);

    /**
     * Display the rendered plot 1:1.
     */
    void viewNormalSize();

    /**
     * Resize to the rendered plot to match the window.
     */
    void viewFitToWindow();

    /**
     * Update the view actions after selecting/deselecting FitToWindow
     */
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
