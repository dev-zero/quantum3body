/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *
 *
 *
 */

#include "quantum3body_simulation_window.hh"

#include <QtGui/QApplication>

int main(int argc, char* argv[])
{
    QApplication app(argc, argv);

    Quantum3BodySimulationWindow window;
    window.show();

    return app.exec();
}

