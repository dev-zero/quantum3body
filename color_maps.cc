/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *                    Christian Reinhardt
 *
 *
 */

#include "color_maps.hh"

QRgb rainbowColorMap(double value)
{
    if(value < 0.0) value = 0.0;
    if(value > 1.0) value = 1.0;

    unsigned int r(0), g(0), b(0);
    if(value < 0.2)
    {
        r = 255*5*(0.2 - value);
        b = 255;
    } else if(value < 0.4)
    {
        b = 255;
        g = 255*5*(value - 0.2);
    } else if(value < 0.6)
    {
        g = 255;
        b = 255*5*(0.6 - value);
    } else if(value < 0.8)
    {
        g = 255;
        r = 255*5*(value - 0.6);
    } else
    {
        r = 255;
        g = 255*5*(1.0 - value);
    }

    return qRgb(r, g, b);
}

