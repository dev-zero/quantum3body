/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *                    Christian Reinhard
 *
 *
 */

#ifndef COLOR_MAPS_HH
#define COLOR_MAPS_HH

#include <QtGui/QRgb>

/**
 * @file
 * This file contains color maps used in several of our
 * programs
 */

/**
 * Transforms the given double value to a color value using
 * all rainbow colors. "red" means a high value, "purple" a low value.
 *
 * @param value This is the value being mapped to the color value.
 *              It should be in the range of 0..1, otherwise it gets
 *              capped.
 */
QRgb rainbowColorMap(double value);

#endif // COLOR_MAPS_HH
