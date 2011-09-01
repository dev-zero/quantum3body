/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *
 *
 *
 */

#include "plot_graphics_item.hh"
#include "color_maps.hh"

#include <QtGui/QPainter>
#include <QtCore/QDebug>

PlotGraphicsItem::PlotGraphicsItem(QGraphicsItem* p) :
    QGraphicsItem(p),
    _sizeX(0),
    _sizeY(0),
    _data(nullptr)
{
}

QRectF PlotGraphicsItem::boundingRect() const
{
    return QRectF(-static_cast<qreal>(_sizeX/2), -static_cast<qreal>(_sizeY/2), _sizeX, _sizeY);
}

void PlotGraphicsItem::paint(QPainter* painter, const QStyleOptionGraphicsItem*, QWidget*)
{
    for (size_t i(0); i < _sizeX; ++i)
    {
        for (size_t j(0); j < _sizeY; ++j)
        {
            painter->setPen(QPen(rainbowColorMap(fabs(_data[j + _sizeY*i]))));
            painter->drawPoint(i, j);
        }
    }
}

