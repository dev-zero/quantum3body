/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *
 *
 *
 */

#ifndef PLOT_GRAPHICS_ITEM_HH
#define PLOT_GRAPHICS_ITEM_HH

#include <QtGui/QGraphicsItem>

#include <complex>

class PlotGraphicsItem :
    public QGraphicsItem
{
public:
    PlotGraphicsItem(QGraphicsItem* p = nullptr);
    
    QRectF boundingRect() const;
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget = nullptr);
    void setData(const std::complex<double>* d, size_t sizeX, size_t sizeY)
    {
        _data = d;
        _sizeX = sizeX;
        _sizeY = sizeY;
    }
private:
    size_t _sizeX, _sizeY;
    const std::complex<double>* _data;
};

#endif // PLOT_GRAPHICS_ITEM_HH
