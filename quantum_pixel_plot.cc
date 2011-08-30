/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *
 *
 *
 */

#include "quantum_pixel_plot.hh"

#include <QtGui/QPainter>
#include <QtCore/QDebug>

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

QuantumPixelPlot::QuantumPixelPlot(QObject *parent)
    : QAbstractTableModel(parent)
{
}

void QuantumPixelPlot::setSpatialData(const std::vector<std::complex<double>>& data, size_t xSize, size_t ySize)
{
    _modelData = data;
    _xSize = xSize;
    _ySize = ySize; 
    reset();
}

int QuantumPixelPlot::rowCount(const QModelIndex & /* parent */) const
{
    return _ySize;
}

int QuantumPixelPlot::columnCount(const QModelIndex & /* parent */) const
{
    return _xSize;
}

QVariant QuantumPixelPlot::data(const QModelIndex &index, int role) const
{
    if (!index.isValid())
        return QVariant();

    size_t idx(index.row() + _ySize*index.column());

    if (role == Qt::DisplayRole)
        return rainbowColorMap(abs(_modelData[idx]));

    if (role == Qt::ToolTipRole)
        return QString("Value: %1 + i%2\nAbsolute: %3")
            .arg(_modelData[idx].real())
            .arg(_modelData[idx].imag())
            .arg(abs(_modelData[idx]));

    return QVariant();
}

QVariant QuantumPixelPlot::headerData(int /* section */,
        Qt::Orientation /* orientation */,
        int role) const
{
    if (role == Qt::SizeHintRole)
        return QSize(1, 1);
    return QVariant();
}

PixelDelegate::PixelDelegate(QObject *parent)
    : QAbstractItemDelegate(parent)
{
    pixelSize = 4;
}

void PixelDelegate::paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const
{
    if (option.state & QStyle::State_Selected)
        painter->fillRect(option.rect, option.palette.highlight());

#ifdef SIMPLE_PIXEL
    painter->fillRect(option.rect, index.model()->data(index, Qt::DisplayRole).value<QRgb>());    
#else
    int size = qMin(option.rect.width(), option.rect.height());
    QRgb color = index.model()->data(index, Qt::DisplayRole).value<QRgb>();
    double radius = size/2.0;

    painter->save();
    painter->setRenderHint(QPainter::Antialiasing, true);
    painter->setPen(Qt::NoPen);

    if (option.state & QStyle::State_Selected)
        painter->setBrush(option.palette.highlightedText());
    else
        painter->setBrush(QBrush(color));
    painter->drawEllipse(QRectF(option.rect.x() + option.rect.width()/2 - radius,
                option.rect.y() + option.rect.height()/2 - radius,
                2*radius, 2*radius));
    painter->restore();
#endif
}

QSize PixelDelegate::sizeHint(const QStyleOptionViewItem & /* option */,
        const QModelIndex & /* index */) const
{
    return QSize(pixelSize, pixelSize);
}

void PixelDelegate::setPixelSize(int size)
{
    pixelSize = size;
}
