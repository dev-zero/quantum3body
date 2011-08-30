/* vim: set sw=4 sts=4 ft=cpp et foldmethod=syntax : */

/*
 * Copyright (c) 2011 Tiziano Müller <tm@dev-zero.ch>
 *                    Christian Reinhardt
 *
 *
 *
 */

#ifndef QUANTUM_PIXEL_PLOT_HH
#define QUANTUM_PIXEL_PLOT_HH

#include <QtCore/QAbstractTableModel>
#include <QtGui/QAbstractItemDelegate>
#include <vector>
#include <complex>

extern QRgb rainbowColorMap(double value);

class QuantumPixelPlot :
    public QAbstractTableModel
{
    Q_OBJECT

public:
    QuantumPixelPlot(QObject* p = nullptr);
    void setSpatialData(const std::complex<double>* data, size_t xSize, size_t ySize);

    int columnCount(const QModelIndex &parent = QModelIndex()) const;
    int rowCount(const QModelIndex &parent = QModelIndex()) const;

    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const;

private:
    const std::complex<double>* _modelData;
    size_t _xSize, _ySize;
};

class PixelDelegate : public QAbstractItemDelegate
{
    Q_OBJECT

public:
    PixelDelegate(QObject *parent = 0);

    void paint(QPainter *painter, const QStyleOptionViewItem &option, const QModelIndex &index) const;

    QSize sizeHint(const QStyleOptionViewItem &option, const QModelIndex &index ) const;

public slots:
    void setPixelSize(int size);

private:
    int pixelSize;
};

#endif // QUANTUM_PIXEL_PLOT_HH
