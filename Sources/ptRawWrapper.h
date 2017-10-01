#ifndef PTRAWWRAPPER_H
#define PTRAWWRAPPER_H

#include "libraw/libraw.h"

#include <QString>

class ptRawWrapper : private LibRaw
{
public:
    ptRawWrapper();

    int open_file(const QString &fileName);
    using LibRaw::unpack;
    using LibRaw::raw2image;

    using LibRaw::unpack_thumb;

    using LibRaw::adjust_bl;
    using LibRaw::subtract_black_internal;
    using LibRaw::scale_colors;

    libraw_data_t &imgdata;
};

#endif // PTRAWWRAPPER_H
