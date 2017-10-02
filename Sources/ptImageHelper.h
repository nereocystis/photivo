/*******************************************************************************
**
** Photivo
**
** Copyright (C) 2012 Michael Munzert <mail@mm-log.com>
**
** This file is part of Photivo.
**
** Photivo is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 3
** as published by the Free Software Foundation.
**
** Photivo is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License
** along with Photivo.  If not, see <http://www.gnu.org/licenses/>.
**
*******************************************************************************/

#ifndef PTIMAGEHELPER_H
#define PTIMAGEHELPER_H

//==============================================================================

#include <QString>

#include <vector>
#include <cstdint>

//==============================================================================

/*!
  \class ptImageHelper

  \brief Collection of functions to work with images on disk.

  This class provides (only static) functions to work with images on disk.
  */

//==============================================================================

namespace Exiv2 {
class ExifData;
class IptcData;
class XmpData;
} // namespace Exiv2

class QImage;

class ptImageHelper
{
public:
  /*! Write a given exif buffer to a file.*/
  static bool WriteExif(const QString              &AFileName,
                        const std::vector<uint8_t> &AExifBuffer,
                        Exiv2::IptcData            &AIptcData,
                        Exiv2::XmpData             &AXmpData);

  /*! Read exif data from file.*/
  static bool ReadExif(const QString        &AFileName,
                       Exiv2::ExifData      &AExifData,
                       std::vector<uint8_t> &AExifBuffer);

  /*! Transfer exif data from one image to another.*/
  static bool TransferExif(const QString ASourceFile,
                           const QString ATargetFile);

  /*! Write a QImage to disk.*/
  static bool DumpImage(QImage       *AImage,
                        const QString AFileName);
private:
  ptImageHelper();
  ~ptImageHelper();
};

//==============================================================================

#endif // PTIMAGEHELPER_H
