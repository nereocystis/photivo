#include "ptRawWrapper.h"

ptRawWrapper::ptRawWrapper() : imgdata(LibRaw::imgdata)
{

}

int ptRawWrapper::open_file(const QString &fileName)
{
    return LibRaw::open_file(fileName.toLocal8Bit().constData());
}
