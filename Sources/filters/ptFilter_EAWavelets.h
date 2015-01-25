/*******************************************************************************
**
** Photivo
**
** Copyright (C) 2015 Bernd Schoeler <brjohn@brother-john.net>
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

#ifndef PTFILTER_EAWavelets_H
#define PTFILTER_EAWavelets_H

#include "ptFilterBase.h"

class ptFilter_EAWavelets: public ptFilterBase {
  Q_OBJECT

public:
  static ptFilterBase *createEAWavelets();

protected:
  QWidget*  doCreateGui() override;
  void      doDefineControls() override;
  bool      doCheckHasActiveCfg() override;
  void      doRunFilter(ptImage *AImage) const override;

private:
  ptFilter_EAWavelets();


private slots:
  void onMasterValueChanged(QString, QVariant ANewValue);
  void onLevelValueChanged(QString AId, QVariant ANewValue);
  
};

#endif // PTFILTER_EAWavelets_H
