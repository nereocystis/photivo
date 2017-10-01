#ifndef PTCOLORPROFILES_H
#define PTCOLORPROFILES_H

#include "ptConstants.h"

#include <tools/finally.h>

#include <lcms2.h>
#include <lcms2_plugin.h>
#include <libraw/libraw_types.h>

#include <cassert>

namespace pt {

cmsHPROFILE MatrixToProfile(cmsMAT3 MatrixXYZToCam) {

  cmsMAT3 MatrixCamToXYZ;
  _cmsMAT3inverse(&MatrixXYZToCam, &MatrixCamToXYZ);

  cmsCIEXYZ RedXYZ = {MatrixCamToXYZ.v[0].n[0], MatrixCamToXYZ.v[1].n[0],
                      MatrixCamToXYZ.v[2].n[0]};
  cmsCIEXYZ GreenXYZ = {MatrixCamToXYZ.v[0].n[1], MatrixCamToXYZ.v[1].n[1],
                        MatrixCamToXYZ.v[2].n[1]};
  cmsCIEXYZ BlueXYZ = {MatrixCamToXYZ.v[0].n[2], MatrixCamToXYZ.v[1].n[2],
                       MatrixCamToXYZ.v[2].n[2]};

  cmsCIExyY RedxyY;
  cmsCIExyY GreenxyY;
  cmsCIExyY BluexyY;

  cmsXYZ2xyY(&RedxyY, &RedXYZ);
  cmsXYZ2xyY(&GreenxyY, &GreenXYZ);
  cmsXYZ2xyY(&BluexyY, &BlueXYZ);

  cmsCIExyYTRIPLE Primaries = {RedxyY, GreenxyY, BluexyY};

  // double gamma  = 0.45;
  // double linear = 0.10;
  // double g = gamma*(1.0-linear)/(1.0-linear*gamma);
  // double a = 1.0/(1.0+linear*(g-1.0));
  // double b = linear*(g-1.0)*a;
  // double c = pow(a*linear+b,g)/linear;
  // double Params[] = {g,a,b,c,linear};
  // Gamma => None (=1)
  cmsToneCurve *Gamma = cmsBuildGamma(NULL, 1.0);
  cmsToneCurve *Gamma3[3];
  Gamma3[0] = Gamma3[1] = Gamma3[2] = Gamma;

  cmsHPROFILE Profile = cmsCreateRGBProfile(cmsD50_xyY(), &Primaries, Gamma3);
  cmsFreeToneCurve(Gamma);

  return Profile;
}

cmsHPROFILE createFlatProfile() {
  cmsMAT3 MatrixXYZToCam;
  for (short j = 0; j < 3; j++) {
    for (short k = 0; k < 3; k++) {
      MatrixXYZToCam.v[j].n[k] = MatrixXYZToRGB[ptSpace_sRGB_D65][j][k];
    }
  }

  cmsHPROFILE Profile = MatrixToProfile(MatrixXYZToCam);
  cmsSetDeviceClass(Profile, cmsSigInputClass);
  cmsSetColorSpace(Profile, cmsSigRgbData);
  cmsSetPCS(Profile, cmsSigXYZData);

  return Profile;
}

cmsHPROFILE createGammaProfile(ptCameraColorGamma ProfileGamma) {

  double Params[6];
  switch (ProfileGamma) {
  case ptCameraColorGamma::sRGB:
    // Parameters for standard (inverse) sRGB in lcms
    Params[0] = 1.0 / 2.4;
    Params[1] = 1.1371189;
    Params[2] = 0.0;
    Params[3] = 12.92;
    Params[4] = 0.0031308;
    Params[5] = -0.055;
    break;
  case ptCameraColorGamma::BT709:
    // Parameters for standard (inverse) BT709 in lcms
    Params[0] = 0.45;
    Params[1] = 1.233405791;
    Params[2] = 0.0;
    Params[3] = 4.5;
    Params[4] = 0.018;
    Params[5] = -0.099;
    break;
  case ptCameraColorGamma::Pure22:
    // Parameters for standard (inverse) 2.2 gamma in lcms
    Params[0] = 1.0 / 2.2;
    Params[1] = 1.0;
    Params[2] = 0.0;
    Params[3] = 0.0;
    Params[4] = 0.0;
    Params[5] = 0.0;
    break;
  default:
    assert(0);
  }

  cmsToneCurve *Gamma = cmsBuildParametricToneCurve(0, 4, Params);
  finally f([&] { cmsFreeToneCurve(Gamma); });

  cmsToneCurve *Gamma4[4];
  Gamma4[0] = Gamma4[1] = Gamma4[2] = Gamma4[3] = Gamma;

  return cmsCreateLinearizationDeviceLink(cmsSigRgbData, Gamma4);
}

cmsHPROFILE createAdobeProfile(const float cam_xyz[4][3]) {
  //  bool HaveFourColorProblem = 0;
  //  for (short j = 9; j < 12; j++) {
  //    if (cam_xyz[j] != 0) {
  //      HaveFourColorProblem = 1;
  //      break;
  //    }
  //  }

  cmsMAT3 MatrixXYZToCam;
  for (short j = 0; j < 3; j++) {
    for (short k = 0; k < 3; k++) {
      MatrixXYZToCam.v[j].n[k] = cam_xyz[j][k];
    }
  }

  cmsHPROFILE Profile = MatrixToProfile(MatrixXYZToCam);
  cmsSetDeviceClass(Profile, cmsSigInputClass);
  cmsSetColorSpace(Profile, cmsSigRgbData);
  cmsSetPCS(Profile, cmsSigXYZData);

  return Profile;
}

} // namespace pt

#endif // PTCOLORPROFILES_H
