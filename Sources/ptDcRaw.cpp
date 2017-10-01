/*******************************************************************************
**
** Photivo
**
** Copyright (C) 2008,2009 Jos De Laender <jos.de_laender@telenet.be>
** Copyright (C) 2009-2012 Michael Munzert <mail@mm-log.com>
** Copyright (C) 2011 Bernd Schoeler <brjohn@brother-john.net>
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
/*
** This is basically the translation into a more or less C++ object of
** dcraw.c -- Dave Coffin's raw photo decoder
** Copyright 1997-2008 by Dave Coffin, dcoffin a cybercom o net
**
*******************************************************************************/

#include "ptDcRaw.h"
#include "ptRawWrapper.h"

#include "ptDefines.h"
#include "ptError.h"
#include "ptConstants.h"
#include "ptCalloc.h"
#include "tools/finally.h"

#include <cassert>

#define NO_JASPER
#ifndef NO_JASPER
#include <jasper/jasper.h>
#endif

// Macro fix for explicit fread returnvalue check.
#define ptfread(ptr,size,n,stream)              \
{                                               \
size_t RV = fread(ptr,size,n,stream);           \
if (RV != (size_t) n) assert(!ferror(stream));  \
}
#define ptfwrite(ptr,size,n,stream)    \
{                                      \
size_t RV = fwrite(ptr,size,n,stream); \
assert(RV == (size_t) n);              \
}
#define ptfscanf(file,format,arg)      \
{                                      \
int RV = fscanf(file,format,arg);      \
assert (RV == 1);                      \
}
#define ptfgets(str,num,file)   \
{                               \
char* RV = fgets(str,num,file); \
assert (RV);                    \
}

inline void VAppend(TImage8RawData &AVector, char* AArray, const int ALength) {
  char* hEnd = AArray + ALength * sizeof(char);
  AVector.insert(AVector.end(), AArray, hEnd);
}

template <typename T> void array_copy(const T &source, T &target) {
  std::copy(std::begin(source), std::end(source), std::begin(target));
}

// The class.
#define CLASS ptDcRaw::
CLASS ptDcRaw()
    : m_rawProcessor(new ptRawWrapper),
      imgdata(m_rawProcessor->imgdata)
{
  //printf("(%s,%d) '%s'\n",__FILE__,__LINE__,__PRETTY_FUNCTION__);

  // This were the original global variables initialized.
  // Now moved into constructor.
  // All m_UserSetting* are obviously ,uh, usersettings
  // that were done via the command line parameters.
  m_UserSetting_ShotSelect=0;
  m_UserSetting_Multiplier[0]=0;
  m_UserSetting_Multiplier[1]=0;
  m_UserSetting_Multiplier[2]=0;
  m_UserSetting_Multiplier[3]=0;
  m_UserSetting_HalfSize=0;
  m_UserSetting_HotpixelReduction=0;
  m_UserSetting_BayerDenoise=0;
  m_UserSetting_CfaLineDn=0;
  m_UserSetting_GreenEquil=0;
  m_UserSetting_CaCorrect=0;
  m_UserSetting_CaRed=0;
  m_UserSetting_CaBlue=0;
  m_UserSetting_AutoWb=0;
  m_UserSetting_CameraWb=0;
  m_UserSetting_CameraMatrix=-1;
  m_UserSetting_UseGreyBox = false;
  m_UserSetting_photivo_ClipMode = ptClipMode_Clip;
  m_UserSetting_photivo_ClipParameter = 0;
  m_UserSetting_Quality = 3;
  m_UserSetting_BlackPoint = -1;
  m_UserSetting_Saturation = -1;
  m_UserSetting_DetailView          = 0;
  m_UserSetting_DetailViewCropX     = 0;
  m_UserSetting_DetailViewCropY     = 0;
  m_UserSetting_DetailViewCropW     = 0;
  m_UserSetting_DetailViewCropH     = 0;
  m_UserSetting_BadPixelsFileName   = NULL;
  m_UserSetting_DarkFrameFileName   = NULL;
  m_UserSetting_AdjustMaximum       = 0;
  m_UserSetting_DenoiseThreshold    = 0;
  m_UserSetting_InterpolationPasses = 0;
  m_UserSetting_MedianPasses        = 0;
  m_UserSetting_ESMedianPasses      = 0;


  // Safety settings to have NULL on uninitialized images.
  m_Image = NULL;
  m_Image_AfterPhase1 = NULL;
  m_Image_AfterPhase2 = NULL;
  m_Image_AfterPhase3 = NULL;
  m_Image_AfterPhase4 = NULL;

  m_Thumb.clear();
  ResetNonUserSettings();
}

////////////////////////////////////////////////////////////////////////////////
//
// Destructor
// Deallocate everything dynamic.
//
////////////////////////////////////////////////////////////////////////////////

CLASS ~ptDcRaw() {

  delete m_rawProcessor;

//printf("(%s,%d) '%s'\n",__FILE__,__LINE__,__PRETTY_FUNCTION__);

  FREE(m_UserSetting_BadPixelsFileName);
  FREE(m_UserSetting_DarkFrameFileName);
  FREE(m_Image);
  FREE(m_Image_AfterPhase1);
  FREE(m_Image_AfterPhase2);
  FREE(m_Image_AfterPhase3);
}

////////////////////////////////////////////////////////////////////////////////
//
// ResetNonUserSettings
// Reset all variables except user settings.
// This is for second entry support.
//
////////////////////////////////////////////////////////////////////////////////

void CLASS ResetNonUserSettings() {
  // Safety settings to have NULL on uninitialized images.
  // And freeing the underlying memory (which was long time a leak !)
  // FREE(NULL) is safe, so the beginsituation is fine too.
  // FREE implies setting of the pointer to NULL
  FREE(m_Image);
  FREE(m_Image_AfterPhase1);
  FREE(m_Image_AfterPhase2);
  FREE(m_Image_AfterPhase3);
  FREE(m_Image_AfterPhase4);

  // This was originally in the identify code, but which is called
  // anyway in the beginning. So this is simply global initialization like
  // anything else.

  m_Tiff_Flip = m_Flip = -1; /* 0 is valid, so -1 is unknown */
  m_Filters = (unsigned)(-1); // hack not to change dcraw.
  m_RawHeight = m_RawWidth = m_Fuji_Width = m_IsFuji = fuji_layout = cr2_slice[0] = 0;
  m_WhiteLevel = m_Height = m_Width = m_TopMargin = m_LeftMargin = 0;
  m_ColorDescriptor[0] = m_Description[0] = m_Artist[0] = 0;
  m_IsoSpeed = m_Shutter = m_Aperture = m_FocalLength = unique_id = 0;
  m_Tiff_NrIFDs = 0;
  memset (white, 0, sizeof white);
  m_Data_Offset = m_Tiff_bps = m_Tiff_Compress = 0;
  m_Kodak_cbpp = zero_after_ff = m_DNG_Version = m_Load_Flags = 0;
  m_TimeStamp = m_ShotOrder = m_Tiff_Samples = m_BlackLevel = m_IsFoveon = 0;
  memset (m_CBlackLevel, 0, sizeof m_CBlackLevel);
  memset (m_CBlackLevel_AfterPhase1, 0, sizeof m_CBlackLevel_AfterPhase1);
  m_MixGreen = m_ProfileLength = data_error = m_ZeroIsBad = 0;
  m_PixelAspect = m_IsRaw = m_RawColor = 1; m_RawColorPhotivo = 0;
  m_TileWidth = m_TileLength = 0;
  m_Raw_Image = 0;
  memset (m_Mask, 0, sizeof m_Mask);
  for (int i=0; i < 4; i++) {
    short c;
    ASSIGN(m_CameraMultipliers[i], i == 1);
    ASSIGN(m_PreMultipliers[i], i < 3);
    ASSIGN(m_D65Multipliers[i], i < 3);
    for (c=0; c<3; c++) m_cmatrix[c][i] = 0;
    for (c=0; c<3; c++) m_MatrixCamRGBToSRGB[c][i] = c == i;
  }
  m_Colors = 3;
  for (int i=0; i < 0x10000; i++) m_Curve[i] = i;

  m_Gamma[0] = 0.45;
  m_Gamma[1] = 4.50;
  m_Gamma[2] = 0;
  m_Gamma[3] = 0;
  m_Gamma[4] = 0;
  m_Gamma[5] = 0;

  m_getbithuff_bitbuf=0;
  m_getbithuff_reset=0;
  m_getbithuff_vbits=0;
  for (int i = 0; i < 0x4000; i++) m_pana_bits_buf[i] = 0;
  m_pana_bits_vbits = 0;
  for (int i = 0; i < 4096; i++) jpeg_buffer[i] = 0;
  for (int i = 0; i < 128; i++) m_sony_decrypt_pad[i] = 0;
  m_sony_decrypt_p = 0;
  for (int i = 0; i < 1024; i++) m_foveon_decoder_huff[i] = 0;

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 3; j++) {
      MatrixXYZToCam[i][j] = 0.0;
    }
  }

  ToCamFunctionInited = 0;
  for (int i = 0; i < 0x20000; i++) ToLABFunctionTable[i] = 0.0;
  ToLABFunctionInited = 0;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 4; j++) {
      MatrixCamToXYZ[i][j] = 0.0;
    }
  }

}

/*
   In order to inline this calculation, I make the risky
   assumption that all filter patterns can be described
   by a repeating pattern of eight rows and two columns

   Do not use the FC or BAYER macros with the Leaf CatchLight,
   because its pattern is 16x16, not 2x8.

   Return values are either 0/1/2/3 = G/M/C/Y or 0/1/2/3 = R/G1/B/G2

  PowerShot 600 PowerShot A50 PowerShot Pro70 Pro90 & G1
  0xe1e4e1e4: 0x1b4e4b1e: 0x1e4b4e1b: 0xb4b4b4b4:

    0 1 2 3 4 5   0 1 2 3 4 5   0 1 2 3 4 5   0 1 2 3 4 5
  0 G M G M G M 0 C Y C Y C Y 0 Y C Y C Y C 0 G M G M G M
  1 C Y C Y C Y 1 M G M G M G 1 M G M G M G 1 Y C Y C Y C
  2 M G M G M G 2 Y C Y C Y C 2 C Y C Y C Y
  3 C Y C Y C Y 3 G M G M G M 3 G M G M G M
      4 C Y C Y C Y 4 Y C Y C Y C
  PowerShot A5  5 G M G M G M 5 G M G M G M
  0x1e4e1e4e: 6 Y C Y C Y C 6 C Y C Y C Y
      7 M G M G M G 7 M G M G M G
    0 1 2 3 4 5
  0 C Y C Y C Y
  1 G M G M G M
  2 C Y C Y C Y
  3 M G M G M G

   All RGB cameras use one of these Bayer grids:

  0x16161616: 0x61616161: 0x49494949: 0x94949494:

    0 1 2 3 4 5   0 1 2 3 4 5   0 1 2 3 4 5   0 1 2 3 4 5
  0 B G B G B G 0 G R G R G R 0 G B G B G B 0 R G R G R G
  1 G R G R G R 1 B G B G B G 1 R G R G R G 1 G B G B G B
  2 B G B G B G 2 G R G R G R 2 G B G B G B 2 R G R G R G
  3 G R G R G R 3 B G B G B G 3 R G R G R G 3 G B G B G B
 */

#define RAW(row,col) \
  m_Raw_Image[(row)*m_RawWidth+(col)]

#define FC(row,col) \
  (m_Filters >> ((((row) << 1 & 14) + ((col) & 1)) << 1) & 3)

#define BAYER(row,col) \
  m_Image[(row)*m_Width + (col)][FC(row,col)]
  // old: m_Image[((row) >> m_Shrink)*m_OutWidth + ((col) >> m_Shrink)][FC(row,col)]

#define BAYER2(row,col) \
  m_Image[(row)*m_Width + (col)][fcol(row,col)]

int CLASS fcol (int row, int col)
{
  static const char filter[16][16] =
  { { 2,1,1,3,2,3,2,0,3,2,3,0,1,2,1,0 },
    { 0,3,0,2,0,1,3,1,0,1,1,2,0,3,3,2 },
    { 2,3,3,2,3,1,1,3,3,1,2,1,2,0,0,3 },
    { 0,1,0,1,0,2,0,2,2,0,3,0,1,3,2,1 },
    { 3,1,1,2,0,1,0,2,1,3,1,3,0,1,3,0 },
    { 2,0,0,3,3,2,3,1,2,0,2,0,3,2,2,1 },
    { 2,3,3,1,2,1,2,1,2,1,1,2,3,0,0,1 },
    { 1,0,0,2,3,0,0,3,0,3,0,3,2,1,2,3 },
    { 2,3,3,1,1,2,1,0,3,2,3,0,2,3,1,3 },
    { 1,0,2,0,3,0,3,2,0,1,1,2,0,1,0,2 },
    { 0,1,1,3,3,2,2,1,1,3,3,0,2,1,3,2 },
    { 2,3,2,0,0,1,3,0,2,0,1,2,3,0,1,0 },
    { 1,3,1,2,3,2,3,2,0,2,0,1,1,0,3,0 },
    { 0,2,0,3,1,0,0,1,1,3,3,2,3,2,2,1 },
    { 2,1,3,2,3,1,2,1,0,3,0,2,0,2,0,2 },
    { 0,3,1,0,0,2,0,3,2,1,3,1,1,3,1,3 } };
  static const char filter2[6][6] =
  { { 1,1,0,1,1,2 },
    { 1,1,2,1,1,0 },
    { 2,0,1,0,2,1 },
    { 1,1,2,1,1,0 },
    { 1,1,0,1,1,2 },
    { 0,2,1,2,0,1 } };

  if (m_Filters == 1) return filter[(row+m_TopMargin)&15][(col+m_LeftMargin)&15];
  if (m_Filters == 2) return filter2[(row+6) % 6][(col+6) % 6];
  return FC(row,col);
}

#ifndef __GLIBC__
char* CLASS my_memmem (char *haystack, size_t haystacklen,
        char *neepte, size_t neeptelen)
{
  char *c;
  for (c = haystack; c <= haystack + haystacklen - neeptelen; c++)
    if (!memcmp (c, neepte, neeptelen))
      return c;
  return 0;
}
#define memmem my_memmem
#endif

float rgb_cam[3][4];
const double xyz_rgb[3][3] = {      /* XYZ from RGB */
  { 0.412453, 0.357580, 0.180423 },
  { 0.212671, 0.715160, 0.072169 },
  { 0.019334, 0.119193, 0.950227 } };
const float d65_white[3] = { 0.950456, 1, 1.088754 };

// Now everything importent is set up, so we can include external demosaicers
#include "dcb/dcb_demosaicing.c"
#include "dcb/dcb_demosaicing_old.c"
#include "vcd/ahd_interpolate_mod.c"
#include "vcd/ahd_partial_interpolate.c"
#include "vcd/refinement.c"
#include "vcd/vcd_interpolate.c"
#include "vcd/es_median_filter.c"
#include "vcd/median_filter_new.c"
#include "perfectraw/lmmse_interpolate.c"
#include "rawtherapee/amaze_interpolate.c"
#include "rawtherapee/cfa_line_dn.c"
#include "rawtherapee/ca_correct.c"
#include "rawtherapee/green_equil.c"

void CLASS merror (void *ptr, const char *where)
{
  if (ptr) return;
  fprintf (stderr,_("%s: Out of memory in %s\n"), m_UserSetting_InputFileName.toLocal8Bit().constData(), where);
  longjmp (m_Failure, 1);
}

void CLASS remove_zeroes()
{
  unsigned row, col, tot, n, r, c;

  for (row=0; row < m_Height; row++)
    for (col=0; col < m_Width; col++)
      if (BAYER(row,col) == 0) {
  tot = n = 0;
  for (r = row-2; r <= row+2; r++)
    for (c = col-2; c <= col+2; c++)
      if (r < m_Height && c < m_Width &&
    FC(r,c) == FC(row,col) && BAYER(r,c))
        tot += (n++,BAYER(r,c));
  if (n) BAYER(row,col) = tot/n;
      }
}

/*
   Seach from the current directory up to the root looking for
   a ".badpixels" file, and fix those pixels now.
 */
void CLASS bad_pixels (const char *cfname)
{
  FILE *fp=0;
  char *cp, line[128];
  int /* len, */ time, row, col, r, c, rad, tot, n;

  if (!m_Filters) return;
  if (cfname)
    fp = fopen (cfname, "r");
/* MASK AWAY IN dcRaw
  else {
    for (len=32 ; ; len *= 2) {
      fname = (char *) MALLOC (len);
      if (!fname) return;
      if (getcwd (fname, len-16)) break;
      FREE (fname);
      if (errno != ERANGE) return;
    }
#if defined(WIN32) || defined(DJGPP)
    if (fname[1] == ':')
      memmove (fname, fname+2, len-2);
    for (cp=fname; *cp; cp++)
      if (*cp == '\\') *cp = '/';
#endif
    cp = fname + strlen(fname);
    if (cp[-1] == '/') cp--;
    while (*fname == '/') {
      strcpy (cp, "/.badpixels");
      if ((fp = fopen (fname, "r"))) break;
      if (cp == fname) break;
      while (*--cp != '/');
    }
    FREE (fname);
  }
*/
  if (!fp) return;
  while (fgets (line, 128, fp)) {
    cp = strchr (line, '#');
    if (cp) *cp = 0;
    if (sscanf (line, "%d %d %d", &col, &row, &time) != 3) continue;
    if ((unsigned) col >= m_Width || (unsigned) row >= m_Height) continue;
    if (time > m_TimeStamp) continue;
    for (tot=n=0, rad=1; rad < 3 && n==0; rad++)
      for (r = row-rad; r <= row+rad; r++)
  for (c = col-rad; c <= col+rad; c++)
    if ((unsigned) r < m_Height && (unsigned) c < m_Width &&
    (r != row || c != col) && fcol(r,c) == fcol(row,col)) {
      tot += BAYER2(r,c);
      n++;
    }
    BAYER2(row,col) = tot/n;
    TRACEKEYVALS("Fixed dead pixel at column","%d",col);
    TRACEKEYVALS("Fixed dead pixel at row","%d",row);
  }
  FCLOSE (fp);
}

void CLASS subtract (const char *fname)
{
  FILE *fp;
  int dim[3]={0,0,0}, comment=0, number=0, error=0, nd=0, c, row, col;
  uint16_t *pixel;

  if (!(fp = fopen (fname, "rb"))) {
    perror (fname);  return;
  }
  if (fgetc(fp) != 'P' || fgetc(fp) != '5') error = 1;
  while (!error && nd < 3 && (c = fgetc(fp)) != EOF) {
    if (c == '#')  comment = 1;
    if (c == '\n') comment = 0;
    if (comment) continue;
    if (isdigit(c)) number = 1;
    if (number) {
      if (isdigit(c)) dim[nd] = dim[nd]*10 + c -'0';
      else if (isspace(c)) {
  number = 0;  nd++;
      } else error = 1;
    }
  }
  if (error || nd < 3) {
    fprintf (stderr,_("%s is not a valid PGM file!\n"), fname);
    FCLOSE (fp);  return;
  } else if (dim[0] != m_Width || dim[1] != m_Height || dim[2] != 65535) {
    fprintf (stderr,_("%s has the wrong dimensions!\n"), fname);
    FCLOSE (fp);  return;
  }
  pixel = (uint16_t*) CALLOC (m_Width, sizeof *pixel);
  merror (pixel, "subtract()");
  for (row=0; row < m_Height; row++) {
    ptfread (pixel, 2, m_Width, fp);
    for (col=0; col < m_Width; col++)
      BAYER(row,col) = MAX (BAYER(row,col) - ntohs(pixel[col]), 0);
  }
  FREE (pixel);
  FCLOSE (fp);
  memset (m_CBlackLevel, 0, sizeof m_CBlackLevel);
  m_BlackLevel = 0;
}

void CLASS gamma_curve (double pwr, double ts, int mode, int imax)
{
  int i;
  double g[6], bnd[2]={0,0}, r;

  g[0] = pwr;
  g[1] = ts;
  g[2] = g[3] = g[4] = 0;
  bnd[g[1] >= 1] = 1;
  if (g[1] && (g[1]-1)*(g[0]-1) <= 0) {
    for (i=0; i < 48; i++) {
      g[2] = (bnd[0] + bnd[1])/2;
      if (g[0]) bnd[(pow(g[2]/g[1],-g[0]) - 1)/g[0] - 1/g[2] > -1] = g[2];
      else  bnd[g[2]/exp(1-1/g[2]) < g[1]] = g[2];
    }
    g[3] = g[2] / g[1];
    if (g[0]) g[4] = g[2] * (1/g[0] - 1);
  }
  if (g[0]) g[5] = 1 / (g[1]*SQR(g[3])/2 - g[4]*(1 - g[3]) +
    (1 - pow(g[3],1+g[0]))*(1 + g[4])/(1 + g[0])) - 1;
  else      g[5] = 1 / (g[1]*SQR(g[3])/2 + 1
    - g[2] - g[3] - g[2]*g[3]*(log(g[3]) - 1)) - 1;
  if (!mode--) {
    memcpy (m_Gamma, g, sizeof m_Gamma);
    return;
  }
  for (i=0; i < 0x10000; i++) {
    m_Curve[i] = 0xffff;
    if ((r = (double) i / imax) < 1)
      m_Curve[i] = 0x10000 * ( mode
  ? (r < g[3] ? r*g[1] : (g[0] ? pow( r,g[0])*(1+g[4])-g[4]    : log(r)*g[2]+1))
  : (r < g[2] ? r/g[1] : (g[0] ? pow((r+g[4])/(1+g[4]),1/g[0]) : exp((r-1)/g[2]))));
  }
}

void CLASS pseudoinverse (double (*in)[3], double (*out)[3], int size)
{
  double work[3][6], num;
  int i, j, k;

  for (i=0; i < 3; i++) {
    for (j=0; j < 6; j++)
      work[i][j] = j == i+3;
    for (j=0; j < 3; j++)
      for (k=0; k < size; k++)
  work[i][j] += in[k][i] * in[k][j];
  }
  for (i=0; i < 3; i++) {
    num = work[i][i];
    for (j=0; j < 6; j++)
      work[i][j] /= num;
    for (k=0; k < 3; k++) {
      if (k==i) continue;
      num = work[k][i];
      for (j=0; j < 6; j++)
  work[k][j] -= work[i][j] * num;
    }
  }
  for (i=0; i < size; i++)
    for (j=0; j < 3; j++)
      for (out[i][j]=k=0; k < 3; k++)
  out[i][j] += work[j][k+3] * in[i][k];
}

void CLASS hat_transform (float *temp, float *base, int st, int size, int sc)
{
  int i;
  for (i=0; i < sc; i++)
    temp[i] = 2*base[st*i] + base[st*(sc-i)] + base[st*(i+sc)];
  for (; i+sc < size; i++)
    temp[i] = 2*base[st*i] + base[st*(i-sc)] + base[st*(i+sc)];
  for (; i < size; i++)
    temp[i] = 2*base[st*i] + base[st*(i-sc)] + base[st*(2*size-2-(i+sc))];
}

////////////////////////////////////////////////////////////////////////////////
//
// ptAdjustMaximum
// This is modelled after the adjust_maximum of libraw.
// It's purpose is false color suppression in the highlights.
// Typically "Threshold" should be between 0.75 and 1.
//
////////////////////////////////////////////////////////////////////////////////

void CLASS ptAdjustMaximum(double Threshold) {
  m_WhiteLevel = LIM((int32_t)(Threshold*(double)m_WhiteLevel),0, 4095);
}

////////////////////////////////////////////////////////////////////////////////
//
// ptAdjustBlacks
//
////////////////////////////////////////////////////////////////////////////////

void CLASS ptAdjustBlacks() {
  auto libraw_int_ptr = m_rawProcessor->imgdata.image;
  m_rawProcessor->imgdata.image = m_Image;
  finally f([&] { m_rawProcessor->imgdata.image = libraw_int_ptr; });
  m_rawProcessor->imgdata.params.user_black = m_UserSetting_BlackPoint;
  m_rawProcessor->imgdata.color.black = m_BlackLevel;
  m_rawProcessor->imgdata.color.maximum = m_WhiteLevel;
  m_rawProcessor->adjust_bl();
  m_rawProcessor->subtract_black_internal();
}

////////////////////////////////////////////////////////////////////////////////
//
// ptScaleColors
// Be aware : ptWaveletDenoising has as side effect that
// m_Blacklevel and m_Whitelevel change.
//
////////////////////////////////////////////////////////////////////////////////

void CLASS ptScaleColors() {
  if (m_UserSetting_Multiplier[0]) {
    array_copy(m_UserSetting_Multiplier,
               m_rawProcessor->imgdata.params.user_mul);
  }

  if (m_UserSetting_UseGreyBox) {
    array_copy(m_UserSetting_GreyBox, m_rawProcessor->imgdata.params.greybox);
  }

  m_rawProcessor->imgdata.params.use_auto_wb = m_UserSetting_AutoWb;
  m_rawProcessor->imgdata.params.use_camera_wb = m_UserSetting_CameraWb;
  m_rawProcessor->imgdata.params.user_sat = m_UserSetting_Saturation;
  m_rawProcessor->imgdata.params.threshold = m_UserSetting_DenoiseThreshold;

  auto libraw_int_ptr = m_rawProcessor->imgdata.image;
  m_rawProcessor->imgdata.image = m_Image;
  finally f([&] { m_rawProcessor->imgdata.image = libraw_int_ptr; });

  m_rawProcessor->scale_colors();

  array_copy(m_rawProcessor->imgdata.color.pre_mul, m_PreMultipliers);
  for (int c = 0; c < 4; c++) {
    ASSIGN(m_Multipliers[c], VALUE(m_PreMultipliers[c]) * 0xffff /
                                 (m_WhiteLevel - m_CBlackLevel[c]));
  }

  m_RawColorPhotivo = m_RawColor;
}

////////////////////////////////////////////////////////////////////////////////
//
// ptHighlight
// Clipping is determined by Clip... FIXME
//
////////////////////////////////////////////////////////////////////////////////

void CLASS ptHighlight(const short  ClipMode,
                       const short  ClipParameter) {

  TRACEKEYVALS("ClipMode","%d",ClipMode);
  TRACEKEYVALS("ClipParameter","%d",ClipParameter);

  TRACEKEYVALS("m_Multi[0]","%f",VALUE(m_Multipliers[0]));
  TRACEKEYVALS("m_Multi[1]","%f",VALUE(m_Multipliers[1]));
  TRACEKEYVALS("m_Multi[2]","%f",VALUE(m_Multipliers[2]));
  TRACEKEYVALS("m_Multi[3]","%f",VALUE(m_Multipliers[3]));
#ifdef TRACE_ORIGIN
  TRACEKEYVALS("Mult File","%s",m_Multipliers[0].File);
  TRACEKEYVALS("Mult Line","%d",m_Multipliers[0].Line);
#endif

  TRACEKEYVALS("m_MinPreMulti","%f",m_MinPreMulti);
  short ColorRGB[4]={0,1,2,1};
#pragma omp parallel for schedule(static) default(shared)
  for (uint16_t Row = 0; Row < m_OutHeight; Row++) {
    for (uint16_t Column = 0; Column < m_OutWidth; Column++) {
      uint32_t Pos = Row*m_OutWidth+Column;

      // Do take into account however that we scaled the image here
      // already !

      // Clip occurs if the sensor reached its max value. Full stop.
      // m_Image[Pos] stands for value after application of the
      // m_Multipliers, which are per construction so that the saturation
      // value of the sensor maps onto 0XFFFF. It is the unclipped pixel.

      short Clipped = 0;
      for (short Color = 0; Color < m_Colors; Color++) {
        if (m_Image[Pos][Color] >=
            (uint16_t)
                ((m_WhiteLevel-m_CBlackLevel[ColorRGB[Color]])*VALUE(m_Multipliers[Color]))) {
          Clipped = 1;
        }
      }

      if (Clipped) {
        // Here it becomes fun. The sensor reached a clipped value.

        // 'ClippedPixel' stands for the pixel which is multiplied
        // by a value that guarantees the saturation value of the
        // least scaled channel maps onto 0XFFFF.
        // That value is obtained by observing that normally the
        // least scaled channel maps on saturation*minimum_multiplier
        // So since maximum_multipliers brings us per definition on 0XFFFF
        // upscaling with maximum_multiplier/minimum_multipier is that
        // value. Via the equivalence with the pre_multipliers where the
        // maximum is per construction 1 , it means that we have to upscale
        // with 1/m_MinPreMulti.
        // That way all saturated pixels are definitely mapped onto 0xFFFF
        uint16_t ClippedPixel[4];
        for (short Color = 0; Color < m_Colors; Color++) {
          // This ensures that the channel with the smallest multiplier
          // is clipped at its saturation level.
          ClippedPixel[Color] =
            MIN(0xFFFF,
                (int32_t)(m_Image[Pos][Color]/m_MinPreMulti));
          // And now we correct it back for the increased exposure.
          // (but clipped stays clipped !)
          ClippedPixel[Color] = (uint16_t)(ClippedPixel[Color]* m_MinPreMulti);
        }

        // From here on there are now different strategies with respect
        // to clipping.

  // Try to remove purple highlights
  if (0)
    if (m_Image[Pos][2]>m_Image[Pos][1] && m_Image[Pos][0]==0xffff) m_Image[Pos][2]=m_Image[Pos][1];

        // Simply use the clipped value.
        if (ptClipMode_Clip == ClipMode) {
          for (short Color = 0; Color < m_Colors; Color++) {
            m_Image[Pos][Color] = ClippedPixel[Color];
          }

        // Or use a value starting from Unclipped->Clipped ,
        // defined by the ClipParameter.
        } else if (ptClipMode_NoClip == ClipMode) {
          for (short Color = 0; Color < m_Colors; Color++) {
            m_Image[Pos][Color] =
              m_Image[Pos][Color] +
               (uint16_t) (ClipParameter/100.0*
                           (ClippedPixel[Color]-m_Image[Pos][Color]));
          }

        // Or restore via Lab, which is basically the same as
        // used in ufraw (and a simplification of Cyril Guyots LCH blending.
        } else if (ptClipMode_Lab == ClipMode) {
          double ClippedLab[3];
          double UnclippedLab[3];
          double FinalLab[3];
          CamToLab(ClippedPixel,ClippedLab);
          CamToLab(m_Image[Pos],UnclippedLab);
          // FIXME / TODO : Clarify and explain.
          // This is a bit bizar in fact (but also in ufraw ...)
          // The result seems to be better when taking the clipped
          // alternatives for ab and unclipped for L.
          // Wouldn't we expect it the other way around ?
          FinalLab[0] = UnclippedLab[0] +
                        ClipParameter/100.0*
                        (ClippedLab[0]-UnclippedLab[0]);
          FinalLab[1] = ClippedLab[1];
          FinalLab[2] = ClippedLab[2];
          LabToCam(FinalLab,m_Image[Pos]);

        // Or restore via HSV, as in ufraw.
        } else if (ptClipMode_HSV == ClipMode) {
          // FIXME / TODO : can this not break in some 4 colour modes ?
          short MaxChannel,MidChannel,MinChannel;
          if (m_Image[Pos][0] > m_Image[Pos][1] &&
              m_Image[Pos][0] > m_Image[Pos][2]) {
            MaxChannel = 0;
      if (m_Image[Pos][1] > m_Image[Pos][2]) {
              MidChannel = 1;
              MinChannel = 2;
            } else {
              MidChannel = 2;
              MinChannel = 1;
            }
          } else if (m_Image[Pos][1] > m_Image[Pos][2]) {
            MaxChannel = 1;
      if (m_Image[Pos][0] > m_Image[Pos][2]) {
              MidChannel = 0;
              MinChannel = 2;
            } else {
              MidChannel = 2;
              MinChannel = 0;
            }
          } else {
            MaxChannel = 2;
      if (m_Image[Pos][0] > m_Image[Pos][1]) {
              MidChannel = 0;
              MinChannel = 1;
            } else {
              MidChannel = 1;
              MinChannel = 0;
            }
          }

          uint16_t UnclippedLuminance = m_Image[Pos][MaxChannel];
          uint16_t ClippedLuminance   = ClippedPixel[MaxChannel];
          double   ClippedSaturation;
          if ( ClippedPixel[MaxChannel]<ClippedPixel[MinChannel] ||
               ClippedPixel[MaxChannel]==0) {
            ClippedSaturation = 0;
    } else {
            ClippedSaturation =
              1.0 - (double)ClippedPixel[MinChannel] / ClippedPixel[MaxChannel];
          }
//warning: variable 'ClippedHue' set but not used [-Wunused-but-set-variable]
//          double ClippedHue;
//    if ( ClippedPixel[MaxChannel]==ClippedPixel[MinChannel] ) {
//            ClippedHue = 0;
//    } else {
//            ClippedHue =
//              ((double)ClippedPixel[MidChannel]-ClippedPixel[MinChannel]) /
//              ((double)ClippedPixel[MaxChannel]-ClippedPixel[MinChannel]);
//          }
          double UnclippedHue;
    if ( m_Image[Pos][MaxChannel]==m_Image[Pos][MinChannel] ) {
            UnclippedHue = 0;
    } else {
            UnclippedHue =
              ((double)m_Image[Pos][MidChannel]-m_Image[Pos][MinChannel]) /
              ((double)m_Image[Pos][MaxChannel]-m_Image[Pos][MinChannel]);
          }
          uint16_t Luminance =
            UnclippedLuminance +
            (uint16_t)(ClipParameter/100.0 *
                       (ClippedLuminance-UnclippedLuminance));
          double Saturation = ClippedSaturation;
          double Hue = UnclippedHue;
    m_Image[Pos][MaxChannel] = Luminance;
    m_Image[Pos][MinChannel] = (uint16_t)(Luminance * (1-Saturation));
    m_Image[Pos][MidChannel] =
            (uint16_t)(Luminance * (1-Saturation + Saturation*Hue));

        } else if (ptClipMode_Blend == ClipMode) {
          // Do nothing at this stage, keep the unclipped image.
          ;
        } else if (ptClipMode_Rebuild == ClipMode) {
          // Do nothing at this stage, keep the unclipped image.
          ;
        } else {
          assert(0); // should not occur.
        }
      }
    }
  }

  if (ptClipMode_Rebuild == ClipMode) {
    ptRebuildHighlights(ClipParameter);
  }
  if (ptClipMode_Blend == ClipMode) {
    ptBlendHighlights();
  }
}

void CLASS pre_interpolate()
{
  int row, col;

  TRACEKEYVALS("pre_interpolate","%s","");

  if (m_Shrink) {
    if (m_UserSetting_HalfSize) {
      // No interpolation will be needed as
      // m_Shrink has caused a 2X2 Bayer area mapped onto one pixel.
      m_Height = m_OutHeight;
      m_Width  = m_OutWidth;
    } else {
      // in photivo we assume that m_Shrink is only set due
      // to m_UserSetting_HalfSize.
      assert(0);
    }
  }
  if (m_Filters > 1000 && m_Colors == 3) {
    if (m_MixGreen) { // 4 color demosaicer will follow
      m_Colors++;
      // Change from dcraw 1.445 to 1.447 wanted "m_MixGreen = !m_UserSetting_HalfSize;"
      // but this doesn't work in Photivo, since we don't run the full pipe of dcraw,
      // since most of the time we start with phase2
    } else {
      // RG1BG2 -> RGB
#pragma omp parallel for schedule(static) default(shared) private(row, col)
      for (row = FC(1,0) >> 1; row < m_Height; row+=2)
        for (col = FC(row,1) & 1; col < m_Width; col+=2)
          m_Image[row*m_Width+col][1] = m_Image[row*m_Width+col][3];
      m_Filters &= ~((m_Filters & 0x55555555) << 1);
    }
  }

  // If m_UserSetting_HalfSize is set no interpolation will
  // be needed. This is registered by m_Filters = 0.
  // (no Bayer array anymore, this means)
  if (m_UserSetting_HalfSize && m_Filters != 2) {
    m_Filters = 0;
  }
}

void CLASS border_interpolate (int border)
{
  unsigned row, col, y, x, f, c, sum[8];

  for (row=0; row < m_Height; row++)
    for (col=0; col < m_Width; col++) {
      if (col==(unsigned) border && row >= (unsigned) border && row < (unsigned) (m_Height-border))
  col = m_Width-border;
      memset (sum, 0, sizeof sum);
      for (y=row-1; y != row+2; y++)
  for (x=col-1; x != col+2; x++)
    if (y < m_Height && x < m_Width) {
    f = fcol(y,x);
      sum[f] += m_Image[y*m_Width+x][f];
      sum[f+4]++;
    }
      f = fcol(row,col);
      for (c=0; c < (unsigned) m_Colors; c++) if (c != f && sum[c+4])
  m_Image[row*m_Width+col][c] = sum[c] / sum[c+4];
    }
}


void CLASS lin_interpolate()
{
  int code[16][16][32], size=16, *ip, sum[4];
  int f, c, i, x, y, row, col, shift, color;
  uint16_t *pix;

  TRACEKEYVALS("Bilinear interpolation","%s","");
  if (m_Filters == 2) size = 6;
  border_interpolate(1);
  for (row=0; row < size; row++)
    for (col=0; col < size; col++) {
      ip = code[row][col]+1;
      f = fcol(row,col);
      memset (sum, 0, sizeof sum);
      for (y=-1; y <= 1; y++)
  for (x=-1; x <= 1; x++) {
    shift = (y==0) + (x==0);
    color = fcol(row+y,col+x);
    if (color == f) continue;
    *ip++ = (m_Width*y + x)*4 + color;
    *ip++ = shift;
    *ip++ = color;
    sum[color] += 1 << shift;
  }
      code[row][col][0] = (ip - code[row][col]) / 3;
      for (c=0; c < m_Colors; c++)
  if (c != f) {
    *ip++ = c;
    *ip++ = 256 / sum[c];
  }
    }
#pragma omp parallel for schedule(static) default(shared) private(row,col,pix,ip,sum,i)
  for (row=1; row < m_Height-1; row++)
    for (col=1; col < m_Width-1; col++) {
      pix = m_Image[row*m_Width+col];
      ip = code[row % size][col % size];
      memset (sum, 0, sizeof sum);
      for (i=*ip++; i--; ip+=3)
  sum[ip[2]] += pix[ip[0]] << ip[1];
      for (i=m_Colors; --i; ip+=2)
  pix[ip[0]] = sum[ip[0]] * ip[1] >> 8;
    }
}

/*
   This algorithm is officially called:

   "Interpolation using a Threshold-based variable number of gradients"

   described in http://scien.stanford.edu/pages/labsite/1999/psych221/projects/99/tingchen/algodep/vargra.html

   I've extended the basic idea to work with non-Bayer filter arrays.
   Gradients are numbered clockwise from NW=0 to W=7.
 */
void CLASS vng_interpolate()
{
  static const int16_t *cp, terms[] = {
    -2,-2,+0,-1,0,0x01, -2,-2,+0,+0,1,0x01, -2,-1,-1,+0,0,0x01,
    -2,-1,+0,-1,0,0x02, -2,-1,+0,+0,0,0x03, -2,-1,+0,+1,1,0x01,
    -2,+0,+0,-1,0,0x06, -2,+0,+0,+0,1,0x02, -2,+0,+0,+1,0,0x03,
    -2,+1,-1,+0,0,0x04, -2,+1,+0,-1,1,0x04, -2,+1,+0,+0,0,0x06,
    -2,+1,+0,+1,0,0x02, -2,+2,+0,+0,1,0x04, -2,+2,+0,+1,0,0x04,
    -1,-2,-1,+0,0,0x80, -1,-2,+0,-1,0,0x01, -1,-2,+1,-1,0,0x01,
    -1,-2,+1,+0,1,0x01, -1,-1,-1,+1,0,0x88, -1,-1,+1,-2,0,0x40,
    -1,-1,+1,-1,0,0x22, -1,-1,+1,+0,0,0x33, -1,-1,+1,+1,1,0x11,
    -1,+0,-1,+2,0,0x08, -1,+0,+0,-1,0,0x44, -1,+0,+0,+1,0,0x11,
    -1,+0,+1,-2,1,0x40, -1,+0,+1,-1,0,0x66, -1,+0,+1,+0,1,0x22,
    -1,+0,+1,+1,0,0x33, -1,+0,+1,+2,1,0x10, -1,+1,+1,-1,1,0x44,
    -1,+1,+1,+0,0,0x66, -1,+1,+1,+1,0,0x22, -1,+1,+1,+2,0,0x10,
    -1,+2,+0,+1,0,0x04, -1,+2,+1,+0,1,0x04, -1,+2,+1,+1,0,0x04,
    +0,-2,+0,+0,1,0x80, +0,-1,+0,+1,1,0x88, +0,-1,+1,-2,0,0x40,
    +0,-1,+1,+0,0,0x11, +0,-1,+2,-2,0,0x40, +0,-1,+2,-1,0,0x20,
    +0,-1,+2,+0,0,0x30, +0,-1,+2,+1,1,0x10, +0,+0,+0,+2,1,0x08,
    +0,+0,+2,-2,1,0x40, +0,+0,+2,-1,0,0x60, +0,+0,+2,+0,1,0x20,
    +0,+0,+2,+1,0,0x30, +0,+0,+2,+2,1,0x10, +0,+1,+1,+0,0,0x44,
    +0,+1,+1,+2,0,0x10, +0,+1,+2,-1,1,0x40, +0,+1,+2,+0,0,0x60,
    +0,+1,+2,+1,0,0x20, +0,+1,+2,+2,0,0x10, +1,-2,+1,+0,0,0x80,
    +1,-1,+1,+1,0,0x88, +1,+0,+1,+2,0,0x08, +1,+0,+2,-1,0,0x40,
    +1,+0,+2,+1,0,0x10
  }, chood[] = { -1,-1, -1,0, -1,+1, 0,+1, +1,+1, +1,0, +1,-1, 0,-1 };
  uint16_t (*brow[5])[4], *pix;
  int prow=8, pcol=2, *ip, *code[16][16], gval[8], gmin, gmax, sum[4];
  int row, col, x, y, x1, x2, y1, y2, t, weight, grads, color, diag;
  int g, diff, thold, num, c;

  lin_interpolate();

  TRACEKEYVALS("VNG interpolation","%s","");

  if (m_Filters == 1) prow = pcol = 16;
  if (m_Filters == 2) prow = pcol =  6;
  ip = (int *) CALLOC (prow*pcol, 1280);
  merror (ip, "vng_interpolate()");
  for (row=0; row < prow; row++)		/* Precalculate for VNG */
    for (col=0; col < pcol; col++) {
      code[row][col] = ip;
      for (cp=terms, t=0; t < 64; t++) {
  y1 = *cp++;  x1 = *cp++;
  y2 = *cp++;  x2 = *cp++;
  weight = *cp++;
  grads = *cp++;
  color = fcol(row+y1,col+x1);
  if (fcol(row+y2,col+x2) != color) continue;
  diag = (fcol(row,col+1) == color && fcol(row+1,col) == color) ? 2:1;
  if (abs(y1-y2) == diag && abs(x1-x2) == diag) continue;
  *ip++ = (y1*m_Width + x1)*4 + color;
  *ip++ = (y2*m_Width + x2)*4 + color;
  *ip++ = weight;
  for (g=0; g < 8; g++)
    if (grads & 1<<g) *ip++ = g;
  *ip++ = -1;
      }
      *ip++ = INT_MAX;
      for (cp=chood, g=0; g < 8; g++) {
  y = *cp++;  x = *cp++;
  *ip++ = (y*m_Width + x) * 4;
  color = fcol(row,col);
  if (fcol(row+y,col+x) != color && fcol(row+y*2,col+x*2) == color)
    *ip++ = (y*m_Width + x) * 8 + color;
  else
    *ip++ = 0;
      }
    }

  brow[4] = (uint16_t (*)[4]) CALLOC (m_Width*3, sizeof **brow);
  merror (brow[4], "vng_interpolate()");
  for (row=0; row < 3; row++)
    brow[row] = brow[4] + row*m_Width;

  for (row=2; row < m_Height-2; row++) {    /* Do VNG interpolation */
    for (col=2; col < m_Width-2; col++) {
      pix = m_Image[row*m_Width+col];
      ip = code[row % prow][col % pcol];
      memset (gval, 0, sizeof gval);
      while ((g = ip[0]) != INT_MAX) {    /* Calculate gradients */
  diff = ABS(pix[g] - pix[ip[1]]) << ip[2];
  gval[ip[3]] += diff;
  ip += 5;
  if ((g = ip[-1]) == -1) continue;
  gval[g] += diff;
  while ((g = *ip++) != -1)
    gval[g] += diff;
      }
      ip++;
      gmin = gmax = gval[0];      /* Choose a threshold */
      for (g=1; g < 8; g++) {
  if (gmin > gval[g]) gmin = gval[g];
  if (gmax < gval[g]) gmax = gval[g];
      }
      if (gmax == 0) {
  memcpy (brow[2][col], pix, sizeof *m_Image);
  continue;
      }
      thold = gmin + (gmax >> 1);
      memset (sum, 0, sizeof sum);
      color = fcol(row,col);
      for (num=g=0; g < 8; g++,ip+=2) {   /* Average the neighbors */
  if (gval[g] <= thold) {
    for (c=0; c < m_Colors; c++)
      if (c == color && ip[1])
        sum[c] += (pix[c] + pix[ip[1]]) >> 1;
      else
        sum[c] += pix[ip[0] + c];
    num++;
  }
      }
      for (c=0; c < m_Colors; c++) {          /* Save to buffer */
  t = pix[color];
  if (c != color)
    t += (sum[c] - sum[color]) / num;
  brow[2][col][c] = CLIP(t);
      }
    }
    if (row > 3)        /* Write buffer to image */
      memcpy (m_Image[(row-2)*m_Width+2], brow[0]+2, (m_Width-4)*sizeof *m_Image);
    for (g=0; g < 4; g++)
      brow[(g-1) & 3] = brow[g];
  }
  memcpy (m_Image[(row-2)*m_Width+2], brow[0]+2, (m_Width-4)*sizeof *m_Image);
  memcpy (m_Image[(row-1)*m_Width+2], brow[1]+2, (m_Width-4)*sizeof *m_Image);
  FREE (brow[4]);
  FREE (code[0][0]);
}

/*
   Patterned Pixel Grouping Interpolation by Alain Desbiolles
*/
void CLASS ppg_interpolate()
{
  const int dir[5] = { 1, m_Width, -1, -m_Width, 1 };
  int row, col, diff[2], guess[2], c, d, i;
  uint16_t (*pix)[4];

  border_interpolate(3);

  TRACEKEYVALS("PPG interpolation","%s","");

#pragma omp parallel private(row,c,col,i,guess,diff,d,pix)
/*  Fill in the green layer with gradients and pattern recognition: */
#pragma omp for
  for (row=3; row < m_Height-3; row++)
    for (col=3+(FC(row,3) & 1), c=FC(row,col); col < m_Width-3; col+=2) {
      pix = m_Image + row*m_Width+col;
      for (i=0; i<2; i++) {
        d = dir[i];
  guess[i] = (pix[-d][1] + pix[0][c] + pix[d][1]) * 2
          - pix[-2*d][c] - pix[2*d][c];
  diff[i] = ( ABS(pix[-2*d][c] - pix[ 0][c]) +
        ABS(pix[ 2*d][c] - pix[ 0][c]) +
        ABS(pix[  -d][1] - pix[ d][1]) ) * 3 +
      ( ABS(pix[ 3*d][1] - pix[ d][1]) +
        ABS(pix[-3*d][1] - pix[-d][1]) ) * 2;
      }
      d = dir[i = diff[0] > diff[1]];
      pix[0][1] = ULIM(guess[i] >> 2, (int32_t)pix[d][1], (int32_t)pix[-d][1]);
    }
/*  Calculate red and blue for each green pixel:    */
#pragma omp for
  for (row=1; row < m_Height-1; row++)
    for (col=1+(FC(row,2) & 1), c=FC(row,col+1); col < m_Width-1; col+=2) {
      pix = m_Image + row*m_Width+col;
      for (i=0; (d=dir[i]) > 0; c=2-c, i++)
  pix[0][c] = CLIP((pix[-d][c] + pix[d][c] + 2*pix[0][1]
      - pix[-d][1] - pix[d][1]) >> 1);
    }
/*  Calculate blue for red pixels and vice versa:   */
#pragma omp for
  for (row=1; row < m_Height-1; row++)
    for (col=1+(FC(row,1) & 1), c=2-FC(row,col); col < m_Width-1; col+=2) {
      pix = m_Image + row*m_Width+col;
      for (i=0; i < 2; i++) {
        d = dir[i]+dir[i+1];
        diff[i] = ABS(pix[-d][c] - pix[d][c]) +
      ABS(pix[-d][1] - pix[0][1]) +
      ABS(pix[ d][1] - pix[0][1]);
  guess[i] = pix[-d][c] + pix[d][c] + 2*pix[0][1]
     - pix[-d][1] - pix[d][1];
      }
      if (diff[0] != diff[1])
  pix[0][c] = CLIP(guess[diff[0] > diff[1]] >> 1);
      else
  pix[0][c] = CLIP((guess[0]+guess[1]) >> 2);
    }
}

/*
   Adaptive Homogeneity-Directed interpolation is based on
   the work of Keigo Hirakawa, Thomas Parks, and Paul Lee.
 */
#define TS 256    /* Tile Size */

void CLASS ahd_interpolate()
{
  int i, j, k, top, left, row, col, tr, tc, c, d, val, hm[2];
  uint16_t (*pix)[4], (*rix)[3];
  static const int dir[4] = { -1, 1, -TS, TS };
  unsigned ldiff[2][4], abdiff[2][4], leps, abeps;
  float r, cbrt[0x10000], xyz[3], xyz_cam[3][4];
  uint16_t (*rgb)[TS][TS][3];
   short (*lab)[TS][TS][3], (*lix)[3];
   char (*homo)[TS][TS], *buffer;

  TRACEKEYVALS("AHD interpolation","%s","");

  for (i=0; i < 0x10000; i++) {
    r = i / 65535.0;
    cbrt[i] = r > 0.008856 ? pow(r,1/3.0) : 7.787*r + 16/116.0;
  }
  for (i=0; i < 3; i++)
    for (j=0; j < m_Colors; j++)
      for (xyz_cam[i][j] = k=0; k < 3; k++)
  xyz_cam[i][j] += MatrixRGBToXYZ[ptSpace_sRGB_D65][i][k] * m_MatrixCamRGBToSRGB[k][j] / D65Reference[i];

  border_interpolate(5);

#pragma omp parallel private(buffer,rgb,lab,homo,top,left,row,c,col,pix,val,d,rix,xyz,lix,tc,tr,ldiff,abdiff,leps,abeps,hm,i,j) firstprivate(cbrt) shared(xyz_cam)
  {
      buffer = (char *) MALLOC (26*TS*TS);    /* 1664 kB */
      merror (buffer, "ahd_interpolate()");
      rgb  = (uint16_t(*)[TS][TS][3]) buffer;
      lab  = (short (*)[TS][TS][3])(buffer + 12*TS*TS);
      homo = (char  (*)[TS][TS])   (buffer + 24*TS*TS);

#pragma omp for schedule(static)
      for (top=2; top < m_Height-5; top += TS-6)
          for (left=2; left < m_Width-5; left += TS-6) {

          /*  Interpolate green horizontally and vertically:    */
          for (row = top; row < top+TS && row < m_Height-2; row++) {
              col = left + (FC(row,left) & 1);
              for (c = FC(row,col); col < left+TS && col < m_Width-2; col+=2) {
                  pix = m_Image + row*m_Width+col;
                  val = ((pix[-1][1] + pix[0][c] + pix[1][1]) * 2
                         - pix[-2][c] - pix[2][c]) >> 2;
                  rgb[0][row-top][col-left][1] = ULIM(val,(int32_t)pix[-1][1],(int32_t)pix[1][1]);
                  val = ((pix[-m_Width][1] + pix[0][c] + pix[m_Width][1]) * 2
                         - pix[-2*m_Width][c] - pix[2*m_Width][c]) >> 2;
                  rgb[1][row-top][col-left][1] = ULIM(val,(int32_t)pix[-m_Width][1],(int32_t)pix[m_Width][1]);
              }
          }
          /*  Interpolate red and blue, and convert to CIELab:    */
          for (d=0; d < 2; d++)
              for (row=top+1; row < top+TS-1 && row < m_Height-3; row++)
                  for (col=left+1; col < left+TS-1 && col < m_Width-3; col++) {
              pix = m_Image + row*m_Width+col;
              rix = &rgb[d][row-top][col-left];
              lix = &lab[d][row-top][col-left];
              if ((c = 2 - FC(row,col)) == 1) {
                  c = FC(row+1,col);
                  val = pix[0][1] + (( pix[-1][2-c] + pix[1][2-c]
                                       - rix[-1][1] - rix[1][1] ) >> 1);
                  rix[0][2-c] = CLIP(val);
                  val = pix[0][1] + (( pix[-m_Width][c] + pix[m_Width][c]
                                       - rix[-TS][1] - rix[TS][1] ) >> 1);
              } else
                  val = rix[0][1] + (( pix[-m_Width-1][c] + pix[-m_Width+1][c]
                                       + pix[+m_Width-1][c] + pix[+m_Width+1][c]
                                       - rix[-TS-1][1] - rix[-TS+1][1]
                                       - rix[+TS-1][1] - rix[+TS+1][1] + 1) >> 2);
              rix[0][c] = CLIP(val);
              c = FC(row,col);
              rix[0][c] = pix[0][c];
              xyz[0] = xyz[1] = xyz[2] = 0.5;
              for (c=0; c < m_Colors; c++) {
                  xyz[0] += xyz_cam[0][c] * rix[0][c];
                  xyz[1] += xyz_cam[1][c] * rix[0][c];
                  xyz[2] += xyz_cam[2][c] * rix[0][c];
              }
              xyz[0] = cbrt[CLIP((int) xyz[0])];
              xyz[1] = cbrt[CLIP((int) xyz[1])];
              xyz[2] = cbrt[CLIP((int) xyz[2])];
              lix[0][0] = (short) (64 * (116 * xyz[1] - 16));
              lix[0][1] = (short) (64 * 500 * (xyz[0] - xyz[1]));
              lix[0][2] = (short) (64 * 200 * (xyz[1] - xyz[2]));
          }
          /*  Build homogeneity maps from the CIELab images:    */
          memset (homo, 0, 2*TS*TS);
          for (row=top+2; row < top+TS-2 && row < m_Height-4; row++) {
              tr = row-top;
              for (col=left+2; col < left+TS-2 && col < m_Width-4; col++) {
                  tc = col-left;
                  for (d=0; d < 2; d++) {
                      lix = &lab[d][tr][tc];
                      for (i=0; i < 4; i++) {
                          ldiff[d][i] = ABS(lix[0][0]-lix[dir[i]][0]);
                          abdiff[d][i] = SQR(lix[0][1]-lix[dir[i]][1])
                                         + SQR(lix[0][2]-lix[dir[i]][2]);
                      }
                  }
                  leps = MIN(MAX(ldiff[0][0],ldiff[0][1]),
                             MAX(ldiff[1][2],ldiff[1][3]));
                  abeps = MIN(MAX(abdiff[0][0],abdiff[0][1]),
                              MAX(abdiff[1][2],abdiff[1][3]));
                  for (d=0; d < 2; d++)
                      for (i=0; i < 4; i++)
                          if (ldiff[d][i] <= leps && abdiff[d][i] <= abeps)
                             homo[d][tr][tc]++;
              }
          }
          /*  Combine the most homogenous pixels for the final result:  */
          for (row=top+3; row < top+TS-3 && row < m_Height-5; row++) {
              tr = row-top;
              for (col=left+3; col < left+TS-3 && col < m_Width-5; col++) {
                  tc = col-left;
                  for (d=0; d < 2; d++)
                      for (hm[d]=0, i=tr-1; i <= tr+1; i++)
                          for (j=tc-1; j <= tc+1; j++)
                              hm[d] += homo[d][i][j];
                  if (hm[0] != hm[1])
                      for (c=0; c<3; c++) m_Image[row*m_Width+col][c] = rgb[hm[1] > hm[0]][tr][tc][c];
                  else
                      for (c=0; c<3; c++) m_Image[row*m_Width+col][c] =
                              (rgb[0][tr][tc][c] + rgb[1][tr][tc][c]) >> 1;
              }
          }
      }
      FREE (buffer);
  }
}
#undef TS

void CLASS fuji_rotate() {
  int i, row, col;
  double step;
  float r, c, fr, fc;
  unsigned ur, uc;
  uint16_t wide, high, (*img)[4], (*pix)[4];

  if (!m_Fuji_Width) return;
  TRACEKEYVALS("Rotating image 45 degrees","%s","");
  TRACEKEYVALS("m_Fuji_Width","%d",m_Fuji_Width);
  TRACEKEYVALS("m_Shrink","%d",m_Shrink);
  TRACEKEYVALS("m_UserSetting_HalfSize","%d",m_UserSetting_HalfSize);
  // XXX JDLA : Was : m_Fuji_Width = (m_Fuji_Width - 1 + m_Shrink) >> m_Shrink;
  m_Fuji_Width =
    (m_Fuji_Width - 1 + m_UserSetting_HalfSize) >> m_UserSetting_HalfSize;
  TRACEKEYVALS("m_Fuji_Width","%d",m_Fuji_Width);
  step = sqrt(0.5);
  wide = (uint16_t) (m_Fuji_Width / step);
  high = (uint16_t) ((m_Height - m_Fuji_Width) / step);
  img = (uint16_t (*)[4]) CALLOC (wide*high, sizeof *img);
  TRACEKEYVALS("m_Width","%d",m_Width);
  TRACEKEYVALS("m_Height","%d",m_Height);
  TRACEKEYVALS("wide","%d",wide);
  TRACEKEYVALS("high","%d",high);
  merror (img, "fuji_rotate()");

  for (row=0; row < high; row++)
    for (col=0; col < wide; col++) {
      ur = (unsigned) (r = m_Fuji_Width + (row-col)*step);
      uc = (unsigned) (c = (row+col)*step);
      if (ur > (unsigned)(m_Height-2) || uc > (unsigned)(m_Width-2)) continue;
      fr = r - ur;
      fc = c - uc;
      pix = m_Image + ur*m_Width + uc;
      for (i=0; i < m_Colors; i++)
  img[row*wide+col][i] = (uint16_t) (
    (pix[    0][i]*(1-fc) + pix[      1][i]*fc) * (1-fr) +
    (pix[m_Width][i]*(1-fc) + pix[m_Width+1][i]*fc) * fr);
    }
  FREE (m_Image);
  m_OutWidth  = m_Width  = wide;
  m_OutHeight = m_Height = high;
  m_Image  = img;
  m_Fuji_Width = 0; // this prevents image rotation when only phase 2 is called

  TRACEKEYVALS("m_Width","%d",m_Width);
  TRACEKEYVALS("m_Height","%d",m_Height);
  TRACEKEYVALS("m_OutWidth","%d",m_OutWidth);
  TRACEKEYVALS("m_OutHeight","%d",m_OutHeight);
}

void CLASS stretch()
{
  uint16_t newdim, (*img)[4], *pix0, *pix1;
  int row, col, c;
  double rc, frac;

  if (m_PixelAspect == 1) return;
  TRACEKEYVALS("Stretching the image","%s","");
  TRACEKEYVALS("m_PixelAspect","%f",m_PixelAspect);
  if (m_PixelAspect < 1) {
    newdim = (uint16_t) (m_Height / m_PixelAspect + 0.5);
    img = (uint16_t (*)[4]) CALLOC (m_Width*newdim, sizeof *img);
    merror (img, "stretch()");
    for (rc=row=0; row < newdim; row++, rc+=m_PixelAspect) {
      frac = rc - (c = (int) rc);
      pix0 = pix1 = m_Image[c*m_Width];
      if (c+1 < m_Height) pix1 += m_Width*4;
      for (col=0; col < m_Width; col++, pix0+=4, pix1+=4)
  for (c=0; c < m_Colors; c++) img[row*m_Width+col][c] =
          (uint16_t)(pix0[c]*(1-frac) + pix1[c]*frac + 0.5);
    }
    m_OutHeight = m_Height = newdim;
  } else {
    newdim = (uint16_t) (m_Width * m_PixelAspect + 0.5);
    img = (uint16_t (*)[4]) CALLOC (m_Height*newdim, sizeof *img);
    merror (img, "stretch()");
    for (rc=col=0; col < newdim; col++, rc+=1/m_PixelAspect) {
      frac = rc - (c = (int) rc);
      pix0 = pix1 = m_Image[c];
      if (c+1 < m_Width) pix1 += 4;
      for (row=0; row < m_Height; row++, pix0+=m_Width*4, pix1+=m_Width*4)
  for (c=0; c < m_Colors; c++) img[row*newdim+col][c] =
         (uint16_t) (pix0[c]*(1-frac) + pix1[c]*frac + 0.5);
    }
    m_OutWidth = m_Width = newdim;
  }
  FREE (m_Image);
  m_Image = img;

  TRACEKEYVALS("m_Width","%d",m_Width);
  TRACEKEYVALS("m_Height","%d",m_Height);
  TRACEKEYVALS("m_OutWidth","%d",m_OutWidth);
  TRACEKEYVALS("m_OutHeight","%d",m_OutHeight);
}

////////////////////////////////////////////////////////////////////////////////
//
// Identify
//
// First method to be called when in an photivo application.
// Hereafter camera, model, and a number of other picture and camera
// parameters are identified.
// Returns 0 on success.
//
////////////////////////////////////////////////////////////////////////////////

short CLASS Identify(const QString &NewInputFile) {

  // This is here to support multiple calls.
  ResetNonUserSettings();

  //if (!NewInputFile.isEmpty()) {
  if (NewInputFile != "") {
    m_UserSetting_InputFileName = NewInputFile;
  }

  auto ret = m_rawProcessor->open_file(m_UserSetting_InputFileName);
  if (ret != LIBRAW_SUCCESS)
  {
      perror (m_UserSetting_InputFileName.toLocal8Bit().constData());
      perror (libraw_strerror(ret));
      return -1;
  }

  m_IsRaw = true;
  const auto &S = m_rawProcessor->imgdata.sizes;
  m_Width = S.width;
  m_Height = S.height;
  m_RawWidth = S.raw_width;
  m_RawHeight = S.raw_height;
  m_TopMargin = S.top_margin;
  m_LeftMargin = S.left_margin;
  m_Flip = S.flip;

  m_Filters = m_rawProcessor->imgdata.idata.filters;
  m_Colors = m_rawProcessor->imgdata.idata.colors;
  m_IsFoveon = m_rawProcessor->imgdata.idata.is_foveon;

  memcpy(m_MatrixCamRGBToSRGB, m_rawProcessor->imgdata.color.rgb_cam, sizeof (float) * 3 * 4);
  array_copy(m_rawProcessor->imgdata.color.pre_mul, m_D65Multipliers);
  // TODO: might be not needed any more: m_CameraMultipliers
  array_copy(m_rawProcessor->imgdata.color.cam_mul, m_CameraMultipliers);

  memcpy(m_CBlackLevel, m_rawProcessor->imgdata.color.cblack, sizeof m_CBlackLevel);
  m_BlackLevel = m_rawProcessor->imgdata.color.black;
  m_WhiteLevel = m_rawProcessor->imgdata.color.maximum;

  m_CameraMake = QString::fromUtf8(m_rawProcessor->imgdata.idata.make);
  m_CameraModel = QString::fromUtf8(m_rawProcessor->imgdata.idata.model);

  /* TODO:
     * bring back intensity scaling
     * bring back camera white balances
     * extra camera profiles
     * use libraw demosaicer
  */

  return !m_IsRaw;
}

////////////////////////////////////////////////////////////////////////////////
//
// RunDcRaw_Phase1
//
// Settings are given via the m_UserSetting* members.
// It will end there where the raw image is loaded.
// Badpixels and darkframes subtracted.
// But before colorscaling and the like so that we can
// intervene in the whitebalances etc.
// Returns 0 on success.
//
////////////////////////////////////////////////////////////////////////////////

short CLASS RunDcRaw_Phase1() {

  // TODO foveon for the moment not in. Need study material
  // to have this +/- right.
  assert (!m_IsFoveon);

  TRACEKEYVALS("PreMult[0]","%f",VALUE(m_PreMultipliers[0]));
  TRACEKEYVALS("PreMult[1]","%f",VALUE(m_PreMultipliers[1]));
  TRACEKEYVALS("PreMult[2]","%f",VALUE(m_PreMultipliers[2]));
  TRACEKEYVALS("PreMult[3]","%f",VALUE(m_PreMultipliers[3]));
#ifdef TRACE_ORIGIN
  TRACEKEYVALS("PreMult File","%s",m_PreMultipliers[0].File);
  TRACEKEYVALS("PreMult Line","%d",m_PreMultipliers[0].Line);
#endif
  TRACEKEYVALS("D65Mult[0]","%f",VALUE(m_D65Multipliers[0]));
  TRACEKEYVALS("D65Mult[1]","%f",VALUE(m_D65Multipliers[1]));
  TRACEKEYVALS("D65Mult[2]","%f",VALUE(m_D65Multipliers[2]));
  TRACEKEYVALS("D65Mult[3]","%f",VALUE(m_D65Multipliers[3]));
#ifdef TRACE_ORIGIN
  TRACEKEYVALS("D65Mult File","%s",m_D65Multipliers[0].File);
  TRACEKEYVALS("D65Mult Line","%d",m_D65Multipliers[0].Line);
#endif

  TRACEKEYVALS("RawWidth","%d",m_RawWidth);
  TRACEKEYVALS("RawHeight","%d",m_RawHeight);
  TRACEKEYVALS("TopMargin","%d",m_TopMargin);
  TRACEKEYVALS("LeftMargin","%d",m_LeftMargin);
  TRACEKEYVALS("BlackLevel","%d",m_BlackLevel);

  // OK for second entry.

  m_OutHeight = m_Height;
  m_OutWidth  = m_Width;

  TRACEKEYVALS("OutWidth","%d",m_OutWidth);
  TRACEKEYVALS("Width","%d",m_Width);

  // Also some reshuffling for second entry problem.
  // not sure how c_matrix comes into play here ...
  int l_CameraMatrix = (m_UserSetting_CameraMatrix < 0) ?
     m_UserSetting_CameraWb : m_UserSetting_CameraMatrix;
  TRACEKEYVALS("US_CameraWB","%d",m_UserSetting_CameraWb);
  TRACEKEYVALS("US_CameraMatr","%d",m_UserSetting_CameraMatrix);
  TRACEKEYVALS("CameraMatrix","%d",l_CameraMatrix);
  TRACEKEYVALS("m_cmatrix[0][0","%f",m_cmatrix[0][0]);
  if (l_CameraMatrix && m_cmatrix[0][0] > 0.25) {
    TRACEKEYVALS("Using CamMatr","%s","Yes");
    memcpy (m_MatrixCamRGBToSRGB, m_cmatrix, sizeof m_cmatrix);
    m_RawColor = 0;
  }

  // Allocation is depending on m_Raw_Image below.
  FREE(m_Image);

  TRACEKEYVALS_QT("CameraMake","%s",m_CameraMake);
  TRACEKEYVALS_QT("CameraModel","%s",m_CameraModel);
  TRACEKEYVALS_QT("InputFile","%s",m_UserSetting_InputFileName);
  TRACEKEYVALS("Filters","%x",m_Filters);
  TRACEKEYVALS("Flip","%x",m_Flip);

  auto ret = m_rawProcessor->unpack();
  if (ret != LIBRAW_SUCCESS) {
      perror ("Cannot unpack!");
      perror (libraw_strerror(ret));
  }
  ret = m_rawProcessor->raw2image();
  if (ret != LIBRAW_SUCCESS) {
      perror ("Cannot read raw!");
      perror (libraw_strerror(ret));
  }
  m_Image = (uint16_t (*)[4]) CALLOC (m_OutHeight*m_OutWidth, sizeof (*m_Image));
  merror (m_Image, "raw image");

  memcpy(m_Image,
         m_rawProcessor->imgdata.image,
         m_OutHeight*m_OutWidth*sizeof(*m_Image));

  if (m_ZeroIsBad) remove_zeroes();
  bad_pixels (m_UserSetting_BadPixelsFileName);
  if (m_UserSetting_DarkFrameFileName)
    subtract (m_UserSetting_DarkFrameFileName);

  // photivo extra.
  // Calculation of the inverse m_MatrixCamRGBToSRGB
  double rgb_cam_transpose[4][3];
  for (short i=0; i<4; i++) for (short j=0; j<3; j++)
    rgb_cam_transpose[i][j] = m_MatrixCamRGBToSRGB[j][i];
  pseudoinverse(rgb_cam_transpose,m_MatrixSRGBToCamRGB,m_Colors);

  // ReportedWidth & Height correct for fuji_rotate stuff and pixelaspects.
  m_ReportedWidth  = m_Width;
  m_ReportedHeight = m_Height;
  if (m_Fuji_Width) {
    m_IsFuji = m_Fuji_Width;
    m_ReportedWidth  = (m_Fuji_Width-1) / sqrt(0.5);
    m_ReportedHeight = (m_Height-m_Fuji_Width+1)/sqrt(0.5);
  }
  if (m_PixelAspect<1)
    m_ReportedHeight = (uint16_t) (m_Height / m_PixelAspect + 0.5);
  if (m_PixelAspect>1)
    m_ReportedWidth = (uint16_t) (m_Width * m_PixelAspect + 0.5);

  // TODO Mike: m_ReportedH/W is never set back to m_H/W, CHECK!

  TRACEKEYVALS("m_Width","%d",m_Width);
  TRACEKEYVALS("m_Height","%d",m_Height);
  TRACEKEYVALS("m_ReportedWidth","%d",m_ReportedWidth);
  TRACEKEYVALS("m_ReportedHeight","%d",m_ReportedHeight);

  TRACEKEYVALS("PreMult[0]","%f",VALUE(m_PreMultipliers[0]));
  TRACEKEYVALS("PreMult[1]","%f",VALUE(m_PreMultipliers[1]));
  TRACEKEYVALS("PreMult[2]","%f",VALUE(m_PreMultipliers[2]));
  TRACEKEYVALS("PreMult[3]","%f",VALUE(m_PreMultipliers[3]));
#ifdef TRACE_ORIGIN
  TRACEKEYVALS("PreMult File","%s",m_PreMultipliers[0].File);
  TRACEKEYVALS("PreMult Line","%d",m_PreMultipliers[0].Line);
#endif
  TRACEKEYVALS("D65Mult[0]","%f",VALUE(m_D65Multipliers[0]));
  TRACEKEYVALS("D65Mult[1]","%f",VALUE(m_D65Multipliers[1]));
  TRACEKEYVALS("D65Mult[2]","%f",VALUE(m_D65Multipliers[2]));
  TRACEKEYVALS("D65Mult[3]","%f",VALUE(m_D65Multipliers[3]));
#ifdef TRACE_ORIGIN
  TRACEKEYVALS("D65Mult File","%s",m_D65Multipliers[0].File);
  TRACEKEYVALS("D65Mult Line","%d",m_D65Multipliers[0].Line);
#endif
  TRACEKEYVALS("Colors","%d",m_Colors);
  TRACEKEYVALS("Filters","%x",m_Filters);
  TRACEKEYVALS("Flip","%x",m_Flip);

  // Cache the image after Phase1.
  FREE(m_Image_AfterPhase1);
  m_Image_AfterPhase1 =
    (uint16_t (*)[4]) CALLOC (m_OutHeight*m_OutWidth, sizeof *m_Image);
  merror (m_Image_AfterPhase1, "main()");
  memcpy(m_Image_AfterPhase1,m_Image,m_OutHeight*m_OutWidth*sizeof(*m_Image));
  // Some other stuff to cache.
  m_Filters_AfterPhase1 = m_Filters;
  m_BlackLevel_AfterPhase1 = m_BlackLevel;
  array_copy(m_CBlackLevel, m_CBlackLevel_AfterPhase1);
  m_WhiteLevel_AfterPhase1 = m_WhiteLevel;
  m_Width_AfterPhase1 = m_Width;
  m_Height_AfterPhase1 = m_Height;
  m_OutWidth_AfterPhase1 = m_OutWidth;
  m_OutHeight_AfterPhase1 = m_OutHeight;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// RunDcRaw_Phase2
//
// Color scaling.
// Interpolation.
// Green mixing.
//
////////////////////////////////////////////////////////////////////////////////

short CLASS RunDcRaw_Phase2(const short NoCache) {

  // Make sure we are starting from the right image.
  if (NoCache) {
    FREE(m_Image_AfterPhase1);
  } else {
    FREE(m_Image);
    m_Width = m_Width_AfterPhase1;
    m_Height = m_Height_AfterPhase1;
    m_OutWidth = m_OutWidth_AfterPhase1;
    m_OutHeight = m_OutHeight_AfterPhase1;
    m_Image =
      (uint16_t (*)[4]) CALLOC (m_OutHeight*m_OutWidth, sizeof *m_Image);
    merror (m_Image, "main()");
    memcpy(m_Image,m_Image_AfterPhase1,m_OutHeight*m_OutWidth*sizeof(*m_Image));
    // Restore some other cached values.
    m_Filters    = m_Filters_AfterPhase1;
    m_BlackLevel = m_BlackLevel_AfterPhase1;
    array_copy(m_CBlackLevel_AfterPhase1, m_CBlackLevel);
    m_WhiteLevel = m_WhiteLevel_AfterPhase1;
  }
  m_Fuji_Width = m_IsFuji;
  unsigned i = m_CBlackLevel[3];
  for (int c=0; c<3; c++) if (i > m_CBlackLevel[c]) i = m_CBlackLevel[c];
  for (int c=0; c<4; c++) m_CBlackLevel[c] -= i;
  m_BlackLevel += i;
  if (m_UserSetting_BlackPoint >= 0) m_BlackLevel = m_UserSetting_BlackPoint;
  for (int c=0; c<4; c++) m_CBlackLevel[c] += m_BlackLevel;
  if (m_UserSetting_Saturation > 0)  m_WhiteLevel = m_UserSetting_Saturation;

  // If m_UserSetting_HalfSize then the BAYER(row,col) macro
  // will map 2X2 pixels of the Bayer array, directly onto one
  // pixel of the image (and the correct color channel) foregoing
  // the need for interpolation and much faster !
  m_Shrink = m_Filters && m_UserSetting_HalfSize;

  TRACEKEYVALS("Shrink","%d",m_Shrink);

  TRACEKEYVALS("BlackLevel","%d",m_BlackLevel);
  TRACEKEYVALS("WhiteLevel","%d",m_WhiteLevel);
  TRACEKEYVALS("Colors","%d",m_Colors);

  TRACEKEYVALS("Phase2 begin Width","%d",m_Width);
  TRACEKEYVALS("Phase2 begin Height","%d",m_Height);
  TRACEKEYVALS("Phase2 begin OutWidth","%d",m_OutWidth);
  TRACEKEYVALS("Phase2 begin OutHeight","%d",m_OutHeight);

  // Crop for detail view
  if (m_UserSetting_DetailView == 1 &&
      m_IsFuji == 0 &&
      m_PixelAspect == 1.0f) {
    ptCrop();
  }

  // Copied from earlier to here also.
  // Enables Phase3 to reenter with a different FourColorRGB setting.
  if ((m_UserSetting_Quality == ptInterpolation_VNG4 || m_UserSetting_HalfSize)
      && m_Filters != 0) { // m_Filters==0 <-> 3 channel RAWs
    m_MixGreen = 1;
  } else {
    m_MixGreen = 0;
  }

  if (!m_IsFoveon && m_UserSetting_AdjustMaximum > 0) {
    ptAdjustMaximum(1-m_UserSetting_AdjustMaximum);
  }

  assert (!m_IsFoveon);
  // if (m_IsFoveon ) foveon_interpolate();

  ptAdjustBlacks();

  if (m_UserSetting_GreenEquil) {
    TRACEKEYVALS("Green equilibration","%s","");
    green_equilibrate((float)m_UserSetting_GreenEquil/100.0f);
  }

  if (m_UserSetting_HotpixelReduction) {
    TRACEKEYVALS("Hotpixel reduction on bayer","%s","");
    ptHotpixelReductionBayer();
  }

  if (m_UserSetting_CfaLineDn !=0) {
    TRACEKEYVALS("Cfa line denoise","%s","");
    cfa_linedn(0.00002*(float)m_UserSetting_CfaLineDn);
  }

  if (!m_IsFoveon) {
    ptScaleColors();
  }

  if (m_UserSetting_CaCorrect !=0) {
    TRACEKEYVALS("CA correction","%s","");
    CA_correct(m_UserSetting_CaRed,m_UserSetting_CaBlue);
  }

  TRACEKEYVALS("Colors","%d",m_Colors);

  // not 1:1 pipe, use FC marco instead of interpolation
  uint16_t    (*TempImage)[4];
  if (m_Shrink && m_Filters != 2) { // -> preinterpolate will set m_Filters = 0 if != 2
    m_OutHeight = (m_Height + 1) / 2;
    m_OutWidth  = (m_Width + 1) / 2;
    TempImage = (uint16_t (*)[4]) CALLOC (m_OutHeight*m_OutWidth, sizeof *TempImage);
    merror (TempImage, "main()");
#pragma omp parallel for schedule(static)
    for (uint16_t row=0; row < m_Height; row++) {
      for (uint16_t col=0; col < m_Width; col++) {
        TempImage[((row) >> m_Shrink)*m_OutWidth + ((col) >> m_Shrink)][FC(row,col)] =
          m_Image[(row)*m_Width + (col)][FC(row,col)];
      }
    }
    FREE(m_Image);
    m_Image = TempImage;
  }

  // RG1BG2 -> RGB (if not 4 colors needed later on)
  pre_interpolate();

  TRACEKEYVALS("Colors","%d",m_Colors);


  // Interpolation/demosaicing according to one of the algorithms.
  if (m_Filters) {
    if (m_UserSetting_BayerDenoise && !m_UserSetting_HalfSize) {
      if (m_UserSetting_BayerDenoise==1) fbdd(0);
      else if (m_UserSetting_BayerDenoise==2) fbdd(1);
    }
    if (m_Filters == 2) // for 3x3 pattern we fix to VNG
      vng_interpolate();
    else if (m_UserSetting_Quality == ptInterpolation_Linear)
      lin_interpolate();
    else if (m_UserSetting_Quality == ptInterpolation_VNG ||
             m_UserSetting_Quality == ptInterpolation_VNG4)
      vng_interpolate();
    else if (m_UserSetting_Quality == ptInterpolation_PPG)
      ppg_interpolate();
    else if (m_UserSetting_Quality == ptInterpolation_DCB)
      dcb(m_UserSetting_InterpolationPasses, 1);
    else if (m_UserSetting_Quality == ptInterpolation_DCBSoft)
      dcb_interpolate_soft_old(m_UserSetting_InterpolationPasses, 1);
    else if (m_UserSetting_Quality == ptInterpolation_DCBSharp)
      dcb_interpolate_sharp_old(m_UserSetting_InterpolationPasses, 1);
    else if (m_UserSetting_Quality == ptInterpolation_AHD)
      ahd_interpolate();
    else if (m_UserSetting_Quality == ptInterpolation_AHD_mod)
      ahd_interpolate_mod();
    else if (m_UserSetting_Quality == ptInterpolation_VCD)
      vcd_interpolate(12);
    else if (m_UserSetting_Quality == ptInterpolation_LMMSE)
      lmmse_interpolate(1);
    else if (m_UserSetting_Quality == ptInterpolation_AMaZE)
      amaze_demosaic();
  } else {
    if (m_UserSetting_HotpixelReduction && m_UserSetting_HalfSize) {
      TRACEKEYVALS("Hotpixel reduction on RGB","%s","");
      ptHotpixelReduction();
    }
    if (m_UserSetting_BayerDenoise && !m_UserSetting_HalfSize) {
      if (m_UserSetting_BayerDenoise==1) fbdd(0);
      else if (m_UserSetting_BayerDenoise==2) fbdd(1);
    }
  }
  // If 1:1 and no interpolation is chosen show the Bayer pattern.

  TRACEKEYVALS("Interpolation type","%d",m_UserSetting_Quality);

  // Additional photivo stuff. Other halvings on request.
  if (m_UserSetting_HalfSize > 1                      ||
      (m_UserSetting_HalfSize == 1 && m_Shrink  == 0) || // 3 channel RAWs
      (m_UserSetting_HalfSize == 1 && m_Filters == 2)) { // 3x3 pattern
    short Factor = m_UserSetting_HalfSize - 1;
    if (m_Shrink == 0 || m_Filters == 2) Factor += 1;

    uint16_t NewHeight = m_Height >> Factor;
    uint16_t NewWidth = m_Width >> Factor;

    short Step = 1 << Factor;
    int Average = 2 * Factor;

    uint16_t (*NewImage)[4] =
      (uint16_t (*)[4]) CALLOC(NewWidth*NewHeight,sizeof(*m_Image));
    ptMemoryError(NewImage,__FILE__,__LINE__);

#pragma omp parallel for schedule(static)
    for (uint16_t Row=0; Row < NewHeight*Step; Row+=Step) {
      for (uint16_t Col=0; Col < NewWidth*Step; Col+=Step) {
        uint32_t  PixelValue[4] = {0,0,0,0};
        for (uint8_t sRow=0; sRow < Step; sRow++) {
          for (uint8_t sCol=0; sCol < Step; sCol++) {
            int32_t index = (Row+sRow)*m_Width+Col+sCol;
            for (short c=0; c < 4; c++) {
              PixelValue[c] += m_Image[index][c];
            }
          }
        }
        for (short c=0; c < 4; c++) {
          NewImage[Row/Step*NewWidth+Col/Step][c]
            = PixelValue[c] >> Average;
        }
      }
    }

    FREE(m_Image);
    m_Height = m_OutHeight = NewHeight;
    m_Width = m_OutWidth = NewWidth;
    m_Image = NewImage;
  }

  // Green mixing
  if (m_MixGreen && m_Colors != 3) {
#pragma omp parallel for schedule(static) default(shared)
    for (uint32_t i=0; i < (uint32_t) m_Height*m_Width; i++)
      m_Image[i][1] = (m_Image[i][1] + m_Image[i][3]) >> 1;
    m_Colors = 3;
  }

  // Median filter.
  if (!m_IsFoveon && m_Colors == 3) {
    if (m_UserSetting_MedianPasses > 0)  median_filter_new();
    if (m_UserSetting_ESMedianPasses > 0 && !m_UserSetting_HalfSize) es_median_filter();
    if (m_UserSetting_EeciRefine == 1)  refinement();
  }

  // Additional cleaning with hotpixel reduction
  if (m_UserSetting_HotpixelReduction && m_UserSetting_Quality != ptInterpolation_Bayer) {
    ptHotpixelReduction();
  }

  // XXX JDLA Additional steps for Fuji
  // And they don't hurt for others as they are early stopped.
  fuji_rotate();
  stretch();

  // Crop for detail view
  if (m_UserSetting_DetailView == 1 &&
      (m_IsFuji != 0 ||
       m_PixelAspect != 1.0f)) {
    ptCrop();
  }

  // Cache the image after Phase2.
  FREE(m_Image_AfterPhase2);
  m_Image_AfterPhase2 =
    (uint16_t (*)[4]) CALLOC (m_OutHeight*m_OutWidth, sizeof *m_Image);
  merror (m_Image_AfterPhase2, "main()");
  memcpy(m_Image_AfterPhase2,m_Image,m_OutHeight*m_OutWidth*sizeof(*m_Image));
  // Some other stuff to cache.
  m_Filters_AfterPhase2 = m_Filters;

  TRACEKEYVALS("BlackLevel","%d",m_BlackLevel);
  TRACEKEYVALS("WhiteLevel","%d",m_WhiteLevel);
  TRACEKEYVALS("Colors","%d",m_Colors);
  TRACEKEYVALS("Filters","%x",m_Filters);

  TRACEKEYVALS("Phase2 end Width","%d",m_Width);
  TRACEKEYVALS("Phase2 end Height","%d",m_Height);
  TRACEKEYVALS("Phase2 end OutWidth","%d",m_OutWidth);
  TRACEKEYVALS("Phase2 end OutHeight","%d",m_OutHeight);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// RunDcRaw_Phase3
//
// Highlights
//
////////////////////////////////////////////////////////////////////////////////

short CLASS RunDcRaw_Phase3(const short NoCache) {

  // Make sure we are starting from the right image.
  if (NoCache) {
    FREE(m_Image_AfterPhase2);
  } else {
    FREE(m_Image);
    m_Image =
      (uint16_t (*)[4]) CALLOC (m_OutHeight*m_OutWidth, sizeof *m_Image);
    merror (m_Image, "main()");
    memcpy(m_Image,m_Image_AfterPhase2,m_OutHeight*m_OutWidth*sizeof(*m_Image));
    // Restore some other cached values.
    m_Filters    = m_Filters_AfterPhase2;
  }

  TRACEKEYVALS("CamMult[0]","%f",VALUE(m_CameraMultipliers[0]));
  TRACEKEYVALS("CamMult[1]","%f",VALUE(m_CameraMultipliers[1]));
  TRACEKEYVALS("CamMult[2]","%f",VALUE(m_CameraMultipliers[2]));
  TRACEKEYVALS("CamMult[3]","%f",VALUE(m_CameraMultipliers[3]));
#ifdef TRACE_ORIGIN
  TRACEKEYVALS("CamMult File","%s",m_CameraMultipliers[0].File);
  TRACEKEYVALS("CamMult Line","%d",m_CameraMultipliers[0].Line);
#endif

  ptHighlight(m_UserSetting_photivo_ClipMode,
              m_UserSetting_photivo_ClipParameter);

  // Cache the image after Phase3.
  FREE(m_Image_AfterPhase3);
  m_Image_AfterPhase3 =
    (uint16_t (*)[4]) CALLOC (m_OutHeight*m_OutWidth, sizeof *m_Image);
  merror (m_Image_AfterPhase3, "main()");
  memcpy(m_Image_AfterPhase3,m_Image,m_OutHeight*m_OutWidth*sizeof(*m_Image));
  // Some other stuff to cache.
  m_Filters_AfterPhase3 = m_Filters;

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// BlendHighlights
// (original from dcraw TODO Refine and analyse)
//
////////////////////////////////////////////////////////////////////////////////

void CLASS ptBlendHighlights() {

  int ClipLevel=INT_MAX;
  int i;
  int j;

  static const float trans[2][4][4] =
  { { { 1,1,1 }, { 1.7320508,-1.7320508,0 }, { -1,-1,2 } },
    { { 1,1,1,1 }, { 1,-1,1,-1 }, { 1,1,-1,-1 }, { 1,-1,-1,1 } } };
  static const float itrans[2][4][4] =
  { { { 1,0.8660254,-0.5 }, { 1,-0.8660254,-0.5 }, { 1,0,1 } },
    { { 1,1,1,1 }, { 1,-1,1,-1 }, { 1,1,-1,-1 }, { 1,-1,-1,1 } } };
  float Cam[2][4], lab[2][4], Sum[2], chratio;

  const int transIdx = m_Colors - 3;
  assert(transIdx > -1);
  if (transIdx > 1) return;

  // to shut up gcc warnings...
  const int localColors = ptMin(static_cast<int>(m_Colors), 4);

  for (short c=0; c<m_Colors; c++) {
    if (ClipLevel > (i = (int)(0xFFFF*VALUE(m_PreMultipliers[c])))) {
      ClipLevel = i;
    }
  }

  for (uint16_t Row=0; Row < m_Height; Row++) {
    for (uint16_t Column=0; Column < m_Width; Column++) {
      short c;
      for (c=0; c<m_Colors; c++) {
        if (m_Image[Row*m_Width+Column][c] > ClipLevel) break;
      }
      if (c == m_Colors) continue; // No clip
      for (c=0; c<m_Colors; c++) {
  Cam[0][c] = m_Image[Row*m_Width+Column][c];
  Cam[1][c] = MIN(Cam[0][c],(float)ClipLevel);
      }
      for (i=0; i < 2; i++) {
  for (c=0; c<m_Colors; c++) {
          for (lab[i][c]=j=0; j < m_Colors; j++) {
      lab[i][c] += trans[transIdx][c][j] * Cam[i][j];
          }
        }
  for (Sum[i]=0,c=1; c < m_Colors; c++) {
    Sum[i] += SQR(lab[i][c]);
        }
      }
      chratio = sqrt(Sum[1]/Sum[0]);
      for (c = 1; c < m_Colors; c++) {
        lab[0][c] *= chratio;
      }
      for (c = 0; (c < localColors); c++) {
        for (Cam[0][c]=j=0; (j < localColors); j++) {
          Cam[0][c] += itrans[transIdx][c][j] * lab[0][j];
        }
      }
      for (c = 0; c < m_Colors; c++) {
        m_Image[Row*m_Width+Column][c] = (uint16_t)(Cam[0][c] / m_Colors);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Hotpixel reduction
//
////////////////////////////////////////////////////////////////////////////////

void CLASS ptHotpixelReduction() {
  uint16_t HotpixelThreshold = 0x00ff;
  uint16_t Threshold = (int32_t) ((1.0-m_UserSetting_HotpixelReduction)*0x2fff);
  uint16_t Width = m_OutWidth;
  uint16_t Height = m_OutHeight;
  // bei 1:2 oder kleiner m_Image hat in jedem Punkt RGBG
  // bei 1:1 m_Image hat Bayer Pattern
  // Länge und Breite OutWidth und OutHeight
  // unterschieptiche pixel in unterschieptichen ebenen ansprechen und werte zur luminanz skalieren
#pragma omp parallel for schedule(static) default(shared)
  for (uint16_t Row=0; Row<Height; Row++) {
    for (uint16_t Col=0; Col<Width; Col++) {
      uint16_t TempValue = 0;
      for (int c=0; c<4; c++) {
        if (m_Image[Row*Width+Col][c] > HotpixelThreshold) {
          if (Row > 1) TempValue = MAX(m_Image[(Row-1)*Width+Col][c],TempValue);
          if (Row < Height-1) TempValue = MAX(m_Image[(Row+1)*Width+Col][c],TempValue);
          if (Col > 1) TempValue = MAX(m_Image[Row*Width+Col-1][c],TempValue);
          if (Col < Width-1) TempValue = MAX(m_Image[Row*Width+Col+1][c],TempValue);
          if (TempValue+Threshold<m_Image[Row*Width+Col][c])
            m_Image[Row*Width+Col][c] = TempValue;
        }
        if (m_Image[Row*Width+Col][c] < 10*HotpixelThreshold) {
          TempValue = 0xffff;
          if (Row > 1) TempValue = MIN(m_Image[(Row-1)*Width+Col][c],TempValue);
          if (Row < Height-1) TempValue = MIN(m_Image[(Row+1)*Width+Col][c],TempValue);
          if (Col > 1) TempValue = MIN(m_Image[Row*Width+Col-1][c],TempValue);
          if (Col < Width-1) TempValue = MIN(m_Image[Row*Width+Col+1][c],TempValue);
          if (TempValue-Threshold>m_Image[Row*Width+Col][c])
            m_Image[Row*Width+Col][c] = TempValue;
        }
      }
    }
  }
}

void CLASS ptHotpixelReductionBayer() {
  uint16_t HotpixelThreshold = 0x00ff;
  uint16_t Threshold = (int32_t) ((1.0-m_UserSetting_HotpixelReduction)*0x2fff);
  uint16_t Width = m_OutWidth;
  uint16_t Height = m_OutHeight;

#pragma omp parallel for schedule(static)
  for (uint16_t Row=0; Row<Height; Row++) {
    for (uint16_t Col=0; Col<Width; Col++) {
      uint16_t TempValue = 0;
      for (int Channel=0; Channel<4; Channel++) {
        if (m_Image[Row*Width+Col][Channel] > HotpixelThreshold) {
          if (Row > 2) TempValue = MAX(m_Image[(Row-2)*Width+Col][Channel],TempValue);
          if (Row < Height-2) TempValue = MAX(m_Image[(Row+2)*Width+Col][Channel],TempValue);
          if (Col > 2) TempValue = MAX(m_Image[Row*Width+Col-2][Channel],TempValue);
          if (Col < Width-2) TempValue = MAX(m_Image[Row*Width+Col+2][Channel],TempValue);
          if (TempValue+Threshold<m_Image[Row*Width+Col][Channel]) {
            for (int c=0; c<4; c++) {
              if (Row > 1) {
                TempValue = MAX(m_Image[(Row-1)*Width+Col][c],TempValue);
                if (Col > 1) TempValue = MAX(m_Image[(Row-1)*Width+Col-1][c],TempValue);
                if (Col < Width-1) TempValue = MAX(m_Image[(Row-1)*Width+Col+1][c],TempValue);
              }
              if (Row < Height-1) {
                TempValue = MAX(m_Image[(Row+1)*Width+Col][c],TempValue);
                if (Col > 1) TempValue = MAX(m_Image[(Row+1)*Width+Col-1][c],TempValue);
                if (Col < Width-1) TempValue = MAX(m_Image[(Row+1)*Width+Col+1][c],TempValue);
              }
              if (Col > 1) TempValue = MAX(m_Image[Row*Width+Col-1][c],TempValue);
              if (Col < Width-1) TempValue = MAX(m_Image[Row*Width+Col+1][c],TempValue);
            }
            if (TempValue+Threshold<m_Image[Row*Width+Col][Channel])
              m_Image[Row*Width+Col][Channel] = TempValue;
          }
        }
      }
      for (int Channel=0; Channel<4; Channel++) {
        if (m_Image[Row*Width+Col][Channel] < 10*HotpixelThreshold) {
          TempValue = 0xffff;
          if (Row > 2) TempValue = MIN(m_Image[(Row-2)*Width+Col][Channel],TempValue);
          if (Row < Height-2) TempValue = MIN(m_Image[(Row+2)*Width+Col][Channel],TempValue);
          if (Col > 2) TempValue = MIN(m_Image[Row*Width+Col-2][Channel],TempValue);
          if (Col < Width-2) TempValue = MIN(m_Image[Row*Width+Col+2][Channel],TempValue);
          if (TempValue-Threshold>m_Image[Row*Width+Col][Channel]) {
            for (int c=0; c<4; c++) {
              if (Row > 1) {
                TempValue = MIN(m_Image[(Row-1)*Width+Col][c],TempValue);
                if (Col > 1) TempValue = MIN(m_Image[(Row-1)*Width+Col-1][c],TempValue);
                if (Col < Width-1) TempValue = MIN(m_Image[(Row-1)*Width+Col+1][c],TempValue);
              }
              if (Row < Height-1) {
                TempValue = MIN(m_Image[(Row+1)*Width+Col][c],TempValue);
                if (Col > 1) TempValue = MIN(m_Image[(Row+1)*Width+Col-1][c],TempValue);
                if (Col < Width-1) TempValue = MIN(m_Image[(Row+1)*Width+Col+1][c],TempValue);
              }
              if (Col > 1) TempValue = MIN(m_Image[Row*Width+Col-1][c],TempValue);
              if (Col < Width-1) TempValue = MIN(m_Image[Row*Width+Col+1][c],TempValue);
            }
            if (TempValue-Threshold>m_Image[Row*Width+Col][Channel])
              m_Image[Row*Width+Col][Channel] = TempValue;
          }
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Crop for detail view
//
////////////////////////////////////////////////////////////////////////////////

inline uint16_t MultOf6(uint16_t AValue) {
  return ptMax(AValue - (AValue % 6), 0);
}

void CLASS ptCrop() {
  assert (m_UserSetting_HalfSize == 0);

  uint16_t (*TempImage)[4];
  uint16_t CropW = m_UserSetting_DetailViewCropW;
  uint16_t CropH = m_UserSetting_DetailViewCropH;
  uint16_t CropX = m_UserSetting_DetailViewCropX;
  uint16_t CropY = m_UserSetting_DetailViewCropY;

  if (m_Filters == 2) { // for 3x3 pattern
    CropW = MultOf6(CropW);
    CropH = MultOf6(CropH);
    CropX = MultOf6(CropX);
    CropY = MultOf6(CropY);
  }

  if (m_Flip & 2) {
    CropX = m_Height - CropX - CropW;
  }
  if (m_Flip & 1) {
    CropY = m_Width - CropY - CropH;
  }
  if (m_Flip & 4) {
    SWAP(CropW, CropH);
    SWAP(CropX, CropY);
  }
  m_OutHeight = CropH;
  m_OutWidth  = CropW;
  TempImage = (uint16_t (*)[4]) CALLOC (m_OutHeight*m_OutWidth, sizeof *TempImage);
  merror (TempImage, "Temp for detail view");

#pragma omp parallel for schedule(static)
  for (uint16_t row=0; row < m_OutHeight; row++) {
    for (uint16_t col=0; col < m_OutWidth; col++) {
      for (short c=0; c<4; c++) {
        TempImage[row *m_OutWidth + col][c] =
          m_Image[(row+CropY)*m_Width + (col+CropX)][c];
      }
    }
  }

  m_Width = m_OutWidth;
  m_Height = m_OutHeight;
  FREE(m_Image);
  m_Image = TempImage;
  TRACEKEYVALS("Phase2 detail view Width","%d",m_Width);
  TRACEKEYVALS("Phase2 detail view Height","%d",m_Height);
  TRACEKEYVALS("Phase2 detail view OutWidth","%d",m_OutWidth);
  TRACEKEYVALS("Phase2 detail view OutHeight","%d",m_OutHeight);
}

////////////////////////////////////////////////////////////////////////////////
//
// MedianFilter
// Repeatepty 3x3 median filter on R-G and B-G channels
// Repeat : m_UserSetting_MedianPasses
//
////////////////////////////////////////////////////////////////////////////////

void CLASS ptMedianFilter() {

  int      Median[9];
  static const uint8_t opt[] =  /* Optimal 9-element median search */
  { 1,2, 4,5, 7,8, 0,1, 3,4, 6,7, 1,2, 4,5, 7,8,
    0,3, 5,8, 4,7, 3,6, 1,4, 2,5, 4,7, 4,2, 6,4, 4,2 };

  for (short Pass=1; Pass <= m_UserSetting_MedianPasses; Pass++) {
    for (short c=0; c < 3; c+=2) {
#pragma omp parallel for schedule(static) default(shared)
      for (int32_t i = 0; i< (int32_t)m_Width*m_Height; i++) {
        m_Image[i][3] = m_Image[i][c];
      }
#pragma omp parallel for schedule(static) default(shared) private(Median)
      for (int32_t n=m_Width; n<m_Width*(m_Height-1); n++) {
        if ((n+1) % m_Width < 2) continue;
        short k=0;
        for (int32_t i = -m_Width; i <= m_Width; i += m_Width) {
          for (int32_t j = i-1; j <= i+1; j++) {
            Median[k++] = m_Image[n+j][3] - m_Image[n+j][1];
          }
        }
        for (unsigned short i=0; i < sizeof opt; i+=2) {
          if (Median[opt[i]] > Median[opt[i+1]]) {
            SWAP (Median[opt[i]] , Median[opt[i+1]]);
          }
        }
        m_Image[n][c] = CLIP(Median[4] + m_Image[n][1]);
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// RebuildHighlights
// (original from dcraw TODO Refine and analyse)
//
////////////////////////////////////////////////////////////////////////////////

#define SCALE (4 >> m_Shrink)
void CLASS ptRebuildHighlights(const short highlight) {
  float *map, sum, wgt, grow;
  int hsat[4], count, spread, change, val, i;
  unsigned high, wide, mrow, mcol, row, col, d, y, x;
  uint16_t *pixel;
  static const signed char dir[8][2] =
    { {-1,-1}, {-1,0}, {-1,1}, {0,1}, {1,1}, {1,0}, {1,-1}, {0,-1} };

  grow = pow (2, 1-highlight);
  for (short c=0;c<m_Colors;c++)  hsat[c] = (int)(32000 * VALUE(m_PreMultipliers[c]));
  short kc = 0;
  for (short c=1; c < m_Colors; c++)
    if (VALUE(m_PreMultipliers[kc]) < VALUE(m_PreMultipliers[c])) kc = c;
  high = m_Height / SCALE;
  wide =  m_Width / SCALE;
  map = (float *) CALLOC (high*wide, sizeof *map);
  merror (map, "recover_highlights()");
  for (short c=0;c<m_Colors;c++) {
    if (c != kc) {
      memset (map, 0, high*wide*sizeof *map);
#pragma omp parallel for schedule(static) private(mrow, mcol, sum, wgt, count, pixel, row, col)
      for (mrow=0; mrow < high; mrow++) {
        for (mcol=0; mcol < wide; mcol++) {
          sum = wgt = count = 0;
          for (row = mrow*SCALE; row < (mrow+1)*SCALE; row++) {
            for (col = mcol*SCALE; col < (mcol+1)*SCALE; col++) {
              pixel = m_Image[row*m_Width+col];
              if (pixel[c] / hsat[c] == 1 && pixel[kc] > 24000) {
                sum += pixel[c];
                wgt += pixel[kc];
                count++;
              }
            }
            if (count == SCALE*SCALE) map[mrow*wide+mcol] = sum / wgt;
          }
        }
      }
      for (spread = (int)(32/grow); spread--; ) {
#pragma omp parallel for schedule(static) private(mrow, mcol, sum, count, x, y, d)
        for (mrow=0; mrow < high; mrow++) {
          for (mcol=0; mcol < wide; mcol++) {
            if (map[mrow*wide+mcol]) continue;
            sum = count = 0;
            for (d=0; d < 8; d++) {
              y = mrow + dir[d][0];
              x = mcol + dir[d][1];
              if (y < high && x < wide && map[y*wide+x] > 0) {
                sum  += (1 + (d & 1)) * map[y*wide+x];
                count += 1 + (d & 1);
              }
            }
            if (count > 3)
              map[mrow*wide+mcol] = - (sum+grow) / (count+grow);
          }
        }
        change = 0;
#pragma omp parallel for schedule(static) private(i) reduction(||:change)
        for (i=0; i < (int)(high*wide); i++) {
          if (map[i] < 0) {
            map[i] = -map[i];
            change = 1;
          }
        }
        if (!change) break;
      }
#pragma omp parallel for schedule(static) private(i)
      for (i=0; i < (int)(high*wide); i++) {
        if (map[i] == 0) map[i] = 1;
      }
#pragma omp parallel for schedule(static) private(mrow, mcol, row, col, pixel, val)
      for (mrow=0; mrow < high; mrow++) {
        for (mcol=0; mcol < wide; mcol++) {
          for (row = mrow*SCALE; row < (mrow+1)*SCALE; row++) {
            for (col = mcol*SCALE; col < (mcol+1)*SCALE; col++) {
              pixel = m_Image[row*m_Width+col];
              if (pixel[c] / hsat[c] > 1) {
                val = (int)(pixel[kc] * map[mrow*wide+mcol]);
                if (pixel[c] < val) pixel[c] = CLIP(val);
              }
            }
          }
        }
      }
    }
  }
  FREE (map);
}
#undef SCALE

////////////////////////////////////////////////////////////////////////////////
//
// LabToCam helper function.
//
// FIXME
// At this moment this is verbatim copies of ptImage.cpp. Room for improvement.
//
////////////////////////////////////////////////////////////////////////////////

void CLASS LabToCam(double Lab[3],uint16_t Cam[4]) {

  if(!ToCamFunctionInited) {
    for(short i=0;i<m_Colors;i++) {
      for (short j=0;j<3;j++) {
        MatrixXYZToCam[i][j] = 0.0;
        for (short k=0;k<3;k++) {
          MatrixXYZToCam[i][j] +=
            m_MatrixSRGBToCamRGB[i][k] * MatrixXYZToRGB[ptSpace_sRGB_D65][k][j];
        }
      }
    }
  ToCamFunctionInited = 1;
  }

  double xyz[3];
  double fx,fy,fz;
  double xr,yr,zr;
  const double epsilon = 216.0/24389.0;
  const double kappa   = 24389.0/27.0;

  const double L = Lab[0];
  const double a = Lab[1];
  const double b = Lab[2];

  const double Tmp1 = (L+16.0)/116.0;

  yr = (L<=kappa*epsilon)?
       (L/kappa):(Tmp1*Tmp1*Tmp1);
  fy = (yr<=epsilon) ? ((kappa*yr+16.0)/116.0) : Tmp1;
  fz = fy - b/200.0;
  fx = a/500.0 + fy;
  const double fz3 = fz*fz*fz;
  const double fx3 = fx*fx*fx;
  zr = (fz3<=epsilon) ? ((116.0*fz-16.0)/kappa) : fz3;
  xr = (fx3<=epsilon) ? ((116.0*fx-16.0)/kappa) : fx3;

  xyz[0] = xr*D65Reference[0]*65535.0 - 0.5;
  xyz[1] = yr*D65Reference[1]*65535.0 - 0.5;
  xyz[2] = zr*D65Reference[2]*65535.0 - 0.5;

  // And finally to RGB via the matrix.
  for (short c=0; c<m_Colors; c++) {
    double Value = 0;
    Value += MatrixXYZToCam[c][0] * xyz[0];
    Value += MatrixXYZToCam[c][1] * xyz[1];
    Value += MatrixXYZToCam[c][2] * xyz[2];
    Cam[c] = (uint16_t) CLIP((int32_t)(Value));
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// CamToLab helper function.
//
////////////////////////////////////////////////////////////////////////////////

void CLASS CamToLab(uint16_t Cam[4], double Lab[3]) {

  // Initialize the lookup table for the RGB->LAB function
  // if this would be the first time.
  if (!ToLABFunctionInited) {
    // Remark : we extend the table well beyond r>1.0 for numerical
    // stability purposes. XYZ>1.0 occurs sometimes and this way
    // it stays stable (srgb->lab->srgb correct within 0.02%)
    for (uint32_t i=0; i<0x20000; i++) {
      double r = (double)(i) / 0xffff;
      ToLABFunctionTable[i] = r > 216.0/24389.0 ?
                              pow(r,1/3.0) :
                              (24389.0/27.0*r + 16.0)/116.0;
    }
    for(short i=0;i<3;i++) {
      for (short j=0;j<m_Colors;j++) {
        MatrixCamToXYZ[i][j] = 0.0;
        for (short k=0;k<3;k++) {
          MatrixCamToXYZ[i][j] +=
            MatrixRGBToXYZ[ptSpace_sRGB_D65][i][k] * m_MatrixCamRGBToSRGB[k][j];
        }
      }
    }
  ToLABFunctionInited = 1;
  }

  // First go to XYZ
  double xyz[3] = {0.5,0.5,0.5};
  for (short Color = 0; Color < m_Colors; Color++) {
    xyz[0] += MatrixCamToXYZ[0][Color] * Cam[Color];
    xyz[1] += MatrixCamToXYZ[1][Color] * Cam[Color];
    xyz[2] += MatrixCamToXYZ[2][Color] * Cam[Color];
  }

  // Reference White
  xyz[0] /= D65Reference[0];
  xyz[1] /= D65Reference[1];
  xyz[2] /= D65Reference[2];

  // Then to Lab
  xyz[0] = ToLABFunctionTable[ (int32_t) MAX(0.0,xyz[0]) ];
  xyz[1] = ToLABFunctionTable[ (int32_t) MAX(0.0,xyz[1]) ];
  xyz[2] = ToLABFunctionTable[ (int32_t) MAX(0.0,xyz[2]) ];

  // L in 0 , a in 1, b in 2
  Lab[0] = 116.0 * xyz[1] - 16.0;
  Lab[1] = 500.0*(xyz[0]-xyz[1]);
  Lab[2] = 200.0*(xyz[1]-xyz[2]);
}

//==============================================================================

TImage8RawData ptDcRaw::thumbnail() {
  m_Thumb.clear();

  if (m_rawProcessor->unpack_thumb() == LIBRAW_SUCCESS) {
    const auto &thumb = m_rawProcessor->imgdata.thumbnail;
    m_Thumb.insert(m_Thumb.end(), &thumb.thumb[0],
                   &thumb.thumb[thumb.tlength]);
    m_ThumbWidth = thumb.twidth;
    m_ThumbHeight = thumb.theight;
  }

  return m_Thumb;
}

