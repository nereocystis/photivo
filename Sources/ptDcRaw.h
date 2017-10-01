/*******************************************************************************
**
** Photivo
**
** Copyright (C) 2008,2009 Jos De Laender <jos.de_laender@telenet.be>
** Copyright (C) 2009,2010 Michael Munzert <mail@mm-log.com>
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
** This is the adaptation of dcraw.c into something that is more or less
** a C++ object.
** It is not implied that this now works correctly with different compilers,
** environments and settings. Lots of testing still needed.
**
*******************************************************************************/

#ifndef DLDCDRAW_H
#define DLDCDRAW_H

// Adaptation of dcraw.c stuff.
#define DCRAW_VERSION "dcraw v9.17"    // Update along with dcraw syncing ...
// $Revision: 1.454 $
// $Date: 2012/12/23 19:25:36 $
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

#include "ptDefines.h"

#include <libraw/libraw_types.h>

#include <QString>

#define _USE_MATH_DEFINES
#include <vector>
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cassert>
#include <csetjmp>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fcntl.h>
#include <sys/types.h>

#ifdef WIN32
#define NO_JPEG
#endif

#ifndef NO_JPEG
#ifdef __cplusplus
  // This hack copes with jpeglib.h that does or doesnt provide the
  // extern internally.
  #pragma push_macro("__cplusplus")
  #undef __cplusplus
  extern "C" {
  #include <jpeglib.h>
  }
  #pragma pop_macro("__cplusplus")
#else
  #include <jpeglib.h>
#endif
#endif

#ifdef LOCALEDIR
#include <libintl.h>
#define _(String) gettext(String)
#else
#define _(String) (String)
#endif
/* jpta away  ?
#ifdef DJGPP
#define fseeko fseek
#define ftello ftell
#else
#define fgetc getc_unlocked
#endif
*/
#ifdef __CYGWIN__
#include <io.h>
#endif
#ifdef WIN32
#include <sys/utime.h>
#include <winsock2.h>
/* jpta : away ?
comment(lib, "ws2_32.lib")
*/
#define snprintf _snprintf
#define strcasecmp stricmp
#define strncasecmp strnicmp
typedef __int64 INT64;
typedef unsigned __int64 UINT64;
#else
#include <unistd.h>
#include <utime.h>
#include <netinet/in.h>
typedef long long INT64;
typedef unsigned long long UINT64;
#endif

#ifdef LJPEG_DECODE
#error Please compile dcraw.c by itself.
#error Do not link it with ljpeg_decode.
#endif

#ifndef LONG_BIT
#define LONG_BIT (8 * sizeof (long))
#endif


////////////////////////////////////////////////////////////////////////////////
//
// TRACE_ORIGIN
//
// Macro to try and analyze some of the origins of variables in the
// dcraw jungle
//
////////////////////////////////////////////////////////////////////////////////

// #define TRACE_ORIGIN

#ifdef TRACE_ORIGIN
#define ASSIGN(x,y) {x.Value=y;strncpy(x.File, __FILE__,1023);x.Line=__LINE__;}
#define VALUE(x) x.Value
#else
#define ASSIGN(x,y) x=y
#define VALUE(x)    x
#endif

class ptRawWrapper;


class ptDcRaw {
  ptRawWrapper *m_rawProcessor;

////////////////////////////////////////////////////////////////////////////////
//
// PUBLIC members
//
////////////////////////////////////////////////////////////////////////////////
public:
  // internal data of libraw
  libraw_data_t &imgdata;

  /*************************************************************************
  m_UserSetting_* :

  Members to be set (or changed from their default) by the user
  in order to define what dcraw will do.

  User settings corresponding to command line optins.
  Those are going to be set when used as embedded in an app.
  Each of them corresponds roughly to a command line parameter of dcraw.
  *************************************************************************/

  // Quality.
  // Please remark that most probably dcraw implementation and
  // manual/documentation differ here.
  //   According to manual :
  //     0 : Use high-speed, low-quality bilinear interpolation.
  //     2 : Use Variable Number of Gradients (VNG) interpolation.
  //     3 : Use Adaptive Homogeneity-Directed (AHD) interpolation.
  //   According to my reading of the code :
  //     0 : Use high-speed, low-quality bilinear interpolation.
  //     1 : Use Variable Number of Gradients (VNG) interpolation.
  //     2 : Use Patterned Pixel Grouping (PPG) interpolation.
  //     3 : Use Adaptive Homogeneity-Directed (AHD) interpolation.
  short       m_UserSetting_Quality;

  // Output half-size color image. Doubles speed (no interpolation needed)
  // 1 : Just half size as in original dcraw.
  // 2 : An additional halving of the image will be done.
  //     (4 pixels to 1 by averaging = bilinear)
  short       m_UserSetting_HalfSize;

  // Automatic hotpixel reduction
  double      m_UserSetting_HotpixelReduction;

  // Bayer denoise
  short       m_UserSetting_BayerDenoise;

  // Green equilibration
  int         m_UserSetting_GreenEquil;

  // CA auto correction
  int         m_UserSetting_CaCorrect;
  float       m_UserSetting_CaRed;
  float       m_UserSetting_CaBlue;

  // CFA line denoise
  int         m_UserSetting_CfaLineDn;

  // Automatic white balance.  The default is to use  a  fixed  color
  // balance based on a white card photographed in sunlight
  short       m_UserSetting_AutoWb;

  // Blackpoint setting.
  // -1 is the default depending on the camera.
  int         m_UserSetting_BlackPoint;

  // Whitepoint setting.
  // -1 is the default depending on the camera.
  int         m_UserSetting_Saturation;

  // Use the white balance specified by the camera.  If this can’t be
  // found, print a warning and revert to the default.
  short       m_UserSetting_CameraWb;

  // <x y width height> Average a grey box for white balance
  // This defines the area that will be used for the white balance.
  bool        m_UserSetting_UseGreyBox;
  unsigned    m_UserSetting_GreyBox[4]; // type like libraw

  // Specify your own raw white balance.  These  multipliers  can  be
  // cut and pasted from the output of dcraw -v.
  float       m_UserSetting_Multiplier[4];

  // Use or not the embedded color matrix of the camera.
  short       m_UserSetting_CameraMatrix;

  // Filename to process.
  QString     m_UserSetting_InputFileName;

  // Detail view
  short       m_UserSetting_DetailView;
  uint16_t    m_UserSetting_DetailViewCropX;
  uint16_t    m_UserSetting_DetailViewCropY;
  uint16_t    m_UserSetting_DetailViewCropW;
  uint16_t    m_UserSetting_DetailViewCropH;

  // Adjust maximum threshold (between 0 and 0.25)
  double      m_UserSetting_AdjustMaximum;

  // Wavelet denoising before color balancing and demosaicing.
  int         m_UserSetting_DenoiseThreshold;

  // Interpolation passes for demosaicing.
  int         m_UserSetting_InterpolationPasses;

  // Median passes after demosaicing.
  int         m_UserSetting_MedianPasses;
  int         m_UserSetting_ESMedianPasses;
  int         m_UserSetting_EeciRefine;

  // photivo in the name is a photivo specific adapation/extension.

  // Clip/Highlight related parameters
  short       m_UserSetting_photivo_ClipMode;
  short       m_UserSetting_photivo_ClipParameter;

  // [0..N-1] Select one raw image or "all" from each file
  unsigned  m_UserSetting_ShotSelect;

  // Bad pixel filename.
  // Default on NULL and not used.
  char*     m_UserSetting_BadPixelsFileName;

  // Dark frame filename.
  // Default on NULL and not used.
  char*     m_UserSetting_DarkFrameFileName;


  /**************************************************************************
    Methods that the user can apply.
  **************************************************************************/

  // Constructor.
  ptDcRaw();

  // Destructor.
  ~ptDcRaw();

  // ResetNonUserSettings. (For second entry from Identify on).
  void ResetNonUserSettings();

  // Identify the input file, the camera parameters etc.
  // See further for stuff that can be read afterwards.
  short  Identify(const QString &NewInputFile = "");

  // Do the raw processing up to the image available.
  short  RunDcRaw_Phase1(); //Load,bad pxs,darkframe.
  short  RunDcRaw_Phase2(const short NoCache=0); // WB,Interpolation
  short  RunDcRaw_Phase3(const short NoCache=0); // Highlights

  // Member functions specifically for photivo or photivo adapted.

  void  ptAdjustBlacks();
  void  ptScaleColors();
  void  ptMedianFilter();
  void  ptHotpixelReduction();
  void  ptHotpixelReductionBayer();
  void  ptHighlight(const short  ClipMode,
                    const short  ClipParameter);
  void  ptAdjustMaximum(double Threshold);

  void  LabToCam(double Lab[3],uint16_t Cam[4]);
  void  CamToLab(uint16_t Cam[4], double Lab[3]);

  void  ptRebuildHighlights(const short Effort);
  void  ptBlendHighlights();
  void  ptCrop();
  TImage8RawData thumbnail();


  /*************************************************************************
  m_*

  Members that are set during the program logic and not by the user.
  Most of them are useful to be *read* by the user. Don't change them.
  The lower ones are questionable if they should be of interest to the user.
  *************************************************************************/

  TImage8RawData m_Thumb;

  // The image !
  uint16_t    (*m_Image)[4];
  uint16_t    (*m_Image_AfterPhase1)[4];  // Cached one.
  uint16_t    (*m_Image_AfterPhase2)[4];  // Cached one.
  uint16_t    (*m_Image_AfterPhase3)[4];  // Cached one.
  uint16_t    (*m_Image_AfterPhase4)[4];  // Cached one.

  // Height and width of output.
  uint16_t    m_Width;
  uint16_t    m_Width_AfterPhase1; // Cached one.
  uint16_t    m_Height;
  uint16_t    m_Height_AfterPhase1; // Cached one.
  uint16_t    m_ReportedWidth; // Take into account already fuji rotate/stretch
  uint16_t    m_ReportedHeight;

  // Camera and photo identifications. Typically available after Identify().
  QString   m_CameraMake;
  QString   m_CameraModel;
  char      m_Artist[64];
  float     m_IsoSpeed;
  float     m_Shutter;
  float     m_Aperture;
  float     m_FocalLength;
  time_t    m_TimeStamp;

  // Nr of raw images found in the file.
  unsigned  m_IsRaw;

  // Not 0 if an embedded ICC profile was found.
  unsigned  m_ProfileLength;

  // Pixel aspect ratio. Whatever this may be.
  float    m_PixelAspect;

  // Dimensions of thumb.
  uint16_t  m_ThumbWidth;
  uint16_t  m_ThumbHeight;

  // Dimensions of the full size image.
  uint16_t  m_RawHeight;
  uint16_t  m_RawWidth;

  // Raw colors.
  short     m_Colors;

  // Daylight multipliers.

  #ifdef TRACE_ORIGIN
  struct {
    float Value;
    char  File[1024];
    int   Line;
  } m_PreMultipliers[4];
  #else
  float     m_PreMultipliers[4];
  #endif

  float     m_MinPreMulti; // Max per construction 1.0 after ptScaleColors.
  //TODO Mike: Max is no longer 1.0, something todo here?

  // Camera multipliers.

  #ifdef TRACE_ORIGIN
  struct {
    float Value;
    char  File[1024];
    int   Line;
  } m_CameraMultipliers[4];
  #else
  float     m_CameraMultipliers[4];
  #endif

  // Multipliers as obtained at the end of 'phase1'. Should be D65 ones (jpta ?)

  #ifdef TRACE_ORIGIN
  struct {
    float Value;
    char  File[1024];
    int   Line;
  } m_D65Multipliers[4];
  #else
  float m_D65Multipliers[4];
  #endif

  // Final multipliers based on pre, white, black

  #ifdef TRACE_ORIGIN
  struct {
    float Value;
    char  File[1024];
    int   Line;
  } m_Multipliers[4];
  #else
  float     m_Multipliers[4];
  #endif


  /*************************************************************************
  From here on it becomes more and more uninteresting
  internal data. If someone can describe it in a useful way, feel free.
  *************************************************************************/
  char      m_ColorDescriptor[5];
  short     m_ByteOrder;  // 0x4949 ("II") means little-endian.
  FILE*     m_OutputFile;
  char      m_Description[512];
  unsigned  m_ShotOrder;
  unsigned  m_Kodak_cbpp;
  unsigned  m_Filters;
  unsigned  m_Filters_AfterPhase1;
  unsigned  m_Filters_AfterPhase2;
  unsigned  m_Filters_AfterPhase3;
  off_t     m_Data_Offset;
  off_t     m_ProfileOffset;
  unsigned  m_Tiff_NrIFDs;
  unsigned  m_Tiff_Samples;
  unsigned  m_Tiff_bps;
  unsigned  m_Tiff_Compress;
  unsigned  m_BlackLevel;
  unsigned  m_CBlackLevel[4];
  unsigned  m_BlackLevel_AfterPhase1;
  unsigned  m_CBlackLevel_AfterPhase1[4];
  unsigned  m_WhiteLevel;
  unsigned  m_WhiteLevel_AfterPhase1;
  unsigned  m_MixGreen;
  unsigned  m_RawColor;
  unsigned  m_RawColorPhotivo; // mirror for Photivo
  unsigned  m_ZeroIsBad;
  unsigned  m_DNG_Version;
  unsigned  m_IsFoveon;
  unsigned  m_TileWidth;
  unsigned  m_TileLength;
  unsigned  m_Load_Flags;
  uint16_t  m_TopMargin;
  short     m_LeftMargin;
  uint16_t  m_Shrink;
  uint16_t  m_OutWidth;
  uint16_t  m_OutWidth_AfterPhase1;
  uint16_t  m_OutHeight;
  uint16_t  m_OutHeight_AfterPhase1;
  uint16_t  m_Fuji_Width;
  uint16_t  m_IsFuji;
  int       m_Flip;
  int       m_Tiff_Flip;
  int       m_Mask[8][4];
  uint16_t  m_Curve[0x10000];
  uint16_t* m_Raw_Image;
  float     m_cmatrix[3][4];
  float     m_MatrixCamRGBToSRGB[3][4];
  double    m_MatrixSRGBToCamRGB[4][3]; // addon photivo
  void      (ptDcRaw::*m_WriteFunction)();
  void      (ptDcRaw::*m_LoadRawFunction)();

  jmp_buf   m_Failure;

  double    m_Gamma[6];

  // TODO : not yet m_ified.

  float     flash_used;
  float     canon_ev;
  unsigned  fuji_layout;
  unsigned  exif_cfa;
  unsigned  unique_id;
  off_t     strip_offset;
  off_t     meta_offset;
  unsigned  zero_after_ff;
  unsigned  data_error;
  uint16_t  white[8][8];
  uint16_t  cr2_slice[3];
  uint16_t  sraw_mul[4];

  unsigned    m_getbithuff_bitbuf;
  int         m_getbithuff_reset;
  int         m_getbithuff_vbits;
  uint8_t     m_pana_bits_buf[0x4000];
  int         m_pana_bits_vbits;
  uint8_t     jpeg_buffer[4096];
  unsigned    m_sony_decrypt_pad[128];
  unsigned    m_sony_decrypt_p;
  unsigned    m_foveon_decoder_huff[1024];
  double      MatrixXYZToCam[4][3];
  short       ToCamFunctionInited;

  // A lookup table used for this purpose.
  // Initialized at the first call.
  double      ToLABFunctionTable[0x20000];
  short       ToLABFunctionInited;
  double      MatrixCamToXYZ[3][4];

  /*************************************************************************
  Internally used member functions. Should be completely uninteresting
  for the user.
  Pseudo automatically generated.
  *************************************************************************/

  void  stretch();
  void  fuji_rotate();
  void  identify();
  float find_green (int bps, int bite, int off0, int off1);
  void  parse_foveon();
  void  parse_redcine();
  char *foveon_gets(int offset,char *str,int len);
  void  parse_cine();
  void  parse_smal(int offset,int fsize);
  void  parse_riff();
  int   parse_jpeg(int offset);
  void  parse_fuji(int offset);
  void  parse_phase_one(int base);
  void  parse_sinar_ia();
  void  parse_rollei();
  void  parse_ciff(int offset,int length);
  void  ciff_block_1030();
  void  parse_external_jpeg();
  int   parse_tiff(int base);
  void  apply_tiff();
  void  parse_minolta(int base);
  void  parse_kodak_ifd(int base);
  void  linear_table(unsigned len);
  void  parse_mos(int offset);
  void  romm_coeff(float romm_cam[3][3]);
  void  parse_exif(int base);
  void  get_timestamp(int reversed);
  void  parse_makernote(int base,int uptag);
  int   parse_tiff_ifd(int base);
  void  tiff_get(unsigned base,unsigned *tag,unsigned *type,unsigned *len,unsigned *save);
  void  ahd_interpolate();
  void  ppg_interpolate();
  void  vng_interpolate();
  void  lin_interpolate();
  void  border_interpolate(int border);
  void  pre_interpolate();
  void  hat_transform(float *temp,float *base,int st,int size,int sc);
  void  pseudoinverse(double(*in)[3],double(*out)[3],int size);
  void  gamma_curve (double pwr, double ts, int mode, int imax);
  void  subtract(const char *fname);
  void  bad_pixels(const char *fname);
  void  canon_sraw_load_raw();
  void  lossless_jpeg_load_raw();
  uint16_t * ljpeg_row(int jrow,struct pt_jhead *jh);
  void  ljpeg_end (struct pt_jhead *jh);
  int   ljpeg_diff (uint16_t *huff);
  int   ljpeg_start(struct pt_jhead *jh,int info_only);
  void  canon_compressed_load_raw();
  int   canon_has_lowbits();
  void  crw_init_tables(unsigned table,uint16_t *huff[2]);
  uint8_t * make_decoder(const uint8_t *source,int level);
  void  init_decoder();
  unsigned  getbits(int nbits);
  void  canon_a5_load_raw();
  unsigned getbithuff(int nbits,uint16_t *huff);
  int   canon_s2is();
  void  remove_zeroes();
  void  canon_600_load_raw();
  void  canon_600_correct();
  void  canon_600_coeff();
  void  canon_600_auto_wb();
  int   canon_600_color(int ratio[2],int mar);
  void  canon_600_fixed_wb(int temp);
  void  canon_load_raw();
  void  lossless_dng_load_raw();
  void  lossy_dng_load_raw();
  void  packed_dng_load_raw();
  void  foveon_dp_load_raw();
  void  foveon_sd_load_raw();
  void  crop_masked_pixels();
  /* void  canon_black (double dark[2],int nblack); */
  void  read_shorts(uint16_t *pixel,int count);
  double  getreal(int type);
  float  int_to_float(int i);
  unsigned  getint(int type);
  unsigned  get4();
  unsigned  sget4(uint8_t *s);
  uint16_t  get2();
  uint16_t  sget2(uint8_t *s);
  void  derror();
  void  merror(void *ptr,const char *where);
  #if !defined(__GLIBC__)
  char *my_memmem(char *haystack,size_t haystacklen,char *neepte,size_t neeptelen);
  #endif
  int  fcol(int row,int col);
  // functions for dcb
  void dcb_pp_old();
  void copy_to_buffer_old(float (*m_Image2)[3]);
  void restore_from_buffer_old(float (*m_Image2)[3]);
  void hid_old();
  void hid2_old();
  void dcb_color_old();
  void dcb_map_old();
  void dcb_map2_old();
  void dcb_correction_old();
  void dcb_correction2_old();
  void dcb_correction3_old();
  void dcb_interpolate_soft_old(int iterations, int dcb_smoothing);
  void dcb_interpolate_sharp_old(int iterations, int dcb_smoothing);
  // functions for new dcb
  void dcb_pp();
  void copy_to_buffer(float (*m_Image2)[3]);
  void restore_from_buffer(float (*m_Image2)[3]);
  void hid();
  void hid2();
  void dcb_color();
  void dcb_color_full();
  void dcb_map();
  void dcb_correction();
  void dcb_correction2();
  void dcb_refinement();
  void fbdd_green();
  void fbdd_green2();
  void fbdd_color();
  void fbdd_correction();
  void fbdd_correction2(double (*m_Image3)[3]);
  void fbdd(int noiserd);
  void dcb(int iterations, int dcb_enhance);
  void lch_to_rgb(double (*m_Image3)[3]);
  void rgb_to_lch(double (*m_Image3)[3]);
  // functions for vcd
  void ahd_interpolate_mod();
  void ahd_partial_interpolate(int threshold_value);
  void es_median_filter();
  void median_filter_new();
  void refinement();
  void vcd_interpolate(int ahd_cutoff);
  // functions for algorithms from perfect raw
  void lmmse_interpolate(int gamma_apply);
  void amaze_demosaic();
  void cfa_linedn(float noise);
  void ddct8x8s(int isgn, float **a);
  void CA_correct(double cared, double cablue);
  void green_equilibrate(float thresh);
};

#endif

////////////////////////////////////////////////////////////////////////////////
