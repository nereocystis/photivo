/*******************************************************************************
**
** Photivo
**
** Copyright (C) 2008,2009 Jos De Laender <jos.de_laender@telenet.be>
** Copyright (C) 2009-2012 Michael Munzert <mail@mm-log.com>
** Copyright (C) 2010-2012 Bernd Schoeler <brjohn@brother-john.net>
** Copyright (C) 2013 Alexander Tzyganenko <tz@fast-report.com>
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
#include "ptDefines.h"
#include "ptInfo.h"
#include "ptCalloc.h"
#include "ptConfirmRequest.h"
#include "ptConstants.h"
#include "ptMessageBox.h"
#include "ptDcRaw.h"
#include "ptProcessor.h"
#include "ptMainWindow.h"
#include "ptViewWindow.h"
#include "ptHistogramWindow.h"
#include "ptGuiOptions.h"
#include "ptSettings.h"
#include "ptError.h"
#include "ptRGBTemperature.h"
#include "ptWhiteBalances.h"
#include "ptCurve.h"
#include "ptTheme.h"
#include "ptWiener.h"
#include "ptParseCli.h"
#include "ptImageHelper.h"
#include "filters/imagespot/ptTuningSpot.h"
#include "qtsingleapplication/qtsingleapplication.h"
#include "filemgmt/ptFileMgrWindow.h"
#include "batch/ptBatchWindow.h"
#include "filters/ptFilterDM.h"
#include "filters/ptFilterBase.h"
#include "ptToolBox.h"
#include "filters/ptFilterUids.h"

#ifdef Q_OS_WIN
  #include "ptEcWin7.h"
  #include "ptWinApi.h"
#endif

#include <wand/magick_wand.h>

#define QT_CLEAN_NAMESPACE

#include <QFileInfo>
#include <QtGui>
#include <QtCore>
#include <QFileDialog>
#include <QColorDialog>
#include <QInputDialog>
#include <QTextCodec>
#include <QDesktopWidget>
#ifdef Q_OS_MAC
  #include <QFileOpenEvent>
#endif

#include <string>
#include <csignal>
#include <exception>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// This is the file where everything is started.
// It starts not only the gui but it is *also* the starting point
// for guiless jobmode. The design is such that the pipe runs
// independent whether or not gui thingies are attached to it.
//
////////////////////////////////////////////////////////////////////////////////

ptDcRaw*     TheDcRaw       = NULL;
ptProcessor* TheProcessor   = NULL;

QStringList FileExtsRaw;
QStringList FileExtsBitmap;
QString     UserDirectory;
QString     ShareDirectory;
// Sidecar related
Exiv2::IptcData IptcData;
Exiv2::XmpData  XmpData;

cmsHPROFILE PreviewColorProfile = NULL;
cmsCIExyY       D65;
cmsCIExyY       D50;
// precalculated color transform
cmsHTRANSFORM ToPreviewTransform = NULL;
void ReadSidecar(const QString& Sidecar);
void SetRatingFromXmp();
void SetTagsFromXmp();

//
// The 'tee' towards the display.
// Visualization :
//
// Raw->Image_AfterDcRaw->Image_AfterLensfun->Image_AfterRGB->Image_AfterLab->
// Image_AfterGREYC->Image_AfterEyeCandy
// -------------------------------------------------------------------------
//                                       |
//                                Somewhere a tee to the preview
//                                       |
//                                  PreviewImage
//
// The pipe changed, Lab and Greyc became LabCC and LabSN


ptImage*  PreviewImage     = NULL;
ptImage*  HistogramImage   = NULL;

// The main windows of the application.
ptMainWindow*      MainWindow      = NULL;
ptViewWindow*      ViewWindow      = NULL;
ptHistogramWindow* HistogramWindow = NULL;
ptFileMgrWindow*   FileMgrWindow   = NULL;
ptBatchWindow*     BatchWindow     = NULL;

// Error dialog for segfaults
ptMessageBox* SegfaultErrorBox;

// Theming
ptTheme* Theme = NULL;

QTranslator appTranslator;

// Gui options and settings.
ptGuiOptions  *GuiOptions = NULL;
ptSettings    *Settings = NULL;

// Lensfun database.
//ptLensfun*  LensfunData = NULL;    // TODO BJ: implement lensfun DB

// Run mode
short NextPhase;
short NextSubPhase;
short ImageSaved;
short ImageCleanUp;

// uint16_t (0,0xffff) to float (0.0, 1.0)
float    ToFloatTable[0x10000];
float    ToFloatABNeutral[0x10000];
uint16_t ToInvertTable[0x10000];
uint16_t ToSRGBTable[0x10000];

// Filter patterns for the filechooser.
QString ChannelMixerFilePattern;
QString CurveFilePattern;
QString JobFilePattern;
QString SettingsFilePattern;
QString ProfilePattern;
QString RawPattern;
QString BitmapPattern;
QString SaveBitmapPattern;

void InitStrings() {
  ChannelMixerFilePattern =
    QCoreApplication::translate("Global Strings","Photivo channelmixer file (*.ptm);;All files (*.*)");
  CurveFilePattern =
    QCoreApplication::translate("Global Strings","Photivo curve file (*.ptc);;All files (*.*)");
  JobFilePattern =
    QCoreApplication::translate("Global Strings","Photivo job file (*.ptj);;All files (*.*)");
  SettingsFilePattern =
    QCoreApplication::translate("Global Strings","Photivo settings file (*.pts);;All files (*.*)");
  ProfilePattern =
    QCoreApplication::translate("Global Strings","ICC colour profiles (*.icc *.icm);;All files (*.*)");

  // QFileDialog has no case insensitive option ...
  RawPattern =
    QCoreApplication::translate("Global Strings","Raw files ("
                                                 "*.arw *.ARW *.Arw "
                                                 "*.bay *.BAY *.Bay "
                                                 "*.bmq *.BMQ *.Bmq "
                                                 "*.cr2 *.CR2 *.Cr2 "
                                                 "*.crw *.CRW *.Crw "
                                                 "*.cs1 *.CS1 *.Cs1 "
                                                 "*.dc2 *.DC2 *.Dc2 "
                                                 "*.dcr *.DCR *.Dcr "
                                                 "*.dng *.DNG *.Dng "
                                                 "*.erf *.ERF *.Erf "
                                                 "*.fff *.FFF *.Fff "
                                                 "*.hdr *.HDR *.Hdr "
                                                 "*.ia  *.IA *.Ia "
                                                 "*.k25 *.K25 "
                                                 "*.kc2 *.KC2 *.Kc2 "
                                                 "*.kdc *.KDC *.Kdc "
                                                 "*.mdc *.MDC *.Mdc "
                                                 "*.mef *.MEF *.Mef "
                                                 "*.mos *.MOS *.Mos "
                                                 "*.mrw *.MRW *.Mrw "
                                                 "*.nef *.NEF *.Nef "
                                                 "*.nrw *.NRW *.Nrw "
                                                 "*.orf *.ORF *.Orf "
                                                 "*.pef *.PEF *.Pef "
                                                 "*.pxn *.PXN *.Pxn "
                                                 "*.qtk *.QTK *.Qtk "
                                                 "*.raf *.RAF *.Raf "
                                                 "*.raw *.RAW *.Raw "
                                                 "*.rdc *.RDC *.Rdc "
                                                 "*.rw2 *.RW2 *.Rw2 "
                                                 "*.sr2 *.SR2 *.Sr2 "
                                                 "*.srf *.SRF *.Srf "
                                                 "*.srw *.SRW *.Srw "
                                                 "*.sti *.STI *.Sti "
                                                 "*.tif *.TIF *.Tif "
                                                 "*.x3f *.X3F *.X3f)"
                                                 ";;Bitmaps ("
                                                 "*.jpeg *.JPEG *.Jpeg "
                                                 "*.jpg *.JPG *.Jpg "
                                                 "*.tiff *.TIFF *.Tiff "
                                                 "*.tif *.TIF *.Tif "
                                                 "*.bmp *.BMP *.Bmp "
                                                 "*.png *.PNG *.Png "
                                                 "*.ppm *.PPm *.Ppm)"
                                                 ";;All files (*.*)");

  BitmapPattern =
    QCoreApplication::translate("Global Strings","Bitmaps ("
                                                 "*.jpeg *.JPEG *.Jpeg "
                                                 "*.jpg *.JPG *.Jpg "
                                                 "*.tiff *.TIFF *.Tiff "
                                                 "*.tif *.TIF *.Tif "
                                                 "*.bmp *.BMP *.Bmp "
                                                 "*.png *.PNG *.Png "
                                                 "*.ppm *.PPm *.Ppm "
                                                 ";;All files (*.*)");

  SaveBitmapPattern =
    QCoreApplication::translate("Global Strings","Jpeg (*.jpg);;"
                                                 "Tiff (*.tiff);;"
                                                 "Png (*.png);;"
                                                 "All files (*.*)");
}

////////////////////////////////////////////////////////////////////////////////
//
// Some function prototypes.
//
////////////////////////////////////////////////////////////////////////////////

ptImageType CheckImageType(QString filename,
                           uint16_t* width, uint16_t* height,
                           ptDcRaw* dcRaw = NULL);
void   RunJob(const QString FileName);
short  ReadJobFile(const QString FileName);
void   WriteOut();
void   UpdatePreviewImage(const ptImage* ForcedImage   = NULL,
                          const short    OnlyHistogram = 0);
void   UpdateCropToolUI();
void   PreCalcTransforms();
void   CB_ZoomFitButton();
void   CB_MenuFileOpen(const short HaveFile);
void   CB_MenuFileExit(const short);
void   ExtFileOpen(const QString file);
void   CB_WritePipeButton();
void   CB_OpenPresetFileButton();
void   CB_OpenSettingsFileButton();
void   CB_CropOrientationButton();
//short  WriteSettingsFile(const QString FileName, const short IsJobFile = 0);
void   SetBackgroundColor(int SetIt);
void   CB_StyleChoice(const QVariant Choice);
void   CB_SliderWidthInput(const QVariant Value);
void Export(const short mode);
void Update(short Phase,
            short SubPhase      = -1,
            short WithIdentify  = 1,
            short ProcessorMode = ptProcessorMode_Preview);
int CalculatePipeSize(const bool NewImage = false);
void CB_OpenSettingsFile(QString SettingsFileName);
void SaveButtonToolTip(const short mode);

int    photivoMain(int Argc, char *Argv[]);
void   CleanupResources();
void copyFolder(QString sourceFolder, QString destFolder);
void CB_PixelReader(const QPointF Point, const ptPixelReading PixelReading);
bool GBusy = false;

// undo-redo & clipboard support
QByteArray ptSettingsToQByteArray();
void ptQByteArrayToSettings(const QByteArray &arr);
void ptAddUndo();
void ptMakeUndo();
void ptMakeRedo();
void ptClearUndoRedo();
void ptMakeFullUndo();
void ptResetSettingsToDefault();
void ptCopySettingsToClipboard();
void ptPasteSettingsFromClipboard();

//==============================================================================

void CreateAllFilters() {
  //                   Filter ID                unique name                         caption postfix
  // Local Edit tab
  GFilterDM->NewFilter("SpotTuning",            Fuid::SpotTuning_Local);
  // RGB tab
  GFilterDM->NewFilter("ChannelMixer",          Fuid::ChannelMixer_RGB);
  GFilterDM->NewFilter("Highlights",            Fuid::Highlights_RGB);
  GFilterDM->NewFilter("ColorIntensity",        Fuid::ColorIntensity_RGB);
  GFilterDM->NewFilter("Brightness",            Fuid::Brightness_RGB);
  GFilterDM->NewFilter("ExposureCorrection",    Fuid::Exposure_RGB);
  GFilterDM->NewFilter("ReinhardBrighten",      Fuid::ReinhardBrighten_RGB);
  GFilterDM->NewFilter("GammaTool",             Fuid::GammaTool_RGB);
  GFilterDM->NewFilter("Normalization",         Fuid::Normalization_RGB);
  GFilterDM->NewFilter("ColorEnhancement",      Fuid::ColorEnhancement_RGB);
  GFilterDM->NewFilter("LMHRecoveryRgb",        Fuid::LMHRecovery_RGB);
  GFilterDM->NewFilter("TextureContrastRgb",    Fuid::TextureContrast_RGB);
  GFilterDM->NewFilter("LocalContrastRgb",      Fuid::LocalContrast1_RGB,           " I");
  GFilterDM->NewFilter("LocalContrastRgb",      Fuid::LocalContrast2_RGB,           " II");
  GFilterDM->NewFilter("SigContrastRgb",        Fuid::SigContrastRgb_RGB);
  GFilterDM->NewFilter("LevelsRgb",             Fuid::Levels_RGB);
  GFilterDM->NewFilter("RgbCurve",              Fuid::RgbCurve_RGB);
  // Lab Color/contrast tab
  GFilterDM->NewFilter("LabTransform",          Fuid::LabTransform_LabCC);
  GFilterDM->NewFilter("ShadowsHighlights",     Fuid::ShadowsHighlights_LabCC);
  GFilterDM->NewFilter("LMHRecoveryLab",        Fuid::LMHRecovery_LabCC);
  GFilterDM->NewFilter("Drc",                   Fuid::Drc_LabCC);
  GFilterDM->NewFilter("SigContrastLab",        Fuid::SigContrastLab_LabCC);
  GFilterDM->NewFilter("TextureCurve",          Fuid::TextureCurve_LabCC);
  GFilterDM->NewFilter("TextureContrastLab",    Fuid::TextureContrast1_LabCC,       " I");
  GFilterDM->NewFilter("TextureContrastLab",    Fuid::TextureContrast2_LabCC,       " II");
  GFilterDM->NewFilter("LocalContrastLab",      Fuid::LocalContrast1_LabCC,         " I");
  GFilterDM->NewFilter("LocalContrastLab",      Fuid::LocalContrast2_LabCC,         " II");
  GFilterDM->NewFilter("LocalContrastStretch",  Fuid::LocalContrastStretch1_LabCC,  " I");
  GFilterDM->NewFilter("LocalContrastStretch",  Fuid::LocalContrastStretch2_LabCC,  " II");
  GFilterDM->NewFilter("Saturation",            Fuid::Saturation_LabCC);
  GFilterDM->NewFilter("ColorBoost",            Fuid::ColorBoost_LabCC);
  GFilterDM->NewFilter("LevelsLab",             Fuid::Levels_LabCC);
  // Lab sharpen/noise tab
  GFilterDM->NewFilter("ImpulseNR",             Fuid::ImpulseNR_LabSN);
  GFilterDM->NewFilter("EAWavelets",            Fuid::EAWavelets_LabSN);
  GFilterDM->NewFilter("GreyCStoration",        Fuid::GreyCStoration_LabSN);
  GFilterDM->NewFilter("Defringe",              Fuid::Defringe_LabSN);
  GFilterDM->NewFilter("WaveletDenoise",        Fuid::WaveletDenoise_LabSN);
  GFilterDM->NewFilter("LumaDenoise",           Fuid::LumaDenoise_LabSN);
  GFilterDM->NewFilter("LumaDenoiseCurve",      Fuid::LumaDenoiseCurve_LabSN);
  GFilterDM->NewFilter("LumaDenoiseCurve",      Fuid::LumaDenoiseCurve2_LabSN,      " II");
  GFilterDM->NewFilter("PyramidDenoise",        Fuid::PyramidDenoise_LabSN);
  GFilterDM->NewFilter("ColorDenoise",          Fuid::ColorDenoise_LabSN);
  GFilterDM->NewFilter("GradientSharpen",       Fuid::GradientSharpen_LabSN);
  GFilterDM->NewFilter("DetailCurve",           Fuid::DetailCurve_LabSN);
  GFilterDM->NewFilter("Wiener",                Fuid::Wiener_LabSN);
  GFilterDM->NewFilter("InvDiffSharpen",        Fuid::InvDiffSharpen_LabSN);
  GFilterDM->NewFilter("UnsharpMask",           Fuid::Usm_LabSN);
  GFilterDM->NewFilter("HighpassSharpen",       Fuid::HighpassSharpen_LabSN);
  GFilterDM->NewFilter("FilmGrain",             Fuid::FilmGrain_LabSN);
  GFilterDM->NewFilter("LabChannelView",        Fuid::ViewLab_LabSN);
  // Lab Eyecandy tab
  GFilterDM->NewFilter("Outline",               Fuid::Outline_LabEyeCandy);
  GFilterDM->NewFilter("LumaByHueCurve",        Fuid::LumaByHueCurve_LabEyeCandy);
  GFilterDM->NewFilter("SatCurve",              Fuid::SatCurve_LabEyeCandy);
  GFilterDM->NewFilter("HueCurve",              Fuid::HueCurve_LabEyeCandy);
  GFilterDM->NewFilter("LCurve",                Fuid::LCurve_LabEyeCandy);
  GFilterDM->NewFilter("ABCurves",              Fuid::ABCurves_LabEyeCandy);
  GFilterDM->NewFilter("ColorContrast",         Fuid::ColorContrast_LabEyeCandy);
  GFilterDM->NewFilter("ToneAdjust",            Fuid::ToneAdjust1_LabEyeCandy,      " I");
  GFilterDM->NewFilter("ToneAdjust",            Fuid::ToneAdjust2_LabEyeCandy,      " II");
  GFilterDM->NewFilter("LumaAdjust",            Fuid::LumaAdjust_LabEyeCandy);
  GFilterDM->NewFilter("SatAdjust",             Fuid::SatAdjust_LabEyeCandy);
  GFilterDM->NewFilter("Tone",                  Fuid::Tone_LabEyeCandy);
  GFilterDM->NewFilter("VignetteLab",           Fuid::Vignette_LabEyeCandy);
  // Eyecandy tab
  GFilterDM->NewFilter("BlackWhite",            Fuid::BlackWhite_EyeCandy);
  GFilterDM->NewFilter("CrossProcessing",       Fuid::CrossProcessing_EyeCandy);
  GFilterDM->NewFilter("SimpleTone",            Fuid::SimpleTone_EyeCandy);
  GFilterDM->NewFilter("ColorTone",             Fuid::ColorTone1_EyeCandy,          " I");
  GFilterDM->NewFilter("ColorTone",             Fuid::ColorTone2_EyeCandy,          " II");
  GFilterDM->NewFilter("SigContrastRgb",        Fuid::SigContrastRgb_EyeCandy);
  GFilterDM->NewFilter("TextureOverlay",        Fuid::TextureOverlay1_EyeCandy,     " I");
  GFilterDM->NewFilter("TextureOverlay",        Fuid::TextureOverlay2_EyeCandy,     " II");
  GFilterDM->NewFilter("GradualOverlay",        Fuid::GradualOverlay1_EyeCandy,     " I");
  GFilterDM->NewFilter("GradualOverlay",        Fuid::GradualOverlay2_EyeCandy,     " II");
  GFilterDM->NewFilter("VignetteRgb",           Fuid::Vignette_EyeCandy);
  GFilterDM->NewFilter("GradualBlur",           Fuid::GradualBlur1_EyeCandy,        " I");
  GFilterDM->NewFilter("GradualBlur",           Fuid::GradualBlur2_EyeCandy,        " II");
  GFilterDM->NewFilter("SoftglowOrton",         Fuid::SoftglowOrton_EyeCandy);
  GFilterDM->NewFilter("ColorIntensity",        Fuid::ColorIntensity_EyeCandy);
  GFilterDM->NewFilter("RToneCurve",            Fuid::RTone_EyeCandy);
  GFilterDM->NewFilter("GToneCurve",            Fuid::GTone_EyeCandy);
  GFilterDM->NewFilter("BToneCurve",            Fuid::BTone_EyeCandy);
  // Output tab
  GFilterDM->NewFilter("RgbCurve",              Fuid::RgbCurve_Out);
  GFilterDM->NewFilter("AfterGammaCurve",       Fuid::AfterGammaCurve_Out);
  GFilterDM->NewFilter("SigContrastRgb",        Fuid::SigContrastRgb_Out);
  GFilterDM->NewFilter("Wiener",                Fuid::Wiener_Out);
}


//==============================================================================


////////////////////////////////////////////////////////////////////////////////
//
// Progress function in the GUI.
// Can later also in job.
//
////////////////////////////////////////////////////////////////////////////////

void ReportProgress(const QString Message) {
  printf("Progress : %s\n",Message.toLocal8Bit().data());
  if (!MainWindow) return;
  MainWindow->StatusLabel->setText(Message);
  MainWindow->StatusLabel->repaint();
  // Workaround to keep the GUI responsive
  // during pipe processing...
  QApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

////////////////////////////////////////////////////////////////////////////////
//
// Main : instantiating toplevel windows, settings and options ..
//
////////////////////////////////////////////////////////////////////////////////

short   InStartup  = 1;

short   JobMode = 0;
QString JobFileName = "";
QString ImageFileToOpen = "";
QString PtsFileToOpen = "";
QString Sidecar = "";
#ifdef Q_OS_MAC
    bool MacGotFileEvent=false;
#endif
QRect DetailViewRect;


#ifndef DLRAW_GIMP_PLUGIN
void SegfaultAbort(int) {
  // prevent infinite recursion if SegfaultAbort() causes another segfault
  std::signal(SIGSEGV, SIG_DFL);
  std::signal(SIGABRT, SIG_DFL);
  SegfaultErrorBox->exec();
//  std::abort();
//  Batch manager needs Photivo to return a value to determine a crash
  exit(EXIT_FAILURE);
}

int main(int Argc, char *Argv[]) {
#ifdef Q_OS_WIN
  WinApi::AttachToParentConsole();
#endif

  int RV = photivoMain(Argc,Argv);
  DestroyMagick();
  DelAndNull(SegfaultErrorBox);
  CleanupResources(); // Not necessary , for debug.

#ifdef Q_OS_WIN
  // Close output channels to a possibly attached console and detach from it.
  // Not sure if this is strictly necessary, but it does not hurt either.
  fclose(stdout);
  fclose(stderr);
  FreeConsole();
#endif

  return RV;
}
#endif


//class QtSingleApplication with redefined event to handle QFileOpenEvent
#ifdef Q_OS_MAC
  class MyQApplication : public QtSingleApplication {
  protected:
    virtual bool event(QEvent *event);
  public:
      bool initialized;
      void macinit();
      MyQApplication(const QString &appId, int & argc, char ** arg);
  };

  bool MyQApplication::event(QEvent *event)
  {
      switch (event->type()) {
      case QEvent::FileOpen:
          ImageFileToOpen=static_cast<QFileOpenEvent *>(event)->file();
          MacGotFileEvent=true;
          if(initialized){
              //skip on first Event
          ExtFileOpen(ImageFileToOpen);
          }
          return true;
      default:
          MacGotFileEvent=false;
          break;
      }
      return QtSingleApplication::event(event);
  }
  MyQApplication::MyQApplication(const QString &appId,int & argc, char ** argv) : QtSingleApplication(appId,argc, argv){
  initialized=false;
  };
  void MyQApplication::macinit(){
      initialized=true;
  };
#endif


#ifdef Q_OS_MAC
MyQApplication* TheApplication;
#else
QtSingleApplication* TheApplication;
#endif

int photivoMain(int Argc, char *Argv[]) {
  QString VerTemp(TOSTRING(APPVERSION));    //also used for the cli syntax error msg below!
  printf("Photivo version %s\n", VerTemp.toLocal8Bit().data());

  InitializeMagick(*Argv);

  //QApplication TheApplication(Argc,Argv);
  QStringList environment = QProcess::systemEnvironment();
  QString user=environment.filter(QRegExp("^USERNAME=|^USER=",Qt::CaseInsensitive)).first();
  if(!user.isEmpty()){
      int l=user.indexOf("=",0);
      user="_"+user.right(user.length()-l-1);
  }
  user="photivo"+user;

#ifdef Q_OS_MAC
  TheApplication = new MyQApplication(user, Argc,Argv);
#else
  TheApplication = new QtSingleApplication(user, Argc,Argv);
#endif
  //close on last window closed, seems to be unneeded
  TheApplication->connect(TheApplication, SIGNAL(lastWindowClosed()), TheApplication, SLOT(quit()));

  // Error handling for segfaults (to prevent silent crashes)
  SegfaultErrorBox = new ptMessageBox(QMessageBox::Critical,
                                      QObject::tr("Photivo crashed"),
                                      QObject::tr("Photivo crashed. You can get help on our flickr forum at\n"
                                      "<a href=\"http://www.flickr.com/groups/photivo/discuss/\">http://www.flickr.com/groups/photivo/discuss/</a>"
                                      "\nWhen you post there make sure to describe your last actions before the crash occurred."),
                                      QMessageBox::Ok);
  SegfaultErrorBox->setTextFormat(Qt::RichText);
  std::signal(SIGSEGV, SegfaultAbort);
  std::signal(SIGABRT, SegfaultAbort);

  // Check for wrong GM quantum depth. We need the 16bit GraphicsMagick.
  ulong QDepth = 0;
  const ulong MinQD = 16;
  MagickGetQuantumDepth(&QDepth);
  if (QDepth < MinQD) {
    QString WrongQDepthMsg = QObject::tr(
        "Fatal error: Wrong GraphicsMagick quantum depth!\n"
        "Found quantum depth %1. Photivo needs at least %2.\n")
        .arg(QDepth).arg(MinQD);
    fprintf(stderr,"%s", WrongQDepthMsg.toLocal8Bit().data());
    ptMessageBox::critical(0, QObject::tr("Photivo: Fatal Error"), WrongQDepthMsg);
    exit(EXIT_FAILURE);
  }


  // Handle cli arguments
  QString PhotivoCliUsageMsg = "<pre>" + QObject::tr(
"Syntax: photivo [inputfile | -i imagefile | -j jobfile |\n"
"                 --load-and-delete imagefile]\n"
"                [--pts ptsfile] [--sidecar sidecarfile] [-h] [--new-instance]\n"
"Options:\n"
"inputfile\n"
"      Specify the image or settings file to load. Works like -i for image files\n"
"      and like --pts for settings files.\n"
"-i imagefile\n"
"      Specify image file to load.\n"
"-j jobfile\n"
"      Specify jobfile for batch processing. Job files are created\n in Photivo\n"
"      and then executed with this option.\n"
"--load-and-delete imagefile\n"
"      Specify temporary file used for Gimp-to-Photivo export. Internal option,\n"
"      not intended for general use. BEWARE! This option deletes imagefile!\n"
"--pts ptsfile\n"
"      Specify settings file to load with the image. Must be used together\n"
"      with -i.\n"
"--sidecar sidecarfile\n"
"      Specify sidecar file to load with the image.\n"
"--new-instance\n"
"      Allow opening another Photivo instance instead of using a currently\n"
"      running Photivo. Job files are always opened in a new instance.\n"
"--no-fmgr or -p\n"
"      Prevent auto-open file manager when Photivo starts.\n"
"--help or -h\n"
"      Display this usage information.\n\n"
"For more documentation visit the wiki: http://photivo.org/photivo/start\n"
  ) + "</pre>";

  ptCliCommands cli = { cliNoAction, "", "", "", false, false };

#ifdef Q_OS_MAC
//Just Skip if engaged by QFileOpenEvent
  if(!MacGotFileEvent) {
#endif

  cli = ParseCli(Argc, Argv);

  // Show help message and exit Photivo
  if (cli.Mode == cliShowHelp) {
  #ifdef Q_OS_WIN
    ptMessageBox::critical(0, QObject::tr("Photivo"), PhotivoCliUsageMsg);
  #else
    fprintf(stderr,"%s",PhotivoCliUsageMsg.toLocal8Bit().data());
  #endif
    exit(EXIT_FAILURE);
  }

  ImageCleanUp = cli.Mode == cliLoadAndDelImage;
  JobMode = cli.Mode == cliProcessJob;
  PtsFileToOpen = cli.PtsFilename;

  if (JobMode) {
    JobFileName = cli.Filename;
  } else {
    ImageFileToOpen = cli.Filename;
  }

  
  if (cli.Sidecar != "" ) {
    Sidecar = cli.Sidecar;
  }

  // QtSingleInstance, add CLI-Switch to skip and allow multiple instances
  // JobMode is always run in a new instance
  // Sent messages are handled by ptMainWindow::OtherInstanceMessage
  if (!JobMode) {
    if (TheApplication->isRunning() && !cli.NewInstance) {
      if (PtsFileToOpen != "") {
        TheApplication->sendMessage("::pts::" + PtsFileToOpen);
      }
      if (ImageFileToOpen != "") {
        if (ImageCleanUp > 0) {
          TheApplication->sendMessage("::tmp::" + ImageFileToOpen);
        } else {
          TheApplication->sendMessage("::img::" + ImageFileToOpen);
        }
      }
      if (Sidecar != ""){
        TheApplication->sendMessage("::sidecar::" + Sidecar);
      }
      TheApplication->activateWindow();
      exit(0);
    }
  }


#ifdef Q_OS_MAC
  } // !MacGotFileEvent
  QDir dir(QApplication::applicationDirPath());
  QApplication::setLibraryPaths(QStringList(dir.absolutePath()));
#endif



  FileExtsRaw << "*.arw" << "*.ARW" << " *.Arw"
              << "*.bay" << "*.BAY" << "*.Bay"
              << "*.bmq" << "*.BMQ" << "*.Bmq"
              << "*.cr2" << "*.CR2" << "*.Cr2"
              << "*.crw" << "*.CRW" << "*.Crw"
              << "*.cs1" << "*.CS1" << "*.Cs1"
              << "*.dc2" << "*.DC2" << "*.Dc2"
              << "*.dcr" << "*.DCR" << "*.Dcr"
              << "*.dng" << "*.DNG" << "*.Dng"
              << "*.erf" << "*.ERF" << "*.Erf"
              << "*.fff" << "*.FFF" << "*.Fff"
              << "*.hdr" << "*.HDR" << "*.Hdr"
              << "*.ia " << "*.IA" << "*.Ia"
              << "*.k25" << "*.K25"
              << "*.kc2" << "*.KC2" << "*.Kc2"
              << "*.kdc" << "*.KDC" << "*.Kdc"
              << "*.mdc" << "*.MDC" << "*.Mdc"
              << "*.mef" << "*.MEF" << "*.Mef"
              << "*.mos" << "*.MOS" << "*.Mos"
              << "*.mrw" << "*.MRW" << "*.Mrw"
              << "*.nef" << "*.NEF" << "*.Nef"
              << "*.nrw" << "*.NRW" << "*.Nrw"
              << "*.orf" << "*.ORF" << "*.Orf"
              << "*.pef" << "*.PEF" << "*.Pef"
              << "*.pxn" << "*.PXN" << "*.Pxn"
              << "*.qtk" << "*.QTK" << "*.Qtk"
              << "*.raf" << "*.RAF" << "*.Raf"
              << "*.raw" << "*.RAW" << "*.Raw"
              << "*.rdc" << "*.RDC" << "*.Rdc"
              << "*.rw2" << "*.RW2" << "*.Rw2"
              << "*.sr2" << "*.SR2" << "*.Sr2"
              << "*.srf" << "*.SRF" << "*.Srf"
              << "*.srw" << "*.SRW" << "*.Srw"
              << "*.sti" << "*.STI" << "*.Sti"
              << "*.tif" << "*.TIF" << "*.Tif"
              << "*.x3f" << "*.X3F" << "*.X3f";

  FileExtsBitmap << "*.jpeg" << "*.JPEG" << "*.Jpeg"
                 << "*.jpg"  << "*.JPG"  << "*.Jpg"
                 << "*.png"  << "*.PNG"  << "*.Png"
                 << "*.tiff" << "*.TIFF" << "*.Tiff"
                 << "*.tif"  << "*.TIF"  << "*.Tif"
                 << "*.bmp"  << "*.BMP"  << "*.Bmp"
                 << "*.png"  << "*.PNG"  << "*.Png"
                 << "*.ppm"  << "*.PPm"  << "*.Ppm";


  // User home folder, where Photivo stores its ini and all Presets, Curves etc
  // %appdata%\Photivo on Windows, ~/.photivo on Linux or the program folder for the
  // portable Windows version.
  short IsPortableProfile = 0;
  QString AppDataFolder = "";
  QString Folder = "";
#ifdef Q_OS_WIN
  IsPortableProfile = QFile::exists("use-portable-profile");
  if (IsPortableProfile != 0) {
    printf("Photivo running in portable mode.\n");
    AppDataFolder = QCoreApplication::applicationDirPath();
    Folder = "";
  } else {
    // WinAPI returns path with native separators "\". We need to change this to "/" for Qt.
    AppDataFolder = WinApi::AppdataFolder();
    // Keeping the leading "/" separate here is important or mkdir will fail.
    Folder = "Photivo/";
  }
#else
  Folder = ".photivo/";
  AppDataFolder = QDir::homePath();
#endif

  UserDirectory = AppDataFolder + "/" + Folder;

  if (IsPortableProfile == 0) {
      QDir home(AppDataFolder);
      if (!home.exists(Folder))
          home.mkdir(Folder);
  }

  QString SettingsFileName = UserDirectory + "photivo.ini";
  // this has to be changed when we move to a different tree structure!
#ifdef __unix__
  QString NewShareDirectory(TOSTRING(PREFIX));
  if (NewShareDirectory.endsWith("/")) NewShareDirectory.chop(1);
  NewShareDirectory.append("/share/photivo/");
#else
  QString NewShareDirectory = QCoreApplication::applicationDirPath().append("/");
#endif
  ShareDirectory = NewShareDirectory;

  QFileInfo SettingsFileInfo(SettingsFileName);
//  short NeedInitialization = 1;
  short FirstStart = 1;
  if (SettingsFileInfo.exists() &&
          SettingsFileInfo.isFile() &&
          SettingsFileInfo.isReadable()) {
      // photivo was initialized
//      NeedInitialization = 0;
      FirstStart = 0;
      printf("Existing settingsfile '%s'\n",SettingsFileName.toLocal8Bit().data());
  } else {
      printf("New settingsfile '%s'\n",SettingsFileName.toLocal8Bit().data());
  }

  printf("User directory: '%s'; \n",UserDirectory.toLocal8Bit().data());
  printf("Share directory: '%s'; \n",NewShareDirectory.toLocal8Bit().data());

  // We need to load the translation before the ptSettings
  QSettings* TempSettings = new QSettings(SettingsFileName, QSettings::IniFormat);

//  if (TempSettings->value("SettingsVersion",0).toInt() < PhotivoSettingsVersion)
//      NeedInitialization = 1;

  // Initialize the user folder if needed
  /* TODO: for testing. Enable the other line below once profile versions are final. */
  if (IsPortableProfile == 0) {
      //if (NeedInitialization == 1 && IsPortableProfile == 0) {
      printf("Initializing/Updating user profile...\n");
      QFile::remove(UserDirectory + "photivo.png");
      QFile::copy(NewShareDirectory + "photivo.png",
              UserDirectory + "photivo.png");
//      QFile::remove(UserDirectory + "photivoLogo.png");
//      QFile::copy(NewShareDirectory + "photivoLogo.png",
//              UserDirectory + "photivoLogo.png");
      QFile::remove(UserDirectory + "photivoPreview.jpg");
      QFile::copy(NewShareDirectory + "photivoPreview.jpg",
              UserDirectory + "photivoPreview.jpg");
      QStringList SourceFolders;
      SourceFolders << NewShareDirectory + "Translations"
          << NewShareDirectory + "Curves"
          << NewShareDirectory + "ChannelMixers"
          << NewShareDirectory + "Presets"
          << NewShareDirectory + "Profiles"
          << NewShareDirectory + "LensfunDatabase"
          << NewShareDirectory + "UISettings";
      QStringList DestFolders;
      DestFolders << UserDirectory + "Translations"
          << UserDirectory + "Curves"
          << UserDirectory + "ChannelMixers"
          << UserDirectory + "Presets"
          << UserDirectory + "Profiles"
          << UserDirectory + "LensfunDatabase"
          << UserDirectory + "UISettings";

      for (int i = 0; i < SourceFolders.size(); i++) {
          copyFolder(SourceFolders.at(i), DestFolders.at(i));
      }
  }


  // Load Translation
  int TranslMode = TempSettings->value("TranslationMode",0).toInt();
  QDir TranslDir(UserDirectory + "Translations");
  QStringList UiLanguages = TranslDir.entryList(QStringList("photivo_*.qm"), QDir::Files|QDir::Readable, QDir::Name).replaceInStrings(".qm", "", Qt::CaseInsensitive);
  UiLanguages.replaceInStrings("photivo_", "", Qt::CaseInsensitive);
  int LangIdx = -1;

  if (TranslMode == 1) {
      LangIdx = UiLanguages.indexOf(TempSettings->value("UiLanguage","").toString());
      if (LangIdx >= 0) {
          QTranslator qtTranslator;
          appTranslator.load("photivo_" + UiLanguages[LangIdx], UserDirectory + "Translations");
          TheApplication->installTranslator(&appTranslator);
          qtTranslator.load("qt_" + UiLanguages[LangIdx], UserDirectory + "Translations");
          TheApplication->installTranslator(&qtTranslator);
          printf("Enabled translation: \"%s\".\n", UiLanguages[LangIdx].toLocal8Bit().data());
      }
  }

  delete TempSettings;

  // Persistent settings.
  // fixed the remember level to 2, since we have settings files now
  short RememberSettingLevel = 2;

  // Load the Settings (are also partly used in JobMode)
  Settings = new ptSettings(RememberSettingLevel, UserDirectory);

  // Set directories for needed files
  Settings->SetValue("UserDirectory", UserDirectory);
  Settings->SetValue("ShareDirectory",NewShareDirectory);
  Settings->SetValue("MainDirectory",QCoreApplication::applicationDirPath().append("/"));
  Settings->SetValue("Sidecar", Sidecar);

  // Set paths once with first start
  if (FirstStart == 1) {
      Settings->SetValue("RawsDirectory", UserDirectory);
      Settings->SetValue("OutputDirectory", UserDirectory);
      Settings->SetValue("UIDirectory", UserDirectory + "UISettings");
      Settings->SetValue("PresetDirectory", UserDirectory + "Presets");
      Settings->SetValue("CurvesDirectory", UserDirectory + "Curves");
      Settings->SetValue("ChannelMixersDirectory", UserDirectory + "ChannelMixers");
      Settings->SetValue("TranslationsDirectory", UserDirectory + "Translations");
      Settings->SetValue("CameraColorProfilesDirectory", UserDirectory + "Profiles/Camera");
      Settings->SetValue("PreviewColorProfilesDirectory", UserDirectory + "Profiles/Preview");
      Settings->SetValue("OutputColorProfilesDirectory", UserDirectory + "Profiles/Output");
      Settings->SetValue("StandardAdobeProfilesDirectory", UserDirectory + "Profiles/Camera/Standard");
      Settings->SetValue("LensfunDatabaseDirectory", UserDirectory + "LensfunDatabase");
      Settings->SetValue("PreviewColorProfile", UserDirectory + "Profiles/Preview/sRGB.icc");
      Settings->SetValue("OutputColorProfile", UserDirectory + "Profiles/Output/sRGB.icc");
      Settings->SetValue("StartupSettingsFile", UserDirectory + "Presets/MakeFancy.pts");
  }

  // Initialize patterns (after translation)
  InitStrings();

  // Init filter data module
  ptFilterDM::CreateInstance();
  CreateAllFilters();

  // Load also the LensfunDatabase.
  printf("Lensfun database: '%s'; \n",Settings->GetString("LensfunDatabaseDirectory").toLocal8Bit().data());
  //  LensfunData = new ptLensfun;    // TODO BJ: implement lensfun DB

  // Instantiate the processor. Spot models are not set here because we
  // do not yet know if we are in GUI or batch mode.
  TheProcessor = new ptProcessor(ReportProgress);

  // First check if we are maybe started as a command line with options.
  // (And thus have to run a batch job)

  if (JobMode) {
    RunJob(JobFileName);
    exit(EXIT_SUCCESS);
  }

  // If falling through to here we are in an interactive non-job mode.
  // Start op the Gui stuff.

  // Start the theme class
  Theme = new ptTheme(TheApplication,
                      (ptTheme::Theme)Settings->GetInt("Style"),
                      (ptTheme::Highlight)Settings->GetInt("StyleHighLight"));
  Theme->setCustomCSS(Settings->GetString("CustomCSSFile"));

  GuiOptions = new ptGuiOptions();

  // Open and keep open the profile for previewing.
  PreviewColorProfile = cmsOpenProfileFromFile(
          Settings->GetString("PreviewColorProfile").toLocal8Bit().data(),
          "r");
  if (!PreviewColorProfile) {
      ptLogError(ptError_FileOpen,
              Settings->GetString("PreviewColorProfile").toLocal8Bit().data());
      return ptError_FileOpen;
  }

  // When loading a file via cli, set file manager directory to that path.
  // Chances are good the user want to work with other files from that dir as well
  if (!JobMode && !ImageCleanUp && (ImageFileToOpen != "")) {
    Settings->SetValue("LastFileMgrLocation", QFileInfo(ImageFileToOpen).absolutePath());
  }

#ifndef PT_WITHOUT_FILEMGR
  Settings->SetValue("PreventFileMgrStartup", int(cli.NoOpenFileMgr));
#endif

  // Construct windows
  MainWindow = new ptMainWindow(QObject::tr("Photivo"));
  QObject::connect(TheApplication, SIGNAL(messageReceived(const QString &)),
                   MainWindow, SLOT(OtherInstanceMessage(const QString &)));
  TheApplication->setActivationWindow(MainWindow);

  ViewWindow =
      new ptViewWindow(MainWindow->ViewFrameCentralWidget, MainWindow);

  ViewWindow->SetPixelReader(CB_PixelReader);

  HistogramWindow =
      new ptHistogramWindow(NULL,MainWindow->HistogramFrameCentralWidget);

#ifndef PT_WITHOUT_FILEMGR
  FileMgrWindow = new ptFileMgrWindow(MainWindow->FileManagerPage);
  MainWindow->FileManagerLayout->addWidget(FileMgrWindow);
  QObject::connect(FileMgrWindow, SIGNAL(fileMgrWindowClosed()),
                   MainWindow, SLOT(CloseFileMgrWindow()));
  QObject::connect(ViewWindow, SIGNAL(openFileMgr()), MainWindow, SLOT(OpenFileMgrWindow()));
#endif
  QObject::connect(ViewWindow, SIGNAL(openBatch()), MainWindow, SLOT(OpenBatchWindow()));

  BatchWindow = new ptBatchWindow(MainWindow->BatchPage);
  MainWindow->BatchLayout->addWidget(BatchWindow);
  QObject::connect(BatchWindow, SIGNAL(BatchWindowClosed()), MainWindow, SLOT(CloseBatchWindow()));

  // Populate Translations combobox
  MainWindow->PopulateTranslationsCombobox(UiLanguages, LangIdx);

  // Theming
  CB_StyleChoice(Settings->GetInt("Style"));

  SetBackgroundColor(Settings->GetInt("BackgroundColor"));

  MainWindow->UpdateToolBoxes();


  //-------------------------------------
  // Initialize main window geometry and position
  // If the screen is not longer present, we move the window to the standard screen

  int      Screen        = Settings->m_IniSettings->value("MainWindowScreen", -1).toInt();
  if (Screen >= qApp->desktop()->screenCount()) {
    Screen = -1;
  }
  QRect    DesktopRect   = qApp->desktop()->screenGeometry(Screen);
  QPoint   UpperLeft     = DesktopRect.topLeft() + QPoint(30, 30);
  QPoint   MainWindowPos = Settings->m_IniSettings->value("MainWindowPos", UpperLeft).toPoint();
  QVariant WinSize       = Settings->m_IniSettings->value("MainWindowSize");
  QSize    MainWindowSize;

  if (!DesktopRect.contains(MainWindowPos)) {
    MainWindowPos = UpperLeft;
  }

  if(!WinSize.isValid() ||
     !(WinSize.toSize().width()  < DesktopRect.width()  + 50 &&
       WinSize.toSize().height() < DesktopRect.height() + 50   )) {
    // ensure a reasonable size if we don’t get one from Settings
    MainWindowSize = QSize(qMin(1200, (int)(DesktopRect.width()  * 0.8)),
                           qMin( 900, (int)(DesktopRect.height() * 0.8)) );
  } else {
    MainWindowSize = WinSize.toSize();
  }

  MainWindow->resize(MainWindowSize);
  MainWindow->move(  MainWindowPos);

  // MainSplitter is the one between tool pane and image pane.
  // ControlSplitter is the one between histogram and tools tabwidget.
  QVariant splitterState = Settings->m_IniSettings->value("MainSplitter");
  if (splitterState.isValid()) {
    MainWindow->MainSplitter->restoreState(splitterState.toByteArray());
  } else {
    // 10000 width for image pane to ensure the 300 for tool pane
    MainWindow->MainSplitter->setSizes(QList<int>() << 300 << 10000);
  }

  splitterState = Settings->m_IniSettings->value("ControlSplitter");
  if (splitterState.isValid()) {
    MainWindow->ControlSplitter->restoreState(splitterState.toByteArray());
  } else {
    MainWindow->ControlSplitter->setSizes(QList<int>() << 100 << 10000);
  }

  if (Settings->m_IniSettings->value("IsMaximized",0).toBool()) {
    MainWindow->showMaximized();
  } else {
    MainWindow->show();
  }
  //-------------------------------------


  // Update the preview image will result in displaying the splash.
  Update(ptProcessorPhase_Preview);

  // Open and keep open the profile for previewing.
  PreviewColorProfile = cmsOpenProfileFromFile(
          Settings->GetString("PreviewColorProfile").toLocal8Bit().data(),
          "r");
  if (!PreviewColorProfile) {
      ptLogError(ptError_FileOpen,
              Settings->GetString("PreviewColorProfile").toLocal8Bit().data());
      assert(PreviewColorProfile);
  }

  try {
    // Start event loops.
    return TheApplication->exec();
  } catch (const exception& E) {
    GInfo->Raise(E.what(), AT);
  } catch (...) {
    GInfo->Raise("Unknown error.", AT);
  }
  return 1;
}

  ////////////////////////////////////////////////////////////////////////////////
  //
  // This is only needed for the gimp integration.
  // Calling over and over photivoMain would otherwise leak like hell
  // (or any comparable place).
  //
  ////////////////////////////////////////////////////////////////////////////////

  void CleanupResources() {
      //printf("(%s,%d) qApp:%p\n",__FILE__,__LINE__,qApp);
      // Not : is done at CB for the exit.
      // delete Settings; // Don't, is done at CB_MenuFileExit
      // Also : do not delete items which are handled by MainWindow, such as
      // ViewWindow or HistogramWindow or CurveWindows
      //  delete LensfunData;    // TODO BJ: implement lensfun DB
      delete TheProcessor;
      delete GuiOptions;
      delete MainWindow;  // Cleans up HistogramWindow and ViewWindow also !
      ViewWindow = NULL;  // needs to be NULL to properly construct MainWindow
      delete PreviewImage;
      delete TheDcRaw;
}
void   ExtFileOpen(const QString file){
    ImageFileToOpen = file;
    CB_MenuFileOpen(1);
}

////////////////////////////////////////////////////////////////////////////////

void ptRemoveFile( const QString FileName) {
  if (QMessageBox::Yes == ptMessageBox::question(
                 MainWindow,
                 QObject::tr("Clean up input file"),
                 "As requested, Photivo will delete the input file " + FileName + ". Proceed?",
                 QMessageBox::No|QMessageBox::Yes)) {
    QFile::remove(FileName);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Hack : called at t=0 after the event loop started.
//
////////////////////////////////////////////////////////////////////////////////

void CB_Event0() {
  // Init Curves : supposed to be in event loop indeed.
  // (f.i. for progress reporting)
  QApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
  if (Settings->GetInt("JobMode") == 0)
    SaveButtonToolTip(Settings->GetInt("SaveButtonMode"));

  // Fill some look up tables
#pragma omp parallel for
  for (uint32_t i=0; i<0x10000; i++) {
    // uint16_t (0,0xffff) to float (0.0, 1.0)
    ToFloatTable[i]     = (float)i*ptInvWP;
    ToFloatABNeutral[i] = (float)i-ptWPHLab;
    ToInvertTable[i]    = ptWP - i;
    // linear RGB to sRGB table
    if (ToFloatTable[i] <= 0.0031308)
      ToSRGBTable[i] = CLIP((int32_t)(12.92*i));
    else
      ToSRGBTable[i] = CLIP((int32_t)((1.055*pow(ToFloatTable[i],1.0/2.4)-0.055)*ptWPf));
  }

  // Init run mode
  NextPhase = ptProcessorPhase_Raw;
  NextSubPhase = ptProcessorPhase_Load;
  ImageSaved = 0;

  PreCalcTransforms();

  if (Settings->GetInt("JobMode") == 0) { // not job mode!
    MainWindow->Settings_2_Form();

    // Load user settings
    if (Settings->GetInt("StartupSettings")) {
      if (ImageCleanUp == 0) {
        CB_OpenSettingsFile(Settings->GetString("StartupSettingsFile"));
      } else { // we got an image from gimp -> neutral display
        CB_OpenSettingsFile(Settings->GetString("PresetDirectory") + "/Neutral_absolute.pts");
      }
      // clean up
      QStringList Temp;
      Temp << "CropX" << "CropY" << "CropW" << "CropH";
      Temp << "RotateW" << "RotateH";
      for (int i = 0; i < Temp.size(); i++) Settings->SetValue(Temp.at(i),0);
    }

    if (PtsFileToOpen != "") {
      CB_OpenSettingsFile(PtsFileToOpen);
    }

    if (Sidecar != "") {
      ReadSidecar(Sidecar);
    }

    if (ImageFileToOpen != "") {
      CB_MenuFileOpen(1);
    }
    MainWindow->UpdateSettings();
    ViewWindow->setFocus();
    HistogramWindow->Init();
  }

  if (Settings->GetStringList("FavouriteTools").isEmpty()) {
    Settings->SetValue(
      "FavouriteTools",
      QStringList({
        "TabWhiteBalance",
        "TabRotation",
        "TabCrop",
        Fuid::ReinhardBrighten_RGB,
        Fuid::TextureContrast2_LabCC,
        Fuid::DetailCurve_LabSN,
        Fuid::Vignette_LabEyeCandy,
        Fuid::SigContrastRgb_Out
      })
    );
  }

#ifndef PT_WITHOUT_FILEMGR
  if (Settings->GetInt("FileMgrIsOpen")) {
    FileMgrWindow->displayThumbnails();
    FileMgrWindow->setFocus(Qt::OtherFocusReason);
  }
#endif

//prepare for further QFileOpenEvent(s)
#ifdef Q_OS_MAC
  TheApplication->macinit();
#endif

  InStartup = 0;
}

////////////////////////////////////////////////////////////////////////////////
//
// Precalculation of color transforms
//
////////////////////////////////////////////////////////////////////////////////

void PreCalcTransforms() {
  cmsToneCurve* Gamma = cmsBuildGamma(NULL, 1.0);
  cmsToneCurve* Gamma3[3];
  Gamma3[0] = Gamma3[1] = Gamma3[2] = Gamma;

  cmsHPROFILE InProfile = 0;

  cmsCIExyY DFromReference;

  switch (Settings->GetInt("WorkColor")) {
    case ptSpace_sRGB_D65 :
    case ptSpace_AdobeRGB_D65 :
      DFromReference = D65;
      break;
    case ptSpace_WideGamutRGB_D50 :
    case ptSpace_ProPhotoRGB_D50 :
      DFromReference = D50;
      break;
    default:
      assert(0);
  }

  ptImage::setCurrentRGB(Settings->GetInt("WorkColor"));

  InProfile = cmsCreateRGBProfile(&DFromReference,
                                  (cmsCIExyYTRIPLE*)&RGBPrimaries[Settings->GetInt("WorkColor")],
                                  Gamma3);

  if (!InProfile) {
    ptLogError(ptError_Profile,"Could not open InProfile profile.");
    return;
  }

  cmsFreeToneCurve(Gamma);

  if (Settings->GetInt("CMQuality") == ptCMQuality_HighResPreCalc) {
    ToPreviewTransform =
      cmsCreateTransform(InProfile,
                         TYPE_RGB_16,
                         PreviewColorProfile,
                         TYPE_RGB_16,
                         Settings->GetInt("PreviewColorProfileIntent"),
                         cmsFLAGS_HIGHRESPRECALC | cmsFLAGS_BLACKPOINTCOMPENSATION);
  } else { // fast sRGB preview also uses the not optimized profile for output
    ToPreviewTransform =
      cmsCreateTransform(InProfile,
                         TYPE_RGB_16,
                         PreviewColorProfile,
                         TYPE_RGB_16,
                         Settings->GetInt("PreviewColorProfileIntent"),
                         cmsFLAGS_NOOPTIMIZE | cmsFLAGS_BLACKPOINTCOMPENSATION);
  }

  cmsCloseProfile(InProfile);
}

////////////////////////////////////////////////////////////////////////////////
//
// Copy folder
//
////////////////////////////////////////////////////////////////////////////////

void copyFolder(QString sourceFolder, QString destFolder)
{
  QDir sourceDir(sourceFolder);
  if(!sourceDir.exists())
    return;
  QDir destDir(destFolder);
  if(!destDir.exists())
  {
    destDir.mkdir(destFolder);
  }
  QStringList files = sourceDir.entryList(QDir::Files);
  for(int i = 0; i< files.count(); i++)
  {
    QString srcName = sourceFolder + "/" + files[i];
    QString destName = destFolder + "/" + files[i];
    QFile::remove(destName);
    QFile::copy(srcName, destName);
  }
  files.clear();
  files = sourceDir.entryList(QDir::AllDirs | QDir::NoDotAndDotDot);
  for(int i = 0; i< files.count(); i++)
  {
    QString srcName = sourceFolder + "/" + files[i];
    QString destName = destFolder + "/" + files[i];
    copyFolder(srcName, destName);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Update
//
////////////////////////////////////////////////////////////////////////////////

void Update(short Phase,
            short SubPhase      /* = -1 */,
            short WithIdentify  /* = 1 */,
            short ProcessorMode /* = ptProcessorMode_Preview */)
{
#ifdef Q_OS_WIN
  ptEcWin7* Win7Taskbar = NULL;
  if (!JobMode) {
    Win7Taskbar = ptEcWin7::GetInstance();
    Win7Taskbar->setProgressState(ptEcWin7::Indeterminate);
  }
#endif

  if (Settings->GetInt("BlockUpdate") == 1) return; // hard block
  if (Settings->GetInt("PipeIsRunning") == 1) {
    // record that we got here with the new settings
    // and make sure the pipe is run after the current
    // pass automatically -> timer
    return;
  } else Settings->SetValue("PipeIsRunning",1);

  if (Phase < ptProcessorPhase_Preview) {
    // main processing
    if (Phase < NextPhase) NextPhase = Phase;
    if (SubPhase > 0 && SubPhase < NextSubPhase) NextSubPhase = SubPhase;
    if (Settings->GetInt("RunMode") == 1) {
      // we're in manual mode!
      MainWindow->UpdateSettings();
    } else {
      try {
        // run processor!
        ImageSaved = 0;
        MainWindow->UpdateSettings();
        if(Settings->GetInt("HaveImage")==1) {
          if (NextPhase < ptProcessorPhase_Output)
            TheProcessor->Run(NextPhase, NextSubPhase, WithIdentify, ProcessorMode);
          UpdatePreviewImage();
        }
        NextPhase = ptProcessorPhase_Output;
        NextSubPhase = ptProcessorPhase_Highlights;
      } catch (std::bad_alloc) {
        // make sure we image gets processed again
        NextPhase = ptProcessorPhase_Raw;
        NextSubPhase = ptProcessorPhase_Load;
        ViewWindow->ShowStatus(ptStatus_Done);
        ReportProgress(QObject::tr("Ready"));
        if (Settings->GetInt("PipeSize") == 0) {
          ptMessageBox::critical(NULL, "Memory error", "Processing aborted, memory error.\nReverting to 1:2 pipe. Reprocess with F5.");

          Settings->SetValue("PipeSize",1);
        } else {
          ptMessageBox::critical(NULL, "Memory error", "Processing aborted, memory error.");
        }
      }
    }
  } else if (Phase == ptProcessorPhase_OnlyHistogram) {
    // only histogram update, don't care about manual mode
    UpdatePreviewImage(NULL,1);
  } else if (Phase == ptProcessorPhase_Preview) {
    // only preview update, don't care about manual mode
    UpdatePreviewImage(NULL,0);
  } else if (Phase == ptProcessorPhase_WriteOut) {
    // write output
    WriteOut();
  } else if (Phase == ptProcessorPhase_ToGimp) {
    // export
    Export(ptExportMode_GimpPipe);
  } else {
    // should not happen!
    assert(0);
  }
  Settings->SetValue("PipeIsRunning",0);

#ifdef Q_OS_WIN
  if (!JobMode) {
    Win7Taskbar->setProgressState(ptEcWin7::NoProgress);
  }
#endif
}

////////////////////////////////////////////////////////////////////////////////
//
// Update
// Overloaded function, expects toolname and will call Update
//
////////////////////////////////////////////////////////////////////////////////

int GetProcessorPhase(const QString GuiName) {
  int Phase = 0;
  QString Tab = MainWindow->ProcessingTabBook->widget(
                  MainWindow->m_GroupBox->value(GuiName)->parentTabIdx() )->objectName();
  if      (Tab == "LocalTab")       Phase = ptProcessorPhase_LocalEdit;
  else if (Tab == "GeometryTab")    Phase = ptProcessorPhase_Geometry;
  else if (Tab == "RGBTab")         Phase = ptProcessorPhase_RGB;
  else if (Tab == "LabCCTab")       Phase = ptProcessorPhase_LabCC;
  else if (Tab == "LabSNTab")       Phase = ptProcessorPhase_LabSN;
  else if (Tab == "LabEyeCandyTab") Phase = ptProcessorPhase_LabEyeCandy;
  else if (Tab == "EyeCandyTab")    Phase = ptProcessorPhase_EyeCandy;
  else if (Tab == "OutTab")         Phase = ptProcessorPhase_Output;
  else                              Phase = ptProcessorPhase_Raw;
  return Phase;
}

//==============================================================================

void Update(const QString GuiName) {
  int Phase = GetProcessorPhase(GuiName);
  // It is assumed that no tool before white balance will use this.
  if (Phase < 2) {
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
  } else {
    Update(Phase);
  }
}


////////////////////////////////////////////////////////////////////////////////
//
// Block tools
//
////////////////////////////////////////////////////////////////////////////////

void BlockTools(const ptBlockToolsMode ANewState, QStringList AExcludeIds = QStringList()) {
  // Set the object name of the widget that is excluded from blocking.
  switch (ANewState) {
    case btmBlockForCrop:
      AExcludeIds << "TabCrop";
      break;
    case btmBlockForSpotRepair:
      AExcludeIds << "TabSpotRepair";
      break;
    default:
      // nothing to do
      break;
  }

  bool EnabledStatus = ANewState == btmUnblock;

  // Handle all necessary widgets outside the processing tabbook
  MainWindow->HistogramFrameCentralWidget->setEnabled(EnabledStatus);
  MainWindow->SearchWidget->setEnabled(EnabledStatus);
  MainWindow->PipeControlWidget->setEnabled(EnabledStatus);
  MainWindow->StatusWidget->setEnabled(EnabledStatus);

  if (MainWindow->m_MovedTools->size() > 0) {
    /* Process moved tools list when UI is not in tab mode (e.g. showing favourites)
      We just cycle through the list of currently visible tools and en/disable them.
    */
    for (QWidget *hToolBox: *MainWindow->m_MovedTools) {
      if (!AExcludeIds.contains(hToolBox->objectName())) {
        if (hToolBox->objectName().contains("-"))  //new-style
          hToolBox->setEnabled(EnabledStatus);
        else                                       // old-style
          ((ptGroupBox*)hToolBox)->SetEnabled(EnabledStatus);
      }
    }

  } else {
    /* Process processing tabbook
      To avoid cycling through all ptGroupBox objects on every tab every time we use the following
      approach: We assume that the tabbook is switched to the appropriate tab when BlockTools() is
      called, i.e. the tab containing the filter that should not be blocked (if there is such a one).
      We en/disable all tabs except the current one completely. Now we only need to cycle through the
      remaining ptGroupBoxes on the current tab.
    */
    int CurrentTab = MainWindow->ProcessingTabBook->currentIndex();

    for (int i = 0; i < MainWindow->ProcessingTabBook->count(); i++) {
      if (i != CurrentTab)
        MainWindow->ProcessingTabBook->setTabEnabled(i, EnabledStatus);
    }

    QList<ptToolBox*> NewToolList =
        MainWindow->ProcessingTabBook->widget(CurrentTab)->findChildren<ptToolBox*>();
    foreach (ptToolBox* Tool, NewToolList) {
      if (!AExcludeIds.contains(Tool->objectName()))
        Tool->setEnabled(EnabledStatus);
    }
    QList<ptGroupBox*> ToolList =
        MainWindow->ProcessingTabBook->widget(CurrentTab)->findChildren<ptGroupBox*>();
    foreach (ptGroupBox* Tool, ToolList) {
      if (!AExcludeIds.contains(Tool->objectName()))
        Tool->SetEnabled(EnabledStatus);
    }
  }

  Settings->SetValue("BlockTools", ANewState);
}

////////////////////////////////////////////////////////////////////////////////
//
// Histogram
//
////////////////////////////////////////////////////////////////////////////////
void HistogramCropDone(const ptStatus ExitStatus, QRect SelectionRect);

void HistogramGetCrop() {
  // Get the crop for the histogram
  if (Settings->GetInt("HistogramCrop")) {
    // Allow to be selected in the view window. And deactivate main.
    BlockTools(btmBlockAll);
    ViewWindow->ShowStatus(QObject::tr("Selection"));
    ViewWindow->StartSimpleRect(HistogramCropDone);
    ViewWindow->setFocus();
  } else {
    ReportProgress(QObject::tr("Updating histogram"));
    Update(ptProcessorPhase_Preview);
    ReportProgress(QObject::tr("Ready"));
  }
}

void HistogramCropDone(const ptStatus ExitStatus, QRect SelectionRect) {
  // Selection is done at this point. Disallow it further and activate main.
  BlockTools(btmUnblock);
  if (ExitStatus == stFailure) {
    Settings->SetValue("HistogramCropX",0);
    Settings->SetValue("HistogramCropY",0);
    Settings->SetValue("HistogramCropW",0);
    Settings->SetValue("HistogramCropH",0);
    Settings->SetValue("HistogramCrop",0);
    return;
  }

  short XScale = 1<<Settings->GetInt("PipeSize");
  short YScale = 1<<Settings->GetInt("PipeSize");

  Settings->SetValue("HistogramCropX", SelectionRect.left() * XScale);
  Settings->SetValue("HistogramCropY", SelectionRect.top() * YScale);
  Settings->SetValue("HistogramCropW", SelectionRect.width() * XScale);
  Settings->SetValue("HistogramCropH", SelectionRect.height() * YScale);

  // Check if the chosen area is large enough
  if (Settings->GetInt("HistogramCropW") < 50 || Settings->GetInt("HistogramCropH") < 50) {
    ptMessageBox::information(0,
      QObject::tr("Selection too small"),
      QObject::tr("Selection rectangle needs to be at least 50x50 pixels in size.\nNo crop, try again."));
    Settings->SetValue("HistogramCropX",0);
    Settings->SetValue("HistogramCropY",0);
    Settings->SetValue("HistogramCropW",0);
    Settings->SetValue("HistogramCropH",0);
    Settings->SetValue("HistogramCrop",0);
  }

  ReportProgress(QObject::tr("Updating histogram"));
  Update(ptProcessorPhase_Preview);
  ReportProgress(QObject::tr("Ready"));
}


////////////////////////////////////////////////////////////////////////////////
//
// Status report in Viewwindow
//
////////////////////////////////////////////////////////////////////////////////

void ViewWindowShowStatus(short State) {
  if (ViewWindow) {
    ViewWindow->ShowStatus(State);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Operations done right before gamma and profile for whatever reason
// Finalrun to allow filters only on the final run
// Resize to forbid resizing (-> histogram only on crop)
//
////////////////////////////////////////////////////////////////////////////////

void BeforeGamma(ptImage* Image, const short FinalRun = 0, const short Resize = 1) {

  if (Settings->GetInt("WebResizeBeforeGamma")==1 && Resize) {
    if (FinalRun == 1) Settings->SetValue("FullOutput",1);
    if (Settings->ToolIsActive("TabWebResize")) {
      ReportProgress(QObject::tr("WebResizing"));
      Image->ptGMResize(Settings->GetInt("WebResizeScale"),
                        0,
                        Settings->GetInt("WebResizeFilter"),
                        Settings->GetInt("WebResizeDimension"));
    }
    if (FinalRun == 1) Settings->SetValue("FullOutput",0);
  }

  // TODO put these curves together as devicelink into lcms

  // BaseCurve.
  auto hFilter = GFilterDM->GetFilterFromName(Fuid::RgbCurve_Out);
  if (hFilter->isActive()) {
    ReportProgress(QObject::tr("Applying base curve"));
    hFilter->runFilter(Image);
  }

  //GammaCompensation
  if (Settings->ToolIsActive("TabGammaCompensation")) {
    ReportProgress(QObject::tr("Applying gamma compensation"));
    ptCurve* CompensationCurve = new ptCurve();
    CompensationCurve->setFromFunc(ptCurve::DeltaGammaTool,Settings->GetDouble("OutputGamma"),
              Settings->GetDouble("OutputLinearity"));
    Image->ApplyCurve(CompensationCurve,7);
    delete CompensationCurve;
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Operations done after gamma and profile for whatever reason
// Finalrun to allow filters only on the final run
// Resize to forbid resizing (-> histogram only on crop)
//
////////////////////////////////////////////////////////////////////////////////

void AfterAll(ptImage* Image, const short FinalRun = 0, const short Resize = 1) {
  // Sigmoidal contrast
  auto hFilter = GFilterDM->GetFilterFromName(Fuid::SigContrastRgb_Out);
  if (hFilter->isActive()) {
    ReportProgress(QObject::tr("Applying RGB Contrast"));
    hFilter->runFilter(Image);
  }

  // After gamma curve
  hFilter = GFilterDM->GetFilterFromName(Fuid::AfterGammaCurve_Out);
  if (hFilter->isActive()) {
    ReportProgress(QObject::tr("Applying after gamma curve"));
    hFilter->runFilter(Image);
  }

  // WebResize for quality reasons done after output profile
  if (Settings->GetInt("WebResizeBeforeGamma")==0 && Resize) {
    if (FinalRun == 1) Settings->SetValue("FullOutput",1);
    if (Settings->ToolIsActive("TabWebResize")) {
      ReportProgress(QObject::tr("WebResizing"));
      Image->ptGMResize(Settings->GetInt("WebResizeScale"),
                        0,
                        Settings->GetInt("WebResizeFilter"),
                        Settings->GetInt("WebResizeDimension"));
    }
    if (FinalRun == 1) Settings->SetValue("FullOutput",0);
  }
}

void EndSharpen(ptImage* Image, cmsHPROFILE Profile, const int Intent) {
  ReportProgress(QObject::tr("Wiener Filter"));
  int InColorSpace = Image->m_ColorSpace;
  cmsHPROFILE InternalProfile = 0;
  if (Profile == NULL || InColorSpace != ptSpace_Profiled) {
    // linear case
    cmsToneCurve* Gamma = cmsBuildGamma(NULL, 1.0);
    cmsToneCurve* Gamma3[3];
    Gamma3[0] = Gamma3[1] = Gamma3[2] = Gamma;

    cmsCIExyY       DFromReference;

    switch (Image->m_ColorSpace) {
      case ptSpace_sRGB_D65 :
      case ptSpace_AdobeRGB_D65 :
        DFromReference = D65;
        break;
      case ptSpace_WideGamutRGB_D50 :
      case ptSpace_ProPhotoRGB_D50 :
        DFromReference = D50;
        break;
      default:
        assert(0);
    }

    InternalProfile = cmsCreateRGBProfile(&DFromReference,
                                    (cmsCIExyYTRIPLE*)&RGBPrimaries[Image->m_ColorSpace],
                                    Gamma3);

    if (!InternalProfile) {
      ptLogError(ptError_Profile,"Could not open InternalProfile profile.");
      return;
    }

    cmsFreeToneCurve(Gamma);
  } else {
    // profiled case
    InternalProfile = Profile;
  }

  cmsHPROFILE LabProfile = 0;
  LabProfile = cmsCreateLab4Profile(NULL);
  // to Lab
  cmsHTRANSFORM Transform;
  Transform = cmsCreateTransform(InternalProfile,
                                 TYPE_RGB_16,
                                 LabProfile,
                                 TYPE_Lab_16,
                                 Intent,
                                 cmsFLAGS_NOOPTIMIZE | cmsFLAGS_BLACKPOINTCOMPENSATION);

  int32_t Size = Image->m_Width*Image->m_Height;
  int32_t Step = 100000;
#pragma omp parallel for schedule(static)
  for (int32_t i = 0; i < Size; i+=Step) {
    int32_t Length = (i+Step)<Size ? Step : Size - i;
    uint16_t* Tile = &(Image->m_Image[i][0]);
    cmsDoTransform(Transform,Tile,Tile,Length);
  }
  Image->m_ColorSpace = ptSpace_Lab;
  // Wiener Filter
  GFilterDM->GetFilterFromName(Fuid::Wiener_Out)->runFilter(Image);

  // to RGB
  Transform = cmsCreateTransform(LabProfile,
                                 TYPE_Lab_16,
                                 InternalProfile,
                                 TYPE_RGB_16,
                                 Intent,
                                 cmsFLAGS_NOOPTIMIZE | cmsFLAGS_BLACKPOINTCOMPENSATION);

#pragma omp parallel for schedule(static)
  for (int32_t i = 0; i < Size; i+=Step) {
    int32_t Length = (i+Step)<Size ? Step : Size - i;
    uint16_t* Tile = &(Image->m_Image[i][0]);
    cmsDoTransform(Transform,Tile,Tile,Length);
  }

  cmsDeleteTransform(Transform);
  cmsCloseProfile(LabProfile);
  if (Profile == NULL || InColorSpace != ptSpace_Profiled)
    cmsCloseProfile(InternalProfile);
  Image->m_ColorSpace = InColorSpace;
}

////////////////////////////////////////////////////////////////////////////////
//
// Determine and update the preview image.
//
// ForcedImage claims that this image should be previewed.
// (somehow of a hack : it is used for the Crop function where
// we want the user to present fast a  new image to be cropped without
// having to calculate through the whole pipe).
//
// OnlyHistogram : hack to sync up the histogram when it was inactive.
//
////////////////////////////////////////////////////////////////////////////////

void UpdatePreviewImage(const ptImage* ForcedImage   /* = NULL  */,
                        const short    OnlyHistogram /* = false */) {

  if (!PreviewImage) PreviewImage = new (ptImage);

  // If we don't have yet a m_Image_AfterGeometry we are in
  // startup condition and show the start screen
  if (!TheProcessor->m_Image_AfterGeometry) {
    MainWindow->ViewFrameStackedWidget->setCurrentWidget(MainWindow->ViewStartPage);
    return;
  }

  MainWindow->ViewFrameStackedWidget->setCurrentWidget(MainWindow->ViewFrameCentralWidget);

  // Fast display of an image, e.g. for cropping
  // -> no histogram
  // -> only exposure for better display
  // -> transfer to preview color space
  if (ForcedImage) {
    PreviewImage->Set(ForcedImage);
    if (Settings->GetDouble("CropExposure") != 0.0) {
      PreviewImage->Expose(pow(2,Settings->GetDouble("CropExposure")), TExposureClipMode::Ratio);
    }
    BeforeGamma(PreviewImage,0,0);

    // Convert from working space to screen space.
    // Using lcms and a standard sRGB or custom profile.
    ptImage* ReturnValue = PreviewImage->lcmsRGBToPreviewRGB(Settings->GetInt("CMQuality") == ptCMQuality_FastSRGB);
    if (!ReturnValue) {
      ptLogError(ptError_lcms,"lcmsRGBToPreviewRGB");
      assert(ReturnValue);
    }
    AfterAll(PreviewImage,0,0);

    ViewWindow->UpdateImage(PreviewImage);
    ReportProgress(QObject::tr("Ready"));
    return;
  }

  ViewWindow->ShowStatus(ptStatus_Updating);
  ReportProgress(QObject::tr("Updating preview image"));

  if (!HistogramImage) HistogramImage = new (ptImage);

  short ActiveTab = MainWindow->GetCurrentTab();

  Settings->SetValue("ShowExposureIndicatorSensor",0);

  // Determine first what is the current image.
  if (ForcedImage) {
    PreviewImage->Set(ForcedImage);
  } else if (Settings->GetInt("PreviewMode") == ptPreviewMode_End) {
    PreviewImage->Set(TheProcessor->m_Image_AfterEyeCandy);
  } else {
    if (!Settings->useRAWHandling()) ActiveTab = MAX(ptLocalTab, ActiveTab);
    switch (ActiveTab) {
      case ptCameraTab:
        Settings->SetValue("ShowExposureIndicatorSensor",1);
        if (Settings->GetInt("ExposureIndicatorSensor")) {
          PreviewImage->Set(TheDcRaw,
                             Settings->GetInt("WorkColor"));
        } else {
          PreviewImage->Set(TheProcessor->m_Image_AfterDcRaw);
        }
        break;
      case ptLocalTab:
        PreviewImage->Set(TheProcessor->m_Image_AfterLocalEdit);
        break;
      case ptGeometryTab:
        PreviewImage->Set(TheProcessor->m_Image_AfterGeometry);
        break;
      case ptRGBTab:
        PreviewImage->Set(TheProcessor->m_Image_AfterRGB);
        break;
      case ptLabCCTab:
        PreviewImage->Set(TheProcessor->m_Image_AfterLabCC);
        break;
      case ptLabSNTab:
        PreviewImage->Set(TheProcessor->m_Image_AfterLabSN);
        break;
      case ptLabEyeCandyTab:
        PreviewImage->Set(TheProcessor->m_Image_AfterLabEyeCandy);
        break;
      case ptEyeCandyTab:
      case ptOutTab:
        PreviewImage->Set(TheProcessor->m_Image_AfterEyeCandy);
        break;
      default:
        // Should not happen.
        assert(0);
    }
  }

  if (Settings->GetInt("HistogramMode")==ptHistogramMode_Linear) {
    HistogramImage->Set(PreviewImage);
    // If we are in a Tab preview mode and in the lab mode
    // we do the conversion to Lab in case it would not have been
    // done yet (due to suppressed for speed in absense of USM or L Curve).
    // This way the histogram is an L Histogram at this point.
    if ( (Settings->GetInt("PreviewMode") == ptPreviewMode_Tab) &&
         (ActiveTab == ptLabCCTab || ActiveTab == ptLabSNTab || ActiveTab == ptLabEyeCandyTab) &&
         (HistogramImage->m_ColorSpace != ptSpace_Lab) ) {
      HistogramImage->RGBToLab();
    }
  } else if (Settings->GetInt("HistogramMode")==ptHistogramMode_Output &&
             !(Settings->GetInt("HistogramCrop") && !Settings->GetInt("WebResize"))) {
    HistogramImage->Set(PreviewImage);


  } else if (Settings->GetInt("HistogramCrop")) {
    HistogramImage->Set(PreviewImage);
  }

  if (PreviewImage->m_ColorSpace == ptSpace_Lab)
    PreviewImage->LabToRGB(Settings->GetInt("WorkColor"));

  uint16_t Width = 0;
  uint16_t Height = 0;
  uint16_t TempCropX = 0;
  uint16_t TempCropY = 0;
  uint16_t TempCropW = 0;
  uint16_t TempCropH = 0;

  if (Settings->GetInt("HistogramMode")==ptHistogramMode_Linear) {
    ReportProgress(QObject::tr("Updating histogram"));

    if (Settings->GetInt("HistogramCrop")) {
      float TmpScaled = TheProcessor->m_ScaleFactor;
      Width = HistogramImage->m_Width;
      Height = HistogramImage->m_Height;
      TempCropX = Settings->GetInt("HistogramCropX")*TmpScaled;
      TempCropY = Settings->GetInt("HistogramCropY")*TmpScaled;
      TempCropW = Settings->GetInt("HistogramCropW")*TmpScaled;
      TempCropH = Settings->GetInt("HistogramCropH")*TmpScaled;
      if ((((TempCropX) + (TempCropW)) > Width) ||
          (((TempCropY) + (TempCropH)) >  Height)) {
        ptMessageBox::information(MainWindow,
               QObject::tr("Histogram selection outside the image"),
               QObject::tr("Histogram selection rectangle too large.\nNo crop, try again."));
        Settings->SetValue("HistogramCropX",0);
        Settings->SetValue("HistogramCropY",0);
        Settings->SetValue("HistogramCropW",0);
        Settings->SetValue("HistogramCropH",0);
        Settings->SetValue("HistogramCrop",0);
      } else {
        HistogramImage->Crop(TempCropX, TempCropY, TempCropW, TempCropH);
      }
    }

    // In case of histogram update only, we're done.
    if (OnlyHistogram) {
      HistogramWindow->UpdateView(HistogramImage);
      ViewWindow->ShowStatus(ptStatus_Done);
      return;
    }
  }

  uint8_t ExposureChannelMask = 0;
  uint16_t OverExposureLevel[3];
  uint16_t UnderExposureLevel[3];
  if (Settings->GetInt("ExposureIndicator")) {
    ReportProgress(QObject::tr("Indicating exposure"));

    if (Settings->GetInt("ExposureIndicatorR")) ExposureChannelMask |= 1;
    if (Settings->GetInt("ExposureIndicatorG")) ExposureChannelMask |= 2;
    if (Settings->GetInt("ExposureIndicatorB")) ExposureChannelMask |= 4;

    if (ActiveTab == ptCameraTab && Settings->GetInt("ExposureIndicatorSensor")){
      OverExposureLevel[0] = CLIP((int32_t)
        ((TheDcRaw->m_WhiteLevel_AfterPhase1-TheDcRaw->m_BlackLevel_AfterPhase1)
        *TheDcRaw->m_Multipliers[0]));
      OverExposureLevel[1] = CLIP((int32_t)
        ((TheDcRaw->m_WhiteLevel_AfterPhase1-TheDcRaw->m_BlackLevel_AfterPhase1)
        *TheDcRaw->m_Multipliers[1]));
      OverExposureLevel[2] = CLIP((int32_t)
        ((TheDcRaw->m_WhiteLevel_AfterPhase1-TheDcRaw->m_BlackLevel_AfterPhase1)
        *TheDcRaw->m_Multipliers[2]));
      UnderExposureLevel[0] = 0x0000;
      UnderExposureLevel[1] = 0x0000;
      UnderExposureLevel[2] = 0x0000;
    } else {
      OverExposureLevel[0] = 0xfff0;
      OverExposureLevel[1] = 0xfff0;
      OverExposureLevel[2] = 0xfff0;
      UnderExposureLevel[0] = 0x000f;
      UnderExposureLevel[1] = 0x000f;
      UnderExposureLevel[2] = 0x000f;
    }
  }

  if (!ForcedImage) {
    BeforeGamma(PreviewImage);
  }

  if (Settings->GetInt("HistogramMode")==ptHistogramMode_Output) {

    ReportProgress(QObject::tr("Updating histogram"));
    if (Settings->GetInt("HistogramCrop") && !Settings->GetInt("WebResize")) {
      HistogramImage->Set(PreviewImage);
    } else {
      BeforeGamma(HistogramImage,0,0);
    }

    ReportProgress(QObject::tr("Converting to output space"));

    cmsHPROFILE OutputColorProfile = NULL;
    OutputColorProfile = cmsOpenProfileFromFile(
                         Settings->GetString("OutputColorProfile").toLocal8Bit().data(), "r");
    if (!OutputColorProfile) {
      ptLogError(ptError_FileOpen,
        Settings->GetString("OutputColorProfile").toLocal8Bit().data());
      assert(OutputColorProfile);
    }

    // Convert from working space to screen space.
    // Using lcms and a standard sRGB or custom profile.

    ptImage* ReturnValue = HistogramImage->lcmsRGBToRGB(
      OutputColorProfile,
      Settings->GetInt("OutputColorProfileIntent"),
      Settings->GetInt("CMQuality"));
    if (!ReturnValue) {
      ptLogError(ptError_lcms,"lcmsRGBToRGB(OutputColorProfile)");
      assert(ReturnValue);
    }

    if (Settings->GetInt("HistogramCrop") && !Settings->GetInt("WebResize")) {
      AfterAll(HistogramImage);
      if (GFilterDM->GetFilterFromName(Fuid::Wiener_Out)->isActive() &&
          Settings->GetInt("PreviewMode") == ptPreviewMode_End)
      {
        EndSharpen(HistogramImage, OutputColorProfile, Settings->GetInt("OutputColorProfileIntent"));
      }
    } else {
      AfterAll(HistogramImage,0,0);
      if (GFilterDM->GetFilterFromName(Fuid::Wiener_Out)->isActive() &&
          Settings->GetInt("PreviewMode") == ptPreviewMode_End)
      {
        EndSharpen(HistogramImage, OutputColorProfile, Settings->GetInt("OutputColorProfileIntent"));
      }
    }

    // Close the output profile.
    cmsCloseProfile(OutputColorProfile);

    if (Settings->GetInt("HistogramCrop")) {
      float TmpScaled = TheProcessor->m_ScaleFactor;
      Width = HistogramImage->m_Width;
      Height = HistogramImage->m_Height;
      TempCropX = Settings->GetInt("HistogramCropX")*TmpScaled;
      TempCropY = Settings->GetInt("HistogramCropY")*TmpScaled;
      TempCropW = Settings->GetInt("HistogramCropW")*TmpScaled;
      TempCropH = Settings->GetInt("HistogramCropH")*TmpScaled;
      if ((((TempCropX) + (TempCropW)) >  Width) ||
          (((TempCropY) + (TempCropH)) >  Height)) {
        ptMessageBox::information(MainWindow,
          QObject::tr("Histogram selection outside the image"),
          QObject::tr("Histogram selection rectangle too large.\nNo crop, try again."));
        Settings->SetValue("HistogramCropX",0);
        Settings->SetValue("HistogramCropY",0);
        Settings->SetValue("HistogramCropW",0);
        Settings->SetValue("HistogramCropH",0);
        Settings->SetValue("HistogramCrop",0);
      } else {
        HistogramImage->Crop(TempCropX, TempCropY, TempCropW, TempCropH);
      }
    }
    if (OnlyHistogram) {
      HistogramWindow->UpdateView(HistogramImage);
      ViewWindow->ShowStatus(ptStatus_Done);
      return;
    }
  }

  ReportProgress(QObject::tr("Converting to screen space"));

  // Convert from working space to screen space.
  // Using lcms and a standard sRGB or custom profile.
  PreviewImage->lcmsRGBToPreviewRGB(Settings->GetInt("CMQuality") == ptCMQuality_FastSRGB);

  if (!ForcedImage) {
    AfterAll(PreviewImage);
    if (GFilterDM->GetFilterFromName(Fuid::Wiener_Out)->isActive() &&
        Settings->GetInt("PreviewMode") == ptPreviewMode_End)
    {
      EndSharpen(PreviewImage, PreviewColorProfile, Settings->GetInt("PreviewColorProfileIntent"));
    }
  }

  if (Settings->GetInt("HistogramMode")==ptHistogramMode_Preview) {
    if (!Settings->GetInt("HistogramCrop") ||
        (Settings->GetInt("HistogramCrop") && !Settings->GetInt("WebResize"))) {
      HistogramImage->Set(PreviewImage);
    } else if (Settings->GetInt("HistogramCrop")) {
      BeforeGamma(HistogramImage,0,0);

      ptImage* ReturnValue = HistogramImage->lcmsRGBToPreviewRGB(Settings->GetInt("CMQuality") == ptCMQuality_FastSRGB);
      if (!ReturnValue) {
        ptLogError(ptError_lcms,"lcmsRGBToPreviewRGB");
        assert(ReturnValue);
      }

      AfterAll(HistogramImage,0,0);
      if (GFilterDM->GetFilterFromName(Fuid::Wiener_Out)->isActive() &&
          Settings->GetInt("PreviewMode") == ptPreviewMode_End)
      {
        EndSharpen(HistogramImage, PreviewColorProfile, Settings->GetInt("PreviewColorProfileIntent"));
      }
    }

    if (Settings->GetInt("HistogramCrop")) {
      float TmpScaled = TheProcessor->m_ScaleFactor;
      Width = HistogramImage->m_Width;
      Height = HistogramImage->m_Height;
      TempCropX = Settings->GetInt("HistogramCropX")*TmpScaled;
      TempCropY = Settings->GetInt("HistogramCropY")*TmpScaled;
      TempCropW = Settings->GetInt("HistogramCropW")*TmpScaled;
      TempCropH = Settings->GetInt("HistogramCropH")*TmpScaled;
      if ((((TempCropX) + (TempCropW)) >  Width) ||
          (((TempCropY) + (TempCropH)) >  Height)) {
        ptMessageBox::information(MainWindow,
          QObject::tr("Histogram selection outside the image"),
          QObject::tr("Histogram selection rectangle too large.\nNo crop, try again."));
        Settings->SetValue("HistogramCropX",0);
        Settings->SetValue("HistogramCropY",0);
        Settings->SetValue("HistogramCropW",0);
        Settings->SetValue("HistogramCropH",0);
        Settings->SetValue("HistogramCrop",0);
      } else {
        HistogramImage->Crop(TempCropX, TempCropY, TempCropW, TempCropH);
      }
    }
  }

  if (!OnlyHistogram) {
    if (Settings->GetInt("HistogramMode")==ptHistogramMode_Output) {
      if (Settings->GetInt("ExposureIndicator"))
        PreviewImage->IndicateExposure(HistogramImage,
                                      Settings->GetInt("ExposureIndicatorOver"),
                                      Settings->GetInt("ExposureIndicatorUnder"),
                                      ExposureChannelMask,
                                      OverExposureLevel,
                                      UnderExposureLevel);
    } else {
      if (Settings->GetInt("ExposureIndicator"))
        PreviewImage->IndicateExposure(Settings->GetInt("ExposureIndicatorOver"),
                                      Settings->GetInt("ExposureIndicatorUnder"),
                                      ExposureChannelMask,
                                      OverExposureLevel,
                                      UnderExposureLevel);
    }

    // Get a special preview
    if (Settings->GetInt("SpecialPreview"))
      PreviewImage->SpecialPreview(Settings->GetInt("SpecialPreview"),
                                   Settings->GetInt("PreviewColorProfileIntent"));

    // Update the ViewWindow and show if needed.
    ViewWindow->UpdateImage(PreviewImage);
  }

  HistogramWindow->UpdateView(HistogramImage);
  ViewWindow->ShowStatus(ptStatus_Done);

  ReportProgress(QObject::tr("Ready"));

  if (!OnlyHistogram) {
    if (Settings->GetInt("WriteBackupSettings")) {
      GFilterDM->WritePresetFile(Settings->GetString("UserDirectory")+"backup.pts");
    }
  }
} // UpdatePreviewImage


////////////////////////////////////////////////////////////////////////////////
//
// RunJob (the non-interactive 'batch' run).
//
////////////////////////////////////////////////////////////////////////////////

void RunJob(const QString JobFileName) {
  Settings->SetValue("JobMode",1);
  // Init
  CB_Event0();

  // Read the gui settings from a file.
  short NextPhase = 1;
  short ReturnValue = GFilterDM->ReadPresetFile(JobFileName, NextPhase);
  // TODO: BJ Implement spot processing
  if (ReturnValue) {
    printf("\nNo valid job file!\n\n");
    return;
  }

  QStringList InputFileNameList = Settings->GetStringList("InputFileNameList");
  assert(InputFileNameList.size());

  do {
    try {
      // Test if we can handle the file
      ptDcRaw* TestDcRaw = new(ptDcRaw);
      Settings->ToDcRaw(TestDcRaw);
      uint16_t InputWidth = 0;
      uint16_t InputHeight = 0;
      ptImageType InputType = CheckImageType(InputFileNameList[0],
                                             &InputWidth, &InputHeight,
                                             TestDcRaw);

      if (InputType <= itUndetermined) {
        // not a supported image type or error reading the file
        QString ErrorMessage = QObject::tr("Cannot decode")
                             + " '"
                             + InputFileNameList[0]
                             + "'" ;
        printf("%s\n",ErrorMessage.toLocal8Bit().data());
        delete TestDcRaw;


      } else {
        // process the image
        if (InputType == itRaw) {
          Settings->SetValue("IsRAW",1);
        } else {
          Settings->SetValue("IsRAW",0);
        }

        if (!Settings->useRAWHandling()) {
          Settings->SetValue("ExposureNormalization", 0.0);
        }

        QFileInfo PathInfo(InputFileNameList[0]);
        if (!Settings->GetString("OutputDirectory").isEmpty()) {
          Settings->SetValue("OutputFileName",
            Settings->GetString("OutputDirectory") + "/" + PathInfo.completeBaseName());
        } else {
          Settings->SetValue("OutputFileName",
            PathInfo.dir().path() + "/" + PathInfo.completeBaseName());
        }
        if (!Settings->GetString("OutputFileNameSuffix").isEmpty()) {
          Settings->SetValue("OutputFileName",
                             Settings->GetString("OutputFileName") +
                             Settings->GetString("OutputFileNameSuffix"));
        }
        else {
          if (!Settings->GetInt("IsRAW")) {
            Settings->SetValue("OutputFileName",
                               Settings->GetString("OutputFileName") + "-new");
          }
        }

        // Here we have the OutputFileName, but extension still to add.
        switch(Settings->GetInt("SaveFormat")) {
          case ptSaveFormat_JPEG :
            Settings->SetValue("OutputFileName",
                               Settings->GetString("OutputFileName") + ".jpg");
            break;
          case ptSaveFormat_PNG :
          case ptSaveFormat_PNG16 :
            Settings->SetValue("OutputFileName",
                               Settings->GetString("OutputFileName") + ".png");
            break;
          case ptSaveFormat_TIFF8 :
          case ptSaveFormat_TIFF16 :
            Settings->SetValue("OutputFileName",
                               Settings->GetString("OutputFileName") + ".tif");
            break;
          default :
            Settings->SetValue("OutputFileName",
                               Settings->GetString("OutputFileName") + ".ppm");
            break;
        }

        // Processing the job.
        delete TheDcRaw;
        delete TheProcessor;
        TheProcessor = new ptProcessor(ReportProgress);
        TheDcRaw = TestDcRaw;
        Settings->SetValue("ImageW",InputWidth);
        Settings->SetValue("ImageH",InputHeight);

        TheProcessor->m_DcRaw = TheDcRaw;

        Settings->SetValue("RunMode",0);
        Settings->SetValue("FullOutput",1);
        TheProcessor->Run(ptProcessorPhase_Raw,ptProcessorPhase_Load,0,1);
        Settings->SetValue("FullOutput",0);
        // And write result.
        Update(ptProcessorPhase_WriteOut);
      }
    } catch (std::bad_alloc) {
      QString ErrorMessage = QObject::tr("Memory error, no conversion for file:")
                           + " '"
                           + InputFileNameList[0]
                           + "'" ;
      printf("%s\n",ErrorMessage.toLocal8Bit().data());
    }

    // Loop over the inputfiles by shifting the next one to index 0
    if (InputFileNameList.size()) {
      InputFileNameList.removeAt(0);
    }

    // Write the changed InputFileNameList back to the settings.
    Settings->SetValue("InputFileNameList",InputFileNameList);

  } while (InputFileNameList.size());
} // RunJob

//==============================================================================
// Prepares the tags list and sets TagsList and DigikamTagsList
void PrepareTags(const QString TagsInput) {

  QString WorkString = TagsInput;
  WorkString.replace("\n",",");
  while (WorkString.contains("  "))
    WorkString.replace("  "," ");
  if (WorkString.startsWith(" ")) WorkString.remove(0,1);
  if (WorkString.endsWith(" ")) WorkString.remove(WorkString.size()-1,1);
  while (WorkString.contains("//"))
    WorkString.replace("//","/");
  while (WorkString.contains("/ /"))
    WorkString.replace("/ /","/");
  if (WorkString.endsWith("/")) WorkString.remove(WorkString.size()-1,1);
  WorkString.replace(" ,",",");
  WorkString.replace("/,",",");
  WorkString.replace(" ,",","); // We need this again!
  WorkString.replace(" /","/");
  WorkString.replace("/ ","/");
  WorkString.replace(",/",",");
  WorkString.replace(", ",",");
  if (WorkString.startsWith("/")) WorkString.remove(0,1);

  QStringList DigikamTags = WorkString.split(",", QString::SkipEmptyParts);

  QStringList Tags;
  QString TempString;
  for (int i = 0; i < DigikamTags.size(); i++) {
    TempString = DigikamTags.at(i);
    TempString.remove(0,TempString.lastIndexOf("/")+1);
    if (TempString.startsWith(" ")) TempString.remove(0,1);
    Tags << TempString;
  }

  Settings->SetValue("DigikamTagsList", DigikamTags);
  Settings->SetValue("TagsList", Tags);
}

//==============================================================================

/*! Write out in one of the output formats (after applying output profile). */
void WriteOut() {
  ptImage* OutImage = NULL;

  if (Settings->GetInt("JobMode") == 1) {
    OutImage = TheProcessor->m_Image_AfterEyeCandy; // Job mode -> no cache
  } else {
    if (!OutImage) OutImage = new(ptImage);
    OutImage->Set(TheProcessor->m_Image_AfterEyeCandy);
  }

  BeforeGamma(OutImage, 1);

  cmsHPROFILE OutputColorProfile = NULL;

  ReportProgress(QObject::tr("Converting to output space"));

  // Prepare and open an output profile.
  OutputColorProfile = cmsOpenProfileFromFile(
    Settings->GetString("OutputColorProfile").toLocal8Bit().data(),
    "r");
  if (!OutputColorProfile) {
    ptLogError(ptError_FileOpen,
   Settings->GetString("OutputColorProfile").toLocal8Bit().data());
    assert(OutputColorProfile);
  }

  // Color space conversion
  ptImage* ReturnValue = OutImage->lcmsRGBToRGB(
    OutputColorProfile,
    Settings->GetInt("OutputColorProfileIntent"),
    Settings->GetInt("CMQuality"));
  if (!ReturnValue) {
    ptLogError(ptError_lcms,"lcmsRGBToRGB(OutputColorProfile)");
    assert(ReturnValue);
  }

  AfterAll(OutImage, 1);
  if (GFilterDM->GetFilterFromName(Fuid::Wiener_Out)->isActive())
    EndSharpen(OutImage, OutputColorProfile, Settings->GetInt("OutputColorProfileIntent"));
  // Close the output profile.
  cmsCloseProfile(OutputColorProfile);

  ReportProgress(QObject::tr("Writing output"));

  bool FileWritten =  OutImage->ptGMCWriteImage(
    Settings->GetString("OutputFileName").toLocal8Bit().data(),
    Settings->GetInt("SaveFormat"),
    Settings->GetInt("SaveQuality"),
    Settings->GetInt("SaveSampling"),
    Settings->GetInt("SaveResolution"),
    Settings->GetString("OutputColorProfile").toLocal8Bit().data(),
    Settings->GetInt("OutputColorProfileIntent"));

  if (!FileWritten) {
    ptMessageBox::warning(MainWindow, QObject::tr("GraphicsMagick Error"), QObject::tr("No output file written."));
  } else {
    ReportProgress(QObject::tr("Writing output (exif)"));

    if (Settings->GetInt("IncludeExif")) {
      if (!ptImageHelper::WriteExif(Settings->GetString("OutputFileName"),
                                    TheProcessor->m_ExifBuffer,
                                    IptcData,
                                    XmpData))
        ptMessageBox::warning(MainWindow, QObject::tr("Exif Error"), QObject::tr("No exif data written."));
    }
  }

  if (Settings->GetInt("JobMode") == 0) delete OutImage;

  if (JobMode == 0) {
    // This variable carries only JobMode, the "JobMode" settings
    // is also used for full image output.
    ReportProgress(QObject::tr("Writing output (settings)"));

    QFileInfo PathInfo(Settings->GetString("OutputFileName"));
    QString SettingsFileName = PathInfo.dir().path() + "/" + PathInfo.completeBaseName() + ".pts";

    GFilterDM->WritePresetFile(SettingsFileName);
  }

  if (FileWritten) {
    QFileInfo hInfo = QFileInfo(Settings->GetString("OutputFileName"));
    ReportProgress(
      QString(QObject::tr("Written %L1 bytes (%L2 MByte)"))
          .arg(hInfo.size())
          .arg((float)hInfo.size()/1024/1024, 0, 'f', 2)
    );
  } else {
    ReportProgress(QObject::tr("Ready"));
  }
}

//==============================================================================

void WritePipe(QString OutputName = "") {

  if (Settings->GetInt("HaveImage")==0) return;

  QStringList InputFileNameList = Settings->GetStringList("InputFileNameList");
  QFileInfo PathInfo(InputFileNameList[0]);
  QString SuggestedFileName = PathInfo.dir().path() + "/" + PathInfo.completeBaseName();
  if (!Settings->GetString("OutputFileNameSuffix").isEmpty())
    SuggestedFileName += Settings->GetString("OutputFileNameSuffix");
  else
    if (!Settings->GetInt("IsRAW")) SuggestedFileName += "-new";
  QString Pattern;

  switch(Settings->GetInt("SaveFormat")) {
    case ptSaveFormat_JPEG :
      SuggestedFileName += ".jpg";
      Pattern = QObject::tr("Jpg images (*.jpg *.jpeg);;All files (*.*)");
      break;
    case ptSaveFormat_PNG :
    case ptSaveFormat_PNG16 :
      SuggestedFileName += ".png";
      Pattern = QObject::tr("PNG images(*.png);;All files (*.*)");
      break;
    case ptSaveFormat_TIFF8 :
    case ptSaveFormat_TIFF16 :
      SuggestedFileName += ".tif";
      Pattern = QObject::tr("Tiff images (*.tif *.tiff);;All files (*.*)");
      break;
    default :
      SuggestedFileName += ".ppm";
      Pattern = QObject::tr("Ppm images (*.ppm);;All files (*.*)");
      break;
  }

  QString FileName = OutputName;
  if (FileName == "")
    FileName = QFileDialog::getSaveFileName(NULL,
                                            QObject::tr("Save File"),
                                            SuggestedFileName,
                                            Pattern);

  if (0 == FileName.size()) return; // Operation cancelled.

  Settings->SetValue("OutputFileName",FileName);

  // Write out (maybe after applying gamma).
  Update(ptProcessorPhase_WriteOut);
  ImageSaved = 1;
}

////////////////////////////////////////////////////////////////////////////////
//
// Utility : Update GuiSettings->m_?Multiplier according
// to the ColorTemperature/GreenIntensity in the GuiSettings.
//
////////////////////////////////////////////////////////////////////////////////

void CalculateMultipliersFromTemperature() {

  // m_D65Multipliers are supposed to be D65
  // (setting Pre to the ratio of the D65 delivers
  // rgbWB = (x,x,x) and 6500 temperature).

  double RefRGB[3];
  TemperatureToRGB(Settings->GetInt("ColorTemperature"),RefRGB);
  RefRGB[1] /= Settings->GetDouble("GreenIntensity");
  if (TheDcRaw->m_RawColorPhotivo) {
   Settings->SetValue("RMultiplier",
                      VALUE(TheDcRaw->m_D65Multipliers[0])/RefRGB[0]);
   Settings->SetValue("GMultiplier",
                      VALUE(TheDcRaw->m_D65Multipliers[1])/RefRGB[1]);
   Settings->SetValue("BMultiplier",
                      VALUE(TheDcRaw->m_D65Multipliers[2])/RefRGB[2]);
  } else {
    // If not raw color, calculate RefRGB in sRGB back to the camera RGB
    for (short c=0; c<TheDcRaw->m_Colors; c++) {
      double InverseMultiplier = 0;
      for (short cc=0; cc<3; cc++) {
        InverseMultiplier += 1/VALUE(TheDcRaw->m_D65Multipliers[c])
                              * TheDcRaw->m_MatrixSRGBToCamRGB[c][cc]
                              * RefRGB[cc];
      }
      switch (c) {
        case 0 :
          Settings->SetValue("RMultiplier",1/InverseMultiplier);
          break;
        case 1 :
          Settings->SetValue("GMultiplier",1/InverseMultiplier);
          break;
        case 2 :
          Settings->SetValue("BMultiplier",1/InverseMultiplier);
          break;
        default :
          assert(0);
      }
    }
  }
}

//==============================================================================
// Handles the display of the pixel values
void CB_PixelReader(const QPointF Point, const ptPixelReading PixelReading) {
  // No display while pipe is running
  if (Settings->GetInt("PipeIsRunning") || Settings->GetInt("HaveImage")==0) {
    HistogramWindow->PixelInfoHide();
    return;
  }

  RGBValue RGB;

  // reset the display
  if (PixelReading == prNone) {
    HistogramWindow->PixelInfoHide();
    return;
  } else if (PixelReading == prLinear) {
    // we use the AfterEyeCandy image as source
    RGB = TheProcessor->m_Image_AfterEyeCandy->GetRGB(Point.x(), Point.y());
  } else {
    // we use the PreviewImage as source
    RGB = PreviewImage->GetRGB(Point.x(), Point.y());
  }

  // We normalize the values to 0..100
  HistogramWindow->PixelInfo(QString::number(RGB.R/655.35, 'f', 2),
                             QString::number(RGB.G/655.35, 'f', 2),
                             QString::number(RGB.B/655.35, 'f', 2));
}

////////////////////////////////////////////////////////////////////////////////
//
// For convenience...
//
////////////////////////////////////////////////////////////////////////////////

void UpdateSettings() {
  MainWindow->UpdateSettings();
}

////////////////////////////////////////////////////////////////////////////////
//
// Callback pertaining to the Tabs structure (switching page)
//
////////////////////////////////////////////////////////////////////////////////

void CB_Tabs(const short) {
  // If we are previewing according to Tab, we now have to update.
  if (Settings->GetInt("PreviewMode") == ptPreviewMode_Tab) {
    Update(ptProcessorPhase_Preview);
  }
}


////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the Menu structure.
// (inclusive 'X' pressing).
//
////////////////////////////////////////////////////////////////////////////////

void CB_MenuFileOpen(const short HaveFile) {
  QStringList OldInputFileNameList = Settings->GetStringList("InputFileNameList");
  QString InputFileName = "";
  if (HaveFile) {
    InputFileName = ImageFileToOpen;
  }

  // Ask user confirmation
  if (Settings->GetInt("HaveImage") == 1) {
    if (!ptConfirmRequest::saveImage(InputFileName)) {
      return;
    }
  }

  // open file dialog
  if (!HaveFile) {
    InputFileName =
      QFileDialog::getOpenFileName(NULL,
                                   QObject::tr("Open Raw"),
                                   Settings->GetString("RawsDirectory"),
                                   RawPattern);
  }

  // no new image file: abort
  if (InputFileName == "") {
    return;
  } else {
    if (!QFile::exists(InputFileName)) {
      ptMessageBox::warning(NULL,
          QObject::tr("File not found"),
          QObject::tr("Input file does not exist.") + "\n\n" + InputFileName);
      return;
    }
  }

  QFileInfo PathInfo(InputFileName);
  Settings->SetValue("RawsDirectory",PathInfo.absolutePath());
  QStringList InputFileNameList;
  InputFileNameList << PathInfo.absoluteFilePath();

  Settings->SetValue("InputFileNameList",InputFileNameList);

  // Feedback to the user
  ReportProgress(QObject::tr("Reading file"));

  // Test if we can handle the file
  ptDcRaw* TestDcRaw = new(ptDcRaw);
  Settings->ToDcRaw(TestDcRaw);
  uint16_t InputWidth = 0;
  uint16_t InputHeight = 0;
  ptImageType InputType = CheckImageType(InputFileNameList[0],
                                         &InputWidth, &InputHeight,
                                         TestDcRaw);

  if (InputType <= itUndetermined) {
    // not a supported image type or error reading the file
    QString ErrorMessage = QObject::tr("Cannot decode")
                         + " '"
                         + InputFileNameList[0]
                         + "'" ;
    ptMessageBox::warning(MainWindow, "Decode error", ErrorMessage);
    Settings->SetValue("InputFileNameList", OldInputFileNameList);
    delete TestDcRaw;
    return;
  }

  // Image type is supported: process the image
  if (InputType == itRaw) {
    Settings->SetValue("IsRAW", 1);
  } else {
    Settings->SetValue("IsRAW", 0);
  }

  Settings->SetEnabled("UseThumbnail", Settings->GetInt("IsRAW") == 1);

  if (!Settings->useRAWHandling()) {
    Settings->SetValue("ExposureNormalization", 0.0);
  }

  if (Settings->GetInt("HaveImage") == 1) {
    // Catch 1:1 pipe size when opening
    if (Settings->GetInt("PipeSize") == 0)
      Settings->SetValue("PipeSize", 1);
  }

  // TODO mike: need to delete the processor here?
  delete TheDcRaw;
  delete TheProcessor;
  // Load user settings
  if (Settings->GetInt("StartupSettings") == 1 &&
      Settings->GetInt("StartupSettingsReset") == 1 &&
      Settings->GetInt("HaveImage") == 1) {
    Settings->SetValue("HaveImage", 0);
    CB_OpenSettingsFile(Settings->GetString("StartupSettingsFile"));
    // clean up
    QStringList Temp;
    Temp << "CropX" << "CropY" << "CropW" << "CropH";
    Temp << "RotateW" << "RotateH";
    for (int i = 0; i < Temp.size(); i++) Settings->SetValue(Temp.at(i), 0);
  }

  if (Settings->GetString("Sidecar") != "") {
    ReadSidecar(Settings->GetString("Sidecar"));
  }

  // load settings file automatically
  QString SettingsFileName = PathInfo.dir().path() + "/" + PathInfo.completeBaseName() + ".pts";
  if (QFile::exists(SettingsFileName)) {
    Settings->SetValue("HaveImage", 0);
    CB_OpenSettingsFile(SettingsFileName);
  }

  // clean up possible detail view cache
  if (Settings->GetInt("DetailViewActive") == 1) {
    Settings->SetValue("DetailViewActive", 0);
    Settings->SetValue("DetailViewCropX",  0);
    Settings->SetValue("DetailViewCropY",  0);
    Settings->SetValue("DetailViewCropW",  0);
    Settings->SetValue("DetailViewCropH",  0);
    Settings->ToDcRaw(TestDcRaw);
  }

  TheDcRaw = TestDcRaw;
  Settings->SetValue("ImageW", InputWidth);
  Settings->SetValue("ImageH", InputHeight);

  if (Settings->GetInt("StartupSwitchAR")) {
    // portrait image
    if ((!Settings->useRAWHandling() &&
         Settings->GetInt("ImageW") < Settings->GetInt("ImageH")) ||
        TheDcRaw->m_Flip & 4) {
      if (Settings->GetInt("AspectRatioW") > Settings->GetInt("AspectRatioH"))
        CB_CropOrientationButton();
    } else { // landscape
      if (Settings->GetInt("AspectRatioW") < Settings->GetInt("AspectRatioH"))
        CB_CropOrientationButton();
    }
  }
  SetRatingFromXmp();
  if (Settings->GetInt("LoadTags"))
    SetTagsFromXmp();

  // reflect RAW or bitmap in GUI
  MainWindow->UpdateToolBoxes();

  TheProcessor = new ptProcessor(ReportProgress);
  TheProcessor->m_DcRaw = TheDcRaw;

  Settings->SetValue("HaveImage", 1);
  short OldRunMode = Settings->GetInt("RunMode");
  Settings->SetValue("RunMode", 0);

  if (Settings->GetInt("AutomaticPipeSize") && Settings->ToolIsActive("TabResize")) {
    CalculatePipeSize(true);
  }

  Update(ptProcessorPhase_Raw, ptProcessorPhase_Load, 0);

  MainWindow->UpdateExifInfo(TheProcessor->m_ExifData);
  Settings->SetValue("PerspectiveFocalLength", Settings->GetDouble("FocalLengthIn35mmFilm"));
  Settings->SetValue("DefishFocalLength",      Settings->GetDouble("FocalLengthIn35mmFilm"));
  Settings->SetValue("LfunFocal",              Settings->GetDouble("FocalLengthIn35mmFilm"));
  double TmpAprt = Settings->GetDouble("ApertureFromExif");
  Settings->SetValue("LfunAperture",           (TmpAprt == 0.0) ? 8.0 : TmpAprt);
  MainWindow->UpdateFilenameInfo(Settings->GetStringList("InputFileNameList"));

  QFileInfo finfo = QFileInfo((Settings->GetStringList("InputFileNameList"))[0]);
  MainWindow->setWindowTitle(
      QString("%1 - %2 - Photivo").arg(finfo.fileName())
                                  .arg(QDir::toNativeSeparators(finfo.absolutePath())) );

  Settings->SetValue("RunMode",OldRunMode);

  ptClearUndoRedo();
}

//==============================================================================

void CB_MenuFileSaveOutput(QString OutputName = "") {
  try {
    if (Settings->GetInt("HaveImage")==0) return;

    QStringList InputFileNameList = Settings->GetStringList("InputFileNameList");
    QFileInfo PathInfo(InputFileNameList[0]);
    QString SuggestedFileName = PathInfo.dir().path() + "/" + PathInfo.completeBaseName();
    if (!Settings->GetString("OutputFileNameSuffix").isEmpty())
      SuggestedFileName += Settings->GetString("OutputFileNameSuffix");
    else
      if (!Settings->GetInt("IsRAW")) SuggestedFileName += "-new";
    QString Pattern;

    switch(Settings->GetInt("SaveFormat")) {
      case ptSaveFormat_JPEG :
        SuggestedFileName += ".jpg";
        Pattern = QObject::tr("Jpg images (*.jpg *.jpeg);;All files (*.*)");
        break;
      case ptSaveFormat_PNG :
      case ptSaveFormat_PNG16 :
        SuggestedFileName += ".png";
        Pattern = QObject::tr("PNG images(*.png);;All files (*.*)");
        break;
      case ptSaveFormat_TIFF8 :
      case ptSaveFormat_TIFF16 :
        SuggestedFileName += ".tif";
        Pattern = QObject::tr("Tiff images (*.tif *.tiff);;All files (*.*)");
        break;
      default :
        SuggestedFileName += ".ppm";
        Pattern = QObject::tr("Ppm images (*.ppm);;All files (*.*)");
        break;
    }

    QString FileName = OutputName;
    if (FileName == "")
      FileName = QFileDialog::getSaveFileName(NULL,
                                              QObject::tr("Save File"),
                                              SuggestedFileName,
                                              Pattern);

    if (0 == FileName.size()) return; // Operation cancelled.

    Settings->SetValue("OutputFileName",FileName);

    short OldRunMode = Settings->GetInt("RunMode");
    Settings->SetValue("RunMode",0);

    // Processing the job.
    delete TheDcRaw;
    delete TheProcessor;
  TheDcRaw = new(ptDcRaw);
    TheProcessor = new ptProcessor(ReportProgress);
    Settings->SetValue("JobMode",1); // Disable caching to save memory
    TheProcessor->m_DcRaw = TheDcRaw;
    Settings->ToDcRaw(TheDcRaw);
    // Run the graphical pipe in full format mode to recreate the image.
    Settings->SetValue("FullOutput",1);
    TheProcessor->Run(ptProcessorPhase_Raw,ptProcessorPhase_Load,1,1);
    Settings->SetValue("FullOutput",0);

    // Write out (maybe after applying gamma).
    Update(ptProcessorPhase_WriteOut);

    delete TheDcRaw;
    delete TheProcessor;
  TheDcRaw = new(ptDcRaw);
    TheProcessor = new ptProcessor(ReportProgress);
    Settings->SetValue("JobMode",0);
    TheProcessor->m_DcRaw = TheDcRaw;
    Settings->ToDcRaw(TheDcRaw);
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Load,0);
    MainWindow->UpdateExifInfo(TheProcessor->m_ExifData);
    Settings->SetValue("RunMode",OldRunMode);
    ImageSaved = 1;
  } catch (std::bad_alloc) {
    ptMessageBox::critical(NULL, "Memory error", "No file written, memory error.");
  }
} // CB_MenuFileSaveOutput

//==============================================================================

void CB_MenuFileExit(const short) {
  if (Settings->GetInt("HaveImage")==1 && ImageSaved == 0) {
    if (!ptConfirmRequest::saveImage()) {
      return;
    }
  }

  // Save current filter config
  if (Settings->GetString("StartupSettingsFile").endsWith("Presets/Latest.pts"))
    GFilterDM->WritePresetFile(Settings->GetString("UserDirectory") + "/Presets/Latest.pts");

  // clean up the input file if we got just a temp file
  if (Settings->GetInt("HaveImage")==1 && ImageCleanUp == 1) {
    QString OldInputFileName = Settings->GetStringList("InputFileNameList")[0];
    ptRemoveFile(OldInputFileName);
    ImageCleanUp--;
  }

  // Delete backup settingsfile
  if (Settings->GetInt("WriteBackupSettings") == 1)
    QFile::remove(Settings->GetString("UserDirectory")+"/backup.pts");

  MainWindow->Form_2_Settings();

  printf("Saving settings ...\n");

#ifndef PT_WITHOUT_FILEMGR
  // this also writes settings.
  delete FileMgrWindow;
#endif

  delete BatchWindow;

  // Store the position of the splitter and main window
  Settings->m_IniSettings->
    setValue("MainSplitter",MainWindow->MainSplitter->saveState());
  Settings->m_IniSettings->
    setValue("ControlSplitter",MainWindow->ControlSplitter->saveState());
  Settings->m_IniSettings->setValue("MainWindowPos",    MainWindow->pos());
  Settings->m_IniSettings->setValue("MainWindowSize",   MainWindow->size());
  Settings->m_IniSettings->setValue("IsMaximized",      MainWindow->windowState() == Qt::WindowMaximized);
  Settings->m_IniSettings->setValue("MainWindowScreen", qApp->desktop()->screenNumber(MainWindow));

  // Store the version of the settings and files
  Settings->m_IniSettings->setValue("SettingsVersion",PhotivoSettingsVersion);

  // Explicitly. The destructor of it cares for persistent settings.
  delete Settings;
#ifdef Q_OS_WIN
  if (!JobMode)
    ptEcWin7::DestroyInstance();
#endif

  ALLOCATED(10000000);
  printf("Exiting Photivo.\n");
  QCoreApplication::exit(EXIT_SUCCESS);
}

////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the gimp
//
////////////////////////////////////////////////////////////////////////////////

void CB_ExportToGimpCheck(const QVariant Check) {
  Settings->SetValue("ExportToGimp",Check);
}

void GimpExport(const short UsePipe) {
  try {
    if (Settings->GetInt("HaveImage")==0) return;

    ReportProgress(QObject::tr("Writing tmp image for gimp"));

    short OldRunMode = Settings->GetInt("RunMode");

    ptImage* ImageForGimp = new ptImage;

    if (UsePipe == 1)
      ImageForGimp->Set(TheProcessor->m_Image_AfterEyeCandy);
    else {
      Settings->SetValue("RunMode",0);

      // Processing the job.
      delete TheDcRaw;
      delete TheProcessor;
    TheDcRaw = new(ptDcRaw);
      TheProcessor = new ptProcessor(ReportProgress);
      Settings->SetValue("JobMode",1); // Disable caching to save memory

      TheProcessor->m_DcRaw = TheDcRaw;
      Settings->ToDcRaw(TheDcRaw);
      // Run the graphical pipe in full format mode to recreate the image.
      Settings->SetValue("FullOutput",1);
      TheProcessor->Run(ptProcessorPhase_Raw,ptProcessorPhase_Load,1,1);
      Settings->SetValue("FullOutput",0);

      ImageForGimp = TheProcessor->m_Image_AfterEyeCandy; // no cache
    }

    BeforeGamma(ImageForGimp);

    cmsHPROFILE OutputColorProfile = NULL;

    ReportProgress(QObject::tr("Converting to output space"));

    // Prepare and open an output profile.
    OutputColorProfile = cmsOpenProfileFromFile(
      Settings->GetString("OutputColorProfile").toLocal8Bit().data(),
      "r");
    if (!OutputColorProfile) {
      ptLogError(ptError_FileOpen,
     Settings->GetString("OutputColorProfile").toLocal8Bit().data());
      assert(OutputColorProfile);
    }

    // Color space conversion
    ptImage* ReturnValue = ImageForGimp->lcmsRGBToRGB(
      OutputColorProfile,
      Settings->GetInt("OutputColorProfileIntent"),
      Settings->GetInt("CMQuality"));
    if (!ReturnValue) {
      ptLogError(ptError_lcms,"lcmsRGBToRGB(OutputColorProfile)");
      assert(ReturnValue);
    }

    AfterAll(ImageForGimp);

    if (GFilterDM->GetFilterFromName(Fuid::Wiener_Out)->isActive())
      EndSharpen(ImageForGimp, OutputColorProfile, Settings->GetInt("OutputColorProfileIntent"));
    // Close the output profile.
    cmsCloseProfile(OutputColorProfile);

    QTemporaryFile ImageFile;
    ImageFile.setFileTemplate(QDir::tempPath()+"/XXXXXX.ppm");
    bool result = ImageFile.open();
    assert (result);
    QString ImageFileName = ImageFile.fileName();
    ImageFile.setAutoRemove(false);
    ImageFile.close();
    printf("(%s,%d) '%s'\n",
           __FILE__,__LINE__,ImageFileName.toLocal8Bit().data());
    ImageForGimp->WriteAsPpm(ImageFileName.toLocal8Bit().data(),16);

    ReportProgress(QObject::tr("Writing tmp exif for gimp"));

    QTemporaryFile ExifFile;
    result = ExifFile.open();
    assert (result);
    QString ExifFileName = ExifFile.fileName();
    ExifFile.setAutoRemove(false);
    printf("(%s,%d) '%s'\n",
           __FILE__,__LINE__,ExifFileName.toLocal8Bit().data());
    QDataStream ExifOut(&ExifFile);
    ExifOut.writeRawData((char *) TheProcessor->m_ExifBuffer.data(),
                         TheProcessor->m_ExifBuffer.size());
    ExifFile.close();

    ReportProgress(QObject::tr("Writing tmp icc for gimp"));

    QTemporaryFile ICCFile;
    result = ICCFile.open();
    assert (result);
    QString ICCFileName = ICCFile.fileName();
    ICCFile.setAutoRemove(false);
    printf("(%s,%d) '%s'\n",
           __FILE__,__LINE__,ICCFileName.toLocal8Bit().data());
    QDataStream ICCOut(&ICCFile);
    FILE* pFile = fopen ( Settings->GetString("OutputColorProfile").toLocal8Bit().data(), "rb" );
    if (pFile==NULL) {
      ptLogError(ptError_FileOpen,Settings->GetString("OutputColorProfile").toLocal8Bit().data());
      exit(EXIT_FAILURE);
    }
    fseek (pFile , 0 , SEEK_END);
    long lSize = ftell (pFile);
    rewind (pFile);

    char* pchBuffer = (char*) CALLOC2(lSize,sizeof(uint8_t));
    ptMemoryError(pchBuffer,__FILE__,__LINE__);

    size_t RV = fread (pchBuffer, 1, lSize, pFile);
    if (RV != (size_t) lSize) {
      ptLogError(ptError_FileOpen,Settings->GetString("OutputColorProfile").toLocal8Bit().data());
      exit(EXIT_FAILURE);
    }
    ICCOut.writeRawData(pchBuffer, lSize);

    FCLOSE (pFile);
    FREE2 (pchBuffer);
    ICCFile.close();

    ReportProgress(QObject::tr("Writing gimp interface file"));

    QTemporaryFile GimpFile;
    GimpFile.setFileTemplate(QDir::tempPath()+"/XXXXXX.ptg");
    result = GimpFile.open();
    assert (result);
    QString GimpFileName = GimpFile.fileName();
    GimpFile.setAutoRemove(false);
    printf("(%s,%d) '%s'\n",
           __FILE__,__LINE__,GimpFileName.toLocal8Bit().data());
    QTextStream Out(&GimpFile);
    Out << ImageFileName << "\n";
    Out << ExifFileName << "\n";
    Out << ICCFileName << "\n";
    GimpFile.close();

    QString GimpExeCommand = Settings->GetString("GimpExecCommand");
    QStringList GimpArguments;
    GimpArguments << GimpFileName;
    QProcess* GimpProcess = new QProcess();
    GimpProcess->startDetached(GimpExeCommand,GimpArguments);

    // clean up
    if (UsePipe == 1) {
      delete ImageForGimp;
    } else {
      delete TheDcRaw;
      delete TheProcessor;
    TheDcRaw = new(ptDcRaw);
      TheProcessor = new ptProcessor(ReportProgress);
      Settings->SetValue("JobMode",0);
      TheProcessor->m_DcRaw = TheDcRaw;
      Settings->ToDcRaw(TheDcRaw);
      Update(ptProcessorPhase_Raw,ptProcessorPhase_Load,0);
      MainWindow->UpdateExifInfo(TheProcessor->m_ExifData);
      Settings->SetValue("RunMode",OldRunMode);
    }
    ReportProgress(QObject::tr("Ready"));
  } catch (std::bad_alloc) {
    ptMessageBox::critical(NULL, "Memory error", "No export, memory error.");
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the Zoom function
//
////////////////////////////////////////////////////////////////////////////////

void CB_ZoomInput(const QVariant Value) {
  ViewWindow->ZoomTo(Value.toFloat() / 100.0);
  // TODO SR: updatesettings doesnt seem to be necessary
//    MainWindow->UpdateSettings(); // To reflect maybe new zoom
}

void CB_ZoomFitButton() {
  ViewWindow->ZoomToFit();
  //MainWindow->UpdateSettings(); // To reflect maybe new zoom
}

void CB_ZoomFullButton() {
  CB_ZoomInput(100);
}

void CB_ZoomStep(int direction) {
  ViewWindow->ZoomStep(direction);
}

void CB_BatchButton() {
  MainWindow->OpenBatchWindow();
}

void CB_FileMgrButton() {
#ifndef PT_WITHOUT_FILEMGR
  MainWindow->OpenFileMgrWindow();
#endif
}

void CB_FullScreenButton(const int State) {
  if (State == 1) {
    MainWindow->raise();
    MainWindow->showFullScreen();
    Settings->SetValue("FullscreenActive", 1);
  } else {
    MainWindow->showNormal();
    Settings->SetValue("FullscreenActive", 0);
  }
  MainWindow->FullScreenButton->setChecked(State);
}

////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the prev/next image buttons
//
////////////////////////////////////////////////////////////////////////////////

QFileInfoList ptGetFilesInTheImageFolder() {
  QStringList InputFileNameList = Settings->GetStringList("InputFileNameList");
  QFileInfo PathInfo(InputFileNameList[0]);
  QString fileName = PathInfo.fileName();
  QDir fileDir = PathInfo.dir();
  fileDir.setFilter(QDir::Files);
  fileDir.setSorting(QDir::Name);

  QStringList fileExts;
  if (Settings->GetInt("FileMgrShowRAWs")) {
    fileExts << FileExtsRaw;
  }
  if (Settings->GetInt("FileMgrShowBitmaps")) {
    fileExts << FileExtsBitmap;
  }
  return fileDir.entryInfoList(fileExts);
}

QString ptGetImageFileName() {
  QStringList InputFileNameList = Settings->GetStringList("InputFileNameList");
  QFileInfo PathInfo(InputFileNameList[0]);
  return PathInfo.fileName();
}

void CB_PreviousImageButton() {
  if (Settings->GetInt("HaveImage")==0) return;

  QString fileName = ptGetImageFileName();
  QFileInfoList files = ptGetFilesInTheImageFolder();

  for (int i = 1; i < files.count(); i++) {
    if (files[i].fileName() == fileName) {
      ImageFileToOpen = files[i - 1].absoluteFilePath();
      CB_MenuFileOpen(1);
      break;
    }
  }
}

void CB_NextImageButton() {
  if (Settings->GetInt("HaveImage")==0) return;

  QString fileName = ptGetImageFileName();
  QFileInfoList files = ptGetFilesInTheImageFolder();

  for (int i = 0; i < files.count() - 1; i++) {
    if (files[i].fileName() == fileName) {
      ImageFileToOpen = files[i + 1].absoluteFilePath();
      CB_MenuFileOpen(1);
      break;
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the Generic Tab
//
////////////////////////////////////////////////////////////////////////////////

void CB_CameraColorChoice(const QVariant Choice) {
  int PreviousChoice = Settings->GetInt("CameraColor");
  Settings->SetValue("CameraColor",Choice);
  if (Choice.toInt() == value(ptCameraColor::Profile)) {
    if (!Settings->GetString("CameraColorProfile").size()) {
      ptMessageBox::warning(MainWindow,
                           QObject::tr("Please load a profile first"),
                           QObject::tr("Please load a profile first"));
      Settings->SetValue("CameraColor",PreviousChoice);
    }
  }

  if (Settings->GetInt("CameraColor") == value(ptCameraColor::Embedded)) {
    // TODO
    ptMessageBox::warning(MainWindow,
                        QObject::tr("Not yet implemented"),
                        QObject::tr("Not yet implemented. Reverting to Adobe."));
    Settings->SetValue("CameraColor",value(ptCameraColor::Adobe_Matrix));
  }

  Update(ptProcessorPhase_Raw,ptProcessorPhase_Highlights);
}

void CB_CameraColorProfileButton() {
  QString ProfileFileName = QFileDialog::getOpenFileName(
                 NULL,
                 QObject::tr("Open Profile"),
                 Settings->GetString("CameraColorProfilesDirectory"),
                 ProfilePattern);

  if (0 == ProfileFileName.size() ) {
    // Canceled just return
    return;
  } else {
    QFileInfo PathInfo(ProfileFileName);
    Settings->SetValue("CameraColorProfilesDirectory",PathInfo.absolutePath());
    Settings->SetValue("CameraColorProfile",PathInfo.absoluteFilePath());
    Settings->SetValue("CameraColor",value(ptCameraColor::Profile));
  }

  Update(ptProcessorPhase_Raw,ptProcessorPhase_Highlights);
}

void CB_CameraColorProfileIntentChoice(const QVariant Choice) {
  Settings->SetValue("CameraColorProfileIntent",Choice);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Highlights);
}

void CB_CameraColorGammaChoice(const QVariant Choice) {
  Settings->SetValue("CameraColorGamma",Choice);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Highlights);
}

void CB_WorkColorChoice(const QVariant Choice) {
  Settings->SetValue("WorkColor",Choice);
  PreCalcTransforms();
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Highlights);
}

void CB_CMQualityChoice(const QVariant Choice) {
  Settings->SetValue("CMQuality",Choice);
  MainWindow->UpdateSettings();
  PreCalcTransforms();
  Update(ptProcessorPhase_Preview);
}

void CB_PreviewColorProfileButton() {
  QString ProfileFileName = QFileDialog::getOpenFileName(
                 NULL,
                 QObject::tr("Open Profile"),
                 Settings->GetString("PreviewColorProfilesDirectory"),
                 ProfilePattern);

  if (0 == ProfileFileName.size() ) {
    // Canceled just return
    return;
  } else {
    QFileInfo PathInfo(ProfileFileName);
    Settings->SetValue("PreviewColorProfilesDirectory",PathInfo.absolutePath());
    Settings->SetValue("PreviewColorProfile",PathInfo.absoluteFilePath());
  }

  // Reflect in gui.
  MainWindow->UpdateSettings();

  // Close old profile and open new one.
  cmsCloseProfile(PreviewColorProfile);
  PreviewColorProfile = cmsOpenProfileFromFile(
          Settings->GetString("PreviewColorProfile").toLocal8Bit().data(),
          "r");
  if (!PreviewColorProfile) {
    ptLogError(ptError_FileOpen,
         Settings->GetString("PreviewColorProfile").toLocal8Bit().data());
    assert(PreviewColorProfile);
  }
  PreCalcTransforms();
  // And update the preview.
  Update(ptProcessorPhase_Preview);
}

void CB_PreviewColorProfileIntentChoice(const QVariant Choice) {
  Settings->SetValue("PreviewColorProfileIntent",Choice);
  PreCalcTransforms();
  Update(ptProcessorPhase_Preview);
}

void CB_OutputColorProfileButton() {
  QString ProfileFileName = QFileDialog::getOpenFileName(
                 NULL,
                 QObject::tr("Open Profile"),
                 Settings->GetString("OutputColorProfilesDirectory"),
                 ProfilePattern);

  if (0 == ProfileFileName.size() ) {
    // Canceled just return
    return;
  } else {
    QFileInfo PathInfo(ProfileFileName);
    Settings->SetValue("OutputColorProfilesDirectory",PathInfo.absolutePath());
    Settings->SetValue("OutputColorProfile",PathInfo.absoluteFilePath());
  }

  // Reflect in gui.
  if (Settings->GetInt("HistogramMode")==ptHistogramMode_Output) {
    if (Settings->GetInt("ExposureIndicator")==1) {
      Update(ptProcessorPhase_Preview);
    } else {
      Update(ptProcessorPhase_OnlyHistogram);
    }
  }
  MainWindow->UpdateSettings();
}

void CB_OutputColorProfileResetButton() {
  Settings->SetValue("OutputColorProfile",
                     (Settings->GetString("UserDirectory") + "Profiles/Output/sRGB.icc").toLocal8Bit().data());
  if (Settings->GetInt("HistogramMode")==ptHistogramMode_Output) {
    if (Settings->GetInt("ExposureIndicator")==1) {
      Update(ptProcessorPhase_Preview);
    } else {
      Update(ptProcessorPhase_OnlyHistogram);
    }
  }
  MainWindow->UpdateSettings();
}

void CB_OutputColorProfileIntentChoice(const QVariant Choice) {
  Settings->SetValue("OutputColorProfileIntent",Choice);
  if (Settings->GetInt("HistogramMode")==ptHistogramMode_Output) {
    if (Settings->GetInt("ExposureIndicator")==1) {
      Update(ptProcessorPhase_Preview);
    } else {
      Update(ptProcessorPhase_OnlyHistogram);
    }
  }
  MainWindow->UpdateSettings();
}

void CB_StyleChoice(const QVariant Choice) {
  Settings->SetValue("Style", Choice);
  ptTheme::Theme newTheme = (ptTheme::Theme)Choice.toInt();

  if (newTheme == ptTheme::thNone) {
    Settings->SetValue("CustomCSSFile","");
    Theme->setCustomCSS("");
  }

  Theme->SwitchTo(newTheme, (ptTheme::Highlight)Settings->GetInt("StyleHighLight"));

#ifdef Q_OS_MAC
//dirty hack to make theming work
  MainWindow->MainSplitter->setStyleSheet("");
#endif

  MainWindow->MainTabBook->setStyle(Theme->style());
  MainWindow->ProcessingTabBook->setStyle(Theme->style());
  MainWindow->BottomContainer->setStyle(Theme->style());
  MainWindow->PipeControlWidget->setStyle(Theme->style());
  MainWindow->MainSplitter->setStyle(Theme->style());
  MainWindow->ControlSplitter->setStyle(Theme->style());
  MainWindow->ViewSplitter->setStyle(Theme->style());
  MainWindow->ViewStartPage->setStyle(Theme->style());

  TheApplication->setPalette(Theme->palette());

  MainWindow->MainTabBook->setStyleSheet(Theme->stylesheet());
  MainWindow->BottomContainer->setStyleSheet(Theme->stylesheet());
  MainWindow->PipeControlWidget->setStyleSheet(Theme->stylesheet());
  MainWindow->StatusWidget->setStyleSheet(Theme->stylesheet());
  MainWindow->SearchWidget->setStyleSheet(Theme->stylesheet());
  MainWindow->ViewStartPageFrame->setStyleSheet(Theme->stylesheet());
#ifdef Q_OS_MAC
//dirty hack to make theming work
  if(Theme->MacStyleFlag){
  MainWindow->MainSplitter->setStyleSheet("background-color:"+Theme->MacBackGround+";");
  }
#endif
  MainWindow->UpdateToolBoxes();
  SetBackgroundColor(Settings->GetInt("BackgroundColor"));
  CB_SliderWidthInput(Settings->GetInt("SliderWidth"));
#ifndef PT_WITHOUT_FILEMGR
  FileMgrWindow->updateTheme();
#endif
  BatchWindow->UpdateTheme();
}

void CB_StyleHighLightChoice(const QVariant Choice) {
  Settings->SetValue("StyleHighLight",Choice);
  CB_StyleChoice(Settings->GetInt("Style"));
}

void CB_LoadStyleButton() {
  QString FileName;

  FileName = QFileDialog::getOpenFileName(NULL,
    QObject::tr("Open Image"),
    Settings->GetString("UserDirectory"),
    QObject::tr("CSS files (*.css *.qss);;All files(*.*)")
  );

  if (FileName.size() == 0) {
    return;
  } else {
    Settings->SetValue("CustomCSSFile", FileName);
    Theme->setCustomCSS(FileName);
    CB_StyleChoice(Settings->GetInt("Style"));
  }
}


void CB_SaveConfirmationCheck(const QVariant Check) {
  Settings->SetValue("SaveConfirmation",Check);
  MainWindow->AutosaveSettingsWidget->setDisabled(Settings->GetInt("SaveConfirmation"));
}

void CB_AutosaveSettingsCheck(const QVariant Check) {
  Settings->SetValue("AutosaveSettings",Check);
}

void CB_ResetSettingsConfirmationCheck(const QVariant Check) {
  Settings->SetValue("ResetSettingsConfirmation",Check);
}

void CB_FullPipeConfirmationCheck(const QVariant Check) {
  Settings->SetValue("FullPipeConfirmation",Check);
}


void CB_WriteBackupSettingsCheck(const QVariant Check) {
  Settings->SetValue("WriteBackupSettings",Check);
}

void CB_MemoryTestInput(const QVariant Value) {
  Settings->SetValue("MemoryTest",0);
  if (Value.toInt()>0)
    if (ptMessageBox::question(MainWindow,
        QObject::tr("Are you sure?"),
        QObject::tr("If you don't stop me, I will waste %1 MB of memory.").arg(Value.toInt()),
        QMessageBox::Ok,QMessageBox::Cancel)==QMessageBox::Ok){
      // allocate orphaned memory for testing
      char (*Test) = (char (*)) CALLOC2(Value.toInt()*1024*1024,1);
      memset(Test, '\0', Value.toInt()*1024*1024);
      ptMessageBox::critical(0,"Feedback","Memory wasted ;-)");
    }
}

void SetDetailViewRect(const ptStatus, QRect rect) {
  DetailViewRect = rect;
}
void CB_PipeSizeChoice(const QVariant Choice) {

  short PreviousPipeSize = Settings->GetInt("PipeSize");

  if (Choice == ptPipeSize_Full &&
      Settings->GetInt("FullPipeConfirmation")==1 &&
      (Settings->GetInt("ImageH") > 2000 ||
       Settings->GetInt("ImageW") > 2000)) {
    ptMessageBox msgBox;
    msgBox.setWindowTitle(QObject::tr("Really switch to 1:1 pipe?"));
    msgBox.setText(QObject::tr("Switching to 1:1 pipe will increase memory usage and processing time greatly.\nAre you sure?"));

    QPushButton *DetailButton = msgBox.addButton(QObject::tr("Detail view"), QMessageBox::ActionRole);
    QPushButton *CancelButton = msgBox.addButton(QMessageBox::Cancel);
    msgBox.addButton(QMessageBox::Ok);

    msgBox.exec();

    if (msgBox.clickedButton() == CancelButton ||
        (msgBox.clickedButton() == DetailButton &&
         Settings->GetInt("HaveImage")==0)) {
      Settings->SetValue("PipeSize", PreviousPipeSize);
      return;
    } else if (msgBox.clickedButton() == DetailButton &&
               Settings->GetInt("HaveImage") == 1) {
      short OldZoomMode = 0;
      if (Settings->GetInt("DetailViewActive") == 0) {
        Settings->SetValue("DetailViewScale", PreviousPipeSize);
        if (TheProcessor->m_Image_DetailPreview == NULL)
          TheProcessor->m_Image_DetailPreview = new ptImage();
        // save a full image if we have several detail views without full view
        if (!Settings->useRAWHandling()) {
          TheProcessor->m_Image_DetailPreview->SetScaled(TheProcessor->m_Image_AfterDcRaw,
                                                         Settings->GetInt("Scaled"));
        } else {
          TheProcessor->m_Image_DetailPreview->Set(TheProcessor->m_Image_AfterDcRaw);
        }
      }
      // display the cache image and start crop mode

      // First : make sure we have Image_AfterDcRaw in the view window.
      // Anything else might have undergone geometric transformations that are
      // impossible to calculate reverse to a spot in dcraw.
      OldZoomMode = Settings->GetInt("ZoomMode");
      ViewWindow->ZoomToFit();
      UpdatePreviewImage(TheProcessor->m_Image_DetailPreview);

      // Allow to be selected in the view window. And deactivate main.
      BlockTools(btmBlockAll);
      ViewWindow->ShowStatus(QObject::tr("Detail view"));
      ViewWindow->StartSimpleRect(SetDetailViewRect);
      while (ViewWindow->interaction() == iaSelectRect) {
        QApplication::processEvents();
      }

      // Selection is done at this point. Disallow it further and activate main.
      BlockTools(btmUnblock);

      if (DetailViewRect.width() >>4 <<4 > 19 &&
          DetailViewRect.height() >>4 <<4 > 19) {
        Settings->SetValue("DetailViewActive",1);
        short CachedPipeSize = Settings->GetInt("DetailViewScale");
        Settings->SetValue("DetailViewCropX", ((DetailViewRect.left() >>4) <<4) << CachedPipeSize);
        Settings->SetValue("DetailViewCropY", ((DetailViewRect.top() >>4) <<4) << CachedPipeSize);
        Settings->SetValue("DetailViewCropW", ((DetailViewRect.width() >>4) <<4) << CachedPipeSize);
        Settings->SetValue("DetailViewCropH", ((DetailViewRect.height() >>4) <<4) << CachedPipeSize);
        Settings->ToDcRaw(TheDcRaw);
      } else {
        ptMessageBox::information(NULL,"No crop","Too small. Please try again!");
        //ViewWindow->Zoom(OldZoom,0);    // TODOSR: re-enable
        Settings->SetValue("ZoomMode",OldZoomMode);
        Update(ptProcessorPhase_Preview);
        if (Settings->GetInt("DetailViewActive")==0) {
          delete TheProcessor->m_Image_DetailPreview;
          TheProcessor->m_Image_DetailPreview = NULL;
        }
        Settings->SetValue("PipeSize",PreviousPipeSize);
        return;
      }

      //ViewWindow->Zoom(OldZoom,0);    // TODOSR: re-enable
      Settings->SetValue("ZoomMode",OldZoomMode);

    } else { // 1:1 full mode
      // clean up possible detail view cache
      if (TheProcessor->m_Image_DetailPreview != NULL) {
        delete TheProcessor->m_Image_DetailPreview;
        TheProcessor->m_Image_DetailPreview = NULL;
        Settings->SetValue("DetailViewActive",0);
        Settings->SetValue("DetailViewCropX", 0);
        Settings->SetValue("DetailViewCropY", 0);
        Settings->SetValue("DetailViewCropW", 0);
        Settings->SetValue("DetailViewCropH", 0);
        Settings->ToDcRaw(TheDcRaw);
      }
    }
  } else {
    // clean up possible detail view cache
    if (TheProcessor->m_Image_DetailPreview != NULL) {
      delete TheProcessor->m_Image_DetailPreview;
      TheProcessor->m_Image_DetailPreview = NULL;
      Settings->SetValue("DetailViewActive",0);
      Settings->SetValue("DetailViewCropX", 0);
      Settings->SetValue("DetailViewCropY", 0);
      Settings->SetValue("DetailViewCropW", 0);
      Settings->SetValue("DetailViewCropH", 0);
      Settings->ToDcRaw(TheDcRaw);
    }
  }

  MainWindow->UpdateToolBoxes();

  Settings->SetValue("PipeSize",Choice);
  short PipeSize = Settings->GetInt("PipeSize");
  short Expansion = PreviousPipeSize-PipeSize;

  if (Settings->GetInt("DetailViewActive")==0) {
    // Following adaptation is needed for the case spot WB is in place.
    if (Expansion > 0) {
      Settings->SetValue("VisualSelectionX",
                         Settings->GetInt("VisualSelectionX")<<Expansion);
      Settings->SetValue("VisualSelectionY",
                         Settings->GetInt("VisualSelectionY")<<Expansion);
      Settings->SetValue("VisualSelectionWidth",
                         Settings->GetInt("VisualSelectionWidth")<<Expansion);
      Settings->SetValue("VisualSelectionHeight",
                         Settings->GetInt("VisualSelectionHeight")<<Expansion);
    } else {
      Expansion = -Expansion;
      Settings->SetValue("VisualSelectionX",
                         Settings->GetInt("VisualSelectionX")>>Expansion);
      Settings->SetValue("VisualSelectionY",
                         Settings->GetInt("VisualSelectionY")>>Expansion);
      Settings->SetValue("VisualSelectionWidth",
                         Settings->GetInt("VisualSelectionWidth")>>Expansion);
      Settings->SetValue("VisualSelectionHeight",
                         Settings->GetInt("VisualSelectionHeight")>>Expansion);
    }
  }

  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
  MainWindow->UpdateExifInfo(TheProcessor->m_ExifData);

  ALLOCATED(10000000);
}


void CB_RunModeCheck(const QVariant Check) {
  Settings->SetValue("RunMode",Check);
  MainWindow->UpdateSettings();
  if (Settings->GetInt("RunMode")==0) {
    Update(ptProcessorPhase_Output);
  }
}

void CB_PreviewModeButton(const QVariant State) {
  Settings->SetValue("PreviewTabMode",State);
  if (Settings->GetInt("PreviewTabMode")) {
    Settings->SetValue("PreviewMode",ptPreviewMode_Tab);
    MainWindow->PreviewModeButton->setChecked(1);
  } else {
    Settings->SetValue("PreviewMode",ptPreviewMode_End);
    MainWindow->PreviewModeButton->setChecked(0);
  }
  Update(ptProcessorPhase_Preview);
}

void CB_RunButton() {
  short OldRunMode = Settings->GetInt("RunMode");
  Settings->SetValue("RunMode",0);
  Update(ptProcessorPhase_Output);
  Settings->SetValue("RunMode",OldRunMode);
  MainWindow->UpdateSettings();
}

void ResetButtonHandler(const short mode) {
  int DoOpen = 1;
  if (mode == ptResetMode_Full) { // full reset
    if ( Settings->GetInt("ResetSettingsConfirmation")==1 ) {
      ptMessageBox msgBox;
      msgBox.setIcon(QMessageBox::Question);
      msgBox.setWindowTitle(QObject::tr("Reset?"));
      msgBox.setText(QObject::tr("Reset to neutral values?\n"));
      msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
      msgBox.setDefaultButton(QMessageBox::Ok);
      if ( msgBox.exec()==QMessageBox::Cancel ) {
        DoOpen = 0;
      }
    }
    if ( DoOpen==1 ) {
      CB_OpenSettingsFile(Settings->GetString("PresetDirectory") + "/Neutral_absolute.pts");
    }
  } else if (mode == ptResetMode_User) { // reset to startup settings
    if ( Settings->GetInt("ResetSettingsConfirmation")==1 ) {
      ptMessageBox msgBox;
      msgBox.setIcon(QMessageBox::Question);
      msgBox.setWindowTitle(QObject::tr("Reset?"));
      msgBox.setText(QObject::tr("Reset to start up settings?\n"));
      msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
      msgBox.setDefaultButton(QMessageBox::Ok);
      if ( msgBox.exec()==QMessageBox::Cancel ) {
        DoOpen = 0;
      }
    }
    if ( DoOpen==1 ) {
      CB_OpenSettingsFile(Settings->GetString("StartupSettingsFile"));
    }
  } else if (mode == ptResetMode_OpenPreset) { // open preset file
    CB_OpenPresetFileButton();
  } else { // open settings file
    CB_OpenSettingsFileButton();
  }
}

void CB_ResetButton() {
  ResetButtonHandler(Settings->GetInt("ResetButtonMode"));
}

void CB_SpecialPreviewChoice(const QVariant Choice) {
  Settings->SetValue("SpecialPreview",Choice);
  Update(ptProcessorPhase_Preview);
}

void CB_GimpExecCommandButton() {

  QString GimpExecCommandString = QFileDialog::getOpenFileName(NULL,
    QObject::tr("Get Gimp command"),
    Settings->GetString("UserDirectory"),
    "All files (*.*)");

  if (0 == GimpExecCommandString.size() ) {
    // Canceled just return
    return;
  } else {
    QFileInfo PathInfo(GimpExecCommandString);
    Settings->SetValue("GimpExecCommand",PathInfo.absoluteFilePath());
  }

  // Reflect in gui.
  MainWindow->UpdateSettings();
}

void CB_RememberSettingLevelChoice(const QVariant Choice) {
  Settings->SetValue("RememberSettingLevel",Choice);
  MainWindow->UpdateSettings();
}

void CB_StartupSettingsCheck(const QVariant State) {
  Settings->SetValue("StartupSettings",State);
}

void CB_StartupSettingsResetCheck(const QVariant State) {
  Settings->SetValue("StartupSettingsReset",State);
}

void CB_StartupSettingsButton() {
  QString StartupSettingsString = QFileDialog::getOpenFileName(NULL,
    QObject::tr("Get preset file"),
    Settings->GetString("PresetDirectory"),
    SettingsFilePattern);

  if (0 == StartupSettingsString.size() ) {
    // Canceled just return
    return;
  } else {
    QFileInfo PathInfo(StartupSettingsString);
    Settings->SetValue("PresetDirectory",PathInfo.absolutePath());
    Settings->SetValue("StartupSettingsFile",PathInfo.absoluteFilePath());
  }

  // Reflect in gui.
  MainWindow->UpdateSettings();
}

void CB_StartupUIModeChoice(const QVariant Choice) {
  Settings->SetValue("StartupUIMode",Choice);
}

void CB_StartupPipeSizeChoice(const QVariant Choice) {
  Settings->SetValue("StartupPipeSize", Choice);
}

void CB_CropInitialZoomChoice(const QVariant Choice) {
  Settings->SetValue("CropInitialZoom",Choice);
}

void CB_InputsAddPowerLawCheck(const QVariant State) {
  Settings->SetValue("InputsAddPowerLaw",State);
  if (Settings->GetInt("InputsAddPowerLaw")) {
    Settings->SetValue("InputPowerFactor",2.2);
  } else {
    Settings->SetValue("InputPowerFactor",1.0);
  }
  Update(ptProcessorPhase_RGB);
}

void CB_ToolBoxModeCheck(const QVariant State) {
  Settings->SetValue("ToolBoxMode",State);
  MainWindow->OnToolBoxesEnabledTriggered(Settings->GetInt("ToolBoxMode"));
}

void CB_TabStatusIndicatorInput(const QVariant Value) {
  Settings->SetValue("TabStatusIndicator",Value);
  MainWindow->UpdateSettings();
}

void CB_PreviewTabModeCheck(const QVariant State) {
  Settings->SetValue("PreviewTabMode",State);
  if (Settings->GetInt("PreviewTabMode")) {
    Settings->SetValue("PreviewMode",ptPreviewMode_Tab);
  } else {
    Settings->SetValue("PreviewMode",ptPreviewMode_End);
  }
  Update(ptProcessorPhase_Preview);
}

void CB_BackgroundColorCheck(const QVariant State) {
  Settings->SetValue("BackgroundColor",State);
  SetBackgroundColor(Settings->GetInt("BackgroundColor"));
}

void CB_BackgroundColorButton() {
  QPixmap Pix(80, 14);
  QColor  Color;
  Color.setRed(Settings->GetInt("BackgroundRed"));
  Color.setGreen(Settings->GetInt("BackgroundGreen"));
  Color.setBlue(Settings->GetInt("BackgroundBlue"));
  QColorDialog Dialog(Color,NULL);
  Dialog.setStyle(Theme->systemStyle());
  Dialog.setPalette(Theme->systemPalette());
  Dialog.exec();
  QColor TestColor = Dialog.selectedColor();
  if (TestColor.isValid()) {
    Color = TestColor;
    Settings->SetValue("BackgroundRed",Color.red());
    Settings->SetValue("BackgroundGreen",Color.green());
    Settings->SetValue("BackgroundBlue",Color.blue());
    Pix.fill(Color);
    MainWindow->BackgroundColorButton->setIcon(Pix);
    if (Settings->GetInt("BackgroundColor")){
      SetBackgroundColor(1);
    }
  }
}

void SetBackgroundColor(int SetIt) {
  if (SetIt) {
    QPalette BGPal;
    BGPal.setColor(QPalette::Background, QColor(Settings->GetInt("BackgroundRed"),
                                                Settings->GetInt("BackgroundGreen"),
                                                Settings->GetInt("BackgroundBlue")));
    ViewWindow->setPalette(BGPal);
    MainWindow->ViewStartPage->setPalette(BGPal);
  } else {
    ViewWindow->setPalette(Theme->palette());
    MainWindow->ViewStartPage->setPalette(Theme->palette());
  }
}

void CB_SliderWidthInput(const QVariant Value) {
  Settings->SetValue("SliderWidth",Value);
  if (Settings->GetInt("SliderWidth") == 0)
//  maximum slider width is equal to the toolbar width
    MainWindow->setStyleSheet("ptSlider { max-width: " + QString("%1").arg(10000) + "px; }");
  else
//  maximum slider width is equal to "SliderWidth"
    MainWindow->setStyleSheet("ptSlider { max-width: " + QString("%1").arg(Settings->GetInt("SliderWidth")) + "px; }");
}

void CB_SaveButtonModeChoice(const QVariant Choice) {
  Settings->SetValue("SaveButtonMode",Choice);
  SaveButtonToolTip(Settings->GetInt("SaveButtonMode"));
}

void CB_ResetButtonModeChoice(const QVariant Value) {
  Settings->SetValue("ResetButtonMode",Value);
}

void CB_SearchBarEnableCheck(const QVariant State) {
  Settings->SetValue("SearchBarEnable",State);
  MainWindow->SearchWidget->setVisible(Settings->GetInt("SearchBarEnable"));
}

void SaveButtonToolTip(const short mode) {
  if (mode==ptOutputMode_Full) {
    MainWindow->WritePipeButton->setToolTip(QObject::tr("Save full size image"));
  } else if (mode==ptOutputMode_Pipe) {
    MainWindow->WritePipeButton->setToolTip(QObject::tr("Save current pipe"));
  } else if (mode==ptOutputMode_Jobfile) {
    MainWindow->WritePipeButton->setToolTip(QObject::tr("Save job file"));
  } else if (mode==ptOutputMode_Settingsfile) {
    MainWindow->WritePipeButton->setToolTip(QObject::tr("Save settings file"));
  } else if (mode==ptOutputMode_Batch) {
    MainWindow->WritePipeButton->setToolTip(QObject::tr("Save and send to batch manager"));
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the Camera Tab
// Partim Generic Corrections.
//
////////////////////////////////////////////////////////////////////////////////

void CB_OpenFileButton() {
  CB_MenuFileOpen(0);
}

void CB_OpenSettingsFile(QString SettingsFileName) {
  if (0 == SettingsFileName.size()) return;
  short NextPhase = 1;
  short ReturnValue = GFilterDM->ReadPresetFile(SettingsFileName, NextPhase);
  if (ReturnValue) {
    ptMessageBox::critical(0,"Error","No valid settings file!\n" + SettingsFileName);
    return;
  }
  if (NextPhase == 1) {
    if (Settings->GetInt("HaveImage")==1) {
      CalculateMultipliersFromTemperature();
    }
    if (Settings->GetInt("AutomaticPipeSize") && Settings->ToolIsActive("TabResize")) {
      if (!CalculatePipeSize())
        Update(ptProcessorPhase_Raw,ptProcessorPhase_Load);
    } else {
      Update(ptProcessorPhase_Raw,ptProcessorPhase_Load);
    }
  } else {
    //TODO: check Nextphase against automatic pipesize
    Update(NextPhase);
  }
  SetRatingFromXmp();
  if (Settings->GetInt("LoadTags"))
    SetTagsFromXmp();
}

void CB_OpenSettingsFileButton() {
  if (!ptConfirmRequest::loadConfig(lcmGeneralCfg)) {
    return;
  }

  ptAddUndo();

  QString SettingsFilePattern =
    QObject::tr("All supported files (*.pts *ptj);;Settings files (*.pts);;Job files (*.ptj);;All files (*.*)");
  QString SettingsFileName =
    QFileDialog::getOpenFileName(NULL,
                                 QObject::tr("Open setting file"),
                                 Settings->GetString("RawsDirectory"),
                                 SettingsFilePattern);
  if (0 == SettingsFileName.size()) return;
  CB_OpenSettingsFile(SettingsFileName);
}

void CB_OpenPresetFileButton() {
  ptAddUndo();

  QString SettingsFilePattern =
    QObject::tr("All supported files (*.pts *ptj);;Settings files (*.pts);;Job files (*.ptj);;All files (*.*)");
  QString SettingsFileName =
    QFileDialog::getOpenFileName(NULL,
                                 QObject::tr("Open preset"),
                                 Settings->GetString("PresetDirectory"),
                                 SettingsFilePattern);
  if (0 == SettingsFileName.size()) return;
  CB_OpenSettingsFile(SettingsFileName);
}

void CB_UseThumbnailCheck(const QVariant Check) {
  Settings->SetValue("UseThumbnail", Check);

  if (Settings->GetInt("IsRAW")) {
    Update(ptProcessorPhase_Raw, ptProcessorPhase_Load);
  }
}


void CB_BadPixelsChoice(const QVariant Choice) {
  Settings->SetValue("HaveBadPixels", Choice);
  short Cancelled = 0;
  if (Choice.toInt() == 1) {
    // Request to load one.
    QString BadPixelsFileName =
      QFileDialog::getOpenFileName(NULL,
                                   QObject::tr("Open 'bad pixels' file"),
                                   NULL);
    if (0 == Settings->GetString("BadPixelsFileName").size() ) {
      // A cancel we'll interprate as no bad pixel
      Settings->SetValue("HaveBadPixels",0);
      Cancelled = 1;
    } else {
      Settings->SetValue("HaveBadPixels",2);
      Settings->SetValue("BadPixelsFileName",BadPixelsFileName);
    }
  }
  if (!Cancelled) {
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Load);
    MainWindow->UpdateExifInfo(TheProcessor->m_ExifData);
  }
}

void CB_DarkFrameChoice(const QVariant Choice) {
  Settings->SetValue("HaveDarkFrame",Choice);
  short Cancelled = 0;
  if (Choice.toInt() == 1) {
    // Request to load one.
    QString DarkFrameFileName =
      QFileDialog::getOpenFileName(NULL,
                                   QObject::tr("Open 'dark frame' file"),
                                   NULL);
    if (0 == Settings->GetString("DarkFrameFileName").size() ) {
      // A cancel we'll interprate as no bad pixel
      Settings->SetValue("HaveDarkFrame",0);
      Cancelled = 1;
    } else {
      Settings->SetValue("HaveDarkFrame",2);
      Settings->SetValue("DarkFrameFileName",DarkFrameFileName);
    }
  }
  if (!Cancelled) {
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Load);
    MainWindow->UpdateExifInfo(TheProcessor->m_ExifData);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the Camera Tab
// Partim White Balance.
//
////////////////////////////////////////////////////////////////////////////////
void SelectSpotWBDone(const ptStatus ExitStatus, const QRect SelectionRect);

void CB_WhiteBalanceChoice(const QVariant Choice) {
  Settings->SetValue("WhiteBalance",Choice);
  if (Settings->GetInt("HaveImage") == 0 || (!Settings->useRAWHandling())) return;

  switch (Choice.toInt()) {
    case ptWhiteBalance_Camera :
    case ptWhiteBalance_Auto :
    case ptWhiteBalance_Manual :
      // In fact all of above just translates to settings into
      // DcRaw handled via GuiSettingsToDcRaw. Nothing more to do.
      Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
      break;

    case ptWhiteBalance_Spot: {
      // First : make sure we have Image_AfterDcRaw in the view window.
      // Anything else might have undergone geometric transformations that are
      // impossible to calculate reverse to a spot in dcraw.
      ViewWindow->SaveZoom();
      ViewWindow->ZoomToFit();
      UpdatePreviewImage(TheProcessor->m_Image_AfterDcRaw);

      // Allow to be selected in the view window. And deactivate main.
      BlockTools(btmBlockAll);
      ViewWindow->ShowStatus(QObject::tr("Spot WB"));
      ReportProgress(QObject::tr("Spot WB"));
      ViewWindow->StartSimpleRect(SelectSpotWBDone);
      ViewWindow->setFocus();
      break;
    }

    default :
      // Here we have presets selected from ptWhiteBalances.
      // GuiSettings->m_WhiteBalance should point
      // to the right index into the array.
      Settings->SetValue("RMultiplier",
                         ptWhiteBalances[Choice.toInt()].m_Multipliers[0]);
      Settings->SetValue("GMultiplier",
                         ptWhiteBalances[Choice.toInt()].m_Multipliers[1]);
      Settings->SetValue("BMultiplier",
                         ptWhiteBalances[Choice.toInt()].m_Multipliers[2]);
      Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
  }
}

void SelectSpotWBDone(const ptStatus ExitStatus, const QRect SelectionRect) {
  // Selection is done at this point. Disallow it further and activate main.
  BlockTools(btmUnblock);

  if (ExitStatus == stSuccess) {
    Settings->SetValue("VisualSelectionX", SelectionRect.left());
    Settings->SetValue("VisualSelectionY", SelectionRect.top());
    Settings->SetValue("VisualSelectionWidth", SelectionRect.width());
    Settings->SetValue("VisualSelectionHeight", SelectionRect.height());

    TRACEKEYVALS("Selection X","%d",
                 Settings->GetInt("VisualSelectionX"));
    TRACEKEYVALS("Selection Y","%d",
                 Settings->GetInt("VisualSelectionY"));
    TRACEKEYVALS("Selection W","%d",
                 Settings->GetInt("VisualSelectionWidth"));
    TRACEKEYVALS("Selection H","%d",
                 Settings->GetInt("VisualSelectionHeight"));
  }

  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
  ViewWindow->RestoreZoom();
}

void CB_SpotWBButton() {
  CB_WhiteBalanceChoice(ptWhiteBalance_Spot);
}

void CB_ColorTemperatureInput(const QVariant Value) {
  Settings->SetValue("ColorTemperature",Value);
  Settings->SetValue("WhiteBalance",ptWhiteBalance_Manual);
  if (!TheDcRaw) return;
  CalculateMultipliersFromTemperature();
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_GreenIntensityInput(const QVariant Value) {
  Settings->SetValue("GreenIntensity",Value);
  Settings->SetValue("WhiteBalance",ptWhiteBalance_Manual);
  if (!TheDcRaw) return;
  CalculateMultipliersFromTemperature();
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_RMultiplierInput(const QVariant Value) {
  Settings->SetValue("RMultiplier",Value);
  Settings->SetValue("WhiteBalance",ptWhiteBalance_Manual);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_GMultiplierInput(const QVariant Value) {
  Settings->SetValue("GMultiplier",Value);
  Settings->SetValue("WhiteBalance",ptWhiteBalance_Manual);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_BMultiplierInput(const QVariant Value) {
  Settings->SetValue("BMultiplier",Value);
  Settings->SetValue("WhiteBalance",ptWhiteBalance_Manual);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_ManualBlackPointCheck(const QVariant State) {
  short OldState = Settings->GetInt("ManualBlackPoint");
  Settings->SetValue("ManualBlackPoint",State);
  if (OldState) {
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
  } else {
    MainWindow->UpdateSettings();
  }
}

void CB_BlackPointInput(const QVariant Value) {
  Settings->SetValue("BlackPoint",Value);
  if (Settings->GetInt("ManualBlackPoint")) {
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
  }
}

void CB_ManualWhitePointCheck(const QVariant State) {
  short OldState = Settings->GetInt("ManualWhitePoint");
  Settings->SetValue("ManualWhitePoint",State);
  if (OldState) {
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
  } else {
    MainWindow->UpdateSettings();
  }
}

void CB_WhitePointInput(const QVariant Value) {
  Settings->SetValue("WhitePoint",Value);
  if (Settings->GetInt("ManualWhitePoint")) {
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the Camera Tab
// Partim Demosaicing.
//
////////////////////////////////////////////////////////////////////////////////

void CB_CaCorrectChoice(const QVariant Choice) {
  Settings->SetValue("CaCorrect",Choice);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_CaRedInput(const QVariant Value) {
  Settings->SetValue("CaRed",Value);
  if (Settings->GetInt("CaCorrect")==ptCACorrect_Manual)
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_CaBlueInput(const QVariant Value) {
  Settings->SetValue("CaBlue",Value);
  if (Settings->GetInt("CaCorrect")==ptCACorrect_Manual)
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_GreenEquilInput(const QVariant Value) {
  Settings->SetValue("GreenEquil",Value);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_CfaLineDenoiseInput(const QVariant Value) {
  Settings->SetValue("CfaLineDenoise",Value);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_AdjustMaximumThresholdInput(const QVariant Value) {
  Settings->SetValue("AdjustMaximumThreshold",Value);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_RawDenoiseThresholdInput(const QVariant Value) {
  Settings->SetValue("RawDenoiseThreshold",Value);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_HotpixelReductionInput(const QVariant Value) {
  Settings->SetValue("HotpixelReduction",Value);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_BayerDenoiseChoice(const QVariant Choice) {
  Settings->SetValue("BayerDenoise",Choice);
  if (!Settings->GetInt("PipeSize")) {
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
  }
}

void CB_InterpolationChoice(const QVariant Choice) {
  Settings->SetValue("Interpolation",Choice);
  MainWindow->UpdateSettings();
  if (!Settings->GetInt("PipeSize")) {
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
  }
}

void CB_InterpolationPassesInput(const QVariant Value) {
  Settings->SetValue("InterpolationPasses",Value);
  if (!Settings->GetInt("PipeSize") &&
      (Settings->GetInt("Interpolation")==ptInterpolation_DCB ||
       Settings->GetInt("Interpolation")==ptInterpolation_DCBSoft ||
       Settings->GetInt("Interpolation")==ptInterpolation_DCBSharp)) {
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
  }
}

void CB_MedianPassesInput(const QVariant Value) {
  Settings->SetValue("MedianPasses",Value);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

void CB_ESMedianPassesInput(const QVariant Value) {
  Settings->SetValue("ESMedianPasses",Value);
  if (!Settings->GetInt("PipeSize")) {
    Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
  }
}

void CB_EeciRefineCheck(const QVariant State) {
  Settings->SetValue("EeciRefine",State);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Demosaic);
}

////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the Camera Tab
// Partim Highlight.
//
////////////////////////////////////////////////////////////////////////////////

void CB_ClipModeChoice(const QVariant Choice) {
  Settings->SetValue("ClipMode",Choice);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Highlights);
}

void CB_ClipParameterInput(const QVariant Value) {
  Settings->SetValue("ClipParameter",Value);
  Update(ptProcessorPhase_Raw,ptProcessorPhase_Highlights);
}


////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to Lensfun (Geometry tab)
//
////////////////////////////////////////////////////////////////////////////////

// General tab

void CB_LfunFocalInput(const QVariant Value) {
  Settings->SetValue("LfunFocal", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunApertureInput(const QVariant Value) {
  Settings->SetValue("LfunAperture", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunDistanceInput(const QVariant Value) {
  Settings->SetValue("LfunDistance", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunScaleInput(const QVariant Value) {
  Settings->SetValue("LfunScale", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunAutoScaleCheck(const QVariant State) {
  Settings->SetValue("LfunAutoScale",State);
  MainWindow->LfunScaleWidget->setEnabled(!State.toBool());
  Update(ptProcessorPhase_Geometry);
}


// CA and Vignette tab

void CB_LfunCAModelChoice(const QVariant Choice) {
  Settings->SetValue("LfunCAModel", Choice);
  MainWindow->UpdateLfunCAUI();
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunCALinearKrInput(const QVariant Value) {
  Settings->SetValue("LfunCALinearKr", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunCALinearKbInput(const QVariant Value) {
  Settings->SetValue("LfunCALinearKb", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunCAPoly3VrInput(const QVariant Value) {
  Settings->SetValue("LfunCAPoly3Vr", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunCAPoly3VbInput(const QVariant Value) {
  Settings->SetValue("LfunCAPoly3Vb", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunCAPoly3CrInput(const QVariant Value) {
  Settings->SetValue("LfunCAPoly3Cr", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunCAPoly3CbInput(const QVariant Value) {
  Settings->SetValue("LfunCAPoly3Cb", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunCAPoly3BrInput(const QVariant Value) {
  Settings->SetValue("LfunCAPoly3Br", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunCAPoly3BbInput(const QVariant Value) {
  Settings->SetValue("LfunCAPoly3Bb", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunVignetteModelChoice(const QVariant Choice) {
  Settings->SetValue("LfunVignetteModel", Choice);
  MainWindow->UpdateLfunVignetteUI();
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunVignettePoly6K1Input(const QVariant Value) {
  Settings->SetValue("LfunVignettePoly6K1", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunVignettePoly6K2Input(const QVariant Value) {
  Settings->SetValue("LfunVignettePoly6K2", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunVignettePoly6K3Input(const QVariant Value) {
  Settings->SetValue("LfunVignettePoly6K3", Value);
  Update(ptProcessorPhase_Geometry);
}


// Lens correction tab

void CB_LfunSrcGeoChoice(const QVariant Choice) {
  Settings->SetValue("LfunSrcGeo", Choice);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunTargetGeoChoice(const QVariant Choice) {
  Settings->SetValue("LfunTargetGeo", Choice);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunDistModelChoice(const QVariant Choice) {
  Settings->SetValue("LfunDistModel", Choice);
  MainWindow->UpdateLfunDistUI();
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunDistPoly3K1Input(const QVariant Value) {
  Settings->SetValue("LfunDistPoly3K1", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunDistPoly5K1Input(const QVariant Value) {
  Settings->SetValue("LfunDistPoly5K1", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunDistPoly5K2Input(const QVariant Value) {
  Settings->SetValue("LfunDistPoly5K2", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunDistFov1OmegaInput(const QVariant Value) {
  Settings->SetValue("LfunDistFov1Omega", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunDistPTLensAInput(const QVariant Value) {
  Settings->SetValue("LfunDistPTLensA", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunDistPTLensBInput(const QVariant Value) {
  Settings->SetValue("LfunDistPTLensB", Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_LfunDistPTLensCInput(const QVariant Value) {
  Settings->SetValue("LfunDistPTLensC", Value);
  Update(ptProcessorPhase_Geometry);
}

////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the Geometry Tab
// Partim Defish
//
////////////////////////////////////////////////////////////////////////////////

void CB_DefishCheck(const QVariant State) {
  Settings->SetValue("Defish",State);
  Update(ptProcessorPhase_Geometry);
}

void CB_DefishFocalLengthInput(const QVariant Value) {
  Settings->SetValue("DefishFocalLength", Value);
  if (Settings->GetInt("Defish"))
    Update(ptProcessorPhase_Geometry);
}

void CB_DefishScaleInput(const QVariant Value) {
  Settings->SetValue("DefishScale", Value);
  if (Settings->GetInt("Defish"))
    Update(ptProcessorPhase_Geometry);
}

void CB_DefishAutoScaleCheck(const QVariant State) {
  Settings->SetValue("DefishAutoScale",State);
  MainWindow->DefishScaleWidget->setEnabled(!State.toBool());
  Update(ptProcessorPhase_Geometry);
}

////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the Geometry Tab
// Partim Rotate
//
////////////////////////////////////////////////////////////////////////////////

void CB_RotateLeftButton() {
  double Value = Settings->GetDouble("Rotate");
  Value -= 90.0;
  if (Value < -180.0) Value += 360.0;
  Settings->SetValue("Rotate",Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_RotateRightButton() {
  double Value = Settings->GetDouble("Rotate");
  Value += 90.0;
  if (Value > 180.0) Value -= 360.0;
  Settings->SetValue("Rotate",Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_RotateAngleButton() {
  if (Settings->GetInt("HaveImage")==0) {
    ptMessageBox::information(MainWindow,
      QObject::tr("No selection"),
      QObject::tr("Open an image first."));
    return;
  }

  // Rerun the part of geometry stage before rotate to get correct preview
  // image in TheProcessor->m_Image_AfterGeometry
  TheProcessor->RunGeometry(ptProcessorStopBefore::Rotate);
  ViewWindow->SaveZoom();
  ViewWindow->ZoomToFit();
  UpdatePreviewImage(TheProcessor->m_Image_AfterGeometry); // Calculate in any case.

  // Allow to be selected in the view window. And deactivate main.
  ViewWindow->ShowStatus(QObject::tr("Get angle"));
  ReportProgress(QObject::tr("Get angle"));

  BlockTools(btmBlockAll);
  ViewWindow->StartLine();
  ViewWindow->setFocus();
}

void RotateAngleDetermined(const ptStatus ExitStatus, double RotateAngle) {
  // Selection is done at this point. Disallow it further and activate main.
  BlockTools(btmUnblock);

  if (ExitStatus == stSuccess) {
    if (RotateAngle < -45.0) {
      RotateAngle += 180.0;
    }
    if (fabs(fabs(RotateAngle) - 90.0) < 45.0) {
      RotateAngle -= 90.0;
    }
    Settings->SetValue("Rotate",RotateAngle);
  }

  ViewWindow->RestoreZoom();
  Update(ptProcessorPhase_Geometry);
}

void CB_RotateInput(const QVariant Value) {
  Settings->SetValue("Rotate",Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_PerspectiveFocalLengthInput(const QVariant Value) {
  Settings->SetValue("PerspectiveFocalLength",Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_PerspectiveTiltInput(const QVariant Value) {
  Settings->SetValue("PerspectiveTilt",Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_PerspectiveTurnInput(const QVariant Value) {
  Settings->SetValue("PerspectiveTurn",Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_PerspectiveScaleXInput(const QVariant Value) {
  Settings->SetValue("PerspectiveScaleX",Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_PerspectiveScaleYInput(const QVariant Value) {
  Settings->SetValue("PerspectiveScaleY",Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_GridCheck(const QVariant State) {
  Settings->SetValue("Grid",State);
  ViewWindow->setGrid(Settings->GetInt("Grid"), Settings->GetInt("GridX"), Settings->GetInt("GridY"));
  Update(ptProcessorPhase_Preview);
}

void CB_GridXInput(const QVariant Value) {
  Settings->SetValue("GridX",Value);
  if (Settings->GetInt("Grid")) {
    ViewWindow->setGrid(Settings->GetInt("Grid"), Settings->GetInt("GridX"), Settings->GetInt("GridY"));
    Update(ptProcessorPhase_Preview);
  }
}

void CB_GridYInput(const QVariant Value) {
  Settings->SetValue("GridY",Value);
  if (Settings->GetInt("Grid")) {
    ViewWindow->setGrid(Settings->GetInt("Grid"), Settings->GetInt("GridX"), Settings->GetInt("GridY"));
    Update(ptProcessorPhase_Preview);
  }
}

void CB_FlipModeChoice(const QVariant Value) {
  Settings->SetValue("FlipMode",Value);
  Update(ptProcessorPhase_Geometry);
}

void CB_GeometryBlockCheck(const QVariant State) {
  Settings->SetValue("GeometryBlock",State);
  Update(ptProcessorPhase_RGB);
}

////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the Geometry Tab
// Partim Crop
//
////////////////////////////////////////////////////////////////////////////////

void CB_FixedAspectRatioCheck(const QVariant Check) {
  Settings->SetValue("FixedAspectRatio", Check);
  MainWindow->UpdateCropToolUI();
  if (ViewWindow->interaction() == iaCrop) {
    ViewWindow->crop()->setAspectRatio(Settings->GetInt("FixedAspectRatio"),
                                       Settings->GetInt("AspectRatioW"),
                                       Settings->GetInt("AspectRatioH"));
    ViewWindow->setFocus();
  }
}

void CB_AspectRatioWChoice(const QVariant Value) {
  Settings->SetValue("AspectRatioW",Value);
  if (ViewWindow->interaction() == iaCrop) {
    ViewWindow->crop()->setAspectRatio(Settings->GetInt("FixedAspectRatio"),
                                       Settings->GetInt("AspectRatioW"),
                                       Settings->GetInt("AspectRatioH"));
    ViewWindow->setFocus();
  }
}

void CB_AspectRatioHChoice(const QVariant Value) {
  Settings->SetValue("AspectRatioH",Value);
  if (ViewWindow->interaction() == iaCrop) {
    ViewWindow->crop()->setAspectRatio(Settings->GetInt("FixedAspectRatio"),
                                       Settings->GetInt("AspectRatioW"),
                                       Settings->GetInt("AspectRatioH"));
    ViewWindow->setFocus();
  }
}

void CB_CropExposureInput(const QVariant Value) {
  Settings->SetValue("CropExposure", Value);

  if (ViewWindow->interaction() == iaCrop) {
    if (GBusy) return;
    else GBusy = true;
    UpdatePreviewImage(TheProcessor->m_Image_AfterGeometry); // Calculate in any case.
    GBusy = false;
  }
}

void CB_CropOrientationButton() {
  int w = Settings->GetInt("AspectRatioW");
  int h = Settings->GetInt("AspectRatioH");

  if (w != h) {
    Settings->SetValue("AspectRatioW", h);
    Settings->SetValue("AspectRatioH", w);
    if (ViewWindow->interaction() == iaCrop) {
      ViewWindow->crop()->flipAspectRatio();
      ViewWindow->setFocus();
    }
  }
}

void CB_CropCenterHorButton() {
  if (ViewWindow->interaction() == iaCrop) {
    ViewWindow->crop()->moveToCenter(1, 0);
    ViewWindow->setFocus();
  }
}

void CB_CropCenterVertButton() {
  if (ViewWindow->interaction() == iaCrop) {
    ViewWindow->crop()->moveToCenter(0, 1);
    ViewWindow->setFocus();
  }
}

void CB_CropGuidelinesChoice(const QVariant Choice) {
  Settings->SetValue("CropGuidelines",Choice);
  if (ViewWindow->interaction() == iaCrop) {
    ViewWindow->crop()->setGuidelines(Choice.toInt());
    ViewWindow->setFocus();
  }
}

void CB_LightsOutChoice(const QVariant Choice) {
  // View window takes care of Settings
  if (ViewWindow->interaction() == iaCrop) {
    ViewWindow->crop()->setLightsOut(Choice.toInt());
    ViewWindow->setFocus();
  }
}


// Prepare and start image crop interaction
void CB_MakeCropButton() {
  if (Settings->GetInt("HaveImage")==0) {
    ptMessageBox::information(MainWindow,
      QObject::tr("No crop possible"),
      QObject::tr("Open an image first."));
    return;
  }

  ViewWindow->ShowStatus(QObject::tr("Prepare"));
  ReportProgress(QObject::tr("Prepare for cropping"));

  // Rerun the part of geometry stage before crop to get correct preview
  // image in TheProcessor->m_Image_AfterGeometry
  TheProcessor->RunGeometry(ptProcessorStopBefore::Crop);
  UpdatePreviewImage(TheProcessor->m_Image_AfterGeometry); // Calculate in any case.

  // Allow to be selected in the view window. And deactivate main.
  ViewWindow->ShowStatus(QObject::tr("Crop"));
  ReportProgress(QObject::tr("Crop"));
  BlockTools(btmBlockForCrop);

  switch (Settings->GetInt("CropInitialZoom")) {
    case ptZoomLevel_Current:
      // nothing to do
      break;

    case ptZoomLevel_Fit:
      ViewWindow->ZoomToFit(0);
      break;

    default:
      ViewWindow->ZoomTo((float)Settings->GetInt("CropInitialZoom") / 100);
      break;
  }

  ViewWindow->StartCrop();          // always start the interaction first,
  MainWindow->UpdateCropToolUI();   // *then* update main window
  ViewWindow->setFocus();
}


// After-crop processing and cleanup.
void CleanupAfterCrop(const ptStatus CropStatus, const QRect CropRect) {
  BlockTools(btmUnblock);

  if (CropStatus == stSuccess) {
    // Account for the pipesize factor.
    int XScale = 1<<Settings->GetInt("PipeSize");
    int YScale = 1<<Settings->GetInt("PipeSize");

    if ((CropRect.width() * XScale < 4) || (CropRect.height() * YScale < 4)) {
      ptMessageBox::information(MainWindow,
          QObject::tr("Crop too small"),
          QObject::tr("Crop rectangle needs to be at least 4x4 pixels in size.\nNo crop, try again."));

      if ((Settings->GetInt("CropW") < 4) || (Settings->GetInt("CropH") < 4)) {
        Settings->SetValue("Crop", 0);
        QCheckBox(MainWindow->CropWidget).setCheckState(Qt::Unchecked);
      }

      if(Settings->GetInt("RunMode")==1) {
        // we're in manual mode!
        Update(ptProcessorPhase_Preview);
      }
    } else {
      Settings->SetValue("Crop",1);
      Settings->SetValue("CropX",CropRect.left() * XScale);
      Settings->SetValue("CropY",CropRect.top() * YScale);
      Settings->SetValue("CropW",CropRect.width() * XScale);
      Settings->SetValue("CropH",CropRect.height() * YScale);
    }

    TRACEKEYVALS("PreviewImageW","%d",PreviewImage->m_Width);
    TRACEKEYVALS("PreviewImageH","%d",PreviewImage->m_Height);
    TRACEKEYVALS("XScale","%d",XScale);
    TRACEKEYVALS("YScale","%d",YScale);
    TRACEKEYVALS("CropX","%d",Settings->GetInt("CropX"));
    TRACEKEYVALS("CropY","%d",Settings->GetInt("CropY"));
    TRACEKEYVALS("CropW","%d",Settings->GetInt("CropW"));
    TRACEKEYVALS("CropH","%d",Settings->GetInt("CropH"));


  // Crop aborted, disable crop checkbox when no previous crop present
  } else {
    if ((Settings->GetInt("CropW") < 4) || (Settings->GetInt("CropH") < 4)) {
      Settings->SetValue("Crop", 0);
    }
  }

  Update(ptProcessorPhase_Geometry);
  MainWindow->UpdateCropToolUI();
}


void CB_ConfirmCropButton() {
  ViewWindow->crop()->stop(stSuccess);    // user confirmed crop
}

void CB_CancelCropButton() {
  ViewWindow->crop()->stop(stFailure);    // user cancelled crop
}


// En/disable an existing crop in the pipe
void CB_CropCheck(const QVariant State) {
  Settings->SetValue("Crop",State);

  if (State.toInt() != 0 &&
      (Settings->GetInt("CropW") <= 4 || Settings->GetInt("CropH") <= 4))
  {
    ptMessageBox::information(MainWindow,
        QObject::tr("No previous crop found"),
        QObject::tr("Set a crop rectangle now."));

    CB_MakeCropButton();
    return;
  }

  Update(ptProcessorPhase_Geometry);
}


////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the Geometry Tab
// Partim Liquid rescale
//
////////////////////////////////////////////////////////////////////////////////

void CB_LqrEnergyChoice(const QVariant Choice) {
  Settings->SetValue("LqrEnergy",Choice);
  Update(ptProcessorPhase_Geometry);
}

void CB_LqrScalingChoice(const QVariant Choice) {
  Settings->SetValue("LqrScaling",Choice);
  MainWindow->UpdateLiquidRescaleUI();
  if (Settings->ToolIsActive("TabLiquidRescale"))
    Update(ptProcessorPhase_Geometry);
}

void CB_LqrHorScaleInput(const QVariant Value) {
  Settings->SetValue("LqrHorScale",Value);
  if (Settings->ToolIsActive("TabLiquidRescale"))
    Update(ptProcessorPhase_Geometry);
}

void CB_LqrVertScaleInput(const QVariant Value) {
  Settings->SetValue("LqrVertScale",Value);
  if (Settings->ToolIsActive("TabLiquidRescale"))
    Update(ptProcessorPhase_Geometry);
}

void CB_LqrWidthInput(const QVariant Value) {
  Settings->SetValue("LqrWidth",Value);
  if (Settings->ToolIsActive("TabLiquidRescale"))
    Update(ptProcessorPhase_Geometry);
}

void CB_LqrHeightInput(const QVariant Value) {
  Settings->SetValue("LqrHeight",Value);
  if (Settings->ToolIsActive("TabLiquidRescale"))
    Update(ptProcessorPhase_Geometry);
}

void CB_LqrVertFirstCheck(const QVariant Check) {
  Settings->SetValue("LqrVertFirst",Check);
  if (Settings->ToolIsActive("TabLiquidRescale"))
    Update(ptProcessorPhase_Geometry);
}


////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to the Geometry Tab
// Partim Resize
//
////////////////////////////////////////////////////////////////////////////////

void CalculatePipeSizeHelper(const short Size, const bool NewImage) {
  if (NewImage) {
    Settings->SetValue("PipeSize",Size);
  } else {
    CB_PipeSizeChoice(Size);
  }
}

// returns 1 if pipe was updated
int CalculatePipeSize(const bool NewImage /* = False */) {
  uint16_t InSize = 0;
  if (Settings->GetInt("HaveImage")==0) return 0;
  if (Settings->GetInt("Crop")==0) {
    InSize = MAX(Settings->GetInt("ImageW"),Settings->GetInt("ImageH"));
  } else {
    InSize = MAX(Settings->GetInt("CropW"),Settings->GetInt("CropH"));
  }
  int s = MAX(0,((int)floor(logf((float)InSize*0.9/(float)Settings->GetInt("ResizeScale"))/logf(2))));
  if (s < Settings->GetInt("PipeSize")) {
    if (Settings->GetInt("RunMode") != 1) {// not manual mode
      ImageSaved = 1; // bad hack to check what happens in the next step
      CalculatePipeSizeHelper(s, NewImage);
      if (ImageSaved == 1 && !NewImage) {
        if (Settings->GetInt("PipeSize")==1) {
          ptMessageBox::information(NULL,"Failure!","Could not run on full size!\nWill stay on half size instead!");
          ImageSaved = 0;
          return 0;
        } else {
          ptMessageBox::information(NULL,"Failure!","Could not run on full size!\nWill run on half size instead!");
          CalculatePipeSizeHelper(1, NewImage);
        }
      }
    } else {
      CalculatePipeSizeHelper(s, NewImage);
    }
    return 1;
  }
  return 0;
}

void CB_ResizeCheck(const QVariant Check) {
  Settings->SetValue("Resize",Check);
  if (Settings->GetInt("AutomaticPipeSize") && Settings->GetInt("Resize")) {
    if (!CalculatePipeSize())
      Update(ptProcessorPhase_Geometry);
  } else {
    Update(ptProcessorPhase_Geometry);
  }
  MainWindow->UpdateExifInfo(TheProcessor->m_ExifData);
}

void CB_ResizeDimensionChoice(const QVariant Choice) {
  Settings->SetValue("ResizeDimension",Choice);
  if (Settings->GetInt("Resize")) {
    if (Settings->GetInt("AutomaticPipeSize")) {
      if (!CalculatePipeSize())
        Update(ptProcessorPhase_Geometry);
    } else {
      Update(ptProcessorPhase_Geometry);
    }
    MainWindow->UpdateExifInfo(TheProcessor->m_ExifData);
  } else {
    MainWindow->UpdateSettings();
  }
}

void CB_ResizeScaleInput(const QVariant Value) {
  Settings->SetValue("ResizeScale",Value);
  if (Settings->GetInt("Resize")) {
    if (Settings->GetInt("AutomaticPipeSize")) {
      if (!CalculatePipeSize())
        Update(ptProcessorPhase_Geometry);
    } else {
      Update(ptProcessorPhase_Geometry);
    }
    MainWindow->UpdateExifInfo(TheProcessor->m_ExifData);
  }
}

void CB_ResizeHeightInput(const QVariant Value) {
  Settings->SetValue("ResizeHeight",Value);
  if (Settings->GetInt("Resize") &&
      Settings->GetInt("ResizeDimension") == ptResizeDimension_WidthHeight) {
    if (Settings->GetInt("AutomaticPipeSize")) {
      if (!CalculatePipeSize())
        Update(ptProcessorPhase_Geometry);
    } else {
      Update(ptProcessorPhase_Geometry);
    }
    MainWindow->UpdateExifInfo(TheProcessor->m_ExifData);
  }
}

void CB_ResizeFilterChoice(const QVariant Choice) {
  Settings->SetValue("ResizeFilter",Choice);
  if (Settings->GetInt("Resize")) {
    if (Settings->GetInt("AutomaticPipeSize")) {
      if (!CalculatePipeSize())
        Update(ptProcessorPhase_Geometry);
    } else {
      Update(ptProcessorPhase_Geometry);
    }
    MainWindow->UpdateExifInfo(TheProcessor->m_ExifData);
  }
}

void CB_AutomaticPipeSizeCheck(const QVariant Check) {
  Settings->SetValue("AutomaticPipeSize",Check);
  if (Settings->GetInt("Resize")) {
    if (Settings->GetInt("AutomaticPipeSize")) {
      if (!CalculatePipeSize())
        Update(ptProcessorPhase_Geometry);
    } else {
      Update(ptProcessorPhase_Geometry);
    }
    MainWindow->UpdateExifInfo(TheProcessor->m_ExifData);
  }
}

////////////////////////////////////////////////////////////////////////////////
//
// Callbacks pertaining to Out Tab
//
////////////////////////////////////////////////////////////////////////////////

void CB_SaveFormatChoice(const QVariant Choice) {
  Settings->SetValue("SaveFormat",Choice);
  MainWindow->UpdateSettings();
}

void CB_EraseExifThumbnailCheck(const QVariant State) {
  Settings->SetValue("EraseExifThumbnail",State);
  TheProcessor->ReadExifBuffer();
}

inline void BlockExport(bool Block) {
  if (!Block) QApplication::processEvents();
  MainWindow->ToGimpButton->setEnabled(!Block);
  QApplication::processEvents();
}

inline void BlockSave(bool Block) {
  if (!Block) QApplication::processEvents();
  MainWindow->WritePipeButton->setEnabled(!Block);
  QApplication::processEvents();
}

//==============================================================================

QString GetCurrentImageBaseName() {
  if (Settings->GetInt("HaveImage")==0) return "";

  QStringList InputFileNameList = Settings->GetStringList("InputFileNameList");
  QFileInfo PathInfo(InputFileNameList[0]);
  return PathInfo.dir().absolutePath() + QDir::separator() + PathInfo.baseName() +
      Settings->GetString("OutputFileNameSuffix") + ".";
}

//==============================================================================

void SaveOutput(const short mode) {
  BlockSave(true);

  if (mode==ptOutputMode_Full) {
    CB_MenuFileSaveOutput();

  } else if (mode==ptOutputMode_Pipe) {
    WritePipe();

  } else if (mode==ptOutputMode_Jobfile) {
    GFilterDM->WriteJobFile();

  } else if (mode==ptOutputMode_Settingsfile) {
    GFilterDM->WritePresetFile(GetCurrentImageBaseName());

  } else if (mode==ptOutputMode_Batch) {
    GFilterDM->SendToBatch(GetCurrentImageBaseName());
  }

  BlockSave(false);
}

//==============================================================================

void Export(const short mode) {
  BlockExport(true);
  if (mode==ptExportMode_GimpPipe && Settings->GetInt("ExportToGimp")) {
    GimpExport(1);
    BlockExport(false);
    return;
  } else if (mode==ptExportMode_GimpFull && Settings->GetInt("ExportToGimp")) {
    GimpExport(0);
    BlockExport(false);
    return;
  }

  // regular export without plugin
  // We generate a temporary name.
  QString ImageFileName = QDir::tempPath() + QDir::separator() +
                          "ptTemp" + QString::number(QDateTime::currentMSecsSinceEpoch()) + ".tif";

  short Temp = Settings->GetInt("SaveFormat");
  Settings->SetValue("SaveFormat", ptSaveFormat_TIFF16);

  if (mode==ptExportMode_GimpPipe && !Settings->GetInt("ExportToGimp")) {
    WritePipe(ImageFileName);
  } else if (mode==ptExportMode_GimpFull && !Settings->GetInt("ExportToGimp")) {
    CB_MenuFileSaveOutput(ImageFileName);
  }

  Settings->SetValue("SaveFormat", Temp);

  QString ExportCommand = Settings->GetString("GimpExecCommand");
  QStringList ExportArguments;
  ExportArguments << ImageFileName;
  QProcess* ExportProcess = new QProcess();
  ExportProcess->startDetached(ExportCommand,ExportArguments);
  BlockExport(false);
}

void CB_WritePipeButton() {
  SaveOutput(Settings->GetInt("SaveButtonMode"));
}

void CB_FileMgrUseThumbMaxRowColCheck(const QVariant checked) {
  Settings->SetValue("FileMgrUseThumbMaxRowCol", checked);
  MainWindow->FileMgrThumbMaxRowColWidget->setEnabled(checked.toBool());
}


////////////////////////////////////////////////////////////////////////////////
//
// Callback dispatcher
//
// For refactoring:
// The control switching the active state of the tool needs always to trigger
// a pipe run (also if it just disabled the tool!) The only exception is the
// external blocked state (blocked or hidden tools!).
//
////////////////////////////////////////////////////////////////////////////////

// For simple switches, just get the value to the settings
void Standard_CB_JustSet (const QString ObjectName, const QVariant Value) {
  QString Temp = ObjectName;
  if (Temp.endsWith("Input") || Temp.endsWith("Check")) Temp.chop(5);
  if (Temp.endsWith("Choice")) Temp.chop(6);
  Settings->SetValue(Temp,Value);
}

// Standard form of call backs, value to settings and pipe run if needed
void Standard_CB_SetAndRun (const QString ObjectName, const QVariant Value) {
  QString Key = ObjectName;
  if (Key.endsWith("Input") || Key.endsWith("Check")) Key.chop(5);
  if (Key.endsWith("Choice")) Key.chop(6);

  QWidget* CurrentControl = Settings->GetGuiWidget(Key);
  if (CurrentControl == NULL) assert(!"Widget not found");
  ptGroupBox* CurrentTool = dynamic_cast<ptGroupBox*>(CurrentControl);

  while (CurrentTool == NULL) {
    CurrentControl = CurrentControl->parentWidget();
    CurrentTool = dynamic_cast<ptGroupBox*>(CurrentControl);
  }
  QString ToolName = CurrentTool->objectName();

  // Save previous state for rerun when disabling a filter
  short PreviousActiveState = Settings->ToolIsActive(ToolName);

  Settings->SetValue(Key,Value);

  short NewActiveState = Settings->ToolIsActive(ToolName);
  CurrentTool->SetActive(NewActiveState);

  if (NewActiveState || PreviousActiveState) {
    Update(ToolName);
  } else {
    MainWindow->UpdateSettings();
  }
}

void CB_InputChanged(const QString ObjectName, const QVariant Value) {
  // No CB processing while in startup phase. Too much
  // noise events otherwise.
  if (InStartup) return;

  if (ObjectName != "ZoomInput" &&
      ObjectName != "PipeSizeChoice" &&
      ObjectName != "RunModeCheck" &&
      ObjectName != "SpecialPreviewChoice" &&
      ObjectName != "UseThumbnailCheck")
  {
    ptAddUndo();
  }

  if (ObjectName == "ZoomInput") {
    CB_ZoomInput(Value);

  #define M_Dispatch(Name)\
  } else if (ObjectName == #Name) { CB_ ## Name (Value);

  #define M_JustSetDispatch(Name)\
  } else if (ObjectName == #Name) { Standard_CB_JustSet(#Name,Value);

  #define M_SetAndRunDispatch(Name)\
  } else if (ObjectName == #Name) { Standard_CB_SetAndRun(#Name,Value);

  M_Dispatch(CameraColorChoice)
  M_Dispatch(CameraColorProfileIntentChoice)
  M_Dispatch(CameraColorGammaChoice)

  M_Dispatch(WorkColorChoice)
  M_Dispatch(CMQualityChoice)

  M_Dispatch(PreviewColorProfileIntentChoice)

  M_Dispatch(StyleChoice)
  M_Dispatch(StyleHighLightChoice)

  M_Dispatch(SaveConfirmationCheck)
  M_Dispatch(AutosaveSettingsCheck)
  M_Dispatch(ResetSettingsConfirmationCheck)
  M_Dispatch(FullPipeConfirmationCheck)

  M_Dispatch(WriteBackupSettingsCheck)

  M_JustSetDispatch(FileMgrThumbnailSizeInput)
  M_JustSetDispatch(FileMgrThumbnailPaddingInput)
  M_Dispatch(FileMgrUseThumbMaxRowColCheck)
  M_JustSetDispatch(FileMgrThumbMaxRowColInput)
  M_JustSetDispatch(FileMgrThumbSaveSizeInput)
  M_JustSetDispatch(FileMgrStartupOpenCheck)
  M_JustSetDispatch(FileMgrThumbCacheSizeInput)

  M_JustSetDispatch(BatchMgrAutosaveCheck)
  M_JustSetDispatch(BatchMgrAutosaveFileChoice)
  M_JustSetDispatch(BatchMgrAutoloadCheck)

  M_Dispatch(MemoryTestInput)

  M_Dispatch(ExportToGimpCheck)

  M_Dispatch(StartupSettingsCheck)
  M_Dispatch(StartupSettingsResetCheck)
  M_Dispatch(StartupUIModeChoice)
  M_Dispatch(StartupPipeSizeChoice)
  M_JustSetDispatch(EscToExitCheck)
  M_JustSetDispatch(LoadTagsCheck)

  M_JustSetDispatch(StartupSwitchARCheck)

  M_Dispatch(CropInitialZoomChoice)

  M_Dispatch(RememberSettingLevelChoice)
  M_Dispatch(InputsAddPowerLawCheck)
  M_Dispatch(ToolBoxModeCheck)
  M_Dispatch(TabStatusIndicatorInput)
  M_Dispatch(PreviewTabModeCheck)
  M_Dispatch(BackgroundColorCheck)
  M_Dispatch(SliderWidthInput)
  M_Dispatch(SaveButtonModeChoice)
  M_Dispatch(ResetButtonModeChoice)
  M_Dispatch(SearchBarEnableCheck)

  M_Dispatch(PipeSizeChoice)
  M_Dispatch(RunModeCheck)
  M_Dispatch(SpecialPreviewChoice)

  M_Dispatch(UseThumbnailCheck)

  M_Dispatch(BadPixelsChoice)
  M_Dispatch(DarkFrameChoice)

  M_Dispatch(WhiteBalanceChoice)
  M_Dispatch(ColorTemperatureInput)
  M_Dispatch(GreenIntensityInput)
  M_Dispatch(RMultiplierInput)
  M_Dispatch(GMultiplierInput)
  M_Dispatch(BMultiplierInput)
  M_Dispatch(ManualBlackPointCheck)
  M_Dispatch(BlackPointInput)
  M_Dispatch(ManualWhitePointCheck)
  M_Dispatch(WhitePointInput)
  M_Dispatch(CaCorrectChoice)
  M_Dispatch(CaRedInput)
  M_Dispatch(CaBlueInput)
  M_Dispatch(GreenEquilInput)
  M_Dispatch(CfaLineDenoiseInput)
  M_Dispatch(AdjustMaximumThresholdInput)
  M_Dispatch(RawDenoiseThresholdInput)
  M_Dispatch(HotpixelReductionInput)
  M_Dispatch(BayerDenoiseChoice)
  M_Dispatch(InterpolationChoice)
  M_Dispatch(InterpolationPassesInput)
  M_Dispatch(MedianPassesInput)
  M_Dispatch(ESMedianPassesInput)
  M_Dispatch(EeciRefineCheck)
  M_Dispatch(ClipModeChoice)
  M_Dispatch(ClipParameterInput)

  M_Dispatch(LfunFocalInput)
  M_Dispatch(LfunApertureInput)
  M_Dispatch(LfunDistanceInput)
  M_Dispatch(LfunScaleInput)
  M_Dispatch(LfunAutoScaleCheck)
  M_Dispatch(LfunCAModelChoice)
  M_Dispatch(LfunCALinearKbInput)
  M_Dispatch(LfunCALinearKrInput)
  M_Dispatch(LfunCAPoly3VrInput)
  M_Dispatch(LfunCAPoly3VbInput)
  M_Dispatch(LfunCAPoly3CrInput)
  M_Dispatch(LfunCAPoly3CbInput)
  M_Dispatch(LfunCAPoly3BrInput)
  M_Dispatch(LfunCAPoly3BbInput)
  M_Dispatch(LfunVignetteModelChoice)
  M_Dispatch(LfunVignettePoly6K1Input)
  M_Dispatch(LfunVignettePoly6K2Input)
  M_Dispatch(LfunVignettePoly6K3Input)
  M_Dispatch(LfunSrcGeoChoice)
  M_Dispatch(LfunTargetGeoChoice)
  M_Dispatch(LfunDistModelChoice)
  M_Dispatch(LfunDistPoly3K1Input)
  M_Dispatch(LfunDistPoly5K1Input)
  M_Dispatch(LfunDistPoly5K2Input)
#if LF_VERSION < (3 << 16)
  M_Dispatch(LfunDistFov1OmegaInput)
#endif
  M_Dispatch(LfunDistPTLensAInput)
  M_Dispatch(LfunDistPTLensBInput)
  M_Dispatch(LfunDistPTLensCInput)

  M_Dispatch(DefishCheck)
  M_Dispatch(DefishFocalLengthInput)
  M_Dispatch(DefishScaleInput)
  M_Dispatch(DefishAutoScaleCheck)

  M_Dispatch(RotateInput)
  M_Dispatch(PerspectiveFocalLengthInput)
  M_Dispatch(PerspectiveTiltInput)
  M_Dispatch(PerspectiveTurnInput)
  M_Dispatch(PerspectiveScaleXInput)
  M_Dispatch(PerspectiveScaleYInput)
  M_Dispatch(GridCheck)
  M_Dispatch(GridXInput)
  M_Dispatch(GridYInput)

  M_Dispatch(CropCheck)
  M_Dispatch(CropGuidelinesChoice)
  M_Dispatch(LightsOutChoice)
  M_Dispatch(FixedAspectRatioCheck)
  M_Dispatch(AspectRatioWChoice)
  M_Dispatch(AspectRatioHChoice)
  M_Dispatch(CropExposureInput)

  M_Dispatch(LqrEnergyChoice)
  M_Dispatch(LqrScalingChoice)
  M_Dispatch(LqrHorScaleInput)
  M_Dispatch(LqrVertScaleInput)
  M_Dispatch(LqrWidthInput)
  M_Dispatch(LqrHeightInput)
  M_Dispatch(LqrVertFirstCheck)

  M_Dispatch(ResizeCheck)
  M_Dispatch(ResizeDimensionChoice)
  M_Dispatch(ResizeScaleInput)
  M_Dispatch(ResizeHeightInput)
  M_Dispatch(ResizeFilterChoice)
  M_Dispatch(AutomaticPipeSizeCheck)

  M_Dispatch(FlipModeChoice)

  M_Dispatch(GeometryBlockCheck)

  M_SetAndRunDispatch(OutputGammaCompensationCheck)
  M_SetAndRunDispatch(OutputGammaInput)
  M_SetAndRunDispatch(OutputLinearityInput)

  M_SetAndRunDispatch(WebResizeChoice)
  M_SetAndRunDispatch(WebResizeDimensionChoice)
  M_SetAndRunDispatch(WebResizeBeforeGammaCheck)
  M_SetAndRunDispatch(WebResizeScaleInput)
  M_SetAndRunDispatch(WebResizeFilterChoice)

  M_Dispatch(OutputColorProfileIntentChoice)

  M_JustSetDispatch(SaveQualityInput)
  M_JustSetDispatch(SaveResolutionInput)
  M_Dispatch(SaveFormatChoice)
  M_JustSetDispatch(SaveSamplingChoice)
  M_JustSetDispatch(IncludeExifCheck)
  M_Dispatch(EraseExifThumbnailCheck)
  M_JustSetDispatch(ImageRatingInput)

  } else {
    fprintf(stderr,"(%s,%d) Unexpected ObjectName '%s'\n",
            __FILE__,__LINE__,ObjectName.toLocal8Bit().data());
    assert(0);
  }
}

//==============================================================================

/*!
 * Reads Sidecar file
 * Places data into \var IptcData and \var XmpData
 */
void ReadSidecar(const QString& Sidecar)
{
  if (Sidecar.size() == 0) return;

  if (Exiv2::ImageFactory::getType(Sidecar.toLocal8Bit().data()) == Exiv2::ImageType::none) {
    return;
  }

  Exiv2::Image::AutoPtr hImage = Exiv2::ImageFactory::open(Sidecar.toLocal8Bit().data());

  if (!hImage.get()) {
    return;
  }
  hImage->readMetadata();
  IptcData = hImage->iptcData();
  XmpData = hImage->xmpData();
}

//==============================================================================

/*!
 * Sets ImageRating value, which is taken from Xmp.xmp.Rating
 * Places data into IptcData and XmpData
 */
void SetRatingFromXmp()
{
  Exiv2::XmpData::const_iterator Pos;
  if ((Pos = XmpData.findKey(Exiv2::XmpKey("Xmp.xmp.Rating"))) != XmpData.end()) {
    std::stringstream str;
    str << *Pos;
    int Rating;
    sscanf(str.str().c_str(),"%d",&Rating);
    Settings->SetValue("ImageRating", Rating);
  }
}

//==============================================================================

/*!
 * Sets Tags value, which is taken from Xmp.digiKam.TagsList
 * Places data into IptcData and XmpData
 */
void SetTagsFromXmp()
{
  Exiv2::XmpData::const_iterator Pos;
  if ((Pos = XmpData.findKey(Exiv2::XmpKey("Xmp.digiKam.TagsList"))) != XmpData.end()) {
    std::stringstream str;
    str << *Pos;
    QString tags = str.str().c_str();
    MainWindow->TagsEditWidget->setPlainText(tags);
  }
}

//==============================================================================

/*! Returns the type of an image file.
  \param filename
    path/filename of the image. An empty string is interpreted as the current global
    input file, i.e. InputFileNameList[0]. Note that \c CheckImageType never changes
    the value of InputFileNameList[0].
  \param width
    A pointer to a variable that holds the width of the image. Pass NULL if you
    do not need to know the width. If the width of the file cannot be determined,
    \c CheckImageType sets this variable to \c 0.
  \param height
    Same as \c width, but for the image height.
  \param dcRaw
    An optional pointer to an existing and properly initialized ptDcRaw object that the
    function should use. Note that the file associated with that dcraw gets changed to
    \c filename.
*/
ptImageType CheckImageType(QString filename,
                           uint16_t* width, uint16_t* height,
                           ptDcRaw* dcRaw /*= NULL*/)
{
  ptImageType result = itUndetermined;
  ptDcRaw* LocalDcRaw = NULL;
  bool UseLocalDcRaw = dcRaw == NULL;

  // Setup dcraw
  if (UseLocalDcRaw) {
    LocalDcRaw = new ptDcRaw;
    Settings->ToDcRaw(LocalDcRaw);
  } else {
    LocalDcRaw = dcRaw;
  }

  // Setup file name
  QStringList InputFileNameList = Settings->GetStringList("InputFileNameList");
  if (filename.isEmpty())
    filename = InputFileNameList[0];

  if (filename != InputFileNameList[0]) {
    LocalDcRaw->m_UserSetting_InputFileName = filename;
  }

  if (LocalDcRaw->Identify() == 0) {
    // we have a raw file
    result = itRaw;
    if (Settings->GetInt("UseThumbnail") == 0) {
      if (width != NULL)  *width  = LocalDcRaw->m_Width;
      if (height != NULL) *height = LocalDcRaw->m_Height;
    } else {
      if (width != NULL)  *width  = LocalDcRaw->m_ThumbWidth;
      if (height != NULL) *height = LocalDcRaw->m_ThumbHeight;
    }

  } else {
    // Not a raw image. We use GraphicsMagick to check for valid Bitmaps.
    MagickWand* image = NewMagickWand();
    MagickPingImage(image, filename.toLocal8Bit().data());

    ExceptionType MagickExcept;
    const char* MagickErrorMsg = MagickGetException(image, &MagickExcept);

    if (MagickExcept == UndefinedException) {
      // image could be pinged without problems: we have a bitmap
      result = itBitmap;
      if (width != NULL) *width = MagickGetImageWidth(image);
      if (height != NULL) *height = MagickGetImageHeight(image);

    } else {
      // not a supported image format
      printf("%s", MagickErrorMsg);
      result = itNotSupported;
      if (width != NULL) *width = 0;
      if (height != NULL) *height = 0;
    }

    DestroyMagickWand(image);
  }

  if (UseLocalDcRaw)
    DelAndNull(LocalDcRaw);

  return result;
}

//==============================================================================
QList<QByteArray> ptUndoBuffer;
QList<QByteArray> ptRedoBuffer;
QByteArray ptClipboard;

////////////////////////////////////////////////////////////////////////////////
//
// ptSettingsToQByteArray
//
// Copies Settings to a byte array
//
////////////////////////////////////////////////////////////////////////////////

QByteArray ptSettingsToQByteArray() {
  ptTempFile* tempFile = new ptTempFile();
  GFilterDM->WritePresetFile(tempFile->fileName());
  QFile file(tempFile->fileName());
  file.open(QIODevice::ReadOnly);
  QByteArray arr = file.readAll();
  file.close();
  delete tempFile;
  return arr;
}


////////////////////////////////////////////////////////////////////////////////
//
// ptQByteArrayToSettings
//
// Loads Settings from a byte array
//
////////////////////////////////////////////////////////////////////////////////

void ptQByteArrayToSettings(const QByteArray &arr) {
  ptTempFile* tempFile = new ptTempFile();
  QFile file(tempFile->fileName());
  file.open(QIODevice::WriteOnly);
  file.write(arr);
  file.close();
  CB_OpenSettingsFile(tempFile->fileName());
  delete tempFile;
}


////////////////////////////////////////////////////////////////////////////////
//
// ptAddUndo
//
// Adds current Settings to undo buffer
//
////////////////////////////////////////////////////////////////////////////////

void ptAddUndo() {
  QByteArray arr = ptSettingsToQByteArray();
  ptUndoBuffer << arr;
  ptRedoBuffer.clear();
}

////////////////////////////////////////////////////////////////////////////////
//
// ptMakeUndo
//
// Loads Settings from undo buffer
//
////////////////////////////////////////////////////////////////////////////////

void ptMakeUndo() {
  if (ptUndoBuffer.isEmpty())
    return;

  QByteArray arr = ptSettingsToQByteArray();
  ptRedoBuffer << arr;
  arr = ptUndoBuffer.last();
  ptUndoBuffer.removeLast();

  ptQByteArrayToSettings(arr);
}

////////////////////////////////////////////////////////////////////////////////
//
// ptMakeRedo
//
// Loads Settings from redo buffer
//
////////////////////////////////////////////////////////////////////////////////

void ptMakeRedo() {
  if (ptRedoBuffer.isEmpty())
    return;

  QByteArray arr = ptSettingsToQByteArray();
  ptUndoBuffer << arr;
  arr = ptRedoBuffer.last();
  ptRedoBuffer.removeLast();

  ptQByteArrayToSettings(arr);
}

////////////////////////////////////////////////////////////////////////////////
//
// ptClearUndoRedo
//
// Clears undo-redo buffers
//
////////////////////////////////////////////////////////////////////////////////

void ptClearUndoRedo() {
  ptUndoBuffer.clear();
  ptRedoBuffer.clear();
}

////////////////////////////////////////////////////////////////////////////////
//
// ptMakeFullUndo
//
// Loads Settings from the first slot of undobuffer
//
////////////////////////////////////////////////////////////////////////////////

void ptMakeFullUndo() {
  if (ptUndoBuffer.isEmpty())
    return;

  QByteArray arr = ptUndoBuffer.first();
  ptQByteArrayToSettings(arr);
  ptClearUndoRedo();
}

////////////////////////////////////////////////////////////////////////////////
//
// ptResetSettingsToDefault
//
// Loads initial Settings
//
////////////////////////////////////////////////////////////////////////////////

void ptResetSettingsToDefault() {
  ptAddUndo();
  CB_OpenSettingsFile(Settings->GetString("StartupSettingsFile"));
}

////////////////////////////////////////////////////////////////////////////////
//
// ptCopySettingsToClipboard
//
// Copies Settings to a "clipboard"
//
////////////////////////////////////////////////////////////////////////////////

void ptCopySettingsToClipboard() {
  ptClipboard = ptSettingsToQByteArray();
}

////////////////////////////////////////////////////////////////////////////////
//
// ptPasteSettingsFromClipboard
//
// Pastes Settings from a "clipboard"
//
////////////////////////////////////////////////////////////////////////////////

void ptPasteSettingsFromClipboard() {
  if (ptClipboard.isEmpty())
    return;

  ptAddUndo();
  ptQByteArrayToSettings(ptClipboard);
}
