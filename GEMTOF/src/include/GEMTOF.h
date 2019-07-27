#ifndef __GEMTOF__
#define __GEMTOF__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH2D.h"
#include "TH1D.h"
#include <iostream>
#include "BHtimetree.h"
#include "SPStree.h"
#include "BMtree.h"
#include "VETOtree.h"
#include "lumigemtree.h"
#include "teletracktree.h"
#include <muserawtree.h>

class GEMTOF:public Plugin
{
 private:
  BHraw *bhraw;
  SPSraw *spsraw;
  BMraw *bmraw;
  VETOraw *vetoraw;
  TeleTracks      * GEM_Tracks;
    double rfoffset = 10000.;


  double alignment[16] = {4,2.5,-4,-8,-4.5,-6,-6,-3,-6.5,-7.5,-7.5,-7.5,-3,-5,-5,-5};

  double rfcut[6] = {6.5,9.5,9.5,12.5,13,16.5}; // e, pi, mu
  double rfcenter = 8.;

  std::string planenames[4]={"Plane A","Plane B","Plane C","Plane D"};
  double track_mx, track_my;
  double track_x[4],track_y[4],track_z[4];


 public:
  GEMTOF(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~GEMTOF();
    MRTRunInfo *theRunInfo;
  // add funtions with return value Long_t here:
    void TOFBM(BHbar *bh, BMbar *bm, int bmbar, int bhplane, int bhbar);
  double reftimeBH=0;
  double reftimeBM=0;
  double reftimeSPS=0;
  double reftimeVETO=0;
  double trigger=0;
    double rfcor[2][16];

  Long_t startup();
  Long_t process();
  Long_t finalize();
  Long_t defineHistograms();

    TH2D *BM0_QCDvTOF[2][16];
  TH2D *BM1_QCDvTOF[2][16];
  TH2D *BM2_QCDvTOF[2][16];
  TH2D *BM3_QCDvTOF[2][16];

  TH2D *BM0_QCDvRF[2][16];
  TH2D *BM1_QCDvRF[2][16];
  TH2D *BM2_QCDvRF[2][16];
  TH2D *BM3_QCDvRF[2][16];

  virtual Long_t cmdline(char * cmd);

  ClassDef(GEMTOF,1);
    };

#endif
