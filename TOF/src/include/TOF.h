#ifndef __TOF__
#define __TOF__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include <iostream>
#include "TH1D.h"
#include "TH2D.h"
#include <map>
#include <iterator>
#include "BHtimetree.h"
#include "SPStree.h"
#include "BMtree.h"
#include "VETOtree.h"

#include <stdbool.h>
#include <string>

#define RFTIME_1290 19.75
#define RES_1290N 0.0245

class TOF:public Plugin
{
protected:


 private:
  BHraw *bhraw;
  SPSraw *spsraw;
  BMraw *bmraw;
  VETOraw *vetoraw;
  double rfoffset = 10000.;
  std::string planenames[4]={"Plane A","Plane B","Plane C","Plane D"};
  std::string wallnames[4]={"LeftFront","RightFront","LeftBack","RightBack"};


  //double rfcut[6] = {6.5,9.5,9.5,12.5,13,16.5}; // e, pi, mu
  //double alignment[16] = {4,2.5,-4,-8,-4.5,-6,-6,-3,-6.5,-7.5,-7.5,-7.5,-3,-5,-5,-5};

  
public:
  Long_t defineHistograms();
  Long_t startup();
  Long_t process();
  Long_t finalize();
  double reftimeBH=0;
  double reftimeBM=0;
  double reftimeSPS=0;
  double reftimeVETO=0;
  double trigger=0;
  double Pi_count_[4]={0,0,0,0};
  TH1D *RF_pi_plane[4];
  TH1D *RF_mu_plane[4];
  TH1D *RF_e_plane[4];
  
  double rfcor[4][16];

  MRTRunInfo *theRunInfo;
  void TOFBM(BHbar *bh, BMbar *bm, int bmbar, int bhplane, int bhbar);
  void TOFSPS(BHbar *bh, SPSbar *sps, int spswall, int spsbar, int bhplane, int bhbar);
  void TOFResolution(double reftimeSPS,double reftimeBH, BHbar *bh1,BHbar *bh2,BHbar *rf);

  // TH2D *BM_QCDvTOF[4][2][16];
  // TH2D *BM_QCDvRF[4][2][16];
  // TH2D *BH_QCDvTOF[4][2][16];
  // TH2D *BH_QCDvRF[4][2][16];


  TOF(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~TOF();
  
  virtual Long_t cmdline(char * cmd);

  ClassDef(TOF,1);
    };

#endif
