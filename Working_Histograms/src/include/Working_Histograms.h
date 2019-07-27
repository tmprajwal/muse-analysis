#ifndef __WORKING_HISTOGRAMS__
#define __WORKING_HISTOGRAMS__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include <iostream>
#include "BHtimetree.h"
#include "StrawTubetree.h"
#include "SPStree.h"
#include "PbGlasstree.h"


class Working_Histograms:public Plugin
{
 private:

  BHraw *bhraw;
  SPSraw *spsraw;
  StrawTube * stt_raw;
  PbGlassraw *pbglassraw;

  double rfoffset = 10000.;

 public:
  
  MRTRunInfo *theRunInfo;
  
  Working_Histograms(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~Working_Histograms();
  // add funtions with return value Long_t here:
  
  Long_t startup();
  Long_t process();
  Long_t finalize();

  virtual Long_t cmdline(char * cmd);

  ClassDef(Working_Histograms,1);
    };

#endif
