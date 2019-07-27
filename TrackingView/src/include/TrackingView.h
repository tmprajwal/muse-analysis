#ifndef __TRACKINGVIEW__
#define __TRACKINGVIEW__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include <iostream>

#include "teletracktree.h"

class TrackingView:public Plugin
{
 private:
  TeleTracks *teletracks;

 public:
  TrackingView(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~TrackingView();
  // add funtions with return value Long_t here:
  
  Long_t defineHistograms();
  Long_t startup();
  Long_t process();
  Long_t finalize();

  virtual Long_t cmdline(char * cmd);

  ClassDef(TrackingView,1);
};

#endif
