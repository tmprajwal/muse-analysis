#ifndef __SCIFI_BEAMPROFILE__
#define __SCIFI_BEAMPROFILE__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include <iostream>
#include "musetdc1190tree.h"

class scifi_beamprofile:public Plugin
{
 private:
  MUSETDC1190 *tdc;

 public:
  scifi_beamprofile(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~scifi_beamprofile();
  // add funtions with return value Long_t here:
  
  Long_t startup();
  Long_t process();
  // Long_t finalize()

  virtual Long_t cmdline(char * cmd);

  void TDCChannel(int plane , int fiber , int *TDCa , int *TDCb );

  ClassDef(scifi_beamprofile,1);
    };

#endif
