#ifndef __GEM_STT__
#define __GEM_STT__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TVector2.h"
#include "TH2D.h"
#include "StrawTubehittree.h"
#include "StrawTubetree.h"
#include <math.h>

#include "teletracktree.h"

#include <iostream>

class GEM_STT:public Plugin
{
 private:
 	TH2D *bothXZ, *bothYZ,*circ,*circXZ,*circYZ;
  StrawTubeHits *STT;
TeleTracks      * GEM_Tracks;
 public:
  GEM_STT(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~GEM_STT();
  // add funtions with return value Long_t here:
      MRTRunInfo *theRunInfo;

   Long_t startup();
   Long_t process();
   Long_t finalize();

  virtual Long_t cmdline(char * cmd);

  ClassDef(GEM_STT,1);
    };

#endif
