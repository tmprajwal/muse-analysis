#ifndef __GEM_BH__
#define __GEM_BH__

#include "TObject.h"
#include "TCanvas.h"
#include "Plugin.h"
#include "TTree.h"
#include <Scinthittree.h>
#include "teletracktree.h"


#include <iostream>

class GEM_BH:public Plugin
{
   protected:
  TH1D *xpos[27];
  TH1D *ypos[27];
  TH1D *xzslope[27];
  TH1D *yzslope[27];

 private:
TeleTracks      * GEM_Tracks;
 ScintHits *Hits;
 TH2D *corrzy;
 TH2D *corrzx;
 TH2D *corry_yp;


 public:
  GEM_BH(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~GEM_BH();
  // add funtions with return value Long_t here:
  MRTRunInfo *theRunInfo;
   Long_t startup();
   Long_t process();
   Long_t finalize();
   void id_to_info(int id, int *plane, int *bar, int *side);
  virtual Long_t cmdline(char * cmd);

  ClassDef(GEM_BH,1);
    };

#endif
