#ifndef __MUSETELETRACKER__
#define __MUSETELETRACKER__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH2I.h"
#include <iostream>

#include "lumigemtree.h"
#include "teletracktree.h"
#include "museadctree.h"
#include "musetdc1190tree.h"
#include<muserawtree.h>


class MUSEteleTracker:public Plugin
{
 private:
  LumiGEM    *clusters;
  TeleTracks *teletracks;


// testbeamanalysistree *tdc1190;


 public:
  MUSEteleTracker(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~MUSEteleTracker();
  // add funtions with return value Long_t here:
  
  Long_t defineHistograms();
  Long_t startup();
  Long_t process();
  Long_t finalize();

  Long_t process_ethan();
  Long_t finalize_ethan();
 // Methods for efficiency (in addition to normal functions above!)
  Long_t startup_efficiency();
  Long_t histos_efficiency(); //make histograms
  //  Long_t get_clusters(); // get cluster info
  Long_t process_efficiency(); // filter and fill histograms
  Long_t finalize_efficiency(); //get the average effifiency


  //Some APVs did not work on DS GEM at DEcember, 2017 data. This flag turn DS GEM OFF by taking in account to define the tracks. 
   bool DSOFF;

  virtual Long_t cmdline(char * cmd);

  ClassDef(MUSEteleTracker,1);



  TH1D *multiplicityUS,*multiplicity4TH,*multiplicityMS,*multiplicityDS;
  TH2D *mult2dUS_MS,*mult2dUS_DS,*mult2dMS_DS;
  TH1D *trackmultiplicityUS,*trackmultiplicityMS,*trackmultiplicityDS;

  TH1D *xresidual;
  TH1D *yresidual;
  //  TH1D *USX_eff,*MSX_eff,*DSX_eff;
  //TH1D *USY_eff,*MSY_eff,*DSY_eff;

  TH1D *tracksprojectedgemUSX_onecut,*tracksprojectedgemUSX_any;
  TH1D *tracksprojectedgemMSX_onecut,*tracksprojectedgemMSX_any;
  TH1D *tracksprojectedgemDSX_onecut,*tracksprojectedgemDSX_any;

  TH1D *tracksprojectedgemUSY_onecut,*tracksprojectedgemUSY_any;
  TH1D *tracksprojectedgemMSY_onecut,*tracksprojectedgemMSY_any;
  TH1D *tracksprojectedgemDSY_onecut,*tracksprojectedgemDSY_any;


  TH2D *h2resmapsx,*h2resmapsy,*singleGEMclust;
};

#endif
