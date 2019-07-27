#ifndef __TRIGGER__
#define __TRIGGER__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TH1D.h"
#include <iostream>
#include "BHtimetree.h"
#include "SPStree.h"
#include "BMtree.h"
#include "VETOtree.h"

class Trigger:public Plugin
{
	protected:
  	double reftime = 0.;
  	double trigger = 0.;
 
 	private:
  SPSraw *spsraw;
  BHraw *bhraw;
  BMraw *bmraw;
  VETOraw *vetoraw; 

 public:
  Trigger(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~Trigger();
  // add funtions with return value Long_t here:
  Long_t defineHistograms();
  Long_t startup();
  Long_t process();
  Long_t finalize();
  double rfoffset = 10000.;

  virtual Long_t cmdline(char * cmd);

  ClassDef(Trigger,1);
    };

#endif
