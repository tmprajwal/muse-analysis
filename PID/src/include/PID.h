#ifndef __PID__
#define __PID__

#include "TObject.h"
#include "Plugin.h"
#include "BHtimetree.h"
#include "ParTypestree.h"

#include "TTree.h"
#include <iostream>
#include <map>
#include <iterator>
class PID:public Plugin
{
 protected:
  double reftime;
  double reftimecfd;
  double rfoffset = 1000.;
  double trigger;

 private:
  BHraw *bhraw;
  ParTypes *Particle;


 public:
  PID(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~PID();
  // add funtions with return value Long_t here:
    void plotplane(BHbar *plane);
  Long_t startup();
  Long_t process();
  Long_t finalize();
  Long_t defineHistograms();

  virtual Long_t cmdline(char * cmd);

  ClassDef(PID,1);
    };

#endif
