#ifndef __LOWLEVEL__
#define __LOWLEVEL__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include <iostream>


class v1190;
class v1290;
class v792;
class TRB3;
class mqdc;
class lowlevel:public Plugin
{
 private:
  v1190 *rawV1190;
  v1290 *rawV1290;
  v792  *rawV792;
  TRB3  *rawTRB3;
  mqdc  *rawMQDC;

 public:
  lowlevel(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~lowlevel();
  // add funtions with return value Long_t here:
  
  Long_t startup();
  Long_t process();
  // Long_t finalize()

  virtual Long_t cmdline(char * cmd);

  ClassDef(lowlevel,1);
    };

#endif
