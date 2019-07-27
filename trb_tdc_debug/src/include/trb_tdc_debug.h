#ifndef __TRB_TDC_DEBUG__
#define __TRB_TDC_DEBUG__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include <iostream>
#include "TRB_TDCtree.h"






class trb_tdc_debug:public Plugin
{
 private:
   TRB_TDC_Board *trb_tdc;
   TTree *tree;


 public:
  trb_tdc_debug(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~trb_tdc_debug();
  // add funtions with return value Long_t here:
  
   Long_t startup();
   Long_t process();
   Long_t finalize();

  virtual Long_t cmdline(char * cmd);

  ClassDef(trb_tdc_debug,1);
    };

#endif
