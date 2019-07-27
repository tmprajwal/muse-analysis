#ifndef __SIPM_HV__
#define __SIPM_HV__

#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TF1.h"
#include <TStyle.h>

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include <iostream>
#include "BHtimetree.h" //header where the Beam Hodoscope data structure is defined;
#include "BMtree.h" // header where the Beam MOnitor data structure is defined;
#include "TPaveText.h"
#include "TCanvas.h"

/*We are plotting exacly the same histograms for BH and BM. But the dimmesions of detector
 * class is different for BH and BM.
 * So it would be smarter to define hispograms as a structure and use it in a future*/
typedef struct bm_bh_hist{
  TH1D *BH[4][16]; //4xSiPM hodoscope planes; 16xbar each;
  TH1D *BM_big[4]; //4xBM big Scintilator bars
  TH1D *BM[32];    //32xSiPm in Beam Monitor plane
} bm_bh_hist;

typedef struct tree_struct_t {
   Int_t plane;
   Int_t bar;
   Int_t up_or_down;
   Float_t pedestal;
   Float_t peak;
   Float_t qdc;
   Float_t gain;
   Float_t shift_dV;
} tree_struct_t;


class SiPM_HV:public Plugin
{
 protected:
/*defining Histograms: each of this structure contains histograms for both BM abd BH*/
  bm_bh_hist MQDC_up;
  bm_bh_hist MQDC_down;
  bm_bh_hist MQDC_mean;

  bm_bh_hist Ped_up;
  bm_bh_hist Ped_down;

  bm_bh_hist Peak_up;
  bm_bh_hist Peak_down;


 private:
  BHraw *bhraw;
  BMraw *bmraw;


 public:
  SiPM_HV(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~SiPM_HV();
  // add funtions with return value Long_t here:
  Long_t defineHistograms();
  Long_t startup();
  Long_t process();
  Long_t finalize();
  void  fit_QDC(double qdc_fit_range_max, double *data_output,TH1D *ped,  TH1D *peak, TH1D *qdc );
  void  fill_gain_tree(TTree *tree, char const * name, int plane, int bar, int up_down, double pedestal, double peak);
  void  fill_gain_hist(int plane, double* data_array, int array_size, char const * D_name, char const* UpDown);
  virtual Long_t cmdline(char * cmd);

  ClassDef(SiPM_HV,1);
    };

#endif
