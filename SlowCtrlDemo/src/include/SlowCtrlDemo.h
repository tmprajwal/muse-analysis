#ifndef __SLOWCTRLDEMO__
#define __SLOWCTRLDEMO__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include <iostream>
#include "slowctrl.h"
#include "muserawtree.h"

#include <TGraph.h>
#include <stdio.h>


class SlowCtrlDemo:public Plugin
{
 private:
  slowctrl::datum * slctrl_fs11_u;
  // slowctrl::datum * slctrl_sipm_c8_l;
  // slowctrl::datum * slctrl_sipm_cl_sum;

  // slowctrl::datum * slctrl_bhc_l[16];
  // slowctrl::datum * slctrl_bhc_r[16];
  // slowctrl::datum * slctrl_bhd_l[16];
  // slowctrl::datum * slctrl_bhd_r[16];

  slowctrl::datum * slctrl_bha_l_sum;
  slowctrl::datum * slctrl_bha_r_sum;
  slowctrl::datum * slctrl_bhb_l_sum;
  slowctrl::datum * slctrl_bhb_r_sum;
  slowctrl::datum * slctrl_bhc_l_sum;
  slowctrl::datum * slctrl_bhc_r_sum;
  slowctrl::datum * slctrl_bhd_l_sum;
  slowctrl::datum * slctrl_bhd_r_sum;

  slowctrl::datum * slctrl_proton;

  // TH1D* bhc_left[16];
  // TH1D* bhc_right[16];
  // TH1D* bhd_left[16];
  // TH1D* bhd_right[16];

  TH1D* bha_left_sum;
  TH1D* bha_right_sum;
  TH1D* bhb_left_sum;
  TH1D* bhb_right_sum;
  TH1D* bhc_left_sum;
  TH1D* bhc_right_sum;
  TH1D* bhd_left_sum;
  TH1D* bhd_right_sum;

  TH1D* proton_scaler;

  // TH1D* rate;

 public:
  SlowCtrlDemo(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~SlowCtrlDemo();
  // add funtions with return value Long_t here:
  
  Long_t defineHistograms();
  Long_t startup();
  Long_t process();
  Long_t finalize();

  Int_t eventnum = 0;

  UInt_t t0_bha_l = 0;
  UInt_t t0_bha_r = 0;
  UInt_t t0_bhb_l = 0;
  UInt_t t0_bhb_r = 0;
  UInt_t t0_bhc_l = 0;
  UInt_t t0_bhc_r = 0;
  UInt_t t0_bhd_l = 0;
  UInt_t t0_bhd_r = 0;

  unsigned long long sum_bha_l = 0;
  unsigned long long sum_bha_r = 0;
  unsigned long long sum_bhb_l = 0;
  unsigned long long sum_bhb_r = 0;
  unsigned long long sum_bhc_l = 0;
  unsigned long long sum_bhc_r = 0;
  unsigned long long sum_bhd_l = 0;
  unsigned long long sum_bhd_r = 0;

  UInt_t t0_proton = 0;
  unsigned long long sum_proton = 0;
  Int_t beamcount = 0;

  // UInt_t initial_time = 0;

  virtual Long_t cmdline(char * cmd);

  ClassDef(SlowCtrlDemo,1);
    };

#endif
