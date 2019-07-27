#include <SlowCtrlDemo.h>

#include <iostream>
#include <cmath>

#include <TTimeStamp.h>
#include <TH1D.h>


SlowCtrlDemo::SlowCtrlDemo(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

SlowCtrlDemo::~SlowCtrlDemo()
{
};


Long_t SlowCtrlDemo::defineHistograms()
{
  /*
  for (int i = 0; i < 16; i++)
  {
    bhc_left[i]  = dH1(TString::Format("Scaler/BH Plane C/Bar %d/Left",i),TString::Format("Bar %d Left Scaler;Time;Rate",i),7200,0,7200);
    bhc_right[i] = dH1(TString::Format("Scaler/BH Plane C/Bar %d/Right",i),TString::Format("Bar %d Right Scaler;Time;Rate",i),7200,0,7200);
    bhd_left[i]  = dH1(TString::Format("Scaler/BH Plane D/Bar %d/Up",i),TString::Format("Bar %d Up Scaler;Time;Rate",i),7200,0,7200);
    bhd_right[i] = dH1(TString::Format("Scaler/BH Plane D/Bar %d/Down",i),TString::Format("Bar %d Down Scaler;Time;Rate",i),7200,0,7200);
  }
  */

  bha_left_sum  = dH1("Scaler/SUM/BH Plane A Left","Plane A Left Scaler SUM;Time(sec);Value",7200,0,7200);
  bha_right_sum = dH1("Scaler/SUM/BH Plane A Right","Plane A Right Scaler SUM;Time(sec);Value",7200,0,7200);
  bhb_left_sum  = dH1("Scaler/SUM/BH Plane B Up","Plane B Up Scaler SUM;Time(sec);Value",7200,0,7200);
  bhb_right_sum = dH1("Scaler/SUM/BH Plane B Down","Plane B Down Scaler SUM;Time(sec);Value",7200,0,7200);
  bhc_left_sum  = dH1("Scaler/SUM/BH Plane C Left","Plane C Left Scaler SUM;Time;Rate",7200,0,7200);
  bhc_right_sum = dH1("Scaler/SUM/BH Plane C Right","Plane C Right Scaler SUM;Time;Rate",7200,0,7200);
  bhd_left_sum  = dH1("Scaler/SUM/BH Plane D Up","Plane D Up Scaler SUM;Time;Rate",7200,0,7200);
  bhd_right_sum = dH1("Scaler/SUM/BH Plane D Down","Plane D Down Scaler SUM;Time;Rate",7200,0,7200);
  
  
  proton_scaler = dH1("Proton Scaler","Primary Proton Current Scaler;Time(sec);Value",7200,0,7200);

  // rate = dH1("Strip Chart","Strip Chart",7200,0,7200);

  return 0;
}

Long_t SlowCtrlDemo::startup()
{

  auto scmanager= (slowctrl::manager*) getMemoryObject("SlowCtrl Manager");

  if (!scmanager)
    {
      debug(0,"Couldn't find SlowCtrl manager. Plugin included?\n");
      exit(-1);
    }

  // slctrl_fs11_u = scmanager->getLastValidByName("MUSE:BL:FS11-U:POS:AV");
/*
  for (int i = 0; i < 16; i++)
  {
    std::string cl_num = std::to_string(i);
    std::string new_cl_num = std::string(2 - cl_num.length(), '0') + cl_num;  
    std::string cl = "MUSE:SIPM:C:" +new_cl_num+":L:SCL";

    std::string cr_num = std::to_string(i);
    std::string new_cr_num = std::string(2 - cr_num.length(), '0') + cr_num;  
    std::string cr = "MUSE:SIPM:C:" +new_cr_num+":R:SCL";

    std::string dl_num = std::to_string(i);
    std::string new_dl_num = std::string(2 - dl_num.length(), '0') + dl_num;  
    std::string dl = "MUSE:SIPM:D:" +new_dl_num+":T:SCL";

    std::string dr_num = std::to_string(i);
    std::string new_dr_num = std::string(2 - dr_num.length(), '0') + dr_num;  
    std::string dr = "MUSE:SIPM:D:" +new_dr_num+":B:SCL";

    slctrl_bhc_l[i] = scmanager->getLastValidByName(cl.c_str());
    slctrl_bhc_r[i] = scmanager->getLastValidByName(cr.c_str());
    slctrl_bhd_l[i] = scmanager->getLastValidByName(dl.c_str());
    slctrl_bhd_r[i] = scmanager->getLastValidByName(dr.c_str());

    // printf(std::string("MUSE:SIPM:C:%02d:R:SCR",i).c_str());
  }
  //slctrl_sipm_c8_l = scmanager->getLastValidByName("MUSE:SIPM:C:08:L:SCL");
*/
  slctrl_bha_l_sum = scmanager->getLastValidByName("MUSE:SIPM:A:L:SUM");
  slctrl_bha_r_sum = scmanager->getLastValidByName("MUSE:SIPM:A:R:SUM");
  slctrl_bhb_l_sum = scmanager->getLastValidByName("MUSE:SIPM:B:T:SUM");
  slctrl_bhb_r_sum = scmanager->getLastValidByName("MUSE:SIPM:B:B:SUM");
  slctrl_bhc_l_sum = scmanager->getLastValidByName("MUSE:SIPM:C:L:SUM");
  slctrl_bhc_r_sum = scmanager->getLastValidByName("MUSE:SIPM:C:R:SUM");
  slctrl_bhd_l_sum = scmanager->getLastValidByName("MUSE:SIPM:D:T:SUM");
  slctrl_bhd_r_sum = scmanager->getLastValidByName("MUSE:SIPM:D:B:SUM");
  

  slctrl_proton = scmanager->getLastValidByName("MUSE:TRIG:MASTER:CNT:INPUT:03:SCL");

  // initial_time = 0;
 
  return 0;
}

Long_t SlowCtrlDemo::process()
{

  // debug(1,"Value of fs11-u: %g status %i at timestamp %lli\n",slctrl_fs11_u->value,slctrl_fs11_u->status,slctrl_fs11_u->timestamp);

  /// make a plot of value vs time. value is rate

  unsigned long long offset = ((30*365+7)*24+1)*60*60*1e6; // 7 leap year, utc+1 is Swiss time
/*
  for (int i = 0; i < 16; i++)
  {
    // unsigned long long mstime = slctrl_bhc_l[i]->timestamp;

    TTimeStamp time_bhc_l((slctrl_bhc_l[i]->timestamp+offset)/1e6);
    TTimeStamp time_bhc_r((slctrl_bhc_r[i]->timestamp+offset)/1e6);
    TTimeStamp time_bhd_l((slctrl_bhd_l[i]->timestamp+offset)/1e6);
    TTimeStamp time_bhd_r((slctrl_bhd_r[i]->timestamp+offset)/1e6);

    if ( eventnum == 0 )
    {
      t0_bhc_l = time_bhc_l.GetTime();
      t0_bhc_r = time_bhc_r.GetTime();
      t0_bhd_l = time_bhd_l.GetTime();
      t0_bhd_r = time_bhd_r.GetTime();
    }
    else
    {
      bhc_left[i]->SetBinContent((int)(time_bhc_l.GetTime()-t0_bhc_l),slctrl_bhc_l[i]->value);
      bhc_right[i]->SetBinContent((int)(time_bhc_r.GetTime()-t0_bhc_r),slctrl_bhc_r[i]->value);
      bhd_left[i]->SetBinContent((int)(time_bhd_l.GetTime()-t0_bhd_l),slctrl_bhd_l[i]->value);
      bhd_right[i]->SetBinContent((int)(time_bhd_r.GetTime()-t0_bhd_r),slctrl_bhd_r[i]->value);
    }

  }
 */

  TTimeStamp time_bha_l((slctrl_bha_l_sum->timestamp+offset)/1e6);
  TTimeStamp time_bha_r((slctrl_bha_r_sum->timestamp+offset)/1e6);
  TTimeStamp time_bhb_l((slctrl_bhb_l_sum->timestamp+offset)/1e6);
  TTimeStamp time_bhb_r((slctrl_bhb_r_sum->timestamp+offset)/1e6);
  TTimeStamp time_bhc_l((slctrl_bhc_l_sum->timestamp+offset)/1e6);
  TTimeStamp time_bhc_r((slctrl_bhc_r_sum->timestamp+offset)/1e6);
  TTimeStamp time_bhd_l((slctrl_bhd_l_sum->timestamp+offset)/1e6);
  TTimeStamp time_bhd_r((slctrl_bhd_r_sum->timestamp+offset)/1e6);

  if ( eventnum == 0 )
  {
    t0_bha_l = time_bha_l.GetTime();
    t0_bha_r = time_bha_r.GetTime();
    t0_bhb_l = time_bhb_l.GetTime();
    t0_bhb_r = time_bhb_r.GetTime();
    t0_bhc_l = time_bhc_l.GetTime();
    t0_bhc_r = time_bhc_r.GetTime();
    t0_bhd_l = time_bhd_l.GetTime();
    t0_bhd_r = time_bhd_r.GetTime();
  }
  else
  {
    bha_left_sum->SetBinContent((int)(time_bha_l.GetTime()-t0_bha_l),slctrl_bha_l_sum->value);
    bha_right_sum->SetBinContent((int)(time_bha_r.GetTime()-t0_bha_r),slctrl_bha_r_sum->value);
    bhb_left_sum->SetBinContent((int)(time_bhb_l.GetTime()-t0_bhb_l),slctrl_bhb_l_sum->value);
    bhb_right_sum->SetBinContent((int)(time_bhb_r.GetTime()-t0_bhb_r),slctrl_bhb_r_sum->value);
    bhc_left_sum->SetBinContent((int)(time_bhc_l.GetTime()-t0_bhc_l),slctrl_bhc_l_sum->value);
    bhc_right_sum->SetBinContent((int)(time_bhc_r.GetTime()-t0_bhc_r),slctrl_bhc_r_sum->value);
    bhd_left_sum->SetBinContent((int)(time_bhd_l.GetTime()-t0_bhd_l),slctrl_bhd_l_sum->value);
    bhd_right_sum->SetBinContent((int)(time_bhd_r.GetTime()-t0_bhd_r),slctrl_bhd_r_sum->value);
  }

  sum_bha_l += slctrl_bha_l_sum->value;
  sum_bha_r += slctrl_bha_r_sum->value;
  sum_bhb_l += slctrl_bhb_l_sum->value;
  sum_bhb_r += slctrl_bhb_r_sum->value;
  sum_bhc_l += slctrl_bhc_l_sum->value;
  sum_bhc_r += slctrl_bhc_r_sum->value;
  sum_bhd_l += slctrl_bhd_l_sum->value;
  sum_bhd_r += slctrl_bhd_r_sum->value;
  

  // Proton Scaler
  TTimeStamp time_proton((slctrl_proton->timestamp+offset)/1e6);
  if ( eventnum == 0 )
  {
    t0_proton = time_proton.GetTime();
    // std::cout << "time t0 " << time_proton.GetTime() << std::endl;
  }
  else
  {
    proton_scaler->SetBinContent((int)(time_proton.GetTime()-t0_proton+1),slctrl_proton->value);
    // std::cout << "time " << time_proton.GetTime()-t0_proton << std::endl;
    // std::cout << "MUSE:TRIG:MASTER:CNT:INPUT:03:SCL time " << slctrl_proton->timestamp+offset << std::endl;
  }

  if ( slctrl_proton->value != 0 )
  {
    sum_proton += slctrl_proton->value;
    beamcount ++;
  }

  eventnum++;
  
  return 0;

}

Long_t SlowCtrlDemo::finalize()
{
  printf("\n ======= summary ======= \n");
  printf(" BH A L average over all events: %.0f kHz \n",sum_bha_l/(double)eventnum/1000.);
  printf(" BH A R average over all events: %.0f kHz \n",sum_bha_r/(double)eventnum/1000.);
  printf(" BH B L average over all events: %.0f kHz \n",sum_bhb_l/(double)eventnum/1000.);
  printf(" BH B R average over all events: %.0f kHz \n",sum_bhb_r/(double)eventnum/1000.);
  // printf("BH C L average over all events: %g\n",sum_bhc_l/(double)eventnum);
  // printf("BH C R average over all events: %g\n",sum_bhc_r/(double)eventnum);
  // printf("BH D L average over all events: %g\n",sum_bhd_l/(double)eventnum);
  // printf("BH D R average over all events: %g\n",sum_bhd_r/(double)eventnum);
  // printf("\n");
  printf(" Primary Proton Current Rate: %.3fk\n",sum_proton/(double)beamcount/1000.);
  printf("\n ==== end of summary ==== \n");

  return 0;
}

Long_t SlowCtrlDemo::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new SlowCtrlDemo(in,out,inf_,outf_,p);
}
}


ClassImp(SlowCtrlDemo);

