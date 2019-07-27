#include <TrackingView.h>

#include<iostream>
#include<cmath>


TrackingView::TrackingView(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
}

TrackingView::~TrackingView()
{
}

Long_t TrackingView::defineHistograms()
{
  return Plugin::ok;
}

Long_t TrackingView::startup()
{
  teletracks=NULL;
  getBranchObject("teletracks", (TObject**)&teletracks);
  if (teletracks==NULL)
    {
      printf(" Cannot find branch >teletracks< in input ROOT file - bailing out!\n");
      return -1;
    };
  printf(" TrackingView (teletracks) @%p\n", teletracks);

  return Plugin::ok;
}

Long_t TrackingView::process()
{
  

  return Plugin::ok;
}

Long_t TrackingView::finalize()
{
  return Plugin::ok;
}


Long_t TrackingView::cmdline(char *cmd)
{
  //add cmdline handling here

  return 0; // 0 = all ok
}


extern "C"{
  Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
  {
    return (Plugin *) new TrackingView(in,out,inf_,outf_,p);
  }
}


ClassImp(TrackingView);
