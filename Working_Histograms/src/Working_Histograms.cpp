#include <Working_Histograms.h>

#include<iostream>
#include<cmath>


Working_Histograms::Working_Histograms(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

Working_Histograms::~Working_Histograms()
{
};

Long_t Working_Histograms::startup()
{

  bhraw = NULL;
  getBranchObject("BH",(TObject **) &bhraw);
  if (!bhraw) {
    debug(0,"Could not find SiPM tree in file\n");
  }

  spsraw = NULL;
  getBranchObject("SPS",(TObject **) &spsraw);
  if (!spsraw) {
    debug(0,"Could not find SPS tree in file\n");
  }

  pbglassraw = NULL;
  getBranchObject("PbGlass",(TObject **) &pbglassraw);
  if(!pbglassraw)
    debug(0,"Could not find PbGlass tree in file\n");

  stt_raw = NULL;
  getBranchObject("StrawTube",(TObject **) &stt_raw);
  if(!stt_raw)
    debug(0,"Could not find STT tree in file\n");


  theRunInfo=(MRTRunInfo *) getFileObject("RunInfo");

  return ok;
}

Long_t Working_Histograms::process()
{

  cd("STT");
  int plane_hits[45]={0};
  for(auto hit:stt_raw->straw)
  {
    if(hit.second.rising)
    {
      int side, plane, straw_in_plane;
      STT_internal_to_logic(hit.first,&side,&plane,&straw_in_plane);

      H2(plane,hit.second.time,"STT Total Drift Time V Plane","STT Total Drift Time V Plane;Plane Number;Time (ns)",20,-0.5,19.5,420,-350,120); 
      H2(plane,straw_in_plane,"STT Wiremap V Plane","STT Wiremap V Plane;plane number;straw number",20,-0.5,19.5,18,-0.5,17.5);
      int planestraw=plane+20*side;
      if (planestraw<45)
        plane_hits[planestraw]++;
    }
  }
  for (int planeside=0;planeside<45;planeside++)
  {
     H2(planeside % 20,plane_hits[planeside],"STT Multiplicity V Plane","STT Multiplicty V Plane;Plane Number;Multiplicity",20,-0.5,19.5,10,-0.5,9.5);
  } 
  cd("..");

  cd("SiPM"); 
  for (int i=0;i<4;i++)
  {
    for (int j=0;j<bhraw->plane[i].size();j++)  
    {
      H1(bhraw->plane[i][j].tdc_trb[0].size(),"Total SiPM Left Multiplicity","Total SiPM Left Multiplicity",10,-0.5,10);
      H1(bhraw->plane[i][j].tdc_trb[1].size(),"Total SiPM Right Multiplicity","Total SiPM Right Multiplicity",10,-0.5,10);
      for(auto hitadc:bhraw->plane[i][j].adc_mqdc[0])
      {
        for(auto hitadc1:bhraw->plane[i][j].adc_mqdc[1])
        {
          H1(sqrt(hitadc*hitadc1),"Total SiPM MQDC Geometric Mean","Total SiPM MQDC Geometric Mean",4096,0,4096);
        }
      }
    }
  }
  cd("..");

  cd("SPS");
  for (int i=0;i<4;i++)
  {
    for (int j=0;j<spsraw->wall[i].size();j++)	
    {
      H1(spsraw->wall[i][j].tdc_trb[0].size(),"Total SPS Left Multiplicity","Total SPS Left Multiplicity",10,-0.5,10);
      H1(spsraw->wall[i][j].tdc_trb[1].size(),"Total SPS Right Multiplicity","Total SPS Right Multiplicity",10,-0.5,10);
      for(auto hitadc:spsraw->wall[i][j].adc_mqdc[0])
      {
        for(auto hitadc1:spsraw->wall[i][j].adc_mqdc[1])
        {
          H1(sqrt(hitadc*hitadc1),"Total SPS MQDC Geometric Mean","Total SPS MQDC Geometric Mean",4096,0,4096);
        }
      } 
    }
  }
  cd("..");

  return ok;
}

Long_t Working_Histograms::finalize()
{

  return ok;
}

Long_t Working_Histograms::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new Working_Histograms(in,out,inf_,outf_,p);
}
}


ClassImp(Working_Histograms);


