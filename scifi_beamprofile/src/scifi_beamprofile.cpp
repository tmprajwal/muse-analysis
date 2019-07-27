#include <scifi_beamprofile.h>

#include<iostream>
#include<cmath>


scifi_beamprofile::scifi_beamprofile(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

scifi_beamprofile::~scifi_beamprofile()
{

};


Long_t scifi_beamprofile::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};



Long_t scifi_beamprofile::startup()
{
  tdc=NULL;
  getBranchObject("v1190", (TObject**)&tdc);
  if (tdc==NULL)
    {
      debug(0,"Could not find v1190 branch\n");
      return -1;
    };

  return ok;
}

Long_t scifi_beamprofile::process()
{
  int TimeOfFlight[128];  //128 TDC channels
  float Position[16];     //16 fibers are connected: 8 Horizontal, 8 Vertical
  for (int i=0;i<128;i++)
    if ((tdc->time[i]>-900) &&(tdc->time[1]>-900))
      {
	H2(i,tdc->time[i],"scifi/raw","raw time",128,-0.5,127.5,512,0,1024*512);
	TimeOfFlight[i] = ((unsigned int) tdc->time[1])-((unsigned int) tdc->time[i]) & 0x7ffff;
	H2(i,TimeOfFlight[i],"scifi/corrected","corrected time",128,-0.5,127.5,1024,0,3073);
      }      

  return ok;
}


void scifi_beamprofile::TDCChannel(int plane , int fiber , int *TDCa , int *TDCb )
{
  int i;
  for (int i  = 17 ; i <= 24 ; i++ )
    if ( fiber == i )
      if ( plane == 1 )	     //Horizontal
	{	*TDCa = i - 1;	*TDCb = i + 7;  }
      else if (plane == 2 )  //Vertical
	{      *TDCa = i + 47;  *TDCb = i + 55; }
}


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new scifi_beamprofile(in,out,inf_,outf_,p);
}
}


ClassImp(scifi_beamprofile);

