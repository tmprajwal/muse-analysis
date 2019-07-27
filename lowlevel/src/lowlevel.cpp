#include <lowlevel.h>

#include<iostream>
#include<cmath>


#include "v1190tree.h"
#include "v1290tree.h"
#include "v792tree.h"
#include "trb3tree.h"
#include "mqdctree.h"



lowlevel::lowlevel(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

lowlevel::~lowlevel()
{
};


Long_t lowlevel::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};





Long_t lowlevel::startup()
{
  // try to grab all debug tree structures.

  rawV1190=NULL;
  getOutBranchObject("V1190",(TObject **) &rawV1190);
  if (!rawV1190)
    getBranchObject("V1190",(TObject **) &rawV1190);

  rawV1290=NULL;
  getOutBranchObject("V1290",(TObject **) &rawV1290);
  if (!rawV1290)
    getBranchObject("V1290",(TObject **) &rawV1290);

  rawV792=NULL;
  getOutBranchObject("V792",(TObject **) &rawV792);
  if (!rawV792)
    getBranchObject("V792",(TObject **) &rawV792);

  rawTRB3=NULL;
  getOutBranchObject("TRB3",(TObject **) &rawTRB3);
  if (!rawTRB3)
    getBranchObject("TRB3",(TObject **) &rawTRB3);

  rawMQDC=NULL;
  getOutBranchObject("MQDC",(TObject **) &rawMQDC);
  if(!rawMQDC)
    getBranchObject("MQDC",(TObject **) &rawMQDC);
  
  debug(0,"\n\n*****************\nFound:\n");
  if (rawV1190) debug(0,"V1190\n");
  if (rawV1290) debug(0,"V1290\n");
  if (rawV792) debug(0,"V792\n");
  if (rawTRB3) debug(0,"TRB3\n");
  if (rawMQDC) debug(0,"MQDC\n");


  debug(0,"\n\n*****************\nNot found:\n");
  if (!rawV1190) debug(0,"V1190\n");
  if (!rawV1290) debug(0,"V1290\n");
  if (!rawV792) debug(0,"V792\n");
  if (!rawTRB3) debug(0,"TRB3\n");
  if (!rawMQDC) debug(0,"MQDC\n");

  
  return Plugin::ok;
}


Long_t lowlevel::process()
{

  if (rawV792)
    {
      cd("V792");
      // Let's draw v792 info
      H1(rawV792->extraHits, "extra hits","Extra Hits",100,-0.5,99.5);
      H1(rawV792->eventCounter,"event counter","Event Counter",1000,-0.5,999.5);
      for(auto unmapped:rawV792->unmapped_channels)
	H1(unmapped.second,TString::Format("Channels/%i",unmapped.first),TString::Format ("V792 Channel %i",unmapped.first),4096,-0.5,4095.5);
      cd("..");
    }
  
  if (rawMQDC)
    {
      cd("MQDC");
      // Let's draw MQDC info
      H1(rawMQDC->extraHits, "extra hits","Extra Hits",100,-0.5,99.5);
      H1(rawMQDC->eventCounter,"event counter","Event Counter",1000,-0.5,999.5);
      for(auto unmapped:rawMQDC->unmapped_channels)
	H1(unmapped.second,TString::Format("Channels/%i",unmapped.first),TString::Format ("MQDC Channel %i",unmapped.first),4096,-0.5,4095.5);
      cd("..");
    } 

  
  if (rawV1190)
    {
      cd("V1190");
      // Let's draw v1190 info
      H1(rawV1190->extraHits, "extra hits","Extra Hits",100,-0.5,99.5);
      H1(rawV1190->status,"status","Status",128,-0.5,127.5);
      int last=-1;
      int count=0;
      for(auto unmapped:rawV1190->unmapped_channels)
	{
	  H1(unmapped.second.time,TString::Format("Channels/%i",unmapped.first),TString::Format ("V1190 Channel %i",unmapped.first),1024,-0.5,524287.5);
	  if (last!= (int)unmapped.first)
	    {
	      if (last>=0)
		H1(count,TString::Format("Multiplicity/%i",last),TString::Format ("V1190 Multiplicity Channel %i",last),32,-0.5,31.5);
	      last=unmapped.first;
	      count=0;
	    }
	  count++;

	}
      if (last>=0)
	H1(count,TString::Format("Multiplicity/%i",last),TString::Format ("V1190 Multiplicity Channel %i",last),32,-0.5,31.5);
      cd("..");
	  
    }
  if (rawV1290)
    {
      cd("V1290");
      // Let's draw v1290 info
      H1(rawV1290->extraHits, "extra hits","Extra Hits",100,-0.5,99.5);
      H1(rawV1290->status,"status","Status",128,-0.5,127.5);
      int last=-1;
      int count=0;
      for(auto unmapped:rawV1290->unmapped_channels)
	{
	  H1(unmapped.second.time,TString::Format("Channels/%i",unmapped.first),TString::Format ("V1290 Channel %i",unmapped.first),1024,-0.5,524287.5);
	  if (last!= (int) unmapped.first)
	    {
	      if (last>=0)
		H1(count,TString::Format("Multiplicity/%i",last),TString::Format ("V1290 Multiplicity Channel %i",last),32,-0.5,31.5);
	      last=unmapped.first;
	      count=0;
	    }
	  count++;
	}
      if (last>=0)
	  H1(count,TString::Format("Multiplicity/%i",last),TString::Format ("V1290 Multiplicity Channel %i",last),32,-0.5,31.5);
      cd("..");
    }

  if (rawTRB3)
    {
      cd("TRB3");

      int last=-1;
      int count=0;
      for(auto unmapped:rawTRB3->unmapped_channels)
	{
	  H1(unmapped.second.time,TString::Format("Channels/0x%x",unmapped.first),TString::Format ("TRB3 Channel %i",unmapped.first),1000,-0.5,1e6);
	  if (last!= (int) unmapped.first)
	    {
	      if (last>=0)
		H1(count,TString::Format("Multiplicity/0x%x",last),TString::Format ("TRB3 Multiplicity Channel 0x%x",last),32,-0.5,31.5);
	      last=unmapped.first;
	      count=0;
	    }
	  count++;
	}
      if (last>=0)
	  H1(count,TString::Format("Multiplicity/0x%x",last),TString::Format ("TRB3 Multiplicity Channel 0x%x",last),32,-0.5,31.5);
      cd("..");
    }
  
  
  return Plugin::ok;
}









extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new lowlevel(in,out,inf_,outf_,p);
}
}


ClassImp(lowlevel);

