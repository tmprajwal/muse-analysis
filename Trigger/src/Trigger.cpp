#include <Trigger.h>

#include<iostream>
#include<cmath>


Trigger::Trigger(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

Trigger::~Trigger()
{
};

template <typename T>
auto findSmallestTime(T &container)
{
  auto best=container.begin();
  for (auto hit=container.begin();hit!=container.end();hit++)
    {
      if (best->second.time > hit->second.time){ //&& !(hit->trailing)
        best=hit;
      }
    }
  return best;
}

///SPECIAL FUNCTION FOR FINDING THE FIRST RISING RF OR TRIG EDGE////
///DONE DUE TO SPECIAL CLASS DEFINITION OF RF IN SIPM///////////////
template <typename T>
auto findReference(T &container)
{
  auto best=container.begin();
  for (auto hit=container.begin();hit!=container.end();hit++)
    {
      if (hit->rising && (best->time< hit->time))
        { //&& !(hit->trailing)
          return hit;
        }
    }
  return container.end();
}


template <typename T>
auto findTrig(T &container)
{
    auto best=container.begin();
  for (auto hit=container.begin();hit!=container.end();hit++)
  {
      if (hit->rising && (best->time< hit->time))
      { //&& !(hit->trailing)
        best=hit;
        //return hit;
      }
  }
  return best;
  //return container.end();
}

template <typename T>
auto findGreatestTime(T &container)
{
  auto best=container.begin();
  for (auto hit=container.begin();hit!=container.end();hit++)
    {
      if (best->second.time< hit->second.time) 
        {
          best=hit;
        }	
    }
  return best;
}


Long_t Trigger::defineHistograms()
{

	return ok;
}

Long_t Trigger::startup()
{
  spsraw = NULL;
  getBranchObject("SPS",(TObject **) &spsraw);
  if (!spsraw) {
    debug(0,"Could not find SPS tree in file\n");
  }

   bhraw = NULL;
  getBranchObject("BH",(TObject **) &bhraw);
  if (!bhraw) {
      debug(0,"Could not find BH tree in file\n");
  }
  bmraw=NULL;
  getBranchObject("BM", (TObject ** ) &bmraw);
  if (!bmraw) {
    debug(0,"Could not find BM tree in file\n");
  }
  vetoraw = NULL;
  getBranchObject("VETO",(TObject **) &vetoraw);
  if(!vetoraw){
    debug(0,"Could not find VETO tree in file\n"); 
  }
	return ok;
}


Long_t Trigger::process()
{



  std::string planenames[]={"Plane A","Plane B","Plane C","Plane D"};



  //SPS Trigger stuff
  cd("SPS");

  H1(spsraw->wall[0][0].tdc_trb[0].size(),"Trigger/Up Trigger Multiplicity","Up Trigger Multiplicity",10,-0.5,9.5);
  H1(spsraw->wall[0][0].tdc_trb[1].size(),"Trigger/Down Trigger Multiplicity","Down Trigger Multiplicity",10,-0.5,9.5);

  auto trigu = findSmallestTime(spsraw->wall[0][0].tdc_trb[0]);
  auto trigd = findSmallestTime(spsraw->wall[0][0].tdc_trb[1]);

  double totalbarup = 0.;
  double totalbardown = 0.;



  if(trigd!=spsraw->wall[0][0].tdc_trb[1].end() && trigd->second.rising)
  {
    for(int i = 0; i<spsraw->wall[1].size(); i++)
    {
    	int totalbackup = 0;
      auto spsup = findSmallestTime(spsraw->wall[1][i].tdc_trb[0]);
      H1(spsraw->wall[1][i].tdc_trb[0].size(),TString::Format("Trigger/up/mult/Front Bar %i Multiplicity if trig",i),TString::Format("Front Bar %i Multiplicity",i),10,-0.5,9.5);
      if(spsup!=spsraw->wall[1][i].tdc_trb[0].end() && spsup->second.rising)
      {
        totalbarup++;
        H1(spsup->second.time-trigd->second.time,TString::Format("Trigger/up/tdiff/Front Bar %i time - trigger",i),TString::Format("Front Bar %i time - trigger",i),1000,-300,100);

        for(int j = 0; j< spsraw->wall[3].size(); j++)
        {
          auto spsbackup = findSmallestTime(spsraw->wall[3][j].tdc_trb[0]);
          if(spsbackup!=spsraw->wall[3][j].tdc_trb[1].end() && spsbackup->second.rising)
          {
	        H1(spsbackup->second.time-trigd->second.time,TString::Format("Trigger/up/tdiff/Back Bar %i time - trigger",j),TString::Format("Back Bar %i time - trigger",j),1000,-300,100);
		  	H1(spsbackup->second.time-spsup->second.time,TString::Format("Trigger/up/tdiff/Back Bar %i time - Front Bar %i",j,i),TString::Format("Back Bar %i time - Front Bar %i",j,i),1000,-300,100);
		  }
		  totalbackup+= spsraw->wall[3][j].tdc_trb[0].size();
          H1(spsraw->wall[3][j].tdc_trb[0].size(),TString::Format("Trigger/up/mult/Back Multiplicity if front bar %i",i),TString::Format("Back Multiplicity if front bar %i",i),10,-0.5,9.5);
          H2(j,spsraw->wall[3][j].tdc_trb[0].size(),TString::Format("Trigger/up/mult/Back Multiplicity v Bar if front bar %i",i),TString::Format("Back Multiplicity v Bar if Front Bar %i Fired;back bar;Multiplicity",i),28,-0.5,27.5,10,-0.5,9.5);
        }
      }
      H1(totalbackup,TString::Format("Trigger/up/mult/Total Back Multiplicity if front bar %i",i),TString::Format("Back Multiplicity if front bar %i",i),10,-0.5,9.5);      
    }
	  H1(totalbarup,"Trigger/up/Front bars fired if trigger","Front bars fired if trigger",10,-0.5,9.5);
  }

  if(trigu!=spsraw->wall[0][0].tdc_trb[0].end() && trigu->second.rising)
  {
    for(int i = 0; i<spsraw->wall[0].size(); i++)
    {
    	int totalbackdown = 0;
      auto spsup = findSmallestTime(spsraw->wall[1][i].tdc_trb[1]);
      H1(spsraw->wall[1][i].tdc_trb[1].size(),TString::Format("Trigger/down/mult/Front Bar %i Multiplicity if trig",i),TString::Format("Front Bar %i Multiplicity",i),10,-0.5,9.5);
      if(spsup!=spsraw->wall[1][i].tdc_trb[1].end() && spsup->second.rising)
      {
        totalbardown++;
      	H1(spsup->second.time-trigu->second.time,TString::Format("Trigger/down/tdiff/Front Bar %i time - trigger",i),TString::Format("Front Bar %i time - trigger",i),1000,-300,100);
        for(int j = 0; j< spsraw->wall[3].size(); j++)
        {
          auto spsbackup = findSmallestTime(spsraw->wall[3][j].tdc_trb[1]);
          if(spsbackup!=spsraw->wall[3][j].tdc_trb[1].end() && spsbackup->second.rising)
		  {
          	H1(spsbackup->second.time-trigu->second.time,TString::Format("Trigger/down/tdiff/Back Bar %i time - trigger",j),TString::Format("Back Bar %i time - trigger",j),1000,-300,100);
   		  	H1(spsbackup->second.time-spsup->second.time,TString::Format("Trigger/down/tdiff/Back Bar %i time - Front Bar %i",j,i),TString::Format("Back Bar %i time - Front Bar %i",j,i),1000,-300,100);
          }
          totalbackdown+= spsraw->wall[3][j].tdc_trb[1].size();
          H1(spsraw->wall[3][j].tdc_trb[1].size(),TString::Format("Trigger/down/mult/Back Multiplicity if front bar %i",i),TString::Format("Back Multiplicity if front bar %i",i),10,-0.5,9.5);
          H2(j,spsraw->wall[3][j].tdc_trb[1].size(),TString::Format("Trigger/down/mult/Back Multiplicity v Bar if front bar %i",i),TString::Format("Back Multiplicity v Bar if Front Bar %i Fired;back bar;Multiplicity",i),28,-0.5,27.5,10,-0.5,9.5);
        }
      }
      H1(totalbackdown,TString::Format("Trigger/down/mult/Total Back Multiplicity if front bar %i",i),TString::Format("Back Multiplicity if front bar %i",i),10,-0.5,9.5);
    }
	  H1(totalbardown,"Trigger/down/Front bars fired if trigger","Front bars fired if trigger",10,-0.5,9.5);
  }



	if(trigd!=spsraw->wall[0][0].tdc_trb[1].end() && trigd->second.rising)
	{
  		if(trigu!=spsraw->wall[0][0].tdc_trb[0].end() && trigu->second.rising)
		{
			H1(trigd->second.time-trigu->second.time,"Timing difference of triggers","Timing difference of triggers",1000,-50,50);
		}
	}

	cd("..");

	cd("BH");


	//BH Trigger stuff
	//We are comparing wiremaps when BM is the trigger and when BH is trigger
	//When BH is trigger we enforce in software that BM fired in same pattern as if BM was trigger
	//I think the end result is that this a poor trigger test.
	//Ethan
	bool bmfired = false;
	for(int j = 0; j < bmraw->plane[1].size(); j++)
	{
		auto bmup = findSmallestTime(bmraw->plane[1][j].tdc_trb[1]);
		auto bmdown = findSmallestTime(bmraw->plane[1][j].tdc_trb[0]);
		if(bmup!=bmraw->plane[1][j].tdc_trb[1].end() && bmup->second.rising)
			if(bmdown!=bmraw->plane[1][j].tdc_trb[0].end() && bmdown->second.rising)
				bmfired = true;
	}
  if(bhraw->trb_reftime.size()==0){
    debug(0,"NO REFERENCE TIME FOR TRB BH -- SKIPPING EVENT\n");
    return ok;
  }
  if(bhraw->trig.size()==0){
    debug(0,"NO TRIGGER TIME FOR bh -- SKIPPING EVENT\n");
    return ok;
  }
  reftime =(double)findReference(bhraw->trb_reftime)->time;
  trigger = (double)findTrig(bhraw->trig)->time;

	if(bmfired)
	{
		for (int i=0;i<4;i++)
		{
			cd (planenames[i].c_str());
			for (int j=0;j<bhraw->plane[i].size();j++)  
			{

				bool vert = false;
				//determine orientation of BH plane
				std::string sipmname[2];
				if((i % 2) == 1)
					vert = true;

				if(vert){
					sipmname[0]={"Down"};
					sipmname[1]={"Up"};
				}
				if(!vert)
				{
					sipmname[0]={"Left"};
					sipmname[1]={"Right"};
				}

				if(bhraw->plane[i][j].tdc_trb[0].size()>0)
				{
					H1(j,TString::Format("Wiremap %s BH %s",planenames[i].c_str(),sipmname[0].c_str()),TString::Format("Wiremap %s BH %s",planenames[i].c_str(),sipmname[0].c_str()),16,-0.5,15.5);
				}
				if(bhraw->plane[i][j].tdc_trb[1].size()>0)
				{
					H1(j,TString::Format("Wiremap %s BH %s",planenames[i].c_str(),sipmname[1].c_str()),TString::Format("Wiremap %s BH %s",planenames[i].c_str(),sipmname[1].c_str()),16,-0.5,15.5);
				}
				if(bhraw->plane[i][j].tdc_trb[0].size()>0 && bhraw->plane[i][j].tdc_trb[1].size()>0)
				{
					H1(j,TString::Format("Wiremap %s BH both",planenames[i].c_str()),TString::Format("Wiremap %s BH both",planenames[i].c_str()),16,-0.5,15.5);
				    auto up = findSmallestTime(bhraw->plane[i][j].tdc_trb[1]);
        			auto down = findSmallestTime(bhraw->plane[i][j].tdc_trb[0]);
			        if(up->second.rising&&down->second.rising&&up!=bhraw->plane[i][j].tdc_trb[1].end()&&down!=bhraw->plane[i][j].tdc_trb[0].end())
			        {
	        	      double alignment[16] = {3.0,1.25,3.8,0.75,3.5,1.75,1.6,2.45,1.2,.25,.1,0,4.5,3.0,2.5,2.5};
			          double RF = fmod((up->second.time+down->second.time)/2-reftime+rfoffset+alignment[j],19.75);
			          H2(j,RF,"RF v Bar","RF v Bar",16,-0.5,15.5,857,0,19.75);
			          H1(RF,TString::Format("RF Bar %i",j),TString::Format("RF Bar %i",j),857,0,19.75);
			          if(RF >14 && RF < 17)
			            H1(j,TString::Format("Wiremap %s BH both Pi",planenames[i].c_str()),TString::Format("Wiremap %s BH both Pi",planenames[i].c_str()),16,-0.5,15.5);
			          if(RF > 8&& RF < 10.5)
			            H1(j,TString::Format("Wiremap %s BH both e",planenames[i].c_str()),TString::Format("Wiremap %s BH both e",planenames[i].c_str()),16,-0.5,15.5);
			          if(RF > 10.5&& RF < 13)
			            H1(j,TString::Format("Wiremap %s BH both Mu",planenames[i].c_str()),TString::Format("Wiremap %s BH both Mu",planenames[i].c_str()),16,-0.5,15.5);

			        }
				}

			}
			cd("..");
		}
	}
	cd("..");

	return ok;
}

Long_t Trigger::finalize()
{
 std::string planenames[]={"Plane A","Plane B","Plane C","Plane D"};

  //Scales the outer BH paddles down by a factor of two to account for extra wide paddles
  cd("BH");
  for(int j = 0; j < 4; j++)
  {
    cd(planenames[j].c_str());
    bool vert = false;
    std::string sipmname[2];
    if((j % 2) == 1)
      vert = true;
    if(vert){
      sipmname[0]={"Down"};
      sipmname[1]={"Up"};
    }
    if(!vert)
    {
      sipmname[0]={"Left"};
      sipmname[1]={"Right"};
    }
    auto down = dH1(TString::Format("Wiremap %s BH %s",planenames[j].c_str(),sipmname[0].c_str()),TString::Format("Wiremap %s BH %s",planenames[j].c_str(),sipmname[0].c_str()),16,-0.5,15.5);
      auto up = dH1(TString::Format("Wiremap %s BH %s",planenames[j].c_str(),sipmname[1].c_str()),TString::Format("Wiremap %s BH %s",planenames[j].c_str(),sipmname[1].c_str()),16,-0.5,15.5);
      auto both = dH1(TString::Format("Wiremap %s BH both",planenames[j].c_str()),TString::Format("Wiremap %s BH both",planenames[j].c_str()),16,-0.5,15.5);
      auto bothe = dH1(TString::Format("Wiremap %s BH both e",planenames[j].c_str()),TString::Format("Wiremap %s BH both e",planenames[j].c_str()),16,-0.5,15.5);
      auto bothmu = dH1(TString::Format("Wiremap %s BH both Mu",planenames[j].c_str()),TString::Format("Wiremap %s BH both Mu",planenames[j].c_str()),16,-0.5,15.5);
      auto bothpi = dH1(TString::Format("Wiremap %s BH both Pi",planenames[j].c_str()),TString::Format("Wiremap %s BH both Pi",planenames[j].c_str()),16,-0.5,15.5);
    for(int i = 1; i < 17; i++)//starts from 1 because bin 0 is the underflow
    {
      //This divides the five outer BH paddles by 2 because they are twice as wide as the inner paddles
      double valuedown = down->GetBinContent(i);
      double valueup = up->GetBinContent(i);
      double valueboth = both->GetBinContent(i);
      double valuebothe = bothe->GetBinContent(i);
      double valuebothmu = bothmu->GetBinContent(i);
      double valuebothpi = bothpi->GetBinContent(i);

      if(i < 6)
       {
        down->SetBinContent(i,valuedown*0.5);
        up->SetBinContent(i,valueup*0.5);
        both->SetBinContent(i,valueboth*0.5);
        bothe->SetBinContent(i,valuebothe*0.5);
        bothmu->SetBinContent(i,valuebothmu*0.5);
        bothpi->SetBinContent(i,valuebothpi*0.5);
       }
      if(i > 11)
      {
        down->SetBinContent(i,valuedown*0.5);
        up->SetBinContent(i,valueup*0.5);
        both->SetBinContent(i,valueboth*0.5);    
        bothe->SetBinContent(i,valuebothe*0.5);
        bothmu->SetBinContent(i,valuebothmu*0.5);
        bothpi->SetBinContent(i,valuebothpi*0.5);      
      }
    }
    cd("..");
  }
  cd("..");
	return ok;
}


Long_t Trigger::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new Trigger(in,out,inf_,outf_,p);
}
}


ClassImp(Trigger);

