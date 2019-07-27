#include <PID.h>

#include <iostream>
#include <cmath>


PID::PID(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

PID::~PID()
{
};
//Various functions for finding different hits in an event
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

template <typename T>
auto findFirstFalling(T &container)
{
  for (auto hit=container.begin();hit!=container.end();hit++)
  {
      if (!hit->second.rising)
      { //&& !(hit->trailing)
        return hit;
      }
    }
return container.end();
}
template <typename T>
auto findFirstRising(T &container)
{
  for (auto hit=container.begin();hit!=container.end();hit++)
  {
      if (hit->second.rising)
      { //&& !(hit->trailing)
        return hit;
      }
    }
return container.end();
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
        //best=hit;
        return hit;
      }
    }
    //return best;
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


void PID::plotplane(BHbar * plane)
{
  auto trbuprising = findFirstRising(plane->tdc_trb[1]);
  auto trbdownrising = findFirstRising(plane->tdc_trb[0]);


  auto trbup = findSmallestTime(plane->tdc_trb[1]);
  auto trbdown = findSmallestTime(plane->tdc_trb[0]);

  auto trbupany = findSmallestTime(plane->tdc_trb[1]);
  auto trbdownany = findSmallestTime(plane->tdc_trb[0]);
  
  auto trbupfalling = findFirstFalling(plane->tdc_trb[1]);
  auto trbdownfalling = findFirstFalling(plane->tdc_trb[0]);
  
  auto trbrf = findSmallestTime(bhraw->plane[0][14].tdc_trb[0]);
  auto trbrf42 = findSmallestTime(bhraw->plane[1][14].tdc_trb[0]);
  auto SAup = findSmallestTime(bhraw->plane[2][0].tdc_trb[0]);
  auto SAdown = findSmallestTime(bhraw->plane[2][0].tdc_trb[1]);

  auto trigbhup = findSmallestTime(bhraw->plane[0][7].tdc_trb[0]);
  auto trigbhdown = findSmallestTime(bhraw->plane[0][7].tdc_trb[1]);

  if(trigbhup!=bhraw->plane[0][7].tdc_trb[0].end() && trigbhdown!=bhraw->plane[0][7].tdc_trb[1].end())
  {
      if(trbuprising!=plane->tdc_trb[1].end() && trbdownrising!=plane->tdc_trb[0].end())
      {
        if(trbuprising->second.rising && trbdownrising->second.rising)
        {
          H1((double)(trbuprising->second.time-trbdownrising->second.time),"Time difference first hit rising","Time difference first hit rising;time (ns);counts",34280,-100,100);
          H1(fmod((double)(trbuprising->second.time+trbdownrising->second.time)/2-reftimecfd+rfoffset,19.75),"RF Spectrum","RF Spectrum",857,0,19.75);
        }
      }

      if(trbup!=plane->tdc_trb[1].end() && trbdown!=plane->tdc_trb[0].end())
      {
        if(trbup->second.rising && trbdown->second.rising)
        {
          H1((double)(trbup->second.time-trbdown->second.time),"Time difference first hit","Time difference first hit;time (ns);counts",34280,-100,100);
          H2((double)(trbup->second.time-trbdown->second.time),fmod(trbup->second.time-reftime+rfoffset,19.75),"Time Difference first hit V right BH RF","Time Difference first hit right BH;Right - Left (ns);Right RF(ns)",350,-10,15,857,0,19.75);
          H2((double)(trbup->second.time-trbdown->second.time),fmod(trbdown->second.time-reftime+rfoffset,19.75),"Time Difference first hit V left BH RF","Time Difference first hit left BH;Right - Left (ns);Left RF (ns)",350,-10,15,857,0,19.75);
          H2((double)(trbup->second.time-trbdown->second.time),plane->tdc_trb[1].size(),"Time Difference first hit V right BH Multiplicity","Time Difference first hit right BH;Right - Left (ns);Right Multiplicity (counts)",70,-10,15,10,-0.5,9.5);
          H2((double)(trbup->second.time-trbdown->second.time),plane->tdc_trb[0].size(),"Time Difference first hit V left BH Multiplicity","Time Difference first hit left BH;Right - Left (ns);Left Multiplicity (counts)",70,-10,15,10,-0.5,9.5);

        }

      }

  }

  if(trbupany!=plane->tdc_trb[1].end() && trbdownany!=plane->tdc_trb[0].end())
  {
    H1((double)(trbupany->second.time-trbdownany->second.time),"Time difference first hit any edge","Time difference first hit any edge;time (ns);counts",350,-10,10);
    H1((double)(trbupany->second.time+trbdownany->second.time)/2,"Time average first hit any edge","Time average first hit any edge;time (ns);counts",350,-1,-1);    
  }
  
  if(trbupfalling!=plane->tdc_trb[1].end() && trbdownfalling!=plane->tdc_trb[0].end())
    if(!trbupfalling->second.rising && !trbdownfalling->second.rising)
    {
      H1((double)(trbupfalling->second.time-trbdownfalling->second.time),"Time difference first hit falling","Time difference first hit falling;time (ns);counts",34280,-100,100);
    }

int numrisingup = 0;
int numfallingup = 0;
int numrisingdown = 0;
int numfallingdown = 0;
//Time difference of up - down for all BHs in the TRB3, all hits in an event
  for(auto trbdown1:plane->tdc_trb[0])
  {
    if(trbdown1.second.rising)
    {
      numrisingdown++;
      H1(trbdown1.second.time - trigger,"TRB time Left","TRB time Left; time(ns);counts",4096,-1,-1);
    }
    if(!trbdown1.second.rising)
      numfallingdown++;
    for(auto trbup1:plane->tdc_trb[1])
    {
      if(trbup1.second.rising && trbdown1.second.rising)
      {
        H1(trbup1.second.time - trbdown1.second.time,"TRB Time Difference for", "TRB Time Difference; Time (ns); counts",3428,-20,20); //
        H1((trbup1.second.time + trbdown1.second.time)/2,"TRB Time Average for", "TRB Time average; Time (ns); counts",3428,-1,-1); //        
      }
    }
  }
  for(auto trbup1:plane->tdc_trb[1])
  {
    if(trbup1.second.rising)
    {
      numrisingup++;
      H1(trbup1.second.time-trigger,"TRB time right","TRB time Right;time (ns);counts",4096,-1,-1);
    }
    if(!trbup1.second.rising)
      numfallingup++;
  }

  //2D plot of multiplicity of rising and falling events for up and down.
  H2(numrisingup,numrisingdown,"Num Rising up v down","Num Rising up v down",10,-0.5,9.5,10,-0.5,9.5);
  H2(numfallingup,numfallingdown,"Num Falling up v down","Num Falling up v down",10,-0.5,9.5,10,-0.5,9.5);

  //Plot multiplicity for all BHs in both TRB3 and V1290
  H1(plane->tdc_trb[0].size(),"TRB Left Multiplicity","TRB Left Multiplicity",10,-0.5,10);
  H1(plane->tdc_trb[1].size(),"TRB Right Multiplicity","TRB Right Multiplicity",10,-0.5,10);


  //Produces QDC plots
  for(auto hitadc:plane->adc_mqdc[0])
  {
    H1(hitadc,"MQDC Left","MQDC Left",4096,0,4096);
    for(auto hitadc1:plane->adc_mqdc[1])
    {
      H1(sqrt(hitadc*hitadc1),"MQDC Geo Mean","MQDC Geo Mean",4096,0,4096);
      if(trbup!=plane->tdc_trb[1].end() && trbdown!=plane->tdc_trb[0].end())
      {
        H2(fmod((double)(trbup->second.time+trbdown->second.time)/2-reftime+rfoffset,19.75),sqrt(hitadc*hitadc1),"RF V QDC","RF V QDC",857,0,19.75,1000,0,1000);
      }
    }
  }
  //RF Plots
if(trbup->second.rising && trbup!=plane->tdc_trb[1].end())
  H1(fmod((double)(trbup->second.time)-reftime+rfoffset,19.75),"RF plot up","RF plot",857,0,19.75);

if(trbdown->second.rising && trbdown!=plane->tdc_trb[0].end())
{
  H1(fmod((double)(trbdown->second.time)-reftime+rfoffset,19.75),"BH PID in RF","BH PID;RF Time (ns);counts",857,0,19.75);

}

if(trbup->second.rising && trbdown->second.rising && trbup!=plane->tdc_trb[1].end() && trbdown!=plane->tdc_trb[0].end())
  H1(fmod((double)(trbup->second.time + trbdown->second.time)/2-reftime+rfoffset,19.75),"RF plot avg","RF plot",857,0,19.75);

if(trbup->second.rising && trbup!=plane->tdc_trb[1].end())
  H1(trbup->second.time-trigger,"trigger plot up","RF plot",857,-1,-1);

if(trbdown->second.rising && trbdown!=plane->tdc_trb[0].end())
  H1(trbdown->second.time-trigger,"trigger plot down","RF plot",857,-1,-1);


if(trbupany!=plane->tdc_trb[1].end())
  H1(trbupany->second.time-trigger,"trigger plot up any edge","RF plot",857,-1,-1);

if(trbdownany!=plane->tdc_trb[0].end())
  H1(trbdownany->second.time-trigger,"trigger plot down any edge","RF plot",857,-1,-1);

if(!trbupfalling->second.rising && trbupfalling!=plane->tdc_trb[1].end())
  H1(trbupfalling->second.time-trigger,"trigger plot up falling edge","RF plot",857,-1,-1);

if(!trbdownfalling->second.rising && trbdownfalling!=plane->tdc_trb[0].end())
  H1(trbdownfalling->second.time-trigger,"trigger plot down falling edge","RF plot",857,-1,-1);

if(trbup->second.rising && trbdown->second.rising && trbup!=plane->tdc_trb[1].end() && trbdown!=plane->tdc_trb[0].end())
{
  H1((trbup->second.time + trbdown->second.time)/2-trigger,"trigger plot avg","RF plot",857,-1,-1);
  H1((trbup->second.time + trbdown->second.time)/2,"avg plot","avg plot",857,-1,-1);
  if(SAup->second.rising && SAdown->second.rising && SAdown!=bhraw->plane[2][0].tdc_trb[1].end() && SAup!=bhraw->plane[2][0].tdc_trb[0].end())
  {
    double TOF = (SAup->second.time+SAdown->second.time)/2-(trbup->second.time+trbdown->second.time)/2;
    double RF = fmod((trbup->second.time+trbdown->second.time)/2 -reftime+rfoffset,19.75);
    H1(TOF-RF,"TOF - RF","TOF - RF",250,20,90);
    H1(TOF,"TOF Between BH and SA","TOF Between BH and SA",250,30,70);
    H2(RF,TOF,"TOF V RF","TOF V RF",857,0,19.75,250,0,70);
    for(auto hitadc:plane->adc_mqdc[0])
      for(auto hitadc1:plane->adc_mqdc[1])
      {
        H2(TOF,sqrt(hitadc*hitadc1),"QDC V TOF","QDC V TOF; TOF (ns);BH QDC ",250,0,70,500,0,500);
      }
  }
}
if(trbup->second.rising && trbup!=plane->tdc_trb[1].end())
  H1(fmod((double)(trbup->second.time)-trigger+rfoffset,19.75),"trigger plot up fmod","RF plot",857,0,19.75);

  if(trbdown->second.rising && trbdown!=plane->tdc_trb[0].end())
  H1(fmod((double)(trbdown->second.time)-trigger+rfoffset,19.75),"trigger plot down fmod","RF plot",857,0,19.75);

  if(trbup->second.rising && trbdown->second.rising && trbup!=plane->tdc_trb[1].end() && trbdown!=plane->tdc_trb[0].end())
  H1(fmod((double)(trbup->second.time + trbdown->second.time)/2-trigger+rfoffset,19.75),"trigger plot avg fmod","RF plot",857,0,19.75);

  for(auto hitadc:plane->adc_mqdc[1])
  {
    H1(hitadc,"MQDC Right","MQDC Right",4096,0,4096);
  }

  return;

}
Long_t PID::defineHistograms()
{

	return ok;
}
Long_t PID::startup()
{


  bhraw = NULL;
  getBranchObject("BH",(TObject **) &bhraw);
  if (!bhraw) {
    debug(0,"Could not find BH tree in file\n");
  }

  Particle = new ParTypes;
  makeBranch("ParticleTypes",(TObject **) &Particle);


	return ok;
}

Long_t PID::process()
{
  //Plug in RF
  if(bhraw->trb_reftime.size()==0){
    debug(0,"NO REFERENCE TIME FOR TRB BH\n");
    return ok;
  }
  if(bhraw->trig.size()==0){
    debug(0,"NO TRIGGER TIME FOR SIPM\n");
    return ok;
  }

  reftime =(double)findReference(bhraw->trb_reftime)->time;
  trigger = (double)findTrig(bhraw->trig)->time;
  std::string planenames[]={"First","Second","Third","IFP"};
  for (int i=0;i<4;i++)
  {
    cd (planenames[i].c_str());
    for (int j=0;j<bhraw->plane[i].size();j++)  
    {
      cd (TString::Format("Paddle %i",j));
      plotplane(&bhraw->plane[i][j]);
      cd("..");
  	}
  	cd("..");
  }

  Particle->clear();
  auto trbRF0 = findSmallestTime(bhraw->plane[1][0].tdc_trb[0]);
  auto trbRF1 = findSmallestTime(bhraw->plane[1][1].tdc_trb[0]);
  auto trbRF2 = findSmallestTime(bhraw->plane[1][2].tdc_trb[0]);
  auto trbRF3 = findSmallestTime(bhraw->plane[1][3].tdc_trb[0]);
  auto trbRF4 = findSmallestTime(bhraw->plane[1][4].tdc_trb[0]);
  auto trbRF5 = findSmallestTime(bhraw->plane[1][5].tdc_trb[0]);
  auto trbRF6 = findSmallestTime(bhraw->plane[1][6].tdc_trb[0]);
  auto trbRF7 = findSmallestTime(bhraw->plane[1][7].tdc_trb[0]);
  auto trbRF8 = findSmallestTime(bhraw->plane[1][8].tdc_trb[0]);
  auto trbRF9 = findSmallestTime(bhraw->plane[1][9].tdc_trb[0]);
  auto trbRF10 = findSmallestTime(bhraw->plane[1][10].tdc_trb[0]);
  auto trbRF11 = findSmallestTime(bhraw->plane[1][11].tdc_trb[0]);
  auto trbRF12 = findSmallestTime(bhraw->plane[1][12].tdc_trb[0]);

  if(trbRF0->second.rising && trbRF0!=bhraw->plane[1][0].tdc_trb[0].end())
  {

    //11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
    ParType type;
    type.id = 0;
    double time_ = fmod((double)(trbRF0->second.time)-reftime+rfoffset,19.75);
    if(time_ >4 && time_ <7)
      type.id = 11;
    else if(time_ > 10 && time_ < 13)
      type.id = 13;
    else if(time_ > 15 && time_ <19)
      type.id = 211;
    Particle->particles.emplace_back(std::move(type));
  }
  else if(trbRF1->second.rising && trbRF1!=bhraw->plane[1][1].tdc_trb[0].end())
  {

    //11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
    ParType type;
    type.id = 0;
    double time_ = fmod((double)(trbRF1->second.time)-reftime+rfoffset,19.75);
    if(time_ >5 && time_ <8)
      type.id = 11;
    else if(time_ > 10 && time_ < 14)
      type.id = 13;
    else if(time_ > 15 && time_ <19.5)
      type.id = 211;
    Particle->particles.emplace_back(std::move(type));
  }
  else if(trbRF2->second.rising && trbRF2!=bhraw->plane[1][2].tdc_trb[0].end())
  {

    //11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
    ParType type;
    type.id = 0;
    double time_ = fmod((double)(trbRF2->second.time)-reftime+rfoffset,19.75);
    if(time_ <1 || time_ >18)
      type.id = 11;
    else if(time_ > 8 && time_ < 11)
      type.id = 13;
    else if(time_ > 13 && time_ <17)
      type.id = 211;
    Particle->particles.emplace_back(std::move(type));
  }
  else if(trbRF3->second.rising && trbRF3!=bhraw->plane[1][3].tdc_trb[0].end())
  {

    //11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
    ParType type;
    type.id = 0;
    double time_ = fmod((double)(trbRF3->second.time)-reftime+rfoffset,19.75);
    if(time_ >1 && time_ <4)
      type.id = 11;
    else if(time_ > 10 && time_ < 14)
      type.id = 13;
    else if(time_ > 16 && time_ <19)
      type.id = 211;
    Particle->particles.emplace_back(std::move(type));
  }
  else if(trbRF4->second.rising && trbRF4!=bhraw->plane[1][4].tdc_trb[0].end())
  {

  	//11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
  	ParType type;
  	type.id = 0;
  	double time_ = fmod((double)(trbRF4->second.time)-reftime+rfoffset,19.75);
  	if(time_ <1 || time_ >18)
  		type.id = 11;
  	else if(time_ > 8 && time_ < 11)
  		type.id = 13;
  	else if(time_ > 13 && time_ <17)
  		type.id = 211;
  	Particle->particles.emplace_back(std::move(type));
  }
  else if(trbRF5->second.rising && trbRF5!=bhraw->plane[1][5].tdc_trb[0].end())
  {

    //11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
    ParType type;
    type.id = 0;
    double time_ = fmod((double)(trbRF5->second.time)-reftime+rfoffset,19.75);
    if(time_ >1 && time_ <4)
      type.id = 11;
    else if(time_ > 10 && time_ < 14)
      type.id = 13;
    else if(time_ > 16 && time_ <19)
      type.id = 211;
    Particle->particles.emplace_back(std::move(type));
  }
  else if(trbRF6->second.rising && trbRF6!=bhraw->plane[1][6].tdc_trb[0].end())
  {

    //11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
    ParType type;
    type.id = 0;
    double time_ = fmod((double)(trbRF6->second.time)-reftime+rfoffset,19.75);
    if(time_ >2 && time_ <5)
      type.id = 11;
    else if(time_ > 11 && time_ < 15)
      type.id = 13;
    else if(time_ > 16 || time_ <1)
      type.id = 211;
    Particle->particles.emplace_back(std::move(type));
  }
  else if(trbRF7->second.rising && trbRF7!=bhraw->plane[1][7].tdc_trb[0].end())
  {

    //11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
    ParType type;
    type.id = 0;
    double time_ = fmod((double)(trbRF7->second.time)-reftime+rfoffset,19.75);
    if(time_ >1 && time_ <3)
      type.id = 11;
    else if(time_ > 10 && time_ < 13)
      type.id = 13;
    else if(time_ > 15 && time_ <19)
      type.id = 211;
    Particle->particles.emplace_back(std::move(type));
  }
  else if(trbRF8->second.rising && trbRF8!=bhraw->plane[1][8].tdc_trb[0].end())
  {

    //11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
    ParType type;
    type.id = 0;
    double time_ = fmod((double)(trbRF8->second.time)-reftime+rfoffset,19.75);
    if(time_ >1 && time_ <3)
      type.id = 11;
    else if(time_ > 10 && time_ < 13)
      type.id = 13;
    else if(time_ > 15 && time_ <19)
      type.id = 211;
    Particle->particles.emplace_back(std::move(type));
  }
  else if(trbRF9->second.rising && trbRF9!=bhraw->plane[1][9].tdc_trb[0].end())
  {

    //11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
    ParType type;
    type.id = 0;
    double time_ = fmod((double)(trbRF9->second.time)-reftime+rfoffset,19.75);
    if(time_ >3 && time_ <6)
      type.id = 11;
    else if(time_ > 12 && time_ < 16)
      type.id = 13;
    else if(time_ > 17 || time_ <2)
      type.id = 211;
    Particle->particles.emplace_back(std::move(type));
  }
  else if(trbRF10->second.rising && trbRF10!=bhraw->plane[1][10].tdc_trb[0].end())
  {

    //11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
    ParType type;
    type.id = 0;
    double time_ = fmod((double)(trbRF10->second.time)-reftime+rfoffset,19.75);
    if(time_ >1 && time_ <4)
      type.id = 11;
    else if(time_ > 10 && time_ < 13)
      type.id = 13;
    else if(time_ > 15 && time_ <19)
      type.id = 211;
    Particle->particles.emplace_back(std::move(type));
  }
  else if(trbRF11->second.rising && trbRF11!=bhraw->plane[1][11].tdc_trb[0].end())
  {

    //11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
    ParType type;
    type.id = 0;
    double time_ = fmod((double)(trbRF11->second.time)-reftime+rfoffset,19.75);
    if(time_ >5 && time_ <8)
      type.id = 11;
    else if(time_ > 14 && time_ < 18)
      type.id = 13;
    else if(time_ > 0 && time_ <4)
      type.id = 211;
    Particle->particles.emplace_back(std::move(type));
  }
  else if(trbRF12->second.rising && trbRF12!=bhraw->plane[1][12].tdc_trb[0].end())
  {

    //11 is electron, 13 is muon, 211 is pion. add minus sign for anti particles (+ charge)
    ParType type;
    type.id = 0;
    double time_ = fmod((double)(trbRF12->second.time)-reftime+rfoffset,19.75);
    if(time_ >1 && time_ <4)
      type.id = 11;
    else if(time_ > 10 && time_ < 14)
      type.id = 13;
    else if(time_ > 16 && time_ <19)
      type.id = 211;
    Particle->particles.emplace_back(std::move(type));
  }  

	return ok;
}

Long_t PID::finalize()
{

	return ok;
}


Long_t PID::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new PID(in,out,inf_,outf_,p);
}
}


ClassImp(PID);

