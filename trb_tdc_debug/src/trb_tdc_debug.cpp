#include <trb_tdc_debug.h>

#include<iostream>
#include<cmath>


trb_tdc_debug::trb_tdc_debug(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

trb_tdc_debug::~trb_tdc_debug()
{
};

//--------------------------------- Extra Functions --------------------------------------//
// define a global variables that will be used in process routing;
int event_counter;


// This find the smallest time in the TDC channel
template <typename T>
auto findSmallestTime(T &container)
{
  auto best=container.begin();
  for (auto hit=container.begin();hit!=container.end();hit++)
    {
      if (best->second.time > hit->second.time){
           best=hit;
      }
    }
  return best;
}

// This find the greatest time in the TDC channel
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







//---------------------------------  Main routing ----------------------------------------//
// startup() routine is running at the beginning of analysis
Long_t trb_tdc_debug::startup()
{
printf("\t\t\tstartup() routine has started:\n");
  //set event counter to 0:
  event_counter = 0;

  //Accessing root branch with data:
  trb_tdc = NULL;
  getBranchObject("trb_tdc_debug",(TObject **) &trb_tdc);
  if (!trb_tdc) {
    debug(0,"Could not find trb_tdc_debug tree in file\n");
  } else{
    debug(0,"Success: trb_tdc_debug tree is found in file\n");
  }


printf("\t\t\tend of startup() routine!\n");
return ok;
}



//process() routine is running for every event
Long_t trb_tdc_debug::process()
{

 std::string tdc_names[]={"TDC_0","TDC_1","TDC_2","TDC_3"};

 //loop over 4 TDCs on trb3, 48 channels each:
  for(int i=0; i<4; i++){
    cd (tdc_names[i].c_str()); //creating TDC folder and entering into the directory

    //loop over all channels:
    for(int j=0; j<48; j++){

      TRB_TDC_Channel * tdc_ch = &trb_tdc->TDC[i][j];
      //time difference between adjesent channels:
      TRB_TDC_Channel * tdc_ch1 = &trb_tdc->TDC[i][j];
      TRB_TDC_Channel * tdc_ch2 = &trb_tdc->TDC[i][(j+1)%48];


       auto small_time = findSmallestTime(tdc_ch-> trb_tdc);
       auto great_time = findGreatestTime(tdc_ch-> trb_tdc);
       double time_diff =  great_time->second.time - small_time->second.time;

       auto small_time1 = findSmallestTime(tdc_ch1-> trb_tdc);
       auto small_time2 = findSmallestTime(tdc_ch2-> trb_tdc);
       double tdc_ch_res = small_time1->second.time - small_time2->second.time;

       //print debug statement:
       if(tdc_ch->trb_tdc.size() != 0 && event_counter < 1){
         debug(0, "TDC: %d, channel: %d, smallest time: %lf (%d), greatest time: %lf (%d), diff = %lf  \n",
               i, j, small_time->second.time, small_time->second.rising, great_time->second.time, great_time->second.rising, time_diff);
         debug(0, "TDC: %d, j: %d,  (j+1)%48: %d\n", i, j, (j+1)%48 );
       }

      //Plotting Some Histograms for each channel:
      cd (TString::Format("Ch_%i",j));
         for (auto time:tdc_ch-> trb_tdc){
           H2(time.second.time,time.second.rising,
              "Timing vs edge", "Timing vs edge; time (ns); edge (0-trailing, 1-leading)",1000,-200,0,5,-0.5,1.5);
         }
         H1(tdc_ch->trb_tdc.size(),"Multiplicity","Multiplicity;number of hits",10,-0.5,9.5);
         H1(time_diff,"Pulse Width","Pulse Width;time, ns",400, 60, 80);
         H1(tdc_ch_res,"Resolution", TString::Format("Ch_%i-Ch_%i;time, ns", j, (j+1)%48), 1000, -5, 10);
      cd("..");

     }

     cd(".."); //exiting from TDC folder
  }

  event_counter ++; // update event counter after each event;
return ok;
}



// finalize() routine is running at the end of analysis
Long_t trb_tdc_debug::finalize()
{
printf("\t\t\tfinalize() routine has started:\n");



printf("\t\t\tend of finalize() routine!\n");
return ok;
}





Long_t trb_tdc_debug::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new trb_tdc_debug(in,out,inf_,outf_,p);
}
}


ClassImp(trb_tdc_debug);
