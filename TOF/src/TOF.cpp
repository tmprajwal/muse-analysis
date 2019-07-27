
#include <TOF.h>

#include <iostream>
#include <cmath>

#include <fstream>
#include <cstdlib>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPad.h>
#include <TStyle.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TChain.h>


TOF::TOF(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

TOF::~TOF()
{
};

template <typename T>
auto findTrig(T &container)
{
	auto best=container.begin();
	for (auto hit=container.begin();hit!=container.end();hit++)
	{
		if (hit->rising && (best->time< hit->time))
      { //&& !(hit->trailing)
      	best=hit;
      }
  }
  return best;
}

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
  auto findSmallestRisingTime(T &container)
  {
  	auto best=container.begin();
  	for (auto hit=container.begin();hit!=container.end();hit++)
  	{
  		if ((best->second.time > hit->second.time)&& (hit->second.rising)) {
  			best=hit;
  		}
  	}
  	return best;
  }

///SPECIAL FUNCTION FOR FINDING THE FIRST RISING RF OR TRIG EDGE////
///DONE DUE TO SPECIAL CLASS DEFINITION OF RF IN bh///////////////
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

void TOF::TOFSPS(BHbar *bh, SPSbar *sps, int spswall, int spsbar, int bhplane, int bhbar) //// To be rewriten 
{
	auto bhleft = findSmallestTime(bh->tdc_trb[0]);
	auto bhright = findSmallestTime(bh->tdc_trb[1]);
	auto spsup = findSmallestTime(sps->tdc_trb[0]);
	auto spsdown = findSmallestTime(sps->tdc_trb[1]);

	if (bhleft!=bh->tdc_trb[0].end() && bhright!=bh->tdc_trb[1].end() && spsup!=sps->tdc_trb[0].end() && spsdown!=sps->tdc_trb[1].end())
	{
		if (bhleft->second.rising && bhright->second.rising && spsup->second.rising && spsdown->second.rising)
		{
			double bhlefttime = bhleft->second.time; 
			double bhrighttime = bhright->second.time; 
			double spsuptime = spsup->second.time; 
			double spsdowntime = spsdown->second.time; 

			double TOF = (spsuptime+spsdowntime)/2-(bhlefttime+bhrighttime)/2;
			double RF = fmod((bhlefttime+bhrighttime)/2-reftimeBH+rfoffset,19.75);
			double Shifted_RF = fmod((bhlefttime+bhrighttime)/2-reftimeBH+rfoffset+rfcor[bhplane][bhbar],19.75);

			H1(TOF,TString::Format("SPS/%s/Bar %d/BH %s/Paddle %d/TOF",wallnames[spswall].c_str(),spsbar,planenames[bhplane].c_str(),bhbar),"TOF;time(ns);counts",5000,0,50);
			H2(RF,TOF,TString::Format("SPS/%s/Bar %d/BH %s/Paddle %d/TOF v RF",wallnames[spswall].c_str(),spsbar,planenames[bhplane].c_str(),bhbar),"TOF v RF;RF Time(ns);TOF(ns)",857,0,19.75,200,10,30);
		}
	}

 //  auto bh1left = findSmallestTime(bh1->tdc_trb[0]);
 //  auto bh1right = findSmallestTime(bh1->tdc_trb[1]);
 //  auto sps1up = findSmallestTime(sps1->tdc_trb[0]);
 //  auto sps1down = findSmallestTime(sps1->tdc_trb[1]);
 //  auto bh2left = findSmallestTime(bh2->tdc_trb[0]);
 //  auto bh2right = findSmallestTime(bh2->tdc_trb[1]);

 //  auto sps1UpLeading = findSmallestTime(sps1Up->tdc_trb[0]);
 //  auto sps1UpTrailing = findSmallestTime(sps1Up->tdc_trb[1]);
 //  auto sps1DownLeading = findSmallestTime(sps1Down->tdc_trb[0]);
 //  auto sps1DownTrailing = findSmallestTime(sps1Down->tdc_trb[1]);


 //  if ((sps1up!=sps1->tdc_trb[0].end()) && (sps1UpLeading!=sps1Up->tdc_trb[0].end()))
 //    H1(sps1up->second.time-sps1UpLeading->second.time,"sps diff","sps diff",100,0,10);

 //  double uptot=0;
 //  double downtot=0;
 //  if ( (sps1UpLeading!=sps1Up->tdc_trb[0].end()) &&(sps1UpTrailing!=sps1Up->tdc_trb[1].end()) &&(sps1DownLeading!=sps1Down->tdc_trb[0].end()) && (sps1DownTrailing!=sps1Down->tdc_trb[1].end()))
 //    {
 //      uptot=sps1UpTrailing->second.time-sps1UpLeading->second.time;
 //      downtot=sps1DownTrailing->second.time-sps1DownLeading->second.time;
	   
 //    }

 //  for(auto hitadc:sps1->adc_mqdc[1])
 //    {
 //      for(auto hitadc1:sps1->adc_mqdc[0])
	// {
	//   // No cuts needed for QDC spectra!
	//   H1(hitadc,"QDC UP","QDC UP ;QDC UP value;counts",6000,-2000,4000);
	//   H1(hitadc1,"QDC DOWN","QDC DOWN ;QDC DOWN value;counts",6000,-2000,4000);
	//   // H1(fmod((sps1up->second.time+sps1down->second.time)/2-reftime,RFTIME_1290) ,"rf","RF; P-DOWN resolution;time [ns];counts",1000,-100,100);

	//   if ( (bh1left!=bh1->tdc_trb[0].end()) && (bh1right!=bh1->tdc_trb[1].end()) )
	//     {
	// 	      if ((sps1up!=sps1->tdc_trb[0].end()) && (sps1down!=sps1->tdc_trb[1].end()))
	// 		{

	// 		  if ((bh2left!=bh2->tdc_trb[0].end()) && (bh2right!=bh2->tdc_trb[1].end()))
	// 		    {

	// 		      double adcavg=sqrt(hitadc*hitadc1);
	// 		      double totavg=sqrt(uptot*downtot);
	// 		      double steffent=(bh1left->second.time+bh1right->second.time+bh2left->second.time+bh2right->second.time)/4. - (sps1up->second.time+sps1down->second.time)/2;
	// 		      double updiff=(bh1left->second.time+bh1right->second.time+bh2left->second.time+bh2right->second.time)/4. - sps1up->second.time;
	// 		      double steffent_totcorr=(steffent+33.3)-(totavg-31)*0.15;
	// 		      double U_minus_D = sps1UpLeading->second.time-sps1DownLeading->second.time;

	// 		      //We want to add up all bh times and divide by four and subtract the avg SPS
	// 		      H1(steffent,"Steffen's T Test","Steffen's T Test;time [ns];counts",15000,-200,200);
	// 		      H1(adcavg,"Geo mean of SPS","Geo mean of SPS; QDC value;counts",4095,-0.5,4094.5);
	// 		      H2(adcavg,steffent,"Steffen's T V Geo mean","Steffen's T V Geo mean;QDC;time[ns]",3000,999.5,2999.5,100,-35,-30);
	// 		      H2(adcavg,totavg,"ADC vs TOT","ADC vs tot",300,499.5,3999.5,100,26,36);
	// 		      H2(totavg,steffent,"TOT vs Steffen's T","ToT vs Steffen's T",100,26,36,150,-35,-30);
	// 		      H2( fmod((sps1UpLeading->second.time+sps1DownLeading->second.time)/2-reftime,RFTIME_1290),totavg,"Tot vs tof","tot vs tof",1000,-100,100,100,26,36);

	// 		      H2( fmod((sps1UpLeading->second.time+sps1DownLeading->second.time)/2-reftime,RFTIME_1290),adcavg,"ADC vs tof","ADC vs tof",1000,-100,100,1000,600,1600);
						  
	// 		      H1(uptot,"Up tot","Up tot",100,25,35);
	// 		      H1(downtot,"Down tot","Down tot",100,25,35);
	// 		      H2(uptot,hitadc,"Up tot vs adc","Up tot vs adc",100,25,35,4096,-0.5,4095.5);

	// 		      H2( fmod((sps1up->second.time+sps1down->second.time)/2-reftime,RFTIME_1290),totavg,"Tot vs tof2","tot vs tof2",1000,-100,100,100,26,36);
	// 		      H2(hitadc,updiff,"test","test",2000,499.5,1999.5,100,-40,-30);
	// 		      H2(totavg,steffent_totcorr,"TOT vs Steffen's T corrected","ToT vs Steffen's T corrected",100, 26,36,100,-1.5,5.5);
	// 		      H1(U_minus_D ,"UP-DOWN resolution","UP-DOWN resolution;time [ns];counts",1000,-5,5);
	// 		      //H1(hitadc,"QDC UP","QDC UP ;QDC UP value;counts",6000,-2000,4000);
	// 		      //H1(hitadc1,"QDC DOWN","QDC DOWN ;QDC DOWN value;counts",6000,-2000,4000);

	// 		    }
	// 		}	
	//     }
	// }	
 //    }
}

void TOF::TOFBM(BHbar *bh, BMbar *bm, int bmbar, int bhplane, int bhbar)
// void TOF::TOFBM(BHbar *bh, BMbar *bm)s
{
	auto bhleft = findSmallestTime(bh->tdc_trb[0]);
	auto bhright = findSmallestTime(bh->tdc_trb[1]);
	auto bmup = findSmallestTime(bm->tdc_trb[0]);
	auto bmdown = findSmallestTime(bm->tdc_trb[1]);

	double bhlefttime = bhleft->second.time;
	double bhrighttime = bhright->second.time;
	double bmdowntime = bmdown->second.time;
	double bmuptime = bmup->second.time;
	
	if(bhleft!=bh->tdc_trb[0].end()  && bhright!=bh->tdc_trb[1].end() && bmup!=bm->tdc_trb[0].end() && bmdown!=bm->tdc_trb[1].end())
	{
	  if(bhright->second.rising && bhleft->second.rising && bmup->second.rising && bmdown->second.rising)
		{
			double TOF = (-(bhlefttime+bhrighttime)/2+(bmuptime+bmdowntime)/2);
			double bmdiff = bmuptime-bmdowntime;
			double bhdiff = bhrighttime-bhlefttime;

			double bmdownmult = bm->tdc_trb[0].size();
			double bmupmult = bm->tdc_trb[1].size();
			double RF = fmod((bhlefttime+bhrighttime)/2-reftimeBH+rfoffset,19.75);
			double RFBM = fmod((bmuptime+bmdowntime)/2-reftimeBM+rfoffset,19.75);
			double shiftedRF = fmod((bhlefttime+bhrighttime)/2-reftimeBH+rfoffset+rfcor[bhplane][bhbar],19.75);
			double time_sub = (-((bhlefttime+bhrighttime)/2- RF)+((bmdowntime+bmuptime)/2)-RFBM);

			H1( time_sub, TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/time_sub",bmbar,planenames[bhplane].c_str(),bhbar),"time_sub;time(ns);counts",2500,-35,35);
			H1(TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/TOF",bmbar,planenames[bhplane].c_str(),bhbar),"TOF;time(ns);counts",2500,10,35);
			H1(RF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/RF",bmbar,planenames[bhplane].c_str(),bhbar),"RF;time(ns)",857,0,19.75);
			H1(RFBM,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/RF BM",bmbar,planenames[bhplane].c_str(),bhbar),"RF;time(ns)",857,0,19.75);
			H1(shiftedRF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/shiftedRF",bmbar,planenames[bhplane].c_str(),bhbar),"Shifted RF;time(ns)",857,0,19.75);
			 H2(RF,TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BH RF V TOF",bmbar,planenames[bhplane].c_str(),bhbar),"BH RF V TOF;RF time (ns);TOF (ns)",857,0,19.75,1000,0,50);
			H2(shiftedRF,TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BH shiftedRF V TOF",bmbar,planenames[bhplane].c_str(),bhbar),"BH Shifted RF V TOF;RF time (ns);TOF (ns)",857,0,19.75,2500,10,35);
			H2(bmdiff,TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/TOF V BM Time Diff",bmbar,planenames[bhplane].c_str(),bhbar),"TOF V BM Time Difference;BM Up Time - Down Time (ns);TOF (ns)",400,-20,20,1000,0,50);
			H2(bhdiff,TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/TOF V BH Time Diff",bmbar,planenames[bhplane].c_str(),bhbar),"TOF V BH Time Difference;BH Up Time - Down Time (ns);TOF (ns)",400,-20,20,1000,0,50);


			// if (RF > rfcut[0] && RF < rfcut[1]) // electron
			// 	H1(TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/TOFee",bmbar,planenames[bhplane].c_str(),bhbar),"TOF Electron;time(ns);counts",1000,0,50);	
			// if (RF > rfcut[2] && RF < rfcut[3]) // pion
			// 	H1(TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/TOFpi",bmbar,planenames[bhplane].c_str(),bhbar),"TOF Pion;time(ns);counts",1000,0,50);
			// if (RF > rfcut[4] && RF < rfcut[5]) // muon
			// 	H1(TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/TOFmu",bmbar,planenames[bhplane].c_str(),bhbar),"TOF Muon;time(ns);counts",1000,0,50);	

			// if ( bm->adc_mqdc[1].size() > 0 && bm->adc_mqdc[0].size() > 0)
			// {	  
			//	auto bmqdcup = bm->adc_mqdc[1][0];
			//	auto bmqdcdown = bm->adc_mqdc[0][0];
			//	H1(sqrt(bmqdcdown*bmqdcup),TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BM QDC GeoMean",bmbar,planenames[bhplane].c2500,10,35_str(),bhbar),"BM QDC Spectrum;Channel",4096,0,4096);
			//	//H2(TOF,sqrt(bmqdcdown*bmqdcup),TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BM QDC V TOF",bmbar,planenames[bhplane].c_str(),bhbar),"BM QDC V TOF;time (ns);channel",1000,0,50,4096,0,4096);
			//	//H2(RF,sqrt(bmqdcdown*bmqdcup),TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BM QDC V RF",bmbar,planenames[bhplane].c_str(),bhbar),"BM QDC V RF;time (ns);channel",857,0,19.75,4096,0,4096);
			//}	

			//if ( bh->adc_mqdc[1].size() > 0 && bh->adc_mqdc[0].size() > 0 )
			//{
			//	auto bhqdcup = bh->adc_mqdc[1][0];
			//	auto bhqdcdown = bh->adc_mqdc[0][0];
			//	H1(sqrt(bhqdcdown*bhqdcup),TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BH QDC GeoMean",bmbar,planenames[bhplane].c_str(),bhbar),"BH QDC Spectrum;Channel",4096,0,4096);
			//	H2(TOF,sqrt(bhqdcdown*bhqdcup),TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BH QDC V TOF",bmbar,planenames[bhplane].c_str(),bhbar),"BH QDC V TOF;time (ns);channel",1000,0,50,4096,0,4096);
			//	H2(RF,sqrt(bhqdcdown*bhqdcup),TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BH QDC V RF",bmbar,planenames[bhplane].c_str(),bhbar),"BH QDC V RF;time (ns);channel",857,0,19.75,4096,0,4096);
		  	
			//}
		}
	
	 }	

	return;
}


Long_t TOF::defineHistograms()
{

 for (int k = 0; k < 4; k ++)
	 {
	 	
			
	RF_pi_plane[k]= dH1(TString::Format("BH Wiremap with BM constraint/%s/both Pi",planenames[k].c_str()),TString::Format("Wiremap %s BH BM constraint Pi",planenames[k].c_str()),16,-0.5,15.5);
	RF_mu_plane[k]= dH1(TString::Format("BH Wiremap with BM constraint/%s/both Mu",planenames[k].c_str()),TString::Format("Wiremap %s BH BM constraint Mu",planenames[k].c_str()),16,-0.5,15.5);
	RF_e_plane[k]= dH1(TString::Format("BH Wiremap with BM constraint/%s/both e",planenames[k].c_str()),TString::Format("Wiremap %s BH BM constraint e",planenames[k].c_str()),16,-0.5,15.5);

}

	// for (int i = 0; i < 4; i ++)
	// {
	// 	for (int j = 2; j < 4; j ++)
	// 	{
	// 		for (int k = 0; k < 16; k ++)
	// 		{
	// 			BM_QCDvTOF[i][j-2][k] = dH2(TString::Format("Beam Monitor/BM Wall 0/BMbar %d/%s/TOF BH paddle %d/BM QDC V TOF",i,planenames[j].c_str(),k),"BM QDC V TOF;TOF(ns);Channel",400,10,40,1996,100,2096);
	// 			BM_QCDvRF[i][j-2][k] = dH2(TString::Format("Beam Monitor/BM Wall 0/BMbar %d/%s/TOF BH paddle %d/BM QDC V RF",i,planenames[j].c_str(),k),"BM QDC V RF;RF(ns);Channel",857,0,19.75,1996,100,2096);
	// 			BH_QCDvTOF[i][j-2][k] = dH2(TString::Format("Beam Monitor/BM Wall 0/BMbar %d/%s/TOF BH paddle %d/BH QDC V TOF",i,planenames[j].c_str(),k),"BH QDC V TOF;TOF(ns);Channel",400,10,40,1996,100,2096);
	// 			BH_QCDvRF[i][j-2][k] = dH2(TString::Format("Beam Monitor/BM Wall 0/BMbar %d/%s/TOF BH paddle %d/BH QDC V RF",i,planenames[j].c_str(),k),"BH QDC V RF;RF(ns);Channel",857,0,19.75,1996,100,2096);
	// 		}

	// 	}
	// }
	
  return ok;
}

Long_t TOF::startup()
{

  
	bhraw = NULL;
	getBranchObject("BH",(TObject **) &bhraw);
	if (!bhraw) {
		debug(0,"Could not find bh tree in file\n");
	}

	spsraw = NULL;
	getBranchObject("SPS",(TObject **) &spsraw);
	if (!spsraw) {
		debug(0,"Could not find SPS tree in file\n");
	}

	bmraw = NULL;
	getBranchObject("BM",(TObject **) &bmraw);
	if (!bmraw)
		debug(0,"Could not find Beam Monitor in file\n");

	vetoraw = NULL;
	getBranchObject("VETO",(TObject **) &vetoraw);
	if (!bmraw)
		debug(0,"Could not find VETO in file\n");


	theRunInfo=(MRTRunInfo *) getFileObject("RunInfo");


	/// align rf
	double rf_center = 8.;
	double rf_peak[4][16] = {{7.38472e+00, 9.32187e+00, 6.68016e+00, 9.56270e+00, 6.67124e+00, 8.43312e+00, 8.41714e+00, 7.30812e+00, 8.61330e+00, 
		                          9.55535e+00, 9.79961e+00, 9.98140e+00, 5.63248e+00, 6.87634e+00, 7.62333e+00, 7.88518e+00},
				 {6.86477e+00, 8.54564e+00,  6.16950e+00, 9.05907e+00, 6.08416e+00,  7.35111e+00, 7.14500e+00, 6.30274e+00, 7.52484e+00, 8.46922e+00,  8.72917e+00, 1.28824e+01, 4.45561e+00, 0,0,0},
                           	  {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		     	          {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};

	for (int i = 0; i < 4; i ++)
		for (int j = 0; j < 16; j ++)
	    		rfcor[i][j] = rf_center - rf_peak[i][j];

	return ok;
}

Long_t TOF::process()
{

	// if  ((bhraw->trb_reftime.size()==0) || (spsraw->trb_reftime.size()==0) )
	if  ((bhraw->trb_reftime.size()==0))
	{
		//debug(0,"No BH reftime found: %i\n",bhraw->trb_reftime.size());
		return ok;
	}
	reftimeBH= findReference(bhraw->trb_reftime)->time;
	if  ((bmraw->trb_reftime.size()==0))
	{
		//debug(0,"No BM reftime found: %i\n",bmraw->trb_reftime.size());
		return ok;
	}
	reftimeBM= findReference(bmraw->trb_reftime)->time;
	//reftimeSPS= findReference(spsraw->trb_reftime)->time;
	// reftimeVETO= findReference(vetoraw->trb_reftime)->time;
	// trigger = findTrig(spsraw->trig)->time;
	//

// 	for (int i = 0; i<vetoraw->plane.size();i++)
// 		if ( vetoraw->plane[i].tdc_trb.size() > 0 )
// 			return ok;
	
	
 	int i = 1; // the plane with Steffen's big bars
	for (int j=0;j<bmraw->plane[i].size();j++)
	{
		if ( bmraw->plane[i][j].tdc_trb[0].size() > 0 && bmraw->plane[i][j].tdc_trb[1].size() > 0 )
		{
	/* for ( int check = j+1; check < bmraw->plane[i].size(); check ++ )
			{
				 // if( bmraw->plane[i][check].tdc_trb[0].size() > 0 && bmraw->plane[i][check].tdc_trb[1].size() > 0 )
					//return ok;
			} */

			for(int k=0; k<4;k++)//only look at last two bh planes
			{
			  	for(int l=0;l<bhraw->plane[k].size();l++)
				{
					if ( bhraw->plane[k][l].tdc_trb[0].size() > 0 && bhraw->plane[k][l].tdc_trb[1].size() > 0 )
					{
					/*	for ( int check = l+1; check < bhraw->plane[k].size(); check ++ )
						{
							  if( bhraw->plane[k][check].tdc_trb[0].size() > 0 && bhraw->plane[k][check].tdc_trb[1].size() > 0 )
								return ok; */
						//}
							TOFBM(&bhraw->plane[k][l],&bmraw->plane[i][j], j, k, l);
					}	
				 }

			  }
		   }
	}

	
	for (int k = 0; k < 4; k++)
	{
		for(int i = 0; i < bhraw->plane[k].size(); i++)
		{
			

	auto bh1 = findSmallestTime(bhraw->plane[k][i].tdc_trb[0]);
	auto bh2 = findSmallestTime(bhraw->plane[k][i].tdc_trb[1]);
	//auto bm1 = findSmallestTime(bm->tdc_trb[0]);
	//auto bm2 = findSmallestTime(bm->tdc_trb[1]);

	//double bh1 = bhleft->second.time;
	//double bh2 = bhright->second.time;
	//double bm1 = bmdown->second.time;
	//double bm2 = bmup->second.time;
    //double RF = fmod((bhlefttime+bhrighttime)/2-reftimeBH+rfoffset,19.75);
			bool vert = false;
			std::string sipmname[2];
			if((k % 2) == 1)
			  vert = true;

			if(vert)
			{
			  sipmname[0]={"Down"};
			  sipmname[1]={"Up"};
			}
			if(!vert)
			{
			  sipmname[0]={"Left"};
			  sipmname[1]={"Right"};
			}

			/// check RF lining up
			
			if (bh1->second.rising && bh1!=bhraw->plane[k][i].tdc_trb[0].end() && bh2->second.rising && bh2!=bhraw->plane[k][i].tdc_trb[1].end())
			{
				double RF = fmod((bh1->second.time+bh2->second.time)/2-reftimeBH+rfoffset,19.75); // not shifted
				double ShiftedRF = fmod((bh1->second.time+bh2->second.time)/2-reftimeBH+rfoffset+rfcor[k][i],19.75);
				H2(i,RF,TString::Format("RF/%s rf per bar",planenames[k].c_str()),TString::Format("RF/%s rf per bar",planenames[k].c_str()),16,-0.5,15.5,857,0,19.75);
				H1(RF,TString::Format("RF/%s rf",planenames[k].c_str()),TString::Format("RF/%s rf",planenames[k].c_str()),857,0,19.75);
				H1(ShiftedRF,TString::Format("ShiftedRF/%s rf",planenames[k].c_str()),TString::Format("ShiftedRF/%s rf",planenames[k].c_str()),857,0,19.75);
				//H1(rf,TString::Format("RF/%s/bar %d",planenames[k].c_str(),i),TString::Format("RF %s bar %d",planenames[k].c_str(),i),857,0,19.75);
			}

			/// Correlations
			
			for(int j = 0; j < bmraw->plane[1].size(); j++) /// Big Bars
			{
				auto bm1 = findSmallestTime(bmraw->plane[1][j].tdc_trb[0]);
				auto bm2 = findSmallestTime(bmraw->plane[1][j].tdc_trb[1]);
				if(bh1->second.rising&&bm1->second.rising&&bh1!=bhraw->plane[k][i].tdc_trb[0].end()&&bm1!=bmraw->plane[1][j].tdc_trb[0].end())
				{
					H2(i,j,TString::Format("Correlations/%s to BM Correlation Down PMT",planenames[k].c_str()),TString::Format("%s to BM Correlation Down PMT;bh paddle;bm paddle",planenames[k].c_str()),16,-0.5,15.5,4,-0.5,3.5);
					if(bh2->second.rising&&bm2->second.rising&&bh2!=bhraw->plane[k][i].tdc_trb[1].end()&&bm2!=bmraw->plane[1][j].tdc_trb[1].end())
						H2(i,j,TString::Format("Correlations/%s to BM Correlation PMT",planenames[k].c_str()),TString::Format("%s to BM Correlation Both PMT;bh paddle;bm paddle",planenames[k].c_str()),16,-0.5,15.5,4,-0.5,3.5);
				}
				if(bh2->second.rising&&bm2->second.rising&&bh2!=bhraw->plane[k][i].tdc_trb[1].end()&&bm2!=bmraw->plane[1][j].tdc_trb[1].end())
					H2(i,j,TString::Format("Correlations/%s to BM Correlation Up PMT",planenames[k].c_str()),TString::Format("%s to BM Correlation Up PMT;bh paddle;bm paddle",planenames[k].c_str()),16,-0.5,15.5,4,-0.5,3.5);
			}

			for(int j = 0; j < bmraw->plane[0].size(); j++) /// Small Bars
			{
				auto bm1 = findSmallestTime(bmraw->plane[0][j].tdc_trb[0]);
				auto bm2 = findSmallestTime(bmraw->plane[0][j].tdc_trb[1]);
				if(bh1->second.rising&&bm1->second.rising&&bh1!=bhraw->plane[k][i].tdc_trb[0].end()&&bm1!=bmraw->plane[0][j].tdc_trb[0].end())
				{
					if(bh2->second.rising&&bm2->second.rising&&bh2!=bhraw->plane[k][i].tdc_trb[1].end()&&bm2!=bmraw->plane[0][j].tdc_trb[1].end())
						H2(i,j,TString::Format("Correlations/%s to BM Hodoscope Correlation PMT",planenames[k].c_str()),TString::Format("%s to BM Hodoscope Correlation Both PMT;bh paddle;bm paddle",planenames[k].c_str()),16,-0.5,15.5,32,-0.5,31.5);
				}					
			}	

			/// Wiremap for beam study
		  	
			if(bhraw->plane[k][i].tdc_trb[0].size()>0)
			  H1(i,TString::Format("BH Wiremap/%s/%s",planenames[k].c_str(),sipmname[0].c_str()),TString::Format("Wiremap %s BH %s",planenames[k].c_str(),sipmname[0].c_str()),16,-0.5,15.5);
			if(bhraw->plane[k][i].tdc_trb[1].size()>0)
			  H1(i,TString::Format("BH Wiremap/%s/%s",planenames[k].c_str(),sipmname[1].c_str()),TString::Format("Wiremap %s BH %s",planenames[k].c_str(),sipmname[1].c_str()),16,-0.5,15.5);		
			if(bhraw->plane[k][i].tdc_trb[0].size()>0 && bhraw->plane[k][i].tdc_trb[1].size()>0)
			{
				H1(i,TString::Format("BH Wiremap/%s/both",planenames[k].c_str()),TString::Format("Wiremap %s BH both",planenames[k].c_str()),16,-0.5,15.5);
				
				for (int j = 0; j < bmraw->plane[1].size(); j++)
				{

					auto bmup = findSmallestTime(bmraw->plane[1][j].tdc_trb[1]);
					auto bmdown = findSmallestTime(bmraw->plane[1][j].tdc_trb[0]);

					if(bh1->second.rising&&bh2->second.rising&&bh2!=bhraw->plane[k][i].tdc_trb[1].end()&&bh1!=bhraw->plane[k][i].tdc_trb[0].end())
					{
						if (bmup->second.rising&&bmdown->second.rising&&bmup!=bmraw->plane[1][j].tdc_trb[1].end()&&bmdown!=bmraw->plane[1][j].tdc_trb[0].end())
						{
							double RF = fmod((bh2->second.time+bh1->second.time)/2-reftimeBH+rfoffset+rfcor[k][i],19.75);
							double TOF = (bmup->second.time+bmdown->second.time)/2 - (bh2->second.time+bh1->second.time)/2;

							H2(i,TOF,TString::Format("BM%d %s TOF per bar",j,planenames[k].c_str()),TString::Format("BM%d %s TOF per bar;Bar Number;TOF(ns)",j,planenames[k].c_str()),16,-0.5,15.5,2500,10,35); 

							 /*if ( i == 12 )
							{
								if(RF >2.5 && RF < 5 && TOF < 24)
									H1(i,TString::Format("BH Wiremap/%s/both Pi",planenames[k].c_str()),TString::Format("Wiremap %s BH Pi",planenames[k].c_str()),16,-0.5,15.5);
								if(RF > 7 && RF < 9 && TOF < 24)
									H1(i,TString::Format("BH Wiremap/%s/both e",planenames[k].c_str()),TString::Format("Wiremap %s BH e",planenames[k].c_str()),16,-0.5,15.5);
								if(RF > 16.2 && RF < 18.2 && TOF < 24)
									H1(i,TString::Format("BH Wiremap/%s/both Mu",planenames[k].c_str()),TString::Format("Wiremap %s BH Mu",planenames[k].c_str()),16,-0.5,15.5);
								if(RF > 18.2 || RF < 1.6 && TOF < 24)
									H1(i,TString::Format("BH Wiremap/%s/Mu Shifted",planenames[k].c_str()),TString::Format("Wiremap %s BH Mu Shifted",planenames[k].c_str()),16,-0.5,15.5);
							} */
							//else
							//if(i==1)
							//{
								if(RF >7.2 && RF < 9.5 && TOF>22.5 && TOF < 23.5)
									RF_pi_plane[k]->Fill(i);

									//H1(i,TString::Format("BH Wiremap/%s/both Pi",planenames[k].c_str()),TString::Format("Wiremap %s BH Pi",planenames[k].c_str()),16,-0.5,15.5);
								if(RF > 12 && RF < 14 && TOF> 20 && TOF < 22.5)
									RF_mu_plane[k]-> Fill(i);
								//	H1(i,TString::Format("BH Wiremap/%s/both e",planenames[k].c_str()),TString::Format("Wiremap %s BH e",planenames[k].c_str()),16,-0.5,15.5);

								if(RF > 6.5 && RF < 8.6 && TOF> 17 && TOF < 19.8)
									RF_e_plane[k]-> Fill(i);
								//	H1(i,TString::Format("BH Wiremap/%s/both Mu",planenames[k].c_str()),TString::Format("Wiremap %s BH Mu",planenames[k].c_str()),16,-0.5,15.5);
								//if(RF > 18.2 || RF < 1.6 && TOF < 22)
									//H1(i,TString::Format("BH Wiremap/%s/Mu Shifted",planenames[k].c_str()),TString::Format("Wiremap %s BH Mu Shifted",planenames[k].c_str()),16,-0.5,15.5);
							

						}
					}
					  
				}
			}	
       
		}	
		
	}
	

	return ok;
}

Long_t TOF::finalize()
{
	for (int k=0; k<4; k++)
	{
	//Pi_count_[k]= RF_pi_plane[k]->GetEntries();

	printf("total number of Pions in Plane%d %f  %f \n",k, RF_pi_plane[k]->GetEntries(), RF_pi_plane[k]->GetMean());
	printf("total number of electrons in Plane%d %f  %f \n",k, RF_e_plane[k]->GetEntries(), RF_e_plane[k]->GetMean());
	printf("total number of muons in Plane%d %f  %f \n",k, RF_mu_plane[k]->GetEntries(), RF_mu_plane[k]->GetMean());

}

	return ok;
}

Long_t TOF::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
	Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
	{
		return (Plugin *) new TOF(in,out,inf_,outf_,p);
	}
}


ClassImp(TOF);







///////Trying out new shit above////////////////////////////////////////////////////////////////

/*
#include <TOF.h>

#include <iostream>
#include <cmath>

#include <fstream>
#include <cstdlib>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TPad.h>
#include <TStyle.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TChain.h>


TOF::TOF(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

TOF::~TOF()
{
};

template <typename T>
auto findTrig(T &container)
{
    auto best=container.begin();
  for (auto hit=container.begin();hit!=container.end();hit++)
  {
      if (hit->rising && (best->time< hit->time))
      { //&& !(hit->trailing)
        best=hit;
      }
    }
return best;
}

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
auto findSmallestRisingTime(T &container) ///find smallest time???
{
	auto best=container.begin();
	for (auto hit=container.begin();hit!=container.end();hit++)
	{
  	  if ((best->second.time > hit->second.time)&& (hit->second.rising)) {
	    best=hit; //this is different too???
  	  }
    }
return best;
}

///SPECIAL FUNCTION FOR FINDING THE FIRST RISING RF OR TRIG EDGE////
///DONE DUE TO SPECIAL CLASS DEFINITION OF RF IN bh///////////////
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

void TOF::plotST(double reftime,BHbar *bh1, SPSbar *sps1, BHbar *bh2, SPSbar *sps1Up, SPSbar *sps1Down) // do I not need this whole section???
{
	auto bh1left = findSmallestTime(bh1->tdc_trb[0]);
	auto bh1right = findSmallestTime(bh1->tdc_trb[1]);
	auto sps1up = findSmallestTime(sps1->tdc_trb[0]);
	auto sps1down = findSmallestTime(sps1->tdc_trb[1]);
	auto bh2left = findSmallestTime(bh2->tdc_trb[0]);
	auto bh2right = findSmallestTime(bh2->tdc_trb[1]);

        auto sps1UpLeading = findSmallestTime(sps1Up->tdc_trb[0]);
        auto sps1UpTrailing = findSmallestTime(sps1Up->tdc_trb[1]);
        auto sps1DownLeading = findSmallestTime(sps1Down->tdc_trb[0]);
        auto sps1DownTrailing = findSmallestTime(sps1Down->tdc_trb[1]);


	if ((sps1up!=sps1->tdc_trb[0].end()) && (sps1UpLeading!=sps1Up->tdc_trb[0].end()))
	    H1(sps1up->second.time-sps1UpLeading->second.time,"sps diff","sps diff",100,0,10);

	double uptot=0;
	double downtot=0;
	if ( (sps1UpLeading!=sps1Up->tdc_trb[0].end()) &&(sps1UpTrailing!=sps1Up->tdc_trb[1].end()) &&(sps1DownLeading!=sps1Down->tdc_trb[0].end()) && (sps1DownTrailing!=sps1Down->tdc_trb[1].end()))
	  {
	    uptot=sps1UpTrailing->second.time-sps1UpLeading->second.time;
	    downtot=sps1DownTrailing->second.time-sps1DownLeading->second.time;
	   
	  }

	for(auto hitadc:sps1->adc_mqdc[1])
	{
		for(auto hitadc1:sps1->adc_mqdc[0])
		{
			// No cuts needed for QDC spectra!
                        H1(hitadc,"QDC UP","QDC UP ;QDC UP value;counts",6000,-2000,4000);
                        H1(hitadc1,"QDC DOWN","QDC DOWN ;QDC DOWN value;counts",6000,-2000,4000);
			// H1(fmod((sps1up->second.time+sps1down->second.time)/2-reftime,RFTIME_1290) ,"rf","RF; P-DOWN resolution;time [ns];counts",1000,-100,100);

			if ( (bh1left!=bh1->tdc_trb[0].end()) && (bh1right!=bh1->tdc_trb[1].end()) )
			{
				if ((sps1up!=sps1->tdc_trb[0].end()) && (sps1down!=sps1->tdc_trb[1].end()))
				{

					if ((bh2left!=bh2->tdc_trb[0].end()) && (bh2right!=bh2->tdc_trb[1].end()))
					{

					  double adcavg=sqrt(hitadc*hitadc1);
					  double totavg=sqrt(uptot*downtot);
					  double steffent=(bh1left->second.time+bh1right->second.time+bh2left->second.time+bh2right->second.time)/4. - (sps1up->second.time+sps1down->second.time)/2;
					  double updiff=(bh1left->second.time+bh1right->second.time+bh2left->second.time+bh2right->second.time)/4. - sps1up->second.time;
					  double steffent_totcorr=(steffent+33.3)-(totavg-31)*0.15;
					  double U_minus_D = sps1UpLeading->second.time-sps1DownLeading->second.time;

					  //We want to add up all bh times and divide by four and subtract the avg SPS
					  H1(steffent,"Steffen's T Test","Steffen's T Test;time [ns];counts",15000,-200,200);
					  H1(adcavg,"Geo mean of SPS","Geo mean of SPS; QDC value;counts",4095,-0.5,4094.5);
					  H2(adcavg,steffent,"Steffen's T V Geo mean","Steffen's T V Geo mean;QDC;time[ns]",3000,999.5,2999.5,100,-35,-30);
					  H2(adcavg,totavg,"ADC vs TOT","ADC vs tot",300,499.5,3999.5,100,26,36);
					  H2(totavg,steffent,"TOT vs Steffen's T","ToT vs Steffen's T",100,26,36,150,-35,-30);
					  H2( fmod((sps1UpLeading->second.time+sps1DownLeading->second.time)/2-reftime,RFTIME_1290),totavg,"Tot vs tof","tot vs tof",1000,-100,100,100,26,36);

					  H2( fmod((sps1UpLeading->second.time+sps1DownLeading->second.time)/2-reftime,RFTIME_1290),adcavg,"ADC vs tof","ADC vs tof",1000,-100,100,1000,600,1600);
					  
					  H1(uptot,"Up tot","Up tot",100,25,35);
					  H1(downtot,"Down tot","Down tot",100,25,35);
					  H2(uptot,hitadc,"Up tot vs adc","Up tot vs adc",100,25,35,4096,-0.5,4095.5);

					  H2( fmod((sps1up->second.time+sps1down->second.time)/2-reftime,RFTIME_1290),totavg,"Tot vs tof2","tot vs tof2",1000,-100,100,100,26,36);
					  H2(hitadc,updiff,"test","test",2000,499.5,1999.5,100,-40,-30);
					  H2(totavg,steffent_totcorr,"TOT vs Steffen's T corrected","ToT vs Steffen's T corrected",100, 26,36,100,-1.5,5.5);
                                          H1(U_minus_D ,"UP-DOWN resolution","UP-DOWN resolution;time [ns];counts",1000,-5,5);
					  //H1(hitadc,"QDC UP","QDC UP ;QDC UP value;counts",6000,-2000,4000);
                                          //H1(hitadc1,"QDC DOWN","QDC DOWN ;QDC DOWN value;counts",6000,-2000,4000);

					}
				}	
			}
		}	
	}
}

void TOF::TOFall(BHbar *bh, SPSbar *sps) //or this one???
{
	auto bhleft = findSmallestTime(bh->tdc_trb[0]);
	auto bhright = findSmallestTime(bh->tdc_trb[1]);
	auto spsup = findSmallestTime(sps->tdc_trb[0]);
	auto spsdown = findSmallestTime(sps->tdc_trb[1]);

	double bhlefttime = bhleft->second.time;
	double bhrighttime = bhright->second.time;
		double spsdowntime = spsdown->second.time;
	double spsuptime = spsup->second.time;

	if(bhleft!=bh->tdc_trb[0].end()  && bhright!=bh->tdc_trb[1].end() && spsup!=sps->tdc_trb[0].end() && spsdown!=sps->tdc_trb[1].end())
	{
		if(bhright->second.rising && bhleft->second.rising && spsup->second.rising && spsdown->second.rising)
		{
			H1((spsuptime+spsdowntime)/2-(bhlefttime+bhrighttime)/2,"TOF","TOF;time(ns);counts",2000,-50,50);
			H2((spsuptime+spsdowntime)/2-(bhlefttime+bhrighttime)/2, fmod((spsuptime+spsdowntime)/2-reftimeSPS+rfoffset,19.75),"TOF V RF SPS","TOF V RF SPS",2000,-50,50,817,0,19.75);
//Ievgen Dec 21. 2017: Added TOF vs RF bh:
			H2((spsuptime+spsdowntime)/2-(bhlefttime+bhrighttime)/2, fmod((bhlefttime+bhrighttime)/2-reftimeBH+rfoffset,19.75),"TOF V RF bh","TOF V RF bh",2000,-50,50,817,0,19.75);
			H2((spsuptime+spsdowntime)/2-(bhlefttime+bhrighttime)/2, spsdowntime-spsuptime,"TOF V dif SPS","TOF V dif SPS",2000,-50,50,817,-40,40);
		}
	}
	if(bhright->second.rising && bhleft->second.rising &&bhleft!=bh->tdc_trb[0].end()  && bhright!=bh->tdc_trb[1].end())
		{H1(fmod((bhrighttime+bhlefttime)/2-reftimeBH+rfoffset,19.75),"RF time BH","RF time BH",817,0,19.75);}
	if(spsup->second.rising && spsdown->second.rising && spsup!=sps->tdc_trb[0].end() && spsdown!=sps->tdc_trb[1].end())
		{H1(fmod((spsuptime+spsdowntime)/2-reftimeSPS+rfoffset,19.75),"RF time SPS","RF time SPS",817,0,19.75);}

	return;
}

void TOF::TOFBM(BHbar *bh, BMbar *bm, int bmbar, int bhplane, int bhbar)
// void TOF::TOFBM(BHbar *bh, BMbar *bm)
{
  
        auto bmqdcup= bm->adc_mqdc[1];
        auto bmqdcdown= bm->adc_mqdc[0];
	auto bhleft = findSmallestTime(bh->tdc_trb[0]);
	auto bhright = findSmallestTime(bh->tdc_trb[1]);
	auto bmup = findSmallestTime(bm->tdc_trb[0]);
	auto bmdown = findSmallestTime(bm->tdc_trb[1]);
	
	double bhlefttime = bhleft->second.time;
	double bhrighttime = bhright->second.time;
	double bmdowntime = bmdown->second.time;
	double bmuptime = bmup->second.time;
	
	// for(auto bmqdcdown:bm->adc_mqdc[1])
	// 	H1(bmqdcdown,"QDC Down","QDC Down;time (ns);channel",4096,0,4096);
	// for(auto bmqdcup:bm->adc_mqdc[0])
	//	H1(bmqdcup,"QDC Up","QDC Up;time (ns);channel",4096,0,4096);
	if(bhleft!=bh->tdc_trb[0].end()  && bhright!=bh->tdc_trb[1].end() && bmup!=bm->tdc_trb[0].end() && bmdown!=bm->tdc_trb[1].end())
	{
		if(bhright->second.rising && bhleft->second.rising && bmup->second.rising && bmdown->second.rising)
		{
			double TOF = (-(bhlefttime+bhrighttime)/2+(bmuptime+bmdowntime)/2);
			

			//	double shiftedTOF;
			//	if (bhplane == 2)
			//	shiftedTOF = TOF-rfcor[bhplane-2][bhbar]+BMrfcor[bmbar];
				//	if (bhplane == 3)

			//double shiftedTOF;
			//if (bhplane == 2)
			//	shiftedTOF = TOF-rfcor[bhplane-2][bhbar]+BMrfcor[bmbar];
			//if (bhplane == 3)

			//	shiftedTOF = TOF-rfcor[bhplane-2][bhbar]+BMrfcor[bmbar]+tofcor[bmbar][bhbar];

			double bmdownmult = bm->tdc_trb[0].size();
			double bmupmult = bm->tdc_trb[1].size();

			//	double RF = fmod((bhlefttime+bhrighttime)/2-reftimeBH+rfoffset+rfcor[bhplane-2][bhbar],19.75);
			double RF = fmod((bhlefttime+bhrighttime)/2-reftimeBH+rfoffset,19.75);

			   // 	 for(auto bmqdcup:bm->adc_mqdc[1])
			   //   {
			   //     for(auto bmqdcdown: bm->adc_mqdc[0])
			   // 	 {
			   // 	   //   H2(RF,bmqdcup,"BM RF V QDC up","BM RF V QDC BM  up;RF time (ns);QDC (channels)",857,0,19.75,4096,0,4096);
			   // 	   //   H2(TOF,bmqdcup,"TOF V QDC up BM","TOF V QDC BM up;time (ns);QDC (channels)",1000,0,50,4096,0,4096);
			   // //   H2(RF, sqrt(bmqdcup*bmqdcdown), "QDC V RF BM", "QDC V RF BM; time(ns);QDC (channels)", 857, 0,19.75, 4096,9,4096);
			   // 	   H2(RF,sqrt(bmqdcdown*bmqdcup),"QDC V RF BM","QDC V RF;time (ns);channel",857,0,19.75,2000,0,2000);
			    	 
			   
			     
			   


			
			H1(TOF,"TOF","TOF;time(ns);counts",1000,0,50);
			//H1(RF,"RF","RF;time(ns)",857,0,19.75);

			 for(auto bmqdcup:bm->adc_mqdc[1])
			     {
			       for(auto bmqdcdown: bm->adc_mqdc[0])
				 {
				   for(auto bhqdcleft:bh->adc_mqdc[0])
				     {
				       for(auto bhqdcright: bh->adc_mqdc[1])
					 {
				   //   H2(RF,bmqdcup,"BM RF V QDC up","BM RF V QDC BM  up;RF time (ns);QDC (channels)",857,0,19.75,4096,0,4096);
				   //   H2(TOF,bmqdcup,"TOF V QDC up BM","TOF V QDC BM up;time (ns);QDC (channels)",1000,0,50,4096,0,4096);
			   //   H2(RF, sqrt(bmqdcup*bmqdcdown), "QDC V RF BM", "QDC V RF BM; time(ns);QDC (channels)", 857, 0,19.75, 4096,9,4096);
					   // H2(RF,sqrt(bmqdcdown*bmqdcup),"QDC V RF BM","QDC V RF;time (ns);channel",857,0,19.75,2000,0,2000);
					   //H1(sqrt(bmqdcdown*bmqdcup),"QDC BM Geometric Mean","QDC BM Geometric Mean;channel;counts",4096,0,4096);
				 
			     
				 
			     
			//	H2(RF,TOF,"BH TOF V RF","BH TOF V RF;RF time (ns);TOF (ns)",857,0,19.75,1000,0,50);
			//	H1(shiftedTOF,"shiftedTOF","TOF with offset;time(ns);counts",1000,0,50);
			// H2(bmuptime-bmdowntime,TOF,"TOF V UP-DOWN","TOF V UP-DOWN Time Diff",200,-10,10,1000,0,50);

			//	TOFperbar[bhplane-2][bmbar]->Fill(bhbar,shiftedTOF);
			// planeTOF->Fill(shiftedTOF);

			//double RF = fmod((bhlefttime+bhrighttime)/2-reftimeBH+rfoffset+rfcor[bhplane-2][bhbar],19.75);

					   
			for(auto bmqdcup:bm->adc_mqdc[0])
			  {
			 // H2(RF,bmqdcup,"BH RF V QDC up","BH RF V QDC up;RF time (ns);QDC (ns)",857,0,19.75,4096,0,4096);
			  H2(TOF,bmqdcup,"TOF V QDC up","TOF V QDC up;time (ns);QDC (ns)",1000,0,50,4096,0,4096);
			  }

			/*
			H1(TOF,"TOF","TOF;time(ns);counts",1000,0,50);
			//H1(RF,"RF","RF;time(ns)",857,0,19.75);
			//H2(RF,TOF,"BH RF V TOF","BH RF V TOF;RF time (ns);TOF (ns)",857,0,19.75,1000,0,50);
			H1(shiftedTOF,"shiftedTOF","TOF with offset;time(ns);counts",1000,0,50);
			 H2(bmuptime-bmdowntime,TOF,"TOF V UP-DOWN","TOF V UP-DOWN Time Diff",200,-10,10,1000,0,50);

			 
			TOFperbar[bhplane-2][bmbar]->Fill(bhbar,shiftedTOF);
			 planeTOF->Fill(shiftedTOF);

			
			 if ( bhplane == 3 && shiftedTOF < 20.8 && shiftedTOF > 18)
			 	TOFee[bmbar][bhbar]->Fill(shiftedTOF);

			*/ ///A change Nov16//
				
				   //if (RF > rfcut[bhplane-2][0] && RF < rfcut[bhplane-2][1])
				   //	H1(TOF,"TOFee","TOF Electron;time(ns);counts",1000,0,50);	
			 // if (RF > rfcut[bhplane-2][2] && RF < rfcut[bhplane-2][3])
			 //	H1(TOF,"TOFpi","TOF Pion;time(ns);counts",1000,0,50);
			 // if (RF > rfcut[bhplane-2][4] && RF < rfcut[bhplane-2][5])
			// //	H1(TOF,"TOFmu","TOF Pion;time(ns);counts",1000,0,50);
/*
				   if (((RF > 4 && RF < 5.5)))
				     //if(sqrt(bmqdcdown*bmqdcup)>650 && sqrt(bmqdcdown*bmqdcup)<900 )
				     //if(sqrt(bhqdcright*bhqdcleft)>130 && sqrt(bhqdcright*bhqdcleft)<200 )
					 {
					 H1(TOF,"TOFee","TOF Electron;time(ns);counts",1000,0,50);
					 //H1(sqrt(bmqdcdown*bmqdcup),"QDCee","Electron Geometric Mean QDC BM;channel;counts",2000,0,2000);
					 }

				   //  if (((RF > 2 && RF < 4.1)||RF>19.2)&&(sqrt(bmqdcdown*bmqdcup)>700 && sqrt(bmqdcdown*bmqdcup)<1100))
				   if ((RF > 19 && RF < 20))
				      //if(sqrt(bmqdcdown*bmqdcup)>750 && sqrt(bmqdcdown*bmqdcup)<1000 )
				      //if(sqrt(bhqdcright*bhqdcleft)>140 && sqrt(bhqdcright*bhqdcleft)<240 )
					 {
				       	H1(TOF,"TOFpi","TOF Pion;time(ns);counts",1000,0,50);
					//H1(sqrt(bmqdcdown*bmqdcup),"QDCpi","Pion Geometric Mean QDC BM;channel;counts",2000,0,2000);
					 }

				   if ((RF > 13 && RF < 14.5))
				     //if(sqrt(bmqdcdown*bmqdcup)>550 && sqrt(bmqdcdown*bmqdcup)<800 )
				     //if(sqrt(bhqdcright*bhqdcleft)>120 && sqrt(bhqdcright*bhqdcleft)<180 )
					 {
				        H1(TOF,"TOFmu","TOF Muon;time(ns);counts",1000,0,50);
					//H1(sqrt(bmqdcdown*bmqdcup),"QDCmu","Muon Geometric Mean QDC BM;channel;counts",2000,0,2000);
					 }

				    
				       

			 // 	for(auto bmqdcup:bm->adc_mqdc[1])
			 // 	{
			 // 	  //	for(auto bmqdcdown:bm->adc_mqdc[0])
			 // 	  //	 {
			 //   // H2(TOF,sqrt(bmqdcdown*bmqdcup),"QDC V TOF","QDC V TOF;time (ns);channel",1000,0,50,4096,0,4096);
			 //   H1(TOF,"TOF","TOF;time(ns);counts",1000,0,50);
			 //   // H1(sqrt(bmqdcdown*bmqdcup), "Geometric mean QDC BM", "Geometric mean QDC;channel; counsts", 4096, 0,4096);
			 //    	H1(RF,"RF","RF;time(ns)",857,0,19.75);
			 // 	//	 H2(RF,sqrt(bmqdcdown*bmqdcup),"QDC V RF","QDC V RF;time (ns);channel",857,0,19.75,4096,0,4096);
			 // 	// }
			 // }

			        
			   // H2(TOF,sqrt(bmqdcdown*bmqdcup),"QDC V TOF","QDC V TOF;time (ns);channel",1000,0,50,4096,0,4096);
				     //  H1(TOF,"TOF","TOF;time(ns);counts",1000,0,50);
			    H1(sqrt(bhqdcleft*bhqdcright), "QDC BH Geometric mean", "QDC BH Geometric mean;channel; counts", 4096, 0,4096);
			   //	H1(RF,"RF","RF;time(ns)",857,0,19.75);
			     H2(RF,sqrt(bhqdcleft*bhqdcright),"QDC V RF BH","QDC V RF BH;time (ns);channel",857,0,19.75,4096,0,4096);

					 	 H1(bhqdcleft,"BH QDC left"," BH QDC left; Channel No.;counts",4096,0,4096);
						 H1(bhqdcright,"BH QDC right"," BH QDC right; Channel No.;counts",4096,0,4096);
					 }
				     }
				 }
			 }
			

			// if (RF > rfcut[bhplane-2][0] && RF < rfcut[bhplane-2][1])
			 	H1(TOF,"TOFee","TOF Electron;time(ns);counts",1000,0,50);	
			// if (RF > rfcut[bhplane-2][2] && RF < rfcut[bhplane-2][3])
			 	H1(TOF,"TOFpi","TOF Pion;time(ns);counts",1000,0,50);
			// if (RF > rfcut[bhplane-2][4] && RF < rfcut[bhplane-2][5])
			 	H1(TOF,"TOFmu","TOF Pion;time(ns);counts",1000,0,50);

			for(auto bmqdcup:bm->adc_mqdc[1])
			{
				for(auto bmqdcdown:bm->adc_mqdc[0])
			  {
			    H2(TOF,sqrt(bmqdcdown*bmqdcup),"QDC V TOF","QDC V TOF;time (ns);channel",1000,0,50,4096,0,4096);
			    //H2(RF,sqrt(bmqdcdown*bmqdcup),"QDC V RF","QDC V RF;time (ns);channel",857,0,19.75,4096,0,4096);
			  }
			}  

		}
	}

	return;
}

void TOF::TOFResolution(double reftimeSPS,double reftimeBH, BHbar *bh1,BHbar *bh2,BHbar *rf) //does this need to be here ???
{
	auto bhleft = findSmallestTime(bh1->tdc_trb[0]);
	auto bhright = findSmallestTime(bh1->tdc_trb[1]);
	auto spsup = findSmallestTime(bh2->tdc_trb[0]);
	auto spsdown = findSmallestTime(bh2->tdc_trb[1]);
	auto rftime = findSmallestTime(rf->tdc_trb[0]);

	double window[6] = {7.,9.,4.,6.5,0.,2.};//e, mu, pi
	double windowtime = fmod((bhleft->second.time+bhright->second.time)/2. - rftime->second.time+1000,19.75);
	H1(fmod((bhleft->second.time+bhright->second.time)/2. - rftime->second.time+1000,19.75),"RF bh","RF bh",806,0,19.75);
	H1(fmod((spsup->second.time+spsdown->second.time)/2. - rftime->second.time+1000,19.75),"RF SA","RF SA",806,0,19.75);
	H1((bhleft->second.time+bhright->second.time)/2. - (spsup->second.time+spsdown->second.time)/2.,"TOF Between bh 50 and SA","TOF Between bh 50 and SA",4122,-100,0);
	return;
}

Long_t TOF::defineHistograms()
{
	for (int i = 0; i < 4; i ++)
	{
		TOFperbar[0][i] = dH2(TString::Format("TOF Plots/TOF per bar BM%d V Plane C", i),TString::Format("TOF per bar BM%d V Plane C;Bar Number;TOF (ns)", i),16,-0.5,15.5,1000,0,50);
		TOFperbar[1][i] = dH2(TString::Format("TOF Plots/TOF per bar BM%d V Plane D", i),TString::Format("TOF per bar BM%d V Plane D;Bar Number;TOF (ns)", i),16,-0.5,15.5,1000,0,50);
	
		// for (int j = 0; j < 16; j ++)
		// 	TOFee[i][j] = dH1(TString::Format("TOF Plots/BM %d/BH %d e TOF",i,j),TString::Format("BH %d e TOF;TOF(ns)",j),400,18,22);
	}

	// planeTOF = dH1("TOF Plots/Plane TOF","Plane TOF (Sum of BM - Sum of BH);TOF (ns)", 1000, 0, 50);

	return ok;
}

Long_t TOF::startup()
{

	bhraw = NULL;
  	getBranchObject("BH",(TObject **) &bhraw);
  	if (!bhraw) {
	  debug(0,"Could not find bh tree in file\n");
  	}

	spsraw = NULL;
  	getBranchObject("SPS",(TObject **) &spsraw);
  	if (!spsraw) {
	  debug(0,"Could not find SPS tree in file\n");
  	}

  	bmraw = NULL;
  	getBranchObject("BM",(TObject **) &bmraw);
  	if (!bmraw)
  		debug(0,"Could not find Beam Monitor in file\n");

	/*
	vetoraw = NULL;
	getBranchObject("VETO",(TObject **) &vetoraw);
	if (!bmraw)
	  debug(0, "Could not find VETO in file\n");
	*/
/*

	double rfpeak[2][16]= {
    {6.48906e+00, 8.06902e+00, 5.56090e+00, 8.21307e+00, 5.36649e+00, 6.89964e+00, 6.95255e+00, 6.23303e+00, 7.33977e+00, 8.09273e+00, 7.99242e+00, 8.35410e+00, 4.22339e+00, 5.28696e+00, 6.21367e+00, 6.35912e+00},
    {6.77723e+00, 8.62393e+00, 5.84697e+00, 8.90379e+00, 5.76910e+00, 7.51783e+00, 7.60129e+00, 6.73340e+00, 7.95756e+00, 8.93531e+00, 9.11514e+00, 9.22474e+00, 4.67324e+00, 6.23261e+00, 6.82969e+00, 7.15897e+00}
  };

  	theRunInfo=(MRTRunInfo *) getFileObject("RunInfo");
	/* **
  	 for (int i = 0; i < 2; i ++)
  	 {
  	 	for (int j = 0; j < 16; j ++)
  	 		rfcor[i][j] = rfcenter - rfpeak[i][j];
  	 }
	 

/*  	for (int i = 0; i < 16; i ++)
  	{
  		rfcor[0][i] = rfcenter - rfpeak[0][i];
  		rfcor[1][i] = rfcenter - rfpeak[1][i] + planedist;
  	}

  	for (int i = 0; i < 4; i ++)
  		BMrfcor[i] = rfcenter - BMrfpeak[i];

  	for (int i = 0; i < 4; i ++)
  	{
  		for (int j = 0; j < 16; j ++)
  		{
  			tofcenter[i] = tofpeak[i][8];
  			tofcor[i][j] = tofcenter[i] - tofpeak[i][j];
  		}
  	}

	return ok;
}

Long_t TOF::process()
{

  	if((bhraw->trb_reftime.size()==0) ||(spsraw->trb_reftime.size()==0) )
    {
      //debug(0,"No reftime found: %i %i\n",bhraw->trb_reftime.size(),spsraw->trb_reftime.size());
        return ok;
    }
	reftimeBH= findReference(bhraw->trb_reftime)->time;
	reftimeBM= findReference(bmraw->trb_reftime)->time;
	reftimeSPS= findReference(spsraw->trb_reftime)->time;
	//reftimeVETO= findReference(vetoraw->trb_reftime)->time;
	//this is where shit is fucked up//
	trigger = findTrig(spsraw->trig)->time;

	/*
	for (int i=0; i<vetoraw->plane.size();i++)
	  if(vetoraw->plane[i].tdc_trb.size()>0)
	    return ok;
	*/  
	
  //reftimebh = 0;
  //reftimeSPS = 0;

  // //  debug(0,"Sizes %i %i\n",bhraw->trb_reftime.size(),spsraw->trb_reftime.size());

  // //H1(reftimebh-reftimeSPS,"reftime diff","reftime diff",1000,-1,-1);
  // //plotST(reftimeSPS,&bhraw->plane[0][0],&spsraw->wall[0][0],&bhraw->plane[0][2], &spsraw->wall[0][1],&spsraw->wall[0][2]);

	// cd("TOF Plots");
	// for(int i = 0; i<4;i++)
	// {
	// 		cd(TString::Format("SPS Wall %i",i));
	// 		for (int j=0;j<spsraw->wall[i].size();j++)
	// 		{
	// 			if(i==1&&j==10)
	// 			{
	// 				for(int k=0; k<2;k++)//only look at first two bh planes
	// 				{
	// 					cd(TString::Format("bh Plane %i",k));
	// 					for(int l=0;l<bhraw->plane[k].size();l++)
	// 					{
	// 						cd(TString::Format("TOF bh paddle %i SPSbar %i",l,j));
	// 						//TOFall(&bhraw->plane[k][l],&spsraw->wall[i][j]);
	// 						cd("..");
	// 					}
	// 					cd("..");
	// 				}
	// 			}
	// 		}
	// 	cd("..");
 //    }
 //    cd("..");


	///Commented out old code in favor of Win's new code////
	/*
   cd("Beam Monitor");
	for(int i = 0; i<2;i++)
	{
			cd(TString::Format("BM Wall %i",i));
			for (int j=0;j<bmraw->plane[i].size();j++)
			{
			  	if(i==1) //the plane with Steffen's big bars
			  	{
					cd(TString::Format("BMbar %i",j));
					for(int k=2; k<4;k++)//only look at last two bh planes
					{
						cd(TString::Format("%s",planenames[k].c_str()));
						for(int l=0;l<bhraw->plane[k].size();l++)
						{
						  // if(l<13) //combined with for loop above, only looks at BH bars from 5 to 12
						  //	{
								cd(TString::Format("TOF BH paddle %i",l));
							       	TOFBM(&bhraw->plane[k][l],&bmraw->plane[i][j], j, k, l);
								// TOFBM(&bhraw->plane[k][l],&bmraw->plane[i][j]);
								cd("..");
							//	}
   //TOFResolution(reftimeSPS,reftimebh,&bhraw->plane[0][0],&bhraw->plane[0][3], &bhraw->plane[0][14]);
						}
						cd("..");
					}
					cd("..");
						}
			}
		cd("..");
    }

   cd("..");
	*/

	///Win's new code starts here

	/*

bool goodevent = true;
 bool goodevent1 = true;
 bool goodevent2 = true;


cd("Beam Monitor");

for(int i = 0; i<2;i++)

{

cd(TString::Format("BM Wall %i",i));

for (int j=0;j<bmraw->plane[i].size();j++)

{

if(i==1) //the plane with Steffen's big bars

{

if ( bmraw->plane[i][j].tdc_trb[0].size() > 0 && bmraw->plane[i][j].tdc_trb[1].size() > 0 )

{

for ( int check = j+1; check < bmraw->plane[i].size(); check ++ )

{

if ( bmraw->plane[i][check].tdc_trb[0].size() > 0 && bmraw->plane[i][check].tdc_trb[1].size() > 0 )

goodevent = false;

}


for(int k=2; k<3;k++)//only look at last two bh planes

{

for(int l=0;l<bhraw->plane[k].size();l++)

{

if ( bhraw->plane[k][l].tdc_trb[0].size() > 0 && bhraw->plane[k][l].tdc_trb[1].size() > 0 )

{

for ( int check = l+1; check < bhraw->plane[k].size(); check ++ )

{

if ( bhraw->plane[k][check].tdc_trb[0].size() > 0 && bhraw->plane[k][check].tdc_trb[1].size() > 0 )

goodevent1 = false;

}
 

 for(int k=3; k<4;k++)//only look at last two bh planes

{

for(int l=0;l<bhraw->plane[k].size();l++)

{

if ( bhraw->plane[k][l].tdc_trb[0].size() > 0 && bhraw->plane[k][l].tdc_trb[1].size() > 0 )

{

for ( int check = l+1; check < bhraw->plane[k].size(); check ++ )

{

if ( bhraw->plane[k][check].tdc_trb[0].size() > 0 && bhraw->plane[k][check].tdc_trb[1].size() > 0 )

goodevent2 = false;

}


 if (goodevent && (goodevent1 || goodevent2))

{

cd(TString::Format("BMbar %i",j));

cd(TString::Format("%s",planenames[k].c_str()));

cd(TString::Format("TOF BH paddle %i",l));

TOFBM(&bhraw->plane[k][l],&bmraw->plane[i][j], j, k, l);

cd("..");

cd("..");

cd("..");

}

}

}
 }
 }
 }

}

}

}

}

cd("..");

}

	
	///my code starts here ////
	

	cd("Beam Monitor");

	bool goodevent1 = true;
	bool goodevent2 = true;
	auto goodcount1 = 0;
	auto goodcount2 = 0;
	

	for(int i = 0; i<2;i++)
	  {
	    cd(TString::Format("BM Wall %i",i));
	    for (int j=0;j<bmraw->plane[i].size();j++)
	      {
		if(i==1) //the plane with Steffen's big bars
		   {
		 if ( bmraw->plane[i][j].tdc_trb[0].size() > 0 && bmraw->plane[i][j].tdc_trb[1].size() > 0 )
		    {
		       goodevent1 = true;
		       goodcount1 = 0; 
			for ( int check = 0; check < bmraw->plane[i].size(); check ++ )
			  {
			    if (bmraw->plane[i][check].tdc_trb[0].size() > 0 && bmraw->plane[i][check].tdc_trb[1].size() > 0 )
			      goodcount1++;
			  }
			if(goodcount1 > 1)
			  goodevent1=true;
			  

			for(int k=2; k<4;k++)//only look at last two bh planes
			  {
			    for(int l=0;l<bhraw->plane[k].size();l++)
			      {
				if ( bhraw->plane[k][l].tdc_trb[0].size() > 0 && bhraw->plane[k][l].tdc_trb[1].size() > 0 )
				{
				  goodevent2 = true;
				  goodcount2 = 0;

				  /* **
				  auto bhleft = findSmallestTime(bhraw->plane[k][l].tdc_trb[0]);
				  auto bhright = findSmallestTime(bhraw->plane[k][l].tdc_trb[1]);
				  double RF = fmod((bhleft->second.time+bhright->second.time)/2-reftimeBH+rfoffset+rfcor[k-2][l],19.75); 
				  H2(l,RF,TString::Format("RF Time per bar %s",planenames[k].c_str()),TString::Format("RF Time per bar, %s;bar number;ns",planenames[k].c_str()),16,-0.5,15.5,857,0,19.75);
    			H1(RF,TString::Format("total RF %s",planenames[k].c_str()),TString::Format("Avg of All %s RF ;time(ns)",planenames[k].c_str()),857,0,19.75);

				  
				    for ( int check = 0; check < bhraw->plane[k].size(); check ++ )
				      {
					if (bhraw->plane[k][check].tdc_trb[0].size() > 0 && bhraw->plane[k][check].tdc_trb[1].size() > 0)
					 goodcount2++;
				      }
				    if(goodcount2 > 1)
				      goodevent2 = true;
				   

				    if (goodevent1 && goodevent2)
				      {
					cd(TString::Format("BMbar %i",j));
					cd(TString::Format("%s",planenames[k].c_str()));
					cd(TString::Format("TOF BH paddle %i",l));
					TOFBM(&bhraw->plane[k][l],&bmraw->plane[i][j], j, k, l);
					cd("..");
					cd("..");
					cd("..");
				      }
			
				
				}
			      }
			  }
		    }
		  }
	      }
	    cd("..");
	  }
	
	////////////

	

   // cd("..");
/*

    for (int k = 2; k < 4; k++)
    {
    	for(int i = 0; i < 16; i++)
    	{
    		auto bh1 = findSmallestTime(bhraw->plane[k][i].tdc_trb[0]);
    		auto bh2 = findSmallestTime(bhraw->plane[k][i].tdc_trb[1]);

    		cd("Correlations");
    		for(int j = 0; j < 4; j++)
    		{
    			auto bm1 = findSmallestTime(bmraw->plane[1][j].tdc_trb[0]);
    			auto bm2 = findSmallestTime(bmraw->plane[1][j].tdc_trb[1]);
    			if(bh1->second.rising&&bm1->second.rising&&bh1!=bhraw->plane[k][i].tdc_trb[0].end()&&bm1!=bmraw->plane[1][j].tdc_trb[0].end())
    			{
    				H2(i,j,TString::Format("%s to BM Correlation Down PMT",planenames[k].c_str()),TString::Format("%s to BM Correlation Down PMT;bh paddle;bm paddle",planenames[k].c_str()),16,-0.5,15.5,4,-0.5,3.5);
    				if(bh2->second.rising&&bm2->second.rising&&bh2!=bhraw->plane[k][i].tdc_trb[1].end()&&bm2!=bmraw->plane[1][j].tdc_trb[1].end())
    					H2(i,j,TString::Format("%s to BM Correlation PMT",planenames[k].c_str()),TString::Format("%s to BM Correlation Both PMT;bh paddle;bm paddle",planenames[k].c_str()),16,-0.5,15.5,4,-0.5,3.5);
    			}
    			if(bh2->second.rising&&bm2->second.rising&&bh2!=bhraw->plane[k][i].tdc_trb[1].end()&&bm2!=bmraw->plane[1][j].tdc_trb[1].end())
    				H2(i,j,TString::Format("%s to BM Correlation Up PMT",planenames[k].c_str()),TString::Format("%s to BM Correlation Up PMT;bh paddle;bm paddle",planenames[k].c_str()),16,-0.5,15.5,4,-0.5,3.5);
    		}

    		for(int j = 0; j < 32; j++)
    		{
    			auto bm1 = findSmallestTime(bmraw->plane[0][j].tdc_trb[0]);
    			auto bm2 = findSmallestTime(bmraw->plane[0][j].tdc_trb[1]);
    			if(bh1->second.rising&&bm1->second.rising&&bh1!=bhraw->plane[k][i].tdc_trb[0].end()&&bm1!=bmraw->plane[0][j].tdc_trb[0].end())
    			{
    				if(bh2->second.rising&&bm2->second.rising&&bh2!=bhraw->plane[k][i].tdc_trb[1].end()&&bm2!=bmraw->plane[0][j].tdc_trb[1].end())
    					H2(i,j,TString::Format("%s to BM Hodoscope Correlation PMT",planenames[k].c_str()),TString::Format("%s to BM Hodoscope Correlation Both PMT;bh paddle;bm paddle",planenames[k].c_str()),16,-0.5,15.5,32,-0.5,31.5);
    			}

    		}
    		cd("..");

    		cd("RF plots");
    		if(bh1->second.rising&&bh1!=bhraw->plane[k][i].tdc_trb[0].end()&&bh2->second.rising&&bh2!=bhraw->plane[k][i].tdc_trb[1].end())
    		{
    			double RF = fmod((bh1->second.time+bh2->second.time)/2-reftimeBH+rfoffset+rfcor[k-2][i],19.75); 
    			H2(i,RF,
			   TString::Format("RF Time per bar %s",planenames[k].c_str()),TString::Format("RF Time per bar, %s;bar number;ns",planenames[k].c_str()),16,-0.5,15.5,857,0,19.75);
    			H1(RF,TString::Format("total RF %s",planenames[k].c_str()),TString::Format("Avg of All %s RF ;time(ns)",planenames[k].c_str()),857,0,19.75);
    		}
    		cd("..");
    	}
    }

    cd ("RF plots");
    for (int j = 0; j < 4; j ++)
    {
    	auto bm1 = findSmallestTime(bmraw->plane[1][j].tdc_trb[0]);
    	auto bm2 = findSmallestTime(bmraw->plane[1][j].tdc_trb[1]);
    	double bmRF = fmod((bm1->second.time+bm2->second.time)/2-reftimeBM+rfoffset+BMrfcor[j],19.75); 
    	if(bm1->second.rising&&bm1!=bmraw->plane[1][j].tdc_trb[0].end()&&bm2->second.rising&&bm2!=bmraw->plane[1][j].tdc_trb[1].end())
    		H2(j,bmRF,"BM RF Time per bar","BM RF Time per bar;bar number;ns",4,-0.5,3.5,857,0,19.75);
    }
    cd ("..");

    cd ("RF plots");
    for (int i = 0; i < 16; i++)
    {
    	auto bh1up = findSmallestTime(bhraw->plane[2][i].tdc_trb[0]);
    	auto bh1dn = findSmallestTime(bhraw->plane[2][i].tdc_trb[1]);

    	if(bh1up->second.rising&&bh1up!=bhraw->plane[2][i].tdc_trb[0].end()&&bh1dn->second.rising&&bh1dn!=bhraw->plane[2][i].tdc_trb[1].end())
    	{
    		for (int j = 0; j < 16; j++)
    		{
    			auto bh2up = findSmallestTime(bhraw->plane[3][j].tdc_trb[0]);
    			auto bh2dn = findSmallestTime(bhraw->plane[3][j].tdc_trb[1]);    			
    			
    			if (bh2up->second.rising&&bh2up!=bhraw->plane[3][j].tdc_trb[0].end()&&bh2dn->second.rising&&bh2dn!=bhraw->plane[3][j].tdc_trb[1].end())
    			{
    				H2(i,j,"BH Correlation","BH Correlation;BH1;BH2",16,-0.5,15.5,16,-0.5,15.5);
    				// if (countrfc==1&&countrfd==1)
    				{
    					double rfC = fmod((bh1up->second.time+bh1dn->second.time)/2-reftimeBH+rfoffset+rfcor[0][i],19.75); 
    					double rfD = fmod((bh2up->second.time+bh2dn->second.time)/2-reftimeBH+rfoffset+rfcor[1][j],19.75);
    					H1(rfC,"rfC","rfC;ns;counts",857,0,19.75);
    					H1(rfD,"rfD","rfD;ns;counts",857,0,19.75);
    					H1(rfC-rfD,"Plane RF Diff","Plane C RF - Plane D RF; ns",100,-5,5);
    					//	 H1((rfcomp[0][i]+rfcomp[1][j])/2.,"Avg RF of C+D","Avg RF of C+D;RF (ns)",857,0,19.75);
					//	H1((rfC+rfD)/2.,"Avg RF of C+D","Avg RF of C+D;RF (ns)",857,0,19.75);//
    					H2(i,rfC-rfD,"RF Diff Per Bar","Plane C RF - Plane D RF Per Bar;Plane D Bar Number; ns",16,-0.5,15.5,100,-5,5);
    					if (i == 8)
					  	{
    						H2(j,rfC-rfD,"C8 - RF Per Bar","Plane C Bar 8 RF - Plane D RF Per Bar;Bar Number; ns",16,-0.5,15.5,100,-5,5);
    					//	H1(rfC-rfD,"C8 - RF","Plane C Bar 8 RF - Plane D RF; ns", 100, -5,5);
    					}
    				}

    				// quick test on plane tof
    				for (int k = 0; k < 4; k ++)
    				{
    					auto bm1 = findSmallestTime(bmraw->plane[1][k].tdc_trb[0]);
    					auto bm2 = findSmallestTime(bmraw->plane[1][k].tdc_trb[1]);
    					if ((bm1->second.rising&&bm1!=bmraw->plane[1][k].tdc_trb[0].end()&&bm2->second.rising&&bm2!=bmraw->plane[1][k].tdc_trb[1].end()))
    					{
					  	double planeCtime = (bh1up->second.time+bh1dn->second.time)/2.+rfcor[0][i];
    						double planeDtime = (bh2up->second.time+bh2dn->second.time)/2.+rfcor[1][j]+tofcor[k][j];
    						double avgCDtime  = (planeCtime+planeDtime)/2.;
    						double BMtime = (bm1->second.time+bm2->second.time)/2+BMrfcor[k];
							H1(BMtime-avgCDtime,"Avg TOF","Avg TOF; ns; counts",1000,0,50);//
    						H2(k,BMtime-avgCDtime,"TOF per BM plane","TOF per BM plane; Bar Number; TOF",4,-0.5,3.5,1000,0,50);
    					}

    				}
    			}
    			
    		}
    	}
    }
    
    cd ("..");
    double offsetsD[16] = {-0.047,1.746,-0.997,2.079,-0.963,0.841,0.888,0.0,1.208,2.147,2.286,2.309,-2.296,-0.7963,-0.216,0.183};
  double offsetsC[16] = {0.593,-0.984,1.521,-1.295,1.558,0.0283,0.0,0.652,-0.397,-1.116,-1.011,-1.386,2.725,1.517,0.940,0.637};
  int i = 1;
    	for(int j = 0; j<16; j++)
    		for(int k = 0; k<16; k++)
    		{
    			auto bm1 = findSmallestTime(bmraw->plane[1][i].tdc_trb[0]);
    			auto bm2 = findSmallestTime(bmraw->plane[1][i].tdc_trb[1]);
    			
    			auto bhc1 = findSmallestTime(bhraw->plane[2][j].tdc_trb[0]);
    			auto bhc2 = findSmallestTime(bhraw->plane[2][j].tdc_trb[1]);

    			auto bhd1 = findSmallestTime(bhraw->plane[3][k].tdc_trb[0]);
    			auto bhd2 = findSmallestTime(bhraw->plane[3][k].tdc_trb[1]);

    			if(bm1->second.rising&&bm2->second.rising&&bm1!=bmraw->plane[1][i].tdc_trb[0].end()&&bm2!=bmraw->plane[1][i].tdc_trb[1].end())
    				if(bhc1->second.rising&&bhc2->second.rising&&bhc1!=bhraw->plane[2][j].tdc_trb[0].end()&&bhc2!=bhraw->plane[2][j].tdc_trb[1].end())
    				{
    					if(bhd1->second.rising&&bhd2->second.rising&&bhd1!=bhraw->plane[3][k].tdc_trb[0].end()&&bhd2!=bhraw->plane[3][k].tdc_trb[1].end())
    						{
    							double BHAVG = ((bhc1->second.time+bhc2->second.time)/2 +offsetsC[j] + (bhd1->second.time+bhd2->second.time)/2 + offsetsD[k])/2;
    							double BHBMTOF = (bm1->second.time+bm2->second.time)/2 - BHAVG;
    							H1(BHAVG-trigger,"BHAVG - Trigger","BHAVG - Trigger",400,-200,-160);
    							double BHTOF = (bhc1->second.time+bhc2->second.time)/2 +offsetsC[j] - (bhd1->second.time+bhd2->second.time)/2 + offsetsD[k];
    							H2(i,BHBMTOF,"TOF Between Avg BH and BM", "TOF Between AVG BH and BM",4,-0.5,3.5,10000,-100,140);
    							H1(BHTOF,"BH TOF","BH TOF",100,-3,3);
    							for(auto hitadc:bmraw->plane[1][2].adc_mqdc[1])
    								for(auto hitadc1:bmraw->plane[1][2].adc_mqdc[0])
    									H2(BHBMTOF,sqrt(hitadc1*hitadc),"QDC V TOF","QDC V TOF",10000,-100,140,2000,0,2000);
    						}
    				}
    		}

    auto c1 = findSmallestTime(bhraw->plane[2][6].tdc_trb[0]);
    auto c2 = findSmallestTime(bhraw->plane[2][6].tdc_trb[1]);
	if(c1->second.rising&&c2->second.rising&&c1!=bhraw->plane[2][6].tdc_trb[0].end()&&c2!=bhraw->plane[2][6].tdc_trb[1].end())
		{
			H1((c1->second.time+c2->second.time)/2 - trigger,"Plane C bar 6 - trigger","Plane C bar 6 - trigger",400,-200,-160);
			for(auto hitadc:bmraw->plane[1][2].adc_mqdc[1])
				for(auto hitadc1:bmraw->plane[1][2].adc_mqdc[0])
				{
				  H2((c1->second.time+c2->second.time)/2 - trigger,sqrt(hitadc*hitadc1),"BM QDC v Plane c bar 6","BM QDC v Plane c bar 6",400,-200,-160,2000,0,2000);
    			}
    	}
   //TOFResolution(reftimeSPS,reftimebh,&bhraw->plane[0][0],&bhraw->plane[0][3], &bhraw->plane[0][14]);

  return ok;
}


Long_t TOF::finalize()
{
	// for ( int i = 0; i < 4; i++ )
	// {
	// 	printf("BM %d e TOF peak:\n", i);
	// 	for ( int j = 0; j < 16; j++ )
	// 	{
	// 		TF1 *f1 = new TF1("f1","gaus",18,22);
	// 		TOFee[i][j]->Fit(f1,"Q");
	// 		printf("%f, ", f1->GetParameter(1));
	// 	}
	// 	printf("\n");
	// }
		
	return ok;
}

Long_t TOF::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new TOF(in,out,inf_,outf_,p);
}
}


ClassImp(TOF);
*/
