#include <GEMTOF.h>

#include<iostream>
#include<cmath>


GEMTOF::GEMTOF(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

GEMTOF::~GEMTOF()
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
void GEMTOF::TOFBM(BHbar *bh, BMbar *bm, int bmbar, int bhplane, int bhbar)
{

	auto bmqdcup= bm->adc_mqdc[1];
	auto bmqdcdown= bm->adc_mqdc[0];
	auto bmtdcdown= bm->tdc_trb[0];
	auto bmtdcup= bm->tdc_trb[1];

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

			double bmdownmult = bm->tdc_trb[0].size();
			double bmupmult = bm->tdc_trb[1].size();
			double RF = fmod((bhlefttime+bhrighttime)/2-reftimeBH+rfoffset+rfcor[bhplane-2][bhbar],19.75);

			H1(TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/TOF",bmbar,planenames[bhplane].c_str(),bhbar),"TOF;time(ns);counts",1000,0,50);
			H1(RF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/RF",bmbar,planenames[bhplane].c_str(),bhbar),"RF;time(ns)",857,0,19.75);
			H2(RF,TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BH RF V TOF",bmbar,planenames[bhplane].c_str(),bhbar),"BH RF V TOF;RF time (ns);TOF (ns)",857,0,19.75,5000,0,50);
			// if (RF > rfcut[0] && RF < rfcut[1]) // electron
			// 	H1(TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/TOFee",bmbar,planenames[bhplane].c_str(),bhbar),"TOF Electron;time(ns);counts",5000,0,50);	
			// if (RF > rfcut[2] && RF < rfcut[3]) // pion
			// 	H1(TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/TOFpi",bmbar,planenames[bhplane].c_str(),bhbar),"TOF Pion;time(ns);counts",5000,0,50);
			// if (RF > rfcut[4] && RF < rfcut[5]) // muon
			// 	H1(TOF,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/TOFmu",bmbar,planenames[bhplane].c_str(),bhbar),"TOF Muon;time(ns);counts",5000,0,50);


			if ( bm->adc_mqdc[1].size() > 0 && bm->adc_mqdc[0].size() > 0) /// H2 with large number of bins will cause memory leak.
			{
				auto bmqdcup = bm->adc_mqdc[1][0];
				auto bmqdcdown = bm->adc_mqdc[0][0];
				H1(sqrt(bmqdcdown*bmqdcup),TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BM QDC GeoMean",bmbar,planenames[bhplane].c_str(),bhbar),"BM QDC Spectrum;Channel",4096,0,4096);
				
				// BM0_QCDvTOF[bhplane][bhbar]->Fill(TOF,sqrt(bmqdcdown*bmqdcup));
				// BM0_QCDvRF[bhplane][bhbar]->Fill(RF,sqrt(bmqdcdown*bmqdcup));
			}

			if ( bh->adc_mqdc[1].size() > 0 && bh->adc_mqdc[0].size() > 0 )
			{
				auto bhqdcup = bh->adc_mqdc[1][0];
				auto bhqdcdown = bh->adc_mqdc[0][0];
				H1(sqrt(bhqdcdown*bhqdcup),TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BH QDC GeoMean",bmbar,planenames[bhplane].c_str(),bhbar),"BH QDC Spectrum;Channel",4096,0,4096);
				// H2(TOF,sqrt(bhqdcdown*bhqdcup),TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BH QDC V TOF",bmbar,planenames[bhplane].c_str(),bhbar),"BH QDC V TOF;time (ns);channel",5000,0,50,4096,0,4096);
				// H2(RF,sqrt(bhqdcdown*bhqdcup),TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/BH QDC V RF",bmbar,planenames[bhplane].c_str(),bhbar),"BH QDC V RF;time (ns);channel",857,0,19.75,4096,0,4096);	
			}

			/// GEM Track 
			H1(track_mx,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/Track/mx",bmbar,planenames[bhplane].c_str(),bhbar),"mx",100,-1.,1.);
			H1(track_my,TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/Track/my",bmbar,planenames[bhplane].c_str(),bhbar),"my",100,-1.,1.);
			for ( int t = 0; t < 4; t ++ )
			{
				H1(track_x[t],TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/Track/x%i",bmbar,planenames[bhplane].c_str(),bhbar,t),"x0",1000,-50,50);
				H1(track_y[t],TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/Track/y%i",bmbar,planenames[bhplane].c_str(),bhbar,t),"y0",1000,-50,50);
				H1(track_z[t],TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/Track/z%i",bmbar,planenames[bhplane].c_str(),bhbar,t),"z0",1000,-50,50);
				H2(track_x[t],track_y[t],TString::Format("Beam Monitor/BM Wall 0/BMbar %i/%s/TOF BH paddle %i/Track/y%i vs x%i",bmbar,planenames[bhplane].c_str(),bhbar,t,t),"y0 vs x0",1000,-50,50,1000,-50,50);
			}	
		}
	}

	return;
}

Long_t GEMTOF::defineHistograms()
{
	// for (int j = 2; j < 4; j ++)
	// {
	// 	for (int k = 0; k < 16; k ++)
	// 	{
	// 		BM0_QCDvTOF[j][k] = new TH2D(TString::Format("Beam Monitor/BM Wall 0/BMbar 0/%s/TOF BH paddle %d/BH QDC V TOF",planenames[j].c_str(),k),"QDC v TOF;TOF(ns);Channel",5000,0,50,4096,0,4096);
	// 		BM1_QCDvTOF[j][k] = new TH2D(TString::Format("Beam Monitor/BM Wall 0/BMbar 1/%s/TOF BH paddle %d/BH QDC V TOF",planenames[j].c_str(),k),"QDC v TOF;TOF(ns);Channel",5000,0,50,4096,0,4096);
	// 		BM2_QCDvTOF[j][k] = new TH2D(TString::Format("Beam Monitor/BM Wall 0/BMbar 2/%s/TOF BH paddle %d/BH QDC V TOF",planenames[j].c_str(),k),"QDC v TOF;TOF(ns);Channel",5000,0,50,4096,0,4096);
	// 		BM3_QCDvTOF[j][k] = new TH2D(TString::Format("Beam Monitor/BM Wall 0/BMbar 3/%s/TOF BH paddle %d/BH QDC V TOF",planenames[j].c_str(),k),"QDC v TOF;TOF(ns);Channel",5000,0,50,4096,0,4096);

	// 		BM0_QCDvRF[j][k] = new TH2D(TString::Format("Beam Monitor/BM Wall 0/BMbar 0/%s/TOF BH paddle %d/BH QDC V RF",planenames[j].c_str(),k),"QDC V RF;RF(ns);Channel",857,0,19.75,4096,0,4096);
	// 		BM1_QCDvRF[j][k] = new TH2D(TString::Format("Beam Monitor/BM Wall 0/BMbar 0/%s/TOF BH paddle %d/BH QDC V RF",planenames[j].c_str(),k),"QDC V RF;RF(ns);Channel",857,0,19.75,4096,0,4096);
	// 		BM2_QCDvRF[j][k] = new TH2D(TString::Format("Beam Monitor/BM Wall 0/BMbar 0/%s/TOF BH paddle %d/BH QDC V RF",planenames[j].c_str(),k),"QDC V RF;RF(ns);Channel",857,0,19.75,4096,0,4096);
	// 		BM3_QCDvRF[j][k] = new TH2D(TString::Format("Beam Monitor/BM Wall 0/BMbar 0/%s/TOF BH paddle %d/BH QDC V RF",planenames[j].c_str(),k),"QDC V RF;RF(ns);Channel",857,0,19.75,4096,0,4096);
	// 	}
	// }
  
  
  return ok;
}

Long_t GEMTOF::startup()
{
  bhraw = NULL;
  getBranchObject("BH",(TObject **) &bhraw);
  if (!bhraw) {
    debug(0,"Could not find bh tree in file\n");
  }
  GEM_Tracks=NULL;
  getBranchObject("teletracks"    ,(TObject **) & GEM_Tracks);
  if(!GEM_Tracks)
  {
  	debug(0,"No Tracks!\n");
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

    for (int i = 0; i < 2; i ++)
  	{
  		for (int j = 0; j < 16; j ++)
  			rfcor[i][j] = alignment[j];
  	}


	return ok;
}

Long_t GEMTOF::process()
{


	if((bhraw->trb_reftime.size()==0) ||(bmraw->trb_reftime.size()==0) )
	{
	  debug(0,"No reftime found: %i %i\n",bhraw->trb_reftime.size(),bmraw->trb_reftime.size());
	  return ok;

	}

	reftimeBH= findReference(bhraw->trb_reftime)->time;
	reftimeBM= findReference(bmraw->trb_reftime)->time;
	reftimeSPS= findReference(spsraw->trb_reftime)->time;
	reftimeVETO= findReference(vetoraw->trb_reftime)->time;
	trigger = 0;

	for (int i = 0; i<vetoraw->plane.size();i++)
			if ( vetoraw->plane[i].tdc_trb.size() > 0 )
				return ok;

	H1(GEM_Tracks->tracks.size(),"track size","track size",30,0,30.);

	for (size_t j = 0 ; j < GEM_Tracks->tracks.size() ; j++)
    {

        if(GEM_Tracks -> tracks.size()>30) continue;
        // if(GEM_Tracks -> tracks.size() != 1) continue;

        track_mx = GEM_Tracks -> tracks[j].mx;
        track_my = GEM_Tracks -> tracks[j].my;
        track_x[0] = GEM_Tracks -> tracks[j].x0;
        track_y[0] = GEM_Tracks -> tracks[j].y0;      
        track_z[0] = GEM_Tracks -> tracks[j].z0;
        track_x[1] = GEM_Tracks -> tracks[j].x1;
        track_y[1] = GEM_Tracks -> tracks[j].y1;
        track_z[1] = GEM_Tracks -> tracks[j].z1;
        track_x[2] = GEM_Tracks -> tracks[j].x2;
        track_y[2] = GEM_Tracks -> tracks[j].y2;
        track_z[2] = GEM_Tracks -> tracks[j].z2;
        //I think number 3 is downstream
        //ETHAN
        track_x[3] = GEM_Tracks -> tracks[j].x3;
        track_y[3] = GEM_Tracks -> tracks[j].y3;
        track_z[3] = GEM_Tracks -> tracks[j].z3;

        if(track_mx!=-1e4 && track_my!=-1e4 && track_x[0]!=-1e4 && track_y[0]!=-1e4)
        {
        	double z = -250;//rough position of the target
            double xtarget = track_mx*(z-track_z[0])+track_x[0];// x= m*z+b, found from Analysis/PSI2018_SUmmer/muse_analysis_Sumer2018_4gem.C
            double ytarget = track_my*(z-track_z[0])+track_y[0];//y=m*z+b
            //So for a variety of Z positions that correspond to the veto entrance and target entrance window you check if
            //xtarget and ytarget are small enough to fit inside the veto or target window. If they are then plot the TOF.
            //If not then skip the event.
			// if ( track_mx > -0.1 && track_mx < 0.1 && track_my > -0.2 && track_my < 0.2 )
			{
				// if ( y0 > 0 && y0 < 35 && y1 > -5 && y1 < 25 && y2 > -5 && y2 < 30 && y3 > -20 && y3 < 35)
				{
					int i = 1; // the plane with Steffen's big bars
					for (int j=0;j<bmraw->plane[i].size();j++)
					{
						if ( bmraw->plane[i][j].tdc_trb[0].size() > 0 && bmraw->plane[i][j].tdc_trb[1].size() > 0 )
						{
							for ( int check = j+1; check < bmraw->plane[i].size(); check ++ )
							{
								if ( bmraw->plane[i][check].tdc_trb[0].size() > 0 && bmraw->plane[i][check].tdc_trb[1].size() > 0 )
									return ok; // check if only one bar gets hit
							}

							for(int k=2; k<4;k++)//only look at last two bh planes
							{
								for(int l=0;l<bhraw->plane[k].size();l++)
								{
									if ( bhraw->plane[k][l].tdc_trb[0].size() > 0 && bhraw->plane[k][l].tdc_trb[1].size() > 0 )
									{
										for ( int check = l+1; check < bhraw->plane[k].size(); check ++ )
										{
											if ( bhraw->plane[k][check].tdc_trb[0].size() > 0 && bhraw->plane[k][check].tdc_trb[1].size() > 0 )
												return ok; // check if only one bar gets hit
										}
										TOFBM(&bhraw->plane[k][l],&bmraw->plane[i][j], j, k, l);
									}	
								}

							}
						}
					}	
				}	
			}
		}

			auto bhleft = findSmallestTime(bhraw->plane[2][11].tdc_trb[0]);
			auto bhright = findSmallestTime(bhraw->plane[2][11].tdc_trb[1]);
			auto bmup = findSmallestTime(bmraw->plane[1][1].tdc_trb[0]);
			auto bmdown = findSmallestTime(bmraw->plane[1][1].tdc_trb[1]);

			double bhlefttime = bhleft->second.time;
			double bhrighttime = bhright->second.time;
			double bmdowntime = bmdown->second.time;
			double bmuptime = bmup->second.time;

			// if(bhleft!=bhraw->plane[2][11].tdc_trb[0].end()  && bhright!=bhraw->plane[2][11].tdc_trb[1].end() && bmup!=bmraw->plane[1][1].tdc_trb[0].end() && bmdown!=bmraw->plane[1][1].tdc_trb[1].end())
			// if(bhleft!=bhraw->plane[2][8].tdc_trb[0].end()  && bhright!=bhraw->plane[2][8].tdc_trb[1].end() && bmup!=bmraw->plane[1][1].tdc_trb[0].end() && bmdown!=bmraw->plane[1][1].tdc_trb[1].end())
			// {
				if(bhright->second.rising && bhleft->second.rising && bmup->second.rising && bmdown->second.rising)
				{
					double TOF = (-(bhlefttime+bhrighttime)/2+(bmuptime+bmdowntime)/2);
					// if (TOF > 21.6 && TOF < 23)
					{
						H1(track_mx,"mx","mx",100,-1.,1.);
						H1(track_my,"my","my",100,-1.,1.);
						H1(track_x[0],"x0","x0",1000,-50,50);
						H1(track_y[0],"y0","y0",1000,-50,50);
						H1(track_x[1],"x1","x1",1000,-50,50);
						H1(track_y[1],"y1","y1",1000,-50,50);
						H1(track_x[2],"x2","x2",1000,-50,50);
						H1(track_y[2],"y2","y2",1000,-50,50);
						H1(track_x[3],"x3","x3",1000,-50,50);
						H1(track_y[3],"y3","y3",1000,-50,50);
					}
				}
			// }
    }


	return ok;
}

Long_t GEMTOF::finalize()
{

	return ok;
}

Long_t GEMTOF::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new GEMTOF(in,out,inf_,outf_,p);
}
}


ClassImp(GEMTOF);

