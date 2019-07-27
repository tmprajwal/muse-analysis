#include <MUSEteleTracker.h>

//#include "cTrack.h"

#include<iostream>
#include<cmath>
#include<vector>
#include<numeric>

#include <fstream>
#include <sstream>
#include <string>
#include "TGraph.h"
#include "TF1.h"


MUSEteleTracker::MUSEteleTracker(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
}

MUSEteleTracker::~MUSEteleTracker()
{
}

Long_t MUSEteleTracker::defineHistograms()
{

	h2resmapsx=dH2(Form("MUSEtele/residualmapsx"), Form("X and Y hitmap weighted by X Residuals on MS GEM"), 50,-50,50,50,-50,50);
	h2resmapsy=dH2(Form("MUSEtele/residualmapsy"), Form("X and Y hitmap weighted by Y Residuals on MS GEM"), 50,-50,50,50,-50,50);

	trackmultiplicityUS=dH1("TeleTracks/US/trackmultiplicityUS","Track Multiplicity US GEM",10,0.5,10.5);
	trackmultiplicityMS=dH1("TeleTracks/MS/trackmultiplicityMS","Track Multiplicity MS GEM",10,0.5,10.5);
	trackmultiplicityDS=dH1("TeleTracks/DS/trackmultiplicityDS","Track Multiplicity DS GEM",10,0.5,10.5);

	return Plugin::ok;
}

Long_t MUSEteleTracker::startup()
{
  // get input branch with GEM clusters:
	clusters=NULL;
	getBranchObject("LumiGEMhits", (TObject**)&clusters);
	if (clusters==NULL)
	{
		printf(" Cannot find branch >LumiGEMhits< in input ROOT file - trying output branch\n");
		getOutBranchObject("LumiGEMhits",(TObject**)&clusters);
		if(clusters==NULL)
		{
			printf("Couldn't find any clusters in any ROOT file :(\n");
			return -1;
		}
	};
	printf(" LumiGEMhits (clusters) @%p\n", clusters);

  // create output branch with tracks:
	teletracks = new TeleTracks();
	makeBranch("teletracks", (TObject**)&teletracks);
	printf(" teletracks %p\n", teletracks);

	return Plugin::ok;
}

int event=0;

//int trks=0;
int trk1=0;
int trk2=0;
//ofstream outf("chk_gem_rotatn.dat");

Long_t MUSEteleTracker::process()
{
	event=event+1;
  //  printf("start a new event %d \n",event);

  // char leftright[2][18] = {"downstream", "upstream"};

      // vector to store all track candidates for this event:
	std::vector <StraightTrack> TrackCands;

	std::vector <int>           whichclusters[4];
	std::vector <double>        chi2;
	std::string GEMnames[] = {"US","4th","MI","DS"};

      // initialize an "empty" straight track:
	StraightTrack aTrack;
      //Target GEMs
	aTrack.x0 = -10000.0;
	aTrack.y0 = -10000.0;
	aTrack.z0 = -10000.0;
	aTrack.x1 = -10000.0;
	aTrack.y1 = -10000.0;
	aTrack.z1 = -10000.0;
	aTrack.x2 = -10000.0;
	aTrack.y2 = -10000.0;
	aTrack.z2 = -10000.0;
	aTrack.x3 = -10000.0;
	aTrack.y3 = -10000.0;
	aTrack.z3 = -10000.0;
      //IFP GEMs
	aTrack.x4 = -10000.0;
	aTrack.y4 = -10000.0;
	aTrack.z4 = -10000.0;
	aTrack.x5 = -10000.0;
	aTrack.y5 = -10000.0;
	aTrack.z5 = -10000.0;

	aTrack.mx = -10000.0;
	aTrack.my = -10000.0;
	aTrack.mxifp = -10000.0;
	aTrack.myifp = -10000.0;
	aTrack.xresidua.push_back(-10000);
	aTrack.yresidua.push_back(-10000);
	aTrack.z.push_back(-10000);
	aTrack.xchi2 = -10000.0;
	aTrack.ychi2 = -10000.0;


	TrackCands.clear();
	chi2.clear();
	teletracks->tracks.clear();

      // loop over all possible combinations of clusters:
	int combmulti=0, combmulti_cut=0;
	int gemmulti[4] = { 0, 0, 0, 0 };
	for (unsigned int g1=0; g1<clusters->hits.size(); g1++)
	{
		if (clusters->hits[g1].GEMid!=0) continue;
		gemmulti[0]++;
		for (unsigned int g2=0; g2<clusters->hits.size(); g2++)
		{
			if (clusters->hits[g2].GEMid!=1) continue;
			gemmulti[1]++;
			for (unsigned int g3=0; g3<clusters->hits.size(); g3++)
			{
				if (clusters->hits[g3].GEMid!=2) continue;
				gemmulti[2]++;
				for (unsigned int g4=0; g4<clusters->hits.size(); g4++)
				{
					if (clusters->hits[g4].GEMid!=3) continue;
					gemmulti[3]++;
				};
			};
		};
	};

      //printf("tele %d: #clusters US: %d MI: %d DS: %d  4th: %d\n", t, gemmulti[0], gemmulti[1], gemmulti[2], gemmulti[3]);
	//I think this is what stops tracking form happening on events with multiple hits
	// if (gemmulti[0]*gemmulti[1]*gemmulti[2]*gemmulti[3]>200) {
	// //	printf("for the new event  %d, %5.3lf, %5.3lf\n",event,aTrack.mx,aTrack.my);
	// 	teletracks->tracks.push_back(aTrack);
	// };


       //if (gemmulti[0]*gemmulti[1]*gemmulti[2]*gemmulti[3]>200) return Plugin::ok; //return Plugin::ok;

       H1(gemmulti[0], Form("MUSEteleTracker/Number of possible clusters US"), Form("Number of possible clusters - US GEM"),
       	11,-0.5,10.5);
       H1(gemmulti[1], Form("MUSEteleTracker/Number of possible clusters 4TH"), Form("Number of possible clusters - 4TH GEM"),
       	11,-0.5,10.5);
       H1(gemmulti[2], Form("MUSEteleTracker/Number of possible clusters MI"), Form("Number of possible clusters - MI GEM"),
       	11,-0.5,10.5);
       H1(gemmulti[3], Form("MUSEteleTracker/Number of possible clusters DS"), Form("Number of possible clusters - DS GEM"),
       	11,-0.5,10.5);


      // The only need is the relative  distances in z here.
      const std::vector<double> zgem = 
     	{ //1846.2, 1946.2, 2046.2, 2146.2,   // OLYMPUS, left sector GEMs
     	-536.0, -474.0 , -412.0, -350.0}; //Ethan and Ron measured for August 2018 Beamtime DS position
     	//	-836.0,-747.0,-712.0,-650};//Rough approximation of middle dowel position
	  
	  double dx, dy;
	  double xsloperec =0.0;
	  double ysloperec =0.0;

	  std::vector<double> x;
	  std::vector<double> y;
	  std::vector<double> xifp;
	  std::vector<double> yifp;
	  std::vector<double> zifp = {0,80};

	  x.resize(4);
	  y.resize(4);




	  bool goodtrack=false;
	  int trks=0;
            DSOFF=false; // Turn OFF DS GEM in track definitions. Used only US and MS GEMs to define the tracks and project the track to the 4TH GEM to get the residuals./////////
      for (unsigned int g1=0; g1<clusters->hits.size(); g1++) // Check for US GEM clusters
      {
	  // printf("check gem 1  %d %d \n",g1,clusters->hits[g1].GEMid);
      	if (clusters->hits[g1].GEMid!=0) continue;
      	x[0]=(clusters->hits[g1].xl*0.4-50.); 
      	y[0]=(clusters->hits[g1].yl*0.4-50.);

	   //printf("gem 1 before %d %d %5.2lf %5.2lf \n",g1,clusters->hits[g1].GEMid,x[0],y[0]);
      	if (fabs(x[0]>45)&& fabs(y[0]>45)) continue;
	  // printf("gem 1 after %d %d  %5.2lf %5.2lf \n",g1,clusters->hits[g1].GEMid,x[0],y[0]);

	  for (unsigned int g2=0; g2<clusters->hits.size(); g2++) // Check for 4th GEM clusters
	  {
	      // printf("check gem 2  %d %d \n",g2,clusters->hits[g2].GEMid);
	  	if (clusters->hits[g2].GEMid!=1) continue;
	  	x[1]=(clusters->hits[g2].xl*0.4-50.); 
	  	y[1]=(clusters->hits[g2].yl*0.4-50.);
	      // printf("gem 2 before %d %d %5.2lf %5.2lf \n",g2,clusters->hits[g2].GEMid,x[1],y[1]);
	  	if (fabs(x[1]>45)&& fabs(y[1]>45)) continue;
	      // printf("gem 2 after %d %d %5.2lf %5.2lf \n",g2,clusters->hits[g2].GEMid,x[1],y[1]);

	      for (unsigned int g3=0; g3<clusters->hits.size(); g3++) // Check for MI GEM clusters
		 //   for (unsigned int g3=gem3clust; g3<gem3clust+1; g3++)
	      {
		      // printf("bool min chi 2 again: %d %f\n",g3,minchi2fit);
		  // printf("check gem 3  %d %d \n",g3,clusters->hits[g3].GEMid);
	      	if (clusters->hits[g3].GEMid!=2) continue;
	      	x[2]=(clusters->hits[g3].xl*0.4-50.); 
	      	y[2]=(clusters->hits[g3].yl*0.4-50.);

		  //  printf("gem 3 before %d %d %5.2lf %5.2lf \n",g3,clusters->hits[g3].GEMid,x[2],y[2]);
	      	if (fabs(x[2]>45)&& fabs(y[2]>45)) continue;

		      for (unsigned int g4=0; g4<clusters->hits.size(); g4++) // Check for DS GEM clusters
		      {
		      	if (clusters->hits[g4].GEMid!=3) continue;
		      	x[3]=(clusters->hits[g4].xl*0.4-50.); 
		      	y[3]=(clusters->hits[g4].yl*0.4-50.);

 ///////////////////////////Need to deactivate DS GEM for tracking because some APVs did not work for some runs. Use the two clusters on US and MS GEMs to define the track and project it to the 4TH GEM. Then calculate the track residuals on the 4TH GEM./////////////////////////	     
		      	if (!DSOFF) {
		      		if (fabs(x[3]>45)&& fabs(y[3]>45)) continue; 
		      	};


		      	double slopezx = -10000;
		      	double slopezy = -10000;
				if ((xsloperec==slopezx)&&(ysloperec==slopezy))continue; 
				goodtrack=true;
				trks=trks+1;

				trackmultiplicityUS->Fill(trks);
				trackmultiplicityMS->Fill(trks);
				trackmultiplicityDS->Fill(trks);

				// printf("gem 3 after sloperec cut %d %5.2lf %5.2lf %5.2lf %5.2lf %f %f %f %f %d %f \n",g3,x[2],y[2],x[0],y[0],dx, dy,slopex,slopey,trks,minchi2fit);
				//This has to be done because of the way they initialize their track at the beginning
				aTrack.xresidua.clear();
				aTrack.yresidua.clear();
				aTrack.z.clear();

				whichclusters[0].push_back(g1);
				whichclusters[0].push_back(g2);
				whichclusters[0].push_back(g3);
				whichclusters[0].push_back(g4);
				//This should be removed from class definition
				aTrack.telescope = 1;		 

		   
		   //      //Quick crap function for tracks between IFP GEMs
	    //   //first GEM
				// for(unsigned int ifp1=0; ifp1<clusters->hits.size(); ifp1++)
				// {
				// 	if (clusters->hits[ifp1].GEMid!=4) continue;
				// 	xifp.push_back(clusters->hits[ifp1].xl*0.4-50.); 
				// 	yifp.push_back(clusters->hits[ifp1].yl*0.4-50.);
				// 	//if (fabs(xifp[0]>40)&& fabs(yifp[0]>40)) continue;
				// 	if(xifp[0]<-40 || xifp[0]>40 || yifp[0]>15 || yifp[0]<-15) continue; //cut on physical region with offset collimator
				// 	//second GEM
				// 	for(unsigned int ifp2=0; ifp2<clusters->hits.size(); ifp2++)
				// 	{
				// 		if (clusters->hits[ifp2].GEMid!=5) continue;
						
				// 		xifp.push_back(clusters->hits[ifp2].xl*0.4-50.); 
				// 		yifp.push_back(clusters->hits[ifp2].yl*0.4-50.);
				// 	  // if (fabs(xifp[1]>)&& fabs(yifp[1]>40)) continue;	
				// 	if(xifp[1]<-40 || xifp[1]>40 || yifp[1]>15 ||yifp[1]<-15) continue; //again cut on physical region
				// 	  //Straight ripped from tracking of normal GEMs below
				// 	  // Calculate track angles using the cluster coordinates on US and MS GEM data:
				// 	if(y.size()!=zgem.size() || x.size()!=zgem.size() )
				// 	{
				// 		debug(1,"Multiplicity > 1 per GEM. Skipping\n");
				// 		//std::cout << "y size x size z size "<< y.size() << " " << x.size() << " " << zgem.size() <<std::endl;
				// 		break;
				// 	}
				// 	double slopexifp = -10000;
				// 	double slopeyifp = -10000;
				// 	slopexifp = (xifp[1]-xifp[0])/(zifp[1]-zifp[0]); // x Slope between US and MI GEMS
				// 	slopeyifp = (yifp[1]-yifp[0])/(zifp[1]-zifp[0]); // y Slope between US and MI GEMS

				// 		//TA DA tracked as good as below
						
				// 	H1(xifp[0],"US X IFP distribution","X IFP Distribution",500,-50,50);	 
				// 	H1(yifp[0],"US Y IFP distribution","Y IFP Distribution",500,-50,50);
				// 	H2(xifp[0],yifp[0],"2D Track Map US - IFP","2D Track Map US - IFP",500,-50,50,500,-50,50);		 
				// 	H1(xifp[1],"DS X IFP distribution","X IFP Distribution",500,-50,50);	 
				// 	H1(yifp[1],"DS Y IFP distribution","Y IFP Distribution",500,-50,50);		 
				// 	H2(xifp[1],yifp[1],"2D Track Map DS - IFP","2D Track Map DS - IFP",500,-50,50,500,-50,50);		 
				// 	H1(slopexifp,"IFP x slope","IFP x slope",100,-1,1);
				// 	H1(slopeyifp,"IFP y slope","IFP y slope",100,-1,1);
				// 	H2(slopexifp,slopeyifp,"2D Slope Map IFP","2D Track Map IFP",100,-1,1,100,-1,1);		 

					//Only do fitting if one hit per GEM
					//At current stage GEMs can have ten clusters per event
					//Too much to deal with now, so only look if total clusters is 4 for Target GEMs and at least 1 in each IFP GEM
					//

					//Begin Least Squares line calculations for target GEMs
					//Averages
					//////////////////////////////////////////////////////////////////////////
					//NOTE: This code only finds the first valid hit for each GEM if there is ONLY ONE HIT IN THE GEM
					// and then attempts to find a least squares fit to those data points.
					//This means it is very sensitive to the accuracy of the cluster finder
					/////////////////////////////////////////////////////////////////////////

					double avgx = std::accumulate(x.begin(),x.end(),0.0)/x.size();
					double avgy = std::accumulate(y.begin(),y.end(),0.0)/y.size();
					double avgz = std::accumulate(zgem.begin(),zgem.end(),0.0)/zgem.size();
					//For z-y and z-x line in target GEMs
					//form of y=mz+b
					//form of x=mz+b
					double numy = 0;
					double denom = 0;
					double numx = 0;
					//Calculating the slope
					for(int i=0; i<y.size(); i++)
					{
						numy += (y[i]-avgy)*(zgem[i]-avgz);
						numx += (x[i]-avgx)*(zgem[i]-avgz);
						denom += pow((zgem[i]-avgz),2);
					}
					//Now slope and intercept and you've got your line
					slopezy = numy/denom;
					slopezx = numx/denom;
					double bzy = avgy-slopezy*avgz;
					double bzx = avgx-slopezx*avgz;
					double xchi2 = 0;
					double ychi2 = 0;

					//residuals calculation
					std::vector<double> xres;
					std::vector<double> yres;
					for(int i=0; i<y.size(); i++)
					{
						double residualx = x[i]-(slopezx*zgem[i]+bzx);
						double residualy = y[i]-(slopezy*zgem[i]+bzy);
						xres.push_back(residualx);
						xchi2 += pow(residualx/0.2,2);//Non reduced chi-sq, chose 0.2 for error as estimate
						ychi2 += pow(residualy/0.2,2);//Non reduced chi-sq, chose 0.2 for error as estimate
						yres.push_back(residualy);
						H1(residualx,TString::Format("%s Individual Hit X Residuals",GEMnames[i].c_str()),"Individual Hit X Residuals;mm",100,-3,3);
						H1(residualy,TString::Format("%s Individual Hit Y Residuals",GEMnames[i].c_str()),"Individual Hit Y Residuals;mm",100,-3,3);
					}
					ychi2 = ychi2/3;//
					xchi2 =xchi2/3;//reduced chi-sq
					H1(xchi2,"X hit reduced chi2","X hit reduced chi2",100,0,10);
					H1(ychi2,"Y hit reduced chi2","Y hit reduced chi2",100,0,10);
					//Ethan's sanity check plots
				    //Turns out I'm insane
					for(int i=0;i<y.size();i++)
					{
						if(x[i]!=-10000)
							H1(x[i],TString::Format("x%i distribution",i),"x distribution",500,-50,50);		  
						if(y[i]!=-10000)
							H1(y[i],TString::Format("y%i distribution",i),"y distribution",500,-50,50);	
						if(x[i]!=-1e4 && y[i]!=-1e4 )
							H2(x[i],y[i],TString::Format("2D Track Map GEM %i",i),TString::Format("2D Track Map GEM %i",i),500,-50,50,500,-50,50);		   
					}
					H1(slopezx,"mx distribution","mx distribution",100,-1,1);		  
					H1(slopezy,"my distribution","my distribution",100,-1,1);		
					H1(bzy,"ZY intercept","ZY intercept",100,-1,-1);
					H1(bzx,"ZX intercept","ZX intercept",100,-1,-1);

					//Since it is a track we should put in the tracks predicted hit positions
					aTrack.x0 = slopezx*zgem[0]+bzx;
					aTrack.y0 = slopezy*zgem[0]+bzy;
					aTrack.z0 = zgem[0];
					aTrack.x1 = slopezx*zgem[1]+bzx;
					aTrack.y1 = slopezy*zgem[1]+bzy;
					aTrack.z1 = zgem[1];
					aTrack.x2 = slopezx*zgem[2]+bzx;
					aTrack.y2 = slopezy*zgem[2]+bzy;
					aTrack.z2 = zgem[2];
					aTrack.x3 = slopezx*zgem[3]+bzx;
					aTrack.y3 = slopezy*zgem[3]+bzy;
					aTrack.z3 = zgem[3];
					aTrack.bx = bzx;
					aTrack.by = bzy;

					// aTrack.x4 = xifp[0];
					// aTrack.y4 = yifp[0];
					// aTrack.z4 = zifp[0];
					// aTrack.x5 = xifp[1];
					// aTrack.y5 = yifp[1];
					// aTrack.z5 = zifp[1];
					// aTrack.bxifp = xifp[0];
					// aTrack.byifp = yifp[0];

					aTrack.mx = slopezx;
					aTrack.my = slopezy;

					//aTrack.mxifp = slopexifp;
					//aTrack.myifp = slopeyifp;
					aTrack.xresidua = xres;
					aTrack.yresidua = yres;
					aTrack.xchi2 = xchi2;
					aTrack.ychi2 = ychi2;
					aTrack.z.push_back(zgem[1]);

					TrackCands.push_back(aTrack);


			  /// creat 2D residual maps ////////////
			  //	  double vec2d=sqrt((dx*dx)+(dy*dy));
			  //   outf<< x[1]<<"\t"<<y[1]<<"\t"<<dx<<"\t"<<dy<<std::endl; 		 

				    h2resmapsx->Fill(x[1],y[1],dx);
				    h2resmapsy->Fill(x[1],y[1],dy);   

				    double thischi2 = 0.0;

			  // printf("before thischi2:  %d,  %d, %d,  %d\n",g1,g2,g3,aTrack.xresidua.size());

				    for (unsigned int j=0; j<aTrack.xresidua.size(); j++){
			    		thischi2 = pow(aTrack.xresidua[j], 2.) + pow(aTrack.yresidua[j], 2.);// THIS IS NOT CHISQ!
			    //  printf("after thischi2:  %d,  %d, %d,  %d, %d, %5.3lf, %5.3lf,%5.3lf\n",g1,g2,g3,j, aTrack.xresidua.size(),aTrack.xresidua[j],aTrack.yresidua[j],thischi2);
					};

			  // printf("before thischi2:  %d, %d,  %d, %d,  %5.3lf,%5.3lf, %5.3lf\n",trks,g1,g2,g3,dx,dy,thischi2);
			  // aTrack.chi_sq=minchi2fit;
			  // if (aTrack.chi_sq <1.0) teletracks->tracks.push_back(aTrack);  // put all the track to the MUSEteletracker tree
			  		teletracks->tracks.push_back(aTrack);  // put all the track to the MUSEteletracker tree
			  // if (aTrack.chi_sq <0.4) printf("use thischi2 cut : %d, %d,  %d, %d,  %5.3lf,%5.3lf, %5.3lf\n",trks,g1,g2,g3,dx,dy,aTrack.chi_sq);

					chi2.push_back(thischi2);
					  
					xsloperec=slopezx;
					ysloperec=slopezy;

		  //TrackCands.clear(); //Ethan Commented this out because why do you clear it here if it is needed right below this to make other plots?
					//   }//IFP GEM 1

	         	//}//IFP GEM 2
			}; //DS GEM
		};//MS GEM
	    };//4TH GEM
	};//US GEM

	if (!goodtrack) {
		//teletracks->tracks.push_back(aTrack);
	//	printf("for everythng after  %d, %5.3lf, %5.3lf\n",event,aTrack.mx,aTrack.my);
	};

	H1(combmulti, Form("MUSEteleTracker/CombinationMulti"),
		Form("combination multiplicity"), 1000, -0.5, 999);
	H1(combmulti_cut, 
		Form("MUSEteleTracker/CombinationMulti_cut"),
		Form("combination multiplicity cut"), 1000, -0.5, 999);
      if (TrackCands.size()==0) return Plugin::ok; //return Plugin::ok;

      // loop over all track candidates and select the best one(s):

      int best       = 0;
      int secondbest = 0;
      double best_ch2       = 10000.000;
      double secondbest_ch2 = 10000.000;

      for (unsigned int i=1; i<TrackCands.size(); i++)
      {
      	if (chi2[i]<chi2[best])
      	{
      		secondbest = best;
      		best = i;	
      		secondbest_ch2=chi2[best];
      	}
      	else {
      		if (((chi2[i]<secondbest_ch2))) 
      		{
      			secondbest_ch2=chi2[i];
      			secondbest =i;
      		};	       
      	};	  
      };	

       ///////////////////////////////////////////////
      //  teletracks->tracks.clear();  // put the best and second best track to the MUSEteletracker tree

   //    aTrack = TrackCands[best];
   //    for (unsigned int i=0; i<TrackCands[best].xresidua.size(); i++)
   //    {
   //    	trk1=trk1+1;
   //    	aTrack.xresidua.push_back(TrackCands[best].xresidua[i]);
   //    	aTrack.yresidua.push_back(TrackCands[best].yresidua[i]);
   //    	aTrack.z.push_back(TrackCands[best].z[i]);	 
	  //    //printf("fill the best tracks  %d %d  %5.3lf, %5.3lf\n",i,trk1,TrackCands[best].xresidua[i],TrackCands[best].yresidua[i]);	  
   //    };
   //     teletracks->tracks.push_back(aTrack);  // put the best track to the MUSEteletracker tree

   //    if (secondbest!=best)  //Hits with only one track gives the same track as the second best track and counts twice. THis get rid of that and fill only the second best track among all the tracks/each event.
   //    {
   //    	trk2=trk2+1;
   //    	aTrack = TrackCands[secondbest];
   //    	aTrack.xresidua = TrackCands[secondbest].xresidua;
   //    	aTrack.yresidua = TrackCands[secondbest].yresidua;
	  //  //printf("fill the second best tracks %d \n",trk2);
	  // // teletracks->tracks.push_back(aTrack);  // put the second best track to the MUSEteletracker tree
   //    };

	///////////////////////////////////////////////////////////////////

      unsigned int thissize = teletracks->tracks.size();
      //printf(" teletracks->tracks.xresidua.size() = %d\n", teletracks->tracks[thissize-1].xresidua.size());

      x[0] = aTrack.x0;
      y[0] = aTrack.y0;
      x[1] = x[0] + aTrack.mx * (zgem[1]-zgem[0]);
      y[1] = y[0] + aTrack.my * (zgem[1]-zgem[0]);
      x[2] = x[0] + aTrack.mx * (zgem[2]-zgem[0]);
      y[2] = y[0] + aTrack.my * (zgem[2]-zgem[0]);
      x[3] = x[0] + aTrack.mx * (zgem[3]-zgem[0]);
      y[3] = y[0] + aTrack.my * (zgem[3]-zgem[0]);
      H2(x[0], y[0], Form("HitmapUS"), Form("Hitmap US GEM For the Best Track"), 100, -50., +50., 100, -50., +50.);
      H2(x[1], y[1], Form("Hitmap4TH"), Form("Hitmap 4TH GEM For the Best Track"), 100, -50., +50., 100, -50., +50.);
      H2(x[2], y[2], Form("HitmapMI"), Form("Hitmap MI GEM For the Best Track"), 100, -50., +50., 100, -50., +50.);
      H2(x[3], y[3], Form("HitmapDS"), Form("Hitmap DS GEM For the Best Track"), 100, -50., +50., 100, -50., +50.);

      //H1(chi2[best], Form("chi2"), Form("chi2 For the Best Track"), 1000., 0., 100.);
      //if (best!=secondbest)
      	//H1(chi2[secondbest], Form("2ndchi2"), Form("2ndchi2 (Chi2 For the Second Best Track"), 1000., 0., 100.);

	   x.clear();
	   y.clear();
	   xifp.clear();
	   yifp.clear();
//If these clears happen, then the output root tree is empty
//duh
//ETHAN
      // aTrack.xresidua.clear();
      // aTrack.yresidua.clear();
      // aTrack.z.clear();
	 
      // TrackCands.clear();
      // chi2.clear();
      // teletracks->tracks.clear();

      return Plugin::ok;
  }

  Long_t MUSEteleTracker::finalize()
  {
  	return Plugin::ok;

  }


  Long_t MUSEteleTracker::cmdline(char *cmd)
  {
  //add cmdline handling here

  return 0; // 0 = all ok
}


extern "C"{
	Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
	{
		return (Plugin *) new MUSEteleTracker(in,out,inf_,outf_,p);
	}
}


ClassImp(MUSEteleTracker);

