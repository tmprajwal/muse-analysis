#include <GEM_BH.h>
#include "TH2.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLegend.h"

#include <iostream>
#include <fstream>
#include <cmath>


GEM_BH::GEM_BH(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

GEM_BH::~GEM_BH()
{
};

Long_t GEM_BH::startup()
{
	Hits=NULL;
	GEM_Tracks=NULL;
    getBranchObject("teletracks"    ,(TObject **) & GEM_Tracks);
	getBranchObject("ScintHits", (TObject **) &Hits);
    theRunInfo=(MRTRunInfo *) getFileObject("RunInfo");
	if(!Hits || !GEM_Tracks)
	{
		debug(0,"Couldn't find Scinthits or GEM_Tracks in tree");
		return 0;
	}

	corrzy = dH2("Z V Y track and BH that fired","Z V Y track and BH that fired;x position (mm);z position (mm)",8*16,-52,60,600,-600,10);
	corrzx = dH2("Z V X track and BH that fired","Z V X track and BH that fired;x position (mm);z position (mm)",8*16,-52,60,600,-600,10);

    corry_yp = dH2("yp - IFP/yp v xpIFP","yp v xpIFP;xpifp;yp",200,-100,100,200,-100,100);

    for(int i=0; i<26; i++)
    {
        xpos[i] = dH1(TString::Format("Beam Profile/X position at z = %i",(-400+25*i)),"X position",400,-75,75);
        ypos[i] = dH1(TString::Format("Beam Profile/Y position at z = %i",(-400+25*i)),"Y position",400,-75,75);
        xzslope[i] = dH1(TString::Format("Beam Profile/ZX slope at z = %i",(-400+25*i)),"ZX slope",50,-0.5,0.5);
        yzslope[i] = dH1(TString::Format("Beam Profile/ZY slope at z = %i",(-400+25*i)),"ZY slope",50,-0.5,0.5);
    }

	return ok;
}

void GEM_BH::id_to_info(int id, int *plane, int *bar, int *side)
{
	*bar = id%16;
	*plane = (id/16)%4;
	*side = (id/(16*4));
	return;
}

Long_t GEM_BH::process()
{
	//if(Hits->hits.size()!=0)
	//	std::cout << "hits.size" << Hits->hits.size() << std::endl;
    double alignment[16] = {3.0,1.25,3.8,0.75,3.5,1.75,1.6,2.45,1.2,.25,.1,0,4.5,3.0,2.5,2.5};

    const int numbars=16;
    //run 395
     // double starte[8] = {8,8,7,6,8,7,6,7};
     // double stope[8] = {11,11,11,10,11,10,8,11};
     // double barcut[8] = {9,11,8,7,10,5,4,6};
     // double startpi[8] = {13,13,12,11,12,11,10,11};
     // double stoppi[8] = {16,16,15,14,16,15,13,15};

    //run380
    // double starte[11] = {7,10,9,10,6,9,9,10,9,8,10};
    // double stope[11] = {11,13,12,13,10,13,12,13,12,10,14};
    // double barcut[11] = {4,9,5,10,12,8,6,3,7,13,11};
    // double startpi[11] = {11,14,13,14,11,14,13,15,12,12,15};
    // double stoppi[11] = {15,17,16,17,13,17,16,18,16,14,18};

    //run258
    // double startmu[10] = {5,3,4,2,3,4,5,4,5,5};
    // double stopmu[10] = {9,6,7,6,5,7,9,8,8,8};
    // double starte[10] = {11,8,9,8,8,9,10,9,10,10};
    // double stope[10] = {14,12,13,11,11,13,14,13,14,13};
    // double barcut[10] = {11,7,5,4,13,6,9,8,3,10};
    // double startpi[10] = {15,13,13,12,12,13,15,13,14,14};
    // double stoppi[10] = {18,16,17,16,15,17,19,17,18,18};

     //run558-560
    // double startmu[numbars] = {14,16,14,16,14,16,16,15,16,17,17,17,12,14,14,15};
    // double stopmu[numbars] = {17,18,16,19,16,18,18,18,18,19,19,20,15,17,18,18};

    // double startpi[numbars] = {0,2,0,3,0,2,2,1,2,4,4,4,18,0,1,1};
    // double stoppi[numbars] = {4,6,3,6,3,5,5,4,5,6,6,6,2,4,4,4};

     //double barcut[numbars] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

    // double starte[numbars] = {5,7,4,7,4,6,6,5,7,8,8,8,3,5,5,6};
    // double stope[numbars] = {8,10,8,10,8,10,9,9,10,10,10,10,6,8,8,8};
    double barcut[numbars] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};

    //Runs for 161+ BH
    double startmu[numbars] = {0,2,0,3,0,1,1,0,2,3,3,3,18,0,0,0};
    double stopmu[numbars] = {3,5,2,6,2,4,5,4,5,6,6,6,1,3,3,4};

    double startpi[numbars] = {9,11,9,12,9,11,11,10,11,12,12,13,8,9,9,9};
    double stoppi[numbars] = {14,15,13,16,13,15,15,14,15,16,16,16,12,13,13,14};

    double starte[numbars] = {4,7,4,7,4,6,6,5,6,7,7,8,3,4,4,4};
    double stope[numbars] = {8,10,7,10,7,9,9,9,10,11,11,11,6,8,8,9};

    //Runs for 117+ BH
    // double startmu[numbars] = {12,14,11,14,11,13,13,12,13,14};//,14,14,10,11,11,11};
    // double stopmu[numbars] = {15,17,15,17,14,16,16,16,17,18};//,18,18,14,15,15,15};

    // double startpi[numbars] = {8,10,7,10,7,9,9,8,10,11};//,11,11,6,8,9,9};
    // double stoppi[numbars] = {11,13,11,13,10,12,12,12};//,13,14,14,14,10,11,11,11};

    // double starte[numbars] = {5,7,4,7,4,6,6,5,6,7};//,7,7,3,5,5,5};
    // double stope[numbars] = {8,10,7,10,7,9,9,8,10,11};//,11,11,6,8,9,9};

    //Runs for 210+ BH
    // double startmu[numbars] = {13,14,12,16,14,15,15,14,16,16,17,17,12,14,14,14};
    // double stopmu[numbars] = {17,19,16,19,16,18,18,18,19,19,19,19,16,17,17,17};

    // double startpi[numbars] = {1,2,0,2,0,1,2,1,2,3,3,3,19,0,0,0};
    // double stoppi[numbars] = {4,6,3,6,3,5,5,4,5,7,6,7,2,4,4,5};

    // double starte[numbars] = {5,7,4,8,4,6,6,5,7,7,8,8,3,5,5,5};
    // double stope[numbars] = {8,10,8,10,7,9,10,8,10,11,11,11,7,8,8,9};

    // double startpi =12.5;
    // double stoppi = 16;

    // double starte = 8;
    // double stope = 11;

    // double startmu = 3;
    // double stopmu = 7;



    for(size_t i = 0; i< Hits->hits.size(); i++)
    {
        int plane, side, bar;
        id_to_info(Hits->hits[i].id,&plane,&bar,&side);
        double RF = Hits->hits[i].rf;
        //H1(RF,TString::Format("RF/RF plane %i bar %i",plane,bar),TString::Format("RF plane %i bar %i;rf time(ns)",plane,bar),857,0,19.75);
        if(plane==3)
        H2(bar,RF+alignment[bar],"RF/RF plane 3","RF/RF plane 3",16,-0.5,15.5,857,0,19.75);
        if(plane==2)
        H2(bar,RF+alignment[bar],"RF/RF plane 2","RF/RF plane 2",16,-0.5,15.5,857,0,19.75);

    }

    //if(GEM_Tracks -> tracks[0].x0!=-1e4) debug(0,"Number of GEM Tracks is %d\n",GEM_Tracks->tracks.size());
    corrzy->Reset();
    corrzx->Reset();
    for (size_t j = 0 ; j < GEM_Tracks->tracks.size() ; j++)
    {

        if(GEM_Tracks -> tracks.size()>30) continue;
        double x0 = GEM_Tracks -> tracks[j].x0;
        double mx = GEM_Tracks -> tracks[j].mx;
        double y0 = GEM_Tracks -> tracks[j].y0;
        double my = GEM_Tracks -> tracks[j].my;
        double z0 = GEM_Tracks -> tracks[j].z0;
        double x1 = GEM_Tracks -> tracks[j].x1;
        double y1 = GEM_Tracks -> tracks[j].y1;
        double z1 = GEM_Tracks -> tracks[j].z1;
        double x2 = GEM_Tracks -> tracks[j].x2;
        double y2 = GEM_Tracks -> tracks[j].y2;
        double z2 = GEM_Tracks -> tracks[j].z2;
        //I think number 3 is downstream
        //ETHAN
        double x3 = GEM_Tracks -> tracks[j].x3;
        double y3 = GEM_Tracks -> tracks[j].y3;
        double z3 = GEM_Tracks -> tracks[j].z3;

        double x4 = GEM_Tracks -> tracks[j].x4;
        double y4 = GEM_Tracks -> tracks[j].y4;
        double z4 = GEM_Tracks -> tracks[j].z4;
        double x5 = GEM_Tracks -> tracks[j].x5;
        double y5 = GEM_Tracks -> tracks[j].y5;
        double z5 = GEM_Tracks -> tracks[j].z5;
        double mxifp = GEM_Tracks->tracks[j].mxifp;
        double myifp = GEM_Tracks->tracks[j].myifp;



        for(int z = z0; z < 1000; z++)
        {

			double xtrack = mxifp*z+x4;
            double ytrack = myifp*z+y4;
            if(mx!=-1e4&&my!=-1e4&&x0!=-1e4&&y0!=-1e4)
            {
                double xtarget = mx*(z-z0)+x0;// x= m*z+b, found from Analysis/PSI2018_SUmmer/muse_analysis_Sumer2018_4gem.C
                double ytarget = my*(z-z0)+y0;//y=m*z+b
                //if(z==-500)
                //std::cout << mx << " " << x0 << std::endl;
                H2(xtarget,z,"Z v X Target Region","Z v X Target Region; X position (mm); Z position (mm)",400,-75,75,1600,-600,1000);
			    H2(ytarget,z,"Z v Y Target Region","Z v Y Target Region; Y position (mm); Z position (mm)",400,-75,75,1600,-600,1000);
                if(z==-200)
                {
                    H1(mx,TString::Format("ZX Slope at z = %i",z),TString::Format("ZX Slope at z = %i",z),50,-0.5,0.5);
                    H1(my,TString::Format("ZY Slope at z = %i",z),TString::Format("ZY Slope at z = %i",z),50,-0.5,0.5);
                    H1(xtarget,TString::Format("X Position at z = %i",z),TString::Format("X Position at z = %i",z),400,-75,75);
                    H1(ytarget,TString::Format("Y Position at z = %i",z),TString::Format("Y Position at z = %i",z),400,-75,75);
                }
                //Used for beam profile measurements which are inputs to Steffen's Simulation
                for(int i = 0; i < 26; i++)
                {
                    if(z==(-500+25*i))
                    {
                        double xtarget = mx*(z-z0)+x0;// x= m*z+b, found from Analysis/PSI2018_SUmmer/muse_analysis_Sumer2018_4gem.C
                        double ytarget = my*(z-z0)+y0;//y=m*z+b
                        xpos[i]->Fill(xtarget);
                        ypos[i]->Fill(ytarget);
                        xzslope[i]->Fill(mx);
                        yzslope[i]->Fill(my);
                    }
                }
            }
            if(mxifp!=-1e4&&myifp!=-1e4&&x4!=-1e4&&y4!=-1e4 && fabs(mxifp)<0.1)
            {
                H2(xtrack,z,"Z V X IFP","Z V X IFP",400,-75,75,1600,-600,1000);
                H2(ytrack,z,"Z V Y IFP","Z V Y IFP",400,-75,75,1600,-600,1000);
            }
		}

    	//A Ron load of correlations between IFP and target
        double xtarg = -(mx*(-250)+x0)/10;//From GEM code this is supposedly the distance from last GEM to target
        double ytarg = -(my*(-250)+y0)/10;// Divide by 10 to go from mm to cm
        double xifp = -(mxifp*(40)+x4)/10;
        double yifp = -(myifp*40 +y4)/10;

        if(mx!=-1e4&&x0!=-1e4&&mxifp!=-1e4&&x4!=-1e4)
        {
            double xangle = -TMath::ATan(mx)*1000;//*180/3.1415;
            double xifpangle = -TMath::ATan(mxifp)*1000;//*180/3.1415;
            double yifpangle = -TMath::ATan(myifp)*1000;//*180/3.1415;

            H2(xangle,xifpangle,"Correlations/XpIFP V XpTarget","XpIFP V XpTarget;xptarget;xpIFP",200,-100,100,200,-100,100);
            H2(xtarg,xifp,"Correlations/X pos IFP V X pos Target","X pos IFP V X pos Target;xtarg (cm);xifp (cm)",50,-5,5,20,-4,4);
            H2(xtarg,xangle,"Correlations/Xp V X pos target","Xp V X pos target;xtarg (cm);xptarget",50,-5,5,200,-100,100);
            H2(xtarg,xifpangle,"Correlations/XpIFP V X pos targ","XpIFP V X pos targ;xtarg (cm);xpifp",50,-5,5,200,-100,100);
            H2(xangle,xifp,"Correlations/X pos IFP V Xptarg","X pos IFP V X slope targ;xptarg;xifp (cm)",200,-30,30,20,-4,4);
            H2(xifp,xifpangle,"Correlations/XpIFP V X pos IFP","XpIFP V X pos IFP;x ifp (cm);xpifp",20,-1,1,80,-40,40);//customized for turtle
            //Correlations to compare with turtle stuff


            H2(xifp,xangle,"Turtle/xp - IFP/xp v x IFP","xp v x IFP;xifp (cm);xp",20,-1,1,50,-100,100);
            H2(yifp,xangle,"Turtle/xp - IFP/xp v y IFP","xp v y IFP;yifp (cm);xp",20,-1.5,1.5,50,-100,100);
            H2(xifpangle,xangle,"Turtle/xp - IFP/xp v xpIFP","xp v xpIFP;xpifp;xp",50,-10,20,50,-100,100);
            H2(yifpangle,xangle,"Turtle/xp - IFP/xp v ypIFP","xp v ypIFP;ypifp;xp",50,-40,10,50,-100,100);
            H2(xifp,xtarg,"Turtle/x - IFP/x targ v x IFP","x targ v x IFP;xifp (cm); xtarg (cm)",20,-1,1,50,-5,5);
            H2(yifp,xtarg,"Turtle/x - IFP/x targ v y IFP","x targ v y IFP;yifp (cm); xtarg (cm)",20,-1.5,1.5,50,-5,5);
            H2(xifpangle,xtarg,"Turtle/x - IFP/x targ v xpIFP","x targ v xpIFP;xpifp; xtarg (cm)",50,-10,20,30,-5,5);
            H2(yifpangle,xtarg,"Turtle/x - IFP/x targ v ypIFP","x targ v ypIFP;ypifp; xtarg (cm)",50,-40,10,30,-5,5);

            H1(xangle,"Turtle/angles/xp","xp",100,-100,100);
            H1(xifpangle,"Turtle/angles/xpifp","xpifp",100,-100,100);
            H1(yifpangle,"Turtle/angles/ypifp","ypifp",100,-100,100);
            H1(xtarg,"Turtle/simple/xtarg","xtarg",50,-5,5);
            H1(xifp,"Turtle/simple/xifp","xifp",20,-4,4);


        }
        if(my!=-1e4&&y0!=-1e4&&myifp!=-1e4&&y4!=-1e4)
        {
            double yangle = -TMath::ATan(my)*1000;//*180/3.1415;
            double xifpangle = -TMath::ATan(mxifp)*1000;//*180/3.1415;
            double yifpangle = -TMath::ATan(myifp)*1000;//*180/3.1415;

            H2(yangle,yifpangle,"Correlations/YpIFP V YpTarget","YpIFP V YpTarget;yptarget;ypIFP",200,-100,100,200,-100,100);
            H2(ytarg,yifp, "Correlations/Y pos IFP V Y pos Target","Y pos IFP V Y pos Target;ytarg (cm);yifp (cm)",50,-5,5,20,-4,4);
            H2(ytarg,yangle,"Correlations/YpTarget V Y pos target","Yptarget V Y pos target;ytarg (cm);yptarget",50,-3,3,200,-100,100);//customized for turtle
            H2(ytarg,yifpangle,"Correlations/YpIFP V Y pos targ","YpIFP V Y pos targ;ytarg (cm);ypifp",50,-5,5,200,-100,100);
            H2(yangle,yifp,"Correlations/Y pos IFP V Yptarg","Y pos IFP V Yptarg;yptarg;yifp",200,-30,30,20,-4,4);
            H2(yifp,yifpangle,"Correlations/YpIFP V Y pos IFP","YpIFP V Y pos IFP;y ifp (cm);ypifp",20,-1.5,1.5,150,-60,60);//customized for turtle

            //Correlations to compare with turtle stuff

            H2(xifp,yangle,"Turtle/yp - IFP/yp v x IFP","yp v x IFP;xifp (cm);yp",20,-1,1,50,-40,40);
            H2(yifp,yangle,"Turtle/yp - IFP/yp v y IFP","yp v y IFP;yifp (cm);yp",20,-1.5,1.5,50,-40,40);
            H2(xifpangle,yangle,"Turtle/yp - IFP/yp v xpIFP","yp v xpIFP;xpifp;yp",50,-40,40,50,-40,40);
            H2(yifpangle,yangle,"Turtle/yp - IFP/yp v ypIFP","yp v ypIFP;ypifp;yp",50,-40,10,50,-40,40);
            H2(xifp,ytarg,"Turtle/y - IFP/y targ v x IFP","y targ v x IFP;xifp (cm); ytarg (cm)",20,-1,1,50,-3,3);
            H2(yifp,ytarg,"Turtle/y - IFP/y targ v y IFP","y targ v y IFP;yifp (cm); ytarg (cm)",20,-1.5,1.5,50,-3,3);
            H2(xifpangle,ytarg,"Turtle/y - IFP/y targ v xpIFP","y targ v xpIFP;xpIFP; ytarg (cm)",100,-10,20,50,-3,3);
            H2(yifpangle,ytarg,"Turtle/y - IFP/y targ v ypIFP","y targ v ypIFP;ypIFP; ytarg (cm)",100,-40,10,50,-3,3);
            H1(yangle,"Turtle/angles/yp","yp",100,-100,100);
            H1(ytarg,"Turtle/simple/ytarg","ytarg",50,-3,3);
            H1(yifp,"Turtle/simple/yifp","yifp",20,-1.5,1.5);


        }

        if(x0!=-1e4 && y0!=-1e4)
    		H2(x0,y0,"US GEM x V y","US GEM x V y",250,-50,50,250,-50,50);
    	if(x1!=-1e4 && y1!=-1e4)
    		H2(x1,y1,"4th GEM x V y","4th GEM x V y",250,-50,50,250,-50,50);
    	if(x2!=-1e4 && y2!=-1e4)
    		H2(x2,y2,"MI GEM x V y","MI GEM x V y",250,-50,50,250,-50,50);
    	if(x3!=-1e4 && y3!=-1e4)
    		H2(x3,y3,"DS GEM x V y","DS GEM x V y",250,-50,50,250,-50,50);
        if(x4!=-1e4 && y4!=-1e4)
            H2(x4,y4,"US-IFP GEM x V y","US-IFP GEM x V y",250,-50,50,250,-50,50);
        if(x5!=-1e4 && y5!=-1e4)
            H2(x5,y5,"DS-IFP GEM x V y","DS-IFP GEM x V y",250,-50,50,250,-50,50);

        for(size_t i = 0; i < Hits->hits.size(); i++)
        {

        	int plane, side, bar;
        	id_to_info(Hits->hits[i].id,&plane,&bar,&side);
        	double zoffsetd = z0-100;//zoffset in mm from Gem to BH plane D, guess from Ievgen
        	double zoffsetc = z0-120;//zoffset in mm from Gem to BH plane C, guess from Ievgen
        	double ypos = -10000;
        	double xpos = -10000;
        	double width = -10000;//bar width
            double RF = Hits->hits[i].rf;
            double timeleft = Hits->hits[i].timeleft;
            double timeright = Hits->hits[i].timeright;
            H1(RF,TString::Format("RF/Hit in GEMs/RF plane %i bar %i",plane,bar),TString::Format("RF plane %i bar %i;rf time(ns)",plane,bar),857,0,19.75);
        	if(plane==3&&bar<5)
        		ypos = -52 + 8*bar; // xpos is the farthest left position of BH
        	else if(plane==3 && bar>4 && bar <8)
        		ypos = -12 +4*(8-bar);
        	else if(plane ==2 && bar >7 && bar<11)
        		ypos = 4*(bar-8);
        	else if(plane==3 && bar >10)
        		ypos = 12 +8*(bar-10);

        	if(plane==3&&bar<5)
        		xpos = -52 + 8*bar; // xpos is the farthest left position of BH
        	else if(plane==3 && bar>4 && bar <8)
        		xpos = -12 +4*(8-bar);
        	else if(plane ==3 && bar >7 && bar<11)
        		xpos = 4*(bar-8);
        	else if(plane==3 && bar >10)
        		xpos = 12 +8*(bar-10);

        	if(xpos!=-1e4)
        		H2(xpos,zoffsetd,"BH Plane D Position","BH Plane D position",8*16,-52,60,150,-500,1000);
        	if(ypos!=-1e4)
        		H2(ypos,zoffsetc,"BH Plane C Position","BH Plane C position",8*16,-52,60,150,-500,1000);

        	for(int z = -600; z < 10; z++)
        	{
        		double xtarg = mx*(z-z0)+x0;// x= m*z+b, found from Analysis/PSI2018_SUmmer/muse_analysis_Sumer2018_4gem.C
        		double ytarg = my*(z-z0)+y0;//y=m*z+b
				//H2(xtarg,z,"Z V X track BH","Z V X track BH",8*16,-52,60,1500,-500,1000);
				if(ypos!=-1e4)
				{
					corrzy->Fill(ypos,zoffsetc);
					corrzy->Fill(ypos+1,zoffsetc);
					corrzy->Fill(ypos+2,zoffsetc);
					corrzy->Fill(ypos+3,zoffsetc);
					corrzy->Fill(ypos+4,zoffsetc);
    				corrzy->Fill(ytarg,z);
				}
				if(xpos!=-1e4 && mx!=-1e4 && x0!=-1e4)
				{
					corrzx->Fill(xpos,zoffsetd);
					corrzx->Fill(xpos+1,zoffsetd);
					corrzx->Fill(xpos+2,zoffsetd);
					corrzx->Fill(xpos+3,zoffsetd);
					corrzx->Fill(xpos+4,zoffsetd);
    				corrzx->Fill(xtarg,z);
				}
					//H2(xpos,zoffset,"Z V X track BH","Z V X track BH",8*16,-52,60,1500,-500,1000);
				// H2(ytarg,z,"Z V Y track","Z V Y track",400,-50,50,1500,-500,1000);
			}

            if(mx!=-1e4 && my !=-1e4 && x0!=-1e4 && y0!=-1e4 && z0!=-1e4)
            {
                for(int z = z0; z < 1000; z++)
                {
                    double xtarg = -(mx*(z-z0)+x0);// x= m*z+b, found from Analysis/PSI2018_SUmmer/muse_analysis_Sumer2018_4gem.C
                    double ytarg = -(my*(z-z0)+y0);//y=m*z+b
                    for(int k=0;k<numbars;k++)
                    {
                        if(plane==3 && bar==barcut[k] && RF > startmu[k] && RF < stopmu[k])
                        {
                            H2(xtarg,z,"mu/Z V X track mu","mu Z vs X track",400,-75,75,1500,-500,1000);
                            H2(ytarg,z,"mu/Z V Y track mu","mu Z vs Y track",400,-75,75,1500,-500,1000);
                            if(z==-250)
                            {
                                double yangle = -TMath::ATan(my)*1000;//*180/3.1415;
                                double xangle = -TMath::ATan(mx)*1000;//*180/3.1415;
                                if(xangle>-8*(xtarg+23))
                                {
                                    H2(xtarg,xangle,"mu/xp V xtarg","mu xp V xtarg;xtarg;xp",400,-75,75,200,-100,100);
                                    H2(ytarg,yangle,"mu/yp V ytarg","mu yp V ytarg;ytarg;yp",400,-75,75,200,-100,100);
                                    H1(xtarg,"mu/mu X distribution at z=-200","mu X distribution at z=-200;x position (mm)",400,-75,75);
                                    H2(xangle,yangle,"mu/yp vs xp at z = -200","yp vs xp at z=-200",100,-1,1,100,-1,1);
                                    H1(ytarg,"mu/mu Y distribution at z=-200","mu Y distribution at z=-200;y position (mm)",400,-75,75);
                                }
                                cd("mu");
                                H1(xangle,"Turtle/angles/x prime","mu x prime;mrad",100,-100,100);
                                H1(yangle,"Turtle/angles/y prime","mu y prime;mrad",100,-100,100);
                                if(mx!=-1e4&&x0!=-1e4&&mxifp!=-1e4&&x4!=-1e4)
                                {
                                    double xifpangle = -TMath::ATan(mxifp)*1000;//*180/3.1415;
                                    double yifpangle = -TMath::ATan(myifp)*1000;//*180/3.1415;

                                    H2(xangle,xifpangle,"Correlations/XpIFP V XpTarget","XpIFP V XpTarget;xptarget;xpIFP",200,-100,100,200,-100,100);
                                    H2(xtarg,xifp,"Correlations/X pos IFP V X pos Target","X pos IFP V X pos Target;xtarg;xifp",100,-50,50,80,-40,40);
                                    H2(xtarg,xangle,"Correlations/Xp V X pos target","Xp V X pos target;xtarg;xptarget",100,-50,50,200,-100,100);
                                    H2(xtarg,xifpangle,"Correlations/XpIFP V X pos targ","XpIFP V X pos targ;xtarg;xpifp",100,-50,50,200,-100,100);
                                    H2(xangle,xifp,"Correlations/X pos IFP V Xptarg","X pos IFP V X slope targ;slopetarg;xifp",200,-100,100,80,-40,40);
                                    H2(xifp,xifpangle,"Correlations/XpIFP V X pos IFP","XpIFP V X pos IFP;x ifp;xpifp",80,-40,40,200,-100,100);
                                    //Correlations to compare with turtle stuff


                                    H2(xifp,xangle,"Turtle/xp - IFP/xp v x IFP","xp v x IFP;xifp (mm);xp",80,-40,40,200,-100,100);
                                    H2(yifp,xangle,"Turtle/xp - IFP/xp v y IFP","xp v y IFP;yifp (mm);xp",80,-40,40,200,-100,100);
                                    H2(xifpangle,xangle,"Turtle/xp - IFP/xp v xpIFP","xp v xpIFP;xpifp;xp",200,-100,100,200,-100,100);
                                    H2(yifpangle,xangle,"Turtle/xp - IFP/xp v ypIFP","xp v ypIFP;ypifp;xp",200,-100,100,200,-100,100);
                                    H2(xifp,xtarg,"Turtle/x - IFP/x targ v x IFP","x targ v x IFP;xifp (mm); xtarg (mm)",80,-40,40,100,-50,50);
                                    H2(yifp,xtarg,"Turtle/x - IFP/x targ v y IFP","x targ v y IFP;yifp (mm); xtarg (mm)",80,-40,40,100,-50,50);
                                    H2(xifpangle,xtarg,"Turtle/x - IFP/x targ v xpIFP","x targ v xpIFP;xpifp; xtarg (mm)",200,-100,100,100,-50,50);
                                    H2(yifpangle,xtarg,"Turtle/x - IFP/x targ v ypIFP","x targ v ypIFP;ypifp; xtarg (mm)",200,-100,100,100,-50,50);

                                }
                                if(my!=-1e4&&y0!=-1e4&&myifp!=-1e4&&y4!=-1e4)
                                {
                                    double xifpangle = -TMath::ATan(mxifp)*1000;//*180/3.1415;
                                    double yifpangle = -TMath::ATan(myifp)*1000;//*180/3.1415;

                                    H2(yangle,yifpangle,"Correlations/YpIFP V YpTarget","YpIFP V YpTarget;yptarget;ypIFP",200,-100,100,200,-100,100);
                                    H2(ytarg,yifp, "Correlations/Y pos IFP V Y pos Target","Y pos IFP V Y pos Target;ytarg;yifp",100,-50,50,80,-40,40);
                                    H2(ytarg,yangle,"Correlations/YpTarget V Y pos target","Yptarget V Y pos target;ytarg;yptarget",100,-50,50,200,-100,100);
                                    H2(ytarg,yifpangle,"Correlations/YpIFP V Y pos targ","YpIFP V Y pos targ;ytarg;ypifp",100,-50,50,200,-100,100);
                                    H2(yangle,yifp,"Correlations/Y pos IFP V Yptarg","Y pos IFP V Yptarg;yptarg;yifp",200,-100,100,80,-40,40);
                                    H2(yifp,yifpangle,"Correlations/YpIFP V Y pos IFP","YpIFP V Y pos IFP;y ifp;ypifp",80,-40,40,200,-100,100);

                                    //Correlations to compare with turtle stuff

                                    H2(xifp,yangle,"Turtle/yp - IFP/yp v x IFP","yp v x IFP;xifp (mm);yp",80,-40,40,1000,-100,100);
                                    H2(yifp,yangle,"Turtle/yp - IFP/yp v y IFP","yp v y IFP;yifp (mm);yp",80,-40,40,1000,-100,100);
                                    H2(xifpangle,yangle,"Turtle/yp - IFP/yp v xpIFP","yp v xpIFP;xpifp;yp",200,-100,100,1000,-100,100);
                                    H2(yifpangle,yangle,"Turtle/yp - IFP/yp v ypIFP","yp v ypIFP;ypifp;yp",200,-100,100,1000,-100,100);
                                    H2(xifp,ytarg,"Turtle/y - IFP/y targ v x IFP","y targ v x IFP;xifp (mm); ytarg (mm)",80,-40,40,100,-50,50);
                                    H2(yifp,ytarg,"Turtle/y - IFP/y targ v y IFP","y targ v y IFP;yifp (mm); ytarg (mm)",80,-40,40,100,-50,50);
                                    H2(xifpangle,ytarg,"Turtle/y - IFP/y targ v xpIFP","y targ v xpIFP;xpIFP; ytarg (mm)",200,-100,100,100,-50,50);
                                    H2(yifpangle,ytarg,"Turtle/y - IFP/y targ v ypIFP","y targ v ypIFP;ypIFP; ytarg (mm)",200,-100,100,100,-50,50);
                                    corry_yp->Fill(xifpangle,yangle);
                                    H1(xifpangle,"Turtle/angles/xifp prime","mu xifp prime;mrad",100,-100,100);
                                    H1(yifpangle,"Turtle/angles/yifp prime","mu yifp prime;mrad",100,-100,100);
                                    H1(xifp*10,"Turtle/simple/xifp","mu xifp;x position (mm)",100,-50,50);
                                    H1(yifp*10,"Turtle/simple/yifp","mu yifp;y position (mm)",100,-50,50);

                                }
                                cd("..");
                            }

                        }
                        if(plane==3 && bar==barcut[k] && RF > starte[k] && RF < stope[k])//&& mxifp > -0.005 && mxifp < 0.00 && xifp >0 && xifp <2 )
                        {
                            H2(xtarg,z,"e/Z V X track e","e Z vs X track",400,-75,75,1500,-500,1000);
                            H2(ytarg,z,"e/Z V Y track e","e Z vs Y track",400,-75,75,1500,-500,1000);
                            if(z==-250)
                            {
                                double yangle = -TMath::ATan(my)*1000;//*180/3.1415;
                                double xangle = -TMath::ATan(mx)*1000;//*180/3.1415;

                                if(xangle>-8*(xtarg+23))
                                {
                                    H2(xtarg,xangle,"e/xp V xtarg","e xp V xtarg;xtarg;xp",400,-75,75,200,-100,100);
                                    H2(ytarg,yangle,"e/yp V ytarg","e yp V ytarg;ytarg;yp",400,-75,75,200,-100,100);
                                    H1(xtarg,"e/e X distribution at z=-200","e X distribution at z=-200;x position (mm)",400,-75,75);
                                    H2(xangle,yangle,"e/yp vs xp at z = -200","yp vs xp at z=-200",200,-30,30,200,-30,30);
                                    H1(ytarg,"e/e Y distribution at z=-200","e Y distribution at z=-200;y position (mm)",400,-75,75);
                                }
                                cd("e");
                                H1(xangle,"Turtle/angles/x prime","e x prime;mrad",100,-100,100);
                                H1(yangle,"Turtle/angles/y prime","e y p;mrad",100,-100,100);
                                if(mx!=-1e4&&x0!=-1e4&&mxifp!=-1e4&&x4!=-1e4)
                                {
                                    double xifpangle = -TMath::ATan(mxifp)*1000;//*180/3.1415;
                                    double yifpangle = -TMath::ATan(myifp)*1000;//*180/3.1415;

                                    H2(xangle,xifpangle,"Correlations/XpIFP V XpTarget","XpIFP V XpTarget;xptarget;xpIFP",200,-100,100,200,-100,100);
                                    H2(xtarg,xifp,"Correlations/X pos IFP V X pos Target","X pos IFP V X pos Target;xtarg;xifp",100,-50,50,80,-40,40);
                                    H2(xtarg,xangle,"Correlations/Xp V X pos target","Xp V X pos target;xtarg;xptarget",100,-50,50,200,-100,100);
                                    H2(xtarg,xifpangle,"Correlations/XpIFP V X pos targ","XpIFP V X pos targ;xtarg;xpifp",100,-50,50,200,-100,100);
                                    H2(xangle,xifp,"Correlations/X pos IFP V Xptarg","X pos IFP V X slope targ;slopetarg;xifp",200,-100,100,80,-40,40);
                                    H2(xifp,xifpangle,"Correlations/XpIFP V X pos IFP","XpIFP V X pos IFP;x ifp;xpifp",80,-40,40,200,-100,100);
                                    //Correlations to compare with turtle stuff


                                    H2(xifp,xangle,"Turtle/xp - IFP/xp v x IFP","xp v x IFP;xifp (mm);xp",80,-40,40,200,-100,100);
                                    H2(yifp,xangle,"Turtle/xp - IFP/xp v y IFP","xp v y IFP;yifp (mm);xp",80,-40,40,200,-100,100);
                                    H2(xifpangle,xangle,"Turtle/xp - IFP/xp v xpIFP","xp v xpIFP;xpifp;xp",200,-100,100,200,-100,100);
                                    H2(yifpangle,xangle,"Turtle/xp - IFP/xp v ypIFP","xp v ypIFP;ypifp;xp",200,-100,100,200,-100,100);
                                    H2(xifp,xtarg,"Turtle/x - IFP/x targ v x IFP","x targ v x IFP;xifp (mm); xtarg (mm)",80,-40,40,100,-50,50);
                                    H2(yifp,xtarg,"Turtle/x - IFP/x targ v y IFP","x targ v y IFP;yifp (mm); xtarg (mm)",80,-40,40,100,-50,50);
                                    H2(xifpangle,xtarg,"Turtle/x - IFP/x targ v xpIFP","x targ v xpIFP;xpifp; xtarg (mm)",200,-100,100,100,-50,50);
                                    H2(yifpangle,xtarg,"Turtle/x - IFP/x targ v ypIFP","x targ v ypIFP;ypifp; xtarg (mm)",200,-100,100,100,-50,50);

                                }
                                if(my!=-1e4&&y0!=-1e4&&myifp!=-1e4&&y4!=-1e4)
                                {
                                    double xifpangle = -TMath::ATan(mxifp)*1000;//*180/3.1415;
                                    double yifpangle = -TMath::ATan(myifp)*1000;//*180/3.1415;

                                    H2(yangle,yifpangle,"Correlations/YpIFP V YpTarget","YpIFP V YpTarget;yptarget;ypIFP",200,-100,100,200,-100,100);
                                    H2(ytarg,yifp, "Correlations/Y pos IFP V Y pos Target","Y pos IFP V Y pos Target;ytarg;yifp",100,-50,50,80,-40,40);
                                    H2(ytarg,yangle,"Correlations/YpTarget V Y pos target","Yptarget V Y pos target;ytarg;yptarget",100,-50,50,200,-100,100);
                                    H2(ytarg,yifpangle,"Correlations/YpIFP V Y pos targ","YpIFP V Y pos targ;ytarg;ypifp",100,-50,50,200,-100,100);
                                    H2(yangle,yifp,"Correlations/Y pos IFP V Yptarg","Y pos IFP V Yptarg;yptarg;yifp",200,-100,100,80,-40,40);
                                    H2(yifp,yifpangle,"Correlations/YpIFP V Y pos IFP","YpIFP V Y pos IFP;y ifp;ypifp",80,-40,40,200,-100,100);

                                    //Correlations to compare with turtle stuff

                                    H2(xifp,yangle,"Turtle/yp - IFP/yp v x IFP","yp v x IFP;xifp (mm);yp",80,-40,40,1000,-100,100);
                                    H2(yifp,yangle,"Turtle/yp - IFP/yp v y IFP","yp v y IFP;yifp (mm);yp",80,-40,40,1000,-100,100);
                                    H2(xifpangle,yangle,"Turtle/yp - IFP/yp v xpIFP","yp v xpIFP;xpifp;yp",200,-100,100,1000,-100,100);
                                    H2(yifpangle,yangle,"Turtle/yp - IFP/yp v ypIFP","yp v ypIFP;ypifp;yp",200,-100,100,1000,-100,100);
                                    H2(xifp,ytarg,"Turtle/y - IFP/y targ v x IFP","y targ v x IFP;xifp (mm); ytarg (mm)",80,-40,40,100,-50,50);
                                    H2(yifp,ytarg,"Turtle/y - IFP/y targ v y IFP","y targ v y IFP;yifp (mm); ytarg (mm)",80,-40,40,100,-50,50);
                                    H2(xifpangle,ytarg,"Turtle/y - IFP/y targ v xpIFP","y targ v xpIFP;xpIFP; ytarg (mm)",200,-100,100,100,-50,50);
                                    H2(yifpangle,ytarg,"Turtle/y - IFP/y targ v ypIFP","y targ v ypIFP;ypIFP; ytarg (mm)",200,-100,100,100,-50,50);
                                    corry_yp->Fill(xifpangle,yangle);
                                    H1(xifpangle,"Turtle/angles/xifp prime","e xifp prime;mrad",100,-100,100);
                                    H1(yifpangle,"Turtle/angles/yifp prime","e yifp prime;mrad",100,-100,100);
                                    H1(xifp*10,"Turtle/simple/xifp","e xifp;x position (mm)",100,-50,50);
                                    H1(yifp*10,"Turtle/simple/yifp","e yifp;y position (mm)",100,-50,50);

                                }
                                cd("..");
                            }

                        }
                        if( plane==3 && bar==barcut[k] && RF > startpi[k] && RF < stoppi[k])//&& mxifp > -.005 && mxifp < 0.00 && xifp >0 && xifp <2)
                        {
                            H2(xtarg,z,"pi/Z V X track pi","pi Z vs X track",400,-75,75,1500,-500,1000);
                            H2(ytarg,z,"pi/Z V Y track pi","pi Z vs Y track",400,-75,75,1500,-500,1000);
                            if(z==-250)
                            {
                                double yangle = -TMath::ATan(my)*1000;//*180/3.1415;
                                double xangle = -TMath::ATan(mx)*1000;//*180/3.1415;

                                if(xangle>-8*(xtarg+23))
                                {
                                    H2(xtarg,xangle,"pi/xp V xtarg","pi xp V xtarg;xtarg;xp",400,-75,75,200,-100,100);
                                    H2(ytarg,yangle,"pi/yp V ytarg","pi yp V ytarg;ytarg;yp",400,-75,75,200,-100,100);
                                    H1(xtarg,"pi/pi X distribution at z=-200","pi X distribution at z=-200;x position (mm)",400,-75,75);
                                    H2(xangle,yangle,"pi/yp vs xp at z = -200","yp vs xp at z=-200",200,-30,30,200,-30,30);
                                    H1(ytarg,"pi/pi Y distribution at z=-200","pi Y distribution at z=-200;y position (mm)",400,-75,75);
                                }
                                cd("pi");
                                H1(xangle,"Turtle/angles/x prime","pi x prime;mrad",100,-100,100);
                                H1(yangle,"Turtle/angles/y prime","pi y prime;mrad",100,-100,100);
                                if(mx!=-1e4&&x0!=-1e4&&mxifp!=-1e4&&x4!=-1e4)
                                {
                                    double xifpangle = -TMath::ATan(mxifp)*1000;//*180/3.1415;
                                    double yifpangle = -TMath::ATan(myifp)*1000;//*180/3.1415;

                                    H2(xangle,xifpangle,"Correlations/XpIFP V XpTarget","XpIFP V XpTarget;xptarget;xpIFP",200,-100,100,200,-100,100);
                                    H2(xtarg,xifp,"Correlations/X pos IFP V X pos Target","X pos IFP V X pos Target;xtarg;xifp",100,-50,50,80,-40,40);
                                    H2(xtarg,xangle,"Correlations/Xp V X pos target","Xp V X pos target;xtarg;xptarget",100,-50,50,200,-100,100);
                                    H2(xtarg,xifpangle,"Correlations/XpIFP V X pos targ","XpIFP V X pos targ;xtarg;xpifp",100,-50,50,200,-100,100);
                                    H2(xangle,xifp,"Correlations/X pos IFP V Xptarg","X pos IFP V X slope targ;slopetarg;xifp",200,-100,100,80,-40,40);
                                    H2(xifp,xifpangle,"Correlations/XpIFP V X pos IFP","XpIFP V X pos IFP;x ifp;xpifp",80,-40,40,200,-100,100);
                                    //Correlations to compare with turtle stuff


                                    H2(xifp,xangle,"Turtle/xp - IFP/xp v x IFP","xp v x IFP;xifp (mm);xp",80,-40,40,200,-100,100);
                                    H2(yifp,xangle,"Turtle/xp - IFP/xp v y IFP","xp v y IFP;yifp (mm);xp",80,-40,40,200,-100,100);
                                    H2(xifpangle,xangle,"Turtle/xp - IFP/xp v xpIFP","xp v xpIFP;xpifp;xp",200,-100,100,200,-100,100);
                                    H2(yifpangle,xangle,"Turtle/xp - IFP/xp v ypIFP","xp v ypIFP;ypifp;xp",200,-100,100,200,-100,100);
                                    H2(xifp,xtarg,"Turtle/x - IFP/x targ v x IFP","x targ v x IFP;xifp (mm); xtarg (mm)",80,-40,40,100,-50,50);
                                    H2(yifp,xtarg,"Turtle/x - IFP/x targ v y IFP","x targ v y IFP;yifp (mm); xtarg (mm)",80,-40,40,100,-50,50);
                                    H2(xifpangle,xtarg,"Turtle/x - IFP/x targ v xpIFP","x targ v xpIFP;xpifp; xtarg (mm)",200,-100,100,100,-50,50);
                                    H2(yifpangle,xtarg,"Turtle/x - IFP/x targ v ypIFP","x targ v ypIFP;ypifp; xtarg (mm)",200,-100,100,100,-50,50);

                                }
                                if(my!=-1e4&&y0!=-1e4&&myifp!=-1e4&&y4!=-1e4)
                                {
                                    double xifpangle = -TMath::ATan(mxifp)*1000;//*180/3.1415;
                                    double yifpangle = -TMath::ATan(myifp)*1000;//*180/3.1415;

                                    H2(yangle,yifpangle,"Correlations/YpIFP V YpTarget","YpIFP V YpTarget;yptarget;ypIFP",200,-100,100,200,-100,100);
                                    H2(ytarg,yifp, "Correlations/Y pos IFP V Y pos Target","Y pos IFP V Y pos Target;ytarg;yifp",100,-50,50,80,-40,40);
                                    H2(ytarg,yangle,"Correlations/YpTarget V Y pos target","Yptarget V Y pos target;ytarg;yptarget",100,-50,50,200,-100,100);
                                    H2(ytarg,yifpangle,"Correlations/YpIFP V Y pos targ","YpIFP V Y pos targ;ytarg;ypifp",100,-50,50,200,-100,100);
                                    H2(yangle,yifp,"Correlations/Y pos IFP V Yptarg","Y pos IFP V Yptarg;yptarg;yifp",200,-100,100,80,-40,40);
                                    H2(yifp,yifpangle,"Correlations/YpIFP V Y pos IFP","YpIFP V Y pos IFP;y ifp;ypifp",80,-40,40,200,-100,100);

                                    //Correlations to compare with turtle stuff

                                    H2(xifp,yangle,"Turtle/yp - IFP/yp v x IFP","yp v x IFP;xifp (mm);yp",80,-40,40,1000,-100,100);
                                    H2(yifp,yangle,"Turtle/yp - IFP/yp v y IFP","yp v y IFP;yifp (mm);yp",80,-40,40,1000,-100,100);
                                    H2(xifpangle,yangle,"Turtle/yp - IFP/yp v xpIFP","yp v xpIFP;xpifp;yp",200,-100,100,1000,-100,100);
                                    H2(yifpangle,yangle,"Turtle/yp - IFP/yp v ypIFP","yp v ypIFP;ypifp;yp",200,-100,100,1000,-100,100);
                                    H2(xifp,ytarg,"Turtle/y - IFP/y targ v x IFP","y targ v x IFP;xifp (mm); ytarg (mm)",80,-40,40,100,-50,50);
                                    H2(yifp,ytarg,"Turtle/y - IFP/y targ v y IFP","y targ v y IFP;yifp (mm); ytarg (mm)",80,-40,40,100,-50,50);
                                    H2(xifpangle,ytarg,"Turtle/y - IFP/y targ v xpIFP","y targ v xpIFP;xpIFP; ytarg (mm)",200,-100,100,100,-50,50);
                                    H2(yifpangle,ytarg,"Turtle/y - IFP/y targ v ypIFP","y targ v ypIFP;ypIFP; ytarg (mm)",200,-100,100,100,-50,50);
                                    corry_yp->Fill(xifpangle,yangle);
                                    H1(xifpangle,"Turtle/angles/xifp prime","pi xifp prime;mrad",100,-100,100);
                                    H1(yifpangle,"Turtle/angles/yifp prime","pi yifp prime;mrad",100,-100,100);
                                    H1(xifp*10,"Turtle/simple/xifp","pi xifp;x position (mm)",100,-50,50);
                                    H1(yifp*10,"Turtle/simple/yifp","pi yifp;y position (mm)",100,-50,50);

                                }
                                cd("..");
                            }
                        }
                    }
                }
            }

            //PID plots
            // for(int k=0;k<numbars;k++)
            // {
            //     if(x0!=-1e4 && y0!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > startmu && RF+alignment[k] < stopmu )
            //         H2(x0,y0,"mu/US GEM x V y cut mu","US GEM x V y",500,-50,50,1000,-50,50);
            //     if(x1!=-1e4 && y1!=-1e4&& plane==3 && bar==barcut[k] && RF+alignment[k] > startmu && RF+alignment[k] < stopmu )
            //         H2(x1,y1,"mu/4th GEM x V y cut mu","4th GEM x V y",500,-50,50,1000,-50,50);
            //     if(x2!=-1e4 && y2!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > startmu && RF+alignment[k] < stopmu )
            //         H2(x2,y2,"mu/MI GEM 2 x V y cut mu","MI GEM x V y",500,-50,50,500,-50,50);
            //     if(x3!=-1e4 && y3!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > startmu && RF+alignment[k] < stopmu )
            //         H2(x3,y3,"mu/DS GEM 3 x V y cut mu","DS GEM x V y",500,-50,50,500,-50,50);

            //     if(x4!=-1e4 && y4!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > startmu && RF+alignment[k] < stopmu )
            //         H2(x4,y4,"mu/US - IFP x V y cut mu","US - IFP x V y",500,-50,50,500,-50,50);
            //     if(x5!=-1e4 && y5!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > startmu && RF+alignment[k] < stopmu )
            //         H2(x5,y5,"mu/DS - IFP x V y cut mu","DS - IFP x V y",500,-50,50,500,-50,50);
            //     if(mxifp!=-1e4 && myifp!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > startmu && RF+alignment[k] < stopmu )
            //         H2(mxifp,myifp,"mu/IFP my V mx cut mu","IFP my V mx",100,-1,1,100,-1,1);


            //     if(x0!=-1e4 && y0!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > starte && RF+alignment[k] < stope )
            //         H2(x0,y0,"e/US GEM x V y cut e","US GEM x V y",500,-50,50,1000,-50,50);
            //     if(x1!=-1e4 && y1!=-1e4&& plane==3 && bar==barcut[k] && RF+alignment[k] > starte && RF+alignment[k] < stope )
            //         H2(x1,y1,"e/4th GEM x V y cut e","4th GEM x V y",500,-50,50,1000,-50,50);
            //     if(x2!=-1e4 && y2!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > starte && RF+alignment[k] < stope )
            //         H2(x2,y2,"e/MI GEM x V y cut e","MI GEM x V y",500,-50,50,500,-50,50);
            //     if(x3!=-1e4 && y3!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > starte && RF+alignment[k] < stope )
            //         H2(x3,y3,"e/DS GEM x V y cut e","DS GEM x V y",500,-50,50,500,-50,50);

            //     if(x4!=-1e4 && y4!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > starte && RF+alignment[k] < stope )
            //         H2(x4,y4,"e/US - IFP x V y cut e","US - IFP x V y",500,-50,50,500,-50,50);
            //     if(x5!=-1e4 && y5!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > starte && RF+alignment[k] < stope )
            //         H2(x5,y5,"e/DS - IFP x V y cut e","DS - IFP x V y",500,-50,50,500,-50,50);
            //     if(mxifp!=-1e4 && myifp!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > starte && RF+alignment[k] < stope )
            //         H2(mxifp,myifp,"e/IFP my V mx cut e","IFP my V mx",100,-1,1,100,-1,1);


            //     if(x0!=-1e4 && y0!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > startpi && RF+alignment[k] < stoppi )
            //         H2(x0,y0,"pi/US GEM x V y cut pi","US GEM x V y",500,-50,50,1000,-50,50);
            //     if(x1!=-1e4 && y1!=-1e4&& plane==3 && bar==barcut[k] && RF+alignment[k] > startpi && RF+alignment[k] < stoppi )
            //         H2(x1,y1,"pi/4th GEM 1 x V y cut pi","4th GEM x V y",500,-50,50,1000,-50,50);
            //     if(x2!=-1e4 && y2!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > startpi && RF+alignment[k] < stoppi )
            //         H2(x2,y2,"pi/MI GEM 2 x V y cut pi","MI GEM x V y",500,-50,50,500,-50,50);
            //     if(x3!=-1e4 && y3!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > startpi && RF+alignment[k] < stoppi )
            //         H2(x3,y3,"pi/DS GEM x V y cut pi","DS GEM x V y",500,-50,50,500,-50,50);

            //     if(x4!=-1e4 && y4!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > startpi && RF+alignment[k] < stoppi )
            //         H2(x4,y4,"pi/US - IFP x V y cut pi","US - IFP x V y",500,-50,50,500,-50,50);
            //     if(x5!=-1e4 && y5!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > startpi && RF+alignment[k] < stoppi )
            //         H2(x5,y5,"pi/DS - IFP x V y cut pi","DS - IFP x V y",500,-50,50,500,-50,50);
            //     if(mxifp!=-1e4 && myifp!=-1e4 && plane==3 && bar==barcut[k] && RF+alignment[k] > startpi && RF+alignment[k] < stoppi )
            //         H2(mxifp,myifp,"pi/IFP my V mx cut pi","IFP my V mx",100,-1,1,100,-1,1);

            // }

					//printf("HERE\n");
			for(int i=0; i<4; i++){
				for(int j=0; j<16; j++){
        	if(x0!=-1e4 && y0!=-1e4 && ((plane==i && bar==j && timeleft <-250 && timeleft>-340)))// || (plane==3 && bar==10)) )
        		H2(x0,y0,TString::Format("GEM-BH/Plane %i/US GEM/x V y cut on BH bar %i", i, j),"US GEM x V y",500,-50,50,1000,-50,50);
        	if(x1!=-1e4 && y1!=-1e4&& ((plane==i && bar==j && timeleft <-250 && timeleft>-340)))// || (plane==3 && bar==10)))
        		H2(x1,y1,TString::Format("GEM-BH/Plane %i/4th GEM/x V y cut on BH bar %i", i, j),"4th GEM x V y",500,-50,50,1000,-50,50);
        	if(x2!=-1e4 && y2!=-1e4 &&((plane==i && bar==j && timeleft <-250 && timeleft>-340)))// || (plane==3 && bar==10)))
        		H2(x2,y2,TString::Format("GEM-BH/Plane %i/MI GEM/x V y cut on BH bar %i", i, j),"MI GEM x V y",500,-50,50,500,-50,50);
        	if(x3!=-1e4 && y3!=-1e4 && ((plane==i && bar==j && timeleft <-250 && timeleft>-340)))// || (plane==3 && bar==10)))
        		H2(x3,y3,TString::Format("GEM-BH/Plane %i/DS GEM/x V y cut on BH bar %i", i, j),"DS GEM x V y",500,-50,50,500,-50,50);
				}
			}

            if(plane==3 && bar==4)
                H1(timeleft,"bar 4 left time","bar 4 left time",1000,-1,-1);
            if(timeleft <-450 && timeleft>-550 && x1!=-1e4)
                std::cout <<"GEM X: " << x1 << " Y: " << y1 << " Bar Number: " << bar << std::endl;

		}


        if(x0!=-1e4 && mx!=-1e4 && y0!=-1e4 && my!=-1e4) {
            debug(1,"x0 = %.2f , mx = %.2f , y0 = %.2f , my = %.2f\n",x0,mx,y0,my);
        }
    }
	return ok;
}

Long_t GEM_BH::finalize()
{
 // std::ofstream outfile;
 // outfile.open("Plots/GEM/BeamProfile.txt", std::ios_base::app);
 // //theRunInfo->runNumber;
 // outfile << "Z pos \t X mean (mm) \t X sig \t Y mean (mm) \t Y sig \t zx slope mean \t zx slope sig \t zy slope mean \t zy slope sig \n";
 //    for(int i=0; i<26; i++)
 //    {
 //        double params[12] ={0.};

 //        TF1 *g1 = new TF1("g1","gaus",-10,10);
 //        TF1 *g2 = new TF1("g2","gaus",-10,10);
 //        TF1 *g3 = new TF1("g3","gaus",-40,40);
 //        TF1 *g4 = new TF1("g4","gaus",-40,40);

 //        xpos[i]->Fit("g1");
 //        ypos[i]->Fit("g2");
 //        xzslope[i]->Fit("g3");
 //        yzslope[i]->Fit("g4");

 //        g1->GetParameters(&params[0]);
 //        g2->GetParameters(&params[3]);
 //        g3->GetParameters(&params[6]);
 //        g4->GetParameters(&params[9]);

 //        outfile << -500+i*25 << "\t" << params[1] << "\t" << params[2] << "\t" << params[4] << "\t" << params[5] << "\t" << params[7] << "\t" << params[8] << "\t" << params[10] << "\t" << params[11] << "\n";

 //    }
 // outfile.close();
   //  auto muxplot = dH1("mu/mu X distribution at z=-200","mu X distribution at z=-200;x position (mm)",400,-75,75);
   //  TCanvas *c1 = new TCanvas();
   //  gStyle->SetOptStat(0);
   //  gStyle->SetOptFit();
   //  c1->Divide(2,2);
   //  c1->cd(1);
   //  muxplot->Draw();
   //  muxplot->Fit("gaus");
   //  c1->cd(3);
   //  muxplot = dH1("mu/mu Y distribution at z=-200","mu Y distribution at z=-200;y position (mm)",400,-75,75);
   //  muxplot->Draw();
   //  muxplot->Fit("gaus");
   //  c1->cd(2);
   //  muxplot = dH1("mu/Turtle/angles/x prime","mu x prime;mrad",100,-100,100);
   //  muxplot->Draw();
   //  muxplot->Fit("gaus");
   //  c1->cd(4);
   //  muxplot = dH1("mu/Turtle/angles/y prime","mu y prime;mrad",100,-100,100);
   //  muxplot->Draw();
   //  muxplot->Fit("gaus");
   //  std::string muxplottitle = "Plots/Mu_targ_distribution_"+std::to_string(theRunInfo->runNumber)+".pdf";
   //  c1->SaveAs(muxplottitle.c_str());


   //  auto pixplot = dH1("pi/pi X distribution at z=-200","pi X distribution at z=-200;x position (mm)",400,-75,75);
   //  TCanvas *c2 = new TCanvas();
   //  gStyle->SetOptStat(0);
   //  gStyle->SetOptFit();
   //  c2->Divide(2,2);
   //  c2->cd(1);
   //  pixplot->Draw();
   //  pixplot->Fit("gaus");
   //  c2->cd(3);
   //  pixplot = dH1("pi/pi Y distribution at z=-200","pi Y distribution at z=-200;y position (mm)",400,-75,75);
   //  pixplot->Draw();
   //  pixplot->Fit("gaus");
   //  c2->cd(2);
   //  pixplot = dH1("pi/Turtle/angles/x prime","pi x prime;mrad",100,-100,100);
   //  pixplot->Draw();
   //  pixplot->Fit("gaus");
   //  c2->cd(4);
   //  pixplot = dH1("pi/Turtle/angles/y prime","pi y prime;mrad",100,-100,100);
   //  pixplot->Draw();
   //  pixplot->Fit("gaus");
   //  std::string pixplottitle = "Plots/Pi_targ_distribution_"+std::to_string(theRunInfo->runNumber)+".pdf";
   //  c2->SaveAs(pixplottitle.c_str());

   // auto explot = dH1("e/e X distribution at z=-200","e X distribution at z=-200;x position (mm)",400,-75,75);
   //  TCanvas *c3 = new TCanvas();
   //  gStyle->SetOptStat(0);
   //  gStyle->SetOptFit();
   //  c3->Divide(2,2);
   //  c3->cd(1);
   //  explot->Draw();
   //  explot->Fit("gaus");
   //  c3->cd(3);
   //  explot = dH1("e/e Y distribution at z=-200","e Y distribution at z=-200;y position (mm)",400,-75,75);
   //  explot->Draw();
   //  explot->Fit("gaus");
   //  c3->cd(2);
   //  explot = dH1("e/Turtle/angles/x prime","e x prime;mrad",100,-100,100);
   //  explot->Draw();
   //  explot->Fit("gaus");
   //  c3->cd(4);
   //  explot = dH1("e/Turtle/angles/y prime","e y prime;mrad",100,-100,100);
   //  explot->Draw();
   //  explot->Fit("gaus");
   //  std::string explottitle = "Plots/e_targ_distribution_"+std::to_string(theRunInfo->runNumber)+".pdf";
   //  c3->SaveAs(explottitle.c_str());

   // auto muxplot = dH1("mu/mu Y distribution at z=-200","mu Y distribution at z=-200;y position (mm)",400,-75,75);
   // auto explot = dH1("e/e Y distribution at z=-200","e Y distribution at z=-200;y position (mm)",400,-75,75);
   // auto pixplot = dH1("pi/pi Y distribution at z=-200","pi Y distribution at z=-200;y position (mm)",400,-75,75);

   // TCanvas *c = new TCanvas();
   // gStyle->SetOptStat(0);
   // muxplot->SetLineColor(kRed+1);
   // explot->SetLineColor(kBlue+1);
   // pixplot->SetLineColor(kGreen+1);
   // explot->SetMarkerStyle(20);
   //  muxplot->SetMarkerStyle(21);
   //  pixplot->SetMarkerStyle(22);
   //  muxplot->SetMarkerColor(kRed+1);
   // explot->SetMarkerColor(kBlue+1);
   // pixplot->SetMarkerColor(kGreen+1);
   //   muxplot->SetLineWidth(2);
   // explot->SetLineWidth(2);
   // pixplot->SetLineWidth(2);
   // pixplot->SetTitle("e/mu/pi distribution comparison");
   // pixplot->SetName("pi");
   // explot->SetName("e");
   // muxplot->SetName("mu");
   // double scalemu = pixplot->Integral()/muxplot->Integral();
   // double scalee = pixplot->Integral()/explot->Integral();
   // muxplot->Rebin(5);
   // explot->Rebin(5);
   // pixplot->Rebin(5);
   // muxplot->Scale(scalemu);
   // explot->Scale(scalee);

   // pixplot->Draw();
   // muxplot->Draw("SAME");
   // explot->Draw("SAME");

   //  auto legend2 = new TLegend(0.1,0.7,0.48,0.9);
   // legend2->SetHeader("Particle Ratios","C"); // option "C" allows to center the header
   // legend2->AddEntry("pi","Pion Y Distribution","ep");
   // legend2->AddEntry("e","Electron Y Distribution","ep");
   //  legend2->AddEntry("mu","Muon Y Distribution","ep");
   // legend2->Draw();

   // c->Update();
   // c->SaveAs("overlap.pdf");


    // auto epplot = dH1("e/Turtle/simple/xifp","e xifp;x position (mm)",80,-40,40);
    // TCanvas *c9 = new TCanvas();
    // c9->Divide(2,2);
    // c9->cd(1);
    // epplot->Draw();
    // c9->cd(3);
    // epplot = dH1("e/Turtle/simple/yifp","e yifp;y position (mm)",80,-40,40);
    // epplot->Draw();
    // c9->cd(2);
    // epplot = dH1("e/Turtle/angles/xifp prime","e xifp prime;mrad",100,-100,100);
    // epplot->Draw();
    // c9->cd(4);
    // epplot = dH1("e/Turtle/angles/yifp prime","e yifp prime;mrad",100,-100,100);
    // epplot->Draw();
    // std::string epplottitle = "Plots/e_ifp_distribution_"+std::to_string(theRunInfo->runNumber)+".pdf";
    // c9->SaveAs(epplottitle.c_str());

    // auto mupplot = dH1("mu/Turtle/simple/xifp","mu xifp;x position (mm)",80,-40,40);
    // TCanvas *c7 = new TCanvas();
    // c7->Divide(2,2);
    // c7->cd(1);
    // mupplot->Draw();
    // c7->cd(3);
    // mupplot = dH1("mu/Turtle/simple/yifp","mu yifp;y position (mm)",80,-40,40);
    // mupplot->Draw();
    // c7->cd(2);
    // mupplot = dH1("mu/Turtle/angles/xifp prime","mu xifp prime;mrad",100,-100,100);
    // mupplot->Draw();
    // c7->cd(4);
    // mupplot = dH1("mu/Turtle/angles/yifp prime","mu yifp prime;mrad",100,-100,100);
    // mupplot->Draw();
    // std::string mupplottitle = "Plots/mu_ifp_distribution_"+std::to_string(theRunInfo->runNumber)+".pdf";
    // c7->SaveAs(mupplottitle.c_str());

    // auto pipplot = dH1("pi/Turtle/simple/xifp","pi xifp;x position (mm)",80,-40,40);
    // TCanvas *c8 = new TCanvas();
    // c8->Divide(2,2);
    // c8->cd(1);
    // pipplot->Draw();
    // c8->cd(3);
    // pipplot = dH1("pi/Turtle/simple/yifp","pi yifp;y position (mm)",80,-40,40);
    // pipplot->Draw();
    // c8->cd(2);
    // pipplot = dH1("pi/Turtle/angles/xifp prime","pi xifp prime;mrad",100,-100,100);
    // pipplot->Draw();
    // c8->cd(4);
    // pipplot = dH1("pi/Turtle/angles/yifp prime","pi yifp prime;mrad",100,-100,100);
    // pipplot->Draw();
    // std::string pipplottitle = "Plots/pi_ifp_distribution_"+std::to_string(theRunInfo->runNumber)+".pdf";
    // c8->SaveAs(pipplottitle.c_str());

	return ok;
}

Long_t GEM_BH::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new GEM_BH(in,out,inf_,outf_,p);
}
}


ClassImp(GEM_BH);
