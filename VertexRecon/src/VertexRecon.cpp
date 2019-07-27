#include <VertexRecon.h>

#include<iostream>
#include<cmath>


VertexRecon::VertexRecon(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

VertexRecon::~VertexRecon()
{
};


Long_t VertexRecon::startup()
{
	GEM_Tracks=NULL;
	STT_Tracks=NULL;
        Hits=NULL;

    getBranchObject("teletracks"    ,(TObject **) & GEM_Tracks);
	getBranchObject("TrackHits", (TObject **) &STT_Tracks);
        getBranchObject("ScintHits", (TObject **) &Hits);
    theRunInfo=(MRTRunInfo *) getFileObject("RunInfo");
	if(!STT_Tracks || !GEM_Tracks)
	{
		debug(0,"Couldn't find STT_Tracks or GEM_Tracks in tree");
		return 0;
	}
	return ok;
}



/*
Email from Steffen Strauch:
Below is a function I used in the analysis of simulated data.  It takes two tracks as inputs.  
Each track is characterized by one point on the path and a direction vector; 
P0 and Q0 for the points, u and v for the directions. 
It returns the vertex as the point of closest approach, the distance of closest approach, and the angle between the two tracks (directional vectors).  
Maybe you like to try this out.

*/
void VertexRecon::find_vertex(TVector3 P0, TVector3 u, TVector3 Q0, TVector3 v, 
		    TVector3 &vertex, double &doca, double &theta) {
	TVector3 w0 = P0 - Q0;
	Double_t a = u * u;
	Double_t b = u * v;
	Double_t c = v * v;
	Double_t d = u * w0;
	Double_t e = v * w0;

	Double_t s = (b*e - c*d) / (a*c - b*b);
	Double_t t = (a*e - b*d) / (a*c - b*b);

  
	TVector3 P = P0 + s * u;
	TVector3 Q = Q0 + t * v;

	vertex = (P + Q) * 0.5;
	doca = (P - Q).Mag();
	theta = TMath::ACos( u.Dot(v) / u.Mag() / v.Mag() );
}

void VertexRecon::id_to_info(int id, int *plane, int *bar, int *side)
{
    *bar = id%16;
    *plane = (id/16)%4;
    *side = (id/(16*4));
    return;
}

//Event by event vertex reconstruction
Long_t VertexRecon::process()
{
        const int numbars=16;

        double barcut[numbars] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
    double startmu[numbars] = {13,14,12,16,14,15,15,14,16,16,17,17,12,14,14,14};
    double stopmu[numbars] = {17,19,16,19,16,18,18,18,19,19,19,19,16,17,17,17};

    double startpi[numbars] = {1,2,0,2,0,1,2,1,2,3,3,3,19,0,0,0};
    double stoppi[numbars] = {4,6,3,6,3,5,5,4,5,7,6,7,2,4,4,5};

    double starte[numbars] = {5,7,4,8,4,6,6,5,7,7,8,8,3,5,5,5};
    double stope[numbars] = {8,10,8,10,7,9,10,8,10,11,11,11,7,8,8,9};
	
	for (size_t m = 0 ; m < GEM_Tracks->tracks.size() ; m++)
 	{

        if(GEM_Tracks -> tracks.size()>30) continue;
        double gemx0 = GEM_Tracks -> tracks[m].x0;//double check your code to verify this is in mm
        double gemmx = TMath::ATan(GEM_Tracks -> tracks[m].mx);
        double gemy0 = GEM_Tracks -> tracks[m].y0;
        double gemmy = TMath::ATan(GEM_Tracks -> tracks[m].my);
        double gemz0 = GEM_Tracks -> tracks[m].z0;
        if(gemx0==-1e4)//no valid track
            continue;
        double z = -1000;
        double xgemtarget = gemmx*(z-gemz0)+gemx0;// x= m*z+b, found from Analysis/PSI2018_SUmmer/muse_analysis_Sumer2018_4gem.C
        double ygemtarget = gemmy*(z-gemz0)+gemy0;//y=m*z+b

        // for(int i =-10000; i<10000;i++)
        // {
        //     xgemtarget = gemmx*(i-gemz0)+gemx0;// x= m*z+b, found from Analysis/PSI2018_SUmmer/muse_analysis_Sumer2018_4gem.C
        //     ygemtarget = gemmy*(i-gemz0)+gemy0;  
        //     H2(xgemtarget,i,"Gem track x","GEM track x",100,-50,50,1000,-10000,10000);
        //     H2(ygemtarget,i,"Gem track y","GEM track y",100,-50,50,1000,-10000,10000);

        // }

        TVector3 P0(gemx0,gemy0,gemz0);//Position of Incoming track
        TVector3 u(gemmx,gemmy,1);//Direction of incoming track
        // std::cout <<"GEM Tracks: " << GEM_Tracks->tracks.size() << std::endl;
        for(size_t i = 0; i < STT_Tracks->tracks.size(); i++)
        {
            if(STT_Tracks->tracks.size()>1)
            std::cout <<"STT Tracks: " << STT_Tracks->tracks.size() << std::endl;
        	double sttx0 = STT_Tracks->tracks[i].x0;//double check your code to verify this is in mm
        	double stty0 = STT_Tracks->tracks[i].y0;
        	double sttz0 = STT_Tracks->tracks[i].z0;
        	double sttmx = M_PI/2-TMath::ATan(STT_Tracks->tracks[i].mx);
        	double sttmy = M_PI/2-TMath::ATan(STT_Tracks->tracks[i].my*TMath::Cos(60*M_PI/180));

            if(sttx0==-100000||stty0==-100000||sttmx==-100000||sttmy==-100000||sttx0>100)//no valid track here
                continue;
           //  z=0;
           //  // for(int j = -10000; j <10000; j++)
           //  // {
           //      double xstttarget = sttmx*(z-sttz0)+sttx0;
           //      double ystttarget = sttmy*(z-sttz0)+stty0;
           //      double xdiff = xstttarget-xgemtarget;
           //      double ydiff = ystttarget-ygemtarget;
           //      double xanglediff = sttmx-gemmx;
           //      double yanglediff = sttmy-gemmy;
           //      H1(sttmx,"STT mx","STT mx",200,-3,3);
           //      H1(sttmy,"STT my","STT my",200,-3,3);
           //      H1(gemmx,"GEM mx","GEM mx",200,-3,3);
           //      H1(gemmy,"GEM my","GEM my",200,-3,3);


           //      H1(xanglediff,"X angle differences between GEM and STT","X angle Differences between GEM and STT;radians",200,-3,3);
           //      H1(yanglediff,"Y angle differences between GEM and STT","Y angle Differences between GEM and STT;radians",200,-3,3);

           //      //Only valid when STT is perpendicular to the beam
           //      std::cout << sttmx << " " << sttmy << " " << gemmy << " " << gemmx << std::endl;
           //      //H2(xstttarget,j,"STT XZ","STT XZ",100,-500,500,100,-10000,10000);
           //      //H2(ystttarget,j,"STT YZ","STT YZ",100,-500,500,100,-10000,10000);

           //      H2(xgemtarget,ygemtarget,"GEM at targ","GEM at targ",100,-50,50,100,-50,50);
           //      H2(xstttarget,ystttarget,"STT at targ","STT at targ",100,-100,100,100,-100,100);
           //      H1(xdiff,"Xdiff at targ","xdiff at targ",100,-100,100);

           //      H2(xdiff,ydiff,"Differences in track position at Target","Differences in track position at target",100,-1,-1,100,-1,-1);
           // // }
            for(size_t l = 0; l < Hits->hits.size(); l++)
            {
            
                int plane, side, bar;
                id_to_info(Hits->hits[l].id,&plane,&bar,&side);
                            double RF = Hits->hits[l].rf;

            	TVector3 Q0(sttx0,stty0,sttz0);
            	TVector3 v(sttmx,sttmy,1);

            	TVector3 vertex(-1000,-1000,-1000);
            	double doca=-1000;
            	double theta=-1000;
                for(int k = 0; k<numbars;k++)
                {
                    if(plane==3 && bar==barcut[k] && RF > starte[k] && RF < stope[k])//&& mxifp > -0.005 && mxifp < 0.00 && xifp >0 && xifp <2 )
                    {
                        cd("e");
                    	find_vertex(P0,u,Q0,v,vertex,doca,theta);
                    	H1(theta*180/M_PI,"Reconstructed Electron Scattering Angle","Reconstructed Electron Scattering Angle; degrees",180,-0.5,179.5);
                    	H1(doca,"DOCA","DOCA; mm",100,0,1000);
                    	H1(vertex.X(),"Reconstruced Electron X vertex","Reconstructed Electron X Vertex; mm",100,-1000,1000);
                    	H1(vertex.Y(),"Reconstruced Electron Y vertex","Reconstructed Electron Y Vertex; mm",100,-1000,1000);
                    	H1(vertex.Z(),"Reconstruced Electron Z vertex","Reconstructed Electron Z Vertex; mm",100,-1000,1000);
                        cd("..");
                    }
                    if(plane==3 && bar==barcut[k] && RF > startmu[k] && RF < stopmu[k])//&& mxifp > -0.005 && mxifp < 0.00 && xifp >0 && xifp <2 )
                    {
                        cd("mu");
                        find_vertex(P0,u,Q0,v,vertex,doca,theta);
                        H1(theta*180/M_PI,"Reconstructed Scattering Angle","Reconstructed Scattering Angle; degrees",180,-0.5,179.5);
                        H1(doca,"DOCA","DOCA; mm",100,0,1000);
                        H1(vertex.X(),"Reconstruced X vertex","Reconstructed X Vertex; mm",100,-1000,1000);
                        H1(vertex.Y(),"Reconstruced Y vertex","Reconstructed Y Vertex; mm",100,-1000,1000);
                        H1(vertex.Z(),"Reconstruced Z vertex","Reconstructed Z Vertex; mm",100,-1000,1000);
                        cd("..");
                    }
                    if(plane==3 && bar==barcut[k] && RF > startpi[k] && RF < stoppi[k])//&& mxifp > -0.005 && mxifp < 0.00 && xifp >0 && xifp <2 )
                    {
                        cd("pi");
                        find_vertex(P0,u,Q0,v,vertex,doca,theta);
                        H1(theta*180/M_PI,"Reconstructed Scattering Angle","Reconstructed Scattering Angle; degrees",180,-0.5,179.5);
                        H1(doca,"DOCA","DOCA; mm",100,0,1000);
                        H1(vertex.X(),"Reconstruced X vertex","Reconstructed X Vertex; mm",100,-1000,1000);
                        H1(vertex.Y(),"Reconstruced Y vertex","Reconstructed Y Vertex; mm",100,-1000,1000);
                        H1(vertex.Z(),"Reconstruced Z vertex","Reconstructed Z Vertex; mm",100,-1000,1000);
                        cd("..");
                    }
                }
            }//end of BH

        }//end of STTs
        break;
    }//end of GEMs
	return ok;
}



Long_t VertexRecon::finalize()
{

	return ok;
}

Long_t VertexRecon::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new VertexRecon(in,out,inf_,outf_,p);
}
}


ClassImp(VertexRecon);

