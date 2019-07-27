#include <GEM_STT.h>

#include<iostream>
#include<cmath>


GEM_STT::GEM_STT(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

GEM_STT::~GEM_STT()
{
};

TVector2 getStrawPos(int id)
{
  // ask gdml or something, for now, make it up.
  
  TVector2 res;
  int side, plane, straw_in_plane;
  STT_internal_to_logic(id,&side,&plane,&straw_in_plane);
  double tube_diameter = 10;
  double distance_between_planes = 8.7;
  double pitch=10.1; //10.1 mm spacing
  double zoffset=0;
  
  //printf("side %i plane %i straw %i\n",side,plane,straw_in_plane);
  if (side==2)
  {
    zoffset = -493;//this far off from Tom's drawing
    double pos=(straw_in_plane)*pitch+(plane % 2)*pitch/2+20;
    res.Set(pos-17*pitch,zoffset+(10-plane)*distance_between_planes);
    res = res.Rotate(0*(M_PI/180));//-60 at scattering position
  }
  if (side==0)
  {
    zoffset = 10.*distance_between_planes+55-493;//
    double pos=straw_in_plane*pitch+(plane % 2)*pitch/2;
    res.Set(pos,zoffset+(10-plane)*distance_between_planes);
    res = res.Rotate(0*(M_PI/180));//-60 at scattering position
  }

  if(side == 1)
  {
      zoffset = 0.;
      double pos = straw_in_plane*pitch +(plane % 2)*pitch/2;//original plus
      res.Set(pos,zoffset+plane*distance_between_planes);
  }
  if (side==3)
  {
    zoffset = 10.;//
    double pos=straw_in_plane*pitch+(plane % 2)*pitch/2;
    res.Set(pos,zoffset+plane*distance_between_planes);
  }
  return res;
};

Long_t GEM_STT::startup()
{

  STT = NULL;
  getOutBranchObject("StrawTubeHits",(TObject **) &STT);
  if(!STT) getBranchObject("StrawTubeHits",(TObject **) &STT); 
  if(!STT) debug(0,"Could not find STT hits in file\n");

  	GEM_Tracks=NULL;
    getBranchObject("teletracks"    ,(TObject **) & GEM_Tracks);

      theRunInfo=(MRTRunInfo *) getFileObject("RunInfo");

    bothXZ=dH2("Tracks XZ","Tracks XZ;Pos X (mm);Pos Z (mm)",501,-1000,1000,4001,-1000,1000);
    bothYZ=dH2("Tracks YZ","Tracks YZ;Pos Y (mm);Pos Z (mm)",501,-1000,1000,4001,-1000,1000);
  circXZ=dH2("circ XZ","circ XZ",501,-1000,1000,4001,-2000,1000);
  circYZ=dH2("circ YZ","circ YZ",501,-1000,1000,4001,-2000,1000);
  circ=dH2("circ","circ",501,-1000,1000,4001,-600,600);



	return ok;
}

Long_t GEM_STT::process()
{
	  const TVector2 offx(365,0);//For GEMs in middle dowel position and STT perpendicular to beam, Tom sees 179 mm separation in drawing
	  const TVector2 offy(285,0);
  int numhitsX = 0;
  int numhitsY = 0;
	for (auto hit:STT->hits)
	{
      //ret=Plugin::stop;
      int side, plane, straw_in_plane;
      STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
      //X coordinate
      if ((side==2||side==0)&& plane <5)
	   {
        if(side==2 && (plane==2||plane==3||plane==4))
          continue;
        if(side==0 && (plane==1||plane==2||plane==3))
          continue;
        if(side==0 && plane == 0 && straw_in_plane==33)
          continue;
      	numhitsX++;
  		}
      if ((side==2||side==0)&& plane >4)
	  {
	    if(side==2 && plane==8 && (straw_in_plane==11||straw_in_plane<6))
	      continue;
	  	numhitsY++;
	  }
	}

	if(numhitsY!=3||numhitsX!=3)
		return ok;

	for (size_t j = 0 ; j < GEM_Tracks->tracks.size() ; j++)
    {
    	bothYZ->Reset();
    	bothXZ->Reset();
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
        double resprimex = 0;
        double resprimey = 0;

        if(x0 == -10000)
        	break;

		for (auto hit:STT->hits)
		{
		      //ret=Plugin::stop;
		      int side, plane, straw_in_plane;
		      STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
		      //X coordinate
		      if ((side==2||side==0)&& plane <5)
			   {
		        if(side==2 && (plane==2||plane==3||plane==4))
		          continue;
		        if(side==0 && (plane==1||plane==2||plane==3))
		          continue;
		        if(side==0 && plane == 0 && straw_in_plane==33)
		          continue;
			      auto pos=getStrawPos(hit.id)-offx;

			      // for (int i =0; i<180;i++)// probably should do -90 to 90 but that means things have to be rewritten more intelligently for the negative angles.
			      // {
			      //  TVector2 d=pos+TVector2(hit.dist,0).Rotate(2*i*M_PI/180);
		       //   	bothXZ->Fill(d.X(),d.Y());
		       //   	circXZ->Fill(d.X(),d.Y());
			      //   circ->Fill(d.X(),d.Y());		         	
			      // } 
         //    	  double linetempX = (pos.X()/mx+pos.Y()- x0)/(mx+1/mx);//this is given by solving equation for two perpendicular intersecting lines 
         //    	  double linetempY = mx*linetempX+x0;
         //    	  resprimex = sqrt(pow(pos.X()-linetempX,2)+pow(pos.Y()-linetempY,2))-hit.dist;
         //    	  H1(resprimex,"X residuals","x residuals",500,-100,500);
			      //   for(int z = -2000; z < 1000; z++)
			      //   {
			      //       if(mx!=-1e4&&my!=-1e4&&x0!=-1e4&&y0!=-1e4)
			      //       {
			      //           double xtarget = mx*(z-z0)+x0;// x= m*z+b, found from Analysis/PSI2018_SUmmer/muse_analysis_Sumer2018_4gem.C
			      //           bothXZ->Fill(xtarget,z);
						   //  circXZ->Fill(xtarget,z);
			      //       }
			      //   }
			      std::vector<double> y;
			      std::vector<double> x;

			      double yposition=-1e6;
			      double xposition=-1e6;
			        double xpos = mx*(pos.Y()-z0)+x0;
			        double resplus = (pos.X()-xpos) +hit.dist;
			        double resminus = (pos.X()-xpos) -hit.dist;
			        double last = 1e9;
			        double temp = 1e9;
			        double ressimp = -1e5;
			        int location = -1;
			        std::vector<double> strawoffset;
			        std::vector<double> chisqX;
			        std::vector<double> resX;
			       for(int kk = 0; kk < 3; kk++)
			       {
			      	//for(int ii = -100; ii < 100; ii++)
			      	//{
			      		yposition = pos.Y();
			      		xposition = pos.X()+kk*1;
			      		y.push_back(pos.Y());
			      		x.push_back(pos.X()+kk*1);
			      		strawoffset.push_back(kk);
			      		temp = sqrt(pow((mx*(yposition-z0)+x0)-xposition,2))-hit.dist;
			      		resX.push_back(temp);
			      		//tempplus = sqrt(pow(xposition -(mx*(yposition-z0)+x0),2))+hit.dist;

			      		if(temp<last)
			      			ressimp = temp;
			      		// else if(tempplus<tempminus && tempplus<last)
			      		// 	ressimp = tempplus;
						chisqX.push_back(pow(temp,2)/pow(0.15,2)/(numhitsX-2));
			      		// last = (tempminus < tempplus) ? tempminus : tempplus;
			      		if(side==0 && plane==0)
			      		H2(xposition,temp,"position V res","position V res",200,-200,200,200,-5,200);
			      		last = temp;	
			      	//}
			      }
			      	double miny = 1e9;
			      	for(int jj = 0; jj < chisqX.size(); jj++)
			      	{
			      		if(miny > chisqX.at(jj))
			      		{
			      			miny = chisqX.at(jj);
			      			location = jj;
			      		}
			      	}
			        // if (resminus<resplus)
			        // 	ressimp = resminus;
			        // else
			        // 	ressimp = resplus;
			        //if(plane==1&&side==2)

		       		//H2(ressimp,mx,Form("Residual V slope Side %i Plane%i",side,plane),"Residual from GEM track projected onto STT; residual (mm)",100,-5,100,100,-1,1);
		       		//H2(ressimp,hit.dist,Form("Residual V hit dist Side %i Plane%i",side,plane),"Residual from GEM track projected onto STT; residual (mm)",100,-5,100,100,0,5);
		       		H1(resX.at(location),Form("Residual from GEM track projected onto STT Side %i Plane%i",side,plane),"Residual from GEM track projected onto STT; residual (mm)",400,-5,100);
		       		H1(resX.at(location),"Total Residual","Total Residual",500,-5,100);
		       		if(side==0&&plane==0)
		       		H2(resX.at(location),strawoffset.at(location),"Residual V offset","Residual V offset",400,-5,100,20,-10,10);
		       		
		       		if(side==0&&plane==0)
		       		H2(x.at(location),pos.X(),"Plane 0 min x","Plane 0 min x",100,-1,-1,100,-1,-1);
		       		
		       		H1(pos.X()-xpos,"Unmodified X residual","Unmodified residual",100,-100,100);
			  }
		      //Y coordinate
		      if ((side==2||side==0)&& plane >4)
		      {
		        if(side==2 && plane==8 && (straw_in_plane==11||straw_in_plane<6))
		          continue;
		        auto pos=getStrawPos(hit.id)-offy;
		   //      for (int i =0; i<180;i++)// probably should do -90 to 90 but that means things have to be rewritten more intelligently for the negative angles.
		   //      {

		   //       TVector2 d=pos+TVector2(hit.dist,0).Rotate(2*i*M_PI/180);
		   //       bothYZ->Fill(d.X(),d.Y());
	    //      	 circYZ->Fill(d.X(),d.Y());
	    //      	 circ->Fill(d.X(),d.Y());
		   //      }
					// double linetempX = (pos.X()/my+pos.Y()- y0)/(my+1/my);//this is given by solving equation for two perpendicular intersecting lines 
     //        		double linetempY = my*linetempX+y0;
     //        		resprimey = sqrt(pow(pos.X()-linetempX,2)+pow(pos.Y()-linetempY,2))-hit.dist;
     //        		H1(resprimey,"Y residuals","y residuals",500,-100,500);
		   //      for(int z = -2000; z < 1000; z++)
		   //      {
		   //          if(mx!=-1e4&&my!=-1e4&&x0!=-1e4&&y0!=-1e4)
		   //          {
		   //              double ytarget = my*(z-z0)+y0;//y=m*z+b
					//     bothYZ->Fill(ytarget,z);
					//     circYZ->Fill(ytarget,z);
		   //      	}
		   //      }
		        double ypos = my*(pos.Y()-z0)+y0;
		        double resplus = (pos.X()-ypos) +hit.dist;
		        double resminus = (pos.X()-ypos) -hit.dist;
		        double ressimp = -1e5;
		        if (resminus<resplus)
		        	ressimp = resminus;
		        else
		        	ressimp = resplus;
		        H1(ressimp,"Y simp res","Y simp res",100,-100,100);
		      }
		}
		//break;
	}



	return ok;
}

Long_t GEM_STT::finalize()
{

	return ok;
}


Long_t GEM_STT::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new GEM_STT(in,out,inf_,outf_,p);
}
}


ClassImp(GEM_STT);

