#include <hough.h>

#include <iostream>
#include <cmath>


hough::hough(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

hough::~hough()
{
};


Long_t hough::startup()
{
  STT = NULL;
  getOutBranchObject("StrawTubeHits",(TObject **) &STT);
  if(!STT) getBranchObject("StrawTubeHits",(TObject **) &STT); 
  if(!STT) debug(0,"Could not find STT hits in file\n");
    SPSHits = NULL;
  getOutBranchObject("SPSHits",(TObject **) &SPSHits);
  if(!SPSHits) getBranchObject("SPSHits",(TObject **) &SPSHits); //Why are we trying again?
  if(!SPSHits) debug(0,"Could not find BH hits in file\n");

    theRunInfo=(MRTRunInfo *) getFileObject("RunInfo");

  Tracks = new TrackHits;
  makeBranch("TrackHits",(TObject **) &Tracks);

  houghplotX=dH2("houghplot XZ","houghplot XZ;theta (degree);rho (mm)",180,-0.5,179.5,601,-600,600);
  houghplotY=dH2("houghplot YZ","houghplot YZ;theta (degree);rho (mm)",180,-0.5,179.5,1201,-600,600);

  circX=dH2("circ XZ","circ XZ",501,-1000,1000,4001,-600,600);
  circY=dH2("circ YZ","circ YZ",501,-1000,1000,4001,-600,600);
  circ=dH2("circ","circ",501,-1000,1000,4001,-600,600);


  lineX=dH2("line XZ","line XZ",151,-500,500,401,-600,600);
  lineY=dH2("line YZ","line YZ",151,-500,500,401,-600,600);

  lineprime=dH2("line prime","line prime",151,100,250,401,-40,40);

  min=dH2("Minimized Tracks XZ","Tracks XZ;Pos X (mm);Pos Z (mm)",501,-1000,1000,4001,-600,600);

  bothX=dH2("Tracks XZ","Tracks XZ;Pos X (mm);Pos Z (mm)",501,-1000,1000,4001,-600,600);
  bothY=dH2("Tracks YZ","Tracks YZ;Pos Y (mm);Pos Z (mm)",501,-1000,1000,4001,-600,600);

  hitresidualX=dH1("Hit Residual XZ","Hit Residual XZ; Residual (mm);counts",201,-10,10);
  hitresidualnormX=dH1("Hit Residual normal XZ","Hit Residual XZ; Residual (mm);counts",201,-10,10);

  hitresidualY=dH1("Hit Residual YZ","Hit Residual YZ; Residual (mm);counts",201,-10,10);
  hitresidualnormY=dH1("Hit Residual normal YZ","Hit Residual YZ; Residual (mm);counts",201,-10,10);


  hitchisq=dH1("Hit Chi Squared","Hit Chi Squared",401,-40,40);

  straws=dH2("Plane","Plane",51,150,200,41,-20,20);
  AngleX=dH1("XZ Angle","Angle XZ",180,-0.5,179.5);
  FinalAnglex=dH1("Final XZ Angle","Angle XZ",360,-179.5,179.5);
  FinalAngley=dH1("Final YZ Angle","Angle YZ",360,-179.5,179.5);


  AngleY=dH1("YZ Angle","Angle YZ",180,-0.5,179.5);
  chiVangle=dH2("Chi sq V Angle","Chi sq V angle;angle (degree);chisq",180,-0.5,179.5,500,0,500);
  chiVangleprime=dH2("Chi sq V Angle prime","Chi sq V angle prime;angle (degree);chisq",180,-0.5,179.5,500,0,500);

  //needed for getting GDML
  // viscotab    = addTab("hough");
  // if(viscotab)
    TEveManager::Create(kFALSE);

  GenerateDet();

  return 0;
  
};


Long_t hough::process()
{

  int ret=0;
  //double angle = 60;//0 for alignment, 60 for scattering
  int LV = 0;
  int LH = 0;
  int RV = 0;
  int RH = 0;

  TrackHit outtrack;
  Tracks->clear();
  outtrack.chisq=-100000;
  outtrack.mx=-100000;
  outtrack.my=-100000;
  outtrack.x0=-100000;
  outtrack.y0=-100000;
  outtrack.z0=-100000;

  outtrack.interceptzx=-100000;
  outtrack.interceptzy=-100000;

  // for(auto hit:STT->hits)
  // {
  //   int side, plane, straw_in_plane;
  //   STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
  //   // if((side==0 && plane <5) || (side==2 && plane <5))
  //   //   LV++;
  //   // if((side==1 && plane <5) || (side==3 && plane <5))
  //   //   RV++;
  //   // if((side==0 && plane >4) || (side==2 && plane >4))
  //   //   LH++;
  //   // if((side==1 && plane >4) || (side==3 && plane >4))
  //   //   RH++;
    
  //   H1(plane,Form("%i side plane multiplicity",side),Form("%i side plane multiplicity",side),10,-0.5,9.5);
  // }
  
  //Don't comment out these commands, code won't run without them.
  houghplotX->Reset();
  houghplotY->Reset();


  //Be careful with ROOT Reset commands below
  //Not sure exactly how it does the reset of the histogram
  //but it will slow this code down by a factor of 5-10
  //At least when the ones below are used
  // circ->Reset();
  // lineX->Reset();
  // bothX->Reset();
  // lineY->Reset();
  // bothY->Reset();
  // min->Reset();

  // if (LV<3||LH<3||LV>10||LH>10)
  // {
  // //std::cout << "Less than 3 hits in STT, skipping event" << std::endl;
  //  Tracks->tracks.emplace_back(std::move(outtrack));
  // return ok;
  // }


  int numhitsX =0;
  int numhitsY =0;
  std::vector<double> angledx;
  std::vector<double> angledy;
  std::vector<double> rhox;
  std::vector<double> rhoy;
  static int maxangles =10;//how many hough maxima to look through

  //this finds the hough transform maxima
  getAngle(&numhitsX,&numhitsY,maxangles,angledx,angledy,rhox,rhoy);

  if (numhitsX<3||numhitsY<3)
  {
    Tracks->tracks.emplace_back(std::move(outtrack));
    return ok;
  }


  int x1,y1,z1, x2,y2,z2;
  double res = 100000;
  double chisq = 100000;
  double xangle=angledx.at(0);
  double xrho = rhox.at(0);
  double yangle=angledy.at(0);
  double yrho = rhoy.at(0);


  AngleX->Fill(xangle);//hopefully a sensible coordinate system
  AngleY->Fill(yangle);//y coordinate, angle

  // for (int i=-400;i<400;i++)
  // {
  //  auto px=TVector2(xrho,i).Rotate(-x1*M_PI/180);
  //  lineX->Fill(px.X(),px.Y());
  //  bothX->Fill(px.X(),px.Y());
  //  auto py=TVector2(yrho,i).Rotate(-x2*M_PI/180);
  //  lineY->Fill(py.X(),py.Y());
  //  bothY->Fill(py.X(),py.Y());
  // }

  double chisqX = 100000;
  double chisqY = 100000;

  std::vector<double> chix;
  // std::vector<double> *par0x = new std::vector<double>;
  // std::vector<double> *par1x = new std::vector<double>;
  std::vector<double> chiy;
  // std::vector<double> *par1y = new std::vector<double>;
  // std::vector<double> *par0y = new std::vector<double>;

  // std::vector<double> *residualX = new std::vector<double>;
  // std::vector<double> *residualY = new std::vector<double>;

  //while loop finds all hough maxima
  for(int i = 0; i<maxangles;i++)
  {
    chisqX=0;
    chisqY = 0;
    //lineX->Reset();
    //lineY->Reset();
    xangle = angledx.at(i);
    yangle = angledy.at(i);
    xrho = rhox.at(i);
    yrho = rhoy.at(i);
    // for (int i=-200;i<200;i++)
    // {
    //   auto px=TVector2(xrho,i).Rotate(-x1*M_PI/180);
    //   lineX->Fill(px.X(),px.Y());
    //   auto py=TVector2(yrho,i).Rotate(-x2*M_PI/180);
    //   lineY->Fill(py.X(),py.Y());
    // }
      
    double xslope = TMath::Tan(xangle*M_PI/180);
    double xint = xrho/TMath::Cos(xangle*M_PI/180);

    double yslope = TMath::Tan(yangle*M_PI/180);
    double yint = yrho/TMath::Cos(yangle*M_PI/180);

    //finding residuals and chi sq
    for (auto hit:STT->hits)
    {
        //ret=Plugin::stop;
        int side, plane, straw_in_plane;
        STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
        if ((side==2||side==0) && plane <5 )//&& hit.dist>0.6 && hit.dist<4.0)
        {
          // if(side==2 && (plane==2||plane==3||plane==4))
          //   continue;
          // if(side==0 && (plane==1||plane==2||plane==3))
          //   continue;
          // if(side==0 && plane == 0 && straw_in_plane==33)
          //   continue;
          res = getResidual(hit.id,hit.dist,xslope,xint);
          hitresidualnormX->Fill(res);
          chisqX = chisqX + (pow(res,2)/pow(0.15,2))/(numhitsX-2);
        }
        if ((side==2||side==0) && plane >4 )//&& hit.dist>0.6 && hit.dist<4.0)
        {
          res = getResidual(hit.id,hit.dist,yslope,yint);
          hitresidualnormY->Fill(res);
          chisqY = chisqY + (pow(res,2)/pow(0.15,2))/(numhitsY-2);
        }
    }
    //fill the output arrays
    chiVangle->Fill(xangle,chisqX);
    chix.push_back(chisqX);
    chiy.push_back(chisqY);
    // par0x->push_back(xslope);
    // par1x->push_back(xint);
    // par0y->push_back(yslope);
    // par1y->push_back(yint);
  }


  int locx = -1;
  int locy = -1;
  //find the min chi-sq and corresponding angle/rho
  locx = getMinLocation(chix);
  locy = getMinLocation(chiy);

  if(locx==-1 || locy==-1)
  {
    Tracks->tracks.emplace_back(std::move(outtrack));
    std::cout << "Chisq too big" << std::endl;
    return ok;
  }


  // H1(angledx->at(locx),"X angle minimized","X angle minimized",180,-0.5,179.5);
  // H1(angledy->at(locy),"Y angle minimized","Y angle minimized",180,-0.5,179.5);
  // double xdiff = par0x->at(locx)-1/TMath::Tan(angledx->at(locx)*M_PI/180);
  // H1(xdiff,"xdiff","xdiff",10,-1,1);
  // double angleprimex = TMath::ATan(par0x->at(locx))*180/M_PI;
  // H1(angleprimex,"X angle minimized from slope","X angle minimized from slope",180,-0.5,179.5);
  // double angleprimey = TMath::ATan(par0y->at(locy))*180/M_PI;
  // H1(angleprimey,"Y angle minimized from slope","Y angle minimized from slope",180,-0.5,179.5);

  // for(int i=-400;i<400;i++)
  // {
  // //auto p=TVector2(dy,i).Rotate(-x*M_PI/180);
  // //both->Fill(p.X(),p.Y());
  //  min->Fill(i,par0x->at(locx)*i+par1x->at(locx));
  // }

// //Do it again but use the starting point of fit as the optimal value found above.
// //Hough Transform isn't perfect for us (only accurate to level of bin size) so have to kind of look around optimal guessed value.

  std::vector<double> chiprimex;
  std::vector<double> par0primex;
  std::vector<double> par1primex;
  std::vector<double> chiprimey;
  std::vector<double> par0primey;
  std::vector<double> par1primey;

  int counter1=0;
  double chisqprimex=100000;
  double resprimex=100000;
  double chisqprimey=100000;
  double resprimey=100000;

  double gran = 0.1;//granularity of search around best hough params
  //chiVangleprime->Reset();
  //lineprime->Reset();
  double xangleprime = angledx.at(locx);
  double yangleprime = angledy.at(locy);
  double xrhoprime = rhox.at(locx);
  double yrhoprime = rhoy.at(locy);

  for(int i = -50;i<50;i++)
  {
    for(int j=-50;j<50;j++)
    {
      double xslope = TMath::Tan(  (xangleprime +gran*i) *M_PI/180);
      double xint = (xrhoprime+gran*j)/(TMath::Cos((xangleprime+gran*i)*M_PI/180));
      double yslope = TMath::Tan(  (yangleprime +gran*i) *M_PI/180);
      double yint = (yrhoprime+gran*j)/(TMath::Cos((yangleprime+gran*i)*M_PI/180));

      chisqprimex = 0;
      chisqprimey = 0;
      for (auto hit:STT->hits)
      {
        //ret=Plugin::stop;
        int side, plane, straw_in_plane;
        STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
        if ((side==2||side==0) && plane <5 )//&& hit.dist>0.6 && hit.dist<4.0)
        {
          // if(side==2 && (plane==2||plane==3||plane==4))
          //   continue;
          // if(side==0 && (plane==1||plane==2||plane==3))
          //   continue;
          // if(side==0 && plane == 0 && straw_in_plane==33)
          //   continue;
          resprimex = getResidual(hit.id,hit.dist,xslope,xint);
          chisqprimex = chisqprimex + pow(resprimex,2)/pow(0.15,2)/(numhitsX-2);
        }
        if ((side==2||side==0) && plane >4 )//&& hit.dist>0.6 && hit.dist<4.0)
        {
          resprimey = getResidual(hit.id,hit.dist,yslope,yint);
          chisqprimey = chisqprimey + pow(resprimey,2)/pow(0.15,2)/(numhitsY-2);
        }
      }

      //std::cout <<" optimizing chi sq: " << chisqprime << std::endl;
      chiprimex.push_back(chisqprimex);
      par0primex.push_back(xslope);
      par1primex.push_back(xint);
      chiprimey.push_back(chisqprimey);
      par0primey.push_back(yslope);
      par1primey.push_back(yint);

      //chiVangleprime->Fill(angledx->at(loc)+gran*i,chisqprime);
    }//end j
  }//end i

  locx =-1;
  locy =-1;

  locx = getMinLocation(chiprimex);
  locy = getMinLocation(chiprimey);
  

  if(locx==-1||locy==-1)
  {
    Tracks->tracks.emplace_back(std::move(outtrack));
    std::cout <<"Chi sq prime too big" << std::endl;
    return ok;
  }


  double anglex = TMath::ATan(par0primex.at(locx))*180/M_PI;
  FinalAnglex->Fill(anglex);
  double angley = (TMath::ATan(par0primey.at(locy))*180/M_PI);
  FinalAngley->Fill(angley);



  chisqX = 0;
  chisqY = 0;
  for (auto hit:STT->hits)
  {
    //ret=Plugin::stop;
    int side, plane, straw_in_plane;
    STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
    if ((side==2||side==0) && plane <5 )//&& hit.dist>0.6 && hit.dist<4.0)
    {
      //   if(side==2 && (plane==2||plane==3||plane==4))
      //     continue;
      // if(side==0 && (plane==1||plane==2||plane==3))
      //     continue;
      // if(side==0 && plane == 0 && straw_in_plane==33)
      //     continue;
      res = getResidual(hit.id,hit.dist,par0primex.at(locx),par1primex.at(locx));
      chisqX = chisqX + pow(res,2)/pow(0.15,2)/(numhitsX-2);
      //H2(res,anglex,"Angle X v residual","Angle X v residual",201,-10,10,180,-0.5,179.5);
      hitresidualX->Fill(res);  
    }
    if ((side==2||side==0) && plane >4 )//&& hit.dist>0.6 && hit.dist<4.0)
    {
      res = getResidual(hit.id,hit.dist,par0primey.at(locy),par1primey.at(locy));
      chisqY = chisqY + pow(res,2)/pow(0.15,2)/(numhitsX-2);
      //H2(res,angley,"Angle Y v residual","Angle Y v residual",201,-10,10,180,-0.5,179.5);
      hitresidualY->Fill(res);  
    }
  }


  // for (int i=-600;i<600;i++)
  // {
  //  //auto px=TVector2(xrho,i).Rotate(-x1*M_PI/180);
  //  bothX->Fill(i,(1/par0primex->at(locx))*i-par1primex->at(locx)/par0primex->at(locx));
  //  // auto py=TVector2(yrho,i).Rotate(-x2*M_PI/180);
  //  bothY->Fill(i,(1/par0primey->at(locy))*i-par1primey->at(locy)/par0primey->at(locy));
  // }

  //solve x=mz+b for x =0 to get z
  double z = -par1primex.at(locx)/par0primex.at(locx);
  //plug that z into y=mz+b to get y
  double y = par0primey.at(locy)*z+par1primey.at(locy);
  TVector3 direction(par0primex.at(locx),par0primey.at(locy),1);

  direction.RotateY(angle*M_PI/180);
  TVector3 position(0,y,z);//we use x = 0 becuase that is the line that divides the experiment between left and right
  position.RotateY(angle*M_PI/180);


  //output tree
  hitchisq->Fill(chisqX);
  outtrack.mx=direction.X();
  outtrack.my=direction.Y();
  outtrack.mz=direction.Z();
  outtrack.x0=position.X();
  outtrack.y0=position.Y();
  outtrack.z0=position.Z();
  Tracks->tracks.emplace_back(std::move(outtrack));


  chix.clear();
  chiy.clear();
  // par0x->clear();
  // par1x->clear();
  angledx.clear();
  // residualX->clear();
  // par0y->clear();
  // par1y->clear();
  angledy.clear();
  // residualY->clear();
  // rhox->clear();
  // rhoy->clear();

  chiprimex.clear();
  par0primex.clear();
  par1primex.clear();
  chiprimey.clear();
  par0primey.clear();
  par1primey.clear();
      return ret;
}

Long_t hough::finalize()
{
  // TF1 *g1 = new TF1("g1","gaus",-0.40,.35);
  // TCanvas *c1 = new TCanvas();

  // hitresidualX->Draw();
  // hitresidualX->Fit("g1");
  // c1->Update();
  // std::string plottitle = "Plots/Tracking/residualX"+std::to_string(theRunInfo->runNumber)+".root";
  // c1->SaveAs(plottitle.c_str());


  return 0; // 0 = all ok
};

Long_t hough::cmdline(char *cmd)
{

  return 0; // 0 = all ok
};




extern "C"{
  Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
  {
    return (Plugin *) new hough(in,out,inf_,outf_,p);
  }
}
