#include <hough.h>


void hough::getAngle(int *numhitsX, int*numhitsY, int maxangles, std::vector<double> &xangle,std::vector<double> &yangle,
	std::vector<double> &xrho, std::vector<double> &yrho)
{
  for (auto hit:STT->hits)
  {
    //ret=Plugin::stop;
    int side, plane, straw_in_plane;
    STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
    //X coordinate
    if ((side==2||side==0)&& plane <5 )//&& hit.dist>0.6 && hit.dist<4.0)
    {
      // if(side==2 && (plane==2||plane==3||plane==4))
      //   continue;
      // if(side==0 && (plane==1||plane==2||plane==3))
      //   continue;
      // if(side==0 && plane == 0 && straw_in_plane==33)
      //   continue;

          *numhitsX = *numhitsX+1;
      auto pos=getLocalPos(hit.id,offsetx,offsety);
      TVector2 xpos(pos.X(),pos.Z());
      for (int i =0; i<180;i++)// probably should do -90 to 90 but that means things have to be rewritten more intelligently for the negative angles.
      {

        TVector2 d=xpos+TVector2(hit.dist,0).Rotate(2*i*M_PI/180);
        //H2(d.X(),d.Y(),"Angle/XZ circle","XZ Circle",501,-1000,1000,4001,-600,600);
        // circX->Fill(d.X(),d.Y());
        // bothX->Fill(d.X(),d.Y());
        // circ->Fill(d.X(),d.Y());
        // lineX->Fill(d.X(),d.Y());
        // min->Fill(d.X(),d.Y());
          
            double x=xpos.Rotate(i*M_PI/180).X();
        for(int j =-3;j<3;j++)
        {
            //double r = gRandom->Gaus(0,0.5);//random smearing. This makes regions with a high density of Hough points stand out in the search for best Hough points
            //double r1 = gRandom->Gaus(0,0.7);//ETHAN added all of the r and r1 nonsense for smearing. If you think it should be removed, go for it.
            //houghplotX->Fill(i,(x+hit.dist));
            houghplotX->Fill(i,(x+hit.dist+j),TMath::Gaus(j,0,3,0));
            //houghplotX->Fill(i,(x+hit.dist/10-1),0.7);

          if (trunc(x+hit.dist+0.5)!=trunc(x-hit.dist+0.5)) 
          {
            houghplotX->Fill(i,x-hit.dist+j,TMath::Gaus(j,0,3,0));
          }
            }
        } 
      }
    //Y coordinate
    if ((side==2||side==0)&& plane >4 )//&& hit.dist/10>0.6 && hit.dist<4.0)
    {
      // if(side==2 && plane==8 && (straw_in_plane==11||straw_in_plane<6))
      //     continue;
      *numhitsY=*numhitsY+1;
      auto pos=getLocalPos(hit.id,offsetx,offsety);
      TVector2 ypos(pos.Y(),pos.Z());
      for (int i =0; i<180;i++)// probably should do -90 to 90 but that means things have to be rewritten more intelligently for the negative angles.
      {

        TVector2 d=ypos+TVector2(hit.dist,0).Rotate(2*i*M_PI/180);
        //H2(d.X(),d.Y(),"Angle/YZ circle","YZ Circle",501,-1000,1000,4001,-600,600);
        // circY->Fill(d.X(),d.Y());
        // bothY->Fill(d.X(),d.Y());
        // circ->Fill(d.X(),d.Y());
        // lineY->Fill(d.X(),d.Y());

        double y=ypos.Rotate(i*M_PI/180).X();
        for(int j =-3;j<3;j++)
        {
          houghplotY->Fill(i,(y+hit.dist+j),TMath::Gaus(j,0,3,0));

          if (trunc(y+hit.dist+0.5)!=trunc(y-hit.dist+0.5)) 
          {
            houghplotY->Fill(i,y-hit.dist+j,TMath::Gaus(j,0,3,0));
          }
        }
      }
    }
  }

  int x1,y1,z1, x2,y2,z2;

  //find the first maxangles number of hough maxima
  for(int i = 0; i < maxangles; i++)
  {
	houghplotX->GetMaximumBin(x1,y1,z1);
	xangle.push_back(houghplotX->GetXaxis()->GetBinCenter(x1));
	xrho.push_back(houghplotX->GetYaxis()->GetBinCenter(y1));
	houghplotX->SetBinContent(x1,y1,-5);

	houghplotY->GetMaximumBin(x2,y2,z2);
	yangle.push_back(houghplotY->GetXaxis()->GetBinCenter(x2));
	yrho.push_back(houghplotY->GetYaxis()->GetBinCenter(y2));
	houghplotY->SetBinContent(x2,y2,-5);
  }



  return;
}