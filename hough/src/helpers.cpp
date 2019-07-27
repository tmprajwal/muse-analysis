#include <hough.h>

void hough::rotAngle(double x)
{
  angle = x;
}

double hough::getResidual(int id,double dist, double slope, double intercept)
{

    auto pos=getLocalPos(id,offsetx,offsety);//this is wire position x,y
    //now need intercept of line perpendicular to fit that passes through the straw
    //the x coordinate of intercept is given by
    int side, plane, straw_in_plane;
  	STT_internal_to_logic(id,&side,&plane,&straw_in_plane);
  	double res=-100;
  	if(plane<5)
  	{
    	double X = (pos.Z()/slope+pos.X()- intercept)/(slope+1/slope);//this is given by solving equation for two perpendicular intersecting lines 
    	double Y = slope*X+intercept;
    	res = sqrt(pow(pos.Z()-X,2)+pow(pos.X()-Y,2))-dist;
    }
    else
    {
		double X = (pos.Z()/slope+pos.Y()- intercept)/(slope+1/slope);//this is given by solving equation for two perpendicular intersecting lines 
    	double Y = slope*X+intercept;
    	res = sqrt(pow(pos.Z()-X,2)+pow(pos.Y()-Y,2))-dist;
    }
	return res;
}

int hough::getMinLocation(std::vector<double> input)
{
  int res = -1;
  int minimum = 10000000000000;//Yes this is stupidly large. Some values of chisq for hough transform are also stupidly large, so we need this.
  for(int k=0;k<input.size();k++)
  {
    if(minimum>input.at(k))
    {
        minimum=input.at(k);
      res=k;
    }
  }
  return res;
}


//this returns the Straw position in the STT local frame
TVector3 hough::getLocalPos(int id, int offsetx[2], int offsety[2])
{  
  TVector3 res;
  int side, plane, straw_in_plane;
  STT_internal_to_logic(id,&side,&plane,&straw_in_plane);
  
  //printf("side %i plane %i straw %i\n",side,plane,straw_in_plane);
  if(side==2)
  {
    double loc[3] = {0,0,0};
    double mas[3] = {0,0,0};
    STTL190[plane][straw_in_plane]->GetMatrix()->LocalToMaster(loc,mas);
    res.SetXYZ(mas[0]*10 + offsetx[1],mas[1]*10 + offsety[1],mas[2]*10);
    res.RotateY(-60*M_PI/180);
    //std::cout << "GDML " << loc[0]*10 << " " <<loc[1]*10 << " "<< loc[2]*10 <<std::endl;
  }
  if(side==0)
  {
    double loc[3] = {0,0,0};
    double mas[3] = {0,0,0};
    STTL160[plane][straw_in_plane]->GetMatrix()->LocalToMaster(loc,mas);
    res.SetXYZ(mas[0]*10+offsetx[0],mas[1]*10+offsety[0],mas[2]*10);
    res.RotateY(-60*M_PI/180);

    //std::cout << "GDML " << mas[0]*10 << " " <<mas[1]*10 << " "<< mas[2]*10 <<std::endl;
  }
  return res;
}


//this returns the straw position in the global frame
TVector2 hough::getStrawPos(int id)//angle is 60 at scattering and 0 for alignment
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
  if(side==2)
  {
    double loc[3] = {0,0,0};
    double mas[3] = {0,0,0};
    STTL190[plane][straw_in_plane]->GetMatrix()->LocalToMaster(loc,mas);
    if(plane<5)
    {
      res.Set(mas[0]*10,mas[2]*10);
    }
    else
      res.Set(mas[1]*10,mas[2]*10);
    //std::cout << "GDML " << mas[0]*10 << " " <<mas[1]*10 << " "<< mas[2]*10 <<std::endl;
  }
  if(side==0)
  {
    double loc[3] = {0,0,0};
    double mas[3] = {0,0,0};
    STTL160[plane][straw_in_plane]->GetMatrix()->LocalToMaster(loc,mas);
    if(plane<5)
    {
      res.Set(mas[0]*10,mas[2]*10);
    }
    else
      res.Set(mas[1]*10,mas[2]*10);
    //std::cout << "GDML " << mas[0]*10 << " " <<mas[1]*10 << " "<< mas[2]*10 <<std::endl;
  }

  // if (side==2)
  // {
  //   zoffset = 493-(10.*distance_between_planes);
  //   double pos=(straw_in_plane)*pitch+(plane % 2)*pitch/2;
  //   res.Set(pos-15*pitch,zoffset+plane*distance_between_planes);
  //   if(plane<5)
  //   {
  //     res = res.Rotate(-angle*(M_PI/180));//-60 at scattering position
  //     std::cout << "HOUG " << res.X() << " " << res.Y() << std::endl;
  //   }
  // }
  // if (side==0)
  // {
  //   zoffset = 493-(10.*distance_between_planes+55+10.*distance_between_planes);//
  //   double pos=straw_in_plane*pitch+(plane % 2)*pitch/2;
  //   res.Set(pos,zoffset+plane*distance_between_planes);
  //   if(plane<5)
  //   res = res.Rotate(-angle*(M_PI/180));//-60 at scattering position
  // }

  // if(side == 1)
  // {
  //     zoffset = 0.;
  //     double pos = straw_in_plane*pitch +(plane % 2)*pitch/2;//original plus
  //     res.Set(pos,zoffset+plane*distance_between_planes);
  //     res = res.Rotate(angle*(M_PI/180));//-60 at scattering position
  // }
  // if (side==3)
  // {
  //   zoffset = 10.;//
  //   double pos=straw_in_plane*pitch+(plane % 2)*pitch/2;
  //   res.Set(pos,zoffset+plane*distance_between_planes);
  //   res = res.Rotate(angle*(M_PI/180));//-60 at scattering position
  // }
  return res;
};