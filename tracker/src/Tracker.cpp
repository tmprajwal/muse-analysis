#include <Tracker.h>
#include "StrawTubetree.h" //for internal_to_logic
#include <TEveEventManager.h>

Tracker::Tracker(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p){
  STTangle = 0;
}
Tracker::~Tracker(){}
Long_t Tracker::setSTTangle(double angle)
{
  STTangle = angle;
  return ok;
}
Long_t Tracker::startup()
{
  new TGeoManager("Geometry", "Geane geometry");
  TGeoManager::Import("genfitGeom.root");

  tab=addTab("GenFit");
  USE_DISPLAY=tab!=NULL;
  //  stopAtEvent(1000);
  // useMC();

  // Setup fitter conditions
  setup_fitter();

  STT = NULL;
  getOutBranchObject("StrawTubeHits",(TObject **) &STT);
  if(!STT) getBranchObject("StrawTubeHits",(TObject **) &STT); //Why are we trying again?
  if(!STT) debug(0,"Could not find STT hits in file\n");
        


  // =========================================================
  // ------------------- GEM initialisation ------------------
  // =========================================================
  if (USE_MC) nGEMS = 10;
  else        nGEMS = 3;

  if(USE_MC)
    // create some GEM positions for now
    for (int i = 0; i < nGEMS; i++) GEMpositions.push_back(100*i);
  else {
    // Attempt to get GEM hits
    GEMs = NULL;
    getOutBranchObject("LumiGEMhits",(TObject **) &GEMs);
    if(!GEMs) getBranchObject("LumiGEMhits",(TObject **) &GEMs); 
    if(!GEMs) debug(0,"Could not find GEM clusters in file\n");
    else debug(0,"Including GEMs into tracker\n");

    // GEMpositions.push_back(1841);
    // GEMpositions.push_back(1941);
    // GEMpositions.push_back(2061);
    GEMpositions.push_back(-71.0);//For June 2017 Beamtime
    GEMpositions.push_back(-65.0);
    GEMpositions.push_back(-56.0);

  }

  if (USE_DISPLAY) display->makeGui(false);
  return ok;
}

Long_t Tracker::setup_fitter()
{
  genfit::MaterialEffects::getInstance()->setNoEffects();                         // Turning off Material effects (for now)
  genfit::FieldManager::getInstance()->init(new genfit::ConstField(0. ,0., 0.));  // No mag field
  if(USE_DISPLAY) display = genfit::EventDisplay::getInstance();
  return ok;
}


std::pair<TVector3,TVector3> STT_getStrawPos(int id, double STTangle)
{
  // ask gdml or something, for now, make it up.
  
  std::pair<TVector3,TVector3> res;
  int side, plane, straw_in_plane;
  STT_internal_to_logic(id,&side,&plane,&straw_in_plane);

  if (side==2)
  {
    double tube_diameter = 10;
    double distance_between_planes = 8.7;
    double pitch=10.1; //10.1 mm spacing
    double zoffset=280;
      
    double pos=(straw_in_plane-8)*pitch-(plane % 2)*pitch/2;
    res.first[0]=-125;
    res.first[1]=pos;
    res.first[2]=zoffset+plane*distance_between_planes;
    res.second[0]=125;
    res.second[1]=pos;
    res.second[2]=zoffset+plane*distance_between_planes;
  }

  if(side ==1)
  {
    double tube_diameter = 10.;
    double distance_between_planes = 8.7;
    double pitch = 10.1;
    double pos = 0;
    double zoffset =0;
    if(plane < 5) //these planes are vertical
    {
      zoffset = 300.;
      pos = (straw_in_plane-8.)*pitch -(plane % 2)*pitch/2;//calculates wire position and puts center wire at 0,0,0
      res.first[1] = -305.; //this is half of the height of these planes, assuming centered on beam roughly
      res.first[0] = pos;
      res.first[2] = zoffset+plane*distance_between_planes;
      res.second[1] = 305.;//half height of planes again
      res.second[0] = pos;
      res.second[2] = zoffset+plane*distance_between_planes;
    }
    if(plane > 4 && plane < 10)//these planes are horizontal
    {
      zoffset = 380.;//distance between front and rear set of planes
      pos = (straw_in_plane-8.)*pitch +(plane % 2)*pitch/2;//makes sure center wire is at beam center
      res.first[1] = pos;
      res.first[0] = -305.; //this is half of the width of these planes, assuming centered on beam roughly
      res.first[2] = zoffset+plane*distance_between_planes;
      res.second[1] = pos;
      res.second[0] = 305.;//half width again
      res.second[2] = zoffset+plane*distance_between_planes;
    }
    else
    {
      printf("This plane should not exist plane %i side %i\n",plane,side);
    }
    
  }
  res.first.RotateY(STTangle*M_PI/180);
  res.second.RotateY(STTangle*M_PI/180);

  return res;
}


Long_t Tracker::process()
{
  if(MAX_EVENT && in->GetReadEvent() > nMAX_EVENT) return ok;
  if (USE_DISPLAY) display->resetEvent();
  // todo add handle RETURN
  process_beamlinetrack("GEMs_HighAmplitudeHits");
  //    process_beamline("GEMs_AllHits");

  process_scatteredtrack("scatteredTrack");

  //process_track();//For June 2017 beamtime with GEMs and STTs in a line

  if (USE_DISPLAY) {
    if (gEve->GetCurrentEvent()) 
      gEve->GetCurrentEvent()->DestroyElements();
    if(scatteredtrack->getNumPoints()>0)// scatteredtrack->getNumPoints()>0)
    {
      display->drawEvent(0);
      std::cout << "Track at Event: " << in->GetReadEvent() << std::endl;
    }
  }

  return ok;

}

Long_t Tracker::process_track()
{
  genfit::AbsKalmanFitter* fitter = new genfit::DAF();//KalmanFitter(20,0.001);
  mastertrack = new genfit::Track(new genfit::RKTrackRep(11), TVector3(0,0,-1700), TVector3(0,0,1));

  TVectorD hitCoords(2);
  const int detId(0); // detector ID
  int planeId(0); // detector plane ID
  int hitId(0); // hit ID
  //  Gather hits in STTs
  // Create a hit in three planes of STTs
  if (STT->hits.size()<3) return ok;
  for (auto hit:STT->hits)
    {   
        auto wire=STT_getStrawPos(hit.id,STTangle);
        int side, plane, straw_in_plane;
        STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
        printf("Which side: %i plane: %i straw: %i distance: %f\n",side, plane, straw_in_plane,hit.dist);
        printf("Wire start %g %g %g Wire end %g %g %g\n\n",wire.first[0],wire.first[1],wire.first[2],wire.second[0],wire.second[1],wire.second[2]);

        auto measurement = new genfit::WireMeasurementNew(hit.dist/10,hit.disterr/10,wire.first,wire.second,side, ++ hitId, NULL);
        measurement->setLeftRightResolution(0);//Sets a preferred side to wires, >0 is right side, <0 is left side.
        measurement->setMaxDistance(5);
        mastertrack->insertPoint(new genfit::TrackPoint(measurement, mastertrack));
    }
  //Gather hits in GEMs
  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
  double detectorResolution(0.1); // resolution of GEMs detectors
  hitCov *= detectorResolution*detectorResolution;

  // Loop over planes
  for (int GEMplane = 0; GEMplane < 3; GEMplane++)
  {
    for (auto hit:GEMs->hits)
    {
      if ((hit.GEMid - 3) == GEMplane)
      {
        //currently this changes coordinate system from native GEMs to a more reasonable world coordinate system
        //however this has "xl" being the vertical axis and "yl" being horizontal
        hitCoords[0] = -(hit.xl)*0.4 + 50; // Turn strip into mm positions (0.4mm per strip) and centers the GEM coord system
        hitCoords[1] = -(hit.yl)*0.4 + 50; // Turn strip into mm positions (0.4mm per strip)
        double z = (GEMpositions[GEMplane]);

        genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
        measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,z), TVector3(1,0,0), TVector3(0,1,0))), planeId++);
        mastertrack->insertPoint(new genfit::TrackPoint(measurement, mastertrack));

        break;
      }
    }
  }

  // Fit the beamline track
  if( mastertrack->getNumPoints() == 0) return ok;  // nothing to fit
  assert(mastertrack.checkConsistency());           // check
  fitter->processTrack(mastertrack);                // do the fit
  assert(mastertrack.checkConsistency());           // check
  // std::cout << beamlinetrack->getFitStatus()->getChi2() << std::endl;

  //Set event display:
  if(USE_DISPLAY) display->addTrack(mastertrack);

  return ok;
}
// ==============================================================================================================
// ------------------------------------------ Beamline Track Subroutines ----------------------------------------
// ==============================================================================================================


Long_t Tracker::process_beamlinetrack(std::string option)
{

  genfit::AbsKalmanFitter* fitter = new genfit::DAF();//KalmanFitter(20,0.001);
  beamlinetrack = new genfit::Track(new genfit::RKTrackRep(11), TVector3(0,0,-1700), TVector3(0,0,1));

  // Grab GEMs, or MC if possible
  if((GEMs) && (!USE_MC)) GatherHits_inGEMs(option);
  if(USE_MC)              CreateHits_inGEMs_MC(option);

  //  GatherHits_inSiPM();

  // Fit the beamline track
  if( beamlinetrack->getNumPoints() == 0) return ok;  // nothing to fit
  assert(beamlinetrack.checkConsistency());           // check
  fitter->processTrack(beamlinetrack);                // do the fit
  assert(beamlinetrack.checkConsistency());           // check
  PlotGEMs(option);
  // std::cout << beamlinetrack->getFitStatus()->getChi2() << std::endl;

  //Set event display:
  if(USE_DISPLAY) display->addTrack(beamlinetrack);

  return ok;
}

// ==============================================================================================================
// ------------------------------------------- Scatter Track Subroutines ----------------------------------------
// ==============================================================================================================


Long_t Tracker::process_scatteredtrack(std::string option)
{

  genfit::AbsKalmanFitter* fitter = new genfit::DAF();
  //genfit::AbsKalmanFitter* fitter = new genfit::KalmanFitter();
  scatteredtrack = new genfit::Track(new genfit::RKTrackRep(11), TVector3(0,0,-1700), TVector3(0,0,1));

  // Grab hits in SPS/STT
   
  if((STT) && (!USE_MC)) GatherHits_inSTTs(option);
  if(USE_MC)              CreateHits_inSTTs_MC(option);
  // Fit the beamline track
  if( scatteredtrack->getNumPoints() == 0)
  {
    debug(1,"No points\n");
    return ok;  // nothing to fit
  }
  assert(scatteredtrack.checkConsistency());           // check
  fitter->processTrack(scatteredtrack);                // do the fit
  assert(scatteredtrack.checkConsistency());           // check
  //std::cout << "Number of points: "<<scatteredtrack->getNumPoints() << std::endl;
  if(scatteredtrack->getNumPoints()==4)
    count4++;
  if(scatteredtrack->getNumPoints()==5)
    count5++;
  if(!beamlinetrack->checkConsistency())
    printf("Track not consistent\n");

  debug(1,"Fit chi^2: %g\n",scatteredtrack->getFitStatus()->getChi2());
  // Plot results
  PlotSTTs("STT");

  //Set event display:
  if(USE_DISPLAY) display->addTrack(scatteredtrack);

  return ok;
}

// ==============================================================================================================
// ------------------------------------------- Real Data Subroutines --------------------------------------------
// ==============================================================================================================





void Tracker::GatherHits_inSTTs(std::string option)
{
  int hitId=0; // hit ID
  // Create a hit in three planes of STTs
  if (STT->hits.size()<3) return;
  for (auto hit:STT->hits)
    {   
        auto wire=STT_getStrawPos(hit.id,STTangle);
        int side, plane, straw_in_plane;
        STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
        printf("Which side: %i plane: %i straw: %i distance: %f\n",side, plane, straw_in_plane,hit.dist);
        printf("Wire start %g %g %g Wire end %g %g %g\n\n",wire.first[0],wire.first[1],wire.first[2],wire.second[0],wire.second[1],wire.second[2]);

        auto measurement = new genfit::WireMeasurementNew(hit.dist,hit.disterr,wire.first,wire.second,side, ++ hitId, NULL);
        measurement->setLeftRightResolution(0);//Sets a preferred side to wires, >0 is right side, <0 is left side.
        measurement->setMaxDistance(5);
        scatteredtrack->insertPoint(new genfit::TrackPoint(measurement, scatteredtrack));
    }
  //debug(0,"Track with %i hits\n",hitId);
}

void Tracker::GatherHits_inGEMs(std::string option)
{
  TVectorD hitCoords(2);
  const int detId(0); // detector ID
  int planeId(0); // detector plane ID
  int hitId(0); // hit ID

  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
  double detectorResolution(0.1); // resolution of GEMs detectors
  hitCov *= detectorResolution*detectorResolution;

  // Loop over planes
  for (int GEMplane = 0; GEMplane < 3; GEMplane++)
  {
    for (auto hit:GEMs->hits)
    {
	    if ((hit.GEMid - 3) == GEMplane)
      {
        //currently this changes coordinate system from native GEMs to a more reasonable world coordinate system
        //however this has "xl" being the vertical axis and "yl" being horizontal
	      hitCoords[0] = -(hit.xl)*0.4 + 50; // Turn strip into mm positions (0.4mm per strip) and centers the GEM coord system
	      hitCoords[1] = -(hit.yl)*0.4 + 50; // Turn strip into mm positions (0.4mm per strip)
	      double z = (GEMpositions[GEMplane]);

	      genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
	      measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,z), TVector3(1,0,0), TVector3(0,1,0))), planeId++);
	      beamlinetrack->insertPoint(new genfit::TrackPoint(measurement, beamlinetrack));

	      if(option == "GEMs_HighAmplitudeHits") break;
      }
    }
  }
}

// ==============================================================================================================
// ----------------------------------------- Simulated Data Subroutines -----------------------------------------
// ==============================================================================================================

void Tracker::CreateHits_inGEMs_MC(std::string option)
{
  TVectorD hitCoords(2);
  const int detId(0); // detector ID
  int planeId(0); // detector plane ID
  int hitId(0); // hit ID

  double detectorResolution(0.1); // resolution of GEMs detectors

  TMatrixDSym hitCov(2);
  hitCov.UnitMatrix();
  hitCov *= detectorResolution*detectorResolution;

  // Loop over planes
  for (int i = 0; i < nGEMS; i++)
    {
      // Create some hits in gaussian smear at 10*i
      hitCoords[0] = gRandom->Gaus(10*i,1); // Turn strip into mm positions (0.4mm per strip)
      hitCoords[1] = gRandom->Gaus(10*i,1); // Turn strip into mm positions (0.4mm per strip)

      double z = (GEMpositions[i]);

      genfit::PlanarMeasurement* measurement = new genfit::PlanarMeasurement(hitCoords, hitCov, detId, ++hitId, NULL);
      measurement->setPlane(genfit::SharedPlanePtr(new genfit::DetPlane(TVector3(0,0,z), TVector3(1,0,0), TVector3(0,1,0))), planeId++);
      beamlinetrack->insertPoint(new genfit::TrackPoint(measurement, beamlinetrack));

    }
}

void Tracker::CreateHits_inSTTs_MC(std::string option)
{
  const int detId(0); // detector ID
  int hitId(0); // hit ID

  double detectorResolution(0.1); // resolution of STT detectors

  double distance_between_planes = 1.0;
  double tube_thickness = 1.0;

  // First wire
  int wire_number = 0;
  int plane_number = 0;

  // Create a hit in three planes of STTs
  for (int plane = 0; plane < 5; plane++)
    {
      TVector3 wirestart, wirestop;
      // start of wire
      wirestart[0] = -125;
      wirestart[1] = 0;
      wirestart[2] = 10+distance_between_planes*plane;

      // end of wire
      wirestop[0] = 125;
      wirestop[1] = 0;
      wirestop[2] = 10+distance_between_planes*plane;

      auto measurement = new genfit::WireMeasurementNew(gRandom->Uniform(0,0.45),detectorResolution,wirestart,wirestop,detId, ++ hitId, NULL);
      measurement->setLeftRightResolution(0); // let DAF decide hit
      measurement->setMaxDistance(tube_thickness/2.0); 

      scatteredtrack->insertPoint(new genfit::TrackPoint(measurement, scatteredtrack));

      wire_number+=5;
    }
}

// ==============================================================================================================
// ---------------------------------------------Plotting Subroutines --------------------------------------------
// ==============================================================================================================

void Tracker::PlotGEMs(std::string option)
{
  if (beamlinetrack->getNumPoints() != nGEMS) return;

  cd(option.c_str());

  // Plot GEM histograms
  for (auto i = 0; i < beamlinetrack->getNumPoints(); i++)
    {
      TVector3 track_at_plane = beamlinetrack->getFittedState(i).getPos();

      double hit_x = beamlinetrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[0];
      double hit_y = beamlinetrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[1];

      // Grab the correct directory
      std::string dir = "GEM_";
      for (int p = 0; p < nGEMS; p++)
        if (track_at_plane.Z() == GEMpositions[p]) dir += std::to_string(p);
      cd(dir.c_str());

      H1(hit_x, "hits_X", "hits_X;hits in X [mm]", 110,-10,100);
      H1(hit_y, "hits_Y", "hits_Y;hits in Y [mm]", 110,-10,100);
      H2(hit_x, hit_y, "hits_XY", "hits_XY;hits in X [mm];hits in Y [mm]", 110,-10,100, 110,-10,100);
      H2(hit_x, hit_y, "../All_hits_XY", "Allhits_XY;hits in X [mm];hits in Y [mm]", 110,-10,100, 110,-10,100);

      H1(track_at_plane.X(), "fitted_track_X", "fitted_track_X;track position in X [mm]", 110,-10,100);
      H1(track_at_plane.Y(), "fitted_track_Y", "fitted_track_Y;track position in Y [mm]", 110,-10,100);
      H2(track_at_plane.X(), track_at_plane.Y(), "fitted_track_XY", "fitted_track_XY;track position in X [mm];track position in Y [mm]", 110,-10,100, 110,-10,100);
      H2(track_at_plane.X(), track_at_plane.Y(), "../All_fitted_track_XY", "All_fitted_track_XY;track position in X [mm];track position in Y [mm]", 110,-10,100, 110,-10,100);

      // Residuals
      double res_x = hit_x - track_at_plane.X();
      double res_y = hit_y - track_at_plane.Y();
      H1(res_x, "residual_X", "residual_X;residual in X [mm]", 200,-5,5);
      H1(res_y, "residual_Y", "residual_Y;residual in Y [mm]", 200,-5,5);
      H2(res_x,res_y, "residual_XY", "residual_XY;residual in X [mm];residual in Y [mm]", 200,-5,5, 200,-5,5);
      H2(res_x,res_y, "../All_residual_XY", "All_residual_XY;residual in X [mm];residual in Y [mm]", 200,-5,5, 200,-5,5);

      cd("..");
    }

  cd("..");

}


void Tracker::PlotSTTs(std::string option)
{
  //    if (scatteredtrack->getNumPoints() != 3) return ok;

  cd(option.c_str());
  TVector3 track_at_tube;
  int i=0;
  double chi2 =0;
  // Plot GEM histograms
  //for (auto i = 0; i < scatteredtrack->getNumPoints(); i++)
    for(auto hit:STT->hits)
    {

      int side, plane, straw_in_plane;
      STT_internal_to_logic(hit.id,&side,&plane,&straw_in_plane);
      //track_at_tube = scatteredtrack->getFittedState(i).getPos();

      try{
        track_at_tube = scatteredtrack->getFittedState(i).getPos();
      }
      catch(const std::exception& e){
        std::cout << "Exception: " << e.what() << std::endl;
        i++;
        continue;
      }

      //        genfit::StateOnPlane reference(scatteredtrack->getCardinalRep());
      //        genfit::SharedPlanePtr firstPlane(scatteredtrack->getPointWithMeasurement(i)->getRawMeasurement(0)->constructPlane(reference));
      //        double residual = reference.getRep()
      //                //extrapolateToPlane(&firstPlane);
      ////        double residual = scatteredtrack->getFittedState(i).extrapolateToMeasurement(scatteredtrack->getPointWithMeasurement()->getRawMeasurement(),false,false);
      //     //   std::cout << "!================ residual : " << residual << std::endl;

      //        TVectorD track_6D = scatteredtrack->getFittedState(i).get6DState();

      //        for (auto var = 0; var < 6; var++)
      //        {
      //            std::string s = "state_" + std::to_string(var);
      //            H1(track_6D[i], s.c_str(), s.c_str(), 1000, -1, -1);
      //        }
      /*
      double hit_drift = scatteredtrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[6];

      double wire_x = (scatteredtrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[3] -
		       scatteredtrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[0])/2.0;

      double wire_y = (scatteredtrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[4] -
		       scatteredtrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[1])/2.0;

      double wire_z = (scatteredtrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[5] -
		       scatteredtrack->getPoint(i)->getRawMeasurement()->getRawHitCoords()[2])/2.0;
      H1(wire_x, "wire_x", "wire_x", 1000,-5,5);
      H1(wire_y, "wire_y", "wire_y", 1000,-5,5);
      H1(wire_z, "wire_z", "wire_z", 1000,0,25);
      
      H1(hit_drift,"hit_drift", "drift;drift distance X [mm]", 1000,-1,1); 
      */
      H1(track_at_tube.X(), "fitted_track_X", "fitted_track_X;track position in X [mm]", 1000,-1,-1);
      H1(track_at_tube.Y(), "fitted_track_Y", "fitted_track_Y;track position in Y [mm]", 1000,-1,-1);
      H1(track_at_tube.Z(), "fitted_track_Z", "fitted_track_Z;track position in Z [mm]", 1000,-1,-1);
      //Attempting to calculate residuals
      auto wire=STT_getStrawPos(hit.id,STTangle);

      auto residual = sqrt(pow(wire.second.Y()-track_at_tube.Y(),2)+pow(wire.second.Z()-track_at_tube.Z(),2)) - fabs(hit.dist);
      /*
      std::cout << "Plane: " << plane << std::endl;
      std::cout << "hit.dist: " << hit.dist << std::endl;
      std::cout << "track X Y Z: " << track_at_tube.X() << " " << track_at_tube.Y()<< " " <<track_at_tube.Z() << std::endl;
      std::cout << "Straw Pos X Y Z: " << wire.second.X() << " " << wire.second.Y() << " " << wire.second.Z() << std::endl;
      std::cout << "Residual in Z = " << residual << " mm" << std::endl;
      */
      if(residual < -0.5)
      {
       H1(residual,"Biased Residual","Residual < -1 mm;residual (mm);counts",1000,-50,50);
       //H1(plane,"Biased Residual Plane","Biased Residual Plane;Plane;counts",5,-0.5,4.5);
       H1(hit.dist,"Biased Residual Distance","Biased Residual Distance;hit distance (mm);counts",100,-1,-1);
      }
      H1(residual,"Hit Residual","Hit Residual;residual (mm);counts",500,-10,10);
      chi2 = chi2 + scatteredtrack->getFitStatus()->getChi2();
      H1(hit.dist,"Hit Distance","Hit Distance;distance (mm);counts",1000,-1,-1);
      i++;//iterate through hits
    
    }

    H1(chi2,"Fitted Chi2","Chi2",2000,-10,100);

  cd("..");

}


Long_t Tracker::finalize()
{
  std::cout << "Ratio of 4 to 5 hits in a track: " << count4/count5 << " tracks with 4 hits: " << count4 << " tracks with 5 hits: " << count5  << std::endl;
  return ok;
};

extern "C"{
  Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
  {
    return (Plugin *) new Tracker(in,out,inf_,outf_,p);
  }
}


ClassImp(Tracker);

