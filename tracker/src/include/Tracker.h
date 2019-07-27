#ifndef __TRACKER__
#define __TRACKER__

#include "TObject.h"
#include "Plugin.h"
#include "TTree.h"
#include "TrackerTree.h"
#include "Hittree.h"
#include "StrawTubehittree.h"
#include "lumigemtree.h"
#include "muserawtree.h"


#include <iostream>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <cmath>


#include <TApplication.h>
#include <TCanvas.h>
#include <TDatabasePDG.h>
#include <TEveManager.h>
#include <TGeoManager.h>
#include <TH1D.h>
#include <TRandom.h>
#include <TStyle.h>
#include <TVector3.h>
#include <vector>

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include "TDatabasePDG.h"
#include <TMath.h>
#include <TString.h>

#include <memory>

#include <EventDisplay.h>
#include <AbsFinitePlane.h>
#include "AbsTrackRep.h"
#include <AbsFitterInfo.h>
#include <TGeoMaterialInterface.h>
#include <HelixTrackModel.h>
#include <MeasurementCreator.h>
#include <AbsMeasurement.h>
#include <AbsTrackRep.h>
#include <ConstField.h>
#include <DetPlane.h>
#include <Exception.h>
#include <FieldManager.h>
#include <KalmanFittedStateOnPlane.h>
#include <AbsKalmanFitter.h>
#include <KalmanFitter.h>
#include <KalmanFitterRefTrack.h>
#include <KalmanFitterInfo.h>
#include <KalmanFitStatus.h>
#include <DAF.h>
#include <GFGbl.h>
#include <MeasuredStateOnPlane.h>
#include <MeasurementOnPlane.h>
#include <PlanarMeasurement.h>
#include <RectangularFinitePlane.h>
#include <ReferenceStateOnPlane.h>
#include <SharedPlanePtr.h>
#include <SpacepointMeasurement.h>
#include <StateOnPlane.h>
#include <Tools.h>
#include <TrackCand.h>
#include <TrackCandHit.h>
#include <Track.h>
#include <TrackPoint.h>
#include <WireMeasurement.h>
#include <WireMeasurementNew.h>


#include <MaterialEffects.h>
#include <RKTools.h>
#include <RKTrackRep.h>
#include <StepLimits.h>
class TGCompositeFrame;
class Tracker:public Plugin
{
 private:
    int nGEMS = 3;
    double STTangle;
    LumiGEM *GEMs;
    StrawTubeHits *STT;

    genfit::Track* beamlinetrack;
    genfit::Track* scatteredtrack;
    genfit::Track* mastertrack;
    genfit::EventDisplay* display;

    std::vector<double> GEMpositions;

    bool USE_MC = false;
    bool USE_DISPLAY = false;
    bool MAX_EVENT = false;
    int nMAX_EVENT;

    
    Long_t process_beamlinetrack(std::string option);
    Long_t process_scatteredtrack(std::string option);
    Long_t process_track();

    // ===============================================
    // ------------- Real Data subroutines------------
    // ===============================================
    void GatherHits_inGEMs(std::string option);
    void GatherHits_inSTTs(std::string option);

    // ===============================================
    // ---------- Simulated Data subroutines----------
    // ===============================================
    void CreateHits_inGEMs_MC(std::string option);
    void CreateHits_inSTTs_MC(std::string option);

    // ===============================================
    // -------------Plotting subroutines--------------
    // ===============================================
    void PlotGEMs(std::string option);
    void PlotSTTs(std::string option);

  // recipe options
    void useMC() {USE_MC = true;}
    void showEventDisplay() {USE_DISPLAY = true;}
    void stopAtEvent(int ev) {MAX_EVENT = true; nMAX_EVENT = ev;}
    TGCompositeFrame *tab;
 public:
    Tracker(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
    virtual ~Tracker();
    double count4 = 0;
    double count5 = 0;
    
    Long_t startup();
    Long_t setup_fitter();
    Long_t process();
    Long_t finalize();
    Long_t defineHistograms() {return ok;}
    Long_t setSTTangle(double angle);
    virtual Long_t cmdline(char * cmd) {return ok;}



  ClassDef(Tracker,1);
    };

#endif
