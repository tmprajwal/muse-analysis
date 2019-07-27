#ifndef __HOUGH__
#define __HOUGH__

#include "TObject.h"
#include "Plugin.h"
#include "TCanvas.h"
#include "TTree.h"
#include <iostream>
#include <map>
#include "TVector2.h"
#include "StrawTubehittree.h"
#include "Scinthittree.h"
#include "SPStree.h"
#include "StrawTubetree.h"
#include "Trackhittree.h"
#include "TH2.h"
#include "TEveManager.h"
#include "TEveGeoNode.h"
#include "TEveScene.h"
#include "TF1.h"
#include "TGeoManager.h"
#include "TGeoNode.h"
#include "TGeoMatrix.h"
#include "TRandom3.h"
#include "TMath.h"

#include "TSystem.h"
#include "TGFrame.h"
#include "TPad.h"
#include "TGraph.h"
#include "TLine.h"
#include "TList.h"
#include "TMarker.h"
#include "TGLabel.h"
#include "TG3DLine.h"
#include "TGMenu.h"
#include "TRootEmbeddedCanvas.h"
#include "TCanvas.h"
#include "TGDockableFrame.h"
#include "TGeoManager.h"
#include "TEveManager.h"
#include "TGedFrame.h"
#include "TGLObject.h"
#include "TGLViewer.h"
#include "TGLSAViewer.h"
#include "TGLEmbeddedViewer.h"
#include "TRootEmbeddedCanvas.h"
#include "TGLViewerBase.h"
#include "TGLCamera.h"
#include "TGLScenePad.h"
#include "TEveGeoNode.h"
#include "TEveViewer.h"
#include "TEveScene.h"
#include "TEvePointSet.h"
#include "TEveLine.h"
#include "TEveElement.h"
#include "TGButton.h"
#include "TExec.h"
#include "TGFrame.h"
#include "TGTextEntry.h"
#include "TGTab.h"
#include "TGIcon.h"
#include "TASImage.h"
#include "TGWindow.h"
#include "TLatex.h"
#include "TError.h"
#include "TVirtualGeoTrack.h"
#include "TEveTrack.h"
#include "TEveTrackPropagator.h"
#include "TEveVector.h"
#include "TPolyLine3D.h"

#include <stdio.h>
#include <math.h>
#include <fstream>
//#include "~/muse/external/cminpack/cminpack.h"
#include "../../../../../../external/cminpack/cminpack.h"
#include <iterator>

class StrawTubeHits;
class TrackHits;
class TGCompositeFrame;
class TRootEmbeddedCanvas;
class TPad;
class TGraph;
class TGLEmbeddedViewer;
class TGLViewer;
class TGLSAViewer;
class TGLScenePad;
class TGeoManager;
class TEveGeoTopNode;
class TGeoNode;
class TGeoVolume;
class TGTextButton;
class TGSplitButton;
class TGPopupMenu;
class TGWindow;
class TGTextEntry;
class TGIcon;
class TASImage;
class TGLabel;
class TGVerticalFrame;
class TString;
class TCanvas;
class TGMainFrame;
class TEveLine;
class TEvePointSet;

//class TH2D;
class hough:public Plugin
{
 private:
  StrawTubeHits *STT;
  TrackHits *Tracks;
   ScintHits *SPSHits;


  TH2D *houghplotX,*houghplotY,*circX,*circY,*min,*circ,*straws,*lineX,*lineY,*bothX,*bothY,*lineprime,*chiVangle,*chiVangleprime;
  TH1D *hitresidualX,*hitresidualnormX,*hitresidualnormY,*hitresidualY,*strawresidual,*hitchisq,*strawchisq,*AngleX,*AngleY,*FinalAnglex,*FinalAngley;

  const char* geofile;
  TGCompositeFrame *tab;



 public:
  hough(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~hough();

  double angle;
  // add funtions with return value Long_t here:
    MRTRunInfo *theRunInfo;
  void GenerateDet();
  Long_t findnodes();
  void printallnodes(TGeoNode *node,int level);
  Long_t startup();
  Long_t process();
  Long_t finalize();
  void rotAngle(double x);
  int offsetx[2] = {0,0};
  int offsety[2] = {0,0};
  void getAngle(int *numhitsX,int *numhitsY, int maxangles, std::vector<double> &xangle,std::vector<double> &yangle,std::vector<double> &xrho, std::vector<double> &yrho);
  TVector2 getStrawPos(int id);
  TVector3 getLocalPos(int id, int offsetx[2], int offsety[2]);

  double getResidual(int id,double dist, double slope, double intercept);

  int getMinLocation(std::vector<double> input);

    TGeoNode* STTL160[10][55];
  TGeoNode* STTL190[10][89];
  TGeoNode* STTR160[10][55];
  TGeoNode* STTR190[10][89];

  TGCompositeFrame *viscotab;

  //static void mychi2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

  virtual Long_t cmdline(char * cmd);

  ClassDef(hough,1);
    };

#endif
