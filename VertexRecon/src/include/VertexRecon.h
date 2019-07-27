#ifndef __VERTEXRECON__
#define __VERTEXRECON__

#include "TObject.h"
#include "Plugin.h"
#include "TVector3.h"
#include "TTree.h"
#include <iostream>
#include <Scinthittree.h>
#include "teletracktree.h"
#include "Trackhittree.h"

class VertexRecon:public Plugin
{
 private:
TeleTracks      * GEM_Tracks;
TrackHits 		* STT_Tracks;
 ScintHits *Hits;


 public:
  VertexRecon(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
  virtual ~VertexRecon();
  // add funtions with return value Long_t here:
    MRTRunInfo *theRunInfo;

  Long_t startup();
  Long_t process();
  Long_t finalize();
     void id_to_info(int id, int *plane, int *bar, int *side);

  void find_vertex(TVector3 P0, TVector3 u, TVector3 Q0, TVector3 v, TVector3 &vertex, double &doca, double &theta);
  virtual Long_t cmdline(char * cmd);

  ClassDef(VertexRecon,1);
    };

#endif
