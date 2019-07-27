#ifndef __SCIFI_ANALYSIS__
#define __SCIFI_ANALYSIS__

#include "TObject.h"
#include "Plugin.h"

#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TVector3.h"

#include "teletracktree.h"
#include "mappedchannelstree.h"
#include <iostream>
#include <vector>
class SciFi;
class SciFiOutput;

class scifi_analysis:public Plugin
{
private:
    SciFi           * rawscifi;
    TeleTracks      * GEM_Tracks;
    SciFiOutput     * scifiout;
    mappedchannels  * mapped;
    v792            * v792tree;
    
public:
    scifi_analysis(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p);
    virtual ~scifi_analysis();
    // add funtions with return value Long_t here:
    
    Long_t startup();
    Long_t process();
    // Long_t finalize()
    
    virtual Long_t cmdline(char * cmd);

    
    double WrapTime(double);
    bool SciFiFired(std::vector<Int_t>, Int_t);

    ClassDef(scifi_analysis,1);
};


#endif
