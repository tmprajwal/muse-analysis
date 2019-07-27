#include <scifi_analysis.h>

#include <iostream>
#include <cmath>
#include "SciFitree.h" //needed for the raw tree
#include "SciFiOutput.h"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
scifi_analysis::scifi_analysis(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
    rawscifi=NULL;
    scifiout=new SciFiOutput();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
scifi_analysis::~scifi_analysis(){delete scifiout;};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Long_t scifi_analysis::startup()
{
    getBranchObject("SciFi"         ,(TObject **) & rawscifi); // rawscifi has now a pointer to this class
    getBranchObject("teletracks"    ,(TObject **) & GEM_Tracks); // rawscifi has now a pointer to this class
    getBranchObject("mappedchannels",(TObject **) & mapped);
    getBranchObject("V792"          ,(TObject **) & v792tree);
    if (!rawscifi || !GEM_Tracks || !mapped || !v792tree){
        debug(0,"Could not find the SciFi / GEM / testbeam / v792 branch in the raw tree. Bailing out\n");
        exit(-1);
    }
    // output
    makeBranch("scifi",(TObject **) &scifiout);
    return Plugin::ok;
};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Long_t scifi_analysis::process()
{ // called every event
    
    float eCut[2] = {5,8} , muCut[2] = {1,3} , piCut[2] = {10,14};
    bool PlotAllHistograms  = true;
    const int Nplane        = 2     , Nfiber        = 16;
    const int FirstFiber    = 13    , LastFiber     = 28;
    const int Nbins         = 210   , tMin          = -150  , tMax          = -90;
    const int DiffTimeMin   = -10   , DiffTimeMax   = 10;
    const int ADCMin        = 5000  , ADCMax        = 4500  , NbinsADC      = 2501;
    double reftime = 0;
    std::vector<Int_t>      U       , V ;
    
    if (rawscifi->reference_time.size()==0)
        return Plugin::ok;
    else
        reftime = rawscifi -> reference_time[rawscifi->reference_time.size()-1].time;
    debug(2, "-------------\nFound reference time %g\n" , reftime);
    
    
    
//    auto sauhits            = testbeamdata->tdc_trb.equal_range(testbeamtime::sau);
//    if (sauhits.first == sauhits.second)    return Plugin::ok;    // only run if there is a hit in SAU
//    testbeamdata -> tdc_trb.end();
//    for (auto hit=sauhits.first ; hit!=sauhits.second ; hit++)
//        debug(0,"hit->second.time = %g\n",hit->second.time);
    

//    double SAU_time         = sauhits.first->second.time - reftime  ,   SAU_wrapped_time = WrapTime( SAU_time );
//    auto sadhits            = testbeamdata->tdc_trb.equal_range(testbeamtime::sad);
//    double SAD_time         = sadhits.first->second.time - reftime  ,   SAD_wrapped_time = WrapTime( SAD_time );
//    double SAavg_time       = ( SAU_time + SAD_time )/2.            ,   SAavg_wr_time    = WrapTime( SAavg_time );
//    if (PlotAllHistograms)  H1( SAU_time ,"RefenceTimes/SA/SAU","time [ns]",Nbins,tMin,tMax);
//    if (PlotAllHistograms)  H1( SAD_time ,"RefenceTimes/SA/SAD","time [ns]",Nbins,tMin,tMax);
//    if (PlotAllHistograms)  H1( SAavg_time ,"RefenceTimes/SA/SA-avg","time [ns]",Nbins,tMin,tMax);
//    debug(1, "SAU = %.6f ns , SAD = %.6f ns, SA avg. = %.6f ns\n" ,SAU_wrapped_time,SAD_wrapped_time,SAavg_wr_time  );
//    H1( SAavg_wr_time ,"RefenceTimes/SA-avg-wrapped-time","time [ns]",Nbins,0,21);

    
    auto sauhits_adc        = mapped->adc_v792.equal_range("sau");
    double SAU_adc          = sauhits_adc.first->second;
    auto sadhits_adc        = mapped->adc_v792.equal_range("sad");
    double SAD_adc          = sadhits_adc.first->second;
    double SAavg_adc        = ( SAU_adc + SAD_adc )/2.;
    if (PlotAllHistograms)  {
        H1( SAU_adc ,"ADCs/SA/SAU","time [ns]",NbinsADC,ADCMin,ADCMax);
        H1( SAD_adc ,"ADCs/SA/SAD","time [ns]",NbinsADC,ADCMin,ADCMax);
        H1( SAavg_adc ,"ADCs/SA/SA-avg","time [ns]",NbinsADC,ADCMin,ADCMax);
        //    H2( SAavg_wr_time , SAavg_adc ,"RefenceTimes/SAavgTDC-vs-SAavgADC" ,"SA avg. TDC [ns] vs. SA avg. ADC ",Nbins,0,21,NbinsADC,ADCMin,ADCMax);
    }
    
    
    // (-) PID using SA-scintillator adc vs TOF
    // ------------------------------------------------------------------------------------------------------------------------------------
//    bool isElectron = (eCut[0]  < SAavg_wr_time && SAavg_wr_time < eCut[1]);
//    bool isMuon     = (muCut[0] < SAavg_wr_time && SAavg_wr_time < muCut[1]);
//    bool isPion     = (piCut[0] < SAavg_wr_time && SAavg_wr_time < piCut[1]);
//    debug(1, "isElectron=%d , isMuon=%d , isPion=%d\n" ,isElectron, isMuon , isPion  );
    char * p_name = "unidentified";
//    if (isElectron)     p_name = "electron";
//    else if (isMuon)    p_name = "muon";
//    else if (isPion)    p_name = "pion";
    
    
    
    
    
    
    
    
    
    
    
    
    
//    // (a) Go through all hits and just draw them
//    // ------------------------------------------------------------------------------------------------------------------------------------
//    if (PlotAllHistograms){
//        for ( int plane = 0 ; plane < Nplane ; plane++ ){
//            for (int side = 0 ; side < 2 ; side++){
//                for (auto channel_time_pair:rawscifi->hits[plane][side]){ // channel_time_pair has as .first the channel nr, and as .second the time;
//                    H1(channel_time_pair.second.time ,Form("SciFi/Plane%i/Side%i/Fiber-%i/Raw",plane,side,channel_time_pair.first) ,"raw time [ns]",Nbins,tMin,tMax);
//                    H1(channel_time_pair.second.time-reftime ,Form("SciFi/Plane%i/Side%i/Fiber-%i/RefTimeSubtracted",plane,side,channel_time_pair.first) ,"time [ns]",Nbins,tMin,tMax);
//                }
//            }
//        }
//    }
    
    
    
    
    
    
    
    
    
    
    // (b) plot (left-right)/2. differnce and (right+left)/2. average timing for all fibers for all planes:
    // ------------------------------------------------------------------------------------------------------------------------------------
    for ( int plane = 0 ; plane < Nplane ; plane++ ){
        for (int fiber = FirstFiber ; fiber <= LastFiber ; fiber++ ){
            auto righthits  = rawscifi->hits[plane][0].equal_range(fiber); // 'auto' automatically get the right type (a complicated std::pair)
            auto lefthits   = rawscifi->hits[plane][1].equal_range(fiber);

            if ((rawscifi -> hits[plane][0].count(fiber) > 0 ) && (rawscifi -> hits[plane][1].count(fiber) > 0 )) { // condition 2 sides of the fiber
                
                for (auto right=righthits.first;right!=righthits.second;right++){       // go through all the right hits
                    for (auto left=lefthits.first;left!=lefthits.second;left++) {       // go through all the left hits
                        
//                        if (PlotAllHistograms) H2( left->second.time - reftime,right->second.time - reftime
//                                                  ,Form("2DScifiCorrelations/Plane-%i/Fiber%i/NoCrossTalkSubrtraction-Left-vs-Right/",plane,fiber)
//                                                  ,"time [ns]",Nbins,tMin,tMax,Nbins,tMin,tMax);

                        if ( ( (right->second.rising == true) && (left->second.rising == true) )
                            && (fabs(left->second.time - right->second.time) < 10) ) {
                            
                            double right_time   = right -> second.time - reftime    , right_wrapped_time    = WrapTime(right_time);
                            double left_time    = left -> second.time - reftime     , left_wrapped_time     = WrapTime(left_time);
                            double avg_time     = (right_time + left_time)/2.       , avg_wrapped_time      = WrapTime(avg_time);
                            double diff_time    = (right_time - left_time)/2.       , diff_wrapped_time     = WrapTime(diff_time);
                            
                            if (PlotAllHistograms){
                                H1( left_time - right_time ,Form("ScifiDiff/Plane%i/Fiber-%i-Left-Right/",plane,fiber) ,"time [ns]",Nbins,DiffTimeMin,DiffTimeMax);
//                                H2(left_time,right_time ,Form("2DScifiCorrelations/Plane-%i/Fiber%i/Left-vs-Right/",plane,fiber)
//                                   ,"time [ns]",Nbins,tMin,tMax,Nbins,tMin,tMax);
//                                H2(left_wrapped_time,right_wrapped_time ,Form("2DScifiCorrelations/Plane-%i/Fiber%i/Wrapped-Left-vs-Right",plane,fiber)
//                                   ,"time [ns]",Nbins,0,21,Nbins,0,21);
                                H1( avg_time ,Form("ScifiAvg-Unwrapped/Plane%i/Fiber-%i/",plane,fiber) ,"time [ns]",Nbins,tMin,tMax);
                                H1( avg_wrapped_time ,Form("ScifiAvg-Wrapped/Plane%i/Fiber-%i/",plane,fiber) ,"time [ns]",Nbins,0,21);
                                H1( diff_time ,Form("ScifiDiff-Unwrapped/Plane%i/Fiber-%i/",plane,fiber) ,"time [ns]",Nbins,2*tMin,-2*tMin);
                                H1( diff_wrapped_time ,Form("ScifiDiff-Wrapped/Plane%i/Fiber-%i/",plane,fiber) ,"time [ns]",Nbins,0,21);
                                H2( avg_wrapped_time , SAavg_adc , Form("PID/SinglePlane/P%iF%ivsSADC",plane,fiber) ,"TOF [ns]",Nbins,0,21,NbinsADC,ADCMin,ADCMax);
                                //                            H1( avg_time ,Form("%s/P%iF%iAvg",p_name,plane,fiber)  ,"time [ns]",Nbins,tMin,tMax);
                                //                            H1( avg_wrapped_time ,Form("%s-Wrapped/P%iF%i",p_name,plane,fiber)  ,"time [ns]",Nbins,0,21);
                                H2( avg_wrapped_time , SAavg_adc , Form("TOFvsADC/p%if%i",plane,fiber) ,"TOF [ns]",Nbins,0,21,NbinsADC,ADCMin,ADCMax);
                            }
                            
                            // output to root tree
                            scifiout -> SciFi_Time[plane][fiber-FirstFiber]     = avg_time;
                            scifiout -> SciFi_ModTime[plane][fiber-FirstFiber]  = avg_wrapped_time;
                        }
                    }
                }
            }
        }
    }
    
    
    
    
    
    
    
    
//    
//    
//    // (c) Coincide U and V planes
//    // ------------------------------------------------------------------------------------------------------------------------------------
//    // WE NOW LOOP OVER ALL U AND ALL V PLANES TO CHECK WHICH FIBER WAS FIRED. A SIMPLER AND FASTER WAY WOULD BE TO IMPLEMENT A VECTOR OF INTEGERS IN THE COOKING STANGE THAT STATES WHICH FIBERS IN EACH PLANE WERE FIRED THIS EVENT.
//    for (int Ufiber = FirstFiber ; Ufiber <= LastFiber ; Ufiber++ ){
//        
//        auto Urighthits  = rawscifi->hits[0][0].equal_range(Ufiber);          // The auto keyword will automatically get the right type, which is some complicated std::pair
//        auto Ulefthits   = rawscifi->hits[0][1].equal_range(Ufiber);
//        size_t NHitsU[2] = {rawscifi->hits[0][0].count(Ufiber),rawscifi->hits[0][1].count(Ufiber)};
//
//        for (int Vfiber = FirstFiber ; Vfiber <= LastFiber ; Vfiber++ ){
//            
//            auto Vrighthits  = rawscifi->hits[1][0].equal_range(Vfiber);          // The auto keyword will automatically get the right type, which is some complicated std::pair
//            auto Vlefthits   = rawscifi->hits[1][1].equal_range(Vfiber);
//            size_t NHitsV[2] = {rawscifi->hits[1][0].count(Ufiber),rawscifi->hits[1][1].count(Ufiber)};
//            
//            if ((NHitsU[0]>0) && (NHitsU[1]>0) && (NHitsV[0]>0) && (NHitsV[1]>0)) { // condition 2 sides of the fiber
//                
//                for (auto Uright=Urighthits.first;Uright!=Urighthits.second;Uright++){       // go through all the Uright hits
//                    for (auto Uleft=Ulefthits.first;Uleft!=Ulefthits.second;Uleft++) {       // go through all the Uleft hits
//                        if (    ( Uright->second.rising==true && Uleft->second.rising==true ) &&  (fabs(Uleft->second.time - Uright->second.time) < 8)  ) {
//                            if (!SciFiFired(U,Ufiber))
//                                U.push_back (Ufiber);
//                            
//                            for (auto Vright=Vrighthits.first;Vright!=Vrighthits.second;Vright++){       // go through all the Vright hits
//                                for (auto Vleft=Vlefthits.first;Vleft!=Vlefthits.second;Vleft++) {       // go through all the Vleft hits
//                                    
//                                    if (    ( Vright->second.rising == true && Vleft->second.rising == true ) &&  (fabs(Vleft->second.time - Vright->second.time) < 8) ){
//                                        if (!SciFiFired(V,Vfiber))
//                                            V.push_back (Vfiber);
//                                        
//                                        
//                                        double Ur_time = Uright->second.time - reftime  , Uright_wrapped_time  = WrapTime( Ur_time );
//                                        double Ul_time = Uleft->second.time - reftime   , Uleft_wrapped_time   = WrapTime( Ul_time );
//                                        double Vr_time = Vright->second.time - reftime  , Vright_wrapped_time  = WrapTime( Vr_time );
//                                        double Vl_time = Vleft->second.time - reftime   , Vleft_wrapped_time   = WrapTime( Vl_time );
//                                        double UV_Coin = ((Ur_time + Ul_time)/2. + (Vr_time + Vl_time)/2.)/2. ,   UV_Coin_wrapped = WrapTime( UV_Coin );
//                                        double UV_Diff = ((Ur_time - Ul_time)/2. + (Vr_time - Vl_time)/2.)/2. ,   UV_Diff_wrapped = WrapTime( UV_Diff );
//                                        double U_V     = ((Ur_time + Ul_time)/2. - (Vr_time + Vl_time)/2.)/2. ,   U_V_wrapped = WrapTime( U_V );
//                                        debug(2,"hit %i: U %d & V %d fired!\n",U.size()-1,Ufiber,Vfiber);
//                                        
//                                        if (PlotAllHistograms){
//                                            //                                        H1( UV_Coin ,Form("UV/Unwrapped/U%iV%/",Ufiber,Vfiber) ,"time [ns]",Nbins,tMin,tMax);
//                                            H1( UV_Coin_wrapped ,Form("UV/Wrapped/U%iV%i/",Ufiber,Vfiber) ,"time [ns]",Nbins,0,21);
//                                            H2( (Uright_wrapped_time + Uleft_wrapped_time)/2., (Vright_wrapped_time + Vleft_wrapped_time)/2.
//                                               ,Form("UV/2D/U%iV%i/",Ufiber,Vfiber) ,"time [ns]",Nbins,0,21,Nbins,0,21);
//                                            H2( (Uright_wrapped_time - Uleft_wrapped_time)/2., (Vright_wrapped_time - Vleft_wrapped_time)/2.
//                                               ,Form("UV-Diff/2D/U%iV%i/",Ufiber,Vfiber) ,"time [ns]",Nbins,0,21,Nbins,0,21);
//                                            H1( UV_Diff_wrapped ,Form("UV-Diff/Wrapped/U%iV%i/",Ufiber,Vfiber) ,"time [ns]",Nbins,0,21);
//                                            H2( UV_Coin_wrapped , SAavg_adc , Form("PID/UV/U%iV%ivsSADC",Ufiber,Vfiber) ,"TOF [ns]",Nbins,0,21,NbinsADC,ADCMin,ADCMax);
//                                            H1( U_V_wrapped ,Form("U-V/U%iV%i/",Ufiber,Vfiber) ,"time [ns]",Nbins,0,21);
//                                            H2( U_V_wrapped , SAavg_adc , Form("PID/UV/U%iV%ivsSADC",Ufiber,Vfiber) ,"TOF [ns]",Nbins,0,21,NbinsADC,ADCMin,ADCMax);
//                                            //                                    float pCuts[6] = {0.5,3.,11,12.5,15.8,17.5};     // pi - mu - e run 3165 (210 MeV/c) U-20 V-22
//                                            //                                    float pCuts[6] = {0,4,4.1,6.5,16,20};     // pi - mu - e run 3166 (153 MeV/c) U-17 V-25
//                                            float pCuts[6] = {4,8,13,15.5,16,20};     // pi - mu - e run 3169 (115 MeV/c) U-20 V-20
//                                            int GoodUFiber = 20 , GoodVFiber = 20;
//                                            if (Ufiber==GoodUFiber && Vfiber==GoodVFiber){
//                                                if (pCuts[0]<UV_Coin_wrapped && UV_Coin_wrapped<pCuts[1])
//                                                    p_name = "pion";
//                                                else if (pCuts[2]<UV_Coin_wrapped && UV_Coin_wrapped<pCuts[3])
//                                                    p_name = "muon";
//                                                else if (pCuts[4]<UV_Coin_wrapped && UV_Coin_wrapped<pCuts[5])
//                                                    p_name = "electron";
//                                                if (p_name == "pion" || p_name == "muon" || p_name == "electron")
//                                                    //                                            H1( UV_Diff_wrapped ,Form("U%dV%d-Diff-%s",Ufiber,Vfiber,p_name) ,"",Nbins,0,21);
//                                                    //                                            H1( UV_Coin_wrapped ,Form("U%dV%d-Avg-%s",Ufiber,Vfiber,p_name) ,"",Nbins,0,21);
//                                                    H1( U_V_wrapped ,Form("U%i(avg)-V%i(avg)-%s",Ufiber,Vfiber,p_name) ,"time [ns]",Nbins,0,21);
//                                            }
//                                            //                                    H1( UV_Coin ,Form("%s/UV/U%iV%iAvg",p_name,Ufiber,Vfiber)  ,"time [ns]",Nbins,tMin,tMax);
//                                            //                                    H1( UV_Coin_wrapped ,Form("%s-Wrapped/UV/U%iV%i",p_name,Ufiber,Vfiber)  ,"time [ns]",Nbins,0,21);
//                                        }
//                                        
//                                        
//                                        //                                    // output to root tree
//                                        //                                    scifiout -> UV_Coincidence_Time[Ufiber-FirstFiber][Vfiber-FirstFiber]   = UV_Coin;
//                                        //                                    scifiout -> UV_ModCoin_Time[Ufiber-FirstFiber][Vfiber-FirstFiber]       = UV_Coin_wrapped;
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }
//    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    // (d) Project GEM tracks to SciFi
    // ------------------------------------------------------------------------------------------------------------------------------------
    Double_t DISTANCE_GEM_SCIFI = 267.5; // [mm]
    std::vector<TVector3>   SFHit   , SFPredHit; // x & y position as determined by the fibers. SciFi_z = z average of U & V
    double x0 , mx , y0 , my;
    // SciFi U ---- 15 mm ---- V -- 7.5 mm ------------------- 26 cm ------------------ GEM1 ------ 12.5 cm ------ GEM2 ------ 12.5 cm ------ GEM3 ------ 10 cm ------ Target
    
    scifiout -> Ntracks = GEM_Tracks->tracks.size();
    if(GEM_Tracks -> tracks[0].x0!=-1e4) debug(1,"Number of GEM Tracks is %d\n",GEM_Tracks->tracks.size());
    for (size_t j = 0 ; j < GEM_Tracks->tracks.size() ; j++){
        x0 = GEM_Tracks -> tracks[j].x0;
        mx = GEM_Tracks -> tracks[j].mx;
        y0 = GEM_Tracks -> tracks[j].y0;
        my = GEM_Tracks -> tracks[j].my;
        if(x0!=-1e4 && mx!=-1e4 && y0!=-1e4 && my!=-1e4) {
            debug(1,"x0 = %.2f , mx = %.2f , y0 = %.2f , my = %.2f\n",x0,mx,y0,my);
            SFPredHit.push_back( TVector3( x0 + mx *(DISTANCE_GEM_SCIFI) , y0 + my * (DISTANCE_GEM_SCIFI) , DISTANCE_GEM_SCIFI + 160 ) );
//            scifiout -> SciFiPredictedHits[j][0] = PredPos.x();
//            scifiout -> SciFiPredictedHits[j][1] = PredPos.y();
//            scifiout -> SciFiPredictedHits[j][2] = PredPos.z();
        }
    }
    for (size_t i = 0 ; i < U.size() ; i++ ){
//        debug(1,"U.at(%d) = %d\n",i,U.at(i));
        for (size_t j = 0 ; j < V.size() ; j++ ){
//            debug(1,"V.at(%d) = %d\n",j,V.at(j));
            SFHit.push_back(TVector3( -sqrt(2.)*( U.at(i) - V.at(j) ) , sqrt(2.)*( (U.at(i) - 21) + (V.at(j) - 21) ) , 427.5 ));
            //        scifiout -> SciFiActualHits[j][0] = HitPos.x();
            //        scifiout -> SciFiActualHits[j][1] = HitPos.y();
            //        scifiout -> SciFiActualHits[j][2] = HitPos.z();
        }
    }
//    scifiout -> NSciFiHits  = SciFi_HitPosition.size();

    
    debug(2,"%d Predicted hits\n" ,SFPredHit.size());
    for (size_t i = 0 ; i < SFPredHit.size() ; i++ ){
        debug(2,"(%g,%g,%g)\n",SFPredHit.at(i).x(),SFPredHit.at(i).y(),SFPredHit.at(i).z());
        for (size_t j = 0 ; j < U.size() ; j++ ){
            if (U.at(j)==17)
                H2(SFPredHit.at(i).x(),SFPredHit.at(i).y(),Form("GEM/GEM+U%d",U.at(j)),"x vs. y [mm]",101,-50,50,101,-50,50);
        }
        H2(SFPredHit.at(i).x(),SFPredHit.at(i).y(),"GEM/GEMprofile","x vs. y [mm]",101,-50,50,101,-50,50);
    }
    debug(2,"%d SciFi hits\n",SFHit.size());
    for (size_t i = 0 ; i < SFHit.size() ; i++ ){
        debug(2,"SciFi hit: (%g,%g,%g)\n",SFHit.at(i).x(),SFHit.at(i).y(),SFHit.at(i).z());
        H2(SFHit.at(i).x(),SFHit.at(i).y(),"BeamProfile/SFHit","x vs. y [mm]",101,-50,50,101,-50,50);
    }
    
//    
//    for (size_t i = 0 ; i < SFHit.size() ; i++ ){
//        for (size_t j = 0 ; j < SFPredHit.size() ; j++ ){
//            double Difference[2];
//            for (int k = 0 ; k < 2 ; k++ ){
//                Difference[k] = scifiout -> SciFiActualHits[i][k] - scifiout -> SciFiPredictedHits[j][k];
//                scifiout -> SciFiHitsDifference[i][j][k] = Difference[k];
//            }
////            debug(1,"Actual Hit %d - Pred Hit %d: Dx = %g mm, Dy = %g mm \n",scifi_hit, pred_hit ,Difference[0] , Difference[1]);
//            H2(Difference[0] , Difference[1] , "BeamProfile/SFHit-SFPredHit","x vs. y [mm]",101,-50,50,101,-50,50);
//        }
//    }
//    
    
    
    
    
    // clear vectors
    U.clear();
    V.clear();
    SFPredHit.clear();
    SFHit.clear();
    debug(1,"--------------\n");

    
    return Plugin::ok;
};





//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
double scifi_analysis::WrapTime(double raw_time){
    return fmod( raw_time + 1005, 19.75 );
};






//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
bool scifi_analysis::SciFiFired(std::vector<Int_t> PlaneHits, Int_t Fiber){
    bool res = false;
    for (size_t i = 0 ; i < PlaneHits.size() ; i++ ){
        if ( Fiber == PlaneHits.at(i) )
            res = true;
    }
    PlaneHits.erase(PlaneHits.begin(),PlaneHits.end());
    return res;
}




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
Long_t scifi_analysis::cmdline(char *cmd)
{
    //add cmdline hanling here
    return 0; // 0 = all ok
};



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
extern "C"{
    Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
    {   return (Plugin *) new scifi_analysis(in,out,inf_,outf_,p);  }

}


ClassImp(scifi_analysis);

