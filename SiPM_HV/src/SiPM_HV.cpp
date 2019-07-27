#include <SiPM_HV.h>

#include<iostream>
#include<cmath>


SiPM_HV::SiPM_HV(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p):Plugin(in,out,inf_,outf_,p)
{
};

SiPM_HV::~SiPM_HV()
{
};



template <typename T>
auto findSmallestTime(T &container)
{
  auto best=container.begin();
  for (auto hit=container.begin();hit!=container.end();hit++)
  {
      if (best->second.time > hit->second.time){ //&& !(hit->trailing)
        best=hit;
      }
    }
return best;
}



//---------------------------------  Main routing ----------------------------------------//
//defining all histograms that we need:
Long_t SiPM_HV::defineHistograms()
{

  debug(0,"\t\t\tDefining histogram for BH SiPM MQDCs \n");
  for(int i=0; i<4; i++){
     //loop over all channels:
     for(int j=0; j<16; j++){
          MQDC_up.BH[i][j] = dH1(TString::Format("Beam Hodoscope/SiPM Plane_%d/Bar_%d/QDC Up",i,j),"QDC UP; QDC channel; counts",4197,-100,4096);
          MQDC_down.BH[i][j] = dH1(TString::Format("Beam Hodoscope/SiPM Plane_%d/Bar_%d/QDC Down",i,j),"QDC Down; QDC channel; counts",4197,-100,4096);

          Ped_up.BH[i][j] = dH1(TString::Format("Beam Hodoscope/SiPM Plane_%d/Bar_%d/Ped. Up",i,j),"QDC Pedestal UP; QDC channel; counts",4197,-100,4096);
          Ped_down.BH[i][j] = dH1(TString::Format("Beam Hodoscope/SiPM Plane_%d/Bar_%d/Ped. Down",i,j),"QDC Pedestal Down; QDC channel; counts",4197,-100,4096);

          Peak_up.BH[i][j] = dH1(TString::Format("Beam Hodoscope/SiPM Plane_%d/Bar_%d/Peak Up",i,j),"QDC Peak UP; QDC channel; counts",4197,-100,4096);
          Peak_down.BH[i][j] = dH1(TString::Format("Beam Hodoscope/SiPM Plane_%d/Bar_%d/Peak Down",i,j),"QDC Peak Down; QDC channel; counts",4197,-100,4096);

 //         MQDC_mean.BH[i][j] = dH1(TString::Format("Beam Hodoscope/SiPM Plane_%d/Bar_%d/G_Mean",i,j),"QDC GeoMean; QDC channel; counts",4197,-100,4096);
     }
  }

  debug(0,"\t\t\tDefining histogram BM Big Bars MQDCs \n");
  for(int i=0; i<4; i++){
     //loop over all channels:
     MQDC_up.BM_big[i] = dH1(TString::Format("Beam Monitor/Big Scintilators/Bar_%d/QDC Up",i),"QDC UP; QDC channel; counts",4197,-100,4096);
     MQDC_down.BM_big[i] = dH1(TString::Format("Beam Monitor/Big Scintilators/Bar_%d/QDC Down",i),"QDC Down; QDC channel; counts",4197,-100,4096);

     Ped_up.BM_big[i] = dH1(TString::Format("Beam Monitor/Big Scintilators/Bar_%d/Ped Up",i),"QDC Pedestal UP; QDC channel; counts",4197,-100,4096);
     Ped_down.BM_big[i] = dH1(TString::Format("Beam Monitor/Big Scintilators/Bar_%d/Ped Down",i),"QDC Pedestal Down; QDC channel; counts",4197,-100,4096);

     Peak_up.BM_big[i] = dH1(TString::Format("Beam Monitor/Big Scintilators/Bar_%d/Peak Up",i),"QDC Peak UP; QDC channel; counts",4197,-100,4096);
     Peak_down.BM_big[i] = dH1(TString::Format("Beam Monitor/Big Scintilators/Bar_%d/Peak Down",i),"QDC Peak Down; QDC channel; counts",4197,-100,4096);


//     MQDC_mean.BM_big[i] = dH1(TString::Format("Beam Monitor/Big Scintilators/Bar_%d/G_Mean",i),"QDC GeoMean; QDC channel; counts",4197,-100,4096);
  }

  debug(0,"\t\t\tDefining histogram BM SiPM MQDCs \n");
  for(int i=0; i<32; i++){
     //loop over all channels:
     MQDC_up.BM[i] = dH1(TString::Format("Beam Monitor/SiPM Plane/Bar_%d/QDC Up",i),"QDC UP; QDC channel; counts",4197,-100,4096);
     MQDC_down.BM[i] = dH1(TString::Format("Beam Monitor/SiPM Plane/Bar_%d/QDC Down",i),"QDC Down; QDC channel; counts",4197,-100,4096);

     Ped_up.BM[i] = dH1(TString::Format("Beam Monitor/SiPM Plane/Bar_%d/Ped Up",i),"QDC Pedestal UP; QDC channel; counts",4197,-100,4096);
     Ped_down.BM[i] = dH1(TString::Format("Beam Monitor/SiPM Plane/Bar_%d/Ped Down",i),"QDC Pedestal Down; QDC channel; counts",4197,-100,4096);

     Peak_up.BM[i] = dH1(TString::Format("Beam Monitor/SiPM Plane/Bar_%d/Peak Up",i),"QDC Peak UP; QDC channel; counts",4197,-100,4096);
     Peak_down.BM[i] = dH1(TString::Format("Beam Monitor/SiPM Plane/Bar_%d/Peak Down",i),"QDC Peak Down; QDC channel; counts",4197,-100,4096);

//     MQDC_mean.BM[i] = dH1(TString::Format("Beam Monitor/SiPM Plane/Bar_%d/G_Mean",i),"QDC GeoMean; QDC channel; counts",4197,-100,4096);
  }

return ok;
}



// startup() routine is running at the beginning of analysis
Long_t SiPM_HV::startup()
{
printf("\t\t\tstartup() routine has started:\n");

//opening SiPM Hodoscope tree in root file
  bhraw = NULL;
  getBranchObject("BH",(TObject **) &bhraw);
  if (!bhraw) {
      debug(0,"Could not find BH tree in file\n");
  }


//opening BM tree: => this we will do later
  bmraw = NULL;
  getBranchObject("BM",(TObject **) &bmraw);
  if (!bmraw) {
      debug(0,"Could not find BM tree in file\n");
  }

printf("\t\t\t End of startup() routine!\n");
return ok;
}






//process() routine is running for every event
Long_t SiPM_HV::process()
{
  //Plug in RF:
  if(bhraw->trb_reftime.size()==0){
    debug(1,"NO REFERENCE TIME FOR TRB SiPM\n");
    return ok;
  }
  //Trigger:
  if(bhraw->trig.size()==0){
    debug(1,"NO TRIGGER TIME FOR SIPM\n");
    return ok;
  }


/*********************** Plots for BH  histograms: **************************************/
    for(int i=0; i<4; i++){
        for(int j=0; j<16; j++){
            auto tdc_up = findSmallestTime(bhraw->plane[i][j].tdc_trb[1]);
            auto tdc_down = findSmallestTime(bhraw->plane[i][j].tdc_trb[0]);

            //filling hist for Up SiPM
            for(auto qdc_up: bhraw-> plane[i][j].adc_mqdc[1]){
                MQDC_up.BH[i][j]  -> Fill(qdc_up);

                //selectring pedistal & peak events based on TDC information:
                if(tdc_up->second.rising && tdc_down->second.rising && qdc_up >=100) {
					Peak_up.BH[i][j]-> Fill(qdc_up);
                } else{
                        Ped_up.BH[i][j]->Fill(qdc_up);
                }//end if statement

            } //end loop over up QDCs

            //filling hist for Down SiPM
            for(auto qdc_down: bhraw-> plane[i][j].adc_mqdc[0]){
                MQDC_down.BH[i][j]-> Fill(qdc_down);

                //selectring pedistal & peak events based on TDC information:
                if(tdc_up->second.rising && tdc_down->second.rising && qdc_down >=100) {
                	Peak_down.BH[i][j]-> Fill(qdc_down);
                } else{
                        Ped_down.BH[i][j]->Fill(qdc_down);
                }//end if statement

            } //end loop over down QDCs


            //filling hist for both QDC_UP * QDC_DOWN;
            for(auto qdc_up: bhraw-> plane[i][j].adc_mqdc[1]){
                for(auto qdc_down: bhraw-> plane[i][j].adc_mqdc[0]){
//                   MQDC_mean.BH[i][j]-> Fill(sqrt(qdc_up*qdc_down));
                    H2(qdc_up,qdc_down,TString::Format("Beam Hodoscope/SiPM Plane_%d/Bar_%d/ QDC UP vs QDC DOWN", i,j),"QDC UP vs QDC DOWN",500,0,500,500,0,500);
                }
            }

        } //end loop over cahnnels
    } //end loop over planes


/*********************** Plots for BM  histograms: **************************************/
/***********For BM in a new mapping plane[0] -> SiPMs; plane[1] -> Big Bars *************/

    for(int j=0; j<32; j++){
        auto tdc_up = findSmallestTime(bmraw->plane[0][j].tdc_trb[1]);
        auto tdc_down = findSmallestTime(bmraw->plane[0][j].tdc_trb[0]);

       //filling hist for Up SiPM
       for(auto qdc_up: bmraw-> plane[0][j].adc_mqdc[1]){
          MQDC_up.BM[j]  -> Fill(qdc_up);

         //selectring pedistal & peak events based on TDC information:
         if(tdc_up->second.rising && tdc_down->second.rising && qdc_up >=100){
                Peak_up.BM[j] -> Fill(qdc_up);
         } else{
                Ped_up.BM[j] -> Fill(qdc_up);
         }//end if statement

       } // end loop over up SiPMs in BM

       //filling hist for Down SiPM
       for(auto qdc_down: bmraw-> plane[0][j].adc_mqdc[0]){
          MQDC_down.BM[j]-> Fill(qdc_down);

         //selectring pedistal & peak events based on TDC information:
         if(tdc_up->second.rising && tdc_down->second.rising && qdc_down >=100){
                Peak_down.BM[j]-> Fill(qdc_down);
         } else{
                Ped_down.BM[j]->Fill(qdc_down);
         }//end if statement

       } //end loop over down SiPM in BM

/*
       //filling hist for Geometric Mean = Sqrt(QDC_UP * QDC_DOWN);
       for(auto qdc_up: bmraw-> plane[0][j].adc_mqdc[1]){
         for(auto qdc_down: bmraw-> plane[0][j].adc_mqdc[0]){
          MQDC_mean.BM[j]-> Fill(sqrt(qdc_up*qdc_down));
         }
       }
*/
   }



/*********************** Plots for BM Big Bars  histograms: *****************************/
    for(int j=0; j<4; j++){
       auto tdc_up = findSmallestTime(bmraw->plane[1][j].tdc_trb[1]);
       auto tdc_down = findSmallestTime(bmraw->plane[1][j].tdc_trb[0]);

       //filling hist for Up Big Bars
       for(auto qdc_up: bmraw-> plane[1][j].adc_mqdc[1]){
          MQDC_up.BM_big[j]  -> Fill(qdc_up);

         //selectring pedistal & peak events based on TDC information:
         if(tdc_up->second.rising && tdc_down->second.rising && qdc_up >=100){
                Peak_up.BM_big[j]-> Fill(qdc_up);
         } else{
                Ped_up.BM_big[j]->Fill(qdc_up);
         }//end if statement

       }//end loop over up Big Bars

       //filling hist for Down Big Bars
       for(auto qdc_down: bmraw-> plane[1][j].adc_mqdc[0]){
          MQDC_down.BM_big[j]-> Fill(qdc_down);

         //selectring pedistal & peak events based on TDC information:
         if(tdc_up->second.rising && tdc_down->second.rising && qdc_down >=100){
                Peak_down.BM_big[j]-> Fill(qdc_down);
         } else{
                Ped_down.BM_big[j]->Fill(qdc_down);
         }//end if statement

       } //end loop over down Big Bars

/*
       //filling hist for Geometric Mean = Sqrt(QDC_UP * QDC_DOWN);
       for(auto qdc_up: bmraw-> plane[1][j].adc_mqdc[1]){
         for(auto qdc_down: bmraw-> plane[1][j].adc_mqdc[0]){
          MQDC_mean.BM_big[j]-> Fill(sqrt(qdc_up*qdc_down));
         }
       }
*/
  }





return ok;
}


/************************************* Fitting Function Declaration ****************************************/

//QDC ped with piecewice  gausian: mean and peak height are the common parameter for both of them.
Double_t fit_Gaus_Com(Double_t *x, Double_t *par) {
  return par[0]*((x[0]<=par[1])*TMath::Gaus(x[0], par[1], par[2]) + (x[0] > par[1])*TMath::Gaus(x[0], par[1], par[3]));
}

//QDC spectra has an exponential background:
Double_t fit_Background(Double_t *x, Double_t *par) {
  return par[0]*TMath::Exp(par[1]*x[0]);
}

Double_t fit_Background_lin(Double_t *x, Double_t *par) {
  return par[0]*x[0] + par[1];
}

Double_t fit_Landau(Double_t *x, Double_t *par) { //
  return par[0]*TMath::Landau(x[0], par[1], par[2]);
}

//QDC spectra will be fitted with 2 assimetric gausians + exponential background:
Double_t fit_Global1(Double_t *x, Double_t *par) {
  return fit_Gaus_Com(x, par) + fit_Gaus_Com(x, &par[4]) + fit_Background(x, &par[8]);
}

//QDC spectra will be fitted with 1 assymetric gausian + landau + linear background:
Double_t fit_Global2(Double_t *x, Double_t *par) {
  return fit_Gaus_Com(x, par) + fit_Landau(x, &par[4]) + fit_Background_lin(x, &par[7]); //changed to landau
}

Double_t fit_Global3(Double_t *x, Double_t *par) {
  return fit_Gaus_Com(x, par) + fit_Landau(x, &par[4]) + fit_Background(x, &par[7]); //changed to landau
}

Double_t fit_Global4(Double_t *x, Double_t *par) { //peak
  return fit_Landau(x, par) + fit_Background(x, &par[3]); //changed to landau
}



// Old Fitting Functions:

//QDC pedestal function
Double_t fit_Ped(Double_t *x, Double_t *par) {
//  return par[0]*TMath::Landau(x[0], par[1], par[2]) ;
  return par[0]*TMath::Gaus(x[0], par[1], par[2]);
}

// QDC Peak function
Double_t fit_Peak(Double_t *x, Double_t *par) {
  return par[0]*TMath::Gaus(x[0], par[1], par[2]);
}

// Sum of background and peak function
Double_t fit_Both(Double_t* x, Double_t *par) {
  return fit_Ped(x, par) + fit_Peak(x,&par[3]);
}


// This function that fits QDCs and returns the array of vaules in the following order:
// {pedestal, peak}
void SiPM_HV:: fit_QDC(double qdc_fit_range_max, double *data_output,TH1D *ped,  TH1D *peak, TH1D *qdc ){ // plane -- detector plane (0-3 for BH, 0-1 for BM); bar -- bar number, det_type = 0(BH)/1(BM)
//defining fitting functions for each bar:
  TF1 *fit_Ped = new TF1("fit_Ped",fit_Gaus_Com,0, 4000, 4);
  //setting names for fitting parameters:
  fit_Ped -> SetParNames("Const_fit", "Mean_fit", "Sigma_1","Sigma_2");
  //Setting Parameters limits for Fitting:
  double max_Ped = ped -> GetMaximum(); //When you are fitting with gausian, you should use this in "SetParameters"
  double mean_Ped = ped -> GetMean();
  double sigma_Ped = ped -> GetStdDev();

  //Setting up contstrains on fitting parameters:
  fit_Ped -> FixParameter(0, max_Ped);
  fit_Ped -> SetParLimits(1, 0.5*mean_Ped, 1.5*mean_Ped);
  fit_Ped -> SetParLimits(2, 0.2*sigma_Ped, 4*sigma_Ped);
  fit_Ped -> SetParLimits(3, 0.2*sigma_Ped, 4*sigma_Ped);
  fit_Ped -> SetLineColor(9); //sets color to violet

  //fitting histograms:
  gStyle->SetOptFit(1);
  ped -> Fit("fit_Ped", "Q", "", mean_Ped - 3*sigma_Ped, mean_Ped + 3*sigma_Ped);  //Fitting only in range from 1 to 100

  /****************  Now we can fit the QDC spectra with fixed pedestal parameters obtained from fitting pedestal separately ***********************/
  /*defining QDC fitting functions for each bar:
      [0]-[3] -> pedestal parameters
      [4]-[6] -> peak parameters
      [7]-[8] -> background
  */
  TF1 *fit_QDC = new TF1("fit_QDC", fit_Global2, 0, 4000, 9); //

  //Setting Parameters limits for Fitting:
  //Up:
  double max_QDC = peak -> GetMaximum();
  double mean_QDC = peak -> GetMean();
  double sigma_QDC = peak -> GetStdDev();

  fit_QDC -> SetParLimits(0, 0.95*fit_Ped->GetParameter(0), fit_Ped->GetParameter(0));
  fit_QDC -> FixParameter(1, fit_Ped->GetParameter(1));
  fit_QDC -> FixParameter(2, fit_Ped->GetParameter(2));
  fit_QDC -> FixParameter(3, fit_Ped->GetParameter(3));

  fit_QDC -> SetParLimits(4, 0.9*5.8*max_QDC, 1.1* 5.8*max_QDC);
  fit_QDC -> SetParLimits(5, 0.8*.9*mean_QDC, 1.2 * .9*mean_QDC); //maybe do a range
  fit_QDC -> SetParLimits(6, 0.4*.4*sigma_QDC, 1.6*.4*sigma_QDC);
  fit_QDC -> SetParLimits(7, -100, 0); // the slope should be a negative
  fit_QDC-> SetParLimits(8, 0, 2*max_QDC); // the value at of line at x=0;
  fit_QDC -> SetLineColor(2); //sets color to red

  //fitting histograms:
  gStyle->SetOptFit(1);
  qdc -> Fit("fit_QDC", "Q", "", 5, qdc_fit_range_max);

  data_output[0] = fit_QDC->GetParameter(1); // pedestal position
  data_output[1] = fit_QDC->GetParameter(5); // peak position
  //printf("data[1] = %lf; data[0] = %lf \n", data_output[1], data_output[0]);

}

void SiPM_HV:: fill_gain_tree(TTree *tree, char const * name, int plane, int bar, int up_down, double pedestal, double peak){
    TString gain_header ("plane/I:bar/I:up_down/I:pedestal/F:peak/F:qdc/F:gain/F:shift_dV/F");
    tree_struct_t tree_struct;
    TBranch* branch = tree->GetBranch(name); //check if branch already exist
    if(!branch){
      printf("Gain branch was not defined\n");
      branch = tree->Branch(name, &tree_struct.plane, gain_header); //if branch doesn't exist, create it
    }
    else{
      printf("Gain branch was defined\n");
      branch->SetAddress(&tree_struct.plane); //if exist, reset the adress to the structure (just in case)
    }
    tree_struct.plane=plane;
    tree_struct.bar = bar;
    tree_struct.up_or_down = up_down;
    tree_struct.pedestal = pedestal;
    tree_struct.peak = peak;
    tree_struct.qdc = peak - pedestal;
    tree_struct.gain = 1.0/(peak - pedestal);
    tree_struct.shift_dV = 0; //not implemented for now!
    tree->Fill();
}


void SiPM_HV:: fill_gain_hist(int plane, double* data_array, int array_size, char  const* D_name, char const *UpDown){
  TH1D* hist = dH1(TString::Format("Gain_%s/Plane_%d_%s", D_name, plane, UpDown),TString::Format("Gain Distribution in %s Plane: %d (%s);Bar_Number",D_name, plane, UpDown),array_size,0, array_size);
       for (int j=0; j<array_size; j++){
          hist->SetBinContent(j+1,data_array[j]);
          hist->GetXaxis()->SetBinLabel(j+1,TString::Format("%d",j));
       }
       hist->Draw();
}



// finalize() routine is running at the end of analysis
Long_t SiPM_HV::finalize()
{
    double PedPeak_BH_Up[4][16], PedPeak_BH_Down[4][16]; //array to store ped peak position differences
    printf("\t\t\tfinalize() routine has started:\n");
/**********************************  Fitting  QDC spectra for BH   ****************************************/
// ----------------------------  BH detector ------------------------------//
    TTree *tree_BH = new TTree("BH_Gain","gain_table");

    for(int i=0; i<4; i++){
//loop over all bars:
        for(int j=0; j<16; j++){
          double qdc_pos_up[2];
          double qdc_pos_down[2];
//Fitting data:
          fit_QDC(1000, qdc_pos_up, Ped_up.BH[i][j], Peak_up.BH[i][j], MQDC_up.BH[i][j]);
          fit_QDC(1000, qdc_pos_down, Ped_down.BH[i][j], Peak_down.BH[i][j], MQDC_down.BH[i][j]);

          PedPeak_BH_Up[i][j] = 1.0/(qdc_pos_up[1] - qdc_pos_up[0]);
          PedPeak_BH_Down[i][j] = 1.0/(qdc_pos_down[1] - qdc_pos_down[0]);
//Filling the tree with gain table:
          fill_gain_tree(tree_BH,  "BH_Gain_Tree", i, j, 1, qdc_pos_up[0], qdc_pos_up[1]);
          fill_gain_tree(tree_BH,  "BH_Gain_Tree", i, j, 0, qdc_pos_down[0], qdc_pos_down[1]);
        } //end for loop over bars
//Plotting histograms of gain
        fill_gain_hist(i, PedPeak_BH_Up[i], 16, "BH", "Up");
        fill_gain_hist(i, PedPeak_BH_Down[i], 16, "BH", "Down");
    }//end for loop over planes

/**********************************  Fitting  QDC spectra for BM   ****************************************/
//-----------------------------  BM detectors: ------------------------------------//
    double PedPeak_BM_Up[32], PedPeak_BM_Down[32]; //array to store ped peak position differences
    TTree *tree_BM = new TTree("BM_Gain","gain_table");
    for(int j=0; j<32; j++){
      double qdc_pos_up[2];
      double qdc_pos_down[2];
//Fitting data:
      fit_QDC(1500, qdc_pos_up, Ped_up.BM[j], Peak_up.BM[j], MQDC_up.BM[j]);
      fit_QDC(1500, qdc_pos_down, Ped_down.BM[j], Peak_down.BM[j], MQDC_down.BM[j]);

      PedPeak_BM_Up[j] = 1.0/(qdc_pos_up[1] - qdc_pos_up[0]);
      PedPeak_BM_Down[j] = 1.0/(qdc_pos_down[1] - qdc_pos_down[0]);
//Filling the tree with gain table:
      fill_gain_tree(tree_BM,  "BM_Gain_Tree", 0, j, 1, qdc_pos_up[0], qdc_pos_up[1]);
      fill_gain_tree(tree_BM,  "BM_Gain_Tree", 0, j, 0, qdc_pos_down[0], qdc_pos_down[1]);
    } //end for loop over bars
//Plotting histograms of gain
    fill_gain_hist(0, PedPeak_BM_Up, 32, "BM", "Up");
    fill_gain_hist(0, PedPeak_BM_Down, 32, "BM", "Down");
// --------------- end fitting QDC BM --------------------------//

// --------------fitting for QDC BM Big Bars---------------------//
    double PedPeak_BM_Big_Up[4], PedPeak_BM_Big_Down[4]; //array to store ped peak position differences
    TTree *tree_BM_Big = new TTree("BM_Big_Gain","gain_table");
    for(int j=0; j<4; j++){
      double qdc_pos_up[2];
      double qdc_pos_down[2];
//Fitting data:
      fit_QDC(3500, qdc_pos_up, Ped_up.BM_big[j], Peak_up.BM_big[j], MQDC_up.BM_big[j]);
      fit_QDC(3500, qdc_pos_down, Ped_down.BM_big[j], Peak_down.BM_big[j], MQDC_down.BM_big[j]);

      PedPeak_BM_Big_Up[j] = 1.0/(qdc_pos_up[1] - qdc_pos_up[0]);
      PedPeak_BM_Big_Down[j] = 1.0/(qdc_pos_down[1] - qdc_pos_down[0]);
//Filling the tree with gain table:
      fill_gain_tree(tree_BM_Big,  "BM_Big_Gain_Tree", 1, j, 1, qdc_pos_up[0], qdc_pos_up[1]);
      fill_gain_tree(tree_BM_Big,  "BM_Big_Gain_Tree", 1, j, 0, qdc_pos_down[0], qdc_pos_down[1]);
    } //end for loop over bars
//Plotting histograms of gain
    fill_gain_hist(1, PedPeak_BM_Big_Up, 4, "BM_Big", "Up");
    fill_gain_hist(1, PedPeak_BM_Big_Down, 4, "BM_Big", "Down");
/**********End of BM_Big fitting******************/


  printf("\t\t\tend of finalize() routine!\n");


return ok;
}




Long_t SiPM_HV::cmdline(char *cmd)
{
  //add cmdline hanling here

  return 0; // 0 = all ok
};


extern "C"{
Plugin *factory(TTree *in, TTree *out,TFile *inf_, TFile * outf_, TObject *p)
{
  return (Plugin *) new SiPM_HV(in,out,inf_,outf_,p);
}
}


ClassImp(SiPM_HV);
