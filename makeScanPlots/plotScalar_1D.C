#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TColor.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TPaletteAxis.h"
#include <iostream>
#include "../CMS_lumi.h"

int mmed(double mh, int code){
    if (code == 800) return ((int)(mh-80000000000))/10000; 
    if (code == 801) return ((int)(mh-80100000000))/10000; 
    if (code == 805) return ((int)(mh-80500000000))/10000; 
    if (code == 806) return ((int)(mh-80600000000))/10000; 
    return -1;
}

int mdm(double mh, int code){
    if (code == 800) return (mh-80000000000)  - ( ((Int_t)(mh-80000000000))/10000 )*10000;
    if (code == 801) return (mh-80100000000)  - ( ((Int_t)(mh-80100000000))/10000 )*10000;
    if (code == 805) return (mh-80500000000)  - ( ((Int_t)(mh-80500000000))/10000 )*10000;
    if (code == 806) return (mh-80600000000)  - ( ((Int_t)(mh-80600000000))/10000 )*10000;
    return -1;
}

int code(double mh){
    return (int)(mh/100000000);
}

/////
static bool addPreliminary = false;
static bool saveOutputFile = true;
static bool addICHEPContours = false;

void plotScalar_1D(string inputFileName, string outputDIR, int dmMass = 1, bool isDMF = false, string coupling = "1",string postfix = "COMB") {
  
  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  // Set the color palette
  bool useNicksPalette = true;
  int ncontours        = 999;

  if (useNicksPalette) {    
    TColor::InitializeColors();
    Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
    Double_t red[9]   = { 243./255., 243./255., 240./255., 240./255., 241./255., 239./255., 186./255., 151./255., 129./255.};
    Double_t green[9] = {   0./255.,  46./255.,  99./255., 149./255., 194./255., 220./255., 183./255., 166./255., 147./255.};
    Double_t blue[9]  = {   6./255.,   8./255.,  36./255.,  91./255., 169./255., 235./255., 246./255., 240./255., 233./255.};
    TColor::CreateGradientColorTable(9, stops, red, green, blue, ncontours);
  }
  else gStyle->SetPalette(70);  
  gStyle->SetNumberContours(ncontours);
  
  // This is where all the plots are made
  TFile *file = TFile::Open(inputFileName.c_str(),"READ");
  TTree *tree = (TTree*)file->Get("limit");

  TGraph* grobs = new TGraph();
  TGraph* grexp = new TGraph();
  TGraph* grexp_1sigma_up   = new TGraphErrors();
  TGraph* grexp_2sigma_up   = new TGraphErrors();
  TGraph* grexp_1sigma_dw   = new TGraphErrors();
  TGraph* grexp_2sigma_dw   = new TGraphErrors();

  double mh;
  double limit;
  float  quantile;
  
  tree->SetBranchAddress("mh",&mh);
  tree->SetBranchAddress("limit",&limit);
  tree->SetBranchAddress("quantileExpected",&quantile);
  
  int expcounter          = 0;
  int exp_up_counter_1s   = 0;
  int exp_down_counter_1s = 0;
  int exp_up_counter_2s   = 0;
  int exp_down_counter_2s = 0;
  int obscounter          = 0;

  double medMin = 100000;
  double medMax = 0;

  cout<<"Loop on the limit tree entries: mass points and quantiles "<<endl;

  vector<float> medMassList;

  for (int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);    
    int c        = code(mh);
    int medmass  = mmed(mh,c);
    int dmmass   = mdm(mh,c);
    
    // for plotting reasons --> smooth plot
    if (medmass < 2* dmmass) continue; // skip off-shell points
    if (dmmass != dmMass) continue; // skip points not belonging to the selected DM mass
    if (medmass > 600) continue;

    // fill expected limit graph
    if(quantile != -1 and not isDMF){
      if (medmass == 30  and dmmass == 1) continue;
      if (medmass == 100 and dmmass == 1) continue;
      if (medmass == 125 and dmmass == 1) continue;
      if (medmass == 315 and dmmass == 1) continue;
    }
    else if (quantile == -1 and not isDMF) {      
      if (medmass == 40  and dmmass == 1) continue;
      if (medmass == 175 and dmmass == 1) continue;
      if (medmass == 315 and dmmass == 1) continue;
    }

    if (quantile == 0.5) {      
      grexp->SetPoint(expcounter, double(medmass), limit);
      expcounter++;
      // find max and min for frame
      if(medmass < medMin)
	medMin = medmass;
      if(medmass > medMax)
	medMax = medmass;      
      medMassList.push_back(medmass);
    }
    else if (quantile == -1){
      grobs->SetPoint(obscounter, double(medmass), limit);
      obscounter++;
    }

    // 1 sigma dw
    else if (quantile < 0.17 && quantile > 0.15 ) {
      grexp_1sigma_dw->SetPoint(exp_down_counter_1s, double(medmass), limit);      
      exp_down_counter_1s++;
    }
    // 1 sigma up
    else if (quantile < 0.85 && quantile > 0.83 ) {
      grexp_1sigma_up->SetPoint(exp_up_counter_1s, double(medmass), limit);      
      exp_up_counter_1s++;
    }

    // 2 sigma dw
    else if (quantile < 0.04 && quantile > 0.02 ) {
      grexp_2sigma_dw->SetPoint(exp_down_counter_2s, double(medmass), limit);      
      exp_down_counter_2s++;
    }
    // 2 sigma up
    else if (quantile < 0.98 && quantile > 0.96 ) {
      grexp_2sigma_up->SetPoint(exp_up_counter_2s, double(medmass), limit);      
      exp_up_counter_2s++;
    }    
  }

  tree->ResetBranchAddresses();

  //// make a spline
  TSpline3 *splineexp = new TSpline3("splineexp",grexp->GetX(),grexp->GetY(),grexp->GetN());
  splineexp->SetLineColor(kBlack);
  splineexp->SetLineStyle(7);
  splineexp->SetLineWidth(2);

  TSpline3 *splineobs = new TSpline3("splineobs",grobs->GetX(),grobs->GetY(),grobs->GetN());
  splineobs->SetLineColor(kBlack);
  splineobs->SetLineWidth(2);

  TGraphAsymmErrors* graph_1sigma_band = new TGraphAsymmErrors();
  TGraphAsymmErrors* graph_2sigma_band = new TGraphAsymmErrors();

  if(exp_up_counter_1s == exp_down_counter_1s and exp_down_counter_1s == expcounter){

    for(int iPoint = 0; iPoint < exp_up_counter_1s; iPoint++){
      double x_central, y_central;
      grexp->GetPoint(iPoint,x_central,y_central);
      graph_1sigma_band->SetPoint(iPoint,x_central,y_central);
      double y_up, y_dw;
      grexp_1sigma_dw->GetPoint(iPoint,x_central,y_dw);
      grexp_1sigma_up->GetPoint(iPoint,x_central,y_up);
      float rangeDw = 0;
      float rangeUp = 0;
      if(iPoint == 0){
	rangeUp = (medMassList.at(iPoint+1)-medMassList.at(iPoint))/2;
      }
      else if(iPoint == exp_up_counter_1s-1){
	rangeDw = (medMassList.at(iPoint)-medMassList.at(iPoint-1))/2;
      }
      else{
	rangeUp = (medMassList.at(iPoint+1)-medMassList.at(iPoint))/2;
	rangeDw = (medMassList.at(iPoint)-medMassList.at(iPoint-1))/2;
      }
      graph_1sigma_band->SetPointError(iPoint,rangeDw,rangeUp,fabs(y_dw-y_central),fabs(y_up-y_central));      
    }
  }
  else {
    cerr<<"Number of expected limits value: mediat, 1-sigma up and 1-sigma down don't match --> skip "<<endl;
    return;
  }

  if(exp_up_counter_2s == exp_down_counter_2s and exp_down_counter_2s == expcounter){

    for(int iPoint = 0; iPoint < exp_up_counter_2s; iPoint++){
      double x_central, y_central;
      grexp->GetPoint(iPoint,x_central,y_central);
      graph_2sigma_band->SetPoint(iPoint,x_central,y_central);
      double y_up, y_dw;
      grexp_2sigma_dw->GetPoint(iPoint,x_central,y_dw);
      grexp_2sigma_up->GetPoint(iPoint,x_central,y_up);
      float rangeDw = 0;
      float rangeUp = 0;
      if(iPoint == 0){
	rangeUp = (medMassList.at(iPoint+1)-medMassList.at(iPoint))/2;
      }
      else if(iPoint == exp_up_counter_1s-1){
	rangeDw = (medMassList.at(iPoint)-medMassList.at(iPoint-1))/2;
      }
      else{
	rangeUp = (medMassList.at(iPoint+1)-medMassList.at(iPoint))/2;
	rangeDw = (medMassList.at(iPoint)-medMassList.at(iPoint-1))/2;
      }
      graph_2sigma_band->SetPointError(iPoint,rangeDw,rangeUp,fabs(y_dw-y_central),fabs(y_up-y_central));      
    }
  }
  else {
    cerr<<"Number of expected limits value: mediat, 2-sigma up and 2-sigma down don't match --> skip "<<endl;
    return;
  }
  
  
  //////////// All the plotting and cosmetics
  TCanvas* canvas = new TCanvas("canvas", "canvas",625,600);
  TH1* frame = canvas->DrawFrame(min(medMin,0.),TMath::MinElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*0.5,
				 medMax,TMath::MaxElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*1.5, "");
  frame->GetXaxis()->SetTitle("m_{med} [GeV]");
  frame->GetYaxis()->SetTitle("95%  CL upper limit on #sigma/#sigma_{theory}");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(0.95);
  frame->Draw();

  if(not addPreliminary)
    CMS_lumi(canvas,"35.9",false,true);
  else
    CMS_lumi(canvas,"35.9",false,false);

  graph_2sigma_band->SetFillColor(kOrange);
  graph_1sigma_band->SetFillColor(kGreen+1);
  graph_2sigma_band->SetLineColor(kOrange);
  graph_1sigma_band->SetLineColor(kGreen+1);

  graph_2sigma_band->Draw("3same");
  graph_1sigma_band->Draw("3same");

  TGraph* graph_obs_ichep = NULL;
  TGraph* graph_exp_ichep = NULL;
  if(addICHEPContours){
    TFile* icheplines = TFile::Open("externalFiles/monojet_S_1D_ICHEP2016.root","READ");
    graph_obs_ichep = (TGraph*) icheplines->Get("obs");
    graph_obs_ichep->SetLineWidth(2);
    graph_obs_ichep->SetLineStyle(1);
    graph_obs_ichep->SetLineColor(kBlue);
    graph_obs_ichep->Draw("Lsame");
    graph_exp_ichep = (TGraph*) icheplines->Get("exp");
    graph_exp_ichep->SetLineWidth(2);
    graph_exp_ichep->SetLineStyle(7);
    graph_exp_ichep->SetLineColor(kBlue);
    graph_exp_ichep->Draw("Lsame");
  }  

  splineexp->Draw("Lsame");
  splineobs->Draw("Lsame");

  TF1* line = new TF1 ("line","1",min(medMin,0.),medMax);
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("L same");

  TLegend *leg = NULL;
  if(not addICHEPContours)
    leg =new TLegend(0.165,0.48,0.44,0.77);
  else
    leg = new TLegend(0.165,0.44,0.52,0.77);  

  leg->AddEntry(splineexp,"Median expected","L");
  leg->AddEntry(graph_1sigma_band,"68% expected","F");
  leg->AddEntry(graph_2sigma_band,"95% expected","F");
  leg->AddEntry(splineobs,"Observed","L");
  if(addICHEPContours){
    leg->AddEntry(graph_obs_ichep,"1703.01651 observed","L");
    leg->AddEntry(graph_exp_ichep,"1703.01651 expected","L");
  }
  leg->AddEntry(line,"#mu = 1","L");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("SAME");
  
  TLatex * tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.030);
  tex->Draw();
  if (coupling == "1")
    tex->DrawLatex(0.175,0.80,("#bf{Scalar med, Dirac DM, m_{DM} = "+to_string(dmMass)+" GeV g_{q} = 1, g_{DM} = 1}").c_str());
  else
    tex->DrawLatex(0.175,0.80,("#bf{Scalar med, Dirac DM, m_{DM} = "+to_string(dmMass)+" GeV g_{q} = 0.25, g_{DM} = 1}").c_str());
  
  gPad->RedrawAxis();
  gPad->Modified(); 
  gPad->Update();
  
  frame->GetYaxis()->SetRangeUser(0.001,7);
  canvas->SaveAs((outputDIR+"/scan_scalar_1D_dmMass_"+to_string(dmMass)+"_g"+string(coupling)+"_"+postfix+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/scan_scalar_1D_dmMass_"+to_string(dmMass)+"_g"+string(coupling)+"_"+postfix+".png").c_str(),"png");

  canvas->SetLogy();
  frame->GetYaxis()->SetRangeUser(TMath::MinElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*0.01,
				  TMath::MaxElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*100);
  canvas->SaveAs((outputDIR+"/scan_scalar_1D_dmMass_"+to_string(dmMass)+"_g"+string(coupling)+"_"+postfix+"_log.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/scan_scalar_1D_dmMass_"+to_string(dmMass)+"_g"+string(coupling)+"_"+postfix+"_log.png").c_str(),"png");

  if(saveOutputFile){
    TFile* outputFile = new TFile((outputDIR+"/limit_scalar_1D.root").c_str(),"RECREATE");
    outputFile->cd();
    splineexp->Write("expected_limit");
    splineobs->Write("observed_limit");
    outputFile->Close();
  }

}
