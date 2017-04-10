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

void plotScalar_1D(string inputFileName, string outputDIR, int dmMass = 1, string coupling = "025",string postfix = "COMB") {
  
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

  int medMin = 100000;
  int medMax = 0;

  cout<<"Loop on the limit tree entries: mass points and quantiles "<<endl;

  vector<float> medMassList;

  for (int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);    
    int c        = code(mh);
    int medmass  = mmed(mh,c);
    int dmmass   = mdm(mh,c);
    
    if (medmass < 2* dmmass) continue; // skip off-shell points
    if (dmmass != dmMass) continue; // skip points not belonging to the selected DM mass

    // for plotting reasons
    if (medmass > 600) continue;
    if (medmass == 90) continue;
    
    // fill expected limit graph
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

    else if (quantile == -1) {      
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

  //cout<<"Found: NexpLim "<<expcounter<<" NObsLim "<<obscounter<<" 1sigmaUp "<<exp_up_counter_1s<<" 1sigmaDw "<<exp_down_counter_1s<<" 2sigmaUp "<<exp_up_counter_2s<<" 2sigmDw "<<exp_up_counter_2s<<" for mDM "<<dmMass<<" medMin "<<medMin<<" medMax "<<medMax<<endl;

  tree->ResetBranchAddresses();

  // Make 1 and 2 sigma brazilian bands
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
  TCanvas* canvas = new TCanvas("canvas", "canvas",600,600);
  TH1* frame = canvas->DrawFrame(medMin,TMath::MinElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*0.5,
				 medMax,TMath::MaxElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*1.5, "");
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle("m_{med} [GeV]");
  frame->GetYaxis()->SetTitle("95%  CL upper limit on #sigma/#sigma_{theory}");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.10);  
  frame->Draw();
  CMS_lumi(canvas,"35.9");

  graph_2sigma_band->SetFillColor(kOrange);
  graph_1sigma_band->SetFillColor(kGreen+1);
  graph_2sigma_band->SetLineColor(kOrange);
  graph_1sigma_band->SetLineColor(kGreen+1);
  
  graph_2sigma_band->Draw("3same");
  graph_1sigma_band->Draw("3same");

  grexp->SetLineColor(kBlack);
  grexp->SetLineStyle(7);
  grexp->SetLineWidth(2);
  grexp->Draw("Csame");

  grobs->SetLineColor(kBlack);
  grobs->SetLineWidth(2);
  grobs->Draw("Csame");

  TF1* line = new TF1 ("line","1",medMin,medMax);
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("L same");

  TLegend *leg = new TLegend(0.175,0.5,0.57,0.77);  
  leg->AddEntry(grobs,"Observed 95% CL","L");
  leg->AddEntry(grexp,"Median expected 95% CL","L");
  leg->AddEntry(graph_1sigma_band,"68% expected","F");
  leg->AddEntry(graph_2sigma_band,"95% expected","F");
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
    tex->DrawLatex(0.175,0.80,("#bf{Scalar med, Dirac DM, m_{med} = "+to_string(dmMass)+" GeV g_{q} = 1, g_{DM} = 1}").c_str());
  else
    tex->DrawLatex(0.175,0.80,("#bf{Scalar med, Dirac DM, m_{med} = "+to_string(dmMass)+" GeV g_{q} = 0.25, g_{DM} = 1}").c_str());
  
  gPad->RedrawAxis();
  gPad->Modified(); 
  gPad->Update();
  
  frame->GetYaxis()->SetRangeUser(0,6);
  canvas->SaveAs((outputDIR+"/scan_scalar_1D_dmMass_"+to_string(dmMass)+"_g"+string(coupling)+"_"+postfix+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/scan_scalar_1D_dmMass_"+to_string(dmMass)+"_g"+string(coupling)+"_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/scan_scalar_1D_dmMass_"+to_string(dmMass)+"_g"+string(coupling)+"_"+postfix+".C").c_str(),"C");

  canvas->SetLogy();
  frame->GetYaxis()->SetRangeUser(TMath::MinElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*0.01,
				  TMath::MaxElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*100);
  canvas->SaveAs((outputDIR+"/scan_scalar_1D_dmMass_"+to_string(dmMass)+"_g"+string(coupling)+"_"+postfix+"_log.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/scan_scalar_1D_dmMass_"+to_string(dmMass)+"_g"+string(coupling)+"_"+postfix+"_log.png").c_str(),"png");
}
