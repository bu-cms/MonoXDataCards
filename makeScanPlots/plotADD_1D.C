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

int dim(int code){
    if (code == 25 or code == 26 or code == 27 or code == 28 or code == 29 or code == 210 or code == 211 or code == 212) return ((int)(2));

    if (code == 34 or code == 35 or code == 36 or code == 37 or code == 38 or code == 39 or code == 310 or code == 311) return ((int)(3));

    if (code == 43 or code == 44 or code == 45 or code == 46 or code == 47 or code == 48 or code == 49) return ((int)(4));

    if (code == 53 or code == 54 or code == 55 or code == 56 or code == 57) return ((int)(5));

    if (code == 63 or code == 64 or code == 65 or code == 66 or code == 67) return ((int)(6));

    return -1;
}

int mdm(int code){
    if (code == 43 or code == 53 or code == 63) return ((int)(3));

    if (code == 34 or code == 44 or code == 54 or code == 64) return ((int)(4));
 
    if (code == 25 or code == 35 or code == 45 or code == 55 or code == 65) return ((int)(5));
 
    if (code == 26 or  code == 36 or code == 46 or code == 56 or code == 66) return ((int)(6));

    if (code == 27 or code == 37 or code == 47 or code == 57 or code == 67)return ((int)(7));

    if (code == 28 or code == 38 or code == 48) return ((int)(8));

    if (code == 29 or code == 39 or code == 49) return ((int)(9));

    if (code == 210 or code == 310) return ((int)(10)); 

    if (code == 211 or code == 311) return ((int)(11));  

    if (code == 212) return ((int)(12));
    return -1;
}

int code(double mh){
    return (int)(mh);
}

void plotADD_1D(string inputFileName, string outputDIR, int dmMass = 1) {
  
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

  int medMin = 100;
  int medMax = 0;

  cout<<"Loop on the limit tree entries: mass points and quantiles "<<dmMass<<endl;
  int phi[]={3,4,5,6,7,8,9,10,11,12};
  int phi_c=10;
  vector<float> medMassList;
  if (dmMass == 2) { int phi_c=8;}
  if (dmMass == 3) { int phi_c=8;}
  if (dmMass == 4) {int phi_c=7;}
  if (dmMass == 5 or dmMass == 6) {int phi_c=7;}
  
  for(int iP = 0; iP < phi_c; iP++){
  for (int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);    
    int c        = code(mh);
    int dd  = dim(c);
    int dmmass   = mdm(c);
    if (dd != dmMass) continue; // skip points not belonging to the selected Dimension

    if(phi[iP]!=dmmass) continue;  // skip points not belonging to the selected DM mass   //yg
    // fill expected limit graph
    if (quantile == 0.5) {

      grexp->SetPoint(expcounter, double(dmmass), limit);
      expcounter++;
      // find max and min for frame
      if(dmmass < medMin)
	medMin = dmmass;

      if(dmmass > medMax)
	medMax = dmmass;

      medMassList.push_back(dmmass);
    }

    else if (quantile == -1) {      
      grobs->SetPoint(obscounter, double(dmmass), limit);
      obscounter++;
    }

    // 1 sigma dw
    else if (quantile < 0.17 && quantile > 0.15 ) {
      grexp_1sigma_dw->SetPoint(exp_down_counter_1s, double(dmmass), limit);      
      exp_down_counter_1s++;
    }
    // 1 sigma up
    else if (quantile < 0.85 && quantile > 0.83 ) {
      grexp_1sigma_up->SetPoint(exp_up_counter_1s, double(dmmass), limit);      
      exp_up_counter_1s++;
    }

    // 2 sigma dw
    else if (quantile < 0.04 && quantile > 0.02 ) {
      grexp_2sigma_dw->SetPoint(exp_down_counter_2s, double(dmmass), limit);      
      exp_down_counter_2s++;
    }
    // 2 sigma up
    else if (quantile < 0.98 && quantile > 0.96 ) {
      grexp_2sigma_up->SetPoint(exp_up_counter_2s, double(dmmass), limit);      
      exp_up_counter_2s++;
    }    
  }
}
  tree->ResetBranchAddresses();
//}
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
      //}

      double x_obs, y_obs;
      grobs->GetPoint(iPoint,x_obs,y_obs);
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
  frame->GetXaxis()->SetTitle("M_{D} [TeV]");
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
  grexp->Draw("Lsame");

  grobs->SetLineColor(kBlack);
  grobs->SetLineWidth(2);
  grobs->Draw("Lsame");

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
      tex->DrawLatex(0.175,0.80,("#bf{ ADD d = "+to_string(dmMass)+"}").c_str());

  gPad->RedrawAxis();
  gPad->Modified(); 
  gPad->Update();
  
  canvas->SaveAs((outputDIR+"/scan_ADD_1D_dmMass_"+to_string(dmMass)+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/scan_ADD_1D_dmMass_"+to_string(dmMass)+".png").c_str(),"pdf");

  canvas->SetLogy();
  frame->GetYaxis()->SetRangeUser(TMath::MinElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*0.1,
				  TMath::MaxElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*200);
  canvas->SaveAs((outputDIR+"/scan_ADD_1D_dmMass_"+to_string(dmMass)+"_log.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/scan_ADD_1D_dmMass_"+to_string(dmMass)+"_log.png").c_str(),"pdf");
}

