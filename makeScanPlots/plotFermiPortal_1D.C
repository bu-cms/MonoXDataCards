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

int mmed(double mh){

 int code = int(mh);
 if (code == 5010) return ((int)(1));
 if (code == 50010 or code == 500100 or code == 500100 or code == 500200 or code == 500300 or code == 500400 or code == 500500 or
     code == 500600 or code == 500700 or code == 500800 or code == 5001000 )
   return ((int)(500));

 if (code == 55010  or code == 550100  or code == 550200 or code == 550300 or code == 550400 or code == 550500 or
     code == 550600 or code == 550700  or code == 550800 or code == 5501000 ) return ((int)(550));

 if (code == 60010  or code == 600100  or code == 600200 or code == 600300 or code == 600400 or code == 600500 or
     code == 600600 or code == 600700  or code == 600800 or code == 6001000 ) return ((int)(600));

 if (code == 65010  or code == 650100  or code == 650200 or code == 650300 or code == 650400 or code == 650500 or
     code == 650600 or code == 650700  or code == 650800 or code == 6501000 ) return ((int)(650));

 if (code == 80010  or code == 800100  or code == 800200 or code == 800300 or code == 800400 or code == 800500 or
     code == 800600 or code == 800700  or code == 800800 or code == 8001000 ) return ((int)(800));

 if (code == 100010   or code == 1000100  or code == 1000200 or code == 1000300 or code == 1000400 or code == 1000500 or
     code == 1000600  or code == 1000700  or code == 1000800 or code == 10001000 ) return ((int)(1000));

 if (code == 120010   or code == 1200100  or code == 1200200 or code == 1200300 or code == 1200400 or code == 1200500 or
     code == 1200600  or code == 1200700  or code == 1200800 or code == 12001000 ) return ((int)(1200));

 if (code == 140010  or code == 1400100  or code == 1400200 or code == 1400300 or code == 1400400 or code == 1400500 or
     code == 1400600 or code == 1400700  or code == 1400800 or code == 14001000 ) return ((int)(1400));

 if (code == 150010  or code == 1500100  or code == 1500200 or code == 1500300 or code == 1500400 or code == 1500500 or
     code == 1500600 or code == 1500700  or code == 1500800 or code == 15001000 ) return ((int)(1500));

 if (code == 160010  or code == 1600100  or code == 1600200 or code == 1600300 or code == 1600400 or code == 1600500 or
     code == 1600600 or code == 1600700  or code == 1600800 or code == 16001000 ) return ((int)(1600));
 if (code == 2000400) return ((int)(2000));
    return -1;
}

int mdm(double mh){

  if (mh == 5010  or mh == 50010  or mh == 55010  or mh == 60010  or mh == 65010 or mh == 80010 or mh == 100010 or
      mh == 120010 or mh == 140010 or mh == 150010 or mh == 160010 )
    return ((int)(10));

  if (mh == 500100  or mh == 550100  or mh == 600100  or mh == 650100 or mh == 800100 or mh == 1000100 or
      mh == 1200100 or mh == 1400100 or mh == 1500100 or mh == 1600100 )
    return ((int)(100));

  if (mh == 500200  or mh == 550200  or mh == 600200  or mh == 650200 or mh == 800200 or mh == 1000200 or
      mh == 1200200 or mh == 1400200 or mh == 1500200 or mh == 1600200 )
    return ((int)(200));

  if (mh == 500300  or mh == 550300  or mh == 600300  or mh == 650300 or mh == 800300 or mh == 1000300 or
      mh == 1200300 or mh == 1400300 or mh == 1500300 or mh == 1600300 )
    return ((int)(300));

  if (mh == 500400  or mh == 550400  or mh == 600400  or mh == 650400 or mh == 800400 or mh == 1000400 or
      mh == 1200400 or mh == 1400400 or mh == 1500400 or mh == 1600400 )
    return ((int)(400));

  if (mh == 500500  or mh == 550500  or mh == 600500  or mh == 650500 or mh == 800500 or mh == 1000500 or
      mh == 1200500 or mh == 1400500 or mh == 1500500 or mh == 1600500 )
    return ((int)(500));

  if (mh == 500600  or mh == 550600  or mh == 600600  or mh == 650600 or mh == 800600 or mh == 1000600 or
      mh == 1200600 or mh == 1400600 or mh == 1500600 or mh == 1600600 )
    return ((int)(600));

  if (mh == 500700  or mh == 550700  or mh == 600700  or mh == 650700 or mh == 800700 or mh == 1000700 or
      mh == 1200700 or mh == 1400700 or mh == 1500700 or mh == 1600700 )
    return ((int)(700));

  if (mh == 500800  or mh == 550800  or mh == 600800  or mh == 650800 or mh == 800800 or mh == 1000800 or
      mh == 1200800 or mh == 1400800 or mh == 1500800 or mh == 1600800 )
    return ((int)(800));

  if (mh == 5001000  or mh == 5501000  or mh == 6001000  or mh == 6501000 or mh == 8001000 or mh == 10001000 or
      mh == 12001000 or mh == 14001000 or mh == 15001000 or mh == 16001000 )
    return ((int)(1000));
  if (mh == 2000400) return ((int)(400));
    return -1;
}

int code(double mh){
    return (int)(mh);
}

void plotFermiPortal_1D(string inputFileName, string outputDIR, int dmMass = 1) {
  
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
  int phi[]={500,550,600,650,800,1000,1200,1400,1500,1600};
  for(int iP = 0; iP < 10; iP++){
  for (int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);    
    int c        = code(mh);
    int medmass  = mmed(mh);
    int dmmass   = mdm(mh);
    if (medmass < 2* dmmass) continue; // skip off-shell points
    if (dmmass != dmMass) continue; // skip points not belonging to the selected DM mass

    if(medmass > 2500) continue;
    if(medmass < 500) continue;

    if(phi[iP]!=medmass) continue;
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
}
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
      }else{
	rangeUp = (medMassList.at(iPoint+1)-medMassList.at(iPoint))/2;
	rangeDw = (medMassList.at(iPoint)-medMassList.at(iPoint-1))/2;
      }

      double x_obs, y_obs;
      grobs->GetPoint(iPoint,x_obs,y_obs);
      graph_1sigma_band->SetPointError(iPoint,rangeDw,rangeUp,fabs(y_dw-y_central),fabs(y_up-y_central));      
      }
    //}
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
  frame->GetXaxis()->SetTitle("m_{#phi} [GeV]");
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
      tex->DrawLatex(0.175,0.80,("#bf{ m_{#chi} = "+to_string(dmMass)+" GeV {#lambda_{u}}= 1}").c_str());

  gPad->RedrawAxis();
  gPad->Modified(); 
  gPad->Update();
  
  canvas->SaveAs((outputDIR+"/scan_fPort_1D_dmMass_"+to_string(dmMass)+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/scan_fPort_1D_dmMass_"+to_string(dmMass)+".png").c_str(),"pdf");

  canvas->SetLogy();
  frame->GetYaxis()->SetRangeUser(TMath::MinElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*0.1,
				  TMath::MaxElement(graph_2sigma_band->GetN(),graph_2sigma_band->GetY())*200);
  canvas->SaveAs((outputDIR+"/scan_fPort_1D_dmMass_"+to_string(dmMass)+"_log.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/scan_fPort_1D_dmMass_"+to_string(dmMass)+"_log.png").c_str(),"pdf");
}

