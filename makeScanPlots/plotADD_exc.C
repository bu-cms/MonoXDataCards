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
#include "Math/Functor.h"
#include "Math/RootFinder.h"
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

  TGraph* grobs = new TGraph();
  TGraph* grexp = new TGraph();
  TGraph* grexp_1sigma_up   = new TGraphErrors();
  TGraph* grexp_2sigma_up   = new TGraphErrors();
  TGraph* grexp_1sigma_dw   = new TGraphErrors();
  TGraph* grexp_2sigma_dw   = new TGraphErrors();

double func_obs(double x) {
    //cout <<b<<end;
    return grobs->Eval(x)-1.;
}
double func_exp(double x) {
    return grexp->Eval(x)-1.;
}
double func_exp_1sigma_up(double x) {
    return grexp_1sigma_up->Eval(x)-1.;
}
double func_exp_2sigma_up(double x) {
    return grexp_2sigma_up->Eval(x)-1.;
}
double func_exp_1sigma_dw(double x) {
    return grexp_1sigma_dw->Eval(x)-1.;
}
double func_exp_2sigma_dw(double x) {
    return grexp_2sigma_dw->Eval(x)-1.;
}

void  setLimit(string inputFileName,int dj) {
  
  
  // This is where all the plots are made
  TFile *file = TFile::Open(inputFileName.c_str(),"READ");
  TTree *tree = (TTree*)file->Get("limit");
  int dmMass=0;
  dmMass=dj;
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

  //cout<<"Loop on the limit tree entries: mass points and quantiles "<<dmMass<<endl;
  int phi[]={3,4,5,6,7,8,9,10,11,12};
  int phi_c=10;
  vector<float> medMassList;
  if (dj == 2) { int phi_c=8;}
  if (dj == 3) { int phi_c=8;}
  if (dj == 4) {int phi_c=7;}
  if (dj == 5 or dj == 6) {int phi_c=7;}
  
  for(int iP = 0; iP < phi_c; iP++){
  for (int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);    
    int c        = code(mh);
    int dd  = dim(c);
    int dmmass   = mdm(c);
    if (dd != dj) continue; // skip points not belonging to the selected Dimension

    if(phi[iP]!=dmmass) continue;  // skip points not belonging to the selected DM mass   //yg

//    cout <<" dd ="<< dd<<" dmmass = "<<dmmass<<" dj "<<dj<<endl;

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
  } //end of tree
}  // end of phi loop
  tree->ResetBranchAddresses();
}

void plotADD_exc(string inputFileName, string outputDIR) {

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


   TH1F *h_obs = new TH1F("h_obs","",5,1.5,6.5);
   TH1F *h_exp = new TH1F("h_exp","",5,1.5,6.5);
   TH1F *h_exp_1sigma_up = new TH1F("h_exp_1sigma_up","",5,1.5,6.5);
   TH1F *h_exp_2sigma_up = new TH1F("h_exp_2sigma_up","",5,1.5,6.5);
   TH1F *h_exp_1sigma_dw = new TH1F("h_exp_1sigma_dw","",5,1.5,6.5);
   TH1F *h_exp_2sigma_dw = new TH1F("h_exp_2sigma_dw","",5,1.5,6.5);
   TH1F *h_obs_8TeV = new TH1F("h_obs_8TeV","",5,1.5,6.5);
   double obs_8TeV[5] = {5.61,4.38,3.86,3.55,3.26};

  int counter_f=0;
  for (int dj=2; dj<7; dj++){
   setLimit(inputFileName,dj);
   ROOT::Math::Functor1D f_ADD(&func_obs);
   ROOT::Math::RootFinder rfbis(ROOT::Math::RootFinder::kGSL_BISECTION );
   if (dj==2) rfbis.SetFunction(f_ADD, 5., 13.);
   if (dj==3) rfbis.SetFunction(f_ADD, 4., 11.);
   if (dj==4) rfbis.SetFunction(f_ADD, 3., 9.);
   if (dj==5 or dj==6) rfbis.SetFunction(f_ADD, 3., 7.);
   rfbis.Solve();
   h_obs->SetBinContent(counter_f+1.5,rfbis.Root());
   h_obs_8TeV->SetBinContent(counter_f+1.5,obs_8TeV[counter_f]);
   cout <<" d= "<<dj<<" obs ="<<rfbis.Root()<<endl;
   ++counter_f;
  }
  counter_f=0;
  for (int dj=2; dj<7; dj++){
   setLimit(inputFileName,dj);
   ROOT::Math::Functor1D f_ADD(&func_exp);
   ROOT::Math::RootFinder rfbis(ROOT::Math::RootFinder::kGSL_BISECTION );
   if (dj==2) rfbis.SetFunction(f_ADD, 5., 13.);
   if (dj==3) rfbis.SetFunction(f_ADD, 4., 11.);
   if (dj==4) rfbis.SetFunction(f_ADD, 3., 9.);
   if (dj==5 or dj==6) rfbis.SetFunction(f_ADD, 3., 7.);
   rfbis.Solve();
   h_exp->SetBinContent(counter_f+1.5,rfbis.Root());
   cout <<" d= "<<dj<<" exp ="<<rfbis.Root()<<endl;
   ++counter_f;
  }
  counter_f=0;
  for (int dj=2; dj<7; dj++){
   setLimit(inputFileName,dj);
   ROOT::Math::Functor1D f_ADD(&func_exp_1sigma_up);
   ROOT::Math::RootFinder rfbis(ROOT::Math::RootFinder::kGSL_BISECTION );
   if (dj==2) rfbis.SetFunction(f_ADD, 5., 13.);
   if (dj==3) rfbis.SetFunction(f_ADD, 4., 11.);
   if (dj==4) rfbis.SetFunction(f_ADD, 3., 9.);
   if (dj==5 or dj==6) rfbis.SetFunction(f_ADD, 3., 7.);
   rfbis.Solve();
   h_exp_1sigma_up->SetBinContent(counter_f+1.5,rfbis.Root());
   cout <<" d= "<<dj<<" h_exp_1sigma_up ="<<rfbis.Root()<<endl;
   ++counter_f;
  }

  counter_f=0;
  for (int dj=2; dj<7; dj++){
   setLimit(inputFileName,dj);
   ROOT::Math::Functor1D f_ADD(&func_exp_2sigma_up);
   ROOT::Math::RootFinder rfbis(ROOT::Math::RootFinder::kGSL_BISECTION );
   if (dj==2) rfbis.SetFunction(f_ADD, 5., 13.);
   if (dj==3) rfbis.SetFunction(f_ADD, 4., 11.);
   if (dj==4) rfbis.SetFunction(f_ADD, 3., 9.);
   if (dj==5 or dj==6) rfbis.SetFunction(f_ADD, 3., 7.);
   rfbis.Solve();
   h_exp_2sigma_up->SetBinContent(counter_f+1.5,rfbis.Root());
   cout <<" d= "<<dj<<" h_exp_2sigma_up ="<<rfbis.Root()<<endl;
   ++counter_f;
  }

  counter_f=0;
  for (int dj=2; dj<7; dj++){
   setLimit(inputFileName,dj);
   ROOT::Math::Functor1D f_ADD(&func_exp_1sigma_dw);
   ROOT::Math::RootFinder rfbis(ROOT::Math::RootFinder::kGSL_BISECTION );
   if (dj==2) rfbis.SetFunction(f_ADD, 5., 13.);
   if (dj==3) rfbis.SetFunction(f_ADD, 4., 11.);
   if (dj==4) rfbis.SetFunction(f_ADD, 3., 9.);
   if (dj==5 or dj==6) rfbis.SetFunction(f_ADD, 3., 7.);
   rfbis.Solve();
   h_exp_1sigma_dw->SetBinContent(counter_f+1.5,rfbis.Root());
   cout <<" d= "<<dj<<" h_exp_1sigma_dw ="<<rfbis.Root()<<endl;
   ++counter_f;
  }

  counter_f=0;
  for (int dj=2; dj<7; dj++){
   setLimit(inputFileName,dj);
   ROOT::Math::Functor1D f_ADD(&func_exp_2sigma_dw);
   ROOT::Math::RootFinder rfbis(ROOT::Math::RootFinder::kGSL_BISECTION );
   if (dj==2) rfbis.SetFunction(f_ADD, 5., 13.);
   if (dj==3) rfbis.SetFunction(f_ADD, 4., 11.);
   if (dj==4) rfbis.SetFunction(f_ADD, 3., 9.);
   if (dj==5 or dj==6) rfbis.SetFunction(f_ADD, 3., 7.);
   rfbis.Solve();
   h_exp_2sigma_dw->SetBinContent(counter_f+1.5,rfbis.Root());
   cout <<" d= "<<dj<<" h_exp_2sigma_dw ="<<rfbis.Root()<<endl;
   ++counter_f;
  }
 
  //////////// All the plotting and cosmetics
  TCanvas* canvas = new TCanvas("canvas", "canvas",600,600);
  TH1* frame = canvas->DrawFrame(0,0,7,17, "");
  h_exp_2sigma_dw->GetXaxis()->SetTitle("Number of Extra Dimensions");
  h_exp_2sigma_dw->GetYaxis()->SetTitle("M_{D} [TeV]");
  h_exp_2sigma_dw->GetXaxis()->SetTitleOffset(1.15);
  h_exp_2sigma_dw->GetYaxis()->SetTitleOffset(1.10);  
  //CMS_lumi(canvas,"35.9");
  CMS_lumi(canvas,"35.9",false,true,false,0.37,0);

  h_exp_2sigma_dw->GetYaxis()->SetRangeUser(0,16);
  h_exp_2sigma_dw->GetXaxis()->SetNdivisions(5);
  h_exp_2sigma_dw->SetLineColor(kYellow);
  h_exp_2sigma_dw->SetFillColor(kYellow);
  h_exp_2sigma_dw->Draw(); 
  h_exp_1sigma_dw->SetLineColor(kGreen);
  h_exp_1sigma_dw->SetFillColor(kGreen);
  h_exp_1sigma_dw->Draw("sames");

  h_exp_2sigma_up->SetLineColor(kGreen);
  h_exp_2sigma_up->SetFillColor(kGreen);
  h_exp_2sigma_up->Draw("sames");
  h_exp_1sigma_up->SetLineColor(kYellow);
  h_exp_1sigma_up->SetFillColor(kYellow);
  h_exp_1sigma_up->Draw("sames");
  h_exp_2sigma_up->SetLineColor(10);
  h_exp_2sigma_up->SetFillColor(10);
  h_exp_2sigma_up->Draw("sames");
  
  h_exp->SetLineColor(kRed);
  h_exp->SetLineWidth(2);
  h_exp->SetLineStyle(2);
  h_exp->Draw("sames");
  h_obs->SetLineWidth(2);
  h_obs->Draw("sames");
  h_obs_8TeV->SetLineWidth(2);
  h_obs_8TeV->SetLineColor(4);
  h_obs_8TeV->Draw("hist sames");

  //TLegend *leg = new TLegend(0.175,0.5,0.57,0.77);  
  TLegend *leg = new TLegend(0.5317726,0.5706806,0.9264214,0.8411867,NULL,"brNDC");
  leg->AddEntry(h_obs,"Observed 95% CL","L");
  leg->AddEntry(h_exp,"Median expected 95% CL","L");
  leg->AddEntry(h_exp_1sigma_dw,"68% expected","F");
  leg->AddEntry(h_exp_1sigma_up,"95% expected","F");
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->Draw("SAME");
  
   TLatex *   tex = new TLatex(0.94,0.95,"35.9 fb^{-1} (13 TeV)");
   tex->SetNDC();
   tex->SetTextAlign(31);
   tex->SetTextFont(42);
   tex->SetTextSize(0.042);
   tex->SetLineWidth(2);
   tex->Draw();
   tex = new TLatex(0.545,0.86,"CMS");
   tex->SetNDC();
   tex->SetTextSize(0.042);
   tex->SetLineWidth(2);
   tex->Draw();
  gPad->RedrawAxis();
  gPad->Modified(); 
  gPad->Update();
  
  canvas->SaveAs((outputDIR+"/Exc_ADD"+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/Exc_ADD"+".png").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/Exc_ADD"+".C").c_str(),"pdf");

}

