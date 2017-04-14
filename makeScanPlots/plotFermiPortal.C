#include "../CMS_lumi.h"

int mmed(double mh){

  int code = int(mh);

  if (code == 50010 or code == 500100 or code == 500100 or code == 500200 or code == 500300 or code == 500400 or code == 500500 or
      code == 500600 or	code == 500700 or code == 500800 or code == 5001000 ) 
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
  

    return -1;
}

int mdm(double mh){

  if (mh == 50010  or mh == 55010  or mh == 60010  or mh == 65010 or mh == 80010 or mh == 100010 or
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


    return -1;
}

/////////
static bool saveOutputFile  = false;
static bool addRelicDensity = true;
static float nbinsX = 800;
static float nbinsY = 500;
static float minX = 0;
static float minY = 1.;
static float maxX = 1600;
static float maxY = 1000;
static float minZ = 0.01;
static float maxZ = 10;

TGraph* relic_gf();

void plotFermiPortal(string inputFileName, string outputDIR) {

  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  // Set the color palette
  bool useNicksPalette = false;
  int ncontours = 999;
  
  if (useNicksPalette) {
    
      TColor::InitializeColors();
      Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
      Double_t red[9]   = { 243./255., 243./255., 240./255., 240./255., 241./255., 239./255., 186./255., 151./255., 129./255.};
      Double_t green[9] = {   0./255.,  46./255.,  99./255., 149./255., 194./255., 220./255., 183./255., 166./255., 147./255.};
      Double_t blue[9]  = {   6./255.,   8./255.,  36./255.,  91./255., 169./255., 235./255., 246./255., 240./255., 233./255.};
      TColor::CreateGradientColorTable(9, stops, red, green, blue, ncontours);
  }
  else 
    gStyle->SetPalette(kBird);
  
  gStyle->SetNumberContours(ncontours);
  
  // This is where all the plots are made
  TFile *file  = TFile::Open(inputFileName.c_str(),"READ");
  TTree *tree  = (TTree*)file->Get("limit");

  TGraph* wm = NULL;  
  if(addRelicDensity)
    wm = relic_gf();
  
  TGraph2D* grexp = new TGraph2D();
  TGraph2D* grexp_up   = new TGraph2D();
  TGraph2D* grexp_down = new TGraph2D();
  TGraph2D* grobs = new TGraph2D();
  TGraph2D* grobu = new TGraph2D();
  TGraph2D* grobd = new TGraph2D();

  double mh;
  double limit;
  float quantile;
  
  tree->SetBranchAddress("mh",&mh);
  tree->SetBranchAddress("limit",&limit);
  tree->SetBranchAddress("quantileExpected",&quantile);
  
  int expcounter       = 0;
  int exp_up_counter   = 0;
  int exp_down_counter = 0;
  int obscounter       = 0;

  double minMedMass = 100000;
  double maxMedMass = -1;
  double minDMMass = 100000;
  double maxDMMass = -1;

  for (int i = 0; i < tree->GetEntries(); i++){

    tree->GetEntry(i);
    ////
    double medmass = mmed(mh);
    double dmmass  = mdm(mh);
    ////
    if(medmass < minMedMass)
      minMedMass = medmass;
    if(medmass > maxMedMass)
      maxMedMass = medmass;
    ////
    if(dmmass < minDMMass)
      minDMMass = dmmass;
    if(dmmass > maxDMMass)
      maxDMMass = dmmass;    

    ////
    if (quantile == 0.5) { // expected limit
      grexp->SetPoint(expcounter, double(medmass), double(dmmass), limit);
      expcounter++;
    }
    
    if (quantile < 0.17 && quantile > 0.14 ) { 
      grexp_down->SetPoint(exp_down_counter, double(medmass), double(dmmass), limit);
      exp_down_counter++;      
    }
    
    if (quantile < 0.85 && quantile > 0.83 ) {
      grexp_up->SetPoint(exp_up_counter, double(medmass), double(dmmass), limit);      
      exp_up_counter++;
    }

    if (quantile == -1) { // observed
      grobs->SetPoint(obscounter, double(medmass), double(dmmass), limit);
      grobu->SetPoint(obscounter, double(medmass), double(dmmass), limit*0.8);
      grobd->SetPoint(obscounter, double(medmass), double(dmmass), limit*1.2);
      obscounter++;      
    }
  }
  
  tree->ResetBranchAddresses();
  
  TH2D* hexp       = new TH2D("hexp", "",      nbinsX, minX, maxX, nbinsY, minY, maxY);
  TH2D* hexp_up    = new TH2D("hexp_up", "",   nbinsX, minX, maxX, nbinsY, minY, maxY);
  TH2D* hexp_down  = new TH2D("hexp_down", "", nbinsX, minX, maxX, nbinsY, minY, maxY);
  TH2D* hobs = new TH2D("hobs", "", nbinsX, minX, maxX, nbinsY, minY, maxY);
  TH2D* hobu = new TH2D("hobu", "", nbinsX, minX, maxX, nbinsY, minY, maxY);
  TH2D* hobd = new TH2D("hobd", "", nbinsX, minX, maxX, nbinsY, minY, maxY);

  // make granularity
  for (int i = 1; i   <= nbinsX; i++) {
    for (int j = 1; j <= nbinsY; j++) {
      hexp_up->SetBinContent(i,j,   grexp_up->Interpolate(hexp_up->GetXaxis()->GetBinCenter(i),hexp_up->GetYaxis()->GetBinCenter(j)));
      hexp_down->SetBinContent(i,j, grexp_down->Interpolate(hexp_down->GetXaxis()->GetBinCenter(i),hexp_down->GetYaxis()->GetBinCenter(j)));
      hexp->SetBinContent(i,j,grexp->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetYaxis()->GetBinCenter(j)));
      hobs->SetBinContent(i,j,grobs->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetYaxis()->GetBinCenter(j)));
      hobu->SetBinContent(i,j,grobu->Interpolate(hobu->GetXaxis()->GetBinCenter(i),hobu->GetYaxis()->GetBinCenter(j)));
      hobd->SetBinContent(i,j,grobd->Interpolate(hobd->GetXaxis()->GetBinCenter(i),hobd->GetYaxis()->GetBinCenter(j)));
    }
  }

  // set lower values for masses below minMedMass
  for(int i = 0; i < nbinsX; i++){
    if(hexp ->GetXaxis()->GetBinCenter(i) < minMedMass){
      for(int j = 0; j < nbinsY; j++){
	hexp->SetBinContent(i,j,hexp ->GetBinContent(hexp->FindBin(minMedMass+hexp->GetXaxis()->GetBinWidth(i),hexp ->GetYaxis()->GetBinCenter(j))));
	hexp_down->SetBinContent(i,j,hexp_down->GetBinContent(hexp_down->FindBin(minMedMass+hexp->GetXaxis()->GetBinWidth(i),hexp_down->GetYaxis()->GetBinCenter(j))));
	hexp_up->SetBinContent(i,j,hexp_up->GetBinContent(hexp_up->FindBin(minMedMass+hexp->GetXaxis()->GetBinWidth(i),hexp_up->GetYaxis()->GetBinCenter(j))));
	hobs->SetBinContent(i,j,hobs ->GetBinContent(hobs->FindBin(minMedMass+hexp->GetXaxis()->GetBinWidth(i),hobs ->GetYaxis()->GetBinCenter(j))));
	hobd->SetBinContent(i,j,hobd->GetBinContent(hobd->FindBin(minMedMass+hexp->GetXaxis()->GetBinWidth(i),hobd->GetYaxis()->GetBinCenter(j))));
	hobu->SetBinContent(i,j,hobu->GetBinContent(hobu->FindBin(minMedMass+hexp->GetXaxis()->GetBinWidth(i),hobu->GetYaxis()->GetBinCenter(j))));
      }
    }
  }

  // set lower values for masses below minDMmass
  for(int i = 0; i < nbinsY; i++){
    if(hexp ->GetYaxis()->GetBinCenter(i) < minDMMass){
      for(int j = 0; j < nbinsX; j++){
	hexp->SetBinContent(j,i,hexp ->GetBinContent(hexp->FindBin(hexp->GetXaxis()->GetBinCenter(j),minDMMass+hexp->GetYaxis()->GetBinWidth(j))));
	hexp_up->SetBinContent(j,i,hexp_up ->GetBinContent(hexp_up->FindBin(hexp_up->GetXaxis()->GetBinCenter(j),minDMMass+hexp->GetYaxis()->GetBinWidth(j))));
	hexp_down->SetBinContent(j,i,hexp_down ->GetBinContent(hexp_down->FindBin(hexp_down->GetXaxis()->GetBinCenter(j),minDMMass+hexp->GetYaxis()->GetBinWidth(j))));
	hobs->SetBinContent(j,i,hobs ->GetBinContent(hobs->FindBin(hobs->GetXaxis()->GetBinCenter(j),minDMMass+hexp->GetYaxis()->GetBinWidth(j))));
	hobu->SetBinContent(j,i,hobu ->GetBinContent(hobu->FindBin(hobu->GetXaxis()->GetBinCenter(j),minDMMass+hexp->GetYaxis()->GetBinWidth(j))));
	hobd->SetBinContent(j,i,hobd ->GetBinContent(hobd->FindBin(hobd->GetXaxis()->GetBinCenter(j),minDMMass+hexp->GetYaxis()->GetBinWidth(j))));
      }
    }
  }

  // set lower values for all off-shell points
  for(int i = 0; i < nbinsX; i++){
    for(int j = 0; j < nbinsY; j++){
      if(hexp->GetXaxis()->GetBinCenter(i) <= hexp->GetYaxis()->GetBinCenter(j) and (hexp->GetXaxis()->GetBinCenter(i) < minMedMass or hexp->GetXaxis()->GetBinCenter(i) > maxMedMass)){
	hexp->SetBinContent(i,j,maxZ);
	hexp_down->SetBinContent(i,j,maxZ);
	hexp_up->SetBinContent(i,j,maxZ);
	hobs->SetBinContent(i,j,maxZ);
	hobd->SetBinContent(i,j,maxZ);
	hobu->SetBinContent(i,j,maxZ);
      }
    }
  }
    
  /// fix remaining bad points
  for(int i = 0; i < nbinsX; i++){
    for(int j = 0; j < nbinsY; j++){
      if(hexp -> GetBinContent(i,j) <= 0) hexp->SetBinContent(i,j,maxZ);
      if(hexp_down -> GetBinContent(i,j) <= 0) hexp_down->SetBinContent(i,j,maxZ);
      if(hexp_up -> GetBinContent(i,j) <= 0) hexp_up->SetBinContent(i,j,maxZ);
      if(hobs -> GetBinContent(i,j) <= 0) hobs->SetBinContent(i,j,maxZ);
      if(hobu -> GetBinContent(i,j) <= 0) hobu->SetBinContent(i,j,maxZ);
      if(hobd -> GetBinContent(i,j) <= 0) hobd->SetBinContent(i,j,maxZ);

      if(hexp -> GetBinContent(i,j) > maxZ) hexp->SetBinContent(i,j,maxZ);
      if(hexp_down -> GetBinContent(i,j) > maxZ) hexp_down->SetBinContent(i,j,maxZ);
      if(hexp_up -> GetBinContent(i,j) > maxZ) hexp_up->SetBinContent(i,j,maxZ);
      if(hobs -> GetBinContent(i,j) > maxZ) hobs->SetBinContent(i,j,maxZ);
      if(hobu -> GetBinContent(i,j) > maxZ) hobu->SetBinContent(i,j,maxZ);
      if(hobd -> GetBinContent(i,j) > maxZ) hobd->SetBinContent(i,j,maxZ);

      if(hexp -> GetBinContent(i,j) < minZ) hexp->SetBinContent(i,j,minZ);
      if(hexp_down -> GetBinContent(i,j) < minZ) hexp_down->SetBinContent(i,j,minZ);
      if(hexp_up -> GetBinContent(i,j) < minZ) hexp_up->SetBinContent(i,j,minZ);
      if(hobs -> GetBinContent(i,j) < minZ) hobs->SetBinContent(i,j,minZ);
      if(hobu -> GetBinContent(i,j) < minZ) hobu->SetBinContent(i,j,minZ);
      if(hobd -> GetBinContent(i,j) < minZ) hobd->SetBinContent(i,j,minZ);
    }
  }

  hexp->Smooth();
  hexp_down->Smooth();
  hexp_up->Smooth();
  hobs->Smooth();
  hobu->Smooth();
  hobd->Smooth();

  ////////////////
  TH2* hexp2 = (TH2*)hexp->Clone("hexp2");
  TH2* hexp2_up = (TH2*)hexp_up->Clone("hexp2_up");
  TH2* hexp2_down = (TH2*)hexp_down->Clone("hexp2_down");
  TH2* hobs2 = (TH2*)hobs->Clone("hobs2");
  TH2* hobu2 = (TH2*)hobu->Clone("hobu2");
  TH2* hobd2 = (TH2*)hobd->Clone("hobd2");


  //////////
  hexp2->SetContour(2);
  hexp2->SetContourLevel(1,1);

  hexp2_up->SetContour(2);
  hexp2_up->SetContourLevel(1,1);

  hexp2_down->SetContour(2);
  hexp2_down->SetContourLevel(1,1);

  hobs2->SetContour(2);
  hobs2->SetContourLevel(1,1);

  hobu2->SetContour(2);
  hobu2->SetContourLevel(1,1);

  hobd2->SetContour(2);
  hobd2->SetContourLevel(1,1);
  
  // All the plotting and cosmetics
  TCanvas* canvas = new TCanvas("canvas", "canvas",625,600);
  canvas->SetLogz();
  
  TH1* frame = canvas->DrawFrame(minX,minY,maxX,maxY, "");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitle("m_{#chi} [GeV]");
  frame->GetXaxis()->SetTitle("m_{#phi} [GeV]");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.20);
  frame->Draw();

  hexp2->SetLineColor(kBlack);
  hexp2_up->SetLineColor(kBlack);
  hexp2_down->SetLineColor(kBlack);    
  hexp2->SetLineWidth(3);
  hexp2_up->SetLineStyle(1);
  hexp2_up->SetLineWidth(1);
  hexp2_down->SetLineStyle(1);
  hexp2_down->SetLineWidth(1);

  hobs2->SetLineColor(kRed);
  hobs2->SetLineWidth(3);
  hobu2->SetLineWidth(1);
  hobd2->SetLineWidth(1);
  hobu2->SetLineColor(kRed);
  hobd2->SetLineColor(kRed);
  
  hobs->SetMinimum(minZ);
  hobs->SetMaximum(maxZ);

  hobs->Draw("COLZ SAME");

  if(addRelicDensity){
    wm->SetFillColor(kBlue);
    wm->SetFillStyle(3005);
    wm->Draw("SAME");
  }

  hexp2_up->Draw("CONT3 SAME");
  hexp2_down->Draw("CONT3 SAME");
  hexp2->Draw("CONT3 SAME");
  hobs2->Draw("CONT3 SAME");
  hobu2->Draw("CONT3 SAME");
  hobd2->Draw("CONT3 SAME");

  CMS_lumi(canvas,"35.9",false,true,false,0,-0.09);

  TLegend *leg = new TLegend(0.1784058,0.5804196,0.5136876,0.8496503,NULL,"brNDC");  //yg
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->AddEntry(hexp2,"Median expected 95% CL","L");
  leg->AddEntry(hexp2_up,"Expected #pm 1 s.d._{experiment}","L");
  leg->AddEntry(hobs2,"Observed 95% CL","L");
  leg->AddEntry(hobu2,"Observed #pm 1 s.d._{theory}","L");
  if(addRelicDensity)
    leg->AddEntry(wm   ,"Thermal Relic","F");
  leg->Draw("SAME");
  
  TLatex *   tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetLineWidth(2);
  tex2->SetTextSize(0.042);
  tex2->SetTextAngle(90);
  tex2->DrawLatex(0.98,0.55,"Observed #sigma_{95% CL}/#sigma_{th}");

  
  gPad->SetRightMargin(0.15);
  gPad->SetLeftMargin(0.13);
  gPad->RedrawAxis();
  gPad->Modified(); 
  gPad->Update();


  canvas->SaveAs((outputDIR+"/scan_fermiportal.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/scan_fermiportal.png").c_str(),"png");

  if(saveOutputFile){
    TFile* outputFile = new TFile((outputDIR+"/fullLikelihood_scan_fermiportal.root").c_str(),"RECREATE");
    outputFile->cd();
    hexp->Write("scan_expected");
    hobs->Write("scan_observed");
    hexp->Write("contour_expected");
    hobs->Write("contour_observed");
    grexp->Write("graph_expected");
    grexp_up->Write("graph_expected_p1s");
    grexp_down->Write("graph_expected_m1s");
    grobs->Write("graph_observed");

    outputFile->Write();
  }
}

/// Relic density
TGraph*relic_gf(){
  
  double *x = new double[1000];
  double *y = new double[1000];

  x[0]=100;
  x[1]=150;
  x[2]=200;
  x[3]=300;
  x[4]=400;
  x[5]=500;
  x[6]=600;
  x[7]=700;
  x[8]=800;
  x[9]=900;
  x[10]=1000;
  x[11]=1025;
  x[12]=1062.966;
  x[13]=1120.97;

  y[0]=1;
  y[1]=5;
  y[2]=10;
  y[3]=25;
  y[4]=45;
  y[5]=75;
  y[6]=110;
  y[7]=160;
  y[8]=220;
  y[9]=300;
  y[10]=410;
  y[11]=440;
  y[12]=490.889;
  y[13]=563.18;
  
	
  TGraph *lrelic = new TGraph(13,x,y);
  lrelic->SetLineColor(kBlue+2);
  lrelic->SetLineWidth(-802);
  lrelic->SetFillStyle(3005);
  lrelic->SetFillColor(kBlue+2);
  return lrelic;
}	

