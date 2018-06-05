#include "../CMS_lumi.h"

int mmed(double mh){

  int code = int(mh);
  if (code == 5010) return ((int)(1));
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

TGraph* produceContour (const int & reduction){
  
  TObjArray *lContoursE = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  std::vector<double> lXE;
  std::vector<double> lYE;
  int lTotalContsE = lContoursE->GetSize();
  for(int i0 = 0; i0 < lTotalContsE; i0++){
    TList * pContLevel = (TList*)lContoursE->At(i0);
    TGraph *pCurv = (TGraph*)pContLevel->First();
    for(int i1 = 0; i1 < pContLevel->GetSize(); i1++){
      for(int i2  = 0; i2 < pCurv->GetN(); i2++) {
	if(i2%reduction != 0) continue; // reduce number of points
	lXE.push_back(pCurv->GetX()[i2]); 
	lYE.push_back(pCurv->GetY()[i2]);
      }
      pCurv->SetLineColor(kRed);                                                                                                            
      pCurv = (TGraph*)pContLevel->After(pCurv);                                                                                        
    }
  }
  if(lXE.size() == 0) {
    lXE.push_back(0); 
    lYE.push_back(0); 
  }
  
  TGraph *lTotalE = new TGraph(lXE.size(),&lXE[0],&lYE[0]);  
  return lTotalE;
}

/////////
static bool saveOutputFile  = false;
static bool addRelicDensity = true;
static bool addTheoreticalLine = false;
//static float nbinsX = 800;
//static float nbinsY = 500;
static float nbinsX = 300;
static float nbinsY = 180;
static float minX = 0;
static float minY = 10.;
static float maxX = 1800;
static float maxY = 1000;
static float minZ = 0.01;
static float maxZ = 10;
static int reductionForContour = 20;

TGraph* relic_gf();
TGraph* theoretical_gf();

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
  
  TGraph* wm2 = NULL;
  if(addTheoreticalLine)
    wm2 = theoretical_gf();

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
  canvas->SetRightMargin(0.15);
  canvas->SetLeftMargin(0.13);
  canvas->SetLogz();
  
  TH1* frame = canvas->DrawFrame(minX,minY,maxX,maxY, "");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitle("m_{#chi} [GeV]");
  frame->GetXaxis()->SetTitle("m_{#phi} [GeV]");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.20);
  frame->Draw();

  hobs->SetMinimum(minZ);
  hobs->SetMaximum(maxZ);

  hexp2_up->GetZaxis()->SetLabelSize(0);
  hexp2_up->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp_up = produceContour(reductionForContour);

  hexp2_down->GetZaxis()->SetLabelSize(0);
  hexp2_down->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp_dw = produceContour(reductionForContour);

  hexp2->GetZaxis()->SetLabelSize(0);
  hexp2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp = produceContour(reductionForContour);

  hobs2->GetZaxis()->SetLabelSize(0);
  hobs2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_obs = produceContour(reductionForContour);

  hobu2->GetZaxis()->SetLabelSize(0);
  hobu2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_obs_up = produceContour(reductionForContour);

  hobd2->GetZaxis()->SetLabelSize(0);
  hobd2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_obs_dw = produceContour(reductionForContour);
  
  frame->Draw();
hobs->Draw("COLZ SAME");

  if(addRelicDensity){
    wm->SetFillColor(kBlue);
    wm->SetFillStyle(3005);
    wm->Draw("SAME");
  }

  if(addTheoreticalLine){
    wm2->SetFillColor(kMagenta+3);
    wm2->SetFillStyle(3005);
    wm2->Draw("SAME");
  }
  contour_exp_up->SetLineColor(kBlack);
  contour_exp->SetLineColor(kBlack);
  contour_exp_dw->SetLineColor(kBlack);
  contour_exp_up->SetLineWidth(1);
  contour_exp->SetLineWidth(3);
  contour_exp_dw->SetLineWidth(1);
  contour_exp_up->SetLineStyle(7);
  contour_exp->SetLineStyle(7);
  contour_exp_dw->SetLineStyle(7);

  contour_exp_up->Draw("Lsame");
  contour_exp_dw->Draw("Lsame");
  contour_exp->Draw("Lsame");

  contour_obs_up->SetLineColor(kRed);
  contour_obs->SetLineColor(kRed);
  contour_obs_dw->SetLineColor(kRed);
  contour_obs_up->SetLineWidth(1);
  contour_obs->SetLineWidth(3);
  contour_obs_dw->SetLineWidth(1);

  //contour_obs_up->Draw("Lsame");
  //contour_obs_dw->Draw("Lsame");
 contour_obs->Draw("Lsame");

  CMS_lumi(canvas,"35.9",false,true,false,0,-0.09);

  TLegend *leg = new TLegend(0.1784058,0.5804196,0.5136876,0.8496503,NULL,"brNDC");  //yg
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);

  leg->AddEntry(contour_exp,"Median expected 95% CL","L");
  leg->AddEntry(contour_exp_up,"Expected #pm 1 s.d._{experiment}","L");
  leg->AddEntry(contour_obs,"Observed 95% CL","L");
  //leg->AddEntry(contour_obs_up,"Observed #pm 1 s.d._{theory}","L");

  if(addRelicDensity)
    leg->AddEntry(wm   ,"#Omega_{c}#timesh^{2} #geq 0.12","F");
  if(addTheoreticalLine)
    leg->AddEntry(wm2   ,"Theoretical Expected (8 TeV)","L");
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
    TFile* outputFile = new TFile((outputDIR+"/scan_fermiportal.root").c_str(),"RECREATE");
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
x[0]=50;     y[0]=3.8;
x[1]=139;    y[1]=11;
x[2]=262;    y[2]=23;
x[3]=336;    y[3]=38;
x[4]=394;    y[4]=53;
x[5]=467;    y[5]=77;
x[6]=541;    y[6]=104;
x[7]=603;    y[7]=134;
x[8]=657;    y[8]=161;
x[9]=727;    y[9]=200;
x[10]=781;   y[10]=235;
x[11]=835;   y[11]=273;
x[12]=850;   y[12]=285;
x[13]=889;   y[13]=320;
x[14]=943;   y[14]=366;
x[15]=997;   y[15]=424;
x[16]=1036;  y[16]=462;
x[17]=1067;  y[17]=501;
x[18]=1121;  y[18]=578;
x[19]=1171;  y[19]=663;
x[20]=1214;  y[20]=744;
x[21]=1241;  y[21]=813;
x[22]=1268;  y[22]=871;
x[23]=1287;  y[23]=952;
x[24]=1303;  y[24]=1002;
x[25]=1314;  y[25]=1048;
x[26]=1326;  y[26]=1114;
x[27]=1338;  y[27]=1141;
x[28]=1349;  y[28]=1172;
x[29]=1361;  y[29]=1206;
x[30]=1384;  y[30]=1249;
x[31]=1407;  y[31]=1284;
x[32]=1427;  y[32]=1311;
x[33]=1454;  y[33]=1353;
x[34]=1488;  y[34]=1399;
x[35]=1519;  y[35]=1446;
x[36]=1550;  y[36]=1488;
x[37]=1562;  y[37]=1500;
  
  TGraph *lrelic = new TGraph(37,x,y);
  lrelic->SetLineColor(kBlue+2);
  lrelic->SetLineWidth(-802);
  lrelic->SetFillStyle(3005);
  lrelic->SetFillColor(kBlue+2);
  return lrelic;
}	

 /// Theorical Line
 TGraph*theoretical_gf(){

   double *x = new double[1000];
   double *y = new double[1000];

x[0]=1.00000000000000e+002;	 y[0]=8.23770491803279e+001;
x[1]=1.06603773584906e+002;	 y[1]=8.72950819672131e+001;
x[2]=1.15408805031447e+002;	 y[2]=9.34426229508197e+001;
x[3]=1.26415094339623e+002;	 y[3]=1.09426229508197e+002;
x[4]=1.35220125786164e+002;	 y[4]=1.22950819672131e+002;
x[5]=1.41823899371069e+002;	 y[5]=1.35245901639344e+002;
x[6]=1.50628930817610e+002;	 y[6]=1.48770491803279e+002;
x[7]=1.61635220125786e+002;	 y[7]=1.58606557377049e+002;
x[8]=1.72641509433962e+002;	 y[8]=1.69672131147541e+002;
x[9]=1.88050314465409e+002;	 y[9]=1.75819672131148e+002;
x[10]=1.99056603773585e+002;	 y[10]=1.79508196721311e+002;
x[11]=2.12264150943396e+002;	 y[11]=1.79508196721311e+002;
x[12]=2.25471698113208e+002;	 y[12]=1.78278688524590e+002;
x[13]=2.43081761006289e+002;	 y[13]=1.75819672131148e+002;
x[14]=2.58490566037736e+002;	 y[14]=1.69672131147541e+002;
x[15]=2.67295597484277e+002;	 y[15]=1.63524590163934e+002;
x[16]=2.76100628930818e+002;	 y[16]=1.52459016393443e+002;
x[17]=2.80503144654088e+002;	 y[17]=1.43852459016393e+002;
x[18]=2.82704402515723e+002;	 y[18]=1.34016393442623e+002;
x[19]=2.87106918238994e+002;	 y[19]=1.22950819672131e+002;
x[20]=2.91509433962264e+002;	 y[20]=1.11885245901639e+002;
x[21]=2.95911949685535e+002;	 y[21]=1.05737704918033e+002;
x[22]=3.00314465408805e+002;	 y[22]=1.02049180327869e+002;
x[23]=3.04716981132075e+002;	 y[23]=1.00819672131148e+002;
x[24]=3.09119496855346e+002;	 y[24]=9.95901639344262e+001;
x[25]=3.13522012578616e+002;	 y[25]=1.02049180327869e+002;
x[26]=3.20125786163522e+002;	 y[26]=1.06967213114754e+002;
x[27]=3.28930817610063e+002;	 y[27]=1.18032786885246e+002;
x[28]=3.35534591194969e+002;	 y[28]=1.31557377049180e+002;
x[29]=3.44339622641509e+002;	 y[29]=1.50000000000000e+002;
x[30]=3.53144654088050e+002;	 y[30]=1.62295081967213e+002;
x[31]=3.59748427672956e+002;	 y[31]=1.70901639344262e+002;
x[32]=3.70754716981132e+002;	 y[32]=1.79508196721311e+002;
x[33]=3.83962264150943e+002;	 y[33]=1.85655737704918e+002;
x[34]=4.03773584905660e+002;	 y[34]=1.89344262295082e+002;
x[35]=4.34591194968553e+002;	 y[35]=1.94262295081967e+002;
x[36]=4.63207547169811e+002;	 y[36]=1.99180327868852e+002;
x[37]=4.80817610062893e+002;	 y[37]=2.04098360655738e+002;
x[38]=4.94025157232704e+002;	 y[38]=2.06557377049180e+002;
x[39]=5.00628930817610e+002;	 y[39]=2.06557377049180e+002;
x[40]=5.05031446540880e+002;	 y[40]=2.02868852459016e+002;
x[41]=5.09433962264151e+002;	 y[41]=1.96721311475410e+002;
x[42]=5.11635220125786e+002;	 y[42]=1.90573770491803e+002;
x[43]=5.16037735849057e+002;	 y[43]=1.81967213114754e+002;
x[44]=5.16037735849057e+002;	 y[44]=1.64754098360656e+002;
x[45]=5.20440251572327e+002;	 y[45]=1.53688524590164e+002;
x[46]=5.20440251572327e+002;	 y[46]=1.41393442622951e+002;
x[47]=5.22641509433962e+002;	 y[47]=1.27868852459016e+002;
x[48]=5.24842767295597e+002;	 y[48]=1.15573770491803e+002;
x[49]=5.27044025157233e+002;	 y[49]=1.02049180327869e+002;
x[50]=5.31446540880503e+002;	 y[50]=9.46721311475410e+001;
x[51]=5.35849056603774e+002;	 y[51]=8.48360655737705e+001;
x[52]=5.42452830188679e+002;	 y[52]=7.13114754098361e+001;
x[53]=5.53459119496855e+002;	 y[53]=6.27049180327869e+001;
x[54]=5.64465408805031e+002;	 y[54]=5.53278688524590e+001;
x[55]=5.75471698113208e+002;	 y[55]=5.04098360655738e+001;
x[56]=5.86477987421384e+002;	 y[56]=4.67213114754098e+001;
x[57]=5.97484276729560e+002;	 y[57]=4.42622950819672e+001;
x[58]=6.06289308176101e+002;	 y[58]=4.42622950819672e+001;
x[59]=6.17295597484277e+002;	 y[59]=4.67213114754098e+001;
x[60]=6.34905660377359e+002;	 y[60]=5.04098360655738e+001;
x[61]=6.48113207547170e+002;	 y[61]=5.40983606557377e+001;
x[62]=6.63522012578616e+002;	 y[62]=6.02459016393443e+001;
x[63]=6.78930817610063e+002;	 y[63]=6.63934426229508e+001;
x[64]=6.87735849056604e+002;	 y[64]=7.00819672131148e+001;
x[65]=6.94339622641509e+002;	 y[65]=7.25409836065574e+001;
x[66]=7.00943396226415e+002;	 y[66]=7.62295081967213e+001;
x[67]=7.05345911949686e+002;	 y[67]=7.62295081967213e+001;
x[68]=7.18553459119497e+002;	 y[68]=7.13114754098361e+001;
x[69]=7.31761006289308e+002;	 y[69]=6.63934426229508e+001;
x[70]=7.42767295597484e+002;	 y[70]=6.02459016393443e+001;
x[71]=7.53773584905660e+002;	 y[71]=5.28688524590164e+001;
x[72]=7.62578616352201e+002;	 y[72]=4.54918032786885e+001;
x[73]=7.71383647798742e+002;	 y[73]=3.44262295081967e+001;
x[74]=7.75786163522013e+002;	 y[74]=2.58196721311475e+001;
x[75]=7.80188679245283e+002;	 y[75]=1.72131147540984e+001;
x[76]=7.82389937106918e+002;	 y[76]=8.60655737704918e+000;
x[77]=7.82389937106918e+002;	 y[77]=-1.22950819672131e+000;


  TGraph *ltheoretical = new TGraph(77,x,y);
  ltheoretical->SetLineColor(kMagenta+3);
  ltheoretical->SetLineWidth(-802);
  ltheoretical->SetFillStyle(3005);
  ltheoretical->SetFillColor(kMagenta+3);
  return ltheoretical;
}
