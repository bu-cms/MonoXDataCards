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
static bool saveOutputFile   = true;
static bool addRelicDensity  = true;
static bool addICHEPContours = false;
//static float nbinsX = 1000;
//static float nbinsY = 600;
static float nbinsX = 425;
static float nbinsY = 250;
static float minX = 0;
static float minY = 1;
static float maxX = 2500;
static float maxY = 1200;
static float minZ = 0.01;
static float maxZ = 10;
static int   reductionForContour = 20;
static bool  addPreliminary = false;
static bool  whiteOut       = true;
static bool  skipPoints     = true;

TGraph* relic_g1();

///////////////
void plotVector(string inputFileName, string outputDIR, bool isDMF = false, string coupling = "025", string energy = "13") {

  if(not isDMF and addICHEPContours) addICHEPContours = false;

  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  if(isDMF) minY = 5;

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

  vector<TGraph*> relicDensity;

  if(addRelicDensity){
    
    if (coupling == "1")
      relicDensity.push_back(relic_g1());
    else if(coupling == "025"){
      TFile* inputFile1 =  TFile::Open("externalFiles/relic_V1.root","READ");
      TList* list = (TList*) inputFile1->Get("mytlist");
      TIter next(list);
      TObject* object = 0;
      while ((object = next())){
	relicDensity.push_back(dynamic_cast<TGraph*>(object));
      }
    }
  }

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
  
  // identify bad limits
  vector<pair<int,int> > goodMassPoint;
  int currentmedmass = -1;
  int currentdmmass  = -1;
  int npoints = 0;

  for(int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);

    int c       = code(mh);
    int medmass = mmed(mh, c);
    int dmmass  = mdm(mh, c);

    if(medmass != currentmedmass or dmmass != currentdmmass){
      if(npoints == 6)
	goodMassPoint.push_back(pair<int,int>(currentmedmass,currentdmmass));
      npoints = 0;
      currentmedmass = medmass;
      currentdmmass  = dmmass;
      npoints++;
    }
    else
      npoints++;
  }

  if(npoints == 6)
    goodMassPoint.push_back(pair<int,int>(currentmedmass,currentdmmass));

  // main loop  
  int expcounter       = 0;
  int exp_up_counter   = 0;
  int exp_down_counter = 0;
  int obscounter       = 0;
  double minmass_exp   = 100000;
  double minmass_obs   = 100000;

  for(int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);
    
    int c       = code(mh);
    int medmass = mmed(mh, c);
    int dmmass  = mdm(mh, c);

    bool isGoodMassPoint = false;
    for(auto mass : goodMassPoint){
      if(medmass == mass.first and dmmass == mass.second){
	isGoodMassPoint = true;
	break;
      }
    }
    if(not isGoodMassPoint){
      cout<<"Bad limit value: medmass "<<medmass<<" dmmass "<<dmmass<<endl;
      continue;
    }

    // remove some point by hand
    if(skipPoints){
      if(not isDMF and coupling == "025"){
	if(medmass == 1800 and (dmmass == 200 or dmmass == 250 or dmmass == 350 or dmmass == 400 or dmmass == 800)) continue;
	if(medmass == 1725 and (dmmass == 200 or dmmass > 700)) continue;
	if(medmass == 1525 and dmmass > 600) continue;
	if(medmass == 1125 and dmmass > 600) continue;
	if(medmass == 600  and dmmass == 350) continue;
	if(medmass == 525  and dmmass == 275) continue;
	if(medmass <= 400  and dmmass >= 400) continue;
	if(quantile == -1 and medmass == 1925 and dmmass == 200) continue;
	if(quantile == -1 and medmass == 1925 and dmmass == 250) continue;
	if(quantile == -1 and medmass == 1800 and dmmass == 300) continue;
	if(quantile == -1 and medmass == 2000 and dmmass == 200) continue;
	
      }
      else if(not isDMF and coupling == "1"){
	if(medmass == 1925  and dmmass == 1000) continue;
	if(medmass == 1800  and dmmass == 1000) continue;
	if(medmass == 1750  and dmmass == 900) continue;
	if(medmass == 1725  and dmmass == 900) continue;
	if(medmass == 1600  and dmmass == 800) continue;
	if(medmass == 1500  and dmmass == 800) continue;
	if(medmass == 1325  and dmmass == 600) continue;
	if(medmass == 1200  and dmmass == 600) continue;
	if(medmass == 1000  and dmmass == 550) continue;
	if(medmass == 925   and dmmass == 500) continue;
      }
      else{
	if(medmass == 1400 and dmmass == 400) continue;
	if(medmass == 1400 and dmmass == 500) continue;
      }      
    }


    if (quantile == 0.5) { // expected limit
      grexp->SetPoint(expcounter, double(medmass), double(dmmass), limit);
      expcounter++;
      if(medmass <= minmass_exp) minmass_exp = medmass;
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
      grobu->SetPoint(obscounter, double(medmass), double(dmmass), limit*0.8);
      grobd->SetPoint(obscounter, double(medmass), double(dmmass), limit*1.2);
      grobs->SetPoint(obscounter, double(medmass), double(dmmass), limit);
      obscounter++;      
      if(medmass <= minmass_obs) minmass_obs = medmass;
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
  
  // fix mass points below min med mass generated                                                                                                                                                 
  for (int i = 1; i   <= nbinsX; i++) {
    for (int j = 1; j <= nbinsY; j++) {
      if(hexp->GetXaxis()->GetBinCenter(i) < minmass_exp){
        hexp_up->SetBinContent(i,j,hexp_up->GetBinContent(hexp_up->FindBin(minmass_exp,hexp_up->GetYaxis()->GetBinCenter(j))));
        hexp_down->SetBinContent(i,j,hexp_down->GetBinContent(hexp_down->FindBin(minmass_exp,hexp_down->GetYaxis()->GetBinCenter(j))));
        hexp->SetBinContent(i,j,hexp->GetBinContent(hexp->FindBin(minmass_exp,hexp->GetYaxis()->GetBinCenter(j))));
      }
      if(hexp->GetXaxis()->GetBinCenter(i) < minmass_obs){
        hobs->SetBinContent(i,j,hobs->GetBinContent(hobs->FindBin(minmass_obs,hobs->GetYaxis()->GetBinCenter(j))));
        hobd->SetBinContent(i,j,hobd->GetBinContent(hobd->FindBin(minmass_obs,hobd->GetYaxis()->GetBinCenter(j))));
        hobu->SetBinContent(i,j,hobu->GetBinContent(hobu->FindBin(minmass_obs,hobu->GetYaxis()->GetBinCenter(j))));
      }
    }
  }

  ////////////
  for(int i = 1; i <= nbinsX; i++){
    for(int j = 1; j <= nbinsY; j++){

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
  double contours[1]; contours[0]=1;
  hexp2->SetContour(1,contours);
  hexp2_up->SetContour(1,contours);
  hexp2_down->SetContour(1,contours);
  hobs2->SetContour(1,contours);
  hobu2->SetContour(1,contours);
  hobd2->SetContour(1,contours);
  
  // All the plotting and cosmetics
  TCanvas* canvas = new TCanvas("canvas", "canvas",625,600);
  canvas->SetRightMargin(0.15);
  canvas->SetLeftMargin(0.14);
  canvas->SetLogz();
  
  TH1* frame = canvas->DrawFrame(minX,minY,maxX,maxY, "");
  frame->GetYaxis()->SetTitle("m_{DM} [GeV]");
  frame->GetXaxis()->SetTitle("m_{med} [GeV]");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.40);
  frame->Draw();

  hobs->GetZaxis()->SetRangeUser(minZ,maxZ);
  hexp->GetZaxis()->SetRangeUser(minZ,maxZ);
  hobu->GetZaxis()->SetRangeUser(minZ,maxZ);
  hobd->GetZaxis()->SetRangeUser(minZ,maxZ);
  hexp_up->GetZaxis()->SetRangeUser(minZ,maxZ);
  hexp_down->GetZaxis()->SetRangeUser(minZ,maxZ);

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

  /// to white out
  if(whiteOut){
    for(int i = 1; i <= nbinsX; i++){
      for(int j = 1; j <= nbinsY; j++){
	if(hexp -> GetXaxis()->GetBinCenter(i) < 2*hexp -> GetYaxis()->GetBinCenter(j) and hexp -> GetBinContent(i,j) == maxZ) hexp->SetBinContent(i,j,minZ*0.01);
	if(hexp_up -> GetXaxis()->GetBinCenter(i) < 2*hexp_up -> GetYaxis()->GetBinCenter(j) and hexp_up -> GetBinContent(i,j) == maxZ) hexp_up->SetBinContent(i,j,minZ*0.01);
	if(hexp_down -> GetXaxis()->GetBinCenter(i) < 2*hexp_down -> GetYaxis()->GetBinCenter(j) and hexp_down -> GetBinContent(i,j) == maxZ) hexp_down->SetBinContent(i,j,minZ*0.01);
	if(hobs -> GetXaxis()->GetBinCenter(i) < 2*hobs -> GetYaxis()->GetBinCenter(j) and hobs -> GetBinContent(i,j) == maxZ) hobs->SetBinContent(i,j,minZ*0.01);
	if(hobu -> GetXaxis()->GetBinCenter(i) < 2*hobu -> GetYaxis()->GetBinCenter(j) and hobu -> GetBinContent(i,j) == maxZ) hobu->SetBinContent(i,j,minZ*0.01);
	if(hobd -> GetXaxis()->GetBinCenter(i) < 2*hobd -> GetYaxis()->GetBinCenter(j) and hobd -> GetBinContent(i,j) == maxZ) hobd->SetBinContent(i,j,minZ*0.01);
      } 
    }

    // more fix by hand
    if(not isDMF){
      for(int i = 1; i <= nbinsX; i++){
	for(int j = 1; j <= nbinsY; j++){
	  
	  // low mass part
	  if(hexp -> GetXaxis()->GetBinCenter(i) < 250 and hexp -> GetYaxis()->GetBinCenter(j) > 220) hexp->SetBinContent(i,j,minZ*0.01);
	  if(hexp_up -> GetXaxis()->GetBinCenter(i) < 250 and hexp_up -> GetYaxis()->GetBinCenter(j) > 220) hexp_up->SetBinContent(i,j,minZ*0.01);
	  if(hexp_down -> GetXaxis()->GetBinCenter(i) < 250 and hexp_down -> GetYaxis()->GetBinCenter(j) > 220) hexp_down->SetBinContent(i,j,minZ*0.01);
	  if(hobs -> GetXaxis()->GetBinCenter(i) < 250 and hobs -> GetYaxis()->GetBinCenter(j) > 220) hobs->SetBinContent(i,j,minZ*0.01);
	  if(hobu -> GetXaxis()->GetBinCenter(i) < 250 and hobu -> GetYaxis()->GetBinCenter(j) > 220) hobu->SetBinContent(i,j,minZ*0.01);
	  if(hobd -> GetXaxis()->GetBinCenter(i) < 250 and hobd -> GetYaxis()->GetBinCenter(j) > 220) hobd->SetBinContent(i,j,minZ*0.01);

	  // low mass part
	  if(hexp -> GetXaxis()->GetBinCenter(i) > 250 and hexp -> GetYaxis()->GetBinCenter(j) > hexp -> GetXaxis()->GetBinCenter(i)/2 + 100) hexp->SetBinContent(i,j,minZ*0.01);
	  if(hobs -> GetXaxis()->GetBinCenter(i) > 250 and hobs -> GetYaxis()->GetBinCenter(j) > hobs -> GetXaxis()->GetBinCenter(i)/2 + 100) hobs->SetBinContent(i,j,minZ*0.01);
	  if(hexp_up -> GetXaxis()->GetBinCenter(i) > 250 and hexp_up -> GetYaxis()->GetBinCenter(j) > hexp_up -> GetXaxis()->GetBinCenter(i)/2+100) hexp_up->SetBinContent(i,j,minZ*0.01);
	  if(hexp_down -> GetXaxis()->GetBinCenter(i) > 250 and hexp_down-> GetYaxis()->GetBinCenter(j) > hexp_down -> GetXaxis()->GetBinCenter(i)/2+100) hexp_down->SetBinContent(i,j,minZ*0.01);
	  if(hobu -> GetXaxis()->GetBinCenter(i) > 250 and hobu -> GetYaxis()->GetBinCenter(j) > hobu -> GetXaxis()->GetBinCenter(i)/2+100) hobu->SetBinContent(i,j,minZ*0.01);
	  if(hobd -> GetXaxis()->GetBinCenter(i) > 250 and hobd -> GetYaxis()->GetBinCenter(j) > hobd -> GetXaxis()->GetBinCenter(i)/2+100) hobd->SetBinContent(i,j,minZ*0.01);    
	} 
      }
    }
    if(isDMF){
      for(int i = 1; i <= nbinsX; i++){
	for(int j = 1; j <= nbinsY; j++){
	  
	  // low mass part
	  if(hexp -> GetXaxis()->GetBinCenter(i) < 250 and hexp -> GetYaxis()->GetBinCenter(j) > 350) hexp->SetBinContent(i,j,minZ*0.01);
	  if(hexp_up -> GetXaxis()->GetBinCenter(i) < 250 and hexp_up -> GetYaxis()->GetBinCenter(j) > 350) hexp_up->SetBinContent(i,j,minZ*0.01);
	  if(hexp_down -> GetXaxis()->GetBinCenter(i) < 250 and hexp_down -> GetYaxis()->GetBinCenter(j) > 350) hexp_down->SetBinContent(i,j,minZ*0.01);
	  if(hobs -> GetXaxis()->GetBinCenter(i) < 250 and hobs -> GetYaxis()->GetBinCenter(j) > 350) hobs->SetBinContent(i,j,minZ*0.01);
	    if(hobu -> GetXaxis()->GetBinCenter(i) < 250 and hobu -> GetYaxis()->GetBinCenter(j) > 350) hobu->SetBinContent(i,j,minZ*0.01);
	    if(hobd -> GetXaxis()->GetBinCenter(i) < 250 and hobd -> GetYaxis()->GetBinCenter(j) > 350) hobd->SetBinContent(i,j,minZ*0.01);
	}
      }      
    }    
  }

  frame->Draw();
  hobs->Draw("COLZ SAME");

  if(addRelicDensity){
    for(auto graph : relicDensity){
      graph->SetLineColor(kBlue+2);
      if(coupling == "025")
	graph->SetLineWidth(802);
      else
	graph->SetLineWidth(-802);
      graph->SetFillStyle(3005);
      graph->SetFillColor(kBlue+2);
      graph->Draw("SAME");
    }
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

  TGraph* graph_obs_ichep = NULL;
  TGraph* graph_exp_ichep = NULL;

  if(addICHEPContours){
    TFile* ichepobs = TFile::Open("externalFiles/monojet_V_MM_ICHEP2016_obs.root","READ");
    graph_obs_ichep = (TGraph*) ichepobs->Get("monojet_obs");
    graph_obs_ichep->SetLineWidth(3);
    graph_obs_ichep->SetLineStyle(1);
    graph_obs_ichep->SetLineColor(kViolet);
    graph_obs_ichep->Draw("Lsame");
    TFile* ichepexp = TFile::Open("externalFiles/monojet_V_MM_ICHEP2016_exp.root","READ");
    graph_exp_ichep = (TGraph*) ichepexp->Get("monojet_exp");
    graph_exp_ichep->SetLineWidth(3);
    graph_exp_ichep->SetLineStyle(7);
    graph_exp_ichep->SetLineColor(kViolet);
    graph_exp_ichep->Draw("Lsame");
  }

  contour_exp_up->Draw("Lsame");
  contour_exp_dw->Draw("Lsame");
  contour_exp->Draw("Lsame");

  contour_obs_up->SetLineColor(kRed);
  contour_obs->SetLineColor(kRed);
  contour_obs_dw->SetLineColor(kRed);
  contour_obs_up->SetLineWidth(1);
  contour_obs->SetLineWidth(3);
  contour_obs_dw->SetLineWidth(1);

  contour_obs_up->Draw("Lsame");
  contour_obs_dw->Draw("Lsame");
  contour_obs->Draw("Lsame");


  if(not addPreliminary)
    CMS_lumi(canvas,"35.9",false,true,false,0,-0.09);
  else
    CMS_lumi(canvas,"35.9",false,false,false,0,-0.09);

  TLegend *leg = NULL;
  if(not addICHEPContours)
    leg = new TLegend(0.165,0.48,0.53,0.77);
  else
    leg = new TLegend(0.165,0.45,0.53,0.78);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->AddEntry(contour_exp,"Median expected 95% CL","L");
  leg->AddEntry(contour_exp_up,"#pm 1 #sigma_{experiment}","L");
  leg->AddEntry(contour_obs,"Observed 95% CL","L");
  leg->AddEntry(contour_obs_up,"Observed #pm theory unc.","L");
  if(addICHEPContours){
    leg->AddEntry(graph_obs_ichep,"1703.01651 observed","L");
    leg->AddEntry(graph_exp_ichep,"1703.01651 expected","L");
  }
  if(addRelicDensity)
    leg->AddEntry(relicDensity.front(),"#Omega_{c}#timesh^{2} #geq 0.12","F");
   

  leg->Draw("SAME");

  TLatex * tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.030);
  tex->Draw();
  if (coupling == "1"){
    tex->DrawLatex(0.175,0.80,"#bf{Vector med, Dirac DM, g_{q} = 1, g_{DM} = 1}");
  }
  else
    tex->DrawLatex(0.175,0.80,"#bf{Vector med, Dirac DM, g_{q} = 0.25, g_{DM} = 1}");
  
  
  TLatex *   tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetLineWidth(2);
  tex2->SetTextSize(0.042);
  tex2->SetTextAngle(90);
  tex2->DrawLatex(0.975,0.55,"Observed #sigma_{95% CL}/#sigma_{th}");
  
  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR+"/scan_vector_g"+string(coupling)+"_"+string(energy)+"TeV_v2.pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_vector_g"+string(coupling)+"_"+string(energy)+"TeV_v2.png").c_str());

  if(saveOutputFile){
    TFile* outputFile = new TFile((outputDIR+"/scan_vector_g"+string(coupling)+"_"+string(energy)+"TeV_v2.root").c_str(),"RECREATE");
    outputFile->cd();
    hobs->Write("Observed_limit");
    hexp->Write("Expected_limit");
    hexp_up->Write("Expected_limit_1sUp");
    hexp_down->Write("Expected_limit_1sDw");
    outputFile->Write();    
  }
}

TGraph*relic_g1(){


  double *x = new double[1000];
  double *y = new double[1000];
  
  
  x[0]=383.945; y[0]=6.075;
  x[1]=386.039; y[1]=8.225;
  x[2]=388.133; y[2]=10.375;
  x[3]=390.227; y[3]=12.525;
  x[4]=392.321; y[4]=14.675;
  x[5]=394.415; y[5]=16.825;
  x[6]=396.509; y[6]=18.975;
  x[7]=398.603; y[7]=21.125;
  x[8]=398.674; y[8]=21.1458;
  x[9]=422.23; y[9]=22.7553;
  x[10]=445.786; y[10]=22.8104;
  x[11]=469.342; y[11]=22.8655;
  x[12]=492.898; y[12]=22.9205;
  x[13]=516.454; y[13]=22.9756;
  x[14]=540.01; y[14]=23.0307;
  x[15]=563.566; y[15]=23.0858;
  x[16]=587.122; y[16]=23.1409;
  x[17]=610.678; y[17]=23.196;
  x[18]=610.984; y[18]=23.275;
  x[19]=619.306; y[19]=25.425;
  x[20]=627.628; y[20]=27.575;
  x[21]=634.234; y[21]=28.0937;
  x[22]=657.79; y[22]=28.5847;
  x[23]=681.346; y[23]=28.6326;
  x[24]=704.902; y[24]=28.6805;
  x[25]=728.458; y[25]=28.7284;
  x[26]=752.014; y[26]=28.7763;
  x[27]=775.57; y[27]=28.8242;
  x[28]=799.126; y[28]=28.8721;
  x[29]=822.682; y[29]=28.92;
  x[30]=846.238; y[30]=28.9679;
  x[31]=869.794; y[31]=29.0158;
  x[32]=893.35; y[32]=29.0637;
  x[33]=916.906; y[33]=29.2578;
  x[34]=938.804; y[34]=29.725;
  x[35]=940.462; y[35]=29.8154;
  x[36]=964.018; y[36]=30.1224;
  x[37]=987.574; y[37]=30.4293;
  x[38]=1011.13; y[38]=30.7363;
  x[39]=1034.69; y[39]=31.0433;
  x[40]=1058.24; y[40]=31.3502;
  x[41]=1081.8; y[41]=31.6572;
  x[42]=1098.51; y[42]=31.875;
  x[43]=1105.35; y[43]=33.4812;
  x[44]=1107.67; y[44]=34.025;
  x[45]=1116.83; y[45]=36.175;
  x[46]=1128.91; y[46]=37.6971;
  x[47]=1144.92; y[47]=38.325;
  x[48]=1152.47; y[48]=38.5117;
  x[49]=1176.02; y[49]=39.0946;
  x[50]=1199.58; y[50]=39.6774;
  x[51]=1214.37; y[51]=40.475;
  x[52]=1223.13; y[52]=40.9475;
  x[53]=1233.71; y[53]=42.625;
  x[54]=1246.69; y[54]=44.6836;
  x[55]=1247.27; y[55]=44.775;
  x[56]=1260.82; y[56]=46.925;
  x[57]=1270.25; y[57]=48.4197;
  x[58]=1274.96; y[58]=49.075;
    x[59]=1293.8; y[59]=50.2066;
  x[60]=1317.36; y[60]=50.7122;
  x[61]=1340.91; y[61]=51.2177;
  x[62]=1341.25; y[62]=51.225;
  x[63]=1364.47; y[63]=51.7233;
  x[64]=1388.03; y[64]=52.2289;
  x[65]=1411.58; y[65]=52.7344;
  x[66]=1435.14; y[66]=53.24;
  x[67]=1441.43; y[67]=53.375;
  x[68]=1458.69; y[68]=53.8498;
  x[69]=1482.25; y[69]=54.4976;
  x[70]=1505.81; y[70]=55.1455;
  x[71]=1517.72; y[71]=55.525;
  x[72]=1529.36; y[72]=56.7616;
  x[73]=1535.06; y[73]=57.675;
  x[74]=1548.49; y[74]=59.825;
  x[75]=1552.92; y[75]=60.5346;
  x[76]=1561.91; y[76]=61.975;
  x[77]=1575.33; y[77]=64.125;
  x[78]=1576.47; y[78]=64.3076;
  x[79]=1588.76; y[79]=66.275;
  x[80]=1600.03; y[80]=68.0806;
  x[81]=1602.18; y[81]=68.425;
  x[82]=1615.6; y[82]=70.575;
  x[83]=1623.59; y[83]=71.3745;
  x[84]=1647.14; y[84]=72.5647;
  x[85]=1651.23; y[85]=72.725;
  x[86]=1670.7; y[86]=73.4892;
  x[87]=1694.25; y[87]=74.4136;
  x[88]=1706.01; y[88]=74.875;
  x[89]=1717.81; y[89]=75.338;
  x[90]=1741.37; y[90]=76.2625;
  x[91]=1760.8; y[91]=77.025;
  x[92]=1764.92; y[92]=77.1869;
  x[93]=1788.48; y[93]=78.1114;
  x[94]=1812.03; y[94]=79.0358;
  x[95]=1813.77; y[95]=79.175;
  x[96]=1835.59; y[96]=81.0713;
  x[97]=1839.33; y[97]=81.325;
  x[98]=1852.62; y[98]=83.475;
  x[99]=1859.15; y[99]=84.5309;
  x[100]=1865.91; y[100]=85.625;
  x[101]=1879.19; y[101]=87.775;
  x[102]=1882.7; y[102]=88.3427;
  x[103]=1892.48; y[103]=89.925;
  x[104]=1905.77; y[104]=92.075;
  x[105]=1906.26; y[105]=92.1546;
  x[106]=1919.05; y[106]=94.225;
  x[107]=1929.81; y[107]=95.7929;
  x[108]=1936.54; y[108]=96.375;
  x[109]=1953.37; y[109]=97.2775;
  x[110]=1976.64; y[110]=98.525;
  x[111]=1976.93; y[111]=98.5404;
  x[112]=2000.48; y[112]=99.8034;
  x[113]=2016.74; y[113]=100.675;
  x[114]=2024.04; y[114]=101.066;
  x[115]=2047.59; y[115]=102.329;
  x[116]=2056.84; y[116]=102.825;
  x[117]=2071.15; y[117]=103.592;
  x[118]=2094.71; y[118]=104.855;
  x[119]=2096.94; y[119]=104.975;
  x[120]=2113.7; y[120]=107.125;
  x[121]=2118.26; y[121]=107.481;
  x[122]=2131.02; y[122]=109.275;
  x[123]=2141.82; y[123]=110.794;
  x[124]=2145.69; y[124]=111.425;
  x[125]=2158.91; y[125]=113.575;
  x[126]=2165.37; y[126]=114.627;
  x[127]=2172.12; y[127]=115.725;
  x[128]=2185.34; y[128]=117.875;
  x[129]=2188.93; y[129]=118.46;
  x[130]=2198.55; y[130]=120.025;
  x[131]=2211.76; y[131]=122.175;
  x[132]=2212.49; y[132]=122.253;
  x[133]=2236.04; y[133]=123.98;
  x[134]=2241.33; y[134]=124.325;
  x[135]=2259.6; y[135]=125.516;
  x[136]=2274.31; y[136]=126.475;
  x[137]=2283.15; y[137]=127.052;
  x[138]=2306.71; y[138]=128.587;
  x[139]=2307.29; y[139]=128.625;
  x[140]=2330.27; y[140]=130.589;
  x[141]=2332.45; y[141]=130.775;
  x[142]=2343.73; y[142]=132.925;
  x[143]=2353.82; y[143]=134.847;
  x[144]=2355.24; y[144]=135.075;
  x[145]=2377.37; y[145]=137.225;
  x[146]=2377.38; y[146]=137.225;
  x[147]=2400.93; y[147]=139.156;
  x[148]=2403.61; y[148]=139.375;
  x[149]=2424.49; y[149]=141.086;
  x[150]=2427.27; y[150]=141.525;
  x[151]=2440.86; y[151]=143.675;
  x[152]=2448.05; y[152]=144.811;
  x[153]=2454.25; y[153]=145.825;
  x[154]=2467.41; y[154]=147.975;
  x[155]=2471.6; y[155]=148.659;
  x[156]=2481.85; y[156]=150.125;
  x[157]=2495.16; y[157]=151.343;
  x[158]=2507.65; y[158]=152.275;
  x[159]=2518.71; y[159]=153.101;
  x[160]=2536.45; y[160]=154.425;
  x[161]=2542.27; y[161]=155.202;
  x[162]=2552.57; y[162]=156.575;
  x[163]=2564.06; y[163]=158.725;
  x[164]=2565.83; y[164]=159.055;
  x[165]=2575.55; y[165]=160.875;
  x[166]=2587.04; y[166]=163.025;
  x[167]=2589.38; y[167]=163.462;
  x[168]=2599.37; y[168]=165.175;
  x[169]=2612.94; y[169]=166.664;
  x[170]=2620.29; y[170]=167.325;
  x[171]=2636.49; y[171]=168.78;
  x[172]=2644.23; y[172]=169.475;
  x[173]=2660.05; y[173]=170.897;
  x[174]=2668.16; y[174]=171.625;
  x[175]=2683.61; y[175]=173.013;
  x[176]=2692.09; y[176]=173.775;
  x[177]=2707.16; y[177]=175.129;
  x[178]=2715.25; y[178]=175.925;
  x[179]=2730.72; y[179]=177.446;
  x[180]=2734.68; y[180]=178.075;
  x[181]=2753.59; y[181]=180.225;
  x[182]=2754.27; y[182]=180.35;
  x[183]=2765.32; y[183]=182.375;
  x[184]=2777.04; y[184]=184.525;
  x[185]=2777.83; y[185]=184.67;
  x[186]=2788.76; y[186]=186.675;
  x[187]=2800.49; y[187]=188.825;
  x[188]=2801.39; y[188]=188.99;
  x[189]=2812.21; y[189]=190.975;
  x[190]=2823.94; y[190]=193.125;
  x[191]=2824.94; y[191]=193.309;
  x[192]=2835.98; y[192]=195.275;
  x[193]=2848.5; y[193]=196.768;
  x[194]=2855.38; y[194]=197.425;
  x[195]=2872.05; y[195]=199.016;
  x[196]=2877.92; y[196]=199.575;
  x[197]=2895.61; y[197]=201.263;
  x[198]=2900.45; y[198]=201.725;
  x[199]=2919.17; y[199]=203.51;
  x[200]=2922.99; y[200]=203.875;
  x[201]=2939.77; y[201]=206.025;
  x[202]=2942.72; y[202]=206.63;
  x[203]=2950.24; y[203]=208.175;
  x[204]=2960.81; y[204]=210.325;
  x[205]=2966.28; y[205]=211.074;
  x[206]=2978.98; y[206]=212.475;
  x[207]=2989.83; y[207]=213.671;
  x[208]=2998.48; y[208]=214.625;
  x[209]=3013.39; y[209]=216.269;
  x[210]=3017.23; y[210]=216.775;
  x[211]=3033.57; y[211]=218.925;
  x[212]=3036.95; y[212]=219.369;
  x[213]=3046.29; y[213]=221.075;
  x[214]=3058.07; y[214]=223.225;
  x[215]=3060.5; y[215]=223.668;
  x[216]=3069.85; y[216]=225.375;
  x[217]=3084.06; y[217]=227.223;
  x[218]=3086.94; y[218]=227.525;
  x[219]=3107.61; y[219]=229.641;
  x[220]=3107.94; y[220]=229.675;
  x[221]=3120.15; y[221]=231.825;
  x[222]=3130.8; y[222]=233.975;
  x[223]=3131.17; y[223]=234.05;
  x[224]=3141.44; y[224]=236.125;
  x[225]=3152.09; y[225]=238.275;
  x[226]=3154.73; y[226]=238.808;
  x[227]=3162.73; y[227]=240.425;
  x[228]=3173.38; y[228]=242.575;
  x[229]=3178.28; y[229]=243.279;
  x[230]=3190.54; y[230]=244.725;
  x[231]=3201.84; y[231]=246.028;
  x[232]=3209.17; y[232]=246.875;
  x[233]=3225.39; y[233]=248.747;
  x[234]=3227.8; y[234]=249.025;
  x[235]=3246.43; y[235]=251.175;
  x[236]=3248.95; y[236]=251.466;
  x[237]=3265.06; y[237]=253.325;
  x[238]=3272.51; y[238]=254.285;
  x[239]=3281.74; y[239]=255.475;
  x[240]=3291.58; y[240]=257.625;
  x[241]=3296.06; y[241]=258.342;
  x[242]=3306.78; y[242]=259.775;
  x[243]=3319.62; y[243]=261.406;
  x[244]=3322.85; y[244]=261.925;
  x[245]=3336.25; y[245]=264.075;
  x[246]=3343.17; y[246]=265.186;
  x[247]=3348.44; y[247]=266.225;
  x[248]=3359.33; y[248]=268.375;
  x[249]=3366.73; y[249]=269.837;
  x[250]=3370.21; y[250]=270.525;
  x[251]=3381.1; y[251]=272.675;
  x[252]=3390.29; y[252]=274.489;
  x[253]=3392.26; y[253]=274.825;
  x[254]=3409.7; y[254]=276.975;
  x[255]=3413.84; y[255]=277.463;
  x[256]=3427.92; y[256]=279.125;
  x[257]=3437.4; y[257]=280.645;
  x[258]=3441.32; y[258]=281.275;
  x[259]=3451.42; y[259]=283.425;
  x[260]=3460.95; y[260]=285.456;
  x[261]=3461.51; y[261]=285.575;
  x[262]=3471.61; y[262]=287.725;
  x[263]=3481.7; y[263]=289.875;
  x[264]=3484.51; y[264]=290.474;
  x[265]=3492.83; y[265]=292.025;
  x[266]=3508.07; y[266]=294.068;
  x[267]=3508.89; y[267]=294.175;
  x[268]=3525.44; y[268]=296.325;
  x[269]=3531.62; y[269]=297.128;
  x[270]=3541.99; y[270]=298.475;
  x[271]=3555.18; y[271]=300.188;
  x[272]=3558.54; y[272]=300.625;
  x[273]=3575.1; y[273]=302.775;
  x[274]=3578.73; y[274]=303.248;
  x[275]=3591.65; y[275]=304.925;
  x[276]=3601.29; y[276]=307.075;
  x[277]=3602.29; y[277]=307.253;
  x[278]=3615.26; y[278]=309.225;
  x[279]=3625.85; y[279]=310.718;
  x[280]=3629.39; y[280]=311.375;
  x[281]=3640.99; y[281]=313.525;
  x[282]=3649.4; y[282]=315.084;
  x[283]=3652.23; y[283]=315.675;
  x[284]=3662.53; y[284]=317.825;
  x[285]=3672.83; y[285]=319.975;
  x[286]=3672.96; y[286]=320.002;
  x[287]=3683.13; y[287]=322.125;
  x[288]=3693.43; y[288]=324.275;
  x[289]=3696.51; y[289]=324.784;
  x[290]=3708.11; y[290]=326.425;
  x[291]=3720.07; y[291]=328.004;
  x[292]=3724.4; y[292]=328.575;
  x[293]=3738.45; y[293]=330.725;
  x[294]=3743.63; y[294]=331.88;
  x[295]=3748.09; y[295]=332.875;
  x[296]=3757.73; y[296]=335.025;
  x[297]=3767.18; y[297]=337.131;
  x[298]=3767.38; y[298]=337.175;
  x[299]=3777.02; y[299]=339.325;
  x[300]=3786.67; y[300]=341.475;
  x[301]=3790.74; y[301]=342.193;
  x[302]=3800.1; y[302]=343.625;
  x[303]=3814.29; y[303]=345.651;
  x[304]=3815.16; y[304]=345.775;
  x[305]=3830.23; y[305]=347.925;
  x[306]=3837.85; y[306]=349.013;
  x[307]=3845.29; y[307]=350.075;
  x[308]=3860.36; y[308]=352.225;
  x[309]=3861.41; y[309]=352.375;
  x[310]=3875.42; y[310]=354.375;
  x[311]=3884.96; y[311]=356.273;
  x[312]=3886.23; y[312]=356.525;
  x[313]=3896.31; y[313]=358.675;
  x[314]=3908.52; y[314]=360.618;
  x[315]=3909.87; y[315]=360.825;
  x[316]=3923.92; y[316]=362.975;
  x[317]=3932.07; y[317]=364.222;
  x[318]=3936.43; y[318]=365.125;
  x[319]=3946.82; y[319]=367.275;
  x[320]=3955.63; y[320]=369.1;
  x[321]=3957.12; y[321]=369.425;
  x[322]=3966.99; y[322]=371.575;
  x[323]=3976.86; y[323]=373.725;
  x[324]=3979.19; y[324]=374.231;
  x[325]=3987.46; y[325]=375.875;
  x[326]=4001.77; y[326]=378.025;

  TGraph *lrelic = new TGraph(778,x,y);
  lrelic->SetLineColor(kBlue+2);
  lrelic->SetLineWidth(-802);
  lrelic->SetFillStyle(3005);
  lrelic->SetFillColor(kBlue+2);
  return lrelic;
}
