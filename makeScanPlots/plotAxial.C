#include "../CMS_lumi.h"

int mmed(double mh, int code){

    if (code == 800) return ((int)(mh-80000000000))/10000; 
    if (code == 801) return ((int)(mh-80100000000))/10000; 
    if (code == 805) return ((int)(mh-80500000000))/10000; 
    if (code == 806) return ((int)(mh-80600000000))/10000; 
    return -1;
}

// 
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

// to produce the contour as a TGprah reducing points
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
static bool addRelicDensity = true;
static bool saveOutputFile  = true;
static bool addICHEPContours = false;
static float nbinsX = 1000;
static float nbinsY = 600;
static float minX = 0;
static float minY = 4;
static float maxX = 2500;
static float maxY = 1200;
static double minZ = 0.01;
static float maxZ = 10;
static int reductionForContour = 20;
static bool whiteOut = true;
static bool addPreliminary = false;
static bool skipPoints = false;

TGraph* relic_g1();

///////////
void plotAxial(string inputFileName, string outputDIR, bool isDMF = false, string coupling = "025", string energy = "13") {

  if(not isDMF and addICHEPContours) addICHEPContours = false;

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

  vector<TGraph*> relicDensity;

  if(addRelicDensity){
    if (coupling == "1")
      relicDensity.push_back(relic_g1());
    else if(coupling == "025"){
      TFile* inputFile1 =  TFile::Open("externalFiles/relic_A1.root","READ");
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

  // identify bad limits files
  int currentmedmass = -1;
  int currentdmmass  = -1;
  int npoints = 0;
  vector<pair<int,int> > goodMassPoint;
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

  /// main loop
  int expcounter       = 0;
  int exp_up_counter   = 0;
  int exp_down_counter = 0;
  int obscounter       = 0;
  double minmass_exp   = 100000;
  double minmass_obs   = 100000;

  for (int i = 0; i < tree->GetEntries(); i++){
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
    if(not isGoodMassPoint){ // printout bad limits
      cout<<"Bad limit value: medmass "<<medmass<<" dmmass "<<dmmass<<endl;
      continue;
    }

    // skip some bad points for a smoother contour
    if(skipPoints){
      if(not isDMF){
	if(medmass == 1800 and (dmmass == 1 or dmmass == 5)) continue;
	if(medmass == 1125 and dmmass >= 600) continue;
	if(medmass == 800  and dmmass >= 400) continue;
	if(medmass == 525  and dmmass == 275) continue;
      }
      else{
	if(medmass == 2500 and dmmass == 325) continue;
	if(medmass == 10 and dmmass == 100) continue;
      }
    }

    // filter out some bad mass points
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
      grobs->SetPoint(obscounter, double(medmass), double(dmmass), limit);
      grobu->SetPoint(obscounter, double(medmass), double(dmmass), limit*0.8);
      grobd->SetPoint(obscounter, double(medmass), double(dmmass), limit*1.2);
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

  // make granularity with linear interpolation on 2D
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
  
  // fix min and max values as well as empty points
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
  hexp2->SetContour(2);
  hexp2->SetContourLevel(1,1);
  
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
          if(hexp -> GetXaxis()->GetBinCenter(i) < 50 and hexp -> GetYaxis()->GetBinCenter(j) > 250) hexp->SetBinContent(i,j,minZ*0.01);
          if(hexp_up -> GetXaxis()->GetBinCenter(i) < 50 and hexp_up -> GetYaxis()->GetBinCenter(j) > 250) hexp_up->SetBinContent(i,j,minZ*0.01);
          if(hexp_down -> GetXaxis()->GetBinCenter(i) < 50 and hexp_down -> GetYaxis()->GetBinCenter(j) > 250) hexp_down->SetBinContent(i,j,minZ*0.01);
          if(hobs -> GetXaxis()->GetBinCenter(i) < 50 and hobs -> GetYaxis()->GetBinCenter(j) > 250) hobs->SetBinContent(i,j,minZ*0.01);
          if(hobu -> GetXaxis()->GetBinCenter(i) < 50 and hobu -> GetYaxis()->GetBinCenter(j) > 250) hobu->SetBinContent(i,j,minZ*0.01);
          if(hobd -> GetXaxis()->GetBinCenter(i) < 50 and hobd -> GetYaxis()->GetBinCenter(j) > 250) hobd->SetBinContent(i,j,minZ*0.01);

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
    TFile* ichepobs = TFile::Open("externalFiles/monojet_AV_MM_ICHEP2016_obs.root","READ");
    graph_obs_ichep = (TGraph*) ichepobs->Get("monojet_obs");
    graph_obs_ichep->SetLineWidth(3);
    graph_obs_ichep->SetLineStyle(1);
    graph_obs_ichep->SetLineColor(kViolet);
    graph_obs_ichep->Draw("Lsame");
    TFile* ichepexp = TFile::Open("externalFiles/monojet_AV_MM_ICHEP2016_exp.root","READ");
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
    leg = new TLegend(0.165,0.49,0.50,0.77);
  else
    leg = new TLegend(0.165,0.47,0.50,0.78);
 
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->AddEntry(contour_exp,"Median expected 95% CL","L");
  leg->AddEntry(contour_exp_up,"#pm 1 #sigma_{experiment}","L");
  leg->AddEntry(contour_obs,"Observed 95% CL","L");
  leg->AddEntry(contour_obs_up,"Observed #pm theory unc.","L");
  if(addICHEPContours){
    leg->AddEntry(graph_obs_ichep,"EXO-16-037 observed","L");
    leg->AddEntry(graph_exp_ichep,"EXO-16-037 expected","L");
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
    tex->DrawLatex(0.175,0.80,"#bf{Axial med, Dirac DM, g_{q} = 1, g_{DM} = 1}");
  }
  else
    tex->DrawLatex(0.175,0.80,"#bf{Axial med, Dirac DM, g_{q} = 0.25, g_{DM} = 1}");
  
  
  TLatex *   tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetLineWidth(2);
  tex2->SetTextSize(0.042);
  tex2->SetTextAngle(90);
  tex2->DrawLatex(0.975,0.55,"Observed #sigma_{95% CL}/#sigma_{th}");  

  canvas->RedrawAxis("sameaxis");
  canvas->SaveAs((outputDIR+"/scan_axial_g"+string(coupling)+"_"+string(energy)+"TeV_v2.pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_axial_g"+string(coupling)+"_"+string(energy)+"TeV_v2.png").c_str());

  if(saveOutputFile){
    TFile* outputFile = new TFile((outputDIR+"/fullLikelihood_scan_axial.root").c_str(),"RECREATE");
    outputFile->cd();
    //hexp->Write("scan_expected");
    //hobs->Write("scan_observed");
    //grexp->Write("graph_expected");
    //grexp_up->Write("graph_expected_p1s");
    //grexp_down->Write("graph_expected_m1s");
    //grobs->Write("graph_observed");
    contour_exp->Write("contour_expected");
    contour_obs->Write("contour_observed");
    outputFile->Write();

  }  
}

TGraph*relic_g1(){

    double *x = new double[1000];
    double *y = new double[1000];

    x[0]=335.032; y[0]=15.975 ;
    x[1]=338.56; y[1]=17.634;
    x[2]=339.179; y[2]=17.925;
    x[3]=342.816; y[3]=19.6352;
    x[4]=343.326; y[4]=19.875;
    x[5]=347.072; y[5]=21.6364;
    x[6]=347.525; y[6]=21.825;
    x[7]=351.328; y[7]=22.6862;
    x[8]=355.584; y[8]=22.9093;
    x[9]=359.84; y[9]=23.1325;
    x[10]=364.096; y[10]=23.3556;
    x[11]=365.872; y[11]=23.775;
    x[12]=368.352; y[12]=24.3608;
    x[13]=371.896; y[13]=25.725;
    x[14]=372.608; y[14]=25.9992;
    x[15]=376.864; y[15]=27.6375;
    x[16]=376.961; y[16]=27.675;
    x[17]=381.12; y[17]=29.2759;
    x[18]=382.027; y[18]=29.625;
    x[19]=385.376; y[19]=30.8199;
    x[20]=389.632; y[20]=31.2927;
    x[21]=393.888; y[21]=31.4943;
    x[22]=395.591; y[22]=31.575;
    x[23]=398.144; y[23]=31.6959;
    x[24]=402.4; y[24]=31.8975;
    x[25]=406.656; y[25]=32.0991;
    x[26]=410.912; y[26]=32.3006;
    x[27]=415.168; y[27]=32.5022;
    x[28]=419.424; y[28]=32.7038;
    x[29]=423.68; y[29]=33.3723;
    x[30]=424.149; y[30]=33.525;
    x[31]=427.936; y[31]=34.7569;
    x[32]=430.167; y[32]=35.475;
    x[33]=432.192; y[33]=35.8886;
    x[34]=436.448; y[34]=36.0717;
    x[35]=440.704; y[35]=36.2548;
    x[36]=444.96; y[36]=36.438;
    x[37]=449.216; y[37]=36.6211;
    x[38]=453.472; y[38]=36.8043;
    x[39]=457.728; y[39]=36.9874;
    x[40]=461.984; y[40]=37.1705;
    x[41]=466.24; y[41]=37.3537;
    x[42]=467.898; y[42]=37.425;
    x[43]=470.496; y[43]=37.5368;
    x[44]=474.752; y[44]=37.7199;
    x[45]=479.008; y[45]=38.1434;
    x[46]=483.264; y[46]=38.755;
    x[47]=487.52; y[47]=38.9465;
    x[48]=491.776; y[48]=39.1142;
    x[49]=496.032; y[49]=39.2818;
    x[50]=498.397; y[50]=39.375;
    x[51]=500.288; y[51]=39.5655;
    x[52]=504.544; y[52]=39.9942;
    x[53]=508.8; y[53]=40.423;
    x[54]=513.056; y[54]=40.8517;
    x[55]=517.312; y[55]=41.2804;
    x[56]=517.755; y[56]=41.325;
    x[57]=521.568; y[57]=42.7704;
    x[58]=522.899; y[58]=43.275;
    x[59]=525.824; y[59]=44.0666;
    x[60]=530.08; y[60]=44.589;
    x[61]=534.336; y[61]=45.1148;
    x[62]=534.665; y[62]=45.225;
    x[63]=538.592; y[63]=46.5422;
    x[64]=540.479; y[64]=47.175;
    x[65]=542.848; y[65]=47.9697;
    x[66]=546.293; y[66]=49.125;
    x[67]=547.104; y[67]=49.3971;
    x[68]=551.36; y[68]=50.8245;
    x[69]=552.107; y[69]=51.075;
    x[70]=555.616; y[70]=52.252;
    x[71]=557.921; y[71]=53.025;
    x[72]=559.872; y[72]=53.513;
    x[73]=564.128; y[73]=54.027;
    x[74]=568.384; y[74]=54.4441;
    x[75]=572.64; y[75]=54.8613;
    x[76]=573.8; y[76]=54.975;
    x[77]=576.896; y[77]=55.2785;
    x[78]=581.152; y[78]=55.6957;
    x[79]=585.408; y[79]=56.1129;
    x[80]=589.664; y[80]=56.5301;
    x[81]=591.05; y[81]=56.925;
    x[82]=593.92; y[82]=57.7425;
    x[83]=597.684; y[83]=58.875;
    x[84]=598.176; y[84]=59.0231;
    x[85]=602.432; y[85]=59.996;
    x[86]=606.688; y[86]=60.4537;
    x[87]=610.747; y[87]=60.825;
    x[88]=610.944; y[88]=60.843;
    x[89]=615.2; y[89]=61.2323;
    x[90]=619.456; y[90]=61.6216;
    x[91]=623.712; y[91]=62.0109;
    x[92]=627.968; y[92]=62.4002;
    x[93]=632.066; y[93]=62.775;
    x[94]=632.224; y[94]=62.7895;
    x[95]=636.48; y[95]=63.1787;
    x[96]=640.736; y[96]=63.568;
    x[97]=644.992; y[97]=63.9573;
    x[98]=649.248; y[98]=64.5576;
    x[99]=651; y[99]=64.725;
    x[100]=653.504; y[100]=65.3346;
    x[101]=657.76; y[101]=66.3709;
    x[102]=659.009; y[102]=66.675;
    x[103]=662.016; y[103]=67.7248;
    x[104]=664.595; y[104]=68.625;
    x[105]=666.272; y[105]=69.2104;
    x[106]=670.181; y[106]=70.575;
    x[107]=670.528; y[107]=70.6961;
    x[108]=674.784; y[108]=72.1817;
    x[109]=675.767; y[109]=72.525;
    x[110]=679.04; y[110]=73.6674;
    x[111]=682.196; y[111]=74.475;
    x[112]=683.296; y[112]=74.683;
    x[113]=687.552; y[113]=75.286;
    x[114]=691.808; y[114]=75.889;
    x[115]=695.591; y[115]=76.425;
    x[116]=696.064; y[116]=76.492;
    x[117]=700.32; y[117]=77.095;
    x[118]=704.576; y[118]=78.1577;
    x[119]=705.255; y[119]=78.375;
    x[120]=708.832; y[120]=79.5202;
    x[121]=711.346; y[121]=80.325;
    x[122]=713.088; y[122]=80.8827;
    x[123]=717.344; y[123]=82.2452;
    x[124]=717.467; y[124]=82.275;
    x[125]=721.6; y[125]=83.0346;
    x[126]=725.856; y[126]=83.6039;
    x[127]=730.112; y[127]=84.1733;
    x[128]=730.498; y[128]=84.225;
    x[129]=734.368; y[129]=84.7427;
    x[130]=738.624; y[130]=85.3121;
    x[131]=742.88; y[131]=85.8815;
    x[132]=745.074; y[132]=86.175;
    x[133]=747.136; y[133]=86.4508;
    x[134]=751.392; y[134]=87.0202;
    x[135]=755.648; y[135]=87.5896;
    x[136]=758.245; y[136]=88.125;
    x[137]=759.904; y[137]=88.4338;
    x[138]=764.16; y[138]=89.1877;
    x[139]=768.416; y[139]=89.7401;
    x[140]=770.996; y[140]=90.075;
    x[141]=772.672; y[141]=90.6762;
    x[142]=776.432; y[142]=92.025;
    x[143]=776.928; y[143]=92.2028;
    x[144]=781.184; y[144]=93.7295;
    x[145]=781.869; y[145]=93.975;
    x[146]=785.44; y[146]=95.2561;
    x[147]=787.305; y[147]=95.925;
    x[148]=789.696; y[148]=96.7827;
    x[149]=792.741; y[149]=97.875;
    x[150]=793.952; y[150]=98.2584;
    x[151]=798.208; y[151]=99.1993;
    x[152]=801.79; y[152]=99.825;
    x[153]=802.464; y[153]=99.9428;
    x[154]=806.72; y[154]=100.686;
    x[155]=810.976; y[155]=101.43;
    x[156]=812.539; y[156]=101.775;
    x[157]=815.232; y[157]=102.37;
    x[158]=819.292; y[158]=103.725;
    x[159]=819.488; y[159]=103.791;
    x[160]=823.744; y[160]=105.211;
    x[161]=825.132; y[161]=105.675;
    x[162]=828; y[162]=106.632;
    x[163]=831.493; y[163]=107.625;
    x[164]=832.256; y[164]=107.795;
    x[165]=836.512; y[165]=108.504;
    x[166]=840.768; y[166]=109.213;
    x[167]=842.94; y[167]=109.575;
    x[168]=845.024; y[168]=109.922;
    x[169]=849.28; y[169]=110.631;
    x[170]=853.536; y[170]=111.34;
    x[171]=854.646; y[171]=111.525;
    x[172]=857.792; y[172]=112.049;
    x[173]=862.048; y[173]=112.758;
    x[174]=866.304; y[174]=113.467;
    x[175]=866.348; y[175]=113.475;
    x[176]=870.56; y[176]=114.313;
    x[177]=874.816; y[177]=115.252;
    x[178]=875.771; y[178]=115.425;
    x[179]=879.072; y[179]=116.631;
    x[180]=881.107; y[180]=117.375;
    x[181]=883.328; y[181]=118.187;
    x[182]=886.443; y[182]=119.325;
    x[183]=887.584; y[183]=119.742;
    x[184]=891.778; y[184]=121.275;
    x[185]=891.84; y[185]=121.298;
    x[186]=896.096; y[186]=122.853;
    x[187]=897.114; y[187]=123.225;
    x[188]=900.352; y[188]=124.408;
    x[189]=902.689; y[189]=125.175;
    x[190]=904.608; y[190]=125.676;
    x[191]=908.864; y[191]=126.529;
    x[192]=911.839; y[192]=127.125;
    x[193]=913.12; y[193]=127.382;
    x[194]=917.376; y[194]=128.234;
    x[195]=921.572; y[195]=129.075;
    x[196]=921.632; y[196]=129.087;
    x[197]=925.888; y[197]=129.94;
    x[198]=929.078; y[198]=131.025;
    x[199]=930.144; y[199]=131.388;
    x[200]=934.4; y[200]=132.852;
    x[201]=934.758; y[201]=132.975;
    x[202]=938.656; y[202]=134.228;
    x[203]=941.792; y[203]=134.925;
    x[204]=942.912; y[204]=135.14;
    x[205]=947.168; y[205]=135.958;
    x[206]=951.424; y[206]=136.775;
    x[207]=951.944; y[207]=136.875;
    x[208]=955.68; y[208]=137.593;
    x[209]=959.936; y[209]=138.41;
    x[210]=962.096; y[210]=138.825;
    x[211]=964.192; y[211]=139.32;
    x[212]=968.448; y[212]=140.324;
    x[213]=970.36; y[213]=140.775;
    x[214]=972.704; y[214]=141.621;
    x[215]=976.015; y[215]=142.725;
    x[216]=976.96; y[216]=142.978;
    x[217]=981.216; y[217]=143.832;
    x[218]=983.901; y[218]=144.675;
    x[219]=985.472; y[219]=145.168;
    x[220]=989.728; y[220]=146.613;
    x[221]=989.763; y[221]=146.625;
    x[222]=993.984; y[222]=148.058;
    x[223]=995.508; y[223]=148.575;
    x[224]=998.24; y[224]=149.503;
    x[225]=1001.25; y[225]=150.525;
    x[226]=1002.5; y[226]=150.947;
    x[227]=1006.75; y[227]=152.392;
    x[228]=1007.01; y[228]=152.475;
    x[229]=1011.01; y[229]=153.501;
    x[230]=1015.26; y[230]=154.317;
    x[231]=1015.83; y[231]=154.425;
    x[232]=1019.52; y[232]=155.133;
    x[233]=1023.78; y[233]=155.95;
    x[234]=1025.99; y[234]=156.375;
    x[235]=1028.03; y[235]=156.766;
    x[236]=1032.29; y[236]=157.582;
    x[237]=1036.16; y[237]=158.325;
    x[238]=1036.54; y[238]=158.398;
    x[239]=1040.8; y[239]=159.574;
    x[240]=1042.99; y[240]=160.275;
    x[241]=1045.06; y[241]=160.911;
    x[242]=1049.31; y[242]=161.88;
    x[243]=1051.19; y[243]=162.225;
    x[244]=1053.57; y[244]=162.66;
    x[245]=1057.82; y[245]=163.44;
    x[246]=1061.83; y[246]=164.175;
    x[247]=1062.08; y[247]=164.217;
    x[248]=1066.34; y[248]=164.944;
    x[249]=1070.59; y[249]=165.672;
    x[250]=1073.25; y[250]=166.125;
    x[251]=1074.85; y[251]=166.386;
    x[252]=1079.1; y[252]=167.079;
    x[253]=1083.36; y[253]=167.771;
    x[254]=1085.38; y[254]=168.075;
    x[255]=1087.62; y[255]=168.346;
    x[256]=1091.87; y[256]=168.543;
    x[257]=1096.13; y[257]=168.964;
    x[258]=1100.38; y[258]=169.611;
    x[259]=1103.11; y[259]=170.025;
    x[260]=1104.64; y[260]=170.258;
    x[261]=1108.9; y[261]=170.905;
    x[262]=1113.15; y[262]=171.552;
    x[263]=1115.93; y[263]=171.975;
    x[264]=1117.41; y[264]=172.199;
    x[265]=1121.66; y[265]=172.846;
    x[266]=1125.92; y[266]=173.493;
    x[267]=1128.83; y[267]=173.925;
    x[268]=1130.18; y[268]=174.088;
    x[269]=1134.43; y[269]=174.274;
    x[270]=1138.69; y[270]=174.46;
    x[271]=1142.94; y[271]=174.646;
    x[272]=1147.2; y[272]=174.832;
    x[273]=1151.46; y[273]=175.162;
    x[274]=1155.71; y[274]=175.768;
    x[275]=1156.46; y[275]=175.875;
    x[276]=1159.97; y[276]=176.374;
    x[277]=1164.22; y[277]=176.981;
    x[278]=1168.48; y[278]=177.587;
    x[279]=1170.15; y[279]=177.825;
    x[280]=1172.74; y[280]=178.194;
    x[281]=1176.99; y[281]=178.642;
    x[282]=1181.25; y[282]=178.831;
    x[283]=1185.5; y[283]=179.007;
    x[284]=1189.76; y[284]=179.184;
    x[285]=1194.02; y[285]=179.361;
    x[286]=1198.27; y[286]=179.537;
    x[287]=1202.53; y[287]=179.714;
    x[288]=1203.57; y[288]=179.775;
    x[289]=1206.78; y[289]=179.963;
    x[290]=1211.04; y[290]=180.534;
    x[291]=1215.3; y[291]=181.106;
    x[292]=1219.55; y[292]=181.677;
    x[293]=1219.91; y[293]=181.725;
    x[294]=1223.81; y[294]=182.175;
    x[295]=1228.06; y[295]=182.384;
    x[296]=1232.32; y[296]=182.553;
    x[297]=1236.58; y[297]=182.721;
    x[298]=1240.83; y[298]=182.89;
    x[299]=1245.09; y[299]=183.058;
    x[300]=1249.34; y[300]=183.226;
    x[301]=1253.6; y[301]=183.395;
    x[302]=1257.86; y[302]=183.563;
    x[303]=1260.52; y[303]=183.675;
    x[304]=1262.11; y[304]=183.741;
    x[305]=1266.37; y[305]=184.281;
    x[306]=1270.62; y[306]=184.821;
    x[307]=1274.88; y[307]=185.168;
    x[308]=1279.14; y[308]=185.35;
    x[309]=1283.39; y[309]=185.511;
    x[310]=1286.4; y[310]=185.625;
    x[311]=1287.65; y[311]=185.672;
    x[312]=1291.9; y[312]=185.834;
    x[313]=1296.16; y[313]=185.995;
    x[314]=1300.42; y[314]=186.156;
    x[315]=1304.67; y[315]=186.317;
    x[316]=1308.93; y[316]=186.478;
    x[317]=1313.18; y[317]=186.64;
    x[318]=1317.44; y[318]=186.801;
    x[319]=1321.7; y[319]=187.267;
    x[320]=1325.84; y[320]=187.575;
    x[321]=1325.95; y[321]=187.582;
    x[322]=1330.21; y[322]=187.737;
    x[323]=1334.46; y[323]=187.891;
    x[324]=1338.72; y[324]=188.046;
    x[325]=1342.98; y[325]=188.2;
    x[326]=1347.23; y[326]=188.355;
    x[327]=1351.49; y[327]=188.509;
    x[328]=1355.74; y[328]=188.664;
    x[329]=1360; y[329]=188.818;
    x[330]=1364.26; y[330]=188.973;
    x[331]=1368.51; y[331]=189.127;
    x[332]=1372.77; y[332]=189.282;
    x[333]=1376.97; y[333]=189.525;
    x[334]=1377.02; y[334]=189.534;
    x[335]=1381.28; y[335]=189.973;
    x[336]=1385.54; y[336]=190.413;
    x[337]=1389.79; y[337]=190.853;
    x[338]=1394.05; y[338]=191.292;
    x[339]=1395.82; y[339]=191.475;
    x[340]=1398.3; y[340]=192.168;
    x[341]=1402.56; y[341]=193.355;
    x[342]=1402.81; y[342]=193.425;
    x[343]=1406.82; y[343]=194.542;
    x[344]=1409.8; y[344]=195.375;
    x[345]=1411.07; y[345]=195.729;
    x[346]=1415.33; y[346]=196.819;
    x[347]=1417.99; y[347]=197.325;
    x[348]=1419.58; y[348]=197.602;
    x[349]=1423.84; y[349]=198.341;
    x[350]=1428.1; y[350]=199.08;
    x[351]=1428.96; y[351]=199.275;
    x[352]=1432.35; y[352]=200.035;
    x[353]=1436.61; y[353]=201.166;
    x[354]=1436.83; y[354]=201.225;
    x[355]=1440.86; y[355]=202.297;
    x[356]=1444.17; y[356]=203.175;
    x[357]=1445.12; y[357]=203.428;
    x[358]=1449.38; y[358]=204.559;
    x[359]=1451.74; y[359]=205.125;
    x[360]=1453.63; y[360]=205.516;
    x[361]=1457.89; y[361]=206.223;
    x[362]=1462.14; y[362]=206.929;
    x[363]=1463.02; y[363]=207.075;
    x[364]=1466.4; y[364]=207.636;
    x[365]=1470.66; y[365]=208.342;
    x[366]=1474.77; y[366]=209.025;
    x[367]=1474.91; y[367]=209.048;
    x[368]=1479.17; y[368]=209.755;
    x[369]=1483.42; y[369]=210.461;
    x[370]=1485.98; y[370]=210.975;
    x[371]=1487.68; y[371]=211.315;
    x[372]=1491.94; y[372]=212.333;
    x[373]=1495.31; y[373]=212.925;
    x[374]=1496.19; y[374]=213.065;
    x[375]=1500.45; y[375]=213.741;
    x[376]=1504.7; y[376]=214.417;
    x[377]=1507.59; y[377]=214.875;
    x[378]=1508.96; y[378]=215.446;
    x[379]=1512.27; y[379]=216.825;
    x[380]=1513.22; y[380]=217.266;
    x[381]=1516.44; y[381]=218.775;
    x[382]=1517.47; y[382]=219.259;
    x[383]=1520.6; y[383]=220.725;
    x[384]=1521.73; y[384]=221.252;
    x[385]=1524.87; y[385]=222.675;
    x[386]=1525.98; y[386]=223.124;
    x[387]=1530.07; y[387]=224.625;
    x[388]=1530.24; y[388]=224.688;
    x[389]=1534.5; y[389]=226.253;
    x[390]=1535.37; y[390]=226.575;
    x[391]=1538.75; y[391]=227.817;
    x[392]=1540.57; y[392]=228.525;
    x[393]=1543.01; y[393]=229.475;
    x[394]=1545.18; y[394]=230.475;
    x[395]=1547.26; y[395]=231.433;
    x[396]=1549.42; y[396]=232.425;
    x[397]=1551.52; y[397]=233.392;
    x[398]=1553.73; y[398]=234.375;
    x[399]=1555.78; y[399]=235.191;
    x[400]=1558.89; y[400]=236.325;
    x[401]=1560.03; y[401]=236.739;
    x[402]=1564.26; y[402]=238.275;
    x[403]=1564.29; y[403]=238.287;
    x[404]=1568.54; y[404]=239.88;
    x[405]=1569.47; y[405]=240.225;
    x[406]=1572.8; y[406]=241.835;
    x[407]=1573.5; y[407]=242.175;
    x[408]=1577.06; y[408]=243.891;
    x[409]=1577.54; y[409]=244.125;
    x[410]=1581.31; y[410]=245.946;
    x[411]=1581.59; y[411]=246.075;
    x[412]=1585.57; y[412]=247.709;
    x[413]=1586.42; y[413]=248.025;
    x[414]=1589.82; y[414]=249.293;
    x[415]=1591.66; y[415]=249.975;
    x[416]=1594.08; y[416]=250.876;
    x[417]=1596.84; y[417]=251.925;
    x[418]=1598.34; y[418]=252.492;
    x[419]=1601.28; y[419]=253.875;
    x[420]=1602.59; y[420]=254.489;
    x[421]=1605.44; y[421]=255.825;
    x[422]=1606.85; y[422]=256.485;
    x[423]=1609.6; y[423]=257.775;
    x[424]=1611.1; y[424]=258.456;
    x[425]=1614.31; y[425]=259.725;
    x[426]=1615.36; y[426]=260.105;
    x[427]=1619.62; y[427]=261.651;
    x[428]=1619.68; y[428]=261.675;
    x[429]=1623.87; y[429]=263.197;
    x[430]=1625.05; y[430]=263.625;
    x[431]=1628.13; y[431]=264.85;
    x[432]=1629.95; y[432]=265.575;
    x[433]=1632.38; y[433]=266.837;
    x[434]=1633.71; y[434]=267.525;
    x[435]=1636.64; y[435]=269.041;
    x[436]=1637.48; y[436]=269.475;
    x[437]=1640.9; y[437]=271.181;
    x[438]=1641.45; y[438]=271.425;
    x[439]=1645.15; y[439]=272.937;
    x[440]=1646.23; y[440]=273.375;
    x[441]=1649.41; y[441]=274.673;
    x[442]=1651.01; y[442]=275.325;
    x[443]=1653.66; y[443]=276.409;
    x[444]=1655.41; y[444]=277.275;
    x[445]=1657.92; y[445]=278.525;
    x[446]=1659.31; y[446]=279.225;
    x[447]=1662.18; y[447]=280.673;
    x[448]=1663.17; y[448]=281.175;
    x[449]=1666.43; y[449]=282.821;
    x[450]=1667.05; y[450]=283.125;
    x[451]=1670.69; y[451]=284.701;
    x[452]=1671.62; y[452]=285.075;
    x[453]=1674.94; y[453]=286.4;
    x[454]=1676.51; y[454]=287.025;
    x[455]=1679.2; y[455]=288.1;
    x[456]=1681.39; y[456]=288.975;
    x[457]=1683.46; y[457]=289.911;
    x[458]=1685.69; y[458]=290.925;
    x[459]=1687.71; y[459]=292.004;
    x[460]=1689.34; y[460]=292.875;
    x[461]=1691.97; y[461]=294.277;
    x[462]=1692.99; y[462]=294.825;
    x[463]=1696.22; y[463]=296.485;
    x[464]=1696.85; y[464]=296.775;
    x[465]=1700.48; y[465]=298.324;
    x[466]=1701.42; y[466]=298.725;
    x[467]=1704.74; y[467]=300.142;
    x[468]=1705.98; y[468]=300.675;
    x[469]=1708.99; y[469]=301.959;
    x[470]=1710.32; y[470]=302.625;
    x[471]=1713.25; y[471]=304.085;
    x[472]=1714.19; y[472]=304.575;
    x[473]=1717.5; y[473]=306.307;
    x[474]=1717.92; y[474]=306.525;
    x[475]=1721.66; y[475]=308.475;
    x[476]=1721.76; y[476]=308.527;
    x[477]=1725.97; y[477]=310.425;
    x[478]=1726.02; y[478]=310.443;
    x[479]=1730.27; y[479]=312.227;
    x[480]=1730.62; y[480]=312.375;
    x[481]=1734.53; y[481]=314.011;
    x[482]=1735.28; y[482]=314.325;
    x[483]=1738.78; y[483]=316.046;
    x[484]=1739.25; y[484]=316.275;
    x[485]=1742.87; y[485]=318.225;
    x[486]=1743.04; y[486]=318.317;
    x[487]=1746.49; y[487]=320.175;
    x[488]=1747.3; y[488]=320.611;
    x[489]=1750.2; y[489]=322.125;
    x[490]=1751.55; y[490]=322.753;
    x[491]=1754.59; y[491]=324.075;
    x[492]=1755.81; y[492]=324.607;
    x[493]=1759.06; y[493]=326.025;
    x[494]=1760.06; y[494]=326.461;
    x[495]=1763.54; y[495]=327.975;
    x[496]=1764.32; y[496]=328.315;
    x[497]=1767.6; y[497]=329.925;
    x[498]=1768.58; y[498]=330.407;
    x[499]=1771.36; y[499]=331.875;
    x[500]=1772.83; y[500]=332.655;
    x[501]=1775.05; y[501]=333.825;
    x[502]=1777.09; y[502]=334.859;
    x[503]=1779.09; y[503]=335.775;
    x[504]=1781.34; y[504]=336.742;
    x[505]=1783.64; y[505]=337.725;
    x[506]=1785.6; y[506]=338.566;
    x[507]=1788.19; y[507]=339.675;
    x[508]=1789.86; y[508]=340.537;
    x[509]=1791.96; y[509]=341.625;
    x[510]=1794.11; y[510]=342.782;
    x[511]=1795.59; y[511]=343.575;
    x[512]=1798.37; y[512]=345.073;
    x[513]=1799.21; y[513]=345.525;
    x[514]=1802.62; y[514]=347.363;
    x[515]=1802.84; y[515]=347.475;
    x[516]=1806.88; y[516]=349.361;
    x[517]=1807.03; y[517]=349.425;
    x[518]=1811.14; y[518]=351.228;
    x[519]=1811.47; y[519]=351.375;
    x[520]=1815.39; y[520]=353.095;
    x[521]=1815.92; y[521]=353.325;
    x[522]=1819.65; y[522]=354.962;
    x[523]=1820.3; y[523]=355.275;
    x[524]=1823.9; y[524]=356.999;
    x[525]=1824.33; y[525]=357.225;
    x[526]=1828.02; y[526]=359.175;
    x[527]=1828.16; y[527]=359.247;
    x[528]=1831.86; y[528]=361.125;
    x[529]=1832.42; y[529]=361.382;
    x[530]=1836.34; y[530]=363.075;
    x[531]=1836.67; y[531]=363.221;
    x[532]=1840.84; y[532]=365.025;
    x[533]=1840.93; y[533]=365.073;
    x[534]=1844.49; y[534]=366.975;
    x[535]=1845.18; y[535]=367.345;
    x[536]=1848.14; y[536]=368.925;
    x[537]=1849.44; y[537]=369.617;
    x[538]=1851.8; y[538]=370.875;
    x[539]=1853.7; y[539]=371.889;
    x[540]=1855.45; y[540]=372.825;
    x[541]=1857.95; y[541]=374.103;
    x[542]=1859.4; y[542]=374.775;
    x[543]=1862.21; y[543]=376.009;
    x[544]=1863.84; y[544]=376.725;
    x[545]=1866.46; y[545]=377.877;
    x[546]=1868.28; y[546]=378.675;
    x[547]=1870.72; y[547]=379.745;
    x[548]=1872.72; y[548]=380.625;
    x[549]=1874.98; y[549]=381.614;
    x[550]=1877.05; y[550]=382.575;
    x[551]=1879.23; y[551]=383.589;
    x[552]=1881.01; y[552]=384.525;
    x[553]=1883.49; y[553]=385.824;
    x[554]=1884.78; y[554]=386.475;
    x[555]=1887.74; y[555]=387.83;
    x[556]=1889.12; y[556]=388.425;
    x[557]=1892; y[557]=389.716;
    x[558]=1893.47; y[558]=390.375;
    x[559]=1896.26; y[559]=391.85;
    x[560]=1897.15; y[560]=392.325;
    x[561]=1900.51; y[561]=394.103;
    x[562]=1900.84; y[562]=394.275;
    x[563]=1904.52; y[563]=396.225;
    x[564]=1904.77; y[564]=396.355;
    x[565]=1908.21; y[565]=398.175;
    x[566]=1909.02; y[566]=398.608;
    x[567]=1912.01; y[567]=400.125;
    x[568]=1913.28; y[568]=400.709;
    x[569]=1916.4; y[569]=402.075;
    x[570]=1917.54; y[570]=402.571;
    x[571]=1920.86; y[571]=404.025;
    x[572]=1921.79; y[572]=404.433;
    x[573]=1925.32; y[573]=405.975;
    x[574]=1926.05; y[574]=406.295;
    x[575]=1929.77; y[575]=407.925;
    x[576]=1930.3; y[576]=408.157;
    x[577]=1934.13; y[577]=409.875;
    x[578]=1934.56; y[578]=410.066;
    x[579]=1937.93; y[579]=411.825;
    x[580]=1938.82; y[580]=412.265;
    x[581]=1942.14; y[581]=413.775;
    x[582]=1943.07; y[582]=414.205;
    x[583]=1946.36; y[583]=415.725;
    x[584]=1947.33; y[584]=416.234;
    x[585]=1950.07; y[585]=417.675;
    x[586]=1951.58; y[586]=418.469;
    x[587]=1953.78; y[587]=419.625;
    x[588]=1955.84; y[588]=420.705;
    x[589]=1957.5; y[589]=421.575;
    x[590]=1960.1; y[590]=422.94;
    x[591]=1961.21; y[591]=423.525;
    x[592]=1964.35; y[592]=425.176;
    x[593]=1964.95; y[593]=425.475;
    x[594]=1968.61; y[594]=427.144;
    x[595]=1969.26; y[595]=427.425;
    x[596]=1972.86; y[596]=428.993;
    x[597]=1973.74; y[597]=429.375;
    x[598]=1977.12; y[598]=430.842;
    x[599]=1978.23; y[599]=431.325;
    x[600]=1981.38; y[600]=432.691;
    x[601]=1982.72; y[601]=433.275;
    x[602]=1985.63; y[602]=434.54;
    x[603]=1987.21; y[603]=435.225;
    x[604]=1989.89; y[604]=436.389;
    x[605]=1991.42; y[605]=437.175;
    x[606]=1994.14; y[606]=438.513;
    x[607]=1995.5; y[607]=439.125;
    x[608]=1998.4; y[608]=440.507;
    x[609]=1999.59; y[609]=441.075;
    x[610]=2002.66; y[610]=442.674;
    x[611]=2003.33; y[611]=443.025;
    x[612]=2006.91; y[612]=444.895;
    x[613]=2007.07; y[613]=444.975;
    x[614]=2010.8; y[614]=446.925;
    x[615]=2011.17; y[615]=447.116;
    x[616]=2014.54; y[616]=448.875;
    x[617]=2015.42; y[617]=449.337;
    x[618]=2018.28; y[618]=450.825;
    x[619]=2019.68; y[619]=451.521;
    x[620]=2022.45; y[620]=452.775;
    x[621]=2023.94; y[621]=453.418;
    x[622]=2026.96; y[622]=454.725;
    x[623]=2028.19; y[623]=455.259;
    x[624]=2031.46; y[624]=456.675;
    x[625]=2032.45; y[625]=457.101;
    x[626]=2035.97; y[626]=458.625;
    x[627]=2036.7; y[627]=458.943;
    x[628]=2040.48; y[628]=460.575;
    x[629]=2040.96; y[629]=460.785;
    x[630]=2044.98; y[630]=462.525;
    x[631]=2045.22; y[631]=462.626;
    x[632]=2049.11; y[632]=464.475;
    x[633]=2049.47; y[633]=464.656;
    x[634]=2053.08; y[634]=466.425;
    x[635]=2053.73; y[635]=466.762;
    x[636]=2056.84; y[636]=468.375;
    x[637]=2057.98; y[637]=468.966;
    x[638]=2060.61; y[638]=470.325;
    x[639]=2062.24; y[639]=471.171;
    x[640]=2064.37; y[640]=472.275;
    x[641]=2066.5; y[641]=473.376;
    x[642]=2068.14; y[642]=474.225;
    x[643]=2070.75; y[643]=475.581;
    x[644]=2071.9; y[644]=476.175;
    x[645]=2075.01; y[645]=477.702;
    x[646]=2075.95; y[646]=478.125;
    x[647]=2079.26; y[647]=479.556;
    x[648]=2080.47; y[648]=480.075;
    x[649]=2083.52; y[649]=481.392;
    x[650]=2084.99; y[650]=482.025;
    x[651]=2087.78; y[651]=483.229;
    x[652]=2089.51; y[652]=483.975;
    x[653]=2092.03; y[653]=485.065;
    x[654]=2094.02; y[654]=485.925;
    x[655]=2096.29; y[655]=486.902;
    x[656]=2098.54; y[656]=487.875;
    x[657]=2100.54; y[657]=488.738;
    x[658]=2102.91; y[658]=489.825;
    x[659]=2104.8; y[659]=490.761;
    x[660]=2106.81; y[660]=491.775;
    x[661]=2109.06; y[661]=492.929;
    x[662]=2110.61; y[662]=493.725;
    x[663]=2113.31; y[663]=495.116;
    x[664]=2114.4; y[664]=495.675;
    x[665]=2117.57; y[665]=497.302;
    x[666]=2118.2; y[666]=497.625;
    x[667]=2121.82; y[667]=499.489;
    x[668]=2121.99; y[668]=499.575;
    x[669]=2125.79; y[669]=501.525;
    x[670]=2126.08; y[670]=501.676;
    x[671]=2129.76; y[671]=503.475;
    x[672]=2130.34; y[672]=503.732;
    x[673]=2134.26; y[673]=505.425;
    x[674]=2134.59; y[674]=505.566;
    x[675]=2138.79; y[675]=507.375;
    x[676]=2138.85; y[676]=507.4;
    x[677]=2143.1; y[677]=509.234;
    x[678]=2143.32; y[678]=509.325;
    x[679]=2147.36; y[679]=511.068;
    x[680]=2147.84; y[680]=511.275;
    x[681]=2151.62; y[681]=512.902;
    x[682]=2152.37; y[682]=513.225;
    x[683]=2155.87; y[683]=514.759;
    
    TGraph *lrelic = new TGraph(683,x,y);
    lrelic->SetLineColor(kBlue+2);    
    lrelic->SetLineWidth(-802);
    lrelic->SetFillStyle(3005);
    lrelic->SetFillColor(kBlue+2);
    return lrelic;

}


