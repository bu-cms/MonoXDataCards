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
      if(i1 > 0) continue;
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
static bool saveOutputFile = true;
static bool addRelicDensity = true;
static float nbinsX = 600;
static float nbinsY = 400;
//static float nbinsX = 80;
//static float nbinsY = 120;
static float minX = 0;
static float minY = 1;
static float maxX = 600;
static float maxY = 400;
static float minZ = 0.1;
static float maxZ = 10;
static int   reductionForContour = 20;
static bool  whiteOut = true;
static bool  addPreliminary = false;
static bool  skipPoints = false;

void plotPseudoScalar(string inputFileName, string outputDIR, string coupling = "1", string energy = "13") {
  
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
    TFile* inputFile1 =  TFile::Open("externalFiles/relic_pseudo_gq1.root","READ");
    TList* list = (TList*) inputFile1->Get("mytlist");
    TIter next(list);
    TObject* object = 0;
    while ((object = next())){
      relicDensity.push_back(dynamic_cast<TGraph*>(object));
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

  int currentmedmass = -1;
  int currentdmmass  = -1;
  int npoints = 0;

  vector<pair<int,int> > goodMassPoint;

  for(int i = 0; i < tree->GetEntries(); i++){
    tree->GetEntry(i);

    int c       = code(mh);
    int medmass = mmed(mh, c);
    int dmmass  = mdm(mh, c);

    if(dmmass > medmass/2) continue;

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

  
  int expcounter       = 0;
  int exp_up_counter   = 0;
  int exp_down_counter = 0;
  int obscounter       = 0;
  int minmass_obs      = 10000;
  int minmass_exp      = 10000;
  double minObs;

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
    if(not isGoodMassPoint){
      cout<<"Bad limit value: medmass "<<medmass<<" dmmass "<<dmmass<<endl;
      continue;
    }


    if (quantile == 0.5) { // expected limit
      grexp->SetPoint(expcounter, double(medmass), double(dmmass), limit);
      expcounter++;
      if(medmass <= minmass_exp and dmmass < medmass/2){// not consider off-shell points
	minmass_exp = medmass;
      }
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
      if(medmass <= minmass_obs and dmmass < medmass/2){
	minObs = limit;
	minmass_obs = medmass;
      }
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

  //////// fix mass points below min med mass generated                                                                                                                            
  for (int i = 1; i   <= nbinsX; i++) {
    for (int j = 1; j <= nbinsY; j++) {
      if(hexp->GetXaxis()->GetBinCenter(i) < minmass_exp and hexp->GetYaxis()->GetBinCenter(j) < hexp->GetXaxis()->GetBinCenter(i)/2){
        hexp_up->SetBinContent(i,j,hexp_up->GetBinContent(hexp_up->FindBin(minmass_exp,hexp_up->GetYaxis()->GetBinCenter(j))));
        hexp_down->SetBinContent(i,j,hexp_down->GetBinContent(hexp_down->FindBin(minmass_exp,hexp_down->GetYaxis()->GetBinCenter(j))));
        hexp->SetBinContent(i,j,hexp->GetBinContent(hexp->FindBin(minmass_exp,hexp->GetYaxis()->GetBinCenter(j))));
      }
      if(hexp->GetXaxis()->GetBinCenter(i) < minmass_obs and hobs->GetYaxis()->GetBinCenter(j) < hobs->GetXaxis()->GetBinCenter(i)/2){
        hobs->SetBinContent(i,j,hobs->GetBinContent(hobs->FindBin(minmass_obs,hobs->GetYaxis()->GetBinCenter(j))));
        hobd->SetBinContent(i,j,hobd->GetBinContent(hobd->FindBin(minmass_obs,hobd->GetYaxis()->GetBinCenter(j))));
        hobu->SetBinContent(i,j,hobu->GetBinContent(hobu->FindBin(minmass_obs,hobu->GetYaxis()->GetBinCenter(j))));
      }
    }
  }
  
  ////////////////  
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

  ////////////////
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
  canvas->SetLeftMargin(0.135);
  canvas->SetLogz();
  
  TH1* frame = canvas->DrawFrame(minX,minY,maxX,maxY, "");
  frame->GetYaxis()->SetTitle("m_{DM} [GeV]");
  frame->GetXaxis()->SetTitle("m_{med} [GeV]");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.30);
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
  TGraph* contour_obs_up = produceContour(1);

  hobd2->GetZaxis()->SetLabelSize(0);
  hobd2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_obs_dw = produceContour(1);

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
  }

  frame->Draw();
  hobs->Draw("COLZ same");

  if(addRelicDensity){
    if(addRelicDensity){
      int iline = 0;
      for(auto graph : relicDensity){
        graph->SetLineColor(kBlue+2);
	if(iline == 0)	  
	  graph->SetLineWidth(802);
	else
	  graph->SetLineWidth(502);
        graph->SetFillStyle(3005);
        graph->SetFillColor(kBlue+2);
        graph->Draw("SAME");
	iline++;
      }
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


  TLegend *leg = new TLegend(0.165,0.48,0.53,0.77);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->AddEntry(contour_exp,"Median expected 95% CL","L");
  leg->AddEntry(contour_exp_up,"#pm 1 #sigma_{experiment}","L");
  leg->AddEntry(contour_obs,"Observed 95% CL","L");
  leg->AddEntry(contour_obs_up,"Observed #pm theory unc.","L");
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
    tex->DrawLatex(0.175,0.80,"#bf{Pseudoscalar med, Dirac DM, g_{q} = 1, g_{DM} = 1}");
  }
  else
    tex->DrawLatex(0.175,0.80,"#bf{Pseudoscalar med, Dirac DM, g_{q} = 0.25, g_{DM} = 1}");
  
  
  TLatex *   tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetLineWidth(2);
  tex2->SetTextSize(0.042);
  tex2->SetTextAngle(90);
  tex2->DrawLatex(0.975,0.55,"Observed #sigma_{95% CL}/#sigma_{th}");
  
  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR+"/scan_pseudoscalar_g"+string(coupling)+"_"+string(energy)+"TeV_v2.pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_pseudoscalar_g"+string(coupling)+"_"+string(energy)+"TeV_v2.png").c_str());

  if(saveOutputFile){
    TFile* outputFile = new TFile((outputDIR+"/scan_pseudoscalar_g"+string(coupling)+"_"+string(energy)+"TeV_v2.root").c_str(),"RECREATE");
    outputFile->cd();
    hobs->Write("Observed_limit");
    hexp->Write("Expected_limit");
    contour_obs->Write("Observed_contour");
    contour_exp->Write("Expected_contour");
    outputFile->Write();
  }
}
