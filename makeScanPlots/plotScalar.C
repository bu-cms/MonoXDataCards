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
static bool saveOutputFile = true;
static bool addRelicDensity = false;
static float nbinsX = 1000;
static float nbinsY = 600;
static float minX = 0;
static float minY = 1;
static float maxX = 600;
static float maxY = 300;
static float minZ = 0.5;
static float maxZ = 10;
static int   reductionForContour = 20;

void plotScalar(string inputFileName, string outputDIR, string coupling = "1", string energy = "13") {

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
  TGraph* dd = NULL;

  if(addRelicDensity){
    TFile* file2 = TFile::Open("");
    TGraph* wm   = (TGraph*)file2->Get("wmap_0");
    TGraph* dd   = (TGraph*)file2->Get("DD_mass");
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
  
  /////////////-----
  int expcounter       = 0;
  int exp_up_counter   = 0;
  int exp_down_counter = 0;
  int obscounter       = 0;
  double minmass_exp = 100000;
  double minmass_obs = 100000;
  double min_exp = 100000;
  double min_obs = 100000;

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
      if(medmass <= 100 and limit > 1.0) continue;
      if(medmass >  100 and limit < 1.0) continue;
      grexp->SetPoint(expcounter, double(medmass), double(dmmass), limit);
      expcounter++;
      if(medmass <= minmass_exp){
        minmass_exp = medmass;
	min_exp = limit;
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
      if(medmass <= 150 and limit > 1.0) continue;
      if(medmass >  150 and limit < 1.0) continue;

      grobs->SetPoint(obscounter, double(medmass), double(dmmass), limit);
      grobu->SetPoint(obscounter, double(medmass), double(dmmass), limit*0.8);
      grobd->SetPoint(obscounter, double(medmass), double(dmmass), limit*1.2);
      obscounter++;      
      if(medmass <= minmass_obs){
        minmass_obs = medmass;
	min_obs = limit;
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

  // fix mass points below min med mass generated                                                                                                                                                     
  for (int i = 1; i   <= nbinsX; i++) {
    for (int j = 1; j <= nbinsY; j++) {
      if(hexp->GetXaxis()->GetBinCenter(i) <= max(minmass_exp,30.) and hexp->GetYaxis()->GetBinCenter(j) < hexp->GetXaxis()->GetBinCenter(i)/2){
	hexp_up->SetBinContent(i,j,min_exp);
	hexp_down->SetBinContent(i,j,min_exp);
	hexp->SetBinContent(i,j,min_exp);
      }
      if(hobs->GetXaxis()->GetBinCenter(i) <= max(minmass_obs,30.) and hobs->GetYaxis()->GetBinCenter(j) < hobs->GetXaxis()->GetBinCenter(i)/2){
	hobs->SetBinContent(i,j,min_obs);
	hobd->SetBinContent(i,j,min_obs);
	hobu->SetBinContent(i,j,min_obs);
      }
    }
  }

  //////////
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
  canvas->SetLeftMargin(0.13);
  canvas->SetLogz();
  
  TH1* frame = canvas->DrawFrame(minX,minY,maxX,maxY, "");
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitle("m_{DM} [GeV]");
  frame->GetXaxis()->SetTitle("m_{med} [GeV]");
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
  TGraph* contour_obs = produceContour(1);

  hobu2->GetZaxis()->SetLabelSize(0);
  hobu2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_obs_up = produceContour(1);

  hobd2->GetZaxis()->SetLabelSize(0);
  hobd2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_obs_dw = produceContour(1);

  frame->Draw();
  hobs->Draw("COLZ SAME");

  if(addRelicDensity){
    wm->SetFillStyle(3005);
    wm->Draw("SAME");
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

  CMS_lumi(canvas,"35.9",false,true,false,0,-0.09);

  TLegend *leg = new TLegend(0.175,0.58,0.45,0.78);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->AddEntry(contour_exp,"Median expected 95% CL","L");
  leg->AddEntry(contour_exp_up,"Expected #pm 1 s.d._{experiment}","L");
  leg->AddEntry(contour_obs,"Observed 95% CL","L");
  leg->AddEntry(contour_obs_up,"Observed #pm 1 s.d._{theory}","L");
  if(addRelicDensity)
    leg->AddEntry(wm   ,"Planck+WMAP Relic","F");
  leg->Draw("SAME");

  TLatex * tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.030);
  tex->Draw();
  if (coupling == "1"){
    tex->DrawLatex(0.175,0.80,"#bf{Scalar med, Dirac DM, g_{q} = 1, g_{DM} = 1}");
  }
  else
    tex->DrawLatex(0.175,0.80,"#bf{Scalar med, Dirac DM, g_{q} = 0.25, g_{DM} = 1}");
    
  TLatex *   tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetLineWidth(2);
  tex2->SetTextSize(0.042);
  tex2->SetTextAngle(90);
  tex2->DrawLatex(0.97,0.55,"Observed #sigma_{95% CL}/#sigma_{th}");
  
  gPad->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR+"/scan_scalar_g"+string(coupling)+"_"+string(energy)+"TeV_v2.pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_scalar_g"+string(coupling)+"_"+string(energy)+"TeV_v2.png").c_str());

  if(saveOutputFile){
    TFile* outputFile = new TFile((outputDIR+"/fullLikelihood_scan_scalar.root").c_str(),"RECREATE");
    outputFile->cd();
    hexp->Write("scan_expected");
    hobs->Write("scan_observed");
    hexp->Write("contour_expected");
    hobs->Write("contour_observed");
    grexp->Write("graph_expected");
    grexp_up->Write("graph_expected_p1s");
    grexp_down->Write("graph_expected_m1s");
    grobs->Write("graph_observed");
    contour_exp->Write("contour_exp");
    contour_obs->Write("contour_obs");

    outputFile->Write();

  }

}
