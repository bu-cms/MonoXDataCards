#include "../CMS_lumi.h"

// to get mass point info
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

// constatins
double const mtop = 173.34;
double const mb   = 4.18;

// Width
double WidthDM(double mMED, double mF){
    double fac = (1- 4*mF*mF/(mMED*mMED));
    double vev = 246.;
    return TMath::Sqrt(fac)*mMED/(8*TMath::Pi());
}
double Width(double mMED, double mF){
    double fac = (1- 4*mF*mF/(mMED*mMED));
    double vev = 246.;
    return TMath::Sqrt(fac)*mF*mF*mMED/(8*TMath::Pi()*vev*vev);
}

bool STOP = false;

// conversion
static float gDM = 1;
static float gq  = 1;

double vecF(double mMED,double mDM){

  if(mMED<=0) return 10;

  if(2*mDM > 0.99*mMED) {
    mMED = 2*mDM;
    STOP=true;
  }
  
  double iGamma = 1;
  if (mMED < 2.*mtop) iGamma = 3*Width(mMED,mb) + WidthDM(mMED,mDM) ;
  else iGamma = 3*Width(mMED,mb) +  3*Width(mMED,mtop) + WidthDM(mMED,mDM) ;

  double  lVal    = 3*mb*mb/(2*TMath::Pi()*246*246);
  double  lDenom  = (mMED*mMED-4.*mDM*mDM)*(mMED*mMED-4.*mDM*mDM);
  lDenom += mMED*mMED*iGamma*iGamma;
  lVal /= lDenom;
  lVal *= mDM *mDM*gDM*gDM*gq*gq;
  lVal *= TMath::Sqrt(1.-(mb*mb)/(mDM*mDM));
  double c = 1.167E-17;

  return c*lVal;
}

void majoranatodirac(TGraph* gr){
  double X;
  double Y;
  double fac = 2;
  for (int p=0;p<gr->GetN();p++){
    gr->GetPoint(p,X,Y);
    gr->SetPoint(p,X,fac*Y);
  }
}

// conversion macro
TGraph * makeOBV(TGraph *Graph1){

  STOP=false;
  TGraph *gr = new TGraph();
  double X;
  double Y;
  int pp = 0;
  Graph1->GetPoint(0,X,Y);
  for (double MDM=mb;MDM<=Y;MDM+=0.1){
    gr->SetPoint(pp,MDM,vecF(X,MDM));
    pp++;
  }

  for (int p =0;p<Graph1->GetN();p++){    
    Graph1->GetPoint(p,X,Y);
    if (!(X>0)) continue; // skip negative mMED 
    if (STOP) continue;   // stop if off-shell      
    if (X < 2*mtop and p%10 !=0) continue;
    gr->SetPoint(pp,Y,vecF(X,Y));    
    pp++;
  }
  gr->GetXaxis()->SetTitle("m_{DM}");
  gr->GetYaxis()->SetTitle("< #sigma v > (cm^{3}/s)");
  gr->SetName(Form("%s_DD",Graph1->GetName()));
  gr->SetLineStyle(Graph1->GetLineStyle());
  gr->SetLineColor(Graph1->GetLineColor());
  gr->SetLineWidth(Graph1->GetLineWidth());

  return gr;
  
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
static float nbinsX = 400;
static float nbinsY = 250;
static float minX = 0;
static float minY = 5; // cannot go below 5 since mDM > mb for the formula used to translate into DM-nucleon xsec
static float maxX = 600;
static float maxY = 300;
static float minZ = 0.1;
static float maxZ = 10;

static float minX_dd  = 5;
static float maxX_dd  = 400;
static double minY_dd = 1e-31;
static double maxY_dd = 1e-22;
static int  reductionForContour = 5;

void plotPseudoScalar_DD(string inputFileName, string outputDirectory, string coupling = "1", string energy = "13") {


  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDirectory).c_str());
  setTDRStyle();
  
  TFile *file = TFile::Open(inputFileName.c_str(),"READ");
  TTree *tree = (TTree*)file->Get("limit");
  
  TFile* file2 = new TFile("FermiLat.root");
  TGraph* grdd = (TGraph*)file2->Get("FermiLAT8Y");   
  
  TGraph2D* grexp = new TGraph2D();
  TGraph2D* grobs = new TGraph2D();
  
  double mh;
  double limit;
  float quantile;
  double minObs = 0;
  double minmass = 100000;

  tree->SetBranchAddress("mh",&mh);
  tree->SetBranchAddress("limit",&limit);
  tree->SetBranchAddress("quantileExpected",&quantile);

  // identify possible bad limits
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

  ////////////
  int expcounter = 0;
  int obscounter = 0;
  
  for (int i = 0; i < tree->GetEntries(); i++){
    
    tree->GetEntry(i);
    
    if (quantile != 0.5 && quantile != -1) continue;
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


    if (quantile == 0.5) {
      expcounter++;
      grexp->SetPoint(expcounter, double(medmass), double(dmmass), limit);
    }
    if (quantile == -1) {
      obscounter++;
      grobs->SetPoint(obscounter, double(medmass), double(dmmass), limit);

      if(medmass <= minmass and dmmass < medmass/2){
        minObs = limit;
        minmass = medmass;
      }     
    }
  }
  
  tree->ResetBranchAddresses();
  
  TH2D* hexp = new TH2D("hexp", "", nbinsX, minX, maxX, nbinsY, minY, maxY);
  TH2D* hobs = new TH2D("hobs", "", nbinsX, minX, maxX, nbinsY, minY, maxY);
  
  // make granularity                                                                                                                                                                           
  for (int i = 1; i   <= nbinsX; i++) {
    for (int j = 1; j <= nbinsY; j++) {
      hexp->SetBinContent(i,j,grexp->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetYaxis()->GetBinCenter(j)));
      hobs->SetBinContent(i,j,grobs->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetYaxis()->GetBinCenter(j)));
      
      if(hexp->GetXaxis()->GetBinCenter(i) <= 30 and hexp->GetYaxis()->GetBinCenter(j) < hexp->GetXaxis()->GetBinCenter(i)/2){
	hexp->SetBinContent(i,j,minObs);
	hobs->SetBinContent(i,j,minObs);
      }
    }
  }
  
  for(int i = 0; i < nbinsX; i++){
    for(int j = 0; j < nbinsY; j++){

      if(hexp -> GetBinContent(i,j) <= 0) hexp->SetBinContent(i,j,maxZ);
      if(hobs -> GetBinContent(i,j) <= 0) hobs->SetBinContent(i,j,maxZ);
      
      if(hexp -> GetBinContent(i,j) > maxZ) hexp->SetBinContent(i,j,maxZ);
      if(hobs -> GetBinContent(i,j) > maxZ) hobs->SetBinContent(i,j,maxZ);
      
      if(hexp -> GetBinContent(i,j) < minZ) hexp->SetBinContent(i,j,minZ);
      if(hobs -> GetBinContent(i,j) < minZ) hobs->SetBinContent(i,j,minZ);
    }
  }
  
  hexp->Smooth();
  hobs->Smooth();

  ////////////////                                                                                                                                                                              
  TH2* hexp2 = (TH2*)hexp->Clone("hexp2");
  TH2* hobs2 = (TH2*)hobs->Clone("hobs2");

  double contours[1]; contours[0]=1;
  hexp2->SetContour(1,contours);
  hobs2->SetContour(1,contours);
  
  // Save expected conoturs as TGraph
  hexp2->Draw("contz list");
  gPad->Update();
  
  TGraph* lTotalE = produceContour(reductionForContour);
  lTotalE->SetLineColor(kBlack);
  lTotalE->SetLineWidth(3);
  lTotalE->SetLineStyle(7);

  /// Save observed contours as TGraph
  hobs2->Draw("contz list");
  gPad->Update();

  TGraph* lTotal = produceContour(reductionForContour);
  lTotal->SetLineColor(kRed);
  lTotal->SetLineWidth(3);
  
  TGraph *DDE_graph = makeOBV(lTotalE);
  TGraph *DD_graph  = makeOBV(lTotal);
  grdd->SetLineColor(kAzure+1);
  grdd->SetLineWidth(3);
  
  
  TCanvas* canvas = new TCanvas("canvas","canvas",750,600);
  canvas->SetLogx();
  canvas->SetLogy();
  
  TH1* frame = canvas->DrawFrame(minX_dd,minY_dd,maxX_dd,maxY_dd,"");
  frame->GetYaxis()->SetTitle("<#sigma v > (cm^{3}/s)");
  frame->GetXaxis()->SetTitle("m_{DM} [GeV]");
  frame->GetXaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.25);
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->CenterTitle();
  frame->Draw();
  
  frame->Draw();
  grdd->Draw("L SAME");
  DDE_graph->Draw("L SAME");
  DD_graph->Draw("L SAME");
  
  gPad->SetRightMargin(0.28);
  gPad->RedrawAxis();
  gPad->Modified();
  gPad->Update();

  TLegend *leg = new TLegend(0.75,0.45,0.97,0.72,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(DDE_graph,"CMS exp. 90% CL","L");
  leg->AddEntry(DD_graph ,"CMS obs. 90% CL","L");
  leg->AddEntry(grdd ,"FermiLAT","L");
  
  leg->Draw("SAME");
  CMS_lumi(canvas,"35.9",false,true,false,0,-0.22);
    
  TLatex * tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.030);
  tex->Draw();
  if (coupling == "1"){
    tex->DrawLatex(0.75,0.82,"#bf{Pseudoscalar med}");
    tex->DrawLatex(0.75,0.79,"#bf{Dirac DM}");
    tex->DrawLatex(0.75,0.75,"#bf{g_{q} = 1, g_{DM} = 1}");
  }
  else{
    tex->DrawLatex(0.75,0.82,"#bf{Pseudoscalar med}");
    tex->DrawLatex(0.75,0.79,"#bf{Dirac DM}");
    tex->DrawLatex(0.75,0.75,"#bf{g_{q} = 0.25, g_{DM} = 1}");
  }
  
  ///////                                                                                                                                                                                          
  canvas->SaveAs((outputDirectory+"/scanDD_pseudo_g"+coupling+"_"+energy+"TeV_v1.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/scanDD_pseudo_g"+coupling+"_"+energy+"TeV_v1.png").c_str(),"png");
  
  if(saveOutputFile){   
    TFile*outfile = new TFile((outputDirectory+"/pseudo_g"+coupling+"_DD.root").c_str(),"RECREATE");
    hobs2->Write("contour_obs");
    hexp2->Write("contour_exp");
    lTotalE->Write("contour_exp_graph");
    lTotal->Write("contour_obs_graph");
    DDE_graph->SetName("expected");
    DD_graph->SetName("observed");
    DDE_graph->Write();
    DD_graph->Write();
    outfile->Write();
    outfile->Close();    
  }


}

