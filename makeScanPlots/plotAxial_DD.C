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

/// formulas for DD limit

double axialF(double mMED,double mDM){  
  if (mMED <= 0) 
    return 10;
  
  mMED/=1000;
  double mR = (0.939*mDM)/(0.939+mDM);
  double c = 2.4E-42;//4.6E-41;
  return c*mR*mR/(mMED*mMED*mMED*mMED);    
}

TGraph * makeOBA(TGraph *Graph1){

  TGraph *gr = new TGraph();
  double X;
  double Y;
  int pp=0;
  Graph1->GetPoint(1,X,Y);
  for (double MDM=1;MDM<=Y;MDM+=0.1){
    
    gr->SetPoint(pp,MDM,axialF(X,MDM));
    pp++;
  }
  for (int p =1;p<Graph1->GetN();p++){
    Graph1->GetPoint(p,X,Y);
    if (!(X >1)) continue;
    if ((X <100)) continue;
    gr->SetPoint(pp,Y,axialF(X,Y));
    pp++;
  }
  gr->GetXaxis()->SetTitle("m_{DM}");
  gr->GetYaxis()->SetTitle("#sigma^{SD}_{DM-proton}");
  gr->SetName(Form("%s_DD",Graph1->GetName()));
  gr->SetLineStyle(Graph1->GetLineStyle());
  gr->SetLineColor(Graph1->GetLineColor());
  gr->SetLineWidth(Graph1->GetLineWidth());
  
  return gr;
}


/////
static float nbinsX = 1500;
static float nbinsY = 1000;
static float minX = 0;
static float minY = 1.;
static float maxX = 2500;
static float maxY = 1000;
static float minZ = 0.01;
static float maxZ = 20;

static float minX_dd = 1;
static float maxX_dd = 1000;
static double minY_dd = 5e-46;
static double maxY_dd = 5e-34;

static bool saveOutputFile = false;

TGraph* Pico2L();
TGraph* Pico60();
TGraph* SuperKtt();
TGraph* IceCubett();

////////

void plotAxial_DD(string inputFileName, string outputDirectory, string coupling = "025", string energy = "13") {

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDirectory).c_str());
  setTDRStyle();

  TFile *file = TFile::Open(inputFileName.c_str(),"READ");
  TTree *tree = (TTree*)file->Get("limit");

  TGraph2D* grexp = new TGraph2D();
  TGraph2D* grobs = new TGraph2D();
  
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
  
  int expcounter = 0;
  int obscounter = 0;

  for (int i = 0; i < tree->GetEntries(); i++){
    
    tree->GetEntry(i);
    
    if (quantile != 0.5 && quantile != -1) continue;

    int c = code(mh);
    int medmass = mmed(mh, c);
    int dmmass = mdm(mh, c);

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

    // filter out some bad mass points                                                                                                                                                           
    if(medmass == 925  and dmmass >= 600) continue;
    if(medmass == 1000 and dmmass >= 600) continue;
    if(medmass == 1125 and dmmass >= 600) continue;
    if(medmass == 1200 and dmmass >= 600) continue;
    if(medmass == 1325 and dmmass >= 600) continue;
    if(medmass == 325) continue;

    if (quantile == 0.5) {
      expcounter++;
      grexp->SetPoint(expcounter, double(medmass), double(dmmass), limit);
    }
    if (quantile == -1) {
      obscounter++;
      grobs->SetPoint(obscounter, double(medmass), double(dmmass), limit);
    }
  }
  tree->ResetBranchAddresses();

  ///
  TH2D* hexp = new TH2D("hexp", "",      nbinsX, minX, maxX, nbinsY, minY, maxY);
  TH2D* hobs = new TH2D("hobs", "", nbinsX, minX, maxX, nbinsY, minY, maxY);
 
  
  // make granularity                                                                                                                                                                                 
  for (int i   = 1; i <= nbinsX; i++) {
    for (int j = 1; j <= nbinsY; j++) {
      hexp->SetBinContent(i,j,grexp->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetYaxis()->GetBinCenter(j)));
      hobs->SetBinContent(i,j,grobs->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetYaxis()->GetBinCenter(j)));
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

  TH2* hexp2 = (TH2*)hexp->Clone("hexp2");
  TH2* hobs2 = (TH2*)hobs->Clone("hobs2");

  double contours[1]; contours[0]=1;
  hexp2->SetContour(1,contours);
  hobs2->SetContour(1,contours);

  hexp2->Draw("contz list");
  gPad->Update();

  /// import the expected contour line
  TObjArray *lContoursE = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  std::vector<double> lXE;
  std::vector<double> lYE;
  int lTotalContsE = lContoursE->GetSize();
  for(int i0 = 0; i0 < lTotalContsE; i0++){ // there can be more than one contour if the are not connected
    TList * pContLevel = (TList*)lContoursE->At(i0);
    TGraph *pCurv = (TGraph*)pContLevel->First();
    for(int i1 = 0; i1 < pContLevel->GetSize(); i1++){
      for(int i2  = 0; i2 < pCurv->GetN(); i2++) {
	lXE.push_back(pCurv->GetX()[i2]); 
	lYE.push_back(pCurv->GetY()[i2]);
      }
      pCurv->SetLineColor(kBlack);                                                                                                            
      pCurv = (TGraph*)pContLevel->After(pCurv);                                                                                        
    }
  }
  if(lXE.size() == 0) {
    lXE.push_back(0); 
    lYE.push_back(0); 
  }

  TGraph *lTotalE = new TGraph(lXE.size(),&lXE[0],&lYE[0]);
  lTotalE->SetLineColor(kRed);
  lTotalE->SetLineWidth(3);

  // observed one
  hobs2->Draw("contz list");
  gPad->Update();

  TObjArray *lContours = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  std::vector<double> lX;
  std::vector<double> lY;
  int lTotalConts = lContours->GetSize();
  for(int i0 = 0; i0 < lTotalConts; i0++){
    TList * pContLevel = (TList*)lContours->At(i0);
    TGraph *pCurv = (TGraph*)pContLevel->First();
    for(int i1 = 0; i1 < pContLevel->GetSize(); i1++){
      for(int i2  = 0; i2 < pCurv->GetN(); i2++) {
	lX.push_back(pCurv->GetX()[i2]);
	lY.push_back(pCurv->GetY()[i2]);
      }
      pCurv->SetLineColor(kBlack);
      pCurv = (TGraph*)pContLevel->After(pCurv);
    }
  }
  if(lX.size() == 0) {
    lX.push_back(0); 
    lY.push_back(0); 
  }

  TGraph *lTotal = new TGraph(lX.size(),&lX[0],&lY[0]);
  lTotal->SetLineColor(kBlack);
  lTotal->SetLineWidth(3);

  // make the DD limits
  TGraph  *lM0 = Pico2L();
  TGraph  *lM1 = Pico60(); //lM1->SetLineStyle(kDashed);
  TGraph  *lM2 = IceCubett();
  TGraph  *lM3 = SuperKtt();

  lM0->SetLineColor(kBlue);
  lM1->SetLineColor(kBlue+2);
  lM2->SetLineColor(kAzure+1);
  lM2->SetLineStyle(7);
  lM3->SetLineColor(kAzure+8);
  lM3->SetLineStyle(7);

  TGraph *DDE_graph = makeOBA(lTotalE);
  TGraph *DD_graph  = makeOBA(lTotal);
  
  TCanvas* canvas = new TCanvas("canvas","canvas",750,600);
  canvas->SetLogx();
  canvas->SetLogy();

  TH1* frame = canvas->DrawFrame(minX_dd,minY_dd,maxX_dd,maxY_dd,"");
  frame->GetYaxis()->SetTitle("#sigma^{SD}_{DM-proton} [cm^{2}]");    
  frame->GetXaxis()->SetTitle("m_{DM} [GeV]");
  frame->GetXaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.25);
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->CenterTitle();
  frame->Draw();

  lM0->Draw("L SAME");
  lM1->Draw("L SAME");  
  lM2->Draw("L SAME");
  lM3->Draw("L SAME");

  DDE_graph->SetLineColor(kBlack);
  DD_graph->SetLineColor(kRed);
  DDE_graph->Draw("L SAME");
  DD_graph->Draw("L SAME");


  canvas->SetLogx();
  canvas->SetLogy();

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
  leg->AddEntry(lM1 ,"PICO-60","L");
  leg->AddEntry(lM0 ,"PICO-2L","L");
  leg->AddEntry(lM2 ,"IceCube #tau^{+}#tau^{-}","L");
  leg->AddEntry(lM3 ,"Super-K #tau^{+}#tau^{-}","L");
    
  leg->Draw("SAME");
  CMS_lumi(canvas,"35.9",false,true,false,0,-0.22);

  canvas->RedrawAxis("samesaxis");

  TLatex * tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.030);
  tex->Draw();
  if (coupling == "1"){
    tex->DrawLatex(0.75,0.82,"#bf{Axial med, Dirac DM,}");
    tex->DrawLatex(0.75,0.78,"#bf{g_{q} = 1, g_{DM} = 1}");
  }
  else{
    tex->DrawLatex(0.75,0.82,"#bf{Axial med, Dirac DM,}");
    tex->DrawLatex(0.75,0.78,"#bf{g_{q} = 0.25, g_{DM} = 1}");
  }
  ///////
  canvas->SaveAs((outputDirectory+"/scanDD_axial_g"+coupling+"_"+energy+"TeV_v1.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/scanDD_axial_g"+coupling+"_"+energy+"TeV_v1.png").c_str(),"png");

  if(saveOutputFile){

    TFile*outfile = new TFile(("axial_g"+coupling+"_DD.root").c_str(),"RECREATE");
    DDE_graph->SetName("expected");
    DD_graph->SetName("observed");
    DDE_graph->Write();
    DD_graph->Write();
    outfile->Write();
    outfile->Close();

  }  
}

TGraph *IceCubett() {
    int i0 = -1;
    double *lX = new double[1000];
    double *lY = new double[1000];
    i0++; lX[i0] = 9.785451025760974;  lY[i0] = 3.57519300867316e-36;
    i0++; lX[i0] = 11.894692473017603; lY[i0] = 7.915478203925451e-37;
    i0++; lX[i0] = 14.936688262634869; lY[i0] = 1.0726470749278954e-37;
    i0++; lX[i0] = 21.83179687403854;  lY[i0] = 6.7206482798073e-39;
    i0++; lX[i0] = 30.8884359647748;   lY[i0] = 1.2932280829152238e-39;
    i0++; lX[i0] = 58.56806967738426;  lY[i0] = 3.925623361155188e-40;
    i0++; lX[i0] = 105.19086393080569; lY[i0] = 1.6337965307944188e-40;
    i0++; lX[i0] = 175.11691609819135; lY[i0] = 5.909811163389418e-41;
    i0++; lX[i0] = 264.4185952824892;  lY[i0] = 3.2560446770685855e-41;
    i0++; lX[i0] = 552.7683532427855;  lY[i0] = 2.73242483102675e-41;
    i0++; lX[i0] = 1206.792640639329;  lY[i0] = 5.706176195201639e-41;
    i0++; lX[i0] = 1802.5485470871574; lY[i0] = 1.1109263788798718e-40;
    i0++; lX[i0] = 2904.74705624372;   lY[i0] = 2.488508258147588e-40;
    i0++; lX[i0] = 4680.903310144558;  lY[i0] = 5.97929177133309e-40;
    i0++; lX[i0] = 6767.920662784252;  lY[i0] = 1.2932280829152238e-39;
    i0++; lX[i0] = 9892.14386559403;   lY[i0] = 2.5177652110127317e-39;
    TGraph *lLimit = new TGraph(i0,lX,lY);
    lLimit->SetLineWidth(3.);
    return lLimit;
}

TGraph *IceCubeWW() {
    int i0 = -1;
    double *lX = new double[1000];
    double *lY = new double[1000];
    i0++; lX[i0] = 20.159400454447134; lY[i0] = 1.2805229753360497e-38;
    i0++; lX[i0] = 35.17779848372143; lY[i0] = 1.218486987237799e-39;
    i0++; lX[i0] = 49.43446640476189; lY[i0] = 2.630985438722169e-40;
    i0++; lX[i0] = 100.68970642844887; lY[i0] = 2.623840449807147e-40;
    i0++; lX[i0] = 246.9092233665341; lY[i0] = 1.2817852845594771e-40;
    i0++; lX[i0] = 502.9126239917861; lY[i0] = 1.4587259517827274e-40;
    i0++; lX[i0] = 982.9633229801225; lY[i0] = 4.183856199578464e-40;
    i0++; lX[i0] = 2993.087429584735; lY[i0] = 4.855589407808413e-39;
    i0++; lX[i0] = 4960.46498138491;  lY[i0] = 1.5084208438229934e-38;
    TGraph *lLimit = new TGraph(i0,lX,lY);
    lLimit->SetLineWidth(3.);
    return lLimit;
}
TGraph *SuperKtt() {
    int i0 = -1;
    double *lX = new double[1000];
    double *lY = new double[1000];
    i0++; lX[i0] = 3.997708395644498; lY[i0] =  2.2090912214681107e-40;
    i0++; lX[i0] = 5.981830941814304; lY[i0] =  1.6808125720951077e-40;
    i0++; lX[i0] = 9.912370243964038; lY[i0] =  1.3142995585251433e-40;
    i0++; lX[i0] = 20.139313905605064; lY[i0] =  1.3507152344698235e-40;
    i0++; lX[i0] = 50.194353560796664; lY[i0] =  1.2108380316531138e-40;
    i0++; lX[i0] = 100.26456553695655; lY[i0] =  1.2108380316531138e-40;
    i0++; lX[i0] = 200.24745234223136; lY[i0] =  1.3881398900233809e-40;
    TGraph *lLimit = new TGraph(i0,lX,lY);
    lLimit->SetLineWidth(3.);
    return lLimit;
}



TGraph*Pico2L() { 
  int i0 = -1;
  double *lX = new double[1000];
  double *lY = new double[1000];

  //i0++;  lX[i0] = 3.3113e+00;   lY[i0] =   8.5638e-33;
  i0++;  lX[i0] = 3.6308e+00;   lY[i0] =   5.4104e-36;
  i0++;  lX[i0] = 3.9811e+00;   lY[i0] =   5.6411e-37;
  i0++;  lX[i0] = 4.3652e+00;   lY[i0] =   1.3643e-37;
  i0++;  lX[i0] = 4.7863e+00;   lY[i0] =   4.8940e-38;
  i0++;  lX[i0] = 5.2481e+00;   lY[i0] =   2.2120e-38;
  i0++;  lX[i0] = 5.7544e+00;   lY[i0] =   1.1705e-38;
  i0++;  lX[i0] = 6.3096e+00;   lY[i0] =   6.9885e-39;
  i0++;  lX[i0] = 6.9183e+00;   lY[i0] =   4.5488e-39;
  i0++;  lX[i0] = 7.5858e+00;   lY[i0] =   3.1742e-39;
  i0++;  lX[i0] = 8.3176e+00;   lY[i0] =   2.3452e-39;
  i0++;  lX[i0] = 9.1201e+00;   lY[i0] =   1.8170e-39;
  i0++;  lX[i0] = 1.0000e+01;   lY[i0] =   1.4649e-39;
  i0++;  lX[i0] = 1.0965e+01;   lY[i0] =   1.2214e-39;
  i0++;  lX[i0] = 1.2023e+01;   lY[i0] =   1.0481e-39;
  i0++;  lX[i0] = 1.3183e+01;   lY[i0] =   9.2175e-40;
  i0++;  lX[i0] = 1.4454e+01;   lY[i0] =   8.2814e-40;
  i0++;  lX[i0] = 1.5849e+01;   lY[i0] =   7.5832e-40;
  i0++;  lX[i0] = 1.7378e+01;   lY[i0] =   7.0569e-40;
  i0++;  lX[i0] = 1.9055e+01;   lY[i0] =   6.6618e-40;
  i0++;  lX[i0] = 2.0893e+01;   lY[i0] =   6.3693e-40;
  i0++;  lX[i0] = 2.2909e+01;   lY[i0] =   6.1596e-40;
  i0++;  lX[i0] = 2.5119e+01;   lY[i0] =   6.0183e-40;
  i0++;  lX[i0] = 3.1623e+01;   lY[i0] =   5.9045e-40;
  i0++;  lX[i0] = 3.9811e+01;   lY[i0] =   6.0684e-40;
  i0++;  lX[i0] = 5.0119e+01;   lY[i0] =   6.4756e-40;
  i0++;  lX[i0] = 6.3096e+01;   lY[i0] =   7.1255e-40;
  i0++;  lX[i0] = 7.9433e+01;   lY[i0] =   8.0411e-40;
  i0++;  lX[i0] = 1.0000e+02;   lY[i0] =   9.2643e-40;
  i0++;  lX[i0] = 1.2589e+02;   lY[i0] =   1.0857e-39;
  i0++;  lX[i0] = 1.5849e+02;   lY[i0] =   1.2900e-39;
  i0++;  lX[i0] = 1.9953e+02;   lY[i0] =   1.5504e-39;
  i0++;  lX[i0] = 2.5119e+02;   lY[i0] =   1.8804e-39;
  i0++;  lX[i0] = 3.1623e+02;   lY[i0] =   2.2976e-39;
  i0++;  lX[i0] = 1.0000e+03;   lY[i0] =   6.7076e-39;
  i0++;  lX[i0] = 3.1623e+03;   lY[i0] =   2.0675e-38;
  i0++;  lX[i0] = 1.0000e+04;   lY[i0] =   6.4850e-38;
    
  TGraph *lLimit = new TGraph(i0,lX,lY);
  lLimit->SetLineWidth(3.);
  return lLimit;
}

TGraph*Pico60() {
  int i0 = -1;
  double *lX = new double[1000];
  double *lY = new double[1000];
  i0++; lX[i0] = 10.867638347999604; lY[i0] = 1.009070636600476*1e-37;
  i0++; lX[i0] = 11.411801474233114; lY[i0] = 6.196648512966848*1e-38;
  i0++; lX[i0] = 13.034473073844428; lY[i0] = 2.7993604566432167*1e-38;
  i0++; lX[i0] = 14.719319205005647; lY[i0] = 1.0367791970603659*1e-38;
  i0++; lX[i0] = 18.28293451100823;  lY[i0] = 3.7711220494242113*1e-39;
  i0++; lX[i0] = 26.51606929269927;  lY[i0] = 1.2308243469132812*1e-39;
  i0++; lX[i0] = 38.859017640523;    lY[i0] = 7.15981763641212*1e-40;
  i0++; lX[i0] = 62.582565374717156; lY[i0] = 5.764800915198607*1e-40;
  i0++; lX[i0] = 110.788609270474;   lY[i0] = 5.869856050894904*1e-40;
  i0++; lX[i0] = 212.95457453334728; lY[i0] = 8.73326162382851*1e-40;
  i0++; lX[i0] = 394.88771974007585; lY[i0] = 1.4221361511653464*1e-39;
  i0++; lX[i0] = 749.8176636967278;  lY[i0] = 2.400999982340401*1e-39;
  i0++; lX[i0] = 1457.7138852764463; lY[i0] = 4.4366873309786245*1e-39;
  i0++; lX[i0] = 2701.945262438203;  lY[i0] = 8.499859846090154*1e-39;
  i0++; lX[i0] = 5707.251556973959;  lY[i0] = 1.7822982971317245*1e-38;
  i0++; lX[i0] = 9508.725458075965;  lY[i0] = 2.8503747409110666*1e-38;
  TGraph *lLimit = new TGraph(i0,lX,lY);
  lLimit->SetLineWidth(3.);
  return lLimit;
}
