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
  if (mMED < 12) return axialF(12,mDM);
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
  Graph1->GetPoint(0,X,Y);
  for (double MDM=1;MDM<=Y;MDM+=0.1){    
    gr->SetPoint(pp,MDM,axialF(X,MDM));
    pp++;
  }
  for (int p =1;p<Graph1->GetN();p++){
    Graph1->GetPoint(p,X,Y);
    if(X < 100) continue;
    gr->SetPoint(pp,Y,axialF(X,Y));
    pp++;
  }

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


/////
static float nbinsX = 1000;
static float nbinsY = 600;
static float minX = 0;
static float minY = 1;
static float maxX = 2500;
static float maxY = 1200;
static float minZ = 0.01;
static float maxZ = 10;

static float minX_dd = 1;
static float maxX_dd = 1200;
static double minY_dd = 1e-45;
static double maxY_dd = 1e-28;

static bool saveOutputFile      = true;
static int  reductionForContour = 20;
static bool addPreliminary = true;

TGraph* Pico2L();
TGraph* Pico60();
TGraph* SuperKtt();
TGraph* SuperKbb();
TGraph* IceCubett();
TGraph* IceCubebb();
TGraph* PicassoFinalGraph();

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
  double minmass = 100000;

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
    if(medmass == 1925 and (dmmass >=10 or dmmass <= 100)) continue;
    if(medmass == 1125 and dmmass >= 600) continue;
    if(medmass == 925  and dmmass >= 600) continue;
    if(medmass == 800  and dmmass >= 400) continue;
    if(medmass == 525  and dmmass >= 275) continue;

    if (quantile == 0.5) {
      grexp->SetPoint(expcounter, double(medmass), double(dmmass), limit);
      expcounter++;
    }
    if (quantile == -1) {
      grobs->SetPoint(obscounter, double(medmass), double(dmmass), limit);
      obscounter++;
      if(medmass <= minmass) minmass = medmass;
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

  // make granularity with linear interpolation on 2D                                                                                                                                               
  for (int i = 1; i   <= nbinsX; i++) {
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
  TGraph* lTotalE = produceContour(reductionForContour);
  lTotalE->SetLineColor(kBlack);
  lTotalE->SetLineStyle(2);
  lTotalE->SetLineWidth(3);

  // observed one
  hobs2->Draw("contz list");
  gPad->Update();

  TGraph* lTotal = produceContour(reductionForContour);
  lTotal->SetLineColor(kRed);
  lTotal->SetLineWidth(3);

  // make the DD limits
  TGraph  *lM0 = Pico60(); //lM1->SetLineStyle(kDashed);
  TGraph  *lM1 = PicassoFinalGraph(); 
  TGraph  *lM2 = IceCubebb();
  TGraph  *lM3 = SuperKbb();
  TGraph  *lM4 = IceCubett();

  lM0->SetLineColor(kBlue);
  lM1->SetLineColor(kBlue+2);
  lM2->SetLineColor(kAzure+1);
  lM2->SetLineStyle(7);
  lM3->SetLineColor(kAzure+8);
  lM3->SetLineStyle(7);
  lM4->SetLineColor(kAzure-7);
  lM4->SetLineStyle(7);

  TGraph *DDE_graph = makeOBA(lTotalE);
  TGraph *DD_graph  = makeOBA(lTotal);

  TCanvas* canvas = new TCanvas("canvas","canvas",600,625);
  canvas->SetLogx();
  canvas->SetLogy();

  TH1* frame = canvas->DrawFrame(minX_dd,minY_dd,maxX_dd,maxY_dd,"");
  frame->GetYaxis()->SetTitle("#sigma^{SD}_{DM-proton} [cm^{2}]");    
  frame->GetXaxis()->SetTitle("m_{DM} [GeV]");
  frame->GetXaxis()->SetLabelSize(0.032);
  frame->GetYaxis()->SetLabelSize(0.032);
  frame->GetXaxis()->SetTitleSize(0.042);
  frame->GetYaxis()->SetTitleSize(0.042);
  frame->GetYaxis()->SetTitleOffset(1.65);
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->CenterTitle();
  frame->Draw();

  lM0->Draw("L SAME");
  lM1->Draw("L SAME");  
  lM2->Draw("L SAME");
  lM3->Draw("L SAME");
  lM4->Draw("L SAME");

  DDE_graph->Draw("C SAME");
  DD_graph->Draw("C SAME");


  canvas->SetLogx();
  canvas->SetLogy();

  gPad->SetLeftMargin(0.15);
  gPad->RedrawAxis();
  gPad->Modified();
  gPad->Update();

  TLegend *leg = new TLegend(0.25,0.52,0.55,0.78,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->AddEntry(DDE_graph,"CMS exp. 90% CL","L");
  leg->AddEntry(DD_graph ,"CMS obs. 90% CL","L");
  leg->AddEntry(lM0 ,"PICO-60","L");
  leg->AddEntry(lM1 ,"Picasso","L");
  leg->AddEntry(lM2 ,"IceCube b#bar{b}","L");
  leg->AddEntry(lM3 ,"Super-K b#bar{b}","L");
  leg->AddEntry(lM4 ,"Super-K t#bar{t}","L");
    
  leg->Draw("SAME");
  if(addPreliminary)
    CMS_lumi(canvas,"35.9",false,false,false,0.05,0);
  else
    CMS_lumi(canvas,"35.9",false,true,false,0.05,0);

  canvas->RedrawAxis("samesaxis");

  TLatex * tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.030);
  tex->Draw();
  if (coupling == "1")
    tex->DrawLatex(0.225,0.81,"#bf{Axial med, Dirac DM, g_{q} = 1, g_{DM} = 1}");
  else
    tex->DrawLatex(0.225,0.81,"#bf{Axial med, Dirac DM, g_{q} = 0.25, g_{DM} = 1}");
  ///////
  canvas->SaveAs((outputDirectory+"/scanDD_axial_g"+coupling+"_"+energy+"TeV_v1.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/scanDD_axial_g"+coupling+"_"+energy+"TeV_v1.png").c_str(),"png");

  if(saveOutputFile){

    TFile*outfile = new TFile((outputDirectory+"/axial_g"+coupling+"_DD.root").c_str(),"RECREATE");
    hobs2->Write("contour_obs");
    hexp2->Write("contour_exp");
    lTotalE->Write("contour_exp_graph");
    lTotal->Write("contour_obs_graph");
    DDE_graph->Write("expected_dd");
    DD_graph->Write("observed_dd");
    outfile->Write();
    outfile->Close();

  } 
   
}

TGraph *IceCubett() {
  int i0 = -1;
  double *lX = new double[1000];
  double *lY = new double[1000];
  i0++; lX[i0] = 174.466541859227;   lY[i0] = 8.183380979327122e-40;
  i0++; lX[i0] = 206.30523955353533; lY[i0] = 4.305154929285435e-40;
  i0++; lX[i0] = 246.49284029453307; lY[i0] = 2.635100144816651e-40;
  i0++; lX[i0] = 313.5434765255382;  lY[i0] = 1.807383107740827e-40;
  i0++; lX[i0] = 402.97353212521045; lY[i0] = 1.4979161553149364e-40;
  i0++; lX[i0] = 557.1360794182834;  lY[i0] = 1.5584397016725806e-40;
  i0++; lX[i0] = 754.2505792137222;  lY[i0] = 1.9588542201429473e-40;
  i0++; lX[i0] = 1098.547008425384;  lY[i0] = 2.6567578744751097e-40;
  i0++; lX[i0] = 1616.5765607149508; lY[i0] = 4.52186814732522e-40;
  i0++; lX[i0] = 2142.9379657265067; lY[i0] = 7.131345681431886e-40;
  i0++; lX[i0] = 2585.8767796981792; lY[i0] = 1.003473137983054e-39;
  i0++; lX[i0] = 3287.669923442391;  lY[i0] = 1.466885365962837e-39;
  i0++; lX[i0] = 4050.7141559436727; lY[i0] = 2.312468131334272e-39;
  i0++; lX[i0] = 5043.409155143994;  lY[i0] = 3.510333100892392e-39;
  i0++; lX[i0] = 6826.42526365628;   lY[i0] = 5.971939160098973e-39;
  i0++; lX[i0] = 8678.867703725591;  lY[i0] = 9.066458245538318e-39;
  i0++; lX[i0] = 9633.750991509341;  lY[i0] = 1.0960884198507952e-38;
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
  
  double *lX = new double[100];
  double *lY = new double[100];
  int i0 = -1;
  i0++; lX[i0] = 3.8931579205245446; lY[i0] = 9.623813407451428e-38;
  i0++; lX[i0] = 3.970282087707827;  lY[i0] = 4.383842019068031e-38;
  i0++; lX[i0] = 4.199820413357338;  lY[i0] = 1.5949524727624943e-38;
  i0++; lX[i0] = 4.482135155497945;  lY[i0] = 6.372032293347673e-39;
  i0++; lX[i0] = 4.826318252605578;  lY[i0] = 2.6926833809088005e-39;
  i0++; lX[i0] = 5.339357267403779;  lY[i0] = 1.1811715112423334e-39;
  i0++; lX[i0] = 5.906064574700994;  lY[i0] = 5.584202759858043e-40;
  i0++; lX[i0] = 6.8346244281270705; lY[i0] = 2.6895089772196782e-40;
  i0++; lX[i0] = 8.348935491889968;  lY[i0] = 1.3699205674505805e-40;
  i0++; lX[i0] = 9.744617813963378;  lY[i0] = 8.415611126903263e-41;
  i0++; lX[i0] = 12.335038810947388; lY[i0] = 5.67553560765733e-41;
  i0++; lX[i0] = 16.04081164050506;  lY[i0] = 4.124854834082995e-41;
  i0++; lX[i0] = 22.215547350883014; lY[i0] = 3.481436069600088e-41;
  i0++; lX[i0] = 32.7678925242856;   lY[i0] = 3.349083563571655e-41;
  i0++; lX[i0] = 45.77331698236027;  lY[i0] = 3.538539580170524e-41;
  i0++; lX[i0] = 64.51149188902427;  lY[i0] = 4.029297802376817e-41;
  i0++; lX[i0] = 87.68299154571334;  lY[i0] = 4.853760102012332e-41;
  i0++; lX[i0] = 121.34641794027017; lY[i0] = 5.957043490614464e-41;
  i0++; lX[i0] = 172.52366371970632; lY[i0] = 7.87886858621807e-41;
  i0++; lX[i0] = 228.20287583832362; lY[i0] = 9.854009733595549e-41;
  i0++; lX[i0] = 310.11274632075424; lY[i0] = 1.3035047088360936e-40;
  i0++; lX[i0] = 421.4383712227371;  lY[i0] = 1.6923184390371026e-40;
  i0++; lX[i0] = 572.6860771785143;  lY[i0] = 2.2809293799712587e-40;
  i0++; lX[i0] = 771.2416717396434;  lY[i0] = 3.017343563226607e-40;
  i0++; lX[i0] = 975.042995639776;   lY[i0] = 3.845657120582829e-40;
  TGraph *lGraph = new TGraph(i0,lX,lY);
  lGraph->SetLineWidth(3);
  lGraph->GetXaxis()->SetTitle("m_{DM}");
  lGraph->GetYaxis()->SetTitle("#sigma_{SD} cm^{2}");
  return lGraph;
}


TGraph *PicassoFinalGraph() {
  double *lX = new double[100];
  double *lY = new double[100];
  int i0 = -1;
  i0++; lX[i0] = 2.2847229213303226; lY[i0] = 2.8119869247447238e-37;
  i0++; lX[i0] = 2.508202006389139; lY[i0] = 1.6579080898297997e-37;
  i0++; lX[i0] = 2.9895074543774984; lY[i0] = 9.989103610795583e-38;
  i0++; lX[i0] = 3.5632885402209293; lY[i0] = 5.88760127871836e-38;
  i0++; lX[i0] = 4.517014276660886; lY[i0] = 3.706050949516562e-38;
  i0++; lX[i0] = 6.3462785304309985; lY[i0] = 2.281188037044168e-38;
  i0++; lX[i0] = 7.878559141104886; lY[i0] = 1.7894099173231214e-38;
  i0++; lX[i0] = 10.619323432682782; lY[i0] = 1.40321207716741e-38;
  i0++; lX[i0] = 13.594509044559207; lY[i0] = 1.201816367022604e-38;
  i0++; lX[i0] = 17.944301377415808; lY[i0] = 1.2005535660497427e-38;
  i0++; lX[i0] = 25.193085761636226; lY[i0] = 1.199011944434211e-38;
  i0++; lX[i0] = 35.728612307775904; lY[i0] = 1.3663870293622583e-38;
  i0++; lX[i0] = 46.66900631153654; lY[i0] = 1.5576118143871882e-38;
  i0++; lX[i0] = 64.83209873748544; lY[i0] = 1.8962948079849698e-38;
  i0++; lX[i0] = 92.87324428048362; lY[i0] = 2.5206846980755147e-38;
  i0++; lX[i0] = 129.0058728576599; lY[i0] = 3.2781417599761544e-38;
  i0++; lX[i0] = 177.35139452114217; lY[i0] = 4.4551490024694034e-38;
  i0++; lX[i0] = 259.30378461665174; lY[i0] = 6.466329909914927e-38;
  i0++; lX[i0] = 356.4906981050408; lY[i0] = 8.596826614474313e-38;
  i0++; lX[i0] = 505.4229981253147; lY[i0] = 1.1941974371322467e-37;
  i0++; lX[i0] = 701.8985085449751; lY[i0] = 1.8116099526783597e-37;
  i0++; lX[i0] = 974.8785668290182; lY[i0] = 2.5167291837581657e-37;
  TGraph *lGraph = new TGraph(i0,lX,lY);
  lGraph->SetLineWidth(3);
  return lGraph;
}

///////
TGraph* IceCubebb() {

  double *lX = new double[100];
  double *lY = new double[100];
  int i0 = -1;
  i0++; lX[i0] = 34.42297015484344; lY[i0] = 1.3987131026472186e-38;
  i0++; lX[i0] = 43.52053661323159; lY[i0] = 1.0353218432956531e-38;
  i0++; lX[i0] = 53.86645302624621; lY[i0] = 7.316807143427178e-39;
  i0++; lX[i0] = 71.82134479782309; lY[i0] = 5.672426068491896e-39;
  i0++; lX[i0] = 98.87867356074929; lY[i0] = 4.008806328898432e-39;
  i0++; lX[i0] = 157.95279137310797; lY[i0] = 3.739937302478771e-39;
  i0++; lX[i0] = 247.01231204435697; lY[i0] = 3.4092850697467935e-39;
  i0++; lX[i0] = 467.75874146675534; lY[i0] = 3.739937302478771e-39;
  i0++; lX[i0] = 975.192358771805; lY[i0] = 2.96730240818886e-39;
  i0++; lX[i0] = 1713.229154611073; lY[i0] = 5.052631065335618e-39;
  // i0++; lX[i0] = 2978.114255518353; lY[i0] = 8.214343584919388e-39;
  //  i0++; lX[i0] = 4319.202043243905; lY[i0] = 1.5702901247293776e-38;
  //  i0++; lX[i0] = 6748.0252680719; lY[i0] = 3.292971255097108e-38;
  //  i0++; lX[i0] = 9580.625677830956; lY[i0] = 6.294988990221836e-38;
  TGraph *lGraph = new TGraph(i0,lX,lY);
  lGraph->SetLineWidth(3);
  return lGraph;
}

TGraph* SuperKbb() {
  double *lX = new double[100];
  double *lY = new double[100];
  int i0 = -1;
  i0++; lX[i0] = 5.959274330558362; lY[i0] = 1.702769172225884e-39;
  i0++; lX[i0] = 9.31957887046634; lY[i0] = 1.5167168884709118e-39;
  i0++; lX[i0] = 14.266563133451648; lY[i0] = 1.4481182276745241e-39;
  i0++; lX[i0] = 19.843960665188213; lY[i0] = 1.4149912974345673e-39;
  i0++; lX[i0] = 32.36945765104957; lY[i0] = 1.825183494319024e-39;
  i0++; lX[i0] = 55.684168678607925; lY[i0] = 2.4658110758226038e-39;
  i0++; lX[i0] = 95.79430086580837; lY[i0] = 3.255088599835043e-39;
  i0++; lX[i0] = 145.05988741704272; lY[i0] = 3.739937302478771e-39;
  i0++; lX[i0] = 295.77687775131477; lY[i0] = 6.669919663030115e-39;
  i0++; lX[i0] = 536.4267393150936; lY[i0] = 1.245883364294993e-38;
  i0++; lX[i0] = 972.9520978888768; lY[i0] = 2.171117945694492e-38;
  i0++; lX[i0] = 1800.5552913529386; lY[i0] = 1.0473708979594488e-37;
  //i0++; lX[i0] = 4533.899422504216; lY[i0] = 9.221978823334266e-37;
  TGraph *lGraph = new TGraph(i0,lX,lY);
  lGraph->SetLineWidth(3);
  return lGraph;
}
