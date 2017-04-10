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

double vecF(double mMED,double mDM){    
    double mR = (0.939*mDM)/(0.939+mDM);
    double c = 6.9e-41*1e12;
    return c*(mR*mR)/(mMED*mMED*mMED*mMED);
}

TGraph * makeOBV(TGraph *Graph1){
    TGraph *gr = new TGraph();
    double X;
    double Y;
    int pp=0;
    Graph1->GetPoint(1,X,Y);
    for (double MDM = 1; MDM < Y; MDM += 0.1){
        gr->SetPoint(pp,MDM,vecF(X,MDM));
        pp++;
    }
    for (int p =0; p < Graph1->GetN(); p++){
        Graph1->GetPoint(p,X,Y);
        if (!(X >0)) continue;
        if (!(Y >0)) continue;
        gr->SetPoint(pp,Y,vecF(X,Y));
        pp++;
    }
    gr->SetName(Form("%s_DD",Graph1->GetName()));
    gr->SetLineStyle(Graph1->GetLineStyle());
    gr->SetLineColor(Graph1->GetLineColor());
    gr->SetLineWidth(Graph1->GetLineWidth());
    
    return gr;
}

/////                                                                                                                                                                                                  
static float nbinsX = 800;
static float nbinsY = 500;
static float minX = 100;
static float minY = 0;
static float maxX = 5000;
static float maxY = 1500;
static float maxZ = 10;

static float minX_dd = 1;
static float maxX_dd = 1400;
static double minY_dd = 5e-47;
static double maxY_dd = 5e-34;

static bool saveOutputFile = false;

TGraph* superCDMS();
TGraph* lux();
TGraph* cdmslite();
TGraph* panda();
TGraph* cresst();

void plotVector_DD (string inputFileName, string outputDirectory, string coupling = "025", string energy = "13") {

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
  
  int expcounter = 0;
  int obscounter = 0;

  for (int i = 0; i < tree->GetEntries(); i++){
    
    tree->GetEntry(i);
    
    if (quantile != 0.5 && quantile != -1) continue;
    int c       = code(mh);
    int medmass = mmed(mh,c);
    int dmmass  = mdm(mh,c);
    //for cosmetic reasons
    if(medmass >  1400 and dmmass >= 150 and dmmass <= 350) continue;
    if(medmass == 2500 and dmmass <= 50) continue;

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
  TH2D* hexp = new TH2D("hexp", "", nbinsX, minX, maxX, nbinsY, minY, maxY);
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
    }
  }

  hexp->Smooth();
  hobs->Smooth();
  
  TH2* hexp2 = (TH2*)hexp->Clone("hexp2");
  TH2* hobs2 = (TH2*)hobs->Clone("hobs2");
  
  hexp2->SetContour(2);
  hexp2->SetContourLevel(1, 1);
  hobs2->SetContour(2);
  hobs2->SetContourLevel(1, 1);
  
  hexp2->Draw("contz list");
  gPad->Update();
  
  TObjArray *lContoursE = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  std::vector<double> lXE;
  std::vector<double> lYE;
  int lTotalContsE = lContoursE->GetSize();
  for(int i0 = 0; i0 < lTotalContsE; i0++){
    TList * pContLevel = (TList*)lContoursE->At(i0);
    TGraph *pCurv = (TGraph*)pContLevel->First();
    for(int i1 = 0; i1 < pContLevel->GetSize(); i1++){
      for(int i2  = 0; i2 < pCurv->GetN(); i2++) {
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
  lTotalE->SetLineColor(1);
  lTotalE->SetLineWidth(3);

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
      pCurv->SetLineColor(kRed);
            pCurv = (TGraph*)pContLevel->After(pCurv);
    }
  }
  if(lX.size() == 0) {
    lX.push_back(0); 
    lY.push_back(0); 
  }

  TGraph *lTotal = new TGraph(lX.size(),&lX[0],&lY[0]);  
  lTotal->SetLineColor(1);
  lTotal->SetLineWidth(3);
  
  TGraph *DDE_graph = makeOBV(lTotalE);
  TGraph *DD_graph  = makeOBV(lTotal);

  TCanvas* canvas = new TCanvas("canvas","canvas",750,600);
  canvas->SetLogx();
  canvas->SetLogy();

  TH1* frame = canvas->DrawFrame(minX_dd,minY_dd,maxX_dd,maxY_dd,"");
  frame->GetYaxis()->SetTitle("#sigma^{SI}_{DM-nucleon} [cm^{2}]");
  frame->GetXaxis()->SetTitle("m_{DM} [GeV]");
  frame->GetXaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetXaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.045);
  frame->GetYaxis()->SetTitleOffset(1.25);
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->CenterTitle();
  frame->Draw();

  TGraph *lM0 = lux();
  TGraph *lM1 = cdmslite();
  TGraph *lM2 = panda();
  TGraph *lM3 = cresst();

  lM0->SetLineColor(kBlue);
  lM1->SetLineColor(kBlue+2);
  lM2->SetLineColor(kAzure+1);
  lM3->SetLineColor(kAzure+8);

  lM0->Draw("L SAME");
  lM1->Draw("L SAME");
  lM2->Draw("L SAME");
  lM3->Draw("L SAME");

  DDE_graph->SetLineColor(kRed);
  DD_graph->SetLineColor(kBlack);
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
  leg->AddEntry(lM0 ,"LUX","L");
  leg->AddEntry(lM1 ,"CDMSLite","L");
  leg->AddEntry(lM2 ,"Panda-X II","L");
  leg->AddEntry(lM3 ,"CRESST-II","L");
  leg->Draw("SAME");

  CMS_lumi(canvas,"35.9",false,true,false,0,-0.22);


  TLatex * tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.030);
  tex->Draw();
  if (coupling == "1"){
    tex->DrawLatex(0.75,0.82,"#bf{Vector med, Dirac DM,}");
    tex->DrawLatex(0.75,0.78,"#bf{g_{q} = 1, g_{DM} = 1}");
  }
  else{
    tex->DrawLatex(0.75,0.82,"#bf{Vector med, Dirac DM,}");
    tex->DrawLatex(0.75,0.78,"#bf{g_{q} = 0.25, g_{DM} = 1}");
  }
    
  ///////                                                                                                                                                                                             
  canvas->SaveAs((outputDirectory+"/scanDD_vector_g"+coupling+"_"+energy+"TeV_v1.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDirectory+"/scanDD_vector_g"+coupling+"_"+energy+"TeV_v1.png").c_str(),"png");
  
  if(saveOutputFile){

    TFile*outfile = new TFile(("vector_g"+coupling+"_DD.root").c_str(),"RECREATE");
    DDE_graph->SetName("expected");
    DD_graph->SetName("observed");
    DDE_graph->Write();
    DD_graph->Write();
    outfile->Write();
    outfile->Close();
  }

}

TGraph *panda(){

  int i0 = -1;
  double *lX = new double[1000];
  double *lY = new double[1000];
  i0++; lX[i0] =4.970 ; lY[i0]= 9.731e-43;
  i0++; lX[i0] =6.285 ; lY[i0]= 9.376e-44;
  i0++; lX[i0] =7.717 ; lY[i0]= 1.870e-44;
  i0++; lX[i0] =10.173 ; lY[i0]= 3.748e-45;
  i0++; lX[i0] =14.394 ; lY[i0]= 9.809e-46;
  i0++; lX[i0] =16.305 ; lY[i0]= 6.979e-46;
  i0++; lX[i0] =19.347 ; lY[i0]= 4.712e-46;
  i0++; lX[i0] =23.526 ; lY[i0]= 3.405e-46;
  i0++; lX[i0] =34.701 ; lY[i0]= 2.511e-46;
  i0++; lX[i0] =42.506 ; lY[i0]= 2.420e-46;
  i0++; lX[i0] =61.783 ; lY[i0]= 2.713e-46;
  i0++; lX[i0] =90.243 ; lY[i0]= 3.430e-46;
  i0++; lX[i0] =126.139 ; lY[i0]= 4.406e-46;
  i0++; lX[i0] =185.600 ; lY[i0]= 5.932e-46;
  i0++; lX[i0] =321.685 ; lY[i0]= 9.791e-46;
  i0++; lX[i0] =457.401 ; lY[i0]= 1.368e-45;
  i0++; lX[i0] =701.572 ; lY[i0]= 2.188e-45;
  i0++; lX[i0] =990.270 ; lY[i0]= 3.323e-45;
  
  TGraph *lLimit = new TGraph(i0,lX,lY);
  lLimit->SetLineWidth(3.);
  return lLimit;
}

TGraph *cdmslite(){
  
  int i0 = -1;
  double *lX = new double[1000];
  double *lY = new double[1000];
  i0++; lX[i0] =1.429 ; lY[i0]= 9.880e-38;
  i0++; lX[i0] =1.473 ; lY[i0]= 3.162e-38;
  i0++; lX[i0] =1.574 ; lY[i0]= 1.075e-38;
  i0++; lX[i0] =1.654 ; lY[i0]= 4.334e-39;
  i0++; lX[i0] =1.731 ; lY[i0]= 1.924e-39;
  i0++; lX[i0] =1.838 ; lY[i0]= 8.338e-40;
  i0++; lX[i0] =1.986 ; lY[i0]= 3.279e-40;
  i0++; lX[i0] =2.169 ; lY[i0]= 1.566e-40;
  i0++; lX[i0] =2.384 ; lY[i0]= 8.338e-41;
  i0++; lX[i0] =2.632 ; lY[i0]= 5.012e-41;
  i0++; lX[i0] =2.857 ; lY[i0]= 3.527e-41;
  i0++; lX[i0] =3.093 ; lY[i0]= 2.701e-41;
  i0++; lX[i0] =3.415 ; lY[i0]= 2.120e-41;
  i0++; lX[i0] =3.827 ; lY[i0]= 1.789e-41;
  i0++; lX[i0] =4.207 ; lY[i0]= 1.566e-41;
  i0++; lX[i0] =4.625 ; lY[i0]= 1.492e-41;
  i0++; lX[i0] =5.106 ; lY[i0]= 1.438e-41;
  i0++; lX[i0] =5.613 ; lY[i0]= 1.321e-41;
  i0++; lX[i0] =6.157 ; lY[i0]= 1.143e-41;
  i0++; lX[i0] =6.768 ; lY[i0]= 9.527e-42;
  i0++; lX[i0] =7.504 ; lY[i0]= 8.040e-42;
  i0++; lX[i0] =8.161 ; lY[i0]= 7.122e-42;
  i0++; lX[i0] =9.009 ; lY[i0]= 6.387e-42;
  i0++; lX[i0] =9.968 ; lY[i0]= 5.796e-42;
  i0++; lX[i0] =11.100 ; lY[i0]= 5.261e-42;
  i0++; lX[i0] =12.307 ; lY[i0]= 5.073e-42;
  i0++; lX[i0] =13.645 ; lY[i0]= 4.775e-42;
  i0++; lX[i0] =14.871 ; lY[i0]= 4.892e-42;
  
  
  TGraph *lLimit = new TGraph(i0,lX,lY);
  lLimit->SetLineWidth(3.);
  return lLimit;
  
}


TGraph *superCDMS() {
  int i0 = -1;
  double *lX = new double[1000];
  double *lY = new double[1000];
  i0++; lX[i0] = 3.5946953342351033; lY[i0]= 9.848667679721256e-40;
  i0++; lX[i0] = 3.75095871431554;   lY[i0] = 4.7369269651270765e-40;
  i0++; lX[i0] = 4.055296786109914;  lY[i0] = 1.731452556791703e-40;
  i0++; lX[i0] = 4.859185478655973;  lY[i0] = 3.335637195211155e-41;
  i0++; lX[i0] = 5.443071002530319;  lY[i0] = 1.4979478319670555e-41;
  i0++; lX[i0] = 6.568479751568511;  lY[i0] = 4.257347230887382e-42;
  i0++; lX[i0] = 7.759712545805536;  lY[i0] = 1.4200965776277166e-42;
  i0++; lX[i0] = 10.195931085528889; lY[i0] = 5.073387251235342e-43;
  i0++; lX[i0] = 13.636666866311272; lY[i0] = 2.930133564061426e-43;
  i0++; lX[i0] = 18.04552501288121;  lY[i0] = 2.0321090962591862e-43;
  i0++; lX[i0] = 23.377104302070727; lY[i0] = 1.941241066769181e-43;
  i0++; lX[i0] = 29.857390063416926; lY[i0] = 2.0321090962591862e-43;
  TGraph *lLimit = new TGraph(i0,lX,lY);
  lLimit->SetLineWidth(3.);
  return lLimit;
}


TGraph *cresst(){

  int i0 = -1;
  double *lX = new double[1000];
  double *lY = new double[1000];
    
  i0++; lX[i0] = 0.502; lY[i0] =  1.929e-36;
  i0++; lX[i0] = 0.535; lY[i0] =  1.164e-36;
  i0++; lX[i0] = 0.566; lY[i0] =  6.032e-37;
  i0++; lX[i0] = 0.603; lY[i0] =  2.755e-37;
  i0++; lX[i0] = 0.638; lY[i0] =  1.428e-37;
  i0++; lX[i0] = 0.675; lY[i0] =  7.402e-38;
  i0++; lX[i0] = 0.720; lY[i0] =  4.036e-38;
  i0++; lX[i0] = 0.767; lY[i0] =  2.257e-38;
  i0++; lX[i0] = 0.811; lY[i0] =  1.468e-38;
  i0++; lX[i0] = 0.879; lY[i0] =  9.317e-39;
  i0++; lX[i0] = 0.942; lY[i0] =  6.708e-39;
  i0++; lX[i0] = 1.004; lY[i0] =  4.953e-39;
  i0++; lX[i0] = 1.106; lY[i0] =  3.390e-39;
  i0++; lX[i0] = 1.242; lY[i0] =  2.206e-39;
  i0++; lX[i0] = 1.405; lY[i0] =  1.629e-39;
  i0++; lX[i0] = 1.517; lY[i0] =  1.331e-39;
  i0++; lX[i0] = 1.727; lY[i0] =  1.087e-39;
  i0++; lX[i0] = 2.000; lY[i0] =  8.442e-40;
  i0++; lX[i0] = 2.159; lY[i0] =  7.073e-40;
  i0++; lX[i0] = 2.339; lY[i0] =  5.223e-40;
  i0++; lX[i0] = 2.483; lY[i0] =  3.760e-40;
  i0++; lX[i0] = 2.992; lY[i0] =  1.237e-40;
  i0++; lX[i0] = 3.285; lY[i0] =  6.914e-41;
  i0++; lX[i0] = 3.523; lY[i0] =  4.978e-41;
  i0++; lX[i0] = 4.273; lY[i0] =  3.079e-41;
  i0++; lX[i0] = 5.166; lY[i0] =  2.274e-41;
  i0++; lX[i0] = 6.184; lY[i0] =  1.637e-41;
  i0++; lX[i0] = 7.703; lY[i0] =  1.209e-41;
  i0++; lX[i0] = 10.393; lY[i0] =  1.065e-41;
  i0++; lX[i0] = 13.654; lY[i0] =  1.149e-41;
  i0++; lX[i0] = 17.466; lY[i0] =  1.304e-41;
  i0++; lX[i0] = 23.330; lY[i0] =  1.637e-41;
  i0++; lX[i0] = 26.476; lY[i0] =  1.480e-41;
  i0++; lX[i0] = 29.748; lY[i0] =  9.628e-42;
     
  TGraph *lLimit = new TGraph(i0,lX,lY);
  lLimit->SetLineWidth(3.);
  return lLimit;

}


TGraph *lux(){

  int i0 = -1;
  double *lX = new double[1000];
  double *lY = new double[1000];
  
  i0++; lX[i0] = 5.230202675339297; lY[i0] =    9.999999999999836e-43;
  i0++; lX[i0] = 5.320175096324759; lY[i0] =    7.906043210907605e-43;
  i0++; lX[i0] = 5.411695265464636; lY[i0] =      6.105402296585314e-43;
  i0++; lX[i0] = 5.50478980785497; lY[i0] =       4.4984326689693795e-43;
  i0++; lX[i0] = 5.599485806609297; lY[i0] =     3.393221771895274e-43;
  i0++; lX[i0] = 5.695810810737687; lY[i0] =     2.6203985288583326e-43;
  i0++; lX[i0] = 5.793792843161313; lY[i0] =     2.0716983998953036e-43;
  i0++; lX[i0] = 5.893460408864902; lY[i0] =     1.5264179671752116e-43;
  i0++; lX[i0] = 5.99484250318941; lY[i0] =      1.2067926406393165e-43;
  i0++; lX[i0] = 6.202868761603334; lY[i0] =     8.891593339164587e-44;
  i0++; lX[i0] = 6.309573444801932; lY[i0] =     6.707035611184258e-44;
  i0++; lX[i0] = 6.418113712446068; lY[i0] =     4.941713361323818e-44;
  i0++; lX[i0] = 6.6408278506348415; lY[i0] =    3.556480306223092e-44;
  i0++; lX[i0] = 6.989473207273485; lY[i0] =     2.500110382617941e-44;
  i0++; lX[i0] = 7.232014017781835; lY[i0] =     1.7575106248547966e-44;
  i0++; lX[i0] = 7.611696816847943; lY[i0] =     1.0481131341546917e-44;
  i0++; lX[i0] = 8.011313071180073; lY[i0] =     7.19685673001147e-45;
  i0++; lX[i0] = 8.576958985908941; lY[i0] =     4.941713361323818e-45;
  i0++; lX[i0] = 9.027251779484576; lY[i0] =     3.641031949310677e-45;
  i0++; lX[i0] = 9.501185073181436; lY[i0] =     2.7464741148160504e-45;
  i0++; lX[i0] = 10.347008713411984; lY[i0] =    1.8858632787726555e-45;
  i0++; lX[i0] = 11.268130079648545; lY[i0] =    1.2354828882567482e-45;
  i0++; lX[i0] = 12.697075549188032; lY[i0] =    7.722449945836333e-46;
  i0++; lX[i0] = 14.803703235666655; lY[i0] =    4.826957437677881e-46;
  i0++; lX[i0] = 16.398903672222477; lY[i0] =    3.4738921120831546e-46;
  i0++; lX[i0] = 19.448624389373624; lY[i0] =    2.4420530945486748e-46;
  i0++; lX[i0] = 22.67543125870801; lY[i0] =     1.930697728883277e-46;
  i0++; lX[i0] = 26.892404165083185; lY[i0] =    1.4563484775012623e-46;
  i0++; lX[i0] = 31.89361179185746; lY[i0] =     1.2354828882567684e-46;
  i0++; lX[i0] = 39.13745601980384; lY[i0] =     1.0481131341546917e-46;
  i0++; lX[i0] = 44.85923383437637; lY[i0] =     1.073030940526174e-46;
  i0++; lX[i0] = 57.937928431613145; lY[i0] =    1.2067926406393363e-46;
  i0++; lX[i0] = 78.75829328777596; lY[i0] =     1.422529313485372e-46;
  i0++; lX[i0] = 100; lY[i0] =                   1.6768329368110306e-46;
  i0++; lX[i0] = 129.1549665014884; lY[i0] =     2.0235896477251638e-46;
  i0++; lX[i0] = 161.21572750178856; lY[i0] =    2.4420530945486748e-46;
  i0++; lX[i0] = 204.696827180752; lY[i0] =      3.0888435964774975e-46;
  i0++; lX[i0] = 264.3761185749101; lY[i0] =     3.906939937054621e-46;
  i0++; lX[i0] = 335.68035509467273; lY[i0] =    5.059197488435864e-46;
  i0++; lX[i0] = 426.2158829015325; lY[i0] =     6.250551925274027e-46;
  i0++; lX[i0] = 532.0175096324763; lY[i0] =     7.906043210907734e-46;
  i0++; lX[i0] = 652.8521141127848; lY[i0] =     1e-45;
  i0++; lX[i0] = 787.5829328777596; lY[i0] =     1.2067926406393363e-45;
  i0++; lX[i0] = 983.0884473994828; lY[i0] =     1.4909716571840752e-45;
  i0++; lX[i0] = 1248.2348288165172; lY[i0] =    1.976598071701651e-45;
  i0++; lX[i0] = 1584.893192461114; lY[i0] =     2.500110382617941e-45;
  i0++; lX[i0] = 1978.3188827841623; lY[i0] =    3.088843596477498e-45;
  i0++; lX[i0] = 2386.5897868585835; lY[i0] =    3.727593720314923e-45;
  i0++; lX[i0] = 2879.116638022353; lY[i0] =     4.498432668969453e-45;
  i0++; lX[i0] = 3593.813663804626; lY[i0] =     5.557736586486892e-45;
  i0++; lX[i0] = 4485.923383437636; lY[i0] =     6.866488450043026e-45;
  i0++; lX[i0] = 5994.8425031894085; lY[i0] =    9.102981779915264e-45;
  i0++; lX[i0] = 7875.829328777596; lY[i0] =     1.2354828882567482e-44;
  i0++; lX[i0] = 11461.96978456578; lY[i0] =     1.7575106248547966e-44;
  i0++; lX[i0] = 16967.959918688968; lY[i0] =   2.6203985288583324e-44;
  i0++; lX[i0] = 26437.61185749096; lY[i0] =    4.0949150623803854e-44;
  i0++; lX[i0] = 35330.366949927295; lY[i0] =   5.557736586486892e-44;
  i0++; lX[i0] = 48026.56010546725; lY[i0] =    7.543120063354607e-44;
  i0++; lX[i0] = 68712.70363478754; lY[i0] =    1.0730309405261566e-43;
  i0++; lX[i0] = 98308.84473994808; lY[i0] =    1.4909716571840507e-43;
    
  TGraph *lLimit = new TGraph(i0,lX,lY);
  lLimit->SetLineWidth(3.);
  return lLimit;
    
}
