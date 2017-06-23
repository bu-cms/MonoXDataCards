#include "../CMS_lumi.h"
#include "externalLibs/RooSplineND.h"

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

TGraph* produceContour (const int & reduction, const bool & applyReduction = true){

  TObjArray *lContoursE = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  std::vector<double> lXE;
  std::vector<double> lYE;
  int lTotalContsE = lContoursE->GetSize();
  for(int i0 = 0; i0 < lTotalContsE; i0++){
    TList * pContLevel = (TList*)lContoursE->At(i0);
    TGraph *pCurv = (TGraph*)pContLevel->First();
    for(int i1 = 0; i1 < pContLevel->GetSize(); i1++){
      for(int i2  = 0; i2 < pCurv->GetN(); i2++) {
        if(applyReduction && i2%reduction != 0) continue; // reduce number of points                                                                                                                 
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
static float nbinsX      = 1000;
static float minX        = 1;
static float minY        = 1;
static float maxX        = 2000;
static float maxY        = 700;
static float nbinsY      = 350;
static float nCoupling   = 100;
static float minCoupling = 0.02;
static float maxCoupling = 1;
static float minZ        = 0.1;
static float maxZ        = 100;
static int   reductionForContour = 12;
static bool  addPreliminary   = false;
static bool  addRelicDensity  = true;
static int   nForInterpolateX = 70;
static int   nForInterpolateY = 55;
//static int   nForInterpolateX = 35;
//static int   nForInterpolateY = 35;
static float minCoupling_spline = 0.01;
static float maxCoupling_spline = 1.0;
static float minX_spline = 1;
static float maxX_spline = 2200;
static float minY_spline = 0.5;
static float maxY_spline = 800;
static float maxXForBrazilianY = 1350;
static float maxXForBrazilianX = 1650;
static float maxXForBrazilianY_DM = 410;
static float maxXForBrazilianX_DM = 520;

void fillLimitGraphs(TTree* tree,
		     TGraph* grexp, TGraph* grexp_up, TGraph* grexp_dw,
		     TGraph* grobs, TGraph* grobs_up, TGraph* grobs_dw,
		     const float & medOverDM = 3, const bool & useDMMass = false,
		     TGraph* grexp_up2 = NULL, TGraph* grexp_dw2 = NULL){

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

  // identify bad mass points
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
  int exp_up2_counter   = 0;
  int exp_down_counter = 0;
  int exp_down2_counter = 0;
  int obscounter       = 0;

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

    // check / use only points for which mediator and DM mass are in this relation
    if(double(medmass) / double(dmmass) != medOverDM) continue;      

    if (quantile == 0.5) { // expected limit
      if(useDMMass){
	grexp->SetPoint(expcounter, double(dmmass), limit);
      }
      else
	grexp->SetPoint(expcounter, double(medmass), limit);
      expcounter++;
    }

    if (quantile < 0.17 && quantile > 0.14 ) {
      if(useDMMass)
	grexp_dw->SetPoint(exp_down_counter, double(dmmass), limit);
      else
	grexp_dw->SetPoint(exp_down_counter, double(medmass), limit);
      exp_down_counter++;
    }

    if(grexp_dw2 and quantile < 0.09 && quantile > 0){
      if(useDMMass)
	grexp_dw2->SetPoint(exp_down2_counter, double(dmmass), limit);
      else
	grexp_dw2->SetPoint(exp_down2_counter, double(medmass), limit);
      exp_down2_counter++;
    }

    if (quantile < 0.85 && quantile > 0.83 ) {
      if(useDMMass)
	grexp_up->SetPoint(exp_up_counter, double(dmmass), limit);
      else
	grexp_up->SetPoint(exp_up_counter, double(medmass),limit);

      exp_up_counter++;
    }    

    if(grexp_up2 and quantile > 0.9 && quantile < 1){
      if(useDMMass)
	grexp_up2->SetPoint(exp_up2_counter, double(dmmass), limit);
      else
	grexp_up2->SetPoint(exp_up2_counter, double(medmass), limit);
      exp_up2_counter++;
    }
    
    if (quantile == -1) { // observed
      if(useDMMass){
	grobs->SetPoint(obscounter, double(dmmass), limit);
	grobs_up->SetPoint(obscounter, double(dmmass), limit*1.2);
	grobs_dw->SetPoint(obscounter, double(dmmass), limit*0.8);
      }
      else{
	grobs->SetPoint(obscounter, double(medmass), limit);
	grobs_up->SetPoint(obscounter, double(medmass), limit*1.2);
	grobs_dw->SetPoint(obscounter, double(medmass), limit*0.8);
      }
      obscounter++;      
    }
  } 
  tree->ResetBranchAddresses();

}

////////---------------
vector<pair<double,TSpline3*> > make1DSpline( const vector<pair<double,TGraph*> > & graph, const string & postfix){

  vector<pair<double,TSpline3*> > spline;
  for(auto gr : graph){
    TString name = Form("spline%s_gq_%f",postfix.c_str(),gr.first);
    name.ReplaceAll(".","p");
    spline.push_back(pair<double,TSpline3*>(gr.first,new TSpline3(name,gr.second->GetX(),gr.second->GetY(),gr.second->GetN())));
  }

  return spline;
}

////////---------------
void fill2DGraph(TGraph2D* graph2D, const vector<pair<double,TSpline3*> > & spline, const vector<pair<double,TGraph*> > & graph){

  int npoints = 0;  
  int ngraph  = 0;
  for(auto gr : graph){
    int ipoint = 0;
    for(int i = 0; i < gr.second->GetN(); i++){
      double x,y;
      gr.second->GetPoint(ipoint,x,y);
      graph2D->SetPoint(npoints,x,gr.first,spline.at(ngraph).second->Eval(x));
      ipoint++;
      npoints++;
    }
    ngraph++;
  }
}


////////---------------
void fill2DGraph(TGraph2D* graph2D, const vector<pair<double,TSpline3*> > & spline,const vector<pair<double,TGraph*> > & graph, float & minx, float & maxx){

  int npoints = 0;  
  int ngraph  = 0;
  for(auto gr : graph){
    int ipoint = 0;
    for(int i = 0; i < gr.second->GetN(); i++){
      double x,y;
      gr.second->GetPoint(ipoint,x,y);
      graph2D->SetPoint(npoints,x,gr.first,spline.at(ngraph).second->Eval(x));
      if(x < minx) minx = x;
      if(x > maxx) maxx = x;
      ipoint++;
      npoints++;
    }
    ngraph++;
  }
}

////////---------------
RooSplineND* makeSplineND (TH2D* histo, RooRealVar & xvar, RooRealVar & yvar, const string & postfix){

  TTree *tree = new TTree("tree","tree"); 
  float  x,y,z;
  tree->Branch("x",&x,"x/F");
  tree->Branch("y",&y,"y/F");
  tree->Branch("z",&z,"z/F");

  for(int i = 1; i <= histo->GetNbinsX(); i++){
    for(int j = 1; j <= histo->GetNbinsY(); j++){
      x = histo->GetXaxis()->GetBinCenter(i);
      y = histo->GetYaxis()->GetBinCenter(j);
      z = histo->GetBinContent(i,j);
      tree->Fill();
    }
  }

  RooArgList list (xvar,yvar);
  RooSplineND *splineND = new RooSplineND(Form("spline_%s",postfix.c_str()),Form("spline_%s",postfix.c_str()),list,tree,"z",1,true);  
  if(tree) delete tree;
  return splineND;
}

////////---------------
void fixBinContent(TH2D* histo, const float & minZ, const float & maxZ){
  ////////////                                                                                                                                                                                         
  for(int i = 1; i <= histo->GetNbinsX(); i++){
    for(int j = 1; j <= histo->GetNbinsY(); j++){
      if(histo->GetBinContent(i,j) <= 0)   histo->SetBinContent(i,j,maxZ);
      if(histo->GetBinContent(i,j) > maxZ) histo->SetBinContent(i,j,maxZ);
      if(histo->GetBinContent(i,j) < minZ) histo->SetBinContent(i,j,minZ);
    }
  }  
}
  

/////// --------
void makeYGraph(TGraph* graph, TGraph* input){
  int effective_points = 0;
  for(int iPoint = 0; iPoint < input->GetN(); iPoint++){
    double x,y;
    input->GetPoint(iPoint,x,y);
    graph->SetPoint(effective_points,y,x);
    effective_points++;
  }
}

/////// --------
void makeXGraph(TGraph* graph, TGraph* input, TGraph* central){
  int effective_points = 0;
  for(int iPoint = 0; iPoint < central->GetN(); iPoint++){
    double x,y;
    central->GetPoint(iPoint,x,y);
    graph->SetPoint(effective_points,x,input->Eval(x));
    effective_points++;
  }
}


////// ---------
void makeUncertaintyBand(TGraphAsymmErrors* graph, TGraph* central, TGraph* band, TGraph* band_up, TGraph* band_dw, const bool & useDMMass, const int & steps = 100){

  int effective_points = 0;

  for(int iPoint = 0; iPoint < central->GetN(); iPoint++){
    ////////////////
    double x,y;
    double x_up,y_up;
    central->GetPoint(iPoint,x,y);

    if(not useDMMass and x < maxXForBrazilianY) continue;
    else if(useDMMass and x < maxXForBrazilianY_DM) continue;

    central->GetPoint(iPoint+1,x_up,y_up);
    //////////////
    graph->SetPoint(effective_points,x,y);
    ////////////
    if(y_up > y and band_up->Eval(y) < band_dw->Eval(y))
      graph->SetPointError(effective_points,x-band_up->Eval(y),band_dw->Eval(y)-x,0,0);
    else if(y_up < y and band_up->Eval(y) < band_dw->Eval(y))
      graph->SetPointError(effective_points,x-band_up->Eval(y),band_dw->Eval(y)-x,0,0);
    else if(y < y_up and band_up->Eval(y) > band_dw->Eval(y))
      graph->SetPointError(effective_points,x-band_dw->Eval(y),band_up->Eval(y)-x,0,0);
    else if(y > y_up and band_up->Eval(y) > band_dw->Eval(y))
      graph->SetPointError(effective_points,x-band_dw->Eval(y),band_up->Eval(y)-x,0,0);
    
    effective_points++;

    for(int nSteps = 0; nSteps < steps; nSteps++){
      graph->SetPoint(effective_points,band->Eval(y+(y_up-y)/nSteps),y+(y_up-y)/nSteps);
      if(y+(y_up-y)/nSteps > y and band_up->Eval(y) < band_dw->Eval(y))
	graph->SetPointError(effective_points,x-band_up->Eval(y),band_dw->Eval(y)-x,0,y+(y_up-y)/nSteps-y);
      else if(y+(y_up-y)/nSteps < y and band_up->Eval(y) < band_dw->Eval(y))
	graph->SetPointError(effective_points,x-band_up->Eval(y),band_dw->Eval(y)-x,y-y+(y_up-y)/nSteps,0);
      else if(y < y+(y_up-y)/nSteps and band_up->Eval(y) > band_dw->Eval(y))
	graph->SetPointError(effective_points,x-band_dw->Eval(y),band_up->Eval(y)-x,0,y+(y_up-y)/nSteps-y);
      else if(y > y+(y_up-y)/nSteps and band_up->Eval(y) > band_dw->Eval(y))
	graph->SetPointError(effective_points,x-band_dw->Eval(y),band_up->Eval(y)-x,y-y+(y_up-y)/nSteps,0);
      effective_points++;
    }
  }
}


////// ---------
void makeUncertaintyBand(TGraphAsymmErrors* graph, TGraph* central, TGraph* band_up, TGraph* band_dw, const bool & useDMMass, const int & steps = 100){

  int effective_points = 0;

  for(int iPoint = 0; iPoint < central->GetN(); iPoint++){
    ////////////////
    double x,y;
    double x_up,y_up;
    central->GetPoint(iPoint,x,y);

    if(not useDMMass and x > maxXForBrazilianX) continue;
    else if(useDMMass and x > maxXForBrazilianX_DM) continue;

    central->GetPoint(iPoint+1,x_up,y_up);
    //////////////
    graph->SetPoint(effective_points,x,y);

    ////////////
    if(x_up > x and band_up->Eval(x) < band_dw->Eval(x))
      graph->SetPointError(effective_points,0,0,y-band_up->Eval(x),band_dw->Eval(x)-y);
    else if(x_up < x and band_up->Eval(x) < band_dw->Eval(x))
      graph->SetPointError(effective_points,0,0,y-band_up->Eval(x),band_dw->Eval(x)-y);
    else if(x < x_up and band_up->Eval(x) > band_dw->Eval(x))
      graph->SetPointError(effective_points,0,0,y-band_dw->Eval(x),band_up->Eval(x)-y);
    else if(x > x_up and band_up->Eval(x) > band_dw->Eval(x))
      graph->SetPointError(effective_points,0,0,y-band_dw->Eval(x),band_up->Eval(x)-y);
    
    effective_points++;

    for(int nSteps = 0; nSteps < steps; nSteps++){
      graph->SetPoint(effective_points,x+(x_up-x)/nSteps,central->Eval(x+(x_up-x)/nSteps));
      if(x+(x_up-x)/nSteps > x and band_up->Eval(x) < band_dw->Eval(x))
	graph->SetPointError(effective_points,0,x+(x_up-x)/nSteps-x,y-band_up->Eval(x),band_dw->Eval(x)-y);
      else if(x+(x_up-x)/nSteps < x and band_up->Eval(x) < band_dw->Eval(x))
	graph->SetPointError(effective_points,x+(x_up-x)/nSteps-x,0,y-band_dw->Eval(x),band_up->Eval(x)-y);      
      else if(x < x+(x_up-x)/nSteps and band_up->Eval(x) > band_dw->Eval(x))
	graph->SetPointError(effective_points,0,x+(x_up-x)/nSteps-x,y-band_dw->Eval(x),band_up->Eval(x)-y);
      else if(x > x+(x_up-x)/nSteps and band_up->Eval(x) > band_dw->Eval(x))
	graph->SetPointError(effective_points,x-x+(x_up-x)/nSteps,0,y-band_dw->Eval(x),band_up->Eval(x)-y);
      effective_points++;
    }
  }
}

  
///////////////
static vector<float> gq_coupling = {1.0,0.75,0.5,0.3,0.25,0.2,0.1,0.05,0.01};
static bool makeCOLZ = false;

void plotVectorCoupling(string outputDIR, bool useDMMass = false, float medOverDM = 3, bool useSplineND = true, string energy = "13") {

  gSystem->Load("externalLibs/RooSplineND_cc.so");

  TGraph* relic_graph_g1 = new TGraph();
  TGraph* relic_graph_g2 = new TGraph();

  if(addRelicDensity){
    TFile* inputRelicDensity = TFile::Open("externalFiles/relic_v6_V_Res_v2.root","READ");
    TTree* relic = (TTree*) inputRelicDensity->Get("relic");
    TTreeReader reader(relic);
    TTreeReaderValue<float> mmed (reader,"med.lMed");
    TTreeReaderValue<float> mdm  (reader,"dm.lMDM");
    TTreeReaderValue<float> gq1  (reader,"gq1.lGQ1");
    TTreeReaderValue<float> gq2  (reader,"gq2.lGQ2");
    int npoints = 0;
    while(reader.Next()){
      if(fabs(*mmed - medOverDM*(*mdm))/(*mdm) > 0.01) continue;
      if(not useDMMass){
        relic_graph_g1->SetPoint(npoints,*mmed,exp(*gq1));
        relic_graph_g2->SetPoint(npoints,*mmed,exp(*gq2));
      }
      else{
        relic_graph_g1->SetPoint(npoints,*mdm,exp(*gq1));
        relic_graph_g2->SetPoint(npoints,*mdm,exp(*gq2));
      }
      npoints++;
    }
    relic_graph_g1->Sort();
    relic_graph_g2->Sort();
  }

  string inputFileName1 = "LimitsForPaperCoupling/higgsCombine_COMB_Vector_gq_1p0.root";
  string inputFileName2 = "LimitsForPaperCoupling/higgsCombine_COMB_Vector_gq_0p75.root";
  string inputFileName3 = "LimitsForPaperCoupling/higgsCombine_COMB_Vector_gq_0p5.root";
  string inputFileName4 = "LimitsForPaperCoupling/higgsCombine_COMB_Vector_gq_0p3.root";
  string inputFileName5 = "LimitsForPaperCoupling/higgsCombine_COMB_Vector_gq_0p25.root";
  string inputFileName6 = "LimitsForPaperCoupling/higgsCombine_COMB_Vector_gq_0p2.root";
  string inputFileName7 = "LimitsForPaperCoupling/higgsCombine_COMB_Vector_gq_0p1.root";
  string inputFileName8 = "LimitsForPaperCoupling/higgsCombine_COMB_Vector_gq_0p05.root";
  string inputFileName9 = "LimitsForPaperCoupling/higgsCombine_COMB_Vector_gq_0p01.root";

  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  int ncontours = 999;
  gStyle->SetPalette(kBird);  
  gStyle->SetNumberContours(ncontours);
  
  // This is where all the plots are made
  vector<TFile*> fileList;
  fileList.push_back(TFile::Open(inputFileName1.c_str(),"READ"));
  fileList.push_back(TFile::Open(inputFileName2.c_str(),"READ"));
  fileList.push_back(TFile::Open(inputFileName3.c_str(),"READ"));
  fileList.push_back(TFile::Open(inputFileName4.c_str(),"READ"));
  fileList.push_back(TFile::Open(inputFileName5.c_str(),"READ"));
  fileList.push_back(TFile::Open(inputFileName6.c_str(),"READ"));
  fileList.push_back(TFile::Open(inputFileName7.c_str(),"READ"));
  fileList.push_back(TFile::Open(inputFileName8.c_str(),"READ"));
  fileList.push_back(TFile::Open(inputFileName9.c_str(),"READ"));

  vector<TTree*> treeList;
  for(auto file: fileList)
    treeList.push_back((TTree*) file->Get("limit"));
  
  vector<pair<double,TGraph*> > grexp_up;
  vector<pair<double,TGraph*> > grexp_dw;
  vector<pair<double,TGraph*> > grexp_up2;
  vector<pair<double,TGraph*> > grexp_dw2;
  vector<pair<double,TGraph*> > grexp;
  vector<pair<double,TGraph*> > grobs;
  vector<pair<double,TGraph*> > grobs_up;
  vector<pair<double,TGraph*> > grobs_dw;

  for(auto gq : gq_coupling){
    grexp.push_back(pair<double,TGraph*>(gq,new TGraph()));
    grexp_up.push_back(pair<double,TGraph*>(gq,new TGraph()));
    grexp_dw.push_back(pair<double,TGraph*>(gq,new TGraph()));
    grexp_up2.push_back(pair<double,TGraph*>(gq,new TGraph()));
    grexp_dw2.push_back(pair<double,TGraph*>(gq,new TGraph()));
    grobs.push_back(pair<double,TGraph*>(gq,new TGraph()));
    grobs_up.push_back(pair<double,TGraph*>(gq,new TGraph()));
    grobs_dw.push_back(pair<double,TGraph*>(gq,new TGraph()));
  }
  
  for(int itree = 0; itree < treeList.size(); itree++){
    if(makeCOLZ)
      fillLimitGraphs(treeList.at(itree),grexp.at(itree).second,grexp_up.at(itree).second,grexp_dw.at(itree).second,
		      grobs.at(itree).second,grobs_up.at(itree).second,grobs_dw.at(itree).second,
		      medOverDM,useDMMass);
    else
      fillLimitGraphs(treeList.at(itree),grexp.at(itree).second,grexp_up.at(itree).second,grexp_dw.at(itree).second,
		      grobs.at(itree).second,grobs_up.at(itree).second,grobs_dw.at(itree).second,
		      medOverDM,useDMMass,grexp_up2.at(itree).second,grexp_dw2.at(itree).second);

  }
  
  //// make a spline to further smooth                                                                                                                                                          
  vector<pair<double,TSpline3*> > splineexp = make1DSpline(grexp,"exp");
  vector<pair<double,TSpline3*> > splineexp_up = make1DSpline(grexp_up,"exp_up");
  vector<pair<double,TSpline3*> > splineexp_dw = make1DSpline(grexp_dw,"exp_dw");
  vector<pair<double,TSpline3*> > splineexp_up2;
  if(not makeCOLZ) 
    splineexp_up2 = make1DSpline(grexp_up2,"exp_up");
  vector<pair<double,TSpline3*> > splineexp_dw2;
  if(not makeCOLZ)
    splineexp_dw2 = make1DSpline(grexp_dw2,"exp_dw");
  vector<pair<double,TSpline3*> > splineobs    = make1DSpline(grobs,"obs");
  vector<pair<double,TSpline3*> > splineobs_up = make1DSpline(grobs_up,"obs_up");
  vector<pair<double,TSpline3*> > splineobs_dw = make1DSpline(grobs_dw,"obs_dw");

  // fill coupling
  TGraph2D* grexp_coupling = new TGraph2D();
  TGraph2D* grexp_coupling_up = new TGraph2D();
  TGraph2D* grexp_coupling_dw = new TGraph2D();
  TGraph2D* grexp_coupling_up2 = new TGraph2D();
  TGraph2D* grexp_coupling_dw2 = new TGraph2D();
  TGraph2D* grobs_coupling = new TGraph2D();
  TGraph2D* grobs_coupling_up = new TGraph2D();
  TGraph2D* grobs_coupling_dw = new TGraph2D();

  fill2DGraph(grexp_coupling,splineexp,grexp);
  fill2DGraph(grexp_coupling_up,splineexp_up,grexp_up);
  fill2DGraph(grexp_coupling_dw,splineexp_dw,grexp_dw);
  fill2DGraph(grobs_coupling_up,splineobs_up,grobs_up);
  fill2DGraph(grobs_coupling_dw,splineobs_dw,grobs_dw);
  if(not makeCOLZ){
    fill2DGraph(grexp_coupling_up2,splineexp_up2,grexp_up2);
    fill2DGraph(grexp_coupling_dw2,splineexp_dw2,grexp_dw2);    
  }
    
  float min_xaxis = 9999;
  float max_xaxis = -1;
  fill2DGraph(grobs_coupling,splineobs,grobs,min_xaxis,max_xaxis);
   
  ///// --> temp histograms 2D with some binning to perform linear interpolation
  TH2D* hexp_coupling = NULL;
  TH2D* hobs_coupling = NULL;
  TH2D* hexp_coupling_up = NULL;
  TH2D* hexp_coupling_dw = NULL;
  TH2D* hexp_coupling_up2 = NULL;
  TH2D* hexp_coupling_dw2 = NULL;
  TH2D* hobs_coupling_up = NULL;
  TH2D* hobs_coupling_dw = NULL;  

  if(useSplineND){
    if(not useDMMass){
      hexp_coupling = new TH2D("hexp_coupling", "",nForInterpolateX,max(min_xaxis,minX),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      hobs_coupling = new TH2D("hobs_coupling", "",nForInterpolateX,max(min_xaxis,minX),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      hexp_coupling_up = new TH2D("hexp_coupling_up", "",nForInterpolateX,max(min_xaxis,minX),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      hobs_coupling_up = new TH2D("hobs_coupling_up", "",nForInterpolateX,max(min_xaxis,minX),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      hexp_coupling_dw = new TH2D("hexp_coupling_dw", "",nForInterpolateX,max(min_xaxis,minX),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      hobs_coupling_dw = new TH2D("hobs_coupling_dw", "",nForInterpolateX,max(min_xaxis,minX),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      if(not makeCOLZ){
	hexp_coupling_up2 = new TH2D("hexp_coupling_up2", "",nForInterpolateX,max(min_xaxis,minX),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
	hexp_coupling_dw2 = new TH2D("hexp_coupling_dw2", "",nForInterpolateX,max(min_xaxis,minX),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      }
    }
    else{
      hexp_coupling = new TH2D("hexp_coupling", "",nForInterpolateX,max(min_xaxis,minY),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      hobs_coupling = new TH2D("hobs_coupling", "",nForInterpolateX,max(min_xaxis,minY),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      hexp_coupling_up = new TH2D("hexp_coupling_up", "",nForInterpolateX,max(min_xaxis,minY),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      hobs_coupling_up = new TH2D("hobs_coupling_up", "",nForInterpolateX,max(min_xaxis,minY),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      hexp_coupling_dw = new TH2D("hexp_coupling_dw", "",nForInterpolateX,max(min_xaxis,minY),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      hobs_coupling_dw = new TH2D("hobs_coupling_dw", "",nForInterpolateX,max(min_xaxis,minY),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);

      if(not makeCOLZ){
	hexp_coupling_up2 = new TH2D("hexp_coupling_up2", "",nForInterpolateX,max(min_xaxis,minY),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
	hexp_coupling_dw2 = new TH2D("hexp_coupling_dw2", "",nForInterpolateX,max(min_xaxis,minY),max_xaxis,nForInterpolateY,minCoupling,maxCoupling);
      }

    }
  }
  else{
    if(not useDMMass){
      hexp_coupling = new TH2D("hexp_coupling", "",nbinsX,max(min_xaxis,minX),max_xaxis,nCoupling,minCoupling,maxCoupling);
      hobs_coupling = new TH2D("hobs_coupling", "",nbinsX,max(min_xaxis,minX),max_xaxis,nCoupling,minCoupling,maxCoupling);
      hexp_coupling_up = new TH2D("hexp_coupling_up", "",nbinsX,max(min_xaxis,minX),max_xaxis,nCoupling,minCoupling,maxCoupling);
      hobs_coupling_up = new TH2D("hobs_coupling_up", "",nbinsX,max(min_xaxis,minX),max_xaxis,nCoupling,minCoupling,maxCoupling);
      hexp_coupling_dw = new TH2D("hexp_coupling_dw", "",nbinsX,max(min_xaxis,minX),max_xaxis,nCoupling,minCoupling,maxCoupling);
      hobs_coupling_dw = new TH2D("hobs_coupling_dw", "",nbinsX,max(min_xaxis,minX),max_xaxis,nCoupling,minCoupling,maxCoupling);     
      if(not makeCOLZ){
	hexp_coupling_up2 = new TH2D("hexp_coupling_up2", "",nbinsX,max(min_xaxis,minX),max_xaxis,nCoupling,minCoupling,maxCoupling);
	hexp_coupling_dw2 = new TH2D("hexp_coupling_dw2", "",nbinsX,max(min_xaxis,minX),max_xaxis,nCoupling,minCoupling,maxCoupling);
      }
    }
    else{
      hexp_coupling = new TH2D("hexp_coupling", "",nbinsX,max(min_xaxis,minY),max_xaxis,nCoupling,minCoupling,maxCoupling);
      hobs_coupling = new TH2D("hobs_coupling", "",nbinsX,max(min_xaxis,minY),max_xaxis,nCoupling,minCoupling,maxCoupling);
      hexp_coupling_up = new TH2D("hexp_coupling_up", "",nbinsX,max(min_xaxis,minY),max_xaxis,nCoupling,minCoupling,maxCoupling);
      hobs_coupling_up = new TH2D("hobs_coupling_up", "",nbinsX,max(min_xaxis,minY),max_xaxis,nCoupling,minCoupling,maxCoupling);
      hexp_coupling_dw = new TH2D("hexp_coupling_dw", "",nbinsX,max(min_xaxis,minY),max_xaxis,nCoupling,minCoupling,maxCoupling);
      hobs_coupling_dw = new TH2D("hobs_coupling_dw", "",nbinsX,max(min_xaxis,minY),max_xaxis,nCoupling,minCoupling,maxCoupling);
      if(not makeCOLZ){
	hexp_coupling_up2 = new TH2D("hexp_coupling_up2", "",nbinsX,max(min_xaxis,minY),max_xaxis,nCoupling,minCoupling,maxCoupling);
	hexp_coupling_dw2 = new TH2D("hexp_coupling_dw2", "",nbinsX,max(min_xaxis,minY),max_xaxis,nCoupling,minCoupling,maxCoupling);
      }
    }
  }

  for(int i = 1; i <= hexp_coupling->GetNbinsX(); i++){
    for(int j = 1; j <= hexp_coupling->GetNbinsY(); j++){
      hexp_coupling->SetBinContent(i,j,grexp_coupling->Interpolate(hexp_coupling->GetXaxis()->GetBinCenter(i),hexp_coupling->GetYaxis()->GetBinCenter(j)));
      hobs_coupling->SetBinContent(i,j,grobs_coupling->Interpolate(hobs_coupling->GetXaxis()->GetBinCenter(i),hobs_coupling->GetYaxis()->GetBinCenter(j)));
      hexp_coupling_up->SetBinContent(i,j,grexp_coupling_up->Interpolate(hexp_coupling_up->GetXaxis()->GetBinCenter(i),hexp_coupling_up->GetYaxis()->GetBinCenter(j)));
      hexp_coupling_dw->SetBinContent(i,j,grexp_coupling_dw->Interpolate(hexp_coupling_dw->GetXaxis()->GetBinCenter(i),hexp_coupling_dw->GetYaxis()->GetBinCenter(j)));
      hobs_coupling_up->SetBinContent(i,j,grobs_coupling_up->Interpolate(hobs_coupling_up->GetXaxis()->GetBinCenter(i),hobs_coupling_up->GetYaxis()->GetBinCenter(j)));
      hobs_coupling_dw->SetBinContent(i,j,grobs_coupling_dw->Interpolate(hobs_coupling_dw->GetXaxis()->GetBinCenter(i),hobs_coupling_dw->GetYaxis()->GetBinCenter(j)));
      if(not makeCOLZ){
	hexp_coupling_up2->SetBinContent(i,j,grexp_coupling_up2->Interpolate(hexp_coupling_up2->GetXaxis()->GetBinCenter(i),hexp_coupling_up2->GetYaxis()->GetBinCenter(j)));
	hexp_coupling_dw2->SetBinContent(i,j,grexp_coupling_dw2->Interpolate(hexp_coupling_dw2->GetXaxis()->GetBinCenter(i),hexp_coupling_dw2->GetYaxis()->GetBinCenter(j)));
      }
    }
  }

  // Smooth more
  hexp_coupling->Smooth();
  hobs_coupling->Smooth();
  hexp_coupling_up->Smooth();
  hobs_coupling_up->Smooth();
  hexp_coupling_dw->Smooth();
  hobs_coupling_dw->Smooth();
  if(not makeCOLZ){
    hexp_coupling_up2->Smooth();
    hexp_coupling_dw2->Smooth();
  }

  ///////-------------
  TFile* outputTemp = new TFile((outputDIR+"/outputTemp_vector.root").c_str(),"RECREATE");
  outputTemp->cd();
    
  RooRealVar* xvar = NULL;
  RooRealVar* yvar = NULL;

  if(useSplineND){    

    if(not useDMMass){
      xvar = new RooRealVar("x","x",(maxX_spline-minX_spline)/2,minX_spline,maxX_spline);
      yvar = new RooRealVar ("y","y",(maxCoupling_spline-minCoupling_spline)/2,minCoupling_spline,maxCoupling_spline);
    }
    else{
      xvar = new RooRealVar("x","x",(maxY_spline-minY_spline)/2,minY_spline,maxY_spline);
      yvar = new RooRealVar ("y","y",(maxCoupling_spline-minCoupling_spline)/2,minCoupling_spline,maxCoupling_spline);
    }
    
    RooSplineND* spline_exp = makeSplineND(hexp_coupling,*xvar,*yvar,"exp");
    RooSplineND* spline_exp_up = makeSplineND(hexp_coupling_up,*xvar,*yvar,"exp_up");
    RooSplineND* spline_exp_dw = makeSplineND(hexp_coupling_dw,*xvar,*yvar,"exp_dw");
    RooSplineND* spline_obs = makeSplineND(hobs_coupling,*xvar,*yvar,"obs");
    RooSplineND* spline_obs_up = makeSplineND(hobs_coupling_up,*xvar,*yvar,"obs_up");
    RooSplineND* spline_obs_dw = makeSplineND(hobs_coupling_dw,*xvar,*yvar,"obs_dw");
    
    RooSplineND* spline_exp_up2 = NULL;
    RooSplineND* spline_exp_dw2 = NULL;
    
    if(not makeCOLZ){
      spline_exp_up2 = makeSplineND(hexp_coupling_up2,*xvar,*yvar,"exp_up2");
      spline_exp_dw2 = makeSplineND(hexp_coupling_dw2,*xvar,*yvar,"exp_dw2");
    }
    
    TH2D* hexp_coupling_ext = NULL;
    TH2D* hexp_coupling_ext_up = NULL;
    TH2D* hexp_coupling_ext_dw = NULL;
    TH2D* hobs_coupling_ext = NULL;
    TH2D* hobs_coupling_ext_up = NULL;
    TH2D* hobs_coupling_ext_dw = NULL;
    TH2D* hexp_coupling_ext_up2 = NULL;
    TH2D* hexp_coupling_ext_dw2 = NULL;

    if(not useDMMass){
      hexp_coupling_ext = new TH2D("hexp_coupling_ext", "",nbinsX,minX,maxX,nCoupling,minCoupling,maxCoupling);
      hexp_coupling_ext_up = new TH2D("hexp_coupling_ext_up", "",nbinsX,minX,maxX,nCoupling,minCoupling,maxCoupling);
      hexp_coupling_ext_dw = new TH2D("hexp_coupling_ext_dw", "",nbinsX,minX,maxX,nCoupling,minCoupling,maxCoupling);
      hobs_coupling_ext = new TH2D("hobs_coupling_ext", "",nbinsX,minX,maxX,nCoupling,minCoupling,maxCoupling);
      hobs_coupling_ext_up = new TH2D("hobs_coupling_ext_up", "",nbinsX,minX,maxX,nCoupling,minCoupling,maxCoupling);
      hobs_coupling_ext_dw = new TH2D("hobs_coupling_ext_dw", "",nbinsX,minX,maxX,nCoupling,minCoupling,maxCoupling);
      if(not makeCOLZ){
	hexp_coupling_ext_up2 = new TH2D("hexp_coupling_ext_up2", "",nbinsX,minX,maxX,nCoupling,minCoupling,maxCoupling);
	hexp_coupling_ext_dw2 = new TH2D("hexp_coupling_ext_dw2", "",nbinsX,minX,maxX,nCoupling,minCoupling,maxCoupling);
      }
    }
    else{
      hexp_coupling_ext = new TH2D("hexp_coupling_ext", "",nbinsY,minY,maxY,nCoupling,minCoupling,maxCoupling);
      hexp_coupling_ext_up = new TH2D("hexp_coupling_ext_up", "",nbinsY,minY,maxY,nCoupling,minCoupling,maxCoupling);
      hexp_coupling_ext_dw = new TH2D("hexp_coupling_ext_dw", "",nbinsY,minY,maxY,nCoupling,minCoupling,maxCoupling);
      hobs_coupling_ext = new TH2D("hobs_coupling_ext", "",nbinsY,minY,maxY,nCoupling,minCoupling,maxCoupling);
      hobs_coupling_ext_up = new TH2D("hobs_coupling_ext_up", "",nbinsY,minY,maxY,nCoupling,minCoupling,maxCoupling);
      hobs_coupling_ext_dw = new TH2D("hobs_coupling_ext_dw", "",nbinsY,minY,maxY,nCoupling,minCoupling,maxCoupling);
      if(not makeCOLZ){
	hexp_coupling_ext_up2 = new TH2D("hexp_coupling_ext_up2", "",nbinsY,minY,maxY,nCoupling,minCoupling,maxCoupling);
	hexp_coupling_ext_dw2 = new TH2D("hexp_coupling_ext_dw2", "",nbinsY,minY,maxY,nCoupling,minCoupling,maxCoupling);
      }
    }

    cout<<"Evaluating the splines "<<endl;
    for(int i = 1; i <= hexp_coupling_ext->GetNbinsX(); i++){
      for(int j = 1; j <= hexp_coupling_ext->GetNbinsY(); j++){
	xvar->setVal(hexp_coupling_ext->GetXaxis()->GetBinCenter(i));
	yvar->setVal(hexp_coupling_ext->GetYaxis()->GetBinCenter(j));
	hexp_coupling_ext->SetBinContent(i,j,spline_exp->getVal());
	hexp_coupling_ext_up->SetBinContent(i,j,spline_exp_up->getVal());
	hexp_coupling_ext_dw->SetBinContent(i,j,spline_exp_dw->getVal());
	hobs_coupling_ext->SetBinContent(i,j,spline_obs->getVal());
	hobs_coupling_ext_up->SetBinContent(i,j,spline_obs_up->getVal());
	hobs_coupling_ext_dw->SetBinContent(i,j,spline_obs_dw->getVal());
	if(not makeCOLZ){
	  hexp_coupling_ext_up2->SetBinContent(i,j,spline_exp_up2->getVal());
	  hexp_coupling_ext_dw2->SetBinContent(i,j,spline_exp_dw2->getVal());
	}
      }
    }
    
    hexp_coupling = hexp_coupling_ext;
    hexp_coupling_up = hexp_coupling_ext_up;
    hexp_coupling_dw = hexp_coupling_ext_dw;
    hobs_coupling    = hobs_coupling_ext;
    hobs_coupling_up = hobs_coupling_ext_up;
    hobs_coupling_dw = hobs_coupling_ext_dw;
    if(not makeCOLZ){
      hexp_coupling_up2 = hexp_coupling_ext_up2;
      hexp_coupling_dw2 = hexp_coupling_ext_dw2;
    }
  }
  
  // extend the relic density line                                                                                                                                                                   
  TGraph* relic_graph_g1_ext = new TGraph();
  TGraph* relic_graph_g2_ext = new TGraph();

  int ipoint = 0;
  for(int i = 1; i <= hexp_coupling->GetNbinsX(); i++){
    relic_graph_g1_ext->SetPoint(ipoint,hexp_coupling->GetXaxis()->GetBinCenter(i),relic_graph_g1->Eval(hexp_coupling->GetXaxis()->GetBinCenter(i)));
    relic_graph_g2_ext->SetPoint(ipoint,hexp_coupling->GetXaxis()->GetBinCenter(i),relic_graph_g2->Eval(hexp_coupling->GetXaxis()->GetBinCenter(i)));
    ipoint++;
  }

  fixBinContent(hexp_coupling,minZ,maxZ);
  fixBinContent(hexp_coupling_up,minZ,maxZ);
  fixBinContent(hexp_coupling_dw,minZ,maxZ);
  fixBinContent(hobs_coupling,minZ,maxZ);
  fixBinContent(hobs_coupling_up,minZ,maxZ);
  fixBinContent(hobs_coupling_dw,minZ,maxZ);
  if(not makeCOLZ){
    fixBinContent(hexp_coupling_up2,minZ,maxZ);
    fixBinContent(hexp_coupling_dw2,minZ,maxZ);
  }


  ////////////////                        
  cout<<"Making contours "<<endl;
  TH2* hexp2 = (TH2*) hexp_coupling->Clone("hexp2");
  TH2* hobs2 = (TH2*) hobs_coupling->Clone("hobs2");
  TH2* hexp2_up = (TH2*) hexp_coupling_up->Clone("hexp2_up");
  TH2* hobs2_up = (TH2*) hobs_coupling_up->Clone("hobs2_up");
  TH2* hexp2_dw = (TH2*) hexp_coupling_dw->Clone("hexp2_dw");
  TH2* hobs2_dw = (TH2*) hobs_coupling_dw->Clone("hobs2_dw");
  TH2* hexp2_up2 = NULL;
  TH2* hexp2_dw2 = NULL;
  if(not makeCOLZ){
    hexp2_up2 = (TH2*) hexp_coupling_up2->Clone("hexp2_up2");
    hexp2_dw2 = (TH2*) hexp_coupling_dw2->Clone("hexp2_dw2");    
  }

  //////////                                                                                                                                                                                           
  double contours[1]; contours[0]=1;
  hexp2->SetContour(1,contours);
  hobs2->SetContour(1,contours);
  hexp2_up->SetContour(1,contours);
  hobs2_up->SetContour(1,contours);
  hexp2_dw->SetContour(1,contours);
  hobs2_dw->SetContour(1,contours);
  if(not makeCOLZ){
    hexp2_up2->SetContour(1,contours);
    hexp2_dw2->SetContour(1,contours);
  }

  // All the plotting and cosmetics                                                                                                                                                                   
  TCanvas* canvas = new TCanvas("canvas", "canvas",650,600);
  if(makeCOLZ){
    canvas->SetRightMargin(0.15);
    canvas->SetLeftMargin(0.13);
  }
  else{
    canvas->SetRightMargin(0.10);
    canvas->SetLeftMargin(0.13);
  }

  canvas->SetLogz();

  TH1* frame = NULL;
  if(not useDMMass)
    frame = canvas->DrawFrame(minX,minCoupling,maxX,maxCoupling,"");
  else
    frame = canvas->DrawFrame(minY,minCoupling,maxY,maxCoupling,"");

  frame->GetYaxis()->CenterTitle();
  if(not useDMMass)
    frame->GetXaxis()->SetTitle("m_{med} [GeV]");
  else
    frame->GetXaxis()->SetTitle("m_{DM} [GeV]");

  frame->GetYaxis()->SetTitle("g_{q}");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.20);
  frame->Draw();

  hexp_coupling->GetZaxis()->SetRangeUser(minZ,maxZ);
  hobs_coupling->GetZaxis()->SetRangeUser(minZ,maxZ);
  hexp_coupling_up->GetZaxis()->SetRangeUser(minZ,maxZ);
  hobs_coupling_up->GetZaxis()->SetRangeUser(minZ,maxZ);
  hexp_coupling_dw->GetZaxis()->SetRangeUser(minZ,maxZ);
  hobs_coupling_dw->GetZaxis()->SetRangeUser(minZ,maxZ);
  if(not makeCOLZ){
    hexp_coupling_up2->GetZaxis()->SetRangeUser(minZ,maxZ);
    hexp_coupling_dw2->GetZaxis()->SetRangeUser(minZ,maxZ);
  }

  hexp2->GetZaxis()->SetLabelSize(0);
  hexp2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp = produceContour(reductionForContour);
  TGraph* contour_exp_back = produceContour(reductionForContour,false);

  hexp2_up->GetZaxis()->SetLabelSize(0);
  hexp2_up->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp_up = NULL;
  if(makeCOLZ)
    contour_exp_up = produceContour(reductionForContour);
  else
    contour_exp_up = produceContour(reductionForContour,false);

  hexp2_dw->GetZaxis()->SetLabelSize(0);
  hexp2_dw->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp_dw = NULL;
  if(makeCOLZ)
    contour_exp_dw = produceContour(reductionForContour);
  else
    contour_exp_dw = produceContour(reductionForContour,false);

  TGraph* contour_exp_up2 = NULL;
  TGraph* contour_exp_dw2 = NULL;

  if(not makeCOLZ){
    hexp2_up2->GetZaxis()->SetLabelSize(0);
    hexp2_up2->Draw("contz list same");
    canvas->Update();
    contour_exp_up2 = produceContour(reductionForContour,false);

    hexp2_dw2->GetZaxis()->SetLabelSize(0);
    hexp2_dw2->Draw("contz list same");
    canvas->Update();
    contour_exp_dw2 = produceContour(reductionForContour,false);
  }

  hobs2->GetZaxis()->SetLabelSize(0);
  hobs2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_obs = produceContour(reductionForContour);

  hobs2_up->GetZaxis()->SetLabelSize(0);
  hobs2_up->Draw("contz list same");
  canvas->Update();
  TGraph* contour_obs_up = produceContour(reductionForContour);

  hobs2_dw->GetZaxis()->SetLabelSize(0);
  hobs2_dw->Draw("contz list same");
  canvas->Update();
  TGraph* contour_obs_dw = produceContour(reductionForContour);
  
  frame->Draw();

  cout<<"Plotting things "<<endl;
  if(makeCOLZ)
    hobs_coupling->Draw("COLZ SAME");
  
  if(addRelicDensity and makeCOLZ){
    relic_graph_g1_ext->SetLineColor(kGreen+3);
    relic_graph_g2_ext->SetLineColor(kGreen+3);
    relic_graph_g1_ext->SetLineWidth(-802);
    relic_graph_g2_ext->SetLineWidth(802);
    relic_graph_g1_ext->SetFillStyle(3005);
    relic_graph_g2_ext->SetFillStyle(3005);
    relic_graph_g1_ext->SetFillColor(kGreen+3);
    relic_graph_g2_ext->SetFillColor(kGreen+3);
    relic_graph_g1_ext->Draw("L SAME");
    relic_graph_g2_ext->Draw("L SAME");                                                                                                                                                           
  }

  TGraphAsymmErrors* graph_1s = NULL;
  TGraphAsymmErrors* graph_2s = NULL;
  TGraphAsymmErrors* graph_1s_y = NULL;
  TGraphAsymmErrors* graph_2s_y = NULL;

  if(makeCOLZ){
    contour_exp_up->SetLineColor(kBlack);
    contour_exp_up->SetLineWidth(2);
    contour_exp_up->SetLineStyle(7);
    contour_exp_up->Draw("Lsame");
    contour_exp_dw->SetLineColor(kBlack);
    contour_exp_dw->SetLineWidth(2);
    contour_exp_dw->SetLineStyle(7);
    contour_exp_dw->Draw("Lsame");
  }
  // make a Brazilian like plot
  else{
    
    graph_1s = new TGraphAsymmErrors();
    graph_2s = new TGraphAsymmErrors();
    graph_1s_y = new TGraphAsymmErrors();
    graph_2s_y = new TGraphAsymmErrors();

    // make band where for each x there is 1 x
    TGraph* band_1_up = new TGraph();
    TGraph* band_1_dw = new TGraph();
    TGraph* band_2_up = new TGraph();
    TGraph* band_2_dw = new TGraph();

    makeXGraph(band_1_up,contour_exp_up,contour_exp_back);
    makeXGraph(band_1_dw,contour_exp_dw,contour_exp_back);
    makeXGraph(band_2_up,contour_exp_up2,contour_exp_back);
    makeXGraph(band_2_dw,contour_exp_dw2,contour_exp_back);

    band_1_up->Write("band_1_up");
    band_1_dw->Write("band_1_dw");
    band_2_up->Write("band_2_up");
    band_2_dw->Write("band_2_dw");
    
    makeUncertaintyBand(graph_1s,contour_exp_back,band_1_up,band_1_dw,useDMMass);
    makeUncertaintyBand(graph_2s,contour_exp_back,band_2_up,band_2_dw,useDMMass);

    TGraph* band_1_y = new TGraph();
    TGraph* band_1_up_y = new TGraph();
    TGraph* band_1_dw_y = new TGraph();
    TGraph* band_2_up_y = new TGraph();
    TGraph* band_2_dw_y = new TGraph();

    makeYGraph(band_1_y,contour_exp_back);
    makeYGraph(band_1_up_y,contour_exp_up);
    makeYGraph(band_1_dw_y,contour_exp_dw);
    makeYGraph(band_2_up_y,contour_exp_up2);
    makeYGraph(band_2_dw_y,contour_exp_dw2);

    band_1_up_y->Write("band_1_up_y");
    band_1_dw_y->Write("band_1_dw_y");
    band_2_up_y->Write("band_2_up_y");
    band_2_dw_y->Write("band_2_dw_y");

    makeUncertaintyBand(graph_1s_y,contour_exp_back,band_1_y,band_1_up_y,band_1_dw_y,useDMMass);
    makeUncertaintyBand(graph_2s_y,contour_exp_back,band_1_y,band_2_up_y,band_2_dw_y,useDMMass);
    
    graph_1s->SetFillStyle(1001);
    graph_1s->SetLineColor(kBlack);
    graph_1s->SetFillColor(kGreen+1);
    graph_2s->SetFillStyle(1001);
    graph_2s->SetFillColor(kOrange);
    graph_2s->SetLineColor(kBlack);

    graph_1s_y->SetFillStyle(1001);
    graph_1s_y->SetLineColor(kBlack);
    graph_1s_y->SetFillColor(kGreen+1);
    graph_2s_y->SetFillStyle(1001);
    graph_2s_y->SetFillColor(kOrange);
    graph_2s_y->SetLineColor(kBlack);

    graph_2s->Draw("2same");
    graph_1s->Draw("2same");
    graph_2s_y->Draw("2same");
    graph_1s_y->Draw("2same");

    graph_2s->Write("graph_2s");
    graph_1s->Write("graph_1s");
    graph_2s_y->Write("graph_2s_y");
    graph_1s_y->Write("graph_1s_y");

  }

  if(addRelicDensity and not makeCOLZ){
    relic_graph_g1_ext->SetLineColor(kBlue);
    relic_graph_g2_ext->SetLineColor(kBlue);
    relic_graph_g1_ext->SetLineWidth(-802);
    relic_graph_g2_ext->SetLineWidth(802);
    relic_graph_g1_ext->SetFillStyle(3005);
    relic_graph_g2_ext->SetFillStyle(3005);
    relic_graph_g1_ext->SetFillColor(kBlue);
    relic_graph_g2_ext->SetFillColor(kBlue);
    relic_graph_g1_ext->Draw("L SAME");
    relic_graph_g2_ext->Draw("L SAME");                                                                                                                                                           
  }


  if(makeCOLZ){
    contour_obs_up->SetLineColor(kRed);
    contour_obs_up->SetLineWidth(2);
    contour_obs_up->SetLineStyle(7);
    contour_obs_up->Draw("Lsame");
    
    contour_obs_dw->SetLineColor(kRed);
    contour_obs_dw->SetLineWidth(2);
    contour_obs_dw->SetLineStyle(7);
    contour_obs_dw->Draw("Lsame");
    
    contour_exp->SetLineColor(kBlack);
    contour_exp->SetLineWidth(3);
    contour_exp->Draw("Lsame");
    
    contour_obs->SetLineColor(kRed);
    contour_obs->SetLineWidth(3);
    contour_obs->Draw("Lsame");
  }
  else{
    contour_obs_up->SetLineColor(kRed);
    contour_obs_up->SetLineWidth(2);
    contour_obs_up->Draw("Lsame");
    
    contour_obs_dw->SetLineColor(kRed);
    contour_obs_dw->SetLineWidth(2);
    contour_obs_dw->Draw("Lsame");
    
    contour_obs->SetLineColor(kBlack);
    contour_obs->SetLineWidth(3);
    contour_obs->Draw("Lsame");
    
    contour_exp->SetLineColor(kBlack);
    contour_exp->SetLineWidth(3);
    contour_exp->SetLineStyle(7);
    contour_exp->Draw("Lsame");
  }

  if(not addPreliminary)
    CMS_lumi(canvas,"35.9",true,true,false,0,-0.09);
  else
    CMS_lumi(canvas,"35.9",true,false,false,0,-0.09);

  
  TLegend* tex = new TLegend(0.15,0.83,0.65,0.89);
  if(makeCOLZ){
    tex->SetFillColor(kWhite);
    tex->SetFillStyle(1001);
  }
  else{
    tex->SetFillColor(0);
    tex->SetFillStyle(0);
  }
  tex->SetBorderSize(0);
  tex->SetTextFont(42);
  tex->SetTextSize(0.027);
  tex->SetHeader("#bf{Vector med, Dirac DM, g_{DM} = 1, m_{med} = 3 #times m_{DM}}");
  tex->Draw("same");

  TLegend *leg = NULL;
  if(makeCOLZ)
    leg = new TLegend(0.15,0.58,0.52,0.81);
  else
    leg = new TLegend(0.15,0.55,0.52,0.81);

  if(makeCOLZ){
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
  }
  else{
    leg->SetFillColor(0);
    leg->SetFillStyle(0);    
  }
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0243902);
  leg->AddEntry(contour_exp,"Median expected 95% CL","L");
  if(makeCOLZ){
    leg->AddEntry(contour_exp_up,"68% expected","L");
  }
  else{
    leg->AddEntry(graph_1s,"68% expected","F");
    leg->AddEntry(graph_2s,"95% expected","F");
  }
  leg->AddEntry(contour_obs,"Observed 95% CL","L");
  leg->AddEntry(contour_obs_up,"Observed #pm theory unc.","L");
  if(addRelicDensity)
    leg->AddEntry(relic_graph_g1_ext,"#Omega_{c}#timesh^{2} #geq 0.12","F");
  leg->Draw("same");



  TLatex *   tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetLineWidth(2);
  tex2->SetTextSize(0.042);
  tex2->SetTextAngle(90);
  if(makeCOLZ){
    tex2->DrawLatex(0.975,0.55,"Observed #sigma_{95% CL}/#sigma_{th}");
  }

  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR+"/scan_vector_mdm_vs_gq_"+string(energy)+"TeV.pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_vs_gq_"+string(energy)+"TeV.png").c_str());

  canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_vs_gq_"+string(energy)+"TeV_log.pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_vs_gq_"+string(energy)+"TeV_log.png").c_str());

  contour_exp_back->Write("contour_exp_back");
  contour_exp_up->Write("contour_exp_up");
  contour_exp_dw->Write("contour_exp_dw");
  contour_exp_up2->Write("contour_exp_up2");
  contour_exp_dw2->Write("contour_exp_dw2");

  outputTemp->Close();
  //  system(("rm "+outputDIR+"/outputTemp_vector.root").c_str());


}

