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
static float nbinsX = 1000;
static float nbinsY = 600;
static float minX = 200;
static float minY = 1;
static float maxX = 2000;
static float maxY = 1000;
static float nCoupling   = 100;
static float minCoupling = 0.02;
static float maxCoupling = 1;
static float minZ = 0.1;
static float maxZ = 100;
static int   reductionForContour = 20;
static bool  addPreliminary = false;
static bool  whiteOut       = true;
static bool  skipPoints     = true;
static bool  addRelicDensity = true;

void fillLimitGraphs(TTree* tree,
		     TGraph2D* grexp, TGraph2D* grexp_up, TGraph2D* grexp_dw,
		     TGraph2D* grobs, TGraph2D* grobs_up, TGraph2D* grobs_dw,
		     const float & coupling){ // in case one needs to produce 2D plots first
  
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
      if(coupling == 0.25){
	if(medmass == 1800 and (dmmass == 200 or dmmass == 250 or dmmass == 350 or dmmass == 400 or dmmass == 800)) continue;
	if(medmass == 1725 and (dmmass == 200 or dmmass > 700)) continue;
	if(medmass == 1525 and dmmass > 600) continue;
	if(medmass == 1125 and dmmass > 600) continue;
	if(medmass == 600  and dmmass == 350) continue;
	if(medmass == 525  and dmmass == 275) continue;
	if(medmass <= 400  and dmmass >= 400) continue;
	if(quantile == -1  and medmass == 1925 and dmmass == 200) continue;
	if(quantile == -1  and medmass == 1925 and dmmass == 250) continue;
	if(quantile == -1  and medmass == 1800 and dmmass == 300) continue;
	if(quantile == -1  and medmass == 2000 and dmmass == 200) continue;
	
      }
      else if(coupling == 1.){
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
      else if(coupling == 0.1){
	if(medmass == 1400 and dmmass == 400) continue;
	if(medmass == 1400 and dmmass == 500) continue;
      }      
    }


    if (quantile == 0.5) { // expected limit
      grexp->SetPoint(expcounter, double(medmass), double(dmmass), limit);
      expcounter++;
    }

    if (quantile < 0.17 && quantile > 0.14 ) {
      grexp_dw->SetPoint(exp_down_counter, double(medmass), double(dmmass), limit);
      exp_down_counter++;
    }

    if (quantile < 0.85 && quantile > 0.83 ) {
      grexp_up->SetPoint(exp_up_counter, double(medmass), double(dmmass), limit);
      exp_up_counter++;
    }    
    
    if (quantile == -1) { // observed
      grobs->SetPoint(obscounter, double(medmass), double(dmmass), limit);
      grobs_up->SetPoint(obscounter, double(medmass), double(dmmass), limit*0.8);
      grobs_dw->SetPoint(obscounter, double(medmass), double(dmmass), limit*1.2);
      obscounter++;      
    }
  } 
  tree->ResetBranchAddresses();
}


void fillLimitGraphs(TTree* tree,
		     TGraph* grexp, TGraph* grexp_up, TGraph* grexp_dw,
		     TGraph* grobs, TGraph* grobs_up, TGraph* grobs_dw,
		     const float & medOverDM = 3, const bool & useDMMass = false){

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
  int exp_down_counter = 0;
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
      if(useDMMass)
	grexp->SetPoint(expcounter, double(dmmass), limit);
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

    if (quantile < 0.85 && quantile > 0.83 ) {
      if(useDMMass)
	grexp_up->SetPoint(exp_up_counter, double(dmmass), limit);
      else
	grexp_up->SetPoint(exp_up_counter, double(medmass),limit);

      exp_up_counter++;
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


///////////////
void plotVectorCoupling(string outputDIR, bool isFull2DGrid = false, bool useDMMass = false, float medOverDM = 3, string energy = "13") {

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
      if(*mmed != medOverDM*(*mdm)) continue;
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
  vector<pair<double,TGraph*> > grexp;
  vector<pair<double,TGraph*> > grobs;
  vector<pair<double,TGraph*> > grobs_up;
  vector<pair<double,TGraph*> > grobs_dw;

  if(isFull2DGrid){

    vector<pair<double,TGraph2D*> > grexp_2d_up; 
    vector<pair<double,TGraph2D*> > grexp_2d_dw; 
    vector<pair<double,TGraph2D*> > grexp_2d; 
    vector<pair<double,TGraph2D*> > grobs_2d; 
    vector<pair<double,TGraph2D*> > grobs_2d_dw; 
    vector<pair<double,TGraph2D*> > grobs_2d_up; 
    
    grexp_2d.push_back(pair<double,TGraph2D*>(1.00,new TGraph2D())); 
    grexp_2d.push_back(pair<double,TGraph2D*>(0.75,new TGraph2D())); 
    grexp_2d.push_back(pair<double,TGraph2D*>(0.50,new TGraph2D())); 
    grexp_2d.push_back(pair<double,TGraph2D*>(0.30,new TGraph2D())); 
    grexp_2d.push_back(pair<double,TGraph2D*>(0.25,new TGraph2D())); 
    grexp_2d.push_back(pair<double,TGraph2D*>(0.20,new TGraph2D())); 
    grexp_2d.push_back(pair<double,TGraph2D*>(0.10,new TGraph2D())); 
    grexp_2d.push_back(pair<double,TGraph2D*>(0.05,new TGraph2D())); 
    grexp_2d.push_back(pair<double,TGraph2D*>(0.01,new TGraph2D())); 

    grexp_2d_up.push_back(pair<double,TGraph2D*>(1.00,new TGraph2D())); 
    grexp_2d_up.push_back(pair<double,TGraph2D*>(0.75,new TGraph2D())); 
    grexp_2d_up.push_back(pair<double,TGraph2D*>(0.50,new TGraph2D())); 
    grexp_2d_up.push_back(pair<double,TGraph2D*>(0.30,new TGraph2D())); 
    grexp_2d_up.push_back(pair<double,TGraph2D*>(0.25,new TGraph2D())); 
    grexp_2d_up.push_back(pair<double,TGraph2D*>(0.20,new TGraph2D())); 
    grexp_2d_up.push_back(pair<double,TGraph2D*>(0.10,new TGraph2D())); 
    grexp_2d_up.push_back(pair<double,TGraph2D*>(0.05,new TGraph2D())); 
    grexp_2d_up.push_back(pair<double,TGraph2D*>(0.01,new TGraph2D())); 

    grexp_2d_dw.push_back(pair<double,TGraph2D*>(1.00,new TGraph2D())); 
    grexp_2d_dw.push_back(pair<double,TGraph2D*>(0.75,new TGraph2D())); 
    grexp_2d_dw.push_back(pair<double,TGraph2D*>(0.50,new TGraph2D())); 
    grexp_2d_dw.push_back(pair<double,TGraph2D*>(0.30,new TGraph2D())); 
    grexp_2d_dw.push_back(pair<double,TGraph2D*>(0.25,new TGraph2D())); 
    grexp_2d_dw.push_back(pair<double,TGraph2D*>(0.20,new TGraph2D())); 
    grexp_2d_dw.push_back(pair<double,TGraph2D*>(0.10,new TGraph2D())); 
    grexp_2d_dw.push_back(pair<double,TGraph2D*>(0.05,new TGraph2D())); 
    grexp_2d_dw.push_back(pair<double,TGraph2D*>(0.01,new TGraph2D())); 

    grobs_2d.push_back(pair<double,TGraph2D*>(1.00,new TGraph2D())); 
    grobs_2d.push_back(pair<double,TGraph2D*>(0.75,new TGraph2D())); 
    grobs_2d.push_back(pair<double,TGraph2D*>(0.50,new TGraph2D())); 
    grobs_2d.push_back(pair<double,TGraph2D*>(0.30,new TGraph2D())); 
    grobs_2d.push_back(pair<double,TGraph2D*>(0.25,new TGraph2D())); 
    grobs_2d.push_back(pair<double,TGraph2D*>(0.20,new TGraph2D())); 
    grobs_2d.push_back(pair<double,TGraph2D*>(0.10,new TGraph2D())); 
    grobs_2d.push_back(pair<double,TGraph2D*>(0.05,new TGraph2D())); 
    grobs_2d.push_back(pair<double,TGraph2D*>(0.01,new TGraph2D())); 

    grobs_2d_up.push_back(pair<double,TGraph2D*>(1.00,new TGraph2D())); 
    grobs_2d_up.push_back(pair<double,TGraph2D*>(0.75,new TGraph2D())); 
    grobs_2d_up.push_back(pair<double,TGraph2D*>(0.50,new TGraph2D())); 
    grobs_2d_up.push_back(pair<double,TGraph2D*>(0.30,new TGraph2D())); 
    grobs_2d_up.push_back(pair<double,TGraph2D*>(0.25,new TGraph2D())); 
    grobs_2d_up.push_back(pair<double,TGraph2D*>(0.20,new TGraph2D())); 
    grobs_2d_up.push_back(pair<double,TGraph2D*>(0.10,new TGraph2D())); 
    grobs_2d_up.push_back(pair<double,TGraph2D*>(0.05,new TGraph2D())); 
    grobs_2d_up.push_back(pair<double,TGraph2D*>(0.01,new TGraph2D()));
 
    grobs_2d_dw.push_back(pair<double,TGraph2D*>(1.00,new TGraph2D())); 
    grobs_2d_dw.push_back(pair<double,TGraph2D*>(0.75,new TGraph2D())); 
    grobs_2d_dw.push_back(pair<double,TGraph2D*>(0.50,new TGraph2D())); 
    grobs_2d_dw.push_back(pair<double,TGraph2D*>(0.30,new TGraph2D())); 
    grobs_2d_dw.push_back(pair<double,TGraph2D*>(0.25,new TGraph2D())); 
    grobs_2d_dw.push_back(pair<double,TGraph2D*>(0.20,new TGraph2D())); 
    grobs_2d_dw.push_back(pair<double,TGraph2D*>(0.10,new TGraph2D())); 
    grobs_2d_dw.push_back(pair<double,TGraph2D*>(0.05,new TGraph2D())); 
    grobs_2d_dw.push_back(pair<double,TGraph2D*>(0.01,new TGraph2D())); 
  
    // fill 2D graphs with limits values
    for(int itree = 0; itree < treeList.size(); itree++)
      fillLimitGraphs(treeList.at(itree),grexp_2d.at(itree).second,grexp_2d_up.at(itree).second,grexp_2d_dw.at(itree).second,
		      grobs_2d.at(itree).second,grobs_2d_up.at(itree).second,grobs_2d_dw.at(itree).second,
		      grexp_2d.at(itree).first);
    
    TH2D* hexp  = new TH2D("hexp", "",      nbinsX, minX, maxX, nbinsY, minY, maxY);
    TH2D* hobs  = new TH2D("hobs", "",      nbinsX, minX, maxX, nbinsY, minY, maxY);
  
    for(auto gr : grexp_2d)
      grexp.push_back(pair<double,TGraph*>(gr.first,new TGraph()));
    for(auto gr : grexp_2d_up)
      grexp_up.push_back(pair<double,TGraph*>(gr.first,new TGraph()));
    for(auto gr : grexp_2d_dw)
      grexp_dw.push_back(pair<double,TGraph*>(gr.first,new TGraph()));
    for(auto gr : grobs_2d)
      grobs.push_back(pair<double,TGraph*>(gr.first,new TGraph()));
    for(auto gr : grobs_2d_up)
      grobs_up.push_back(pair<double,TGraph*>(gr.first,new TGraph()));
    for(auto gr : grobs_2d_dw)
      grobs_dw.push_back(pair<double,TGraph*>(gr.first,new TGraph()));
  
    for(int igraph = 0; igraph < grexp.size(); igraph++){
      int ipoint_exp = 0;
      int ipoint_obs = 0;

      for(int i = 1; i <= nbinsX; i++){	

	if(not useDMMass){
	  
	  /// fix bad limits for coupling 1.0
	  if(grexp.at(igraph).first == 1.0){ 	  
	    if((hexp->GetXaxis()->GetBinCenter(i) < 700 or hexp->GetXaxis()->GetBinCenter(i) > 780) and (hexp->GetXaxis()->GetBinCenter(i) < 450 or hexp->GetXaxis()->GetBinCenter(i) > 500)){	    

	      grexp.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						grexp_2d.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));	    
	      grexp_up.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						   grexp_2d_up.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));	    
	      grexp_dw.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						   grexp_2d_dw.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));	    
	      ipoint_exp++;
	    }
	    
	    if((hobs->GetXaxis()->GetBinCenter(i) < 720 or hobs->GetXaxis()->GetBinCenter(i) > 800) and (hobs->GetXaxis()->GetBinCenter(i) < 450 or hobs->GetXaxis()->GetBinCenter(i) > 500)){

	      grobs.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),
						grobs_2d.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	      grobs_up.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),
						   grobs_2d_up.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	      grobs_dw.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),
						   grobs_2d_dw.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	      ipoint_obs++;
	    }
	  }
	  
	  // skip some spikes for gq  = 0.25
	  else if(grexp.at(igraph).first == 0.25){
	    if((hexp->GetXaxis()->GetBinCenter(i) < 600 or hexp->GetXaxis()->GetBinCenter(i) > 800) and (hexp->GetXaxis()->GetBinCenter(i) < 60 or hexp->GetXaxis()->GetBinCenter(i) > 100) and 
	       (hexp->GetXaxis()->GetBinCenter(i) < 30 or hexp->GetXaxis()->GetBinCenter(i) > 50) and (hexp->GetXaxis()->GetBinCenter(i) < 120 or hexp->GetXaxis()->GetBinCenter(i) > 300)){
	      
	      grexp.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						grexp_2d.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	      grexp_up.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						grexp_2d_up.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	      grexp_dw.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						grexp_2d_dw.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	      ipoint_exp++;
	    }
	    
	    
	    if((hobs->GetXaxis()->GetBinCenter(i) < 600 or hobs->GetXaxis()->GetBinCenter(i) > 800) and (hobs->GetXaxis()->GetBinCenter(i) < 60 or hobs->GetXaxis()->GetBinCenter(i) > 100) and 
	       (hobs->GetXaxis()->GetBinCenter(i) < 30 or hobs->GetXaxis()->GetBinCenter(i) > 50) and (hobs->GetXaxis()->GetBinCenter(i) < 120 or hobs->GetXaxis()->GetBinCenter(i) > 300)){

	      grobs.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),grobs_2d.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	      grobs_up.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),
						   grobs_2d_up.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	      grobs_dw.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),
						   grobs_2d_dw.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	      ipoint_obs++;
	    }
	  }
	  else if(grexp.at(igraph).first == 0.1){
	    
	    if(hexp->GetXaxis()->GetBinCenter(i) < 700 or hexp->GetXaxis()->GetBinCenter(i) > 780){

	      grexp.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						grexp_2d.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	      grexp_up.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						   grexp_2d_up.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	      grexp_dw.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						   grexp_2d_dw.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	      ipoint_exp++;
	    }
	    
	    if(hobs->GetXaxis()->GetBinCenter(i) < 700 or hobs->GetXaxis()->GetBinCenter(i) > 780){

	      grobs.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),
						grobs_2d.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	      grobs_up.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),
						   grobs_2d_up.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	      grobs_dw.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),
						   grobs_2d_dw.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	      ipoint_obs++;
	    }
	  }
	  else{
	    grexp.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
					      grexp_2d.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	    grexp_up.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						 grexp_2d_up.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	    grexp_dw.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						 grexp_2d_dw.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
  
	    ipoint_exp++;

	    grobs.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),
					      grobs_2d.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	    grobs_up.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),
						 grobs_2d_up.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	    grobs_dw.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),
						 grobs_2d_dw.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
       
	    ipoint_obs++;
	  }
	}
	else{

	  grexp.at(igraph).second->SetPoint(i,hexp->GetXaxis()->GetBinCenter(i)/3,
					    grexp_2d.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	  grexp_up.at(igraph).second->SetPoint(i,hexp->GetXaxis()->GetBinCenter(i)/3,
					       grexp_2d_up.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	  grexp_dw.at(igraph).second->SetPoint(i,hexp->GetXaxis()->GetBinCenter(i)/3,
					       grexp_2d_dw.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));

	  grobs.at(igraph).second->SetPoint(i,hobs->GetXaxis()->GetBinCenter(i)/3,grobs_2d.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));

	  grobs_up.at(igraph).second->SetPoint(i,hobs->GetXaxis()->GetBinCenter(i)/3,
					       grobs_2d_up.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	  grobs_dw.at(igraph).second->SetPoint(i,hobs->GetXaxis()->GetBinCenter(i)/3,
					       grobs_2d_dw.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	}
      }
    }    
  }
  else{

    grexp.push_back(pair<double,TGraph*>(1.00,new TGraph()));
    grexp.push_back(pair<double,TGraph*>(0.75,new TGraph()));
    grexp.push_back(pair<double,TGraph*>(0.50,new TGraph()));
    grexp.push_back(pair<double,TGraph*>(0.30,new TGraph()));
    grexp.push_back(pair<double,TGraph*>(0.25,new TGraph()));
    grexp.push_back(pair<double,TGraph*>(0.20,new TGraph()));
    grexp.push_back(pair<double,TGraph*>(0.10,new TGraph()));
    grexp.push_back(pair<double,TGraph*>(0.05,new TGraph()));
    grexp.push_back(pair<double,TGraph*>(0.01,new TGraph()));

    grexp_up.push_back(pair<double,TGraph*>(1.00,new TGraph()));
    grexp_up.push_back(pair<double,TGraph*>(0.75,new TGraph()));
    grexp_up.push_back(pair<double,TGraph*>(0.50,new TGraph()));
    grexp_up.push_back(pair<double,TGraph*>(0.30,new TGraph()));
    grexp_up.push_back(pair<double,TGraph*>(0.25,new TGraph()));
    grexp_up.push_back(pair<double,TGraph*>(0.20,new TGraph()));
    grexp_up.push_back(pair<double,TGraph*>(0.10,new TGraph()));
    grexp_up.push_back(pair<double,TGraph*>(0.05,new TGraph()));
    grexp_up.push_back(pair<double,TGraph*>(0.01,new TGraph()));

    grexp_dw.push_back(pair<double,TGraph*>(1.00,new TGraph()));
    grexp_dw.push_back(pair<double,TGraph*>(0.75,new TGraph()));
    grexp_dw.push_back(pair<double,TGraph*>(0.50,new TGraph()));
    grexp_dw.push_back(pair<double,TGraph*>(0.30,new TGraph()));
    grexp_dw.push_back(pair<double,TGraph*>(0.25,new TGraph()));
    grexp_dw.push_back(pair<double,TGraph*>(0.20,new TGraph()));
    grexp_dw.push_back(pair<double,TGraph*>(0.10,new TGraph()));
    grexp_dw.push_back(pair<double,TGraph*>(0.05,new TGraph()));
    grexp_dw.push_back(pair<double,TGraph*>(0.01,new TGraph()));

    grobs.push_back(pair<double,TGraph*>(1.00,new TGraph()));
    grobs.push_back(pair<double,TGraph*>(0.75,new TGraph()));
    grobs.push_back(pair<double,TGraph*>(0.50,new TGraph()));
    grobs.push_back(pair<double,TGraph*>(0.30,new TGraph()));
    grobs.push_back(pair<double,TGraph*>(0.25,new TGraph()));
    grobs.push_back(pair<double,TGraph*>(0.20,new TGraph()));
    grobs.push_back(pair<double,TGraph*>(0.10,new TGraph()));
    grobs.push_back(pair<double,TGraph*>(0.05,new TGraph()));
    grobs.push_back(pair<double,TGraph*>(0.01,new TGraph()));

    grobs_up.push_back(pair<double,TGraph*>(1.00,new TGraph()));
    grobs_up.push_back(pair<double,TGraph*>(0.75,new TGraph()));
    grobs_up.push_back(pair<double,TGraph*>(0.50,new TGraph()));
    grobs_up.push_back(pair<double,TGraph*>(0.30,new TGraph()));
    grobs_up.push_back(pair<double,TGraph*>(0.25,new TGraph()));
    grobs_up.push_back(pair<double,TGraph*>(0.20,new TGraph()));
    grobs_up.push_back(pair<double,TGraph*>(0.10,new TGraph()));
    grobs_up.push_back(pair<double,TGraph*>(0.05,new TGraph()));
    grobs_up.push_back(pair<double,TGraph*>(0.01,new TGraph()));

    grobs_dw.push_back(pair<double,TGraph*>(1.00,new TGraph()));
    grobs_dw.push_back(pair<double,TGraph*>(0.75,new TGraph()));
    grobs_dw.push_back(pair<double,TGraph*>(0.50,new TGraph()));
    grobs_dw.push_back(pair<double,TGraph*>(0.30,new TGraph()));
    grobs_dw.push_back(pair<double,TGraph*>(0.25,new TGraph()));
    grobs_dw.push_back(pair<double,TGraph*>(0.20,new TGraph()));
    grobs_dw.push_back(pair<double,TGraph*>(0.10,new TGraph()));
    grobs_dw.push_back(pair<double,TGraph*>(0.05,new TGraph()));
    grobs_dw.push_back(pair<double,TGraph*>(0.01,new TGraph()));

    for(int itree = 0; itree < treeList.size(); itree++)
      fillLimitGraphs(treeList.at(itree),grexp.at(itree).second,grexp_up.at(itree).second,grexp_dw.at(itree).second,
		      grobs.at(itree).second,grobs_up.at(itree).second,grobs_dw.at(itree).second,
		      medOverDM,useDMMass);
    
  }
  
  //// make a spline to further smooth                                                                                                                                                          
  vector<pair<double,TSpline3*> > splineexp;
  vector<pair<double,TSpline3*> > splineexp_up;
  vector<pair<double,TSpline3*> > splineexp_dw;
  vector<pair<double,TSpline3*> > splineobs;
  vector<pair<double,TSpline3*> > splineobs_up;
  vector<pair<double,TSpline3*> > splineobs_dw;
  for(auto gr : grexp){
    TString name = Form("splineexp_gq_%f",gr.first);
    name.ReplaceAll(".","p");
    splineexp.push_back(pair<double,TSpline3*>(gr.first,new TSpline3(name,gr.second->GetX(),gr.second->GetY(),gr.second->GetN())));
  }
  for(auto gr : grexp_up){
    TString name = Form("splineexp_up_gq_%f",gr.first);
    name.ReplaceAll(".","p");
    splineexp_up.push_back(pair<double,TSpline3*>(gr.first,new TSpline3(name,gr.second->GetX(),gr.second->GetY(),gr.second->GetN())));
  }
  for(auto gr : grexp_dw){
    TString name = Form("splineexp_dw_gq_%f",gr.first);
    name.ReplaceAll(".","p");
    splineexp_dw.push_back(pair<double,TSpline3*>(gr.first,new TSpline3(name,gr.second->GetX(),gr.second->GetY(),gr.second->GetN())));
  }

  for(auto gr : grobs){
    TString name = Form("splineobs_gq_%f",gr.first);
    name.ReplaceAll(".","p");
    splineobs.push_back(pair<double,TSpline3*>(gr.first,new TSpline3(name,gr.second->GetX(),gr.second->GetY(),gr.second->GetN())));
  }
  for(auto gr : grobs_up){
    TString name = Form("splineobs_up_gq_%f",gr.first);
    name.ReplaceAll(".","p");
    splineobs_up.push_back(pair<double,TSpline3*>(gr.first,new TSpline3(name,gr.second->GetX(),gr.second->GetY(),gr.second->GetN())));
  }
  for(auto gr : grobs_dw){
    TString name = Form("splineobs_dw_gq_%f",gr.first);
    name.ReplaceAll(".","p");
    splineobs_dw.push_back(pair<double,TSpline3*>(gr.first,new TSpline3(name,gr.second->GetX(),gr.second->GetY(),gr.second->GetN())));
  }

  // fill coupling
  TGraph2D* grexp_coupling = new TGraph2D();
  TGraph2D* grexp_coupling_up = new TGraph2D();
  TGraph2D* grexp_coupling_dw = new TGraph2D();
  TGraph2D* grobs_coupling = new TGraph2D();
  TGraph2D* grobs_coupling_up = new TGraph2D();
  TGraph2D* grobs_coupling_dw = new TGraph2D();

  int npoints;  
  int ngraph = 0;
  for(auto gr : grexp){
    int ipoint = 0;
    for(int i = 0; i < gr.second->GetN(); i++){
      double x,y;
      gr.second->GetPoint(ipoint,x,y);
      grexp_coupling->SetPoint(npoints,x,gr.first,splineexp.at(ngraph).second->Eval(x));
      ipoint++;
      npoints++;
    }
    ngraph++;
  }

  ngraph = 0;
  npoints = 0;
  for(auto gr : grexp_up){
    int ipoint = 0;
    for(int i = 0; i < gr.second->GetN(); i++){
      double x,y;
      gr.second->GetPoint(ipoint,x,y);
      grexp_coupling_up->SetPoint(npoints,x,gr.first,splineexp_up.at(ngraph).second->Eval(x));
      ipoint++;
      npoints++;
    }
    ngraph++;
  }
  
  ngraph = 0;
  npoints = 0;
  for(auto gr : grexp_dw){
    int ipoint = 0;
    for(int i = 0; i < gr.second->GetN(); i++){
      double x,y;
      gr.second->GetPoint(ipoint,x,y);
      grexp_coupling_dw->SetPoint(npoints,x,gr.first,splineexp_dw.at(ngraph).second->Eval(x));
      ipoint++;
      npoints++;
    }
    ngraph++;
  }

  float min_xaxis = 9999;
  float max_xaxis = -1;
  ngraph = 0;
  npoints = 0;
  for(auto gr : grobs){
    int ipoint = 0;
    for(int i = 0; i < gr.second->GetN(); i++){
      double x,y;
      gr.second->GetPoint(ipoint,x,y);
      grobs_coupling->SetPoint(npoints,x,gr.first,splineobs.at(ngraph).second->Eval(x));
      ipoint++;
      if(x < min_xaxis) min_xaxis = x;
      if(x > max_xaxis) max_xaxis = x;
      npoints++;
    }
    ngraph++;
  }

  ngraph = 0;
  npoints = 0;
  for(auto gr : grobs_up){
    int ipoint = 0;
    for(int i = 0; i < gr.second->GetN(); i++){
      double x,y;
      gr.second->GetPoint(ipoint,x,y);
      grobs_coupling_up->SetPoint(npoints,x,gr.first,splineobs_up.at(ngraph).second->Eval(x));
      ipoint++;
      npoints++;
    }
    ngraph++;
  }

  ngraph = 0;
  npoints = 0;
  for(auto gr : grobs_dw){
    int ipoint = 0;
    for(int i = 0; i < gr.second->GetN(); i++){
      double x,y;
      gr.second->GetPoint(ipoint,x,y);
      grobs_coupling_dw->SetPoint(npoints,x,gr.first,splineobs_dw.at(ngraph).second->Eval(x));
      ipoint++;
      npoints++;
    }
    ngraph++;
  }

  /////
  TH2D* hexp_coupling = new TH2D("hexp_coupling", "",nbinsX,min_xaxis,max_xaxis,nCoupling,minCoupling,maxCoupling);
  TH2D* hobs_coupling = new TH2D("hobs_coupling", "",nbinsX,min_xaxis,max_xaxis,nCoupling,minCoupling,maxCoupling);
  TH2D* hexp_coupling_up = new TH2D("hexp_coupling_up", "",nbinsX,min_xaxis,max_xaxis,nCoupling,minCoupling,maxCoupling);
  TH2D* hobs_coupling_up = new TH2D("hobs_coupling_up", "",nbinsX,min_xaxis,max_xaxis,nCoupling,minCoupling,maxCoupling);
  TH2D* hexp_coupling_dw = new TH2D("hexp_coupling_dw", "",nbinsX,min_xaxis,max_xaxis,nCoupling,minCoupling,maxCoupling);
  TH2D* hobs_coupling_dw = new TH2D("hobs_coupling_dw", "",nbinsX,min_xaxis,max_xaxis,nCoupling,minCoupling,maxCoupling);

  for(int i = 1; i < hexp_coupling->GetNbinsX(); i++){
    for(int j = 1; j < hexp_coupling->GetNbinsY(); j++){
      hexp_coupling->SetBinContent(i,j,grexp_coupling->Interpolate(hexp_coupling->GetXaxis()->GetBinCenter(i),hexp_coupling->GetYaxis()->GetBinCenter(j)));
      hobs_coupling->SetBinContent(i,j,grobs_coupling->Interpolate(hobs_coupling->GetXaxis()->GetBinCenter(i),hobs_coupling->GetYaxis()->GetBinCenter(j)));
      hexp_coupling_up->SetBinContent(i,j,grexp_coupling_up->Interpolate(hexp_coupling_up->GetXaxis()->GetBinCenter(i),hexp_coupling_up->GetYaxis()->GetBinCenter(j)));
      hexp_coupling_dw->SetBinContent(i,j,grexp_coupling_dw->Interpolate(hexp_coupling_dw->GetXaxis()->GetBinCenter(i),hexp_coupling_dw->GetYaxis()->GetBinCenter(j)));
      hobs_coupling_up->SetBinContent(i,j,grobs_coupling_up->Interpolate(hobs_coupling_up->GetXaxis()->GetBinCenter(i),hobs_coupling_up->GetYaxis()->GetBinCenter(j)));
      hobs_coupling_dw->SetBinContent(i,j,grobs_coupling_dw->Interpolate(hobs_coupling_dw->GetXaxis()->GetBinCenter(i),hobs_coupling_dw->GetYaxis()->GetBinCenter(j)));
    }
  }

  // contour
  hexp_coupling->Smooth();
  hobs_coupling->Smooth();
  hexp_coupling_up->Smooth();
  hobs_coupling_up->Smooth();
  hexp_coupling_dw->Smooth();
  hobs_coupling_dw->Smooth();

  // extend the relic density line                                                                                                                                                                   
  TGraph* relic_graph_g1_ext = new TGraph();
  TGraph* relic_graph_g2_ext = new TGraph();

  int ipoint = 0;
  for(int i = 1; i <= hexp_coupling->GetNbinsX(); i++){
    relic_graph_g1_ext->SetPoint(ipoint,hexp_coupling->GetXaxis()->GetBinCenter(i),relic_graph_g1->Eval(hexp_coupling->GetXaxis()->GetBinCenter(i)));
    relic_graph_g2_ext->SetPoint(ipoint,hexp_coupling->GetXaxis()->GetBinCenter(i),relic_graph_g2->Eval(hexp_coupling->GetXaxis()->GetBinCenter(i)));
    ipoint++;
  }

  ////////////                                                                                                                                                                                         
  for(int i = 1; i <= hexp_coupling->GetNbinsX(); i++){
    for(int j = 1; j <= hexp_coupling->GetNbinsY(); j++){

      if(hexp_coupling -> GetBinContent(i,j) <= 0) hexp_coupling->SetBinContent(i,j,maxZ);
      if(hobs_coupling -> GetBinContent(i,j) <= 0) hobs_coupling->SetBinContent(i,j,maxZ);
      if(hexp_coupling_up -> GetBinContent(i,j) <= 0) hexp_coupling_up->SetBinContent(i,j,maxZ);
      if(hobs_coupling_up -> GetBinContent(i,j) <= 0) hobs_coupling_up->SetBinContent(i,j,maxZ);
      if(hexp_coupling_dw -> GetBinContent(i,j) <= 0) hexp_coupling_dw->SetBinContent(i,j,maxZ);
      if(hobs_coupling_dw -> GetBinContent(i,j) <= 0) hobs_coupling_dw->SetBinContent(i,j,maxZ);

      if(hexp_coupling -> GetBinContent(i,j) > maxZ) hexp_coupling->SetBinContent(i,j,maxZ);
      if(hobs_coupling -> GetBinContent(i,j) > maxZ) hobs_coupling->SetBinContent(i,j,maxZ);
      if(hexp_coupling_up -> GetBinContent(i,j) > maxZ) hexp_coupling_up->SetBinContent(i,j,maxZ);
      if(hobs_coupling_up -> GetBinContent(i,j) > maxZ) hobs_coupling_up->SetBinContent(i,j,maxZ);
      if(hexp_coupling_dw -> GetBinContent(i,j) > maxZ) hexp_coupling_dw->SetBinContent(i,j,maxZ);
      if(hobs_coupling_dw -> GetBinContent(i,j) > maxZ) hobs_coupling_dw->SetBinContent(i,j,maxZ);

      if(hexp_coupling -> GetBinContent(i,j) < minZ) hexp_coupling->SetBinContent(i,j,minZ);
      if(hobs_coupling -> GetBinContent(i,j) < minZ) hobs_coupling->SetBinContent(i,j,minZ);
      if(hexp_coupling_up -> GetBinContent(i,j) < minZ) hexp_coupling_up->SetBinContent(i,j,minZ);
      if(hobs_coupling_up -> GetBinContent(i,j) < minZ) hobs_coupling_up->SetBinContent(i,j,minZ);
      if(hexp_coupling_dw -> GetBinContent(i,j) < minZ) hexp_coupling_dw->SetBinContent(i,j,minZ);
      if(hobs_coupling_dw -> GetBinContent(i,j) < minZ) hobs_coupling_dw->SetBinContent(i,j,minZ);
    }
  }

  ////////////////                                                                                                                                                                                    
  TH2* hexp2 = (TH2*) hexp_coupling->Clone("hexp2");
  TH2* hobs2 = (TH2*) hobs_coupling->Clone("hobs2");
  TH2* hexp2_up = (TH2*) hexp_coupling_up->Clone("hexp2");
  TH2* hobs2_up = (TH2*) hobs_coupling_up->Clone("hobs2");
  TH2* hexp2_dw = (TH2*) hexp_coupling_dw->Clone("hexp2");
  TH2* hobs2_dw = (TH2*) hobs_coupling_dw->Clone("hobs2");

  //////////                                                                                                                                                                                           
  double contours[1]; contours[0]=1;
  hexp2->SetContour(1,contours);
  hobs2->SetContour(1,contours);
  hexp2_up->SetContour(1,contours);
  hobs2_up->SetContour(1,contours);
  hexp2_dw->SetContour(1,contours);
  hobs2_dw->SetContour(1,contours);

  // All the plotting and cosmetics                                                                                                                                                                   
  TCanvas* canvas = new TCanvas("canvas", "canvas",650,600);
  canvas->SetRightMargin(0.15);
  canvas->SetLeftMargin(0.13);
  canvas->SetLogz();

  TH1* frame = canvas->DrawFrame(minX,minCoupling,maxX,maxCoupling,"");
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

  hexp2->GetZaxis()->SetLabelSize(0);
  hexp2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp = produceContour(reductionForContour);

  hexp2_up->GetZaxis()->SetLabelSize(0);
  hexp2_up->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp_up = produceContour(reductionForContour);

  hexp2_dw->GetZaxis()->SetLabelSize(0);
  hexp2_dw->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp_dw = produceContour(reductionForContour);

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

  hobs_coupling->Draw("COLZ SAME");

  if(addRelicDensity){
    relic_graph_g1_ext->SetLineColor(kGreen+3);
    relic_graph_g2_ext->SetLineColor(kGreen+3);
    relic_graph_g1_ext->SetLineWidth(-802);
    relic_graph_g2_ext->SetLineWidth(802);
    relic_graph_g1_ext->SetFillStyle(3005);
    relic_graph_g2_ext->SetFillStyle(3005);
    relic_graph_g1_ext->SetFillColor(kGreen+3);
    relic_graph_g2_ext->SetFillColor(kGreen+3);
    relic_graph_g1_ext->Draw("L SAME");
    //    relic_graph_g2_ext->Draw("L SAME");                                                                                                                                                           
  }

  contour_exp_up->SetLineColor(kBlack);
  contour_exp_up->SetLineWidth(2);
  contour_exp_up->SetLineStyle(7);
  contour_exp_up->Draw("Lsame");

  contour_exp_dw->SetLineColor(kBlack);
  contour_exp_dw->SetLineWidth(2);
  contour_exp_dw->SetLineStyle(7);
  contour_exp_dw->Draw("Lsame");
  
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

  if(not addPreliminary)
    CMS_lumi(canvas,"35.9",true,true,false,0,-0.09);
  else
    CMS_lumi(canvas,"35.9",true,false,false,0,-0.09);

  
  TLegend* tex = new TLegend(0.15,0.83,0.65,0.89);
  tex->SetFillColor(kWhite);
  tex->SetFillStyle(1001);
  tex->SetBorderSize(0);
  tex->SetTextFont(42);
  tex->SetTextSize(0.027);
  tex->SetHeader("#bf{Axial med, Dirac DM, g_{DM} = 1, m_{med} = 3 #times m_{DM}}");
  tex->Draw("same");

  TLegend *leg = new TLegend(0.15,0.58,0.52,0.81);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(1001);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->SetTextSize(0.0243902);
  leg->AddEntry(contour_exp,"Median expected 95% CL","L");
  leg->AddEntry(contour_exp_up,"68% expected","L");
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
  tex2->DrawLatex(0.975,0.55,"Observed #sigma_{95% CL}/#sigma_{th}");

  canvas->RedrawAxis("sameaxis");

  canvas->SaveAs((outputDIR+"/scan_vector_mdm_vs_gq_"+string(energy)+"TeV.pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_vs_gq_"+string(energy)+"TeV.png").c_str());

  canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_vs_gq_"+string(energy)+"TeV_log.pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_vs_gq_"+string(energy)+"TeV_log.png").c_str());

}
