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
static float minX = 1;
static float minY = 1;
static float maxX = 2500;
static float maxY = 1200;
static float nCoupling   = 25;
static float minCoupling = 0.1;
static float maxCoupling = 1;
static float minZ = 0.01;
static float maxZ = 10;
static int   reductionForContour = 20;
static bool  addPreliminary = false;
static bool  whiteOut       = true;
static bool  skipPoints     = true;

void fillLimitGraphs(TTree* tree,TGraph2D* grexp, TGraph2D* grobs, const float & coupling){
  
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
    
    if (quantile == -1) { // observed
      grobs->SetPoint(obscounter, double(medmass), double(dmmass), limit);
      obscounter++;      
    }
  }
  
  tree->ResetBranchAddresses();
}


///////////////
void plotVectorCoupling(string outputDIR, bool useDMMass = false, string energy = "13") {

  string inputFileName1 = "LimitsForPaperCoupling/higgsCombine_COMB_Vector_1p0.root";
  string inputFileName2 = "LimitsForPaperCoupling/higgsCombine_COMB_Vector_0p25.root";
  string inputFileName3 = "LimitsForPaperCoupling/higgsCombine_COMB_Vector_0p1.root";

  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  int ncontours = 999;
  gStyle->SetPalette(kBird);  
  gStyle->SetNumberContours(ncontours);
  
  // This is where all the plots are made
  TFile *file1  = TFile::Open(inputFileName1.c_str(),"READ");
  TTree *tree1  = (TTree*)file1->Get("limit");
  TFile *file2  = TFile::Open(inputFileName2.c_str(),"READ");
  TTree *tree2  = (TTree*)file2->Get("limit");
  TFile *file3  = TFile::Open(inputFileName3.c_str(),"READ");
  TTree *tree3  = (TTree*)file3->Get("limit");

  vector<pair<double,TGraph2D*> > grexp; 
  grexp.push_back(pair<double,TGraph2D*>(1.00,new TGraph2D())); 
  grexp.push_back(pair<double,TGraph2D*>(0.25,new TGraph2D())); 
  grexp.push_back(pair<double,TGraph2D*>(0.10,new TGraph2D())); 
  vector<pair<double,TGraph2D*> > grobs; 
  grobs.push_back(pair<double,TGraph2D*>(1.00,new TGraph2D())); 
  grobs.push_back(pair<double,TGraph2D*>(0.25,new TGraph2D())); 
  grobs.push_back(pair<double,TGraph2D*>(0.10,new TGraph2D())); 
  
  // fill 2D graphs with limits values
  fillLimitGraphs(tree1,grexp.at(0).second,grobs.at(0).second,grexp.at(0).first);
  fillLimitGraphs(tree2,grexp.at(1).second,grobs.at(1).second,grexp.at(1).first);
  fillLimitGraphs(tree3,grexp.at(2).second,grobs.at(2).second,grexp.at(2).first);

  TH2D* hexp  = new TH2D("hexp", "",      nbinsX, minX, maxX, nbinsY, minY, maxY);
  TH2D* hobs  = new TH2D("hobs", "",      nbinsX, minX, maxX, nbinsY, minY, maxY);

  // project in 1 dimention for each coupling
  vector<pair<double,TGraph*> > grexp_project;
  vector<pair<double,TGraph*> > grobs_project;
  
  for(auto gr : grexp)
    grexp_project.push_back(pair<double,TGraph*>(gr.first,new TGraph()));
  for(auto gr : grobs)
    grobs_project.push_back(pair<double,TGraph*>(gr.first,new TGraph()));
  
  for(int igraph = 0; igraph < grexp_project.size(); igraph++){
    int ipoint_exp = 0;
    int ipoint_obs = 0;
    for(int i = 1; i   <= nbinsX; i++){

      if(not useDMMass){

	/// fix bad limits for coupling 1.0
	if(grexp_project.at(igraph).first == 1.0){ 

	  if((hexp->GetXaxis()->GetBinCenter(i) < 700 or hexp->GetXaxis()->GetBinCenter(i) > 780) and (hexp->GetXaxis()->GetBinCenter(i) < 450 or hexp->GetXaxis()->GetBinCenter(i) > 500)){

	    grexp_project.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						      grexp.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));	    
	    ipoint_exp++;
	  }

	  if((hobs->GetXaxis()->GetBinCenter(i) < 720 or hobs->GetXaxis()->GetBinCenter(i) > 800) and (hobs->GetXaxis()->GetBinCenter(i) < 450 or hobs->GetXaxis()->GetBinCenter(i) > 500)){
            grobs_project.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),
						      grobs.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
            ipoint_obs++;
          }
	}

	// skip some spikes for gq  = 0.25
	else if(grexp_project.at(igraph).first == 0.25){
	  if((hexp->GetXaxis()->GetBinCenter(i) < 600 or hexp->GetXaxis()->GetBinCenter(i) > 800) and (hexp->GetXaxis()->GetBinCenter(i) < 60 or hexp->GetXaxis()->GetBinCenter(i) > 100) and (hexp->GetXaxis()->GetBinCenter(i) < 30 or hexp->GetXaxis()->GetBinCenter(i) > 50) and (hexp->GetXaxis()->GetBinCenter(i) < 120 or hexp->GetXaxis()->GetBinCenter(i) > 300)){

	    grexp_project.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),
						      grexp.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	    ipoint_exp++;
	  }
	  
	  
	  if((hobs->GetXaxis()->GetBinCenter(i) < 600 or hobs->GetXaxis()->GetBinCenter(i) > 800) and (hobs->GetXaxis()->GetBinCenter(i) < 60 or hobs->GetXaxis()->GetBinCenter(i) > 100) and (hobs->GetXaxis()->GetBinCenter(i) < 30 or hobs->GetXaxis()->GetBinCenter(i) > 50) and (hobs->GetXaxis()->GetBinCenter(i) < 120 or hobs->GetXaxis()->GetBinCenter(i) > 300)){
	    grobs_project.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),grobs.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
	    ipoint_obs++;
	  }
	}
	else if(grexp_project.at(igraph).first == 0.1){

	  if(hexp->GetXaxis()->GetBinCenter(i) < 700 or hexp->GetXaxis()->GetBinCenter(i) > 780){
	    grexp_project.at(igraph).second->SetPoint(ipoint_exp,hexp->GetXaxis()->GetBinCenter(i),grexp.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	    ipoint_exp++;
	  }

	  if(hobs->GetXaxis()->GetBinCenter(i) < 700 or hobs->GetXaxis()->GetBinCenter(i) > 780){
            grobs_project.at(igraph).second->SetPoint(ipoint_obs,hobs->GetXaxis()->GetBinCenter(i),grobs.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
            ipoint_obs++;
	  }
	}
      }
      else{
	grexp_project.at(igraph).second->SetPoint(i,hexp->GetXaxis()->GetBinCenter(i)/3,grexp.at(igraph).second->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetXaxis()->GetBinCenter(i)/3.));
	grobs_project.at(igraph).second->SetPoint(i,hobs->GetXaxis()->GetBinCenter(i)/3,grobs.at(igraph).second->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetXaxis()->GetBinCenter(i)/3.));
      }
    }
  }

  TFile* output = new TFile("output.root","RECREATE");
  output->cd();
  grobs_project.at(0).second->Write("gq1");
  grobs_project.at(1).second->Write("g025");
  grobs_project.at(2).second->Write("gq01");
  
  //// make a spline                                                                                                                                                                                   
  vector<pair<double,TSpline3*> > splineexp;
  vector<pair<double,TSpline3*> > splineobs;
  for(auto gr : grexp_project){
    TString name = Form("splineexp_gq_%f",gr.first);
    name.ReplaceAll(".","p");
    splineexp.push_back(pair<double,TSpline3*>(gr.first,new TSpline3(name,gr.second->GetX(),gr.second->GetY(),gr.second->GetN())));
  }

  for(auto gr : grobs_project){
    TString name = Form("splineobs_gq_%f",gr.first);
    name.ReplaceAll(".","p");
    splineobs.push_back(pair<double,TSpline3*>(gr.first,new TSpline3(name,gr.second->GetX(),gr.second->GetY(),gr.second->GetN())));
  }

  // fill coupling
  TGraph2D* grexp_coupling = new TGraph2D();
  TGraph2D* grobs_coupling = new TGraph2D();

  int npoints;  
  int ngraph = 0;
  for(auto gr : grexp_project){
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

  float min_xaxis = 9999;
  float max_xaxis = -1;
  int npointsX;  
  ngraph = 0;
  for(auto gr : grobs_project){
    int ipoint = 0;
    for(int i = 0; i < gr.second->GetN(); i++){
      double x,y;
      gr.second->GetPoint(ipoint,x,y);
      grobs_coupling->SetPoint(npoints,x,gr.first,splineobs.at(ngraph).second->Eval(x));
      ipoint++;
      npoints++;
      if(x < min_xaxis) min_xaxis = x;
      if(x > max_xaxis) max_xaxis = x;
    }
    npointsX = ipoint;
    ngraph++;
  }

  /////
  TH2D* hexp_coupling = new TH2D("hexp_coupling", "",npoints,min_xaxis,max_xaxis,nCoupling,minCoupling,maxCoupling);
  TH2D* hobs_coupling = new TH2D("hobs_coupling", "",npoints,min_xaxis,max_xaxis,nCoupling,minCoupling,maxCoupling);

  for(int i = 1; i < hexp_coupling->GetNbinsX(); i++){
    for(int j = 1; j < hexp_coupling->GetNbinsY(); j++){
      hexp_coupling->SetBinContent(i,j,grexp_coupling->Interpolate(hexp_coupling->GetXaxis()->GetBinCenter(i),hexp_coupling->GetYaxis()->GetBinCenter(j)));
      hobs_coupling->SetBinContent(i,j,grobs_coupling->Interpolate(hobs_coupling->GetXaxis()->GetBinCenter(i),hobs_coupling->GetYaxis()->GetBinCenter(j)));
    }
  }

  // contour
  hexp_coupling->Smooth();
  hobs_coupling->Smooth();

  ////////////                                                                                                                                                                                         
  for(int i = 1; i <= hexp_coupling->GetNbinsX(); i++){
    for(int j = 1; j <= hexp_coupling->GetNbinsY(); j++){

      if(hexp_coupling -> GetBinContent(i,j) <= 0) hexp_coupling->SetBinContent(i,j,maxZ);
      if(hobs_coupling -> GetBinContent(i,j) <= 0) hobs_coupling->SetBinContent(i,j,maxZ);
      if(hexp_coupling -> GetBinContent(i,j) > maxZ) hexp_coupling->SetBinContent(i,j,maxZ);
      if(hobs_coupling -> GetBinContent(i,j) > maxZ) hobs_coupling->SetBinContent(i,j,maxZ);
      if(hexp_coupling -> GetBinContent(i,j) < minZ) hexp_coupling->SetBinContent(i,j,minZ);
      if(hobs_coupling -> GetBinContent(i,j) < minZ) hobs_coupling->SetBinContent(i,j,minZ);
    }
  }


  ////////////////                                                                                                                                                                                    
  TH2* hexp2 = (TH2*) hexp_coupling->Clone("hexp2");
  TH2* hobs2 = (TH2*) hobs_coupling->Clone("hobs2");

  //////////                                                                                                                                                                                           
  double contours[1]; contours[0]=1;
  hexp2->SetContour(1,contours);
  hobs2->SetContour(1,contours);

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
  
  hexp2->GetZaxis()->SetLabelSize(0);
  hexp2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp = produceContour(reductionForContour);

  hobs2->GetZaxis()->SetLabelSize(0);
  hobs2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_obs = produceContour(reductionForContour);
  
  frame->Draw();
  hobs_coupling->Draw("COLZ SAME");
  
  contour_exp->SetLineColor(kBlack);
  contour_exp->SetLineWidth(3);
  contour_exp->Draw("Lsame");

  contour_obs->SetLineColor(kRed);
  contour_obs->SetLineWidth(3);
  contour_obs->Draw("Lsame");

  if(not addPreliminary)
    CMS_lumi(canvas,"35.9",false,true,false,0,-0.09);
  else
    CMS_lumi(canvas,"35.9",false,false,false,0,-0.09);
  
  TLegend *leg = new TLegend(0.175,0.48,0.50,0.75);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->AddEntry(contour_exp,"Median expected 95% CL","L");
  leg->AddEntry(contour_obs,"Observed 95% CL","L");
  leg->Draw("same");

  TLatex * tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.030);
  tex->DrawLatex(0.175,0.80,"#bf{Vector med, Dirac DM, g_{DM} = 1, m_{med} = 3m_{DM}}");

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
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_vs_gq_"+string(energy)+"TeV.root").c_str());

  canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_vs_gq_"+string(energy)+"TeV_log.pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_vs_gq_"+string(energy)+"TeV_log.png").c_str());
  
}
