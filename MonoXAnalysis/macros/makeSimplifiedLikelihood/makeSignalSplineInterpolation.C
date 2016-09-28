#include "../CMS_lumi.h"

//////////////
vector<double> mediatorPtBinning  = {150,200,250,300,350,400,450,500,600,700,850,1000,1250};
vector<double> mediatorEtaBinning = {0.,1,2,3,5};
vector<double> mediatorMassBinning = {};
vector<double> darkMatterMassBinning = {};
//////////////
static double minMedMass = 100.;
//////////////
static float massStep = 25;
//////////////
static float ptStep  = 50;
static float etaStep = 1;

// use the small trees produced by makeTreesForInterpolation/makeSmallGenTree.C and used for gen-interpolation to make a spline
void makeSignalSplineInterpolation (string inputFileName, string outputDIR, string category){
  
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  if(category != "monojet" and category != "monoV"){
    cerr<<"Problem with category definition --> should be monojet or monoV --> return "<<endl;
    return;
  }
    
  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TTree* inputTree  = (TTree*) inputFile->Get("tree");

  gSystem->Load("libHiggsAnalysisCombinedLimit.so"); 
  // for each mass point we need to store the correct efficiency after storing the related histograms for numerator and denominator
  TTreeReader reader(inputTree);
  TTreeReaderValue<int> id (reader,"id");
  TTreeReaderValue<double> genMediatorPt (reader,"genMediatorPt");
  TTreeReaderValue<double> genMediatorEta (reader,"genMediatorEta");
  TTreeReaderValue<double> genMediatorMass (reader,"genMediatorMass");
  TTreeReaderValue<double> genX1Pt (reader,"genX1Pt");
  TTreeReaderValue<double> genX1Eta (reader,"genX1Eta");
  TTreeReaderValue<double> genX1Mass (reader,"genX1Mass");
  TTreeReaderValue<double> genX2Pt (reader,"genX2Pt");
  TTreeReaderValue<double> genX2Eta (reader,"genX2Eta");
  TTreeReaderValue<double> genX2Mass (reader,"genX2Mass");
  TTreeReaderValue<double> pfMetPt (reader,"pfMetPt");
  TTreeReaderValue<double> weightPU (reader,"weightPU");
  TTreeReaderValue<double> weightTurnOn (reader,"weightTurnOn");
  TTreeReaderValue<double> genWeight (reader,"genWeight");

  // make map
  // first loop to identify mediator points)
  while(reader.Next()){
    if(std::find(mediatorMassBinning.begin(),mediatorMassBinning.end(),*genMediatorMass) == mediatorMassBinning.end())
      mediatorMassBinning.push_back(*genMediatorMass);
    if(std::find(darkMatterMassBinning.begin(),darkMatterMassBinning.end(),*genX1Mass) == darkMatterMassBinning.end())
      darkMatterMassBinning.push_back(*genX1Mass);
  }

  std::sort(mediatorMassBinning.begin(),mediatorMassBinning.end());
  std::sort(darkMatterMassBinning.begin(),darkMatterMassBinning.end());

  TH2F* massDenominator = new TH2F("massDenominator","",mediatorMassBinning.size()-1,&mediatorMassBinning[0],darkMatterMassBinning.size()-1,&darkMatterMassBinning[0]);
  TH2F* massNumerator = new TH2F("massNumerator","",mediatorMassBinning.size()-1,&mediatorMassBinning[0],darkMatterMassBinning.size()-1,&darkMatterMassBinning[0]);
  TH2F* mediatorDenominator = new TH2F("mediatorDenominator","",mediatorPtBinning.size()-1,&mediatorPtBinning[0],mediatorEtaBinning.size()-1,&mediatorEtaBinning[0]);
  TH2F* mediatorNumerator = new TH2F("mediatorNumerator","",mediatorPtBinning.size()-1,&mediatorPtBinning[0],mediatorEtaBinning.size()-1,&mediatorEtaBinning[0]);
  massDenominator->Sumw2();
  massNumerator->Sumw2();
  mediatorDenominator->Sumw2();
  mediatorNumerator->Sumw2();

  reader.SetEntry(0);
  while(reader.Next()){

    massDenominator->Fill(*genMediatorMass,*genX1Mass,(*weightPU)*(*genWeight));
    mediatorDenominator->Fill(*genMediatorPt,fabs(*genMediatorEta),(*weightPU)*(*genWeight));
    if(category == "monojet" and *id == 1){
      massNumerator->Fill(*genMediatorMass,*genX1Mass,(*weightPU)*(*genWeight)*(*weightTurnOn));
      if(*genMediatorMass > minMedMass)
	mediatorNumerator->Fill(*genMediatorPt,fabs(*genMediatorEta),(*weightPU)*(*genWeight)*(*weightTurnOn));
    }
    else if(category == "monoV" and *id == 2){
      massNumerator->Fill(*genMediatorMass,*genX1Mass,(*weightPU)*(*genWeight)*(*weightTurnOn));
      if(*genMediatorMass > minMedMass)
	mediatorNumerator->Fill(*genMediatorPt,fabs(*genMediatorEta),(*weightPU)*(*genWeight)*(*weightTurnOn));
    }    
  }

  // calculate the efficiency
  TH2F* efficiencyMass = (TH2F*) massNumerator->Clone("efficiencyMass");
  efficiencyMass->Reset();
  efficiencyMass->Divide(massNumerator,massDenominator,1,1,"B");
  TH2F* efficiencyMediator = (TH2F*) mediatorNumerator->Clone("efficiencyMediator");
  efficiencyMediator->Reset();
  efficiencyMediator->Divide(mediatorNumerator,mediatorDenominator,1,1,"B");

  TCanvas* canvas = new TCanvas("canvas","",650,600);
  canvas->cd();
  canvas->SetLeftMargin(0.12);
  canvas->SetRightMargin(0.18);
  canvas->SetLogx();
  efficiencyMass->GetXaxis()->SetTitle("Mediator Mass [GeV]");
  efficiencyMass->GetYaxis()->SetTitle("DM Mass [GeV]");
  efficiencyMass->GetYaxis()->SetTitleOffset(1.1);
  efficiencyMass->GetZaxis()->SetTitle("Signal Efficiency");
  efficiencyMass->GetZaxis()->SetTitleOffset(1.1);
  efficiencyMass->Draw("colz text");
  canvas->SaveAs((outputDIR+"/efficiencyMass_2D.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/efficiencyMass_2D.pdf").c_str(),"pdf");

  canvas->SetLogx(0);
  efficiencyMediator->GetXaxis()->SetTitle("Mediator p_{T} [GeV]");
  efficiencyMediator->GetYaxis()->SetTitle("Mediator |#eta|");
  efficiencyMediator->GetZaxis()->SetTitle("Signal Efficiency");
  efficiencyMediator->Draw("colz text");
  canvas->SaveAs((outputDIR+"/efficiencyMediator_2D.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/efficiencyMediator_2D.pdf").c_str(),"pdf");
 

  TGraph2D *graphEfficiencyMass = new TGraph2D();
  TGraph2D *graphEfficiencyMediator = new TGraph2D();
  graphEfficiencyMass->SetMarkerSize(0.8);
  graphEfficiencyMass->SetMarkerStyle(20);
  graphEfficiencyMediator->SetMarkerSize(0.8);
  graphEfficiencyMediator->SetMarkerStyle(20);
 
  int nPoint = 0;
  for(int binX = 0; binX < efficiencyMass->GetNbinsX()+1; binX++){
    for(int binY = 0; binY < efficiencyMass->GetNbinsY()+1; binY++){
      if(efficiencyMass->GetBinContent(binX+1,binY+1) == 0) continue;
      graphEfficiencyMass->SetPoint(nPoint,efficiencyMass->GetXaxis()->GetBinLowEdge(binX+1),efficiencyMass->GetYaxis()->GetBinLowEdge(binY+1),efficiencyMass->GetBinContent(binX+1,binY+1));
      nPoint++;
    }
  }
  nPoint = 0;
  for(int binX = 0; binX < efficiencyMediator->GetNbinsX()+1; binX++){
    for(int binY = 0; binY < efficiencyMediator->GetNbinsY()+1; binY++){
      if(efficiencyMediator->GetBinContent(binX+1,binY+1) == 0) continue;
      graphEfficiencyMediator->SetPoint(nPoint,efficiencyMediator->GetXaxis()->GetBinCenter(binX+1),efficiencyMediator->GetYaxis()->GetBinCenter(binY+1),efficiencyMediator->GetBinContent(binX+1,binY+1));
      nPoint++;
    }
  }

  // Create the spline
  RooRealVar* mmed  = new RooRealVar("mmed","",(mediatorMassBinning.back()+mediatorMassBinning.front())/2,mediatorMassBinning.front(),mediatorMassBinning.back());
  RooRealVar* mdm   = new RooRealVar("mdm","",(darkMatterMassBinning.front()+darkMatterMassBinning.back())/2,darkMatterMassBinning.front(),darkMatterMassBinning.back());
  RooRealVar* bpt   = new RooRealVar("bpt","",(mediatorPtBinning.front()+mediatorPtBinning.back())/2,mediatorPtBinning.front(),mediatorPtBinning.back());
  RooRealVar* beta  = new RooRealVar("beta","",(mediatorEtaBinning.front()+mediatorEtaBinning.back())/2,mediatorEtaBinning.front(),mediatorEtaBinning.back());
  
  //Convert histograms in TTree
  TTree* treeMass = new TTree("treeMass","treeMass");
  TTree* treeMediator = new TTree("treeMediator","treeMediator");
  float medMass = 0, dmMass = 0, bosonpt = 0, bosoneta = 0, effm = 0, effb = 0;
  treeMass->Branch("mmed",&medMass,"mmed/F");
  treeMass->Branch("mdm",&dmMass,"mdm/F");
  treeMediator->Branch("bpt",&bosonpt,"bpt/F");
  treeMediator->Branch("beta",&bosoneta,"beta/F");
  treeMass->Branch("effm",&effm,"effm/F");
  treeMediator->Branch("effb",&effb,"effb/F");

  for(int binX = 0; binX < efficiencyMass->GetNbinsX()+1; binX++){
    for(int binY = 0; binY < efficiencyMass->GetNbinsY()+1; binY++){
      if(efficiencyMass->GetBinContent(binX+1,binY+1) == 0) continue;
      medMass = 0; dmMass = 0; effm = 0; 
      medMass = efficiencyMass->GetXaxis()->GetBinLowEdge(binX+1);
      dmMass  = efficiencyMass->GetYaxis()->GetBinLowEdge(binY+1);
      effm = efficiencyMass->GetBinContent(binX+1,binY+1);
      if(effm != 0)
	treeMass->Fill();
    }
  }

  for(int binX = 0; binX < efficiencyMediator->GetNbinsX()+1; binX++){
    for(int binY = 0; binY < efficiencyMediator->GetNbinsY()+1; binY++){
      if(efficiencyMediator->GetBinContent(binX+1,binY+1) == 0) continue;
      bosonpt  = 0; bosoneta = 0; effb = 0; 
      bosonpt  = efficiencyMediator->GetXaxis()->GetBinCenter(binX+1);
      bosoneta = efficiencyMediator->GetYaxis()->GetBinCenter(binY+1);
      effb = efficiencyMediator->GetBinContent(binX+1,binY+1);
      if(effb != 0)
	treeMediator->Fill();
    }
  }

  RooArgList list1 (*mmed,*mdm);
  RooArgList list2 (*bpt,*beta);
  
  RooSplineND *spline_mmed_mdm = new RooSplineND("spline_mmed_mdm","spline_mmed_mdm",list1,treeMass,"effm",1,true);  
  RooSplineND *spline_bpt_beta = new RooSplineND("spline_bpt_beta","spline_bpt_beta",list2,treeMediator,"effb",1,true);  
  TGraph2D* splineMassGraph     = new TGraph2D();
  TGraph2D* splineMediatorGraph = new TGraph2D();

  // Create interpolated graphs
  nPoint = 0;
  for(float iBinX = mediatorMassBinning.front(); iBinX < mediatorMassBinning.back(); iBinX=iBinX+massStep){
    for(float iBinY = darkMatterMassBinning.front(); iBinY < darkMatterMassBinning.back(); iBinY=iBinY+massStep){
      mmed->setVal(iBinX);
      mdm->setVal(iBinY);
      splineMassGraph->SetPoint(nPoint,iBinX,iBinY,spline_mmed_mdm->getVal());
      nPoint++;
    }
  }

  nPoint = 0;
  for(float iBinX = mediatorPtBinning.front(); iBinX < mediatorPtBinning.back(); iBinX=iBinX+ptStep){
    for(float iBinY = mediatorEtaBinning.front(); iBinY < mediatorEtaBinning.back(); iBinY=iBinY+etaStep){
      bpt->setVal(iBinX);
      beta->setVal(iBinY);
      splineMediatorGraph->SetPoint(nPoint,iBinX,iBinY,spline_bpt_beta->getVal());
      nPoint++;
    }
  }

  canvas->SetLogx();  
  graphEfficiencyMass->Draw("P0");
  splineMassGraph->Draw("surf1 same");
  splineMassGraph->Draw("P0");
  canvas->SaveAs((outputDIR+"/efficiencyMass_2D_Graph.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/efficiencyMass_2D_Graph.pdf").c_str(),"pdf");
  canvas->SetLogx(0);
  graphEfficiencyMediator->Draw("P0");
  splineMediatorGraph->Draw("surf1 same");
  graphEfficiencyMediator->Draw("P0same");
  canvas->SaveAs((outputDIR+"/efficiencyMediator_2D_Graph.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/efficiencyMediator_2D_Graph.pdf").c_str(),"pdf");    
  
}
