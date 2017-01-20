#include <fstream>
#include "TROOT.h"
#include "../CMS_lumi.h"
#include "triggerUtils.h"

/// loop on the events to build the efficiency for each pt-eta bin
void makeTrigger(TTree* tree, // tree 
		 TH2F*  histoNum,
		 TH2F*  histoDen,
		 const bool & isMC,
		 const bool & isSingleMuon,
		 const string & outputDIR,
		 const string & triggerFlag,
		 const float & ptMin,
		 const float & ptMax,
		 const float & etaMin,
		 const float & etaMax,
		 string triggerFlag_2 = "",
		 int    nbins = 60,
		 float  xmin  = 60,
		 float  xmax  = 120){

  // passing, fail and all probes
  TH1F hpass("hpass",  "", nbins, xmin, xmax);
  TH1F hfail("hfail",  "", nbins, xmin, xmax);
  TH1F hall ("hall",   "", nbins, xmin, xmax);

  // passing and failing fractions
  TH1F hp("hp" , "", 1, xmin, xmax);
  TH1F ha("ha" , "", 1, xmin, xmax);
  TH1F hr("hr" , "", 1, xmin, xmax);

  hpass.Sumw2();
  hfail.Sumw2();
  hall.Sumw2();
  hp.Sumw2();
  ha.Sumw2();

  // pileup-re-weight from external file                                                                                                                                        
  TFile* pufile = new TFile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_35p9fb.root");
  TH1* puhist   = (TH1*) pufile->Get("puhist");

  TTreeReader reader(tree);
  TTreeReaderValue<float> mass (reader, "mass");
  TTreeReaderValue<float> eta  (reader, "eta");
  TTreeReaderValue<float> phi  (reader, "phi");
  TTreeReaderValue<float> pt   (reader, "pt");
  TTreeReaderValue<float> nvtx (reader, "nvtx");
  TTreeReaderValue<float> wgt  (reader, "wgt");
  TTreeReaderValue<int>   id   (reader, "tightid");  // single electorn and muon trigger efficiency are useful in the single and double lepton CR were we ask at least one tight lepton
  TTreeReaderValue<int>   hlt  (reader,triggerFlag.c_str());
  string hltsafe = "wgt";
  if(not isSingleMuon)
    hltsafe = "hltsafeid";
  TTreeReaderValue<int>   hltsafeid (reader,hltsafe.c_str());

  if(triggerFlag_2 == "")
    triggerFlag_2 = triggerFlag;
  TTreeReaderValue<int>   hlt2 (reader,triggerFlag_2.c_str());
  
        
  string mcTrueVar;
  string mcMassVar;
  if(isMC){
    mcTrueVar = "mcTrue";
    mcMassVar = "mcMass";
  }
  else{// just place holder to run the code
    mcTrueVar = "tightid";
    mcMassVar = "wgt";
  }

  TTreeReaderValue<int>    mcTrue  (reader, mcTrueVar.c_str());
  TTreeReaderValue<float>  mcMass  (reader, mcMassVar.c_str());
  
  // weight sum                                                                                                                                                                 
  double wgtsum = 0;
  if(isMC){
    while(reader.Next()){
      wgtsum += *wgt;
    }
  }
 
  reader.SetEntry(0);
  // loop on the event                                                                                                                                                          
  while(reader.Next()){
    Float_t puwgt = 0.;
    if (*nvtx <= 40) puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));    
    // check the probe lepton information --> if it falls inside the bins                                                                                                       
    if(*pt < ptMin or *pt > ptMax) continue;
    if(fabs(*eta) < etaMin or fabs(*eta) > etaMax) continue;    
    // if not matched to a genLepton skip --> we want to extract the true templateds                                                                                            
    if(isMC and not *mcTrue) continue;
    if(*id <= 0) continue;
    if(not isSingleMuon and *hltsafeid <= 0) continue;
    if((*hlt == 0 and *hlt2 == 0) and isMC)      
      hfail.Fill(*mass, puwgt*(*wgt)/wgtsum);
    else if((*hlt == 0 and *hlt2 == 0) and not isMC)
      hfail.Fill(*mass);    
    else if((*hlt == 1 || *hlt2 == 1) and isMC){ 
      hpass.Fill(*mass, puwgt*(*wgt)/wgtsum);
      hp.Fill(*mass, puwgt*(*wgt)/wgtsum);
    }
    else if((*hlt == 1 or *hlt2 == 1) and not isMC){      
      hpass.Fill(*mass);
      hp.Fill(*mass);
    }

    if(isMC){
      ha.Fill(*mass, puwgt*(*wgt)/wgtsum);
      hall.Fill(*mass, puwgt*(*wgt)/wgtsum);
    }
    else{
      ha.Fill(*mass);
      hall.Fill(*mass);
    }
  }
  hr.Divide(&hp,&ha,1,1,"B");
  // useful just for MC negative weights
  for (int j = 1; j <= nbins; j++) {
    if (hpass.GetBinContent(j) < 0.) hpass.SetBinContent(j, 0.0);
    if (hfail.GetBinContent(j) < 0.) hfail.SetBinContent(j, 0.0);
    if (hp.GetBinContent(j)    < 0.) hpass.SetBinContent(j, 0.0);
    if (ha.GetBinContent(j)    < 0.) hfail.SetBinContent(j, 0.0);
    if (hall.GetBinContent(j)  < 0.) hall.SetBinContent(j, 0.0);
  }

  histoNum->SetBinContent(histoNum->GetXaxis()->FindBin((etaMax+etaMin)/2),histoNum->GetYaxis()->FindBin((ptMax+ptMin)/2),hp.GetBinContent(1));
  histoNum->SetBinError(histoNum->GetXaxis()->FindBin((etaMax+etaMin)/2),histoNum->GetYaxis()->FindBin((ptMax+ptMin)/2),hp.GetBinError(1));
  histoDen->SetBinContent(histoDen->GetXaxis()->FindBin((etaMax+etaMin)/2),histoDen->GetYaxis()->FindBin((ptMax+ptMin)/2),ha.GetBinContent(1));
  histoDen->SetBinError(histoDen->GetXaxis()->FindBin((etaMax+etaMin)/2),histoDen->GetYaxis()->FindBin((ptMax+ptMin)/2),ha.GetBinError(1));
  
  if(triggerFlag_2 == "")
    std::cout << "Efficiency for "+triggerFlag+" -- pT [" << ptMin << ", " << ptMax << "], eta [" << etaMin << ", " << etaMax << "]  :  " << hr.GetBinContent(1) << " +/- " << hr.GetBinError(1) << std::endl; 
  else
    std::cout << "Efficiency for "+triggerFlag+" || "+triggerFlag_2+" -- pT [" << ptMin << ", " << ptMax << "], eta [" << etaMin << ", " << etaMax << "]  :  " << hr.GetBinContent(1) << " +/- " << hr.GetBinError(1) << std::endl; 

  pufile->Close();
  
}

// plotting the 2D map and 1D projections
void plotTriggerEfficiency(TCanvas*      canvas,
			   TEfficiency*  efficiency,
			   bool     isMuon,
			   string   outputDIR,
			   bool     isHighPt = false,
			   float    lumi   = 0.59,
			   bool     isLogz = false,
			   bool     doFit  = false){


  gPad->SetBottomMargin(0.10);
  gPad->SetRightMargin(0.20);

  TH2* histo = efficiency->CreateHistogram();  
  histo->SetName(efficiency->GetName());
  canvas->cd();
  if(isMuon){
    histo->GetYaxis()->SetTitle("Muon p_{T} [GeV]");
    histo->GetXaxis()->SetTitle("Muon |#eta|");
  }
  else{
    histo->GetYaxis()->SetTitle("Electron p_{T} [GeV]");
    histo->GetXaxis()->SetTitle("Electron |#eta|");
  }
  histo->GetZaxis()->SetTitle("Trigger Efficiency");
  gStyle->SetPaintTextFormat(".2f");
  histo->Draw("colz text");
  CMS_lumi(canvas,string(Form("%.2f",lumi)),true);
  if(isLogz)
    canvas->SetLogz();
  canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/"+string(histo->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+string(histo->GetName())+".pdf").c_str(),"pdf");
  canvas->SetLogy(0);
  
  gPad->SetRightMargin(0.06);
  canvas->SetLogz(0);
  // project by hand
  TGraphAsymmErrors* projection_pt  = new TGraphAsymmErrors();
  TF1  fitfunc ("fitfunc", ErfCB, 0, 500, 5);
  if(not isMuon){
    if(not isHighPt)
      fitfunc.SetParameters(22., 1., 2., 2.,1.);
    else
      fitfunc.SetParameters(100., 3., 5., 5.,1.);
  }
  else
    fitfunc.SetParameters(16., 5., 10., 2.,1.);
  
  for(int iBinX = 0; iBinX < histo->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < histo->GetNbinsY(); iBinY++){      
      int globalBin = histo->GetBin(iBinX+1,iBinY+1);
      projection_pt->SetPoint(iBinY+1,histo->GetYaxis()->GetBinCenter(iBinY+1),efficiency->GetEfficiency(globalBin));
      projection_pt->SetPointError(iBinY+1,histo->GetYaxis()->GetBinWidth(iBinY+1)/2,histo->GetYaxis()->GetBinWidth(iBinY+1)/2,efficiency->GetEfficiencyErrorLow(globalBin),efficiency->GetEfficiencyErrorUp(globalBin));
    }
    if(not isMuon)
      projection_pt->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
    else
      projection_pt->GetXaxis()->SetTitle("Muon p_{T} [GeV]");
    projection_pt->GetYaxis()->SetTitle("Trigger Efficiency");
    projection_pt->SetMarkerSize(1);
    projection_pt->SetMarkerStyle(20);
    projection_pt->SetMarkerColor(kBlack);
    projection_pt->SetLineColor(kBlack);
    projection_pt->Draw("AE1P");
    fitfunc.SetLineColor(kBlue);
    fitfunc.SetLineWidth(2);
    if(doFit)
      projection_pt->Fit(&fitfunc);
    CMS_lumi(canvas,string(Form("%.2f",lumi)),true);
    canvas->SaveAs((outputDIR+"/"+string(histo->GetName())+string(Form("_projection_pt_eta_%.1f_%.1f.png",histo->GetXaxis()->GetBinLowEdge(iBinX+1),
								       histo->GetXaxis()->GetBinLowEdge(iBinX+2)))).c_str(),"png");
    canvas->SaveAs((outputDIR+"/"+string(histo->GetName())+string(Form("_projection_pt_eta_%.1f_%.1f.pdf",histo->GetXaxis()->GetBinLowEdge(iBinX+1),
								       histo->GetXaxis()->GetBinLowEdge(iBinX+2)))).c_str(),"pdf");
  }
}
  
// Main function
void makeSingleLeptonTriggerEfficiency(string inputDIR, // where trees are located
				       string outputDIR, // to store plots and files
				       bool   isMC,
				       bool   isSingleMuon,
				       float  lumi = 0.59,
				       bool   doFit = false){
  
  
  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);
  
  // take the input tree --> valid for both data and MC
  TChain* tree = NULL;
  if(isSingleMuon)
    tree = new TChain("muontnptree/fitter_tree");
  else
    tree = new TChain("electrontnptree/fitter_tree");

  tree->Add((inputDIR+"/*root").c_str());

  // select binning for efficiency measurement
  vector<float> binningPt;
  vector<float> binningHighPt;
  vector<float> binningEta;
  if(isSingleMuon){
    binningPt      = {20.,25,30,40,45,50,60,70,85,100,125,150,175,200,250,300};
    binningHighPt  = {90,100,105,110,115,120,125,135,150,175,200};
    binningEta     = {-2.4,-1.2,1.2,2.4};
  }
  else{
    binningPt      = {20.,25,30,40,45,50,60,70,85,115,130,150,175,200,250,300};
    binningHighPt  = {90,100,105,110,115,120,125,150,175,200,250};
    binningEta     = {0,0.5,1,1.444,1.56,2.1,2.5};
  }
  
  //prepare the histograms for muons
  TH2F* Passing_mu20      = new TH2F("Passing_muIso20","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_muIsoTk20 = new TH2F("Passing_muIsoTk20","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_mu22      = new TH2F("Passing_muIso22","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_muIsoTk22 = new TH2F("Passing_muIsoTk22","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_mu        = new TH2F("Passing_muIso","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_muIsoTk   = new TH2F("Passing_muIsoTk","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);

  TH2F* Total_mu20      = new TH2F("Total_muIso20","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_muIsoTk20 = new TH2F("Total_muIsoTk20","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_mu22      = new TH2F("Total_muIso22","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_muIsoTk22 = new TH2F("Total_muIsoTk22","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_mu        = new TH2F("Total_muIso","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_muIsoTk   = new TH2F("Total_muIsoTk","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);

  //prepare the histograms for electrons
  TH2F* Passing_ele242p1wploose = new TH2F("Passing_ele242p1wploose","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_ele252p1wptight = new TH2F("Passing_ele252p1wptight","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_ele272p1wploose = new TH2F("Passing_ele272p1wploose","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_ele272p1wptight = new TH2F("Passing_ele272p1wptight","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_ele27wptight    = new TH2F("Passing_ele27wptight","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_ele105          = new TH2F("Passing_ele105","",binningEta.size()-1,&binningEta[0],binningHighPt.size()-1,&binningHighPt[0]);
  TH2F* Passing_ele             = new TH2F("Passing_ele","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);

  TH2F* Total_ele242p1wploose = new TH2F("Total_ele242p1wploose","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_ele252p1wptight = new TH2F("Total_ele252p1wptight","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_ele272p1wploose = new TH2F("Total_ele272p1wploose","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_ele272p1wptight = new TH2F("Total_ele272p1wptight","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_ele27wptight    = new TH2F("Total_ele27wptight","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_ele105          = new TH2F("Total_ele105","",binningEta.size()-1,&binningEta[0],binningHighPt.size()-1,&binningHighPt[0]);
  TH2F* Total_ele             = new TH2F("Total_ele","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  
  // start looping on the pt and eta bins
  for(size_t ipt = 0; ipt < binningPt.size()-1; ipt++){
    for(size_t ieta = 0; ieta < binningEta.size()-1; ieta++){
      if(isSingleMuon){
	makeTrigger(tree,Passing_mu20,Total_mu20,isMC,isSingleMuon,outputDIR,"hltmu20",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1)); // mu20
	makeTrigger(tree,Passing_muIsoTk20,Total_muIsoTk20,isMC,isSingleMuon,outputDIR,"hlttkmu20",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1)); 
	makeTrigger(tree,Passing_mu22,Total_mu22,isMC,isSingleMuon,outputDIR,"hltmu22",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1)); 
	makeTrigger(tree,Passing_muIsoTk22,Total_muIsoTk22,isMC,isSingleMuon,outputDIR,"hlttkmu22",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1)); 
	makeTrigger(tree,Passing_mu,Total_mu,isMC,isSingleMuon,outputDIR,"hltmu",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1)); 
	makeTrigger(tree,Passing_muIsoTk,Total_muIsoTk,isMC,isSingleMuon,outputDIR,"hlttkmu",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1)); 
      }
      else{
	makeTrigger(tree,Passing_ele242p1wploose,Total_ele242p1wploose,isMC,isSingleMuon,outputDIR,"hltele24eta2p1wpl",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1));
	makeTrigger(tree,Passing_ele252p1wptight,Total_ele252p1wptight,isMC,isSingleMuon,outputDIR,"hltele25eta2p1wpt",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1));
	makeTrigger(tree,Passing_ele272p1wploose,Total_ele272p1wploose,isMC,isSingleMuon,outputDIR,"hltele27eta2p1wpl",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1));
	makeTrigger(tree,Passing_ele272p1wptight,Total_ele272p1wptight,isMC,isSingleMuon,outputDIR,"hltele27eta2p1wpt",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1));
	makeTrigger(tree,Passing_ele27wptight,Total_ele27wptight,isMC,isSingleMuon,outputDIR,"hltele27wpt",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1));
	makeTrigger(tree,Passing_ele,Total_ele,isMC,isSingleMuon,outputDIR,"hltele",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1));
      }
    }
  }

  
  TEfficiency* trgeff_mu20      = NULL;
  TEfficiency* trgeff_muIsoTk20 = NULL;
  TEfficiency* trgeff_mu22      = NULL;
  TEfficiency* trgeff_muIsoTk22 = NULL;
  TEfficiency* trgeff_mu        = NULL;
  TEfficiency* trgeff_muIsoTk   = NULL;

  TEfficiency* trgeff_ele242p1wploose = NULL;
  TEfficiency* trgeff_ele252p1wptight = NULL;
  TEfficiency* trgeff_ele272p1wploose = NULL;
  TEfficiency* trgeff_ele272p1wptight = NULL;
  TEfficiency* trgeff_ele27wptight    = NULL;
  TEfficiency* trgeff_ele105          = NULL;
  TEfficiency* trgeff_ele             = NULL;

  if(isSingleMuon){
    trgeff_mu20 = new TEfficiency(*Passing_mu20,*Total_mu20);
    trgeff_mu20->SetName("trgeff_mu20");
    trgeff_muIsoTk20 = new TEfficiency(*Passing_muIsoTk20,*Total_muIsoTk20);
    trgeff_muIsoTk20->SetName("trgeff_muIsoTk20");
    trgeff_mu22 = new TEfficiency(*Passing_mu22,*Total_mu22);
    trgeff_mu22->SetName("trgeff_mu22");
    trgeff_muIsoTk22 = new TEfficiency(*Passing_muIsoTk22,*Total_muIsoTk22);
    trgeff_muIsoTk22->SetName("trgeff_muIsoTk22");
    trgeff_mu = new TEfficiency(*Passing_mu,*Total_mu);
    trgeff_mu->SetName("trgeff_mu");
    trgeff_muIsoTk = new TEfficiency(*Passing_muIsoTk,*Total_muIsoTk);
    trgeff_muIsoTk->SetName("trgeff_muIsoTk");
  }
  else{
    trgeff_ele242p1wploose = new TEfficiency(*Passing_ele242p1wploose,*Total_ele242p1wploose);
    trgeff_ele242p1wploose->SetName("trgeff_ele242p1wploose");
    trgeff_ele272p1wploose = new TEfficiency(*Passing_ele272p1wploose,*Total_ele272p1wploose);
    trgeff_ele272p1wploose->SetName("trgeff_ele272p1wploose");
    trgeff_ele252p1wptight = new TEfficiency(*Passing_ele252p1wptight,*Total_ele252p1wptight);
    trgeff_ele252p1wptight->SetName("trgeff_ele252p1wptight");
    trgeff_ele272p1wptight = new TEfficiency(*Passing_ele272p1wptight,*Total_ele272p1wptight);
    trgeff_ele272p1wptight->SetName("trgeff_ele272p1wptight");
    trgeff_ele27wptight = new TEfficiency(*Passing_ele27wptight,*Total_ele27wptight);
    trgeff_ele27wptight->SetName("trgeff_ele27wptight");
    trgeff_ele105       =  new TEfficiency(*Passing_ele105,*Total_ele105);
    trgeff_ele105->SetName("trgeff_ele105");
    trgeff_ele          =  new TEfficiency(*Passing_ele,*Total_ele);
    trgeff_ele->SetName("trgeff_ele");
  }

  // set the output name
  string name = "triggerEfficiency";
  if(isMC)
    name+="_MC";
  else
    name+="_DATA";
  if(isSingleMuon)
    name += "_SingleMuon";
  else
    name += "_SingleElectron";

    // Plotting part
  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();

  if(isSingleMuon){
    plotTriggerEfficiency(canvas,trgeff_mu20,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_muIsoTk20,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_mu22,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_muIsoTk22,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_mu,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_muIsoTk,isSingleMuon,outputDIR,false,lumi,doFit);
  }
  else{ 
    plotTriggerEfficiency(canvas,trgeff_ele242p1wploose,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele252p1wptight,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele272p1wploose,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele272p1wptight,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele105,isSingleMuon,outputDIR,true,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele27wptight,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele,isSingleMuon,outputDIR,false,lumi,doFit);
  }

  
  TFile* triggerEfficiency = new TFile((outputDIR+"/"+name+".root").c_str(),"RECREATE");
  triggerEfficiency->cd();
  if(isSingleMuon){
    if(trgeff_mu20)
      trgeff_mu20->Write();
    if(trgeff_muIsoTk20)
      trgeff_muIsoTk20->Write();
    if(trgeff_mu22)
      trgeff_mu22->Write();
    if(trgeff_muIsoTk22)
      trgeff_muIsoTk22->Write();
    if(trgeff_mu)
      trgeff_mu->Write();
    if(trgeff_muIsoTk)
      trgeff_muIsoTk->Write();
  }
  else{
    if(trgeff_ele242p1wploose)
      trgeff_ele242p1wploose->Write();
    if(trgeff_ele252p1wptight)
      trgeff_ele252p1wptight->Write();
    if(trgeff_ele272p1wploose)
      trgeff_ele272p1wploose->Write();
    if(trgeff_ele272p1wptight)
      trgeff_ele272p1wptight->Write();
    if(trgeff_ele27wptight)
      trgeff_ele27wptight->Write();
    if(trgeff_ele105)
      trgeff_ele105->Write();
    if(trgeff_ele)
      trgeff_ele->Write();
  }
  triggerEfficiency->Close();
}

