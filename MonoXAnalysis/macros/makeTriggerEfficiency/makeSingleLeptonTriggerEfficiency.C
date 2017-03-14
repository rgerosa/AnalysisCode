#include <fstream>
#include "TROOT.h"
#include "../CMS_lumi.h"
#include "triggerUtils.h"

static float muonTagPtCut = 25;
static float electronTagPtCut = 35;
static float muonTagEtaCut = 2.4;
static float electronTagEtaCut = 2.5;


/// loop on the events to build the efficiency for each pt-eta bin
void makeTrigger(TTree* tree, // tree 
		 TH2F*  histoNum,
		 TH2F*  histoDen,
		 const bool & isMC,
		 const bool & isSingleMuon,
		 const string & outputDIR,
		 const string & triggerFlag,
		 const vector<float> & ptBin,
		 const vector<float> & etaBin,
		 const bool & isabseta = false,
		 string triggerFlagAlt = "",		 
		 int    nbins = 60,
		 float  xmin  = 60,
		 float  xmax  = 120){

  // passing, fail and all probes
  vector<TH1F>  hpass;
  vector<TH1F>  hfail;
  vector<TH1F>  hall;
  vector<TH1F>  hp;
  vector<TH1F>  ha;
  vector<TH1F>  hr;
  vector<float> ptMin;

  for(size_t ptbin = 0 ; ptbin <  ptBin.size()-1; ptbin++){
    for(size_t etabin = 0 ; etabin <  etaBin.size()-1; etabin++){
      ptMin.push_back(ptBin.at(ptbin));
      hpass.push_back(TH1F(Form("hpass_pt_%.1f_%.1f_eta_%.2f_%.2f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1)),"",nbins, xmin, xmax));
      hfail.push_back(TH1F(Form("hfail_pt_%.1f_%.1f_eta_%.2f_%.2f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1)),"",nbins, xmin, xmax));
      hall.push_back(TH1F(Form("hall_pt_%.1f_%.1f_eta_%.2f_%.2f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1)),"",nbins, xmin, xmax));
      hp.push_back(TH1F(Form("hp_pt_%.1f_%.1f_eta_%.2f_%.2f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1)),"",1, xmin, xmax));
      ha.push_back(TH1F(Form("ha_pt_%.1f_%.1f_eta_%.2f_%.2f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1)),"",1, xmin, xmax));
      hr.push_back(TH1F(Form("hr_pt_%.1f_%.1f_eta_%.2f_%.2f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1)),"",1, xmin, xmax));
      
      hpass.back().Sumw2();
      hfail.back().Sumw2();
      hall.back().Sumw2();

      hp.back().Sumw2();
      ha.back().Sumw2();
      hr.back().Sumw2();
    }
  }


  // pileup-re-weight from external file                                                                                                                                        
  TFile* pufile = new TFile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt_36.40_summer16.root");
  TH1* puhist   = (TH1*) pufile->Get("puhist");
  
  TTreeReader reader(tree);
  TTreeReaderValue<float> mass (reader, "mass");
  TTreeReaderValue<float> probe_mass  (reader, "probe_mass");
  TTreeReaderValue<float> eta  (reader, "eta");
  TTreeReaderValue<float> phi  (reader, "phi");
  TTreeReaderValue<float> pt   (reader, "pt");
  TTreeReaderValue<float> tag_mass (reader, "tag_mass");
  TTreeReaderValue<float> tag_eta  (reader, "tag_eta");
  TTreeReaderValue<float> tag_phi  (reader, "tag_phi");
  TTreeReaderValue<float> tag_pt   (reader, "tag_pt");
  TTreeReaderValue<float> nvtx (reader, "nvtx");
  TTreeReaderValue<float> wgt  (reader, "wgt");
  TTreeReaderValue<int>   id   (reader, "tightid");  // tight lepton id used offline
  TTreeReaderValue<int>   hlt  (reader,triggerFlag.c_str());
  if(triggerFlagAlt == "")
    triggerFlagAlt = triggerFlag;
  TTreeReaderValue<int>   hltAlt  (reader,triggerFlagAlt.c_str());
  

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

  // start loopin on events
  string currentTree ;
  while(reader.Next()){

    if(currentTree != reader.GetTree()->GetCurrentFile()->GetName()){
      cout<<"current tree "<<reader.GetTree()->GetCurrentFile()->GetName()<<endl;
      currentTree = reader.GetTree()->GetCurrentFile()->GetName();
    }

    ///////// better selections for tag lepton
    if(isSingleMuon and *tag_pt < muonTagPtCut) continue;
    else if(not isSingleMuon and *tag_pt < electronTagPtCut) continue;
    if(isSingleMuon and fabs(*tag_eta) > muonTagEtaCut) continue;
    else if(not isSingleMuon and fabs(*tag_eta) > electronTagEtaCut) continue;
    
    // ask for a tight id on the probe leg
    if(*id <= 0) continue;
    // ask for dR(tag,probe)
    float dphi = fabs(*phi-*tag_phi);
    if(dphi > TMath::Pi())
      dphi = 2*TMath::Pi()-dphi;
    float dR = sqrt(dphi*dphi+fabs(*eta-*tag_eta)*fabs(*eta-*tag_eta));
    if(dR < 0.3) continue;
    
    if(isSingleMuon and (TString(currentTree).Contains("Run2016B") or TString(currentTree).Contains("Run2016C") or TString(currentTree).Contains("Run2016D") or
			 TString(currentTree).Contains("Run2016E") or TString(currentTree).Contains("Run2016F"))){
      if(((*eta > 1.2 and *tag_eta > 1.2) or (*eta < -1.2 and *tag_eta < -1.2)) and dphi < 1.2) continue;
    }
    
    // find the right bin                                                                                                                                                                          
    size_t ptbin = 0;
    size_t etabin = 0;
    for( ; ptbin <  ptBin.size()-1; ptbin++){
      if( *pt <= ptBin.at(ptbin) or *pt > ptBin.at(ptbin+1)) continue;
      else break;
    }
    for( ; etabin <  etaBin.size()-1; etabin++){
      if(isabseta and ( fabs(*eta) <= etaBin.at(etabin) or fabs(*eta) > etaBin.at(etabin+1))) continue;
      else if(not isabseta and (*eta <= etaBin.at(etabin) or *eta > etaBin.at(etabin+1))) continue;
      else break;
    }
    
    // remove boundaries                                                                                                                                                                       
    if(*pt >= ptBin.at(ptBin.size()-1))
      continue;
    if(*pt < ptBin.at(0))
      continue;
    if(isabseta and fabs(*eta) >= etaBin.at(etaBin.size()-1))
      continue;
    else if(not isabseta and *eta >= etaBin.at(etaBin.size()-1))
      continue;
    else if(not isabseta and *eta <= etaBin.at(0))
      continue;

    TString bin = Form("pt_%.1f_%.1f_eta_%.2f_%.2f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1));
    size_t binHisto = 0;
    for( ; binHisto < hpass.size()-1; binHisto++){
      if(TString(hpass.at(binHisto).GetName()).Contains(bin))
        break;
    }

    // truth matching in mc
    if(isMC and not *mcTrue) continue;

    // pu weight
    Float_t puwgt = 1.;
    if (*nvtx <= 60 and isMC) puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));    
    
    // failing trigger
    if((*hlt == 0 and *hltAlt == 0) and isMC)      
      hfail.at(binHisto).Fill(*mass, puwgt*(*wgt)/wgtsum);
    else if((*hlt == 0 and *hltAlt == 0) and not isMC)
      hfail.at(binHisto).Fill(*mass);    
    // passing trigger
    else if((*hlt == 1 or *hltAlt == 1) and isMC){ 
      hpass.at(binHisto).Fill(*mass, puwgt*(*wgt)/wgtsum);
      hp.at(binHisto).Fill(*mass, puwgt*(*wgt)/wgtsum);
    }
    else if((*hlt == 1 or *hltAlt == 1) and not isMC){      
      hpass.at(binHisto).Fill(*mass);
      hp.at(binHisto).Fill(*mass);
    }
    
    // All events (passing + failing trigger)
    if(isMC){
      ha.at(binHisto).Fill(*mass, puwgt*(*wgt)/wgtsum);
      hall.at(binHisto).Fill(*mass, puwgt*(*wgt)/wgtsum);
    }
    else{
      ha.at(binHisto).Fill(*mass);
      hall.at(binHisto).Fill(*mass);
    }
  }
  
  // Calculate efficiency making ratios
  for(size_t isize = 0; isize < hp.size(); isize++){
    // useful just for MC negative weights
    for (int j = 1; j <= nbins; j++) {
      if (hpass.at(isize).GetBinContent(j) < 0.) hpass.at(isize).SetBinContent(j, 0.0);
      if (hfail.at(isize).GetBinContent(j) < 0.) hfail.at(isize).SetBinContent(j, 0.0);
      if (hall.at(isize).GetBinContent(j)  < 0.) hall.at(isize).SetBinContent(j, 0.0);
    }
    if (hp.at(isize).GetBinContent(1)  < 0.) hp.at(isize).SetBinContent(1, 0.0);
    if (ha.at(isize).GetBinContent(1)  < 0.) ha.at(isize).SetBinContent(1, 0.0);

    hr.at(isize).Divide(&hp.at(isize),&ha.at(isize),1,1,"B");    
  }

  /// Fill summary histogram
  for(size_t etabin = 1 ; etabin <  histoNum->GetNbinsX()+1; etabin++){
    for(size_t ptbin = 1 ; ptbin < histoNum->GetNbinsY()+1; ptbin++){

      TString name  = Form("pt_%.1f_%.1f_eta_%.2f_%.2f",histoNum->GetYaxis()->GetBinLowEdge(ptbin),histoNum->GetYaxis()->GetBinLowEdge(ptbin+1),
			   histoNum->GetXaxis()->GetBinLowEdge(etabin),histoNum->GetXaxis()->GetBinLowEdge(etabin+1));
      size_t binHisto = 0;
      for( ; binHisto < hp.size()-1; binHisto++){
	if(TString(hp.at(binHisto).GetName()).Contains(name))
	  break;
      }
      
      histoNum->SetBinContent(etabin,ptbin,hp.at(binHisto).GetBinContent(1));
      histoNum->SetBinError(etabin,ptbin,hp.at(binHisto).GetBinError(1));
      histoDen->SetBinContent(etabin,ptbin,ha.at(binHisto).GetBinContent(1));
      histoDen->SetBinError(etabin,ptbin,ha.at(binHisto).GetBinError(1));

      std::cout << "Efficiency for "+triggerFlag+" -- pT [" << histoNum->GetYaxis()->GetBinLowEdge(ptbin) << ", " << histoNum->GetYaxis()->GetBinLowEdge(ptbin+1) << "], eta [" << histoNum->GetXaxis()->GetBinLowEdge(etabin) << ", " << histoNum->GetXaxis()->GetBinLowEdge(etabin+1) << "]  :  " << hr.at(binHisto).GetBinContent(1) << " +/- " << hr.at(binHisto).GetBinError(1) << std::endl; 
      
    }  
  }


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
    if(isMuon)
      projection_pt->GetYaxis()->SetRangeUser(0.5,1.05);
    else
      projection_pt->GetYaxis()->SetRangeUser(0.5,1.05);

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
  bool isabseta = false;

  if(isSingleMuon){
    binningPt      = {20.,22.,24.,26.,28.,30,33.,36.,40,50,60,70,85,100,125,150,175,200,250,300};
    binningHighPt  = {90,100,105,110,115,120,125,135,150,175,200};
    binningEta     = {-2.4,-2.1,-1.2,-0.9,-0.45,0.,0.45,0.9,1.2,2.1,2.4};
  }
  else{
    binningPt      = {25,30,40,45,50,60,70,85,115,130,150,175,200,250,300};
    binningHighPt  = {90,100,105,110,115,120,125,150,175,200,250};
    binningEta     = {-2.5,-2.1,-1.56,-1.44,-1,-0.5,0,0.5,1,1.444,1.56,2.1,2.5};
  }
  
  
  //prepare the histograms for muons
  TH2F* Passing_muIso24   = new TH2F("Passing_muIso24","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_mu50      = new TH2F("Passing_mu50","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_mu        = new TH2F("Passing_mu","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);

  Passing_muIso24->Sumw2();
  Passing_mu50->Sumw2();
  Passing_mu->Sumw2();

  TH2F* Total_muIso24   = new TH2F("Total_muIso24","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_mu50      = new TH2F("Total_mu50","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_mu        = new TH2F("Total_mu","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);

  Total_muIso24->Sumw2();
  Total_mu50->Sumw2();
  Total_mu->Sumw2();

  //prepare the histograms for electrons
  TH2F* Passing_ele242p1wploose = new TH2F("Passing_ele242p1wploose","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_ele252p1wptight = new TH2F("Passing_ele252p1wptight","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_ele272p1wploose = new TH2F("Passing_ele272p1wploose","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_ele272p1wptight = new TH2F("Passing_ele272p1wptight","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_ele27wptight    = new TH2F("Passing_ele27wptight","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Passing_ele105          = new TH2F("Passing_ele105","",binningEta.size()-1,&binningEta[0],binningHighPt.size()-1,&binningHighPt[0]);
  TH2F* Passing_ele115          = new TH2F("Passing_ele115","",binningEta.size()-1,&binningEta[0],binningHighPt.size()-1,&binningHighPt[0]);
  TH2F* Passing_ele             = new TH2F("Passing_ele","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);

  Passing_ele242p1wploose->Sumw2();
  Passing_ele252p1wptight->Sumw2();
  Passing_ele272p1wploose->Sumw2();
  Passing_ele272p1wptight->Sumw2();
  Passing_ele27wptight->Sumw2();
  Passing_ele105->Sumw2();
  Passing_ele115->Sumw2();
  Passing_ele->Sumw2();

  TH2F* Total_ele242p1wploose = new TH2F("Total_ele242p1wploose","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_ele252p1wptight = new TH2F("Total_ele252p1wptight","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_ele272p1wploose = new TH2F("Total_ele272p1wploose","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_ele272p1wptight = new TH2F("Total_ele272p1wptight","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_ele27wptight    = new TH2F("Total_ele27wptight","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);
  TH2F* Total_ele105          = new TH2F("Total_ele105","",binningEta.size()-1,&binningEta[0],binningHighPt.size()-1,&binningHighPt[0]);
  TH2F* Total_ele115          = new TH2F("Total_ele115","",binningEta.size()-1,&binningEta[0],binningHighPt.size()-1,&binningHighPt[0]);
  TH2F* Total_ele             = new TH2F("Total_ele","",binningEta.size()-1,&binningEta[0],binningPt.size()-1,&binningPt[0]);

  Total_ele242p1wploose->Sumw2();
  Total_ele252p1wptight->Sumw2();
  Total_ele272p1wploose->Sumw2();
  Total_ele272p1wptight->Sumw2();
  Total_ele27wptight->Sumw2();
  Total_ele105->Sumw2();
  Total_ele115->Sumw2();
  Total_ele->Sumw2();
  
  // start looping on the pt and eta bins
  if(isSingleMuon){
    makeTrigger(tree,Passing_muIso24,Total_muIso24,isMC,isSingleMuon,outputDIR,"hltmu24",binningPt,binningEta,isabseta,"hlttkmu24"); 
    makeTrigger(tree,Passing_mu50,Total_mu50,isMC,isSingleMuon,outputDIR,"hltmu50",binningPt,binningEta,isabseta,"hlttkmu50"); 
    makeTrigger(tree,Passing_mu,Total_mu,isMC,isSingleMuon,outputDIR,"hltmu",binningPt,binningEta,isabseta,"hlttkmu"); 
  }
  else{
    makeTrigger(tree,Passing_ele242p1wploose,Total_ele242p1wploose,isMC,isSingleMuon,outputDIR,"hltele24eta2p1wpl",binningPt,binningEta,isabseta);
    makeTrigger(tree,Passing_ele252p1wptight,Total_ele252p1wptight,isMC,isSingleMuon,outputDIR,"hltele25eta2p1wpt",binningPt,binningEta,isabseta);
    makeTrigger(tree,Passing_ele272p1wploose,Total_ele272p1wploose,isMC,isSingleMuon,outputDIR,"hltele27eta2p1wpl",binningPt,binningEta,isabseta);
    makeTrigger(tree,Passing_ele272p1wptight,Total_ele272p1wptight,isMC,isSingleMuon,outputDIR,"hltele27eta2p1wpt",binningPt,binningEta,isabseta);
    makeTrigger(tree,Passing_ele27wptight,Total_ele27wptight,isMC,isSingleMuon,outputDIR,"hltele27wpt",binningPt,binningEta,isabseta);
    makeTrigger(tree,Passing_ele105,Total_ele105,isMC,isSingleMuon,outputDIR,"hltele105",binningPt,binningEta,isabseta);
    makeTrigger(tree,Passing_ele115,Total_ele115,isMC,isSingleMuon,outputDIR,"hltele115",binningPt,binningEta,isabseta);
    makeTrigger(tree,Passing_ele,Total_ele,isMC,isSingleMuon,outputDIR,"hltele",binningPt,binningEta,isabseta);
  }

  
  TEfficiency* trgeff_muIso24 = NULL;
  TEfficiency* trgeff_mu50    = NULL;
  TEfficiency* trgeff_mu      = NULL;

  TEfficiency* trgeff_ele242p1wploose = NULL;
  TEfficiency* trgeff_ele252p1wptight = NULL;
  TEfficiency* trgeff_ele272p1wploose = NULL;
  TEfficiency* trgeff_ele272p1wptight = NULL;
  TEfficiency* trgeff_ele27wptight    = NULL;
  TEfficiency* trgeff_ele105          = NULL;
  TEfficiency* trgeff_ele115          = NULL;
  TEfficiency* trgeff_ele             = NULL;

  if(isSingleMuon){
    trgeff_muIso24   = new TEfficiency(*Passing_muIso24,*Total_muIso24);
    trgeff_muIso24->SetName("trgeff_mu24");
    trgeff_mu50   = new TEfficiency(*Passing_mu50,*Total_mu50);
    trgeff_mu50->SetName("trgeff_mu50");
    trgeff_mu   = new TEfficiency(*Passing_mu,*Total_mu);
    trgeff_mu->SetName("trgeff_mu");
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
    trgeff_ele115       =  new TEfficiency(*Passing_ele115,*Total_ele115);
    trgeff_ele115->SetName("trgeff_ele115");
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
    plotTriggerEfficiency(canvas,trgeff_muIso24,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_mu50,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_mu,isSingleMuon,outputDIR,false,lumi,doFit);
  }
  else{ 
    plotTriggerEfficiency(canvas,trgeff_ele242p1wploose,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele252p1wptight,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele272p1wploose,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele272p1wptight,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele105,isSingleMuon,outputDIR,true,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele115,isSingleMuon,outputDIR,true,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele27wptight,isSingleMuon,outputDIR,false,lumi,doFit);
    plotTriggerEfficiency(canvas,trgeff_ele,isSingleMuon,outputDIR,false,lumi,doFit);
  }

  
  TFile* triggerEfficiency = new TFile((outputDIR+"/"+name+".root").c_str(),"RECREATE");
  triggerEfficiency->cd();
  if(isSingleMuon){
    if(trgeff_muIso24)
      trgeff_muIso24->Write();
    if(trgeff_mu50)
      trgeff_mu50->Write();
    if(trgeff_mu)
      trgeff_mu->Write();
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
    if(trgeff_ele115)
      trgeff_ele115->Write();
    if(trgeff_ele)
      trgeff_ele->Write();
  }
  triggerEfficiency->Close();
}

