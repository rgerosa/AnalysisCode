#include <fstream>
#include "TROOT.h"
#include "../CMS_lumi.h"
#include "triggerUtils.h"


void makeTrigger(TTree* tree, // tree 
		 TH2F*  trigeff,
		 const bool & isMC,
		 const bool & isSingleMuon,
		 const string & outputDIR,
		 const string & triggerFlag,
		 const float & ptMin,
		 const float & ptMax,
		 const float & etaMin,
		 const float & etaMax,
		 int   nbins = 60,
		 float xmin  = 65,
		 float xmax  = 115){

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
  TFile* pufile = new TFile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt.root");
  TH1* puhist = (TH1*)pufile->Get("puhist");

  TTreeReader reader(tree);
  TTreeReaderValue<float> mass (reader, "mass");
  TTreeReaderValue<float> eta  (reader, "eta");
  TTreeReaderValue<float> phi  (reader, "phi");
  TTreeReaderValue<float> pt   (reader, "pt");
  TTreeReaderValue<float> nvtx (reader, "nvtx");
  TTreeReaderValue<float> wgt  (reader, "wgt");
  TTreeReaderValue<int>   id   (reader, "tightid");  // single electorn and muon trigger efficiency are useful in the single and double lepton CR were we ask at least one tight lepton
  TTreeReaderValue<int>   hlt (reader,triggerFlag.c_str());
  TTreeReaderValue<int>   hlt22 (reader,"hltmu22");
  
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
    if(not isSingleMuon and fabs(*eta) > 1.442 and fabs(*eta) < 1.566) continue;
    // if not matched to a genLepton skip --> we want to extract the true templateds                                                                                            
    if(isMC and not *mcTrue) continue;
    if(*id <= 0) continue;

    if(*hlt == 0 and isMC)      
      hfail.Fill(*mass, puwgt*(*wgt)/wgtsum);
    else if(*hlt == 0 and not isMC)
      hfail.Fill(*mass);
    else if(*hlt == 1 and isMC){ 
      hpass.Fill(*mass, puwgt*(*wgt)/wgtsum);
      hp.Fill(*mass, puwgt*(*wgt)/wgtsum);
    }
    else if(*hlt == 1 and not isMC and *hlt22 == 0){      
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

  // useful just for MC negative weights
  for (int j = 1; j <= nbins; j++) {
    if (hpass.GetBinContent(j) < 0.) hpass.SetBinContent(j, 0.0);
    if (hfail.GetBinContent(j) < 0.) hfail.SetBinContent(j, 0.0);
    if (hp.GetBinContent(j) < 0.) hpass.SetBinContent(j, 0.0);
    if (ha.GetBinContent(j) < 0.) hfail.SetBinContent(j, 0.0);
    if (hall.GetBinContent(j)  < 0.) hall.SetBinContent(j, 0.0);
  }
  // identical to TGraphAsymErrors
  hr.Divide(&hp, &ha, 1.0, 1.0, "B");
  trigeff->SetBinContent(trigeff->GetXaxis()->FindBin((ptMax+ptMin)/2),trigeff->GetYaxis()->FindBin((etaMax+etaMin)/2),hr.GetBinContent(1));
  trigeff->SetBinError(trigeff->GetXaxis()->FindBin((ptMax+ptMin)/2),trigeff->GetYaxis()->FindBin((etaMax+etaMin)/2),hr.GetBinError(1));

  std::cout << "Efficiency for "+triggerFlag+" -- pT [" << ptMin << ", " << ptMax << "], eta [" << etaMin << ", " << etaMax << "]  :  " << hr.GetBinContent(1) << " +/- " << hr.GetBinError(1) << std::endl; 
}

void plotTriggerEfficiency(TCanvas* canvas,
			   TH2F*    histo,
			   bool     isMuon,
			   string   outputDIR,
			   float    lumi   = 0.59,
			   bool     isLogz = false){


  gPad->SetBottomMargin(0.10);
  gPad->SetRightMargin(0.20);

  canvas->cd();
  if(isMuon){
    histo->GetXaxis()->SetTitle("Muon p_{T} [GeV]");
    histo->GetYaxis()->SetTitle("Muon |#eta|");
  }
  else{
    histo->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
    histo->GetYaxis()->SetTitle("Electron |#eta|");
  }
  histo->GetZaxis()->SetTitle("Trigger Efficiency");
  histo->Draw("colz");
  CMS_lumi(canvas,string(Form("%.2f",lumi)),true);
  if(isLogz)
    canvas->SetLogz();
  canvas->SaveAs((outputDIR+"/"+string(histo->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+string(histo->GetName())+".pdf").c_str(),"pdf");

  gPad->SetRightMargin(0.06);

  canvas->SetLogz(0);
  // project by hand
  TH1F projection_pt ("projection_pt","",histo->GetXaxis()->GetXbins()->GetSize()-1,histo->GetXaxis()->GetXbins()->GetArray());
  TF1  fitfunc ("fitfunc", ErfCB, 0, 500, 5);
  if(not isMuon)
    fitfunc.SetParameters(22., 5., 10., 2., 1.);
  else
    fitfunc.SetParameters(16., 5., 10., 2., 1.);
  
  for(int iBinY = 0; iBinY < histo->GetNbinsY(); iBinY++){
    for(int iBinX = 0; iBinX < histo->GetNbinsX(); iBinX++){      
      projection_pt.SetBinContent(iBinX+1,histo->GetBinContent(iBinX+1,iBinY+1));
      projection_pt.SetBinError(iBinX+1,histo->GetBinError(iBinX+1,iBinY+1));
    }
    projection_pt.SetMarkerSize(1);
    projection_pt.SetMarkerStyle(20);
    projection_pt.SetMarkerColor(kBlack);
    projection_pt.SetLineColor(kBlack);
    projection_pt.Draw("E1P");
    fitfunc.SetLineColor(kBlue);
    fitfunc.SetLineWidth(2);
    projection_pt.Fit(&fitfunc);
    CMS_lumi(canvas,string(Form("%.2f",lumi)),true);
    canvas->SaveAs((outputDIR+"/"+string(histo->GetName())+string(Form("_projection_pt_eta_%d_%d.png",int(projection_pt.GetYaxis()->GetBinLowEdge(iBinY+1)),
								       int(projection_pt.GetYaxis()->GetBinLowEdge(iBinY+2))))).c_str(),"png");
    canvas->SaveAs((outputDIR+"/"+string(histo->GetName())+string(Form("_projection_pt_eta_%d_%d.pdf",int(projection_pt.GetYaxis()->GetBinLowEdge(iBinY+1)),
								       int(projection_pt.GetYaxis()->GetBinLowEdge(iBinY+2))))).c_str(),"pdf");
  }
}
  

void makeTriggerEfficiency(string inputDIR, // where trees are located
			   string outputDIR, // to store plots and files
			   bool   isMC,
			   bool   isSingleMuon){
  
  
  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  
  // take the input tree --> valid for both data and MC
  TChain* tree = NULL;
  if(isSingleMuon)
    tree = new TChain("muontnptree/fitter_tree");
  else
    tree = new TChain("electrontnptree/fitter_tree");

  tree->Add((inputDIR+"/*root").c_str());

  // select binning for efficiency measurement
  vector<float> binningPt;
  vector<float> binningEta;
  if(isSingleMuon){
    binningPt  = {10,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,32.5,35,37.5,40,45,50,60,70,85,100,125,150,200};
    binningEta = {0,1.2,2.4};
  }
  else{
    binningPt  = {10,15,17.5,18.5,20,20.5,21,21.5,22,22.5,23,23.5,24,25,26,27,28,29,30,32.5,35,40,45,50,60,70,85,100,125,150,200};
    binningEta = {0,1.5,2.5};
  }
    
  //prepare the histograms for muons
  TH2F* trigeff_mu20      = new TH2F("trigeff_muIso20","",binningPt.size()-1,&binningPt[0],binningEta.size()-1,&binningEta[0]);
  TH2F* trigeff_muIsoTk20 = new TH2F("trigeff_muIsoTk20","",binningPt.size()-1,&binningPt[0],binningEta.size()-1,&binningEta[0]);
  TH2F* trigeff_mu22      = new TH2F("trigeff_muIso22","",binningPt.size()-1,&binningPt[0],binningEta.size()-1,&binningEta[0]);
  TH2F* trigeff_muIsoTk22 = new TH2F("trigeff_muIsoTk22","",binningPt.size()-1,&binningPt[0],binningEta.size()-1,&binningEta[0]);
  TH2F* trigeff_mu        = new TH2F("trigeff_muIso","",binningPt.size()-1,&binningPt[0],binningEta.size()-1,&binningEta[0]);
  TH2F* trigeff_muIsoTk   = new TH2F("trigeff_muIsoTk","",binningPt.size()-1,&binningPt[0],binningEta.size()-1,&binningEta[0]);

  //prepare the histograms for electrons
  TH2F* trigeff_ele23wploose = new TH2F("trigeff_ele23wploose","",binningPt.size()-1,&binningPt[0],binningEta.size()-1,&binningEta[0]);
  TH2F* trigeff_ele27wploose = new TH2F("trigeff_ele27wploose","",binningPt.size()-1,&binningPt[0],binningEta.size()-1,&binningEta[0]);
  TH2F* trigeff_ele25wptight = new TH2F("trigeff_ele25wptight","",binningPt.size()-1,&binningPt[0],binningEta.size()-1,&binningEta[0]);
  TH2F* trigeff_ele105       = new TH2F("trigeff_ele105","",binningPt.size()-1,&binningPt[0],binningEta.size()-1,&binningEta[0]);
  TH2F* trigeff_ele          = new TH2F("trigeff_ele","",binningPt.size()-1,&binningPt[0],binningEta.size()-1,&binningEta[0]);
  
  // start looping on the pt and eta bins
  for(size_t ipt = 0; ipt < binningPt.size()-1; ipt++){
    for(size_t ieta = 0; ieta < binningEta.size()-1; ieta++){
      if(isSingleMuon){
	makeTrigger(tree,trigeff_mu20,isMC,isSingleMuon,outputDIR,"hltmu20",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1)); // mu20
	//	makeTrigger(tree,trigeff_muIsoTk20,isMC,isSingleMuon,outputDIR,"hlttkmu20",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1)); 
	//	makeTrigger(tree,trigeff_mu22,isMC,isSingleMuon,outputDIR,"hltmu22",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1)); 
	//	makeTrigger(tree,trigeff_muIsoTk22,isMC,isSingleMuon,outputDIR,"hlttkmu22",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1)); 
	//	makeTrigger(tree,trigeff_mu,isMC,isSingleMuon,outputDIR,"hltmu",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1)); 
	//	makeTrigger(tree,trigeff_muIsoTk,isMC,isSingleMuon,outputDIR,"hlttkmu",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1)); 
      }
      else{
	makeTrigger(tree,trigeff_ele23wploose,isMC,isSingleMuon,outputDIR,"hltele23",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1));
	//makeTrigger(tree,trigeff_ele27wploose,isMC,isSingleMuon,outputDIR,"hltele27",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1));
	//	makeTrigger(tree,trigeff_ele25wptight,isMC,isSingleMuon,outputDIR,"hltele25WPTight",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1));
	//	makeTrigger(tree,trigeff_ele105,isMC,isSingleMuon,outputDIR,"hltele105",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1));
	//	makeTrigger(tree,trigeff_ele,isMC,isSingleMuon,outputDIR,"hltele",binningPt.at(ipt),binningPt.at(ipt+1),binningEta.at(ieta),binningEta.at(ieta+1));
      }
    }
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
  TCanvas* canvas = new TCanvas("canvas","",625,600);
  canvas->cd();

  if(isSingleMuon){
    plotTriggerEfficiency(canvas,trigeff_mu20,isSingleMuon,outputDIR);
    //    plotTriggerEfficiency(canvas,trigeff_muIsoTk20,isSingleMuon,outputDIR);
    //    plotTriggerEfficiency(canvas,trigeff_mu22,isSingleMuon,outputDIR);
    //    plotTriggerEfficiency(canvas,trigeff_muIsoTk22,isSingleMuon,outputDIR);
    //    plotTriggerEfficiency(canvas,trigeff_mu,isSingleMuon,outputDIR);
    //    plotTriggerEfficiency(canvas,trigeff_muIsoTk,isSingleMuon,outputDIR);
  }
  else{
    plotTriggerEfficiency(canvas,trigeff_ele23wploose,isSingleMuon,outputDIR);
    //plotTriggerEfficiency(canvas,trigeff_ele27wploose,isSingleMuon,outputDIR);
    //    plotTriggerEfficiency(canvas,trigeff_ele25wptight,isSingleMuon,outputDIR);
    //    plotTriggerEfficiency(canvas,trigeff_ele105,isSingleMuon,outputDIR);
    //    plotTriggerEfficiency(canvas,trigeff_ele,isSingleMuon,outputDIR);
  }

  
  TFile* triggerEfficiency = new TFile((outputDIR+"/"+name+".root").c_str(),"RECREATE");
  triggerEfficiency->cd();
  if(isSingleMuon){
    trigeff_mu20->Write();
    trigeff_muIsoTk20->Write();
    trigeff_mu22->Write();
    trigeff_muIsoTk22->Write();
    trigeff_mu->Write();
    trigeff_muIsoTk->Write();
  }
  else{
    trigeff_ele23wploose->Write();
    trigeff_ele27wploose->Write();
    trigeff_ele25wptight->Write();
    trigeff_ele105->Write();
    trigeff_ele->Write();
  }
  triggerEfficiency->Close();
}

