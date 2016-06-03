#include "../CMS_lumi.h"

/// make the templates
void maketemplate(const string & inputDIR, 
		  const bool   & isMuon, 
		  const string & outputDIR, 
		  const string & idName,
		  const float  & ptMin,
		  const float  & ptMax,
		  const float  & etaMin,
		  const float  & etaMax,
		  bool  smooth = true,
		  int   nbins = 60, 
		  float xmin  = 65, 
		  float xmax  = 115) {

  TChain* inputMC = NULL;
  if(isMuon)
    inputMC = new TChain("muontnptree/fitter_tree");
  else
    inputMC = new TChain("electrontnptree/fitter_tree");

  inputMC->Add((inputDIR+"/*root").c_str());
  
  // pileup-re-weight from external file
  TFile* pufile = new TFile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt.root");
  TH1* puhist = (TH1*)pufile->Get("puhist");

  TTreeReader reader(inputMC);
  TTreeReaderValue<float> mass (reader, "mass");
  TTreeReaderValue<float> eta  (reader, "eta");
  TTreeReaderValue<float> phi  (reader, "phi");
  TTreeReaderValue<float> pt   (reader, "pt");
  TTreeReaderValue<float> nvtx (reader, "nvtx");
  TTreeReaderValue<float> wgt  (reader, "wgt");
  TTreeReaderValue<int> id   (reader, idName.c_str());
  TTreeReaderValue<int> mcTrue  (reader, "mcTrue");
  TTreeReaderValue<float>    mcMass  (reader, "mcMass");

  TH1F hpass("hpass", "", nbins, xmin, xmax);
  TH1F hfail("hfail", "", nbins, xmin, xmax);
  TH1F hp("hp" , "", 1, xmin, xmax);
  TH1F ha("ha" , "", 1, xmin, xmax);
  TH1F hr("hr" , "", 1, xmin, xmax);
  hpass.Sumw2();
  hfail.Sumw2();
  hp.Sumw2();
  ha.Sumw2();
  hr.Sumw2();
  
  // weight sum
  double wgtsum = 0;
  while(reader.Next()){
    wgtsum += *wgt;
  }

  reader.SetEntry(0);
  // loop on the event
  while(reader.Next()){
    
    Float_t puwgt = 0.;
    if (*nvtx <= 40) puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));

    // check the probe lepton information --> if it falls inside the bins
    if(*pt < ptMin or *pt > ptMax) continue;
    if(fabs(*eta) < etaMin or fabs(*eta) > etaMax) continue;
    if(not isMuon and fabs(*eta) > 1.442 and fabs(*eta) < 1.566) continue;
    // if not matched to a genLepton skip --> we want to extract the true templateds
    if(not *mcTrue) continue;
    if (*id == 0) hfail.Fill(*mass, puwgt*(*wgt)/wgtsum);
    if (*id >  0) hpass.Fill(*mass, puwgt*(*wgt)/wgtsum);
    if (*id >  0) hp.Fill(*mass, puwgt*(*wgt)/wgtsum);
    ha.Fill(*mass, puwgt*(*wgt)/wgtsum);
  }

  for (int j = 1; j <= nbins; j++) {
    if (hpass.GetBinContent(j) < 0.) hpass.SetBinContent(j, 0.0);
    if (hfail.GetBinContent(j) < 0.) hfail.SetBinContent(j, 0.0);
  }
  
  // efficiency value with Binomial uncertainty
  hr.Divide(&hp, &ha, 1.0, 1.0, "B");
  // smooth templates two times
  if(smooth){// amcnlo as lower stat
    if(ptMin >= 50){
      hpass.Rebin(2);
      hfail.Rebin(2);
    }
  }
  
  // Build the RooHistPdf  
  string fileName = "";
  if(isMuon)
    fileName = "template_TaP_muon_"+idName+"_pt_"+string(Form("%.1f",ptMin))+"_"+string(Form("%.1f",ptMax))+"_eta_"+string(Form("%.1f",etaMin))+"_"+string(Form("%.1f",etaMax));
  else
    fileName = "template_TaP_electron_"+idName+"_pt_"+string(Form("%.1f",ptMin))+"_"+string(Form("%.1f",ptMax))+"_eta_"+string(Form("%.1f",etaMin))+"_"+string(Form("%.1f",etaMax));

  TFile outfile((outputDIR+"/"+fileName+".root").c_str(), "RECREATE");
  hr.Write("efficiency");

  RooRealVar m("mass", "", xmin, xmax);
  RooDataHist dhpass("dhpass", "", RooArgList(m), RooFit::Import(hpass), 0);
  RooDataHist dhfail("dhfail", "", RooArgList(m), RooFit::Import(hfail), 0);
  RooHistPdf signalPassMC("signalPassMC", "", RooArgSet(m), dhpass);
  RooHistPdf signalFailMC("signalFailMC", "", RooArgSet(m), dhfail);

  RooWorkspace w("w", "");
  w.import(dhpass);
  w.import(dhfail);
  w.import(signalPassMC);
  w.import(signalFailMC);
  w.Write();
  
  cout << "Efficiency -- pT [" << ptMin << ", " << ptMax << "], eta [" << etaMin << ", " << etaMax << "]   :   " << hr.GetBinContent(1) << " +/- " << hr.GetBinError(1) << endl;

  TCanvas* c1 = new TCanvas("c1","",600,650);
  c1->cd();
  TH1F* hpassForPlot = (TH1F*) hpass.Clone("hpassForPlot");
  hpassForPlot->Scale(1./hpassForPlot->Integral());
  TH1F* hfailForPlot = (TH1F*) hfail.Clone("hfailForPlot");
  hfailForPlot->Scale(1./hfailForPlot->Integral());
  hfailForPlot->SetLineColor(kRed);
  hpassForPlot->SetLineColor(kBlue);
  hfailForPlot->SetLineWidth(2);
  hpassForPlot->SetLineWidth(2);  
  hpassForPlot->Draw("hist");
  hfailForPlot->Draw("hist same");

  TLegend* leg = new TLegend(0.2,0.7,0.4,0.9);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->AddEntry(hpassForPlot,"Passing probes","L");
  leg->AddEntry(hfailForPlot,"Failing probes","L");
  leg->Draw("same");

  c1->SaveAs((outputDIR+"/"+fileName+".png").c_str(),"png");
  c1->SaveAs((outputDIR+"/"+fileName+".pdf").c_str(),"pdf");
}

void makeTnPTemplates(string inputDIR, bool isMuon, string outputDIR) {

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());

  vector<float> ptBinMuon      = {10.,20.,30.,40.,50.,70.,200};
  vector<float> ptBinElectron  = {10.,20.,30.,40.,50.,70.,200};
  vector<float> etaBinMuon     = {0.,1.2,2.4};
  vector<float> etaBinElectron = {0.,1.5,2.5};

  setTDRStyle();
  
  if(isMuon){
    for(size_t ipt = 0; ipt < ptBinMuon.size()-1; ipt++){
      for(size_t ieta = 0; ieta < etaBinMuon.size()-1; ieta++){
	string idName = "looseid";
	maketemplate(inputDIR, isMuon, outputDIR, idName, ptBinMuon.at(ipt), ptBinMuon.at(ipt+1), etaBinMuon.at(ieta), etaBinMuon.at(ieta+1),true);    
	idName = "tightid";
	maketemplate(inputDIR, isMuon, outputDIR, idName, ptBinMuon.at(ipt), ptBinMuon.at(ipt+1), etaBinMuon.at(ieta), etaBinMuon.at(ieta+1),true);

      }
    }
  }
  else{
    for(size_t ipt = 0; ipt < ptBinElectron.size()-1; ipt++){
      for(size_t ieta = 0; ieta < etaBinElectron.size()-1; ieta++){
	string idName = "vetoid";
	maketemplate(inputDIR, isMuon, outputDIR, idName, ptBinElectron.at(ipt), ptBinElectron.at(ipt+1), etaBinElectron.at(ieta), etaBinElectron.at(ieta+1));    	
	idName = "tightid";
	maketemplate(inputDIR, isMuon, outputDIR, idName, ptBinElectron.at(ipt), ptBinElectron.at(ipt+1), etaBinElectron.at(ieta), etaBinElectron.at(ieta+1));
      }
    }
  }
}
