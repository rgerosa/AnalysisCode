#include "../makeTriggerEfficiency/triggerUtils.h"
#include "TnPBinning.h"
#include "../CMS_lumi.h"

vector<TH1*> projectionMC;
vector<TH1*> projectionMC_Alternative;
vector<TH1*> projectionDATA_RooCMSShape;
vector<TH1*> projectionDATA_Exp;
vector<TH1*> projectionDATA_Alternative;
vector<TH1*> projectionDATA_Analytical;
vector<TH1*> projectionSF_RooCMSShape;
vector<TH1*> projectionSF_Exp;
vector<TH1*> projectionSF_Alternative;
vector<TH1*> projectionSF_Analytical;
map<string,TFile*> tagAndProbeFits_RooCMSShape;
map<string,TFile*> tagAndProbeFits_Exp;
map<string,TFile*> tagAndProbeFits_Alternative;
map<string,TFile*> tagAndProbeFits_Analytical;
bool addTurnOnFits_;

void fillEfficiencyMC(TH2F* efficiency, const string & directory, const string & leptonType, const string & typeID, const int & ipu){

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  TFile* inputFile = NULL;
  vector<float> ptBin;
  vector<float> etaBin;
  vector<float> puBin;
  
  if(leptonType == "muon"){
    ptBin  = ptBinMuon;
    etaBin = etaBinMuon;
    puBin  = nvtxBinMuon;
  }
  else if(leptonType == "electron"){
    ptBin  = ptBinElectron;
    etaBin = etaBinElectron;
    puBin  = nvtxBinElectron;
  }
  else if(leptonType == "photon"){
    ptBin = ptBinPhoton;
    etaBin = etaBinPhoton;
    puBin = nvtxBinPhoton;
  }
  
  for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++){
    for(size_t ieta = 0; ieta < etaBin.size()-1; ieta++){
      system(("ls "+directory+" | grep root | grep "+leptonType+" | grep "+typeID+" | grep _pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.1f",etaBin.at(ieta)))+"_"+string(Form("%.1f",etaBin.at(ieta+1)))+"_pu_"+string(Form("%.1f_%.1f",puBin.at(ipu),puBin.at(ipu+1)))+" > file"+leptonType).c_str());
      ifstream file;
      file.open(("file"+leptonType).c_str());
      int nFiles = 0;
      string name;
      if(file.is_open()){
	while(!file.eof()){
	  string line;
	  getline(file,line);	    
	  if(line == "" or line =="\n") continue;
	  nFiles++;
	  name = line;
	}
      }
      file.close();
      system(("rm file"+leptonType).c_str());
      if(nFiles != 1){
	cerr<<"file not found -->check please "<<endl;
      }
      //	cout<<"Opening: "<<directory+"/"+name<<endl;
      inputFile = TFile::Open((directory+"/"+name).c_str(),"READ");
      if(inputFile->TestBit(TFile::kRecovered)) continue;
      TH1F* eff = (TH1F*) inputFile->Get("efficiency");
      efficiency->SetBinContent(ieta+1,ipt+1,eff->GetBinContent(1));
      efficiency->SetBinError(ieta+1,ipt+1,eff->GetBinError(1));
      }
  }
  return;
}

// Fill efficiency from Data
void fillEfficiencyData(TH2F* efficiency, const string & directory, const string & leptonType, const string & typeID, const int & ipu, const string & postfix = "RooCMSShape"){

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  TFile* inputFile = NULL;
  vector<float> ptBin;
  vector<float> etaBin;
  vector<float> puBin;
  
  if(leptonType == "muon"){
    ptBin  = ptBinMuon;
    etaBin = etaBinMuon;
    puBin  = nvtxBinMuon;
  }
  else if(leptonType == "electron"){
    ptBin = ptBinElectron;
    etaBin = etaBinElectron;
    puBin  = nvtxBinElectron;
  }
  else if(leptonType == "photon"){
    ptBin = ptBinPhoton;
    etaBin = etaBinPhoton;
    puBin  = nvtxBinPhoton;
  }
   
  for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++){
    for(size_t ieta = 0; ieta < etaBin.size()-1; ieta++){
      if(postfix != "Alternative" and postfix != "Analytical")
        system(("ls "+directory+" | grep root | grep "+leptonType+" | grep "+typeID+" | grep _pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.1f",etaBin.at(ieta)))+"_"+string(Form("%.1f",etaBin.at(ieta+1)))+"_pu_"+Form("%.1f_%.1f",puBin.at(ipu),puBin.at(ipu+1))+" | grep "+postfix+" | grep -v Alternative | grep -v Analytical > file"+leptonType).c_str());
      else
        system(("ls "+directory+" | grep root | grep "+leptonType+" | grep "+typeID+" | grep _pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.1f",etaBin.at(ieta)))+"_"+string(Form("%.1f",etaBin.at(ieta+1)))+"_pu_"+Form("%.1f_%.1f",puBin.at(ipu),puBin.at(ipu+1))+" | grep "+postfix+" > file"+leptonType).c_str());
      
      ifstream file;
      file.open(("file"+leptonType).c_str());
      int nFiles = 0;
      string name;
      if(file.is_open()){
	while(!file.eof()){
	  string line;
	  getline(file,line);	    
	  if(line == "" or line =="\n") continue;
	  nFiles++;
	  name = line;
	}
      }
      file.close();
      system(("rm file"+leptonType).c_str());
      if(nFiles != 1){
	cerr<<" find not found -->check please "<<endl;
      }
      //      cout<<"Opening: "<<directory+"/"+name<<endl;
      inputFile = TFile::Open((directory+"/"+name).c_str(),"READ");
      if(inputFile->TestBit(TFile::kRecovered)) continue;
      if(postfix == "RooCMSShape" and not TString(inputFile->GetName()).Contains("Alternative"))
	tagAndProbeFits_RooCMSShape["pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.1f",etaBin.at(ieta)))+"_"+string(Form("%.1f",etaBin.at(ieta+1)))+"_pu_"+string(Form("%.1f",puBin.at(ipu)))+"_"+string(Form("%.1f",puBin.at(ipu+1)))] = inputFile;
      else if(postfix == "Exp" and not TString(inputFile->GetName()).Contains("Alternative"))
	tagAndProbeFits_Exp["pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.1f",etaBin.at(ieta)))+"_"+string(Form("%.1f",etaBin.at(ieta+1)))+"_pu_"+string(Form("%.1f",puBin.at(ipu)))+"_"+string(Form("%.1f",puBin.at(ipu+1)))] = inputFile;
	else if(postfix == "Alternative")
	  tagAndProbeFits_Alternative["pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.1f",etaBin.at(ieta)))+"_"+string(Form("%.1f",etaBin.at(ieta+1)))+"_pu_"+string(Form("%.1f",puBin.at(ipu)))+"_"+string(Form("%.1f",puBin.at(ipu+1)))] = inputFile;
	else if(postfix == "Analytical")
	  tagAndProbeFits_Analytical["pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.1f",etaBin.at(ieta)))+"_"+string(Form("%.1f",etaBin.at(ieta+1)))+"_pu_"+string(Form("%.1f",puBin.at(ipu)))+"_"+string(Form("%.1f",puBin.at(ipu+1)))] = inputFile;
      
      RooWorkspace* workspace = (RooWorkspace*) inputFile->FindObjectAny("w");      
      RooFitResult* fitResult = (RooFitResult*) workspace->obj("fitresults");
      if(fitResult){
	RooArgList parameter    = fitResult->floatParsFinal();
	RooRealVar* eff         = dynamic_cast<RooRealVar*>(parameter.find("efficiency"));	
	efficiency->SetBinContent(ieta+1,ipt+1,eff->getVal());
	efficiency->SetBinError(ieta+1,ipt+1,(eff->getErrorHi()+fabs(eff->getErrorLo()))/2);
      }
    }
  }
  return;
}

// plot efficiency
void plotEfficiency(TCanvas* canvas, TH2* histoEfficiency, const string & outputDIR, const string & leptonType, const float & lumi, const bool & isMC, const bool & isScaleFactor, const int & ipu){

  vector<float> puBin;
  
  if(leptonType == "muon")
    puBin  = nvtxBinMuon;
  else if(leptonType == "electron")
    puBin  = nvtxBinElectron;
  else if(leptonType == "photon")
    puBin  = nvtxBinPhoton;

  canvas->Clear();
  canvas->SetRightMargin(0.20);
  canvas->SetLeftMargin(0.12);
  if(not isMC and not isScaleFactor){
    for(int iBinX = 0; iBinX < histoEfficiency->GetNbinsX(); iBinX++){
      for(int iBinY = 0; iBinY < histoEfficiency->GetNbinsY(); iBinY++){
	if(histoEfficiency->GetBinContent(iBinX+1,iBinY+1)/projectionMC.at(iBinX)->GetBinContent(iBinY+1) < 0.5 or histoEfficiency->GetBinContent(iBinX+1,iBinY+1)/projectionMC.at(iBinX)->GetBinContent(iBinY+1) > 1.5)
	  histoEfficiency->SetBinContent(iBinX+1,iBinY+1,histoEfficiency->GetBinContent(iBinX+1,iBinY+2)*0.9);
      }
    }
  }
  
  if(leptonType == "muon"){
    histoEfficiency->GetYaxis()->SetTitle("Muon p_{T} (GeV)");
    histoEfficiency->GetYaxis()->SetTitleOffset(1.1);
    histoEfficiency->GetXaxis()->SetTitle("Muon |#eta|");
    histoEfficiency->GetXaxis()->SetTitleOffset(1.1);
  }
  else if(leptonType == "electron"){
    histoEfficiency->GetYaxis()->SetTitle("Electron p_{T} (GeV)");
    histoEfficiency->GetYaxis()->SetTitleOffset(1.1);
    histoEfficiency->GetXaxis()->SetTitle("Electron |#eta|");
    histoEfficiency->GetXaxis()->SetTitleOffset(1.1);
  }
  else if(leptonType == "photon"){
    histoEfficiency->GetYaxis()->SetTitle("Photon p_{T} (GeV)");
    histoEfficiency->GetYaxis()->SetTitleOffset(1.1);
    histoEfficiency->GetXaxis()->SetTitle("Photon |#eta|");
    histoEfficiency->GetXaxis()->SetTitleOffset(1.1);
  }
  if(not isScaleFactor and leptonType != "photon")
    histoEfficiency->GetZaxis()->SetTitle("Efficiency Lepton ID");
  else if(not isScaleFactor and leptonType == "photon")
    histoEfficiency->GetZaxis()->SetTitle("Efficiency Photon ID");
  else
    histoEfficiency->GetZaxis()->SetTitle("ScaleFactor");

  histoEfficiency->GetZaxis()->SetTitleOffset(1.25);
  gStyle->SetPaintTextFormat(".3f");
  histoEfficiency->Draw("colz text");
  canvas->SetLogy();
  CMS_lumi(canvas,string(Form("%.2f",lumi)),true);
  //  canvas->SaveAs((outputDIR+"/"+string(histoEfficiency->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+string(histoEfficiency->GetName())+".pdf").c_str(),"pdf");
  canvas->SetLogy(0);
  // project by hand as a function of pt                                                                                                                                       
  canvas->SetRightMargin(0.06);
  canvas->SetLeftMargin(0.15);
  string postfix;
  if(isMC)
    postfix = "MC";
  else
    postfix = "Data";

  for(int iBinX = 0; iBinX < histoEfficiency->GetNbinsX(); iBinX++){
    TH1F* projection_pt = new TH1F((string(histoEfficiency->GetName())+"pt_projection_eta_"+string(Form("%.1f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+1)))+"_"+string(Form("%.1f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))).c_str(),"",histoEfficiency->GetYaxis()->GetXbins()->GetSize()-1,histoEfficiency->GetYaxis()->GetXbins()->GetArray());
    
    for(int iBinY = 0; iBinY < histoEfficiency->GetNbinsY(); iBinY++){
      projection_pt->SetBinContent(iBinY+1,histoEfficiency->GetBinContent(iBinX+1,iBinY+1));
      projection_pt->SetBinError(iBinY+1,histoEfficiency->GetBinError(iBinX+1,iBinY+1));
    }

    // style options
    if(leptonType == "electron")
      projection_pt->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
    else if(leptonType == "muon")
      projection_pt->GetXaxis()->SetTitle("Muon p_{T} [GeV]");
    else if(leptonType == "photon")
      projection_pt->GetXaxis()->SetTitle("Photon p_{T} [GeV]");

    if(not isScaleFactor and leptonType != "photon")
      projection_pt->GetYaxis()->SetTitle("Efficiency Lepton ID");
    else if(not isScaleFactor and leptonType == "photon")
      projection_pt->GetYaxis()->SetTitle("Efficiency Photon ID");
    else
      projection_pt->GetYaxis()->SetTitle("Scale Factor");

    projection_pt->GetYaxis()->SetTitleOffset(1.25);

    projection_pt->SetMarkerSize(1);
    projection_pt->SetMarkerStyle(20);
    projection_pt->SetMarkerColor(kBlack);
    projection_pt->SetLineColor(kBlack);
    if(not isScaleFactor)
      projection_pt->GetYaxis()->SetRangeUser(projection_pt->GetMinimum()*0.9,1.05);
    else 
      projection_pt->GetYaxis()->SetRangeUser(projection_pt->GetMinimum()*0.95,1.03);

    projection_pt->Draw("E1P");

    if(not isScaleFactor){
      // add fits
      if(addTurnOnFits_){
	TF1*  fitfunc = new TF1(("fitfunc_"+to_string(iBinX)).c_str(), ErfCB, 0, 200, 5);
	if(leptonType == "electron" or leptonType == "photon")
	  fitfunc->SetParameters(22., 5., 10., 2., 1.);
	else
	  fitfunc->SetParameters(16., 5., 10., 2., 1.);        
	fitfunc->SetLineColor(kBlue);
	fitfunc->SetLineWidth(2);
	projection_pt->Fit(fitfunc);
      }
    }    
    // store canvas
    CMS_lumi(canvas,string(Form("%.2f",lumi)),true);
    //    canvas->SaveAs((outputDIR+"/"+string(histoEfficiency->GetName())+"_vs_pt_eta_"+string(Form("%.1f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+1)))+"_"+string(Form("%.1f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))+".png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/"+string(histoEfficiency->GetName())+"_vs_pt_eta_"+string(Form("%.1f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+1)))+"_"+string(Form("%.1f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))+".pdf").c_str(),"pdf");
    // store projections in eta bins
    if(not isScaleFactor){
      if(isMC){
	if(TString(histoEfficiency->GetName()).Contains("Alternative"))
	  projectionMC_Alternative.push_back(projection_pt);
	else
	  projectionMC.push_back(projection_pt);
      }
      else{
	if(TString(histoEfficiency->GetName()).Contains("RooCMSShape") and not TString(histoEfficiency->GetName()).Contains("Alternative") and not TString(histoEfficiency->GetName()).Contains("Analytical"))
	  projectionDATA_RooCMSShape.push_back(projection_pt);
	else if(TString(histoEfficiency->GetName()).Contains("Exp"))
	  projectionDATA_Exp.push_back(projection_pt);
	else if(TString(histoEfficiency->GetName()).Contains("Alternative"))
	  projectionDATA_Alternative.push_back(projection_pt);
	else if(TString(histoEfficiency->GetName()).Contains("Analytical"))
	  projectionDATA_Analytical.push_back(projection_pt);	
      }
    }

    else{
      if(TString(histoEfficiency->GetName()).Contains("Alternative")){
	projectionSF_Alternative.push_back(projection_pt);
	continue;
      }
      if(TString(histoEfficiency->GetName()).Contains("Analytical")){
	projectionSF_Analytical.push_back(projection_pt);
	continue;
      }
      if(TString(histoEfficiency->GetName()).Contains("RooCMSShape")){
	projectionSF_RooCMSShape.push_back(projection_pt);
	continue;
      }
      if(TString(histoEfficiency->GetName()).Contains("Exp")){
	projectionSF_Exp.push_back(projection_pt);
	continue;
      }
    }
  }
}

/// make tag and probe fit
void makeTagAndProbeFits(const map<string,TFile*> & tagAndProbeFits, const string & outputDIR, const string & typeID, const float & lumi){

  int iFile = 0;
  TCanvas* canvas = new TCanvas("canvasTagAndProbe","",1100,600);
  canvas->cd();
  TPad* pad1 = new TPad("pad1TagAndProbe","",0,0,0.5,1);
  pad1->Draw();
  TPad* pad2 = new TPad("pad2TagAndProbe","",0.5,0,1,1);
  pad2->Draw();
  
  for(auto imap : tagAndProbeFits){
    iFile++;    
    RooWorkspace* workspace = (RooWorkspace*) imap.second->FindObjectAny("w");
    if(workspace == NULL or not workspace) continue;
    RooRealVar* mass = (RooRealVar*) workspace->obj("mass");
    mass->setBins(40);
    // for better plotting
    RooDataHist* histPass = (RooDataHist*) workspace->data("passDataHist");
    RooDataHist* histFail = (RooDataHist*) workspace->data("failDataHist");
    // total pdf pass and fail
    RooAddPdf* pdfPass = (RooAddPdf*) workspace->obj("pdfPass");
    RooAddPdf* pdfFail = (RooAddPdf*) workspace->obj("pdfFail");
    RooAbsPdf* backgroundFail = (RooAddPdf*) workspace->obj("backgroundFail");
    RooAbsPdf* backgroundPass = (RooAddPdf*) workspace->obj("backgroundPass");    
    RooFormulaVar* nBkgFail = (RooFormulaVar*) workspace->obj("nBkgFail");
    RooFormulaVar* nSignalFail = (RooFormulaVar*) workspace->obj("nSigFail");
    RooFormulaVar* nBkgPass = (RooFormulaVar*) workspace->obj("nBkgPass");
    RooFormulaVar* nSignalPass = (RooFormulaVar*) workspace->obj("nSigPass");

    canvas->cd();
    pad1->cd();
    RooPlot* framePass = mass->frame();
    RooFitResult* fitResult = (RooFitResult*) workspace->obj("fitresults");
    RooArgList parlist = fitResult->floatParsFinal();
    framePass->SetTitle("");
    framePass->GetXaxis()->SetTitle("mass (GeV)");
    framePass->GetYaxis()->SetTitle("Events / GeV");
    framePass->GetYaxis()->SetTitleOffset(1.1);
    histPass->plotOn(framePass,RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::LineColor(kBlack),RooFit::DrawOption("EP"),RooFit::DataError(RooAbsData::Poisson),RooFit::Name(histPass->GetName()));
    backgroundPass->plotOn(framePass,RooFit::LineColor(kRed),RooFit::DrawOption("L"),RooFit::Normalization(nBkgPass->getVal(),RooAbsReal::NumEvent),RooFit::Name(backgroundPass->GetName()));
    pdfPass->plotOn(framePass,RooFit::LineColor(kBlue),RooFit::DrawOption("L"),RooFit::Normalization(nBkgPass->getVal()+nSignalPass->getVal(),RooAbsReal::NumEvent),RooFit::Name(pdfPass->GetName()));
    histPass->plotOn(framePass,RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::LineColor(kBlack),RooFit::DrawOption("EP"),RooFit::DataError(RooAbsData::Poisson));
    float chi2 = framePass->chiSquare(parlist.getSize());
    framePass->Draw();
    TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);    
    leg->AddEntry((TObject*)0,"Passing Probes","");
    for(int obj = 0; obj < int(framePass->numItems()); obj++){
      string name = framePass->nameOf(obj);
      TObject* theObj = framePass->getObject(obj);
      if(name == string(backgroundPass->GetName()))
	leg->AddEntry(theObj,"Background Pdf","L");
      else if(name == string(pdfPass->GetName()))
	leg->AddEntry(theObj,"Total S+B Pdf","L");
      else if(name == string(histPass->GetName()))
	leg->AddEntry(theObj,"Data","PE");
    }
    leg->AddEntry((TObject*)0,Form("#chi^{2} = %.2f",chi2),"");
    leg->Draw("same");
    CMS_lumi(pad1,string(Form("%.2f",lumi)),true,false);

    canvas->cd();
    pad2->cd();
    RooPlot* frameFail = mass->frame();
    frameFail->SetTitle("");
    frameFail->GetXaxis()->SetTitle("mass (GeV)");
    frameFail->GetYaxis()->SetTitle("Events / GeV");
    frameFail->GetYaxis()->SetTitleOffset(1.1);
    histFail->plotOn(frameFail,RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::LineColor(kBlack),RooFit::DrawOption("EP"),RooFit::DataError(RooAbsData::Poisson),RooFit::Name(histFail->GetName()));
    backgroundFail->plotOn(frameFail,RooFit::LineColor(kRed),RooFit::DrawOption("L"),RooFit::Normalization(nBkgFail->getVal(),RooAbsReal::NumEvent),RooFit::Name(backgroundFail->GetName()));
    pdfFail->plotOn(frameFail,RooFit::LineColor(kBlue),RooFit::DrawOption("L"),RooFit::Normalization(nBkgFail->getVal()+nSignalFail->getVal(),RooAbsReal::NumEvent),RooFit::Name(pdfFail->GetName()));
    histFail->plotOn(frameFail,RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::LineColor(kBlack),RooFit::DrawOption("EP"),RooFit::DataError(RooAbsData::Poisson));
    chi2 = frameFail->chiSquare(parlist.getSize());
    frameFail->Draw();
    TLegend* leg2 = new TLegend(0.6,0.7,0.9,0.9);
    leg2->SetFillStyle(0);
    leg2->SetFillColor(0);
    leg2->SetBorderSize(0);
    leg2->AddEntry((TObject*)0,"Failing Probes","");
    for(int obj = 0; obj < int(frameFail->numItems()); obj++){
      string name = frameFail->nameOf(obj);
      TObject* theObj = frameFail->getObject(obj);
      if(name == string(backgroundFail->GetName()))
	leg2->AddEntry(theObj,"Background Pdf","L");
      else if(name == string(pdfFail->GetName()))
	leg2->AddEntry(theObj,"Total S+B Pdf","L");
      else if(name == string(histFail->GetName()))
	leg2->AddEntry(theObj,"Data","PE");
    }
    leg2->AddEntry((TObject*)0,Form("#chi^{2} = %.2f",chi2),"");
    leg2->Draw("same");
    CMS_lumi(pad2,string(Form("%.2f",lumi)),true,false);    
    //    canvas->SaveAs((outputDIR+"/"+imap.first+".png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/"+imap.first+".pdf").c_str(),"pdf");
    
    delete histPass;
    delete histFail;
    delete framePass;
    delete frameFail;
    delete nBkgFail;
    delete nSignalFail;
    delete nBkgPass;
    delete nSignalPass;
    delete backgroundFail;
    delete backgroundPass;
    delete pdfPass;
    delete pdfFail;
    delete mass;
  }    
}


/// mkae scale factor for lepton ID
void makeLeptonIDScaleFactors(string inputTagAndProbeFitDIR, // direcory with root files with tag and probe results
			      string inputTagAndProbeTemplateDIR, // input directory with templates and efficiencies extracted from the MC
			      string leptonType, // muons or electrons
			      string typeID, // id type: looseid, tightid, vetoid ..
			      string outputDIR, // where plots and root files with 2D scale factors are stored,
			      string inputTagAndProbeTemplateDIR_Alternative = "",
			      float  lumi = 0.86,
			      bool   addTurnOnFits = false,
			      bool   addAnalyticalFits = false
			      ){


  gSystem->Load("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeTagAndProbe/PDFs/RooErfExpPdf_cc.so");
  gSystem->Load("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeTagAndProbe/PDFs/RooCMSShape_cc.so");
 
  addTurnOnFits_ = addTurnOnFits;

  if(leptonType != "muon" and leptonType!= "electron" and leptonType!= "photon"){
    cerr<<"Wrong lepton type --> electron, muon or photon --> return "<<endl;
    return;
  }

  // start wit the MC
  system(("mkdir -p "+outputDIR).c_str());

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  TGaxis::SetMaxDigits(3);

  if(typeID != "looseid" and typeID != "tightid" and typeID != "vetoid" and typeID != "mediumid"){
    cerr<<"typeID is not valid --> please check --> return"<<endl;
    return;
  }

  vector<float> puBins;
  if(leptonType == "muon")
    puBins = nvtxBinMuon;
  else if(leptonType == "electron")
    puBins = nvtxBinElectron;
  else if(leptonType == "photon")
    puBins = nvtxBinPhoton;

  /// fill the vector for tag and probe templates
  cout<<"Run MC efficiency "<<endl;
  vector<TH2F*> histoEfficiencyMC;
  TCanvas* canvas = new TCanvas("canvas","",625,600);
  canvas->cd();
  for(size_t ipu = 0; ipu < puBins.size()-1; ipu++){
    if(leptonType == "muon")
      histoEfficiencyMC.push_back(new TH2F(("efficiencyMC_muon_"+typeID+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
					   etaBinMuon.size()-1,&etaBinMuon[0],ptBinMuon.size()-1,&ptBinMuon[0])); 
    else if(leptonType == "electron")
      histoEfficiencyMC.push_back(new TH2F(("efficiencyMC_electron_"+typeID+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
					   etaBinElectron.size()-1,&etaBinElectron[0],ptBinElectron.size()-1,&ptBinElectron[0])); 
    else if(leptonType == "photon")
      histoEfficiencyMC.push_back(new TH2F(("efficiencyMC_photon_"+typeID+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
					   etaBinPhoton.size()-1,&etaBinPhoton[0],ptBinPhoton.size()-1,&ptBinPhoton[0])); 
    // fill histo
    histoEfficiencyMC.back()->Sumw2();
    fillEfficiencyMC(histoEfficiencyMC.back(),inputTagAndProbeTemplateDIR,leptonType,typeID,ipu);
    plotEfficiency(canvas,histoEfficiencyMC.back(),outputDIR,leptonType,lumi,true,false,ipu);
  }

  vector<TH2F*> histoEfficiencyMC_Alternative;
  if(inputTagAndProbeTemplateDIR_Alternative != ""){ //means that there is also an alternative set of templates
    for(size_t ipu = 0; ipu < puBins.size()-1; ipu++){
      if(leptonType == "muon")
	histoEfficiencyMC_Alternative.push_back(new TH2F(("efficiencyMC_Alternative_muon_"+typeID+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
							 etaBinMuon.size()-1,&etaBinMuon[0],ptBinMuon.size()-1,&ptBinMuon[0]));
    else if(leptonType == "electron")
      histoEfficiencyMC_Alternative.push_back(new TH2F(("efficiencyMC_Alternative_electron_"+typeID+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
						       etaBinElectron.size()-1,&etaBinElectron[0],ptBinElectron.size()-1,&ptBinElectron[0]));
    else if(leptonType == "photon")
      histoEfficiencyMC_Alternative.push_back(new TH2F(("efficiencyMC_Alternative_photon_"+typeID+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
						       etaBinPhoton.size()-1,&etaBinPhoton[0],ptBinPhoton.size()-1,&ptBinPhoton[0]));

      histoEfficiencyMC_Alternative.back()->Sumw2();
      fillEfficiencyMC(histoEfficiencyMC_Alternative.back(),inputTagAndProbeTemplateDIR_Alternative,leptonType,typeID,ipu);
      plotEfficiency(canvas,histoEfficiencyMC_Alternative.back(),outputDIR,leptonType,lumi,true,false,ipu);  
    }
  }
  
  cout<<"End MC efficiency "<<endl;
  
  // switch to data --> make also the passing and failing plots
  cout<<"Run Data efficiency "<<endl;
  vector<TH2F*> histoEfficiencyDATA_RooCMSShape;
  vector<TH2F*> histoEfficiencyDATA_Exp;
  vector<TH2F*> histoEfficiencyDATA_Alternative;
  vector<TH2F*> histoEfficiencyDATA_Analytical;

  for(size_t ipu = 0; ipu < puBins.size()-1; ipu++){
    if(leptonType == "muon"){
      histoEfficiencyDATA_RooCMSShape.push_back( new TH2F(("efficiencyDATA_muon_"+typeID+"_RooCMSShape"+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
							  etaBinMuon.size()-1,&etaBinMuon[0],ptBinMuon.size()-1,&ptBinMuon[0]));
      histoEfficiencyDATA_Exp.push_back( new TH2F(("efficiencyDATA_muon_"+typeID+"_Exp"+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
						  etaBinMuon.size()-1,&etaBinMuon[0],ptBinMuon.size()-1,&ptBinMuon[0]));
    }
    else if(leptonType == "electron"){
      histoEfficiencyDATA_RooCMSShape .push_back( new TH2F(("efficiencyDATA_electron_"+typeID+"_RooCMSShape"+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
							   etaBinElectron.size()-1,&etaBinElectron[0],ptBinElectron.size()-1,&ptBinElectron[0]));
      histoEfficiencyDATA_Exp .push_back( new TH2F(("efficiencyDATA_electron_"+typeID+"_Exp"+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
						   etaBinElectron.size()-1,&etaBinElectron[0],ptBinElectron.size()-1,&ptBinElectron[0]));
    }
    else if(leptonType == "photon"){
      histoEfficiencyDATA_RooCMSShape .push_back( new TH2F(("efficiencyDATA_photon_"+typeID+"_RooCMSShape"+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
							   etaBinPhoton.size()-1,&etaBinPhoton[0],ptBinPhoton.size()-1,&ptBinPhoton[0]));
      histoEfficiencyDATA_Exp .push_back( new TH2F(("efficiencyDATA_photon_"+typeID+"_Exp"+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
						   etaBinPhoton.size()-1,&etaBinPhoton[0],ptBinPhoton.size()-1,&ptBinPhoton[0]));
    }
    
  
    histoEfficiencyDATA_RooCMSShape.back()->Sumw2();
    histoEfficiencyDATA_Exp.back()->Sumw2();
    // fill histos
    fillEfficiencyData(histoEfficiencyDATA_RooCMSShape.back(),inputTagAndProbeFitDIR,leptonType,typeID,ipu,"RooCMSShape");  
    fillEfficiencyData(histoEfficiencyDATA_Exp.back(),inputTagAndProbeFitDIR,leptonType,typeID,ipu,"Exp");
    //plot histos
    plotEfficiency(canvas,histoEfficiencyDATA_RooCMSShape.back(),outputDIR,leptonType,lumi,false,false,ipu);
    plotEfficiency(canvas,histoEfficiencyDATA_Exp.back(),outputDIR,leptonType,lumi,false,false,ipu);
  }

  if(inputTagAndProbeTemplateDIR_Alternative != ""){    
    for(size_t ipu = 0; ipu < puBins.size()-1; ipu++){
      if(leptonType == "muon")
	histoEfficiencyDATA_Alternative.push_back( new TH2F(("efficiencyDATA_muon_"+typeID+"_Alternative"+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
							    etaBinMuon.size()-1,&etaBinMuon[0],ptBinMuon.size()-1,&ptBinMuon[0]));
      else if(leptonType == "electron")
	histoEfficiencyDATA_Alternative.push_back( new TH2F(("efficiencyDATA_electron_"+typeID+"_Alternative"+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
							    etaBinElectron.size()-1,&etaBinElectron[0],ptBinElectron.size()-1,&ptBinElectron[0]));
      else if(leptonType == "photon")
	histoEfficiencyDATA_Alternative.push_back( new TH2F(("efficiencyDATA_photon_"+typeID+"_Alternative"+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
							    etaBinPhoton.size()-1,&etaBinPhoton[0],ptBinPhoton.size()-1,&ptBinPhoton[0]));
    
      histoEfficiencyDATA_Alternative.back()->Sumw2();
      // fill histos
      fillEfficiencyData(histoEfficiencyDATA_Alternative.back(),inputTagAndProbeFitDIR,leptonType,typeID,ipu,"Alternative");
      plotEfficiency(canvas,histoEfficiencyDATA_Alternative.back(),outputDIR,leptonType,lumi,false,false,ipu);
    }
  }
  
  if(addAnalyticalFits){
    for(size_t ipu = 0; ipu < puBins.size()-1; ipu++){
      if(leptonType == "muon")
	histoEfficiencyDATA_Analytical.push_back(new TH2F(("efficiencyDATA_muon_"+typeID+"_Analytical"+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
							  etaBinMuon.size()-1,&etaBinMuon[0],ptBinMuon.size()-1,&ptBinMuon[0]));
      else if(leptonType == "electron")
	histoEfficiencyDATA_Analytical.push_back(new TH2F(("efficiencyDATA_electron_"+typeID+"_Analytical"+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
							  etaBinElectron.size()-1,&etaBinElectron[0],ptBinElectron.size()-1,&ptBinElectron[0]));
      else if(leptonType == "photon")
	histoEfficiencyDATA_Analytical.push_back(new TH2F(("efficiencyDATA_photon_"+typeID+"_Analytical"+"_pu_"+Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))).c_str(),"",
							  etaBinPhoton.size()-1,&etaBinPhoton[0],ptBinPhoton.size()-1,&ptBinPhoton[0]));
      
      histoEfficiencyDATA_Analytical.back()->Sumw2();
      // fill histos
      fillEfficiencyData(histoEfficiencyDATA_Analytical.back(),inputTagAndProbeFitDIR,leptonType,typeID,ipu,"Analytical");
      plotEfficiency(canvas,histoEfficiencyDATA_Analytical.back(),outputDIR,leptonType,lumi,false,false,ipu);
    }
  }

  cout<<"End Data efficiency "<<endl;
  
  // calculate the scale factor and store it
  cout<<"Calculate SF "<<endl;
  vector<TH2F*> histoEfficiencySF_RooCMSShape;
  vector<TH2F*> histoEfficiencySF_Exp;
  vector<TH2F*> histoEfficiencySF_Alternative;
  vector<TH2F*> histoEfficiencySF_Analytical;

  for(size_t ipu = 0; ipu < puBins.size()-1; ipu++){
    histoEfficiencySF_RooCMSShape.push_back((TH2F*) histoEfficiencyDATA_RooCMSShape.at(ipu)->Clone(Form("histoEfficiencySF_RooCMSShape_pu_%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))));
    histoEfficiencySF_RooCMSShape.back()->Reset();
    histoEfficiencySF_RooCMSShape.back()->SetName(Form("histoEfficiencySF_RooCMSShape_pu_%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1)));
    histoEfficiencySF_RooCMSShape.back()->Add(histoEfficiencyDATA_RooCMSShape.at(ipu));
    histoEfficiencySF_RooCMSShape.back()->Divide(histoEfficiencyMC.at(ipu));
    histoEfficiencySF_Exp.push_back((TH2F*) histoEfficiencyDATA_Exp.at(ipu)->Clone(Form("histoEfficiencySF_Exp_pu_%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))));
    histoEfficiencySF_Exp.back()->Reset();
    histoEfficiencySF_Exp.back()->SetName(Form("histoEfficiencySF_Exp_pu_%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1)));
    histoEfficiencySF_Exp.back()->Add(histoEfficiencyDATA_Exp.at(ipu));
    histoEfficiencySF_Exp.back()->Divide(histoEfficiencyMC.at(ipu));
    
    plotEfficiency(canvas,histoEfficiencySF_RooCMSShape.back(),outputDIR,leptonType,lumi,false,true,ipu);
    plotEfficiency(canvas,histoEfficiencySF_Exp.back(),outputDIR,leptonType,lumi,false,true,ipu);
  
    if(histoEfficiencyDATA_Alternative.size() != 0 and histoEfficiencyMC_Alternative.size() != 0){
      histoEfficiencySF_Alternative.push_back((TH2F*) histoEfficiencyDATA_Alternative.at(ipu)->Clone(Form("histoEfficiencySF_Alternative_pu_%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))));
      histoEfficiencySF_Alternative.back()->Reset();
      histoEfficiencySF_Alternative.back()->SetName(Form("histoEfficiencySF_Alternative_pu_%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1)));
      histoEfficiencySF_Alternative.back()->Add(histoEfficiencyDATA_Alternative.at(ipu));
      histoEfficiencySF_Alternative.back()->Divide(histoEfficiencyMC_Alternative.at(ipu));
      plotEfficiency(canvas,histoEfficiencySF_Alternative.back(),outputDIR,leptonType,lumi,false,true,ipu);
    }

    if(addAnalyticalFits){
      histoEfficiencySF_Analytical.push_back((TH2F*) histoEfficiencyDATA_Analytical.at(ipu)->Clone(Form("histoEfficiencySF_Analytical_pu_%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1))));
      histoEfficiencySF_Analytical.back()->Reset();
      histoEfficiencySF_Analytical.back()->SetName(Form("histoEfficiencySF_Analytical_pu_%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1)));
      histoEfficiencySF_Analytical.back()->Add(histoEfficiencyDATA_Analytical.at(ipu));
      histoEfficiencySF_Analytical.back()->Divide(histoEfficiencyMC.at(ipu));
      plotEfficiency(canvas,histoEfficiencySF_Analytical.back(),outputDIR,leptonType,lumi,false,true,ipu);
    }
  }

  TCanvas* can = new TCanvas("can","",600,650);
  can->cd();
  TPad *   pad1 = new TPad("pad1","", 0, 0.3,1.0,1.0);
  pad1->Draw();
  TPad *   pad2 = new TPad("pad2","", 0, 0.15,1.0,0.3);
  pad2->Draw();
  TPad *   pad3 = new TPad("pad3","", 0, 0.0,1.0,0.15);
  pad3->Draw();

  vector<TH1F*> sysUnc;
  vector<TH1F*> sysUnc_alt;
  vector<TH1F*> sysUnc_alt2;

  for(size_t iproj = 0; iproj < projectionMC.size(); iproj++){

    projectionMC.at(iproj)->SetLineColor(kRed);
    projectionMC.at(iproj)->SetMarkerColor(kRed);
    projectionMC.at(iproj)->SetMarkerSize(0.75);
    TF1* funz = NULL;
    if(addTurnOnFits_){
      funz = (TF1*) projectionMC.at(iproj)->GetListOfFunctions()->At(0);
      funz->SetLineColor(kRed);
    }
    pad1->cd();
    projectionMC.at(iproj)->Draw("E1P");

    if(projectionMC_Alternative.size() > iproj){
      projectionMC_Alternative.at(iproj)->SetLineColor(kCyan+1);
      projectionMC_Alternative.at(iproj)->SetMarkerColor(kCyan+1);
      projectionMC_Alternative.at(iproj)->SetMarkerSize(0.75);
      projectionMC_Alternative.at(iproj)->SetMarkerStyle(24);
      TF1* funz = NULL;
      if(addTurnOnFits_){
	funz = (TF1*) projectionMC_Alternative.at(iproj)->GetListOfFunctions()->At(0);
	funz->SetLineColor(kCyan);
      }
      pad1->cd();
      projectionMC_Alternative.at(iproj)->Draw("E1Psame");
    }

    if(projectionDATA_RooCMSShape.size() < iproj) continue;

    projectionDATA_RooCMSShape.at(iproj)->SetLineColor(kBlue);
    projectionDATA_RooCMSShape.at(iproj)->SetMarkerColor(kBlue);
    projectionDATA_RooCMSShape.at(iproj)->SetMarkerSize(0.75);
    if(addTurnOnFits_){
      funz = (TF1*) projectionDATA_RooCMSShape.at(iproj)->GetListOfFunctions()->At(0);
      funz->SetLineColor(kBlue);
    }
    projectionDATA_RooCMSShape.at(iproj)->Draw("E1Psame");

    if(projectionDATA_Exp.size() < iproj) continue;

    projectionDATA_Exp.at(iproj)->SetLineColor(kBlack);
    projectionDATA_Exp.at(iproj)->SetMarkerColor(kBlack);
    projectionDATA_Exp.at(iproj)->SetMarkerSize(0.75);
    if(addTurnOnFits_){
      funz = (TF1*) projectionDATA_Exp.at(iproj)->GetListOfFunctions()->At(0);
      funz->SetLineColor(kBlack);
    }
    projectionDATA_Exp.at(iproj)->Draw("E1Psame");

    if(projectionDATA_Alternative.size() >iproj){
      projectionDATA_Alternative.at(iproj)->SetLineColor(kGreen+1);
      projectionDATA_Alternative.at(iproj)->SetMarkerColor(kGreen+1);
      projectionDATA_Alternative.at(iproj)->SetMarkerSize(0.75);
      projectionDATA_Alternative.at(iproj)->SetMarkerStyle(24);
      if(addTurnOnFits_){
	funz = (TF1*) projectionDATA_Alternative.at(iproj)->GetListOfFunctions()->At(0);
	funz->SetLineColor(kGreen+1);
      }
      projectionDATA_Alternative.at(iproj)->Draw("E1Psame");
    }

    if(projectionDATA_Analytical.size() >iproj){
      projectionDATA_Analytical.at(iproj)->SetLineColor(kGray+1);
      projectionDATA_Analytical.at(iproj)->SetMarkerColor(kGray+1);
      projectionDATA_Analytical.at(iproj)->SetMarkerSize(0.75);
      projectionDATA_Analytical.at(iproj)->SetMarkerStyle(24);
      if(addTurnOnFits_){
	funz = (TF1*) projectionDATA_Analytical.at(iproj)->GetListOfFunctions()->At(0);
	funz->SetLineColor(kGray+1);
      }
      projectionDATA_Analytical.at(iproj)->Draw("E1Psame");
    }
  
    CMS_lumi(pad1,string(Form("%.2f",lumi)),true);
    
    TLegend* leg = new TLegend(0.65,0.25,0.90,0.42);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(projectionMC.at(iproj),"MC DY #rightarrow ll","PE");
    if(projectionMC_Alternative.size() >iproj)
       leg->AddEntry(projectionMC_Alternative.at(iproj),"MC Alternative DY #rightarrow ll","PE");
    leg->AddEntry(projectionDATA_RooCMSShape.at(iproj),"DATA : Bkg RooCMS","PE");
    leg->AddEntry(projectionDATA_Exp.at(iproj),"DATA : Bkg Exp","PE");
    if(projectionDATA_Alternative.size() >iproj)
      leg->AddEntry(projectionDATA_Alternative.at(iproj),"DATA: Alternative signal template","PE");
    if(projectionDATA_Analytical.size() >iproj)
      leg->AddEntry(projectionDATA_Analytical.at(iproj),"DATA: Analytical fit","PE");
    leg->Draw("same");

    pad2->cd();    
    pad2->SetGridy();
    projectionSF_RooCMSShape.at(iproj)->SetLineColor(kBlue);
    projectionSF_RooCMSShape.at(iproj)->SetMarkerColor(kBlue);
    projectionSF_RooCMSShape.at(iproj)->SetMarkerSize(0.75);
    projectionSF_Exp.at(iproj)->SetLineColor(kBlack);
    projectionSF_Exp.at(iproj)->SetMarkerColor(kBlack);
    projectionSF_Exp.at(iproj)->SetMarkerSize(0.75);
    projectionSF_RooCMSShape.at(iproj)->GetYaxis()->SetTitle("SF");
    projectionSF_RooCMSShape.at(iproj)->GetXaxis()->SetTitle("");
    projectionSF_RooCMSShape.at(iproj)->GetYaxis()->CenterTitle();
    projectionSF_RooCMSShape.at(iproj)->SetNdivisions(510);
    projectionSF_RooCMSShape.at(iproj)->GetXaxis()->SetLabelSize(0);
    projectionSF_RooCMSShape.at(iproj)->GetYaxis()->SetLabelSize(0.15);
    projectionSF_RooCMSShape.at(iproj)->GetYaxis()->SetTitleSize(0.20);
    projectionSF_RooCMSShape.at(iproj)->GetYaxis()->SetTitleOffset(0.25);

    if(projectionSF_Alternative.size() > iproj){
      projectionSF_Alternative.at(iproj)->SetLineColor(kGreen+1);
      projectionSF_Alternative.at(iproj)->SetMarkerColor(kGreen+1);
      projectionSF_Alternative.at(iproj)->SetMarkerSize(0.75);
      projectionSF_Alternative.at(iproj)->SetMarkerStyle(24);
    }

    if(projectionSF_Analytical.size() > iproj){
      projectionSF_Analytical.at(iproj)->SetLineColor(kGray+1);
      projectionSF_Analytical.at(iproj)->SetMarkerColor(kGray+1);
      projectionSF_Analytical.at(iproj)->SetMarkerSize(0.75);
      projectionSF_Analytical.at(iproj)->SetMarkerStyle(24);
    }

    projectionSF_RooCMSShape.at(iproj)->Draw("E1P");
    projectionSF_Exp.at(iproj)->Draw("E1Psame");
    if(projectionSF_Alternative.size() > iproj)
      projectionSF_Alternative.at(iproj)->Draw("E1Psame");
    if(projectionSF_Analytical.size() > iproj)
      projectionSF_Analytical.at(iproj)->Draw("E1Psame");

    pad3->cd();
    pad3->SetGridy();
    TH1F* sysError = (TH1F*) projectionSF_RooCMSShape.at(iproj)->Clone("Uncertainty");
    TH1F* uncPlot  = (TH1F*) projectionSF_RooCMSShape.at(iproj)->Clone("UncPlot");
    sysError->Divide(projectionSF_Exp.at(iproj));
    uncPlot->Divide(projectionSF_Exp.at(iproj));

    TH1F* sysError_alt = (TH1F*) projectionSF_RooCMSShape.at(iproj)->Clone("Uncertainty_alt");
    if(projectionSF_Alternative.size() > iproj)
      sysError_alt->Divide(projectionSF_Alternative.at(iproj));    

    TH1F* sysError_alt2 = (TH1F*) projectionSF_RooCMSShape.at(iproj)->Clone("Uncertainty_alt2");
    if(projectionSF_Analytical.size() > iproj)
      sysError_alt2->Divide(projectionSF_Analytical.at(iproj));    

    float max = 0;
    for(int iBin = 0; iBin < projectionSF_RooCMSShape.at(iproj)->GetNbinsX(); iBin++){
      // stotal error fit + shape sys
      sysError->SetBinError(iBin+1,fabs(1-sysError->GetBinContent(iBin+1)));
      sysError_alt->SetBinError(iBin+1,fabs(1-sysError_alt->GetBinContent(iBin+1)));
      sysError_alt2->SetBinError(iBin+1,fabs(1-sysError_alt2->GetBinContent(iBin+1)));
      
      if(projectionSF_Alternative.size() > iproj)
	uncPlot->SetBinError(iBin+1,sqrt(fabs(1-sysError->GetBinContent(iBin+1))*fabs(1-sysError->GetBinContent(iBin+1))+
					 fabs(1-sysError_alt->GetBinContent(iBin+1))*fabs(1-sysError_alt->GetBinContent(iBin+1))+
					 projectionSF_RooCMSShape.at(iproj)->GetBinError(iBin+1)*projectionSF_RooCMSShape.at(iproj)->GetBinError(iBin+1)));
      else
	uncPlot->SetBinError(iBin+1,sqrt(fabs(1-sysError->GetBinContent(iBin+1))*fabs(1-sysError->GetBinContent(iBin+1))+projectionSF_RooCMSShape.at(iproj)->GetBinError(iBin+1)*projectionSF_RooCMSShape.at(iproj)->GetBinError(iBin+1)));


      if(projectionSF_Analytical.size() > iproj)
	uncPlot->SetBinError(iBin+1,sqrt(uncPlot->GetBinError(iBin+1)*uncPlot->GetBinError(iBin+1)+
					 fabs(1-sysError_alt2->GetBinContent(iBin+1))*fabs(1-sysError_alt2->GetBinContent(iBin+1))));
      
      
      if(uncPlot->GetBinError(iBin+1) > max)
	max = uncPlot->GetBinError(iBin+1);
      if(uncPlot->GetBinError(iBin+1) < 0.001)
	uncPlot->SetBinError(iBin+1,0.001);
      uncPlot->SetBinContent(iBin+1,1);
    }
    uncPlot->GetYaxis()->SetTitle("Uncertainty");
    uncPlot->GetYaxis()->CenterTitle();
    uncPlot->GetXaxis()->SetTitle("");

    uncPlot->SetFillColor(kGray);
    uncPlot->GetYaxis()->SetRangeUser(1-fabs(max)-0.005,1+fabs(max)+0.005);
    TLine* line = new TLine(uncPlot->GetBinLowEdge(1),1,uncPlot->GetBinLowEdge(sysError->GetNbinsX()+1),1);
    line->SetLineColor(kBlack);
    line->SetLineStyle(7);
    line->SetLineWidth(2);
    uncPlot->SetMarkerSize(0);
    uncPlot->Draw("E2");
    line->Draw("Lsame");
    pad3->RedrawAxis("sameaxis");    
    sysUnc.push_back(sysError); // store histogram with systematic uncertainty
    sysUnc_alt.push_back(sysError_alt); // store histogram with systematic uncertainty
    sysUnc_alt2.push_back(sysError_alt2); // store histogram with systematic uncertainty
    //    can->SaveAs((outputDIR+"/summaryScaleFactor_vs_pt_etaBin_"+string(Form("%d",int(iproj)))+".png").c_str(),"png");
    can->SaveAs((outputDIR+"/summaryScaleFactor_vs_pt_etaBin_"+string(Form("%d",int(iproj)))+".pdf").c_str(),"pdf");    
  }      

  
  // use RooCMSShape as central value + uncertaint from alternative background description --> total 2D hist SF stored in root file.
  cout<<"Scale factor Nominal "<<endl;
  for(size_t ipu = 0; ipu < histoEfficiencySF_RooCMSShape.size(); ipu++){
    for(int iBinX = 0; iBinX < histoEfficiencySF_RooCMSShape.at(ipu)->GetNbinsX(); iBinX++){
      for(int iBinY = 0; iBinY < histoEfficiencySF_RooCMSShape.at(ipu)->GetNbinsY(); iBinY++){      
	histoEfficiencySF_RooCMSShape.at(ipu)->SetBinError(iBinX+1,iBinY+1,sqrt(histoEfficiencySF_RooCMSShape.at(ipu)->GetBinError(iBinX+1,iBinY+1)*histoEfficiencySF_RooCMSShape.at(ipu)->GetBinError(iBinX+1,iBinY+1)+sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)+sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1+sysUnc_alt2.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc_alt2.at(iBinX+ipu)->GetBinError(iBinY+1))));    
	cout<<"Scale Factor etaBin ["<<histoEfficiencySF_RooCMSShape.at(ipu)->GetYaxis()->GetBinLowEdge(iBinY+1)<<" - "<<histoEfficiencySF_RooCMSShape.at(ipu)->GetYaxis()->GetBinLowEdge(iBinY+2)<<"]"<<" ptBin ["<<histoEfficiencySF_RooCMSShape.at(ipu)->GetXaxis()->GetBinLowEdge(iBinX+1)<<" - "<<histoEfficiencySF_RooCMSShape.at(ipu)->GetXaxis()->GetBinLowEdge(iBinX+2)<<"]"<<" : "<<histoEfficiencySF_RooCMSShape.at(ipu)->GetBinContent(iBinX+1,iBinY+1)<<" err "<<histoEfficiencySF_RooCMSShape.at(ipu)->GetBinError(iBinX+1,iBinY+1)<<endl;
      }
    }
  }

  cout<<"Scale factor Exp "<<endl;
  for(size_t ipu = 0; ipu < histoEfficiencySF_Exp.size(); ipu++){    
    for(int iBinX = 0; iBinX < histoEfficiencySF_Exp.at(ipu)->GetNbinsX(); iBinX++){
      for(int iBinY = 0; iBinY < histoEfficiencySF_Exp.at(ipu)->GetNbinsY(); iBinY++){      
	histoEfficiencySF_Exp.at(ipu)->SetBinError(iBinX+1,iBinY+1,sqrt(histoEfficiencySF_Exp.at(ipu)->GetBinError(iBinX+1,iBinY+1)*histoEfficiencySF_Exp.at(ipu)->GetBinError(iBinX+1,iBinY+1)+sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)+sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)+sysUnc_alt2.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc_alt2.at(iBinX+ipu)->GetBinError(iBinY+1)));      
      }
    }
  }

  cout<<"Scale factor Alternative "<<endl;
  if(histoEfficiencySF_Alternative.size() != 0){
    for(size_t ipu = 0; ipu < histoEfficiencySF_Alternative.size(); ipu++){
      for(int iBinX = 0; iBinX < histoEfficiencySF_Alternative.at(ipu)->GetNbinsX(); iBinX++){
	for(int iBinY = 0; iBinY < histoEfficiencySF_Alternative.at(ipu)->GetNbinsY(); iBinY++){      
	  histoEfficiencySF_Alternative.at(ipu)->SetBinError(iBinX+1,iBinY+1,sqrt(histoEfficiencySF_Alternative.at(ipu)->GetBinError(iBinX+1,iBinY+1)*histoEfficiencySF_Alternative.at(ipu)->GetBinError(iBinX+1,iBinY+1)+sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)+sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)+sysUnc_alt2.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc_alt2.at(iBinX+ipu)->GetBinError(iBinY+1)));      
	}
      }
    }
  }

  cout<<"Scale factor Analytical "<<endl;
  if(histoEfficiencySF_Analytical.size()!=0){
    for(size_t ipu = 0; ipu < histoEfficiencySF_Analytical.size(); ipu++){
      for(int iBinX = 0; iBinX < histoEfficiencySF_Analytical.at(ipu)->GetNbinsX(); iBinX++){
	for(int iBinY = 0; iBinY < histoEfficiencySF_Analytical.at(ipu)->GetNbinsY(); iBinY++){      
	  histoEfficiencySF_Analytical.at(ipu)->SetBinError(iBinX+1,iBinY+1,sqrt(histoEfficiencySF_Analytical.at(ipu)->GetBinError(iBinX+1,iBinY+1)*histoEfficiencySF_Analytical.at(ipu)->GetBinError(iBinX+1,iBinY+1)+sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)+sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)+sysUnc_alt2.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc_alt2.at(iBinX+ipu)->GetBinError(iBinY+1)));      
	}
      }
    }
  }
  
  cout<<"Output plots "<<endl;  
  TFile* outputScaleFactor = new TFile((outputDIR+"/scaleFactor_"+leptonType+"_"+typeID+".root").c_str(),"RECREATE");
  outputScaleFactor->cd();

  for(int ipu = 0; ipu < puBins.size()-1; ipu++){
    histoEfficiencySF_RooCMSShape.at(ipu)->Write(("scaleFactor_"+leptonType+"_"+typeID+"_RooCMSShape"+string(Form("_pu_%d_%d",int(puBins.at(ipu)),int(puBins.at(ipu+1))))).c_str());
    histoEfficiencySF_Exp.at(ipu)->Write(("scaleFactor_"+leptonType+"_"+typeID+"_Exp"+string(Form("_pu_%d_%d",int(puBins.at(ipu)),int(puBins.at(ipu+1))))).c_str());
    if(histoEfficiencySF_Alternative.size() != 0)
      histoEfficiencySF_Alternative.at(ipu)->Write(("scaleFactor_"+leptonType+"_"+typeID+"_Alternative"+string(Form("_pu_%d_%d",int(puBins.at(ipu)),int(puBins.at(ipu+1))))).c_str());
    if(histoEfficiencySF_Analytical.size() != 0)
      histoEfficiencySF_Analytical.at(ipu)->Write(("scaleFactor_"+leptonType+"_"+typeID+"_Analytical"+string(Form("_pu_%d_%d",int(puBins.at(ipu)),int(puBins.at(ipu+1))))).c_str());
  }
  outputScaleFactor->Close();
  // store fits  
  system(("mkdir -p "+outputDIR+"/RooCMSShape").c_str());
  makeTagAndProbeFits(tagAndProbeFits_RooCMSShape,outputDIR+"/RooCMSShape",typeID,lumi);
  system(("mkdir -p "+outputDIR+"/Exp").c_str());
  makeTagAndProbeFits(tagAndProbeFits_Exp,outputDIR+"/Exp",typeID,lumi);
  system(("mkdir -p "+outputDIR+"/Alternative").c_str());
  if(tagAndProbeFits_Alternative.size() != 0)
    makeTagAndProbeFits(tagAndProbeFits_Alternative,outputDIR+"/Alternative",typeID,lumi);
  system(("mkdir -p "+outputDIR+"/Analytical").c_str());
  if(tagAndProbeFits_Analytical.size() != 0)
    makeTagAndProbeFits(tagAndProbeFits_Analytical,outputDIR+"/Analytical",typeID,lumi);  
}
