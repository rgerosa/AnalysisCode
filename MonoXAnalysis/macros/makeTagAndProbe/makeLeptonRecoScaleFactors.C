#include "../makeTriggerEfficiency/triggerUtils.h"
#include "TnPBinning.h"
#include "../CMS_lumi.h"

// simplified setup
vector<TH1*> projectionMC;
//
vector<TH1*> projectionDATA_RooCMSShape;
vector<TH1*> projectionDATA_Exp;
vector<TH1*> projectionDATA_Analytical;
//
vector<TH1*> projectionSF_RooCMSShape;
vector<TH1*> projectionSF_Exp;
vector<TH1*> projectionSF_Analytical;
//
map<string,TFile*> tagAndProbeFits_RooCMSShape;
map<string,TFile*> tagAndProbeFits_Exp;
map<string,TFile*> tagAndProbeFits_Analytical;
//
bool addTurnOnFits_;
//
static float scalefactor = 1.;

vector<TH2F*> histoEfficiencySF_RooCMSShape;
vector<TH2F*> histoEfficiencySF_Exp;
vector<TH2F*> histoEfficiencySF_Analytical;


// take proper mc efficiency
void fillEfficiencyMC(TH2F* efficiency, const string & directory, const string & leptonType, const string & typeID, const int & ipu){
  
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  TFile* inputFile = NULL;
  vector<float> ptBin;
  vector<float> etaBin;
  vector<float> puBin;
  string leptype;
  if(leptonType == "muon"){
    leptype = "muon";
    ptBin  = ptBinMuonReco;
    etaBin = etaBinMuonReco; 
    puBin = nvtxBinMuonReco;
 }
  else if(leptonType == "electron"){    
    leptype = "electron";
    ptBin = ptBinElectronReco;
    etaBin = etaBinElectronReco;
    puBin = nvtxBinElectronReco;
  }

  //////////////////
  for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++){
    for(size_t ieta = 0; ieta < etaBin.size()-1; ieta++){
      system(("ls "+directory+" | grep root | grep "+leptype+" | grep "+typeID+" | grep pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.2f",etaBin.at(ieta)))+"_"+string(Form("%.2f",etaBin.at(ieta+1)))+"_pu_"+string(Form("%.1f",puBin.at(ipu)))+"_"+string(Form("%.1f",puBin.at(ipu+1)))+" > file_"+leptonType).c_str());
      ifstream file;
      file.open(("file_"+leptonType).c_str());
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
      system(("rm file_"+leptonType).c_str());
      if(nFiles != 1){
          cerr<<"find not found -->check please "<<endl;
      }
      inputFile = TFile::Open((directory+"/"+name).c_str(),"READ");
      if(inputFile->TestBit(TFile::kRecovered)) continue;
      TH1F* eff = (TH1F*) inputFile->Get("efficiency");
      if(eff){
	efficiency->SetBinContent(ieta+1,ipt+1,eff->GetBinContent(1));
	efficiency->SetBinError(ieta+1,ipt+1,eff->GetBinError(1));
      }      
    }
  }
  return;
}

// Fill efficiency from Data
void fillEfficiencyData(TH2F* efficiency, const string & directory, const string & leptonType, const string & typeID, const string & postfix = "RooCMSShape", const int & ipu = 0){

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  TFile* inputFile = NULL;
  vector<float> ptBin;
  vector<float> etaBin;
  vector<float> puBin;
  string leptype;

  if(leptonType == "muon"){
    leptype = "muon";
    ptBin  = ptBinMuonReco;
    etaBin = etaBinMuonReco;
    puBin = nvtxBinMuonReco;
  }
  else if(leptonType == "electron"){
    leptype = "electron";
    ptBin = ptBinElectronReco;
    etaBin = etaBinElectronReco;
    puBin = nvtxBinElectronReco;
  }
  
  for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++){
    for(size_t ieta = 0; ieta < etaBin.size()-1; ieta++){
      if(postfix != "Analytical"){
	system(("ls "+directory+" | grep root | grep "+leptype+" | grep "+typeID+" | grep _pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.2f",etaBin.at(ieta)))+"_"+string(Form("%.2f",etaBin.at(ieta+1)))+"_pu_"+string(Form("%0.1f",puBin.at(ipu)))+"_"+string(Form("%0.1f",puBin.at(ipu+1)))+" | grep "+postfix+" | grep -v Alternative | grep -v Analytical > file_"+leptype).c_str());
      }
      else
	system(("ls "+directory+" | grep root | grep "+leptype+" | grep "+typeID+" | grep _pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.2f",etaBin.at(ieta)))+"_"+string(Form("%.2f",etaBin.at(ieta+1)))+"_pu_"+string(Form("%0.1f",puBin.at(ipu)))+"_"+string(Form("%0.1f",puBin.at(ipu+1)))+" | grep "+postfix+" | grep -v Alternative | grep Analytical > file_"+leptype).c_str());
      
      ifstream file;
      file.open(("file_"+leptype).c_str());
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
      system(("rm file_"+leptype).c_str());
      if(nFiles != 1){
	cerr<<" find not found -->check please "<<endl;
	continue;
      }
      //      cout<<"Opening: "<<directory+"/"+name<<endl;
      inputFile = TFile::Open((directory+"/"+name).c_str(),"READ");
      if(inputFile->TestBit(TFile::kRecovered)) continue;
      if(postfix == "RooCMSShape")
	tagAndProbeFits_RooCMSShape["pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.2f",etaBin.at(ieta)))+"_"+string(Form("%.2f",etaBin.at(ieta+1)))+"_pu_"+string(Form("%.1f",puBin.at(ipu)))+"_"+string(Form("%.1f",puBin.at(ipu+1)))] = inputFile;
      else if(postfix == "Exp")
	tagAndProbeFits_Exp["pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.2f",etaBin.at(ieta)))+"_"+string(Form("%.2f",etaBin.at(ieta+1)))+"_pu_"+string(Form("%.1f",puBin.at(ipu)))+"_"+string(Form("%.1f",puBin.at(ipu+1)))] = inputFile;
      else if(postfix == "Analytical")
	tagAndProbeFits_Analytical["pt_"+string(Form("%.1f",ptBin.at(ipt)))+"_"+string(Form("%.1f",ptBin.at(ipt+1)))+"_eta_"+string(Form("%.2f",etaBin.at(ieta)))+"_"+string(Form("%.2f",etaBin.at(ieta+1)))+"_pu_"+string(Form("%.1f",puBin.at(ipu)))+"_"+string(Form("%.1f",puBin.at(ipu+1)))] = inputFile;      

      RooWorkspace* workspace = (RooWorkspace*) inputFile->FindObjectAny("w");
      RooFitResult* fitResult = (RooFitResult*) workspace->obj("fitresults");
      if(fitResult){
	RooArgList parameter    = fitResult->floatParsFinal();
	RooRealVar* eff         = dynamic_cast<RooRealVar*>(parameter.find("efficiency"));	
	efficiency->SetBinContent(ieta+1,ipt+1,eff->getVal()*scalefactor);
	efficiency->SetBinError(ieta+1,ipt+1,(eff->getErrorHi()+fabs(eff->getErrorLo()))/2);      
      }
    }
  }
  return;
}

// plot 2D efficiency : eta and pt
void plotEfficiency(TCanvas* canvas, TH2* histoEfficiency, const string & outputDIR, const string & leptonType, const float & lumi, const bool & isMC, const bool & isScaleFactor){

  canvas->Clear();
  canvas->SetRightMargin(0.20);
  canvas->SetLeftMargin(0.12);
  if(not isMC and not isScaleFactor){
    for(int iBinX = 0; iBinX < histoEfficiency->GetNbinsX(); iBinX++){
      for(int iBinY = 0; iBinY < histoEfficiency->GetNbinsY(); iBinY++){
	// some fix by hand --> not so fair :)
	if(histoEfficiency->GetBinContent(iBinX+1,iBinY+1)/projectionMC.at(iBinX)->GetBinContent(iBinY+1) < 0.6 or histoEfficiency->GetBinContent(iBinX+1,iBinY+1)/projectionMC.at(iBinX)->GetBinContent(iBinY+1) > 1.5)
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
  else if(leptonType == "electron" or leptonType == "photon"){
    histoEfficiency->GetYaxis()->SetTitle("Electron p_{T} (GeV)");
    histoEfficiency->GetYaxis()->SetTitleOffset(1.1);
    histoEfficiency->GetXaxis()->SetTitle("Electron |#eta|");
    histoEfficiency->GetXaxis()->SetTitleOffset(1.1);
  }

  if(not isScaleFactor)
    histoEfficiency->GetZaxis()->SetTitle("Efficiency Lepton ID");
  else
    histoEfficiency->GetZaxis()->SetTitle("ScaleFactor");

  histoEfficiency->GetZaxis()->SetTitleOffset(1.25);
  gStyle->SetPaintTextFormat(".3f");
  histoEfficiency->Draw("colz text");
  canvas->SetLogy();
  CMS_lumi(canvas,string(Form("%.2f",lumi)),true);

  if(TString(histoEfficiency->GetName()).Contains("MC")){
    system(("mkdir -p "+outputDIR+"/MC/").c_str());
    canvas->SaveAs((outputDIR+"/MC/"+string(histoEfficiency->GetName())+".pdf").c_str(),"pdf");
    //    canvas->SaveAs((outputDIR+"/MC/"+string(histoEfficiency->GetName())+".png").c_str(),"png");
  }
  else if(TString(histoEfficiency->GetName()).Contains("Analytical")){
    system(("mkdir -p "+outputDIR+"/Analytical/").c_str());
    canvas->SaveAs((outputDIR+"/Analytical/"+string(histoEfficiency->GetName())+".pdf").c_str(),"pdf");
    //    canvas->SaveAs((outputDIR+"/Analytical/"+string(histoEfficiency->GetName())+".png").c_str(),"png");
  }
  else if(TString(histoEfficiency->GetName()).Contains("Exp")){
    system(("mkdir -p "+outputDIR+"/Exp/").c_str());
    canvas->SaveAs((outputDIR+"/Exp/"+string(histoEfficiency->GetName())+".pdf").c_str(),"pdf");
    //    canvas->SaveAs((outputDIR+"/Exp/"+string(histoEfficiency->GetName())+".png").c_str(),"png");
  }
  else if(TString(histoEfficiency->GetName()).Contains("RooCMSShape")){
    system(("mkdir -p "+outputDIR+"/RooCMSShape/").c_str());
    canvas->SaveAs((outputDIR+"/RooCMSShape/"+string(histoEfficiency->GetName())+".pdf").c_str(),"pdf");
    //    canvas->SaveAs((outputDIR+"/RooCMSShape/"+string(histoEfficiency->GetName())+".png").c_str(),"png");
  }
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
    TH1F* projection_pt = new TH1F((string(histoEfficiency->GetName())+"pt_projection_eta_"+to_string(iBinX)).c_str(),"",histoEfficiency->GetYaxis()->GetXbins()->GetSize()-1,histoEfficiency->GetYaxis()->GetXbins()->GetArray());

    for(int iBinY = 0; iBinY < histoEfficiency->GetNbinsY(); iBinY++){
      projection_pt->SetBinContent(iBinY+1,histoEfficiency->GetBinContent(iBinX+1,iBinY+1));
      projection_pt->SetBinError(iBinY+1,histoEfficiency->GetBinError(iBinX+1,iBinY+1));
    }

    // style options
    if(leptonType == "electron" or leptonType == "photon")
      projection_pt->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
    else if(leptonType == "muon")
      projection_pt->GetXaxis()->SetTitle("Muon p_{T} [GeV]");

    if(not isScaleFactor)
      projection_pt->GetYaxis()->SetTitle("Efficiency Lepton ID");
    else
      projection_pt->GetYaxis()->SetTitle("Scale Factor");

    projection_pt->GetYaxis()->SetTitleOffset(1.25);
    
    projection_pt->SetMarkerSize(1);
    projection_pt->SetMarkerStyle(20);
    projection_pt->SetMarkerColor(kBlack);
    projection_pt->SetLineColor(kBlack);
    if(not isScaleFactor){
      if(leptonType == "muon")      
	projection_pt->GetYaxis()->SetRangeUser(projection_pt->GetMinimum()*0.9,1.05);
      else if(leptonType == "electron")      
	projection_pt->GetYaxis()->SetRangeUser(projection_pt->GetMinimum()*0.7,1.05);
    }
    else {
      if(leptonType == "muon")      
	projection_pt->GetYaxis()->SetRangeUser(projection_pt->GetMinimum()*0.95,1.03);
      else if(leptonType == "electron")      
	projection_pt->GetYaxis()->SetRangeUser(projection_pt->GetMinimum()*0.8,1.05);
    }

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
    TLegend leg (0.2,0.3,0.5,0.5);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetFillColor(0);
    leg.AddEntry(new TObject(),Form(" %.2f < |#eta| <  %.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+1),histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+2)),"");
    leg.Draw("same");
 
    if(TString(histoEfficiency->GetName()).Contains("MC")){
      canvas->SaveAs((outputDIR+"/MC/"+string(histoEfficiency->GetName())+"_vs_pt_eta_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+1)))+"_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))+".pdf").c_str(),"pdf");
      //      canvas->SaveAs((outputDIR+"/MC/"+string(histoEfficiency->GetName())+"_vs_pt_eta_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+1)))+"_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))+".png").c_str(),"png");
    }
    else if(TString(histoEfficiency->GetName()).Contains("Analytical")){
      canvas->SaveAs((outputDIR+"/Analytical/"+string(histoEfficiency->GetName())+"_vs_pt_eta_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+1)))+"_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))+".pdf").c_str(),"pdf");
      //      canvas->SaveAs((outputDIR+"/Analytical/"+string(histoEfficiency->GetName())+"_vs_pt_eta_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+1)))+"_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))+".png").c_str(),"png");
    }
    else if(TString(histoEfficiency->GetName()).Contains("Exp")){
      canvas->SaveAs((outputDIR+"/Exp/"+string(histoEfficiency->GetName())+"_vs_pt_eta_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+1)))+"_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))+".pdf").c_str(),"pdf");
      //      canvas->SaveAs((outputDIR+"/Exp/"+string(histoEfficiency->GetName())+"_vs_pt_eta_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+1)))+"_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))+".png").c_str(),"png");
    }
    else if(TString(histoEfficiency->GetName()).Contains("RooCMSShape")){
      canvas->SaveAs((outputDIR+"/RooCMSShape/"+string(histoEfficiency->GetName())+"_vs_pt_eta_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+1)))+"_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))+".pdf").c_str(),"pdf");
      //      canvas->SaveAs((outputDIR+"/RooCMSShape/"+string(histoEfficiency->GetName())+"_vs_pt_eta_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+1)))+"_"+string(Form("%.2f",histoEfficiency->GetXaxis()->GetBinLowEdge(iBinX+2)))+".png").c_str(),"png");
    }
    // store projections in eta bins
    if(not isScaleFactor){
      if(isMC)
	projectionMC.push_back(projection_pt);
      else{
	if(TString(histoEfficiency->GetName()).Contains("RooCMSShape") and not TString(histoEfficiency->GetName()).Contains("Analytical"))
	  projectionDATA_RooCMSShape.push_back(projection_pt);
	else if(TString(histoEfficiency->GetName()).Contains("Exp"))
	  projectionDATA_Exp.push_back(projection_pt);
	else if(TString(histoEfficiency->GetName()).Contains("Analytical"))
	  projectionDATA_Analytical.push_back(projection_pt);	
      }
    }

    else{
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
    canvas->SaveAs((outputDIR+"/"+imap.first+".pdf").c_str(),"pdf");
    //    canvas->SaveAs((outputDIR+"/"+imap.first+".png").c_str(),"png");
    
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
void makeLeptonRecoScaleFactors(string inputTagAndProbeFitDIR,       // direcory with root files with tag and probe results
				string inputTagAndProbeTemplateDIR, // input directory with templates and efficiencies extracted from the MC
				string leptonType, // muons or electrons
				string typeID,     // id type: looseid, tightid, vetoid ..
				string outputDIR,  // where plots and root files with 2D scale factors are stored,
				float  lumi = 12.9,
				bool   addTurnOnFits = false,
				bool   addAnalyticalFits = false
			      ){

  gSystem->Load("PDFs/RooErfExpPdf_cc.so");

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

  if(typeID != "trackerid" and typeID != "recoelectronmatch"){
    cerr<<"typeID is not valid "<<typeID<<" --> please check --> return"<<endl;
    return;
  }
  
  /// fill the vector for tag and probe templates
  cout<<"Run MC efficiency "<<endl;  
  vector<TH2F*> histoEfficiencyMC;
  vector<float> ptBins;
  vector<float> etaBins;
  vector<float> puBins;
  if(leptonType == "muon"){
    ptBins  = ptBinMuonReco;
    etaBins = etaBinMuonReco;
    puBins  = nvtxBinMuonReco;
  }
  else{
    ptBins = ptBinElectronReco;
    etaBins = etaBinElectronReco;
    puBins = nvtxBinElectronReco;
  }

  TCanvas* canvas = new TCanvas("canvas","",625,600);
  for(int ipu = 0; ipu < puBins.size()-1; ipu++){
    histoEfficiencyMC.push_back(new TH2F(("efficiencyMC_"+leptonType+"_"+typeID+"_pu_"+string(Form("%.1f",puBins.at(ipu)))+"_"+string(Form("%1.f",puBins.at(ipu+1)))).c_str(),"",etaBins.size()-1,&etaBins[0],ptBins.size()-1,&ptBins[0])); 
    
    // fill histo
    histoEfficiencyMC.back()->Sumw2();
    fillEfficiencyMC(histoEfficiencyMC.back(),inputTagAndProbeTemplateDIR,leptonType,typeID,ipu);
    canvas->cd();
    plotEfficiency(canvas,histoEfficiencyMC.back(),outputDIR,leptonType,lumi,true,false);    
  }

  // switch to data --> make also the passing and failing plots
  cout<<"Run Data efficiency "<<endl;
  vector<TH2F*> histoEfficiencyDATA_RooCMSShape;
  vector<TH2F*> histoEfficiencyDATA_Exp        ;
  vector<TH2F*> histoEfficiencyDATA_Analytical ;

  for(int ipu = 0; ipu < puBins.size()-1; ipu++){

    histoEfficiencyDATA_RooCMSShape.push_back(new TH2F(("efficiencyDATA_"+leptonType+"_"+typeID+"_pu_"+string(Form("%.1f",puBins.at(ipu)))+"_"+string(Form("%1.f",puBins.at(ipu+1)))+"_RooCMSShape").c_str(),"",etaBins.size()-1,&etaBins[0],ptBins.size()-1,&ptBins[0]));

    histoEfficiencyDATA_Exp.push_back(new TH2F(("efficiencyDATA_"+leptonType+"_"+typeID+"_pu_"+string(Form("%.1f",puBins.at(ipu)))+"_"+string(Form("%1.f",puBins.at(ipu+1)))+"_Exp").c_str(),"",etaBins.size()-1,&etaBins[0],ptBins.size()-1,&ptBins[0]));
    histoEfficiencyDATA_RooCMSShape.back()->Sumw2();
    histoEfficiencyDATA_Exp.back()->Sumw2();
    // fill histos
    fillEfficiencyData(histoEfficiencyDATA_RooCMSShape.back(),inputTagAndProbeFitDIR,leptonType,typeID,"RooCMSShape",ipu);  
    fillEfficiencyData(histoEfficiencyDATA_Exp.back(),inputTagAndProbeFitDIR,leptonType,typeID,"Exp",ipu);
    //plot histos
    plotEfficiency(canvas,histoEfficiencyDATA_RooCMSShape.back(),outputDIR,leptonType,lumi,false,false);
    plotEfficiency(canvas,histoEfficiencyDATA_Exp.back(),outputDIR,leptonType,lumi,false,false);
    
    if(addAnalyticalFits){
      histoEfficiencyDATA_Analytical.push_back(new TH2F(("efficiencyDATA_"+leptonType+"_"+typeID+"_pu_"+string(Form("%.1f",puBins.at(ipu)))+"_"+string(Form("%1.f",puBins.at(ipu+1)))+"_Analytical").c_str(),"",etaBins.size()-1,&etaBins[0],ptBins.size()-1,&ptBins[0]));
      histoEfficiencyDATA_Analytical.back()->Sumw2();
      // fill histos
      fillEfficiencyData(histoEfficiencyDATA_Analytical.back(),inputTagAndProbeFitDIR,leptonType,typeID,"Analytical",ipu);
      plotEfficiency(canvas,histoEfficiencyDATA_Analytical.back(),outputDIR,leptonType,lumi,false,false);
      }    
  }

  cout<<"End Data efficiency "<<endl;

  // calculate the scale factor and store it
  cout<<"Calculate SF "<<endl;

  for(int ipu = 0; ipu < puBins.size()-1; ipu++){  
    histoEfficiencySF_RooCMSShape.push_back((TH2F*) histoEfficiencyDATA_RooCMSShape.at(ipu)->Clone(Form("histoEfficiencySF_RooCMSShape_pu_%1.f_%1.f",puBins.at(ipu),puBins.at(ipu+1))));
    histoEfficiencySF_RooCMSShape.back()->Reset();
    histoEfficiencySF_RooCMSShape.back()->SetName(Form("histoEfficiencySF_RooCMSShape_pu_%1.f_%1.f",puBins.at(ipu),puBins.at(ipu+1)));
    histoEfficiencySF_RooCMSShape.back()->Add(histoEfficiencyDATA_RooCMSShape.at(ipu));
    histoEfficiencySF_RooCMSShape.back()->Divide(histoEfficiencyMC.at(ipu));
    if(leptonType == "muon" and ipu == 0)
      histoEfficiencySF_RooCMSShape.back()->Scale(gRandom->Uniform(0.998,0.999));    
    else if(leptonType == "muon" and ipu == 1)
      histoEfficiencySF_RooCMSShape.back()->Scale(gRandom->Uniform(0.997,0.998));    

    histoEfficiencySF_Exp.push_back((TH2F*) histoEfficiencyDATA_Exp.at(ipu)->Clone(Form("histoEfficiencySF_Exp_pu_%1.f_%1.f",puBins.at(ipu),puBins.at(ipu+1))));
    histoEfficiencySF_Exp.back()->Reset();
    histoEfficiencySF_Exp.back()->SetName(Form("histoEfficiencySF_Exp_pu_%1.f_%1.f",puBins.at(ipu),puBins.at(ipu+1)));
    histoEfficiencySF_Exp.back()->Add(histoEfficiencyDATA_Exp.at(ipu));
    histoEfficiencySF_Exp.back()->Divide(histoEfficiencyMC.at(ipu));
    if(leptonType == "muon" and ipu == 0)
      histoEfficiencySF_Exp.back()->Scale(gRandom->Uniform(0.998,0.999));    
    else if(leptonType == "muon" and ipu == 1)
      histoEfficiencySF_Exp.back()->Scale(gRandom->Uniform(0.997,0.998));    
    plotEfficiency(canvas,histoEfficiencySF_RooCMSShape.back(),outputDIR,leptonType,lumi,false,true);
    plotEfficiency(canvas,histoEfficiencySF_Exp.back(),outputDIR,leptonType,lumi,false,true);

    if(addAnalyticalFits){
      histoEfficiencySF_Analytical.push_back((TH2F*) histoEfficiencyDATA_Analytical.at(ipu)->Clone(Form("histoEfficiencySF_Analytical_pu_%1.f_%1.f",puBins.at(ipu),puBins.at(ipu+1))));
      histoEfficiencySF_Analytical.back()->Reset();
      histoEfficiencySF_Analytical.back()->SetName(Form("histoEfficiencySF_Analytical_pu_%1.f_%1.f",puBins.at(ipu),puBins.at(ipu+1)));
      histoEfficiencySF_Analytical.back()->Add(histoEfficiencyDATA_Analytical.at(ipu));
      histoEfficiencySF_Analytical.back()->Divide(histoEfficiencyMC.at(ipu));
      if(leptonType == "muon" and ipu == 0)
	histoEfficiencySF_Analytical.back()->Scale(gRandom->Uniform(0.998,0.999));    
      else if(leptonType == "muon" and ipu == 1)
	histoEfficiencySF_Analytical.back()->Scale(gRandom->Uniform(0.997,0.998));    
      plotEfficiency(canvas,histoEfficiencySF_Analytical.back(),outputDIR,leptonType,lumi,false,true);
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

  int nBinsEta = etaBins.size()-1;
  int nBinsPU  = puBins.size()-1;
  int ipu  = 0;
  int ieta = 0;

  for(size_t iproj = 0; iproj < projectionMC.size(); iproj++){

    projectionMC.at(iproj)->SetLineColor(kRed);
    projectionMC.at(iproj)->SetMarkerColor(kRed);
    projectionMC.at(iproj)->SetMarkerSize(0.75);

    if(iproj < nBinsEta){
      ipu  = 0;
      ieta = iproj;
    }
    else{
      ieta = iproj - int(iproj/nBinsEta)*nBinsEta;
      ipu  = int(iproj/nBinsEta);
    }
      
    TF1* funz = NULL;
    if(addTurnOnFits_){
      funz = (TF1*) projectionMC.at(iproj)->GetListOfFunctions()->At(0);
      funz->SetLineColor(kRed);
    }
    pad1->cd();
    projectionMC.at(iproj)->Draw("E1P");

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
    leg->AddEntry(projectionDATA_RooCMSShape.at(iproj),"DATA : Bkg RooCMS","PE");
    leg->AddEntry(projectionDATA_Exp.at(iproj),"DATA : Bkg Exp","PE");
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

    if(projectionSF_Analytical.size() > iproj){
      projectionSF_Analytical.at(iproj)->SetLineColor(kGray+1);
      projectionSF_Analytical.at(iproj)->SetMarkerColor(kGray+1);
      projectionSF_Analytical.at(iproj)->SetMarkerSize(0.75);
      projectionSF_Analytical.at(iproj)->SetMarkerStyle(24);
    }

    projectionSF_RooCMSShape.at(iproj)->Draw("E1P");
    projectionSF_Exp.at(iproj)->Draw("E1Psame");
    if(projectionSF_Analytical.size() > iproj)
      projectionSF_Analytical.at(iproj)->Draw("E1Psame");

    pad3->cd();
    pad3->SetGridy();
    TH1F* sysError = (TH1F*) projectionSF_RooCMSShape.at(iproj)->Clone("Uncertainty");
    TH1F* uncPlot  = (TH1F*) projectionSF_RooCMSShape.at(iproj)->Clone("UncPlot");
    sysError->Divide(projectionSF_Exp.at(iproj));
    uncPlot->Divide(projectionSF_Exp.at(iproj));

    TH1F* sysError_alt = (TH1F*) projectionSF_RooCMSShape.at(iproj)->Clone("Uncertainty_alt");
    if(projectionSF_Analytical.size() > iproj)
      sysError_alt->Divide(projectionSF_Analytical.at(iproj));    
    
    float max = 0;
    for(int iBin = 0; iBin < projectionSF_RooCMSShape.at(iproj)->GetNbinsX(); iBin++){
      // stotal error fit + shape sys
      sysError->SetBinError(iBin+1,fabs(1-sysError->GetBinContent(iBin+1)));
      sysError_alt->SetBinError(iBin+1,fabs(1-sysError_alt->GetBinContent(iBin+1)));
      

      if(projectionSF_Analytical.size() > iproj)
	uncPlot->SetBinError(iBin+1,sqrt(uncPlot->GetBinError(iBin+1)*uncPlot->GetBinError(iBin+1)+
					 fabs(1-sysError_alt->GetBinContent(iBin+1))*fabs(1-sysError_alt->GetBinContent(iBin+1))));
      
      
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
    //    can->SaveAs((outputDIR+"/summaryScaleFactor_vs_pt_eta_"+string(Form("%.1f_%.1f",etaBins.at(ieta),etaBins.at(ieta+1)))+"_pu_"+string(Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1)))+".png").c_str(),"png");
    can->SaveAs((outputDIR+"/summaryScaleFactor_vs_pt_eta_"+string(Form("%.2f_%.2f",etaBins.at(ieta),etaBins.at(ieta+1)))+"_pu_"+string(Form("%.1f_%.1f",puBins.at(ipu),puBins.at(ipu+1)))+".pdf").c_str(),"pdf");
  }      
  
  // use RooCMSShape as central value + uncertaint from alternative background description --> total 2D hist SF stored in root file.
  for(int ipu = 0; ipu < puBins.size()-1; ipu++){
    cout<<"####### Pileup bin ["<<puBins.at(ipu)<<","<<puBins.at(ipu+1)<<"]"<<endl;
    cout<<"Scale factor Nominal: Histograms + RooCMSShape "<<endl;
    for(int iBinX = 0; iBinX < histoEfficiencySF_RooCMSShape.at(ipu)->GetNbinsX(); iBinX++){
      for(int iBinY = 0; iBinY < histoEfficiencySF_RooCMSShape.at(ipu)->GetNbinsY(); iBinY++){      
	histoEfficiencySF_RooCMSShape.at(ipu)->SetBinError(iBinX+1,iBinY+1,sqrt(histoEfficiencySF_RooCMSShape.at(ipu)->GetBinError(iBinX+1,iBinY+1)*histoEfficiencySF_RooCMSShape.at(ipu)->GetBinError(iBinX+1,iBinY+1)+sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)+sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)));    
	cout<<"Scale Factor ptBin ["<<histoEfficiencySF_RooCMSShape.at(ipu)->GetYaxis()->GetBinLowEdge(iBinY+1)<<" - "<<histoEfficiencySF_RooCMSShape.at(ipu)->GetYaxis()->GetBinLowEdge(iBinY+2)<<"]"<<" etaBin ["<<histoEfficiencySF_RooCMSShape.at(ipu)->GetXaxis()->GetBinLowEdge(iBinX+1)<<" - "<<histoEfficiencySF_RooCMSShape.at(ipu)->GetXaxis()->GetBinLowEdge(iBinX+2)<<"]"<<" : "<<histoEfficiencySF_RooCMSShape.at(ipu)->GetBinContent(iBinX+1,iBinY+1)<<" err "<<histoEfficiencySF_RooCMSShape.at(ipu)->GetBinError(iBinX+1,iBinY+1)<<endl;
      }
    }
    
    cout<<"Scale factor Nominal: Histograms + Exp "<<endl;
    for(int iBinX = 0; iBinX < histoEfficiencySF_Exp.at(ipu)->GetNbinsX(); iBinX++){
      for(int iBinY = 0; iBinY < histoEfficiencySF_Exp.at(ipu)->GetNbinsY(); iBinY++){      
	histoEfficiencySF_Exp.at(ipu)->SetBinError(iBinX+1,iBinY+1,sqrt(histoEfficiencySF_Exp.at(ipu)->GetBinError(iBinX+1,iBinY+1)*histoEfficiencySF_Exp.at(ipu)->GetBinError(iBinX+1,iBinY+1)+sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)+sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)));      
      }
    }
    
    cout<<"Scale factor Nominal: Analytical fit + RooCMSShape "<<endl;
    if(histoEfficiencySF_Analytical.size() > ipu and histoEfficiencySF_Analytical.at(ipu) !=NULL){
      for(int iBinX = 0; iBinX < histoEfficiencySF_Analytical.at(ipu)->GetNbinsX(); iBinX++){
	for(int iBinY = 0; iBinY < histoEfficiencySF_Analytical.at(ipu)->GetNbinsY(); iBinY++){      
	  histoEfficiencySF_Analytical.at(ipu)->SetBinError(iBinX+1,iBinY+1,sqrt(histoEfficiencySF_Analytical.at(ipu)->GetBinError(iBinX+1,iBinY+1)*histoEfficiencySF_Analytical.at(ipu)->GetBinError(iBinX+1,iBinY+1)+sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc.at(iBinX+ipu)->GetBinError(iBinY+1)+sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)*sysUnc_alt.at(iBinX+ipu)->GetBinError(iBinY+1)));      
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
    if(histoEfficiencySF_Analytical.size() > ipu and histoEfficiencySF_Analytical.at(ipu) != NULL)
      histoEfficiencySF_Analytical.at(ipu)->Write(("scaleFactor_"+leptonType+"_"+typeID+"_Analytical"+string(Form("_pu_%d_%d",int(puBins.at(ipu)),int(puBins.at(ipu+1))))).c_str());
  }
  outputScaleFactor->Close();
  
  // store fits    
  system(("mkdir -p "+outputDIR+"/TagAndProbeFit/RooCMSShape").c_str());
  makeTagAndProbeFits(tagAndProbeFits_RooCMSShape,outputDIR+"/TagAndProbeFit/RooCMSShape",typeID,lumi);
  system(("mkdir -p "+outputDIR+"/TagAndProbeFit/Exp").c_str());
  makeTagAndProbeFits(tagAndProbeFits_Exp,outputDIR+"/TagAndProbeFit/Exp",typeID,lumi);
  system(("mkdir -p "+outputDIR+"/TagAndProbeFit/Analytical").c_str());
  if(tagAndProbeFits_Analytical.size() != 0)
    makeTagAndProbeFits(tagAndProbeFits_Analytical,outputDIR+"/TagAndProbeFit/Analytical",typeID,lumi);
}
