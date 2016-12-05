#include "../CMS_lumi.h"
#include "TnPBinning.h"
#include "PDFs/RooErfExpPdf.h"

// bin in abs eta instead of eta
static bool  isabseta     = true;
static float luminosity   = 35.9;
// reco efficiency or id
static bool  isRecoEff    = false;
// normalize plots to a.u.
static bool  plotRescaled = true;
// perform a continous fit as well
static bool  simFit       = true;
// which kind of PDF you want to use
static bool  usePolynomial  = false;
static bool  useCBShape     = false;
static bool  useRooCMSShape = false;
static bool  useErfExpPdf   = false;
// RooCMSShape for 10-20 pt, or high pt, RooErfExp in the middle pt range
static bool  useMixedModel  = true;
static float ptThrehsoldForRebin = 70.;

/// make the templates
void maketemplate(const string & inputDIR, 
		  const string & leptonType, 
		  const string & outputDIR, 
		  const string & idName,
		  const float  & ptMin,
		  const float  & ptMax,
		  const float  & etaMin,
		  const float  & etaMax,
		  const float  & nvtxMin,
		  const float  & nvtxMax,
		  bool  smooth = true,
		  int   nbins = 60, 
		  float xmin  = 60, 
		  float xmax  = 120) {


  TChain* inputMC = NULL;
  if(leptonType == "muon")
    inputMC = new TChain("muontnptree/fitter_tree");
  else if(leptonType == "electron" and not isRecoEff)
    inputMC = new TChain("electrontnptree/fitter_tree");
  else if(leptonType == "photon" or (leptonType == "electron" and isRecoEff))
    inputMC = new TChain("photontnptree/fitter_tree");

  inputMC->Add((inputDIR+"/*root").c_str());

  map<string,float> crossSection;
  crossSection["100to200"] = 148.;
  crossSection["200to400"] = 40.94;
  crossSection["400to600"] = 5.49;
  crossSection["600toInf"] = 2.19;

  // pileup-re-weight from external file
  TFile* pufile = new TFile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_35p9fb.root");
  TH1* puhist = (TH1*)pufile->Get("puhist");
  
  TTreeReader reader(inputMC);
  TTreeReaderValue<float> mass (reader, "mass");
  TTreeReaderValue<float> eta  (reader, "eta");
  TTreeReaderValue<float> phi  (reader, "phi");
  TTreeReaderValue<float> pt   (reader, "pt");
  TTreeReaderValue<float> nvtx (reader, "nvtx");
  TTreeReaderValue<float> wgt  (reader, "wgt");
  TTreeReaderValue<int>   mcTrue (reader, "mcTrue");
  TTreeReaderValue<int>   id   (reader, idName.c_str());
  string idtrigger = "";
  string trackerid = "";
  string standaloneid = "";
  if(leptonType == "muon"){
    trackerid    = "trackerid";
    standaloneid = "standaloneid";
    idtrigger    = idName;
  }
  else if((leptonType == "electron" and isRecoEff) or leptonType == "photon"){
    standaloneid = "recoelectronmatch";
    trackerid = "recoelectronmatch";
    idtrigger = idName;
  }
  else{
    trackerid = idName ;
    standaloneid = idName;
    if(leptonType == "electron")
      idtrigger = "hltsafeid";
    else
      idtrigger = idName;
	
  }

  TTreeReaderValue<int>   isStandAloneMuon (reader, standaloneid.c_str());
  TTreeReaderValue<int>   isTrackerMuon    (reader, trackerid.c_str());
  TTreeReaderValue<int>   idtrig           (reader,idtrigger.c_str());    
  TTreeReaderValue<float> mcMass (reader, "mcMass");
  
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
  string currentTree ;
  
  while(reader.Next()){
   
    // probe definition
    if(leptonType == "muon"){
      if(isRecoEff == true and not *isStandAloneMuon) continue;
      else if(isRecoEff == false and not (*isStandAloneMuon || *isTrackerMuon)) continue;
    }

    // in case HT binned samples
    float lumiwgt = 1.;
    float xsec = 1;
    for(auto entry : crossSection){      
      if(TString(reader.GetTree()->GetCurrentFile()->GetName()).Contains(entry.first)){
	lumiwgt = luminosity*1000*entry.second;      
	xsec = entry.second;
      }
    }
    
    if(currentTree != reader.GetTree()->GetCurrentFile()->GetName()){
      cout<<"current tree "<<reader.GetTree()->GetCurrentFile()->GetName()<<" xsec "<<xsec<<endl;
      currentTree = reader.GetTree()->GetCurrentFile()->GetName();
    }

    Float_t puwgt = 0.;
    if (*nvtx <= 40) puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
    else puwgt = 1;

    // check the probe lepton information --> if it falls inside the bins
    if(*pt < ptMin or *pt > ptMax) continue;
    if(isabseta and (fabs(*eta) < etaMin or fabs(*eta) > etaMax)) continue;
    else if(not isabseta and (*eta < etaMin or *eta > etaMax)) continue;
    if(*nvtx <= nvtxMin or *nvtx > nvtxMax) continue;
    // if not matched to a genLepton skip --> we want to extract the true templateds

    if(not *mcTrue) continue;

    if (*id == 0 or *idtrig == 0) hfail.Fill(*mass, lumiwgt*puwgt*(*wgt)/wgtsum);
    if (*id >  0 and *idtrig > 0) hpass.Fill(*mass, lumiwgt*puwgt*(*wgt)/wgtsum);
    if (*id >  0 and *idtrig > 0) hp.Fill(*mass, lumiwgt*puwgt*(*wgt)/wgtsum);
    ha.Fill(*mass, lumiwgt*puwgt*(*wgt)/wgtsum);
  }

  for (int j = 1; j <= nbins; j++) {
    if (hpass.GetBinContent(j) < 0.) hpass.SetBinContent(j, 0.0);
    if (hfail.GetBinContent(j) < 0.) hfail.SetBinContent(j, 0.0);
  }
  
  // efficiency value with Binomial uncertainty
  hr.Divide(&hp, &ha, 1.0, 1.0, "B");
  // smooth templates two times
  if(smooth){// amcnlo as lower stat
    if(ptMin >= ptThrehsoldForRebin){
      hpass.Rebin(2);
      hfail.Rebin(2);
    }
  }
  
  // Build the RooHistPdf  
  string fileName = "template_TaP_"+leptonType+"_"+idName+"_pt_"+string(Form("%.1f",ptMin))+"_"+string(Form("%.1f",ptMax))+"_eta_"+string(Form("%.1f",etaMin))+"_"+string(Form("%.1f",etaMax))+"_pu_"+string(Form("%.1f",nvtxMin))+"_"+string(Form("%.1f",nvtxMax));
  TFile outfile((outputDIR+"/"+fileName+".root").c_str(), "RECREATE");
  hr.Write("efficiency");

  /// lineshape observable
  RooRealVar m("mass", "", xmin, xmax);
  // passing and failing sample binnned data
  RooDataHist dhpass("dhpass", "", RooArgList(m), RooFit::Import(hpass), 0);
  RooDataHist dhfail("dhfail", "", RooArgList(m), RooFit::Import(hfail), 0);
  // passing and failing histogram pdf
  RooHistPdf signalPassMC("signalPassMC", "", RooArgSet(m), dhpass);
  RooHistPdf signalFailMC("signalFailMC", "", RooArgSet(m), dhfail);
  
  // Resolution for passing sample -> simple CB
  RooRealVar  zPeakP1    ("zPeakP1","",91.,85.,101);
  RooRealVar  zWidthP1   ("zWidthP1","",2.1,0.,3.4);
  RooRealVar  sigmaP1    ("sigmaP1","",2.,0.,10.);
  RooVoigtian zPeakPass  ("zPeakPass","",m,zPeakP1,zWidthP1,sigmaP1);
  RooRealVar  meanP2     ("meanP2","",62.,60.,88.);
  RooRealVar  sigmaP2    ("sigmaP2","",2.,1.,10.);
  RooRealVar  alphaP2    ("alphaP2","",0.1,-10,10.);
  RooRealVar  nP2        ("nP2","",1.,0.,10.);
  RooCBShape  zTailPass  ("zTailPass","",m,meanP2,sigmaP2,alphaP2,nP2);
  RooRealVar  fracPass   ("fracPass","",0.9,0.,1.);
  //RooAddPdf   signalPassAna ("signalPassAna","",RooArgList(zPeakPass,zTailPass),fracPass);
  RooAddPdf   signalPassAna ("signalPass","",RooArgList(zPeakPass,zTailPass),fracPass); // we will trust MC signal for passing sample in the analytical fit
  // create a CB for the resolution
  RooRealVar zPeakF1     ("zPeakF1","",91.,85.,101);
  RooRealVar zWidthF1    ("zWidthF1","",2.1,0.,3.4);
  RooRealVar sigmaF1     ("sigmaF1","",2.,0.,10.);
  RooVoigtian zPeakFail  ("zPeakFail","",m,zPeakF1,zWidthF1,sigmaF1); 
  // poly
  RooRealVar   a0        ("a0","a0",0.5,-10.,10.) ;
  RooRealVar   a1        ("a1","a1",0.5,-10.,10.) ;
  RooChebychev zPol      ("zRooChebychev","",m,RooArgList(a0,a1));
  // CB
  RooRealVar  meanF2     ("meanF2","",65.,60.,85.);
  RooRealVar  sigmaF2    ("sigmaF2","",2.,1.,20.);
  RooRealVar  alphaF2    ("alphaF2","",0.1,-10,10.);
  RooRealVar  nF2        ("nF2","",1.,-10.,10.);
  RooCBShape  zCB        ("zCrystalBall","",m,meanF2,sigmaF2,alphaF2,nF2);  
  // RooCMSShape
  RooRealVar alphaFail ("alphaFail","alphaFail",125,30,300);
  RooRealVar betaFail  ("betaFail","betaFail",0.01,-1,1);
  RooRealVar gammaFail ("gammaFail","gammaFail",0.05,-0.5,0.5);
  RooRealVar peakFail  ("peakFail","peakFail",91.2,60,120);
  RooCMSShape zRooCMSShape ("zRooCMSShape","",m,alphaFail,betaFail,gammaFail,peakFail);

  // RooErfExp
  RooRealVar   cExp ("cExp","",   -0.005,-1,0.);
  RooRealVar   offset ("offset","",90,60,120);
  RooRealVar   width ("width","",  10,0.,25);
  RooErfExpPdf erfExp ("zRooErfExp","",m,cExp,offset,width);

  RooRealVar  fracFail   ("fracFail","",0.9,0.,1.);
  RooAddPdf*  signalFailAna = NULL;  // failing signal will be convoluted with a gaussian.
  if(usePolynomial)
    signalFailAna = new RooAddPdf("signalFailAna","",RooArgList(zPeakFail,zPol),fracFail);
  else if(useCBShape)
    signalFailAna = new RooAddPdf("signalFailAna","",RooArgList(zPeakFail,zCB),fracFail);
  else if(useRooCMSShape or (useMixedModel and (ptMin < 20 )))
    signalFailAna = new RooAddPdf("signalFailAna","",RooArgList(zPeakFail,zRooCMSShape),fracFail);
  else if(useErfExpPdf or (useMixedModel and  (ptMin >= 20 )))
    signalFailAna = new RooAddPdf("signalFailAna","",RooArgList(zPeakFail,erfExp),fracFail);

  /// prepare for simultaneous fit --> define efficiency, total pass and fial, as well as extended pdfs
  RooRealVar efficiency ("efficiency", "Efficiency", dhpass.sumEntries()/(dhfail.sumEntries()+dhpass.sumEntries()),0, 1); 
  RooRealVar    nTot    ("nTot","",dhfail.sumEntries()+dhpass.sumEntries(),0,(dhfail.sumEntries()+dhpass.sumEntries())*100);  
  RooFormulaVar nPass   ("nPass","@0*@1",RooArgList(efficiency,nTot));
  RooFormulaVar nFail   ("nFail","(1-@0)*@1",RooArgList(efficiency,nTot));  
  RooRealVar    nSignal ("nSignal","",dhpass.sumEntries(),0,dhpass.sumEntries()*100);
  RooRealVar    nBack   ("nBack","",dhfail.sumEntries(),0,dhfail.sumEntries()*100);

  RooExtendPdf  *pdfPass, *pdfFail;
  if(simFit){
    pdfPass = new RooExtendPdf("pdfPassAna","",signalPassAna,nPass);
    pdfFail = new RooExtendPdf("pdfFailAna","",*signalFailAna,nFail);
  }
  else{
    pdfPass = new RooExtendPdf("pdfPassAna","",signalPassAna,nSignal);
    pdfFail = new RooExtendPdf("pdfFailAna","",*signalFailAna,nBack);
  }

  //combined category and pdf
  RooCategory category("category","category") ;
  category.defineType("Pass");
  category.defineType("Fail");

  RooSimultaneous simPdf("simPdfAna","simultaneous pdf",category) ;
  simPdf.addPdf(*pdfPass,"Pass");
  simPdf.addPdf(*pdfFail,"Fail");

  // combined dataset
  std::map <string,RooDataHist*> histoMap;
  histoMap ["Pass"] = &dhpass;
  histoMap ["Fail"] = &dhfail;
  RooDataHist combData = RooDataHist( "combData", "combined data", RooArgList(m),category,histoMap);

  if(simFit){
    RooAbsReal* simNLL = simPdf.createNLL(combData,RooFit::Extended(true));
    RooMinuit minuit(*simNLL);
    cout<<" Minimize with Migrad "<<endl;
    minuit.migrad();
    cout<<" Minimize with Hesse "<<endl;
    minuit.hesse();
    cout<<" Final fit "<<endl;
    simPdf.fitTo(combData,RooFit::Extended(true),RooFit::SumW2Error(kTRUE));
  }
  else{
    pdfPass->fitTo(dhpass,RooFit::SumW2Error(kTRUE));
    pdfFail->fitTo(dhfail,RooFit::SumW2Error(kTRUE));    
  }
  RooWorkspace w("w", "");
  // RooDataHist used for efficiency
  w.import(dhpass);
  w.import(dhfail);
  w.import(signalPassMC);
  w.import(signalFailMC);
  // fix pdf parameters before storing the workspace
  RooArgList passParam = RooArgList(*(signalPassAna.getParameters(dhpass)));
  RooArgList failParam = RooArgList(*(signalFailAna->getParameters(dhfail)));
  for(int i = 0; i < passParam.getSize(); i++)
    ((RooRealVar*) passParam.at(i))->setConstant(kTRUE);  
  for(int i = 0; i < failParam.getSize(); i++)
    ((RooRealVar*) failParam.at(i))->setConstant(kTRUE);
  

  // RooHistPdf used for tag and probe fit
  if(simFit)
    w.import(simPdf);
  else{
    w.import(signalPassAna);
    w.import(*signalFailAna);
  }
  w.Write();

  cout << "Efficiency -- pT [" << ptMin << ", " << ptMax << "], eta [" << etaMin << ", " << etaMax << "]   :  pu [" << nvtxMin << " ,  "<<nvtxMax<<" ] " << hr.GetBinContent(1) << " +/- " << hr.GetBinError(1) << endl;

  // plotting part
  TCanvas* c1 = new TCanvas("c1","",600,650);
  c1->cd();
  RooPlot* frame = m.frame();
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("mass (GeV)");
  frame->GetYaxis()->SetTitle("a.u.");
  frame->GetYaxis()->SetTitleOffset(1.1);
  if(plotRescaled){    
    dhpass.plotOn(frame,RooFit::LineColor(kRed),RooFit::MarkerColor(kRed),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::DataError(RooAbsData::SumW2),RooFit::DrawOption("EP"),RooFit::Rescale(1./hpass.Integral()));
    dhfail.plotOn(frame,RooFit::LineColor(kBlack),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::DataError(RooAbsData::SumW2),RooFit::DrawOption("EP"),RooFit::Rescale(1./hfail.Integral()));
    
    hpass.Scale(1./hpass.Integral());
    hfail.Scale(1./hfail.Integral());
    frame->SetMaximum(max(hpass.GetMaximum(),hfail.GetMaximum())*1.2);
    pdfPass->plotOn(frame,RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::DrawOption("L"),RooFit::Normalization(1,RooAbsReal::NumEvent));
    pdfFail->plotOn(frame,RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::DrawOption("L"),RooFit::Normalization(1,RooAbsReal::NumEvent));
  }
  else{
    dhpass.plotOn(frame,RooFit::LineColor(kRed),RooFit::MarkerColor(kRed),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::DataError(RooAbsData::SumW2),RooFit::DrawOption("EP"),RooFit::Binning(40,m.getMin(),m.getMax()));
    dhfail.plotOn(frame,RooFit::LineColor(kBlack),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::DataError(RooAbsData::SumW2),RooFit::DrawOption("EP"),RooFit::Binning(40,m.getMin(),m.getMax()));
    pdfPass->plotOn(frame,RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::DrawOption("L"));
    pdfFail->plotOn(frame,RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::DrawOption("L"));
  }
  frame->Draw();  
  TLegend* leg = new TLegend(0.2,0.75,0.5,0.92);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  hpass.SetLineColor(kRed);
  hfail.SetLineColor(kBlack);
  hfail.SetLineWidth(2);
  hpass.SetLineWidth(2);
  leg->AddEntry(&hpass,"Pass Probes","L");
  leg->AddEntry(&hfail,"Fail Probes","L");
  leg->Draw("same");
  CMS_lumi(c1,"",true);
  c1->SaveAs((outputDIR+"/"+fileName+".png").c_str(),"png");
  c1->SaveAs((outputDIR+"/"+fileName+".pdf").c_str(),"pdf");

  if(frame) delete frame;
  if(leg) delete leg;
  if(c1) delete c1;
}

void makeTnPTemplates(string inputDIR, string leptonType, string outputDIR, bool isRecoEfficiency = false) {

  TGaxis::SetMaxDigits(3);

  gSystem->Load("PDFs/RooErfExpPdf_cc.so");

  isRecoEff = isRecoEfficiency;

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  if(leptonType != "muon" and leptonType!= "electron" and leptonType!= "photon"){
    cerr<<"Problem with the lepton type : muon or electron or photons --> return "<<endl;
    return;
  }
    

  if(not isRecoEff){ //sequence for lepon id template maker

    if(leptonType == "muon"){
      for(size_t ipt = 0; ipt < ptBinMuon.size()-1; ipt++){
	for(size_t ieta = 0; ieta < etaBinMuon.size()-1; ieta++){
	  for(size_t invtx = 0; invtx < nvtxBinMuon.size()-1; invtx++){
	    string idName = "looseid";
	    maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinMuon.at(ipt),ptBinMuon.at(ipt+1),etaBinMuon.at(ieta),etaBinMuon.at(ieta+1),nvtxBinMuon.at(invtx),nvtxBinMuon.at(invtx+1),true);    
	    idName = "tightid";
	    maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinMuon.at(ipt),ptBinMuon.at(ipt+1),etaBinMuon.at(ieta),etaBinMuon.at(ieta+1),nvtxBinMuon.at(invtx),nvtxBinMuon.at(invtx+1),true);    
	  }
	}
      }
    }
    else if(leptonType == "electron"){
      for(size_t ipt = 0; ipt < ptBinElectron.size()-1; ipt++){
	for(size_t ieta = 0; ieta < etaBinElectron.size()-1; ieta++){
	  for(size_t invtx = 0; invtx < nvtxBinElectron.size()-1; invtx++){
	    string idName = "vetoid";
	    maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinElectron.at(ipt),ptBinElectron.at(ipt+1),etaBinElectron.at(ieta),etaBinElectron.at(ieta+1),nvtxBinMuon.at(invtx),nvtxBinMuon.at(invtx+1));    	
	    idName = "tightid";
	    maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinElectron.at(ipt),ptBinElectron.at(ipt+1),etaBinElectron.at(ieta),etaBinElectron.at(ieta+1),nvtxBinMuon.at(invtx),nvtxBinMuon.at(invtx+1));    	
	  }
	} 
      }
    }
    else if (leptonType == "photon"){
      for(size_t ipt = 0; ipt < ptBinPhoton.size()-1; ipt++){
	for(size_t ieta = 0; ieta < etaBinPhoton.size()-1; ieta++){
	  for(size_t invtx = 0; invtx < nvtxBinPhoton.size()-1; invtx++){
	    string idName = "mediumid";
	    maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinPhoton.at(ipt),ptBinPhoton.at(ipt+1),etaBinPhoton.at(ieta),etaBinPhoton.at(ieta+1),nvtxBinPhoton.at(invtx),nvtxBinPhoton.at(invtx+1));
	  }
	}     
      }
    }
  }
  // template for reco efficiency
  else{

    if (leptonType == "muon"){
      for(size_t ipt = 0; ipt < ptBinMuonReco.size()-1; ipt++){
	for(size_t ieta = 0; ieta < etaBinMuonReco.size()-1; ieta++){
	  for(size_t invtx = 0; invtx < nvtxBinMuonReco.size()-1; invtx++){
	    string idName = "trackerid";
	    maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinMuonReco.at(ipt),ptBinMuonReco.at(ipt+1),etaBinMuonReco.at(ieta),etaBinMuonReco.at(ieta+1),nvtxBinMuonReco.at(invtx),nvtxBinMuonReco.at(invtx+1),true);    
	  }
	}
      }
    }
    else if(leptonType == "electron"){
      for(size_t ipt = 0; ipt < ptBinElectronReco.size()-1; ipt++){
	for(size_t ieta = 0; ieta < etaBinElectronReco.size()-1; ieta++){
	  for(size_t invtx = 0; invtx < nvtxBinElectronReco.size()-1; invtx++){
	    string idName = "recoelectronmatch";
	    leptonType = "photon";
	    maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinElectronReco.at(ipt),ptBinElectronReco.at(ipt+1),etaBinElectronReco.at(ieta),etaBinElectronReco.at(ieta+1),nvtxBinElectronReco.at(invtx),nvtxBinElectronReco.at(invtx+1));    
	  }
	}
      }
    }
  }
}
