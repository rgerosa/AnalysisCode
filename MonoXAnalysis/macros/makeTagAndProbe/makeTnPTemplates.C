#include "../CMS_lumi.h"
#include "TnPBinning.h"
#include "PDFs/RooErfExpPdf.h"

// bin in abs eta instead of eta
static bool  isabseta     = false;
// luminosity value
static float luminosity   = 35.9;
// if summer or spring 16 for pu-weight
static bool  isSUMMER16 = true;
// reco efficiency or id
static bool  isRecoEff    = false;
// normalize plots to a.u.
static bool  plotRescaled = true;
// apply hltsafe id for electrons
static bool  useHLTSafeID = false;
// perform a continous fit as well
static bool  simFit       = true;
// which kind of PDF you want to use for simulate signal/background componend 
static bool  usePolynomial  = false;
static bool  useCBShape     = false;
static bool  useRooCMSShape = false;
static bool  useErfExpPdf   = false;
// RooCMSShape for 10-20 pt, or high pt, RooErfExp in the middle pt range
static bool  useMixedModel  = true;
// make less bins due to low MC stat
static float ptThrehsoldForRebin = 100.;

/// make the templates
void maketemplate(const string & inputDIR, 
		  const string & leptonType, 
		  const string & outputDIR, 
		  const string & idName,
		  const vector<float>  & ptBin,
		  const vector<float>  & etaBin,
		  const vector<float>  & nvtxBin,
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

  // in case of exclusive HT bins --> set xsections
  map<string,float> crossSection;
  crossSection["100to200"] = 148.;
  crossSection["200to400"] = 40.94;
  crossSection["400to600"] = 5.49;
  crossSection["600toInf"] = 2.19;

  // pileup-re-weight from external file
  TFile* pufile = NULL;
  if(not isSUMMER16)
    pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_35p9fb.root");
  else
    pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt_36.40_summer16.root");
  TH1* puhist = (TH1*)pufile->Get("puhist");

  // set branches
  TTreeReader reader(inputMC);
  TTreeReaderValue<float> mass (reader, "mass");
  TTreeReaderValue<float> eta  (reader, "eta");
  TTreeReaderValue<float> phi  (reader, "phi");
  TTreeReaderValue<float> pt   (reader, "pt");
  TTreeReaderValue<float> nvtx (reader, "nvtx");
  TTreeReaderValue<float> wgt  (reader, "wgt");
  TTreeReaderValue<int>   mcTrue (reader, "mcTrue");
  TTreeReaderValue<int>   id   (reader, idName.c_str());

  // some info usefult to measure reco-efficiency
  string idtrigger    = "";
  string trackerid    = "";
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
    trackerid     = idName ;
    standaloneid  = idName;
    if(leptonType == "electron" and useHLTSafeID)
      idtrigger   = "hltsafeid";
    else
      idtrigger = idName;    
  }

  TTreeReaderValue<int>   isStandAloneMuon (reader, standaloneid.c_str());
  TTreeReaderValue<int>   isTrackerMuon    (reader, trackerid.c_str());
  TTreeReaderValue<int>   idtrig           (reader,idtrigger.c_str());    
  TTreeReaderValue<float> mcMass (reader, "mcMass");
  
  // histograms
  vector<TH1F> hpass;
  vector<TH1F> hfail;
  vector<TH1F> hp;
  vector<TH1F> ha;
  vector<TH1F> hr;  
  vector<float> ptMin;

  for(size_t ptbin = 0 ; ptbin <  ptBin.size()-1; ptbin++){
    for(size_t etabin = 0 ; etabin <  etaBin.size()-1; etabin++){
      for(size_t nvtxbin = 0 ; nvtxbin <  nvtxBin.size()-1; nvtxbin++){
	ptMin.push_back(ptBin.at(ptbin));
	hpass.push_back(TH1F(Form("hpass_pt_%.1f_%.1f_eta_%.2f_%.2f_pu_%.1f_%.1f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1),
				  nvtxBin.at(nvtxbin),nvtxBin.at(nvtxbin+1)), "", nbins, xmin, xmax));	
	hfail.push_back(TH1F(Form("hfail_pt_%.1f_%.1f_eta_%.2f_%.2f_pu_%.1f_%.1f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1),
				  nvtxBin.at(nvtxbin),nvtxBin.at(nvtxbin+1)), "", nbins, xmin, xmax));
	hp.push_back(TH1F(Form("hp_pt_%.1f_%.1f_eta_%.2f_%.2f_pu_%.1f_%.1f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1),
			       nvtxBin.at(nvtxbin),nvtxBin.at(nvtxbin+1)), "",1,xmin, xmax));
	ha.push_back(TH1F(Form("ha_pt_%.1f_%.1f_eta_%.2f_%.2f_pu_%.1f_%.1f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1),
			       nvtxBin.at(nvtxbin),nvtxBin.at(nvtxbin+1)), "",1,xmin, xmax));
	hr.push_back(TH1F(Form("hr_pt_%.1f_%.1f_eta_%.2f_%.2f_pu_%.1f_%.1f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1),
			       nvtxBin.at(nvtxbin),nvtxBin.at(nvtxbin+1)), "",1,xmin, xmax));
	hpass.back().Sumw2();
	hfail.back().Sumw2();
	hp.back().Sumw2();
	ha.back().Sumw2();
	hr.back().Sumw2();
      }
    }
  }

  // weight sum
  cout<<"Calculate sum of weights "<<endl;
  double wgtsum = 0;
  while(reader.Next()){
    wgtsum += *wgt;
  }
  reader.SetEntry(0);

  // loop on the event
  string currentTree ;  
  cout<<"Start looping on the event "<<endl;
  while(reader.Next()){
   
    // probe definition
    if(leptonType == "muon"){ // apply probe selection for muons
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

    // read pu-wgt
    Float_t puwgt = 0.;
    if (*nvtx <= 40) puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
    else puwgt = 1;

    // find the right bin
    size_t ptbin = 0;
    size_t etabin = 0;
    size_t nvtxbin = 0;

    for( ; ptbin <  ptBin.size()-1; ptbin++){
      if( *pt <= ptBin.at(ptbin) or *pt > ptBin.at(ptbin+1)) continue;
      else break;
    }
    for( ; etabin <  etaBin.size()-1; etabin++){
      if(isabseta and ( fabs(*eta) <= etaBin.at(etabin) or fabs(*eta) > etaBin.at(etabin+1))) continue;
      else if(not isabseta and (*eta <= etaBin.at(etabin) or *eta > etaBin.at(etabin+1))) continue;
      else break;
    }
    for( ; nvtxbin <  nvtxBin.size()-1; nvtxbin++){
      if( *nvtx <= nvtxBin.at(nvtxbin) or *nvtx > nvtxBin.at(nvtxbin+1)) continue;
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
    if(*nvtx >= nvtxBin.at(nvtxBin.size()-1))
      continue;
    if(*nvtx < nvtxBin.at(0))
      continue;

    TString bin = Form("pt_%.1f_%.1f_eta_%.2f_%.2f_pu_%.1f_%.1f",ptBin.at(ptbin),ptBin.at(ptbin+1),etaBin.at(etabin),etaBin.at(etabin+1),nvtxBin.at(nvtxbin),nvtxBin.at(nvtxbin+1));
    size_t binHisto = 0;
    for( ; binHisto < hpass.size()-1; binHisto++){
      if(TString(hpass.at(binHisto).GetName()).Contains(bin))
    	break;
    }
      
    // if not matched to a genLepton skip --> we want to extract the true templateds
    if(not *mcTrue) continue;

    // fill passing and failing samples
    if (*id == 0) hfail.at(binHisto).Fill(*mass, lumiwgt*puwgt*(*wgt)/wgtsum);
    if (*id >  0) hpass.at(binHisto).Fill(*mass, lumiwgt*puwgt*(*wgt)/wgtsum);
    if (*id >  0) hp.at(binHisto).Fill(*mass, lumiwgt*puwgt*(*wgt)/wgtsum);
    ha.at(binHisto).Fill(*mass, lumiwgt*puwgt*(*wgt)/wgtsum);
  }

  // fix in case of problems
  for(size_t size = 0; size < hpass.size(); size++){
    for (int j = 1; j <= nbins; j++) {
      if (hpass.at(size).GetBinContent(j) < 0.) hpass.at(size).SetBinContent(j, 0.0);
      if (hfail.at(size).GetBinContent(j) < 0.) hfail.at(size).SetBinContent(j, 0.0);
    }
  
    // efficiency value with Binomial uncertainty --> matching efficiency for MC
    hr.at(size).Divide(&hp.at(size), &ha.at(size), 1.0, 1.0, "B");
    // smooth templates two times
    if(smooth){// amcnlo as lower stat
      hpass.at(size).Rebin(2);
      hfail.at(size).Rebin(2);
    }
      
    // Build the RooHistPdf
    TString name (hpass.at(size).GetName());
    TString fileName = Form("template_TaP_%s_%s",leptonType.c_str(),idName.c_str())+name.ReplaceAll("hpass_","");
    cout << "Efficiency -- "<<fileName<<" = "<<hr.at(size).GetBinContent(1)<<" +/- " <<hr.at(size).GetBinError(1)<<endl;
    TFile outfile((outputDIR+"/"+string(fileName)+".root").c_str(), "RECREATE");
    hr.at(size).Write("efficiency");

    /// lineshape observable
    RooRealVar m("mass", "", xmin, xmax);
    // passing and failing sample binnned data
    RooDataHist dhpass("dhpass", "", RooArgList(m), RooFit::Import(hpass.at(size)), 0);
    RooDataHist dhfail("dhfail", "", RooArgList(m), RooFit::Import(hfail.at(size)), 0);
    // passing and failing histogram pdf
    RooHistPdf signalPassMC("signalPassMC", "", RooArgSet(m), dhpass);
    RooHistPdf signalFailMC("signalFailMC", "", RooArgSet(m), dhfail);
    
    RooRealVar  zPeakP1    ("zPeakP1","",91.,85.,101);
    RooRealVar  zWidthP1   ("zWidthP1","",2.1,0.,3.4);
    RooRealVar  sigmaP1    ("sigmaP1","",2.,0.,10.);
    RooVoigtian zPeakPass  ("zPeakPass","",m,zPeakP1,zWidthP1,sigmaP1); // Z-peak as Voigtian in pass sample
    
    RooRealVar  meanP2     ("meanP2","",62.,60.,88.);
    RooRealVar  sigmaP2    ("sigmaP2","",2.,1.,10.);
    RooRealVar  alphaP2    ("alphaP2","",0.1,-10,10.);
    RooRealVar  nP2        ("nP2","",1.,0.,10.);
    RooCBShape  zTailPass  ("zTailPass","",m,meanP2,sigmaP2,alphaP2,nP2); //CB as resolution model in pass sample
    
    RooRealVar  fracPass   ("fracPass","",0.9,0.,1.);
    RooAddPdf   signalPassAna ("signalPass","",RooArgList(zPeakPass,zTailPass),fracPass); // we will trust MC signal for passing sample in the analytical fit
    
    // create a CB for the resolution
    RooRealVar zPeakF1     ("zPeakF1","",91.,85.,101);
    RooRealVar zWidthF1    ("zWidthF1","",2.1,0.,3.4);
    RooRealVar sigmaF1     ("sigmaF1","",2.,0.,10.);
    RooVoigtian zPeakFail  ("zPeakFail","",m,zPeakF1,zWidthF1,sigmaF1); // Z-peak as Voigtian in fail sample 
    
    // poly --> fail sample
    RooRealVar   a0        ("a0","a0",0.5,-10.,10.) ;
    RooRealVar   a1        ("a1","a1",0.5,-10.,10.) ;
    RooChebychev zPol      ("zRooChebychev","",m,RooArgList(a0,a1));
    
    // CB --> fail sample
    RooRealVar  meanF2     ("meanF2","",65.,60.,85.);
    RooRealVar  sigmaF2    ("sigmaF2","",2.,1.,20.);
    RooRealVar  alphaF2    ("alphaF2","",0.1,-10,10.);
    RooRealVar  nF2        ("nF2","",1.,-10.,10.);
    RooCBShape  zCB        ("zCrystalBall","",m,meanF2,sigmaF2,alphaF2,nF2);  
    
    // RooCMSShape --> fail sample
    RooRealVar alphaFail ("alphaFail","alphaFail",125,30,300);
    RooRealVar betaFail  ("betaFail","betaFail",0.01,-1,1);
    RooRealVar gammaFail ("gammaFail","gammaFail",0.05,-0.5,0.5);
    RooRealVar peakFail  ("peakFail","peakFail",91.2,60,120);
    RooCMSShape zRooCMSShape ("zRooCMSShape","",m,alphaFail,betaFail,gammaFail,peakFail);
    
    // RooErfExp --> fail sample
    RooRealVar   cExp ("cExp","",   -0.005,-1,0.);
    RooRealVar   offset ("offset","",90,60,120);
    RooRealVar   width ("width","",  10,0.,25);
    RooErfExpPdf erfExp ("zRooErfExp","",m,cExp,offset,width);
    
    RooRealVar  fracFail   ("fracFail","",0.9,0.,1.);
    RooAddPdf*  signalFailAna = NULL;  // failing signal will be convoluted with a gaussian.

    if(usePolynomial) // make fail sample as pol2 + Voig
      signalFailAna = new RooAddPdf("signalFailAna","",RooArgList(zPeakFail,zPol),fracFail);
    else if(useCBShape) // make fail sample as CB + Voig
      signalFailAna = new RooAddPdf("signalFailAna","",RooArgList(zPeakFail,zCB),fracFail);    
    else if(useRooCMSShape or (useMixedModel and (ptMin.at(size) < 20 ))) // make fail sample as RooCMS + Voig
      signalFailAna = new RooAddPdf("signalFailAna","",RooArgList(zPeakFail,zRooCMSShape),fracFail);
    else if(useErfExpPdf or (useMixedModel and  (ptMin.at(size) >= 20 )))  // make fail sample as RooExpPdf + Voig 
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
    
    // Store info in workspace
    RooWorkspace w("w", "");
    // RooDataHist used for efficiency
    w.import(dhpass);
    w.import(dhfail);
    w.import(signalPassMC);
    w.import(signalFailMC);
    // fix pdf parameters before storing the workspace --> set all of them to constant
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

    // plotting part
    TCanvas* c1 = new TCanvas("c1","",600,650);
    c1->cd();
    RooPlot* frame = m.frame();
    frame->SetTitle("");
    frame->GetXaxis()->SetTitle("mass (GeV)");
    frame->GetYaxis()->SetTitle("a.u.");
    frame->GetYaxis()->SetTitleOffset(1.1);
    
    if(plotRescaled){ // plot shape rescaled to 1
      dhpass.plotOn(frame,RooFit::LineColor(kRed),RooFit::MarkerColor(kRed),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::DataError(RooAbsData::SumW2),RooFit::DrawOption("EP"),RooFit::Rescale(1./hpass.at(size).Integral()));
      dhfail.plotOn(frame,RooFit::LineColor(kBlack),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::DataError(RooAbsData::SumW2),RooFit::DrawOption("EP"),RooFit::Rescale(1./hfail.at(size).Integral()));    
      hpass.at(size).Scale(1./hpass.at(size).Integral());
      hfail.at(size).Scale(1./hfail.at(size).Integral());
      frame->SetMaximum(max(hpass.at(size).GetMaximum(),hfail.at(size).GetMaximum())*1.2);
      pdfPass->plotOn(frame,RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::DrawOption("L"),RooFit::Normalization(1,RooAbsReal::NumEvent));
      pdfFail->plotOn(frame,RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::DrawOption("L"),RooFit::Normalization(1,RooAbsReal::NumEvent));
    }
    else{ // use proper MC normalization (scaled to lumi if xsec is properly set)
      dhpass.plotOn(frame,RooFit::LineColor(kRed),RooFit::MarkerColor(kRed),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::DataError(RooAbsData::SumW2),RooFit::DrawOption("EP"),RooFit::Binning(40,m.getMin(),m.getMax()));
      dhfail.plotOn(frame,RooFit::LineColor(kBlack),RooFit::MarkerColor(kBlack),RooFit::MarkerSize(1),RooFit::MarkerStyle(20),RooFit::DataError(RooAbsData::SumW2),RooFit::DrawOption("EP"),RooFit::Binning(40,m.getMin(),m.getMax()));
      pdfPass->plotOn(frame,RooFit::LineColor(kRed),RooFit::LineWidth(2),RooFit::DrawOption("L"));
      pdfFail->plotOn(frame,RooFit::LineColor(kBlack),RooFit::LineWidth(2),RooFit::DrawOption("L"));
    }

    frame->Draw();  
    TLegend* leg = new TLegend(0.2,0.75,0.5,0.92);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    hpass.at(size).SetLineColor(kRed);
    hfail.at(size).SetLineColor(kBlack);
    hfail.at(size).SetLineWidth(2);
    hpass.at(size).SetLineWidth(2);
    leg->AddEntry(&hpass.at(size),"Pass Probes","L");
    leg->AddEntry(&hfail.at(size),"Fail Probes","L");
    leg->Draw("same");
    CMS_lumi(c1,"",true);
    c1->SaveAs((outputDIR+"/"+string(fileName)+".png").c_str(),"png");
    c1->SaveAs((outputDIR+"/"+string(fileName)+".pdf").c_str(),"pdf");
    
    if(frame) delete frame;
    if(leg) delete leg;
    if(c1) delete c1;
  }
  if(pufile) pufile->Close();
}

/// main function
void makeTnPTemplates(string inputDIR, string leptonType, string outputDIR, bool isRecoEfficiency = false, bool isSummer16 = true) {

  // Basic style options
  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  TGaxis::SetMaxDigits(3);

  // Load custom pdf
  gSystem->Load("PDFs/RooErfExpPdf_cc.so");

  // Parse few options
  isRecoEff = isRecoEfficiency;
  isSUMMER16 = isSummer16;

  // Basic checks
  if(leptonType != "muon" and leptonType!= "electron" and leptonType!= "photon"){
    cerr<<"Problem with the lepton type : muon or electron or photons --> return "<<endl;
    return;
  }
    

  ////sequence for lepon id template maker
  if(not isRecoEff){ 
    if(leptonType == "muon"){
	    string idName = "looseid"; // loose id templates
	    maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinMuon,etaBinMuon,nvtxBinMuon,true);    
	    idName = "tightid"; // tight id templates
	    maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinMuon,etaBinMuon,nvtxBinMuon,true);    
    }
  
    else if(leptonType == "electron"){
      string idName = "vetoid"; // vetoid templates
      maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinElectron,etaBinElectron,nvtxBinElectron);    	
      idName = "tightid";// tight id templates
      maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinElectron,etaBinElectron,nvtxBinElectron);    	
    }
    else if (leptonType == "photon"){
      string idName = "mediumid"; // medium id templates
      maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinPhoton,etaBinPhoton,nvtxBinPhoton);
    }
  }
  // template for reco efficiency
  else{
    if (leptonType == "muon"){
      string idName = "trackerid"; // muon reco efficiency
      maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinMuonReco,etaBinMuonReco,nvtxBinMuonReco,true);    
    }
    else if(leptonType == "electron"){
      string idName = "recoelectronmatch";
      leptonType = "photon"; // electron reco efficiecny
      maketemplate(inputDIR,leptonType,outputDIR,idName,ptBinElectronReco,etaBinElectronReco,nvtxBinElectronReco);    
    }
  }
}
