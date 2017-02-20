#include "TnPBinning.h"
#include "PDFs/RooCMSShape.h"

class TagAndProbeBin {
  
public:
  TagAndProbeBin();
  TagAndProbeBin(double  ptmin, double  ptmax, double  etamin, double  etamax, double  nvtxmin, double  nvtxmax){
    ptMin = ptmin;
    ptMax = ptmax;
    etaMin = etamin;
    etaMax = etamax;
    nvtxMin = nvtxmin;
    nvtxMax = nvtxmax;
  }

  double ptMin;
  double ptMax;
  double etaMin;
  double etaMax;
  double nvtxMin;
  double nvtxMax;

};

// only valid for analytical fits
bool fixBackgroundParameters = true;
bool floatSignalShapeParameters = true;

// for the time being only binned fit can be done with the same binning adopted in the templateFile
void makeTagAndProbeFits(string inputDIR, // input directory where download files
			 string outputDIR, // output dir typically on afs
			 string templateFileName, // file where templates and pdfs fits on MC are places
			 string typeID, // id to be considered --> branch in the input tree
			 int leptonPID, // pdgId of the lepton			 
			 TagAndProbeBin selectionBin, // binning to be used for selection
			 string backgroundType, // two options: RooCMSShape or Exp
			 string signalType = "", //Nominal, Alternative template or analytical fit
			 bool isEOSDIR = true, // tell wheter the directory is on a local pc or EOS
			 bool isRecoEff = false, // tell if is a reco efficiency measirement or id
			 bool isAbsEta = true, // apply the selection on absolute eta
			 bool doAnalyticalFit = false // if analytical fit need to be perfomred
			 ){

  gSystem->Load("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeTagAndProbe/PDFs/RooCMSShape_cc.so");
  gSystem->Load("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/macros/makeTagAndProbe/PDFs/RooErfExpPdf_cc.so");
  gROOT->SetBatch(kTRUE);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);

  cout<<"############ Tag and probe analysis ###########"<<endl;
  cout<<"######### Template file "<<templateFileName<<endl;
  cout<<"######### Lepton id "<<typeID<<endl;
  cout<<"######### Lepton pdg id "<<leptonPID<<endl;
  cout<<"######### backgroundType "<<backgroundType<<endl;
  cout<<"######### signalType "<<signalType<<endl;
  cout<<"######### isRecoEff "<<isRecoEff<<endl;
  cout<<"######### [ptMin,ptMax] "<<selectionBin.ptMin<<" "<<selectionBin.ptMax<<endl;
  cout<<"######### [etaMin,etaMax] "<<selectionBin.etaMin<<" "<<selectionBin.etaMax<<endl;
  cout<<"######### [nvtxMin,nvtxMax] "<<selectionBin.nvtxMin<<" "<<selectionBin.nvtxMax<<endl;
  cout<<"###############################################"<<endl;

  ///////////////////
  string dirName;
  if (leptonPID == 13)
    dirName = "muontnptree";
  else if(leptonPID == 11 and not isRecoEff)
    dirName = "electrontnptree";
  else if(leptonPID == 11 and isRecoEff)
    dirName = "photontnptree";
  else if(leptonPID == 22)
    dirName = "photontnptree";
  else{    
    cerr<<"Lepton type not recognized --> return"<<endl;
    return;
  }

  /////////////
  TChain* tree = new TChain((dirName+"/fitter_tree").c_str(),(dirName+"/fitter_tree").c_str());  
  if(not isEOSDIR)
    tree->Add((inputDIR+"/*root").c_str());
  else{
    system(("/afs/cern.ch/project/eos/installation/cms/bin/eos.select ls "+inputDIR+" | grep root > file.temp").c_str());
    ifstream infile;
    infile.open("file.temp");
    if(infile.is_open()){
      string line;
      while(!infile.eof()){
	getline(infile,line);
	if(TString(line).Contains("root") and line != ""){	  
	  tree->Add(("root://eoscms.cern.ch//eos/cms/"+inputDIR+"/"+line).c_str());
	}
      }
    }
    infile.close();
    system("rm -r file.temp");
  }
  
  // checks
  if(leptonPID == 11){
    if(typeID != "vetoid" and typeID != "tightid" and typeID != "recoelectronmatch"){
      cerr<<"ID not recognized for electrons -> exit"<<endl;
      return;
    }
  }      
  
  else if(leptonPID == 13){
    if(typeID != "looseid" and typeID != "tightid" and typeID != "trackerid"){
      cerr<<"ID not recognized for muons -> exit"<<endl;
      return;
    }
  }      
  else if(leptonPID == 22){
    if(typeID != "mediumid"){
      cerr<<"ID not recognized for photons -> exit"<<endl;
      return;
    }
  }      

  // make output file name
  string postfix = "";
  if(leptonPID == 11 and not isRecoEff)
    postfix = "electron";
  else if(leptonPID == 11 and isRecoEff)
    postfix = "electronReco";
  else if (leptonPID == 13 and not isRecoEff)
    postfix = "muon";
  else if (leptonPID == 13 and isRecoEff)
    postfix = "muonReco";
  else if(leptonPID == 22)
    postfix = "photon";

  postfix += "_"+typeID+"_pt_"+Form("%.1f",selectionBin.ptMin)+"_"+Form("%.1f",selectionBin.ptMax)+"_eta_"+Form("%.2f",selectionBin.etaMin)+"_"+Form("%.2f",selectionBin.etaMax)+"_pu_"+Form("%.1f",selectionBin.nvtxMin)+"_"+Form("%.1f",selectionBin.nvtxMax);

  ///////////////////
  string outputFileName;
  outputFileName = "efficiency-data"+postfix;  
  if(backgroundType == "RooCMSShape")
    outputFileName +="_RooCMSShape";
  else if(backgroundType == "Exponential")
    outputFileName +="_Exponential";
  else if(signalType != "")
    outputFileName += "_Alternative";
  if(doAnalyticalFit)
    outputFileName += "_Analytical";

  /////////////////////////////////
  TTreeReader reader(tree);
  TTreeReaderValue<float> m    (reader, "mass");
  TTreeReaderValue<float> eta  (reader, "eta");
  TTreeReaderValue<float> phi  (reader, "phi");
  TTreeReaderValue<float> pt   (reader, "pt");
  TTreeReaderValue<float> nvtx (reader, "nvtx");
  TTreeReaderValue<int>   id   (reader, typeID.c_str());
 
  // open template file and load the binning
  system(("mkdir -p "+outputDIR).c_str());
  TFile* outputFile = new TFile((outputDIR+"/"+outputFileName+".root").c_str(),"RECREATE");
  outputFile->cd();

  TFile* templateFile = TFile::Open(templateFileName.c_str());
  RooWorkspace* workspace = (RooWorkspace*) templateFile->Get("w");
  RooRealVar* mass = (RooRealVar*) workspace->obj("mass");
  TH1F* passingEvents = new TH1F("passingEvents","passingEvents", mass->getBinning().numBins(), mass->getBinning().array());
  TH1F* failingEvents = new TH1F("failingEvents","failingEvents", mass->getBinning().numBins(), mass->getBinning().array());
  
  cout<<"######### Loop on data "<<endl;
  long int nEvents = 0;
  long int nPart   = 100000;
  long int nTotal  = tree->GetEntries();
  while(reader.Next()){    
    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;
    if(*pt <= selectionBin.ptMin) continue;
    if(*pt > selectionBin.ptMax) continue;
    if(isAbsEta and fabs(*eta) <= selectionBin.etaMin) continue;
    if(isAbsEta and fabs(*eta) > selectionBin.etaMax) continue;
    if(not isAbsEta and *eta <= selectionBin.etaMin) continue;
    if(not isAbsEta and *eta > selectionBin.etaMax) continue;
    if(*nvtx <= selectionBin.nvtxMin) continue;
    if(*nvtx > selectionBin.nvtxMax) continue;
    if(*id)
      passingEvents->Fill(*m);
    else
      failingEvents->Fill(*m);    
  }
  cout<<endl;

  ///////////////////////////
  cout<<"####### Passing events "<<passingEvents->Integral()<<endl;
  cout<<"####### Failing events "<<failingEvents->Integral()<<endl;

  RooDataHist* passDataHist = new RooDataHist("passDataHist","",RooArgList(*mass), RooFit::Import(*passingEvents), 0);
  RooDataHist* failDataHist = new RooDataHist("failDataHist","",RooArgList(*mass), RooFit::Import(*failingEvents), 0);

  cout<<"######### Building Sim pdf"<<endl;

  RooAbsPdf* signalPass = NULL;
  RooAbsPdf* signalFail = NULL;

  RooHistPdf* signalPassMC = (RooHistPdf*) workspace->pdf("signalPassMC");
  RooHistPdf* signalFailMC = (RooHistPdf*) workspace->pdf("signalFailMC");
  RooRealVar* mP = new RooRealVar("mP","",0.,-2.,2.);
  RooRealVar* mF = new RooRealVar("mF","",0.,-3.,3.);
  RooRealVar* sP = new RooRealVar("sP","",0.5,0.,2.);
  RooRealVar* sF = new RooRealVar("sF","",0.5,0.,2.);
  
  RooGaussian* signalPassSmear = new RooGaussian("signalPassSmear","",*mass,*mP,*sP);
  RooGaussian* signalFailSmear = new RooGaussian("signalFailSmear","",*mass,*mF,*sF);
  if(not doAnalyticalFit){
    signalPass = new RooFFTConvPdf("signalPass","",*mass,*signalPassMC,*signalPassSmear);
    signalFail = new RooFFTConvPdf("signalFail","",*mass,*signalFailMC,*signalFailSmear);
  }
  else{
    signalPass = new RooFFTConvPdf("signalPassPreFit","",*mass,*signalPassMC,*signalPassSmear);
    signalFail = new RooFFTConvPdf("signalFailPreFit","",*mass,*signalFailMC,*signalFailSmear);
  }
  signalPass->Print();
  signalFail->Print();

  /////////////////////////////////
  RooAbsPdf* backgroundPass = NULL;
  RooAbsPdf* backgroundFail = NULL;
  if(backgroundType == "RooCMSShape"){
    RooRealVar* alphaBkgPass = new RooRealVar("alphaBkgPass","",105.,80.,150.);
    RooRealVar* alphaBkgFail = new RooRealVar("alphaBkgFail","",125.,80.,150.);
    RooRealVar* betaBkgPass = new RooRealVar("betaBkgPass","",0.01,-1.,1.);
    RooRealVar* betaBkgFail = new RooRealVar("betaBkgFail","",0.01,-1.,1.);
    RooRealVar* gammaBkgPass = new RooRealVar("gammaBkgPass","",0.05,-0.1,0.5);
    RooRealVar* gammaBkgFail = new RooRealVar("gammaBkgFail","",0.05,-0.1,0.5);
    RooRealVar* peakBkgPass = new RooRealVar("peakBkgPass","",91.2,90,92.5);
    RooRealVar* peakBkgFail = new RooRealVar("peakBkgFail","",91.2,90,92.5);
    backgroundPass = new RooCMSShape("backgroundPass","",*mass,*alphaBkgPass,*betaBkgPass,*gammaBkgPass,*peakBkgPass);
    backgroundFail = new RooCMSShape("backgroundFail","",*mass,*alphaBkgFail,*betaBkgFail,*gammaBkgFail,*peakBkgFail);    
  }
  else if(backgroundType == "Exponential"){
    RooRealVar* cBkgP = new RooRealVar("cBkgP","",-0.05,-2,2);
    RooRealVar* cBkgF = new RooRealVar("cBkgF","",-0.05,-2,2);
    backgroundPass = new RooExponential("backgroundPass","",*mass,*cBkgP);
    backgroundFail = new RooExponential("backgroundFail","",*mass,*cBkgF);
  }
  else{
    cerr<<"######### Not recognized background type --> exit "<<endl;
    return;
  }

  backgroundPass->Print();
  backgroundFail->Print();

  ///////////////////////
  RooRealVar* efficiency = new RooRealVar("efficiency","",0.8,0,1);
  RooRealVar* effBkg     = new RooRealVar("effBkg","",0.8,0,1);
  RooRealVar* numTot     = new RooRealVar("numTot","",passingEvents->Integral()+failingEvents->Integral(),0.,(passingEvents->Integral()+failingEvents->Integral())*10.);
  RooRealVar* fSigAll    = new RooRealVar("fSigAll","",0.9,0.,1.);
  RooFormulaVar* nSigPass = new RooFormulaVar("nSigPass","@0*@1*@2",RooArgList(*numTot,*fSigAll,*efficiency));
  RooFormulaVar* nSigFail = new RooFormulaVar("nSigFail","@0*@1*(1-@2)",RooArgList(*numTot,*fSigAll,*efficiency));
  RooFormulaVar* nBkgPass = new RooFormulaVar("nBkgPass","@0*(1-@1)*@2",RooArgList(*numTot,*fSigAll,*effBkg));
  RooFormulaVar* nBkgFail = new RooFormulaVar("nBkgFail","@0*(1-@1)*(1-@2)",RooArgList(*numTot,*fSigAll,*effBkg));

  efficiency->Print();
  effBkg->Print();
  numTot->Print();
  fSigAll->Print();
  nSigPass->Print();
  nSigFail->Print();
  nBkgPass->Print();
  nBkgFail->Print();

  RooExtendPdf* signalPass_extend = new RooExtendPdf("signalPass_extend","",*signalPass,*nSigPass); 
  RooExtendPdf* signalFail_extend = new RooExtendPdf("signalFail_extend","",*signalFail,*nSigFail); 
  RooExtendPdf* backgroundPass_extend = new RooExtendPdf("backgroundPass_extend","",*backgroundPass,*nBkgPass); 
  RooExtendPdf* backgroundFail_extend = new RooExtendPdf("backgroundFail_extend","",*backgroundFail,*nBkgFail); 

  signalPass_extend->Print();
  signalFail_extend->Print();
  backgroundPass_extend->Print();
  backgroundFail_extend->Print();

  RooAddPdf* pdfPass = new RooAddPdf("pdfPass","",RooArgList(*signalPass_extend,*backgroundPass_extend));
  RooAddPdf* pdfFail = new RooAddPdf("pdfFail","",RooArgList(*signalFail_extend,*backgroundFail_extend));
  pdfPass->Print();
  pdfFail->Print();

  ////////////////////////
  RooCategory category("category","category") ;
  category.defineType("Pass");
  category.defineType("Fail");

  RooSimultaneous* simPdf =  new RooSimultaneous("simPdfAna","simultaneous pdf",category) ;
  simPdf->addPdf(*pdfPass,"Pass");
  simPdf->addPdf(*pdfFail,"Fail");
  simPdf->Print();

  // combined dataset                                                                                                                                                                                 
  std::map <string,RooDataHist*> histoMap;
  histoMap ["Pass"] = passDataHist;
  histoMap ["Fail"] = failDataHist;
  RooDataHist* combData = new RooDataHist( "combData", "combined data", RooArgList(*mass),category,histoMap);  
  combData->Print();

  // make the fit
  cout<<"######### Start the fit"<<endl;
  RooAbsReal* simNLL = simPdf->createNLL(*combData,RooFit::Extended(true));
  simNLL->Print();
  RooMinimizer minimizer(*simNLL);   
  minimizer.Print();
  cout<<"######### Start a scan fit "<<endl;
  minimizer.minimize("Minuit2","Scan");  
  cout<<"######### Start a migrad fit "<<endl;
  minimizer.minimize("Minuit2","migrad");
  cout<<"######### Start a migrad improved fit "<<endl;
  minimizer.minimize("Minuit2","hesse");
  std::cout<<"######### Minimizer Minuit2 --> Extended + Minos on efficiency  ###########"<<std::endl;
  RooFitResult* result = simPdf->fitTo(*combData, RooFit::Save(true), RooFit::Extended(true), RooFit::Strategy(2),RooFit::Minos(*efficiency),RooFit::Strategy(2));

  // re-do the fit fixing background yields as well to the fit with templates
  if(doAnalyticalFit){
    signalPass = (RooAbsPdf*) workspace->pdf("signalPass");
    if(floatSignalShapeParameters){
      RooArgList paramPdfPassSig = RooArgList(*signalPass->getParameters(combData));
      for(int iparam = 0; iparam < paramPdfPassSig.getSize(); iparam++){
	((RooRealVar*) paramPdfPassSig.at(iparam))->setConstant(kFALSE);
      }
    }
    RooAbsPdf*   signalFailAna = (RooAbsPdf*) workspace->pdf("signalFailAna");
    RooGaussian* signalFailSmear_ana = new RooGaussian("signalFailSmear_ana","",*mass,*mF,*sF);
    signalFail = new RooFFTConvPdf("signalFail","",*mass,*signalFailAna,*signalFailSmear_ana);   
    if(floatSignalShapeParameters){
      RooArgList paramPdfPassFail = RooArgList(*signalFailAna->getParameters(combData));
      for(int iparam = 0; iparam < paramPdfPassFail.getSize(); iparam++){
	((RooRealVar*) paramPdfPassFail.at(iparam))->setConstant(kFALSE);
      }
    }
    RooExtendPdf* signalPass_extend_ana = new RooExtendPdf("signalPass_extend_ana","",*signalPass,*nSigPass); 
    RooExtendPdf* signalFail_extend_ana = new RooExtendPdf("signalFail_extend_ana","",*signalFail,*nSigFail); 
    RooExtendPdf* backgroundPass_extend_ana = new RooExtendPdf("backgroundPass_extend_ana","",*backgroundPass,*nBkgPass); 
    RooExtendPdf* backgroundFail_extend_ana = new RooExtendPdf("backgroundFail_extend_ana","",*backgroundFail,*nBkgFail);
    // fix background from initial fit
    if(fixBackgroundParameters){
      fSigAll->setConstant(kTRUE);
      effBkg->setConstant(kTRUE);
      RooArgList paramPdfPass = RooArgList(*backgroundPass->getParameters(combData));
      for(int iparam = 0; iparam < paramPdfPass.getSize(); iparam++){
	((RooRealVar*) paramPdfPass.at(iparam))->setConstant(kTRUE);
      }
      RooArgList paramPdfFail = RooArgList(*backgroundFail->getParameters(combData));
      for(int iparam = 0; iparam < paramPdfFail.getSize(); iparam++)
	((RooRealVar*) paramPdfFail.at(iparam))->setConstant(kTRUE);
    }

    signalPass_extend_ana->Print();
    signalFail_extend_ana->Print();
    backgroundPass_extend_ana->Print();
    backgroundFail_extend_ana->Print();
    ///
    pdfPass = new RooAddPdf("pdfPass","",RooArgList(*signalPass_extend_ana,*backgroundPass_extend_ana));
    pdfFail = new RooAddPdf("pdfFail","",RooArgList(*signalFail_extend_ana,*backgroundFail_extend_ana));
    pdfPass->Print();
    pdfFail->Print();

    simPdf =  new RooSimultaneous("simPdf","simultaneous pdf",category) ;
    simPdf->addPdf(*pdfPass,"Pass");
    simPdf->addPdf(*pdfFail,"Fail");
    simPdf->Print();

    // make the fit                                                                                                                                                                                
    cout<<"######### Start the fit"<<endl;
    simNLL = simPdf->createNLL(*combData,RooFit::Extended(true));
    simNLL->Print();
    RooMinimizer minimizer(*simNLL);
    minimizer.Print();
    cout<<"######### Start a scan fit "<<endl;
    minimizer.minimize("Minuit2","Scan");
    cout<<"######### Start a migrad fit "<<endl;
    minimizer.minimize("Minuit2","migrad");
    cout<<"######### Start a migrad improved fit "<<endl;
    minimizer.minimize("Minuit2","hesse");
    std::cout<<"######### Minimizer Minuit2 --> Extended + Minos on efficiency  ###########"<<std::endl;
    result = simPdf->fitTo(*combData, RooFit::Save(true), RooFit::Extended(true), RooFit::Strategy(2),RooFit::Minos(*efficiency),RooFit::Strategy(2));
  }

  cout<<"Store in output file"<<endl;
  outputFile->cd();
  RooWorkspace w("w", "");
  w.import(*passDataHist);
  w.import(*failDataHist);
  w.import(*simPdf);
  w.import(*result,"fitresults");
  w.saveSnapshot("finalState",w.components());
  RooRealVar* e = (RooRealVar*) result->floatParsFinal().find("efficiency");
  efficiency->setVal(e->getVal());
  Double_t errLo = e->getErrorLo(), errHi = e->getErrorHi();
  if (errLo == 0 && e->getVal() < 0.5) errLo = e->getMin()-e->getVal();
  if (errHi == 0 && e->getVal() > 0.5) errHi = e->getMax()-e->getVal();
  efficiency->setAsymError(errLo, errHi);
  
  if (passingEvents->Integral() * failingEvents->Integral() == 0) {
    RooRealVar* nTot = (RooRealVar*) result->floatParsFinal().find("numTot");
    RooRealVar* fSig = (RooRealVar*) result->floatParsFinal().find("fSigAll");
    double cerr = ROOT::Math::beta_quantile( 1-(1.0-.68540158589942957)/2, 1, nTot->getVal() * fSig->getVal() ); 
    if (passingEvents == 0) {
      efficiency->setVal(0);
      efficiency->setAsymError(0,cerr);
    } else {
      efficiency->setVal(1);
      efficiency->setAsymError(-cerr,0);
    }
  }
  w.Write();
  outputFile->Close();
  cout<<"Close application "<<endl;
}
