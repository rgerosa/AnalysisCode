#include "makePhotonPurityUtils.h"
#include <fstream>

// to decide which data to run
vector<string> RunEra = {"Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G","Run2016H"};
// photon pt bins
static vector<float> ptBins = {175,200,225,250,280,320,360,400,500,650,1000};
static int nBinPhotonIso  = 35;
static float photonIsoMax = 12;
static float photonIsoMin = 0;
static bool debug = true;

void makePhotonPurityFit(string inputDirectory, // directory with dataFiles
			 float  lumi, // luminosity
			 string outputDIR,
			 bool   addSystematics = false,
			 string inputDirectorySignalMC = "",
			 string inputDirectoryBackgroundMC = ""
			 ){

  system(("mkdir -p "+outputDIR).c_str());

  //from twiki https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Selection_implementation_details
  photonID mediumID (0.0396,0.01022,0.441,2.725,0.0148,0.000017,2.571,0.0047); // set wp for medium id

  // bins for purity and histograms
  vector<fitPurity> dataHisto;
  vector<fitPurity> signalTemplateRND04;
  vector<fitPurity> signalTemplateRND08;
  vector<fitPurity> backgroundTemplateData;

  for(size_t ibin = 0; ibin < ptBins.size()-1; ibin++){

    dataHisto.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
				  new TH1F(Form("dataHisto_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso,photonIsoMin,photonIsoMax)));  
    signalTemplateRND04.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
					    new TH1F(Form("signalTemplateRND04_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso,photonIsoMin,photonIsoMax)));  
    signalTemplateRND08.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
					    new TH1F(Form("signalTemplateRND08_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso,photonIsoMin,photonIsoMax)));  
    backgroundTemplateData.push_back(fitPurity(ptBins.at(ibin),ptBins.at(ibin+1),
					       new TH1F(Form("backgroundTemplateData_pt_%d_%d",int(ptBins.at(ibin)),int(ptBins.at(ibin+1))),"",nBinPhotonIso,photonIsoMin,photonIsoMax)));  
    
    dataHisto.back().phHisto->Sumw2();
    signalTemplateRND04.back().phHisto->Sumw2();
    signalTemplateRND08.back().phHisto->Sumw2();
    backgroundTemplateData.back().phHisto->Sumw2();
  }
  
  // style
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  // add specific files to the chain --> data
  TChain* chain_data = new TChain("tree/tree");
  system(("ls "+inputDirectory+" | grep root > file.list").c_str());
  ifstream infile("file.list");
  string line;
  if(infile.is_open()){
    while(getline(infile,line)){
      if(line == "") continue;
      if(not TString(line).Contains("root")) continue;
      bool found = false;
      for(auto era : RunEra){
	if(TString(line).Contains(era.c_str()))
	  found = true;
      }
      if(found)
	chain_data->Add((inputDirectory+"/"+line).c_str());
    }
  }
  
  // set all branches
  TTreeReader reader (chain_data);
  TTreeReaderValue<float> phpt  (reader,"pPuritypt");
  TTreeReaderValue<float> pheta (reader,"pPurityeta");
  TTreeReaderValue<float> phphi (reader,"pPurityphi");
  TTreeReaderValue<float> phElVeto (reader,"phPurityElectronVeto");
  TTreeReaderValue<float> phPHIso (reader,"phPurityPHiso");
  TTreeReaderValue<float> phCHIso (reader,"phPurityCHiso");
  TTreeReaderValue<float> phNHIso (reader,"phPurityNHiso");
  TTreeReaderValue<float> phHoE   (reader,"phPurityhoe");
  TTreeReaderValue<float> phSieie (reader,"phPuritysieie");
  TTreeReaderValue<float> phPHIsoRND04 (reader,"phPurityRND04PHiso");
  TTreeReaderValue<float> phPHIsoRND08 (reader,"phPurityRND08PHiso");
  TTreeReaderValue<float> phEAEgamma (reader,"phPurityEAEGamma");
  TTreeReaderValue<float> rho (reader,"rho");
  TTreeReaderValue<float> xsec (reader,"xsec");
  TTreeReaderValue<float> wgt (reader,"wgt");
  TTreeReaderValue<double> wgtsum (reader,"wgtsum");
  TTreeReaderValue<unsigned int> nphotons (reader,"nphotons");
  TTreeReaderValue<unsigned int> nphotonsPurity (reader,"nphotonsPurity");
  TTreeReaderValue<unsigned int> nelectrons (reader,"nelectrons");
  TTreeReaderValue<unsigned int> nmuons (reader,"nmuons");
  TTreeReaderValue<unsigned int> ntausraw (reader,"ntausrawold");
  TTreeReaderValue<unsigned int> nbjets (reader,"nbjetslowpt");
  TTreeReaderValue<unsigned int> njets (reader,"njets");
  TTreeReaderValue<UChar_t> hltp165     (reader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175     (reader,"hltphoton175");
  TTreeReaderValue<unsigned int> run         (reader,"run");
  TTreeReaderValue<unsigned int> lumisection (reader,"lumi");
  TTreeReaderValue<unsigned int> event       (reader,"event");
  TTreeReaderValue<unsigned int> nvtx        (reader,"nvtx");
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
  TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
  TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
  TTreeReaderValue<vector<float> > chfrac  (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac  (reader,"combinejetNHfrac");
  TTreeReaderValue<float> pmet        (reader,"t1phmet");
  TTreeReaderValue<float> pmetphi     (reader,"t1phmetphi");
  TTreeReaderValue<float> metpf       (reader,"pfmet");
  TTreeReaderValue<float> metcalo     (reader,"calomet");
  TTreeReaderValue<float> jpmdphi (reader,"incjetphmetdphimin4");
  TTreeReaderValue<float> wzpt   (reader,"wzpt");
  TTreeReaderValue<float> wzeta  (reader,"wzeta");
  TTreeReaderValue<float> wzphi  (reader,"wzphi");
  TTreeReaderValue<int>   wzid  (reader,"wzid");
  TTreeReaderValue<int>   ismatch  (reader,"ismatch");

  // output file
  TFile* outputFile = new TFile((outputDIR+"/PhotonPurityFitResult.root").c_str(),"RECREATE");
  outputFile->cd();

  
  // loop on data events
  cout<<"Number of evedns in data "<<chain_data->GetEntries()<<endl;
  long int nTotal = chain_data->GetEntries();
  long int nEvents = 0;
  while(reader.Next()){
    cout.flush();
    if(nEvents % 100000 == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;
    
    if(not (*hltp165 or *hltp175)) continue;
    if(*fhbhe == 0 or *fhbiso == 0 or *fcsct == 0 or *feeb == 0 or *fetp == 0 or *fvtx == 0 or *fbadmu == 0 or *fbadch == 0 or *fcsc == 0) continue;
    if(fabs(*pheta) > 1.4442) continue;
    if(*nmuons != 0) continue;
    if(*nelectrons != 0) continue;
    if(*nbjets != 0) continue;
    if(*ntausraw != 0) continue;
    if(*njets < 1) continue;
    if(jetpt->at(0) < 100) continue;
    if(fabs(jeteta->at(0)) > 2.5) continue;
    if(chfrac->at(0) < 0.1) continue;
    if(nhfrac->at(0) > 0.8) continue;
    if(*jpmdphi < 0.5) continue;

    // apply photon id
    if(*nphotonsPurity < 1) continue;
    if(*phpt < 175) continue;
    if(*phElVeto == 0) continue;
    if(*phCHIso > mediumID.chadiso) continue; // already corrected for effective area
    if(*phNHIso > mediumID.nhadiso0+mediumID.nhadiso1*(*phpt) + mediumID.nhadiso2*(*phpt)*(*phpt)) continue; // already corrected for effective area
    if(*phHoE > mediumID.HoE) continue;

    // select the histogram given the pt of the photon
    unsigned int bin = 0;
    for( ; bin < dataHisto.size(); bin++){
      if(*phpt >= dataHisto.at(bin).ptMin and *phpt < dataHisto.at(bin).ptMax) break;	
    }

    // distributions in data without effective area correction
    if(*phSieie < mediumID.sigmaieie){
      dataHisto.at(bin).phHisto->Fill(max(0.,double(*phPHIso-*rho*(*phEAEgamma))));
      signalTemplateRND04.at(bin).phHisto->Fill(max(0.,double(*phPHIsoRND04-*rho*(*phEAEgamma))));
      signalTemplateRND08.at(bin).phHisto->Fill(max(0.,double(*phPHIsoRND08-*rho*(*phEAEgamma))));
      dataHisto.at(bin).ptMean += *phpt;
      signalTemplateRND04.at(bin).ptMean += *phpt;
      signalTemplateRND08.at(bin).ptMean += *phpt;      
    }
    else{
      backgroundTemplateData.at(bin).phHisto->Fill(max(0.,double(*phPHIso-*rho*(*phEAEgamma))));
      backgroundTemplateData.at(bin).ptMean += *phpt;
    }					           
  }

  for(auto bin: dataHisto)
    bin.ptMean = bin.ptMean/bin.phHisto->Integral();
  for(auto bin: signalTemplateRND04)
    bin.ptMean = bin.ptMean/bin.phHisto->Integral();
  for(auto bin: signalTemplateRND08)
    bin.ptMean = bin.ptMean/bin.phHisto->Integral();
  for(auto bin: backgroundTemplateData)
    bin.ptMean = bin.ptMean/bin.phHisto->Integral();

  // Build Model for fit
  vector<RooRealVar*>   observable; // one observable for each pt-bin
  vector<RooWorkspace*> worksapceRND04;
  vector<RooWorkspace*> worksapceRND08;
  vector<RooDataHist*>  dataHisto_roo; // conversion from TH1F* to RooDataHist*
  vector<RooDataHist*>  signalTemplateRND04_roo;
  vector<RooDataHist*>  signalTemplateRND08_roo;
  vector<RooDataHist*>  backgroundTemplateData_roo;
  vector<RooHistPdf*>   signalTemplateRND04Pdf; //Conversion from histogram to binned Pdfs
  vector<RooHistPdf*>   signalTemplateRND08Pdf;
  vector<RooHistPdf*>   backgroundTemplateDataRND04Pdf;
  vector<RooHistPdf*>   backgroundTemplateDataRND08Pdf;
  vector<RooAbsPdf*>    totalPdfRND04;  // total Pdf used in the fit
  vector<RooAbsPdf*>    totalPdfRND08;
  vector<RooAbsReal*>   nllRND04;
  vector<RooAbsReal*>   nllRND08;
  vector<RooFitResult*> fitResultRND04; // fit result to extract parameters
  vector<RooFitResult*> fitResultRND08; 

  TGraphAsymmErrors* photonPurityRND04 = new TGraphAsymmErrors();
  TGraphAsymmErrors* photonPurityRND08 = new TGraphAsymmErrors();

  TCanvas* canvas = new TCanvas("canvas","canvas",600,700);
  canvas->cd();

  for(size_t isize = 0; isize < dataHisto.size(); isize++){

    int ptMin = int(dataHisto.at(isize).ptMin);
    int ptMax = int(dataHisto.at(isize).ptMax);
    
    // create observable
    observable.push_back(new RooRealVar(Form("phiso_pt_%d_%d",ptMin,ptMax),"",0,dataHisto.at(isize).phHisto->GetXaxis()->GetBinLowEdge(1),
					dataHisto.at(isize).phHisto->GetXaxis()->GetBinLowEdge(dataHisto.at(isize).phHisto->GetNbinsX())));
    //create workspace
    worksapceRND04.push_back(new RooWorkspace(Form("wsRND04_pt_%d_%d",ptMin,ptMax),Form("wsRND04_pt_%d_%d",ptMin,ptMax)));
    worksapceRND08.push_back(new RooWorkspace(Form("wsRND08_pt_%d_%d",ptMin,ptMax),Form("wsRND08_pt_%d_%d",ptMin,ptMax)));

    // data-histogram
    RooArgList observable_list (*observable.back());
    dataHisto_roo.push_back(new RooDataHist(Form("RooDataHisto_pt_%d_%d",ptMin,ptMax),"",observable_list, dataHisto.at(isize).phHisto));
    signalTemplateRND04_roo.push_back(new RooDataHist(Form("RooSignalTemplateDataRND04_pt_%d_%d",ptMin,ptMax),"",observable_list,signalTemplateRND04.at(isize).phHisto));
    signalTemplateRND08_roo.push_back(new RooDataHist(Form("RooSignalTemplateDataRND08_pt_%d_%d",ptMin,ptMax),"",observable_list,signalTemplateRND08.at(isize).phHisto));
    backgroundTemplateData_roo.push_back(new RooDataHist(Form("RooBackgroundTemplateData_pt_%d_%d",ptMin,ptMax),"",observable_list, backgroundTemplateData.at(isize).phHisto));

    // integral in data wrt neutral iso selection
    double dataIntegralErr = 0;
    double dataIntegral =  dataHisto.at(isize).phHisto->IntegralAndError(dataHisto.at(isize).phHisto->FindBin(dataHisto.at(isize).phHisto->GetBinLowEdge(0)),
									dataHisto.at(isize).phHisto->FindBin(mediumID.phiso0+mediumID.phiso1*dataHisto.at(isize).ptMean),dataIntegralErr);
    
    // Pdfs
    signalTemplateRND04Pdf.push_back(new RooHistPdf(Form("signalTemplateRND04Pdf_pt_%d_%d",ptMin,ptMax),"",observable_list,*signalTemplateRND04_roo.back()));
    signalTemplateRND08Pdf.push_back(new RooHistPdf(Form("signalTemplateRND08Pdf_pt_%d_%d",ptMin,ptMax),"",observable_list,*signalTemplateRND08_roo.back()));
    backgroundTemplateDataRND04Pdf.push_back(new RooHistPdf(Form("backgroundTemplateDataRND04Pdf_pt_%d_%d",ptMin,ptMax),"",observable_list,*backgroundTemplateData_roo.back()));
    backgroundTemplateDataRND08Pdf.push_back(new RooHistPdf(Form("backgroundTemplateDataRND08Pdf_pt_%d_%d",ptMin,ptMax),"",observable_list,*backgroundTemplateData_roo.back()));

    // make total Pdf
    RooRealVar sigNormRND04 (Form("sigNormRND04_pt_%d_%d",ptMin,ptMax),"",signalTemplateRND04.at(isize).phHisto->Integral(),0,signalTemplateRND04.at(isize).phHisto->Integral()*10);   
    RooRealVar sigNormRND08 (Form("sigNormRND08_pt_%d_%d",ptMin,ptMax),"",signalTemplateRND08.at(isize).phHisto->Integral(),0,signalTemplateRND08.at(isize).phHisto->Integral()*10);   
    RooRealVar bkgNormDataRND04  (Form("bkgNormDataRND04_pt_%d_%d",ptMin,ptMax),"",backgroundTemplateData.at(isize).phHisto->Integral(),0,backgroundTemplateData.at(isize).phHisto->Integral()*10);
    RooRealVar bkgNormDataRND08  (Form("bkgNormDataRND08_pt_%d_%d",ptMin,ptMax),"",backgroundTemplateData.at(isize).phHisto->Integral(),0,backgroundTemplateData.at(isize).phHisto->Integral()*10);

    RooExtendPdf signalExtendRND04Pdf (Form("signalExtendRND04Pdf_pt_%d_%d",ptMin,ptMax),"",*signalTemplateRND04Pdf.back(),sigNormRND04);    
    RooExtendPdf signalExtendRND08Pdf (Form("signalExtendRND08Pdf_pt_%d_%d",ptMin,ptMax),"",*signalTemplateRND08Pdf.back(),sigNormRND08);
    RooExtendPdf bkgExtendDataRND04Pdf (Form("bkgExdendDataRND04Pdf_pt_%d_%d",ptMin,ptMax),"",*backgroundTemplateDataRND04Pdf.back(),bkgNormDataRND04);
    RooExtendPdf bkgExtendDataRND08Pdf (Form("bkgExdendDataRND04Pdf_pt_%d_%d",ptMin,ptMax),"",*backgroundTemplateDataRND08Pdf.back(),bkgNormDataRND08);

    totalPdfRND04.push_back(new RooAddPdf(Form("totalPdfRND04_%d_%d",ptMin,ptMax),"",signalExtendRND04Pdf,bkgExtendDataRND04Pdf));
    totalPdfRND08.push_back(new RooAddPdf(Form("totalPdfRND08_%d_%d",ptMin,ptMax),"",signalExtendRND08Pdf,bkgExtendDataRND08Pdf));
    
    // create NLL
    if(not debug){
      nllRND04.push_back(totalPdfRND04.back()->createNLL(*dataHisto_roo.back(),RooFit::Extended(),RooFit::Verbose(-1),RooFit::SumW2Error(kTRUE)));
      nllRND08.push_back(totalPdfRND08.back()->createNLL(*dataHisto_roo.back(),RooFit::Extended(),RooFit::Verbose(-1),RooFit::SumW2Error(kTRUE)));
    }
    else{
      nllRND04.push_back(totalPdfRND04.back()->createNLL(*dataHisto_roo.back(),RooFit::Extended(),RooFit::SumW2Error(kTRUE)));
      nllRND08.push_back(totalPdfRND08.back()->createNLL(*dataHisto_roo.back(),RooFit::Extended(),RooFit::SumW2Error(kTRUE)));
    }

    //make fits    
    RooMinimizer mfitRND04(*nllRND04.back());
    RooMinimizer mfitRND08(*nllRND08.back());
    if(not debug){
      mfitRND04.setVerbose(kFALSE);
      mfitRND04.setPrintLevel(-1);
      mfitRND08.setVerbose(kFALSE);
      mfitRND08.setPrintLevel(-1);
    }

    cout<<"Minimize random cone R = 0.4 for ptMin "<<ptMin<<" ptMax "<<ptMax<<endl;
    mfitRND04.minimize("Minuit2","minimize");
    cout<<"Minimize random cone R = 0.4 hesse"<<endl;
    mfitRND04.minimize("Minuit2","hesse");
    cout<<"Estimate minos errors for all parameters"<<endl;
    mfitRND04.minos();
    fitResultRND04.push_back(mfitRND04.save(Form("fitResult_pt_%d_%d",ptMin,ptMax)));
    fitResultRND04.back()->Print();

    cout<<"Minimize random cone R = 0.8 for ptMin "<<ptMin<<" ptMax "<<ptMax<<endl;
    mfitRND08.minimize("Minuit2","minimize");
    cout<<"Minimize random cone R = 0.8 hesse"<<endl;
    mfitRND08.minimize("Minuit2","hesse");
    cout<<"Estimate minos errors for all parameters"<<endl;
    mfitRND08.minos();
    fitResultRND08.push_back(mfitRND08.save(Form("fitResult_pt_%d_%d",ptMin,ptMax)));
    fitResultRND08.back()->Print();
    
    // import data in the workspace
    worksapceRND04.back()->import(*dataHisto_roo.back());
    worksapceRND04.back()->import(*totalPdfRND04.back());
    worksapceRND04.back()->import(*fitResultRND04.back());

    // import data in the workspace
    worksapceRND08.back()->import(*dataHisto_roo.back());
    worksapceRND08.back()->import(*totalPdfRND08.back());
    worksapceRND08.back()->import(*fitResultRND08.back());

    // calculate the photon purity
    cout<<"####### Purity RND 04 "<<endl;
    cout<<"Photon pt bin "<<ptMin<<" "<<ptMax<<" observed event "<<dataHisto_roo.back()->sumEntries()<<" signal events RND04: "<<sigNormRND04.getVal()<<" pm "<<sigNormRND04.getError()
	<<" background events "<<bkgNormDataRND04.getVal()<<" pm "<<bkgNormDataRND04.getError()<<endl;
    
    // integrate pdfs
    observable.back()->setRange("isolated",dataHisto.at(isize).phHisto->GetBinLowEdge(0),mediumID.phiso0+mediumID.phiso1*dataHisto.at(isize).ptMean);
    RooRealVar* int_sigRND04 = (RooRealVar*) signalExtendRND04Pdf.createIntegral(*observable.back(),RooFit::NormSet(*observable.back()),RooFit::Range("isolated"));
    RooRealVar* int_bkgRND04 = (RooRealVar*) bkgExtendDataRND04Pdf.createIntegral(*observable.back(),RooFit::NormSet(*observable.back()),RooFit::Range("isolated"));

    cout<<"Integral for isolation < 2.571+0.0047*pt: observed rate "<<dataIntegral<<" pm "<<dataIntegralErr<<" signal events RND04: "<<int_sigRND04->getVal()<<" pm "<<int_sigRND04->getError()<<" background events "<<int_bkgRND04->getVal()<<" pm "<<int_bkgRND04->getError()<<endl;

    // fill purity information
    double purityRND04 = int_sigRND04->getVal()/(int_sigRND04->getVal()+int_bkgRND04->getVal());
    double purityErrRND04_lo =  sqrt(int_sigRND04->getErrorLo()*int_sigRND04->getErrorLo()*(int_bkgRND04->getVal()*int_bkgRND04->getVal()/pow((int_sigRND04->getVal()+int_bkgRND04->getVal()),4))+
				    int_bkgRND04->getErrorLo()*int_bkgRND04->getErrorLo()*(int_sigRND04->getVal()*int_sigRND04->getVal()/pow((int_sigRND04->getVal()+int_bkgRND04->getVal()),4)));
    double purityErrRND04_hi =  sqrt(int_sigRND04->getErrorHi()*int_sigRND04->getErrorHi()*(int_bkgRND04->getVal()*int_bkgRND04->getVal()/pow((int_sigRND04->getVal()+int_bkgRND04->getVal()),4))+
				    int_bkgRND04->getErrorHi()*int_bkgRND04->getErrorHi()*(int_sigRND04->getVal()*int_sigRND04->getVal()/pow((int_sigRND04->getVal()+int_bkgRND04->getVal()),4)));

    cout<<"Purity value "<<purityRND04<<" - "<<purityErrRND04_lo<<" + "<<purityErrRND04_hi<<endl;

    photonPurityRND04->SetPoint(isize,(ptMax+ptMin)/2,purityRND04);
    photonPurityRND04->SetPointError(isize,(ptMax-ptMin)/2,(ptMax-ptMin)/2,purityErrRND04_lo,purityErrRND04_hi);

    RooRealVar phPurityRND04 (Form("photonPurityRND04_pt_%d_%d",ptMin,ptMax),"",purityRND04,0,1);
    phPurityRND04.setAsymError(purityErrRND04_lo,purityErrRND04_hi);
    worksapceRND04.back()->import(phPurityRND04);

    /////////////
    cout<<"####### Purity RND 08 "<<endl;
    cout<<"Photon pt bin "<<ptMin<<" "<<ptMax<<" observed event "<<dataHisto_roo.back()->sumEntries()<<" signal events RND08: "<<sigNormRND08.getVal()<<" pm "<<sigNormRND08.getError()
	<<" background events "<<bkgNormDataRND08.getVal()<<" pm "<<bkgNormDataRND08.getError()<<endl;
    
    // integrate pdfs
    observable.back()->setRange("isolated",dataHisto.at(isize).phHisto->GetBinLowEdge(0),mediumID.phiso0+mediumID.phiso1*dataHisto.at(isize).ptMean);
    RooRealVar* int_sigRND08 = (RooRealVar*) signalExtendRND08Pdf.createIntegral(*observable.back(),RooFit::NormSet(*observable.back()),RooFit::Range("isolated"));
    RooRealVar* int_bkgRND08 = (RooRealVar*) bkgExtendDataRND08Pdf.createIntegral(*observable.back(),RooFit::NormSet(*observable.back()),RooFit::Range("isolated"));

    cout<<"Integral for isolation < 2.571+0.0087*pt: observed rate "<<dataIntegral<<" pm "<<dataIntegralErr<<" signal events RND08: "<<int_sigRND08->getVal()<<" pm "<<int_sigRND08->getError()<<" background events "<<int_bkgRND08->getVal()<<" pm "<<int_bkgRND08->getError()<<endl;

    // fill purity information
    double purityRND08 = int_sigRND08->getVal()/(int_sigRND08->getVal()+int_bkgRND08->getVal());
    double purityErrRND08_lo =  sqrt(int_sigRND08->getErrorLo()*int_sigRND08->getErrorLo()*(int_bkgRND08->getVal()*int_bkgRND08->getVal()/pow((int_sigRND08->getVal()+int_bkgRND08->getVal()),4))+
				    int_bkgRND08->getErrorLo()*int_bkgRND08->getErrorLo()*(int_sigRND08->getVal()*int_sigRND08->getVal()/pow((int_sigRND08->getVal()+int_bkgRND08->getVal()),4)));
    double purityErrRND08_hi =  sqrt(int_sigRND08->getErrorHi()*int_sigRND08->getErrorHi()*(int_bkgRND08->getVal()*int_bkgRND08->getVal()/pow((int_sigRND08->getVal()+int_bkgRND08->getVal()),4))+
				    int_bkgRND08->getErrorHi()*int_bkgRND08->getErrorHi()*(int_sigRND08->getVal()*int_sigRND08->getVal()/pow((int_sigRND08->getVal()+int_bkgRND08->getVal()),4)));

    cout<<"Purity value "<<purityRND08<<" - "<<purityErrRND08_lo<<" + "<<purityErrRND08_hi<<endl;

    photonPurityRND08->SetPoint(isize,(ptMax+ptMin)/2,purityRND08);
    photonPurityRND08->SetPointError(isize,(ptMax-ptMin)/2,(ptMax-ptMin)/2,purityErrRND08_lo,purityErrRND08_hi);

    RooRealVar phPurityRND08 (Form("photonPurityRND08_pt_%d_%d",ptMin,ptMax),"",purityRND08,0,1);
    phPurityRND08.setAsymError(purityErrRND08_lo,purityErrRND08_hi);
    worksapceRND08.back()->import(phPurityRND08);

    // save in the output file
    worksapceRND04.back()->writeToFile(outputFile->GetName());
    // save in the output file
    worksapceRND08.back()->writeToFile(outputFile->GetName());

    // plot fit result:
    plotFitResult(canvas,dataHisto.at(isize).phHisto,signalExtendRND04Pdf,bkgExtendDataRND04Pdf,*observable.back(),outputDIR,ptMin,ptMax,"RND04");
    plotFitResult(canvas,dataHisto.at(isize).phHisto,signalExtendRND08Pdf,bkgExtendDataRND08Pdf,*observable.back(),outputDIR,ptMin,ptMax,"RND08");
  }
  
   
  // plot purity result (nominal results) 
  TCanvas* canvas2 = new TCanvas("canvas2","",600,650);
  canvas2->cd();
  canvas2->SetTickx(1);
  canvas2->SetTicky(1);
  canvas2->cd();
  canvas2->SetRightMargin(0.06);

  photonPurityRND04->SetLineColor(kBlack);
  photonPurityRND04->SetMarkerColor(kBlack);
  photonPurityRND04->SetMarkerStyle(20);
  photonPurityRND04->SetMarkerSize(1);
  photonPurityRND04->GetXaxis()->SetTitle("photon p_{T} [GeV]");
  photonPurityRND04->GetYaxis()->SetTitle("purity");
  photonPurityRND04->GetYaxis()->SetRangeUser(0.5,1.2);
  photonPurityRND04->Draw("AP");

  photonPurityRND08->SetLineColor(kRed);
  photonPurityRND08->SetMarkerColor(kRed);
  photonPurityRND08->SetMarkerStyle(20);
  photonPurityRND08->SetMarkerSize(1);
  photonPurityRND08->Draw("Psame");

  CMS_lumi(canvas2,"36.2");
  
  TLegend leg (0.6,0.3,0.9,0.6);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(photonPurityRND04,"Purity with #DeltaR = 0.4","PE");
  leg.AddEntry(photonPurityRND08,"Purity with #DeltaR = 0.8","PE");
  leg.Draw("same");

  canvas2->SaveAs((outputDIR+"/photonPurityResultData.png").c_str(),"png");
  canvas2->SaveAs((outputDIR+"/photonPurityResultData.pdf").c_str(),"pdf");

  if(addSystematics){ // implement sys uncertainties using MC

  }

  outputFile->Close();
}
