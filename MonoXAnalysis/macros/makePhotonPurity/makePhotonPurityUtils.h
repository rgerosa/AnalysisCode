#include "../CMS_lumi.h"

static float deltaRMatching = 0.3;

// to define samples
enum class Sample {data, gjets, qcd};

// photon id values --> class to handle it                                                                                                                                                           
class photonID {

 public:
 photonID(const float & HoE, const float & sigmaieie, const float&  sigmaieie_sideband, const float & chadiso, const float & nhadiso0, const float & nhadiso1, const float & nhadiso2, const float & phiso0, const float & phiso1):
  HoE(HoE),
    sigmaieie(sigmaieie),
    sigmaieie_sideband(sigmaieie_sideband),
    chadiso(chadiso),
    nhadiso0(nhadiso0),
    nhadiso1(nhadiso1),
    nhadiso2(nhadiso2),
    phiso0(phiso0),
    phiso1(phiso1){
    };
  ~photonID(){};

  float HoE;
  float sigmaieie;
  float sigmaieie_sideband;
  float chadiso;
  float nhadiso0;
  float nhadiso1;
  float nhadiso2;
  float phiso0;
  float phiso1;
};

// class for bin fir                                                                                                                                                                                 
class fitPurity {
 public:
 fitPurity(const float & ptMin, const float & ptMax, TH1F* histo):
  ptMin(ptMin),
    ptMax(ptMax),
    phHisto(histo){
      ptMean = 0;
  }
    ~fitPurity(){};

  float ptMin;
  float ptMax;
  float ptMean;
  TH1F* phHisto;
};

/////////////
void fillHistograms(TTree* chain,const Sample & sample, 
		    vector<fitPurity>  & dataHisto, 
		    vector<fitPurity>  & signalTemplateRND04, 
		    vector<fitPurity>  & signalTemplateRND08,
		    vector<fitPurity>  & backgroundTemplate,
		    const photonID & mediumID,
		    const vector<TH1*> & khists,
		    const float & lumi = 36.2){


  /// apply pileup and trigger turn on
  TFile* pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_35p9fb.root");
  TH1*   puhist = (TH1*) pufile->Get("puhist");
  ///
  TFile* triggerfile_SinglePhoton  = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/Monojet/photonTriggerEfficiency_photon.root");
  TEfficiency* triggerphoton       = (TEfficiency*) triggerfile_SinglePhoton->Get("eff_recoil");
  TGraphAsymmErrors* triggerphoton_graph = triggerphoton->CreateGraph();

  /// set all branches                                                                                                                                                                                
  TTreeReader reader (chain);
  TTreeReaderValue<float> phpt  (reader,"phPuritypt");
  TTreeReaderArray<float> pheta (reader,"pPurityheta.phPurityeta");
  TTreeReaderValue<float> phphi (reader,"phPurityphi");
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
  TTreeReaderValue<unsigned int> nphotons (reader,"nphotons");
  TTreeReaderValue<unsigned int> nphotonsPurity (reader,"nphotonsPurity");
  TTreeReaderValue<unsigned int> nelectrons (reader,"nelectrons");
  TTreeReaderValue<unsigned int> nmuons (reader,"nmuons");
  TTreeReaderValue<unsigned int> ntausraw (reader,"ntaus");
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

  // loop on data events                                                                                                                                                                        
  if(sample == Sample::data) 
    cout<<"Number of evedns in data "<<chain->GetEntries()<<endl;  
  else if(sample == Sample::qcd)
    cout<<"Number of evedns in QCD MC "<<chain->GetEntries()<<endl;  
  else if(sample == Sample::gjets)
    cout<<"Number of evedns in Gjets MC "<<chain->GetEntries()<<endl;  


  // calculate sum of weights in case of MC sample
  vector<double> wgtsum;
  string currentFile = "";
  if(sample == Sample::qcd or sample == Sample::gjets){
    while(reader.Next()){
      if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName()  != currentFile){
	currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() ;
	wgtsum.push_back(0);
	wgtsum.back() += *wgt;
      }
      else
	wgtsum.back() += *wgt;
    }
    reader.SetEntry(0);
  }
  
  long int nTotal = chain->GetEntries();
  long int nEvents = 0;

  currentFile = "";
  int ifile = 0;
  while(reader.Next()){
    
    if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      ifile ++;
    }
    else if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile == ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
    }
      
    cout.flush();
    if(nEvents % 100000 == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    if(not (*hltp165 or *hltp175)) continue;
    if(*fhbhe == 0 or *fhbiso == 0 or *fcsct == 0 or *feeb == 0 or *fetp == 0 or *fvtx == 0 or *fbadmu == 0 or *fbadch == 0) continue;
    if(fabs(pheta[0]) > 1.4442) continue;
    if(*nmuons     != 0) continue;
    if(*nelectrons != 0) continue;
    if(*nbjets     != 0) continue;
    if(*ntausraw   != 0) continue;
    if(*njets < 1) continue;
    if(jetpt->at(0) < 100) continue;
    if(fabs(jeteta->at(0)) > 2.5) continue;
    if(chfrac->at(0) < 0.1) continue;
    if(nhfrac->at(0) > 0.8) continue;
    
    // apply photon id: except sigma-ieta-ieta (later) and electron veto as well as photon isolation                                                                                      
    if(*nphotonsPurity < 1) continue;
    if(*phpt < 175) continue;
    if(*phCHIso > mediumID.chadiso) continue; // already corrected for effective area                                                                                                          
    if(*phNHIso > mediumID.nhadiso0+mediumID.nhadiso1*(*phpt) + mediumID.nhadiso2*(*phpt)*(*phpt)) continue; // already corrected for effective area                                           
    if(*phHoE   > mediumID.HoE)  continue;    

    // select the histogram given the pt of the photon                                                                                                                                           
    unsigned int bin = 0;
    if(signalTemplateRND04.size() != 0){
      for( ; bin < signalTemplateRND04.size()-1; bin++){
	if(*phpt >= signalTemplateRND04.at(bin).ptMin and *phpt < signalTemplateRND04.at(bin).ptMax) break;
      }
      if(*phpt >  signalTemplateRND04.back().ptMax)
	bin = signalTemplateRND04.size()-1;
    }
    else if(backgroundTemplate.size() != 0){
      for( ; bin < backgroundTemplate.size()-1; bin++){
	if(*phpt >= backgroundTemplate.at(bin).ptMin and *phpt < backgroundTemplate.at(bin).ptMax) break;
      }
      if(*phpt >  backgroundTemplate.back().ptMax)
	bin = backgroundTemplate.size()-1;
    }
    
    double evtwgt = 1.;
    if(sample == Sample::gjets or sample == Sample::qcd){
      evtwgt *= lumi*(*xsec)*(*wgt);
      if (*nvtx <= 60)
	evtwgt *= puhist->GetBinContent(puhist->FindBin(*nvtx));      
      evtwgt *= triggerphoton_graph->Eval(min(double(*pmet),triggerphoton_graph->GetXaxis()->GetXmax()));
    }

    // apply k-factor
    double kwgt = 1.0;      

    // distributions in data without effective area correction                                                                                                                                      
    if(*phSieie < mediumID.sigmaieie){
      // adding track veto --> in order to add more stat in the background template
      if(*phElVeto == 0) continue;
      if(sample == Sample::data){ // fill also data histogram
	dataHisto.at(bin).phHisto->Fill(max(0.,double(*phPHIso-*rho*(*phEAEgamma))));
	signalTemplateRND04.at(bin).phHisto->Fill(max(0.,double(*phPHIsoRND04-*rho*(*phEAEgamma))));
	signalTemplateRND08.at(bin).phHisto->Fill(max(0.,double(*phPHIsoRND08-*rho*(*phEAEgamma))));
	dataHisto.at(bin).ptMean += *phpt;
	signalTemplateRND04.at(bin).ptMean += *phpt;
	signalTemplateRND08.at(bin).ptMean += *phpt;
      }
      else if(sample == Sample::gjets){// don't fill data histogram

	double genpt = *wzpt;
	for (size_t i = 0; i < khists.size(); i++) {
	  if (khists[i]) {
	    if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) 
	      genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	    if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) 
	      genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
	    kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
	  }
	}

	// matching condition
	if(*wzpt <= 0) continue;
	float deltaEta_gen = fabs(pheta[0]-*wzeta);
	float deltaPhi_gen = fabs(*phphi-*wzphi);
	if(deltaPhi_gen > TMath::Pi()) deltaPhi_gen = 2*TMath::Pi() - deltaPhi_gen;
	if(sqrt(deltaEta_gen*deltaEta_gen+deltaPhi_gen*deltaPhi_gen) > deltaRMatching) continue;
	if(not *ismatch) continue;
	// fill histograms
	signalTemplateRND04.at(bin).phHisto->Fill(max(0.,double(*phPHIso-*rho*(*phEAEgamma))),evtwgt*kwgt/(wgtsum.at(ifile)));
	signalTemplateRND04.at(bin).ptMean += *phpt;
      }
    }
    else if(*phSieie > mediumID.sigmaieie and *phSieie < mediumID.sigmaieie_sideband){
      if(sample == Sample::data){
	backgroundTemplate.at(bin).phHisto->Fill(max(0.,double(*phPHIso-*rho*(*phEAEgamma))));
	backgroundTemplate.at(bin).ptMean += *phpt;
      }
      else if(sample == Sample::qcd){
	if(*wzpt <= 0) continue;
	// matching condition                                                                                                                                                                         
        float deltaEta_gen = fabs(pheta[0]-*wzeta);
        float deltaPhi_gen = fabs(*phphi-*wzphi);
        if(deltaPhi_gen > TMath::Pi()) deltaPhi_gen = 2*TMath::Pi() - deltaPhi_gen;
        if(sqrt(deltaEta_gen*deltaEta_gen+deltaPhi_gen*deltaPhi_gen) < deltaRMatching) continue;
	////
	backgroundTemplate.at(bin).phHisto->Fill(max(0.,double(*phPHIso-*rho*(*phEAEgamma))),evtwgt*kwgt/(wgtsum.at(ifile)));
	backgroundTemplate.at(bin).ptMean += *phpt;
      }
    }
  }

  cout<<endl;
  if(pufile) pufile->Close();
  if(triggerfile_SinglePhoton) triggerfile_SinglePhoton->Close();
}


/////////////
void makePurityFit(RooWorkspace* ws, 
		   const fitPurity & dataTemplate, const fitPurity & signalTemplate, const fitPurity & backgroundTemplate, 
		   const photonID & mediumID,
		   const bool & debug){

  // create observable                                                                                                                                                                               
  RooRealVar observable ("photoniso","",0,dataTemplate.phHisto->GetXaxis()->GetBinLowEdge(1),dataTemplate.phHisto->GetXaxis()->GetBinLowEdge(dataTemplate.phHisto->GetNbinsX()+1));
  if(debug)
    observable.Print();

  // data-histogram                                                                                                                                                                                  
  RooArgList observable_list (observable);
  RooDataHist RooDataHisto      ("data","",observable_list, dataTemplate.phHisto);
  RooDataHist RooSignalTemplate ("signalTemplateData","",observable_list,signalTemplate.phHisto);
  RooDataHist RooBackgroundTemplate ("backgroundTemplateData","",observable_list, backgroundTemplate.phHisto);

  if(debug){
    RooDataHisto.Print();
    RooSignalTemplate.Print();
    RooBackgroundTemplate.Print();
  }

  // integral in data wrt neutral iso selection                                                                                                                                                      
  double dataIntegralErr = 0;
  double dataIntegral =  dataTemplate.phHisto->IntegralAndError(dataTemplate.phHisto->FindBin(dataTemplate.phHisto->GetBinLowEdge(0)),
								       dataTemplate.phHisto->FindBin(mediumID.phiso0+mediumID.phiso1*dataTemplate.ptMean),dataIntegralErr);

  if(debug)
    cout<<"Data integral : "<<dataIntegral<<" error "<<dataIntegralErr<<endl;


  // Pdfs                                                                                                                                                                                            
  RooHistPdf signalTemplatePdf ("signalTemplatePdf","",observable_list,RooSignalTemplate);
  RooHistPdf backgroundTemplatePdf ("backgroundTemplatePdf","",observable_list,RooBackgroundTemplate);
  if(debug){
    signalTemplatePdf.Print();
    backgroundTemplatePdf.Print();
  }

  // make total Pdf                                                                                                                                                                                  
  RooRealVar signalNorm     ("signalNorm","",signalTemplate.phHisto->Integral()*0.9,0,signalTemplate.phHisto->Integral()*10);
  RooRealVar backgroundNorm ("backgroundNorm","",signalTemplate.phHisto->Integral()*0.1,0,signalTemplate.phHisto->Integral()*10);
  if(debug){
    signalNorm.Print();
    backgroundNorm.Print();
  }
    
  RooExtendPdf signalExtendPdf ("signalExtendPdf","",signalTemplatePdf,signalNorm);
  RooExtendPdf backgroundExtendPdf ("backgroundExdendPdf","",backgroundTemplatePdf,backgroundNorm);

  if(debug){
    signalExtendPdf.Print();
    backgroundExtendPdf.Print();
  }

  // total pdf
  RooAddPdf totalPdf ("totalPdf","",signalExtendPdf,backgroundExtendPdf);
  if(debug)
    totalPdf.Print();

  RooAbsReal* nll = NULL;
  if(not debug)
    nll = totalPdf.createNLL(RooDataHisto,RooFit::Extended(kTRUE),RooFit::Verbose(-1));
  else
    nll = totalPdf.createNLL(RooDataHisto,RooFit::Extended(kTRUE));

  if(debug)
    nll->Print();
  
  //make fits                                                                                                                                                                                        
  RooMinimizer mfit(*nll);
  if(not debug){
    mfit.setVerbose(kFALSE);
    mfit.setPrintLevel(-1);
    mfit.setVerbose(kFALSE);
    mfit.setPrintLevel(-1);
  }

  cout<<"Minimize for ptMin "<<dataTemplate.ptMin<<" ptMax "<<dataTemplate.ptMax<<endl;
  mfit.minimize("Minuit2","minimize");
  cout<<"Minimize hesse"<<endl;
  mfit.minimize("Minuit2","hesse");
  cout<<"Estimate minos errors for all parameters"<<endl;
  mfit.minos();
  RooFitResult* fitResult = mfit.save("fitResult");
  if(debug)
    fitResult->Print();

  // import data in the workspace                                                                                                                                                                    
  ws->import(RooDataHisto);
  ws->import(RooSignalTemplate);
  ws->import(RooBackgroundTemplate);
  ws->import(totalPdf);
  ws->import(*nll);
  ws->import(*fitResult);

  // calculate the photon purity                                                                                                                                                                     
  if(debug)
    cout<<"####### Purity measurement for photon pt bin "<<dataTemplate.ptMin<<" "<<dataTemplate.ptMax<<" observed event "<<RooDataHisto.sumEntries()<<" signal events: "<<signalNorm.getVal()<<" pm "<<signalNorm.getError()<<" background events "<<backgroundNorm.getVal()<<" pm "<<backgroundNorm.getError()<<endl;

  // integrate pdfs                                                                                                                                                                                  
  observable.setRange("isolated",dataTemplate.phHisto->GetBinLowEdge(0),mediumID.phiso0+mediumID.phiso1*dataTemplate.ptMean);
  RooRealVar* int_sig = (RooRealVar*) signalExtendPdf.createIntegral(observable,RooFit::NormSet(observable),RooFit::Range("isolated"));
  RooRealVar* int_bkg = (RooRealVar*) backgroundExtendPdf.createIntegral(observable,RooFit::NormSet(observable),RooFit::Range("isolated"));
  
  if(debug)
    cout<<"Integral for isolation < 2.571+0.0047*pt: observed rate "<<dataIntegral<<" pm "<<dataIntegralErr<<" signal events : "<<int_sig->getVal()<<" pm "<<int_sig->getError()<<" background events "<<int_bkg->getVal()<<" pm "<<int_bkg->getError()<<endl;

  // fill purity information                                                                                                                                                                         
  double purity       = int_sig->getVal()/(int_sig->getVal()+int_bkg->getVal());
  double purityErr_lo =  sqrt(int_sig->getErrorLo()*int_sig->getErrorLo()*(int_bkg->getVal()*int_bkg->getVal()/pow((int_sig->getVal()+int_bkg->getVal()),4))+
				   int_bkg->getErrorLo()*int_bkg->getErrorLo()*(int_sig->getVal()*int_sig->getVal()/pow((int_sig->getVal()+int_bkg->getVal()),4)));
  double purityErr_hi =  sqrt(int_sig->getErrorHi()*int_sig->getErrorHi()*(int_bkg->getVal()*int_bkg->getVal()/pow((int_sig->getVal()+int_bkg->getVal()),4))+
				   int_bkg->getErrorHi()*int_bkg->getErrorHi()*(int_sig->getVal()*int_sig->getVal()/pow((int_sig->getVal()+int_bkg->getVal()),4)));

  if(debug)
    cout<<"Purity value: "<<purity<<" - "<<purityErr_lo<<" + "<<purityErr_hi<<endl;

  RooRealVar phPurity ("photonPurity","",purity,0,1);
  phPurity.setAsymError(purityErr_lo,purityErr_hi);
  ws->import(phPurity);
  
}



////////////
void plotFitResult(TCanvas* canvas, TH1F* data, const RooAbsPdf & signal, const RooAbsPdf & background, const RooRealVar & x, 
		   const string & outputDIR, const int & ptMin, const int & ptMax, const string & postfix){

  //create histograms from pdfs
  TH1F* signal_hist = (TH1F*) signal.createHistogram(Form("signal%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax),x,RooFit::Binning(data->GetNbinsX(),data->GetXaxis()->GetBinLowEdge(1),data->GetXaxis()->GetBinLowEdge(data->GetNbinsX()+1)));
  signal_hist->Sumw2();
  TH1F* background_hist = (TH1F*) background.createHistogram(Form("background%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax),x,RooFit::Binning(data->GetNbinsX(),data->GetXaxis()->GetBinLowEdge(1),data->GetXaxis()->GetBinLowEdge(data->GetNbinsX()+1)));
  background_hist->Sumw2();

  TH1F* totalHist = (TH1F*) signal_hist->Clone(Form("totalHist%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax));
  totalHist->Add(background_hist);
  
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);
  
  data->SetLineColor(kBlack);
  data->SetMarkerColor(kBlack);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(1);

  data->GetXaxis()->SetLabelSize(0);
  data->GetYaxis()->SetTitle("Events / GeV");
  data->GetYaxis()->SetRangeUser(data->GetMinimum()*0.1,data->GetMaximum()*100);
  data->Draw("EP");

  background_hist->SetLineColor(kGreen+1);
  background_hist->SetLineWidth(2);
  background_hist->Draw("hist same");

  signal_hist->SetLineColor(kBlue);
  signal_hist->SetLineWidth(2);
  signal_hist->Draw("hist same");

  totalHist->SetLineColor(kRed);
  totalHist->SetLineWidth(2);
  totalHist->Draw("hist same");

  TLegend leg (0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(data,"Data","PE");
  leg.AddEntry(totalHist,"Total Fit","L");
  leg.AddEntry(signal_hist,"Signal Component","L");
  leg.AddEntry(background_hist,"Background Component","L");
  leg.Draw("same");
  canvas->RedrawAxis("sameaxis");

  pad2->Draw();
  pad2->cd();

  TH1* frame =  (TH1*) data->Clone(Form("frame_pt_%d_%d",ptMin,ptMax));
  frame->Reset("ICES");  
  frame->GetYaxis()->SetRangeUser(0.4,1.6);
  frame->GetYaxis()->SetNdivisions(5);
  frame->GetYaxis()->SetTitle("Data/Exp.");
  frame->GetXaxis()->SetTitle("Photon Isolation [GeV]");

  TH1F* ratio = (TH1F*) data->Clone(Form("ratio%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax));
  TH1F* denFit_noerr = (TH1F*) totalHist->Clone(Form("denFit_noerr%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax));
  TH1F* denFit = (TH1F*) totalHist->Clone(Form("denFit%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax));
  for(int iBin = 1; iBin < denFit->GetNbinsX()+1; iBin++)
    denFit_noerr->SetBinError(iBin,0.);

  ratio->Divide(denFit_noerr);
  denFit->Divide(denFit_noerr);
  denFit->SetFillColor(kGray);
  frame->Draw();
  ratio->Draw("EPsame");
  denFit->Draw("E2same");

  TF1* line = new TF1(Form("line%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax),"1",ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1));
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  ratio->Draw("EPsame");
  pad2->RedrawAxis("sameaxis");
    

  CMS_lumi(canvas,"36.2");

  canvas->SetLogy();

  canvas->SaveAs((outputDIR+"/photonPurity"+postfix+"_pt_"+to_string(ptMin)+"_"+to_string(ptMax)+".png").c_str());
  canvas->SaveAs((outputDIR+"/photonPurity"+postfix+"_pt_"+to_string(ptMin)+"_"+to_string(ptMax)+".pdf").c_str());

}
