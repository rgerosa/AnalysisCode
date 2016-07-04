// pt binning to study ID observables
#include "../CMS_lumi.h"

static vector<float> ptBin = {175.,1000};
const float luminosity = 2.61;
const int nvtxMin = 0;
const int nvtxMax = 50;

vector<TH1F*> photonPt_data_HoE;
vector<TH1F*> photonPt_data_SIetaIeta;
vector<TH1F*> photonPt_data_Eveto;
vector<TH1F*> photonPt_data_CHIso;
vector<TH1F*> photonPt_data_NHIso;
vector<TH1F*> photonPt_data_EMIso;

vector<TH1F*> photonPt_mc_HoE;
vector<TH1F*> photonPt_mc_SIetaIeta;
vector<TH1F*> photonPt_mc_Eveto;
vector<TH1F*> photonPt_mc_CHIso;
vector<TH1F*> photonPt_mc_NHIso;
vector<TH1F*> photonPt_mc_EMIso;

vector<TH1F*> photonEta_data;
vector<TH1F*> photonSCEta_data;
vector<TH1F*> photonSCEnergy_data;
vector<TH1F*> photonSCRawEnergy_data;
vector<TH1F*> photonHoverE_data;
vector<TH1F*> photonSigmaIetaIeta_data;
vector<TH1F*> photonElectronVeto_data;
vector<TH1F*> photonCHIso_data;
vector<TH1F*> photonNHIso_data;
vector<TH1F*> photonEMIso_data;

vector<TH1F*> photonEta_mc;
vector<TH1F*> photonSCEta_mc;
vector<TH1F*> photonSCEnergy_mc;
vector<TH1F*> photonSCRawEnergy_mc;
vector<TH1F*> photonHoverE_mc;
vector<TH1F*> photonSigmaIetaIeta_mc;
vector<TH1F*> photonElectronVeto_mc;
vector<TH1F*> photonCHIso_mc;
vector<TH1F*> photonNHIso_mc;
vector<TH1F*> photonEMIso_mc;

float getEffectiveArea(string isotype, float abseta, bool isData){

  if(isotype == "CH" and isData){
    return 0.;
  }
  else if(isotype == "NH" and isData){
    if(abseta <= 1) return 0.0599;
    else if(abseta > 1 and abseta <= 1.479) return 0.0819;
    else return 0.;
  }
  else if(isotype == "EM" and isData){
    if(abseta <=1 ) return 0.1271;
    else if(abseta > 1 and abseta <= 1.479) return 0.1101;
    else return 0.;
  }
  if(isotype == "CH" and not isData){
    return 0.;
  }
  else if(isotype == "NH" and not isData){
    if(abseta <= 1) return 0.0599;
    else if(abseta > 1 and abseta <= 1.479) return 0.0819;
    else return 0.;
  }
  else if(isotype == "EM" and not isData){
    if(abseta <=1 ) return 0.1271;
    else if(abseta > 1 and abseta <= 1.479) return 0.1101;
    else return 0.;
  }
  return 0.;
}

void plotHistogram(TH1* histo_data, TH1* histo_mc, TCanvas* canvas, const string & plotDIR, const string & observable){
  
  canvas->cd();
  canvas->Clear();
  TPad* pad1 = new TPad("pad1","pad1",0.,0.2,1.,1.);
  pad1->Draw();
  canvas->cd();
  TPad* pad2 = new TPad("pad2","pad2",0.,0.,1.,0.2);
  pad2->Draw();
  canvas->cd();
  pad1->cd();
  histo_data->GetXaxis()->SetTitle(observable.c_str());
  histo_data->GetYaxis()->SetTitle("Events");
  histo_data->SetLineColor(kBlack);
  histo_data->SetMarkerColor(kBlack);
  histo_data->SetMarkerStyle(20);
  histo_data->SetMarkerSize(1);
  histo_data->GetYaxis()->SetRangeUser(1,max(histo_data->GetMaximum(),histo_mc->GetMaximum())*1000);
  histo_mc->SetFillColor(kAzure-3);
  histo_mc->SetLineColor(kRed);
  // calculate impurity by hand assuming to be ~ 7%
  TH1F* histo_qcd = (TH1F*) histo_data->Clone("histo_qcd");
  histo_qcd->Scale(0.08);
  histo_qcd->SetFillColor(kGray-2);
  histo_qcd->SetLineColor(kBlack);
  histo_mc->Add(histo_qcd);

  histo_mc->SetLineWidth(2);
  histo_data->Draw("EP");
  histo_mc->Draw("hist same");
  histo_qcd->Draw("hist same");
  histo_data->Draw("EPsame");
  CMS_lumi(pad1,Form("%.2f",luminosity),true);
  TLegend* leg = new TLegend(0.75,0.75,0.9,0.9);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histo_data,"Data","PE");
  leg->AddEntry(histo_mc,"#gamma+jets","F");
  leg->AddEntry(histo_qcd,"QCD","F");
  leg->Draw("same");
  pad1->SetLogy();
  pad1->RedrawAxis("sameaxis");

  canvas->cd();
  pad2->cd();
  TH1* ratio = (TH1*) histo_data->Clone("ratio");
  ratio->GetXaxis()->SetLabelSize(0);
  ratio->Divide(histo_mc);
  ratio->GetXaxis()->SetTitle("");
  ratio->GetYaxis()->SetTitle("Data/MC");
  ratio->GetYaxis()->SetRangeUser(0.,2.);
  ratio->Draw("EP");
  TF1* line = new TF1("line","1",ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1));
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  canvas->cd();		      
  canvas->SaveAs((plotDIR+"/"+string(histo_data->GetName())+".png").c_str(),"png");
  canvas->SaveAs((plotDIR+"/"+string(histo_data->GetName())+".pdf").c_str(),"pdf");
}

void fillHistograms(TTree* tree, const float & ptMin, const float & ptMax, const size_t & ipt, const bool & isData, const double & wgtsum = 1, const bool & addRecoil = false){

  TTreeReader reader(tree);
  TTreeReaderValue<vector<double> > phoPt (reader,"photonPt");
  TTreeReaderValue<vector<double> > phoEta(reader,"photonEta");
  TTreeReaderValue<vector<double> > phoPhi(reader,"photonPhi");
  TTreeReaderValue<vector<double> > phoSCEta       (reader,"photonSCEta");
  TTreeReaderValue<vector<double> > phoSCEnergy    (reader,"photonSCEnergy");
  TTreeReaderValue<vector<double> > phoSCRawEnergy (reader,"photonSCRawEnergy");
  TTreeReaderValue<vector<double> > phoHoverE      (reader,"photonHOverE");
  TTreeReaderValue<vector<double> > phoSigmaIetaIeta (reader,"photonSigmaIetaIeta");
  TTreeReaderValue<vector<double> > phoElectronVeto  (reader,"photonElectronVeto");
  TTreeReaderValue<vector<double> > phoCHIso  (reader,"photonChargedIso");
  TTreeReaderValue<vector<double> > phoNHIso  (reader,"photonNeutralIso");
  TTreeReaderValue<vector<double> > phoEMIso  (reader,"photonEMIso");
  TTreeReaderValue<UChar_t> hltp165    (reader,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175    (reader,"hltphoton175");
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<unsigned int> nbjets    (reader,"nbjetslowpt");
  TTreeReaderValue<unsigned int> nvtx      (reader,"nvtx");
  TTreeReaderValue<unsigned int> ntaus     (reader,"ntausraw");
  TTreeReaderValue<unsigned int> nmuon     (reader,"nmuons");
  TTreeReaderValue<unsigned int> nelectron (reader,"nelectrons");
  TTreeReaderValue<unsigned int> nphotons (reader,"nphotons");
  TTreeReaderValue<double> xsec         (reader,"xsec");
  TTreeReaderValue<double> wgt          (reader,"wgt");
  TTreeReaderValue<double> rho          (reader,"rho");
  TTreeReaderValue<vector<double> > jetpt        (reader,"combinejetpt");
  TTreeReaderValue<vector<double> > jeteta       (reader,"combinejeteta");
  TTreeReaderValue<vector<double> > jetphi       (reader,"combinejetphi");
  TTreeReaderValue<vector<double> > jetCHfrac    (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<double> > jetNHfrac    (reader,"combinejetNHfrac");
  TTreeReaderValue<double> wzpt     (reader,"wzpt");
  TTreeReaderValue<double> t1pfmet  (reader,"t1pfmet");
  TTreeReaderValue<double> t1pfmetphi  (reader,"t1pfmetphi");
  TTreeReaderValue<double> t1phmet     (reader,"t1phmet");
  TTreeReaderValue<double> phpt   (reader,"phpt");
  TTreeReaderValue<double> pheta  (reader,"pheta");
  TTreeReaderValue<int>    phidm  (reader,"phidm");

  TFile triggerfile_SinglePhoton("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/photonTriggerEfficiency_jetHT.root");
  TEfficiency*  triggerphoton            = (TEfficiency*)triggerfile_SinglePhoton.Get("efficiency_photon_pfht");
  TGraphAsymmErrors* triggerphoton_graph = triggerphoton->CreateGraph();
  
  TFile pufile ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt_2.60.root");
  TH1* puhist = (TH1*) pufile.Get("puhist");  

  TFile kffile ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/uncertainties_EWK_24bins.root");
  TH1*  alohist  = (TH1*) kffile.Get("GJets_LO/inv_pt_G");
  TH1* aewkhist  = (TH1*) kffile.Get("EWKcorr/photon");
  aewkhist->Divide(alohist);

  vector<TH1*> khist; khist.push_back(aewkhist);
  while(reader.Next()){    
    //basic selections for trigger, filters, bjet, taus and lepton vetoes    
    unsigned int hlt = *hltp165 + *hltp175;
    // trigger for data
    if(hlt == 0) continue;
    // met filters
    if(*fhbhe   == 0) continue;
    if(*fhbiso  == 0) continue;
    if(*fcsc    == 0) continue;
    if(*fcsct   == 0) continue;
    if(*feeb    == 0) continue;
    if(*fetp    == 0) continue;
    if(*fvtx    == 0) continue;    
    // no bjets
    if(*nbjets  >= 1) continue;    
    // no tau-jet
    if(*ntaus   >= 1) continue;
    // no muons
    if(*nmuon   >= 1)   continue;
    // no electrons
    if(*nelectron >= 1) continue;
    // mono-jet cut to select gamma+jets
    if(jetpt->size() == 0) continue;
    if(jetpt->at(0)  < 100)    continue;
    if(fabs(jeteta->at(0)) > 2.5) continue;
    if(jetCHfrac->at(0) < 0.1)  continue;
    if(jetNHfrac->at(0) > 0.8)  continue;
    // photon pt, eta requirements ... only one reco photon with pT > 35 GeV
    if(phoPt->size() == 0)      continue;
    if(phoPt->at(0)  < 175.)    continue;    
    if(fabs(phoSCEta->at(0)) > 1.479) continue;
    if(phoPt->size() != 1) continue;
    // Select the right pt bin
    if(phoPt->at(0) < ptMin or phoPt->at(0) > ptMax) continue;
    // pileup selection
    if(*nvtx <= nvtxMin or *nvtx > nvtxMax) continue;
    

    if(addRecoil){
      // recoil using the leading reco photon
      double metx = *t1pfmet*cos(*t1pfmetphi);
      double mety = *t1pfmet*sin(*t1pfmetphi);
      metx += phoPt->at(0)*cos(phoPhi->at(0));
      mety += phoPt->at(0)*sin(phoPhi->at(0));
      
      double recoil = sqrt(metx*metx+mety*mety);
      if(recoil < 200) continue;
      double recoilphi = atan(mety/metx); 
      
      double minDeltaPhi = 99;
      for(int ijet = 0; ijet < jetpt->size() and ijet < 4; ijet++){
	float deltaPhi = fabs(jetphi->at(ijet)-recoilphi);
	if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi()-deltaPhi;      
	// photon - jet dPhi
	float phoJetDphi = fabs(phoPhi->at(0)-jetphi->at(ijet));
	if(phoJetDphi > TMath::Pi()) phoJetDphi = 2*TMath::Pi()-phoJetDphi;	
	if(sqrt(phoJetDphi*phoJetDphi+fabs(jeteta->at(ijet)-phoEta->at(0))*fabs(jeteta->at(ijet)-phoEta->at(0))) < 0.4) continue;
	// dphi recoul jet
	if(deltaPhi < minDeltaPhi)
	  minDeltaPhi = deltaPhi;
      }
      //if(minDeltaPhi < 0.5) continue;
    }


    // in case of data
    if(isData){            
      // EMIso + energy and eta
      if(phoHoverE->at(0) < 0.05 and
	 phoSigmaIetaIeta->at(0) < 0.0102 and
	 phoElectronVeto->at(0) and
	 max(0.,phoCHIso->at(0)-*rho*getEffectiveArea("CH",fabs(phoSCEta->at(0)),true)) < 3.32 and
	 max(0.,phoNHIso->at(0)-*rho*getEffectiveArea("NH",fabs(phoSCEta->at(0)),true)) < 1.92+phoPt->at(0)*0.014+phoPt->at(0)*phoPt->at(0)*1.9e-05){

	photonEMIso_data.at(ipt)->Fill(max(0.,phoEMIso->at(0)-*rho*getEffectiveArea("EM",fabs(phoSCEta->at(0)),true)));
	photonEta_data.at(ipt)->Fill(phoEta->at(0));
	photonSCEta_data.at(ipt)->Fill(phoSCEta->at(0));
	photonSCEnergy_data.at(ipt)->Fill(phoSCEnergy->at(0));
	photonSCRawEnergy_data.at(ipt)->Fill(phoSCRawEnergy->at(0));
	photonPt_data_EMIso.at(ipt)->Fill(phoPt->at(0));
      }
      // NH iso
      if(phoHoverE->at(0) < 0.05 and
	 phoSigmaIetaIeta->at(0) < 0.0102 and
	 phoElectronVeto->at(0) and 
	 max(0.,phoCHIso->at(0)-*rho*getEffectiveArea("CH",fabs(phoSCEta->at(0)),true)) < 3.32 and
	 max(0.,phoEMIso->at(0)-*rho*getEffectiveArea("EM",fabs(phoSCEta->at(0)),true)) < 0.81 + phoPt->at(0)*0.0053){
	
	photonNHIso_data.at(ipt)->Fill(max(0.,phoNHIso->at(0)-*rho*getEffectiveArea("NH",fabs(phoSCEta->at(0)),true)));
	photonPt_data_NHIso.at(ipt)->Fill(phoPt->at(0));
      }
      // CH iso
      if(phoHoverE->at(0) < 0.05 and
	 phoSigmaIetaIeta->at(0) < 0.0102 and
	 phoElectronVeto->at(0) and
	 max(0.,phoNHIso->at(0)-*rho*getEffectiveArea("NH",fabs(phoSCEta->at(0)),true)) < 1.92+phoPt->at(0)*0.014+phoPt->at(0)*phoPt->at(0)*1.9e-05 and
	 max(0.,phoEMIso->at(0)-*rho*getEffectiveArea("EM",fabs(phoSCEta->at(0)),true)) < 0.81 + phoPt->at(0)*0.0053){
	photonCHIso_data.at(ipt)->Fill(max(0.,phoCHIso->at(0)-*rho*getEffectiveArea("CH",fabs(phoSCEta->at(0)),true)));
	photonPt_data_CHIso.at(ipt)->Fill(phoPt->at(0));
      }
      
      
      // HoverE
      if(phoSigmaIetaIeta->at(0) < 0.0102 and
	 phoElectronVeto->at(0) and
	 max(0.,phoNHIso->at(0)-*rho*getEffectiveArea("NH",fabs(phoSCEta->at(0)),true)) < 1.92+phoPt->at(0)*0.014+phoPt->at(0)*phoPt->at(0)*1.9e-05 and
	 max(0.,phoEMIso->at(0)-*rho*getEffectiveArea("EM",fabs(phoSCEta->at(0)),true)) < 0.81 + phoPt->at(0)*0.0053 and
	 max(0.,phoCHIso->at(0)-*rho*getEffectiveArea("CH",fabs(phoSCEta->at(0)),true)) < 3.32 ){
	photonHoverE_data.at(ipt)->Fill(phoHoverE->at(0));
	photonPt_data_HoE.at(ipt)->Fill(phoPt->at(0));
      }
      
      // sigma ieta ieta
      if(phoHoverE->at(0) < 0.05 and
	 phoElectronVeto->at(0) and
	 max(0.,phoNHIso->at(0)-*rho*getEffectiveArea("NH",fabs(phoSCEta->at(0)),true)) < 1.92+phoPt->at(0)*0.014+phoPt->at(0)*phoPt->at(0)*1.9e-05 and
	 max(0.,phoEMIso->at(0)-*rho*getEffectiveArea("EM",fabs(phoSCEta->at(0)),true)) < 0.81 + phoPt->at(0)*0.0053 and
	 max(0.,phoCHIso->at(0)-*rho*getEffectiveArea("CH",fabs(phoSCEta->at(0)),true)) < 3.32 ){
	photonSigmaIetaIeta_data.at(ipt)->Fill(phoSigmaIetaIeta->at(0));
	photonPt_data_SIetaIeta.at(ipt)->Fill(phoPt->at(0));
      }
	
      //electorn veto
      if(phoHoverE->at(0) < 0.05 and
	 phoSigmaIetaIeta->at(0) < 0.0102  and
	 max(0.,phoNHIso->at(0)-*rho*getEffectiveArea("NH",fabs(phoSCEta->at(0)),true)) < 1.92+phoPt->at(0)*0.014+phoPt->at(0)*phoPt->at(0)*1.9e-05 and
	 max(0.,phoEMIso->at(0)-*rho*getEffectiveArea("EM",fabs(phoSCEta->at(0)),true)) < 0.81 + phoPt->at(0)*0.0053 and
	 max(0.,phoCHIso->at(0)-*rho*getEffectiveArea("CH",fabs(phoSCEta->at(0)),true)) < 3.32 ){
	
	photonElectronVeto_data.at(ipt)->Fill(phoElectronVeto->at(0));
	photonPt_data_Eveto.at(ipt)->Fill(phoPt->at(0));
      }
    }
    else{
      
      double pwgt  = 1;
      pwgt *= triggerphoton_graph->Eval(min(phoPt->at(0),triggerphoton_graph->GetXaxis()->GetXmax()));
      if(*nvtx < 40)
      	pwgt *= puhist->GetBinContent(puhist->FindBin(*nvtx));
      //Gen level info --> NLO re-weight                                                                                                                                    
      Double_t kwgt = 1.0;
      double genpt = *wzpt;
      if (*wzpt < 150. ) genpt = 150.;
      if (*wzpt > 1000.) genpt = 999.;
      for (size_t i = 0; i < khist.size(); i++) {
      	if (khist[i]) {
	  kwgt *= khist[i]->GetBinContent(khist[i]->FindBin(genpt));
      	}
      }      
      pwgt *= (*wgt)*(*xsec)*luminosity*kwgt/(wgtsum);

      // EMIso + energy and eta
      if(phoHoverE->at(0) < 0.05 and
	 phoSigmaIetaIeta->at(0) < 0.0102 and
	 phoElectronVeto->at(0) and
	 max(0.,phoCHIso->at(0)-*rho*getEffectiveArea("CH",fabs(phoSCEta->at(0)),true)) < 3.32 and
	 max(0.,phoNHIso->at(0)-*rho*getEffectiveArea("NH",fabs(phoSCEta->at(0)),true)) < 1.92+phoPt->at(0)*0.014+phoPt->at(0)*phoPt->at(0)*1.9e-05){

	photonEMIso_mc.at(ipt)->Fill(max(0.,phoEMIso->at(0)-*rho*getEffectiveArea("EM",fabs(phoSCEta->at(0)),true)),pwgt);
	photonEta_mc.at(ipt)->Fill(phoEta->at(0),pwgt);
	photonSCEta_mc.at(ipt)->Fill(phoSCEta->at(0),pwgt);
	photonSCEnergy_mc.at(ipt)->Fill(phoSCEnergy->at(0),pwgt);
	photonSCRawEnergy_mc.at(ipt)->Fill(phoSCRawEnergy->at(0),pwgt);
	photonPt_mc_EMIso.at(ipt)->Fill(phoPt->at(0),pwgt);

      }
      // NH iso
      if(phoHoverE->at(0) < 0.05 and
	 phoSigmaIetaIeta->at(0) < 0.0102 and
	 photonElectronVeto_mc.at(ipt) and 
	 max(0.,phoCHIso->at(0)-*rho*getEffectiveArea("CH",fabs(phoSCEta->at(0)),true)) < 3.32 and
	 max(0.,phoEMIso->at(0)-*rho*getEffectiveArea("EM",fabs(phoSCEta->at(0)),true)) < 0.81 + phoPt->at(0)*0.0053){
	photonNHIso_mc.at(ipt)->Fill(max(0.,phoNHIso->at(0)-*rho*getEffectiveArea("NH",fabs(phoSCEta->at(0)),true)),pwgt);
	photonPt_mc_NHIso.at(ipt)->Fill(phoPt->at(0),pwgt);
      }
      // CH iso
      if(phoHoverE->at(0) < 0.05 and
	 phoSigmaIetaIeta->at(0) < 0.0102 and
	 photonElectronVeto_mc.at(ipt) and
	 max(0.,phoNHIso->at(0)-*rho*getEffectiveArea("NH",fabs(phoSCEta->at(0)),true)) < 1.92+phoPt->at(0)*0.014+phoPt->at(0)*phoPt->at(0)*1.9e-05 and
	 max(0.,phoEMIso->at(0)-*rho*getEffectiveArea("EM",fabs(phoSCEta->at(0)),true)) < 0.81 + phoPt->at(0)*0.0053){
	photonCHIso_mc.at(ipt)->Fill(max(0.,phoCHIso->at(0)-*rho*getEffectiveArea("CH",fabs(phoSCEta->at(0)),true)),pwgt);
	photonPt_mc_CHIso.at(ipt)->Fill(phoPt->at(0),pwgt);
      }
      
      
      // HoverE
      if(phoSigmaIetaIeta->at(0) < 0.0102 and
	 photonElectronVeto_mc.at(ipt) and
	 max(0.,phoNHIso->at(0)-*rho*getEffectiveArea("NH",fabs(phoSCEta->at(0)),true)) < 1.92+phoPt->at(0)*0.014+phoPt->at(0)*phoPt->at(0)*1.9e-05 and
	 max(0.,phoEMIso->at(0)-*rho*getEffectiveArea("EM",fabs(phoSCEta->at(0)),true)) < 0.81 + phoPt->at(0)*0.0053 and
	 max(0.,phoCHIso->at(0)-*rho*getEffectiveArea("CH",fabs(phoSCEta->at(0)),true)) < 3.32 ){
	photonHoverE_mc.at(ipt)->Fill(phoHoverE->at(0),pwgt);
	photonPt_mc_HoE.at(ipt)->Fill(phoPt->at(0),pwgt);
      }

      // sigma ieta ieta
      if(phoHoverE->at(0) < 0.05 and
	 photonElectronVeto_mc.at(ipt) and
	 max(0.,phoNHIso->at(0)-*rho*getEffectiveArea("NH",fabs(phoSCEta->at(0)),true)) < 1.92+phoPt->at(0)*0.014+phoPt->at(0)*phoPt->at(0)*1.9e-05 and
	 max(0.,phoEMIso->at(0)-*rho*getEffectiveArea("EM",fabs(phoSCEta->at(0)),true)) < 0.81 + phoPt->at(0)*0.0053 and
	 max(0.,phoCHIso->at(0)-*rho*getEffectiveArea("CH",fabs(phoSCEta->at(0)),true)) < 3.32 ){
	photonSigmaIetaIeta_mc.at(ipt)->Fill(phoSigmaIetaIeta->at(0),pwgt);
	photonPt_mc_SIetaIeta.at(ipt)->Fill(phoPt->at(0),pwgt);
      }
      
      //electorn veto
      if(phoHoverE->at(0) < 0.05 and
	 phoSigmaIetaIeta->at(0) < 0.0102  and
	 max(0.,phoNHIso->at(0)-*rho*getEffectiveArea("NH",fabs(phoSCEta->at(0)),true)) < 1.92+phoPt->at(0)*0.014+phoPt->at(0)*phoPt->at(0)*1.9e-05 and
	 max(0.,phoEMIso->at(0)-*rho*getEffectiveArea("EM",fabs(phoSCEta->at(0)),true)) < 0.81 + phoPt->at(0)*0.0053 and
	 max(0.,phoCHIso->at(0)-*rho*getEffectiveArea("CH",fabs(phoSCEta->at(0)),true)) < 3.32 ){
	photonElectronVeto_mc.at(ipt)->Fill(phoElectronVeto->at(0),pwgt);
	photonPt_mc_Eveto.at(ipt)->Fill(phoPt->at(0),pwgt);
      }
    }     
  }    
  triggerfile_SinglePhoton.Close();
  pufile.Close();
  kffile.Close();
}


void makePhotonIDVariableComparison(string inputDIRData, string inputDIRMC, string ouputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+ouputDIR).c_str());
  cout<<"Start chain files"<<endl;
  TChain* chain_data   = new TChain("tree/tree");
  chain_data->Add((inputDIRData+"/*root").c_str());
  TChain* chain_mc_1   = new TChain("tree/tree");
  chain_mc_1->Add((inputDIRMC+"/tree_GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root").c_str());
  TChain* chain_mc_2   = new TChain("tree/tree");
  chain_mc_2->Add((inputDIRMC+"/tree_GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root").c_str());
  TChain* chain_mc_3   = new TChain("tree/tree");
  chain_mc_3->Add((inputDIRMC+"/tree_GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root").c_str());
  TChain* chain_mc_4   = new TChain("tree/tree");
  chain_mc_4->Add((inputDIRMC+"/tree_GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root").c_str());
  TChain* chain_mc_5   = new TChain("tree/tree");
  chain_mc_5->Add((inputDIRMC+"/tree_GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root").c_str());

  TChain* chain_genmc_1  = new TChain("gentree/gentree");
  chain_genmc_1->Add((inputDIRMC+"/tree_GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root").c_str());
  TChain* chain_genmc_2  = new TChain("gentree/gentree");
  chain_genmc_2->Add((inputDIRMC+"/tree_GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root").c_str());
  TChain* chain_genmc_3  = new TChain("gentree/gentree");
  chain_genmc_3->Add((inputDIRMC+"/tree_GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root").c_str());
  TChain* chain_genmc_4  = new TChain("gentree/gentree");
  chain_genmc_4->Add((inputDIRMC+"/tree_GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root").c_str());
  TChain* chain_genmc_5  = new TChain("gentree/gentree");
  chain_genmc_5->Add((inputDIRMC+"/tree_GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root").c_str());
  cout<<"End to chain files "<<endl;
  
  cout<<"Calculate mc weight "<<endl;
  TTreeReader treeReader_genmc_1(chain_genmc_1);
  TTreeReaderValue<double> wgtsign_genmc_1 (treeReader_genmc_1,"wgtsign");
  double weightsum_genmc_1 = 0.;
  while(treeReader_genmc_1.Next()){
    weightsum_genmc_1 += (*wgtsign_genmc_1);
  }

  TTreeReader treeReader_genmc_2(chain_genmc_2);
  TTreeReaderValue<double> wgtsign_genmc_2 (treeReader_genmc_2,"wgtsign");
  double weightsum_genmc_2 = 0.;
  while(treeReader_genmc_2.Next()){
    weightsum_genmc_2 += (*wgtsign_genmc_2);
  }

  TTreeReader treeReader_genmc_3(chain_genmc_3);
  TTreeReaderValue<double> wgtsign_genmc_3 (treeReader_genmc_3,"wgtsign");
  double weightsum_genmc_3 = 0.;
  while(treeReader_genmc_3.Next()){
    weightsum_genmc_3 += (*wgtsign_genmc_3);
  }

  TTreeReader treeReader_genmc_4(chain_genmc_4);
  TTreeReaderValue<double> wgtsign_genmc_4 (treeReader_genmc_4,"wgtsign");
  double weightsum_genmc_4 = 0.;
  while(treeReader_genmc_4.Next()){
    weightsum_genmc_4 += (*wgtsign_genmc_4);
  }

  TTreeReader treeReader_genmc_5(chain_genmc_5);
  TTreeReaderValue<double> wgtsign_genmc_5 (treeReader_genmc_5,"wgtsign");
  double weightsum_genmc_5 = 0.;
  while(treeReader_genmc_5.Next()){
    weightsum_genmc_5 += (*wgtsign_genmc_5);
  }

  cout<<"Create histogram "<<endl;
  // create histograms
  for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++){
    photonEta_data.push_back(new TH1F(Form("photonEta_data_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",20,-1.5,1.5));
    photonSCEta_data.push_back(new TH1F(Form("photonSCEta_data_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",20,-1.5,1.5));
    photonSCEnergy_data.push_back(new TH1F(Form("photonSCEnergy_data_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",50,100,1300));
    photonSCRawEnergy_data.push_back(new TH1F(Form("photonSCRawEnergy_data_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",50,100,1300));
    photonHoverE_data.push_back(new TH1F(Form("photonHoverE_data_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",30,0,0.1));
    photonSigmaIetaIeta_data.push_back(new TH1F(Form("photonSigmaIetaIeta_data_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,0,0.03));
    photonElectronVeto_data.push_back(new TH1F(Form("photonElectronVeto_data_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",2,0,2));
    photonCHIso_data.push_back(new TH1F(Form("photonCHIso_data_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",30,0,6));
    photonNHIso_data.push_back(new TH1F(Form("photonNHIso_data_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",30,0,6));
    photonEMIso_data.push_back(new TH1F(Form("photonEMIso_data_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",30,0,6));

    photonPt_data_HoE.push_back(new TH1F(Form("photonPt_data_HoE_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,175.,1000.));
    photonPt_data_SIetaIeta.push_back(new TH1F(Form("photonPt_data_SIetaIeta_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,175.,1000.));
    photonPt_data_Eveto.push_back(new TH1F(Form("photonPt_data_Eveto_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,175.,1000.));
    photonPt_data_CHIso.push_back(new TH1F(Form("photonPt_data_CHIso_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,175.,1000.));
    photonPt_data_NHIso.push_back(new TH1F(Form("photonPt_data_NHIso_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,175.,1000.));
    photonPt_data_EMIso.push_back(new TH1F(Form("photonPt_data_EMIso_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,175.,1000.));

    photonEta_mc.push_back(new TH1F(Form("photonEta_mc_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",20,-1.5,1.5));
    photonSCEta_mc.push_back(new TH1F(Form("photonSCEta_mc_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",20,-1.5,1.5));
    photonSCEnergy_mc.push_back(new TH1F(Form("photonSCEnergy_mc_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",50,100,1300));
    photonSCRawEnergy_mc.push_back(new TH1F(Form("photonSCRawEnergy_mc_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",50,100,1300));
    photonHoverE_mc.push_back(new TH1F(Form("photonHoverE_mc_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",30,0,0.1));
    photonSigmaIetaIeta_mc.push_back(new TH1F(Form("photonSigmaIetaIeta_mc_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,0,0.03));
    photonElectronVeto_mc.push_back(new TH1F(Form("photonElectronVeto_mc_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",2,0,2));
    photonCHIso_mc.push_back(new TH1F(Form("photonCHIso_mc_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",30,0,6));
    photonNHIso_mc.push_back(new TH1F(Form("photonNHIso_mc_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",30,0,6));
    photonEMIso_mc.push_back(new TH1F(Form("photonEMIso_mc_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",30,0,6));

    photonPt_mc_HoE.push_back(new TH1F(Form("photonPt_mc_HoE_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,175.,1000.));
    photonPt_mc_SIetaIeta.push_back(new TH1F(Form("photonPt_mc_SIetaIeta_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,175.,1000.));
    photonPt_mc_Eveto.push_back(new TH1F(Form("photonPt_mc_Eveto_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,175.,1000.));
    photonPt_mc_CHIso.push_back(new TH1F(Form("photonPt_mc_CHIso_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,175.,1000.));
    photonPt_mc_NHIso.push_back(new TH1F(Form("photonPt_mc_NHIso_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,175.,1000.));
    photonPt_mc_EMIso.push_back(new TH1F(Form("photonPt_mc_EMIso_%d_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",40,175.,1000.));

  }
  
  // loop on the event
  cout<<"Start pt bin analysis "<<endl;
  for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++){
    cout<<"Pt Bin Data : "<<ptBin.at(ipt)<<" : "<<ptBin.at(ipt+1)<<endl;
    fillHistograms(chain_data,ptBin.at(ipt),ptBin.at(ipt+1),ipt,true,1,true);
    cout<<"Pt Bin MC   : "<<ptBin.at(ipt)<<" : "<<ptBin.at(ipt+1)<<endl;
    fillHistograms(chain_mc_1,ptBin.at(ipt),ptBin.at(ipt+1),ipt,false,weightsum_genmc_1,true);
    fillHistograms(chain_mc_2,ptBin.at(ipt),ptBin.at(ipt+1),ipt,false,weightsum_genmc_2,true);
    fillHistograms(chain_mc_3,ptBin.at(ipt),ptBin.at(ipt+1),ipt,false,weightsum_genmc_3,true);
    fillHistograms(chain_mc_4,ptBin.at(ipt),ptBin.at(ipt+1),ipt,false,weightsum_genmc_4,true);
    fillHistograms(chain_mc_5,ptBin.at(ipt),ptBin.at(ipt+1),ipt,false,weightsum_genmc_5,true);
  }
  
  TCanvas* canvas = new TCanvas("canvas","",600,700);
  cout<<"Start plotting "<<endl;
  // plot histograms
  for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++){
    plotHistogram(photonEta_data.at(ipt),photonEta_mc.at(ipt),canvas,ouputDIR,"#eta^{#gamma}");
    plotHistogram(photonSCEta_data.at(ipt),photonSCEta_mc.at(ipt),canvas,ouputDIR,"#eta_{SC}^{#gamma}");
    plotHistogram(photonSCEnergy_data.at(ipt),photonSCEnergy_mc.at(ipt),canvas,ouputDIR,"#E_{SC}^{#gamma} [GeV]");
    plotHistogram(photonSCRawEnergy_data.at(ipt),photonSCRawEnergy_mc.at(ipt),canvas,ouputDIR,"E_{SC,raw}^{#gamma} [GeV]");
    plotHistogram(photonHoverE_data.at(ipt),photonHoverE_mc.at(ipt),canvas,ouputDIR,"H/E");
    plotHistogram(photonSigmaIetaIeta_data.at(ipt),photonSigmaIetaIeta_mc.at(ipt),canvas,ouputDIR,"#sigma_{i#eta i#eta}");
    plotHistogram(photonElectronVeto_data.at(ipt),photonElectronVeto_mc.at(ipt),canvas,ouputDIR,"Electron Veto");
    plotHistogram(photonCHIso_data.at(ipt),photonCHIso_mc.at(ipt),canvas,ouputDIR,"CH Iso [GeV]");
    plotHistogram(photonNHIso_data.at(ipt),photonNHIso_mc.at(ipt),canvas,ouputDIR,"NH Iso [GeV]");
    plotHistogram(photonEMIso_data.at(ipt),photonEMIso_mc.at(ipt),canvas,ouputDIR,"EM Iso [GeV]");

    // Plot Pt
    plotHistogram(photonPt_data_HoE.at(ipt),photonPt_mc_HoE.at(ipt),canvas,ouputDIR,"p_{T}^{#gamma} (GeV)");
    plotHistogram(photonPt_data_SIetaIeta.at(ipt),photonPt_mc_SIetaIeta.at(ipt),canvas,ouputDIR,"p_{T}^{#gamma} (GeV)");
    plotHistogram(photonPt_data_Eveto.at(ipt),photonPt_mc_Eveto.at(ipt),canvas,ouputDIR,"p_{T}^{#gamma} (GeV)");
    plotHistogram(photonPt_data_CHIso.at(ipt),photonPt_mc_CHIso.at(ipt),canvas,ouputDIR,"p_{T}^{#gamma} (GeV)");
    plotHistogram(photonPt_data_EMIso.at(ipt),photonPt_mc_EMIso.at(ipt),canvas,ouputDIR,"p_{T}^{#gamma} (GeV)");
    plotHistogram(photonPt_data_NHIso.at(ipt),photonPt_mc_NHIso.at(ipt),canvas,ouputDIR,"p_{T}^{#gamma} (GeV)");
  }
  cout<<"End plotting "<<endl;
}

