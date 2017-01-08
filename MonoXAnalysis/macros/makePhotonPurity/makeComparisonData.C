#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

//////////////
//////////////
//////////////

static float luminosity = 36.2;

TH1F* photonpt_histo  = new TH1F("photonpt","", 50,-100,100);
TH1F* photoneta_histo = new TH1F("photoneta","",50,-0.15,0.15);
TH1F* photonphi_histo = new TH1F("photonphi","",50,-0.06,0.06);
TH1F* photonPHIso_histo = new TH1F("photonPHIso","",120,-5.0,5.0);
TH1F* photonCHIso_histo = new TH1F("photonCHIso","",120,-7.5,7.5);
TH1F* photonNHIso_histo = new TH1F("photonNHIso","",120,-5.0,5.0);
TH1F* photonHoE_histo = new TH1F("photonHoE","",50,-0.1,0.1);
TH1F* photonNum_histo = new TH1F("nphotons","",3,-1.5,1.5);
TH1F* photonEveto_histo = new TH1F("photonEveto","",3,-1.5,1.5);
TH1F* photonSieie_histo = new TH1F("photonSieie","",50,-0.025,0.025);

TH2F* photonpt_histo_2D  = new TH2F("photonpt_2D","", 25,175,1000,25,175,1000);
TH2F* photoneta_histo_2D = new TH2F("photoneta_2D","",25,-1.5,1.5,25,-1.5,1.5);
TH2F* photonphi_histo_2D = new TH2F("photonphi_2D","",25,-3.14,3.14,25,-3.14,3.14);
TH2F* photonPHIso_histo_2D = new TH2F("photonPHIso_2D","",60,0,25,60,0,25);
TH2F* photonCHIso_histo_2D = new TH2F("photonCHIso_2D","",60,0,35,60,0,35);
TH2F* photonNHIso_histo_2D = new TH2F("photonNHIso_2D","",60,0,25,60,0,25);
TH2F* photonHoE_histo_2D = new TH2F("photonHoE_2D","",25,0,0.06,25,0,0.06);
TH2F* photonSieie_histo_2D = new TH2F("photonSieie_2D","",50,0,0.03,50,0,0.03);

TH1F* photonpt_histo_only1  = new TH1F("photonpt_only1","", 30,175,1000);
TH1F* photoneta_histo_only1 = new TH1F("photoneta_only1","",30,-1.5,1.5);
TH1F* photonphi_histo_only1 = new TH1F("photonphi_only1","",30,-3.14,3.14);
TH1F* photonPHIso_histo_only1 = new TH1F("photonPHIso_only1","",50,0,25);
TH1F* photonCHIso_histo_only1 = new TH1F("photonCHIso_only1","",50,0,35);
TH1F* photonNHIso_histo_only1 = new TH1F("photonNHIso_only1","",50,0,25);
TH1F* photonHoE_histo_only1 = new TH1F("photonHoE_only1","",30,0,0.06);
TH1F* photonNum_histo_only1 = new TH1F("nphotons_only1","",3,-1.5,1.5);
TH1F* photonEveto_histo_only1 = new TH1F("photonEveto_only1","",3,-1.5,1.5);
TH1F* photonSieie_histo_only1 = new TH1F("photonSieie_only1","",50,0,0.03);

TH1F* photonpt_histo_only2  = new TH1F("photonpt_only2","", 30,175,1000);
TH1F* photoneta_histo_only2 = new TH1F("photoneta_only2","",30,-1.5,1.5);
TH1F* photonphi_histo_only2 = new TH1F("photonphi_only2","",30,-3.14,3.14);
TH1F* photonPHIso_histo_only2 = new TH1F("photonPHIso_only2","",50,0,25);
TH1F* photonCHIso_histo_only2 = new TH1F("photonCHIso_only2","",50,0,35);
TH1F* photonNHIso_histo_only2 = new TH1F("photonNHIso_only2","",50,0,25);
TH1F* photonHoE_histo_only2 = new TH1F("photonHoE_only2","",30,0,0.06);
TH1F* photonNum_histo_only2 = new TH1F("nphotons_only2","",3,-1.5,1.5);
TH1F* photonEveto_histo_only2 = new TH1F("photonEveto_only2","",3,-1.5,1.5);
TH1F* photonSieie_histo_only2 = new TH1F("photonSieie_only2","",50,0,0.03);

//////////////
void drawPlot(TCanvas* pad1, TH1* histo_1, string observable, string outputDIR, bool logScale = true){

  gStyle->SetOptStat(111111111);
 
  pad1->cd();
  histo_1->GetYaxis()->SetLabelSize(0.035);
  histo_1->GetYaxis()->SetTitleSize(0.042);
  histo_1->GetXaxis()->SetLabelSize(0.035);
  histo_1->GetXaxis()->SetTitleSize(0.042);
  histo_1->GetYaxis()->SetTitleOffset(1.35);
  histo_1->GetYaxis()->SetTitle("Events");
  histo_1->GetXaxis()->SetTitle(observable.c_str());
  if(histo_1->GetMinimum() != 0 and logScale)
    histo_1->GetYaxis()->SetRangeUser(histo_1->GetMinimum(),histo_1->GetMaximum()*100);
  else if(logScale)
    histo_1->GetYaxis()->SetRangeUser(0.001,histo_1->GetMaximum()*100);
  else
    histo_1->GetYaxis()->SetRangeUser(0,histo_1->GetMaximum()*1.5);

  histo_1->SetLineColor(kBlack);
  histo_1->SetLineWidth(2);  
  
  histo_1->Draw("HIST");
  CMS_lumi(pad1,Form("%.1f",luminosity));  
  histo_1->Draw("HIST SAME");

  if(logScale)
    pad1->SetLogy();

  pad1->SaveAs((outputDIR+"/"+string(histo_1->GetName())+"_comparison.png").c_str(),"png");
  pad1->SaveAs((outputDIR+"/"+string(histo_1->GetName())+"_comparison.pdf").c_str(),"pdf");

  return ;
}

//////////////
void drawPlot2D(TCanvas* pad1, TH2* histo_1, string observable, string outputDIR, bool logScale = true){

  gStyle->SetOptStat(0);
  pad1->SetLogy(0);

  pad1->SetRightMargin(0.15);
  pad1->SetLeftMargin(0.12);

  pad1->cd();
  histo_1->GetYaxis()->SetLabelSize(0.035);
  histo_1->GetYaxis()->SetTitleSize(0.042);
  histo_1->GetXaxis()->SetLabelSize(0.035);
  histo_1->GetXaxis()->SetTitleSize(0.042);
  histo_1->GetYaxis()->SetTitleOffset(1.35);
  histo_1->GetYaxis()->SetTitle((observable+" Prompt").c_str());
  histo_1->GetXaxis()->SetTitle((observable+" ReReco").c_str());

  histo_1->Draw("COLZ");
  CMS_lumi(pad1,Form("%.1f",luminosity));  

  if(logScale)
    pad1->SetLogz();

  pad1->SaveAs((outputDIR+"/"+string(histo_1->GetName())+"_comparison.png").c_str(),"png");
  pad1->SaveAs((outputDIR+"/"+string(histo_1->GetName())+"_comparison.pdf").c_str(),"pdf");

  return ;
}

//////////////
void makeComparisonData(string inputDIR_1, string inputDIR_2, string outputDIR, bool useOnlyICHEP, bool doOppositeMatching, bool applyPhotonID){

  if(useOnlyICHEP)
    luminosity = 12.9;
  else
    luminosity = 36.2;
    
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  gROOT->SetBatch();

  photonpt_histo->Sumw2();
  photoneta_histo->Sumw2();
  photonphi_histo->Sumw2();
  photonPHIso_histo->Sumw2();
  photonCHIso_histo->Sumw2();
  photonNHIso_histo->Sumw2();
  photonHoE_histo->Sumw2();
  photonNum_histo->Sumw2();
  photonEveto_histo->Sumw2();
  photonSieie_histo->Sumw2();
  
  photonpt_histo_2D->Sumw2();
  photoneta_histo_2D->Sumw2();
  photonphi_histo_2D->Sumw2();
  photonPHIso_histo_2D->Sumw2();
  photonCHIso_histo_2D->Sumw2();
  photonNHIso_histo_2D->Sumw2();
  photonHoE_histo_2D->Sumw2();
  photonSieie_histo_2D->Sumw2();

  TChain* chain_1 = new TChain("tree/tree");
  TChain* chain_2 = new TChain("tree/tree");

  chain_1->Add((inputDIR_1+"/*root").c_str());
  chain_2->Add((inputDIR_2+"/*root").c_str());

  ///////////
  cout<<"Build first chain with index "<<endl;
  unsigned int run = 0, lumi = 0, event = 0;
  unsigned int njets = 0, nmuons = 0, nelectrons = 0, nphotons = 0, ntausraw = 0, ntaus = 0, nincjets = 0, nbjets = 0;
  UChar_t hltp165 = 0, hltp175 = 0, fhbhe = 0, fhbiso = 0, fcsct = 0, feeb = 0, fetp = 0, fvtx = 0, fcsc = 0;
  float t1met = 0, t1metphi = 0;
  float phPuritypt = 0, phPurityeta = 0, phPurityphi = 0, phPurityElectronVeto = 0, phPurityPHiso = 0;
  float phPurityCHiso = 0, phPurityNHiso = 0, phPurityhoe = 0, phPuritysieie = 0, phPurityRND04PHiso = 0, phPurityRND08PHiso = 0;
  float phPurityEAEGamma = 0, rho = 0;
  float t1phmet = 0;
  vector<float> *jeteta = 0, *jetpt = 0, *jetphi = 0,*jetm = 0;

  chain_1->SetBranchStatus("*",kFALSE);
  chain_1->SetBranchStatus("run",kTRUE);
  chain_1->SetBranchStatus("lumi",kTRUE);
  chain_1->SetBranchStatus("event",kTRUE);
  chain_1->SetBranchStatus("hltphoton*",kTRUE);
  chain_1->SetBranchStatus("flag*",kTRUE);
  chain_1->SetBranchStatus("njets",kTRUE);
  chain_1->SetBranchStatus("nmuons",kTRUE);
  chain_1->SetBranchStatus("nelectrons",kTRUE);
  chain_1->SetBranchStatus("nphotons*",kTRUE);
  chain_1->SetBranchStatus("ntaus*",kTRUE);
  chain_1->SetBranchStatus("nbjetslowpt",kTRUE);
  chain_1->SetBranchStatus("combinejeteta",kTRUE);
  chain_1->SetBranchStatus("combinejetpt", kTRUE);
  chain_1->SetBranchStatus("combinejetphi",kTRUE);
  chain_1->SetBranchStatus("combinejetm",kTRUE);
  chain_1->SetBranchStatus("t1pfmet*",kTRUE);
  chain_1->SetBranchStatus("t1phmet*",kTRUE);
  chain_1->SetBranchStatus("rho",kTRUE);
  chain_1->SetBranchStatus("phPurity*",kTRUE);

  chain_1->SetBranchAddress("run",&run);
  chain_1->SetBranchAddress("lumi",&lumi);
  chain_1->SetBranchAddress("event",&event);
  chain_1->SetBranchAddress("hltphoton165",&hltp165);
  chain_1->SetBranchAddress("hltphoton175",&hltp175);
  chain_1->SetBranchAddress("flaghbhenoise",&fhbhe);
  chain_1->SetBranchAddress("flaghbheiso",&fhbiso);
  chain_1->SetBranchAddress("flagcsctight",&fcsct);
  chain_1->SetBranchAddress("flageebadsc",&feeb);
  chain_1->SetBranchAddress("flagecaltp",&fetp);
  chain_1->SetBranchAddress("flaggoodvertices",&fvtx);
  chain_1->SetBranchAddress("flagglobaltighthalo",&fcsc);
  chain_1->SetBranchAddress("njets",&njets);
  chain_1->SetBranchAddress("nmuons",&nmuons);
  chain_1->SetBranchAddress("nelectrons",&nelectrons);
  chain_1->SetBranchAddress("nphotonsPurity",&nphotons);
  chain_1->SetBranchAddress("nelectrons",&nelectrons);
  chain_1->SetBranchAddress("ntaus",&ntaus);
  chain_1->SetBranchAddress("ntausrawold",&ntausraw);
  chain_1->SetBranchAddress("nbjetslowpt",&nbjets);
  chain_1->SetBranchAddress("combinejeteta",&jeteta);
  chain_1->SetBranchAddress("combinejetpt", &jetpt);
  chain_1->SetBranchAddress("combinejetphi",&jetphi);
  chain_1->SetBranchAddress("combinejetm",&jetm);
  chain_1->SetBranchAddress("phPuritypt",&phPuritypt);
  chain_1->SetBranchAddress("pPurityheta",&phPurityeta);
  chain_1->SetBranchAddress("phPurityphi",&phPurityphi);
  chain_1->SetBranchAddress("phPurityElectronVeto",&phPurityElectronVeto);
  chain_1->SetBranchAddress("phPurityPHiso",&phPurityPHiso);
  chain_1->SetBranchAddress("phPurityCHiso",&phPurityCHiso);
  chain_1->SetBranchAddress("phPurityNHiso",&phPurityNHiso);
  chain_1->SetBranchAddress("phPurityhoe",&phPurityhoe);
  chain_1->SetBranchAddress("phPuritysieie",&phPuritysieie);
  chain_1->SetBranchAddress("phPurityRND04PHiso",&phPurityRND04PHiso);
  chain_1->SetBranchAddress("phPurityRND08PHiso",&phPurityRND08PHiso);
  chain_1->SetBranchAddress("phPurityEAEGamma",&phPurityEAEGamma);
  chain_1->SetBranchAddress("rho",&rho);
  chain_1->SetBranchAddress("t1pfmet",&t1met);
  chain_1->SetBranchAddress("t1phmet",&t1phmet);
  chain_1->SetBranchAddress("t1pfmetphi",&t1metphi);

  chain_1->BuildIndex("run","event");

  ///
  cout<<"Build second chain with index "<<endl;
  unsigned int run_alt = 0, lumi_alt = 0, event_alt = 0;
  unsigned int njets_alt = 0, nmuons_alt = 0, nelectrons_alt = 0, nphotons_alt = 0, ntausraw_alt = 0, ntaus_alt = 0, nincjets_alt = 0, nbjets_alt = 0;
  UChar_t hltp165_alt = 0, hltp175_alt = 0, fhbhe_alt = 0, fhbiso_alt = 0, fcsct_alt = 0, feeb_alt = 0, fetp_alt = 0, fvtx_alt = 0, fcsc_alt = 0;
  double t1met_alt = 0, t1metphi_alt = 0;
  double phPuritypt_alt = 0, phPurityeta_alt = 0, phPurityphi_alt = 0, phPurityElectronVeto_alt = 0, phPurityPHiso_alt = 0;
  double phPurityCHiso_alt = 0, phPurityNHiso_alt = 0, phPurityhoe_alt = 0, phPuritysieie_alt = 0, phPurityRND04PHiso_alt = 0, phPurityRND08PHiso_alt = 0;
  double phPurityEAEGamma_alt = 0, rho_alt = 0;
  vector<double> *jeteta_alt = 0, *jetpt_alt = 0, *jetphi_alt = 0,*jetm_alt = 0;
  double t1phmet_alt = 0;

  chain_2->SetBranchStatus("*",kFALSE);
  chain_2->SetBranchStatus("run",kTRUE);
  chain_2->SetBranchStatus("lumi",kTRUE);
  chain_2->SetBranchStatus("event",kTRUE);
  chain_2->SetBranchStatus("hltphoton*",kTRUE);
  chain_2->SetBranchStatus("flag*",kTRUE);
  chain_2->SetBranchStatus("njets",kTRUE);
  chain_2->SetBranchStatus("nmuons",kTRUE);
  chain_2->SetBranchStatus("nelectrons",kTRUE);
  chain_2->SetBranchStatus("nphotons*",kTRUE);
  chain_2->SetBranchStatus("ntaus*",kTRUE);
  chain_2->SetBranchStatus("nbjetslowpt",kTRUE);
  chain_2->SetBranchStatus("combinejeteta",kTRUE);
  chain_2->SetBranchStatus("combinejetpt", kTRUE);
  chain_2->SetBranchStatus("combinejetphi",kTRUE);
  chain_2->SetBranchStatus("combinejetm",kTRUE);
  chain_2->SetBranchStatus("t1pfmet*",kTRUE);
  chain_2->SetBranchStatus("t1phmet*",kTRUE);
  chain_2->SetBranchStatus("rho",kTRUE);
  chain_2->SetBranchStatus("phPurity*",kTRUE);

  chain_2->SetBranchAddress("run",&run_alt);
  chain_2->SetBranchAddress("lumi",&lumi_alt);
  chain_2->SetBranchAddress("event",&event_alt);
  chain_2->SetBranchAddress("hltphoton165",&hltp165_alt);
  chain_2->SetBranchAddress("hltphoton175",&hltp175_alt);
  chain_2->SetBranchAddress("flaghbhenoise",&fhbhe_alt);
  chain_2->SetBranchAddress("flaghbheiso",&fhbiso_alt);
  chain_2->SetBranchAddress("flagcsctight",&fcsct_alt);
  chain_2->SetBranchAddress("flageebadsc",&feeb_alt);
  chain_2->SetBranchAddress("flagecaltp",&fetp_alt);
  chain_2->SetBranchAddress("flaggoodvertices",&fvtx_alt);
  chain_2->SetBranchAddress("flagglobaltighthalo",&fcsc_alt);
  chain_2->SetBranchAddress("njets",&njets_alt);
  chain_2->SetBranchAddress("nmuons",&nmuons_alt);
  chain_2->SetBranchAddress("nelectrons",&nelectrons_alt);
  chain_2->SetBranchAddress("nphotonsPurity",&nphotons_alt);
  chain_2->SetBranchAddress("nelectrons",&nelectrons_alt);
  chain_2->SetBranchAddress("ntaus",&ntaus_alt);
  chain_2->SetBranchAddress("ntausraw",&ntausraw_alt);
  chain_2->SetBranchAddress("nbjetslowpt",&nbjets_alt);
  chain_2->SetBranchAddress("combinejeteta",&jeteta_alt);
  chain_2->SetBranchAddress("combinejetpt", &jetpt_alt);
  chain_2->SetBranchAddress("combinejetphi",&jetphi_alt);
  chain_2->SetBranchAddress("combinejetm",&jetm_alt);
  chain_2->SetBranchAddress("phPuritypt",&phPuritypt_alt);
  chain_2->SetBranchAddress("pPurityheta",&phPurityeta_alt);
  chain_2->SetBranchAddress("phPurityphi",&phPurityphi_alt);
  chain_2->SetBranchAddress("phPurityElectronVeto",&phPurityElectronVeto_alt);
  chain_2->SetBranchAddress("phPurityPHiso",&phPurityPHiso_alt);
  chain_2->SetBranchAddress("phPurityCHiso",&phPurityCHiso_alt);
  chain_2->SetBranchAddress("phPurityNHiso",&phPurityNHiso_alt);
  chain_2->SetBranchAddress("phPurityhoe",&phPurityhoe_alt);
  chain_2->SetBranchAddress("phPuritysieie",&phPuritysieie_alt);
  chain_2->SetBranchAddress("phPurityRND04PHiso",&phPurityRND04PHiso_alt);
  chain_2->SetBranchAddress("phPurityRND08PHiso",&phPurityRND08PHiso_alt);
  chain_2->SetBranchAddress("phPurityEAEGamma",&phPurityEAEGamma_alt);
  chain_2->SetBranchAddress("rho",&rho_alt);
  chain_2->SetBranchAddress("t1pfmet",&t1met_alt);
  chain_2->SetBranchAddress("t1phmet",&t1phmet_alt);
  chain_2->SetBranchAddress("t1pfmetphi",&t1metphi_alt);
  
  chain_2->BuildIndex("run","event");

  cout<<"Number of events in chain 1 "<<chain_1->GetEntries()<<endl;
  cout<<"Number of events in chain 2 "<<chain_2->GetEntries()<<endl;

  long int nEvents = 0;
  int nPart        = 10000;
  long int nTotal  = chain_1->GetEntries();
  long int notMatchedEvents = 0;
  long int effectiveEvents  = 0;

  for(long int entry = 0; entry < chain_1->GetEntries(); entry++){
    
    int entry_status = chain_1->GetEntry(entry);
    if(entry_status <= 0) continue;
    
    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;
    
    // to select only ICHEP data
    if(useOnlyICHEP and run > 276242) continue;  

    //////// -- event selection --> everything except for Photon ID selections, i.e. H/E, Sietaieta, Isolation ..etc
    // trigger
    int hlt = 0;
    hlt = hltp165+hltp175;
    if(not hlt) continue;     
    /// apply met filters
    if(not fhbhe or not fhbiso or not fcsct or not feeb or not feeb or not fetp or not fvtx or not fcsc) continue;    
    /// vetoes
    if(ntausraw != 0) continue;
    if(nbjets   != 0) continue;
    if(nmuons   != 0) continue;
    if(nelectrons != 0) continue;
    if(nphotons == 0) continue;
    // photon candidate
    if(phPuritypt < 175) continue;
    if(fabs(phPurityeta) > 1.5) continue;

    if(applyPhotonID){
      if(phPurityhoe   > 0.0396) continue;
      if(phPuritysieie > 0.01022) continue;
      if(phPurityElectronVeto == 0) continue;
      if(phPurityNHiso > 2.725 + 0.0148*phPuritypt + 0.000017*phPuritypt*phPuritypt) continue;
    }

    // jets
    if(njets < 1) continue;
    // check overlap with leading jet --> cleaning is done wrt loose photons by default in the trees
    int ijet = 0;
    for(size_t i = 0 ; i < jetpt->size(); i++){
      float deta = fabs(jeteta->at(i)-phPurityeta);
      float dphi = fabs(phPurityphi-jetphi->at(i));
      if(dphi > TMath::Pi())
	dphi = 2*TMath::Pi()-dphi;
      if(sqrt(deta*deta+dphi*dphi) > 0.4){
	ijet = i;
	break;
      }
    }

    if(jetpt->at(ijet) < 100) continue;
    if(fabs(jeteta->at(ijet)) > 2.5) continue;

    // jet-met dphi
    int njet = 0;
    float mindphi = 99;
    for(size_t i= 0 ; i < jetpt->size(); i++){
      // check pt
      if(jetpt->at(i) < 30) continue;
      // check if overlaps with the photon candidate
      float deta = fabs(jeteta->at(i)-phPurityeta);
      float dphi = fabs(phPurityphi-jetphi->at(i));
      if(dphi > TMath::Pi())
	dphi = 2*TMath::Pi()-dphi;
      if(sqrt(deta*deta+dphi*dphi) < 0.4) continue;
      njet++;
      if(njet > 4) continue;
      // calculate px and py
      float metx = t1met*cos(t1metphi)+phPuritypt*cos(phPurityphi);
      float mety = t1met*sin(t1metphi)+phPuritypt*sin(phPurityphi);      
      TLorentzVector met;
      met.SetPxPyPzE(metx,mety,0.,sqrt(metx*metx+mety*mety));
      float dphitemp = fabs(met.Phi()-jetphi->at(ijet));
      if(dphitemp > TMath::Pi())
	dphitemp = 2*TMath::Pi()-dphitemp;
      if(dphitemp < mindphi)
	mindphi = dphitemp; 
    }
    if(mindphi < 0.5) continue;

    // access to the other tree
    entry_status = chain_2->GetEntryWithIndex(run,event);
    effectiveEvents++;

    // events in chain_1 but not in 2
    if(entry_status <= 0){
      notMatchedEvents++;
      photonpt_histo_only1->Fill(phPuritypt);
      photoneta_histo_only1->Fill(phPurityeta);
      photonphi_histo_only1->Fill(phPurityphi);
      photonPHIso_histo_only1->Fill(phPurityPHiso);
      photonCHIso_histo_only1->Fill(phPurityCHiso);
      photonNHIso_histo_only1->Fill(phPurityNHiso);
      photonHoE_histo_only1->Fill(phPurityhoe);
      photonNum_histo_only1->Fill(nphotons);
      photonEveto_histo_only1->Fill(phPurityElectronVeto);
      photonSieie_histo_only1->Fill(phPuritysieie);
    }
    else{
      // common events
      if(run != run_alt or lumi != lumi_alt or event != event_alt){
	cerr<<"Problem in associating: [run,lumi,event] = ["<<run<<","<<lumi<<","<<event<<"] with ["<<run_alt<<","<<lumi_alt<<","<<event_alt<<"]"<<endl;
	continue;
      }
      
      //////////////
      photonpt_histo->Fill(phPuritypt-phPuritypt_alt);
      photoneta_histo->Fill(phPurityeta-phPurityeta_alt);
      photonphi_histo->Fill(phPurityphi-phPurityphi_alt);
      photonPHIso_histo->Fill(phPurityPHiso-phPurityPHiso_alt);
      photonNHIso_histo->Fill(phPurityNHiso-phPurityNHiso_alt);
      photonCHIso_histo->Fill(phPurityCHiso-phPurityCHiso_alt);
      photonHoE_histo->Fill(phPurityhoe-phPurityhoe_alt);
      photonNum_histo->Fill(nphotons-nphotons_alt);
      photonEveto_histo->Fill(phPurityElectronVeto-phPurityElectronVeto_alt);
      photonSieie_histo->Fill(phPuritysieie-phPuritysieie_alt);
      //////////////
      photonpt_histo_2D->Fill(phPuritypt,phPuritypt_alt);
      photoneta_histo_2D->Fill(phPurityeta,phPurityeta_alt);
      photonphi_histo_2D->Fill(phPurityphi,phPurityphi_alt);
      photonPHIso_histo_2D->Fill(phPurityPHiso,phPurityPHiso_alt);
      photonNHIso_histo_2D->Fill(phPurityNHiso,phPurityNHiso_alt);
      photonCHIso_histo_2D->Fill(phPurityCHiso,phPurityCHiso_alt);
      photonHoE_histo_2D->Fill(phPurityhoe,phPurityhoe_alt);
      photonSieie_histo_2D->Fill(phPuritysieie,phPuritysieie_alt);
    }
  }
  cout<<endl;
  cout<<"Not matched events after selections "<<notMatchedEvents<<" event passing selections "<<effectiveEvents<<" i.e. "<<double(notMatchedEvents)/effectiveEvents*100 << " %"<<endl;

  if(doOppositeMatching){

    nEvents = 0;
    nPart   = 10000;
    nTotal  = chain_2->GetEntries();
    notMatchedEvents = 0;
    effectiveEvents  = 0;

    for(long int entry = 0; entry < chain_2->GetEntries(); entry++){

      int entry_status = chain_2->GetEntry(entry);
      if(entry_status <= 0) continue;
      
      cout.flush();
      if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
      nEvents++;

      // to select only ICHEP data
      if(useOnlyICHEP and run_alt > 276242) continue;  

      //////// -- event selection --> everything except for Photon ID selections, i.e. H/E, Sietaieta, Isolation ..etc
      // trigger
      int hlt = 0;
      hlt = hltp165_alt+hltp175_alt;
      if(not hlt) continue;     
      /// apply met filters
      if(not fhbhe_alt or not fhbiso_alt or not fcsct_alt or not feeb_alt or not feeb_alt or not fetp_alt or not fvtx_alt or not fcsc_alt) continue;      
      /// vetoes
      if(ntausraw_alt != 0) continue;
      if(nbjets_alt   != 0) continue;
      if(nmuons_alt   != 0) continue;
      if(nelectrons_alt != 0) continue;
      if(nphotons_alt == 0) continue;
      // photon candidate
      if(phPuritypt_alt < 175) continue;
      if(fabs(phPurityeta_alt) > 1.5) continue;

      if(applyPhotonID){
	if(phPurityhoe_alt   > 0.0396) continue;
	if(phPuritysieie_alt > 0.01022) continue;
	if(phPurityElectronVeto_alt == 0) continue;
	if(phPurityNHiso_alt > 2.725 + 0.0148*phPuritypt_alt + 0.000017*phPuritypt_alt*phPuritypt_alt) continue;
      }
      
      // jets
      if(njets_alt < 1) continue;
      // check overlap with leading jet --> cleaning is done wrt loose photons by default in the trees
      int ijet = 0;
      for(size_t i = 0 ; i < jetpt_alt->size(); i++){
	float deta = fabs(jeteta_alt->at(i)-phPurityeta_alt);
	float dphi = fabs(phPurityphi_alt-jetphi_alt->at(i));
	if(dphi > TMath::Pi())
	  dphi = 2*TMath::Pi()-dphi;
	if(sqrt(deta*deta+dphi*dphi) > 0.4){
	  ijet = i;
	  break;
	}
      }

      if(jetpt_alt->at(ijet) < 100) continue;
      if(fabs(jeteta_alt->at(ijet)) > 2.5) continue;
      // jet-met dphi
      int njet = 0;
      float mindphi = 99;
      for(size_t i= 0 ; i < jetpt_alt->size(); i++){
	// check pt
	if(jetpt_alt->at(i) < 30) continue;
	// check if overlaps with the photon candidate
	float deta = fabs(jeteta_alt->at(i)-phPurityeta_alt);
	float dphi = fabs(phPurityphi_alt-jetphi_alt->at(i));
	if(dphi > TMath::Pi())
	  dphi = 2*TMath::Pi()-dphi;
	if(sqrt(deta*deta+dphi*dphi) < 0.4) continue;
	njet++;
	if(njet > 4) continue;
	// calculate px and py
	float metx = t1met_alt*cos(t1metphi_alt)+phPuritypt_alt*cos(phPurityphi_alt);
	float mety = t1met_alt*sin(t1metphi_alt)+phPuritypt_alt*sin(phPurityphi_alt);      
	TLorentzVector met;
	met.SetPxPyPzE(metx,mety,0.,sqrt(metx*metx+mety*mety));
	float dphitemp = fabs(met.Phi()-jetphi_alt->at(ijet));
	if(dphitemp > TMath::Pi())
	  dphitemp = 2*TMath::Pi()-dphitemp;
	if(dphitemp < mindphi)
	  mindphi = dphitemp; 
      }    
      if(mindphi < 0.5) continue;
      
      // access to the other tree
      effectiveEvents++;
      entry_status = chain_1->GetEntryWithIndex(run_alt,event_alt);
      
      // events in chain_1 but not in 2
      if(entry_status <= 0){
	notMatchedEvents++;
	photonpt_histo_only2->Fill(phPuritypt_alt);
	photoneta_histo_only2->Fill(phPurityeta_alt);
	photonphi_histo_only2->Fill(phPurityphi_alt);
	photonPHIso_histo_only2->Fill(phPurityPHiso_alt);
	photonCHIso_histo_only2->Fill(phPurityCHiso_alt);
	photonNHIso_histo_only2->Fill(phPurityNHiso_alt);
	photonHoE_histo_only2->Fill(phPurityhoe_alt);
	photonNum_histo_only2->Fill(nphotons_alt);
	photonEveto_histo_only2->Fill(phPurityElectronVeto_alt);
	photonSieie_histo_only2->Fill(phPuritysieie_alt);
      }
    }
    
    cout<<endl;
    cout<<"Not matched events after selections "<<notMatchedEvents<<" event passing selections "<<effectiveEvents<<" i.e. "<<double(notMatchedEvents)/effectiveEvents*100<< " %"<<endl;    
  }
  
  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();

  drawPlot(canvas,photonpt_histo,"#Delta p_{T} [GeV]",outputDIR);  
  drawPlot(canvas,photoneta_histo,"#Delta #eta",outputDIR);  
  drawPlot(canvas,photonphi_histo,"#Delta #phi",outputDIR);  
  drawPlot(canvas,photonPHIso_histo,"#Delta Photon Iso [GeV]",outputDIR);  
  drawPlot(canvas,photonCHIso_histo,"#Delta Charged Iso [GeV]",outputDIR);  
  drawPlot(canvas,photonNHIso_histo,"#Delta Neutral Iso [GeV]",outputDIR);  
  drawPlot(canvas,photonHoE_histo,"#Delta H/E",outputDIR);    
  drawPlot(canvas,photonNum_histo,"#Delta N_{ph}",outputDIR);  
  drawPlot(canvas,photonEveto_histo,"#Delta E-veto",outputDIR);  
  drawPlot(canvas,photonSieie_histo,"#Delta #sigma_{i#eta-i#eta} [GeV]",outputDIR);  

  drawPlot2D(canvas,photonpt_histo_2D,"p_{T} [GeV]",outputDIR);  
  drawPlot2D(canvas,photoneta_histo_2D,"#eta",outputDIR);  
  drawPlot2D(canvas,photonphi_histo_2D,"#phi",outputDIR);  
  drawPlot2D(canvas,photonPHIso_histo_2D,"Photon Iso [GeV]",outputDIR);  
  drawPlot2D(canvas,photonCHIso_histo_2D,"Charged Iso [GeV]",outputDIR);  
  drawPlot2D(canvas,photonNHIso_histo_2D,"Neutral Iso [GeV]",outputDIR);  
  drawPlot2D(canvas,photonHoE_histo_2D,"H/E",outputDIR);  
  drawPlot2D(canvas,photonSieie_histo_2D,"#sigma_{i#eta-i#eta} [GeV]",outputDIR);  

  drawPlot(canvas,photonpt_histo_only1,"p_{T} [GeV]",outputDIR);  
  drawPlot(canvas,photoneta_histo_only1,"#eta",outputDIR);  
  drawPlot(canvas,photonphi_histo_only1,"#phi",outputDIR);  
  drawPlot(canvas,photonPHIso_histo_only1,"Photon Iso [GeV]",outputDIR);  
  drawPlot(canvas,photonCHIso_histo_only1,"Charged Iso [GeV]",outputDIR);  
  drawPlot(canvas,photonNHIso_histo_only1,"Neutral Iso [GeV]",outputDIR);  
  drawPlot(canvas,photonHoE_histo_only1,"H/E",outputDIR);  
  drawPlot(canvas,photonSieie_histo_only1,"#sigma_{i#eta-i#eta} [GeV]",outputDIR);  

  if(doOppositeMatching){
    drawPlot(canvas,photonpt_histo_only2,"p_{T} [GeV]",outputDIR);  
    drawPlot(canvas,photoneta_histo_only2,"#eta",outputDIR);  
    drawPlot(canvas,photonphi_histo_only2,"#phi",outputDIR);  
    drawPlot(canvas,photonPHIso_histo_only2,"Photon Iso [GeV]",outputDIR);  
    drawPlot(canvas,photonCHIso_histo_only2,"Charged Iso [GeV]",outputDIR);  
    drawPlot(canvas,photonNHIso_histo_only2,"Neutral Iso [GeV]",outputDIR);  
    drawPlot(canvas,photonHoE_histo_only2,"H/E",outputDIR);  
    drawPlot(canvas,photonSieie_histo_only2,"#sigma_{i#eta-i#eta} [GeV]",outputDIR);  
  }

}
