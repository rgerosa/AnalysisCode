#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

//////////////
//////////////
//////////////

TH1F* photonpt_histo  = new TH1F("photonpt","", 50,-100,100);
TH1F* photoneta_histo = new TH1F("photoneta","",50,-0.15,0.15);
TH1F* photonphi_histo = new TH1F("photonphi","",50,-0.06,0.06);
TH1F* photonPHIso_histo = new TH1F("photonPHIso","",80,-2.0,2.0);
TH1F* photonCHIso_histo = new TH1F("photonCHIso","",80,-5.0,5.0);
TH1F* photonNHIso_histo = new TH1F("photonNHIso","",80,-2.0,2.0);
TH1F* photonHoE_histo = new TH1F("photonHoE","",50,-0.1,0.1);
TH1F* photonNum_histo = new TH1F("nphotons","",3,-1.5,1.5);
TH1F* photonEveto_histo = new TH1F("photonEveto","",3,-1.5,1.5);
TH1F* photonSieie_histo = new TH1F("photonSieie","",50,-0.025,0.025);

TH2F* photonpt_histo_2D  = new TH2F("photonpt_2D","", 25,175,1000,25,175,1000);
TH2F* photoneta_histo_2D = new TH2F("photoneta_2D","",25,-1.5,1.5,25,-1.5,1.5);
TH2F* photonphi_histo_2D = new TH2F("photonphi_2D","",25,-3.14,3.14,25,-3.14,3.14);
TH2F* photonPHIso_histo_2D = new TH2F("photonPHIso_2D","",40,0,15,40,0,15);
TH2F* photonCHIso_histo_2D = new TH2F("photonCHIso_2D","",40,0,25,40,0,25);
TH2F* photonNHIso_histo_2D = new TH2F("photonNHIso_2D","",25,0,15,25,0,15);
TH2F* photonHoE_histo_2D = new TH2F("photonHoE_2D","",25,0,0.1,25,0,0.1);
TH2F* photonSieie_histo_2D = new TH2F("photonSieie_2D","",50,0,0.03,50,0,0.03);

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
  CMS_lumi(pad1,"12.9");  
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
  CMS_lumi(pad1,"12.9");  

  if(logScale)
    pad1->SetLogz();

  pad1->SaveAs((outputDIR+"/"+string(histo_1->GetName())+"_comparison.png").c_str(),"png");
  pad1->SaveAs((outputDIR+"/"+string(histo_1->GetName())+"_comparison.pdf").c_str(),"pdf");

  return ;
}

//////////////
void makeComparisonData(string inputDIR_1, string inputDIR_2, bool useOnlyICHEP, string outputDIR){


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
  UChar_t      hltp165 = 0, hltp175 = 0;
  UChar_t      fhbhe = 0, fhbiso = 0, fcsct = 0, feeb = 0, fetp = 0, fvtx = 0, fcsc = 0;
  unsigned int njets = 0, nmuons = 0, nelectrons = 0, nphotons = 0, ntausraw = 0, ntaus = 0, nincjets = 0, nbjets = 0;
  vector<float> *jeteta = 0, *jetpt = 0, *jetphi = 0,*jetm = 0,*chfrac = 0,*nhfrac = 0;
  float jpmdphi = 0;
  float phPuritypt = 0,phPurityeta = 0,phPurityphi = 0,phPurityElectronVeto = 0,phPurityPHiso = 0;
  float phPurityCHiso = 0,phPurityNHiso = 0,phPurityhoe = 0,phPuritysieie = 0,phPurityRND04PHiso = 0,phPurityRND08PHiso = 0;
  float phPurityEAEGamma = 0,rho = 0;

  chain_1->SetBranchStatus("*",kFALSE);
  chain_1->SetBranchStatus("run",kTRUE);
  chain_1->SetBranchStatus("lumi",kTRUE);
  chain_1->SetBranchStatus("event",kTRUE);
  chain_1->SetBranchStatus("hlt*",kTRUE);
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
  chain_1->SetBranchStatus("combinejetCHfrac",kTRUE);
  chain_1->SetBranchStatus("combinejetNHfrac",kTRUE);
  chain_1->SetBranchStatus("t1*",kTRUE);
  chain_1->SetBranchStatus("incjet*metdphimin4",kTRUE);
  chain_1->SetBranchStatus("ph*",kTRUE);
  chain_1->SetBranchStatus("rho*",kTRUE);
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
  chain_1->SetBranchAddress("combinejetCHfrac",&chfrac);
  chain_1->SetBranchAddress("combinejetNHfrac",&nhfrac);
  chain_1->SetBranchAddress("phPuritypt",&phPuritypt);
  chain_1->SetBranchAddress("phPurityeta",&phPurityeta);
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
  chain_1->SetBranchAddress("incjetphmetdphimin4",&jpmdphi);

  chain_1->BuildIndex("run","event");

  ///
  cout<<"Build second chain with index "<<endl;
  unsigned int run_alt = 0, lumi_alt = 0, event_alt = 0;
  UChar_t hltp165_alt = 0,hltp175_alt = 0;
  UChar_t fhbhe_alt = 0,fhbiso_alt = 0,fcsct_alt = 0,feeb_alt = 0,fetp_alt = 0,fvtx_alt = 0,fcsc_alt = 0;
  unsigned int njets_alt = 0,nmuons_alt = 0,nelectrons_alt = 0,nphotons_alt = 0,ntausraw_alt = 0,ntaus_alt = 0,nincjets_alt = 0,nbjets_alt = 0;
  vector<double> *jeteta_alt = 0, *jetpt_alt = 0, *jetphi_alt = 0,*jetm_alt = 0,*chfrac_alt = 0,*nhfrac_alt = 0;
  double jpmdphi_alt = 0;
  double phPuritypt_alt = 0,phPurityeta_alt = 0,phPurityphi_alt = 0,phPurityElectronVeto_alt = 0,phPurityPHiso_alt = 0;
  double phPurityCHiso_alt = 0,phPurityNHiso_alt = 0,phPurityhoe_alt = 0,phPuritysieie_alt = 0,phPurityRND04PHiso_alt = 0,phPurityRND08PHiso_alt = 0;
  double phPurityEAEGamma_alt = 0,rho_alt = 0;

  chain_2->SetBranchStatus("*",kFALSE);
  chain_2->SetBranchStatus("run",kTRUE);
  chain_2->SetBranchStatus("lumi",kTRUE);
  chain_2->SetBranchStatus("event",kTRUE);
  chain_2->SetBranchStatus("hlt*",kTRUE);
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
  chain_2->SetBranchStatus("combinejetCHfrac",kTRUE);
  chain_2->SetBranchStatus("combinejetNHfrac",kTRUE);
  chain_2->SetBranchStatus("t1*",kTRUE);
  chain_2->SetBranchStatus("incjet*metdphimin4",kTRUE);
  chain_2->SetBranchStatus("ph*",kTRUE);
  chain_2->SetBranchStatus("rho*",kTRUE);
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
  chain_2->SetBranchAddress("combinejetCHfrac",&chfrac_alt);
  chain_2->SetBranchAddress("combinejetNHfrac",&nhfrac_alt);
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
  chain_2->SetBranchAddress("incjetphmetdphimin4",&jpmdphi_alt);

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
    // jets
    if(njets < 1) continue;
    if(jetpt->at(0) < 100) continue;
    if(fabs(jeteta->at(0)) > 2.5) continue;
    // photon candidate
    if(phPuritypt < 175) continue;
    if(fabs(phPurityeta) > 1.5) continue;
    if(jpmdphi < 0.5) continue;
             
    // access to the other tree
    effectiveEvents++;
    entry_status = chain_2->GetEntryWithIndex(run,event);
    
    if(entry_status <= 0){
      notMatchedEvents++;
      continue;
    }

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
  cout<<endl;
  cout<<"Not matched events after selections "<<notMatchedEvents<<" event passing selections "<<effectiveEvents<<" i.e. "<<double(notMatchedEvents)/effectiveEvents*100<<endl;

  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();

  drawPlot(canvas,photonpt_histo,"#Delta p_{T} [GeV]",outputDIR);  
  drawPlot(canvas,photoneta_histo,"#Delta #eta",outputDIR);  
  drawPlot(canvas,photonphi_histo,"#Delta #phi",outputDIR);  
  drawPlot(canvas,photonPHIso_histo,"#Delta Photon Iso [GeV]",outputDIR);  
  drawPlot(canvas,photonCHIso_histo,"#Delta Charged Iso [GeV]",outputDIR);  
  drawPlot(canvas,photonNHIso_histo,"#Delta Neutral Iso [GeV]",outputDIR);  
  drawPlot(canvas,photonHoE_histo,"#Delta H/E",outputDIR);  
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
}
