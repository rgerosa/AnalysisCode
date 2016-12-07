#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static float recoilSelection = 400;

TH1F* bosonpt_histo  = new TH1F("bosonpt","",50,-100,100);
TH1F* bosoneta_histo = new TH1F("bosoneta","",50,-0.15,0.15);
TH1F* bosonphi_histo = new TH1F("bosonphi","",50,-0.1,0.1);
TH1F* recoil_histo   = new TH1F("recoil","",50,-100,100);
TH1F* recoilphi_histo = new TH1F("recoilphi","",50,-0.25,0.25);
TH1F* pfmet_histo     = new TH1F("pfmet","",50,-100,100);
TH1F* pfmetphi_histo  = new TH1F("pfmetphi","",50,-0.25,0.25);
TH1F* jetpt_histo  = new TH1F("jetpt","",50,-100,100);


TH2F* bosonpt_histo_2D  = new TH2F("bosonpt_2d","",50,200,1100,50,200,1100);
TH2F* bosoneta_histo_2D = new TH2F("bosoneta_2d","",50,-2.5,2.5,50,-2.5,2.5);
TH2F* bosonphi_histo_2D = new TH2F("bosonphi_2d","",50,-3.14,3.14,50,-3.14,3.14);
TH2F* recoil_histo_2D   = new TH2F("recoil_2d","",50,300,1200,50,300,1200);
TH2F* recoilphi_histo_2D = new TH2F("recoilphi_2d","",50,-3.14,3.14,50,-3.14,3.14);
TH2F* pfmet_histo_2D     = new TH2F("pfmet_2d","",50,0,300,50,0,300);
TH2F* pfmetphi_histo_2D  = new TH2F("pfmetphi_2d","",50,-3.14,3.14,50,-3.14,3.14);
TH2F* jetpt_histo_2D     = new TH2F("jetpt_2d","",50,100,900,50,100,900);


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



void makeComparisonData(string inputDIR_1, string inputDIR_2, Sample sample, Category category, string outputDIR){

  bosonpt_histo->Sumw2();
  bosoneta_histo->Sumw2();
  bosonphi_histo->Sumw2();
  recoil_histo->Sumw2();
  recoilphi_histo->Sumw2();
  pfmet_histo->Sumw2();
  pfmetphi_histo->Sumw2();
  jetpt_histo->Sumw2();

  bosonpt_histo_2D->Sumw2();
  bosoneta_histo_2D->Sumw2();
  bosonphi_histo_2D->Sumw2();
  recoil_histo_2D->Sumw2();
  recoilphi_histo_2D->Sumw2();
  pfmet_histo_2D->Sumw2();
  pfmetphi_histo_2D->Sumw2();
  jetpt_histo_2D->Sumw2();

  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  gROOT->SetBatch();

  TChain* chain_1 = new TChain("tree/tree");
  TChain* chain_2 = new TChain("tree/tree");

  if(sample == Sample::sig){
    chain_1->Add((inputDIR_1+"/sigfilter/*root").c_str());
    chain_2->Add((inputDIR_2+"/sigfilter/*root").c_str());
  }
  else if(sample == Sample::wmn){
    chain_1->Add((inputDIR_1+"/wmnfilter/*root").c_str());
    chain_2->Add((inputDIR_2+"/wmnfilter/*root").c_str());
  }
  else if(sample == Sample::wen){
    chain_1->Add((inputDIR_1+"/wenfilter/*root").c_str());
    chain_2->Add((inputDIR_2+"/wenfilter/*root").c_str());
  }
  else if(sample == Sample::zmm){
    chain_1->Add((inputDIR_1+"/zmmfilter/*root").c_str());
    chain_2->Add((inputDIR_2+"/zmmfilter/*root").c_str());
  }
  else if(sample == Sample::zee){
    chain_1->Add((inputDIR_1+"/zeefilter/*root").c_str());
    chain_2->Add((inputDIR_2+"/zeefilter/*root").c_str());
  }
  else if(sample == Sample::gam){
    chain_1->Add((inputDIR_1+"/gamfilter/*root").c_str());
    chain_2->Add((inputDIR_2+"/gamfilter/*root").c_str());
  }

  ///////////
  cout<<"Build first chain with index "<<endl;
  unsigned int run=0;
  unsigned int lumi=0;
  unsigned int event=0;
  UChar_t hltm90=0,hltm100=0,hltm110=0,hltm120=0,hltmwm120=0,hltmwm170=0,hltmwm300=0,hltmwm90=0,hlte=0,hltenoiso=0,hltp165=0,hltp175=0;
  UChar_t fhbhe=0,fhbiso=0,fcsct=0,feeb=0,fetp=0,fvtx=0,fbadmu=0,fbadch=0,fcsc=0;
  unsigned int njets=0,nmuons=0,nelectrons=0,nphotons=0,ntausraw=0,nincjets=0,nbjets=0;
  vector<float> *jeteta=0, *jetpt=0, *jetphi=0,*jetm=0,*chfrac=0,*nhfrac=0;
  float incjetmumetdphimin4=0,incjetelmetdphimin4=0,incjetphmetdphimin4=0;
  float met=0,metphi=0,mmet=0,mmetphi=0,emet=0,emetphi=0,pmet=0,pmetphi=0,pfmet=0,pfmetphi=0,metcalo=0;
  float jmmdphi=0,jemdphi=0,jpmdphi=0;
  int   mu1pid=0,mu2pid=0,mu1id=0,mu2id=0;
  float mu1pt=0,mu2pt=0,mu1eta=0,mu2eta=0,mu1phi=0,mu2phi=0;
  int   el1pid=0,el2pid=0,el1id=0,el2id=0;
  float el1pt=0,el2pt=0,el1eta=0,el2eta=0,el1phi=0,el2phi=0;
  int   phidm=0;
  float phpt=0,pheta=0,phphi=0;
  float wmt=0,wemt=0,zmass=0,zeemass=0,zpt=0,zeept=0,zeeeta=0,zeta=0,zphi=0,zeephi=0;
  vector<float> *boostedJetpt=0,*boostedJeteta=0,*boostedJetphi=0,*prunedJetm=0,*boostedJettau2=0,*boostedJettau1=0;

  chain_1->SetBranchStatus("*",kFALSE);
  chain_1->SetBranchStatus("run",kTRUE);
  chain_1->SetBranchStatus("lumi",kTRUE);
  chain_1->SetBranchStatus("event",kTRUE);
  chain_1->SetBranchStatus("hlt*",kTRUE);
  chain_1->SetBranchStatus("flag*",kTRUE);
  chain_1->SetBranchStatus("njets",kTRUE);
  chain_1->SetBranchStatus("nmuons",kTRUE);
  chain_1->SetBranchStatus("nelectrons",kTRUE);
  chain_1->SetBranchStatus("nphotons",kTRUE);
  chain_1->SetBranchStatus("ntausrawold",kTRUE);
  chain_1->SetBranchStatus("njetsinc",kTRUE);
  chain_1->SetBranchStatus("nbjetslowpt",kTRUE);
  chain_1->SetBranchStatus("combinejeteta",kTRUE);
  chain_1->SetBranchStatus("combinejetpt", kTRUE);
  chain_1->SetBranchStatus("combinejetphi",kTRUE);
  chain_1->SetBranchStatus("combinejetm",kTRUE);
  chain_1->SetBranchStatus("combinejetCHfrac",kTRUE);
  chain_1->SetBranchStatus("combinejetNHfrac",kTRUE);
  chain_1->SetBranchStatus("t1*",kTRUE);
  chain_1->SetBranchStatus("pfmet",kTRUE);
  chain_1->SetBranchStatus("incjet*metdphimin4",kTRUE);
  chain_1->SetBranchStatus("mu*",kTRUE);
  chain_1->SetBranchStatus("el*",kTRUE);
  chain_1->SetBranchStatus("ph*",kTRUE);
  chain_1->SetBranchStatus("wmt",kTRUE);
  chain_1->SetBranchStatus("wemt",kTRUE);
  chain_1->SetBranchStatus("zmass",kTRUE);
  chain_1->SetBranchStatus("zeemass",kTRUE);
  chain_1->SetBranchStatus("zpt",kTRUE);
  chain_1->SetBranchStatus("zeept",kTRUE);
  chain_1->SetBranchStatus("boostedJet*",kTRUE);
  chain_1->SetBranchStatus("prunedJetm",kTRUE);

  chain_1->SetBranchAddress("run",&run);
  chain_1->SetBranchAddress("lumi",&lumi);
  chain_1->SetBranchAddress("event",&event);
  chain_1->SetBranchAddress("hltmet90",&hltm90);
  chain_1->SetBranchAddress("hltmet100",&hltm100);
  chain_1->SetBranchAddress("hltmet110",&hltm110);
  chain_1->SetBranchAddress("hltmet120",&hltm120);
  chain_1->SetBranchAddress("hltmetwithmu170",&hltmwm170);
  chain_1->SetBranchAddress("hltmetwithmu300",&hltmwm300);
  chain_1->SetBranchAddress("hltmetwithmu90",&hltmwm90);
  chain_1->SetBranchAddress("hltsingleel",&hlte);
  chain_1->SetBranchAddress("hltelnoiso",&hltenoiso);
  chain_1->SetBranchAddress("hltphoton165",&hltp165);
  chain_1->SetBranchAddress("hltphoton175",&hltp175);
  chain_1->SetBranchAddress("flaghbhenoise",&fhbhe);
  chain_1->SetBranchAddress("flaghbheiso",&fhbiso);
  chain_1->SetBranchAddress("flagcsctight",&fcsct);
  chain_1->SetBranchAddress("flageebadsc",&feeb);
  chain_1->SetBranchAddress("flagecaltp",&fetp);
  chain_1->SetBranchAddress("flaggoodvertices",&fvtx);
  chain_1->SetBranchAddress("flagbadpfmu",&fbadmu);
  chain_1->SetBranchAddress("flagbadchpf",&fbadch);
  chain_1->SetBranchAddress("flagglobaltighthalo",&fcsc);
  chain_1->SetBranchAddress("njets",&njets);
  chain_1->SetBranchAddress("nmuons",&nmuons);
  chain_1->SetBranchAddress("nelectrons",&nelectrons);
  chain_1->SetBranchAddress("nphotons",&nphotons);
  chain_1->SetBranchAddress("nelectrons",&nelectrons);
  chain_1->SetBranchAddress("ntausrawold",&ntausraw);
  chain_1->SetBranchAddress("njetsinc",&nincjets);
  chain_1->SetBranchAddress("nbjetslowpt",&nbjets);
  chain_1->SetBranchAddress("combinejeteta",&jeteta);
  chain_1->SetBranchAddress("combinejetpt", &jetpt);
  chain_1->SetBranchAddress("combinejetphi",&jetphi);
  chain_1->SetBranchAddress("combinejetm",&jetm);
  chain_1->SetBranchAddress("combinejetCHfrac",&chfrac);
  chain_1->SetBranchAddress("combinejetNHfrac",&nhfrac);
  chain_1->SetBranchAddress("t1pfmet",&met);
  chain_1->SetBranchAddress("t1pfmetphi",&metphi);
  chain_1->SetBranchAddress("t1mumet",&mmet);
  chain_1->SetBranchAddress("t1mumetphi",&mmetphi);
  chain_1->SetBranchAddress("t1elmet",&emet);
  chain_1->SetBranchAddress("t1elmetphi",&emetphi);
  chain_1->SetBranchAddress("t1phmet",&pmet);
  chain_1->SetBranchAddress("t1phmetphi",&pmetphi);
  chain_1->SetBranchAddress("pfmet",&pfmet);
  chain_1->SetBranchAddress("pfmetphi",&pfmetphi);
  chain_1->SetBranchAddress("calomet",&metcalo);
  chain_1->SetBranchAddress("incjetmumetdphimin4",&incjetmumetdphimin4);
  chain_1->SetBranchAddress("incjetelmetdphimin4",&incjetelmetdphimin4);
  chain_1->SetBranchAddress("incjetphmetdphimin4",&incjetphmetdphimin4);
  chain_1->SetBranchAddress("mu1pid",&mu1pid);
  chain_1->SetBranchAddress("mu2pid",&mu2pid);
  chain_1->SetBranchAddress("mu1id",&mu1id);
  chain_1->SetBranchAddress("mu2id",&mu2id);
  chain_1->SetBranchAddress("mu1pt",&mu1pt);
  chain_1->SetBranchAddress("mu2pt",&mu2pt);
  chain_1->SetBranchAddress("mu1eta",&mu1eta);
  chain_1->SetBranchAddress("mu2eta",&mu2eta);
  chain_1->SetBranchAddress("mu1phi",&mu1phi);
  chain_1->SetBranchAddress("mu2phi",&mu2phi);
  chain_1->SetBranchAddress("el1pid",&el1pid);
  chain_1->SetBranchAddress("el2pid",&el2pid);
  chain_1->SetBranchAddress("el1id",&el1id);
  chain_1->SetBranchAddress("el2id",&el2id);
  chain_1->SetBranchAddress("el1pt",&el1pt);
  chain_1->SetBranchAddress("el2pt",&el2pt);
  chain_1->SetBranchAddress("el1eta",&el1eta);
  chain_1->SetBranchAddress("el2eta",&el2eta);
  chain_1->SetBranchAddress("el1phi",&el1phi);
  chain_1->SetBranchAddress("el2phi",&el2phi);
  chain_1->SetBranchAddress("phidm",&phidm);
  chain_1->SetBranchAddress("phpt",&phpt);
  chain_1->SetBranchAddress("pheta",&pheta);
  chain_1->SetBranchAddress("phphi",&phphi);
  chain_1->SetBranchAddress("wmt",&wmt);
  chain_1->SetBranchAddress("wemt",&wemt);
  chain_1->SetBranchAddress("zmass",&zmass);
  chain_1->SetBranchAddress("zeemass",&zeemass);
  chain_1->SetBranchAddress("zpt",&zpt);
  chain_1->SetBranchAddress("zeept",&zeept);
  chain_1->SetBranchAddress("zeta",&zeta);
  chain_1->SetBranchAddress("zeeeta",&zeeeta);
  chain_1->SetBranchAddress("zphi",&zphi);
  chain_1->SetBranchAddress("zeephi",&zeephi);
  chain_1->SetBranchAddress("boostedJetpt",&boostedJetpt);
  chain_1->SetBranchAddress("boostedJeteta",&boostedJeteta);
  chain_1->SetBranchAddress("boostedJetphi",&boostedJetphi);
  chain_1->SetBranchAddress("boostedJettau2",&boostedJettau2);
  chain_1->SetBranchAddress("boostedJettau1",&boostedJettau1);
  chain_1->SetBranchAddress("prunedJetm",&prunedJetm);

  chain_1->BuildIndex("run","event");

  ///
  cout<<"Build second chain with index "<<endl;
  unsigned int run_alt=0;
  unsigned int lumi_alt=0;
  unsigned int event_alt=0;
  double met_alt=0,metphi_alt=0,mmet_alt=0,mmetphi_alt=0,emet_alt=0,emetphi_alt=0,pmet_alt=0,pmetphi_alt=0,pfmet_alt=0,pfmetphi_alt=0,metcalo_alt=0;
  double mu1pt_alt=0,mu2pt_alt=0,mu1eta_alt=0,mu2eta_alt=0,mu1phi_alt=0,mu2phi_alt=0;
  double el1pt_alt=0,el2pt_alt=0,el1eta_alt=0,el2eta_alt=0,el1phi_alt=0,el2phi_alt=0;
  double phpt_alt=0,pheta_alt=0,phphi_alt=0;
  double zpt_alt=0, zeept_alt=0, zeta_alt=0, zeeeta_alt=0, zphi_alt=0, zeephi_alt=0;
  vector<float> *jetpt_alt=0;

  chain_2->SetBranchStatus("*",kFALSE);
  chain_2->SetBranchStatus("run",kTRUE);
  chain_2->SetBranchStatus("lumi",kTRUE);
  chain_2->SetBranchStatus("event",kTRUE);
  chain_2->SetBranchStatus("t1*",kTRUE);
  chain_2->SetBranchStatus("pfmet",kTRUE);
  chain_2->SetBranchStatus("mu*",kTRUE);
  chain_2->SetBranchStatus("el*",kTRUE);
  chain_2->SetBranchStatus("ph*",kTRUE);
  chain_2->SetBranchStatus("z*",kTRUE);
  chain_2->SetBranchStatus("combinejetpt",kTRUE);


  chain_2->SetBranchAddress("run",&run_alt);
  chain_2->SetBranchAddress("lumi",&lumi_alt);
  chain_2->SetBranchAddress("event",&event_alt);
  chain_2->SetBranchAddress("t1pfmet",&met_alt);
  chain_2->SetBranchAddress("t1pfmetphi",&metphi_alt);
  chain_2->SetBranchAddress("t1mumet",&mmet_alt);
  chain_2->SetBranchAddress("t1mumetphi",&mmetphi_alt);
  chain_2->SetBranchAddress("t1elmet",&emet_alt);
  chain_2->SetBranchAddress("t1elmetphi",&emetphi_alt);
  chain_2->SetBranchAddress("t1phmet",&pmet_alt);
  chain_2->SetBranchAddress("t1phmetphi",&pmetphi_alt);
  chain_2->SetBranchAddress("pfmet",&pfmet_alt);
  chain_2->SetBranchAddress("pfmetphi",&pfmetphi_alt);
  chain_2->SetBranchAddress("calomet",&metcalo_alt);
  chain_2->SetBranchAddress("mu1pt",&mu1pt_alt);
  chain_2->SetBranchAddress("mu2pt",&mu2pt_alt);
  chain_2->SetBranchAddress("mu1eta",&mu1eta_alt);
  chain_2->SetBranchAddress("mu2eta",&mu2eta_alt);
  chain_2->SetBranchAddress("mu1phi",&mu1phi_alt);
  chain_2->SetBranchAddress("mu2phi",&mu2phi_alt);
  chain_2->SetBranchAddress("el1pt",&el1pt_alt);
  chain_2->SetBranchAddress("el2pt",&el2pt_alt);
  chain_2->SetBranchAddress("el1eta",&el1eta_alt);
  chain_2->SetBranchAddress("el2eta",&el2eta_alt);
  chain_2->SetBranchAddress("el1phi",&el1phi_alt);
  chain_2->SetBranchAddress("el2phi",&el2phi_alt);
  chain_2->SetBranchAddress("phpt",&phpt_alt);
  chain_2->SetBranchAddress("pheta",&pheta_alt);
  chain_2->SetBranchAddress("phphi",&phphi_alt);
  chain_2->SetBranchAddress("zpt",&zpt_alt);
  chain_2->SetBranchAddress("zeept",&zeept_alt);
  chain_2->SetBranchAddress("zeta",&zeta_alt);
  chain_2->SetBranchAddress("zeeeta",&zeeeta_alt);
  chain_2->SetBranchAddress("zphi",&zphi_alt);
  chain_2->SetBranchAddress("zeephi",&zeephi_alt);
  chain_2->SetBranchAddress("combinejetpt",&jetpt_alt);
 

  chain_2->BuildIndex("run","event");

  cout<<"Number of events in chain 1 "<<chain_1->GetEntries()<<endl;
  cout<<"Number of events in chain 2 "<<chain_2->GetEntries()<<endl;
  long int nEvents = 0;
  int nPart        = 10000;
  long int nTotal  = chain_1->GetEntries();
  for(long int entry = 0; entry < chain_1->GetEntries(); entry++){

    int entry_status = chain_1->GetEntry(entry);
    if(entry_status <= 0) continue;

    //if(nEvents > nTotal/10) break;

    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;
    ////////
    int hlt = 0;
    if(sample == Sample::sig or sample == Sample::wmn or sample == Sample::zmm)
      hlt = hltm90+hltm100+hltm110+hltm120+hltmwm120+hltmwm170+hltmwm300+hltmwm90;
    else if(sample == Sample::wen or sample == Sample::zee)
      hlt = hlte+hltenoiso;
    else if(sample == Sample::gam)
      hlt = hltp165+hltp175;
    
    if(not hlt) continue;
     
    /// apply met filters
    if(not fhbhe or not fhbiso or not fcsct or not feeb or not feeb or not fetp or not fvtx or not fbadmu or not fbadch or not fcsc) continue;

    /// jets
    if(njets < 1 and (category == Category::monojet or category == Category::monoV)) continue;
    else if(category == Category::VBF and nincjets < 2) continue;
    // vetos
    if(ntausraw != 0 ) continue;
    if(nbjets != 0) continue;
    
    // tag objets
    if(sample == Sample::wen and (el1pt < 40 or el1id != 1 or wemt > 160 or pfmet < 50 or nmuons != 0 or nphotons != 0)) continue;
    else if(sample == Sample::wmn and (mu1pt < 20 or mu1id != 1 or wmt > 160 or nelectrons != 0 or nphotons != 0)) continue;
    else if(sample == Sample::gam and (phpt < 120 or fabs(pheta) > 1.442 or phidm != 1 or nmuons != 0 or nelectrons != 0)) continue;
    else if(sample == Sample::zmm and (not ((mu1pt > 20 and mu1id == 1 and fabs(mu1eta) < 2.4) or (mu2pt > 20 and mu2id == 1 and fabs(mu1eta) < 2.4)) or 
				       zmass < 60 or zmass > 120 or mu1id == mu2id or nelectrons != 0 or nphotons != 0)) continue;
    else if(sample == Sample::zee and (not ((el1pt > 40 and el1id == 1 and fabs(el1eta) < 2.5) or (el2pt > 40 and el2id == 1 and fabs(el1eta) < 2.45)) or 
				       zeemass < 60 or zeemass > 120 or el1id == el2id or nmuons != 0 or nphotons != 0)) continue;
    // min-dphi
    if((sample == Sample::wmn or sample == Sample::zmm or sample == Sample::sig) and incjetmumetdphimin4 < 0.5) continue;
    else if((sample == Sample::wen or sample == Sample::zee) and incjetelmetdphimin4 < 0.5) continue;
    else if(sample == Sample::gam and incjetphmetdphimin4 < 0.5) continue;

    // recoil
    if((sample == Sample::wmn or sample == Sample::zmm or sample == Sample::sig) and mmet < recoilSelection) continue;
    else if((sample == Sample::wen or sample == Sample::zee) and emet < recoilSelection) continue;
    else if(sample == Sample::gam and pmet < recoilSelection) continue;

    //
    Double_t metden = 0.0;
    if (sample == Sample::sig || sample == Sample::qcd) {metden = mmet;}
    else if (sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu){ metden = mmet;}
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel){ metden = emet;}
    else if (sample == Sample::qcdgam || sample == Sample::gam)  { metden = pmet;}    
    if(fabs(met-metcalo)/metden > 0.5) continue;

    if(category == Category::monojet){

      if(jetpt->at(0) < 100) continue;
      if(fabs(jeteta->at(0)) > 2.5) continue;
      if(chfrac->at(0) < 0.1) continue;
      if(nhfrac->at(0) > 0.8) continue;
      bool goodMonoV = false;

      if(boostedJetpt->size() != 0 and boostedJetpt->at(0) > 250 and fabs(boostedJeteta->at(0)) < 2.4 and boostedJettau2->at(0)/boostedJettau1->at(0) < 0.6 and 
	 prunedJetm->at(0) > 65 and prunedJetm->at(0) < 105 and metden > 250)
	goodMonoV = true;
      if(goodMonoV) continue;     

     
      // access to the other tree
      entry_status = chain_2->GetEntryWithIndex(run,event);

      if(entry_status <= 0) continue;

      if(run != run_alt or lumi != lumi_alt or event != event_alt){
	cerr<<"Problem in associating: [run,lumi,event] = ["<<run<<","<<lumi<<","<<event<<"] with ["<<run_alt<<","<<lumi_alt<<","<<event_alt<<"]"<<endl;
	continue;
      }

      if(sample == Sample::zmm){
	bosonpt_histo->Fill(zpt-zpt_alt);
	bosoneta_histo->Fill(zeta-zeta_alt);
	bosonphi_histo->Fill(zphi-zphi_alt);
	recoil_histo->Fill(mmet-mmet_alt);
	recoilphi_histo->Fill(mmetphi-mmetphi_alt);
	pfmet_histo->Fill(pfmet-pfmet_alt);
	pfmetphi_histo->Fill(pfmetphi-pfmetphi_alt);
      }
      else if(sample == Sample::zee){
	bosonpt_histo->Fill(zeept-zeept_alt);
	bosoneta_histo->Fill(zeeeta-zeeeta_alt);
	bosonphi_histo->Fill(zeephi-zeephi_alt);
	recoil_histo->Fill(emet-emet_alt);
	recoilphi_histo->Fill(emetphi-emetphi_alt);
	pfmet_histo->Fill(pfmet-pfmet_alt);
	pfmetphi_histo->Fill(pfmetphi-pfmetphi_alt);
      }
      else if(sample == Sample::wen){
	bosonpt_histo->Fill(el1pt-el1pt_alt);
	bosoneta_histo->Fill(el1eta-el1eta_alt);
	bosonphi_histo->Fill(el1phi-el1phi_alt);
	recoil_histo->Fill(emet-emet_alt);
	recoilphi_histo->Fill(emetphi-emetphi_alt);
	pfmet_histo->Fill(pfmet-pfmet_alt);
	pfmetphi_histo->Fill(pfmetphi-pfmetphi_alt);
      }
      else if(sample == Sample::wmn){
	bosonpt_histo->Fill(mu1pt-mu1pt_alt);
	bosoneta_histo->Fill(mu1eta-mu1eta_alt);
	bosonphi_histo->Fill(mu1phi-mu1phi_alt);
	recoil_histo->Fill(mmet-mmet_alt);
	recoilphi_histo->Fill(mmetphi-mmetphi_alt);
	pfmet_histo->Fill(pfmet-pfmet_alt);
	pfmetphi_histo->Fill(pfmetphi-pfmetphi_alt);
      }
      else if(sample == Sample::gam){	
	bosonpt_histo->Fill(phpt-phpt_alt);
	bosoneta_histo->Fill(pheta-pheta_alt);
	bosonphi_histo->Fill(phphi-phphi_alt);
	recoil_histo->Fill(pmet-pmet_alt);
	recoilphi_histo->Fill(pmetphi-pmetphi_alt);
	pfmet_histo->Fill(pfmet-pfmet_alt);
	pfmetphi_histo->Fill(pfmetphi-pfmetphi_alt);
	jetpt_histo->Fill(jetpt->at(0)-jetpt_alt->at(0));

	bosonpt_histo_2D->Fill(phpt,phpt_alt);
	bosoneta_histo_2D->Fill(pheta,pheta_alt);
	bosonphi_histo_2D->Fill(phphi,phphi_alt);
	recoil_histo_2D->Fill(pmet,pmet_alt);
	recoilphi_histo_2D->Fill(pmetphi,pmetphi_alt);
	pfmet_histo_2D->Fill(pfmet,pfmet_alt);
	pfmetphi_histo_2D->Fill(pfmetphi,pfmetphi_alt);
	jetpt_histo_2D->Fill(jetpt->at(0),jetpt_alt->at(0));

      }
      else if(sample == Sample::sig){
	recoil_histo->Fill(met-met_alt);
	recoilphi_histo->Fill(metphi-metphi_alt);
	pfmet_histo->Fill(pfmet-pfmet_alt);
	pfmetphi_histo->Fill(pfmetphi-pfmetphi_alt);
      }
    }
    //    else if(category == Category::monoV){
    //    }
    //    else if(category == Category::VBF){
    //    }
  }


  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();

  drawPlot(canvas,bosonpt_histo,"boson p_{T}",outputDIR);  
  drawPlot(canvas,bosoneta_histo,"boson #eta",outputDIR);  
  drawPlot(canvas,bosonphi_histo,"boson #phi",outputDIR);  
  drawPlot(canvas,recoil_histo,"Recoil",outputDIR);  
  drawPlot(canvas,recoilphi_histo,"Recoil #phi",outputDIR);  
  drawPlot(canvas,pfmet_histo,"Met",outputDIR);  
  drawPlot(canvas,pfmetphi_histo,"Met #phi",outputDIR);  

  drawPlot2D(canvas,bosonpt_histo_2D,"boson p_{T}",outputDIR);  
  drawPlot2D(canvas,bosoneta_histo_2D,"boson #eta",outputDIR);  
  drawPlot2D(canvas,bosonphi_histo_2D,"boson #phi",outputDIR);  
  drawPlot2D(canvas,recoil_histo_2D,"Recoil",outputDIR);  
  drawPlot2D(canvas,recoilphi_histo_2D,"Recoil #phi",outputDIR);  
  drawPlot2D(canvas,pfmet_histo_2D,"Met",outputDIR);  
  drawPlot2D(canvas,pfmetphi_histo_2D,"Met #phi",outputDIR);  

  

}
