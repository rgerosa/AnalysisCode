#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static float recoilSelection = 600;

TH1F* bosonpt_histo   = new TH1F("bosonpt","",50,-100,100);
TH1F* bosoneta_histo  = new TH1F("bosoneta","",50,-0.07,0.07);
TH1F* bosonphi_histo  = new TH1F("bosonphi","",50,-0.01,0.01);
TH1F* recoil_histo    = new TH1F("recoil","",50,-100,100);
TH1F* recoilphi_histo = new TH1F("recoilphi","",50,-0.25,0.25);
TH1F* pfmet_histo     = new TH1F("pfmet","",50,-100,100);
TH1F* pfmetphi_histo  = new TH1F("pfmetphi","",50,-0.25,0.25);
TH1F* jetpt_histo     = new TH1F("jetpt","",50,-30,30);
TH1F* jetpt2_histo    = new TH1F("jetpt2","",50,-30,30);
TH1F* njet_histo      = new TH1F("njet","",10,-5,5);
TH1F* njetinc_histo   = new TH1F("njetinc","",10,-5,5);
TH1F* jeteta_histo    = new TH1F("jeteta","",50,-0.1,0.1);
TH1F* jeteta2_histo   = new TH1F("jeteta2","",50,-0.1,0.1);
TH1F* jetmetdphi_histo   = new TH1F("jetmetdphi","",50,-0.25,0.25);

TH2F* bosonpt_histo_2D  = new TH2F("bosonpt_2d","",50,200,1100,50,200,1100);
TH2F* bosoneta_histo_2D = new TH2F("bosoneta_2d","",50,-2.5,2.5,50,-2.5,2.5);
TH2F* bosonphi_histo_2D = new TH2F("bosonphi_2d","",50,-3.14,3.14,50,-3.14,3.14);
TH2F* recoil_histo_2D   = new TH2F("recoil_2d","",50,300,1200,50,300,1200);
TH2F* recoilphi_histo_2D = new TH2F("recoilphi_2d","",50,-3.14,3.14,50,-3.14,3.14);
TH2F* pfmet_histo_2D     = new TH2F("pfmet_2d","",50,0,300,50,0,300);
TH2F* pfmetphi_histo_2D  = new TH2F("pfmetphi_2d","",50,-3.14,3.14,50,-3.14,3.14);
TH2F* jetpt_histo_2D     = new TH2F("jetpt_2d","",50,100,900,50,100,900);
TH2F* jetpt2_histo_2D    = new TH2F("jetpt2_2d","",50,100,900,50,100,900);
TH2F* jeteta_histo_2D    = new TH2F("jeteta_2d","",30,-2.5,2.5,30,-2.5,2.5);
TH2F* jeteta2_histo_2D   = new TH2F("jeteta2_2d","",30,-2.5,2.5,30,-2.5,2.5);
TH2F* njet_histo_2D      = new TH2F("njet_2d","",5,0,5,5,0,5);
TH2F* njetinc_histo_2D   = new TH2F("njetinc_2d","",5,0,5,5,0,5);
TH2F* jetmetdphi_histo_2D = new TH2F("jetmetdphi_2d","",50,0,3.14,50,0,3.14);

  
TH1F* bosonpt_histo_only1   = new TH1F("bosonpt_only1","",40,120,1000);
TH1F* bosoneta_histo_only1  = new TH1F("bosoneta_only1","",30,-1.5,1.5);
TH1F* bosonphi_histo_only1  = new TH1F("bosonphi_only1","",30,-3.14,3.14);
TH1F* recoil_histo_only1    = new TH1F("recoil_only1","",40,200,1200);
TH1F* recoilphi_histo_only1 = new TH1F("recoilphi_only1","",30,-3.14,3.14);
TH1F* pfmet_histo_only1     = new TH1F("pfmet_only1","",40,0,400);
TH1F* pfmetphi_histo_only1  = new TH1F("pfmetphi_only1","",30,-3.14,3.14);
TH1F* jetpt_histo_only1     = new TH1F("jetpt_only1","",40,100,1000);
TH1F* jetpt2_histo_only1    = new TH1F("jetpt2_only1","",40,30,500);
TH1F* njet_histo_only1      = new TH1F("njet_only1","",7,1,8);
TH1F* njetinc_histo_only1      = new TH1F("njetinc_only1","",7,1,8);
TH1F* jeteta_histo_only1    = new TH1F("jeteta_only1","",30,-2.5,2.5);
TH1F* jeteta2_histo_only1   = new TH1F("jeteta2_only1","",50,-2.5,2.5);
TH1F* jetmetdphi_histo_only1   = new TH1F("jetmetdphi_only1","",50,0.5,3.14);

TH1F* bosonpt_histo_only2   = new TH1F("bosonpt_only2","",40,120,1000);
TH1F* bosoneta_histo_only2  = new TH1F("bosoneta_only2","",30,-1.5,1.5);
TH1F* bosonphi_histo_only2  = new TH1F("bosonphi_only2","",30,-3.14,3.14);
TH1F* recoil_histo_only2    = new TH1F("recoil_only2","",40,200,1200);
TH1F* recoilphi_histo_only2 = new TH1F("recoilphi_only2","",30,-3.14,3.14);
TH1F* pfmet_histo_only2     = new TH1F("pfmet_only2","",40,0,400);
TH1F* pfmetphi_histo_only2  = new TH1F("pfmetphi_only2","",30,-3.14,3.14);
TH1F* jetpt_histo_only2     = new TH1F("jetpt_only2","",40,100,1000);
TH1F* jetpt2_histo_only2    = new TH1F("jetpt2_only2","",40,30,500);
TH1F* njet_histo_only2      = new TH1F("njet_only2","",7,1,7);
TH1F* njetinc_histo_only2   = new TH1F("njetinc_only2","",7,1,7);
TH1F* jeteta_histo_only2    = new TH1F("jeteta_only2","",30,-2.5,2.5);
TH1F* jeteta2_histo_only2   = new TH1F("jeteta2_only2","",50,-2.5,2.5);
TH1F* jetmetdphi_histo_only2   = new TH1F("jetmetdphi_only2","",50,0.5,3.14);
  


class eventID{

public:
  eventID(unsigned int & eventid, unsigned int & runid){
    eventid_ = eventid;
    runid_ = runid;
  }

 bool operator== (const eventID & a) const{
    if(eventid_ == a.eventid_ and runid_ == a.runid_) return true;
    else return false;
  }

 bool operator< (const eventID & a) const{
    if(runid_ < a.runid_) return true;
    else if(runid_ > a.runid_) return false;
    else if(runid_ == a.runid_){
      if(eventid_ <= a.eventid_) return true;
      else return false;
    }
    else return false;
  }
  
  unsigned int eventid_;
  unsigned int runid_;

};


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



void makeComparisonData(string inputDIR_1, string inputDIR_2, Sample sample, Category category, string outputDIR, bool useOnlyICHEP){

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

  TFile* outputFile = new TFile((outputDIR+"/outputTreeReduced.root").c_str(),"RECREATE");

  TTree* tree_1 = new TTree("tree_1","tree_1");

  float bosonpt = 0, bosoneta = 0, bosonphi = 0, recoil = 0, recoilphi = 0, pfmet = 0, pfmetphi = 0, leadingjetpt = 0, leadingjeteta = 0, secondjetpt = 0, secondjeteta = 0, jetmetdphi = 0;
  unsigned int njetcentral = 0, njetinclusive = 0, runnumber = 0, eventnumber = 0;

  tree_1->Branch("bosonpt",&bosonpt,"bosonpt/F");
  tree_1->Branch("bosoneta",&bosoneta,"bosoneta/F");
  tree_1->Branch("bosonphi",&bosoneta,"bosonphi/F");
  tree_1->Branch("recoil",&recoil,"recoil/F");
  tree_1->Branch("recoilphi",&recoilphi,"recoilphi/F");
  tree_1->Branch("pfmet",&pfmet,"pfmet/F");
  tree_1->Branch("pfmetphi",&pfmetphi,"pfmetphi/F");
  tree_1->Branch("leadingjetpt",&leadingjetpt,"leadingjetpt/F");
  tree_1->Branch("leadingjeteta",&leadingjeteta,"leadingjeteta/F");
  tree_1->Branch("secondjetpt",&secondjetpt,"secondjetpt/F");
  tree_1->Branch("secondjeteta",&secondjeteta,"secondjeteta/F");
  tree_1->Branch("jetmetdphi",&jetmetdphi,"jetmetdphi/F");
  tree_1->Branch("njetcentral",&njetcentral,"njetcentral/i");
  tree_1->Branch("njetinclusive",&njetinclusive,"njetinclusive/i");
  tree_1->Branch("runnumber",&runnumber,"runnumber/i");
  tree_1->Branch("eventnumber",&eventnumber,"eventnumber/i");

  ///////////
  cout<<"Build first chain with index "<<endl;
  TTreeReader reader_1(chain_1);
  // general info                                                                                                                                                                                     
  TTreeReaderValue<unsigned int> run    (reader_1,"run");
  TTreeReaderValue<unsigned int> lumi   (reader_1,"lumi");
  TTreeReaderValue<unsigned int> event  (reader_1,"event");
  TTreeReaderValue<UChar_t> hltm90      (reader_1,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100     (reader_1,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110     (reader_1,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120     (reader_1,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120   (reader_1,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170   (reader_1,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300   (reader_1,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90    (reader_1,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hlte        (reader_1,"hltsingleel");
  TTreeReaderValue<UChar_t> hltenoiso   (reader_1,"hltelnoiso");
  TTreeReaderValue<UChar_t> hltm        (reader_1,"hltsinglemu");
  TTreeReaderValue<UChar_t> hltp165     (reader_1,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175     (reader_1,"hltphoton175");
  TTreeReaderValue<UChar_t> fhbhe  (reader_1,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader_1,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsct  (reader_1,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader_1,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader_1,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader_1,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader_1,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (reader_1,"flagbadchpf");
  TTreeReaderValue<UChar_t> fcsc   (reader_1,"flagglobaltighthalo");
  TTreeReaderValue<unsigned int> njets      (reader_1,"njets");
  TTreeReaderValue<unsigned int> nincjets   (reader_1,"njetsinc");
  TTreeReaderValue<unsigned int> ntausraw   (reader_1,"ntausrawold");
  TTreeReaderValue<unsigned int> nbjets     (reader_1,"nbjetslowpt");
  TTreeReaderValue<unsigned int> nmuons     (reader_1,"nmuons");
  TTreeReaderValue<unsigned int> nphotons     (reader_1,"nphotons");
  TTreeReaderValue<unsigned int> nelectrons     (reader_1,"nelectrons");
  TTreeReaderValue<vector<float> > boostedJetpt    (reader_1,"boostedJetpt");
  TTreeReaderValue<vector<float> > boostedJeteta   (reader_1,"boostedJeteta");
  TTreeReaderValue<vector<float> > prunedJetm      (reader_1,"prunedJetm");
  TTreeReaderValue<vector<float> > boostedJettau2  (reader_1,"boostedJettau2");
  TTreeReaderValue<vector<float> > boostedJettau1  (reader_1,"boostedJettau1");
  TTreeReaderValue<vector<float> > jeteta  (reader_1,"combinejeteta");
  TTreeReaderValue<vector<float> > jetpt   (reader_1,"combinejetpt");
  TTreeReaderValue<vector<float> > chfrac  (reader_1,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac  (reader_1,"combinejetNHfrac");
  TTreeReaderValue<float> met         (reader_1,"t1pfmet");
  TTreeReaderValue<float> metcalo     (reader_1,"calomet");
  TTreeReaderValue<float> metphi      (reader_1,"t1pfmetphi");
  TTreeReaderValue<float> mmet        (reader_1,"t1mumet");
  TTreeReaderValue<float> mmetphi     (reader_1,"t1mumetphi");
  TTreeReaderValue<float> emet        (reader_1,"t1elmet");
  TTreeReaderValue<float> emetphi     (reader_1,"t1elmetphi");
  TTreeReaderValue<float> pmet        (reader_1,"t1phmet");
  TTreeReaderValue<float> pmetphi     (reader_1,"t1phmetphi");
  TTreeReaderValue<float> jmmdphi (reader_1,"incjetmumetdphimin4");
  TTreeReaderValue<float> jemdphi (reader_1,"incjetelmetdphimin4");
  TTreeReaderValue<float> jpmdphi (reader_1,"incjetphmetdphimin4");
  TTreeReaderValue<int>   mu1pid (reader_1,"mu1pid");
  TTreeReaderValue<int>   mu2pid (reader_1,"mu2pid");
  TTreeReaderValue<int>   mu1id  (reader_1,"mu1id");
  TTreeReaderValue<int>   mu2id  (reader_1,"mu2id");
  TTreeReaderValue<float> mu1pt  (reader_1,"mu1pt");
  TTreeReaderValue<float> mu2pt  (reader_1,"mu2pt");
  TTreeReaderValue<float> mu1eta (reader_1,"mu1eta");
  TTreeReaderValue<float> mu2eta (reader_1,"mu2eta");
  TTreeReaderValue<float> mu1phi (reader_1,"mu1phi");
  TTreeReaderValue<float> mu2phi (reader_1,"mu2phi");
  TTreeReaderValue<int>   el1pid (reader_1,"el1pid");
  TTreeReaderValue<int>   el2pid (reader_1,"el2pid");
  TTreeReaderValue<int>   el1id  (reader_1,"el1id");
  TTreeReaderValue<int>   el2id  (reader_1,"el2id");
  TTreeReaderValue<float> el1pt  (reader_1,"el1pt");
  TTreeReaderValue<float> el2pt  (reader_1,"el2pt");
  TTreeReaderValue<float> el1eta (reader_1,"el1eta");
  TTreeReaderValue<float> el2eta (reader_1,"el2eta");
  TTreeReaderValue<float> el1phi (reader_1,"el1phi");
  TTreeReaderValue<float> el2phi (reader_1,"el2phi");
  TTreeReaderValue<int>   phidm  (reader_1,"phidm");
  TTreeReaderValue<float> phpt   (reader_1,"phpt");
  TTreeReaderValue<float> pheta  (reader_1,"pheta");
  TTreeReaderValue<float> phphi  (reader_1,"phphi");
  TTreeReaderValue<float> wmt    (reader_1,"wmt");
  TTreeReaderValue<float> wemt   (reader_1,"wemt");
  TTreeReaderValue<float> zmass  (reader_1,"zmass");
  TTreeReaderValue<float> zeemass(reader_1,"zeemass");
  TTreeReaderValue<float> zmmpt  (reader_1,"zpt");
  TTreeReaderValue<float> zeept  (reader_1,"zeept");
  TTreeReaderValue<float> zeeeta (reader_1,"zeeeta");
  TTreeReaderValue<float> zeephi (reader_1,"zeephi");
  TTreeReaderValue<float> zmmeta (reader_1,"zeta");
  TTreeReaderValue<float> zmmphi (reader_1,"zphi");

  cout<<"Number of events in chain 1 "<<chain_1->GetEntries()<<endl;

  long int nEvents = 0;
  int nPart        = 10000;
  long int nTotal  = chain_1->GetEntries();

  vector<eventID> eventPassing_chain1;

  while(reader_1.Next()){

    bosonpt = 0; bosoneta = 0; bosonphi = 0; recoil = 0; recoilphi = 0; pfmet = 0; pfmetphi = 0; leadingjetpt = 0; leadingjeteta = 0; secondjetpt = 0; secondjeteta = 0;
    njetcentral = 0; njetinclusive = 0; runnumber = 0; eventnumber = 0;
    

    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events chain 1 "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    if(useOnlyICHEP and *run > 276811) continue;

    // recoil
    if((sample == Sample::wmn or sample == Sample::zmm or sample == Sample::sig) and *mmet < recoilSelection) continue;
    else if((sample == Sample::wen or sample == Sample::zee) and *emet < recoilSelection) continue;
    else if(sample == Sample::gam and *pmet < recoilSelection) continue;

    /// jets
    if(*njets < 1 and (category == Category::monojet or category == Category::monoV)) continue;
    else if(category == Category::VBF and *nincjets < 2) continue;

    // vetos
    if(*ntausraw != 0) continue;
    if(*nbjets != 0) continue;

    /// apply met filters
    if(not *fhbhe or not *fhbiso or not *fcsct or not *feeb or not *feeb or not *fetp or not *fvtx or not *fbadmu or not *fbadch or not *fcsc) continue;

        ////////
    int hlt = 0;
    if(sample == Sample::sig or sample == Sample::wmn or sample == Sample::zmm)
      hlt = *hltm90+*hltm100+*hltm110+*hltm120+*hltmwm120+*hltmwm170+*hltmwm300+*hltmwm90;
    else if(sample == Sample::wen or sample == Sample::zee)
      hlt = *hlte+*hltenoiso;
    else if(sample == Sample::gam)
      hlt = *hltp165+*hltp175;
    
    if(not hlt) continue;
     
    // tag objets
    if(sample == Sample::wen and (*el1pt < 40 or *el1id != 1 or *wemt > 160 or *met < 50 or *nmuons != 0 or *nphotons != 0 or *nelectrons != 1)) continue;
    else if(sample == Sample::wmn and (*mu1pt < 20 or *mu1id != 1 or *wmt > 160 or *nelectrons != 0 or *nphotons != 0 or *nmuons != 1)) continue;
    else if(sample == Sample::gam and (*phpt < 175 or fabs(*pheta) > 1.442 or *phidm != 1 or *nmuons != 0 or *nelectrons != 0 or *nphotons != 1)) continue;
    else if(sample == Sample::zmm and (not ((*mu1pt > 20 and *mu1id == 1 and fabs(*mu1eta) < 2.4) or (*mu2pt > 20 and *mu2id == 1 and fabs(*mu1eta) < 2.4)) or 
				       *zmass < 60 or *zmass > 120 or *mu1id == *mu2id or *nelectrons != 0 or *nphotons != 0 or *nmuons != 2)) continue;
    else if(sample == Sample::zee and (not ((*el1pt > 40 and *el1id == 1 and fabs(*el1eta) < 2.5) or (*el2pt > 40 and *el2id == 1 and fabs(*el1eta) < 2.45)) or 
				       *zeemass < 60 or *zeemass > 120 or *el1id == *el2id or *nmuons != 0 or *nphotons != 0 or *nelectrons != 2)) continue;
    // min-dphi
    if((sample == Sample::wmn or sample == Sample::zmm or sample == Sample::sig) and *jmmdphi < 0.5) continue;
    else if((sample == Sample::wen or sample == Sample::zee) and *jemdphi < 0.5) continue;
    else if(sample == Sample::gam and *jpmdphi < 0.5) continue;


    //
    Double_t metden = 0.0;    
    if (sample == Sample::sig || sample == Sample::qcd) {metden = *mmet;}
    else if (sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu){ metden = *mmet;}
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel){ metden = *emet;}
    else if (sample == Sample::qcdgam || sample == Sample::gam)  { metden = *pmet;}    
    if(fabs(*met-*metcalo)/metden > 0.5) continue;

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
    }

    eventPassing_chain1.push_back(eventID(*event,*run));

    if (sample == Sample::sig || sample == Sample::qcd){ 
      bosonpt = *mmet; bosoneta = 0; bosonphi = *mmetphi; recoil = *mmet; recoilphi = *mmetphi; pfmet = *mmet; pfmetphi = *mmetphi; 
      leadingjetpt = jetpt->at(0); leadingjeteta =  jeteta->at(0); njetcentral = *njets; njetinclusive = *nincjets; runnumber = *run; eventnumber = *event;
      if(jetpt->size() > 1){
	secondjetpt = jetpt->at(1);
	secondjeteta = jeteta->at(1);
      }
      jetmetdphi = *jmmdphi;
    }
    else if (sample == Sample::zmm){
      TLorentzVector mu1, mu2;
      mu1.SetPtEtaPhiM(*mu1pt,*mu1eta,*mu1phi,0);
      mu2.SetPtEtaPhiM(*mu2pt,*mu2eta,*mu2phi,0);
      bosonpt = (mu1+mu2).Pt(); bosoneta = (mu1+mu2).Eta(); bosonphi = (mu1+mu2).Phi(); recoil = *mmet; recoilphi = *mmetphi; pfmet = *met; pfmetphi = *metphi; 
      leadingjetpt = jetpt->at(0); leadingjeteta =  jeteta->at(0); njetcentral = *njets; njetinclusive = *nincjets; runnumber = *run; eventnumber = *event;
      if(jetpt->size() > 1){
	secondjetpt = jetpt->at(1);
	secondjeteta = jeteta->at(1);
      }      
      jetmetdphi = *jmmdphi;
    }
    else if (sample == Sample::zee){
      TLorentzVector el1, el2;
      el1.SetPtEtaPhiM(*el1pt,*el1eta,*el1phi,0);
      el2.SetPtEtaPhiM(*el2pt,*el2eta,*el2phi,0);
      bosonpt = (el1+el2).Pt(); bosoneta = (el1+el2).Eta(); bosonphi = (el1+el2).Phi(); recoil = *emet; recoilphi = *emetphi; pfmet = *met; pfmetphi = *metphi; 
      leadingjetpt = jetpt->at(0); leadingjeteta =  jeteta->at(0); njetcentral = *njets; njetinclusive = *nincjets; runnumber = *run; eventnumber = *event;
      if(jetpt->size() > 1){
	secondjetpt = jetpt->at(1);
	secondjeteta = jeteta->at(1);
      }      
      jetmetdphi = *jemdphi;
    }
    else if (sample == Sample::wen){
      bosonpt = *el1pt; bosoneta = *el1eta; bosonphi = *el1phi; recoil = *emet; recoilphi = *emetphi; pfmet = *met; pfmetphi = *metphi; 
      leadingjetpt = jetpt->at(0); leadingjeteta =  jeteta->at(0); njetcentral = *njets; njetinclusive = *nincjets; runnumber = *run; eventnumber = *event;
      if(jetpt->size() > 1){
	secondjetpt = jetpt->at(1);
	secondjeteta = jeteta->at(1);
      }      
      jetmetdphi = *jemdphi;
    }
    else if (sample == Sample::wmn){
      bosonpt = *mu1pt; bosoneta = *mu1eta; bosonphi = *mu1phi; recoil = *mmet; recoilphi = *mmetphi; pfmet = *met; pfmetphi = *metphi; 
      leadingjetpt = jetpt->at(0); leadingjeteta =  jeteta->at(0); njetcentral = *njets; njetinclusive = *nincjets; runnumber = *run; eventnumber = *event;
      if(jetpt->size() > 1){
	secondjetpt = jetpt->at(1);
	secondjeteta = jeteta->at(1);
      }      
      jetmetdphi = *jmmdphi;
    }
    else if (sample == Sample::gam){
      bosonpt = *phpt; bosoneta = *pheta; bosonphi = *phphi; recoil = *pmet; recoilphi = *pmetphi; pfmet = *met; pfmetphi = *metphi; 
      leadingjetpt = jetpt->at(0); leadingjeteta =  jeteta->at(0); njetcentral = *njets; njetinclusive = *nincjets; runnumber = *run; eventnumber = *event;
      if(jetpt->size() > 1){
	secondjetpt = jetpt->at(1);
	secondjeteta = jeteta->at(1);
      }      
      jetmetdphi = *jpmdphi;
    }
    tree_1->Fill();
  }
  cout<<endl;
  cout<<"Number of events passing selections in chain_1 is "<<eventPassing_chain1.size()<<endl;
  tree_1->Write();

  ///////////
  TTreeReader reader_2(chain_2);
  // general info                                                                                                                                                                                     
  TTreeReaderValue<unsigned int> run_alt    (reader_2,"run");
  TTreeReaderValue<unsigned int> lumi_alt   (reader_2,"lumi");
  TTreeReaderValue<unsigned int> event_alt  (reader_2,"event");
  TTreeReaderValue<UChar_t> hltm90_alt      (reader_2,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100_alt     (reader_2,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110_alt     (reader_2,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120_alt     (reader_2,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120_alt   (reader_2,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170_alt   (reader_2,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300_alt   (reader_2,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90_alt    (reader_2,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hlte_alt        (reader_2,"hltsingleel");
  TTreeReaderValue<UChar_t> hltenoiso_alt   (reader_2,"hltelnoiso");
  TTreeReaderValue<UChar_t> hltm_alt        (reader_2,"hltsinglemu");
  TTreeReaderValue<UChar_t> hltp165_alt     (reader_2,"hltphoton165");
  TTreeReaderValue<UChar_t> hltp175_alt     (reader_2,"hltphoton175");
  TTreeReaderValue<UChar_t> fhbhe_alt  (reader_2,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso_alt (reader_2,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsct_alt  (reader_2,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb_alt   (reader_2,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp_alt   (reader_2,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx_alt   (reader_2,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu_alt (reader_2,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch_alt (reader_2,"flagbadchpf");
  TTreeReaderValue<UChar_t> fcsc_alt   (reader_2,"flagglobaltighthalo");
  TTreeReaderValue<unsigned int> njets_alt      (reader_2,"njets");
  TTreeReaderValue<unsigned int> nincjets_alt   (reader_2,"njetsinc");
  TTreeReaderValue<unsigned int> ntausraw_alt   (reader_2,"ntausraw");
  TTreeReaderValue<unsigned int> nbjets_alt     (reader_2,"nbjetslowpt");
  TTreeReaderValue<unsigned int> nmuons_alt     (reader_2,"nmuons");
  TTreeReaderValue<unsigned int> nphotons_alt     (reader_2,"nphotons");
  TTreeReaderValue<unsigned int> nelectrons_alt     (reader_2,"nelectrons");
  TTreeReaderValue<vector<double> > boostedJetpt_alt    (reader_2,"boostedJetpt");
  TTreeReaderValue<vector<double> > boostedJeteta_alt   (reader_2,"boostedJeteta");
  TTreeReaderValue<vector<double> > prunedJetm_alt      (reader_2,"prunedJetm");
  TTreeReaderValue<vector<double> > boostedJettau2_alt  (reader_2,"boostedJettau2");
  TTreeReaderValue<vector<double> > boostedJettau1_alt  (reader_2,"boostedJettau1");
  TTreeReaderValue<vector<double> > jeteta_alt  (reader_2,"combinejeteta");
  TTreeReaderValue<vector<double> > jetpt_alt   (reader_2,"combinejetpt");
  TTreeReaderValue<vector<double> > chfrac_alt  (reader_2,"combinejetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac_alt  (reader_2,"combinejetNHfrac");
  TTreeReaderValue<double> met_alt         (reader_2,"t1pfmet");
  TTreeReaderValue<double> metcalo_alt     (reader_2,"calomet");
  TTreeReaderValue<double> metphi_alt      (reader_2,"t1pfmetphi");
  TTreeReaderValue<double> mmet_alt        (reader_2,"t1mumet");
  TTreeReaderValue<double> mmetphi_alt     (reader_2,"t1mumetphi");
  TTreeReaderValue<double> emet_alt        (reader_2,"t1elmet");
  TTreeReaderValue<double> emetphi_alt     (reader_2,"t1elmetphi");
  TTreeReaderValue<double> pmet_alt        (reader_2,"t1phmet");
  TTreeReaderValue<double> pmetphi_alt     (reader_2,"t1phmetphi");
  TTreeReaderValue<double> jmmdphi_alt (reader_2,"incjetmumetdphimin4");
  TTreeReaderValue<double> jemdphi_alt (reader_2,"incjetelmetdphimin4");
  TTreeReaderValue<double> jpmdphi_alt (reader_2,"incjetphmetdphimin4");
  TTreeReaderValue<int>   mu1pid_alt (reader_2,"mu1pid");
  TTreeReaderValue<int>   mu2pid_alt (reader_2,"mu2pid");
  TTreeReaderValue<int>   mu1id_alt  (reader_2,"mu1id");
  TTreeReaderValue<int>   mu2id_alt  (reader_2,"mu2id");
  TTreeReaderValue<double> mu1pt_alt  (reader_2,"mu1pt");
  TTreeReaderValue<double> mu2pt_alt  (reader_2,"mu2pt");
  TTreeReaderValue<double> mu1eta_alt (reader_2,"mu1eta");
  TTreeReaderValue<double> mu2eta_alt (reader_2,"mu2eta");
  TTreeReaderValue<double> mu1phi_alt (reader_2,"mu1phi");
  TTreeReaderValue<double> mu2phi_alt (reader_2,"mu2phi");
  TTreeReaderValue<int>   el1pid_alt (reader_2,"el1pid");
  TTreeReaderValue<int>   el2pid_alt (reader_2,"el2pid");
  TTreeReaderValue<int>   el1id_alt  (reader_2,"el1id");
  TTreeReaderValue<int>   el2id_alt  (reader_2,"el2id");
  TTreeReaderValue<double> el1pt_alt  (reader_2,"el1pt");
  TTreeReaderValue<double> el2pt_alt  (reader_2,"el2pt");
  TTreeReaderValue<double> el1eta_alt (reader_2,"el1eta");
  TTreeReaderValue<double> el2eta_alt (reader_2,"el2eta");
  TTreeReaderValue<double> el1phi_alt (reader_2,"el1phi");
  TTreeReaderValue<double> el2phi_alt (reader_2,"el2phi");
  TTreeReaderValue<int>    phidm_alt  (reader_2,"phidm");
  TTreeReaderValue<double> phpt_alt   (reader_2,"phpt");
  TTreeReaderValue<double> pheta_alt  (reader_2,"pheta");
  TTreeReaderValue<double> phphi_alt  (reader_2,"phphi");
  TTreeReaderValue<double> wmt_alt    (reader_2,"wmt");
  TTreeReaderValue<double> wemt_alt   (reader_2,"wemt");
  TTreeReaderValue<double> zmass_alt  (reader_2,"zmass");
  TTreeReaderValue<double> zeemass_alt (reader_2,"zeemass");
  TTreeReaderValue<double> zmmpt_alt  (reader_2,"zpt");
  TTreeReaderValue<double> zeept_alt  (reader_2,"zeept");
  TTreeReaderValue<double> zeeeta_alt (reader_2,"zeeeta");
  TTreeReaderValue<double> zeephi_alt (reader_2,"zeephi");
  TTreeReaderValue<double> zmmeta_alt (reader_2,"zeta");
  TTreeReaderValue<double> zmmphi_alt (reader_2,"zphi");


  TTree* tree_2 = new TTree("tree_2","tree_2");

  bosonpt = 0, bosoneta = 0, bosonphi = 0, recoil = 0, recoilphi = 0, pfmet = 0, pfmetphi = 0, leadingjetpt = 0, leadingjeteta = 0, secondjetpt = 0, secondjeteta = 0; jetmetdphi = 0;
  njetcentral = 0, njetinclusive = 0, runnumber = 0, eventnumber = 0;

  tree_2->Branch("bosonpt",&bosonpt,"bosonpt/F");
  tree_2->Branch("bosoneta",&bosoneta,"bosoneta/F");
  tree_2->Branch("bosonphi",&bosoneta,"bosonphi/F");
  tree_2->Branch("recoil",&recoil,"recoil/F");
  tree_2->Branch("recoilphi",&recoilphi,"recoilphi/F");
  tree_2->Branch("pfmet",&pfmet,"pfmet/F");
  tree_2->Branch("pfmetphi",&pfmetphi,"pfmetphi/F");
  tree_2->Branch("leadingjetpt",&leadingjetpt,"leadingjetpt/F");
  tree_2->Branch("leadingjeteta",&leadingjeteta,"leadingjeteta/F");
  tree_2->Branch("secondjetpt",&secondjetpt,"secondjetpt/F");
  tree_2->Branch("secondjeteta",&secondjeteta,"secondjeteta/F");
  tree_2->Branch("jetmetdphi",&jetmetdphi,"jetmetdphi/F");
  tree_2->Branch("njetcentral",&njetcentral,"njetcentral/i");
  tree_2->Branch("njetinclusive",&njetinclusive,"njetinclusive/i");
  tree_2->Branch("runnumber",&runnumber,"runnumber/i");
  tree_2->Branch("eventnumber",&eventnumber,"eventnumber/i");

  cout<<"Number of events in chain 2 "<<chain_2->GetEntries()<<endl;
  nEvents = 0;
  nPart        = 10000;
  nTotal  = chain_2->GetEntries();

  vector<eventID> eventPassing_chain2;

  while(reader_2.Next()){

    bosonpt = 0, bosoneta = 0, bosonphi = 0, recoil = 0, recoilphi = 0, pfmet = 0, pfmetphi = 0, leadingjetpt = 0, leadingjeteta = 0, secondjetpt = 0, secondjeteta = 0; jetmetdphi = 0;
    njetcentral = 0, njetinclusive = 0, runnumber = 0, eventnumber = 0;

    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events chain 1 "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    if(useOnlyICHEP and *run_alt > 276811) continue;

    // recoil
    if((sample == Sample::wmn or sample == Sample::zmm or sample == Sample::sig) and *mmet_alt < recoilSelection) continue;
    else if((sample == Sample::wen or sample == Sample::zee) and *emet_alt < recoilSelection) continue;
    else if(sample == Sample::gam and *pmet_alt < recoilSelection) continue;

    /// jets
    if(*njets_alt < 1 and (category == Category::monojet or category == Category::monoV)) continue;
    else if(category == Category::VBF and *nincjets_alt < 2) continue;

    // vetos
    if(*ntausraw_alt != 0) continue;
    if(*nbjets_alt != 0) continue;

    /// apply met filters
    if(not *fhbhe_alt or not *fhbiso_alt or not *fcsct_alt or not *feeb_alt or not *feeb_alt or not *fetp_alt or not *fvtx_alt or not *fbadmu_alt or not *fbadch_alt or not *fcsc_alt) continue;

        ////////
    int hlt = 0;
    if(sample == Sample::sig or sample == Sample::wmn or sample == Sample::zmm)
      hlt = *hltm90_alt+*hltm100_alt+*hltm110_alt+*hltm120_alt+*hltmwm120_alt+*hltmwm170_alt+*hltmwm300_alt+*hltmwm90_alt;
    else if(sample == Sample::wen or sample == Sample::zee)
      hlt = *hlte_alt+*hltenoiso_alt;
    else if(sample == Sample::gam)
      hlt = *hltp165_alt+*hltp175_alt;
    
    if(not hlt) continue;
     
    // tag objets
    if(sample == Sample::wen and (*el1pt_alt < 40 or *el1id_alt != 1 or *wemt_alt > 160 or *met_alt < 50 or *nmuons_alt != 0 or *nphotons_alt != 0 or *nelectrons_alt != 1)) continue;
    else if(sample == Sample::wmn and (*mu1pt_alt < 20 or *mu1id_alt != 1 or *wmt_alt > 160 or *nelectrons_alt != 0 or *nphotons_alt != 0 or *nmuons_alt != 1)) continue;
    else if(sample == Sample::gam and (*phpt_alt < 175 or fabs(*pheta_alt) > 1.442 or *phidm_alt != 1 or *nmuons_alt != 0 or *nelectrons_alt != 0 or *nphotons_alt != 1)) continue;
    else if(sample == Sample::zmm and (not ((*mu1pt_alt > 20 and *mu1id_alt == 1 and fabs(*mu1eta_alt) < 2.4) or (*mu2pt_alt > 20 and *mu2id_alt == 1 and fabs(*mu1eta_alt) < 2.4)) or 
				       *zmass_alt < 60 or *zmass_alt > 120 or *mu1id_alt == *mu2id_alt or *nelectrons_alt != 0 or *nphotons_alt != 0 or *nmuons_alt != 2)) continue;
    else if(sample == Sample::zee and (not ((*el1pt_alt > 40 and *el1id_alt == 1 and fabs(*el1eta_alt) < 2.5) or (*el2pt_alt > 40 and *el2id_alt == 1 and fabs(*el1eta_alt) < 2.45)) or 
				       *zeemass_alt < 60 or *zeemass_alt > 120 or *el1id_alt == *el2id_alt or *nmuons_alt != 0 or *nphotons_alt != 0 or *nelectrons_alt != 2)) continue;
    // min-dphi
    if((sample == Sample::wmn or sample == Sample::zmm or sample == Sample::sig) and *jmmdphi_alt < 0.5) continue;
    else if((sample == Sample::wen or sample == Sample::zee) and *jemdphi_alt < 0.5) continue;
    else if(sample == Sample::gam and *jpmdphi_alt < 0.5) continue;


    //
    Double_t metden = 0.0;
    if (sample == Sample::sig || sample == Sample::qcd) {metden = *mmet_alt;}
    else if (sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu){ metden = *mmet_alt;}
    else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel){ metden = *emet_alt;}
    else if (sample == Sample::qcdgam || sample == Sample::gam)  { metden = *pmet_alt;}    
    if(fabs(*met_alt-*metcalo_alt)/metden > 0.5) continue;

    if(category == Category::monojet){

      if(jetpt_alt->at(0) < 100) continue;
      if(fabs(jeteta_alt->at(0)) > 2.5) continue;
      if(chfrac_alt->at(0) < 0.1) continue;
      if(nhfrac_alt->at(0) > 0.8) continue;
      bool goodMonoV = false;

      if(boostedJetpt_alt->size() != 0 and boostedJetpt_alt->at(0) > 250 and fabs(boostedJeteta_alt->at(0)) < 2.4 and boostedJettau2_alt->at(0)/boostedJettau1_alt->at(0) < 0.6 and 
	 prunedJetm_alt->at(0) > 65 and prunedJetm_alt->at(0) < 105 and metden > 250)
	goodMonoV = true;
      if(goodMonoV) continue;     
    }

    eventPassing_chain2.push_back(eventID(*event_alt,*run_alt));

    if (sample == Sample::sig || sample == Sample::qcd){ 
      bosonpt = *mmet_alt; bosoneta = 0; bosonphi = *mmetphi_alt; recoil = *mmet_alt; recoilphi = *mmetphi_alt; pfmet = *mmet_alt; pfmetphi = *mmetphi_alt; 
      leadingjetpt = jetpt_alt->at(0); leadingjeteta =  jeteta_alt->at(0); njetcentral = *njets_alt; njetinclusive = *nincjets_alt; runnumber = *run_alt; eventnumber = *event_alt;
      if(jetpt_alt->size() > 1){
	secondjetpt = jetpt_alt->at(1);
	secondjeteta = jeteta_alt->at(1);
      }
      jetmetdphi = *jmmdphi_alt;
    }
    else if (sample == Sample::zmm){
      TLorentzVector mu1, mu2;
      mu1.SetPtEtaPhiM(*mu1pt_alt,*mu1eta_alt,*mu1phi_alt,0);
      mu2.SetPtEtaPhiM(*mu2pt_alt,*mu2eta_alt,*mu2phi_alt,0);
      bosonpt = (mu1+mu2).Pt(); bosoneta = (mu1+mu2).Eta(); bosonphi = (mu1+mu2).Phi(); recoil = *mmet_alt; recoilphi = *mmetphi_alt; pfmet = *met_alt; pfmetphi = *metphi_alt; 
      leadingjetpt = jetpt_alt->at(0); leadingjeteta =  jeteta_alt->at(0); njetcentral = *njets_alt; njetinclusive = *nincjets_alt; runnumber = *run_alt; eventnumber = *event_alt;
      if(jetpt_alt->size() > 1){
	secondjetpt = jetpt_alt->at(1);
	secondjeteta = jeteta_alt->at(1);
      }      
      jetmetdphi = *jmmdphi_alt;
    }
    else if (sample == Sample::zee){
      TLorentzVector el1, el2;
      el1.SetPtEtaPhiM(*el1pt_alt,*el1eta_alt,*el1phi_alt,0);
      el2.SetPtEtaPhiM(*el2pt_alt,*el2eta_alt,*el2phi_alt,0);
      bosonpt = (el1+el2).Pt(); bosoneta = (el1+el2).Eta(); bosonphi = (el1+el2).Phi(); recoil = *emet_alt; recoilphi = *emetphi_alt; pfmet = *met_alt; pfmetphi = *metphi_alt; 
      leadingjetpt = jetpt_alt->at(0); leadingjeteta =  jeteta_alt->at(0); njetcentral = *njets_alt; njetinclusive = *nincjets_alt; runnumber = *run_alt; eventnumber = *event_alt;
      if(jetpt_alt->size() > 1){
	secondjetpt = jetpt_alt->at(1);
	secondjeteta = jeteta_alt->at(1);
      }      
      jetmetdphi = *jemdphi_alt;
    }
    else if (sample == Sample::wen){
      bosonpt = *el1pt_alt; bosoneta = *el1eta_alt; bosonphi = *el1phi_alt; recoil = *emet_alt; recoilphi = *emetphi_alt; pfmet = *met_alt; pfmetphi = *metphi_alt; 
      leadingjetpt = jetpt_alt->at(0); leadingjeteta =  jeteta_alt->at(0); njetcentral = *njets_alt; njetinclusive = *nincjets_alt; runnumber = *run_alt; eventnumber = *event_alt;
      if(jetpt_alt->size() > 1){
	secondjetpt = jetpt_alt->at(1);
	secondjeteta = jeteta_alt->at(1);
      }      
      jetmetdphi = *jemdphi_alt;
    }
    else if (sample == Sample::wmn){
      bosonpt = *mu1pt_alt; bosoneta = *mu1eta_alt; bosonphi = *mu1phi_alt; recoil = *mmet_alt; recoilphi = *mmetphi_alt; pfmet = *met_alt; pfmetphi = *metphi_alt; 
      leadingjetpt = jetpt_alt->at(0); leadingjeteta =  jeteta_alt->at(0); njetcentral = *njets_alt; njetinclusive = *nincjets_alt; runnumber = *run_alt; eventnumber = *event_alt;
      if(jetpt_alt->size() > 1){
	secondjetpt = jetpt_alt->at(1);
	secondjeteta = jeteta_alt->at(1);
      }      
      jetmetdphi = *jmmdphi_alt;
    }
    else if (sample == Sample::gam){
      bosonpt = *phpt_alt; bosoneta = *pheta_alt; bosonphi = *phphi_alt; recoil = *pmet_alt; recoilphi = *pmetphi_alt; pfmet = *met_alt; pfmetphi = *metphi_alt; 
      leadingjetpt = jetpt_alt->at(0); leadingjeteta =  jeteta_alt->at(0); njetcentral = *njets_alt; njetinclusive = *nincjets_alt; runnumber = *run_alt; eventnumber = *event_alt;
      if(jetpt_alt->size() > 1){
	secondjetpt = jetpt_alt->at(1);
	secondjeteta = jeteta_alt->at(1);
      }      
      jetmetdphi = *jpmdphi_alt;
    }
    tree_2->Fill();

  }
  cout<<endl;
  cout<<"Number of events passing selections in chain_2 is "<<eventPassing_chain2.size()<<endl;
  tree_2->Write();

  
  std::sort(eventPassing_chain1.begin(),eventPassing_chain1.end());
  std::sort(eventPassing_chain2.begin(),eventPassing_chain2.end());

  vector<eventID> commonEvents;
  vector<eventID> eventIn1Only;
  vector<eventID> eventIn2Only;

  for(auto event : eventPassing_chain1){
    if(std::find(eventPassing_chain2.begin(),eventPassing_chain2.end(),event) != eventPassing_chain2.end())
      commonEvents.push_back(event);
    else
      eventIn1Only.push_back(event);
  }

  for(auto event : eventPassing_chain2){
    if(std::find(eventPassing_chain1.begin(),eventPassing_chain1.end(),event) != eventPassing_chain1.end()){
      if(std::find(commonEvents.begin(),commonEvents.end(),event) == commonEvents.end()) cout<<"Problem in finding common events please check "<<endl;
    }
    else
      eventIn2Only.push_back(event);
  }

  cout<<"Common events among the two datasets "<<commonEvents.size()<<" i.e. "<<100*double(commonEvents.size())/eventPassing_chain1.size()<<" % of dataset 1 and "<<100*double(commonEvents.size())/eventPassing_chain2.size()<<" % of dataset 2"<<endl;
  cout<<"Events in 1 but non in 2 "<<eventIn1Only.size()<<" i.e. "<<100*double(eventIn1Only.size())/eventPassing_chain1.size()<<" % of dataset 1"<<endl;
  cout<<"Events in 2 but non in 1 "<<eventIn2Only.size()<<" i.e. "<<100*double(eventIn2Only.size())/eventPassing_chain2.size()<<" % of dataset 2"<<endl;


  // Loop over common events and fill diff histograms
  bosonpt_histo->Sumw2();
  bosoneta_histo->Sumw2();
  bosonphi_histo->Sumw2();
  recoil_histo->Sumw2();
  recoilphi_histo->Sumw2();
  pfmet_histo->Sumw2();
  pfmetphi_histo->Sumw2();
  jetpt_histo->Sumw2();
  jetpt2_histo->Sumw2();
  njet_histo->Sumw2();
  njetinc_histo->Sumw2();
  jeteta_histo->Sumw2();
  jeteta2_histo->Sumw2();
  jetmetdphi_histo->Sumw2();

  bosonpt_histo_2D->Sumw2();
  bosoneta_histo_2D->Sumw2();
  bosonphi_histo_2D->Sumw2();
  recoil_histo_2D->Sumw2();
  recoilphi_histo_2D->Sumw2();
  pfmet_histo_2D->Sumw2();
  pfmetphi_histo_2D->Sumw2();
  jetpt_histo_2D->Sumw2();
  jetpt2_histo_2D->Sumw2();
  jeteta_histo_2D->Sumw2();
  jeteta2_histo_2D->Sumw2();
  njet_histo_2D->Sumw2();
  njetinc_histo_2D->Sumw2();
  jetmetdphi_histo_2D->Sumw2();

  bosonpt_histo_only1->Sumw2();
  bosoneta_histo_only1->Sumw2();
  bosonphi_histo_only1->Sumw2();
  recoil_histo_only1->Sumw2();
  recoilphi_histo_only1->Sumw2();
  pfmet_histo_only1->Sumw2();
  pfmetphi_histo_only1->Sumw2();
  jetpt_histo_only1->Sumw2();
  jetpt2_histo_only1->Sumw2();
  njet_histo_only1->Sumw2();
  njetinc_histo_only1->Sumw2();
  jeteta_histo_only1->Sumw2();
  jeteta2_histo_only1->Sumw2();
  jetmetdphi_histo_only1->Sumw2();

  bosonpt_histo_only2->Sumw2();
  bosoneta_histo_only2->Sumw2();
  bosonphi_histo_only2->Sumw2();
  recoil_histo_only2->Sumw2();
  recoilphi_histo_only2->Sumw2();
  pfmet_histo_only2->Sumw2();
  pfmetphi_histo_only2->Sumw2();
  jetpt_histo_only2->Sumw2();
  jetpt2_histo_only2->Sumw2();
  njet_histo_only2->Sumw2();
  njetinc_histo_only2->Sumw2();
  jeteta_histo_only2->Sumw2();
  jeteta2_histo_only2->Sumw2();
  jetmetdphi_histo_only2->Sumw2();

  // loop on tree_1 and tree_2 just checking common events
  TTreeReader myreader1 ("tree_1",outputFile);
  TTreeReaderValue<float> bosonpt_1 (myreader1,"bosonpt");
  TTreeReaderValue<float> bosoneta_1 (myreader1,"bosoneta");
  TTreeReaderValue<float> bosonphi_1 (myreader1,"bosonphi");
  TTreeReaderValue<float> recoil_1 (myreader1,"recoil");
  TTreeReaderValue<float> recoilphi_1 (myreader1,"recoilphi");
  TTreeReaderValue<float> pfmet_1 (myreader1,"pfmet");
  TTreeReaderValue<float> pfmetphi_1 (myreader1,"pfmetphi");
  TTreeReaderValue<float> leadingjetpt_1 (myreader1,"leadingjetpt");
  TTreeReaderValue<float> leadingjeteta_1 (myreader1,"leadingjeteta");
  TTreeReaderValue<float> secondjetpt_1 (myreader1,"secondjetpt");
  TTreeReaderValue<float> secondjeteta_1 (myreader1,"secondjeteta");
  TTreeReaderValue<float> jetmetdphi_1 (myreader1,"jetmetdphi");
  TTreeReaderValue<unsigned int> njetcentral_1   (myreader1,"njetcentral");
  TTreeReaderValue<unsigned int> njetinclusive_1 (myreader1,"njetinclusive");
  TTreeReaderValue<unsigned int> runnumber_1     (myreader1,"runnumber");
  TTreeReaderValue<unsigned int> eventnumber_1   (myreader1,"eventnumber");

  float bosonpt_2, bosoneta_2, bosonphi_2, recoil_2, recoilphi_2, pfmet_2, leadingjetpt_2, leadingjeteta_2, secondjetpt_2, secondjeteta_2, jetmetdphi_2, pfmetphi_2;
  unsigned int njetcentral_2, njetinclusive_2, runnumber_2, eventnumber_2;

  tree_2->SetBranchAddress("bosonpt",&bosonpt_2);
  tree_2->SetBranchAddress("bosoneta",&bosoneta_2);
  tree_2->SetBranchAddress("bosonphi",&bosonphi_2);
  tree_2->SetBranchAddress("recoil",&recoil_2);
  tree_2->SetBranchAddress("recoilphi",&recoilphi_2);
  tree_2->SetBranchAddress("pfmet",&pfmet_2);
  tree_2->SetBranchAddress("pfmetphi",&pfmetphi_2);
  tree_2->SetBranchAddress("leadingjetpt",&leadingjetpt_2);
  tree_2->SetBranchAddress("leadingjeteta",&leadingjeteta_2);
  tree_2->SetBranchAddress("secondjetpt",&secondjetpt_2);
  tree_2->SetBranchAddress("secondjeteta",&secondjeteta_2);
  tree_2->SetBranchAddress("jetmetdphi",&jetmetdphi_2);
  tree_2->SetBranchAddress("njetcentral",&njetcentral_2);
  tree_2->SetBranchAddress("njetinclusive",&njetinclusive_2);
  tree_2->SetBranchAddress("runnumber",&runnumber_2);
  tree_2->SetBranchAddress("eventnumber",&eventnumber_2);

  while(myreader1.Next()){
    eventID event_tmp(*eventnumber_1,*runnumber_1);
    auto iterator = std::find(commonEvents.begin(),commonEvents.end(),event_tmp);
    if(iterator != commonEvents.end()){
      // find a good event
      for(long int iEvent = 0; iEvent < tree_2->GetEntries(); iEvent++){
	tree_2->GetEntry(iEvent);
	if(event_tmp.eventid_ != eventnumber_2) continue;
	if(event_tmp.runid_ != runnumber_2) continue;
	// when found
	bosonpt_histo_2D->Fill(*bosonpt_1,bosonpt_2);
	bosoneta_histo_2D->Fill(*bosoneta_1,bosoneta_2);
	bosonphi_histo_2D->Fill(*bosonphi_1,bosonphi_2);
	recoil_histo_2D->Fill(*recoil_1,recoil_2);
	recoilphi_histo_2D->Fill(*recoilphi_1,recoilphi_2);
	pfmet_histo_2D->Fill(*pfmet_1,pfmet_2);
	pfmetphi_histo_2D->Fill(*pfmetphi_1,pfmetphi_2);
	jetpt_histo_2D->Fill(*leadingjetpt_1,leadingjetpt_2);
	jetpt2_histo_2D->Fill(*secondjetpt_1,secondjetpt_2);
	jeteta_histo_2D->Fill(*leadingjeteta_1,leadingjeteta_2);
	jeteta2_histo_2D->Fill(*secondjeteta_1,secondjeteta_2);
	jetmetdphi_histo_2D->Fill(*jetmetdphi_1,jetmetdphi_2);	
	njet_histo_2D->Fill(float(*njetcentral_1),float(njetcentral_2));
	njetinc_histo_2D->Fill(float(*njetinclusive_1),float(njetinclusive_2));

	bosonpt_histo->Fill(*bosonpt_1-bosonpt_2);
	bosoneta_histo->Fill(*bosoneta_1-bosoneta_2);
	bosonphi_histo->Fill(*bosonphi_1-bosonphi_2);
	recoil_histo->Fill(*recoil_1-recoil_2);
	recoilphi_histo->Fill(*recoilphi_1-recoilphi_2);
	pfmet_histo->Fill(*pfmet_1-pfmet_2);
	pfmetphi_histo->Fill(*pfmetphi_1-pfmetphi_2);
	jetpt_histo->Fill(*leadingjetpt_1-leadingjetpt_2);
	jetpt2_histo->Fill(*secondjetpt_1-secondjetpt_2);
	jeteta_histo->Fill(*leadingjeteta_1-leadingjeteta_2);
	jeteta2_histo->Fill(*secondjeteta_1-secondjeteta_2);
	jetmetdphi_histo->Fill(*jetmetdphi_1-jetmetdphi_2);	
	njet_histo->Fill(float(*njetcentral_1)-float(njetcentral_2));
	njetinc_histo->Fill(float(*njetinclusive_1)-float(njetinclusive_2));

	break;
      }
    }
    else{
      auto iterator = std::find(eventIn1Only.begin(),eventIn1Only.end(),event_tmp);
      if(iterator != eventIn1Only.end()){
	bosonpt_histo_only1->Fill(*bosonpt_1);
	bosoneta_histo_only1->Fill(*bosoneta_1);
	bosonphi_histo_only1->Fill(*bosonphi_1);
	recoil_histo_only1->Fill(*recoil_1);
	recoilphi_histo_only1->Fill(*recoilphi_1);
	pfmet_histo_only1->Fill(*pfmet_1);
	pfmetphi_histo_only1->Fill(*pfmetphi_1);
	jetpt_histo_only1->Fill(*leadingjetpt_1);
	jetpt2_histo_only1->Fill(*secondjetpt_1);
	jeteta_histo_only1->Fill(*leadingjeteta_1);
	jeteta2_histo_only1->Fill(*secondjeteta_1);
	jetmetdphi_histo_only1->Fill(*jetmetdphi_1);
	njet_histo_only1->Fill(float(*njetcentral_1));
	njetinc_histo_only1->Fill(float(*njetinclusive_1));
      }
      else{cout<<" Event not found --> is a problem "<<endl;}
    }
  }

  for(long int iEvent = 0; iEvent < tree_2->GetEntries(); iEvent++){
    tree_2->GetEntry(iEvent);
    eventID event_tmp(eventnumber_2,runnumber_2);
    auto iterator = std::find(eventIn2Only.begin(),eventIn2Only.end(),event_tmp);
    if(iterator != eventIn2Only.end()){
	bosonpt_histo_only2->Fill(bosonpt_2);
	bosoneta_histo_only2->Fill(bosoneta_2);
	bosonphi_histo_only2->Fill(bosonphi_2);
	recoil_histo_only2->Fill(recoil_2);
	recoilphi_histo_only2->Fill(recoilphi_2);
	pfmet_histo_only2->Fill(pfmet_2);
	pfmetphi_histo_only2->Fill(pfmetphi_2);
	jetpt_histo_only2->Fill(leadingjetpt_2);
	jetpt2_histo_only2->Fill(secondjetpt_2);
	jeteta_histo_only2->Fill(leadingjeteta_2);
	jeteta2_histo_only2->Fill(secondjeteta_2);
	jetmetdphi_histo_only2->Fill(jetmetdphi_2);
	njet_histo_only2->Fill(float(njetcentral_2));
	njetinc_histo_only2->Fill(float(njetinclusive_2));
    }    
  }

  // Draw histograms
  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();

  drawPlot(canvas,bosonpt_histo,"boson p_{T}",outputDIR);  
  drawPlot(canvas,bosoneta_histo,"boson #eta",outputDIR);  
  drawPlot(canvas,bosonphi_histo,"boson #phi",outputDIR);  
  drawPlot(canvas,recoil_histo,"Recoil",outputDIR);  
  drawPlot(canvas,recoilphi_histo,"Recoil #phi",outputDIR);  
  drawPlot(canvas,pfmet_histo,"Met",outputDIR);  
  drawPlot(canvas,pfmetphi_histo,"Met #phi",outputDIR);  
  drawPlot(canvas,jetpt_histo,"p_{T}^{j1}",outputDIR);  
  drawPlot(canvas,jetpt2_histo,"p_{T}^{j2}",outputDIR);  
  drawPlot(canvas,jeteta_histo,"#eta_{j1}",outputDIR);  
  drawPlot(canvas,jeteta2_histo,"#eta_{j2}",outputDIR);  
  drawPlot(canvas,njet_histo,"N_{jet}",outputDIR);  
  drawPlot(canvas,njetinc_histo,"N_{jet}",outputDIR);  
  drawPlot(canvas,jetmetdphi_histo,"#Delta#phi(jet,met)",outputDIR);  
  
  drawPlot2D(canvas,bosonpt_histo_2D,"boson p_{T}",outputDIR);  
  drawPlot2D(canvas,bosoneta_histo_2D,"boson #eta",outputDIR);  
  drawPlot2D(canvas,bosonphi_histo_2D,"boson #phi",outputDIR);  
  drawPlot2D(canvas,recoil_histo_2D,"Recoil",outputDIR);  
  drawPlot2D(canvas,recoilphi_histo_2D,"Recoil #phi",outputDIR);  
  drawPlot2D(canvas,pfmet_histo_2D,"Met",outputDIR);  
  drawPlot2D(canvas,pfmetphi_histo_2D,"Met #phi",outputDIR);  
  drawPlot2D(canvas,jetpt_histo_2D,"p_{T}^{j1}",outputDIR);  
  drawPlot2D(canvas,jetpt2_histo_2D,"p_{T}^{j2}",outputDIR);  
  drawPlot2D(canvas,jeteta_histo_2D,"#eta_{j1}",outputDIR);  
  drawPlot2D(canvas,jeteta2_histo_2D,"#eta_{j2}",outputDIR);  
  drawPlot2D(canvas,njet_histo_2D,"N_{jet}",outputDIR);  
  drawPlot2D(canvas,njetinc_histo_2D,"N_{jet}",outputDIR);  
  drawPlot2D(canvas,jetmetdphi_histo_2D,"#Delta#phi(jet,met)",outputDIR);  

  /////////////////////
  drawPlot(canvas,bosonpt_histo_only1,"boson p_{T}",outputDIR);  
  drawPlot(canvas,bosoneta_histo_only1,"boson #eta",outputDIR);  
  drawPlot(canvas,bosonphi_histo_only1,"boson #phi",outputDIR);  
  drawPlot(canvas,recoil_histo_only1,"Recoil",outputDIR);  
  drawPlot(canvas,recoilphi_histo_only1,"Recoil #phi",outputDIR);  
  drawPlot(canvas,pfmet_histo_only1,"Met",outputDIR);  
  drawPlot(canvas,pfmetphi_histo_only1,"Met #phi",outputDIR);  
  drawPlot(canvas,jetpt_histo_only1,"p_{T}^{j1}",outputDIR);  
  drawPlot(canvas,jetpt2_histo_only1,"p_{T}^{j2}",outputDIR);  
  drawPlot(canvas,jeteta_histo_only1,"#eta_{j1}",outputDIR);  
  drawPlot(canvas,jeteta2_histo_only1,"#eta_{j2}",outputDIR);  
  drawPlot(canvas,njet_histo_only1,"N_{jet}",outputDIR);  
  drawPlot(canvas,njetinc_histo_only1,"N_{jet}",outputDIR);  
  drawPlot(canvas,jetmetdphi_histo_only1,"#Delta#phi(jet,met)",outputDIR);  

  drawPlot(canvas,bosonpt_histo_only2,"boson p_{T}",outputDIR);  
  drawPlot(canvas,bosoneta_histo_only2,"boson #eta",outputDIR);  
  drawPlot(canvas,bosonphi_histo_only2,"boson #phi",outputDIR);  
  drawPlot(canvas,recoil_histo_only2,"Recoil",outputDIR);  
  drawPlot(canvas,recoilphi_histo_only2,"Recoil #phi",outputDIR);  
  drawPlot(canvas,pfmet_histo_only2,"Met",outputDIR);  
  drawPlot(canvas,pfmetphi_histo_only2,"Met #phi",outputDIR);  
  drawPlot(canvas,jetpt_histo_only2,"p_{T}^{j1}",outputDIR);  
  drawPlot(canvas,jetpt2_histo_only2,"p_{T}^{j2}",outputDIR);  
  drawPlot(canvas,jeteta_histo_only2,"#eta_{j1}",outputDIR);  
  drawPlot(canvas,jeteta2_histo_only2,"#eta_{j2}",outputDIR);  
  drawPlot(canvas,njet_histo_only2,"N_{jet}",outputDIR);  
  drawPlot(canvas,jetmetdphi_histo_only2,"#Delta#phi(jet,met)",outputDIR);  

  /// Make a larger check --> take eventIn1Only and see if these events are in chain_2 and if they fail a specific cut somewhere
  cout<<"Event found in chain_1 but not in chain_2 --> check what selection are failing"<<endl;
  reader_2.SetEntry(0);
  TH1F* selections = new TH1F("selections","",12,0,12);
  selections->GetXaxis()->SetBinLabel(1,"Found events "); selections->SetBinContent(1,0);
  selections->GetXaxis()->SetBinLabel(2,"Recoil       "); selections->SetBinContent(2,0);
  selections->GetXaxis()->SetBinLabel(3,"Njet         "); selections->SetBinContent(3,0);
  selections->GetXaxis()->SetBinLabel(4,"Tau-veto     "); selections->SetBinContent(4,0);
  selections->GetXaxis()->SetBinLabel(5,"B-veto       "); selections->SetBinContent(5,0);
  selections->GetXaxis()->SetBinLabel(6,"Met filters  "); selections->SetBinContent(6,0);
  selections->GetXaxis()->SetBinLabel(7,"Trigger      "); selections->SetBinContent(7,0);
  selections->GetXaxis()->SetBinLabel(8,"Photon-id    "); selections->SetBinContent(8,0);
  selections->GetXaxis()->SetBinLabel(9,"min-dphi     "); selections->SetBinContent(9,0);
  selections->GetXaxis()->SetBinLabel(10,"metcalo      "); selections->SetBinContent(10,0);
  selections->GetXaxis()->SetBinLabel(11,"mono-jet id  "); selections->SetBinContent(11,0);
  selections->GetXaxis()->SetBinLabel(12,"mono-V veto  "); selections->SetBinContent(12,0);

  nEvents=0;
  nTotal = chain_2->GetEntries();
  while(reader_2.Next()){

    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events chain 2 "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    eventID event_tmp(*event_alt,*run_alt);
    auto iterator = std::find(eventIn1Only.begin(),eventIn1Only.end(),event_tmp);
    if(iterator != eventIn1Only.end()){
      // found event in the trees
      selections->SetBinContent(1,selections->GetBinContent(1)+1);      

      // recoil                                                                                                                                                                                   
      if((sample == Sample::wmn or sample == Sample::zmm or sample == Sample::sig) and *mmet_alt < recoilSelection) continue;
      else if((sample == Sample::wen or sample == Sample::zee) and *emet_alt < recoilSelection) continue;
      else if(sample == Sample::gam and *pmet_alt < recoilSelection) continue;
      selections->SetBinContent(2,selections->GetBinContent(2)+1);

      /// jets                                                                                                                                                                                 
      if(*njets_alt < 1 and (category == Category::monojet or category == Category::monoV)) continue;
      else if(category == Category::VBF and *nincjets_alt < 2) continue;
      selections->SetBinContent(3,selections->GetBinContent(3)+1);

      // vetos                                                                                                                                                                                     
      if(*ntausraw_alt != 0) continue;
      selections->SetBinContent(4,selections->GetBinContent(4)+1);
      if(*nbjets_alt != 0) continue;
      selections->SetBinContent(5,selections->GetBinContent(5)+1);

      /// apply met filters                                                                                                                                                                         
      if(not *fhbhe_alt or not *fhbiso_alt or not *fcsct_alt or not *feeb_alt or not *feeb_alt or not *fetp_alt or not *fvtx_alt or not *fbadmu_alt or not *fbadch_alt or not *fcsc_alt) continue;
      selections->SetBinContent(6,selections->GetBinContent(6)+1);

      ////////                                                                                                                                                                                       
      int hlt = 0;
      if(sample == Sample::sig or sample == Sample::wmn or sample == Sample::zmm)
	hlt = *hltm90_alt+*hltm100_alt+*hltm110_alt+*hltm120_alt+*hltmwm120_alt+*hltmwm170_alt+*hltmwm300_alt+*hltmwm90_alt;
      else if(sample == Sample::wen or sample == Sample::zee)
	hlt = *hlte_alt+*hltenoiso_alt;
      else if(sample == Sample::gam)
	hlt = *hltp165_alt+*hltp175_alt;

      if(not hlt) continue;
      selections->SetBinContent(7,selections->GetBinContent(7)+1);

      // tag objets                                                                                                                                                                   
      if(sample == Sample::wen and (*el1pt_alt < 40 or *el1id_alt != 1 or *wemt_alt > 160 or *met_alt < 50 or *nmuons_alt != 0 or *nphotons_alt != 0 or *nelectrons_alt != 1)) continue;
      else if(sample == Sample::wmn and (*mu1pt_alt < 20 or *mu1id_alt != 1 or *wmt_alt > 160 or *nelectrons_alt != 0 or *nphotons_alt != 0 or *nmuons_alt != 1)) continue;
      else if(sample == Sample::gam and (*phpt_alt < 175 or fabs(*pheta_alt) > 1.442 or *phidm_alt != 1 or *nmuons_alt != 0 or *nelectrons_alt != 0 or *nphotons_alt != 1)) continue;
      else if(sample == Sample::zmm and (not ((*mu1pt_alt > 20 and *mu1id_alt == 1 and fabs(*mu1eta_alt) < 2.4) or (*mu2pt_alt > 20 and *mu2id_alt == 1 and fabs(*mu1eta_alt) < 2.4)) or
					 *zmass_alt < 60 or *zmass_alt > 120 or *mu1id_alt == *mu2id_alt or *nelectrons_alt != 0 or *nphotons_alt != 0 or *nmuons_alt != 2)) continue;
      else if(sample == Sample::zee and (not ((*el1pt_alt > 40 and *el1id_alt == 1 and fabs(*el1eta_alt) < 2.5) or (*el2pt_alt > 40 and *el2id_alt == 1 and fabs(*el1eta_alt) < 2.45)) or
					 *zeemass_alt < 60 or *zeemass_alt > 120 or *el1id_alt == *el2id_alt or *nmuons_alt != 0 or *nphotons_alt != 0 or *nelectrons_alt != 2)) continue;

      selections->SetBinContent(8,selections->GetBinContent(8)+1);

      // min-dphi                                                                                                                                                                                     
      if((sample == Sample::wmn or sample == Sample::zmm or sample == Sample::sig) and *jmmdphi_alt < 0.5) continue;
      else if((sample == Sample::wen or sample == Sample::zee) and *jemdphi_alt < 0.5) continue;
      else if(sample == Sample::gam and *jpmdphi_alt < 0.5) continue;

      selections->SetBinContent(9,selections->GetBinContent(9)+1);

      //                                                                                                                                                                                        
      Double_t metden = 0.0;
      if (sample == Sample::sig || sample == Sample::qcd) {metden = *mmet_alt;}
      else if (sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu){ metden = *mmet_alt;}
      else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel){ metden = *emet_alt;}
      else if (sample == Sample::qcdgam || sample == Sample::gam)  { metden = *pmet_alt;}
      if(fabs(*met_alt-*metcalo_alt)/metden > 0.5) continue;

      selections->SetBinContent(10,selections->GetBinContent(10)+1);


      if(category == Category::monojet){

	if(jetpt_alt->at(0) < 100) continue;
	if(fabs(jeteta_alt->at(0)) > 2.5) continue;
	if(chfrac_alt->at(0) < 0.1) continue;
	if(nhfrac_alt->at(0) > 0.8) continue;
	selections->SetBinContent(11,selections->GetBinContent(11)+1);
	
	bool goodMonoV = false;

	if(boostedJetpt_alt->size() != 0 and boostedJetpt_alt->at(0) > 250 and fabs(boostedJeteta_alt->at(0)) < 2.4 and boostedJettau2_alt->at(0)/boostedJettau1_alt->at(0) < 0.6 and
	   prunedJetm_alt->at(0) > 65 and prunedJetm_alt->at(0) < 105 and metden > 250)
	  goodMonoV = true;
	if(goodMonoV) continue;

	selections->SetBinContent(12,selections->GetBinContent(12)+1);
      }
    }
  }
  cout<<endl;
  cout<<"Total events = "<<eventIn1Only.size()<<endl;
  for(int iBin = 1; iBin <= selections->GetNbinsX(); iBin++){
    if(iBin != 1)
      cout<<selections->GetXaxis()->GetBinLabel(iBin)<<" = "<<selections->GetBinContent(iBin)<<" reduction "<<fabs(selections->GetBinContent(iBin)-selections->GetBinContent(iBin-1))/eventIn1Only.size()<<endl;
    else
      cout<<selections->GetXaxis()->GetBinLabel(iBin)<<" = "<<selections->GetBinContent(iBin)<<endl;
  }
  
  /// Make a larger check --> take eventIn2Only and see if these events are in chain_1 and if they fail a specific cut somewhere
  cout<<"Event found in chain_2 but not in chain_1 --> check what selection are failing"<<endl;
  reader_1.SetEntry(0);
  selections->Reset();
  selections->GetXaxis()->SetBinLabel(1,"Found events "); selections->SetBinContent(1,0);
  selections->GetXaxis()->SetBinLabel(2,"Recoil       "); selections->SetBinContent(2,0);
  selections->GetXaxis()->SetBinLabel(3,"Njet         "); selections->SetBinContent(3,0);
  selections->GetXaxis()->SetBinLabel(4,"Tau-veto     "); selections->SetBinContent(4,0);
  selections->GetXaxis()->SetBinLabel(5,"B-veto       "); selections->SetBinContent(5,0);
  selections->GetXaxis()->SetBinLabel(6,"Met filters  "); selections->SetBinContent(6,0);
  selections->GetXaxis()->SetBinLabel(7,"Trigger      "); selections->SetBinContent(7,0);
  selections->GetXaxis()->SetBinLabel(8,"Photon-id    "); selections->SetBinContent(8,0);
  selections->GetXaxis()->SetBinLabel(9,"min-dphi     "); selections->SetBinContent(9,0);
  selections->GetXaxis()->SetBinLabel(10,"metcalo      "); selections->SetBinContent(10,0);
  selections->GetXaxis()->SetBinLabel(11,"mono-jet id  "); selections->SetBinContent(11,0);
  selections->GetXaxis()->SetBinLabel(12,"mono-V veto  "); selections->SetBinContent(12,0);
  nEvents=0;
  nTotal = chain_1->GetEntries();
  while(reader_1.Next()){

    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events chain 1 "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    eventID event_tmp(*event,*run);
    auto iterator = std::find(eventIn2Only.begin(),eventIn2Only.end(),event_tmp);
    if(iterator != eventIn2Only.end()){
      // found event in the trees
      selections->SetBinContent(1,selections->GetBinContent(1)+1);      

      // recoil                                                                                                                                                                                   
      if((sample == Sample::wmn or sample == Sample::zmm or sample == Sample::sig) and *mmet < recoilSelection) continue;
      else if((sample == Sample::wen or sample == Sample::zee) and *emet < recoilSelection) continue;
      else if(sample == Sample::gam and *pmet < recoilSelection) continue;
      selections->SetBinContent(2,selections->GetBinContent(2)+1);

      /// jets                                                                                                                                                                                 
      if(*njets < 1 and (category == Category::monojet or category == Category::monoV)) continue;
      else if(category == Category::VBF and *nincjets < 2) continue;
      selections->SetBinContent(3,selections->GetBinContent(3)+1);

      // vetos                                                                                                                                                                                     
      if(*ntausraw != 0) continue;
      selections->SetBinContent(4,selections->GetBinContent(4)+1);
      if(*nbjets != 0) continue;
      selections->SetBinContent(5,selections->GetBinContent(5)+1);

      /// apply met filters                                                                                                                                                                         
      if(not *fhbhe or not *fhbiso or not *fcsct or not *feeb or not *feeb or not *fetp or not *fvtx or not *fbadmu or not *fbadch or not *fcsc) continue;
      selections->SetBinContent(6,selections->GetBinContent(6)+1);

      ////////                                                                                                                                                                                       
      int hlt = 0;
      if(sample == Sample::sig or sample == Sample::wmn or sample == Sample::zmm)
	hlt = *hltm90+*hltm100+*hltm110+*hltm120+*hltmwm120+*hltmwm170+*hltmwm300+*hltmwm90;
      else if(sample == Sample::wen or sample == Sample::zee)
	hlt = *hlte+*hltenoiso;
      else if(sample == Sample::gam)
	hlt = *hltp165+*hltp175;

      if(not hlt) continue;
      selections->SetBinContent(7,selections->GetBinContent(7)+1);

      // tag objets                                                                                                                                                                   
      if(sample == Sample::wen and (*el1pt < 40 or *el1id != 1 or *wemt > 160 or *met < 50 or *nmuons != 0 or *nphotons != 0 or *nelectrons != 1)) continue;
      else if(sample == Sample::wmn and (*mu1pt < 20 or *mu1id != 1 or *wmt > 160 or *nelectrons != 0 or *nphotons != 0 or *nmuons != 1)) continue;
      else if(sample == Sample::gam and (*phpt < 175 or fabs(*pheta) > 1.442 or *phidm != 1 or *nmuons != 0 or *nelectrons != 0 or *nphotons != 1)) continue;
      else if(sample == Sample::zmm and (not ((*mu1pt > 20 and *mu1id == 1 and fabs(*mu1eta) < 2.4) or (*mu2pt > 20 and *mu2id == 1 and fabs(*mu1eta) < 2.4)) or
					 *zmass < 60 or *zmass > 120 or *mu1id == *mu2id or *nelectrons != 0 or *nphotons != 0 or *nmuons != 2)) continue;
      else if(sample == Sample::zee and (not ((*el1pt > 40 and *el1id == 1 and fabs(*el1eta) < 2.5) or (*el2pt > 40 and *el2id == 1 and fabs(*el1eta) < 2.45)) or
					 *zeemass < 60 or *zeemass > 120 or *el1id == *el2id or *nmuons != 0 or *nphotons != 0 or *nelectrons != 2)) continue;

      selections->SetBinContent(8,selections->GetBinContent(8)+1);

      // min-dphi                                                                                                                                                                                     
      if((sample == Sample::wmn or sample == Sample::zmm or sample == Sample::sig) and *jmmdphi < 0.5) continue;
      else if((sample == Sample::wen or sample == Sample::zee) and *jemdphi < 0.5) continue;
      else if(sample == Sample::gam and *jpmdphi < 0.5) continue;

      selections->SetBinContent(9,selections->GetBinContent(9)+1);

      //                                                                                                                                                                                        
      Double_t metden = 0.0;
      if (sample == Sample::sig || sample == Sample::qcd) {metden = *mmet;}
      else if (sample == Sample::zmm || sample == Sample::wmn || sample == Sample::topmu){ metden = *mmet;}
      else if (sample == Sample::zee || sample == Sample::wen || sample == Sample::topel){ metden = *emet;}
      else if (sample == Sample::qcdgam || sample == Sample::gam)  { metden = *pmet;}
      if(fabs(*met-*metcalo)/metden > 0.5) continue;

      selections->SetBinContent(10,selections->GetBinContent(10)+1);


      if(category == Category::monojet){

	if(jetpt->at(0) < 100) continue;
	if(fabs(jeteta->at(0)) > 2.5) continue;
	if(chfrac->at(0) < 0.1) continue;
	if(nhfrac->at(0) > 0.8) continue;
	selections->SetBinContent(11,selections->GetBinContent(11)+1);
	
	bool goodMonoV = false;

	if(boostedJetpt->size() != 0 and boostedJetpt->at(0) > 250 and fabs(boostedJeteta->at(0)) < 2.4 and boostedJettau2->at(0)/boostedJettau1->at(0) < 0.6 and
	   prunedJetm->at(0) > 65 and prunedJetm->at(0) < 105 and metden > 250)
	  goodMonoV = true;
	if(goodMonoV) continue;

	selections->SetBinContent(12,selections->GetBinContent(12)+1);
      }
    }
  }
  cout<<endl;
  cout<<"Total events = "<<eventIn2Only.size()<<endl;
  for(int iBin = 1; iBin <= selections->GetNbinsX(); iBin++){
    if(iBin != 1)
      cout<<selections->GetXaxis()->GetBinLabel(iBin)<<" = "<<selections->GetBinContent(iBin)<<" reduction "<<fabs(selections->GetBinContent(iBin)-selections->GetBinContent(iBin-1))/eventIn2Only.size()<<endl;
    else
      cout<<selections->GetXaxis()->GetBinLabel(iBin)<<" = "<<selections->GetBinContent(iBin)<<endl;
  }
  
  
  outputFile->Close();

}
