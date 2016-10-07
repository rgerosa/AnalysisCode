#include "../CMS_lumi.h"

static bool doVBF = true;

void fillHisto(TH1* histo, double val){ // embed the overflow
  
  if(val < histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1))
    histo->Fill(val);
  else
    histo->Fill(histo->GetXaxis()->GetBinCenter(histo->GetNbinsX()));

}

void fillHisto(TH2* histo, double valx, double valy){ // Embed- the overflow
  
  double x = 0;
  if(valx < histo->GetXaxis()->GetBinLowEdge(histo->GetNbinsX()+1))
    x = valx;
  else
    x =histo->GetXaxis()->GetBinCenter(histo->GetNbinsX());
  double y = 0;
  if(valy < histo->GetYaxis()->GetBinLowEdge(histo->GetNbinsY()+1))
    y = valy;
  else
    y = histo->GetYaxis()->GetBinCenter(histo->GetNbinsY());

  histo->Fill(x,y);
}


void plotDistribution(TCanvas* canvas, TH1* qqH_m1, TH1* qqH_m2, TH1* qqH_m3,  const string & label, const string & outputDIR){
  
  // normalize to 1
  canvas->cd();
  qqH_m1->Scale(1./qqH_m1->Integral());
  qqH_m2->Scale(1./qqH_m2->Integral());
  qqH_m3->Scale(1./qqH_m3->Integral());
  qqH_m1->GetXaxis()->SetTitle(label.c_str());
  qqH_m1->GetYaxis()->SetTitle("a.u");
  qqH_m1->SetLineColor(kBlack);
  qqH_m1->SetLineWidth(2);
  qqH_m1->Draw("hist");

  qqH_m2->SetLineColor(kBlue);
  qqH_m3->SetLineColor(kRed);
  qqH_m2->Draw("hist same");
  qqH_m3->Draw("hist same");

  qqH_m1->GetYaxis()->SetRangeUser(0.,max(qqH_m1->GetMaximum(),max(qqH_m2->GetMaximum(),qqH_m3->GetMaximum()))*1.2);
  CMS_lumi(canvas,"");
  TLegend leg(0.7,0.7,0.9,0.9);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  if(doVBF){
    leg.AddEntry(qqH_m1,"VBF H = 125 GeV","L");
    leg.AddEntry(qqH_m2,"VBF H = 400 GeV","L");
    leg.AddEntry(qqH_m3,"VBF H = 800 GeV","L");
  }
  else{
    leg.AddEntry(qqH_m1,"ggH = 125 GeV","L");
    leg.AddEntry(qqH_m2,"ggH = 400 GeV","L");
    leg.AddEntry(qqH_m3,"ggH = 800 GeV","L");
  }

  leg.Draw("same");
  
  canvas->SaveAs((outputDIR+"/"+string(qqH_m1->GetName())+".pdf").c_str(),"pdf");
  qqH_m1->GetYaxis()->SetRangeUser(max(0.0001,min(qqH_m1->GetMinimum(),min(qqH_m2->GetMinimum(),qqH_m3->GetMinimum()))*0.8),
				   max(qqH_m1->GetMaximum(),max(qqH_m2->GetMaximum(),qqH_m3->GetMaximum()))*1000);
  canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/"+string(qqH_m1->GetName())+"_log.pdf").c_str(),"pdf");
  canvas->SetLogy(0);
  
}


void plotDistribution(TCanvas* canvas, TH2* qqH_m1, TH2* qqH_m2, TH2* qqH_m3, const string & labelX, const string & labelY, const string & outputDIR){
  
  system(("mkdir -p "+outputDIR+"/correlation/").c_str());
  // normalize to 1
  canvas->cd();
  canvas->SetRightMargin(0.18);
  qqH_m1->Scale(1./qqH_m1->Integral());
  qqH_m2->Scale(1./qqH_m2->Integral());
  qqH_m3->Scale(1./qqH_m3->Integral());

  TProfile* qqHProfile_m1 = qqH_m1->ProfileX(Form("%s_pfx",qqH_m1->GetName()));
  qqHProfile_m1->SetMarkerColor(kBlack);
  qqHProfile_m1->SetMarkerStyle(20);
  qqHProfile_m1->SetMarkerSize(1);
  qqHProfile_m1->GetXaxis()->SetTitle(labelX.c_str());
  qqHProfile_m1->GetYaxis()->SetTitle(labelY.c_str());
  qqHProfile_m1->Draw("EP");

  TProfile* qqHProfile_m2 = qqH_m2->ProfileX(Form("%s_pfx",qqH_m2->GetName()));
  qqHProfile_m2->SetMarkerColor(kBlue);
  qqHProfile_m2->SetMarkerStyle(20);
  qqHProfile_m2->SetMarkerSize(1);
  qqHProfile_m2->Draw("EPsame");

  TProfile* qqHProfile_m3 = qqH_m3->ProfileX(Form("%s_pfx",qqH_m3->GetName()));
  qqHProfile_m3->SetMarkerColor(kRed);
  qqHProfile_m3->SetMarkerStyle(20);
  qqHProfile_m3->SetMarkerSize(1);
  qqHProfile_m3->Draw("EPsame");

  CMS_lumi(canvas,"");
  TLegend leg(0.6,0.6,0.9,0.9);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  if(doVBF){
    leg.AddEntry(qqHProfile_m1,"VBF H = 125 GeV","P");
    leg.AddEntry(qqHProfile_m2,"VBF H = 400 GeV","P");
    leg.AddEntry(qqHProfile_m3,"VBF H = 800 GeV","P");
  }
  else{
    leg.AddEntry(qqHProfile_m1,"ggH = 125 GeV","P");
    leg.AddEntry(qqHProfile_m2,"ggH = 400 GeV","P");
    leg.AddEntry(qqHProfile_m3,"ggH = 800 GeV","P");
  }
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/correlation/"+string(qqH_m1->GetName())+".pdf").c_str(),"pdf");
  canvas->SetLogz();
  canvas->SaveAs((outputDIR+"/correlation/"+string(qqH_m1->GetName())+"_log.pdf").c_str(),"pdf");
  canvas->SetLogz(0);
  
}

void makeVBFSignalStudies(string outputPlot, bool doVBFPlots = true){

  doVBF = doVBFPlots;

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputPlot).c_str());

  TChain* vbftree_m1 = new TChain("tree/tree");
  TChain* vbftree_m2 = new TChain("tree/tree");
  TChain* vbftree_m3 = new TChain("tree/tree");
  if(doVBFPlots){
    vbftree_m1->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/*VBF*125*root");
    vbftree_m2->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/*VBF*400*root");
    vbftree_m3->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/*VBF*800*root");
  }
  else{
    vbftree_m1->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/*GluGlu*125*root");
    vbftree_m2->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/*GluGlu*400*root");
    vbftree_m3->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/*GluGlu*800*root");
  }

  TTreeReader myReader (vbftree_m1);
  TTreeReaderValue<unsigned int> njetsinc  (myReader,"njetsinc");
  TTreeReaderValue<unsigned int> nphotons  (myReader,"nphotons");
  TTreeReaderValue<unsigned int> nelectrons   (myReader,"nelectrons");
  TTreeReaderValue<unsigned int> ntaus   (myReader,"ntausraw");
  TTreeReaderValue<unsigned int> nmuons (myReader,"nmuons");
  TTreeReaderValue<unsigned int> nbjets  (myReader,"nbjetslowpt");
  TTreeReaderValue<UChar_t> fhbhe   (myReader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso  (myReader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc    (myReader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb    (myReader,"flageebadsc");
  TTreeReaderValue<vector<double> > jetpt    (myReader,"combinejetpt");
  TTreeReaderValue<vector<double> > jeteta  (myReader,"combinejeteta");
  TTreeReaderValue<vector<double> > jetphi  (myReader,"combinejetphi");
  TTreeReaderValue<vector<double> > jetbtag  (myReader,"combinejetbtag");
  TTreeReaderValue<vector<double> > jetm     (myReader,"combinejetm");
  TTreeReaderValue<vector<double> > chfrac   (myReader,"combinejetCHfrac");
  TTreeReaderValue<vector<double> > nhfrac   (myReader,"combinejetNHfrac");
  TTreeReaderValue<double > mediatorMass     (myReader,"dmmass");
  TTreeReaderValue<double > mediatorPt       (myReader,"dmpt");
  TTreeReaderValue<double > mediatorEta      (myReader,"dmeta");
  TTreeReaderValue<double > mediatorPhi      (myReader,"dmphi");
  TTreeReaderValue<double > x1Mass     (myReader,"dmX1mass");
  TTreeReaderValue<double > x1Pt       (myReader,"dmX1pt");
  TTreeReaderValue<double > x1Eta      (myReader,"dmX1eta");
  TTreeReaderValue<double > x1Phi      (myReader,"dmX1phi");
  TTreeReaderValue<double > x2Mass     (myReader,"dmX2mass");
  TTreeReaderValue<double > x2Pt       (myReader,"dmX2pt");
  TTreeReaderValue<double > x2Eta      (myReader,"dmX2eta");
  TTreeReaderValue<double > x2Phi      (myReader,"dmX2phi");
  TTreeReaderValue<double> mmetphi     (myReader,"t1mumetphi");
  TTreeReaderValue<double> mmet        (myReader,"t1mumet");
  TTreeReaderValue<double> genmet      (myReader,"genmet");
  TTreeReaderValue<double> genmetphi   (myReader,"genmetphi");
  TTreeReaderValue<double> jmmdphi4    (myReader,"incjetmumetdphimin4");
  TTreeReaderValue<double> jmmdphi     (myReader,"incjetmumetdphimin");  
  TTreeReaderValue<UChar_t> hltm90     (myReader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm120    (myReader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120  (myReader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm170  (myReader,"hltmetwithmu170");
  TTreeReaderValue<UChar_t> hltmwm300  (myReader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90   (myReader,"hltmetwithmu90");
  
  TH1F* qqH_m1_bosonPt = new TH1F("qqH_m1_bosonPt","qqH_m1_bosonPt",25,150,1250);
  TH1F* qqH_m1_bosonEta = new TH1F("qqH_m1_bosonEta","qqH_m1_bosonEta",25,-3.5,3.5);
  TH1F* qqH_m1_bosonPhi = new TH1F("qqH_m1_bosonPhi","qqH_m1_bosonPhi",25,0.,3.14);
  TH1F* qqH_m1_bosonPz    = new TH1F("qqH_m1_bosonPz","qqH_m1_bosonPz",25,150,1250);
  TH1F* qqH_m1_bosonPzOPt = new TH1F("qqH_m1_bosonPzOPt","qqH_m1_bosonPzOPt",20,0,20);
  TH1F* qqH_m1_jetpt1   = new TH1F("qqH_m1_jetpt1","qqH_m1_jetpt1",25,50,500);
  TH1F* qqH_m1_jetpt2   = new TH1F("qqH_m1_jetpt2","qqH_m1_jetpt2",25,50,500);
  TH1F* qqH_m1_jetpz1   = new TH1F("qqH_m1_jetpz1","qqH_m1_jetpz1",25,50,500);  
  TH1F* qqH_m1_jetpz2   = new TH1F("qqH_m1_jetpz2","qqH_m1_jetpz2",25,50,500);
  TH1F* qqH_m1_jetpz1Opt   = new TH1F("qqH_m1_jetpz1Opt","qqH_m1_jetpz1",20,0,20);  
  TH1F* qqH_m1_jetpz2Opt   = new TH1F("qqH_m1_jetpz2OPt","qqH_m1_jetpz2",20,0,20);
  TH1F* qqH_m1_jeteta1  = new TH1F("qqH_m1_jeteta1","qqH_m1_jeteta1",25,-3.5,3.5);
  TH1F* qqH_m1_jeteta2  = new TH1F("qqH_m1_jeteta2","qqH_m1_jeteta2",25,-3.5,3.5);
  TH1F* qqH_m1_jetphi1  = new TH1F("qqH_m1_jetphi1","qqH_m1_jetphi1",25,0,3.14);
  TH1F* qqH_m1_jetphi2  = new TH1F("qqH_m1_jetphi2","qqH_m1_jetphi2",25,0,3.14);
  TH1F* qqH_m1_jetDeltaPhi  = new TH1F("qqH_m1_jetDeltaPhi","qqH_m1_jetDeltaPhi",25,0,3.14);
  TH1F* qqH_m1_jetDeltaEta  = new TH1F("qqH_m1_jetDeltaEta","qqH_m1_jetDeltaEta",25,0,9);
  TH1F* qqH_m1_Mjj     = new TH1F("qqH_m1_Mjj","qqH_m1_Mjj",40,400,2500);
  TH1F* qqH_m1_jetmetDphi4 = new TH1F("qqH_m1_jetmetdphi4","qqH_m1_jetmetdphi4",25,0,3.14);
  TH1F* qqH_m1_jetmetDphi  = new TH1F("qqH_m1_jetmetdphi","qqH_m1_jetmetdphi",25,0,3.14);
  TH1F* qqH_m1_jetmediatorDphi = new TH1F("qqH_m1_jetmediatorDphi","qqH_m1_jetmediatorDphi",25,0,3.14);
  TH1F* qqH_m1_jetmediatorDeta = new TH1F("qqH_m1_jetmediatorDeta","qqH_m1_jetmediatorDeta",25,0,9);
  TH1F* qqH_m1_MT = new TH1F("qqH_m1_MT","qqH_m1_MT",30,0,600);
  TH1F* qqH_m1_MetMT = new TH1F("qqH_m1_MetMT","qqH_m1_MetMT",30,100,1000);

  qqH_m1_bosonPt->Sumw2();
  qqH_m1_bosonEta->Sumw2();
  qqH_m1_bosonPhi->Sumw2();
  qqH_m1_bosonPz->Sumw2();
  qqH_m1_bosonPzOPt->Sumw2();
  qqH_m1_jetpt1->Sumw2();
  qqH_m1_jetpt2->Sumw2();
  qqH_m1_jetpz1->Sumw2();
  qqH_m1_jetpz2->Sumw2();
  qqH_m1_jetpz1Opt->Sumw2();
  qqH_m1_jetpz2Opt->Sumw2();
  qqH_m1_jeteta1->Sumw2();
  qqH_m1_jeteta2->Sumw2();
  qqH_m1_jetphi1->Sumw2();
  qqH_m1_jetphi2->Sumw2();
  qqH_m1_jetDeltaPhi->Sumw2();
  qqH_m1_jetDeltaEta->Sumw2();
  qqH_m1_Mjj->Sumw2();
  qqH_m1_jetmetDphi->Sumw2();
  qqH_m1_jetmetDphi4->Sumw2();
  qqH_m1_jetmediatorDphi->Sumw2();
  qqH_m1_jetmediatorDeta->Sumw2();
  qqH_m1_MT->Sumw2();
  qqH_m1_MetMT->Sumw2();

  TH2F* qqH_m1_mjj_detajj = new TH2F("qqH_m1_mjj_detajj","qqH_m1_mjj_detajj",15,400,2500,10,2,8);
  TH2F* qqH_m1_mjj_jetmetDphi = new TH2F("qqH_m1_mjj_jetmetDphi","qqH_m1_mjj_jetmetDphi",15,400,2500,10,0.5,3.14);
  TH2F* qqH_m1_mjj_jetDphi = new TH2F("qqH_m1_mjj_jetDphi","qqH_m1_mjj_jetDphi",15,400,2500,10,0.5,3.14);
  TH2F* qqH_m1_mjj_MT = new TH2F("qqH_m1_mjj_MT","qqH_m1_mjj_MT",15,400,2500,10,0,500);
  TH2F* qqH_m1_detajj_jetmetDphi = new TH2F("qqH_m1_detajj_jetmetDphi","qqH_m1_detajj_jetmetDphi",10,2,8,10,0.5,3.14);
  TH2F* qqH_m1_detajj_jetDphi = new TH2F("qqH_m1_detajj_jetDphi","qqH_m1_detajj_jetDphi",10,2,8,10,0.5,3.14);
  TH2F* qqH_m1_detajj_MT = new TH2F("qqH_m1_detajj_MT","qqH_m1_detajj_MT",10,2,8,10,0,500);

  qqH_m1_mjj_detajj->Sumw2();
  qqH_m1_mjj_jetmetDphi->Sumw2();
  qqH_m1_mjj_jetDphi->Sumw2();
  qqH_m1_mjj_MT->Sumw2();
  qqH_m1_detajj_jetmetDphi->Sumw2();
  qqH_m1_detajj_jetDphi->Sumw2();
  qqH_m1_detajj_MT->Sumw2();
  
  ///////////////

  TH1F* qqH_m2_bosonPt = new TH1F("qqH_m2_bosonPt","qqH_m2_bosonPt",25,150,1250);
  TH1F* qqH_m2_bosonEta = new TH1F("qqH_m2_bosonEta","qqH_m2_bosonEta",25,-3.5,3.5);
  TH1F* qqH_m2_bosonPhi = new TH1F("qqH_m2_bosonPhi","qqH_m2_bosonPhi",25,0.,3.14);
  TH1F* qqH_m2_bosonPz    = new TH1F("qqH_m2_bosonPz","qqH_m2_bosonPz",25,150,1250);
  TH1F* qqH_m2_bosonPzOPt = new TH1F("qqH_m2_bosonPzOPt","qqH_m2_bosonPzOPt",20,0,20);
  TH1F* qqH_m2_jetpt1   = new TH1F("qqH_m2_jetpt1","qqH_m2_jetpt1",25,50,500);
  TH1F* qqH_m2_jetpt2   = new TH1F("qqH_m2_jetpt2","qqH_m2_jetpt2",25,50,500);
  TH1F* qqH_m2_jetpz1   = new TH1F("qqH_m2_jetpz1","qqH_m2_jetpz1",25,50,500);  
  TH1F* qqH_m2_jetpz2   = new TH1F("qqH_m2_jetpz2","qqH_m2_jetpz2",25,50,500);
  TH1F* qqH_m2_jetpz1Opt   = new TH1F("qqH_m2_jetpz1Opt","qqH_m2_jetpz1",20,0,20);  
  TH1F* qqH_m2_jetpz2Opt   = new TH1F("qqH_m2_jetpz2OPt","qqH_m2_jetpz2",20,0,20);
  TH1F* qqH_m2_jeteta1  = new TH1F("qqH_m2_jeteta1","qqH_m2_jeteta1",25,-3.5,3.5);
  TH1F* qqH_m2_jeteta2  = new TH1F("qqH_m2_jeteta2","qqH_m2_jeteta2",25,-3.5,3.5);
  TH1F* qqH_m2_jetphi1  = new TH1F("qqH_m2_jetphi1","qqH_m2_jetphi1",25,0,3.14);
  TH1F* qqH_m2_jetphi2  = new TH1F("qqH_m2_jetphi2","qqH_m2_jetphi2",25,0,3.14);
  TH1F* qqH_m2_jetDeltaPhi  = new TH1F("qqH_m2_jetDeltaPhi","qqH_m2_jetDeltaPhi",25,0,3.14);
  TH1F* qqH_m2_jetDeltaEta  = new TH1F("qqH_m2_jetDeltaEta","qqH_m2_jetDeltaEta",25,0,9);
  TH1F* qqH_m2_Mjj     = new TH1F("qqH_m2_Mjj","qqH_m2_Mjj",40,400,2500);
  TH1F* qqH_m2_jetmetDphi4 = new TH1F("qqH_m2_jetmetdphi4","qqH_m2_jetmetdphi4",25,0,3.14);
  TH1F* qqH_m2_jetmetDphi  = new TH1F("qqH_m2_jetmetdphi","qqH_m2_jetmetdphi",25,0,3.14);
  TH1F* qqH_m2_jetmediatorDphi = new TH1F("qqH_m2_jetmediatorDphi","qqH_m2_jetmediatorDphi",25,0,3.14);
  TH1F* qqH_m2_jetmediatorDeta = new TH1F("qqH_m2_jetmediatorDeta","qqH_m2_jetmediatorDeta",25,0,9);
  TH1F* qqH_m2_MT = new TH1F("qqH_m2_MT","qqH_m2_MT",30,0,600);
  TH1F* qqH_m2_MetMT = new TH1F("qqH_m2_MetMT","qqH_m2_MetMT",30,100,1000);

  qqH_m2_bosonPt->Sumw2();
  qqH_m2_bosonEta->Sumw2();
  qqH_m2_bosonPhi->Sumw2();
  qqH_m2_bosonPz->Sumw2();
  qqH_m2_bosonPzOPt->Sumw2();
  qqH_m2_jetpt1->Sumw2();
  qqH_m2_jetpt2->Sumw2();
  qqH_m2_jetpz1->Sumw2();
  qqH_m2_jetpz2->Sumw2();
  qqH_m2_jetpz1Opt->Sumw2();
  qqH_m2_jetpz2Opt->Sumw2();
  qqH_m2_jeteta1->Sumw2();
  qqH_m2_jeteta2->Sumw2();
  qqH_m2_jetphi1->Sumw2();
  qqH_m2_jetphi2->Sumw2();
  qqH_m2_jetDeltaPhi->Sumw2();
  qqH_m2_jetDeltaEta->Sumw2();
  qqH_m2_Mjj->Sumw2();
  qqH_m2_jetmetDphi->Sumw2();
  qqH_m2_jetmetDphi4->Sumw2();
  qqH_m2_jetmediatorDphi->Sumw2();
  qqH_m2_jetmediatorDeta->Sumw2();
  qqH_m2_MT->Sumw2();
  qqH_m2_MetMT->Sumw2();

  TH2F* qqH_m2_mjj_detajj = new TH2F("qqH_m2_mjj_detajj","qqH_m2_mjj_detajj",15,400,2500,10,2,8);
  TH2F* qqH_m2_mjj_jetmetDphi = new TH2F("qqH_m2_mjj_jetmetDphi","qqH_m2_mjj_jetmetDphi",15,400,2500,10,0.5,3.14);
  TH2F* qqH_m2_mjj_jetDphi = new TH2F("qqH_m2_mjj_jetDphi","qqH_m2_mjj_jetDphi",15,400,2500,10,0.5,3.14);
  TH2F* qqH_m2_mjj_MT = new TH2F("qqH_m2_mjj_MT","qqH_m2_mjj_MT",15,400,2500,10,0,500);
  TH2F* qqH_m2_detajj_jetmetDphi = new TH2F("qqH_m2_detajj_jetmetDphi","qqH_m2_detajj_jetmetDphi",10,2,8,10,0.5,3.14);
  TH2F* qqH_m2_detajj_jetDphi = new TH2F("qqH_m2_detajj_jetDphi","qqH_m2_detajj_jetDphi",10,2,8,10,0.5,3.14);
  TH2F* qqH_m2_detajj_MT = new TH2F("qqH_m2_detajj_MT","qqH_m2_detajj_MT",10,2,8,10,0,500);

  qqH_m2_mjj_detajj->Sumw2();
  qqH_m2_mjj_jetmetDphi->Sumw2();
  qqH_m2_mjj_jetDphi->Sumw2();
  qqH_m2_mjj_MT->Sumw2();
  qqH_m2_detajj_jetmetDphi->Sumw2();
  qqH_m2_detajj_jetDphi->Sumw2();
  qqH_m2_detajj_MT->Sumw2();


  ///////////////

  TH1F* qqH_m3_bosonPt = new TH1F("qqH_m3_bosonPt","qqH_m3_bosonPt",25,150,1250);
  TH1F* qqH_m3_bosonEta = new TH1F("qqH_m3_bosonEta","qqH_m3_bosonEta",25,-3.5,3.5);
  TH1F* qqH_m3_bosonPhi = new TH1F("qqH_m3_bosonPhi","qqH_m3_bosonPhi",25,0.,3.14);
  TH1F* qqH_m3_bosonPz    = new TH1F("qqH_m3_bosonPz","qqH_m3_bosonPz",25,150,1250);
  TH1F* qqH_m3_bosonPzOPt = new TH1F("qqH_m3_bosonPzOPt","qqH_m3_bosonPzOPt",20,0,20);
  TH1F* qqH_m3_jetpt1   = new TH1F("qqH_m3_jetpt1","qqH_m3_jetpt1",25,50,500);
  TH1F* qqH_m3_jetpt2   = new TH1F("qqH_m3_jetpt2","qqH_m3_jetpt2",25,50,500);
  TH1F* qqH_m3_jetpz1   = new TH1F("qqH_m3_jetpz1","qqH_m3_jetpz1",25,50,500);  
  TH1F* qqH_m3_jetpz2   = new TH1F("qqH_m3_jetpz2","qqH_m3_jetpz2",25,50,500);
  TH1F* qqH_m3_jetpz1Opt   = new TH1F("qqH_m3_jetpz1Opt","qqH_m3_jetpz1",20,0,20);  
  TH1F* qqH_m3_jetpz2Opt   = new TH1F("qqH_m3_jetpz2OPt","qqH_m3_jetpz2",20,0,20);
  TH1F* qqH_m3_jeteta1  = new TH1F("qqH_m3_jeteta1","qqH_m3_jeteta1",25,-3.5,3.5);
  TH1F* qqH_m3_jeteta2  = new TH1F("qqH_m3_jeteta2","qqH_m3_jeteta2",25,-3.5,3.5);
  TH1F* qqH_m3_jetphi1  = new TH1F("qqH_m3_jetphi1","qqH_m3_jetphi1",25,0,3.14);
  TH1F* qqH_m3_jetphi2  = new TH1F("qqH_m3_jetphi2","qqH_m3_jetphi2",25,0,3.14);
  TH1F* qqH_m3_jetDeltaPhi  = new TH1F("qqH_m3_jetDeltaPhi","qqH_m3_jetDeltaPhi",25,0,3.14);
  TH1F* qqH_m3_jetDeltaEta  = new TH1F("qqH_m3_jetDeltaEta","qqH_m3_jetDeltaEta",25,0,9);
  TH1F* qqH_m3_Mjj     = new TH1F("qqH_m3_Mjj","qqH_m3_Mjj",40,400,2500);
  TH1F* qqH_m3_jetmetDphi4 = new TH1F("qqH_m3_jetmetdphi4","qqH_m3_jetmetdphi4",25,0,3.14);
  TH1F* qqH_m3_jetmetDphi  = new TH1F("qqH_m3_jetmetdphi","qqH_m3_jetmetdphi",25,0,3.14);
  TH1F* qqH_m3_jetmediatorDphi = new TH1F("qqH_m3_jetmediatorDphi","qqH_m3_jetmediatorDphi",25,0,3.14);
  TH1F* qqH_m3_jetmediatorDeta = new TH1F("qqH_m3_jetmediatorDeta","qqH_m3_jetmediatorDeta",25,0,9);
  TH1F* qqH_m3_MT = new TH1F("qqH_m3_MT","qqH_m3_MT",30,0,600);
  TH1F* qqH_m3_MetMT = new TH1F("qqH_m3_MetMT","qqH_m3_MetMT",30,100,1000);

  qqH_m3_bosonPt->Sumw2();
  qqH_m3_bosonEta->Sumw2();
  qqH_m3_bosonPhi->Sumw2();
  qqH_m3_bosonPz->Sumw2();
  qqH_m3_bosonPzOPt->Sumw2();
  qqH_m3_jetpt1->Sumw2();
  qqH_m3_jetpt2->Sumw2();
  qqH_m3_jetpz1->Sumw2();
  qqH_m3_jetpz2->Sumw2();
  qqH_m3_jetpz1Opt->Sumw2();
  qqH_m3_jetpz2Opt->Sumw2();
  qqH_m3_jeteta1->Sumw2();
  qqH_m3_jeteta2->Sumw2();
  qqH_m3_jetphi1->Sumw2();
  qqH_m3_jetphi2->Sumw2();
  qqH_m3_jetDeltaPhi->Sumw2();
  qqH_m3_jetDeltaEta->Sumw2();
  qqH_m3_Mjj->Sumw2();
  qqH_m3_jetmetDphi->Sumw2();
  qqH_m3_jetmetDphi4->Sumw2();
  qqH_m3_jetmediatorDphi->Sumw2();
  qqH_m3_jetmediatorDeta->Sumw2();
  qqH_m3_MT->Sumw2();
  qqH_m3_MetMT->Sumw2();

  TH2F* qqH_m3_mjj_detajj = new TH2F("qqH_m3_mjj_detajj","qqH_m3_mjj_detajj",15,400,2500,10,2,8);
  TH2F* qqH_m3_mjj_jetmetDphi = new TH2F("qqH_m3_mjj_jetmetDphi","qqH_m3_mjj_jetmetDphi",15,400,2500,10,0.5,3.14);
  TH2F* qqH_m3_mjj_jetDphi = new TH2F("qqH_m3_mjj_jetDphi","qqH_m3_mjj_jetDphi",15,400,2500,10,0.5,3.14);
  TH2F* qqH_m3_mjj_MT = new TH2F("qqH_m3_mjj_MT","qqH_m3_mjj_MT",15,400,2500,10,0,500);
  TH2F* qqH_m3_detajj_jetmetDphi = new TH2F("qqH_m3_detajj_jetmetDphi","qqH_m3_detajj_jetmetDphi",10,2,8,10,0.5,3.14);
  TH2F* qqH_m3_detajj_jetDphi = new TH2F("qqH_m3_detajj_jetDphi","qqH_m3_detajj_jetDphi",10,2,8,10,0.5,3.14);
  TH2F* qqH_m3_detajj_MT = new TH2F("qqH_m3_detajj_MT","qqH_m3_detajj_MT",10,2,8,10,0,500);

  qqH_m3_mjj_detajj->Sumw2();
  qqH_m3_mjj_jetmetDphi->Sumw2();
  qqH_m3_mjj_jetDphi->Sumw2();
  qqH_m3_mjj_MT->Sumw2();
  qqH_m3_detajj_jetmetDphi->Sumw2();
  qqH_m3_detajj_jetDphi->Sumw2();
  qqH_m3_detajj_MT->Sumw2();
  

  ////////////

  while(myReader.Next()){
    // met filters
    if (*fhbhe  == 0 or *fhbiso == 0 or *feeb == 0 or *fcsc  == 0) continue;
    // number of jets
    if (*njetsinc  < 2) continue;
    // b-veto
    if (*nbjets > 0) continue;
    // vetos
    if (*nmuons > 0)     continue;
    if (*nelectrons > 0) continue;
    if (*ntaus > 0)      continue;
    if (*nphotons > 0)   continue;
    if (*hltm90 == 0 and *hltm120 == 0 and *hltmwm120 == 0 and *hltmwm170 == 0 and *hltmwm300 == 0 and *hltmwm90 == 0 ) continue;
    // relaxed vbf jet pt cuts
    if (jetpt->at(0) < 50) continue; 
    if (jetpt->at(1) < 50) continue; 
    if (fabs(jeteta->at(0)) > 4.7) continue;
    if (fabs(jeteta->at(1)) > 4.7) continue;    
    if (fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) < 0.1) continue;
    if (fabs(jeteta->at(0)) < 2.5 and nhfrac->at(0) > 0.8) continue;
    // relaxed met cut
    if (*mmet < 150) continue;
    if (*jmmdphi4 < 0.5) continue;

    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
    // relaxed VBF cuts
    if (fabs(jeteta->at(0)-jeteta->at(1)) < 2) continue;
    if ((jet1+jet2).M() < 450) continue;
    
    TLorentzVector mediator;
    mediator.SetPtEtaPhiM(*mediatorPt,*mediatorEta,*mediatorPhi,*mediatorMass);
    TLorentzVector met;
    met.SetPxPyPzE(*mmet*TMath::Cos(*mmetphi),*mmet*TMath::Sin(*mmetphi),0.,*mmet);
    // fill some hisotgrams
    fillHisto(qqH_m1_bosonPt,*mediatorPt);
    fillHisto(qqH_m1_bosonEta,*mediatorEta);
    fillHisto(qqH_m1_bosonPhi,fabs(*mediatorPhi));
    fillHisto(qqH_m1_bosonPz,mediator.Pz());
    fillHisto(qqH_m1_bosonPzOPt,mediator.Pz()/mediator.Pt());
    fillHisto(qqH_m1_jetpt1,jetpt->at(0));
    fillHisto(qqH_m1_jetpt2,jetpt->at(1));
    fillHisto(qqH_m1_jetpz1,jet1.Pz());
    fillHisto(qqH_m1_jetpz2,jet2.Pz());
    fillHisto(qqH_m1_jetpz1Opt,jet1.Pz()/jet1.Pt());
    fillHisto(qqH_m1_jetpz2Opt,jet1.Pz()/jet1.Pt());
    fillHisto(qqH_m1_jeteta1,jeteta->at(0));
    fillHisto(qqH_m1_jeteta2,jeteta->at(1));
    fillHisto(qqH_m1_jetphi1,fabs(jetphi->at(0)));
    fillHisto(qqH_m1_jetphi2,fabs(jetphi->at(1)));
    fillHisto(qqH_m1_jetDeltaPhi,fabs(jet1.DeltaPhi(jet2)));
    fillHisto(qqH_m1_jetDeltaEta,fabs(jet1.Eta()-jet2.Eta()));
    fillHisto(qqH_m1_Mjj,(jet1+jet2).M());
    fillHisto(qqH_m1_jetmetDphi,*jmmdphi);
    fillHisto(qqH_m1_jetmetDphi4,*jmmdphi4);
    fillHisto(qqH_m1_jetmediatorDphi,fabs(mediator.DeltaPhi(jet1+jet2)));
    fillHisto(qqH_m1_jetmediatorDeta,fabs(mediator.Eta()-(jet1+jet2).Eta()));
    fillHisto(qqH_m1_MT,sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));    
    fillHisto(qqH_m1_MetMT,sqrt(2*(jet1+jet2).Pt()*met.Pt()*(1-cos((jet1+jet2).DeltaPhi(met)))));

    // 2D part
    fillHisto(qqH_m1_mjj_detajj,(jet1+jet2).M(),fabs(jet1.Eta()-jet2.Eta()));
    fillHisto(qqH_m1_mjj_jetmetDphi,(jet1+jet2).M(),*jmmdphi);
    fillHisto(qqH_m1_mjj_jetDphi,(jet1+jet2).M(),fabs(jet1.DeltaPhi(jet2)));
    fillHisto(qqH_m1_mjj_MT,(jet1+jet2).M(),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));    
    fillHisto(qqH_m1_detajj_jetmetDphi,fabs(jet1.Eta()-jet2.Eta()),*jmmdphi);
    fillHisto(qqH_m1_detajj_jetDphi,fabs(jet1.Eta()-jet2.Eta()),fabs(jet1.DeltaPhi(jet2)));
    fillHisto(qqH_m1_detajj_MT,fabs(jet1.Eta()-jet2.Eta()),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));

  }



  ////////////
  myReader.SetTree(vbftree_m2);
  myReader.SetEntry(0);
  
  while(myReader.Next()){
    // met filters
    if (*fhbhe  == 0 or *fhbiso == 0 or *feeb == 0 or *fcsc  == 0) continue;
    // number of jets
    if (*njetsinc  < 2) continue;
    // b-veto
    if (*nbjets > 0) continue;
    // vetos
    if (*nmuons > 0)     continue;
    if (*nelectrons > 0) continue;
    if (*ntaus > 0)      continue;
    if (*nphotons > 0)   continue;
    if (*hltm90 == 0 and *hltm120 == 0 and *hltmwm120 == 0 and *hltmwm170 == 0 and *hltmwm300 == 0 and *hltmwm90 == 0 ) continue;
    // relaxed vbf jet pt cuts
    if (jetpt->at(0) < 50) continue; 
    if (jetpt->at(1) < 50) continue; 
    if (fabs(jeteta->at(0)) > 4.7) continue;
    if (fabs(jeteta->at(1)) > 4.7) continue;    
    if (fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) < 0.1) continue;
    if (fabs(jeteta->at(0)) < 2.5 and nhfrac->at(0) > 0.8) continue;
    // relaxed met cut
    if (*mmet < 150) continue;
    if (*jmmdphi4 < 0.5) continue;

    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
    // relaxed VBF cuts
    if (fabs(jeteta->at(0)-jeteta->at(1)) < 2) continue;
    if ((jet1+jet2).M() < 450) continue;
    
    TLorentzVector mediator;
    mediator.SetPtEtaPhiM(*mediatorPt,*mediatorEta,*mediatorPhi,*mediatorMass);
    TLorentzVector met;
    met.SetPxPyPzE(*mmet*TMath::Cos(*mmetphi),*mmet*TMath::Sin(*mmetphi),0.,*mmet);
    // fill some hisotgrams
    fillHisto(qqH_m2_bosonPt,*mediatorPt);
    fillHisto(qqH_m2_bosonEta,*mediatorEta);
    fillHisto(qqH_m2_bosonPhi,fabs(*mediatorPhi));
    fillHisto(qqH_m2_bosonPz,mediator.Pz());
    fillHisto(qqH_m2_bosonPzOPt,mediator.Pz()/mediator.Pt());
    fillHisto(qqH_m2_jetpt1,jetpt->at(0));
    fillHisto(qqH_m2_jetpt2,jetpt->at(1));
    fillHisto(qqH_m2_jetpz1,jet1.Pz());
    fillHisto(qqH_m2_jetpz2,jet2.Pz());
    fillHisto(qqH_m2_jetpz1Opt,jet1.Pz()/jet1.Pt());
    fillHisto(qqH_m2_jetpz2Opt,jet1.Pz()/jet1.Pt());
    fillHisto(qqH_m2_jeteta1,jeteta->at(0));
    fillHisto(qqH_m2_jeteta2,jeteta->at(1));
    fillHisto(qqH_m2_jetphi1,fabs(jetphi->at(0)));
    fillHisto(qqH_m2_jetphi2,fabs(jetphi->at(1)));
    fillHisto(qqH_m2_jetDeltaPhi,fabs(jet1.DeltaPhi(jet2)));
    fillHisto(qqH_m2_jetDeltaEta,fabs(jet1.Eta()-jet2.Eta()));
    fillHisto(qqH_m2_Mjj,(jet1+jet2).M());
    fillHisto(qqH_m2_jetmetDphi,*jmmdphi);
    fillHisto(qqH_m2_jetmetDphi4,*jmmdphi4);
    fillHisto(qqH_m2_jetmediatorDphi,fabs(mediator.DeltaPhi(jet1+jet2)));
    fillHisto(qqH_m2_jetmediatorDeta,fabs(mediator.Eta()-(jet1+jet2).Eta()));
    fillHisto(qqH_m2_MT,sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));    
    fillHisto(qqH_m2_MetMT,sqrt(2*(jet1+jet2).Pt()*met.Pt()*(1-cos((jet1+jet2).DeltaPhi(met)))));

    // 2D part
    fillHisto(qqH_m2_mjj_detajj,(jet1+jet2).M(),fabs(jet1.Eta()-jet2.Eta()));
    fillHisto(qqH_m2_mjj_jetmetDphi,(jet1+jet2).M(),*jmmdphi);
    fillHisto(qqH_m2_mjj_jetDphi,(jet1+jet2).M(),fabs(jet1.DeltaPhi(jet2)));
    fillHisto(qqH_m2_mjj_MT,(jet1+jet2).M(),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));    
    fillHisto(qqH_m2_detajj_jetmetDphi,fabs(jet1.Eta()-jet2.Eta()),*jmmdphi);
    fillHisto(qqH_m2_detajj_jetDphi,fabs(jet1.Eta()-jet2.Eta()),fabs(jet1.DeltaPhi(jet2)));
    fillHisto(qqH_m2_detajj_MT,fabs(jet1.Eta()-jet2.Eta()),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));

  }




  ////////////
  myReader.SetTree(vbftree_m3);
  myReader.SetEntry(0);
  
  while(myReader.Next()){
    // met filters
    if (*fhbhe  == 0 or *fhbiso == 0 or *feeb == 0 or *fcsc  == 0) continue;
    // number of jets
    if (*njetsinc  < 2) continue;
    // b-veto
    if (*nbjets > 0) continue;
    // vetos
    if (*nmuons > 0)     continue;
    if (*nelectrons > 0) continue;
    if (*ntaus > 0)      continue;
    if (*nphotons > 0)   continue;
    if (*hltm90 == 0 and *hltm120 == 0 and *hltmwm120 == 0 and *hltmwm170 == 0 and *hltmwm300 == 0 and *hltmwm90 == 0 ) continue;
    // relaxed vbf jet pt cuts
    if (jetpt->at(0) < 50) continue; 
    if (jetpt->at(1) < 50) continue; 
    if (fabs(jeteta->at(0)) > 4.7) continue;
    if (fabs(jeteta->at(1)) > 4.7) continue;    
    if (fabs(jeteta->at(0)) < 2.5 and chfrac->at(0) < 0.1) continue;
    if (fabs(jeteta->at(0)) < 2.5 and nhfrac->at(0) > 0.8) continue;
    // relaxed met cut
    if (*mmet < 150) continue;
    if (*jmmdphi4 < 0.5) continue;

    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
    // relaxed VBF cuts
    if (fabs(jeteta->at(0)-jeteta->at(1)) < 2) continue;
    if ((jet1+jet2).M() < 450) continue;
    
    TLorentzVector mediator;
    mediator.SetPtEtaPhiM(*mediatorPt,*mediatorEta,*mediatorPhi,*mediatorMass);
    TLorentzVector met;
    met.SetPxPyPzE(*mmet*TMath::Cos(*mmetphi),*mmet*TMath::Sin(*mmetphi),0.,*mmet);
    // fill some hisotgrams
    fillHisto(qqH_m3_bosonPt,*mediatorPt);
    fillHisto(qqH_m3_bosonEta,*mediatorEta);
    fillHisto(qqH_m3_bosonPhi,fabs(*mediatorPhi));
    fillHisto(qqH_m3_bosonPz,mediator.Pz());
    fillHisto(qqH_m3_bosonPzOPt,mediator.Pz()/mediator.Pt());
    fillHisto(qqH_m3_jetpt1,jetpt->at(0));
    fillHisto(qqH_m3_jetpt2,jetpt->at(1));
    fillHisto(qqH_m3_jetpz1,jet1.Pz());
    fillHisto(qqH_m3_jetpz2,jet2.Pz());
    fillHisto(qqH_m3_jetpz1Opt,jet1.Pz()/jet1.Pt());
    fillHisto(qqH_m3_jetpz2Opt,jet1.Pz()/jet1.Pt());
    fillHisto(qqH_m3_jeteta1,jeteta->at(0));
    fillHisto(qqH_m3_jeteta2,jeteta->at(1));
    fillHisto(qqH_m3_jetphi1,fabs(jetphi->at(0)));
    fillHisto(qqH_m3_jetphi2,fabs(jetphi->at(1)));
    fillHisto(qqH_m3_jetDeltaPhi,fabs(jet1.DeltaPhi(jet2)));
    fillHisto(qqH_m3_jetDeltaEta,fabs(jet1.Eta()-jet2.Eta()));
    fillHisto(qqH_m3_Mjj,(jet1+jet2).M());
    fillHisto(qqH_m3_jetmetDphi,*jmmdphi);
    fillHisto(qqH_m3_jetmetDphi4,*jmmdphi4);
    fillHisto(qqH_m3_jetmediatorDphi,fabs(mediator.DeltaPhi(jet1+jet2)));
    fillHisto(qqH_m3_jetmediatorDeta,fabs(mediator.Eta()-(jet1+jet2).Eta()));
    fillHisto(qqH_m3_MT,sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));    
    fillHisto(qqH_m3_MetMT,sqrt(2*(jet1+jet2).Pt()*met.Pt()*(1-cos((jet1+jet2).DeltaPhi(met)))));

    // 2D part
    fillHisto(qqH_m3_mjj_detajj,(jet1+jet2).M(),fabs(jet1.Eta()-jet2.Eta()));
    fillHisto(qqH_m3_mjj_jetmetDphi,(jet1+jet2).M(),*jmmdphi);
    fillHisto(qqH_m3_mjj_jetDphi,(jet1+jet2).M(),fabs(jet1.DeltaPhi(jet2)));
    fillHisto(qqH_m3_mjj_MT,(jet1+jet2).M(),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));    
    fillHisto(qqH_m3_detajj_jetmetDphi,fabs(jet1.Eta()-jet2.Eta()),*jmmdphi);
    fillHisto(qqH_m3_detajj_jetDphi,fabs(jet1.Eta()-jet2.Eta()),fabs(jet1.DeltaPhi(jet2)));
    fillHisto(qqH_m3_detajj_MT,fabs(jet1.Eta()-jet2.Eta()),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));

  }

  
  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();
  /////////////
  plotDistribution(canvas,qqH_m1_bosonPt,qqH_m2_bosonPt,qqH_m3_bosonPt,"Mediator p_{T} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_m1_bosonEta,qqH_m2_bosonEta,qqH_m3_bosonEta,"Mediator #eta",outputPlot);
  plotDistribution(canvas,qqH_m1_bosonPhi,qqH_m2_bosonPhi,qqH_m3_bosonPhi,"Mediator #phi",outputPlot);
  plotDistribution(canvas,qqH_m1_bosonPz,qqH_m2_bosonPz,qqH_m3_bosonPz,"Mediator p_{z} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_m1_bosonPzOPt,qqH_m2_bosonPzOPt,qqH_m3_bosonPzOPt,"Mediator p_{z}/p_{T}",outputPlot);
  plotDistribution(canvas,qqH_m1_jetpt1,qqH_m2_jetpt1,qqH_m3_jetpt1,"p^{j1}_{T} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_m1_jetpt2,qqH_m2_jetpt2,qqH_m3_jetpt2,"p^{j2}_{T} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_m1_jetpz1,qqH_m2_jetpz1,qqH_m3_jetpz1,"p^{j1}_{z} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_m1_jetpz2,qqH_m2_jetpz2,qqH_m3_jetpz2,"p^{j2}_{z} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_m1_jetpz1Opt,qqH_m2_jetpz1Opt,qqH_m3_jetpz1Opt,"p^{j1}_{z}/p^{j1}_{T}",outputPlot);
  plotDistribution(canvas,qqH_m1_jetpz2Opt,qqH_m2_jetpz2Opt,qqH_m3_jetpz2Opt,"p^{j1}_{z}/p^{j1}_{T}",outputPlot);
  plotDistribution(canvas,qqH_m1_jeteta1,qqH_m2_jeteta1,qqH_m3_jeteta1,"#eta_{j1}",outputPlot);
  plotDistribution(canvas,qqH_m1_jeteta2,qqH_m2_jeteta2,qqH_m3_jeteta2,"#eta_{j2}",outputPlot);
  plotDistribution(canvas,qqH_m1_jetphi1,qqH_m2_jetphi1,qqH_m3_jetphi1,"#phi_{j1}",outputPlot);
  plotDistribution(canvas,qqH_m1_jetphi2,qqH_m2_jetphi2,qqH_m3_jetphi2,"#phi_{j2}",outputPlot);
  plotDistribution(canvas,qqH_m1_jetDeltaPhi,qqH_m2_jetDeltaPhi,qqH_m3_jetDeltaPhi,"#Delta#phi_{j1,j2}",outputPlot);
  plotDistribution(canvas,qqH_m1_jetDeltaEta,qqH_m2_jetDeltaEta,qqH_m3_jetDeltaEta,"#Delta#eta_{j1,j2}",outputPlot);
  plotDistribution(canvas,qqH_m1_Mjj,qqH_m2_Mjj,qqH_m3_Mjj,"M_{j1,j2}",outputPlot);
  plotDistribution(canvas,qqH_m1_jetmetDphi,qqH_m2_jetmetDphi,qqH_m3_jetmetDphi,"#Delta#phi_{jet,met}",outputPlot);
  plotDistribution(canvas,qqH_m1_jetmetDphi4,qqH_m2_jetmetDphi4,qqH_m3_jetmetDphi4,"#Delta#phi_{4-jet,met}",outputPlot);
  plotDistribution(canvas,qqH_m1_jetmediatorDeta,qqH_m2_jetmediatorDeta,qqH_m3_jetmediatorDeta,"#Delta#eta_{jet,med}",outputPlot);
  plotDistribution(canvas,qqH_m1_jetmediatorDphi,qqH_m2_jetmediatorDphi,qqH_m3_jetmediatorDphi,"#Delta#phi_{jet,med}",outputPlot);
  plotDistribution(canvas,qqH_m1_MetMT,qqH_m2_MetMT,qqH_m3_MetMT,"m_{T}(jj,met) [GeV]",outputPlot);
  plotDistribution(canvas,qqH_m1_MT,qqH_m2_MT,qqH_m3_MT,"m_{T}(jj) [GeV]",outputPlot);
  /////////////
  plotDistribution(canvas,qqH_m1_mjj_detajj,qqH_m2_mjj_detajj,qqH_m3_mjj_detajj,"m_{jj} [GeV]","#Delta#eta_{jj}",outputPlot);
  plotDistribution(canvas,qqH_m1_mjj_jetmetDphi,qqH_m2_mjj_jetmetDphi,qqH_m3_mjj_jetmetDphi,"m_{jj} [GeV]","#Delta#phi_{jet,met}",outputPlot);
  plotDistribution(canvas,qqH_m1_mjj_jetDphi,qqH_m2_mjj_jetDphi,qqH_m3_mjj_jetDphi,"m_{jj} [GeV]","#Delta#phi_{jj}",outputPlot);
  plotDistribution(canvas,qqH_m1_mjj_MT,qqH_m2_mjj_MT,qqH_m3_mjj_MT,"m_{jj} [GeV]","m_{T}(jj) [GeV]",outputPlot);
  plotDistribution(canvas,qqH_m1_detajj_jetmetDphi,qqH_m2_detajj_jetmetDphi,qqH_m3_detajj_jetmetDphi,"#Delta#eta_{jj}","#Delta#phi_{jet,met}",outputPlot);
  plotDistribution(canvas,qqH_m1_detajj_jetDphi,qqH_m2_detajj_jetDphi,qqH_m3_detajj_jetDphi,"#Delta#eta_{jj}","#Delta#phi_{jj}",outputPlot);
  plotDistribution(canvas,qqH_m1_detajj_MT,qqH_m2_detajj_MT,qqH_m3_detajj_MT,"#Delta#eta_{jj}","m_{T}(jj) [GeV]",outputPlot);
  
}
