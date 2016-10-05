#include "../CMS_lumi.h"

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


void plotDistribution(TCanvas* canvas, TH1* qqH, TH1* ggH, const string & label, const string & outputDIR){
  
  // normalize to 1
  canvas->cd();
  qqH->Scale(1./qqH->Integral());
  ggH->Scale(1./ggH->Integral());
  qqH->GetXaxis()->SetTitle(label.c_str());
  qqH->GetYaxis()->SetTitle("a.u");
  qqH->SetLineColor(kBlack);
  qqH->SetLineWidth(2);
  qqH->Draw("hist");
  
  ggH->SetLineColor(kRed);
  ggH->SetLineWidth(2);
  ggH->Draw("hist same");

  qqH->GetYaxis()->SetRangeUser(0.,max(qqH->GetMaximum(),ggH->GetMaximum())*1.2);
  CMS_lumi(canvas,"");
  TLegend leg(0.7,0.7,0.9,0.9);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.AddEntry(qqH,"VBF H = 125 GeV","L");
  leg.AddEntry(ggH,"ggH H = 125 GeV","L");
  leg.Draw("same");
  
  //  canvas->SaveAs((outputDIR+"/"+string(qqH->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+string(qqH->GetName())+".pdf").c_str(),"pdf");

  qqH->GetYaxis()->SetRangeUser(max(0.0001,min(qqH->GetMinimum(),ggH->GetMinimum())*0.8),max(qqH->GetMaximum(),ggH->GetMaximum())*1000);
  canvas->SetLogy();
  //  canvas->SaveAs((outputDIR+"/"+string(qqH->GetName())+"_log.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+string(qqH->GetName())+"_log.pdf").c_str(),"pdf");
  canvas->SetLogy(0);
  
}



void plotDistribution(TCanvas* canvas, TH2* qqH, TH2* ggH, const string & labelX, const string & labelY, const string & outputDIR){
  
  system(("mkdir -p "+outputDIR+"/correlation/").c_str());
  // normalize to 1
  canvas->cd();
  canvas->SetRightMargin(0.18);
  qqH->Scale(1./qqH->Integral());
  ggH->Scale(1./ggH->Integral());

  TGraph2D* qqHGraph = new TGraph2D();
  qqHGraph->SetNpx(300);
  qqHGraph->SetNpy(300);
  int nPoint = 0;
  for(int iBinX = 0; iBinX < qqH->GetNbinsX() ; iBinX++){
    for(int iBinY = 0; iBinY < qqH->GetNbinsY() ; iBinY++){
      qqHGraph->SetPoint(nPoint,qqH->GetXaxis()->GetBinCenter(iBinX+1),qqH->GetYaxis()->GetBinCenter(iBinY+1),qqH->GetBinContent(iBinX+1,iBinY+1));      
      nPoint++;
    }
  }

  TH2* qqHPlot = qqHGraph->GetHistogram();
  qqHPlot->GetXaxis()->SetTitle(labelX.c_str());
  qqHPlot->GetYaxis()->SetTitle(labelY.c_str());
  qqHPlot->GetZaxis()->SetTitle("a.u");   
  qqHPlot->Draw("colz");

  TProfile* qqHProfile = qqH->ProfileX(Form("%s_pfx",qqH->GetName()));
  qqHProfile->SetMarkerColor(kBlack);
  qqHProfile->SetMarkerStyle(20);
  qqHProfile->SetMarkerSize(1);
  qqHProfile->Draw("EPsame");

  CMS_lumi(canvas,"");
  TLegend leg(0.6,0.6,0.9,0.9);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.SetBorderSize(0);
  leg.AddEntry((TObject*)0,"VBF H = 125 GeV","");
  leg.AddEntry((TObject*)0,Form("Correlation = %.2f",qqHPlot->GetCorrelationFactor()),"");
  leg.Draw("same");

  //  canvas->SaveAs((outputDIR+"/correlation/"+string(qqH->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/correlation/"+string(qqH->GetName())+".pdf").c_str(),"pdf");
  canvas->SetLogz();
  //  canvas->SaveAs((outputDIR+"/correlation/"+string(qqH->GetName())+"_log.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/correlation/"+string(qqH->GetName())+"_log.pdf").c_str(),"pdf");
  canvas->SetLogz(0);

  TGraph2D* ggHGraph = new TGraph2D();
  ggHGraph->SetNpx(300);
  ggHGraph->SetNpy(300);
  nPoint = 0;
  for(int iBinX = 0; iBinX < ggH->GetNbinsX() ; iBinX++){
    for(int iBinY = 0; iBinY < ggH->GetNbinsY() ; iBinY++){
      ggHGraph->SetPoint(nPoint,ggH->GetXaxis()->GetBinCenter(iBinX+1),ggH->GetYaxis()->GetBinCenter(iBinY+1),ggH->GetBinContent(iBinX+1,iBinY+1));      
      nPoint++;
    }
  }

  TH2* ggHPlot = ggHGraph->GetHistogram();
  ggHPlot->GetXaxis()->SetTitle(labelX.c_str());
  ggHPlot->GetYaxis()->SetTitle(labelY.c_str());
  ggHPlot->GetZaxis()->SetTitle("a.u");   
  ggHPlot->Draw("colz");
  CMS_lumi(canvas,"");

  TProfile* ggHProfile = ggH->ProfileX(Form("%s_pfx",ggH->GetName()));
  ggHProfile->SetMarkerColor(kBlack);
  ggHProfile->SetMarkerStyle(20);
  ggHProfile->SetMarkerSize(1);
  ggHProfile->Draw("EPsame");
 
  leg.Clear();
  leg.AddEntry((TObject*)0,"ggH H = 125 GeV","");
  leg.AddEntry((TObject*)0,Form("Correlation = %.2f",ggHPlot->GetCorrelationFactor()),"");
  leg.Draw();

  //  canvas->SaveAs((outputDIR+"/correlation/"+string(ggH->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/correlation/"+string(ggH->GetName())+".pdf").c_str(),"pdf");
  canvas->SetLogz();
  //  canvas->SaveAs((outputDIR+"/correlation/"+string(ggH->GetName())+"_log.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/correlation/"+string(ggH->GetName())+"_log.pdf").c_str(),"pdf");
  canvas->SetLogz(0);
  
}

void makeVBFSignalStudies(string outputPlot){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputPlot).c_str());

  TChain* vbftree = new TChain("tree/tree");
  TChain* ggHtree = new TChain("tree/tree");
  vbftree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/*VBF*125*root");
  vbftree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/*VBF*150*root");
  vbftree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/*VBF*200*root");
  ggHtree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/*GluGlu*125*root");
  ggHtree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/*GluGlu*150*root");
  ggHtree->Add("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/MetCut/Production_30_09_2016/HiggsInvisible/sigfilter/*GluGlu*200*root");
  
  
  TTreeReader myReader (vbftree);
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
  
  TH1F* qqH_bosonPt = new TH1F("qqH_bosonPt","qqH_bosonPt",25,150,1250);
  TH1F* qqH_bosonEta = new TH1F("qqH_bosonEta","qqH_bosonEta",25,-3.5,3.5);
  TH1F* qqH_bosonPhi = new TH1F("qqH_bosonPhi","qqH_bosonPhi",25,0.,3.14);
  TH1F* qqH_bosonPz    = new TH1F("qqH_bosonPz","qqH_bosonPz",25,150,1250);
  TH1F* qqH_bosonPzOPt = new TH1F("qqH_bosonPzOPt","qqH_bosonPzOPt",20,0,20);
  TH1F* qqH_jetpt1   = new TH1F("qqH_jetpt1","qqH_jetpt1",25,50,500);
  TH1F* qqH_jetpt2   = new TH1F("qqH_jetpt2","qqH_jetpt2",25,50,500);
  TH1F* qqH_jetpz1   = new TH1F("qqH_jetpz1","qqH_jetpz1",25,50,500);  
  TH1F* qqH_jetpz2   = new TH1F("qqH_jetpz2","qqH_jetpz2",25,50,500);
  TH1F* qqH_jetpz1Opt   = new TH1F("qqH_jetpz1Opt","qqH_jetpz1",20,0,20);  
  TH1F* qqH_jetpz2Opt   = new TH1F("qqH_jetpz2OPt","qqH_jetpz2",20,0,20);
  TH1F* qqH_jeteta1  = new TH1F("qqH_jeteta1","qqH_jeteta1",25,-3.5,3.5);
  TH1F* qqH_jeteta2  = new TH1F("qqH_jeteta2","qqH_jeteta2",25,-3.5,3.5);
  TH1F* qqH_jetphi1  = new TH1F("qqH_jetphi1","qqH_jetphi1",25,0,3.14);
  TH1F* qqH_jetphi2  = new TH1F("qqH_jetphi2","qqH_jetphi2",25,0,3.14);
  TH1F* qqH_jetDeltaPhi  = new TH1F("qqH_jetDeltaPhi","qqH_jetDeltaPhi",25,0,3.14);
  TH1F* qqH_jetDeltaEta  = new TH1F("qqH_jetDeltaEta","qqH_jetDeltaEta",25,0,9);
  TH1F* qqH_Mjj     = new TH1F("qqH_Mjj","qqH_Mjj",40,400,2500);
  TH1F* qqH_jetmetDphi4 = new TH1F("qqH_jetmetdphi4","qqH_jetmetdphi4",25,0,3.14);
  TH1F* qqH_jetmetDphi  = new TH1F("qqH_jetmetdphi","qqH_jetmetdphi",25,0,3.14);
  TH1F* qqH_jetmediatorDphi = new TH1F("qqH_jetmediatorDphi","qqH_jetmediatorDphi",25,0,3.14);
  TH1F* qqH_jetmediatorDeta = new TH1F("qqH_jetmediatorDeta","qqH_jetmediatorDeta",25,0,9);
  TH1F* qqH_MT = new TH1F("qqH_MT","qqH_MT",30,0,600);
  TH1F* qqH_MetMT = new TH1F("qqH_MetMT","qqH_MetMT",30,100,1000);

  qqH_bosonPt->Sumw2();
  qqH_bosonEta->Sumw2();
  qqH_bosonPhi->Sumw2();
  qqH_bosonPz->Sumw2();
  qqH_bosonPzOPt->Sumw2();
  qqH_jetpt1->Sumw2();
  qqH_jetpt2->Sumw2();
  qqH_jetpz1->Sumw2();
  qqH_jetpz2->Sumw2();
  qqH_jetpz1Opt->Sumw2();
  qqH_jetpz2Opt->Sumw2();
  qqH_jeteta1->Sumw2();
  qqH_jeteta2->Sumw2();
  qqH_jetphi1->Sumw2();
  qqH_jetphi2->Sumw2();
  qqH_jetDeltaPhi->Sumw2();
  qqH_jetDeltaEta->Sumw2();
  qqH_Mjj->Sumw2();
  qqH_jetmetDphi->Sumw2();
  qqH_jetmetDphi4->Sumw2();
  qqH_jetmediatorDphi->Sumw2();
  qqH_jetmediatorDeta->Sumw2();
  qqH_MT->Sumw2();
  qqH_MetMT->Sumw2();

  TH2F* qqH_mjj_detajj = new TH2F("qqH_mjj_detajj","qqH_mjj_detajj",15,400,2500,10,2,8);
  TH2F* qqH_mjj_jetmetDphi = new TH2F("qqH_mjj_jetmetDphi","qqH_mjj_jetmetDphi",15,400,2500,10,0.5,3.14);
  TH2F* qqH_mjj_jetDphi = new TH2F("qqH_mjj_jetDphi","qqH_mjj_jetDphi",15,400,2500,10,0.5,3.14);
  TH2F* qqH_mjj_MT = new TH2F("qqH_mjj_MT","qqH_mjj_MT",15,400,2500,10,0,500);
  TH2F* qqH_detajj_jetmetDphi = new TH2F("qqH_detajj_jetmetDphi","qqH_detajj_jetmetDphi",10,2,8,10,0.5,3.14);
  TH2F* qqH_detajj_jetDphi = new TH2F("qqH_detajj_jetDphi","qqH_detajj_jetDphi",10,2,8,10,0.5,3.14);
  TH2F* qqH_detajj_MT = new TH2F("qqH_detajj_MT","qqH_detajj_MT",10,2,8,10,0,500);

  qqH_mjj_detajj->Sumw2();
  qqH_mjj_jetmetDphi->Sumw2();
  qqH_mjj_jetDphi->Sumw2();
  qqH_mjj_MT->Sumw2();
  qqH_detajj_jetmetDphi->Sumw2();
  qqH_detajj_jetDphi->Sumw2();
  qqH_detajj_MT->Sumw2();
  
  TH1F* ggH_bosonPt = new TH1F("ggH_bosonPt","ggH_bosonPt",25,150,1250);
  TH1F* ggH_bosonEta = new TH1F("ggH_bosonEta","ggH_bosonEta",25,-3.5,3.5);
  TH1F* ggH_bosonPhi = new TH1F("ggH_bosonPhi","ggH_bosonPhi",25,0.,3.14);
  TH1F* ggH_bosonPz    = new TH1F("ggH_bosonPz","ggH_bosonPz",25,150,1250);
  TH1F* ggH_bosonPzOPt = new TH1F("ggH_bosonPzOPt","ggH_bosonPzOPt",20,0,20);
  TH1F* ggH_jetpt1   = new TH1F("ggH_jetpt1","ggH_jetpt1",25,50,500);
  TH1F* ggH_jetpt2   = new TH1F("ggH_jetpt2","ggH_jetpt2",25,50,500);
  TH1F* ggH_jetpz1   = new TH1F("ggH_jetpz1","ggH_jetpz1",25,50,500);  
  TH1F* ggH_jetpz2   = new TH1F("ggH_jetpz2","ggH_jetpz2",25,50,500);
  TH1F* ggH_jetpz1Opt   = new TH1F("ggH_jetpz1Opt","ggH_jetpz1",20,0,20);  
  TH1F* ggH_jetpz2Opt   = new TH1F("ggH_jetpz2OPt","ggH_jetpz2",20,0,20);
  TH1F* ggH_jeteta1  = new TH1F("ggH_jeteta1","ggH_jeteta1",25,-3.5,3.5);
  TH1F* ggH_jeteta2  = new TH1F("ggH_jeteta2","ggH_jeteta2",25,-3.5,3.5);
  TH1F* ggH_jetphi1  = new TH1F("ggH_jetphi1","ggH_jetphi1",25,0,3.14);
  TH1F* ggH_jetphi2  = new TH1F("ggH_jetphi2","ggH_jetphi2",25,0,3.14);
  TH1F* ggH_jetDeltaPhi  = new TH1F("ggH_jetDeltaPhi","ggH_jetDeltaPhi",25,0,3.14);
  TH1F* ggH_jetDeltaEta  = new TH1F("ggH_jetDeltaEta","ggH_jetDeltaEta",25,0,9);
  TH1F* ggH_Mjj     = new TH1F("ggH_Mjj","ggH_Mjj",40,400,2500);
  TH1F* ggH_jetmetDphi4 = new TH1F("ggH_jetmetdphi4","ggH_jetmetdphi4",25,0,3.14);
  TH1F* ggH_jetmetDphi  = new TH1F("ggH_jetmetdphi","ggH_jetmetdphi",25,0,3.14);
  TH1F* ggH_jetmediatorDphi = new TH1F("ggH_jetmediatorDphi","ggH_jetmediatorDphi",25,0,3.14);
  TH1F* ggH_jetmediatorDeta = new TH1F("ggH_jetmediatorDeta","ggH_jetmediatorDeta",25,0,9);
  TH1F* ggH_MT = new TH1F("ggH_MT","ggH_MT",30,0,600);
  TH1F* ggH_MetMT = new TH1F("ggH_MetMT","ggH_MetMT",30,100,1000);

  ggH_bosonPt->Sumw2();
  ggH_bosonEta->Sumw2();
  ggH_bosonPhi->Sumw2();
  ggH_bosonPz->Sumw2();
  ggH_bosonPzOPt->Sumw2();
  ggH_jetpt1->Sumw2();
  ggH_jetpt2->Sumw2();
  ggH_jetpz1->Sumw2();
  ggH_jetpz2->Sumw2();
  ggH_jetpz1Opt->Sumw2();
  ggH_jetpz2Opt->Sumw2();
  ggH_jeteta1->Sumw2();
  ggH_jeteta2->Sumw2();
  ggH_jetphi1->Sumw2();
  ggH_jetphi2->Sumw2();
  ggH_jetDeltaPhi->Sumw2();
  ggH_jetDeltaEta->Sumw2();
  ggH_Mjj->Sumw2();
  ggH_jetmetDphi->Sumw2();
  ggH_jetmetDphi4->Sumw2();
  ggH_jetmediatorDphi->Sumw2();
  ggH_jetmediatorDeta->Sumw2();
  ggH_MT->Sumw2();
  ggH_MetMT->Sumw2();

  TH2F* ggH_mjj_detajj = new TH2F("ggH_mjj_detajj","ggH_mjj_detajj",15,400,2500,10,2,8);
  TH2F* ggH_mjj_jetmetDphi = new TH2F("ggH_mjj_jetmetDphi","ggH_mjj_jetmetDphi",15,400,2500,10,0.5,3.14);
  TH2F* ggH_mjj_jetDphi = new TH2F("ggH_mjj_jetDphi","ggH_mjj_jetDphi",15,400,2500,10,0.5,3.14);
  TH2F* ggH_mjj_MT = new TH2F("ggH_mjj_MT","ggH_mjj_MT",15,400,2500,10,0,500);
  TH2F* ggH_detajj_jetmetDphi = new TH2F("ggH_detajj_jetmetDphi","ggH_detajj_jetmetDphi",10,2,8,10,0.5,3.14);
  TH2F* ggH_detajj_jetDphi = new TH2F("ggH_detajj_jetDphi","ggH_detajj_jetDphi",10,2,8,10,0.5,3.14);
  TH2F* ggH_detajj_MT = new TH2F("ggH_detajj_MT","ggH_detajj_MT",10,2,8,10,0,500);

  ggH_mjj_detajj->Sumw2();
  ggH_mjj_jetmetDphi->Sumw2();
  ggH_mjj_jetDphi->Sumw2();
  ggH_mjj_MT->Sumw2();
  ggH_detajj_jetmetDphi->Sumw2();
  ggH_detajj_jetDphi->Sumw2();
  ggH_detajj_MT->Sumw2();

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
    fillHisto(qqH_bosonPt,*mediatorPt);
    fillHisto(qqH_bosonEta,*mediatorEta);
    fillHisto(qqH_bosonPhi,fabs(*mediatorPhi));
    fillHisto(qqH_bosonPz,mediator.Pz());
    fillHisto(qqH_bosonPzOPt,mediator.Pz()/mediator.Pt());
    fillHisto(qqH_jetpt1,jetpt->at(0));
    fillHisto(qqH_jetpt2,jetpt->at(1));
    fillHisto(qqH_jetpz1,jet1.Pz());
    fillHisto(qqH_jetpz2,jet2.Pz());
    fillHisto(qqH_jetpz1Opt,jet1.Pz()/jet1.Pt());
    fillHisto(qqH_jetpz2Opt,jet1.Pz()/jet1.Pt());
    fillHisto(qqH_jeteta1,jeteta->at(0));
    fillHisto(qqH_jeteta2,jeteta->at(1));
    fillHisto(qqH_jetphi1,fabs(jetphi->at(0)));
    fillHisto(qqH_jetphi2,fabs(jetphi->at(1)));
    fillHisto(qqH_jetDeltaPhi,fabs(jet1.DeltaPhi(jet2)));
    fillHisto(qqH_jetDeltaEta,fabs(jet1.Eta()-jet2.Eta()));
    fillHisto(qqH_Mjj,(jet1+jet2).M());
    fillHisto(qqH_jetmetDphi,*jmmdphi);
    fillHisto(qqH_jetmetDphi4,*jmmdphi4);
    fillHisto(qqH_jetmediatorDphi,fabs(mediator.DeltaPhi(jet1+jet2)));
    fillHisto(qqH_jetmediatorDeta,fabs(mediator.Eta()-(jet1+jet2).Eta()));
    fillHisto(qqH_MT,sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));    
    fillHisto(qqH_MetMT,sqrt(2*(jet1+jet2).Pt()*met.Pt()*(1-cos((jet1+jet2).DeltaPhi(met)))));

    // 2D part
    fillHisto(qqH_mjj_detajj,(jet1+jet2).M(),fabs(jet1.Eta()-jet2.Eta()));
    fillHisto(qqH_mjj_jetmetDphi,(jet1+jet2).M(),*jmmdphi);
    fillHisto(qqH_mjj_jetDphi,(jet1+jet2).M(),fabs(jet1.DeltaPhi(jet2)));
    fillHisto(qqH_mjj_MT,(jet1+jet2).M(),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));    
    fillHisto(qqH_detajj_jetmetDphi,fabs(jet1.Eta()-jet2.Eta()),*jmmdphi);
    fillHisto(qqH_detajj_jetDphi,fabs(jet1.Eta()-jet2.Eta()),fabs(jet1.DeltaPhi(jet2)));
    fillHisto(qqH_detajj_MT,fabs(jet1.Eta()-jet2.Eta()),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));

  }

  myReader.SetTree(ggHtree);
  myReader.SetEntry(0);

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
    fillHisto(ggH_bosonPt,*mediatorPt);
    fillHisto(ggH_bosonEta,*mediatorEta);
    fillHisto(ggH_bosonPhi,fabs(*mediatorPhi));
    fillHisto(ggH_bosonPz,mediator.Pz());
    fillHisto(ggH_bosonPzOPt,mediator.Pz()/mediator.Pt());
    fillHisto(ggH_jetpt1,jetpt->at(0));
    fillHisto(ggH_jetpt2,jetpt->at(1));
    fillHisto(ggH_jetpz1,jet1.Pz());
    fillHisto(ggH_jetpz2,jet2.Pz());
    fillHisto(ggH_jetpz1Opt,jet1.Pz()/jet1.Pt());
    fillHisto(ggH_jetpz2Opt,jet1.Pz()/jet1.Pt());
    fillHisto(ggH_jeteta1,jeteta->at(0));
    fillHisto(ggH_jeteta2,jeteta->at(1));
    fillHisto(ggH_jetphi1,fabs(jetphi->at(0)));
    fillHisto(ggH_jetphi2,fabs(jetphi->at(1)));
    fillHisto(ggH_jetDeltaPhi,fabs(jet1.DeltaPhi(jet2)));
    fillHisto(ggH_jetDeltaEta,fabs(jet1.Eta()-jet2.Eta()));
    fillHisto(ggH_Mjj,(jet1+jet2).M());
    fillHisto(ggH_jetmetDphi,*jmmdphi);
    fillHisto(ggH_jetmetDphi4,*jmmdphi4);
    fillHisto(ggH_jetmediatorDphi,fabs(mediator.DeltaPhi(jet1+jet2)));
    fillHisto(ggH_jetmediatorDeta,fabs(mediator.Eta()-(jet1+jet2).Eta()));
    fillHisto(ggH_MT,sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));    
    fillHisto(ggH_MetMT,sqrt(2*(jet1+jet2).Pt()*met.Pt()*(1-cos((jet1+jet2).DeltaPhi(met)))));

    fillHisto(ggH_mjj_detajj,(jet1+jet2).M(),fabs(jet1.Eta()-jet2.Eta()));
    fillHisto(ggH_mjj_jetmetDphi,(jet1+jet2).M(),*jmmdphi);
    fillHisto(ggH_mjj_jetDphi,(jet1+jet2).M(),fabs(jet1.DeltaPhi(jet2)));
    fillHisto(ggH_mjj_MT,(jet1+jet2).M(),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));    
    fillHisto(ggH_detajj_jetmetDphi,fabs(jet1.Eta()-jet2.Eta()),*jmmdphi);
    fillHisto(ggH_detajj_jetDphi,fabs(jet1.Eta()-jet2.Eta()),fabs(jet1.DeltaPhi(jet2)));
    fillHisto(ggH_detajj_MT,fabs(jet1.Eta()-jet2.Eta()),sqrt(2*jet1.Pt()*jet2.Pt()*(1-cos(jet1.DeltaPhi(jet2)))));
  
  }
  
  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();
  /////////////
  plotDistribution(canvas,qqH_bosonPt,ggH_bosonPt,"Mediator p_{T} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_bosonEta,ggH_bosonEta,"Mediator #eta",outputPlot);
  plotDistribution(canvas,qqH_bosonPhi,ggH_bosonPhi,"Mediator #phi",outputPlot);
  plotDistribution(canvas,qqH_bosonPz,ggH_bosonPz,"Mediator p_{z} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_bosonPzOPt,ggH_bosonPzOPt,"Mediator p_{z}/p_{T}",outputPlot);
  plotDistribution(canvas,qqH_jetpt1,ggH_jetpt1,"p^{j1}_{T} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_jetpt2,ggH_jetpt2,"p^{j2}_{T} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_jetpz1,ggH_jetpz1,"p^{j1}_{z} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_jetpz2,ggH_jetpz2,"p^{j2}_{z} [GeV]",outputPlot);
  plotDistribution(canvas,qqH_jetpz1Opt,ggH_jetpz1Opt,"p^{j1}_{z}/p^{j1}_{T}",outputPlot);
  plotDistribution(canvas,qqH_jetpz2Opt,ggH_jetpz2Opt,"p^{j1}_{z}/p^{j1}_{T}",outputPlot);
  plotDistribution(canvas,qqH_jeteta1,ggH_jeteta1,"#eta_{j1}",outputPlot);
  plotDistribution(canvas,qqH_jeteta2,ggH_jeteta2,"#eta_{j2}",outputPlot);
  plotDistribution(canvas,qqH_jetphi1,ggH_jetphi1,"#phi_{j1}",outputPlot);
  plotDistribution(canvas,qqH_jetphi2,ggH_jetphi2,"#phi_{j2}",outputPlot);
  plotDistribution(canvas,qqH_jetDeltaPhi,ggH_jetDeltaPhi,"#Delta#phi_{j1,j2}",outputPlot);
  plotDistribution(canvas,qqH_jetDeltaEta,ggH_jetDeltaEta,"#Delta#eta_{j1,j2}",outputPlot);
  plotDistribution(canvas,qqH_Mjj,ggH_Mjj,"M_{j1,j2}",outputPlot);
  plotDistribution(canvas,qqH_jetmetDphi,ggH_jetmetDphi,"#Delta#phi_{jet,met}",outputPlot);
  plotDistribution(canvas,qqH_jetmetDphi4,ggH_jetmetDphi4,"#Delta#phi_{4-jet,met}",outputPlot);
  plotDistribution(canvas,qqH_jetmediatorDeta,ggH_jetmediatorDeta,"#Delta#eta_{jet,med}",outputPlot);
  plotDistribution(canvas,qqH_jetmediatorDphi,ggH_jetmediatorDphi,"#Delta#phi_{jet,med}",outputPlot);
  plotDistribution(canvas,qqH_MetMT,ggH_MetMT,"m_{T}(jj,met) [GeV]",outputPlot);
  plotDistribution(canvas,qqH_MT,ggH_MT,"m_{T}(jj) [GeV]",outputPlot);
  /////////////
  plotDistribution(canvas,qqH_mjj_detajj,ggH_mjj_detajj,"m_{jj} [GeV]","#Delta#eta_{jj}",outputPlot);
  plotDistribution(canvas,qqH_mjj_jetmetDphi,ggH_mjj_jetmetDphi,"m_{jj} [GeV]","#Delta#phi_{jet,met}",outputPlot);
  plotDistribution(canvas,qqH_mjj_jetDphi,ggH_mjj_jetDphi,"m_{jj} [GeV]","#Delta#phi_{jj}",outputPlot);
  plotDistribution(canvas,qqH_mjj_MT,ggH_mjj_MT,"m_{jj} [GeV]","m_{T}(jj) [GeV]",outputPlot);
  plotDistribution(canvas,qqH_detajj_jetmetDphi,ggH_detajj_jetmetDphi,"#Delta#eta_{jj}","#Delta#phi_{jet,met}",outputPlot);
  plotDistribution(canvas,qqH_detajj_jetDphi,ggH_detajj_jetDphi,"#Delta#eta_{jj}","#Delta#phi_{jj}",outputPlot);
  plotDistribution(canvas,qqH_detajj_MT,ggH_detajj_MT,"#Delta#eta_{jj}","m_{T}(jj) [GeV]",outputPlot);
  
}
