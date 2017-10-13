#include "../CMS_lumi.h"

vector<float> bosonPt      {150.,200.,250.,300.,400.,500.,650.,1000};
vector<float> bosonPt_tail {150.,200.,250.,300.,400.,500.,650.,1000.};
vector<float> mjj_2D       {200.,600.,1000.,1400,2000.,5000.};

vector<float> jetpt2_2D    {40. ,80. ,120.,180.,250.,600.};
vector<float> mjj_bin      {200.,600.,1000.,1500.,2000.,2500.,5000.};

static float mjj            = 1000;
static float mjjrelaxed     = 200;
static float detajj         = 3.0;
static float detajjrelaxed  = 1.0;
static float leadingJetVBF  = 80;
static float trailingJetVBF = 40;
static float dphijj         = 1.5;
static float dphijjrelaxed  = 1.5;

enum class Sample {znn, zll, wjet, gam};
float lumi_ = 1;

void calculateSumWgt(TChain* tree, vector<double> & sumwgt_vec, const Sample & sample){

  TTreeReader reader (tree);
  TTreeReaderValue<float> xsec  (reader,"xsec");  
  TTreeReaderValue<float> wgt   (reader,"wgt");  
  TTreeReaderValue<int>   wzid  (reader,"wzid");  

  double sumwgt = 0;
  string currentFile = "";
  while(reader.Next()){

    //check if file name switched or not                                                                                                                                                               
    if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != ""){
      cout<<"Tree name "<<currentFile<<" wgt "<<sumwgt<<endl;
      sumwgt_vec.push_back(sumwgt);
      sumwgt = 0;
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
    }
    else if(currentFile == ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      sumwgt = 0;
    }
    
    // filter away bad events
    if((sample == Sample::znn or sample == Sample::zll) and fabs(*wzid) != 23) continue;
    else if(sample == Sample::wjet and fabs(*wzid) != 24) continue;
    else if(sample == Sample::gam and fabs(*wzid) != 22) continue;
    sumwgt += *wgt;
    
  }
  cout<<"Tree name "<<dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName()<<" wgt "<<sumwgt<<endl;
  sumwgt_vec.push_back(sumwgt);
  
}

void plotDistributions (const Sample & sample, TH1F* histo1, TH1F* histo2, TH1F* histo3 = NULL, const string & outputDIR = "", const string & xAxisLabel = "", const string & postfix = "" ){ 

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
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

  canvas->cd();

  if(TString(postfix).Contains("mjj")) 
    histo1->GetXaxis()->SetNdivisions(506);

  histo1->GetYaxis()->SetTitleOffset(1.1);
  histo1->GetXaxis()->SetTitleOffset(1.1);
  if(histo3)
    histo1->GetYaxis()->SetRangeUser(min(histo1->GetMinimum(),min(histo2->GetMinimum(),histo3->GetMinimum()))*0.7,
				     max(histo1->GetMaximum(),max(histo2->GetMaximum(),histo3->GetMaximum()))*1.3);
  else
    histo1->GetYaxis()->SetRangeUser(min(histo1->GetMinimum(),histo2->GetMinimum())*0.9,
				     max(histo1->GetMaximum(),histo2->GetMaximum())*1.1);
  
  histo1->SetLineColor(kBlack);
  histo1->SetMarkerColor(kBlack);
  histo1->SetLineWidth(2);
  histo1->SetMarkerSize(1);

  histo2->SetLineColor(kBlue);
  histo2->SetMarkerColor(kBlue);
  histo2->SetLineWidth(2);
  histo2->SetMarkerSize(1);

  if(histo3){
    histo3->SetLineColor(kRed);
    histo3->SetMarkerColor(kRed);
    histo3->SetLineWidth(2);
    histo3->SetMarkerSize(1);
  }


  TH1F* histo1_band = (TH1F*) histo1->Clone("histo1_band");
  histo1_band->SetFillColor(kBlack);
  histo1_band->SetFillStyle(3001);

  TH1F* histo2_band = (TH1F*) histo2->Clone("histo2_band");
  histo2_band->SetFillColor(kBlue);
  histo2_band->SetFillStyle(3001);

  TH1F* histo3_band = NULL;
  if(histo3){
    histo3_band = (TH1F*) histo3->Clone("histo3_band");
    histo3_band->SetFillColor(kRed);
    histo3_band->SetFillStyle(3001);
  }
  
  histo1->GetYaxis()->SetTitle("K-factor");
  histo1->GetXaxis()->SetLabelSize(0);
  histo1->GetXaxis()->SetTitleSize(0);
  histo1->Draw("hist");
  histo1_band->Draw("E2 same");
  histo2_band->Draw("E2 same");
  if(histo3_band)
    histo3_band->Draw("E2 same");

  histo1->Draw("hist same");
  histo2->Draw("hist same");
  if(histo3)
    histo3->Draw("hist same");

  TLegend leg (0.65,0.65,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  if(histo3_band){
    leg.AddEntry(histo1_band,"K-factor monojet","FL");
    leg.AddEntry(histo2_band,"K-factor VBF loose","FL");
    leg.AddEntry(histo3_band,"K-factor VBF tight","FL");
  }
  else{
    leg.AddEntry(histo1_band,"K-factor VBF loose","FL");
    leg.AddEntry(histo2_band,"K-factor VBF tight","FL");
  }
  leg.Draw("same");

  TLatex* latex = new TLatex ();
  latex->SetNDC();
  latex->SetTextSize(0.6*gPad->GetTopMargin());
  latex->SetTextFont(62);
  latex->SetTextAlign(11);
  if(sample == Sample::znn)
    latex->DrawLatex(0.15,0.95,"CMS Simulation Z+jets");
  else if(sample == Sample::zll)
    latex->DrawLatex(0.15,0.95,"CMS Simulation DY+jets");
  else if(sample == Sample::wjet)
    latex->DrawLatex(0.15,0.95,"CMS Simulation W+jets");
  else if(sample == Sample::gam)
    latex->DrawLatex(0.15,0.95,"CMS Simulation #gamma+jets");

  pad2->Draw();
  pad2->cd();

  //ratio plot
  TH1* ratio_1 = (TH1*) histo2->Clone("ratio_2");
  ratio_1->Divide(histo1);
  TH1* ratio_2 = NULL;
  if(histo3){
    ratio_2 = (TH1*) histo3->Clone("ratio_3");
    ratio_2->Divide(histo1);
  }

  ratio_1->SetLineColor(kBlue);
  ratio_1->SetMarkerColor(kBlue);
  ratio_1->GetYaxis()->SetTitle("Ratio");
  if(sample == Sample::gam)
    ratio_1->GetYaxis()->SetRangeUser(0.85,1.5);
  else
    ratio_1->GetYaxis()->SetRangeUser(0.85,1.25);
  
  ratio_1->GetXaxis()->SetTitle(xAxisLabel.c_str());
  ratio_1->GetYaxis()->CenterTitle();
  ratio_1->GetYaxis()->SetTitleOffset(1.5);
  ratio_1->GetYaxis()->SetLabelSize(0.04);
  ratio_1->GetYaxis()->SetTitleSize(0.04);
  ratio_1->GetXaxis()->SetLabelSize(0.04);
  ratio_1->GetXaxis()->SetTitleSize(0.05);  
  ratio_1->GetYaxis()->SetNdivisions(5);
  ratio_1->Draw("hist");

  if(ratio_2){
    ratio_2->SetLineColor(kRed);
    ratio_2->SetMarkerColor(kRed);
    ratio_2->Draw("hist same");
  }

  TH1F* ratio_1_band = (TH1F*) ratio_1->Clone("ratio_3_band");
  ratio_1_band->SetFillColor(kBlue);
  ratio_1_band->SetFillStyle(3001);
  ratio_1_band->Draw("E2 same");

  if(ratio_2){
    TH1F* ratio_2_band = (TH1F*) ratio_2->Clone("ratio_2_band");
    ratio_2_band->SetFillColor(kOrange+1);
    ratio_2_band->SetFillStyle(3001);
    ratio_2_band->Draw("E2 same");
  }

  if(ratio_2)
    ratio_2->Draw("hist same");
  ratio_1->Draw("hist same");
  
  if(sample == Sample::znn){
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_znn.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_znn.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::zll){
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_zll.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_zll.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::wjet){
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_wjet.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_wjet.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::gam){
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_gam.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_gam.pdf").c_str(),"pdf");
  }
}


/// -----
void plotDistributions(const Sample & sample, vector<TH1F*> & histos, const string & outputDIR, const vector<TString> & labels, const string & xAxisLabel,
		       const string & postfix){

  TCanvas* canvas = new TCanvas("canvas", "canvas", 700, 700);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();

  histos.at(0)->Draw("hist");
  int icolor = 1;
  double min = 10000;
  double max = -10000;

  if(TString(postfix).Contains("mjj")) 
    histos.at(0)->GetXaxis()->SetNdivisions(506);

  for(auto hist : histos){
    if(icolor == 3) icolor++;
    if(icolor == 5) icolor++;
    if(icolor == 6) icolor++;
    hist->GetXaxis()->SetTitle(xAxisLabel.c_str());
    hist->SetLineColor(icolor);
    hist->SetMarkerColor(icolor);
    hist->SetMarkerStyle(20);
    hist->SetMarkerSize(1);
    hist->SetLineWidth(2);
    hist->GetYaxis()->SetTitle("k-factor");
    hist->Draw("hist same");
    hist->Draw("P same");
    icolor++;
    if(hist->GetMinimum() < min) min = hist->GetMinimum() ;
    if(hist->GetMaximum() > max) max = hist->GetMaximum() ;

  }

  histos.at(0)->GetYaxis()->SetRangeUser(0.3,2.8);
  
  TLegend leg (0.7,0.7,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  int iobj = 0;
  for(auto lab : labels){
    leg.AddEntry(histos.at(iobj),lab,"EP");
    iobj++;
  }
  leg.Draw("same");

  TLatex* latex = new TLatex ();
  latex->SetNDC();
  latex->SetTextSize(0.6*gPad->GetTopMargin());
  latex->SetTextFont(62);
  latex->SetTextAlign(11);
  if(sample == Sample::znn)
    latex->DrawLatex(0.15,0.95,"CMS Simulation Z+jets");
  else if(sample == Sample::zll)
    latex->DrawLatex(0.15,0.95,"CMS Simulation DY+jets");
  else if(sample == Sample::wjet)
    latex->DrawLatex(0.15,0.95,"CMS Simulation W+jets");
  else if(sample == Sample::gam)
    latex->DrawLatex(0.15,0.95,"CMS Simulation #gamma+jets");
  
  if(sample == Sample::znn){
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_znn.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_znn.pdf").c_str(),"pdf");
    }
  else if(sample == Sample::zll){
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_zll.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_zll.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::wjet){
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_wjet.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_wjet.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::gam){
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_gam.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/kfactor_"+postfix+"_gam.pdf").c_str(),"pdf");
  }

}

///// ----
void makeKFactorVBF(string inputDIR_LO, string inputDIR_NLO, string outputDIR, Sample sample){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  float scale_lo = 1;
  float scale_nlo = 1;
  if(sample == Sample::znn or sample == Sample::zll)
    scale_nlo = 3;

  system(("mkdir -p "+outputDIR).c_str());

  TChain* tree_LO = new TChain("gentree/tree");
  tree_LO->Add((inputDIR_LO+"/*root").c_str());

  TChain* tree_NLO = new TChain("gentree/tree");
  tree_NLO->Add((inputDIR_NLO+"/*root").c_str());

  /// -----
  string postfix;
  if(sample == Sample::znn)       postfix = "znn";
  else if(sample == Sample::zll)  postfix = "zll";
  else if(sample == Sample::wjet) postfix = "wjet";
  else if(sample == Sample::gam)  postfix = "gam";

  // histograms
  TH1F* bosonPt_LO_monojet  = new TH1F("bosonPt_LO_monojet","",bosonPt.size()-1,&bosonPt[0]);
  TH1F* bosonPt_NLO_monojet = new TH1F("bosonPt_NLO_monojet","",bosonPt.size()-1,&bosonPt[0]);
  TH1F* bosonPt_LO_vbf      = new TH1F("bosonPt_LO_vbf","", bosonPt.size()-1,&bosonPt[0]);
  TH1F* bosonPt_NLO_vbf     = new TH1F("bosonPt_NLO_vbf","", bosonPt.size()-1,&bosonPt[0]);
  TH1F* bosonPt_LO_vbf_tight      = new TH1F("bosonPt_LO_vbf_tight","", bosonPt.size()-1,&bosonPt[0]);
  TH1F* bosonPt_NLO_vbf_tight     = new TH1F("bosonPt_NLO_vbf_tight","", bosonPt.size()-1,&bosonPt[0]);

  TH1F* mjj_LO_vbf      = new TH1F("mjj_LO_vbf","", mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* mjj_NLO_vbf     = new TH1F("mjj_NLO_vbf","", mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* mjj_LO_vbf_tight      = new TH1F("mjj_LO_vbf_tight","", mjj_bin.size()-1,&mjj_bin[0]);
  TH1F* mjj_NLO_vbf_tight     = new TH1F("mjj_NLO_vbf_tight","", mjj_bin.size()-1,&mjj_bin[0]);

  bosonPt_LO_monojet->Sumw2();
  bosonPt_NLO_monojet->Sumw2();
  bosonPt_LO_vbf->Sumw2();
  bosonPt_NLO_vbf->Sumw2();
  bosonPt_LO_vbf_tight->Sumw2();
  bosonPt_NLO_vbf_tight->Sumw2();

  mjj_LO_vbf->Sumw2();
  mjj_NLO_vbf->Sumw2();
  mjj_LO_vbf_tight->Sumw2();
  mjj_NLO_vbf_tight->Sumw2();

  vector<TH1F*> bosonpt_mjj_LO_vbf;
  vector<TH1F*> bosonpt_mjj_NLO_vbf;
  for(size_t iBin = 0; iBin < mjj_2D.size()-1; iBin++){
    if(mjj_2D.at(iBin) < 1500){ // use a more granular bin in boson-pt
      bosonpt_mjj_LO_vbf.push_back(new TH1F(Form("bosonpt_LO_mjj_%d_%d",int(mjj_2D.at(iBin)),int(mjj_2D.at(iBin+1))),"",bosonPt.size()-1,&bosonPt[0]));
      bosonpt_mjj_NLO_vbf.push_back(new TH1F(Form("bosonpt_NLO_mjj_%d_%d",int(mjj_2D.at(iBin)),int(mjj_2D.at(iBin+1))),"",bosonPt.size()-1,&bosonPt[0]));
    }
    else{
      bosonpt_mjj_LO_vbf.push_back(new TH1F(Form("bosonpt_LO_mjj_%d_%d",int(mjj_2D.at(iBin)),int(mjj_2D.at(iBin+1))),"",bosonPt_tail.size()-1,&bosonPt_tail[0]));
      bosonpt_mjj_NLO_vbf.push_back(new TH1F(Form("bosonpt_NLO_mjj_%d_%d",int(mjj_2D.at(iBin)),int(mjj_2D.at(iBin+1))),"",bosonPt_tail.size()-1,&bosonPt_tail[0]));
    }
    bosonpt_mjj_LO_vbf.back()->Sumw2();
    bosonpt_mjj_NLO_vbf.back()->Sumw2();    
  }

  vector<TH1F*> bosonpt_jetpt2_LO_vbf;
  vector<TH1F*> bosonpt_jetpt2_NLO_vbf;
  for(size_t iBin = 0; iBin < jetpt2_2D.size()-1; iBin++){
    if(jetpt2_2D.at(iBin) < 1500){ // use a more granular bin in boson-pt
      bosonpt_jetpt2_LO_vbf.push_back(new TH1F(Form("bosonpt_LO_jetpt2_%d_%d",int(jetpt2_2D.at(iBin)),int(jetpt2_2D.at(iBin+1))),"",bosonPt.size()-1,&bosonPt[0]));
      bosonpt_jetpt2_NLO_vbf.push_back(new TH1F(Form("bosonpt_NLO_jetpt2_%d_%d",int(jetpt2_2D.at(iBin)),int(jetpt2_2D.at(iBin+1))),"",bosonPt.size()-1,&bosonPt[0]));
    }
    else{
      bosonpt_jetpt2_LO_vbf.push_back(new TH1F(Form("bosonpt_LO_jetpt2_%d_%d",int(jetpt2_2D.at(iBin)),int(jetpt2_2D.at(iBin+1))),"",bosonPt_tail.size()-1,&bosonPt_tail[0]));
      bosonpt_jetpt2_NLO_vbf.push_back(new TH1F(Form("bosonpt_NLO_jetpt2_%d_%d",int(jetpt2_2D.at(iBin)),int(jetpt2_2D.at(iBin+1))),"",bosonPt_tail.size()-1,&bosonPt_tail[0]));
    }
    bosonpt_jetpt2_LO_vbf.back()->Sumw2();
    bosonpt_jetpt2_NLO_vbf.back()->Sumw2();    
  }
  
  // outputFile

  TFile* output = new TFile((outputDIR+"/kfactor_"+postfix+".root").c_str(),"RECREATE");
  output->cd();

  // calculate sumwgt
  vector<double> sumwgt_lo;
  calculateSumWgt(tree_LO,sumwgt_lo,sample);
  vector<double> sumwgt_nlo;
  calculateSumWgt(tree_NLO,sumwgt_nlo,sample);


  // Loop on LO trees
  int ifile = 0;
  TTreeReader reader (tree_LO);
  TTreeReaderValue<float> xsec  (reader,"xsec");  
  TTreeReaderValue<float> wgt   (reader,"wgt");  
  TTreeReaderValue<int>   wzid  (reader,"wzid");  
  TTreeReaderValue<int>   l1id  (reader,"l1id");  
  TTreeReaderValue<int>   l2id  (reader,"l2id");  
  TTreeReaderValue<unsigned int>   njetsinc  (reader,"njetsinc");  
  TTreeReaderValue<unsigned int>   njets  (reader,"njets");  
  TTreeReaderValue<float> wzmass  (reader,"wzmass");  
  TTreeReaderValue<float> wzpt    (reader,"wzpt");  
  TTreeReaderValue<float> wzpt_lhe (reader,"mvpt");  
  TTreeReaderValue<float> wzeta   (reader,"wzeta");  
  TTreeReaderValue<float> wzphi   (reader,"wzphi");  
  
  TTreeReaderValue<float> l1pt   (reader,"l1pt");  
  TTreeReaderValue<float> l1eta  (reader,"l1eta");  
  TTreeReaderValue<float> l1phi  (reader,"l1phi");  
  TTreeReaderValue<float> l2pt   (reader,"l2pt");  
  TTreeReaderValue<float> l2eta  (reader,"l2eta");  
  TTreeReaderValue<float> l2phi  (reader,"l2phi");  
  
  TTreeReaderValue<vector<float> > jetpt  (reader,"jetpt");  
  TTreeReaderValue<vector<float> > jeteta (reader,"jeteta");  
  TTreeReaderValue<vector<float> > jetphi (reader,"jetphi");  
  TTreeReaderValue<vector<float> > jetmass (reader,"jetmass");  

  cout<<"Start loop on LO trees "<<endl;
  string currentFile = "";
  int nEvents = 0;
  while(reader.Next()){

    //check if file name switched or not                                                                                                                                                               
    if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      cout<<"Looping on: "<<currentFile<<endl;
      ifile++;
    }
    else if(currentFile == ""){
      currentFile =  dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      cout<<"Looping on: "<<currentFile<<endl;
      ifile= 0;
    }

    if(*wzpt < bosonPt.front()) continue;
    if(*wzpt_lhe < bosonPt.front()) continue;

    /// --- filter away bad events
    if((sample == Sample::znn or sample == Sample::zll) and fabs(*wzid) != 23) continue;
    else if(sample == Sample::wjet and fabs(*wzid) != 24) continue;
    else if(sample == Sample::gam  and fabs(*wzid) != 22) continue;

    /// ---
    if(sample == Sample::wjet and ((fabs(*l1id) == 15 and fabs(*l2id) == 16) or (fabs(*l2id) == 15 and fabs(*l1id) == 16))) continue; // skip taus 
    if(sample == Sample::zll  and ((fabs(*l1id) == 15 and fabs(*l2id) == 15))) continue; // skip taus

    // for photons only look at the central region                                                                                                                                                  
    if(sample == Sample::gam and fabs(*wzeta) > 1.449) continue;

    /////
    if(sample == Sample::zll and fabs(*l1id) == 11 and fabs(*l1eta) > 2.5) continue;
    if(sample == Sample::zll and fabs(*l2id) == 11 and fabs(*l1eta) > 2.5) continue;
    if(sample == Sample::zll and fabs(*l1eta) > 2.4) continue;
    if(sample == Sample::zll and fabs(*l1eta) > 2.4) continue;

    ////
    if(sample == Sample::wjet and fabs(*l1id) == 11 and fabs(*l1eta) > 2.5) continue;
    if(sample == Sample::wjet and fabs(*l2id) == 11 and fabs(*l1eta) > 2.5) continue;
    if(sample == Sample::wjet and fabs(*l1id) == 13 and fabs(*l1eta) > 2.4) continue;
    if(sample == Sample::wjet and fabs(*l2id) == 13 and fabs(*l1eta) > 2.4) continue;

    ///
    if(sample == Sample::wjet and fabs(*l1id) == 11 and *l1pt < 40) continue;
    if(sample == Sample::wjet and fabs(*l1id) == 13 and *l1pt < 20) continue;
    
    ///
    TLorentzVector lepton1, lepton2;
    lepton1.SetPtEtaPhiM(*l1pt,*l1eta,*l1phi,0.);
    lepton2.SetPtEtaPhiM(*l2pt,*l2eta,*l2phi,0.);
    ///
    if(sample == Sample::zll and max(*l1pt,*l2pt) < 20) continue;
    if(sample == Sample::zll and min(*l1pt,*l2pt) < 10) continue;
    if(sample == Sample::zll and (lepton1+lepton2).M() < 60) continue;
    if(sample == Sample::zll and (lepton1+lepton2).M() > 120) continue;

    ///
    vector<TLorentzVector> jets;
    for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
      TLorentzVector jet4V; jet4V.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetmass->at(ijet));
      if(sample == Sample::wjet or sample == Sample::znn or sample == Sample::zll){
	float dphi1 = fabs(jetphi->at(ijet)-*l1phi);
	float dphi2 = fabs(jetphi->at(ijet)-*l2phi);
	if(dphi1 > TMath::Pi())
	  dphi1 = 2*TMath::Pi()-dphi1;
	if(dphi2 > TMath::Pi())
	  dphi2 = 2*TMath::Pi()-dphi2;
	if(sqrt(dphi1*dphi1+fabs(jeteta->at(ijet)-*l1eta)*fabs(jeteta->at(ijet)-*l1eta)) > 0.4 and
	   sqrt(dphi2*dphi2+fabs(jeteta->at(ijet)-*l2eta)*fabs(jeteta->at(ijet)-*l2eta)) > 0.4
	   ){ // check cleaning
	  jets.push_back(jet4V);
	}
      }
      else if(sample == Sample::gam){
	float dphi1 = fabs(jetphi->at(ijet)-*wzphi);
	if(dphi1 > TMath::Pi())
	  dphi1 = 2*TMath::Pi()-dphi1;
	if(sqrt(dphi1*dphi1+fabs(jeteta->at(ijet)-*wzeta)*fabs(jeteta->at(ijet)-*wzeta)) > 0.4)
	  jets.push_back(jet4V);
        }
    }
    
    if(jets.size() < 1) continue;
    
    // calculate min-dphi at gen level where met is boson 4V
    float mindphi   = 100;
    for(size_t ijet = 0; ijet < jets.size(); ijet++){
      if(ijet > 3) break; // limiting min dphi to first 4 leading jets
      float dphi = fabs(*wzphi-jets.at(ijet).Phi());
      if(dphi > TMath::Pi())
	dphi = 2*TMath::Pi()-dphi;
      if(dphi < mindphi)
	mindphi = dphi;	
    }
    
    
    // apply monojet-like selections
    if(jets.size() >= 1 and 
       jets.at(0).Pt() > 100 and 
       fabs(jets.at(0).Eta()) < 2.5 and 
       *wzpt >= 150 and 
       mindphi > 0.5)
      
      bosonPt_LO_monojet->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_lo/sumwgt_lo.at(ifile));
    
    if(jets.size() >= 2 and 
       jets.at(0).Pt() > leadingJetVBF and 
       jets.at(1).Pt() > trailingJetVBF and 
       fabs(jets.at(0).Eta()) < 4.7 and  
       fabs(jets.at(1).Eta()) < 4.7 and 
       mindphi > 0.5){      
      
      if(fabs(jets.at(1).Eta()) > 3 and fabs(jets.at(0).Eta()) > 3) continue;
      if((jets.at(0)+jets.at(1)).M() < mjjrelaxed) continue;
      if(jets.at(0).Eta()*jets.at(1).Eta() > 0) continue; 
      if(fabs(jets.at(0).Eta()-jets.at(1).Eta()) < detajjrelaxed) continue;
      float deltaPhi = fabs(jets.at(0).Phi()-jets.at(1).Phi());
      if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi()-deltaPhi;
      if(deltaPhi < dphijjrelaxed){
	
	bosonPt_LO_vbf->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_lo/sumwgt_lo.at(ifile));      
	
	if(*wzpt > 250) // above 250 GeV
	  mjj_LO_vbf->Fill((jets.at(0)+jets.at(1)).M(),lumi_*(*wgt)*(*xsec)*scale_lo/sumwgt_lo.at(ifile));      
	
	for(size_t iBin = 0; iBin < mjj_2D.size()-1; iBin++){
	  if((jets.at(0)+jets.at(1)).M() >= mjj_2D.at(iBin) and (jets.at(0)+jets.at(1)).M() < mjj_2D.at(iBin+1)){
	    bosonpt_mjj_LO_vbf.at(iBin)->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_lo/sumwgt_lo.at(ifile));
	    break;
	  }
	}
	
	for(size_t iBin = 0; iBin < jetpt2_2D.size()-1; iBin++){
	  if(jets.at(1).Pt() >= jetpt2_2D.at(iBin) and jets.at(1).Pt() < jetpt2_2D.at(iBin+1)){
	    bosonpt_jetpt2_LO_vbf.at(iBin)->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_lo/sumwgt_lo.at(ifile));
	    break;
	  }
	}
      }
      
      if(fabs(jets.at(0).Eta()-jets.at(1).Eta()) < detajj) continue;
      if(deltaPhi < dphijj){
	if((jets.at(0)+jets.at(1)).M() > mjj) // above the Mjj cut
	  bosonPt_LO_vbf_tight->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_lo/sumwgt_lo.at(ifile));      	
	if(*wzpt > 250) // above 250 GeV
	  mjj_LO_vbf_tight->Fill((jets.at(0)+jets.at(1)).M(),lumi_*(*wgt)*(*xsec)*scale_lo/sumwgt_lo.at(ifile));      	
      }
    }
  }

  //// -----
  reader.SetTree(tree_NLO);
  reader.SetEntry(0);
    
  cout<<"Start Loop on NLO file "<<endl;
  ifile = 0;
  currentFile = "";

  while(reader.Next()){

    //check if file name switched or not                                                                                                                                                               
    if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      cout<<"Looping on: "<<currentFile<<endl;
      ifile++;
    }
    else if(currentFile == ""){
      currentFile =  dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      cout<<"Looping on: "<<currentFile<<endl;
      ifile= 0;
    }
    
    if(*wzpt < bosonPt.front()) continue;
    if(*wzpt_lhe < bosonPt.front()) continue;
    
    // filter away bad events with no matching
    if((sample == Sample::znn or sample == Sample::zll) and fabs(*wzid) != 23) continue;
    else if(sample == Sample::wjet and fabs(*wzid) != 24) continue;
    else if(sample == Sample::gam  and fabs(*wzid) != 22) continue;

    /// ---
    if(sample == Sample::wjet and ((fabs(*l1id) == 15 and fabs(*l2id) == 16) or (fabs(*l2id) == 15 and fabs(*l1id) == 16))) continue; // skip taus 
    if(sample == Sample::zll  and ((fabs(*l1id) == 15 and fabs(*l2id) == 15))) continue; // skip taus
    

    // for photons only look at the central region                                                                                                                                                  
    if(sample == Sample::gam and fabs(*wzeta) > 1.449) continue;

    /////
    if(sample == Sample::zll and fabs(*l1id) == 11 and fabs(*l1eta) > 2.5) continue;
    if(sample == Sample::zll and fabs(*l2id) == 11 and fabs(*l1eta) > 2.5) continue;
    if(sample == Sample::zll and fabs(*l1eta) > 2.4) continue;
    if(sample == Sample::zll and fabs(*l1eta) > 2.4) continue;

    ////
    if(sample == Sample::wjet and fabs(*l1id) == 11 and fabs(*l1eta) > 2.5) continue;
    if(sample == Sample::wjet and fabs(*l2id) == 11 and fabs(*l1eta) > 2.5) continue;
    if(sample == Sample::wjet and fabs(*l1id) == 13 and fabs(*l1eta) > 2.4) continue;
    if(sample == Sample::wjet and fabs(*l2id) == 13 and fabs(*l1eta) > 2.4) continue;

    ///
    if(sample == Sample::wjet and fabs(*l1id) == 11 and *l1pt < 40) continue;
    if(sample == Sample::wjet and fabs(*l1id) == 13 and *l1pt < 20) continue;


    TLorentzVector lepton1, lepton2;
    lepton1.SetPtEtaPhiM(*l1pt,*l1eta,*l1phi,0.);
    lepton2.SetPtEtaPhiM(*l2pt,*l2eta,*l2phi,0.);
    
    if(sample == Sample::zll and max(*l1pt,*l2pt) < 20) continue;
    if(sample == Sample::zll and min(*l1pt,*l2pt) < 10) continue;
    if(sample == Sample::zll and (lepton1+lepton2).M() < 60) continue;
    if(sample == Sample::zll and (lepton1+lepton2).M() > 120) continue;
    
    // avoid events with large weight
    if(sample == Sample::wjet and *wzpt_lhe >= 400 and fabs(*wgt) > 100) continue;
    if(sample == Sample::wjet and *wzpt_lhe >= 600 and fabs(*wgt) > 10)  continue;
    if(sample == Sample::znn  and *wzpt_lhe >= 400 and fabs(*wgt) > 10)  continue;
    if(sample == Sample::zll  and *wzpt_lhe >= 400 and fabs(*wgt) > 10)  continue;
    if(sample == Sample::znn  and *wzpt_lhe >= 650 and fabs(*wgt) > 1)   continue;
    if(sample == Sample::zll  and *wzpt_lhe >= 650 and fabs(*wgt) > 1)   continue;
    
    vector<TLorentzVector> jets;
    for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
      TLorentzVector jet4V; jet4V.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetmass->at(ijet));
      if(sample == Sample::wjet or sample == Sample::znn or sample == Sample::zll){
	float dphi1 = fabs(jetphi->at(ijet)-*l1phi);
	float dphi2 = fabs(jetphi->at(ijet)-*l2phi);
	if(dphi1 > TMath::Pi())
	  dphi1 = 2*TMath::Pi()-dphi1;
	if(dphi2 > TMath::Pi())
	  dphi2 = 2*TMath::Pi()-dphi2;
	if(sqrt(dphi1*dphi1+fabs(jeteta->at(ijet)-*l1eta)*fabs(jeteta->at(ijet)-*l1eta)) > 0.4 and
	   sqrt(dphi2*dphi2+fabs(jeteta->at(ijet)-*l2eta)*fabs(jeteta->at(ijet)-*l2eta)) > 0.4){ // check cleaning
	  jets.push_back(jet4V);
	}
      }
      else if(sample == Sample::gam){
	float dphi1 = fabs(jetphi->at(ijet)-*wzphi);
	if(dphi1 > TMath::Pi())
	  dphi1 = 2*TMath::Pi()-dphi1;
	if(sqrt(dphi1*dphi1+fabs(jeteta->at(ijet)-*wzeta)*fabs(jeteta->at(ijet)-*wzeta)) > 0.4)
	  jets.push_back(jet4V);
      }     
    }
    
    if(jets.size() < 1) continue;     
      // calculate min-dphi at gen level where met is boson 4V
      float mindphi   = 100;
      for(size_t ijet = 0; ijet < jets.size(); ijet++){
	if(ijet > 3) break; // limiting min dphi to first 4 leading jets
	float dphi = fabs(*wzphi-jets.at(ijet).Phi());
	if(dphi > TMath::Pi())
	  dphi = 2*TMath::Pi()-dphi;
	if(dphi < mindphi)
	  mindphi = dphi;	
      }
      
      
      // apply monojet-like selections
      if(jets.size() >= 1 and 
	 jets.at(0).Pt() > 100 and 
	 fabs(jets.at(0).Eta()) < 2.5 and 
	 *wzpt >= 150 and 
	 mindphi > 0.5){
	bosonPt_NLO_monojet->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_nlo/sumwgt_nlo.at(ifile));
      }
      
      if(jets.size() >= 2 and 
	 jets.at(0).Pt() > leadingJetVBF and 
	 jets.at(1).Pt() > trailingJetVBF and 
	 fabs(jets.at(0).Eta()) < 4.7 and  
	 fabs(jets.at(1).Eta()) < 4.7 and 
	 mindphi > 0.5){      

	if(fabs(jets.at(1).Eta()) > 3 and fabs(jets.at(0).Eta()) > 3) continue;
	if((jets.at(0)+jets.at(1)).M() < mjjrelaxed) continue;
	if(jets.at(0).Eta()*jets.at(1).Eta() > 0) continue; 
	if(fabs(jets.at(0).Eta()-jets.at(1).Eta()) < detajjrelaxed) continue;
	float deltaPhi = fabs(jets.at(0).Phi()-jets.at(1).Phi());
	if(deltaPhi > TMath::Pi()) deltaPhi = 2*TMath::Pi()-deltaPhi;
	if(deltaPhi < dphijjrelaxed){
	  
	  bosonPt_NLO_vbf->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_nlo/sumwgt_nlo.at(ifile));      

	  if(*wzpt > 250) // above 250 GeV
	    mjj_NLO_vbf->Fill((jets.at(0)+jets.at(1)).M(),lumi_*(*wgt)*(*xsec)*scale_nlo/sumwgt_nlo.at(ifile));      
	  
	  for(size_t iBin = 0; iBin < mjj_2D.size()-1; iBin++){
	    if((jets.at(0)+jets.at(1)).M() >= mjj_2D.at(iBin) and (jets.at(0)+jets.at(1)).M() < mjj_2D.at(iBin+1)){
	      bosonpt_mjj_NLO_vbf.at(iBin)->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_nlo/sumwgt_nlo.at(ifile));
	      break;
	    }
	  }
	  
	  for(size_t iBin = 0; iBin < jetpt2_2D.size()-1; iBin++){
	    if(jets.at(1).Pt() >= jetpt2_2D.at(iBin) and jets.at(1).Pt() < jetpt2_2D.at(iBin+1)){
	      bosonpt_jetpt2_NLO_vbf.at(iBin)->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_nlo/sumwgt_nlo.at(ifile));
	      break;
	    }
	  }
	}
	
	if(fabs(jets.at(0).Eta()-jets.at(1).Eta()) < detajj) continue;
	if(deltaPhi < dphijj){
	  if((jets.at(0)+jets.at(1)).M() > mjj)
	    bosonPt_NLO_vbf_tight->Fill(*wzpt,lumi_*(*wgt)*(*xsec)*scale_nlo/sumwgt_nlo.at(ifile));	
	  if(*wzpt > 250) 
	    mjj_NLO_vbf_tight->Fill((jets.at(0)+jets.at(1)).M(),lumi_*(*wgt)*(*xsec)*scale_nlo/sumwgt_nlo.at(ifile));	
	}
      }
  }
  
  
  ///// -----
  TH1F* kfactor_monojet = (TH1F*) bosonPt_NLO_monojet->Clone("kfactor_monojet");
  kfactor_monojet->Divide(bosonPt_LO_monojet);
  TH1F* kfactor_vbf_tight = (TH1F*) bosonPt_NLO_vbf_tight->Clone("kfactor_vbf_tight");
  kfactor_vbf_tight->Divide(bosonPt_LO_vbf_tight);
  TH1F* kfactor_vbf = (TH1F*) bosonPt_NLO_vbf->Clone("kfactor_vbf");
  kfactor_vbf->Divide(bosonPt_LO_vbf);
  plotDistributions(sample,kfactor_monojet,kfactor_vbf,kfactor_vbf_tight,outputDIR,"boson p_{T} [GeV]","bosonpt");

  TH1F* kfactor_vbf_mjj = (TH1F*) mjj_NLO_vbf->Clone("kfactor_vbf_mjj");
  kfactor_vbf_mjj->Divide(mjj_LO_vbf);
  TH1F* kfactor_vbf_tight_mjj = (TH1F*) mjj_NLO_vbf_tight->Clone("kfactor_vbf_tight_mjj");
  kfactor_vbf_tight_mjj->Divide(mjj_LO_vbf_tight);
  plotDistributions(sample,kfactor_vbf_mjj,kfactor_vbf_tight_mjj,NULL,outputDIR,"M_{jj} [GeV]","mjj");

  ///// ---- 
  vector<TH1F*> kfactor_vbf_vs_mjj;
  vector<TString> label_vs_mjj;
  
  ///// ---- 
  for(size_t ihist = 0; ihist < bosonpt_mjj_NLO_vbf.size(); ihist++){
    kfactor_vbf_vs_mjj.push_back((TH1F*) bosonpt_mjj_NLO_vbf.at(ihist)->Clone(Form("kfactor_vbf_mjj_%d_%d",int(mjj_2D.at(ihist)),int(mjj_2D.at(ihist+1)))));
    kfactor_vbf_vs_mjj.back()->Divide(bosonpt_mjj_LO_vbf.at(ihist));
    label_vs_mjj.push_back(Form("%d < M_{jj} < %d",int(mjj_2D.at(ihist)),int(mjj_2D.at(ihist+1))));
  }

  plotDistributions(sample,kfactor_vbf_vs_mjj,outputDIR,label_vs_mjj,"boson p_{T} [GeV]","bosonpt_mjj");


  ///// ---- 
  vector<TH1F*> kfactor_vbf_vs_jetpt2;
  vector<TString> label_vs_jetpt2;
  
  ///// ---- 
  for(size_t ihist = 0; ihist < bosonpt_jetpt2_NLO_vbf.size(); ihist++){
    kfactor_vbf_vs_jetpt2.push_back((TH1F*) bosonpt_jetpt2_NLO_vbf.at(ihist)->Clone(Form("kfactor_vbf_jetpt2_%d_%d",int(jetpt2_2D.at(ihist)),int(jetpt2_2D.at(ihist+1)))));
    kfactor_vbf_vs_jetpt2.back()->Divide(bosonpt_jetpt2_LO_vbf.at(ihist));
    label_vs_jetpt2.push_back(Form("%d < M_{jj} < %d",int(jetpt2_2D.at(ihist)),int(jetpt2_2D.at(ihist+1))));
  }

  plotDistributions(sample,kfactor_vbf_vs_jetpt2,outputDIR,label_vs_jetpt2,"p_{T}^{j2} [GeV]","bosonpt_pt2");
  
  bosonPt_NLO_monojet->Write();
  bosonPt_NLO_vbf->Write();
  bosonPt_NLO_vbf->Write();
  bosonPt_LO_monojet->Write();
  bosonPt_LO_vbf->Write();
  bosonPt_LO_vbf->Write();

  mjj_NLO_vbf->Write();
  mjj_NLO_vbf->Write();
  mjj_LO_vbf->Write();
  mjj_LO_vbf->Write();

  for(auto hist: kfactor_vbf_vs_mjj)
    hist->Write();

  for(auto hist: kfactor_vbf_vs_jetpt2)
    hist->Write();

  output->Close();
}
