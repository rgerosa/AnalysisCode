#include "../CMS_lumi.h"

enum class Sample {zjet,wjet};
enum class Category {monojet, VBFrelaxed, VBF};

static float luminosity = 35.9;

////////////
void calculateSumWg(TChain* chain, vector<double> & sumwgt_vec, Sample sample){

  TTreeReader reader (chain);
  TTreeReaderValue<float> xsec  (reader,"xsec");
  TTreeReaderValue<float> wgt   (reader,"wgt");
  TTreeReaderValue<int>   wzid  (reader,"wzid");
  
  double sumwgt = 0;
  string currentFile = "";
  while(reader.Next()){    

    //check if file name switched or not
    if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != ""){
      sumwgt_vec.push_back(sumwgt);
      sumwgt = 0;
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
    }
    else if(currentFile == ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      sumwgt = 0;
    }

    // filter away bad events with no matching                                                                                                                                                        
    if(sample == Sample::zjet and fabs(*wzid) != 23)
      continue;
    else if(sample == Sample::wjet and fabs(*wzid) != 24)
      continue;
    sumwgt += *wgt;
  }    
  sumwgt_vec.push_back(sumwgt);
}

//////////
void fillHistograms(TChain* chain, vector<TH1F*> & histogram, vector<TH1*> & khists, const Sample & sample, const Category & category, vector<double> & sumwgt){
  
  TTreeReader reader (chain);
  TTreeReaderValue<float> xsec  (reader,"xsec");
  TTreeReaderValue<float> wgt   (reader,"wgt");
  TTreeReaderValue<int>   wzid  (reader,"wzid");
  TTreeReaderValue<int>   l1id  (reader,"l1id");
  TTreeReaderValue<int>   l2id  (reader,"l2id");
  TTreeReaderValue<unsigned int>   njetsinc  (reader,"njetsinc");
  TTreeReaderValue<unsigned int>   njets  (reader,"njets");
  TTreeReaderValue<float> wzmass  (reader,"wzmass");
  TTreeReaderValue<float> wzpt    (reader,"wzpt");
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
  
  ///////////////////////////////
  int ifile = 0;
  string currentFile = "";
  while(reader.Next()){

    //check if file name switched or not                                                                                                                                                               
    if(dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName() != currentFile and currentFile != ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      ifile++;
      cout<<currentFile<<" "<<ifile<<endl;
    }
    else if(currentFile == ""){
      currentFile = dynamic_cast<TChain*>(reader.GetTree())->GetFile()->GetName();
      ifile = 0;
      cout<<currentFile<<" "<<ifile<<endl;
    }

    // filter away bad events with no matching                                                                                                                                                        
    if(sample == Sample::zjet and fabs(*wzid) != 23)
      continue;
    else if(sample == Sample::wjet and fabs(*wzid) != 24)
      continue;
    
    //////////////////
    if(*wzpt < 250) continue; // --> V-boson cut
    if(sample == Sample::wjet and (fabs(*l1id) == 15 or fabs(*l1id) == 16 or fabs(*l2id) == 15 or fabs(*l2id) == 16)) continue; // skip taus 

    Double_t kwgt = 1.0;
    double genpt  = *wzpt;
    for (size_t i = 0; i < khists.size(); i++) {
      if (khists[i]) {
	if(genpt <= khists[i]->GetXaxis()->GetBinLowEdge(1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(1) + 1;
	if(genpt >= khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)) genpt = khists[i]->GetXaxis()->GetBinLowEdge(khists[i]->GetNbinsX()+1)-1;
	kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(genpt));
      }
    }

    float mindphi = 100;
    for(size_t ijet = 0; ijet < jetphi->size(); ijet++){
      if(ijet > 3) break; // limiting min dphi to first 4 leading jets                                                                                                                                
      float dphi = fabs(*wzphi-jetphi->at(ijet));
      if(dphi > TMath::Pi())
	dphi = 2*TMath::Pi()-dphi;
      if(dphi < mindphi)
	mindphi = dphi;
    }

    TLorentzVector jet1, jet2;
    ///////////////
    if(category == Category::VBFrelaxed){
      if(*njetsinc < 2) continue;
      if(jetpt->size() < 2) continue;
      if(fabs(jeteta->at(0)) > 4.7) continue;
      if(fabs(jeteta->at(1)) > 4.7) continue;
      if(fabs(jetpt->at(0)) < 80) continue;
      if(fabs(jetpt->at(1)) < 40) continue;
      if(fabs(jeteta->at(0)-jeteta->at(1)) < 1.0) continue;
      if(jeteta->at(0)*jeteta->at(1) > 0) continue;
      if(fabs(mindphi) < 0.5) continue;
      jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetmass->at(0));
      jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetmass->at(1));          
      if((jet1+jet2).M() < 200) continue;
      if(fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) > 3) continue;
    }
    else if(category == Category::VBF){
      if(*njetsinc < 2) continue;
      if(jetpt->size() < 2) continue;
      if(fabs(jeteta->at(0)) > 4.7) continue;
      if(fabs(jeteta->at(1)) > 4.7) continue;
      if(fabs(jetpt->at(0)) < 80) continue;
      if(fabs(jetpt->at(1)) < 40) continue;
      if(fabs(jeteta->at(0)-jeteta->at(1)) < 4.0) continue;
      if(jeteta->at(0)*jeteta->at(1)) continue;
      if(fabs(mindphi) < 0.5) continue;
      jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetmass->at(0));
      jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetmass->at(1));          
      if((jet1+jet2).M() < 1300) continue;
      if(fabs(jet1.Eta()) > 3 and fabs(jet2.Eta()) > 3) continue;
    }

    for(auto histo : histogram){
      double fillvar = 0.;
      TString name (histo->GetName());
      if(name.Contains("bosonpt"))
	fillvar = *wzpt;
      else if(name.Contains("jetpt1"))
	fillvar = jetpt->at(0);
      else if(name.Contains("jetpt2"))
	fillvar = jetpt->at(1);
      else if(name.Contains("detajj"))
	fillvar = fabs(jeteta->at(0)-jeteta->at(1));
      else if(name.Contains("mjj"))
	fillvar = (jet1+jet2).M();
      else if(name.Contains("dphijj"))
	fillvar = fabs(jet1.DeltaPhi(jet2));      	

      histo->Fill(fillvar,*xsec*luminosity*(*wgt)*kwgt/(sumwgt.at(ifile)));

    }
  }              
}

////////
void drawHistogram(TH1F* histogram_lo, TH1F* histogram_nlo, TH1F* histogram_reweight, const Sample & sample, const string & outputDIR){

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);

  // normalize to a.u.
  histogram_lo->Scale(1./histogram_lo->Integral());
  histogram_nlo->Scale(1./histogram_nlo->Integral());
  histogram_reweight->Scale(1./histogram_reweight->Integral());

  //
  TH1* ratio_1 =  (TH1*) histogram_nlo->Clone("ratio_1");
  TH1* ratio_2 =  (TH1*) histogram_reweight->Clone("ratio_2");
  ratio_1->Divide(histogram_lo);
  ratio_2->Divide(histogram_lo);
  histogram_lo->GetXaxis()->SetLabelSize(0);
  histogram_lo->GetXaxis()->SetTitleSize(0);
  
  TString name (histogram_lo->GetName());
  string postfix = "";
  if(name.Contains("bosonpt")){
    ratio_1->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
    postfix = "bosonpt";
  }
  else if(name.Contains("jetpt1")){
    ratio_1->GetXaxis()->SetTitle("p_{T}^{j1} [GeV]");
    postfix = "jetpt1";
  }
  else if(name.Contains("jetpt2")){
    ratio_1->GetXaxis()->SetTitle("p_{T}^{j2} [GeV]");
    postfix = "jetpt2";
  }
  else if(name.Contains("detajj")){
    ratio_1->GetXaxis()->SetTitle("#Delta#eta_{jj}");
    postfix = "detajj";
  }
  else if(name.Contains("dphijj")){
    ratio_1->GetXaxis()->SetTitle("#Delta#phi_{jj}");
    postfix = "dphijj";
  }
  else if(name.Contains("mjj")){
    ratio_1->GetXaxis()->SetTitle("M_{jj} [GeV]");
    postfix = "mjj";
  }
  
  if(name.Contains("detajj") or name.Contains("dphijj"))
    histogram_lo->GetYaxis()->SetRangeUser(0,max(histogram_lo->GetMaximum(),max(histogram_nlo->GetMaximum(),histogram_reweight->GetMaximum()))*1.5);
  else
    histogram_lo->GetYaxis()->SetRangeUser(0.0001,max(histogram_lo->GetMaximum(),max(histogram_nlo->GetMaximum(),histogram_reweight->GetMaximum()))*100);
  
  histogram_lo->GetYaxis()->SetTitle("a.u.");
  histogram_lo->GetYaxis()->SetTitleOffset(1.1);
  histogram_lo->SetLineColor(kBlack);
  histogram_lo->SetLineWidth(2);
  histogram_nlo->SetLineColor(kRed);
  histogram_nlo->SetLineWidth(2);
  histogram_reweight->SetLineColor(kBlue);
  histogram_reweight->SetLineWidth(2);

  histogram_lo->Draw("hist");
  histogram_nlo->Draw("hist same");
  histogram_reweight->Draw("hist same");
  
  ///
  CMS_lumi(canvas,"");

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  if(sample == Sample::zjet){
    postfix += "_zjet";
    leg->AddEntry(histogram_lo,"Z+jets LO","L");
    leg->AddEntry(histogram_nlo,"Z+jets NLO","L");
    leg->AddEntry(histogram_reweight,"Z+jets Re-weight","L");
  }      
  else if(sample == Sample::wjet){
    postfix += "_wjet";
    leg->AddEntry(histogram_lo,"W+jets LO","L");
    leg->AddEntry(histogram_nlo,"W+jets NLO","L");
    leg->AddEntry(histogram_reweight,"W+jets Re-weight","L");
  }
  leg->Draw("same");

  if(name.Contains("detajj") or name.Contains("dphijj"))
    canvas->SetLogy(0);
  else
    canvas->SetLogy();

  pad2->Draw();
  pad2->cd();
  ratio_1->GetYaxis()->SetTitle("NLO/LO");
  ratio_1->GetYaxis()->SetTitleOffset(1.20);
  ratio_1->GetYaxis()->SetTitleSize(0.04);
  ratio_1->GetYaxis()->SetLabelSize(0.03);
  ratio_1->GetYaxis()->SetNdivisions(505);
  ratio_1->GetXaxis()->SetTitleOffset(1.10);
  ratio_1->SetLineColor(kRed);
  ratio_1->SetLineWidth(2);
  ratio_2->SetLineColor(kBlue);
  ratio_2->SetLineWidth(2);
  ratio_1->Draw("hist");
  ratio_2->Draw("hist same");

  canvas->SaveAs((outputDIR+"/"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/"+postfix+".root").c_str(),"root");

  if(canvas) delete canvas;
}

/////////////
void makeVJetsShapeComparisons(string inputDIR_LO, string inputDIR_NLO, string kfactorFile, Sample sample, Category category, string outputDIR){

  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  TFile* kffile = TFile::Open(kfactorFile.c_str());
  TH1* hist_nloqcd    = NULL;
  TH1* hist_loqcd     = NULL;

  if(sample == Sample::zjet){
    hist_nloqcd    = (TH1*) kffile->Get("ZJets_012j_NLO/nominal");
    hist_loqcd     = (TH1*) kffile->Get("ZJets_LO/inv_pt");
  }
  else if(sample == Sample::wjet){
    hist_nloqcd    = (TH1*) kffile->Get("WJets_012j_NLO/nominal");
    hist_loqcd     = (TH1*) kffile->Get("WJets_LO/inv_pt");
  }
 
  hist_nloqcd->Divide(hist_loqcd);

  vector<TH1*> khists; khists.push_back(hist_nloqcd);
  vector<TH1*> ehists;

  TChain* chain_lo = new TChain("gentree/tree");
  TChain* chain_nlo = new TChain("gentree/tree");
  
  chain_lo->Add((inputDIR_LO+"/*root").c_str());
  chain_nlo->Add((inputDIR_NLO+"/*root").c_str());

  vector<TH1F*> histogram_lo;
  vector<TH1F*> histogram_nlo;
  vector<TH1F*> histogram_reweight;

  histogram_lo.push_back(new TH1F("histogram_lo_bosonpt","",25,200,1200));
  histogram_lo.push_back(new TH1F("histogram_lo_jetpt1","",25,80,700));
  histogram_lo.push_back(new TH1F("histogram_lo_jetpt2","",25,40,400));
  histogram_lo.push_back(new TH1F("histogram_lo_dphijj","",25,0,3.14));
  histogram_lo.push_back(new TH1F("histogram_lo_detajj","",20,1,7));
  histogram_lo.push_back(new TH1F("histogram_lo_mjj","",25,200,3500));
  
  for(auto hist : histogram_lo)
    hist->Sumw2();

  histogram_nlo.push_back(new TH1F("histogram_nlo_bosonpt","",25,200,1200));
  histogram_nlo.push_back(new TH1F("histogram_nlo_jetpt1","",25,80,700));
  histogram_nlo.push_back(new TH1F("histogram_nlo_jetpt2","",25,40,400));
  histogram_nlo.push_back(new TH1F("histogram_nlo_dphijj","",25,0,3.14));
  histogram_nlo.push_back(new TH1F("histogram_nlo_detajj","",20,1,7));
  histogram_nlo.push_back(new TH1F("histogram_nlo_mjj","",25,200,3500));

  for(auto hist : histogram_nlo)
    hist->Sumw2();

  histogram_reweight.push_back(new TH1F("histogram_reweight_bosonpt","",25,200,1200));
  histogram_reweight.push_back(new TH1F("histogram_reweight_jetpt1","",25,80,700));
  histogram_reweight.push_back(new TH1F("histogram_reweight_jetpt2","",25,40,400));
  histogram_reweight.push_back(new TH1F("histogram_reweight_dphijj","",25,0,3.14));
  histogram_reweight.push_back(new TH1F("histogram_reweight_detajj","",20,1,7));
  histogram_reweight.push_back(new TH1F("histogram_reweight_mjj","",25,200,3500));

  for(auto hist : histogram_reweight)
    hist->Sumw2();

  vector<double> sumwgt_lo;
  vector<double> sumwgt_nlo;
  cout<<"Calculate sumwgt for LO samples "<<endl;
  calculateSumWg(chain_lo,sumwgt_lo,sample);
  cout<<"Calculate sumwgt for NLO samples "<<endl;
  calculateSumWg(chain_nlo,sumwgt_nlo,sample);
  
  fillHistograms(chain_lo, histogram_lo, ehists, sample, category, sumwgt_lo);
  fillHistograms(chain_nlo, histogram_nlo, ehists, sample, category, sumwgt_nlo);
  fillHistograms(chain_lo, histogram_reweight, khists, sample, category, sumwgt_lo);

  for(size_t ihist = 0; ihist < histogram_lo.size(); ihist++)
    drawHistogram(histogram_lo.at(ihist),histogram_nlo.at(ihist),histogram_reweight.at(ihist),sample,outputDIR);

  string postfix;
  if(sample == Sample::zjet)
    postfix = "_zjet";
  else if(sample == Sample::wjet)
    postfix = "_wjet";
    

  TFile* outputFile = new TFile((outputDIR+"/outputFile_"+postfix+".root").c_str(),"RECREATE");
  outputFile->cd();
  for(auto hist : histogram_lo)
    hist->Write();
  for(auto hist : histogram_nlo)
    hist->Write();
  for(auto hist : histogram_reweight)
    hist->Write();
  
  outputFile->Close();

  
}
