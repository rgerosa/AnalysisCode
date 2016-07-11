// pt binning to study ID observables
#include "../CMS_lumi.h"

static vector<float> ptBin  = {10.,30.,50,100,150,500};
static vector<float> etaBin = {-2.4,-1.2,0,1.2,2.4};
const int nvtxMin = 0;
const int nvtxMax = 50;
const float luminosity = 4.30;

vector<TH1F*> muonChi2_data;
vector<TH1F*> muonGlobalId_data;
vector<TH1F*> muonPfId_data;
vector<TH1F*> muonNvalidHit_data;
vector<TH1F*> muonNpixelHit_data;
vector<TH1F*> muonNtrackerLayer_data;
vector<TH1F*> muonNStation_data;
vector<TH1F*> muonDxy_data;
vector<TH1F*> muonDz_data;
vector<TH1F*> muonNhiso_data;
vector<TH1F*> muonChiso_data;
vector<TH1F*> muonPuiso_data;
vector<TH1F*> muonEMIso_data;
vector<TH1F*> muonZMass_data;
vector<TH1F*> muonPt_data;
vector<TH1F*> muonEta_data;

vector<TH1F*> muonChi2_mc;
vector<TH1F*> muonPfId_mc;
vector<TH1F*> muonGlobalId_mc;
vector<TH1F*> muonNvalidHit_mc;
vector<TH1F*> muonNpixelHit_mc;
vector<TH1F*> muonNStation_mc;
vector<TH1F*> muonNtrackerLayer_mc;
vector<TH1F*> muonDxy_mc;
vector<TH1F*> muonDz_mc;
vector<TH1F*> muonNhiso_mc;
vector<TH1F*> muonChiso_mc;
vector<TH1F*> muonPuiso_mc;
vector<TH1F*> muonEMIso_mc;
vector<TH1F*> muonZMass_mc;
vector<TH1F*> muonPt_mc;
vector<TH1F*> muonEta_mc;

vector<TH1F*> muonEfficiency_num_data;
vector<TH1F*> muonEfficiency_num_mc;
vector<TH1F*> muonEfficiency_den_data;
vector<TH1F*> muonEfficiency_den_mc;

void plotEfficiency(TH1* histo_num_data, TH1* histo_den_data, TH1* histo_num_mc, TH1* histo_den_mc, TCanvas* canvas, const string & plotDIR){

  canvas->cd();
  canvas->Clear();

  TGraphAsymmErrors* eff_data = new TGraphAsymmErrors();
  TGraphAsymmErrors* eff_mc = new TGraphAsymmErrors();
  
  histo_num_data->GetYaxis()->SetRangeUser(0.85,1.10);
  histo_num_data->Draw("AXIS");

  eff_data->BayesDivide(histo_num_data,histo_den_data);
  eff_mc->BayesDivide(histo_num_mc,histo_den_mc);

  eff_data->SetLineColor(kBlack);
  eff_data->SetMarkerColor(kBlack);
  eff_data->SetMarkerStyle(20);
  eff_data->SetMarkerSize(1);


  eff_mc->SetLineColor(kRed);
  eff_mc->SetMarkerColor(kRed);
  eff_mc->SetFillColor(kRed);
  eff_mc->SetFillStyle(3001);
  eff_mc->SetMarkerStyle(20);
  eff_mc->SetMarkerSize(1);

  eff_data->Draw("P0same");
  eff_mc->Draw("P0same");

  CMS_lumi(canvas,Form("%.2f",luminosity),true);
  
  TLegend* leg = new TLegend(0.7,0.78,0.9,0.92);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(eff_data,"Data","PE");
  leg->AddEntry(eff_mc,"MC","PE");
  leg->Draw("same");

  //  canvas->SaveAs((plotDIR+"/"+string(histo_num_data->GetName())+".png").c_str(),"png");
  canvas->SaveAs((plotDIR+"/"+string(histo_num_data->GetName())+".pdf").c_str(),"pdf");

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
  histo_data->GetYaxis()->SetTitle("A.U.");
  histo_data->SetLineColor(kBlack);
  histo_data->SetMarkerColor(kBlack);
  histo_data->SetMarkerStyle(20);
  histo_data->SetMarkerSize(1);
  histo_data->GetYaxis()->SetRangeUser(1,max(histo_data->GetMaximum(),histo_mc->GetMaximum())*1000);

  histo_data->Scale(1./histo_data->Integral());
  histo_mc->Scale(1./histo_mc->Integral());

  histo_mc->SetFillColor(kCyan+1);
  histo_mc->SetFillStyle(3001);
  histo_mc->SetLineColor(kRed);
  histo_mc->SetLineWidth(2);
  histo_data->Draw("EP");
  histo_mc->Draw("hist same");
  histo_data->Draw("EPsame");
  CMS_lumi(pad1,Form("%.2f",luminosity),true);
  TLegend* leg = new TLegend(0.75,0.75,0.9,0.9);
  leg->SetFillColor(kWhite);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histo_data,"Data","PE");
  leg->AddEntry(histo_mc,"Z+jets","F");
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
  ratio->GetYaxis()->SetRangeUser(0.5,1.5);
  ratio->GetYaxis()->SetTitleSize(ratio->GetYaxis()->GetTitleSize()*1.5);
  ratio->GetYaxis()->SetLabelSize(ratio->GetYaxis()->GetLabelSize()*1.5);
  ratio->Draw("EP0");
  TF1* line = new TF1("line","1",ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1));
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  canvas->cd();		      
  //  canvas->SaveAs((plotDIR+"/"+string(histo_data->GetName())+".png").c_str(),"png");
  canvas->SaveAs((plotDIR+"/"+string(histo_data->GetName())+".pdf").c_str(),"pdf");
}

void fillHistograms(TTree* tree, const float & ptMin, const float & ptMax, const float & etaMin, const float & etaMax, const size_t & ivec, const bool & isData, const double & wgtsum = 1){

  TTreeReader reader(tree);
  TTreeReaderValue<float> pt   (reader,"pt");
  TTreeReaderValue<float> eta  (reader,"eta");
  TTreeReaderValue<float> mass (reader,"mass");
  TTreeReaderValue<float> nvtx (reader,"nvtx");
  TTreeReaderValue<float> wgt  (reader,"wgt");
  TTreeReaderValue<float> nstation  (reader,"nstation");
  TTreeReaderValue<float> chi2  (reader,"chi2");
  TTreeReaderValue<float> nvalidhit (reader,"nvalidhit");
  TTreeReaderValue<float> npixelhit (reader,"npixelhit");
  TTreeReaderValue<float> ntrackerlayerhit (reader,"ntrackerlayerhit");
  TTreeReaderValue<float> dxy (reader,"dxy");
  TTreeReaderValue<float> dz (reader,"dz");
  TTreeReaderValue<float> nhiso (reader,"nhiso");
  TTreeReaderValue<float> emiso (reader,"emiso");
  TTreeReaderValue<float> puiso (reader,"puiso");
  TTreeReaderValue<float> chiso (reader,"chiso");
  TTreeReaderValue<int> pfid (reader,"pfid");
  TTreeReaderValue<int> globalid (reader,"globalid");

  TTreeReaderValue<int> hltmu (reader,"hltmu");
  TTreeReaderValue<int> hlttkmu (reader,"hlttkmu");
    
  TFile pufile ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_4p0fb.root");
  TH1* puhist = (TH1*) pufile.Get("puhist");  

  TFile triggerfile_SinglMu("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/triggerEfficiency_DATA_SingleMuon.root");
  TH2F* triggermu_eff = (TH2F*) triggerfile_SinglMu.Get("trigeff_muIso");

  float totalMuons = 0;

  while(reader.Next()){    
    //basic selections for trigger, filters, bjet, taus and lepton vetoes    
    unsigned int hlt = *hltmu + *hlttkmu;
    // trigger for data
    if(hlt == 0) continue;
    if(*pt < ptMin or *pt > ptMax) continue;
    if(*eta < etaMin or *eta > etaMax) continue;

    // to ensure a background free sample
    // super loose isolation cut
    if((*chiso+max(0.,*nhiso+*emiso -0.5*(*puiso)))/(*pt) > 0.35) continue;

    // in case of data
    if(isData){            
      // to ensure a background free sample
      if(*globalid == 1){
	// no need for N-1 efficiency
	muonChi2_data.at(ivec)->Fill(*chi2);
	muonPfId_data.at(ivec)->Fill(*pfid);
	muonNvalidHit_data.at(ivec)->Fill(*nvalidhit);
	muonNpixelHit_data.at(ivec)->Fill(*npixelhit);
	muonNtrackerLayer_data.at(ivec)->Fill(*ntrackerlayerhit);
	muonNStation_data.at(ivec)->Fill(*nstation);
	muonDxy_data.at(ivec)->Fill(*dxy);
	muonDz_data.at(ivec)->Fill(*dz);
	muonNhiso_data.at(ivec)->Fill(*nhiso);
	muonChiso_data.at(ivec)->Fill(*chiso);
	muonEMIso_data.at(ivec)->Fill(*puiso);
	muonPuiso_data.at(ivec)->Fill(*emiso);
	muonZMass_data.at(ivec)->Fill(*mass);
      }
      if((*chiso+max(0.,*nhiso+*emiso -0.5*(*puiso)))/(*pt) < 0.25){	
	muonGlobalId_data.at(ivec)->Fill(*globalid);
      }

      // efficiency studies
      for(int iBin = 0; iBin < muonEfficiency_den_data.at(ivec)->GetNbinsX()+1; iBin++)
	muonEfficiency_den_data.at(ivec)->SetBinContent(iBin+1,muonEfficiency_den_data.at(ivec)->GetBinContent(iBin+1)+1);
      if(*globalid == 1)
	muonEfficiency_num_data.at(ivec)->SetBinContent(1,muonEfficiency_num_data.at(ivec)->GetBinContent(1)+1);
      if(*globalid == 1 and *pfid == 1)
	muonEfficiency_num_data.at(ivec)->SetBinContent(2,muonEfficiency_num_data.at(ivec)->GetBinContent(2)+1);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10)
	muonEfficiency_num_data.at(ivec)->SetBinContent(3,muonEfficiency_num_data.at(ivec)->GetBinContent(3)+1);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0)
	muonEfficiency_num_data.at(ivec)->SetBinContent(4,muonEfficiency_num_data.at(ivec)->GetBinContent(4)+1);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0 and *nstation > 1)
	muonEfficiency_num_data.at(ivec)->SetBinContent(5,muonEfficiency_num_data.at(ivec)->GetBinContent(5)+1);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0 and *nstation > 1 and fabs(*dxy) < 0.2)
	muonEfficiency_num_data.at(ivec)->SetBinContent(6,muonEfficiency_num_data.at(ivec)->GetBinContent(6)+1);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0 and *nstation > 1 and fabs(*dz) < 0.5)
	muonEfficiency_num_data.at(ivec)->SetBinContent(7,muonEfficiency_num_data.at(ivec)->GetBinContent(7)+1);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0 and *nstation > 1 and fabs(*dz) < 0.5 and *npixelhit > 0)
	muonEfficiency_num_data.at(ivec)->SetBinContent(8,muonEfficiency_num_data.at(ivec)->GetBinContent(8)+1);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0 and *nstation > 1 and fabs(*dz) < 0.5 and *npixelhit > 0 and *ntrackerlayerhit > 5)
	muonEfficiency_num_data.at(ivec)->SetBinContent(9,muonEfficiency_num_data.at(ivec)->GetBinContent(9)+1);
      float iso = (*chiso + max(0.,*nhiso+*emiso-0.5*(*puiso)))/(*pt);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0 and *nstation > 1 and fabs(*dz) < 0.5 and *npixelhit > 0 and *ntrackerlayerhit > 5 and iso < 0.15)
	muonEfficiency_num_data.at(ivec)->SetBinContent(10,muonEfficiency_num_data.at(ivec)->GetBinContent(10)+1);

    }
    else{
      
      double pwgt  = 1;
      pwgt *= triggermu_eff->GetBinContent(triggermu_eff->FindBin(*pt,*eta));
      if(*nvtx < 50)
      	pwgt *= puhist->GetBinContent(puhist->FindBin(*nvtx));
      // we assume that NLO DYJets simulation is considered, otherwise kfactors
      pwgt *= (*wgt)/(wgtsum);
      // no need for N-1 efficiency
      if(*globalid == 1){
	muonChi2_mc.at(ivec)->Fill(*chi2,pwgt);
	muonPfId_mc.at(ivec)->Fill(*pfid,pwgt);
	muonNvalidHit_mc.at(ivec)->Fill(*nvalidhit,pwgt);
	muonNpixelHit_mc.at(ivec)->Fill(*npixelhit,pwgt);
	muonNtrackerLayer_mc.at(ivec)->Fill(*ntrackerlayerhit,pwgt);
	muonNStation_mc.at(ivec)->Fill(*nstation,pwgt);
	muonDxy_mc.at(ivec)->Fill(*dxy,pwgt);
	muonDz_mc.at(ivec)->Fill(*dz,pwgt);
	muonNhiso_mc.at(ivec)->Fill(*nhiso,pwgt);
	muonChiso_mc.at(ivec)->Fill(*chiso,pwgt);
	muonEMIso_mc.at(ivec)->Fill(*puiso,pwgt);
	muonPuiso_mc.at(ivec)->Fill(*emiso,pwgt);
	muonZMass_mc.at(ivec)->Fill(*mass,pwgt);
      }
      if((*chiso+max(0.,*nhiso+*emiso -0.5*(*puiso)))/(*pt) < 0.25){	
      	muonGlobalId_mc.at(ivec)->Fill(*globalid,pwgt);
      }

      for(int iBin = 0; iBin < muonEfficiency_den_mc.at(ivec)->GetNbinsX()+1; iBin++)
	muonEfficiency_den_mc.at(ivec)->SetBinContent(iBin+1,muonEfficiency_den_mc.at(ivec)->GetBinContent(iBin+1)+pwgt);
      if(*globalid == 1)
	muonEfficiency_num_mc.at(ivec)->SetBinContent(1,muonEfficiency_num_mc.at(ivec)->GetBinContent(1)+pwgt);
      if(*globalid == 1 and *pfid == 1)
	muonEfficiency_num_mc.at(ivec)->SetBinContent(2,muonEfficiency_num_mc.at(ivec)->GetBinContent(2)+pwgt);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10)
	muonEfficiency_num_mc.at(ivec)->SetBinContent(3,muonEfficiency_num_mc.at(ivec)->GetBinContent(3)+pwgt);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0)
	muonEfficiency_num_mc.at(ivec)->SetBinContent(4,muonEfficiency_num_mc.at(ivec)->GetBinContent(4)+pwgt);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0 and *nstation > 1)
	muonEfficiency_num_mc.at(ivec)->SetBinContent(5,muonEfficiency_num_mc.at(ivec)->GetBinContent(5)+pwgt);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0 and *nstation > 1 and fabs(*dxy) < 0.2)
	muonEfficiency_num_mc.at(ivec)->SetBinContent(6,muonEfficiency_num_mc.at(ivec)->GetBinContent(6)+pwgt);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0 and *nstation > 1 and fabs(*dz) < 0.5)
	muonEfficiency_num_mc.at(ivec)->SetBinContent(7,muonEfficiency_num_mc.at(ivec)->GetBinContent(7)+pwgt);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0 and *nstation > 1 and fabs(*dz) < 0.5 and *npixelhit > 0)
	muonEfficiency_num_mc.at(ivec)->SetBinContent(8,muonEfficiency_num_mc.at(ivec)->GetBinContent(8)+pwgt);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0 and *nstation > 1 and fabs(*dz) < 0.5 and *npixelhit > 0 and *ntrackerlayerhit > 5)
	muonEfficiency_num_mc.at(ivec)->SetBinContent(9,muonEfficiency_num_mc.at(ivec)->GetBinContent(9)+pwgt);
      float iso = (*chiso + max(0.,*nhiso+*emiso-0.5*(*puiso)))/(*pt);
      if(*globalid == 1 and *pfid == 1 and *chi2 < 10 and *nvalidhit > 0 and *nstation > 1 and fabs(*dz) < 0.5 and *npixelhit > 0 and *ntrackerlayerhit > 5 and iso < 0.15)
	muonEfficiency_num_mc.at(ivec)->SetBinContent(10,muonEfficiency_num_mc.at(ivec)->GetBinContent(10)+pwgt);
    }
  }    
  triggerfile_SinglMu.Close();
  pufile.Close();
}




void fillPtHistograms(TTree* tree, const float & etaMin, const float & etaMax, const size_t & ivec, const bool & isData, const double & wgtsum = 1){

  TTreeReader reader(tree);
  TTreeReaderValue<float> pt   (reader,"pt");
  TTreeReaderValue<float> eta  (reader,"eta");
  TTreeReaderValue<float> nvtx (reader,"nvtx");
  TTreeReaderValue<float> wgt  (reader,"wgt");
  TTreeReaderValue<int> globalid (reader,"globalid");
  TTreeReaderValue<int> hltmu (reader,"hltmu");
  TTreeReaderValue<int> hlttkmu (reader,"hlttkmu");
  TTreeReaderValue<int> tightid (reader,"tightid");
  TTreeReaderValue<float> nhiso (reader,"nhiso");
  TTreeReaderValue<float> emiso (reader,"emiso");
  TTreeReaderValue<float> puiso (reader,"puiso");
  TTreeReaderValue<float> chiso (reader,"chiso");
    
  TFile pufile ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_4p0fb.root");
  TH1* puhist = (TH1*) pufile.Get("puhist");  

  TFile triggerfile_SinglMu("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/triggerEfficiency_DATA_SingleMuon.root");
  TH2F* triggermu_eff = (TH2F*) triggerfile_SinglMu.Get("trigeff_muIso");

  while(reader.Next()){    
    //basic selections for trigger, filters, bjet, taus and lepton vetoes    
    unsigned int hlt = *hltmu + *hlttkmu;
    // trigger for data
    if(hlt == 0) continue;
    if(*eta < etaMin or *eta > etaMax) continue;

    // to ensure a background free sample
    // super loose isolation cut
    if((*chiso+max(0.,*nhiso+*emiso -0.5*(*puiso)))/(*pt) > 0.35) continue;
    
    // in case of data
    if(isData){            
      // to ensure a background free sample
      if(*globalid == 1 and *tightid == 1)
	// no need for N-1 efficiency
	muonPt_data.at(ivec)->Fill(*pt);
    }
    else{
      
      double pwgt  = 1;
      pwgt *= triggermu_eff->GetBinContent(triggermu_eff->FindBin(*pt,*eta));
      if(*nvtx < 50)
      	pwgt *= puhist->GetBinContent(puhist->FindBin(*nvtx));
      // we assume that NLO DYJets simulation is considered, otherwise kfactors
      pwgt *= (*wgt)/(wgtsum);
      // no need for N-1 efficiency
      if(*globalid == 1 and *tightid == 1)
	muonPt_mc.at(ivec)->Fill(*pt,pwgt);
    }
  }    
  triggerfile_SinglMu.Close();
  pufile.Close();
}

void fillEtaHistograms(TTree* tree, const float & ptMin, const float & ptMax, const size_t & ivec, const bool & isData, const double & wgtsum = 1){

  TTreeReader reader(tree);
  TTreeReaderValue<float> pt   (reader,"pt");
  TTreeReaderValue<float> eta  (reader,"eta");
  TTreeReaderValue<float> nvtx (reader,"nvtx");
  TTreeReaderValue<float> wgt  (reader,"wgt");
  TTreeReaderValue<int> globalid (reader,"globalid");
  TTreeReaderValue<int> hltmu (reader,"hltmu");
  TTreeReaderValue<int> hlttkmu (reader,"hlttkmu");
  TTreeReaderValue<int> tightid (reader,"tightid");
  TTreeReaderValue<float> nhiso (reader,"nhiso");
  TTreeReaderValue<float> emiso (reader,"emiso");
  TTreeReaderValue<float> puiso (reader,"puiso");
  TTreeReaderValue<float> chiso (reader,"chiso");
    
  TFile pufile ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/puwrt_4p0fb.root");
  TH1* puhist = (TH1*) pufile.Get("puhist");  

  TFile triggerfile_SinglMu("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/triggerEfficiency_DATA_SingleMuon.root");
  TH2F* triggermu_eff = (TH2F*) triggerfile_SinglMu.Get("trigeff_muIso");

  while(reader.Next()){    
    //basic selections for trigger, filters, bjet, taus and lepton vetoes    
    unsigned int hlt = *hltmu + *hlttkmu;
    // trigger for data
    if(hlt == 0) continue;
    if(*pt < ptMin or *pt > ptMax) continue;

    // to ensure a background free sample
    // super loose isolation cut
    if((*chiso+max(0.,*nhiso+*emiso -0.5*(*puiso)))/(*pt) > 0.35) continue;
    
    // in case of data
    if(isData){            
      // to ensure a background free sample
      if(*globalid == 1 and *tightid == 1)
	// no need for N-1 efficiency
	muonEta_data.at(ivec)->Fill(*eta);
    }
    else{
      
      double pwgt  = 1;
      pwgt *= triggermu_eff->GetBinContent(triggermu_eff->FindBin(*pt,*eta));
      if(*nvtx < 50)
      	pwgt *= puhist->GetBinContent(puhist->FindBin(*nvtx));
      // we assume that NLO DYJets simulation is considered, otherwise kfactors
      pwgt *= (*wgt)/(wgtsum);
      // no need for N-1 efficiency
      if(*globalid == 1 and *tightid == 1)
	muonEta_mc.at(ivec)->Fill(*eta,pwgt);
    }
  }    
  triggerfile_SinglMu.Close();
  pufile.Close();
}



void makeMuonIDVariableComparison(string inputDIRData, string inputDIRMC, string ouputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+ouputDIR).c_str());
  cout<<"Start chain files"<<endl;
  TChain* chain_data   = new TChain("muontnptree/fitter_tree");
  chain_data->Add((inputDIRData+"/*root").c_str());
  TChain* chain_mc     = new TChain("muontnptree/fitter_tree");
  chain_mc->Add((inputDIRMC+"/*root").c_str());
  cout<<"End to chain files "<<endl;
  
  cout<<"Calculate mc weight "<<endl;
  TTreeReader treeReader_gen(chain_mc);
  TTreeReaderValue<float> wgt (treeReader_gen,"wgt");
  double weightsum = 0.;
  while(treeReader_gen.Next()){
    weightsum += (*wgt);
  }
  cout<<"Create histogram "<<endl;

  // create histograms
  for(size_t ieta = 0; ieta < etaBin.size()-1; ieta++){
    for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++){

      muonChi2_data.push_back(new TH1F(Form("trackchi2_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					    int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",50,0,20));
      muonChi2_data.back()->Sumw2();
      muonPfId_data.push_back(new TH1F(Form("pfid_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					    int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",2,0,2));
      muonPfId_data.back()->Sumw2();
      muonGlobalId_data.push_back(new TH1F(Form("globalid_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",2,0,2));
      muonGlobalId_data.back()->Sumw2();
      muonNvalidHit_data.push_back(new TH1F(Form("nvalidhit_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						 int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",15,0,15));
      muonNvalidHit_data.back()->Sumw2();
      muonNpixelHit_data.push_back(new TH1F(Form("npixelhit_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						 int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",10,0,10));
      muonNpixelHit_data.back()->Sumw2();
      muonNtrackerLayer_data.push_back(new TH1F(Form("ntrackerlayer_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						     int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",30,0,30));
      muonNtrackerLayer_data.back()->Sumw2();
      muonNStation_data.push_back(new TH1F(Form("nstation_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",8,0,8));
      muonNStation_data.back()->Sumw2();
      muonDxy_data.push_back(new TH1F(Form("dxy_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					   int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",40,-0.25,0.25));
      muonDxy_data.back()->Sumw2();
      muonDz_data.push_back(new TH1F(Form("dz_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					  int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",40,-0.5,0.5));
      muonDz_data.back()->Sumw2();
      muonNhiso_data.push_back(new TH1F(Form("nhiso_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					     int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",40,0,15));
      muonNhiso_data.back()->Sumw2();
      muonChiso_data.push_back(new TH1F(Form("chiso_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",					     
					     int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",40,0,15));					
      muonChiso_data.back()->Sumw2();
      muonPuiso_data.push_back(new TH1F(Form("puiso_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					     int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",40,0,15));
      muonPuiso_data.back()->Sumw2();
      muonEMIso_data.push_back(new TH1F(Form("emiso_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					     int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",40,0,25));
      muonEMIso_data.back()->Sumw2();
      muonZMass_data.push_back(new TH1F(Form("zmass_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					     int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",25,60,120));
      muonZMass_data.back()->Sumw2();


      muonChi2_mc.push_back(new TH1F(Form("trackchi2_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					    int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",50,0,20));
      muonChi2_mc.back()->Sumw2();
      muonPfId_mc.push_back(new TH1F(Form("pfid_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					    int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",2,0,2));
      muonPfId_mc.back()->Sumw2();
      muonGlobalId_mc.push_back(new TH1F(Form("globalid_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",2,0,2));
      muonGlobalId_mc.back()->Sumw2();
      muonNvalidHit_mc.push_back(new TH1F(Form("nvalidhit_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						 int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",15,0,15));
      muonNvalidHit_mc.back()->Sumw2();
      muonNpixelHit_mc.push_back(new TH1F(Form("npixelhit_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						 int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",10,0,10));
      muonNpixelHit_mc.back()->Sumw2();
      muonNtrackerLayer_mc.push_back(new TH1F(Form("ntrackerlayer_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						     int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",30,0,30));
      muonNtrackerLayer_mc.back()->Sumw2();
      muonNStation_mc.push_back(new TH1F(Form("nstation_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",8,0,8));
      muonNStation_mc.back()->Sumw2();
      muonDxy_mc.push_back(new TH1F(Form("dxy_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					   int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",40,-0.25,0.25));
      muonDxy_mc.back()->Sumw2();
      muonDz_mc.push_back(new TH1F(Form("dz_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",40,-0.5,0.5));
      muonDz_mc.back()->Sumw2();
      muonNhiso_mc.push_back(new TH1F(Form("nhiso_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					     int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",40,0,15));
      muonNhiso_mc.back()->Sumw2();
      muonChiso_mc.push_back(new TH1F(Form("chiso_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",					     
					     int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",40,0,15));					
      muonChiso_mc.back()->Sumw2();
      muonPuiso_mc.push_back(new TH1F(Form("puiso_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					     int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",40,0,15));
      muonPuiso_mc.back()->Sumw2();
      muonEMIso_mc.push_back(new TH1F(Form("emiso_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					   int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",40,0,25));
      muonEMIso_mc.back()->Sumw2();

      muonZMass_mc.push_back(new TH1F(Form("zmass_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
					     int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",25,60,120));
      muonZMass_mc.back()->Sumw2();

      /// create efficiency histograms
      muonEfficiency_num_data.push_back(new TH1F(Form("mueff_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						      int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",10,0,10));
      muonEfficiency_num_mc.push_back(new TH1F(Form("mueff_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						    int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",10,0,10));

      muonEfficiency_num_data.back()->Sumw2();
      muonEfficiency_num_mc.back()->Sumw2();

      muonEfficiency_den_data.push_back(new TH1F(Form("mueff_den_data_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						      int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",10,0,10));
      muonEfficiency_den_mc.push_back(new TH1F(Form("mueff_den_mc_ptMin_%d_ptMax_%d_etaMin_%.1f_etaMax_%.1f",
						    int(ptBin.at(ipt)),int(ptBin.at(ipt+1)),etaBin.at(ieta),etaBin.at(ieta+1)),"",10,0,10));

      muonEfficiency_den_data.back()->Sumw2();
      muonEfficiency_den_mc.back()->Sumw2();

      muonEfficiency_num_data.back()->GetXaxis()->SetBinLabel(1,"Global ID");
      muonEfficiency_num_data.back()->GetXaxis()->SetBinLabel(2,"PF ID");
      muonEfficiency_num_data.back()->GetXaxis()->SetBinLabel(3,"#chi^{2}");
      muonEfficiency_num_data.back()->GetXaxis()->SetBinLabel(4,"N_{valid}");
      muonEfficiency_num_data.back()->GetXaxis()->SetBinLabel(5,"N_{station}");
      muonEfficiency_num_data.back()->GetXaxis()->SetBinLabel(6,"dxy");
      muonEfficiency_num_data.back()->GetXaxis()->SetBinLabel(7,"dz");
      muonEfficiency_num_data.back()->GetXaxis()->SetBinLabel(8,"N_{pixel}");
      muonEfficiency_num_data.back()->GetXaxis()->SetBinLabel(9,"N_{layer}");
      muonEfficiency_num_data.back()->GetXaxis()->SetBinLabel(10,"Iso");

      muonEfficiency_num_mc.back()->GetXaxis()->SetBinLabel(1,"Global ID");
      muonEfficiency_num_mc.back()->GetXaxis()->SetBinLabel(2,"PF ID");
      muonEfficiency_num_mc.back()->GetXaxis()->SetBinLabel(3,"#chi^{2}");
      muonEfficiency_num_mc.back()->GetXaxis()->SetBinLabel(4,"N_{valid}");
      muonEfficiency_num_mc.back()->GetXaxis()->SetBinLabel(5,"N_{station}");
      muonEfficiency_num_mc.back()->GetXaxis()->SetBinLabel(6,"dxy");
      muonEfficiency_num_mc.back()->GetXaxis()->SetBinLabel(7,"dz");
      muonEfficiency_num_mc.back()->GetXaxis()->SetBinLabel(8,"N_{pixel}");
      muonEfficiency_num_mc.back()->GetXaxis()->SetBinLabel(9,"N_{layer}");
      muonEfficiency_num_mc.back()->GetXaxis()->SetBinLabel(10,"Iso");
      
    }
  }

  // create pt distributions per eta bins
  for(size_t ieta = 0; ieta < etaBin.size()-1; ieta++){
    muonPt_data.push_back(new TH1F(Form("muonpt_data_etaMin_%.1f_etaMax_%1.f",etaBin.at(ieta),etaBin.at(ieta+1)),"",60,10,350));
    muonPt_data.back()->Sumw2();
    muonPt_mc.push_back(new TH1F(Form("muonpt_mc_etaMin_%.1f_etaMax_%1.f",etaBin.at(ieta),etaBin.at(ieta+1)),"",60,10,350));
    muonPt_mc.back()->Sumw2();
  }
  // create eta distributions per pt bin
  for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++){
    muonEta_data.push_back(new TH1F(Form("muoneta_data_ptMin_%d_ptMax_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",25,-2.4,2.4));
    muonEta_data.back()->Sumw2();
    muonEta_mc.push_back(new TH1F(Form("muoneta_mc_ptMin_%d_ptMax_%d",int(ptBin.at(ipt)),int(ptBin.at(ipt+1))),"",25,-2.4,2.4));
    muonEta_mc.back()->Sumw2();
  }

  size_t ipos = 0;
  // loop on the event
  cout<<"Start filling ID variables "<<endl;
  for(size_t ieta = 0; ieta < etaBin.size()-1; ieta++){
    cout<<"Eta bin "<<etaBin.at(ieta)<<" : "<<etaBin.at(ieta+1)<<endl;
    for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++){
      cout<<"Pt Bin Data : "<<ptBin.at(ipt)<<" : "<<ptBin.at(ipt+1)<<endl;
      fillHistograms(chain_data,ptBin.at(ipt),ptBin.at(ipt+1),etaBin.at(ieta),etaBin.at(ieta+1),ipos,true,1);
      cout<<"Pt Bin MC   : "<<ptBin.at(ipt)<<" : "<<ptBin.at(ipt+1)<<endl;
      fillHistograms(chain_mc,ptBin.at(ipt),ptBin.at(ipt+1),etaBin.at(ieta),etaBin.at(ieta+1),ipos,false,weightsum);
      ipos++;
    }
  }

  cout<<"Start filling pt spectrum "<<endl;
  // Fill inclusive pt histograms
  for(size_t ieta = 0; ieta < etaBin.size()-1; ieta++){
    cout<<"Eta bin Data : "<<etaBin.at(ieta)<<" : "<<etaBin.at(ieta+1)<<endl;
    fillPtHistograms(chain_data,etaBin.at(ieta),etaBin.at(ieta+1),ieta,true,1);
    cout<<"Eta bin MC   : "<<etaBin.at(ieta)<<" : "<<etaBin.at(ieta+1)<<endl;
    fillPtHistograms(chain_mc,etaBin.at(ieta),etaBin.at(ieta+1),ieta,false,weightsum);
  }

  cout<<"Start filling eta spectrum "<<endl;
  // Fill inclusive pt histograms
  for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++){
    cout<<"Pt bin Data : "<<ptBin.at(ipt)<<" : "<<ptBin.at(ipt+1)<<endl;
    fillEtaHistograms(chain_data,ptBin.at(ipt),ptBin.at(ipt+1),ipt,true,1);
    cout<<"Pt bin MC   : "<<ptBin.at(ipt)<<" : "<<ptBin.at(ipt+1)<<endl;
    fillEtaHistograms(chain_mc,ptBin.at(ipt),ptBin.at(ipt+1),ipt,false,weightsum);
  }


  TCanvas* canvas = new TCanvas("canvas","",600,700);
  cout<<"Start plotting "<<endl;
  // plot histograms
  ipos = 0;

  for(size_t ieta = 0; ieta < etaBin.size()-1; ieta++){
    for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++){
      plotHistogram(muonChi2_data.at(ipos),muonChi2_mc.at(ipos),canvas,ouputDIR,"#chi^{2} track");
      plotHistogram(muonPfId_data.at(ipos),muonPfId_mc.at(ipos),canvas,ouputDIR,"PF ID");
      plotHistogram(muonGlobalId_data.at(ipos),muonGlobalId_mc.at(ipos),canvas,ouputDIR,"Global ID");
      plotHistogram(muonNvalidHit_data.at(ipos),muonNvalidHit_mc.at(ipos),canvas,ouputDIR,"N_{valid}^{hit}");
      plotHistogram(muonNpixelHit_data.at(ipos),muonNpixelHit_mc.at(ipos),canvas,ouputDIR,"N_{pixel}^{hit}");
      plotHistogram(muonNtrackerLayer_data.at(ipos),muonNtrackerLayer_mc.at(ipos),canvas,ouputDIR,"N_{tracker}^{layer}");
      plotHistogram(muonNStation_data.at(ipos), muonNStation_mc.at(ipos),canvas,ouputDIR,"N_{station}");
      plotHistogram(muonDxy_data.at(ipos), muonDxy_mc.at(ipos),canvas,ouputDIR,"d_{xy}");
      plotHistogram(muonDz_data.at(ipos), muonDz_mc.at(ipos),canvas,ouputDIR,"d_{z}");
      plotHistogram(muonNhiso_data.at(ipos), muonNhiso_mc.at(ipos),canvas,ouputDIR,"Iso_{NH} [GeV]");
      plotHistogram(muonChiso_data.at(ipos), muonChiso_mc.at(ipos),canvas,ouputDIR,"Iso_{CH} [GeV]");
      plotHistogram(muonPuiso_data.at(ipos), muonPuiso_mc.at(ipos),canvas,ouputDIR,"Iso_{PU} [GeV]");
      plotHistogram(muonEMIso_data.at(ipos), muonEMIso_mc.at(ipos),canvas,ouputDIR,"Iso_{EM} [GeV]");      
      plotHistogram(muonZMass_data.at(ipos), muonZMass_mc.at(ipos),canvas,ouputDIR,"M_{#mu#mu} [GeV]");      
      plotEfficiency(muonEfficiency_num_data.at(ipos),muonEfficiency_den_data.at(ipos),muonEfficiency_num_mc.at(ipos),muonEfficiency_den_mc.at(ipos),canvas,ouputDIR);
      ipos++;
    }
  }

  for(size_t ieta = 0; ieta < etaBin.size()-1; ieta++)
    plotHistogram(muonPt_data.at(ieta), muonPt_mc.at(ieta),canvas,ouputDIR,"p_{T}^{muon} [GeV]");

  for(size_t ipt = 0; ipt < ptBin.size()-1; ipt++)
    plotHistogram(muonEta_data.at(ipt), muonEta_mc.at(ipt),canvas,ouputDIR,"#eta^{muon}");

  cout<<"End plotting "<<endl;
}

