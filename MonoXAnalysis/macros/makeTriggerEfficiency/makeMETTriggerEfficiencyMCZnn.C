#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

// recoil binning for monojet                                                                                                                                                                    
vector <float> bins_monojet_recoil = {50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 265., 280., 300, 320., 340., 360., 380., 400., 430., 460., 490., 520, 550., 580., 610., 650., 700., 740., 800., 900., 1000.,1250};

enum class Sample {sig,wmn,wnmu};

void calculateSumWeight(vector<TTree*> gentree, vector<double> & wgtsum){

  for(auto tree: gentree){
    wgtsum.push_back(0);
  }

  int itree = 0;
  for(auto tree: gentree){
    TTreeReader reader(tree);
    TTreeReaderValue<float> wgt (reader,"wgt");
    //////////////////                                                                                                                                                                           
    long int nTotal = tree->GetEntries();
    long int nEvents = 0;
    long int nPart = 100000;
    cout<<"Looping on itree "<<itree<<" of "<<gentree.size()<<" Total number of events: "<<nTotal<<endl;
    while(reader.Next()){
      if(nEvents > 500000) break;
      cout.flush();
      if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
      nEvents++;
      wgtsum.at(itree) += *wgt;
    }
    cout<<endl;
    cout<<"Sum of weigths for this tree "<<wgtsum.at(itree)<<endl;
    itree++;
  }
}

//////////
void makeTriggerAnalysis(vector<TTree*> trees, TH1F* hnum, TH1F* hden, const Sample & sample, vector<double> wgtsum, float luminosity){
  
  cout<<"Loop on trees "<<endl;
  TFile* pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt_36.40_summer16.root");
  TH1* puhist = (TH1*)pufile->Get("puhist");

  TFile* singleMuon = TFile::Open("triggerEfficiency_tagAndProbe_singleMu_v3/triggerEfficiency_DATA_SingleMuon.root");
  TEfficiency* singleMu_eff = (TEfficiency*) singleMuon->Get("trgeff_mu");
  TH2F* singleMu_efficiency = (TH2F*) singleMu_eff->CreateHistogram();

  int itree = 0;
  for(auto tree : trees){    

    TTreeReader reader(tree);
    TTreeReaderValue<float> xsec         (reader,"xsec");
    TTreeReaderValue<float> wgt          (reader,"wgt");
    TTreeReaderValue<UChar_t> hltm90     (reader,"hltmet90");
    TTreeReaderValue<UChar_t> hltm100    (reader,"hltmet100");
    TTreeReaderValue<UChar_t> hltm110    (reader,"hltmet110");
    TTreeReaderValue<UChar_t> hltm120    (reader,"hltmet120");
    TTreeReaderValue<UChar_t> hltmwm120  (reader,"hltmetwithmu120");
    TTreeReaderValue<UChar_t> hltmwm100  (reader,"hltmetwithmu100");
    TTreeReaderValue<UChar_t> hltmwm110  (reader,"hltmetwithmu110");
    TTreeReaderValue<UChar_t> hltmwm170  (reader,"hltmetwithmu170");
    TTreeReaderValue<UChar_t> hltmwm300  (reader,"hltmetwithmu300");
    TTreeReaderValue<UChar_t> hltmwm90   (reader,"hltmetwithmu90");
    TTreeReaderValue<UChar_t> hltjm      (reader,"hltjetmet");
    TTreeReaderValue<UChar_t> hltsinglemu (reader,"hltsinglemu");
    TTreeReaderValue<UChar_t> hltdoublemu (reader,"hltdoublemu");
    TTreeReaderValue<float>   mu1pt      (reader,"mu1pt");
    TTreeReaderValue<float>   mu1eta     (reader,"mu1eta");
    TTreeReaderValue<float>   mu1phi     (reader,"mu1phi");
    TTreeReaderValue<int>     mu1id      (reader,"mu1id");
    TTreeReaderValue<int>     mu1pid     (reader,"mu1pid");
    TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
    TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
    TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
    TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
    TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
    TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
    TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
    TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
    TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");
    TTreeReaderValue<unsigned int> nvtx        (reader,"nvtx");
    TTreeReaderValue<unsigned int> ntaus       (reader,"ntaus");
    TTreeReaderValue<unsigned int> nmuons      (reader,"nmuons");
    TTreeReaderValue<unsigned int> nelectrons  (reader,"nelectrons");
    TTreeReaderValue<unsigned int> nphotons    (reader,"nphotons");
    TTreeReaderValue<unsigned int> nincjets    (reader,"njetsinc");
    TTreeReaderValue<unsigned int> nbjets      (reader,"nbjetslowpt");
    TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
    TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
    TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
    TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
    TTreeReaderValue<vector<float> > jetchfrac  (reader,"combinejetCHfrac");
    TTreeReaderValue<vector<float> > jetnhfrac  (reader,"combinejetNHfrac");
    TTreeReaderValue<float> met         (reader,"t1pfmet");
    TTreeReaderValue<float> metphi      (reader,"t1pfmetphi");
    TTreeReaderValue<float> mmet        (reader,"t1mumet");
    TTreeReaderValue<float> mmetphi     (reader,"t1mumetphi");
    TTreeReaderValue<float> metpf       (reader,"pfmet");
    TTreeReaderValue<float> metcalo     (reader,"calomet");
    TTreeReaderValue<float> jmmdphi (reader,"incjetmumetdphimin4");

    //////////////////                                                                                                                                                                              
    long int nTotal = tree->GetEntries();
    cout<<"Looping on itree "<<itree<<" of "<<trees.size()<<" Total number of events: "<<nTotal<<endl;
    long int nEvents = 0;    
    long int nPart = 100000;

    while(reader.Next()){
      cout.flush();
      if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
      nEvents++;
      if(nEvents > 500000) break;

      if(*nbjets   != 0)    continue;
      if(*ntaus    != 0)    continue;
      if(*nphotons  != 0)   continue;
      if(*nelectrons != 0 ) continue;

      if(not *fcsc)  continue;
      if(not *fcsct) continue;
      if(not *feeb)  continue;
      if(not *fetp)  continue;
      if(not *fvtx)  continue;
      if(not *fbadmu) continue;
      if(not *fbadch) continue;
      if(not *fhbhe)  continue;
      if(not *fhbiso) continue;

      // apply calo met cleaning
      if(fabs(*metpf-*metcalo)/(*mmet) > 0.5) continue;
      // ask single muon trigger
      if(sample == Sample::wmn and not *hltsinglemu) continue;
      else if(sample == Sample::wnmu and not (*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm)) continue;
      
      // apply jet pt selections
      if(*nincjets < 1) continue;
      if(jetpt->size() == 0) continue;
      if(jetpt->at(0) < 100) continue;
      if(fabs(jeteta->at(0)) > 2.5) continue;
      if(jetchfrac->at(0) < 0.1) continue;
      if(jetnhfrac->at(0) > 0.8) continue;
      // apply jet-met dphi--> not for gen level analysis
      if(*jmmdphi < 0.5) continue;

      // apply standard reco-level cuts
      if(sample == Sample::wmn){
	if(*mu1pt < 20) continue;
	if(fabs(*mu1eta) > 2.4) continue;
	if(*mu1id  !=1) continue;
	if(*nmuons !=1) continue;
	float dphi = fabs(*mu1phi-*metphi);
	if(dphi > TMath::Pi())
	  dphi = 2*TMath::Pi()-dphi;
	float mtw = sqrt(2*(*mu1pt)*(*met)*(1-cos(dphi)));
	if(mtw > 160) continue;
      }
      else if(sample == Sample::wnmu){
	if(*mu1pt < 20) continue;
        if(fabs(*mu1eta) > 2.4) continue;
        if(*mu1id  !=1) continue;
        if(*nmuons !=1) continue;
        float dphi = fabs(*mu1phi-*metphi);
        if(dphi > TMath::Pi())
          dphi = 2*TMath::Pi()-dphi;
        float mtw = sqrt(2*(*mu1pt)*(*met)*(1-cos(dphi)));
        if(mtw > 160) continue;
      }
      else if(sample == Sample::sig){
	  if(*nmuons != 0) continue;
      }
      
      // pileup re-weight
      double puwgt = 1;
      if(*nvtx < 60)
	puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));

      double mueff = singleMu_efficiency->GetBinContent(singleMu_efficiency->FindBin(*mu1eta,min(double(*mu1pt),singleMu_efficiency->GetYaxis()->GetXmax()-1)));     
      hden->Fill(*mmet,luminosity*(*xsec)*(*wgt)*puwgt*mueff/wgtsum.at(itree));

      if(sample == Sample::wmn and (*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm))
	hnum->Fill(*mmet,luminosity*(*xsec)*(*wgt)*puwgt*mueff/wgtsum.at(itree));     
      else if(sample == Sample::wnmu and not *hltsinglemu)
	hnum->Fill(*mmet,luminosity*(*xsec)*(*wgt)*puwgt/wgtsum.at(itree));
    }
    
    cout<<endl;
    itree++;
  }
  
  if(pufile) pufile->Close();
  if(singleMuon) singleMuon->Close();
}

void makeMETTriggerEfficiencyMC(string inputDIR, string outputDIR, float luminosity = 35.9, bool applyGenSelection = false){

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TFile* inputZvvTurnOn = TFile::Open("triggerEfficiencyMC_withTrigDen_v2/metTriggerEfficiencyMC_monojet_recoil.root");
  TF1* fitfunc_monojet_recoil_zvv = (TF1*) inputZvvTurnOn->Get("func_zvv");
  TEfficiency* eff_monojet_recoil_zvv = (TEfficiency*) inputZvvTurnOn->Get("efficiency_zvv");

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
 
  // input tree                                                                                                                                                                                      
  vector<TFile*> filelist_wjet;
  vector<TTree*> tree_wjet;
  vector<TTree*> gentree_wjet;

  cout<<"Read list of files for Wjet process "<<endl;
  system(("ls "+inputDIR+" | grep WJets > list_dir.txt").c_str());
  ifstream file_dir_wjet("list_dir.txt");
  if(file_dir_wjet.is_open()){
    string line;
    while(!file_dir_wjet.eof()){
      getline(file_dir_wjet,line);
      if(line == "") continue;
      system(("find "+inputDIR+"/"+line+" -name  \"*.root\" > list.txt").c_str());
      ifstream file_wjet("list.txt");
      if(file_wjet.is_open()){
	string line2;
	while(!file_wjet.eof()){
	  getline(file_wjet,line2);
	  if(TString(line2).Contains("failed")) continue;	  
	  if(line == "" or not TString(line2).Contains("root")) continue;
	  cout<<"Open Wjet file with name: "<<line2<<endl;
	  filelist_wjet.push_back(TFile::Open(line2.c_str()));
	  tree_wjet.push_back((TTree*) filelist_wjet.back()->Get("tree/tree"));
	  gentree_wjet.push_back((TTree*) filelist_wjet.back()->Get("gentree/gentree"));
	}
      }
      system("rm list.txt");
    }
  }
  system("rm list_dir.txt");

  
  // sum of weights
  vector<double> wgtsum_wjet;
  cout<<"######### Calculate Weights for : W+jets "<<endl;
  calculateSumWeight(gentree_wjet,wgtsum_wjet);

  /////// start analysis
  TF1 *fitfunc_monojet_recoil_wmn = new TF1("fitfunc_monojet_recoil_wmn",ErfCB,bins_monojet_recoil.front(), bins_monojet_recoil.back(),5);
  fitfunc_monojet_recoil_wmn->SetParameters(120., 25., 30., 4., 1.);  
  TH1F* hnum_monojet_recoil_wmn = new TH1F("hnum_monojet_recoil_wmn", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  TH1F* hden_monojet_recoil_wmn = new TH1F("hden_monojet_recoil_wmn", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  hnum_monojet_recoil_wmn->Sumw2();
  hden_monojet_recoil_wmn->Sumw2();

  /////// start analysis  
  TF1 *fitfunc_monojet_recoil_wmn_2 = new TF1("fitfunc_monojet_recoil_wmn_2",ErfCB,bins_monojet_recoil.front(), bins_monojet_recoil.back(),5);
  fitfunc_monojet_recoil_wmn_2->SetParameters(120., 25., 30., 4., 1.);  
  TH1F* hnum_monojet_recoil_wmn_2 = new TH1F("hnum_monojet_recoil_wmn_2", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  TH1F* hden_monojet_recoil_wmn_2 = new TH1F("hden_monojet_recoil_wmn_2", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  hnum_monojet_recoil_wmn_2->Sumw2();
  hden_monojet_recoil_wmn_2->Sumw2();

  ////
  cout<<"######### Loop on W+jets trees for Wmn selection "<<endl;
  makeTriggerAnalysis(tree_wjet,hnum_monojet_recoil_wmn,hden_monojet_recoil_wmn,Sample::wmn,wgtsum_wjet,luminosity);
  cout<<"######### Loop on W+jets trees for Wmn selection "<<endl;
  makeTriggerAnalysis(tree_wjet,hnum_monojet_recoil_wmn_2,hden_monojet_recoil_wmn_2,Sample::wnmu,wgtsum_wjet,luminosity);

  // Make efficiencies
  eff_monojet_recoil_zvv->SetMarkerColor(kBlack);
  eff_monojet_recoil_zvv->SetLineColor(kBlack);
  eff_monojet_recoil_zvv->SetMarkerStyle(20);
  eff_monojet_recoil_zvv->SetMarkerSize(1);
  fitfunc_monojet_recoil_zvv->SetLineColor(kBlack);
  fitfunc_monojet_recoil_zvv->SetLineWidth(2);
  fitfunc_monojet_recoil_zvv->SetLineStyle(7);

  TEfficiency* eff_monojet_recoil_wmn = new TEfficiency(*hnum_monojet_recoil_wmn,*hden_monojet_recoil_wmn);
  eff_monojet_recoil_wmn->SetMarkerColor(kBlue);
  eff_monojet_recoil_wmn->SetLineColor(kBlue);
  eff_monojet_recoil_wmn->SetMarkerStyle(20);
  eff_monojet_recoil_wmn->SetMarkerSize(1);
  fitfunc_monojet_recoil_wmn->SetLineColor(kBlue);
  fitfunc_monojet_recoil_wmn->SetLineWidth(2);
  fitfunc_monojet_recoil_wmn->SetLineStyle(7);
  eff_monojet_recoil_wmn->Fit(fitfunc_monojet_recoil_wmn,"RE");

  TEfficiency* eff_monojet_recoil_wmn_2 = new TEfficiency(*hnum_monojet_recoil_wmn_2,*hden_monojet_recoil_wmn_2);
  eff_monojet_recoil_wmn_2->SetMarkerColor(kBlue);
  eff_monojet_recoil_wmn_2->SetLineColor(kBlue);
  eff_monojet_recoil_wmn_2->SetMarkerStyle(20);
  eff_monojet_recoil_wmn_2->SetMarkerSize(1);
  fitfunc_monojet_recoil_wmn_2->SetLineColor(kBlue);
  fitfunc_monojet_recoil_wmn_2->SetLineWidth(2);
  fitfunc_monojet_recoil_wmn_2->SetLineStyle(7);
  eff_monojet_recoil_wmn_2->Fit(fitfunc_monojet_recoil_wmn_2,"RE");


  TGraphAsymmErrors* graph_zvv = eff_monojet_recoil_zvv->CreateGraph();
  TGraphAsymmErrors* graph_wmn = eff_monojet_recoil_wmn->CreateGraph();
  TGraphAsymmErrors* graph_wmn_2 = eff_monojet_recoil_wmn_2->CreateGraph();
  TGraphAsymmErrors* graph_eff = (TGraphAsymmErrors*) graph_wmn->Clone("graph_eff");
  for(int ipoint = 0; ipoint < graph_wmn->GetN(); ipoint++){
    graph_eff->SetPoint(0,0,0);
    double x1,y1;
    graph_wmn->GetPoint(ipoint,x1,y1);
    double x2,y2;
    graph_wmn_2->GetPoint(ipoint,x2,y2);
    graph_eff->SetPoint(ipoint,x1,y1/(1-y2));
  }
  
  

  // Plotting final result for MC turn ons
  TH1* frame = canvas->DrawFrame(bins_monojet_recoil.front(),0.,bins_monojet_recoil.back(), 1.1, "");
  frame->GetXaxis()->SetTitle("Recoil [GeV]");
  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetYaxis()->SetLabelSize(0.8*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetTitleSize(0.8*frame->GetYaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleOffset(1.0);

  canvas->SetRightMargin(0.075);
  canvas->SetTopMargin(0.06);
  canvas->Draw();
  canvas->cd();
  frame->Draw();
  CMS_lumi(canvas,string(Form("%.2f",luminosity)),true);
  
  eff_monojet_recoil_zvv->Draw("Lsame");
  graph_eff->Draw("Lsame");
  
  TLegend leg (0.6,0.3,0.9,0.5);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(eff_monojet_recoil_zvv,"Z #rightarrow #nu#nu","EPL");
  leg.AddEntry(graph_eff,"Estimated","EPL");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/metTriggerEfficiencyMC_monojet_recoil.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/metTriggerEfficiencyMC_monojet_recoil.pdf").c_str(),"pdf");

  // Plotting final result for MC turn ons
  TH1* frame2 = canvas->DrawFrame(150,0.85,600,1.05, "");
  frame2->GetXaxis()->SetTitle("Recoil [GeV]");
  frame2->GetYaxis()->SetTitle("Trigger Efficiency");
  frame2->GetYaxis()->SetLabelSize(0.8*frame->GetYaxis()->GetLabelSize());
  frame2->GetXaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame2->GetYaxis()->SetTitleSize(0.8*frame->GetYaxis()->GetTitleSize());
  frame2->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame2->GetXaxis()->SetTitleOffset(1.0);
  
  canvas->cd();
  frame2->Draw();
  CMS_lumi(canvas,string(Form("%.2f",luminosity)),true);
  
  eff_monojet_recoil_zvv->Draw("Lsame");
  graph_eff->Draw("Lsame");
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/metTriggerEfficiencyMC_monojet_recoil_zoom.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/metTriggerEfficiencyMC_monojet_recoil_zoom.pdf").c_str(),"pdf");

  TFile* outputFile = new TFile((outputDIR+"/metTriggerEfficiencyComparison.root").c_str(),"RECREATE");
  outputFile->cd();

  fitfunc_monojet_recoil_zvv->Write("func_zvv");
  fitfunc_monojet_recoil_wmn->Write("func_wmn");
  fitfunc_monojet_recoil_wmn_2->Write("func_wnmu");
  eff_monojet_recoil_zvv->Write("efficiency_zvv");
  eff_monojet_recoil_wmn->Write("efficiency_wmn");
  eff_monojet_recoil_wmn_2->Write("efficiency_wnmu");
  outputFile->Close();
}
