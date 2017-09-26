#include <iostream>
#include <sstream>
#include <cmath>
#include "../triggerUtils.h"
#include "../../CMS_lumi.h"

// recoil binning
vector <float> bins_monojet_recoil = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 225., 250., 275., 300., 325., 350., 400., 450., 500., 550., 650., 800., 1000., 1250};
vector <float> bins_monoV_recoil   = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 180., 200., 250., 300., 350., 400., 450., 500., 550., 650., 800., 1000., 1500};

// eras
vector<string> RunEra = {"Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G","Run2016H"};

// Options
static bool  drawUncertaintyBand   = false;
static bool  useDoubleMuonTriggers = true;
static bool  applyJetSelections    = true;

// Samples
enum class Sample {wmn,wen,zmm,zee};

/// plotting result
void plotTurnOn(TCanvas* canvas, TEfficiency* eff, TF1* fitfunc, const string & axisLabel, const string & postfix, const string & outputDIR);

/// confidence interval for the band
void GetConfidenceIntervals(TF1* funz, TH1F* obj, Double_t cl, const TFitResultPtr & fitResult);

/// main function
void makeMETTriggerEfficiencyMonoJet_Data(string inputDIR, string outputDIR, Sample sample, bool doubleLepton = false) {

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();

  // input tree
  TChain* tree = new TChain("tree/tree");
  // use only a subset of directories
  if(sample == Sample::zmm or sample == Sample::wmn)
    system(("ls "+inputDIR+"  | grep SingleMu > list_dir.txt").c_str());
  else if(sample == Sample::wen)
    system(("ls "+inputDIR+"  | grep SingleEle > list_dir.txt").c_str());
  ifstream dirlist ("list_dir.txt");
  string dirname;
  if(dirlist.is_open()){
    while(not dirlist.eof()){      
      getline(dirlist,dirname);
      bool found = false;
      for(auto era : RunEra){
	if(dirname.find(era) != string::npos)
	  found = true;
      }      	
      if(found == false) continue;
      system(("find "+inputDIR+"/"+dirname+" -name  \"*.root\" > list.txt").c_str());
      ifstream file("list.txt");                                                                                                                                                                   
      if(file.is_open()){                                                                                                                                                                        
	string line;                                                                                                                                                                               
	while(!file.eof()){                                                                                                                                                                       
	  getline(file,line);                                         
	  if(TString(line).Contains("failed")) continue;
	  if(line == "" or not TString(line).Contains("root")) continue;	  
	  cout<<"adding following file: "<<line<<endl;
	  tree->Add(line.c_str());
	}
      }
      system("rm list.txt");
    }
  }
  system("rm list_dir.txt");

  ////
  string postfix = "recoil_monojet";
  if(sample == Sample::wmn) postfix += "_Wmn";
  else if(sample == Sample::wen) postfix += "_Wen";
  else if(sample == Sample::zmm) postfix += "_Zmm";
  else if(sample == Sample::zee) postfix += "_Zee";

  TFile* outputFile = new TFile((outputDIR+"/triggerEfficiency_"+postfix+".root").c_str(),"RECREATE");
  outputFile->cd();

  // fitting function for the turn-on  as a function of recoil
  TF1 *fitfunc_monojet_recoil = new TF1("fitfunc_monojet_recoil","[0]*1./((1.+[1]*exp(-[2]*(x-[3])))^(1./[4]))",bins_monojet_recoil.front(), bins_monojet_recoil.back());
  fitfunc_monojet_recoil->SetParameters(1.,0.01,0.02,50,0.01);
  fitfunc_monojet_recoil->SetParLimits(0,0.,1.01);
  fitfunc_monojet_recoil->SetParLimits(1,-100.,100.);
  fitfunc_monojet_recoil->SetParLimits(3,0.,500.);
  fitfunc_monojet_recoil->SetParLimits(4,0.,100);

  TH1F* hnum_monojet_recoil = new TH1F("hnum_monojet_recoil", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  TH1F* hden_monojet_recoil = new TH1F("hden_monojet_recoil", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  hnum_monojet_recoil->Sumw2();
  hden_monojet_recoil->Sumw2();
  
  // define numerator as event with tight muon + trigger requirement
  // define denominator as an event with a tight muon passing single muon trigger
  TTreeReader reader(tree);
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
  TTreeReaderValue<UChar_t> hlte      (reader,"hltsingleel");
  TTreeReaderValue<float>   mu1pt     (reader,"mu1pt");
  TTreeReaderValue<float>   mu1eta    (reader,"mu1eta");
  TTreeReaderValue<float>   mu1phi    (reader,"mu1phi");
  TTreeReaderValue<int>     mu1id     (reader,"mu1id");
  TTreeReaderValue<int>     mu1pid     (reader,"mu1pid");
  TTreeReaderValue<float>   mu2pt     (reader,"mu2pt");
  TTreeReaderValue<float>   mu2eta    (reader,"mu2eta");
  TTreeReaderValue<float>   mu2phi    (reader,"mu2phi");
  TTreeReaderValue<int>     mu2id     (reader,"mu2id");
  TTreeReaderValue<int>     mu2pid     (reader,"mu2pid");
  TTreeReaderValue<float>   el1pt     (reader,"el1pt");
  TTreeReaderValue<float>   el1eta    (reader,"el1eta");
  TTreeReaderValue<float>   el1phi    (reader,"el1phi");
  TTreeReaderValue<int>     el1id     (reader,"el1id");
  TTreeReaderValue<int>     el1pid    (reader,"el1pid");
  TTreeReaderValue<float>   el2pt     (reader,"el2pt");
  TTreeReaderValue<float>   el2eta    (reader,"el2eta");
  TTreeReaderValue<float>   el2phi    (reader,"el2phi");
  TTreeReaderValue<int>     el2id     (reader,"el2id");
  TTreeReaderValue<int>     el2pid    (reader,"el2pid");
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");                                                                                                                                       
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");                                                                                                                                       
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
  TTreeReaderValue<float> emet        (reader,"t1elmet");
  TTreeReaderValue<float> emetphi     (reader,"t1elmetphi");  
  TTreeReaderValue<float> metpf       (reader,"pfmet");
  TTreeReaderValue<float> metcalo     (reader,"calomet");
  TTreeReaderValue<float> jmmdphi (reader,"incjetmumetdphimin4");
  TTreeReaderValue<float> jemdphi (reader,"incjetelmetdphimin4");
  TTreeReaderValue<float> zmass (reader,"zmass");
  TTreeReaderValue<float> zeemass (reader,"zeemass");
 
  //////////////////
  long int nTotal = tree->GetEntries();
  cout<<"Total number of events: "<<nTotal<<endl;
  long int nEvents = 0;

  TH1D* efficiencyMonojetSelections = new TH1D("efficiencyMonojetSelections","efficiencyMonojetSelections",14,0,13);
  efficiencyMonojetSelections->Sumw2();
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(1,"Total");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(2,"B-veto");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(3,"Tau-veto");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(4,"Photon-veto");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(5,"MET filters");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(6,"HLT single muon");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(7,"Lepton pT and eta");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(8,"Lepton ID");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(9,"Lepton Veto");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(10,"Calo-PF MET");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(11,"Jet met dphi");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(12,"Jet cuts");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(13,"Trigger");

  long int nPart = 100000;
  while(reader.Next()){

    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    // define the denominator
    efficiencyMonojetSelections->SetBinContent(1,efficiencyMonojetSelections->GetBinContent(1)+1);
    if(*nbjets != 0) continue;
    efficiencyMonojetSelections->SetBinContent(2,efficiencyMonojetSelections->GetBinContent(2)+1);
    ///
    if(*ntaus != 0) continue;
    ///
    efficiencyMonojetSelections->SetBinContent(3,efficiencyMonojetSelections->GetBinContent(3)+1);
    if(*nphotons  != 0) continue;
    efficiencyMonojetSelections->SetBinContent(4,efficiencyMonojetSelections->GetBinContent(4)+1);
    //
    if(not *fcsc)  continue;
    if(not *feeb)  continue;
    if(not *fetp)  continue;
    if(not *fvtx)  continue;
    if(not *fbadmu) continue;
    if(not *fbadch) continue; 
    if(not *fhbhe)  continue;
    if(not *fhbiso) continue;
    efficiencyMonojetSelections->SetBinContent(5,efficiencyMonojetSelections->GetBinContent(5)+1);
    
    if(sample == Sample::wmn or sample == Sample::zmm){
      // trigger selection for the denominator
      if(not doubleLepton){
	if(not *hltsinglemu) continue;    
      }
      else if(doubleLepton and not useDoubleMuonTriggers){
	if(not *hltsinglemu) continue;
      }
      else if(doubleLepton and useDoubleMuonTriggers){
	if(not *hltdoublemu) continue;
	if(not *hltsinglemu) continue;
      }
    }
    else if(sample == Sample::wen or sample == Sample::zee){
      if(not *hlte) continue;
    }

    efficiencyMonojetSelections->SetBinContent(6,efficiencyMonojetSelections->GetBinContent(6)+1);

    if(sample == Sample::wmn or sample == Sample::zmm){
      if(*mu1pt < 20) continue;
      if(fabs(*mu1eta) > 2.4) continue;
      efficiencyMonojetSelections->SetBinContent(7,efficiencyMonojetSelections->GetBinContent(7)+1);
    }
    else if(sample == Sample::wen or sample == Sample::zee){
      if(*el1pt < 40) continue;
      if(fabs(*el1eta) > 2.5) continue;
    }
    ////
    if(sample == Sample::wmn){
      if(*mu1id != 1) continue;
      if(*nmuons != 1) continue;
      float dphi = fabs(*mu1phi-*metphi);
      if(dphi > TMath::Pi())
	dphi = 2*TMath::Pi()-dphi;
      float mtw = sqrt(2*(*mu1pt)*(*met)*(1-cos(dphi)));
      if(mtw > 160) continue;
    }
    else if (sample == Sample::zmm){
	if(*nmuons != 2) continue;
	if(not ((*mu1pt > 20 and *mu1id == 1) or (*mu2pt > 20 and *mu2id == 1))) continue;
	if(*mu1pid == *mu2pid) continue;
	if(*zmass < 60 or *zmass > 120) continue;
      }
      else if (sample == Sample::wen){
	if(*el1id != 1) continue;
	if(*nelectrons != 1) continue;
	float dphi = fabs(*el1phi-*metphi);
	if(dphi > TMath::Pi())
	  dphi = 2*TMath::Pi()-dphi;
	float mtw = sqrt(2*(*el1pt)*(*met)*(1-cos(dphi)));
	if(mtw > 160) continue;
	if(*met < 50) continue;
      }
      else if (sample == Sample::zee){
	if(*nelectrons != 2) continue;
	if(not ((*el1pt > 40 and *el1id == 1) or (*el2pt > 40 and *el2id == 1))) continue;
	if(*el1pid == *el2pid) continue;
	if(*zeemass < 60 or *zeemass > 120) continue;
      }

      efficiencyMonojetSelections->SetBinContent(8,efficiencyMonojetSelections->GetBinContent(8)+1);

      if((sample == Sample::wmn or sample == Sample::zmm) and *nelectrons > 0 ) continue;
      else if((sample == Sample::wen or sample == Sample::zee) and *nmuons > 0 ) continue;
      efficiencyMonojetSelections->SetBinContent(9,efficiencyMonojetSelections->GetBinContent(9)+1);

      if((sample == Sample::wmn or sample == Sample::zmm) and fabs(*metpf-*metcalo)/(*mmet) > 0.5) continue;
      else if((sample == Sample::wen or sample == Sample::zee) and fabs(*metpf-*metcalo)/(*emet) > 0.5) continue;
      efficiencyMonojetSelections->SetBinContent(10,efficiencyMonojetSelections->GetBinContent(10)+1);

      if((sample == Sample::wmn or sample == Sample::zmm) and applyJetSelections and *jmmdphi < 0.5) continue;      
      else if((sample == Sample::wen or sample == Sample::zee) and  applyJetSelections and *jemdphi < 0.5) continue;
      efficiencyMonojetSelections->SetBinContent(11,efficiencyMonojetSelections->GetBinContent(11)+1);

      float metVal = *mmet;
      if(sample == Sample::wen or sample == Sample::zee)
	metVal = *met;

      // MONOJET
      if(applyJetSelections and 
	 jetpt->at(0) > 100 and 
	 fabs(jeteta->at(0)) < 2.5 and jetchfrac->at(0) > 0.1 and jetnhfrac->at(0) < 0.8){
	
	// monojet vs recoil and met
	hden_monojet_recoil->Fill(metVal);
	efficiencyMonojetSelections->SetBinContent(12,efficiencyMonojetSelections->GetBinContent(12)+1);
	// numerator monojet
	if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm){
	  hnum_monojet_recoil->Fill(metVal);	
	  efficiencyMonojetSelections->SetBinContent(13,efficiencyMonojetSelections->GetBinContent(13)+1);	 
	}
      }
      else if(not applyJetSelections){
	// monojet vs recoil and met
	hden_monojet_recoil->Fill(metVal);
	efficiencyMonojetSelections->SetBinContent(12,efficiencyMonojetSelections->GetBinContent(12)+1);
	// numerator monojet
	if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm){
	  hnum_monojet_recoil->Fill(metVal);	
	  efficiencyMonojetSelections->SetBinContent(13,efficiencyMonojetSelections->GetBinContent(13)+1);	 
	}	
      }      
  }
    
  cout<<endl;
  
  //////////////
  for(int iBin = 0; iBin < efficiencyMonojetSelections->GetNbinsX(); iBin++){
    efficiencyMonojetSelections->SetBinContent(iBin+1,efficiencyMonojetSelections->GetBinContent(iBin+1)/double(nEvents));
    cout<<"Monojet selection efficiency: "<<efficiencyMonojetSelections->GetXaxis()->GetBinLabel(iBin+1)<<" eff "<<efficiencyMonojetSelections->GetBinContent(iBin+1)<<endl;
  }


  TEfficiency* eff_monojet_recoil = new TEfficiency(*hnum_monojet_recoil,*hden_monojet_recoil);
  eff_monojet_recoil->SetMarkerColor(kBlack);
  eff_monojet_recoil->SetLineColor(kBlack);
  eff_monojet_recoil->SetMarkerStyle(20);
  eff_monojet_recoil->SetMarkerSize(1);
  fitfunc_monojet_recoil->SetLineColor(kBlue);
  fitfunc_monojet_recoil->SetLineWidth(2);
  
  plotTurnOn(canvas,eff_monojet_recoil,fitfunc_monojet_recoil,"Recoil [GeV]",postfix,outputDIR);

  outputFile->cd();
  eff_monojet_recoil->Write();
  fitfunc_monojet_recoil->Write();
  outputFile->Close();
  
}

void plotTurnOn(TCanvas* canvas, TEfficiency* eff, TF1* fitfunc, const string & axisLabel, const string & postfix, const string & outputDIR){
  
  
  TH1* frame = canvas->DrawFrame(fitfunc->GetXaxis()->GetXmin(),0.,fitfunc->GetXaxis()->GetXmax(), 1.1, "");
  frame->GetXaxis()->SetTitle(axisLabel.c_str());  
  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetYaxis()->SetLabelSize(0.85*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetLabelSize(0.85*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetTitleSize(0.95*frame->GetYaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleSize(0.95*frame->GetXaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleOffset(1.1);
  frame->GetYaxis()->SetTitleOffset(1.2);

  // make the fit
  TGraphAsymmErrors* graph = eff->CreateGraph();
  TFitResultPtr fitResult = graph->Fit(fitfunc,"RS");
  int npoints       = 350;                                                                                                                                                                   
  TH1F* error_band  = new TH1F(Form("%s_error_band",fitfunc->GetName()),"",npoints,fitfunc->GetXaxis()->GetXmin(),fitfunc->GetXaxis()->GetXmax());
  if(drawUncertaintyBand)
    (TVirtualFitter::GetFitter())->GetConfidenceIntervals(error_band,0.68);
  
  canvas->SetRightMargin(0.075);
  canvas->SetTopMargin(0.06);
  canvas->Draw();
  canvas->cd();
  frame->Draw();
  if(drawUncertaintyBand){    
    error_band->SetFillColor(kBlue);
    error_band->SetFillStyle(3001);
    error_band->Draw("e4 same");  
  }
  eff->Draw("E1PSAME");
  fitfunc->Draw("SAME");
  
  canvas->RedrawAxis();
  CMS_lumi(canvas,"35.9",true);

  canvas->SaveAs((outputDIR+"/metTriggerEfficiency_"+string(postfix)+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/metTriggerEfficiency_"+string(postfix)+".pdf").c_str(),"pdf");

}

