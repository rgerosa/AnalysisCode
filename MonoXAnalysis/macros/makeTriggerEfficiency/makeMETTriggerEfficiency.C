#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

 const bool  drawUncertaintyBand   = false;


/// plotting result
void plotTurnOn(TFile* outputFile, TCanvas* canvas, TEfficiency* eff, TF1* fitfunc, const string & axisLabel, const TString & postfix, const string & ouputDIR, const float  & luminosity, 
    const bool & singleMuon, const TString & banner = "");

/// confidence interval for the band
void GetConfidenceIntervals(TF1* funz, TH1F* obj, Double_t cl, const TFitResultPtr & fitResult);

/// main function
void makeMETTriggerEfficiency(string inputDIR, string ouputDIR, float luminosity = 0.81, bool isMuon = true, bool doubleLepton = false) {


 // recoil binning
 vector <float> bins_monojet_recoil   = {0.,50.,60.,70.,80.,85.,95.,100., 110., 120., 130., 140., 150., 160., 170., 180., 190., 200., 225., 250., 275., 300., 325., 350., 400., 450., 500., 550., 650., 800., 1000., 1250};

 // eras
 vector<string> RunEra = {"Run2017B"};

 const bool  useDoubleMuonTriggers = true;
 const bool  applyJetSelections    = true;


  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+ouputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();

  // input tree
  TChain* tree = new TChain("tree/tree");
  // use only a subset of directories
  if(isMuon)
    system(("ls "+inputDIR+"  | grep SingleMu > list_dir.txt").c_str());
  else
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
      //system("rm list.txt");
    }
  }
  //system("rm list_dir.txt");

  // fitting function for the turn-on  as a function of recoil
  TF1 *fitfunc_monojet_recoil = new TF1("fitfunc_monojet_recoil",ErfCB,bins_monojet_recoil.front(), bins_monojet_recoil.back(),5);
  fitfunc_monojet_recoil->SetParameters(120., 25., 30., 4., 1.);

  TH1F* hnum_monojet_recoil = new TH1F("hnum_monojet_recoil", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  TH1F* hden_monojet_recoil = new TH1F("hden_monojet_recoil", "", bins_monojet_recoil.size()-1, &bins_monojet_recoil[0]);
  hnum_monojet_recoil->Sumw2();
  hden_monojet_recoil->Sumw2();

  // define numerator as event with tight muon + trigger requirement
  // define denominator as an event with a tight muon passing single muon trigger
  TTreeReader reader(tree);
  TTreeReaderValue<unsigned int> run    (reader,"run");
  TTreeReaderValue<unsigned int> lumi   (reader,"lumi");
  TTreeReaderValue<unsigned int> event  (reader,"event");
  TTreeReaderValue<UChar_t> hltm90     (reader,"hltmet90");
  TTreeReaderValue<UChar_t> hltm100    (reader,"hltmet100");
  TTreeReaderValue<UChar_t> hltm110    (reader,"hltmet110");
  TTreeReaderValue<UChar_t> hltm120    (reader,"hltmet120");
  TTreeReaderValue<UChar_t> hltmwm120  (reader,"hltmetwithmu120");
  TTreeReaderValue<UChar_t> hltmwm130  (reader,"hltmetwithmu130");
  TTreeReaderValue<UChar_t> hltmwm140  (reader,"hltmetwithmu140");
  //TTreeReaderValue<UChar_t> hltmwm170  (reader,"hltmetwithmu170");
  //TTreeReaderValue<UChar_t> hltmwm300  (reader,"hltmetwithmu300");
  TTreeReaderValue<UChar_t> hltmwm90   (reader,"hltmetwithmu90");
  TTreeReaderValue<UChar_t> hltjm      (reader,"hltjetmet");
  TTreeReaderValue<UChar_t> hltsinglemu (reader,"hltsinglemu");
  TTreeReaderValue<UChar_t> hltdoublemu (reader,"hltdoublemu");
  TTreeReaderValue<UChar_t> hlte      (reader,"hltsingleel");
  TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");                                                                                                                                       
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");                                                                                                                                       
  TTreeReaderValue<unsigned int> ntaus       (reader,"ntausold");
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
  TTreeReaderValue<float> metpf       (reader,"pfmet");
  TTreeReaderValue<float> metcalo     (reader,"calomet");
  TTreeReaderValue<float> jmmdphi (reader,"incjetmumetdphimin4");
  TTreeReaderValue<float> jemdphi (reader,"incjetelmetdphimin4");
 
    TString lepton;
    if (isMuon) lepton = "mu";
    else        lepton = "el";

    TTreeReaderValue<int>     lept1id     (reader,lepton+"1id");
    TTreeReaderValue<unsigned int> nleptons;
    if(isMuon)      nleptons = TTreeReaderValue<unsigned int>(reader,"nmuons");
    else            nleptons = TTreeReaderValue<unsigned int>(reader,"nelectrons");
    TTreeReaderValue<float>   lep1pt     (reader,lepton+"1pt");
    TTreeReaderValue<float>   lep1phi    (reader,lepton+"1phi");
    TTreeReaderValue<float>   lep1eta    (reader,lepton+"1eta");
    TTreeReaderValue<int>     lep1id     (reader,lepton+"1id");
    TTreeReaderValue<int>     lep1pid     (reader,lepton+"1pid");
    TTreeReaderValue<float>   lep2pt     (reader,lepton+"2pt");
    TTreeReaderValue<float>   lep2eta    (reader,lepton+"2eta");
    TTreeReaderValue<float>   lep2phi    (reader,lepton+"2phi");
    TTreeReaderValue<int>     lep2id     (reader,lepton+"2id");
    TTreeReaderValue<int>     lep2pid     (reader,lepton+"2pid");
    TTreeReaderValue<UChar_t> hltsinglelepton (reader,"hltsingle"+lepton);
    TTreeReaderValue<UChar_t> hltdoublelepton (reader,"hltdouble"+lepton);
    TTreeReaderValue<float> lepmet        (reader,"t1"+lepton+"met");
    TTreeReaderValue<float> lepmetphi     (reader,"t1"+lepton+"metphi");

  //////////////////
  long int nTotal = tree->GetEntries();
  cout<<"Total number of events: "<<nTotal<<endl;
  long int nEvents = 0;

  TH1D* efficiencyMonojetSelections = new TH1D("efficiencyMonojetSelections","efficiencyMonojetSelections",16,0,17);
  efficiencyMonojetSelections->Sumw2();
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(1,"Total");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(2,"B-veto");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(3,"Tau-veto");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(4,"Photon-veto");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(5,"MET filters");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(6,"HLT single lepton");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(7,"Lepton pT and eta");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(8,"Lepton ID");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(9,"Lepton Veto");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(10,"Calo-PF MET");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(11,"Jet met dphi");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(12,"Jet cuts");
  efficiencyMonojetSelections->GetXaxis()->SetBinLabel(13,"Trigger");

    std::map<TString, double> lepton_cuts;
    lepton_cuts["el_pt"] = 40.;
    lepton_cuts["el_eta"] = 2.5;
    lepton_cuts["mu_pt"] = 20.;
    lepton_cuts["mu_eta"] = 2.4;

    double lepton_met_recoil = 0;

    long int nPart = 100000;
    while(reader.Next()){
    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;

    lepton_met_recoil = *lepmet;
    if(!isMuon) 
      lepton_met_recoil = *met;
    //hden_monojet_recoil->Fill(*lepmet);
    // define the denominator
    efficiencyMonojetSelections->SetBinContent(1,efficiencyMonojetSelections->GetBinContent(1)+1);

    if(*nbjets != 0) continue;
    efficiencyMonojetSelections->SetBinContent(2,efficiencyMonojetSelections->GetBinContent(2)+1);

    if(*ntaus != 0) continue;
    efficiencyMonojetSelections->SetBinContent(3,efficiencyMonojetSelections->GetBinContent(3)+1);

    if(*nphotons  != 0) continue;
    efficiencyMonojetSelections->SetBinContent(4,efficiencyMonojetSelections->GetBinContent(4)+1);

    if(not *fcsc)  continue;
    if(not *feeb)  continue;
    if(not *fetp)  continue;
    if(not *fvtx)  continue;
    if(not *fbadmu) continue;
    if(not *fbadch) continue; 
    if(not *fhbhe)  continue;
    if(not *fhbiso) continue;

    efficiencyMonojetSelections->SetBinContent(5,efficiencyMonojetSelections->GetBinContent(5)+1);

      // trigger selection for the denominator
      if(not doubleLepton and not *hltsinglelepton) continue;
      else if(doubleLepton and not useDoubleMuonTriggers and not *hltsinglelepton) continue;
      else if(doubleLepton and useDoubleMuonTriggers and (not *hltdoublelepton or not *hltsinglelepton) ) continue;
            
      efficiencyMonojetSelections->SetBinContent(6,efficiencyMonojetSelections->GetBinContent(6)+1);

      if(*lep1pt < lepton_cuts[(lepton+"_pt")] or fabs(*lep1eta) > lepton_cuts[(lepton+"_eta")] ) continue;

      efficiencyMonojetSelections->SetBinContent(7,efficiencyMonojetSelections->GetBinContent(7)+1);



      if(not doubleLepton  and *lept1id != 1) continue;      
      if(not doubleLepton  and *nleptons != 1) continue;      
      else if(doubleLepton and *nleptons != 2) continue;    

      if(doubleLepton){
          if(*lep1pid == *lep2pid) continue; //opposite charge
          TLorentzVector lept1, lept2;
          lept1.SetPtEtaPhiM(*lep1pt,*lep1eta,*lep1phi,0.);
          lept2.SetPtEtaPhiM(*lep2pt,*lep2eta,*lep2phi,0.);
          if((*lep1pt > lepton_cuts[(lepton+"_pt")] and *lep1id != 1) or (*lep2pt > 20 and *lep2id != 1)) continue;
          if((lept1+lept2).M() < 60 or (lept1+lept2).M() > 120) continue;
      }

      efficiencyMonojetSelections->SetBinContent(8,efficiencyMonojetSelections->GetBinContent(8)+1);

      if(isMuon and *nelectrons > 0 ) continue;
      if(!isMuon and *nmuons > 0 ) continue;
      efficiencyMonojetSelections->SetBinContent(9,efficiencyMonojetSelections->GetBinContent(9)+1);

      if(fabs(*metpf-*metcalo)/(*lepmet) > 0.5) continue;
      efficiencyMonojetSelections->SetBinContent(10,efficiencyMonojetSelections->GetBinContent(10)+1);

      if(applyJetSelections and *jmmdphi < 0.5 and isMuon) continue;      
      if(applyJetSelections and *jemdphi < 0.5 and !isMuon) continue;      
      efficiencyMonojetSelections->SetBinContent(11,efficiencyMonojetSelections->GetBinContent(11)+1);

      // transverse mass cut
      if(not doubleLepton){
          float dphi = fabs(*lep1phi-*lepmetphi);
          if(dphi > TMath::Pi())
            dphi = 2*TMath::Pi()-dphi;
          float mtw = sqrt(2*(*lep1pt)*(lepton_met_recoil)*(1-cos(dphi)));
          if(mtw > 160) continue;
          //if(*met < 50) continue;  ///electron?
      }
      
       // MONOJET
       if(applyJetSelections and (jetpt->at(0) < 100. or fabs(jeteta->at(0)) > 2.5 or jetchfrac->at(0) < 0.1 or jetnhfrac->at(0) > 0.8 or jetchfrac->at(0) > 0.997)) continue;
  
       // monojet vs recoil and met
       hden_monojet_recoil->Fill(lepton_met_recoil);

       efficiencyMonojetSelections->SetBinContent(12,efficiencyMonojetSelections->GetBinContent(12)+1);
       // numerator monojet
       //if(*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm90 or *hltmwm100 or *hltmwm110 or *hltmwm170 or *hltmwm300 or *hltjm ){
       if(*hltm110 or *hltm120 or *hltmwm120 or *hltmwm130 or *hltmwm140 or *hltjm ){
         hnum_monojet_recoil->Fill(lepton_met_recoil); 
         efficiencyMonojetSelections->SetBinContent(13,efficiencyMonojetSelections->GetBinContent(13)+1);   
       }
      
       if(lepton_met_recoil > 300 and not (*hltm90 or *hltm100 or *hltm110 or *hltm120 or *hltmwm120 or *hltmwm130 or *hltmwm140 or *hltjm   )){
         cout<<"Monojet analysis: run "<<*run<<" lumi "<<*lumi<<" event "<<*event<<" t1lepmet "<<lepton_met_recoil<<" lepton pt "<<*lep1pt<<" lepton eta "<<*lep1eta<<" jet pt "<<jetpt->at(0)<<" eta "<<jeteta->at(0)<<" njets "<<*nincjets<<" pfmet "<<*met<<" calomet "<<*metcalo<<endl;
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
  TFile* outputFile = new TFile((ouputDIR+"/metTriggerEfficiency_"+lepton+"_recoil_monojet.root"),"RECREATE");

  if(isMuon)
    plotTurnOn(outputFile,canvas,eff_monojet_recoil,fitfunc_monojet_recoil,"Recoil [GeV]","recoil_monojet",ouputDIR,luminosity,isMuon);
  else
    plotTurnOn(outputFile,canvas,eff_monojet_recoil,fitfunc_monojet_recoil,"#slash{E}_{T} [GeV]","recoil_monojet",ouputDIR,luminosity,isMuon);

  outputFile->cd();
  efficiencyMonojetSelections->Write("selection");
  hnum_monojet_recoil->Write();
  hden_monojet_recoil->Write();
  outputFile->Close();

}

void plotTurnOn(TFile *outputFile,TCanvas* canvas, TEfficiency* eff, TF1* fitfunc, const string & axisLabel, const TString & postfix, const string & ouputDIR, const float  & luminosity, const bool & singleMuon, const TString & banner){


  TH1* frame = canvas->DrawFrame(fitfunc->GetXaxis()->GetXmin(),0.,fitfunc->GetXaxis()->GetXmax(), 1.1, "");
  frame->GetXaxis()->SetTitle(axisLabel.c_str());  
  frame->GetYaxis()->SetTitle("Trigger Efficiency");
  frame->GetYaxis()->SetLabelSize(0.8*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetTitleSize(0.8*frame->GetYaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleOffset(1.0);

  // make the fit
  TGraphAsymmErrors* graph = eff->CreateGraph();
  TFitResultPtr fitResult = graph->Fit(fitfunc,"RS","",100,2000);
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
  CMS_lumi(canvas,string(Form("%.2f",luminosity)),true);

  TLegend leg (0.6,0.3,0.9,0.4);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  leg.AddEntry((TObject*)0,banner,"");
  leg.Draw("same");

  TString lepton = "";
  if(!singleMuon) lepton = "ele_";

  canvas->SaveAs((ouputDIR+"/metTriggerEfficiency_"+lepton+string(postfix)+".png"),"png");
  canvas->SaveAs((ouputDIR+"/metTriggerEfficiency_"+lepton+string(postfix)+".pdf"),"pdf");

  outputFile->cd();
  // efficiency
  eff->Write("efficiency");
  fitfunc->Write("efficiency_func");

}

