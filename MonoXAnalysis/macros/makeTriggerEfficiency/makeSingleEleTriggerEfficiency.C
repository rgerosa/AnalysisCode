#include <iostream>
#include <sstream>
#include <cmath>
#include "triggerUtils.h"
#include "../CMS_lumi.h"

vector<string> RunEra = {"Run2016B","Run2016C","Run2016D","Run2016E","Run2016F","Run2016G","Run2016H"};


// calculation from the jetHT dataset
void makeSingleElectronTriggerEfficiency(string inputDIR, string outputDIR, float lumi = 12.9, bool doFit = false) {

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1410065408);

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 600);
  canvas->cd();
  
  TF1 *fitfunc = new TF1("fitfunc", ErfCB, 20, 1200, 5);
  fitfunc->SetParameters(30., 5., 5., 4., 1.);

  vector<float> binsPt  = {100,150,200,250,300,350,400,500,600,700,800,1000,1200,1500};
  vector<float> binsEta = {0,0.75,1.5,2.5};  
  fitfunc->SetRange(binsPt.front(),binsPt.back());
  
  TChain* tree = new TChain("tree/tree");

  system(("ls "+inputDIR+"  | grep root > list_dir.txt").c_str());
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
  
  TH2F* hnum = new TH2F("hnum", "", binsEta.size()-1, &binsEta[0],binsPt.size()-1, &binsPt[0]);
  TH2F* hden = new TH2F("hden", "", binsEta.size()-1, &binsEta[0],binsPt.size()-1, &binsPt[0]);
  hnum->Sumw2();
  hden->Sumw2();

  TH2F* hnum_recover = new TH2F("hnum_recover", "", binsEta.size()-1, &binsEta[0],binsPt.size()-1, &binsPt[0]);
  TH2F* hden_recover = new TH2F("hden_recover", "", binsEta.size()-1, &binsEta[0],binsPt.size()-1, &binsPt[0]);
  hnum_recover->Sumw2();
  hden_recover->Sumw2();

  TTreeReader reader(tree);
  TTreeReaderValue<unsigned int> run    (reader,"run");
  TTreeReaderValue<unsigned int> event  (reader,"event");
  TTreeReaderValue<UChar_t> hlte        (reader,"hltsingleel27");
  TTreeReaderValue<UChar_t> hltenoiso   (reader,"hltelnoiso");
  TTreeReaderValue<UChar_t> hltPFHT400  (reader,"hltPFHT400");
  TTreeReaderValue<UChar_t> hltPFHT650  (reader,"hltPFHT650");
  TTreeReaderValue<UChar_t> hltPFHT800  (reader,"hltPFHT800");
  TTreeReaderValue<UChar_t> hltEcalPFHT800  (reader,"hltEcalHT800");
  TTreeReaderValue<float>   el1pt    (reader,"el1pt");
  TTreeReaderValue<float>   el1eta   (reader,"el1eta");
  TTreeReaderValue<float>   el1phi   (reader,"el1phi");
  TTreeReaderValue<int>     el1id    (reader,"el1id");
  TTreeReaderValue<int>     el1idt   (reader,"el1idt");
  TTreeReaderValue<UChar_t> fhbhe    (reader,"flaghbhenoise");
  TTreeReaderValue<UChar_t> fhbiso   (reader,"flaghbheiso");
  TTreeReaderValue<UChar_t> fcsc   (reader,"flagglobaltighthalo");
  TTreeReaderValue<UChar_t> fcsct  (reader,"flagcsctight");
  TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
  TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
  TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
  TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
  TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");
  TTreeReaderValue<unsigned int> ntausraw    (reader,"ntausrawold");
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
  TTreeReaderValue<float> emet        (reader,"t1elmet");
  TTreeReaderValue<float> emetphi     (reader,"t1elmetphi");
  TTreeReaderValue<float> metpf       (reader,"pfmet");
  TTreeReaderValue<float> metcalo     (reader,"calomet");
  TTreeReaderValue<float> jemdphi (reader,"incjetelmetdphimin4");
  //TTreeReaderValue<float> wemt   (reader,"wemt");
 
  long int nTotal = tree->GetEntries();
  cout<<"Total number of events: "<<nTotal<<endl;
  long int nEvents = 0;

  long int nPart = 100000;
  while(reader.Next()){
    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/nTotal*100<<" % ";
    nEvents++;
    if(*nmuons   != 0) continue;
    if(*ntausraw != 0) continue;
    if(*nphotons != 0) continue;
    if(*nbjets   != 0) continue;
    if(*nelectrons != 1) continue;
    // to kill QCD in a sample of Wen events
    if(*met < 50) continue;
    if(fabs(*el1eta) > 2.5) continue;
    if(*el1id  != 1)  continue;
    if(*el1idt != 1)  continue;
    if(*fhbhe == 0 or *fhbiso == 0 or *fcsc == 0 or *fcsct == 0 or *feeb == 0 or *fetp == 0 or *fvtx == 0 or *fbadmu == 0 or *fbadch == 0) continue;
    //if(*wemt > 160)    continue; 
    if(*jemdphi < 0.5) continue;
    if(fabs(*el1eta) > 1.445 and fabs(*el1eta) < 1.55) continue; //exclude gap
    if(fabs(*el1eta) > 2.2)  continue;
    if(jetpt->size() == 0)   continue;
    if(jetpt->at(0)  <  100) continue;
    if(fabs(jeteta->at(0)) > 2.5) continue;
    if(jetchfrac->at(0) < 0.1) continue;
    if(jetnhfrac->at(0) > 0.8) continue;

    if(*hltPFHT400 or *hltPFHT650 or *hltPFHT800){
      hden->Fill(fabs(*el1eta),*el1pt);
      if(*hlte or *hltenoiso){
	hnum->Fill(fabs(*el1eta),*el1pt);
      }
    }
    if(*hltPFHT400 or *hltPFHT650 or *hltPFHT800){
      hden_recover->Fill(fabs(*el1eta),*el1pt);
      if(*hlte or *hltenoiso or *hltEcalPFHT800)
	hnum_recover->Fill(fabs(*el1eta),*el1pt);
    }    
  }
  
  TEfficiency* efficiency         = new TEfficiency(*hnum,*hden);
  TEfficiency* efficiency_recover = new TEfficiency(*hnum_recover,*hden_recover);
  TH2* histoEff         = efficiency->CreateHistogram();
  TH2* histoEff_recover = efficiency_recover->CreateHistogram();
  // in order to plot the 2D histo
  fitfunc->SetLineColor(kBlue);
  fitfunc->SetLineWidth(2);
  
  gPad->SetBottomMargin(0.10);
  gPad->SetRightMargin(0.17);
  TH1* frame = canvas->DrawFrame(binsEta.front(),binsPt.front(), binsEta.back(), binsPt.back(), "");  
  frame->GetYaxis()->SetTitle("Electron p_{T} [GeV]");
  frame->GetXaxis()->SetTitle("Electron #eta");
  frame->GetZaxis()->SetTitle("Trigger Efficiency");
  frame->GetYaxis()->SetLabelSize(0.8*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame->GetZaxis()->SetLabelSize(0.8*frame->GetXaxis()->GetLabelSize());
  frame->GetYaxis()->SetTitleSize(0.8*frame->GetYaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame->GetZaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());
  frame->GetXaxis()->SetTitleOffset(1.0);  
  canvas->SetTopMargin(0.06);
  canvas->Draw();
  canvas->cd();
  canvas->SetLogy();
  gStyle->SetPaintTextFormat(".2f");
  efficiency->Draw("colztext same");
  canvas->RedrawAxis();
  CMS_lumi(canvas,string(Form("%.1f",lumi)),true);
  canvas->SaveAs((outputDIR+"/electronTriggerEff.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/electronTriggerEff.pdf").c_str(),"pdf");
  canvas->SetLogy(0);

  canvas->Draw();
  canvas->cd();
  canvas->SetLogy();
  gStyle->SetPaintTextFormat(".2f");
  efficiency_recover->Draw("colztext same");
  canvas->RedrawAxis();
  CMS_lumi(canvas,string(Form("%.1f",lumi)),true);
  canvas->SaveAs((outputDIR+"/electronTriggerEff_recover.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/electronTriggerEff_recover.pdf").c_str(),"pdf");
  canvas->SetLogy(0);

  TCanvas* canvas2 = new TCanvas("canvas2", "canvas2", 600, 625);
  canvas2->cd();
  gPad->SetRightMargin(0.06);
  for(int iBinX = 0; iBinX < histoEff->GetNbinsX(); iBinX++){
    TGraphAsymmErrors* projection_pt = new TGraphAsymmErrors();
    TGraphAsymmErrors* projection_pt_recover = new TGraphAsymmErrors();
    for(int iBinY = 0; iBinY < histoEff->GetNbinsY(); iBinY++){
      int globalBin = histoEff->GetBin(iBinX+1,iBinY+1);
      projection_pt->SetPoint(iBinY+1,histoEff->GetYaxis()->GetBinCenter(iBinY+1),efficiency->GetEfficiency(globalBin));
      projection_pt->SetPointError(iBinY+1,histoEff->GetYaxis()->GetBinWidth(iBinY+1)/2,histoEff->GetYaxis()->GetBinWidth(iBinY+1)/2,efficiency->GetEfficiencyErrorLow(globalBin),efficiency->GetEfficiencyErrorUp(globalBin));
    }

    for(int iBinY = 0; iBinY < histoEff_recover->GetNbinsY(); iBinY++){
      int globalBin = histoEff_recover->GetBin(iBinX+1,iBinY+1);
      projection_pt_recover->SetPoint(iBinY+1,histoEff_recover->GetYaxis()->GetBinCenter(iBinY+1),efficiency_recover->GetEfficiency(globalBin));
      projection_pt_recover->SetPointError(iBinY+1,histoEff_recover->GetYaxis()->GetBinWidth(iBinY+1)/2,histoEff_recover->GetYaxis()->GetBinWidth(iBinY+1)/2,efficiency_recover->GetEfficiencyErrorLow(globalBin),efficiency_recover->GetEfficiencyErrorUp(globalBin));
    }

    projection_pt->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
    projection_pt->GetXaxis()->SetRangeUser(binsPt.front(),binsPt.back());
    projection_pt->GetYaxis()->SetTitle("Trigger Efficiency");
    projection_pt->GetYaxis()->SetRangeUser(0.8,1.1);
    projection_pt->SetMarkerSize(1);
    projection_pt->SetMarkerStyle(20);
    projection_pt->SetMarkerColor(kRed);
    projection_pt->SetLineColor(kRed);
    projection_pt->Draw("AE1P");
    projection_pt_recover->SetMarkerSize(1);
    projection_pt_recover->SetMarkerStyle(20);
    projection_pt_recover->SetMarkerColor(kBlue);
    projection_pt_recover->SetLineColor(kBlue);
    projection_pt_recover->Draw("E1Psame");

    TLegend* leg = new TLegend(0.6,0.25,0.9,0.45);
    leg->SetFillStyle(0);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->AddEntry(projection_pt,"HLT_Ele27 || HLT_Ele_105","PE");
    leg->AddEntry(projection_pt_recover,"HLT_Ele27 || HLT_Ele_105 || Ecal_HT800","PE");
    leg->Draw("same");
    
    if(doFit)
      projection_pt->Fit(fitfunc);
    CMS_lumi(canvas2,string(Form("%.1f",lumi)),true);
    canvas2->SaveAs((outputDIR+"/"+string(histoEff->GetName())+string(Form("_projection_pt_eta_%.1f_%.1f.png",histoEff->GetXaxis()->GetBinLowEdge(iBinX+1),
									   histoEff->GetXaxis()->GetBinLowEdge(iBinX+2)))).c_str(),"png");
    canvas2->SaveAs((outputDIR+"/"+string(histoEff->GetName())+string(Form("_projection_pt_eta_%.1f_%.1f.pdf",histoEff->GetXaxis()->GetBinLowEdge(iBinX+1),
										 histoEff->GetXaxis()->GetBinLowEdge(iBinX+2)))).c_str(),"pdf");
  }
  
  TFile* outputFile = new TFile((outputDIR+"/triggerEfficiency_SingleElectron_jetHT.root").c_str(),"RECREATE");
  outputFile->cd();
  efficiency->Write("efficiency");
  efficiency_recover->Write("efficiency_recover");
}

