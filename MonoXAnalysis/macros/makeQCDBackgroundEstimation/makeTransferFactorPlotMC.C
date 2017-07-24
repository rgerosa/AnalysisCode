#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void makePlots(TCanvas* canvas,TH1* histoA, TH1* histoB, TH1* histoC, TH1* histoD, string xAxisTitle, string yAxisTitle, string postfix, string outputDIR){

  histoA->GetXaxis()->SetTitle(xAxisTitle.c_str());
  histoB->GetXaxis()->SetTitle(xAxisTitle.c_str());
  histoC->GetXaxis()->SetTitle(xAxisTitle.c_str());
  histoD->GetXaxis()->SetTitle(xAxisTitle.c_str());
  histoA->GetYaxis()->SetTitle(yAxisTitle.c_str());
  histoB->GetYaxis()->SetTitle(yAxisTitle.c_str());
  histoC->GetYaxis()->SetTitle(yAxisTitle.c_str());
  histoD->GetYaxis()->SetTitle(yAxisTitle.c_str());
  
  histoA->SetLineColor(kBlack);
  histoB->SetLineColor(kBlue);
  histoC->SetLineColor(kRed);
  histoD->SetLineColor(kCyan+1);
  histoA->SetLineWidth(2);
  histoB->SetLineWidth(2);
  histoC->SetLineWidth(2);
  histoD->SetLineWidth(2);

  histoA->Draw("hist");
  histoB->Draw("hist same");
  histoC->Draw("hist same");
  histoD->Draw("hist same");

  TLegend* leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histoA,"E_{T}^{miss} [150,250] GeV, #Delta#phi < 0.5","L");
  leg->AddEntry(histoB,"E_{T}^{miss} > 250 GeV, #Delta#phi < 0.5","L");
  leg->AddEntry(histoC,"E_{T}^{miss} [150,250] GeV, #Delta#phi > 0.5","L");
  leg->AddEntry(histoD,"E_{T}^{miss} > 250 GeV, #Delta#phi > 0.5","L");
  leg->Draw("same");

  CMS_lumi(canvas,"35.9");

  histoA->GetYaxis()->SetRangeUser(min(histoA->GetMinimum(),min(histoB->GetMinimum(),min(histoC->GetMinimum(),histoD->GetMinimum())))*0.1,max(histoA->GetMaximum(),max(histoB->GetMaximum(),max(histoC->GetMaximum(),histoD->GetMaximum())))*100);
  
  if(histoA->GetMinimum() == 0) histoA->SetMinimum(0.001);
  
  canvas->SetLogy();  
  canvas->SaveAs((outputDIR+"/distributions_"+postfix+"_qcd_mc.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/distributions_"+postfix+"_qcd_mc.pdf").c_str(),"pdf");
}

//////
void makeRatioPlots(TCanvas* canvas,TH1* ratioDA, TH1* ratioDB, TH1* ratioDC, string xAxisTitle, string yAxisTitle, string postfix, string outputDIR){

  ratioDA->GetXaxis()->SetTitle(xAxisTitle.c_str());
  ratioDB->GetXaxis()->SetTitle(xAxisTitle.c_str());
  ratioDC->GetXaxis()->SetTitle(xAxisTitle.c_str());
  ratioDA->GetYaxis()->SetTitle(yAxisTitle.c_str());
  ratioDB->GetYaxis()->SetTitle(yAxisTitle.c_str());
  ratioDC->GetYaxis()->SetTitle(yAxisTitle.c_str());
  
  ratioDA->SetLineColor(kBlack);
  ratioDB->SetLineColor(kBlue);
  ratioDC->SetLineColor(kRed);
  ratioDA->SetLineWidth(2);
  ratioDB->SetLineWidth(2);
  ratioDC->SetLineWidth(2);

  ratioDA->Draw("hist");
  ratioDB->Draw("hist same");
  ratioDC->Draw("hist same");

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(ratioDA,"Ratio D/A","L");
  leg->AddEntry(ratioDB,"Ratio D/B","L");
  leg->AddEntry(ratioDC,"Ratio D/C","L");
  leg->Draw("same");

  CMS_lumi(canvas,"");

  ratioDA->GetYaxis()->SetRangeUser(min(ratioDA->GetMinimum(),min(ratioDB->GetMinimum(),ratioDC->GetMinimum()))*0.7,
				    max(ratioDA->GetMaximum(),max(ratioDB->GetMaximum(),ratioDC->GetMaximum()))*2.0);

  if(ratioDA->GetMinimum() == 0) ratioDA->SetMinimum(0.001);
  
  canvas->SaveAs((outputDIR+"/ratios_"+postfix+"_qcd_mc.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/ratios_"+postfix+"_qcd_mc.pdf").c_str(),"pdf");
  
}

void makeRatioFit(TCanvas* canvas, TH1* ratio, const string & outputDIR, const string & xAxisTitle, const string & yAxisTitle){

  canvas->cd();
  canvas->SetRightMargin(0.06);
  ratio->GetXaxis()->SetTitle(xAxisTitle.c_str());
  ratio->GetYaxis()->SetTitle(yAxisTitle.c_str());
  ratio->SetMarkerColor(kBlack);
  ratio->SetMarkerSize(1.0);
  ratio->SetMarkerStyle(20);
  ratio->Draw("EP");

  canvas->SetLogy(0);

  CMS_lumi(canvas,"");
  
  ratio->GetYaxis()->SetRangeUser(0,ratio->GetMaximum()*2.0);
  ratio->GetYaxis()->SetTitleOffset(1.10);
  ratio->GetXaxis()->SetTitleOffset(1.10);
  ratio->GetYaxis()->SetLabelSize(0.035);

  TF1* fitExp = new TF1("fitExp","[0]*TMath::Exp([1]*x)",ratio->GetXaxis()->GetXmin(),ratio->GetXaxis()->GetXmax());
  TFitResultPtr result = ratio->Fit(fitExp,"RMS");
  TF1* func = ratio->GetFunction("fitExp");
  
  func->SetLineColor(kBlue);
  func->SetLineWidth(2);
  ratio->Draw("EPsame");

  canvas->SaveAs((outputDIR+"/"+Form("%s",ratio->GetName())+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+Form("%s",ratio->GetName())+".pdf").c_str(),"pdf");
  
}


///////////
static bool doSmooth = true;
void makeTransferFactorPlotMC (string inputPath, string outputDIR, float lumi, Category category, bool isEOS){

  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  initializeBinning();
  gROOT->SetBatch(kTRUE);

  gStyle->SetOptStat(1000000000);
  gStyle->SetOptFit(1);  

  
  TChain* chain = new TChain("tree/tree");  

  if(isEOS)
    system(("/afs/cern.ch/project/eos/installation/cms/bin/eos.select find "+inputPath+" -name \"*.root\" | grep -v failed > file.temp").c_str());
  else
    system(("find "+inputPath+" -name \"*.root\" | grep -v failed > file.temp").c_str());

  ifstream infile;
  string line;
  infile.open("file.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "" or not TString(line).Contains(".root")) continue;
      if(isEOS)
	chain->Add(("root://eoscms.cern.ch//"+line).c_str());
      else
	chain->Add(line.c_str());
    }
  }
  infile.close();
  system("rm file.temp");

  TTreeReader reader(chain);
  TTreeReaderValue<UChar_t> hltPFHT350 (reader,"hltPFHT350");
  TTreeReaderValue<UChar_t> hltPFHT400 (reader,"hltPFHT400");
  TTreeReaderValue<UChar_t> hltPFHT475 (reader,"hltPFHT475");
  TTreeReaderValue<UChar_t> hltPFHT600 (reader,"hltPFHT600");
  TTreeReaderValue<UChar_t> hltPFHT650 (reader,"hltPFHT650");
  TTreeReaderValue<UChar_t> hltPFHT800 (reader,"hltPFHT800");
  TTreeReaderValue<UChar_t> hltPFHT900 (reader,"hltPFHT900");

  TTreeReaderValue<unsigned int> lumisection (reader,"lumi");
  TTreeReaderValue<unsigned int> event       (reader,"event");
  TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
  TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
  TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
  TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
  TTreeReaderValue<vector<float> > chfrac  (reader,"combinejetCHfrac");
  TTreeReaderValue<vector<float> > nhfrac  (reader,"combinejetNHfrac");
  TTreeReaderValue<float> xsec        (reader,"xsec");
  TTreeReaderValue<float> met         (reader,"t1pfmet");
  TTreeReaderValue<float> jmmdphi     (reader,"incjetmumetdphimin4");
  TTreeReaderValue<float> metcalo     (reader,"calomet");
  TTreeReaderValue<float> wgtpileup   (reader,"wgtpileup");
  TTreeReaderValue<float> wgtbtag     (reader,"wgtbtag");
  TTreeReaderValue<double> wgtsum      (reader,"wgtsum");
  TTreeReaderValue<float> wgt         (reader,"wgt");

  //// output file
  TFile* outputFile = new TFile((outputDIR+"/distributions_qcd_mc.root").c_str(),"RECREATE");
  outputFile->cd();

  vector<double> bins_mjj = {200.,500.,800.,1300.,2200.,5000.};
  
  TH1F* histoQCD_mjj_inclusive = new TH1F("histoQCD_mjj_inclusive","",bins_mjj.size()-1,&bins_mjj[0]);
  TH1F* histoQCD_mjj_regionA = new TH1F("histoQCD_mjj_regionA","",bins_mjj.size()-1,&bins_mjj[0]);
  TH1F* histoQCD_mjj_regionB = new TH1F("histoQCD_mjj_regionB","",bins_mjj.size()-1,&bins_mjj[0]);
  TH1F* histoQCD_mjj_regionC = new TH1F("histoQCD_mjj_regionC","",bins_mjj.size()-1,&bins_mjj[0]);
  TH1F* histoQCD_mjj_regionD = new TH1F("histoQCD_mjj_regionD","",bins_mjj.size()-1,&bins_mjj[0]);

  histoQCD_mjj_inclusive->Sumw2();
  histoQCD_mjj_regionA->Sumw2();
  histoQCD_mjj_regionB->Sumw2();
  histoQCD_mjj_regionC->Sumw2();
  histoQCD_mjj_regionD->Sumw2();

  vector<double> bins_ht = selectBinning("HT",category);
  TH1F* histoQCD_ht_regionA = new TH1F("histoQCD_ht_regionA","",bins_ht.size()-1,&bins_ht[0]);
  TH1F* histoQCD_ht_regionB = new TH1F("histoQCD_ht_regionB","",bins_ht.size()-1,&bins_ht[0]);
  TH1F* histoQCD_ht_regionC = new TH1F("histoQCD_ht_regionC","",bins_ht.size()-1,&bins_ht[0]);
  TH1F* histoQCD_ht_regionD = new TH1F("histoQCD_ht_regionD","",bins_ht.size()-1,&bins_ht[0]);

  histoQCD_ht_regionA->Sumw2();
  histoQCD_ht_regionB->Sumw2();
  histoQCD_ht_regionC->Sumw2();
  histoQCD_ht_regionD->Sumw2();

  long int totalEntries = chain->GetEntries();
  long int nPart = 100000;
  long int nEvents = 0;
  double sumwgt  = 0;

  while(reader.Next()){

    // require one of the trigger to be accpted                                                                                                                                                      
    int passHT = *hltPFHT350+*hltPFHT400+*hltPFHT475+*hltPFHT600+*hltPFHT650+*hltPFHT800+*hltPFHT900;
    if( passHT == 0) continue;

    if(jetpt->size() < 2) continue;
    
    TLorentzVector jet1, jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));
    
    histoQCD_mjj_inclusive->Fill((jet1+jet2).M(),*xsec*lumi*(*wgt)*(*wgtpileup)*(*wgtbtag)/(*wgtsum));
    
    cout.flush();
    if(nEvents % nPart == 0) cout<<"\r"<<"Analyzing events "<<double(nEvents)/totalEntries*100<<" % ";
    nEvents++;

    if(fabs(*met-*metcalo)/(*met) > 0.5) continue;

    if(category == Category::VBFrelaxed){
      if(jetpt->at(0) < 80) continue;
      if(jetpt->at(1) < 40) continue;
      if(fabs(jeteta->at(0)) > 4.7) continue;
      if(fabs(jeteta->at(1)) > 4.7) continue;
      if(fabs(jeteta->at(0)-jeteta->at(1)) < 1.0) continue;
      if(fabs(jet1.DeltaPhi(jet2)) > 1.3) continue;
      if((jet1+jet2).M() < 200) continue;
      if(fabs(jeteta->at(0)) < 2.4 and chfrac->at(0) < 0.1) continue;
      if(fabs(jeteta->at(0)) < 2.4 and nhfrac->at(0) > 0.8) continue;
      if(jeteta->at(0)*jeteta->at(1) > 0 ) continue;
   
    }
    else if(category == Category::VBF){
      if(jetpt->at(0) < 80) continue;
      if(jetpt->at(1) < 40) continue;
      if(fabs(jeteta->at(0)) > 4.7) continue;
      if(fabs(jeteta->at(1)) > 4.7) continue;
      if(fabs(jeteta->at(0)-jeteta->at(1)) < 4.0) continue;
      if(fabs(jet1.DeltaPhi(jet2)) > 1.5) continue;
      if((jet1+jet2).M() < 1300) continue;
      if(fabs(jeteta->at(0)) < 2.4 and chfrac->at(0) < 0.1) continue;
      if(fabs(jeteta->at(0)) < 2.4 and nhfrac->at(0) > 0.8) continue;
      if(jeteta->at(0)*jeteta->at(1) > 0 ) continue;
    }

    float ht = 0;
    for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
      if(jetpt->at(ijet) > 30 and fabs(jeteta->at(ijet)) < 4.7)
	ht += jetpt->at(ijet);
    }

    if(category == Category::VBF or category == Category::VBFrelaxed){
      if(*met > 100 and *met < 160 and *jmmdphi < 0.5){
	histoQCD_mjj_regionA->Fill((jet1+jet2).M(),*xsec*lumi*(*wgt)*(*wgtpileup)*(*wgtbtag)/(*wgtsum));
	histoQCD_ht_regionA->Fill(ht,*xsec*lumi*(*wgt)*(*wgtpileup)*(*wgtbtag)/(*wgtsum));
      }
      
      else if(*met > 250 and *jmmdphi < 0.5){
	histoQCD_mjj_regionB->Fill((jet1+jet2).M(),*xsec*lumi*(*wgt)*(*wgtpileup)*(*wgtbtag)/(*wgtsum));
	histoQCD_ht_regionB->Fill(ht,*xsec*lumi*(*wgt)*(*wgtpileup)*(*wgtbtag)/(*wgtsum));
      }
      else if(*met > 100 and *met < 160 and *jmmdphi > 0.5){
	histoQCD_mjj_regionC->Fill((jet1+jet2).M(),*xsec*lumi*(*wgt)*(*wgtpileup)*(*wgtbtag)/(*wgtsum));
	histoQCD_ht_regionC->Fill(ht,*xsec*lumi*(*wgt)*(*wgtpileup)*(*wgtbtag)/(*wgtsum));
      }
      else if(*met > 250 and *jmmdphi > 0.5){
	histoQCD_mjj_regionD->Fill((jet1+jet2).M(),*xsec*lumi*(*wgt)*(*wgtpileup)*(*wgtbtag)/(*wgtsum));
	histoQCD_ht_regionD->Fill(ht,*xsec*lumi*(*wgt)*(*wgtpileup)*(*wgtbtag)/(*wgtsum));
      }      
    }
  }  
  cout<<endl;

  cout<<"Expected rate in region A: "<<histoQCD_mjj_regionA->Integral()<<endl;
  cout<<"Expected rate in region B: "<<histoQCD_mjj_regionB->Integral()<<endl;
  cout<<"Expected rate in region C: "<<histoQCD_mjj_regionC->Integral()<<endl;
  cout<<"Expected rate in region D: "<<histoQCD_mjj_regionD->Integral()<<endl;
  
  cout<<"Region A: efficiency "<<histoQCD_mjj_regionA->Integral()/(histoQCD_mjj_inclusive->Integral())*100<<"%"<<endl;
  cout<<"Region B: efficiency "<<histoQCD_mjj_regionB->Integral()/(histoQCD_mjj_inclusive->Integral())*100<<"%"<<endl;
  cout<<"Region C: efficiency "<<histoQCD_mjj_regionC->Integral()/(histoQCD_mjj_inclusive->Integral())*100<<"%"<<endl;
  cout<<"Region D: efficiency "<<histoQCD_mjj_regionD->Integral()/(histoQCD_mjj_inclusive->Integral())*100<<"%"<<endl;

  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();

  histoQCD_ht_regionA->GetXaxis()->SetNdivisions(505);
  histoQCD_ht_regionB->GetXaxis()->SetNdivisions(505);
  histoQCD_ht_regionC->GetXaxis()->SetNdivisions(505);
  histoQCD_ht_regionD->GetXaxis()->SetNdivisions(505);
  histoQCD_mjj_regionA->GetXaxis()->SetNdivisions(505);
  histoQCD_mjj_regionB->GetXaxis()->SetNdivisions(505);
  histoQCD_mjj_regionC->GetXaxis()->SetNdivisions(505);
  histoQCD_mjj_regionD->GetXaxis()->SetNdivisions(505);

  if(doSmooth){
    histoQCD_ht_regionA->Smooth();
    histoQCD_ht_regionB->Smooth();
    histoQCD_ht_regionC->Smooth();
    histoQCD_ht_regionD->Smooth();
    histoQCD_mjj_regionA->Smooth();
    histoQCD_mjj_regionB->Smooth();
    histoQCD_mjj_regionC->Smooth();
    histoQCD_mjj_regionD->Smooth();
  }

  // ratios D/A
  TH1* ratio_mjj_DoverA = (TH1*) histoQCD_mjj_regionD->Clone("ratio_mjj_DoverA");
  ratio_mjj_DoverA->Divide(histoQCD_mjj_regionA);
  // ratios D/B
  TH1* ratio_mjj_DoverB = (TH1*) histoQCD_mjj_regionD->Clone("ratio_mjj_DoverB");
  ratio_mjj_DoverB->Divide(histoQCD_mjj_regionB);

  // ratios D/C
  TH1* ratio_mjj_DoverC = (TH1*) histoQCD_mjj_regionD->Clone("ratio_mjj_DoverC");
  ratio_mjj_DoverC->Divide(histoQCD_mjj_regionC);

  // ratios B/A
  TH1* ratio_mjj_BoverA = (TH1*) histoQCD_mjj_regionB->Clone("ratio_mjj_BoverA");
  ratio_mjj_BoverA->Divide(histoQCD_mjj_regionA);

  // ratios C/A
  TH1* ratio_mjj_CoverA = (TH1*) histoQCD_mjj_regionC->Clone("ratio_mjj_CoverA");
  ratio_mjj_CoverA->Divide(histoQCD_mjj_regionA);

  outputFile->cd();  
  histoQCD_ht_regionA->Write("histoQCD_ht_regionA");
  histoQCD_ht_regionB->Write("histoQCD_ht_regionB");
  histoQCD_ht_regionC->Write("histoQCD_ht_regionC");
  histoQCD_ht_regionD->Write("histoQCD_ht_regionD");

  histoQCD_mjj_regionA->Write("histoQCD_mjj_regionA");
  histoQCD_mjj_regionB->Write("histoQCD_mjj_regionB");
  histoQCD_mjj_regionC->Write("histoQCD_mjj_regionC");
  histoQCD_mjj_regionD->Write("histoQCD_mjj_regionD");

  ratio_mjj_DoverA->Write("ratio_mjj_DoverA");
  ratio_mjj_DoverB->Write("ratio_mjj_DoverB");
  ratio_mjj_DoverC->Write("ratio_mjj_DoverC");

  ratio_mjj_BoverA->Write("ratio_mjj_BoverA");
  ratio_mjj_CoverA->Write("ratio_mjj_CoverA");

  makePlots(canvas,histoQCD_mjj_regionA,histoQCD_mjj_regionB,histoQCD_mjj_regionC,histoQCD_mjj_regionD,"M_{jj} [GeV]","Events","mjj",outputDIR);
  makePlots(canvas,histoQCD_ht_regionA,histoQCD_ht_regionB,histoQCD_ht_regionC,histoQCD_ht_regionD,"H_{T} [GeV]","Events","ht",outputDIR);

  makeRatioFit(canvas,ratio_mjj_DoverA,outputDIR,"M_{jj} [GeV]","Ratio");
  makeRatioFit(canvas,ratio_mjj_DoverB,outputDIR,"M_{jj} [GeV]","Ratio");
  makeRatioFit(canvas,ratio_mjj_DoverC,outputDIR,"M_{jj} [GeV]","Ratio");
  makeRatioFit(canvas,ratio_mjj_BoverA,outputDIR,"M_{jj} [GeV]","Ratio");
  makeRatioFit(canvas,ratio_mjj_CoverA,outputDIR,"M_{jj} [GeV]","Ratio");
  
  outputFile->Close();

}
