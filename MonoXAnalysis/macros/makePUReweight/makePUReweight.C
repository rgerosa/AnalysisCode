#include "../CMS_lumi.h"

static int nbins     = 40;
static float nvtxMin = 0.;
static float nvtxMax = 40.;

void makeHisto(TH1F* hist, TChain* chain, bool isData,string controlregion, float lumi){

  TTreeReader reader(chain);

  // scale factors
  TFile sffile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF/leptonIDsfs.root");
  TH2*  msflhist = (TH2*)sffile.Get("muon_loose_SF");
  TH2*  msfthist = (TH2*)sffile.Get("muon_tight_SF");
  TH2*  esflhist = (TH2*)sffile.Get("electron_veto_SF");
  TH2*  esfthist = (TH2*)sffile.Get("electron_tight_SF");

  // trigger
  TFile trefile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF/leptonTrigsfs.root");
  TH2*  trehist = (TH2*)trefile.Get("hltel27_SF");

  const char* wgtsumvar;
  const char* wgtpuvar;
  const char* wgtbtagvar;

  if (not isData)   {
    wgtsumvar  = "wgtsum";
    wgtpuvar   = "wgtpileup";
    wgtbtagvar = "wgtbtag";
  }
  else {
    wgtsumvar  = "wgt";
    wgtpuvar   = "wgt";
    wgtbtagvar = "wgt";
  }

  TTreeReaderValue<unsigned char>   fcsc   (reader, "flagcsctight");
  TTreeReaderValue<unsigned char>   fhbhe  (reader, "flaghbhenoise");
  TTreeReaderValue<unsigned char>   fhbhei (reader, "flaghbheiso");
  TTreeReaderValue<unsigned char>   feesc  (reader, "flageebadsc");
  TTreeReaderValue<unsigned char>   fecal  (reader, "flagecaltp");
  TTreeReaderValue<unsigned char>   fnvtx  (reader, "flaggoodvertices");

  TTreeReaderValue<unsigned char>   hmnm90 (reader, "hltmet90");
  TTreeReaderValue<unsigned char>   hmnm120(reader, "hltmet120");
  TTreeReaderValue<unsigned char>   hmwm90 (reader, "hltmetwithmu90");
  TTreeReaderValue<unsigned char>   hmwm120(reader, "hltmetwithmu120");
  TTreeReaderValue<unsigned char>   hmwm170(reader, "hltmetwithmu170");
  TTreeReaderValue<unsigned char>   hmwm300(reader, "hltmetwithmu300");
  TTreeReaderValue<unsigned char>   hsmu   (reader, "hltsinglemu");
  TTreeReaderValue<unsigned char>   hsele  (reader, "hltsingleel");
  TTreeReaderValue<unsigned char>   hph175  (reader, "hltphoton175");
  TTreeReaderValue<unsigned char>   hph165  (reader, "hltphoton165");

  TTreeReaderValue<unsigned int>    nvtx   (reader, "nvtx");
  TTreeReaderValue<unsigned int>    nbjets (reader, "nbjetslowpt");
  TTreeReaderValue<double>          wgtsum (reader, wgtsumvar);
  TTreeReaderValue<double>          wgtpu  (reader, wgtpuvar);
  TTreeReaderValue<double>          wgtbtag (reader, wgtbtagvar);
  TTreeReaderValue<double>          wgt    (reader, "wgt");
  TTreeReaderValue<double>          xsec   (reader, "xsec");
  TTreeReaderValue<int>             mu1pid (reader, "mu1pid");
  TTreeReaderValue<int>             mu2pid (reader, "mu2pid");
  TTreeReaderValue<int>             mu1id  (reader, "mu1id");
  TTreeReaderValue<int>             mu2id  (reader, "mu2id");
  TTreeReaderValue<double>          mu1pt  (reader, "mu1pt");
  TTreeReaderValue<double>          mu1eta (reader, "mu1eta");
  TTreeReaderValue<double>          mu2pt  (reader, "mu2pt");
  TTreeReaderValue<double>          mu2eta (reader, "mu2eta");

  TTreeReaderValue<int>             el1pid (reader, "el1pid");
  TTreeReaderValue<int>             el2pid (reader, "el2pid");
  TTreeReaderValue<int>             el1id  (reader, "el1id");
  TTreeReaderValue<int>             el2id  (reader, "el2id");
  TTreeReaderValue<double>          el1pt  (reader, "el1pt");
  TTreeReaderValue<double>          el1eta (reader, "el1eta");
  TTreeReaderValue<double>          el2pt  (reader, "el2pt");
  TTreeReaderValue<double>          el2eta (reader, "el2eta");

  TTreeReaderValue<double>          t1pfmet(reader, "t1pfmet");
  TTreeReaderValue<double>          t1mumet(reader, "t1mumet");
  TTreeReaderValue<double>          t1elmet(reader, "t1elmet");
  TTreeReaderValue<double>          t1phmet(reader, "t1phmet");

  while (reader.Next()) {
    double weight = 1.0;
    double kfact  = 1.0;
    double puwgt  = 1.0;
    double effsf  = 1.0;
    double trgsf  = 1.0;

    if (isData && *fcsc    == 0) continue;
    if (isData && *fhbhe   == 0) continue;
    if (isData && *fhbhei  == 0) continue;
    if (isData && *feesc   == 0) continue;
    if (isData && *fecal   == 0) continue;
    if (isData && *fnvtx   == 0) continue;

    if (*nbjets > 0) continue;

    if(controlregion == "wmn"){
      unsigned char hlt = (*hmnm90) + (*hmnm120) + (*hmwm90) +(*hmwm120) + (*hmwm170) +(*hmwm300) + (*hsmu);
      if (isData and hlt == 0) continue;
      effsf *= msfthist->GetBinContent(msfthist->FindBin(min(999., *mu1pt), fabs(*mu1eta)));
    }
    
    if (controlregion == "zmm") {
      unsigned char hlt = (*hmnm90) + (*hmnm120) + (*hmwm90) +(*hmwm120) + (*hmwm170) +(*hmwm300) +(*hsmu);
      if (isData and hlt == 0) continue;
      bool istight = false;
      if (*mu1id == 1 && *mu1pt > 20.) istight = true;
      if (*mu2id == 1 && *mu2pt > 20.) istight = true;
      if (!istight) continue;
      if (*mu1id == 1) effsf *= msfthist->GetBinContent(msfthist->FindBin(min(999., *mu1pt), fabs(*mu1eta)));
      else             effsf *= msflhist->GetBinContent(msflhist->FindBin(min(999., *mu1pt), fabs(*mu1eta)));
      if (*mu2id == 1) effsf *= msfthist->GetBinContent(msfthist->FindBin(min(999., *mu2pt), fabs(*mu2eta)));
      else             effsf *= msflhist->GetBinContent(msflhist->FindBin(min(999., *mu2pt), fabs(*mu2eta)));
    }

    if (controlregion == "wen") {
      unsigned char hlt = *hsele;
      if (isData and hlt == 0) continue;
      if (*t1pfmet < 50) continue;
      effsf *= esfthist->GetBinContent(esfthist->FindBin(min(999., *el1pt), fabs(*el1eta)));
      trgsf *= trehist ->GetBinContent(trehist ->FindBin(min(999., *el1pt), fabs(*el1eta)));
    }
    
    if (controlregion == "zee") {
      unsigned char hlt = (*hsele) + (*hph165) + (*hph175);
      if (isData and hlt == 0) continue;
      bool istight = false;
      if (*el1id == 1 && *el1pt > 40.) istight = true;
      if (*el2id == 1 && *el2pt > 40.) istight = true;
      if (!istight) continue;
      if (*el1id == 1) effsf *= esfthist->GetBinContent(esfthist->FindBin(min(999., *el1pt), fabs(*el1eta)));
      else             effsf *= esflhist->GetBinContent(esflhist->FindBin(min(999., *el1pt), fabs(*el1eta)));
      if (*el2id == 1) effsf *= esfthist->GetBinContent(esfthist->FindBin(min(999., *el2pt), fabs(*el2eta)));
      else             effsf *= esflhist->GetBinContent(esflhist->FindBin(min(999., *el2pt), fabs(*el2eta)));
    }

    if (not isData) weight = (*xsec)*lumi*(*wgt)*(kfact)*(puwgt)*(trgsf)*(effsf)*(*wgtbtag)/(*wgtsum);
    hist->Fill(*nvtx, weight);   
  }

  sffile .Close();
  trefile.Close();
}

// k-factors not needed sicne done with inclusive samples
void makePUReweight(string baseDIR, string controlregion, string outputDIR, float lumi){

  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  system(("mkdir -p "+outputDIR).c_str());

  string fileNameData;

  if( controlregion == "zmm" or controlregion == "wmn")
    fileNameData = baseDIR+"/SingleMuon/"+controlregion+"filter";
  else if(controlregion == "zee" or controlregion == "wen")
    fileNameData = baseDIR+"/SingleElectron/"+controlregion+"filter";
  else{
    cerr<<"Control region is not adapt for nvtx re-weight calulation --> exit"<<endl;
    return;
  }

  TChain* chainData    = new TChain("tree/tree");
  TChain* chainZLL     = new TChain("tree/tree");
  TChain* chainWLN     = new TChain("tree/tree");
  TChain* chainTop     = new TChain("tree/tree");
  TChain* chainDiBoson = new TChain("tree/tree");

  chainData->Add((fileNameData+"/*root").c_str());
  chainZLL->Add((baseDIR+"/DYJets/"+controlregion+"filter/*root").c_str());
  chainWLN->Add((baseDIR+"/WJets/"+controlregion+"filter/*root").c_str());
  chainTop->Add((baseDIR+"/Top/"+controlregion+"filter/*root").c_str());
  chainDiBoson->Add((baseDIR+"/DiBoson/"+controlregion+"filter/*root").c_str());
  

  TH1F*  dathist = new TH1F("dathist","",nbins,nvtxMin,nvtxMax);
  TH1F*  zllhist = new TH1F("zllhist","",nbins,nvtxMin,nvtxMax);
  TH1F*  wlnhist = new TH1F("wlnhist","",nbins,nvtxMin,nvtxMax);
  TH1F*  tophist = new TH1F("tophist","",nbins,nvtxMin,nvtxMax);
  TH1F*  dibhist = new TH1F("dibhist","",nbins,nvtxMin,nvtxMax);

  dathist->Sumw2();
  zllhist->Sumw2();
  wlnhist->Sumw2();
  tophist->Sumw2();
  dibhist->Sumw2();
  
  makeHisto(dathist,chainData,true,controlregion,lumi);
  makeHisto(zllhist,chainZLL,false,controlregion,lumi);
  makeHisto(wlnhist,chainWLN,false,controlregion,lumi);  
  makeHisto(tophist,chainTop,false,controlregion,lumi);
  makeHisto(dibhist,chainDiBoson,false,controlregion,lumi);

  TCanvas* nvtx = new TCanvas("nvtx","");
  nvtx->cd();
  dathist->SetMarkerColor(kBlack);
  dathist->SetLineColor(kBlack);
  dathist->SetMarkerSize(1);
  dathist->SetMarkerStyle(20);
  dathist->GetXaxis()->SetTitle("N_{PV}");
  dathist->GetYaxis()->SetTitle("Events");
  
  dathist->Draw("EP");

  TH1F* totalMC = (TH1F*) zllhist->Clone("totalMC");

  totalMC->Add(wlnhist);
  totalMC->Add(tophist);
  totalMC->Add(dibhist);

  cout<<"Data entries "<<dathist->Integral()<<" total MC "<<totalMC->Integral()<<" zll "<<zllhist->Integral()<<" wln "<<wlnhist->Integral()<<" top "<<tophist->Integral()<<" diboson "<<dibhist->Integral()<<endl;
    
  totalMC->SetFillColor(kAzure-9);
  totalMC->SetFillStyle(3001);
  totalMC->SetLineColor(kBlack);  
  totalMC->Draw("hist same");
  
  dathist->Draw("EP same");

  nvtx->SaveAs((outputDIR+"/nvtx_"+controlregion+".pdf").c_str(),"pdf");
  nvtx->SaveAs((outputDIR+"/nvtx_"+controlregion+".png").c_str(),"pdf");
  
  TFile* output = new TFile((outputDIR+"/purwt"+controlregion+".root").c_str(),"RECREATE");
  output->cd();
  dathist->Scale(1./dathist->Integral());
  totalMC->Scale(1./totalMC->Integral());

  nvtx->SaveAs((outputDIR+"/nvtx_shape_"+controlregion+".pdf").c_str(),"pdf");
  nvtx->SaveAs((outputDIR+"/nvtx_shape_"+controlregion+".png").c_str(),"pdf");

  TH1F* puhist = (TH1F*) dathist->Clone("puhist");
  puhist->Divide(totalMC);
  puhist->Write("puhist");
  output->Close();
}
