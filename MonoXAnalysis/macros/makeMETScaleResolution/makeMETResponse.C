#include "../CMS_lumi.h"

// met binning used in this study
vector<double> ZPT_bins  = {0.0, 5., 10., 15., 20., 30., 40., 60., 80., 100., 125, 150., 175., 200., 400.};
vector<double> Nvtx_bins = {0.0, 5.5, 10.5, 15.5, 20.5, 25.5, 30.5, 50.};
vector<double> ZEta_bins = {0.0, 0.5, 1.0,  1.5, 2.0, 2.5};

class analysisBin {

public:
  analysisBin(){};
  analysisBin(const string & met, const string & observable, const string & category, const vector<double> & bins = {}):
    met_(met),
    observable_(observable),
    category_(category),
    bins_(bins)
  {

    if(met_ != "pfmet" and met_ != "t1pfmet" and met_ != "puppit1pfmet" and met_ != "puppipfmet")
      cout<<"Un-recognized met option --> please check"<<endl;
    if(observable_ != "zpt" and observable_ != "nvtx" and observable_ != "zeta")
      cout<<"Un-recognized observable option --> please check"<<endl;
    if(category_ != "zmm" and category_ !=  "zee")
      cout<<"Un-recognized category option --> please check"<<endl;

    if(bins_.empty()){
      if(observable_ == "zpt")
        bins_ = ZPT_bins;
      else if(observable_ == "zeta")
        bins_ = ZEta_bins;
      else if(observable_ == "nvtx")
        bins_ = Nvtx_bins;
    }
  }
  ~analysisBin(){};
  vector<double> bins_;
  string met_;
  string category_;
  string observable_;
  size_t iBin_;
};


Double_t dphi(Double_t phi1, Double_t phi2) {
  double result = phi1 - phi2;
  while (result > M_PI) result -= 2*M_PI;
  while (result <= -M_PI) result += 2*M_PI;
  return result;
}

// Loop on the event in order to make the response in a specific met bin --> re-loop for each bin -> to be changed                                                           
double drawplot(TTree* tree,
                TH1* hist,
                const bool & isMC,
                const double & xmin, const double & xmax,
                const analysisBin & binning,
                const float & lumi) {


  TFile* pufile = new TFile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt.root");
  TH1*   puhist = (TH1*)pufile->Get("puhist");

  TFile sffile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF/leptonIDsfs.root");
  TH2*  msflhist = (TH2*)sffile.Get("muon_loose_SF");
  TH2*  msfthist = (TH2*)sffile.Get("muon_tight_SF");
  TH2*  esflhist = (TH2*)sffile.Get("electron_veto_SF");
  TH2*  esfthist = (TH2*)sffile.Get("electron_tight_SF");

  TFile trefile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF/leptonTrigsfs.root");
  TH2*  trehist = (TH2*)trefile.Get("hltel27_SF");

  TTreeReader reader(tree);
  const char* wgtsumvar;
  const char* wgtpuvar;
  const char* wgtbtagvar;
  
  if (isMC)   {
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
  TTreeReaderValue<unsigned char>   hph165 (reader, "hltphoton165");
  TTreeReaderValue<unsigned char>   hph175 (reader, "hltphoton175");
  TTreeReaderValue<unsigned char>   hsmu   (reader, "hltsinglemu");
  TTreeReaderValue<unsigned char>   hsele  (reader, "hltsingleel");
  TTreeReaderValue<unsigned int>    njets  (reader, "njets");
  TTreeReaderValue<unsigned int>    nvtx   (reader, "nvtx");
  TTreeReaderValue<unsigned int>    nbjets (reader, "nbjetslowpt");
  TTreeReaderValue<double>          wgtsum (reader, wgtsumvar);
  TTreeReaderValue<double>          wgtpu  (reader, wgtpuvar);
  TTreeReaderValue<double>          wgtbtag (reader, wgtbtagvar);
  TTreeReaderValue<double>          wgt    (reader, "wgt");
  TTreeReaderValue<double>          xsec   (reader, "xsec");

  TString metVar = Form("%s",binning.met_.c_str());
  if(binning.category_ == "zmm")
    metVar.ReplaceAll("pf","mu");
  else if(binning.category_ == "zee")
    metVar.ReplaceAll("pf","el");

  TTreeReaderValue<double>          met    (reader, metVar);
  TTreeReaderValue<double>          metphi (reader, metVar+"phi");
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

  TTreeReaderValue<double>          zpt   (reader, "zpt");
  TTreeReaderValue<double>          zeta  (reader, "zeta");
  TTreeReaderValue<double>          zphi  (reader, "zphi");
  TTreeReaderValue<double>          zmass (reader, "zmass");
  TTreeReaderArray<double>          zeept (reader,"zeept.zeeept");
  TTreeReaderValue<double>          zeeeta  (reader, "zeeeta");
  TTreeReaderValue<double>          zeephi  (reader, "zeephi");
  TTreeReaderValue<double>          zeemass (reader, "zeemass");

  double yield  = 0.0;
  double zptsum = 0.0;

  // loop on events                                                                                                                                                           
  while(reader.Next()){

    double weight = 1.0;
    double kfact  = 1.0;
    double puwgt  = 0.0;
    double effsf  = 1.0;
    double trgsf  = 1.0;

    if (*nvtx <= 40) puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));

    if (!isMC && *fcsc    == 0) continue;
    if (!isMC && *fhbhe   == 0) continue;
    if (!isMC && *fhbhei  == 0) continue;
    if (!isMC && *feesc   == 0) continue;
    if (!isMC && *fecal   == 0) continue;
    if (!isMC && *fnvtx   == 0) continue;

    double obs = 0;
    if(binning.observable_ == "nvtx")
      obs = *nvtx;
    else if(binning.observable_ == "zpt"  and binning.category_ == "zmm")
      obs = *zpt;
    else if(binning.observable_ == "zpt"  and binning.category_ == "zee")
      obs = zeept[0];
    else if(binning.observable_ == "zeta" and binning.category_ == "zmm")
      obs = fabs(*zeta);
    else if(binning.observable_ == "zeta" and binning.category_ == "zee")
      obs = fabs(*zeeeta);

    if(*nbjets > 1) continue;

    if(binning.category_ == "zmm"){
      unsigned char hlt = (*hmnm90) + (*hmnm120) + (*hmwm90) + (*hmwm120) + (*hmwm170) + (*hmwm300) +(*hsmu);
      if(not isMC and hlt == 0) continue;
      if(obs < binning.bins_.at(binning.iBin_)) continue;
      if(obs >= binning.bins_.at(binning.iBin_+1)) continue;
      bool istight = false;
      if (*mu1id == 1 && *mu1pt > 20.) istight = true;
      if (*mu2id == 1 && *mu2pt > 20.) istight = true;
      if (!istight) continue;
      if (*mu1id == 1) effsf *= msfthist->GetBinContent(msfthist->FindBin(min(999., *mu1pt), fabs(*mu1eta)));
      else             effsf *= msflhist->GetBinContent(msflhist->FindBin(min(999., *mu1pt), fabs(*mu1eta)));
      if (*mu2id == 1) effsf *= msfthist->GetBinContent(msfthist->FindBin(min(999., *mu2pt), fabs(*mu2eta)));
      else             effsf *= msflhist->GetBinContent(msflhist->FindBin(min(999., *mu2pt), fabs(*mu2eta)));
    }

    if(binning.category_ == "zee"){
      unsigned char hlt = (*hsele) + (*hph165) + (*hph175);
      if (not isMC and hlt == 0) continue;
      bool istight = false;
      if(obs < binning.bins_.at(binning.iBin_)) continue;
      if(obs >= binning.bins_.at(binning.iBin_+1)) continue;
      if (*el1id == 1 && *el1pt > 40.) istight = true;
      if (*el2id == 1 && *el2pt > 40.) istight = true;
      if (!istight) continue;
      if (*el1id == 1) effsf *= esfthist->GetBinContent(esfthist->FindBin(min(999., *el1pt), fabs(*el1eta)));
      else             effsf *= esflhist->GetBinContent(esflhist->FindBin(min(999., *el1pt), fabs(*el1eta)));
      if (*el2id == 1) effsf *= esfthist->GetBinContent(esfthist->FindBin(min(999., *el2pt), fabs(*el2eta)));
      else             effsf *= esflhist->GetBinContent(esflhist->FindBin(min(999., *el2pt), fabs(*el2eta)));
    }

    // make x and y projections                                                                                                                                               
    double metx = 0.0;
    double mety = 0.0;
    double uphi = 0.0;
    double mphi = 0.0;
    double u1 = 0.0;
    double u2 = 0.0;
    metx += *met * cos(*metphi);
    mety += *met * sin(*metphi);
    uphi = atan2(-sin(*zphi)  , -cos(*zphi));
    mphi = atan2(-sin(*metphi), -cos(*metphi));
    u1 = *met * cos(dphi(mphi, uphi));
    u2 = *met * sin(dphi(mphi, uphi));
    double fillvar = u1;

    if (fillvar >= xmax) continue;
    if (fillvar <  xmin) continue;
    // being at response pleteau
    if(binning.observable_ != "zpt"){
      if(binning.category_ == "zmm" and *zpt <= 50) continue;
      if(binning.category_ == "zee" and zeept[0] <= 50) continue;
    }

    double evtwgt = 1.0;
    if (isMC) weight = (*xsec)*(lumi)*(*wgt)*(kfact)*(puwgt)*(trgsf)*(effsf)*(*wgtbtag)/(*wgtsum);
    hist->Fill(fillvar, weight);
    // count the number of events accouring to the weight and the total zpt sum in the bin                                                                                     
    yield  += evtwgt;
    if(binning.category_ == "zmm")
      zptsum += *zpt*evtwgt;
    else if(binning.category_ == "zee")
      zptsum += zeept[0]*evtwgt;
  }

  cout << hist->GetName() << " integral : " << yield << ", <Obs> = " << zptsum/yield << " hist GetMean "<<hist->GetMean()<<endl;
  return zptsum/yield;
}


void getresponse(TTree* tree,
		 const char* histname, bool isMC,
		 const analysisBin & binning,
		 double & val, double & err,
		 const float  & lumi,
		 const string & outputDIR
		 ) {


  // define response histogram
  int nbins   = 135;
  double xmin = -ZPT_bins.back()*1.25;
  double xmax = ZPT_bins.back()*1.25;
  TH1F  hist(histname, "", nbins, xmin, xmax);

  double zptavg = drawplot(tree, &hist, isMC, xmin, xmax, binning, lumi);

  double fitrangemin = hist.GetMean() - 3.*hist.GetStdDev();
  double fitrangemax = hist.GetMean() + 3.*hist.GetStdDev();

  TF1 fitfunc("fitfunc", "gaus(0)", xmin, xmax);
  hist.Fit(&fitfunc, "", "", fitrangemin, fitrangemax);
  // Take the mean value of u parallel and divide by zptavg (mean zpt in the bin considered)
  val = fitfunc.GetParameter(1)/zptavg;
  err = fitfunc.GetParError(1)/zptavg;

  TH1F* h = (TH1F*)hist.Clone("h");
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  hist.Draw("EP");
  if(isMC)
    c->SaveAs((outputDIR+"/histMC_"+binning.category_+"_"+binning.observable_+"_"+to_string(binning.bins_.at(binning.iBin_))+"_"+to_string(binning.bins_.at(binning.iBin_+1))+"_metScale.png").c_str(),"png");
  else
    c->SaveAs((outputDIR+"/histDATA_"+binning.category_+"_"+binning.observable_+"_"+to_string(binning.bins_.at(binning.iBin_))+"_"+to_string(binning.bins_.at(binning.iBin_+1))+"_metScale.png").c_str(),"png");

}

// main function for met response plots
void makeMETResponse(string baseDIR,   // directory with ntuples                                                                                                               
		     string category,  // "zmm" or "zee" at the moment                                                                                                         
		     string met,       // can be "t1pfmet", "pfmet", "puppit1pfmet", "puppipfmet"                                                                              
		     string observable, // can be "zpt","zeta" ot "nvtx"                                                                                                       
		     string outputDIR,
		     float  lumi = 0.218){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputDIR).c_str());

  if(met != "t1pfmet" and observable!= "pfmet" and met != "puppit1pfmet" and met != "puppipfmet"){
    cerr<<"Not a good observable --> return"<<endl;
    return;
  }

  if(observable != "zpt" and observable != "zeta" and observable != "nvtx"){
    cerr<<"Not a good observable --> return"<<endl;
    return;
  }

  // select the binning                                                                                                                                                        
  vector<double> bins;
  if(observable == "zpt")
    bins = ZPT_bins;
  else if(observable == "zeta")
    bins = ZEta_bins;
  else if(observable == "nvtx")
    bins = Nvtx_bins;
  analysisBin binning (met,observable,category,bins);

  TChain* zlltree = new TChain("tree/tree");
  TChain* dattree = new TChain("tree/tree");

  if( category != "zmm" and category != "zee"){
    cerr<<"Problem with category --> atm only zmm and zee can be used --> exit"<<endl;
    return;
  }

  zlltree->Add((baseDIR+"/DYJets/"+category+"filter/*root").c_str());
  if(category == "zmm")
    dattree->Add((baseDIR+"/SingleMuon/"+category+"filter/*root").c_str());
  else if(category == "zee")
    dattree->Add((baseDIR+"/SingleElectron/"+category+"filter/*root").c_str());

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
  TPad *pad1      = new TPad("pad1","pad1",0,0.3,1,1.0);
  canvas->cd();
  TPad *pad2      = new TPad("pad2","pad2",0,0.1,1,0.3);
    
  double val = 0.0;
  double err = 0.0;

  std::vector<double> zllvalvector;
  std::vector<double> zllerrvector;
  std::vector<double> datvalvector;
  std::vector<double> daterrvector;

  /// loop on the met bins                                                                                                                                                  
  for (int i = 0; i < binning.bins_.size()-1; i++) {
    binning.iBin_ = i;
    getresponse(zlltree, "zllhist", true, binning, val, err, lumi, outputDIR);
    zllvalvector.push_back(val);
    zllerrvector.push_back(err);
    getresponse(dattree, "dathist", false, binning, val, err, lumi, outputDIR);
    datvalvector.push_back(val);
    daterrvector.push_back(err);
  }

  TH1F* zllres = new TH1F("zllres", "", bins.size()-1, &bins[0]);
  TH1F* datres = new TH1F("datres", "", bins.size()-1, &bins[0]);

  for (int i = 1; i <= bins.size()-1; i++) {
    zllres->SetBinContent(i, zllvalvector[i-1]);
    zllres->SetBinError  (i, zllerrvector[i-1]);
    datres->SetBinContent(i, datvalvector[i-1]);
    datres->SetBinError  (i, daterrvector[i-1]);
  }
 

  zllres->SetLineColor(kRed); 
  zllres->SetMarkerColor(kRed); 
  zllres->SetMarkerSize(1); 
  zllres->SetMarkerStyle(20);
  datres->SetMarkerColor(kBlack);
  datres->SetMarkerSize(1);
  datres->SetLineColor(kBlack);
  datres->SetMarkerStyle(20);

  canvas->cd();
  TH1* frame = canvas->DrawFrame(bins.front(), 0.4, bins.back(), 1.1, "");
  if(observable == "zpt")
    frame->GetXaxis()->SetTitle("Z p_{T} [GeV]");
  else if(observable == "zeta")
    frame->GetXaxis()->SetTitle("Z |#eta|");
  else if(observable == "nvtx")
    frame->GetXaxis()->SetTitle("N_{PV}");
  frame->GetYaxis()->SetTitle("E_{T}^{miss} Response");
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetLabelSize(0.9*frame->GetYaxis()->GetLabelSize());
  frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());

  pad1->SetTopMargin(0.06);
  pad1->SetRightMargin(0.075);
  pad1->SetBottomMargin(0.035);
  pad1->Draw();
  pad1->cd();
  frame ->Draw();
  zllres->Draw("PE SAME");
  datres->Draw("PE SAME");

  TLegend* leg = new TLegend(0.6, 0.3, 0.9, 0.6);
  leg->SetFillColor(0);
  leg->AddEntry(datres, "Data");
  if(category == "zmm")
    leg->AddEntry(zllres, "Z(#mu#mu) MC");
  else if(category == "zee")
    leg->AddEntry(zllres, "Z(ee) MC");

  leg->Draw("SAME");

  pad1->RedrawAxis();  
  TString lumi_ = Form("%.2f",lumi);
  CMS_lumi(pad1,string(lumi_),true);

  canvas->cd();
  pad2->SetTopMargin(0.08);
  pad2->SetRightMargin(0.075);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();

  TH1* dahist = (TH1*)datres->Clone("dahist");
  TH1* mchist = (TH1*)zllres->Clone("mchist");
  dahist->Divide(mchist);

  dahist->GetXaxis()->SetLabelSize(3.0*dahist->GetXaxis()->GetLabelSize());
  dahist->GetYaxis()->SetLabelSize(3.0*dahist->GetYaxis()->GetLabelSize());
  dahist->GetYaxis()->SetRangeUser(0.9, 1.1);
  dahist->GetYaxis()->SetNdivisions(504);
  dahist->SetMarkerSize(0);
  dahist->SetMarkerStyle(20);
  dahist->SetMarkerSize(1.0);
  dahist->Draw("PE");

  pad1->cd();
  pad1->RedrawAxis();

  canvas->SaveAs((outputDIR+"/metResponse_"+category+"_"+observable+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/metResponse_"+category+"_"+observable+".png").c_str(),"png");
  
}
