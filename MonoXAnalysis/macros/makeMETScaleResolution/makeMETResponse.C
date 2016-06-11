#include "../CMS_lumi.h"

// met binning used in this study
vector<double> ZPT_bins      = {0.0, 5., 10., 15., 20., 30., 40., 60., 80., 100., 125, 155.,190.,250.};
vector<double> GammaPT_bins  = {175.,200,250,300,500};
vector<double> Nvtx_bins     = {0.0, 5.5, 10.5, 15.5, 20.5, 25.5, 30.5, 50.};
vector<double> ZEta_bins     = {0.0, 0.5, 1.0,  1.5, 2.0, 2.5};
vector<double> GammaEta_bins = {0.0, 0.5, 1.0,  1.5};

TH1F* bosonPtWeight = NULL;
bool  metNoHF_        ;
bool  isVoigModel_    ;

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
    if(observable_ != "zpt" and observable_ != "nvtx" and observable_ != "zeta" and observable != "phpt" and observable != "pheta")
      cout<<"Un-recognized observable option --> please check"<<endl;
    if(category_ != "zmm" and category_ !=  "zee" and category_ != "gam")
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


void getBosonPt(TTree* tree, TH1* hist, const bool & isMC, const float & lumi, const string & category){

  TFile* pufile = new TFile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt.root");
  TH1*   puhist = (TH1*)pufile->Get("puhist");

  TTreeReader reader(tree);
  const char* wgtsumvar = "";
  const char* wgtpuvar = "";
  const char* wgtbtagvar = "";

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
  TTreeReaderValue<double>          wgtsum (reader, wgtsumvar);
  TTreeReaderValue<double>          wgtpu  (reader, wgtpuvar);
  TTreeReaderValue<double>          wgtbtag (reader, wgtbtagvar);
  TTreeReaderValue<double>          wgt    (reader, "wgt");
  TTreeReaderValue<double>          xsec   (reader, "xsec");
  TTreeReaderValue<double>          zpt   (reader, "zpt");
  TTreeReaderValue<double>          zeept (reader, "zeept");  
  TTreeReaderValue<double>          phpt  (reader, "phpt");  
  TTreeReaderValue<unsigned int>    nvtx   (reader, "nvtx");

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


    double fillvar = 0;
    if(category == "zmm")
      fillvar = *zpt;
    else if(category == "zee")
      fillvar = *zeept;
    else if(category == "gam")
      fillvar = *phpt;

    double evtwgt = 1.0;
    if (isMC) weight = (*xsec)*(lumi)*(*wgt)*(kfact)*(puwgt)*(trgsf)*(effsf)*(*wgtbtag)/(*wgtsum);
    hist->Fill(fillvar, weight);
  }
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

  // electron and muon ID scale factor files                                                                                                                                    
  TFile sffile_eleTight("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/scaleFactor_electron_tightid.root");
  TFile sffile_eleVeto("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/scaleFactor_electron_vetoid.root");
  TFile sffile_muTight("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/scaleFactor_muon_tightid.root");
  TFile sffile_muLoose("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF_2016/scaleFactor_muon_looseid.root");

  TH2*  msfloose = (TH2*)sffile_muLoose.Get("scaleFactor_muon_looseid_RooCMSShape");
  TH2*  msftight = (TH2*)sffile_muTight.Get("scaleFactor_muon_tightid_RooCMSShape");
  TH2*  esfveto  = (TH2*)sffile_eleVeto.Get("scaleFactor_electron_vetoid_RooCMSShape");
  TH2*  esftight = (TH2*)sffile_eleTight.Get("scaleFactor_electron_tightid_RooCMSShape");

  // Photon ID scale factor                                                                                                                                                     
  TFile sffile_phoLoose("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF_2016/scaleFactor_photon_looseid.root");
  TFile sffile_phoMedium("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF_2016/scaleFactor_photon_mediumid.root");
  TH2*  psfloose  = (TH2*)sffile_phoLoose.Get("scaleFactor_photon_looseid_RooCMSShape");
  TH2*  psfmedium = (TH2*)sffile_phoMedium.Get("scaleFactor_photon_mediumid_RooCMSShape");

  // Photon Purity                                                                                                                                                              
  TFile purityfile_photon ("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/PhotonSFandEffandPurity_Lumi2p1fb_0202.root");
  TH2*  purhist = (TH2*) purityfile_photon.Get("PhotonPurity");

  // trigger files used for 2016                                                                                                                                                
  TFile triggerfile_SinglEle("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/triggerEfficiency_DATA_SingleElectron.root");
  TFile triggerfile_SingleMu("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/triggerEfficiency_DATA_SingleMuon.root");
  TFile triggerfile_MET("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/metTriggerEfficiency.root");
  TFile triggerfile_SinglePhoton("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/photonTriggerEfficiency.root");

  TH2*  triggerelhist = (TH2*) triggerfile_SinglEle.Get("trigeff_ele27wptight");
  TH2*  triggermuhist = (TH2*) triggerfile_SingleMu.Get("trigeff_muIso");
  TF1*  triggermet = (TF1*) triggerfile_MET.Get("efficiency_func");
  TF1*  triggerphoton = (TF1*)triggerfile_SinglePhoton.Get("efficiency_func");


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
  else if(binning.category_ == "gam")
    metVar.ReplaceAll("pf","ph");

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

  TTreeReaderValue<int>             phidm  (reader, "phidm");
  TTreeReaderValue<int>             phidt  (reader, "phidt");
  TTreeReaderValue<double>          phpt   (reader, "phpt");
  TTreeReaderValue<double>          pheta  (reader, "pheta");
  TTreeReaderValue<double>          phphi  (reader, "phphi");

  TTreeReaderValue<double>          zpt   (reader, "zpt");
  TTreeReaderValue<double>          zeta  (reader, "zeta");
  TTreeReaderValue<double>          zphi  (reader, "zphi");
  TTreeReaderValue<double>          zmass (reader, "zmass");
  TTreeReaderValue<double>          zeept (reader,"zeept");
  TTreeReaderValue<double>          zeeeta  (reader, "zeeeta");
  TTreeReaderValue<double>          zeephi  (reader, "zeephi");
  TTreeReaderValue<double>          zeemass (reader, "zeemass");

  TTreeReaderValue<double>          pfmetchargedhadron(reader,"pfmetchargedhadron");
  TTreeReaderValue<double>          pfmetchargedhadronphi(reader,"pfmetchargedhadronphi");
  TTreeReaderValue<double>          pfmetneutralhadron(reader,"pfmetneutralhadron");
  TTreeReaderValue<double>          pfmetneutralhadronphi(reader,"pfmetneutralhadronphi");
  TTreeReaderValue<double>          pfmetelectrons(reader,"pfmetelectrons");
  TTreeReaderValue<double>          pfmetelectronsphi(reader,"pfmetelectronsphi");
  TTreeReaderValue<double>          pfmetmuons(reader,"pfmetmuons");
  TTreeReaderValue<double>          pfmetmuonsphi(reader,"pfmetmuonsphi");
  TTreeReaderValue<double>          pfmetphotons(reader,"pfmetphotons");
  TTreeReaderValue<double>          pfmetphotonsphi(reader,"pfmetphotonsphi");

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
      obs = *zeept;
    else if(binning.observable_ == "phpt" and binning.category_ == "gam")
      obs = *phpt;
    else if(binning.observable_ == "zeta" and binning.category_ == "zmm")
      obs = fabs(*zeta);
    else if(binning.observable_ == "zeta" and binning.category_ == "zee")
      obs = fabs(*zeeeta);
    else if(binning.observable_ == "pheta" and binning.category_ == "gam")
      obs = fabs(*pheta);

    if(*nbjets > 1) continue;

    if (binning.category_ == "gam") {
      unsigned char hlt = (*hph165) + (*hph175);
      if (not isMC and hlt == 0) continue;
      if (*phpt < 175. || fabs(*pheta) > 1.4442) continue;
      trgsf *= triggerphoton->Eval(min(*phpt,triggerphoton->GetXaxis()->GetXmax()));
      effsf *= psfmedium->GetBinContent(psfmedium->FindBin(min(*phpt,psfmedium->GetXaxis()->GetBinLowEdge(psfmedium->GetNbinsX()+1)-1),*pheta));
      if(obs < binning.bins_.at(binning.iBin_)) continue;
      if(obs >= binning.bins_.at(binning.iBin_+1)) continue;
    }


    if(binning.category_ == "zmm"){
      unsigned char hlt = (*hmnm90) + (*hmnm120) + (*hmwm90) + (*hmwm120) + (*hmwm170) + (*hmwm300) +(*hsmu);
      if(not isMC and hlt == 0) continue;
      if(obs < binning.bins_.at(binning.iBin_)) continue;
      if(obs >= binning.bins_.at(binning.iBin_+1)) continue;
      bool istight = false;
      if (*mu1id == 1 && *mu1pt > 20.) istight = true;
      if (*mu2id == 1 && *mu2pt > 20.) istight = true;
      if (!istight) continue;
      if (*mu1id == 1 ) effsf *= msftight->GetBinContent(msftight->FindBin(min(*mu1pt,msftight->GetXaxis()->GetBinLowEdge(msftight->GetNbinsX()+1)-1),*mu1eta));
      else              effsf *= msfloose->GetBinContent(msfloose->FindBin(min(*mu1pt,msfloose->GetXaxis()->GetBinLowEdge(msfloose->GetNbinsX()+1)-1),*mu1eta));
      if (*mu2id == 1 ) effsf *= msftight->GetBinContent(msftight->FindBin(min(*mu2pt,msftight->GetXaxis()->GetBinLowEdge(msftight->GetNbinsX()+1)-1),*mu2eta));
      else              effsf *= msfloose->GetBinContent(msfloose->FindBin(min(*mu2pt,msfloose->GetXaxis()->GetBinLowEdge(msfloose->GetNbinsX()+1)-1),*mu2eta));

      if (*mu1id == 1 and *mu2id == 1)  // both tight efficiency is eff1+eff2-eff1*eff2 that at plateau is ~1                                                                    
        trgsf *= 1;
      else if(*mu1id == 1 and *mu2id != 1)
        trgsf *= triggermuhist->GetBinContent(triggermuhist->FindBin(min(*mu1pt,triggermuhist->GetXaxis()->GetBinLowEdge(triggermuhist->GetNbinsX()+1)-1),*mu1eta));
      else if(*mu1id != 1 and *mu2id == 1)
        trgsf *= triggermuhist->GetBinContent(triggermuhist->FindBin(min(*mu2pt,triggermuhist->GetXaxis()->GetBinLowEdge(triggermuhist->GetNbinsX()+1)-1),*mu2eta));
    }

    if(binning.category_ == "zee"){
      unsigned char hlt = (*hsele) + (*hph165) + (*hph175);
      if (not isMC and hlt == 0) continue;
      if(obs < binning.bins_.at(binning.iBin_)) continue;
      if(obs >= binning.bins_.at(binning.iBin_+1)) continue;
      bool istight = false;
      if (*el1id == 1 && *el1pt > 40.) istight = true;
      if (*el2id == 1 && *el2pt > 40.) istight = true;
      if (!istight) continue;
      if (*el1id == 1 and *el2id == 1)  // both tight efficiency is eff1+eff2-eff1*eff2 that at plateau is ~1                                                                    
        trgsf *= 1;
      else if(*el1id == 1 and *el2id != 1)
        trgsf *= triggerelhist->GetBinContent(triggerelhist->FindBin(min(*el1pt,triggerelhist->GetXaxis()->GetBinLowEdge(triggerelhist->GetNbinsX()+1)-1),*el1eta));
      else if(*el1id != 1 and *el2id == 1)
        trgsf *= triggerelhist->GetBinContent(triggerelhist->FindBin(min(*el2pt,triggerelhist->GetXaxis()->GetBinLowEdge(triggerelhist->GetNbinsX()+1)-1),*el2eta));
    }

    // make x and y projections                                                                                                                                               
    double metx = 0.0;
    double mety = 0.0;
    double uphi = 0.0;
    double mphi = 0.0;
    double u1 = 0.0;
    double u2 = 0.0;

    if(binning.category_ == "zmm")
      uphi = atan2(-sin(*zphi)  , -cos(*zphi));
    else if(binning.category_ == "zee")
      uphi = atan2(-sin(*zeephi)  , -cos(*zeephi));
    else if(binning.category_ == "gam")
      uphi = atan2(-sin(*phphi)  , -cos(*phphi));

    if(not metNoHF_){
      mphi = atan2(-sin(*metphi), -cos(*metphi));
      u1 = *met * cos(dphi(mphi, uphi));
      u2 = *met * sin(dphi(mphi, uphi));
    }
    else{
      if(binning.category_ == "zmm"){
        metx = *pfmetchargedhadron*cos(*pfmetchargedhadronphi)+*pfmetneutralhadron*cos(*pfmetneutralhadronphi)+*pfmetelectrons*cos(*pfmetelectronsphi)+*pfmetphotons*cos(*pfmetphotonsphi);
        mety = *pfmetchargedhadron*sin(*pfmetchargedhadronphi)+*pfmetneutralhadron*sin(*pfmetneutralhadronphi)+*pfmetelectrons*sin(*pfmetelectronsphi)+*pfmetphotons*sin(*pfmetphotonsphi);

        u1 = -sqrt(metx*metx+mety*mety)*cos(dphi(uphi,atan2(mety,metx)));
        u2 = -sqrt(metx*metx+mety*mety)*sin(dphi(uphi,atan2(mety,metx)));
      }
      else if(binning.category_ == "zee"){
        metx = *pfmetchargedhadron*cos(*pfmetchargedhadronphi)+*pfmetneutralhadron*cos(*pfmetneutralhadronphi)+*pfmetmuons*cos(*pfmetmuonsphi)+*pfmetphotons*cos(*pfmetphotonsphi);
        mety = *pfmetchargedhadron*sin(*pfmetchargedhadronphi)+*pfmetneutralhadron*sin(*pfmetneutralhadronphi)+*pfmetmuons*sin(*pfmetmuonsphi)+*pfmetphotons*sin(*pfmetphotonsphi);

        u1 = -sqrt(metx*metx+mety*mety)*cos(dphi(uphi,atan2(mety,metx)));
        u2 = -sqrt(metx*metx+mety*mety)*sin(dphi(uphi,atan2(mety,metx)));
      }
      else if(binning.category_ == "gam"){
        metx = *pfmetchargedhadron*cos(*pfmetchargedhadronphi)+*pfmetneutralhadron*cos(*pfmetneutralhadronphi)+*pfmetmuons*cos(*pfmetmuonsphi)+*pfmetelectrons*cos(*pfmetelectronsphi);
        mety = *pfmetchargedhadron*sin(*pfmetchargedhadronphi)+*pfmetneutralhadron*sin(*pfmetneutralhadronphi)+*pfmetmuons*sin(*pfmetmuonsphi)+*pfmetelectrons*sin(*pfmetelectronsphi);

        u1 = -sqrt(metx*metx+mety*mety)*cos(dphi(uphi,atan2(mety,metx)));
        u2 = -sqrt(metx*metx+mety*mety)*sin(dphi(uphi,atan2(mety,metx)));
      }
    }

    double fillvar = u1;
    if (fillvar >= xmax) continue;
    if (fillvar <  xmin) continue;
    // being at response pleteau
    if(binning.observable_ != "zpt" and binning.category_ != "gam"){
      if(binning.category_ == "zmm" and *zpt   <= 50) continue;
      if(binning.category_ == "zee" and *zeept <= 50) continue;
    }

    double evtwgt = 1.0;    
    double zptwgt = 1.0;
    if(binning.category_      == "zmm" and isMC and bosonPtWeight != NULL) zptwgt = bosonPtWeight->GetBinContent(bosonPtWeight->FindBin(*zpt));
    else if(binning.category_ == "zee" and isMC and bosonPtWeight != NULL) zptwgt = bosonPtWeight->GetBinContent(bosonPtWeight->FindBin(*zeept));
    else if(binning.category_ == "gam" and isMC and bosonPtWeight != NULL) zptwgt = bosonPtWeight->GetBinContent(bosonPtWeight->FindBin(*phpt));

    if (isMC) weight = (*xsec)*(lumi)*(*wgt)*(kfact)*(puwgt)*(trgsf)*(effsf)*(*wgtbtag)*zptwgt/(*wgtsum);
    hist->Fill(fillvar, weight);
 
    // count the number of events accouring to the weight and the total zpt sum in the bin                                                                                     
    yield  += evtwgt;
    if(binning.category_ == "zmm")
      zptsum += *zpt*evtwgt;
    else if(binning.category_ == "zee")
      zptsum += *zeept*evtwgt;
    else if(binning.category_ == "gam")
      zptsum += *phpt*evtwgt;
  }

  cout << hist->GetName() << " integral : " << hist->Integral() << ", <Obs> = " << zptsum/yield << " hist GetMean "<<hist->GetMean()<<endl;
  return zptsum/yield;
}


TH1F* getresponse(TTree* tree,
		 const string & histname, bool isMC,
		 const analysisBin & binning,
		 double & val, double & err,
		 const float  & lumi,
		 const string & outputDIR
		 ) {


  // define response histogram
  int nbins   = 150;  
  double xmin = -ZPT_bins.back()*1.25; // integrated over z-pt in case  the analysis is done vs NPV or Zeta                                                                    
  double xmax = ZPT_bins.back()*1.25;

  if(binning.category_ == "gam"){
    xmin = -GammaPT_bins.back()*1.25;
    xmax =  GammaPT_bins.back()*1.25;
  }

  if(binning.bins_.at(binning.iBin_) >= 100){
    xmin = ZPT_bins.at(binning.iBin_)*0.3;
    xmax = ZPT_bins.at(binning.iBin_+1)*2;
    if(binning.category_ == "gam"){
      xmin = GammaPT_bins.at(binning.iBin_)*0.3;
      xmax = GammaPT_bins.at(binning.iBin_+1)*2;      
    }
  }
  else{
    xmin = -50;
    xmax = 150;
  }

  TH1F*  hist =  new TH1F(histname.c_str(), "", nbins, xmin, xmax);
  double zptavg = drawplot(tree, hist, isMC, xmin, xmax, binning, lumi);

  double fitrangemin = 0;
  double fitrangemax = 0;

  TF1* fitfunc = 0;
  if(isVoigModel_){
   fitrangemin = hist->GetMean() - 5.*hist->GetStdDev();
   fitrangemax = hist->GetMean() + 5.*hist->GetStdDev();
   fitfunc = new TF1("fitfunc","[0]*TMath::Voigt(x-[1],[2],[3],4)",xmin,xmax);
   fitfunc->SetParameter(0,hist->Integral());
   fitfunc->SetParameter(1,hist->GetMean());
   fitfunc->SetParameter(3,hist->GetRMS());
   fitfunc->SetParameter(2,hist->GetRMS());
  }
  else{
    fitrangemin = hist->GetMean() - 3.*hist->GetStdDev();
    fitrangemax = hist->GetMean() + 3.*hist->GetStdDev();
    fitfunc = new TF1("fitfunc", "gaus(0)", xmin, xmax);
  }
 if(isMC){
   fitfunc->SetLineColor(kBlue);
   fitfunc->SetLineWidth(2);
 }
 else{
   fitfunc->SetLineColor(kRed);
   fitfunc->SetLineWidth(2);
 }
 
 hist->Fit(fitfunc, "QME", "", fitrangemin, fitrangemax);
 // take the resolution from fitting u-perpendicular as default                                                                                                               
 val = fitfunc->GetParameter(1)/zptavg;
 err = fitfunc->GetParError(1)/zptavg;

 return hist;

}

// main function for met response plots
void makeMETResponse(string baseDIR,   // directory with ntuples                                                                                                               
		     string category,  // "zmm" or "zee" at the moment                                                                                                         
		     string met,       // can be "t1pfmet", "pfmet", "puppit1pfmet", "puppipfmet"                                                                              
		     string observable, // can be "zpt","zeta" ot "nvtx"                                                                                                       
		     string outputDIR,
		     float  lumi = 0.590,
		     bool   doBosonPtRewight = true,
		     bool   isVoigModel = false,
		     bool   metNoHF = false){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  isVoigModel_ = isVoigModel;
  metNoHF_ = metNoHF;
  
  system(("mkdir -p "+outputDIR).c_str());

  if(met != "t1pfmet" and met != "pfmet" and met != "puppit1pfmet" and met != "puppipfmet"){
    cerr<<"Not a good observable --> return"<<endl;
    return;
  }

  if(observable != "zpt" and observable != "zeta" and observable != "nvtx" and observable != "phpt" and observable != "pheta"){
    cerr<<"Not a good observable --> return"<<endl;
    return;
  }

  // select the binning                                                                                                                                                        
  vector<double> bins;
  if(observable == "zpt")
    bins = ZPT_bins;
  else if(observable == "phpt")
    bins = GammaPT_bins;
  else if(observable == "zeta")
    bins = ZEta_bins;
  else if(observable == "pheta")
    bins = GammaEta_bins;
  else if(observable == "nvtx")
    bins = Nvtx_bins;
  analysisBin binning (met,observable,category,bins);

  TChain* zlltree = new TChain("tree/tree");
  TChain* dattree = new TChain("tree/tree");

  if( category != "zmm" and category != "zee" and category != "gam"){
    cerr<<"Problem with category --> atm only zmm and zee and gam can be used --> exit"<<endl;
    return;
  }

  if(category != "gam")
    zlltree->Add((baseDIR+"/DYJets/"+category+"filter/*root").c_str());
  else
    zlltree->Add((baseDIR+"/PhotonJets/"+category+"filter/*root").c_str());

  if(category == "zmm")
    dattree->Add((baseDIR+"/SingleMuon/"+category+"filter/*root").c_str());
  else if(category == "zee")
    dattree->Add((baseDIR+"/SingleElectron/"+category+"filter/*root").c_str());
  else if(category == "gam")
    dattree->Add((baseDIR+"/SinglePhoton/"+category+"filter/*root").c_str());

  // add also minor backgrounds to the zlltree
  zlltree->Add((baseDIR+"/Top/"+category+"filter/*root").c_str());
  zlltree->Add((baseDIR+"/WJets/"+category+"filter/*root").c_str());
  zlltree->Add((baseDIR+"/DiBoson/"+category+"filter/*root").c_str());

  // boson pt re-weight
  TH1F* bosonPt_data = NULL;
  TH1F* bosonPt_MC = NULL;
  if(category != "gam"){
    bosonPt_data = new TH1F("bosonPt_data","",70,0,ZPT_bins.back());
    bosonPt_MC   = new TH1F("bosonPt_MC",""  ,70,0,ZPT_bins.back());
  }
  else{
    bosonPt_data = new TH1F("bosonPt_data","",50,GammaPT_bins.front(),GammaPT_bins.back());
    bosonPt_MC   = new TH1F("bosonPt_MC",""  ,50,GammaPT_bins.front(),GammaPT_bins.back());
  }
  bosonPt_data->Sumw2();
  bosonPt_MC->Sumw2();

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);

  if(doBosonPtRewight){
    getBosonPt(dattree,bosonPt_data,false,lumi,category);
    getBosonPt(zlltree,bosonPt_MC,true,lumi,category);
    
    bosonPt_data->Scale(1./bosonPt_data->Integral());
    bosonPt_MC->Scale(1./bosonPt_MC->Integral());
    bosonPtWeight = (TH1F*) bosonPt_data->Clone("bosonPtWeight");
    bosonPtWeight->Divide(bosonPt_MC);
    canvas->cd();
    bosonPtWeight->Draw("EP");
    canvas->SaveAs((outputDIR+"/bosonPtWeight.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/bosonPtWeight.pdf").c_str(),"pdf");
  }


  TPad *pad1      = new TPad("pad1","pad1",0,0.3,1,1.0);
  TPad *pad2      = new TPad("pad2","pad2",0,0.1,1,0.3);
    
  double val = 0.0;
  double err = 0.0;

  std::vector<double> zllvalvector;
  std::vector<double> zllerrvector;
  std::vector<double> datvalvector;
  std::vector<double> daterrvector;

  TH1F* zllhist = NULL;
  TH1F* dathist = NULL;
  /// loop on the met bins                                                                                                                                                  
  for (int i = 0; i < binning.bins_.size()-1; i++) {
    // plotting
    TCanvas* can = new TCanvas(("can_"+to_string(binning.bins_.at(i))+"_"+to_string(binning.bins_.at(i+1))).c_str(),"",600,650);
    can->cd();
    binning.iBin_ = i;
    zllhist = getresponse(zlltree, "zllhist_"+to_string(binning.bins_.at(i))+"_"+to_string(binning.bins_.at(i+1)), true, binning, val, err, lumi, outputDIR);
    zllvalvector.push_back(val);
    zllerrvector.push_back(err);
    dathist = getresponse(dattree, "dathist_"+to_string(binning.bins_.at(i))+"_"+to_string(binning.bins_.at(i+1)), false, binning, val, err, lumi, outputDIR);
    datvalvector.push_back(val);
    daterrvector.push_back(err);

    zllhist->SetFillColor(kBlue);
    zllhist->SetFillStyle(3001);
    zllhist->SetLineColor(kBlack);
    zllhist->GetXaxis()->SetTitle("u_{#parallel} (GeV)");
    zllhist->GetYaxis()->SetTitle("Events");
    zllhist->GetYaxis()->SetRangeUser(0,max(zllhist->GetMaximum(),dathist->GetMaximum())*1.2);
    zllhist->Draw("hist"); 
    zllhist->GetListOfFunctions()->At(0)->Draw("Lsame");
    double chi2_mc = ((TF1*) zllhist->GetListOfFunctions()->At(0))->GetChisquare()/((TF1*) zllhist->GetListOfFunctions()->At(0))->GetNDF();
    TString lumi_ = Form("%.2f",lumi);
    CMS_lumi(can,string(lumi_),true);
    
    dathist->SetMarkerColor(kBlack);
    dathist->SetMarkerStyle(20);
    dathist->SetMarkerSize(1);
    dathist->SetLineColor(1);
    dathist->Draw("EPsame");
    dathist->GetListOfFunctions()->At(0)->Draw("Lsame");
    double chi2_data = ((TF1*) dathist->GetListOfFunctions()->At(0))->GetChisquare()/((TF1*) dathist->GetListOfFunctions()->At(0))->GetNDF();

    TLegend* leg = new TLegend(0.65,0.65,0.9,0.9);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->AddEntry(dathist,"Data","EP");
    if(category != "gam")
      leg->AddEntry(zllhist,"DY MC","F");
    else
      leg->AddEntry(zllhist,"#gamma+jets MC","F");

    leg->AddEntry((TObject*)0,Form("#chi^{2} (MC) = %.2f",chi2_mc),"");
    leg->AddEntry((TObject*)0,Form("#chi^{2} (DATA) = %.2f",chi2_data),"");
    leg->Draw("same");
    can->SaveAs((outputDIR+"/dataVsMC_"+category+"_"+observable+"_"+to_string(binning.bins_.at(binning.iBin_))+"_"+to_string(binning.bins_.at(binning.iBin_+1))+".pdf").c_str(),"pdf");
    can->SaveAs((outputDIR+"/dataVsMC_"+category+"_"+observable+"_"+to_string(binning.bins_.at(binning.iBin_))+"_"+to_string(binning.bins_.at(binning.iBin_+1))+".png").c_str(),"png");
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
  else if(observable == "phpt")
    frame->GetXaxis()->SetTitle("#gamma p_{T} [GeV]");
  else if(observable == "zeta")
    frame->GetXaxis()->SetTitle("Z |#eta|");
  else if(observable == "pheta")
    frame->GetXaxis()->SetTitle("#gamma |#eta|");
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
  zllres->SetLineWidth(2);
  zllres->Draw("hist SAME");
  datres->Draw("PE SAME");

  TLegend* leg = new TLegend(0.6, 0.3, 0.9, 0.6);
  leg->SetFillColor(0);
  leg->AddEntry(datres, "Data","PE");
  if(category == "zmm")
    leg->AddEntry(zllres, "Z(#mu#mu) MC","L");
  else if(category == "zee")
    leg->AddEntry(zllres, "Z(ee) MC","L");
  else if(category == "gam")
    leg->AddEntry(zllres, "#gamma MC","L");

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
  dahist->GetYaxis()->SetTitle("Data/MC");  
  dahist->GetYaxis()->SetTitleSize(dahist->GetYaxis()->GetTitleSize()*2);
  dahist->Draw("PE");

  pad1->cd();
  pad1->RedrawAxis();

  canvas->SaveAs((outputDIR+"/metResponse_"+category+"_"+observable+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/metResponse_"+category+"_"+observable+".png").c_str(),"png");
  
}
