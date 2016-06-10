#ifndef MAKEHIST_H
#define MAKEHIST_H

enum class Sample { qcd, sig, gam, wmn, zmm, wen, zee, topmu, topel};

vector<string> split(const string &s, char delim) {
    vector<string> parts;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        parts.push_back(item);
    }
    return parts;
}

template <typename T>
struct is_vector : std::false_type { };

template <typename... T> struct is_vector<std::vector<T...> > : std::true_type { };

template<typename T>
double getValueFromVar(std::true_type , TTreeReaderValue<T>& var, size_t i) {
    return double((*var)[i]);
}

template<typename T>
double getValueFromVar(std::false_type, TTreeReaderValue<T>& var, size_t) {
    return double(*var);
}

template<typename T>
void templatedMakeHist(TTree* tree, TH1* hist, const char* varstr, bool isMC, Sample chan, double lumi, vector<TH1*> khists, size_t index, bool isInclusive) {

  TTreeReader reader(tree);

  // pileup re-weight after 864 pb-1
  TFile pufile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt.root");
  TH1*  puhist = (TH1*)pufile.Get("puhist");
  
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

  // Branches 
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
    TTreeReaderValue<unsigned char>   hsmu   (reader, "hltsinglemu");
    TTreeReaderValue<unsigned char>   hsele  (reader, "hltsingleel");
    TTreeReaderValue<unsigned char>   hph165 (reader, "hltphoton165");
    TTreeReaderValue<unsigned char>   hph175 (reader, "hltphoton175");

    TTreeReaderValue<unsigned int>    njets  (reader, "njets");
    TTreeReaderValue<unsigned int>    nvtx   (reader, "nvtx");
    TTreeReaderValue<unsigned int>    njetsin(reader, "njetsinc");
    TTreeReaderValue<unsigned int>    nbjets (reader, "nbjetslowpt");

    TTreeReaderValue<vector<double> > jetpt  (reader, "centraljetpt");
    TTreeReaderValue<vector<double> > chfrac (reader, "centraljetCHfrac");
    TTreeReaderValue<vector<double> > nhfrac (reader, "centraljetNHfrac");

    TTreeReaderValue<T>               var    (reader, varstr);
    TTreeReaderValue<double>          wgtsum (reader, wgtsumvar);
    TTreeReaderValue<double>          wgtpu  (reader, wgtpuvar);
    TTreeReaderValue<double>          wgtbtag (reader, wgtbtagvar);
    TTreeReaderValue<double>          wgt    (reader, "wgt");
    TTreeReaderValue<double>          xsec   (reader, "xsec");
    TTreeReaderValue<double>          wzpt   (reader, "wzpt");

    TTreeReaderValue<double>          t1pfmet(reader, "t1pfmet");
    TTreeReaderValue<double>          t1mumet(reader, "t1mumet");
    TTreeReaderValue<double>          t1elmet(reader, "t1elmet");
    TTreeReaderValue<double>          t1phmet(reader, "t1phmet");

    TTreeReaderValue<double>          jetmetm(reader, "incjetmumetdphimin4");
    TTreeReaderValue<double>          jetmete(reader, "incjetelmetdphimin4");
    TTreeReaderValue<double>          jetmetp(reader, "incjetphmetdphimin4");

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

    TTreeReaderValue<double>          wmt (reader, "wmt");
    TTreeReaderValue<double>          wemt(reader, "wemt");



    Double_t max = hist->GetBinLowEdge(hist->GetNbinsX()) + hist->GetBinWidth(hist->GetNbinsX());
    Double_t mid = hist->GetBinLowEdge(hist->GetNbinsX()) + hist->GetBinWidth(hist->GetNbinsX())/ 2.0;

    // event loop
    while (reader.Next()) {

        double weight = 1.0;
        double kfact  = 1.0;
        double puwgt  = 1.0;
        double effsf  = 1.0;
        double trgsf  = 1.0;

        if (!isMC && *fcsc    == 0) continue;
        if (!isMC && *fhbhe   == 0) continue;
        if (!isMC && *fhbhei  == 0) continue;
        if (!isMC && *feesc   == 0) continue;
        if (!isMC && *fecal   == 0) continue;
        if (!isMC && *fnvtx   == 0) continue;

	if (not isInclusive and *njets < 1) continue;
        if ((chan != Sample::topmu and chan != Sample::topel) and *nbjets > 0) continue;	
	if (not isInclusive and (*chfrac)[0] < 0.1 ) continue;
	if (not isInclusive and (*nhfrac)[0] > 0.8 ) continue;	
	if (not isInclusive and (*jetpt )[0] < 100.) continue;
        double jetmetdphi = *jetmetm;
        if (chan == Sample::qcd || chan == Sample::gam) jetmetdphi = *jetmetp; 
        if (chan == Sample::wen || chan == Sample::zee || chan == Sample::topel) jetmetdphi = *jetmete; 
        if (not isInclusive and jetmetdphi < 0.5) continue;

        double recoil = *t1mumet;
        if (chan == Sample::qcd || chan == Sample::gam) recoil = *t1phmet;
        if (chan == Sample::wen || chan == Sample::zee || chan == Sample::topel) recoil = *t1elmet; 
	if (not isInclusive and recoil < 200.) continue;

        if (chan == Sample::qcd) {
            unsigned char hlt = (*hph165) + (*hph175);
            if (not isMC and hlt == 0) continue;
	    if (*phpt < 175. || fabs(*pheta) > 1.4442) continue;
	    trgsf *= triggerphoton->Eval(min(*phpt,triggerphoton->GetXaxis()->GetXmax()));
	    weight *= (1.0 - purhist->GetBinContent(purhist->FindBin(min(*phpt,purhist->GetXaxis()->GetBinLowEdge(purhist->GetNbinsX()+1)-1), fabs(*pheta))));
        }
	
        if (chan == Sample::sig) {
            unsigned char hlt = (*hmnm90) + (*hmnm120) + (*hmwm90) + (*hmwm120) + (*hmwm170) + (*hmwm300);      
            if (not isMC and hlt == 0) continue; 
	    if (not isInclusive) 
	      trgsf *= triggermet->Eval(min(*t1mumet,triggermet->GetXaxis()->GetXmax()));
        }
	
        if (chan == Sample::gam) {
	  unsigned char hlt = (*hph165) + (*hph175);      
	  if (not isMC and hlt == 0) continue; 
	  if (*phpt < 175. || fabs(*pheta) > 1.4442 || *phidm !=1) continue;
	  trgsf *= triggerphoton->Eval(min(*phpt,triggerphoton->GetXaxis()->GetXmax()));
	  effsf *= psfmedium->GetBinContent(psfmedium->FindBin(min(*phpt,psfmedium->GetXaxis()->GetBinLowEdge(psfmedium->GetNbinsX()+1)-1),*pheta));
        }

        if (chan == Sample::wmn) {
	  unsigned char hlt = (*hmnm90) + (*hmnm120) + (*hmwm90) +(*hmwm120) + (*hmwm170) +(*hmwm300);
	  if(isInclusive)
	    hlt += (*hsmu);
	  if (not isMC and hlt == 0) continue; 	  
	  if (*mu1pt < 20 or *mu1id != 1) continue;
	  effsf *= msftight->GetBinContent(msftight->FindBin(min(*mu1pt,msftight->GetXaxis()->GetBinLowEdge(msftight->GetNbinsX()+1)-1),*mu1eta));
	  if(not isInclusive)
	    trgsf *= triggermet->Eval(*t1mumet);
	  else
	    trgsf *= triggermuhist->GetBinContent(triggermuhist->FindBin(min(*mu1pt,triggermuhist->GetXaxis()->GetBinLowEdge(triggermuhist->GetNbinsX()+1)-1),*mu1eta));
	  if(isInclusive and *wmt < 60) continue;
        }        

        if (chan == Sample::zmm) {
	  unsigned char hlt = (*hmnm90) + (*hmnm120) + (*hmwm90) +(*hmwm120) + (*hmwm170) +(*hmwm300);      
	  if(isInclusive) hlt += (*hsmu);
	  if (not isMC and hlt == 0) continue;
	  bool istight = false;
	  if (*mu1id == 1 && *mu1pt > 20.) istight = true;
	  if (*mu2id == 1 && *mu2pt > 20.) istight = true;
	  if (!istight) continue;	 
	  if (*mu1id == 1 ) effsf *= msftight->GetBinContent(msftight->FindBin(min(*mu1pt,msftight->GetXaxis()->GetBinLowEdge(msftight->GetNbinsX()+1)-1),*mu1eta));
	  else              effsf *= msfloose->GetBinContent(msfloose->FindBin(min(*mu1pt,msfloose->GetXaxis()->GetBinLowEdge(msfloose->GetNbinsX()+1)-1),*mu1eta));
	  if (*mu2id == 1 ) effsf *= msftight->GetBinContent(msftight->FindBin(min(*mu2pt,msftight->GetXaxis()->GetBinLowEdge(msftight->GetNbinsX()+1)-1),*mu2eta));
	  else              effsf *= msfloose->GetBinContent(msfloose->FindBin(min(*mu2pt,msfloose->GetXaxis()->GetBinLowEdge(msfloose->GetNbinsX()+1)-1),*mu2eta));

	  if(not isInclusive) trgsf *= triggermet->Eval(*t1mumet);
	  else{
	    if (*mu1id == 1 and *mu2id == 1)  // both tight efficiency is eff1+eff2-eff1*eff2 that at plateau is ~1
	      trgsf *= 1;
	    else if(*mu1id == 1 and *mu2id != 1) 
	      trgsf *= triggermuhist->GetBinContent(triggermuhist->FindBin(min(*mu1pt,triggermuhist->GetXaxis()->GetBinLowEdge(triggermuhist->GetNbinsX()+1)-1),*mu1eta));
	    else if(*mu1id != 1 and *mu2id == 1) 
	      trgsf *= triggermuhist->GetBinContent(triggermuhist->FindBin(min(*mu2pt,triggermuhist->GetXaxis()->GetBinLowEdge(triggermuhist->GetNbinsX()+1)-1),*mu2eta));
	  }
        }
	
        if (chan == Sample::wen) {
	  unsigned char hlt = *hsele;
	  if (not isMC and hlt == 0) continue;
	  if (*el1pt < 40 or *el1id != 1) continue;
	  if (isInclusive and *wemt < 60) continue;
	  effsf *= esftight->GetBinContent(esftight->FindBin(min(*el1pt,esftight->GetXaxis()->GetBinLowEdge(esftight->GetNbinsX()+1)-1),*el1eta));
	  trgsf *= triggerelhist ->GetBinContent(triggerelhist->FindBin(min(*el1pt,triggerelhist->GetXaxis()->GetBinLowEdge(triggerelhist->GetNbinsX()+1)-1),*el1eta));
        }
	
        if (chan == Sample::zee) {	  
	  unsigned char hlt = (*hsele) + (*hph165) + (*hph175);
	  if (not isMC and hlt == 0) continue;
	  bool istight = false;
	  if (*el1id == 1 && *el1pt > 40.) istight = true;
	  if (*el2id == 1 && *el2pt > 40.) istight = true;
	  if (!istight) continue;
	  
	  if (*el1id == 1 ) effsf *= esftight->GetBinContent(esftight->FindBin(min(*el1pt,esftight->GetXaxis()->GetBinLowEdge(esftight->GetNbinsX()+1)-1),*el1eta));
          else              effsf *= esfveto->GetBinContent(esfveto->FindBin(min(*el1pt,esfveto->GetXaxis()->GetBinLowEdge(esfveto->GetNbinsX()+1)-1),*el1eta));
          if (*el2id == 1 ) effsf *= esftight->GetBinContent(esftight->FindBin(min(*el2pt,esftight->GetXaxis()->GetBinLowEdge(esftight->GetNbinsX()+1)-1),*el2eta));
          else              effsf *= esfveto->GetBinContent(esfveto->FindBin(min(*el2pt,esfveto->GetXaxis()->GetBinLowEdge(esfveto->GetNbinsX()+1)-1),*el2eta));
	  
	  if (*el1id == 1 and *el2id == 1)  // both tight efficiency is eff1+eff2-eff1*eff2 that at plateau is ~1
	    trgsf *= 1;
	  else if(*el1id == 1 and *el2id != 1) 
	    trgsf *= triggerelhist->GetBinContent(triggerelhist->FindBin(min(*el1pt,triggerelhist->GetXaxis()->GetBinLowEdge(triggerelhist->GetNbinsX()+1)-1),*el1eta));
	  else if(*el1id != 1 and *el2id == 1) 
	    trgsf *= triggerelhist->GetBinContent(triggerelhist->FindBin(min(*el2pt,triggerelhist->GetXaxis()->GetBinLowEdge(triggerelhist->GetNbinsX()+1)-1),*el2eta)); 
        }
	
	if (chan == Sample::topmu){
	  unsigned char hlt = (*hmnm90) + (*hmnm120) + (*hmwm90) +(*hmwm120) + (*hmwm170) +(*hmwm300);
          if(isInclusive) hlt += (*hsmu);
          if (not isMC and hlt == 0) continue;
	  if (*mu1pt < 20 or *mu1id != 1) continue;
	  if(*nbjets < 1) continue;	  

	  effsf *= msftight->GetBinContent(msftight->FindBin(min(*mu1pt,msftight->GetXaxis()->GetBinLowEdge(msftight->GetNbinsX()+1)-1),*mu1eta));	  
          if(not isInclusive)
	    trgsf *= triggermet->Eval(*t1mumet);
	  else
	    trgsf *= triggermuhist->GetBinContent(triggermuhist->FindBin(min(*mu1pt,triggermuhist->GetXaxis()->GetBinLowEdge(triggermuhist->GetNbinsX()+1)-1),*mu1eta));
	}

	if (chan == Sample::topel){
	  unsigned char hlt = *hsele;
	  if (not isMC and hlt == 0) continue;	  
	  if (*el1pt < 40 or *el1id != 1) continue;
	  if(*nbjets < 1) continue;	  
	  effsf *= esftight->GetBinContent(esftight->FindBin(min(*el1pt,esftight->GetXaxis()->GetBinLowEdge(esftight->GetNbinsX()+1)-1),*el1eta));	  
	  trgsf *= triggerelhist->GetBinContent(triggerelhist->FindBin(min(*el1pt,triggerelhist->GetXaxis()->GetBinLowEdge(triggerelhist->GetNbinsX()+1)-1),*el1eta));
	}

	if ((*nvtx) <= 40) puwgt = puhist->GetBinContent(puhist->FindBin(*nvtx));
	else puwgt = 1.0;

        double genpt = *wzpt;
        if (*wzpt < 100. ) genpt = 100.;
        if (*wzpt > 1000.) genpt = 999.;
        for (size_t i = 0; i < khists.size(); i++) {
            if (khists[i]) {
                kfact = khists[i]->GetBinContent(khists[i]->FindBin(genpt));
            }
        }

	if (isMC) weight *= (*xsec)*(lumi)*(*wgt)*(kfact)*(puwgt)*(trgsf)*(*wgtbtag)*(effsf)/(*wgtsum);	
        is_vector<T> isv;
        Double_t fillval = getValueFromVar(isv, var, index); 
        if (fillval >= max) fillval = mid;
        hist->Fill(fillval, weight);
    }

    pufile .Close();
    sffile_eleTight.Close();
    sffile_eleVeto.Close();
    sffile_muTight.Close();
    sffile_muLoose.Close();
    sffile_phoLoose.Close();
    sffile_phoMedium.Close();
    purityfile_photon.Close();
    triggerfile_SinglEle.Close();
    triggerfile_SingleMu.Close();
    triggerfile_MET.Close();
    triggerfile_SinglePhoton.Close();
}


void makeGenericHist(string sfilename, TH1* hist, const char* varname, const char* cut, bool isMC, Sample chan, double lumi, vector<TH1*> khists, bool isInclusive) {
    string varstr = varname;
    vector<string> varstrvec = split(varstr, '[');
    if (varstrvec.size() == 0) return;
    size_t index = 0;
    if (varstrvec.size() > 1) {
        vector<string> idxstrvec = split(varstrvec[1], ']');
        if (idxstrvec.size() > 0) {
            index = atoi(idxstrvec[0].c_str());
        }
    }

    const char* filename = sfilename.c_str();
    TFile* file = new TFile(filename);
    if(file == 0 or file == NULL or file->IsZombie()){
      cerr<<"File not found --> name :"<<file->GetName()<<" --> return "<<endl;
      return;
    }
    TTree* tree = (TTree*)file->Get("tree/tree");
    TTree* cuttree = NULL;
    if (string(cut) == "") cuttree = tree;
    else {
        gROOT->cd();
        cuttree = tree->CopyTree(cut);
    }

    const char* branchType = "";
    TObjArray* branchList = tree->GetListOfBranches();
    for(Int_t i = 0; i < branchList->GetEntries(); ++i ) {
        TBranch *branch = static_cast<TBranch*>(branchList->At(i));
        if (!branch || !branch->GetLeaf(branch->GetName())) continue;
        if (varstrvec[0] == string(branch->GetName())) branchType = branch->GetLeaf(branch->GetName())->GetTypeName();
    }

    string branchT(branchType);
    if     (branchT == "") return;
    else if(branchT == "int"            || branchT == "Int_t"                              ) templatedMakeHist<Int_t>           (cuttree, hist, varstrvec[0].c_str(), isMC, chan, lumi, khists, index, isInclusive);
    else if(branchT == "unsigned"       || branchT == "unsigned int" || branchT == "UInt_t") templatedMakeHist<UInt_t>          (cuttree, hist, varstrvec[0].c_str(), isMC, chan, lumi, khists, index, isInclusive);
    else if(branchT == "unsigned char"  || branchT == "UChar_t"                            ) templatedMakeHist<UChar_t>         (cuttree, hist, varstrvec[0].c_str(), isMC, chan, lumi, khists, index, isInclusive);
    else if(branchT == "float"          || branchT == "Float_t"                            ) templatedMakeHist<Float_t>         (cuttree, hist, varstrvec[0].c_str(), isMC, chan, lumi, khists, index, isInclusive);
    else if(branchT == "double"         || branchT == "Double_t"                           ) templatedMakeHist<Double_t>        (cuttree, hist, varstrvec[0].c_str(), isMC, chan, lumi, khists, index, isInclusive);
    else if(branchT == "vector<double>" || branchT == "vector<Double_t>"                   ) templatedMakeHist<vector<double> > (cuttree, hist, varstrvec[0].c_str(), isMC, chan, lumi, khists, index, isInclusive);

    file->Close();
}

#endif
