#ifndef MAKEHIST_H
#define MAKEHIST_H

enum class Sample { qcd, sig, gam, wmn, zmm, wen, zee};

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
void templatedMakeHist(TTree* tree, TH1* hist, const char* varstr, bool isMC, Sample chan, double lumi, vector<TH1*> khists, size_t index) {

    TTreeReader reader(tree);

    TFile pufile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt.root");
    TH1*  puhist = (TH1*)pufile.Get("puhist");

    TFile sffile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/leptonSF/leptonIDsfs.root");
    TH2*  msflhist = (TH2*)sffile.Get("muon_loose_SF");
    TH2*  msfthist = (TH2*)sffile.Get("muon_tight_SF");
    TH2*  esflhist = (TH2*)sffile.Get("electron_veto_SF");
    TH2*  esfthist = (TH2*)sffile.Get("electron_tight_SF");

    TFile psffile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/PhotonSFandEffandPurity_Lumi2p1fb_2211.root");
    TH2*  psfhist = (TH2*)psffile.Get("PhotonSF");
    TH2*  purhist = (TH2*)psffile.Get("PhotonPurity");

    TFile trefile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF/leptonTrigsfs.root");
    TH2*  trehist = (TH2*)trefile.Get("hltel27_SF");

    TFile trmfile("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF/mettrigSF.root");
    TH1*  trmhist = (TH1*)trmfile.Get("mettrigSF");

    const char* wgtsumvar;
    const char* wgtpuvar;

    if (isMC)   {
        wgtsumvar = "wgtsum";
        wgtpuvar  = "wgtpileup";
    }
    else {
        wgtsumvar = "wgt";
        wgtpuvar  = "wgt";
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
    TTreeReaderValue<double>          wgt    (reader, "wgt");
    TTreeReaderValue<double>          xsec   (reader, "xsec");
    TTreeReaderValue<double>          wzpt   (reader, "wzpt");

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
    TTreeReaderValue<double>          phpt   (reader, "phpt");
    TTreeReaderValue<double>          pheta  (reader, "pheta");

    Double_t max = hist->GetBinLowEdge(hist->GetNbinsX()) + hist->GetBinWidth(hist->GetNbinsX());
    Double_t mid = hist->GetBinLowEdge(hist->GetNbinsX()) + hist->GetBinWidth(hist->GetNbinsX())/ 2.0;

    while (reader.Next()) {
        double weight = 1.0;
        double kfact  = 1.0;
        double puwgt  = 0.0;
        double effsf  = 1.0;
        double trgsf  = 1.0;

        if (!isMC && *fcsc    == 0) continue;
        if (!isMC && *fhbhe   == 0) continue;
        if (!isMC && *fhbhei  == 0) continue;
        if (!isMC && *feesc   == 0) continue;
        if (!isMC && *fecal   == 0) continue;
        if (!isMC && *fnvtx   == 0) continue;

        if (*njets  < 1) continue;
        if (*nbjets > 0) continue;
        if ((*chfrac)[0] < 0.1 ) continue;
        if ((*nhfrac)[0] > 0.8 ) continue;
        if ((*jetpt )[0] < 100.) continue;

        double jetmetdphi = *jetmetm;
        if (chan == Sample::qcd || chan == Sample::gam) jetmetdphi = *jetmetp; 
        if (chan == Sample::wen || chan == Sample::zee) jetmetdphi = *jetmete; 
        if (jetmetdphi < 0.5) continue;

        double recoil = *t1mumet;
        if (chan == Sample::qcd || chan == Sample::gam) recoil = *t1phmet;
        if (chan == Sample::wen || chan == Sample::zee) recoil = *t1elmet; 
        if (recoil < 200.) continue;

        if (chan == Sample::qcd) {
            unsigned char hlt = (*hph165) + (*hph175);
            if (hlt == 0) continue;
            if (*phpt < 175. || fabs(*pheta) > 1.4442) continue;
        }

        if (chan == Sample::sig) {
            unsigned char hlt = (*hmnm90) + (*hmnm120) + (*hmwm90) + (*hmwm120) + (*hmwm170) + (*hmwm300);      
            if (hlt == 0) continue; 
            trgsf *= trmhist->GetBinContent(trmhist->FindBin(min(999., recoil)));
        }

        if (chan == Sample::gam) {
            unsigned char hlt = (*hph165) + (*hph175);      
            if (hlt == 0) continue; 
            if (*phpt < 175. || fabs(*pheta) > 1.4442) continue;
        }

        if (chan == Sample::wmn) {
            unsigned char hlt = (*hmnm90) + (*hmnm120);      
            if (hlt == 0) continue; 
            effsf *= msfthist->GetBinContent(msfthist->FindBin(min(999., *mu1pt), fabs(*mu1eta)));
            trgsf *= trmhist ->GetBinContent(trmhist ->FindBin(min(999., recoil)));
        }        

        if (chan == Sample::zmm) {
            unsigned char hlt = (*hmnm90) + (*hmnm120);
            if (hlt == 0) continue;
            bool istight = false;
            if (*mu1id == 1 && *mu1pt > 20.) istight = true;
            if (*mu2id == 1 && *mu2pt > 20.) istight = true;
            if (!istight) continue;
            if (*mu1id == 1) effsf *= msfthist->GetBinContent(msfthist->FindBin(min(999., *mu1pt), fabs(*mu1eta)));
            else             effsf *= msflhist->GetBinContent(msflhist->FindBin(min(999., *mu1pt), fabs(*mu1eta)));
            if (*mu2id == 1) effsf *= msfthist->GetBinContent(msfthist->FindBin(min(999., *mu2pt), fabs(*mu2eta)));
            else             effsf *= msflhist->GetBinContent(msflhist->FindBin(min(999., *mu2pt), fabs(*mu2eta)));
            trgsf *= trmhist->GetBinContent(trmhist->FindBin(min(999., recoil)));
        }

        if (chan == Sample::wen) {
            unsigned char hlt = *hsele;
            if (hlt == 0) continue;
            effsf *= esfthist->GetBinContent(esfthist->FindBin(min(999., *el1pt), fabs(*el1eta)));
            trgsf *= trehist ->GetBinContent(trehist ->FindBin(min(999., *el1pt), fabs(*el1eta)));
        }

        if (chan == Sample::zee) {
            unsigned char hlt = (*hsele) + (*hph165) + (*hph175);
            if (hlt == 0) continue;
            bool istight = false;
            if (*el1id == 1 && *el1pt > 40.) istight = true;
            if (*el2id == 1 && *el2pt > 40.) istight = true;
            if (!istight) continue;
            if (*el1id == 1) effsf *= esfthist->GetBinContent(esfthist->FindBin(min(999., *el1pt), fabs(*el1eta)));
            else             effsf *= esflhist->GetBinContent(esflhist->FindBin(min(999., *el1pt), fabs(*el1eta)));
            if (*el2id == 1) effsf *= esfthist->GetBinContent(esfthist->FindBin(min(999., *el2pt), fabs(*el2eta)));
            else             effsf *= esflhist->GetBinContent(esflhist->FindBin(min(999., *el2pt), fabs(*el2eta)));
        }

        //if ((*nvtx) <= 35) puwgt = puhist->GetBinContent(*nvtx);
        puwgt = 1.0;
        double genpt = *wzpt;
        if (*wzpt < 100. ) genpt = 100.;
        if (*wzpt > 1000.) genpt = 999.;
        for (size_t i = 0; i < khists.size(); i++) {
            if (khists[i]) {
                kfact = khists[i]->GetBinContent(khists[i]->FindBin(genpt));
            }
        }

        if (isMC) weight = (*xsec)*(lumi)*(*wgt)*(kfact)*(puwgt)*(trgsf)*(effsf)/(*wgtsum);
        if (chan == Sample::qcd) weight *= (1.0 - purhist->GetBinContent(purhist->FindBin(min(999., *phpt), fabs(*pheta))));
        is_vector<T> isv;
        Double_t fillval = getValueFromVar(isv, var, index); 
        if (fillval >= max) fillval = mid;
        hist->Fill(fillval, weight);
    }

    pufile .Close();
    sffile .Close();
    psffile.Close();
    trefile.Close();
    trmfile.Close();
}


void makeGenericHist(string sfilename, TH1* hist, const char* varname, const char* cut, bool isMC, Sample chan, double lumi, vector<TH1*> khists) {
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
    else if(branchT == "int"            || branchT == "Int_t"                              ) templatedMakeHist<Int_t>           (cuttree, hist, varstrvec[0].c_str(), isMC, chan, lumi, khists, index);
    else if(branchT == "unsigned"       || branchT == "unsigned int" || branchT == "UInt_t") templatedMakeHist<UInt_t>          (cuttree, hist, varstrvec[0].c_str(), isMC, chan, lumi, khists, index);
    else if(branchT == "unsigned char"  || branchT == "UChar_t"                            ) templatedMakeHist<UChar_t>         (cuttree, hist, varstrvec[0].c_str(), isMC, chan, lumi, khists, index);
    else if(branchT == "float"          || branchT == "Float_t"                            ) templatedMakeHist<Float_t>         (cuttree, hist, varstrvec[0].c_str(), isMC, chan, lumi, khists, index);
    else if(branchT == "double"         || branchT == "Double_t"                           ) templatedMakeHist<Double_t>        (cuttree, hist, varstrvec[0].c_str(), isMC, chan, lumi, khists, index);
    else if(branchT == "vector<double>" || branchT == "vector<Double_t>"                   ) templatedMakeHist<vector<double> > (cuttree, hist, varstrvec[0].c_str(), isMC, chan, lumi, khists, index);

    file->Close();
}

#endif
