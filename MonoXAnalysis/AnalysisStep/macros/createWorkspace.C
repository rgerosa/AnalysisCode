#include <vector>

bool connectWZ = true;

void addTemplate(string procname, RooArgList& varlist, RooWorkspace& ws, TH1F* hist) {
    RooDataHist rhist(procname.c_str(), "", varlist, hist);
    ws.import(rhist);
}

void makeBinList(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* hist, RooArgList& binlist, bool setConst = false) {
    for (int i = 1; i <= hist->GetNbinsX(); i++) {
        stringstream binss;
        binss << procname << "_bin" << i;
        RooRealVar* binvar;
        if (!setConst) binvar = new RooRealVar(binss.str().c_str(), "", hist->GetBinContent(i), 0., hist->GetBinContent(i)*2.0);
        else           binvar = new RooRealVar(binss.str().c_str(), "", hist->GetBinContent(i));
        binlist.add(*binvar);
    }

    stringstream normss;
    normss << procname << "_norm";

    RooParametricHist phist(procname.c_str(), "", var, binlist, *hist);
    RooAddition norm(normss.str().c_str(), "", binlist);

    ws.import(phist);
    ws.import(norm,RooFit::RecycleConflictNodes());

}

void makeConnectedBinList(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* rhist, vector<TH1*> syst, const RooArgList& srbinlist, RooArgList* crbinlist=NULL) {
    if (crbinlist == NULL) crbinlist = new RooArgList();

    for (int i = 1; i <= rhist->GetNbinsX(); i++) {
        stringstream rbinss;
        rbinss << "r_" << procname << "_bin" << i;
        RooRealVar* rbinvar = new RooRealVar(rbinss.str().c_str(), "", rhist->GetBinContent(i));

        stringstream rerrbinss;
        rerrbinss << procname << "_bin" << i << "_Runc";
        RooRealVar* rerrbinvar = new RooRealVar(rerrbinss.str().c_str(), "", 0., -5., 5.);

        stringstream binss;
        binss << procname << "_bin" << i;

        RooArgList fobinlist;
        fobinlist.add(srbinlist[i-1]);
        fobinlist.add(*rbinvar);
        fobinlist.add(*rerrbinvar);

        stringstream formss;
        formss << "@0/";
        formss << "(";
        formss << "@1";
        formss << "*(1+" << rhist->GetBinError(i)/rhist->GetBinContent(i) << "*@2)";
        for (int j = 0; j < syst.size(); j++) {
            stringstream systbinss;
            systbinss << procname << "_bin" << i << "_" << syst[j]->GetName();
            RooRealVar* systbinvar = new RooRealVar(systbinss.str().c_str(), "", 0., -5., 5.);
            fobinlist.add(*systbinvar);
            formss << "*(1+" << syst[j]->GetBinContent(i) << "*@" << j+3 << ")";
        }
        formss << ")";

        RooFormulaVar* binvar = new RooFormulaVar(binss.str().c_str(), "", formss.str().c_str(), RooArgList(fobinlist));
        crbinlist->add(*binvar);
    }

    stringstream normss;
    normss << procname << "_norm";

    RooParametricHist phist(procname.c_str(), "", var, *crbinlist, *rhist);
    RooAddition norm(normss.str().c_str(),"", *crbinlist);

    ws.import(phist);
    ws.import(norm, RooFit::RecycleConflictNodes());
}

void createWorkspace(){
    gSystem->Load("libHiggsAnalysisCombinedLimit.so");
    
    TFile *outfile = new TFile("workspace.root","RECREATE");
    RooWorkspace wspace("w","w");

    RooRealVar met("met","E_{T}^{miss}",200,1000);
    RooArgList vars(met);

    // Templates
    TFile* templatesfile = new TFile("templates.root");
    TFile* datafile      = new TFile("data.root");

    // ---------------------------- SIGNAL REGION -------------------------------------------------------------------//
    // Data
    addTemplate("data_obs_SR", vars, wspace, (TH1F*)templatesfile->Get("datahist"));

    // Signal shape
    addTemplate("DM_SR", vars, wspace, (TH1F*)templatesfile->Get("sig1hist"));
    
    // Znunu background
    TH1F* znn_SR_hist = (TH1F*)templatesfile->Get("zinvhist");
    RooArgList znn_SR_bins;
    makeBinList("Znunu_SR", met, wspace, znn_SR_hist, znn_SR_bins);

    // WJets background
    TH1F* wln_SR_hist = (TH1F*)templatesfile->Get("wjethist");
    RooArgList wln_SR_bins;
    vector<TH1*> wln_SR_syst;
    wln_SR_syst.push_back((TH1F*)templatesfile->Get("ZW_Theory"));
    if (!connectWZ) makeBinList("WJets_SR", met, wspace, wln_SR_hist, wln_SR_bins);
    else   makeConnectedBinList("WJets_SR", met, wspace, (TH1F*)templatesfile->Get("zwjcorhist"), wln_SR_syst, znn_SR_bins, &wln_SR_bins);

    // Other MC backgrounds
    addTemplate("ZJets_SR"     , vars, wspace, (TH1F*)templatesfile->Get("zjethist"));
    addTemplate("Top_SR"       , vars, wspace, (TH1F*)templatesfile->Get("tbkghist"));
    addTemplate("QCD_SR"       , vars, wspace, (TH1F*)templatesfile->Get("qbkghist"));
    addTemplate("Dibosons_SR"  , vars, wspace, (TH1F*)templatesfile->Get("dbkghist"));

    // ---------------------------- CONTROL REGION (Dimuon) -----------------------------------------------------------------//
    addTemplate("data_obs_ZM", vars, wspace, (TH1F*)templatesfile->Get("datahistzmm"));
    vector<TH1*>   znn_ZM_syst;
    makeConnectedBinList("Znunu_ZM", met, wspace, (TH1F*)templatesfile->Get("zmmcorhist"), znn_ZM_syst, znn_SR_bins);

    // Other MC backgrounds in dimuon control region
    addTemplate("WJets_ZM"     , vars, wspace, (TH1F*)templatesfile->Get("lbkghistzmm"));
    addTemplate("Top_ZM"       , vars, wspace, (TH1F*)templatesfile->Get("tbkghistzmm"));
    addTemplate("QCD_ZM"       , vars, wspace, (TH1F*)templatesfile->Get("qbkghistzmm"));
    addTemplate("Dibosons_ZM"  , vars, wspace, (TH1F*)templatesfile->Get("dbkghistzmm"));

    // ---------------------------- CONTROL REGION (Dielectron) -----------------------------------------------------------------//
    addTemplate("data_obs_ZE"  , vars, wspace, (TH1F*)templatesfile->Get("datahistzee"));
    vector<TH1*> znn_ZE_syst;
    makeConnectedBinList("Znunu_ZE", met, wspace, (TH1F*)templatesfile->Get("zeecorhist"), znn_ZE_syst, znn_SR_bins);

    // Other MC backgrounds in dielectron control region
    addTemplate("WJets_ZE"     , vars, wspace, (TH1F*)templatesfile->Get("lbkghistzee"));
    addTemplate("Top_ZE"       , vars, wspace, (TH1F*)templatesfile->Get("tbkghistzee"));
    addTemplate("QCD_ZE"       , vars, wspace, (TH1F*)templatesfile->Get("qbkghistzee"));
    addTemplate("Dibosons_ZE"  , vars, wspace, (TH1F*)templatesfile->Get("dbkghistzee"));

    // ---------------------------- CONTROL REGION (Photon+Jets) -----------------------------------------------------------------//
    addTemplate("data_obs_GJ"  , vars, wspace, (TH1F*)templatesfile->Get("datahistgam"));
    vector<TH1*> znn_GJ_syst;
    znn_GJ_syst.push_back((TH1F*)templatesfile->Get("ZG_Theory"));
    makeConnectedBinList("Znunu_GJ", met, wspace, (TH1F*)templatesfile->Get("gamcorhist"), znn_GJ_syst, znn_SR_bins);

    // Other MC backgrounds photon+jets control region
    addTemplate("QCD_GJ"     , vars, wspace, (TH1F*)templatesfile->Get("qbkghistgam"));

    // ---------------------------- CONTROL REGION (Single muon) -----------------------------------------------------------------//
    addTemplate("data_obs_WM"  , vars, wspace, (TH1F*)templatesfile->Get("datahistwmn"));
    vector<TH1*> wln_WM_syst;
    makeConnectedBinList("WJets_WM", met, wspace, (TH1F*)templatesfile->Get("wmncorhist"), wln_WM_syst, wln_SR_bins);

    // Other MC backgrounds in single muon control region
    addTemplate("ZJets_WM"     , vars, wspace, (TH1F*)templatesfile->Get("lbkghistwmn"));
    addTemplate("Top_WM"       , vars, wspace, (TH1F*)templatesfile->Get("tbkghistwmn"));
    addTemplate("QCD_WM"       , vars, wspace, (TH1F*)templatesfile->Get("qbkghistwmn"));
    addTemplate("Dibosons_WM"  , vars, wspace, (TH1F*)templatesfile->Get("dbkghistwmn"));

    // ---------------------------- CONTROL REGION (Single electron) -----------------------------------------------------------------//
    addTemplate("data_obs_WE"  , vars, wspace, (TH1F*)templatesfile->Get("datahistwen"));
    vector<TH1*> wln_WE_syst;
    makeConnectedBinList("WJets_WE", met, wspace, (TH1F*)templatesfile->Get("wencorhist"), wln_WE_syst, wln_SR_bins);

    // Other MC backgrounds in single electron control region
    addTemplate("ZJets_WE"     , vars, wspace, (TH1F*)templatesfile->Get("lbkghistwen"));
    addTemplate("Top_WE"       , vars, wspace, (TH1F*)templatesfile->Get("tbkghistwen"));
    addTemplate("QCD_WE"       , vars, wspace, (TH1F*)templatesfile->Get("qbkghistwen"));
    addTemplate("Dibosons_WE"  , vars, wspace, (TH1F*)templatesfile->Get("dbkghistwen"));

    // ---------------------------- Write out the workspace -----------------------------------------------------------------//
    outfile->cd();
    wspace.Write();
    outfile->Close();
}
