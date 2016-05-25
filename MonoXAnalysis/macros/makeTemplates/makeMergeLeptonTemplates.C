void makeMergeLeptonTemplates(string templateFile, string zmmFile, string zeeFile, string wmnFile, string wenFile, string outputFile, vector<string> obs, 
			      bool doShapeSys = false){

  vector<string> shapeSys = {"bUp","bDw","metJetUp","metJetDw","metResUp","metResDw","metUncUp","metUncDw"};

  TFile* inputTemplate = TFile::Open(templateFile.c_str());
  TFile* zmmTemplate = TFile::Open(zmmFile.c_str());
  TFile* zeeTemplate = TFile::Open(zeeFile.c_str());
  TFile* wmnTemplate = TFile::Open(wmnFile.c_str());
  TFile* wenTemplate = TFile::Open(wenFile.c_str());
    
  system(("cp "+templateFile+" "+outputFile).c_str());
  TFile* output = TFile::Open(outputFile.c_str(),"UPDATE");

  // make central value
  for(auto observable : obs){
    TH1F* zmm_num = (TH1F*) zmmTemplate->Get(("nhist_"+observable).c_str());
    TH1F* zmm_den = (TH1F*) zmmTemplate->Get(("dhist_"+observable).c_str());
    TH1F* zee_num = (TH1F*) zeeTemplate->Get(("nhist_"+observable).c_str());
    TH1F* zee_den = (TH1F*) zeeTemplate->Get(("dhist_"+observable).c_str());
    TH1F* wmn_num = (TH1F*) wmnTemplate->Get(("nhist_"+observable).c_str());
    TH1F* wmn_den = (TH1F*) wmnTemplate->Get(("dhist_"+observable).c_str());
    TH1F* wen_num = (TH1F*) wenTemplate->Get(("nhist_"+observable).c_str());
    TH1F* wen_den = (TH1F*) wenTemplate->Get(("dhist_"+observable).c_str());
    
    output->cd();
    if(not output->GetDirectory("TF_ZL"))
      output->mkdir("TF_ZL");
    output->cd("TF_ZL");
    
    // cloning TFs
    TH1F* zll_tf = (TH1F*) zmm_num->Clone(("zllcorhist_"+observable).c_str());
    zmm_den->Add(zee_den);
    zll_tf->Divide(zmm_den);
    zll_tf->Write();
    
    output->cd();
    if(not output->GetDirectory("TF_WL"))
      output->mkdir("TF_WL");
    output->cd("TF_WL");
    TH1F* wln_tf = (TH1F*) wmn_num->Clone(("wlncorhist_"+observable).c_str());
    wmn_den->Add(wen_den);
    wln_tf->Divide(wmn_den);
    wln_tf->Write();

    //Cloning and adding templates 
    output->cd();
    if(not output->GetDirectory("ZL"))
      output->mkdir("ZL");
    output->cd("ZL");
    
    TH1F* data_zmm = (TH1F*) inputTemplate->FindObjectAny(("datahistzmm_"+observable).c_str());
    TH1F* data_zee = (TH1F*) inputTemplate->FindObjectAny(("datahistzee_"+observable).c_str());
    TH1F* data_zll = (TH1F*) data_zmm->Clone(TString(data_zmm->GetName()).ReplaceAll("zmm","zll"));
    data_zll->Add(data_zee);
    data_zll->Write();
    
    TH1F* vl_zmm = (TH1F*) inputTemplate->FindObjectAny(("vlbkghistzmm_"+observable).c_str());
    TH1F* vl_zee = (TH1F*) inputTemplate->FindObjectAny(("vlbkghistzee_"+observable).c_str());
    TH1F* vl_zll = (TH1F*) vl_zmm->Clone(TString(vl_zmm->GetName()).ReplaceAll("zmm","zll"));
    vl_zll->Add(vl_zee);
    vl_zll->Write();
    
    TH1F* vll_zmm = (TH1F*) inputTemplate->FindObjectAny(("vllbkghistzmm_"+observable).c_str());
    TH1F* vll_zee = (TH1F*) inputTemplate->FindObjectAny(("vllbkghistzee_"+observable).c_str());
    TH1F* vll_zll = (TH1F*) vll_zmm->Clone(TString(vll_zmm->GetName()).ReplaceAll("zmm","zll"));
    vll_zll->Add(vll_zee);
    vll_zll->Write();
    
    TH1F* qc_zmm = (TH1F*) inputTemplate->FindObjectAny(("qbkghistzmm_"+observable).c_str());
    TH1F* qc_zee = (TH1F*) inputTemplate->FindObjectAny(("qbkghistzee_"+observable).c_str());
    TH1F* qc_zll = (TH1F*) qc_zmm->Clone(TString(qc_zmm->GetName()).ReplaceAll("zmm","zll"));
    qc_zll->Add(qc_zee);
    qc_zll->Write();
    
    TH1F* db_zmm = (TH1F*) inputTemplate->FindObjectAny(("dbkghistzmm_"+observable).c_str());
    TH1F* db_zee = (TH1F*) inputTemplate->FindObjectAny(("dbkghistzee_"+observable).c_str());
    TH1F* db_zll = (TH1F*) db_zmm->Clone(TString(db_zmm->GetName()).ReplaceAll("zmm","zll"));
    db_zll->Add(db_zee);
    db_zll->Write();
    
    TH1F* gm_zmm = (TH1F*) inputTemplate->FindObjectAny(("gbkghistzmm_"+observable).c_str());
    TH1F* gm_zee = (TH1F*) inputTemplate->FindObjectAny(("gbkghistzee_"+observable).c_str());
    TH1F* gm_zll = (TH1F*) gm_zmm->Clone(TString(gm_zmm->GetName()).ReplaceAll("zmm","zll"));
    gm_zll->Add(gm_zee);
    gm_zll->Write();
    
    TH1F* tt_zmm = (TH1F*) inputTemplate->FindObjectAny(("tbkghistzmm_"+observable).c_str());
    TH1F* tt_zee = (TH1F*) inputTemplate->FindObjectAny(("tbkghistzee_"+observable).c_str());
    TH1F* tt_zll = (TH1F*) tt_zmm->Clone(TString(tt_zmm->GetName()).ReplaceAll("zmm","zll"));
    tt_zll->Add(tt_zee);
    tt_zll->Write();
    
    output->cd();
    if(not output->GetDirectory("ZL/shapeSys"))
      output->mkdir("ZL/shapeSys");
    output->cd("ZL/shapeSys");
    for(auto sys : shapeSys){
      TH1F* vl_zmm = (TH1F*) inputTemplate->FindObjectAny(("vlbkghistzmm_"+sys+"_"+observable).c_str());
      TH1F* vl_zee = (TH1F*) inputTemplate->FindObjectAny(("vlbkghistzee_"+sys+"_"+observable).c_str());
      TH1F* vl_zll = (TH1F*) vl_zmm->Clone(TString(vl_zmm->GetName()).ReplaceAll("zmm","zll"));
      vl_zll->Add(vl_zee);
      vl_zll->Write();
      
      TH1F* db_zmm = (TH1F*) inputTemplate->FindObjectAny(("dbkghistzmm_"+sys+"_"+observable).c_str());
      TH1F* db_zee = (TH1F*) inputTemplate->FindObjectAny(("dbkghistzee_"+sys+"_"+observable).c_str());
      TH1F* db_zll = (TH1F*) db_zmm->Clone(TString(db_zmm->GetName()).ReplaceAll("zmm","zll"));
      db_zll->Add(db_zee);
      db_zll->Write();
      
      TH1F* gm_zmm = (TH1F*) inputTemplate->FindObjectAny(("gbkghistzmm_"+sys+"_"+observable).c_str());
      TH1F* gm_zee = (TH1F*) inputTemplate->FindObjectAny(("gbkghistzee_"+sys+"_"+observable).c_str());
      TH1F* gm_zll = (TH1F*) gm_zmm->Clone(TString(gm_zmm->GetName()).ReplaceAll("zmm","zll"));
      gm_zll->Add(gm_zee);
      gm_zll->Write();
      
      TH1F* tt_zmm = (TH1F*) inputTemplate->FindObjectAny(("tbkghistzmm_"+sys+"_"+observable).c_str());
      TH1F* tt_zee = (TH1F*) inputTemplate->FindObjectAny(("tbkghistzee_"+sys+"_"+observable).c_str());
      TH1F* tt_zll = (TH1F*) tt_zmm->Clone(TString(tt_zmm->GetName()).ReplaceAll("zmm","zll"));
      tt_zll->Add(tt_zee);
      tt_zll->Write();
    }
    
    
    //Cloning and adding templates 
    output->cd();
    if(not output->GetDirectory("WL"))
      output->mkdir("WL");
    output->cd("WL");
    
    TH1F* data_wmn = (TH1F*) inputTemplate->FindObjectAny(("datahistwmn_"+observable).c_str());
    TH1F* data_wen = (TH1F*) inputTemplate->FindObjectAny(("datahistwen_"+observable).c_str());
    TH1F* data_wln = (TH1F*) data_wmn->Clone(TString(data_wmn->GetName()).ReplaceAll("wmn","wln"));
    data_wln->Add(data_wen);
    data_wln->Write();
    
    TH1F* vl_wmn = (TH1F*) inputTemplate->FindObjectAny(("vlbkghistwmn_"+observable).c_str());
    TH1F* vl_wen = (TH1F*) inputTemplate->FindObjectAny(("vlbkghistwen_"+observable).c_str());
    TH1F* vl_wln = (TH1F*) vl_wmn->Clone(TString(vl_wmn->GetName()).ReplaceAll("wmn","wln"));
    vl_wln->Add(vl_wen);
    vl_wln->Write();
    
    TH1F* vll_wmn = (TH1F*) inputTemplate->FindObjectAny(("vllbkghistwmn_"+observable).c_str());
    TH1F* vll_wen = (TH1F*) inputTemplate->FindObjectAny(("vllbkghistwen_"+observable).c_str());
    TH1F* vll_wln = (TH1F*) vll_wmn->Clone(TString(vll_wmn->GetName()).ReplaceAll("wmn","wln"));
    vll_wln->Add(vll_wen);
    vll_wln->Write();

    TH1F* qc_wmn = (TH1F*) inputTemplate->FindObjectAny(("qbkghistwmn_"+observable).c_str());
    TH1F* qc_wen = (TH1F*) inputTemplate->FindObjectAny(("qbkghistwen_"+observable).c_str());
    TH1F* qc_wln = (TH1F*) qc_wmn->Clone(TString(qc_wmn->GetName()).ReplaceAll("wmn","wln"));
    qc_wln->Add(qc_wen);
    qc_wln->Write();
    
    TH1F* db_wmn = (TH1F*) inputTemplate->FindObjectAny(("dbkghistwmn_"+observable).c_str());
    TH1F* db_wen = (TH1F*) inputTemplate->FindObjectAny(("dbkghistwen_"+observable).c_str());
    TH1F* db_wln = (TH1F*) db_wmn->Clone(TString(db_wmn->GetName()).ReplaceAll("wmn","wln"));
    db_wln->Add(db_wen);
    db_wln->Write();
  
    TH1F* gm_wmn = (TH1F*) inputTemplate->FindObjectAny(("gbkghistwmn_"+observable).c_str());
    TH1F* gm_wen = (TH1F*) inputTemplate->FindObjectAny(("gbkghistwen_"+observable).c_str());
    TH1F* gm_wln = (TH1F*) gm_wmn->Clone(TString(gm_wmn->GetName()).ReplaceAll("wmn","wln"));
    gm_wln->Add(gm_wen);
    gm_wln->Write();
    
    TH1F* tt_wmn = (TH1F*) inputTemplate->FindObjectAny(("tbkghistwmn_"+observable).c_str());
    TH1F* tt_wen = (TH1F*) inputTemplate->FindObjectAny(("tbkghistwen_"+observable).c_str());
    TH1F* tt_wln = (TH1F*) tt_wmn->Clone(TString(tt_wmn->GetName()).ReplaceAll("wmn","wln"));
    tt_wln->Add(tt_wen);
    tt_wln->Write();
    
    output->cd();
    if(not output->GetDirectory("WL/shapeSys"))
      output->mkdir("WL/shapeSys");
    output->cd("WL/shapeSys");
    for(auto sys : shapeSys){
      TH1F* vll_wmn = (TH1F*) inputTemplate->FindObjectAny(("vllbkghistwmn_"+sys+"_"+observable).c_str());
      TH1F* vll_wen = (TH1F*) inputTemplate->FindObjectAny(("vllbkghistwen_"+sys+"_"+observable).c_str());
      TH1F* vll_wln = (TH1F*) vll_wmn->Clone(TString(vll_wmn->GetName()).ReplaceAll("wmn","wln"));
      vll_wln->Add(vll_wen);
      vll_wln->Write();
      
      TH1F* db_wmn = (TH1F*) inputTemplate->FindObjectAny(("dbkghistwmn_"+sys+"_"+observable).c_str());
      TH1F* db_wen = (TH1F*) inputTemplate->FindObjectAny(("dbkghistwen_"+sys+"_"+observable).c_str());
      TH1F* db_wln = (TH1F*) db_wmn->Clone(TString(db_wmn->GetName()).ReplaceAll("wmn","wln"));
      db_wln->Add(db_wen);
      db_wln->Write();

      TH1F* gm_wmn = (TH1F*) inputTemplate->FindObjectAny(("gbkghistwmn_"+sys+"_"+observable).c_str());
      TH1F* gm_wen = (TH1F*) inputTemplate->FindObjectAny(("gbkghistwen_"+sys+"_"+observable).c_str());
      TH1F* gm_wln = (TH1F*) gm_wmn->Clone(TString(gm_wmn->GetName()).ReplaceAll("wmn","wln"));
      gm_wln->Add(gm_wen);
      gm_wln->Write();
    
      TH1F* tt_wmn = (TH1F*) inputTemplate->FindObjectAny(("tbkghistwmn_"+sys+"_"+observable).c_str());
      TH1F* tt_wen = (TH1F*) inputTemplate->FindObjectAny(("tbkghistwen_"+sys+"_"+observable).c_str());
      TH1F* tt_wln = (TH1F*) tt_wmn->Clone(TString(tt_wmn->GetName()).ReplaceAll("wmn","wln"));
      tt_wln->Add(tt_wen);
      tt_wln->Write();
    }
    
  }

  
  output->Close();
  
}
