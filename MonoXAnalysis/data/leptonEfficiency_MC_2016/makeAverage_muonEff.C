void makeAverage(){

  TFile* inputBCDEF_id = TFile::Open("~/Downloads/EfficienciesAndSF_BCDEF.root");
  TFile* inputGH_id = TFile::Open("~/Downloads/EfficienciesAndSF_GH.root");
  TFile* inputBCDEF_iso = TFile::Open("~/Downloads/EfficienciesAndSF_BCDEF-2.root");
  TFile* inputGH_iso = TFile::Open("~/Downloads/EfficienciesAndSF_GH-2.root");

  TH2F* efficiency_BCDEF_loose_id = (TH2F*) inputBCDEF_id->Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/abseta_pt_MC");
  TH2F* efficiency_GH_loose_id = (TH2F*) inputGH_id->Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/efficienciesMC/abseta_pt_MC");
  TH2F* efficiency_BCDEF_tight_id = (TH2F*) inputBCDEF_id->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/abseta_pt_MC");
  TH2F* efficiency_GH_tight_id = (TH2F*) inputGH_id->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/efficienciesMC/abseta_pt_MC");

  TH2F* efficiency_BCDEF_loose_iso = (TH2F*) inputBCDEF_iso->Get("LooseISO_LooseID_pt_eta/efficienciesMC/abseta_pt_MC");
  TH2F* efficiency_GH_loose_iso = (TH2F*) inputGH_iso->Get("LooseISO_LooseID_pt_eta/efficienciesMC/abseta_pt_MC");
  TH2F* efficiency_BCDEF_tight_iso = (TH2F*) inputBCDEF_iso->Get("TightISO_TightID_pt_eta/efficienciesMC/abseta_pt_MC");
  TH2F* efficiency_GH_tight_iso = (TH2F*) inputGH_iso->Get("TightISO_TightID_pt_eta/efficienciesMC/abseta_pt_MC");
  
  float lumi_BCDEF = 19.72;
  float lumi_GH = 16.14;
  float lumi_total = 35.9;

  if(efficiency_BCDEF_loose_id->GetNbinsX() != efficiency_GH_loose_id->GetNbinsX())
    cout<<"Problem i.e. different x-axis binning for loose id "<<endl;
  if(efficiency_BCDEF_loose_id->GetNbinsY() != efficiency_GH_loose_id->GetNbinsY())
    cout<<"Problem i.e. different y-axis binning for loose id "<<endl;

  TH2F* efficiency_MC_loose_id = (TH2F*) efficiency_BCDEF_loose_id->Clone("efficiency_MC_loose_id");

  for(int iBinX = 0; iBinX < efficiency_BCDEF_loose_id->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < efficiency_BCDEF_loose_id->GetNbinsY(); iBinY++){
      efficiency_MC_loose_id->SetBinContent(iBinX+1,iBinY+1,
					    (efficiency_BCDEF_loose_id->GetBinContent(iBinX+1,iBinY+1)*(lumi_BCDEF/lumi_total)+efficiency_GH_loose_id->GetBinContent(iBinX+1,iBinY+1)*(lumi_GH/lumi_total))/(lumi_BCDEF/lumi_total+lumi_GH/lumi_total));            
      cout<<"Loose id : BCDEF "<<efficiency_BCDEF_loose_id->GetBinContent(iBinX+1,iBinY+1)<<" GH "<<efficiency_GH_loose_id->GetBinContent(iBinX+1,iBinY+1)<<" mean "<<efficiency_MC_loose_id->GetBinContent(iBinX+1,iBinY+1)<<endl;
    }
  }

  if(efficiency_BCDEF_tight_id->GetNbinsX() != efficiency_GH_tight_id->GetNbinsX())
    cout<<"Problem i.e. different x-axis binning for tight id "<<endl;
  if(efficiency_BCDEF_tight_id->GetNbinsY() != efficiency_GH_tight_id->GetNbinsY())
    cout<<"Problem i.e. different y-axis binning for tight id "<<endl;

  TH2F* efficiency_MC_tight_id = (TH2F*) efficiency_BCDEF_tight_id->Clone("efficiency_MC_tight_id");

  for(int iBinX = 0; iBinX < efficiency_BCDEF_tight_id->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < efficiency_BCDEF_tight_id->GetNbinsY(); iBinY++){
      efficiency_MC_tight_id->SetBinContent(iBinX+1,iBinY+1,
					    (efficiency_BCDEF_tight_id->GetBinContent(iBinX+1,iBinY+1)*(lumi_BCDEF/lumi_total)+efficiency_GH_tight_id->GetBinContent(iBinX+1,iBinY+1)*(lumi_GH/lumi_total))/(lumi_BCDEF/lumi_total+lumi_GH/lumi_total));
      
      cout<<"Tight id : BCDEF "<<efficiency_BCDEF_tight_id->GetBinContent(iBinX+1,iBinY+1)<<" GH "<<efficiency_GH_tight_id->GetBinContent(iBinX+1,iBinY+1)<<" mean "<<efficiency_MC_tight_id->GetBinContent(iBinX+1,iBinY+1)<<endl;
    }
  }
  ///////////

  if(efficiency_BCDEF_loose_iso->GetNbinsX() != efficiency_GH_loose_iso->GetNbinsX())
    cout<<"Problem i.e. different x-axis binning for loose iso "<<endl;
  if(efficiency_BCDEF_loose_iso->GetNbinsY() != efficiency_GH_loose_iso->GetNbinsY())
    cout<<"Problem i.e. different y-axis binning for loose iso "<<endl;

  TH2F* efficiency_MC_loose_iso = (TH2F*) efficiency_BCDEF_loose_iso->Clone("efficiency_MC_loose_iso");

  for(int iBinX = 0; iBinX < efficiency_BCDEF_loose_iso->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < efficiency_BCDEF_loose_iso->GetNbinsY(); iBinY++){
      efficiency_MC_loose_iso->SetBinContent(iBinX+1,iBinY+1,
					    (efficiency_BCDEF_loose_iso->GetBinContent(iBinX+1,iBinY+1)*(lumi_BCDEF/lumi_total)+efficiency_GH_loose_iso->GetBinContent(iBinX+1,iBinY+1)*(lumi_GH/lumi_total))/(lumi_BCDEF/lumi_total+lumi_GH/lumi_total));
    
      cout<<"Loose iso : BCDEF "<<efficiency_BCDEF_loose_iso->GetBinContent(iBinX+1,iBinY+1)<<" GH "<<efficiency_GH_loose_iso->GetBinContent(iBinX+1,iBinY+1)<<" mean "<<efficiency_MC_loose_iso->GetBinContent(iBinX+1,iBinY+1)<<endl;
    }
  }

  if(efficiency_BCDEF_tight_iso->GetNbinsX() != efficiency_GH_tight_iso->GetNbinsX())
    cout<<"Problem i.e. different x-axis binning for tight iso "<<endl;
  if(efficiency_BCDEF_tight_iso->GetNbinsY() != efficiency_GH_tight_iso->GetNbinsY())
    cout<<"Problem i.e. different y-axis binning for tight iso "<<endl;

  TH2F* efficiency_MC_tight_iso = (TH2F*) efficiency_BCDEF_tight_iso->Clone("efficiency_MC_tight_iso");

  for(int iBinX = 0; iBinX < efficiency_BCDEF_tight_iso->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < efficiency_BCDEF_tight_iso->GetNbinsY(); iBinY++){
      efficiency_MC_tight_iso->SetBinContent(iBinX+1,iBinY+1,
					    (efficiency_BCDEF_tight_iso->GetBinContent(iBinX+1,iBinY+1)*(lumi_BCDEF/lumi_total)+efficiency_GH_tight_iso->GetBinContent(iBinX+1,iBinY+1)*(lumi_GH/lumi_total))/(lumi_BCDEF/lumi_total+lumi_GH/lumi_total));
      
      cout<<"Tight iso : BCDEF "<<efficiency_BCDEF_tight_iso->GetBinContent(iBinX+1,iBinY+1)<<" GH "<<efficiency_GH_tight_iso->GetBinContent(iBinX+1,iBinY+1)<<" mean "<<efficiency_MC_tight_iso->GetBinContent(iBinX+1,iBinY+1)<<endl;
    }
  }

  TFile* outputFile = new TFile("muon_efficiency_MC.root","RECREATE");
  outputFile->cd();

  efficiency_MC_loose_id->Write("efficiency_MC_loose_id");
  efficiency_MC_loose_iso->Write("efficiency_MC_loose_iso");
  efficiency_MC_tight_id->Write("efficiency_MC_tight_id");
  efficiency_MC_tight_iso->Write("efficiency_MC_tight_iso");

  outputFile->Close();
}
