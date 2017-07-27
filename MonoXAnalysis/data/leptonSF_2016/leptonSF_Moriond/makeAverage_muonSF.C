void makeAverage_muonSF(){

  TFile* inputBCDEF_id  = TFile::Open("~/EfficienciesAndSF_BCDEF.root");
  TFile* inputGH_id     = TFile::Open("~/EfficienciesAndSF_GH.root");
  TFile* inputBCDEF_iso = TFile::Open("~/EfficienciesAndSF_BCDEF-2.root");
  TFile* inputGH_iso    = TFile::Open("~/EfficienciesAndSF_GH-2.root");

  TH2F* scalefactor_BCDEF_loose_id = (TH2F*) inputBCDEF_id->Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
  TH2F* scalefactor_GH_loose_id = (TH2F*) inputGH_id->Get("MC_NUM_LooseID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
  TH2F* scalefactor_BCDEF_tight_id = (TH2F*) inputBCDEF_id->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");
  TH2F* scalefactor_GH_tight_id = (TH2F*) inputGH_id->Get("MC_NUM_TightID_DEN_genTracks_PAR_pt_eta/abseta_pt_ratio");

  TH2F* scalefactor_BCDEF_loose_iso = (TH2F*) inputBCDEF_iso->Get("LooseISO_LooseID_pt_eta/abseta_pt_ratio");
  TH2F* scalefactor_GH_loose_iso = (TH2F*) inputGH_iso->Get("LooseISO_LooseID_pt_eta/abseta_pt_ratio");
  TH2F* scalefactor_BCDEF_tight_iso = (TH2F*) inputBCDEF_iso->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
  TH2F* scalefactor_GH_tight_iso = (TH2F*) inputGH_iso->Get("TightISO_TightID_pt_eta/abseta_pt_ratio");
  
  float lumi_BCDEF = 19.72;
  float lumi_GH    = 16.14;
  float lumi_total = 35.9;

  if(scalefactor_BCDEF_loose_id->GetNbinsX() != scalefactor_GH_loose_id->GetNbinsX())
    cout<<"Problem i.e. different x-axis binning for loose id "<<endl;
  if(scalefactor_BCDEF_loose_id->GetNbinsY() != scalefactor_GH_loose_id->GetNbinsY())
    cout<<"Problem i.e. different y-axis binning for loose id "<<endl;

  TH2F* scalefactor_loose_id = (TH2F*) scalefactor_BCDEF_loose_id->Clone("scalefactors_MuonLooseId_Muon");

  for(int iBinX = 0; iBinX < scalefactor_BCDEF_loose_id->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < scalefactor_BCDEF_loose_id->GetNbinsY(); iBinY++){
      scalefactor_loose_id->SetBinContent(iBinX+1,iBinY+1,
					    (scalefactor_BCDEF_loose_id->GetBinContent(iBinX+1,iBinY+1)*(lumi_BCDEF/lumi_total)+scalefactor_GH_loose_id->GetBinContent(iBinX+1,iBinY+1)*(lumi_GH/lumi_total))/(lumi_BCDEF/lumi_total+lumi_GH/lumi_total));            
      cout<<"Loose id : BCDEF "<<scalefactor_BCDEF_loose_id->GetBinContent(iBinX+1,iBinY+1)<<" GH "<<scalefactor_GH_loose_id->GetBinContent(iBinX+1,iBinY+1)<<" mean "<<scalefactor_loose_id->GetBinContent(iBinX+1,iBinY+1)<<endl;
    }
  }

  if(scalefactor_BCDEF_tight_id->GetNbinsX() != scalefactor_GH_tight_id->GetNbinsX())
    cout<<"Problem i.e. different x-axis binning for tight id "<<endl;
  if(scalefactor_BCDEF_tight_id->GetNbinsY() != scalefactor_GH_tight_id->GetNbinsY())
    cout<<"Problem i.e. different y-axis binning for tight id "<<endl;

  TH2F* scalefactor_tight_id = (TH2F*) scalefactor_BCDEF_tight_id->Clone("scalefactors_TightId_Muon");

  for(int iBinX = 0; iBinX < scalefactor_BCDEF_tight_id->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < scalefactor_BCDEF_tight_id->GetNbinsY(); iBinY++){
      scalefactor_tight_id->SetBinContent(iBinX+1,iBinY+1,
					    (scalefactor_BCDEF_tight_id->GetBinContent(iBinX+1,iBinY+1)*(lumi_BCDEF/lumi_total)+scalefactor_GH_tight_id->GetBinContent(iBinX+1,iBinY+1)*(lumi_GH/lumi_total))/(lumi_BCDEF/lumi_total+lumi_GH/lumi_total));
      
      cout<<"Tight id : BCDEF "<<scalefactor_BCDEF_tight_id->GetBinContent(iBinX+1,iBinY+1)<<" GH "<<scalefactor_GH_tight_id->GetBinContent(iBinX+1,iBinY+1)<<" mean "<<scalefactor_tight_id->GetBinContent(iBinX+1,iBinY+1)<<endl;
    }
  }
  ///////////

  if(scalefactor_BCDEF_loose_iso->GetNbinsX() != scalefactor_GH_loose_iso->GetNbinsX())
    cout<<"Problem i.e. different x-axis binning for loose iso "<<endl;
  if(scalefactor_BCDEF_loose_iso->GetNbinsY() != scalefactor_GH_loose_iso->GetNbinsY())
    cout<<"Problem i.e. different y-axis binning for loose iso "<<endl;

  TH2F* scalefactor_loose_iso = (TH2F*) scalefactor_BCDEF_loose_iso->Clone("scalefactors_Iso_MuonLooseId");

  for(int iBinX = 0; iBinX < scalefactor_BCDEF_loose_iso->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < scalefactor_BCDEF_loose_iso->GetNbinsY(); iBinY++){
      scalefactor_loose_iso->SetBinContent(iBinX+1,iBinY+1,
					    (scalefactor_BCDEF_loose_iso->GetBinContent(iBinX+1,iBinY+1)*(lumi_BCDEF/lumi_total)+scalefactor_GH_loose_iso->GetBinContent(iBinX+1,iBinY+1)*(lumi_GH/lumi_total))/(lumi_BCDEF/lumi_total+lumi_GH/lumi_total));
    
      cout<<"Loose iso : BCDEF "<<scalefactor_BCDEF_loose_iso->GetBinContent(iBinX+1,iBinY+1)<<" GH "<<scalefactor_GH_loose_iso->GetBinContent(iBinX+1,iBinY+1)<<" mean "<<scalefactor_loose_iso->GetBinContent(iBinX+1,iBinY+1)<<endl;
    }
  }

  if(scalefactor_BCDEF_tight_iso->GetNbinsX() != scalefactor_GH_tight_iso->GetNbinsX())
    cout<<"Problem i.e. different x-axis binning for tight iso "<<endl;
  if(scalefactor_BCDEF_tight_iso->GetNbinsY() != scalefactor_GH_tight_iso->GetNbinsY())
    cout<<"Problem i.e. different y-axis binning for tight iso "<<endl;

  TH2F* scalefactor_tight_iso = (TH2F*) scalefactor_BCDEF_tight_iso->Clone("scalefactors_Iso_MuonTightId");

  for(int iBinX = 0; iBinX < scalefactor_BCDEF_tight_iso->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < scalefactor_BCDEF_tight_iso->GetNbinsY(); iBinY++){
      scalefactor_tight_iso->SetBinContent(iBinX+1,iBinY+1,
					    (scalefactor_BCDEF_tight_iso->GetBinContent(iBinX+1,iBinY+1)*(lumi_BCDEF/lumi_total)+scalefactor_GH_tight_iso->GetBinContent(iBinX+1,iBinY+1)*(lumi_GH/lumi_total))/(lumi_BCDEF/lumi_total+lumi_GH/lumi_total));
      
      cout<<"Tight iso : BCDEF "<<scalefactor_BCDEF_tight_iso->GetBinContent(iBinX+1,iBinY+1)<<" GH "<<scalefactor_GH_tight_iso->GetBinContent(iBinX+1,iBinY+1)<<" mean "<<scalefactor_tight_iso->GetBinContent(iBinX+1,iBinY+1)<<endl;
    }
  }

  TFile* outputFile = new TFile("muon_scalefactors.root","RECREATE");
  outputFile->cd();

  scalefactor_loose_id->Write();
  scalefactor_loose_iso->Write();
  scalefactor_tight_id->Write();
  scalefactor_tight_iso->Write();

  outputFile->Close();
}
