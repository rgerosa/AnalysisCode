void makeWJetsTauEfficiency(string outputFileName){

  string inputDIR_HT100to200 = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/ProductionMC_21_07_2017/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/";
  string inputDIR_HT200to400 = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/ProductionMC_21_07_2017/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/";
  string inputDIR_HT400to600 = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/ProductionMC_21_07_2017/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/";
  string inputDIR_HT600to800 = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/ProductionMC_21_07_2017/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/";
  string inputDIR_HT800to1200 = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/ProductionMC_21_07_2017/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/";
  string inputDIR_HT1200to2500 = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/ProductionMC_21_07_2017/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/";
  string inputDIR_HT2500toInf = "/eos/cms/store/group/phys_exotica/monojet/rgerosa/ProductionMC_21_07_2017/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/";

  system(("find "+inputDIR_HT100to200+" -name \"*root\" | grep -v failed > file_ht100to200.temp").c_str());
  system(("find "+inputDIR_HT200to400+" -name \"*root\" | grep -v failed > file_ht200to400.temp").c_str());
  system(("find "+inputDIR_HT400to600+" -name \"*root\" | grep -v failed > file_ht400to600.temp").c_str());
  system(("find "+inputDIR_HT600to800+" -name \"*root\" | grep -v failed > file_ht600to800.temp").c_str());
  system(("find "+inputDIR_HT800to1200+" -name \"*root\" | grep -v failed > file_ht800to1200.temp").c_str());
  system(("find "+inputDIR_HT1200to2500+" -name \"*root\" | grep -v failed > file_ht1200to2500.temp").c_str());
  system(("find "+inputDIR_HT2500toInf+" -name \"*root\" | grep -v failed > file_ht2500toInf.temp").c_str());

  TH1F* tauDenom_HT100to200 =  0;
  TH1F* tauNum_HT100to200 =  0;
  TH1F* tauDenom_HT200to400 = 0;
  TH1F* tauNum_HT200to400 = 0;
  TH1F* tauDenom_HT400to600 = 0;
  TH1F* tauNum_HT400to600 = 0;
  TH1F* tauDenom_HT600to800 = 0;
  TH1F* tauNum_HT600to800 = 0;
  TH1F* tauDenom_HT800to1200 = 0;
  TH1F* tauNum_HT800to1200 = 0;
  TH1F* tauDenom_HT1200to2500 = 0;
  TH1F* tauNum_HT1200to2500 = 0;
  TH1F* tauDenom_HT2500toInf = 0;
  TH1F* tauNum_HT2500toInf = 0;


  cout<<"Fill histogram HT 100-200"<<endl;  
  ifstream infile;
  string line;
  infile.open("file_ht100to200.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "" or not TString(line).Contains(".root")) continue;
      TFile* temp = TFile::Open(line.c_str(),"READ");
      if(tauDenom_HT100to200 == 0){
	tauDenom_HT100to200 = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"))->Clone("tauDenom_HT100to200");
	tauDenom_HT100to200->SetDirectory(0);
      }
      else
	tauDenom_HT100to200->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"));
      if(tauNum_HT100to200 == 0){
	tauNum_HT100to200 = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"))->Clone("tauNum_HT100to200");
	tauNum_HT100to200->SetDirectory(0);
      }
      else
	tauNum_HT100to200->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"));
      temp->Close();
    }
  }
  infile.close();

  cout<<"Fill histogram HT 200-400"<<endl;
  infile.open("file_ht200to400.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "" or not TString(line).Contains(".root")) continue;
      TFile* temp = TFile::Open(line.c_str(),"READ");
      if(tauDenom_HT200to400 == 0){
	tauDenom_HT200to400 = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"))->Clone("tauDenom_HT200to400");
	tauDenom_HT200to400->SetDirectory(0);
      }
      else
	tauDenom_HT200to400->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"));
      if(tauNum_HT200to400 == 0){
	tauNum_HT200to400 = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"))->Clone("tauNum_HT200to400");
	tauNum_HT200to400->SetDirectory(0);
      }
      else
	tauNum_HT200to400->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"));
      temp->Close();
    }
  }
  infile.close();
  
  cout<<"Fill histogram HT 400-600"<<endl;
  infile.open("file_ht400to600.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "" or not TString(line).Contains(".root")) continue;
      TFile* temp = TFile::Open(line.c_str(),"READ");
      if(tauDenom_HT400to600 == 0){
	tauDenom_HT400to600 = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"))->Clone("tauDenom_HT400to600");
	tauDenom_HT400to600->SetDirectory(0);
      }
      else
	tauDenom_HT400to600->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"));
      if(tauNum_HT400to600 == 0){
	tauNum_HT400to600 = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"))->Clone("tauNum_HT400to600");
	tauNum_HT400to600->SetDirectory(0);
      }
      else
	tauNum_HT400to600->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"));
      temp->Close();
    }
  }
  infile.close();

  cout<<"Fill histogram HT 600-800"<<endl;
  infile.open("file_ht600to800.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "" or not TString(line).Contains(".root")) continue;
      TFile* temp = TFile::Open(line.c_str(),"READ");
      if(tauDenom_HT600to800 == 0){
	tauDenom_HT600to800 = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"))->Clone("tauDenom_HT600to800");
	tauDenom_HT600to800->SetDirectory(0);
      }
      else
	tauDenom_HT600to800->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"));
      if(tauNum_HT600to800 == 0){
	tauNum_HT600to800 = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"))->Clone("tauNum_HT600to800");
	tauNum_HT600to800->SetDirectory(0);
      }
      else
	tauNum_HT600to800->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"));
      temp->Close();
    }
  }
  infile.close();

  cout<<"Fill histogram HT 800-1200"<<endl;
  infile.open("file_ht800to1200.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "" or not TString(line).Contains(".root")) continue;
      TFile* temp = TFile::Open(line.c_str(),"READ");
      if(tauDenom_HT800to1200 == 0){
	tauDenom_HT800to1200 = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"))->Clone("tauDenom_HT800to1200");
	tauDenom_HT800to1200->SetDirectory(0);
      }
      else
	tauDenom_HT800to1200->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"));
      if(tauNum_HT800to1200 == 0){
	tauNum_HT800to1200 = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"))->Clone("tauNum_HT800to1200");
	tauNum_HT800to1200->SetDirectory(0);
      }
      else
	tauNum_HT800to1200->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"));
      temp->Close();
    }
  }
  infile.close();

  cout<<"Fill histogram HT 1200-2500"<<endl;
  infile.open("file_ht1200to2500.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "" or not TString(line).Contains(".root")) continue;
      TFile* temp = TFile::Open(line.c_str(),"READ");
      if(tauDenom_HT1200to2500 == 0){
	tauDenom_HT1200to2500 = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"))->Clone("tauDenom_HT1200to2500");
	tauDenom_HT1200to2500->SetDirectory(0);
      }
      else
	tauDenom_HT1200to2500->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"));
      if(tauNum_HT1200to2500 == 0){
	tauNum_HT1200to2500 = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"))->Clone("tauNum_HT1200to2500");
	tauNum_HT1200to2500->SetDirectory(0);
      }
      else
	tauNum_HT1200to2500->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"));
      temp->Close();
    }
  }
  infile.close();

  cout<<"Fill histogram HT 2500-Inf"<<endl;
  infile.open("file_ht2500toInf.temp",ifstream::in);
  if(infile.is_open()){
    while(!infile.eof()){
      getline(infile,line);
      if(line == "" or not TString(line).Contains(".root")) continue;
      TFile* temp = TFile::Open(line.c_str(),"READ");
      if(tauDenom_HT2500toInf == 0){
	tauDenom_HT2500toInf = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"))->Clone("tauDenom_HT2500toInf");
	tauDenom_HT2500toInf->SetDirectory(0);
      }
      else
	tauDenom_HT2500toInf->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Denom_tau"));
      if(tauNum_HT2500toInf == 0){
	tauNum_HT2500toInf = (TH1F*)((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"))->Clone("tauNum_HT2500toInf");
	tauNum_HT2500toInf->SetDirectory(0);
      }
      else
	tauNum_HT2500toInf->Add((TH1F*) temp->Get("taueff/eff_VLooseTauOldDM_Num_tau"));
      temp->Close();
    }
  }
  infile.close();


  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();

  TH1F* efficiency_tau_HT100to200 = (TH1F*) tauNum_HT100to200->Clone("efficiency_tau_HT100to200");
  efficiency_tau_HT100to200->Divide(tauDenom_HT100to200);
  TH1F* efficiency_tau_HT200to400 = (TH1F*) tauNum_HT200to400->Clone("efficiency_tau_HT200to400");
  efficiency_tau_HT200to400->Divide(tauDenom_HT200to400);
  TH1F* efficiency_tau_HT400to600 = (TH1F*) tauNum_HT400to600->Clone("efficiency_tau_HT400to600");
  efficiency_tau_HT400to600->Divide(tauDenom_HT400to600);
  TH1F* efficiency_tau_HT600to800 = (TH1F*) tauNum_HT600to800->Clone("efficiency_tau_HT600to800");
  efficiency_tau_HT600to800->Divide(tauDenom_HT600to800);
  TH1F* efficiency_tau_HT800to1200 = (TH1F*) tauNum_HT800to1200->Clone("efficiency_tau_HT800to1200");
  efficiency_tau_HT800to1200->Divide(tauDenom_HT800to1200);
  TH1F* efficiency_tau_HT1200to2500 = (TH1F*) tauNum_HT1200to2500->Clone("efficiency_tau_HT1200to2500");
  efficiency_tau_HT1200to2500->Divide(tauDenom_HT1200to2500);
  TH1F* efficiency_tau_HT2500toInf = (TH1F*) tauNum_HT2500toInf->Clone("efficiency_tau_HT2500ToInf");
  efficiency_tau_HT2500toInf->Divide(tauDenom_HT2500toInf);

  // take the sum of weights:
  double sumwgt_ht100to200, sumwgt_ht200to400, sumwgt_ht400to600, sumwgt_ht600to800, sumwgt_ht800to1200, sumwgt_ht1200to2500, sumwgt_ht2500toInf;
  float  xsec_ht100to200, xsec_ht200to400, xsec_ht400to600, xsec_ht600to800, xsec_ht800to1200, xsec_ht1200to2500, xsec_ht2500toInf;

  TFile* wjet_sig_ht100to200 = TFile::Open("/eos/cms/store/group/phys_exotica/monojet/rgerosa/SkimmedProductionMC_21_07_2017/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/sigfilter/sig_WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","READ");
  TFile* wjet_sig_ht200to400 = TFile::Open("/eos/cms/store/group/phys_exotica/monojet/rgerosa/SkimmedProductionMC_21_07_2017/WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/sigfilter/sig_WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","READ");
  TFile* wjet_sig_ht400to600 = TFile::Open("/eos/cms/store/group/phys_exotica/monojet/rgerosa/SkimmedProductionMC_21_07_2017/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/sigfilter/sig_WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","READ");
  TFile* wjet_sig_ht600to800 = TFile::Open("/eos/cms/store/group/phys_exotica/monojet/rgerosa/SkimmedProductionMC_21_07_2017/WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/sigfilter/sig_WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","READ");
  TFile* wjet_sig_ht800to1200 = TFile::Open("/eos/cms/store/group/phys_exotica/monojet/rgerosa/SkimmedProductionMC_21_07_2017/WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/sigfilter/sig_WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","READ");
  TFile* wjet_sig_ht1200to2500 = TFile::Open("/eos/cms/store/group/phys_exotica/monojet/rgerosa/SkimmedProductionMC_21_07_2017/WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/sigfilter/sig_WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","READ");
  TFile* wjet_sig_ht2500toInf = TFile::Open("/eos/cms/store/group/phys_exotica/monojet/rgerosa/SkimmedProductionMC_21_07_2017/WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/sigfilter/sig_WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root","READ");

  TTree* tree_sig_ht100to200 = (TTree*) wjet_sig_ht100to200->Get("tree/tree");
  tree_sig_ht100to200->SetBranchAddress("wgtsum",&sumwgt_ht100to200);
  tree_sig_ht100to200->SetBranchAddress("xsec",&xsec_ht100to200);
  tree_sig_ht100to200->GetEntry(0);

  TTree* tree_sig_ht200to400 = (TTree*) wjet_sig_ht200to400->Get("tree/tree");
  tree_sig_ht200to400->SetBranchAddress("wgtsum",&sumwgt_ht200to400);
  tree_sig_ht200to400->SetBranchAddress("xsec",&xsec_ht200to400);
  tree_sig_ht200to400->GetEntry(0);

  TTree* tree_sig_ht400to600 = (TTree*) wjet_sig_ht400to600->Get("tree/tree");
  tree_sig_ht400to600->SetBranchAddress("wgtsum",&sumwgt_ht400to600);
  tree_sig_ht400to600->SetBranchAddress("xsec",&xsec_ht400to600);
  tree_sig_ht400to600->GetEntry(0);

  TTree* tree_sig_ht600to800 = (TTree*) wjet_sig_ht600to800->Get("tree/tree");
  tree_sig_ht600to800->SetBranchAddress("wgtsum",&sumwgt_ht600to800);
  tree_sig_ht600to800->SetBranchAddress("xsec",&xsec_ht600to800);
  tree_sig_ht600to800->GetEntry(0);

  TTree* tree_sig_ht800to1200 = (TTree*) wjet_sig_ht800to1200->Get("tree/tree");
  tree_sig_ht800to1200->SetBranchAddress("wgtsum",&sumwgt_ht800to1200);
  tree_sig_ht800to1200->SetBranchAddress("xsec",&xsec_ht800to1200);
  tree_sig_ht800to1200->GetEntry(0);

  TTree* tree_sig_ht1200to2500 = (TTree*) wjet_sig_ht1200to2500->Get("tree/tree");
  tree_sig_ht1200to2500->SetBranchAddress("wgtsum",&sumwgt_ht1200to2500);
  tree_sig_ht1200to2500->SetBranchAddress("xsec",&xsec_ht1200to2500);
  tree_sig_ht1200to2500->GetEntry(0);

  TTree* tree_sig_ht2500toInf = (TTree*) wjet_sig_ht2500toInf->Get("tree/tree");
  tree_sig_ht2500toInf->SetBranchAddress("wgtsum",&sumwgt_ht2500toInf);
  tree_sig_ht2500toInf->SetBranchAddress("xsec",&xsec_ht2500toInf);
  tree_sig_ht2500toInf->GetEntry(0);

  outputFile->cd();

  efficiency_tau_HT100to200->Write();
  efficiency_tau_HT200to400->Write();
  efficiency_tau_HT400to600->Write();
  efficiency_tau_HT600to800->Write();
  efficiency_tau_HT800to1200->Write();
  efficiency_tau_HT1200to2500->Write();
  efficiency_tau_HT2500toInf->Write();

  TH2F* efficiency_tau_MC = (TH2F*) efficiency_tau_HT100to200->Clone("efficiency_tau_MC");

  for(int iBinX = 0; iBinX < efficiency_tau_HT100to200->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < efficiency_tau_HT100to200->GetNbinsY(); iBinY++){
      efficiency_tau_MC->SetBinContent(iBinX+1,iBinY+1,(
							efficiency_tau_HT100to200->GetBinContent(iBinX+1,iBinY+1)*xsec_ht100to200/sumwgt_ht100to200+
							efficiency_tau_HT200to400->GetBinContent(iBinX+1,iBinY+1)*xsec_ht200to400/sumwgt_ht200to400+
							efficiency_tau_HT400to600->GetBinContent(iBinX+1,iBinY+1)*xsec_ht400to600/sumwgt_ht400to600+
							efficiency_tau_HT600to800->GetBinContent(iBinX+1,iBinY+1)*xsec_ht600to800/sumwgt_ht600to800+
							efficiency_tau_HT800to1200->GetBinContent(iBinX+1,iBinY+1)*xsec_ht800to1200/sumwgt_ht800to1200+
							efficiency_tau_HT1200to2500->GetBinContent(iBinX+1,iBinY+1)*xsec_ht1200to2500/sumwgt_ht1200to2500+
							efficiency_tau_HT2500toInf->GetBinContent(iBinX+1,iBinY+1)*xsec_ht2500toInf/sumwgt_ht2500toInf)/
				       (xsec_ht100to200/sumwgt_ht100to200+xsec_ht200to400/sumwgt_ht200to400+xsec_ht400to600/sumwgt_ht400to600+xsec_ht600to800/sumwgt_ht600to800+xsec_ht800to1200/sumwgt_ht800to1200+xsec_ht1200to2500/sumwgt_ht1200to2500+xsec_ht2500toInf/sumwgt_ht2500toInf));
    }
  }

  efficiency_tau_MC->Write();
  system("rm file_ht*.temp");

  outputFile->Close();

}
