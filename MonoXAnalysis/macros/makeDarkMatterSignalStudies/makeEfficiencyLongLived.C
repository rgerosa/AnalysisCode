#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void plotEfficiency(TCanvas* canvas,vector<TH1F*> eff, const vector<string> & label_legend, const string & label, const string & outputDIR, const float & lumi, const string & postfix){
  
  canvas->cd();
  int icolor = 1;
  TLegend leg (0.18,0.25,0.4,0.5);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry((TObject*)(0),label.c_str(),"");
  for(int ieff = 0; ieff < eff.size(); ieff++){
    eff.at(ieff)->GetXaxis()->SetTitleSize(0);
    eff.at(ieff)->GetYaxis()->SetTitle("Efficiency");
    if(icolor == 3) icolor++;
    if(icolor == 5) icolor++;
    eff.at(ieff)->SetLineColor(icolor);
    eff.at(ieff)->SetMarkerColor(icolor);
    eff.at(ieff)->SetMarkerStyle(20);
    eff.at(ieff)->SetMarkerSize(1);
    if(ieff == 0)
      eff.at(ieff)->Draw("EPL");
    else
      eff.at(ieff)->Draw("EPLsame");
    
    leg.AddEntry(eff.at(ieff),label_legend.at(ieff).c_str(),"EP");
    icolor++;
  }
  
  leg.Draw("same");
  eff.at(0)->GetYaxis()->SetRangeUser(0.001,10);
  eff.at(0)->GetXaxis()->SetNdivisions(-414);
  canvas->SetLogy();

  CMS_lumi(canvas,Form("%.1f",lumi));
  
  canvas->SaveAs((outputDIR+"/eff_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/eff_"+postfix+".pdf").c_str(),"pdf");

}

void makeEfficiencyLongLived(string inputFileDIR, 
			     string outputDIR, 
			     string outputFileName, 
			     string observable, 
			     Category category,
			     float lumi = 35.9){

  if(category != Category::monojet and category != Category::monoV){
    cerr<<"If category is different from monojet and monoV the analysis cannot be run --> check"<<endl;
    return;
  }
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  initializeBinning();

  vector<TFile*> inputFiles;
  vector<TTree*> inputTrees;
  vector<TTree*> inputGenTrees;
  
  system(("ls "+inputFileDIR+" | grep root > file.list").c_str());
  ifstream infile ("file.list");
  if(infile.is_open()){
    std::string line;
    while(std::getline(infile,line)){
      if(line == "" or line == "\n") continue;
      inputFiles.push_back(TFile::Open((inputFileDIR+"/"+line).c_str(),"READ"));
      inputTrees.push_back((TTree*) inputFiles.back()->Get("tree/tree"));
      inputGenTrees.push_back((TTree*) inputFiles.back()->Get("gentree/gentree"));
    }
  }
  infile.close();
  system("rm file.list");
  
  vector<double> sumOfWeights ;
  cout<<"Calculate sum of weights "<<endl;
  int ifile = 0;
  for(auto tree : inputGenTrees){
    cout<<"Weights for "<<inputFiles.at(ifile)->GetName()<<endl;
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");
    sumOfWeights.push_back(0);
    while(reader.Next()){
      sumOfWeights.back() += *wgt;
    }
    ifile++;
  }

  cout<<"Calculate Efficiency "<<endl;
  vector<TH1F*> event_inclusive;
  vector<TH1F*> event_leptonVeto;
  vector<TH1F*> event_tauVeto;
  vector<TH1F*> event_bVeto;
  vector<TH1F*> event_photonVeto;
  vector<TH1F*> event_caloMet;
  vector<TH1F*> event_jetCut;
  vector<TH1F*> event_jetCleaning;
  vector<TH1F*> event_jetMetDphi;
  vector<TH1F*> event_metCut;
  vector<TH1F*> event_monoVVeto;
  vector<TH1F*> event_expWeight;

  vector<double> bins = selectBinning(observable,category);
  vector<TH1F*> final_templates;
  vector<TH1F*> mediator_pt;
  vector<TH1F*> ratio_x2_pt;
  vector<TH1F*> response_mediator_pt;
  vector<TH1F*> ratio_x1_pt;
  vector<TH1F*> ratio_x1x2_pt;
  vector<TH1F*> response_x1_pt;


  // load all the weights for the SR selection                                                                                                                                                       
  TFile* pufile = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/npvWeight/purwt_36.40_summer16.root");
  TH1* puhist = (TH1*) pufile->Get("puhist");
  //////
  TFile* triggerfile_MET = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/triggerSF_2016/trigger_MORIOND/Monojet/metTriggerEfficiency_recoil_monojet.root");
  ///// ---                                                                                                                                                                                           
  TEfficiency* triggermet = (TEfficiency*) triggerfile_MET->Get("efficiency");
  TGraphAsymmErrors* triggermet_graph = triggermet->CreateGraph();

  // start loop on trees
  ifile = 0;
  cout<<"Start loop on trees "<<endl;
  for(auto tree: inputTrees){

    cout<<"Event loop: "<<inputFiles.at(ifile)->GetName()<<endl;

    TTreeReader reader (tree);

    TTreeReaderValue<unsigned int> nvtx (reader,"nvtx");
    TTreeReaderValue<float>  xsec       (reader,"xsec");
    TTreeReaderValue<float>  wgt        (reader,"wgt");

    TTreeReaderValue<UChar_t> hltm90    (reader,"hltmet90");
    TTreeReaderValue<UChar_t> hltm100   (reader,"hltmet100");
    TTreeReaderValue<UChar_t> hltm110   (reader,"hltmet110");
    TTreeReaderValue<UChar_t> hltm120   (reader,"hltmet120");
    TTreeReaderValue<UChar_t> hltmwm90  (reader,"hltmetwithmu90");
    TTreeReaderValue<UChar_t> hltmwm120 (reader,"hltmetwithmu120");
    TTreeReaderValue<UChar_t> hltmwm170 (reader,"hltmetwithmu170");
    TTreeReaderValue<UChar_t> hltmwm300 (reader,"hltmetwithmu300");

    TTreeReaderValue<UChar_t> fhbhe  (reader,"flaghbhenoise");
    TTreeReaderValue<UChar_t> fhbiso (reader,"flaghbheiso");
    TTreeReaderValue<UChar_t> fcsc  (reader,"flagglobaltighthalo");
    TTreeReaderValue<UChar_t> feeb   (reader,"flageebadsc");
    TTreeReaderValue<UChar_t> fetp   (reader,"flagecaltp");
    TTreeReaderValue<UChar_t> fvtx   (reader,"flaggoodvertices");
    TTreeReaderValue<UChar_t> fbadmu (reader,"flagbadpfmu");
    TTreeReaderValue<UChar_t> fbadch (reader,"flagbadchpf");
    
    TTreeReaderValue<unsigned int> njets      (reader,"njets");
    TTreeReaderValue<unsigned int> nphotons   (reader,"nphotons");
    TTreeReaderValue<unsigned int> nelectrons (reader,"nelectrons");
    TTreeReaderValue<unsigned int> ntaus      (reader,"ntausold");
    TTreeReaderValue<unsigned int> nmuons     (reader,"nmuons");
    TTreeReaderValue<unsigned int> nbjets     (reader,"nbjetslowpt");
    
    TTreeReaderValue<vector<float> > jetpt   (reader,"combinejetpt");
    TTreeReaderValue<vector<float> > jeteta  (reader,"combinejeteta");
    TTreeReaderValue<vector<float> > jetphi  (reader,"combinejetphi");
    TTreeReaderValue<vector<float> > jetbtag (reader,"combinejetbtag");
    TTreeReaderValue<vector<float> > jetm    (reader,"combinejetm");
    TTreeReaderValue<vector<float> > chfrac  (reader,"combinejetCHfrac");
    TTreeReaderValue<vector<float> > nhfrac  (reader,"combinejetNHfrac");

    TTreeReaderValue<vector<float> > boostedJetpt    (reader,"boostedJetpt");
    TTreeReaderValue<vector<float> > boostedJeteta   (reader,"boostedJeteta");
    TTreeReaderValue<vector<float> > boostedJetphi   (reader,"boostedJetphi");
    TTreeReaderValue<vector<float> > boostedJetm     (reader,"boostedJetm");
    TTreeReaderValue<vector<float> > prunedJetm      (reader,"prunedJetm");
    TTreeReaderValue<vector<float> > boostedJettau2  (reader,"boostedJettau2");
    TTreeReaderValue<vector<float> > boostedJettau1  (reader,"boostedJettau1");

    TTreeReaderValue<float> mmet       (reader,"t1mumet");
    TTreeReaderValue<float> jmmdphi    (reader,"incjetmumetdphimin4");
    TTreeReaderValue<float> metcalo    (reader,"calomet");

    TTreeReaderValue<float> dmpt       (reader,"dmpt");
    TTreeReaderValue<float> dmeta      (reader,"dmeta");
    TTreeReaderValue<float> dmphi      (reader,"dmphi");
    TTreeReaderValue<float> dmmass     (reader,"dmmass");

    TTreeReaderValue<float> dmX1mass (reader,"dmX1mass");
    TTreeReaderValue<float> dmX1eta (reader,"dmX1eta");
    TTreeReaderValue<float> dmX1phi (reader,"dmX1phi");
    TTreeReaderValue<float> dmX1pt (reader,"dmX1pt");

    TTreeReaderValue<float> dmX2mass (reader,"dmX2mass");
    TTreeReaderValue<float> dmX2eta (reader,"dmX2eta");
    TTreeReaderValue<float> dmX2phi (reader,"dmX2phi");
    TTreeReaderValue<float> dmX2pt (reader,"dmX2pt");

    TTreeReaderValue<float> dmX1mass_ll (reader,"dmX1mass_ll");
    TTreeReaderValue<float> dmX1eta_ll (reader,"dmX1eta_ll");
    TTreeReaderValue<float> dmX1phi_ll (reader,"dmX1phi_ll");
    TTreeReaderValue<float> dmX1pt_ll (reader,"dmX1pt_ll");

    TTreeReaderValue<float> dmX2mass_ll (reader,"dmX2mass_ll");
    TTreeReaderValue<float> dmX2eta_ll (reader,"dmX2eta_ll");
    TTreeReaderValue<float> dmX2phi_ll (reader,"dmX2phi_ll");
    TTreeReaderValue<float> dmX2pt_ll (reader,"dmX2pt_ll");


    stringstream fileName(TString(inputFiles.at(ifile)->GetName()).ReplaceAll(".root","").Data());
    std::string segment;
    std::vector<std::string> seglist;
    std::vector<std::string> seglist_v2;

    while(std::getline(fileName, segment,'/')){
      seglist.push_back(segment);
    }

    stringstream fileName2(seglist.back());
    while(std::getline(fileName2,segment,'_')){
      seglist_v2.push_back(segment);
    }
    
    std::vector<std::string> Mmed;
    std::vector<std::string> Mchi2;
    std::vector<std::string> ctau;
    
    for(auto name : seglist_v2){
      TString name_tmp (name);
      if(name_tmp.Contains("Mphi")){
	stringstream fileName_tmp(name_tmp.Data());
	while(std::getline(fileName_tmp, segment,'-')){
	  Mmed.push_back(segment);
	}
      }
      else if(name_tmp.Contains("Mchi2")){
	stringstream fileName_tmp(name_tmp.Data());
	while(std::getline(fileName_tmp, segment,'-')){
	  Mchi2.push_back(segment);
	}
      }
      else if(name_tmp.Contains("ctau")){
	stringstream fileName_tmp(name_tmp.Data());
        while(std::getline(fileName_tmp, segment,'-')){
          ctau.push_back(segment);
        }
      }	
    }

    event_inclusive.push_back(new TH1F(Form("event_inclusive_mchi2_%s_ctau_%s",Mchi2.back().c_str(),ctau.back().c_str()),"",1,0,10000)); //we are interested in a inclusive efficiency
    event_leptonVeto.push_back(new TH1F(Form("event_leptonVeto_mchi2_%s_ctau_%s",Mchi2.back().c_str(),ctau.back().c_str()),"",1,0,10000));
    event_tauVeto.push_back(new TH1F(Form("event_tauVeto_mchi2_%s_ctau_%s",Mchi2.back().c_str(),ctau.back().c_str()),"",1,0,10000));
    event_bVeto.push_back(new TH1F(Form("event_bVeto_mchi2_%s_ctau_%s",Mchi2.back().c_str(),ctau.back().c_str()),"",1,0,10000));
    event_photonVeto.push_back(new TH1F(Form("event_photonVeto_mchi2_%s_ctau_%s",Mchi2.back().c_str(),ctau.back().c_str()),"",1,0,10000));
    event_caloMet.push_back(new TH1F(Form("event_caloMet_mchi2_%s_ctau_%s",Mchi2.back().c_str(),ctau.back().c_str()),"",1,0,10000));
    event_jetCut.push_back(new TH1F(Form("event_jetCut_mchi2_%s_ctau_%s",Mchi2.back().c_str(),ctau.back().c_str()),"",1,0,10000));
    event_jetCleaning.push_back(new TH1F(Form("event_jetCleaning_mchi2_%s_ctau_%s",Mchi2.back().c_str(),ctau.back().c_str()),"",1,0,10000));
    event_jetMetDphi.push_back(new TH1F(Form("event_jetMetDphi_mchi2_%s_ctau_%s",Mchi2.back().c_str(),ctau.back().c_str()),"",1,0,10000));
    event_metCut.push_back(new TH1F(Form("event_metCut_mchi2_%s_ctau_%s",Mchi2.back().c_str(),ctau.back().c_str()),"",1,0,10000));
    event_monoVVeto.push_back(new TH1F(Form("event_monoVVeto_mchi2_%s_ctau_%s",Mchi2.back().c_str(),ctau.back().c_str()),"",1,0,10000));
    event_expWeight.push_back(new TH1F(Form("event_expWeight_mchi2_%s_ctau_%s",Mchi2.back().c_str(),ctau.back().c_str()),"",1,0,10000));

    event_inclusive.back()->Sumw2();
    event_leptonVeto.back()->Sumw2();
    event_tauVeto.back()->Sumw2();
    event_bVeto.back()->Sumw2();
    event_photonVeto.back()->Sumw2();
    event_caloMet.back()->Sumw2();
    event_jetCut.back()->Sumw2();
    event_jetCleaning.back()->Sumw2();
    event_jetMetDphi.back()->Sumw2();
    event_metCut.back()->Sumw2();
    event_monoVVeto.back()->Sumw2();
    event_expWeight.back()->Sumw2();

    final_templates.push_back(new TH1F(Form("MonoJ_mMED_%s_mchi2_%s_ctau_%s",Mmed.back().c_str(),Mchi2.back().c_str(),ctau.back().c_str()),"",bins.size()-1,&bins[0]));
    final_templates.back()->Sumw2();


    mediator_pt.push_back(new TH1F(Form("mediator_pt_mMED_%s_mchi2_%s_ctau_%s",Mmed.back().c_str(),Mchi2.back().c_str(),ctau.back().c_str()),"",75,0,1200));
    mediator_pt.back()->Sumw2();
    ratio_x2_pt.push_back(new TH1F(Form("ratio_x2_pt_mMED_%s_mchi2_%s_ctau_%s",Mmed.back().c_str(),Mchi2.back().c_str(),ctau.back().c_str()),"",75,0,5));
    ratio_x2_pt.back()->Sumw2();
    response_mediator_pt.push_back(new TH1F(Form("response_mediator_pt_mMED_%s_mchi2_%s_ctau_%s",Mmed.back().c_str(),Mchi2.back().c_str(),ctau.back().c_str()),"",100,-400,400));
    response_mediator_pt.back()->Sumw2();
    ratio_x1_pt.push_back(new TH1F(Form("ratio_x1_pt_mMED_%s_mchi2_%s_ctau_%s",Mmed.back().c_str(),Mchi2.back().c_str(),ctau.back().c_str()),"",75,0,5));
    ratio_x1_pt.back()->Sumw2();
    ratio_x1x2_pt.push_back(new TH1F(Form("ratio_x1x2_pt_mMED_%s_mchi2_%s_ctau_%s",Mmed.back().c_str(),Mchi2.back().c_str(),ctau.back().c_str()),"",75,0,1.5));
    ratio_x1x2_pt.back()->Sumw2();
    response_x1_pt.push_back(new TH1F(Form("response_x1_pt_mMED_%s_mchi2_%s_ctau_%s",Mmed.back().c_str(),Mchi2.back().c_str(),ctau.back().c_str()),"",100,-400,400));
    response_x1_pt.back()->Sumw2();

    ///////////////////////////
    while(reader.Next()){ // start applying selections

      event_inclusive.back()->Fill(*mmet,*xsec*lumi*(*wgt)/(sumOfWeights.at(ifile)));

      if(*nmuons != 0) continue;
      if(*nelectrons != 0) continue;

      event_leptonVeto.back()->Fill(*mmet,*xsec*lumi*(*wgt)/(sumOfWeights.at(ifile)));

      if(*ntaus != 0) continue;

      event_tauVeto.back()->Fill(*mmet,*xsec*lumi*(*wgt)/(sumOfWeights.at(ifile)));

      if(*nphotons != 0) continue;

      event_photonVeto.back()->Fill(*mmet,*xsec*lumi*(*wgt)/(sumOfWeights.at(ifile)));

      if(fabs(*mmet-*metcalo)/(*mmet) > 0.5) continue;

      event_caloMet.back()->Fill(*mmet,*xsec*lumi*(*wgt)/(sumOfWeights.at(ifile)));
      
      if(*njets < 1) continue;
      if(fabs(jeteta->at(0)) > 2.5) continue;
      if(jetpt->at(0) < 100) continue;
      
      event_jetCut.back()->Fill(*mmet,*xsec*lumi*(*wgt)/(sumOfWeights.at(ifile)));

      double puwgt   = puhist->GetBinContent(puhist->FindBin(*nvtx));
      double trgwgt  = triggermet_graph->Eval(min(double(*mmet),triggermet_graph->GetXaxis()->GetXmax()));

      mediator_pt.back()->Fill(*dmpt,puwgt*trgwgt*(*xsec)*lumi*(*wgt)/(sumOfWeights.at(ifile)));
      ratio_x2_pt.back()->Fill(*dmX1pt/(*dmX2pt),puwgt*trgwgt*(*xsec)*lumi*(*wgt)/(sumOfWeights.at(ifile)));
      response_mediator_pt.back()->Fill(*mmet-*dmpt,puwgt*trgwgt*(*xsec)*lumi*(*wgt)/(sumOfWeights.at(ifile)));
      ratio_x1x2_pt.back()->Fill(*dmX1pt_ll/(*dmX1pt),puwgt*trgwgt*(*xsec)*lumi*(*wgt)/(sumOfWeights.at(ifile)));
      ratio_x1x2_pt.back()->Fill(*dmX2pt_ll/(*dmX2pt),puwgt*trgwgt*(*xsec)*lumi*(*wgt)/(sumOfWeights.at(ifile)));      
      ratio_x1_pt.back()->Fill(*dmX1pt_ll/(*dmX2pt_ll),puwgt*trgwgt*(*xsec)*lumi*(*wgt)/(sumOfWeights.at(ifile)));

      TLorentzVector x1vec, x2vec;
      x1vec.SetPtEtaPhiM(*dmX1pt_ll,*dmX1eta_ll,*dmX1phi_ll,*dmX1mass_ll);
      x2vec.SetPtEtaPhiM(*dmX2pt_ll,*dmX2eta_ll,*dmX2phi_ll,*dmX2mass_ll);
      response_x1_pt.back()->Fill(*mmet-(x1vec+x2vec).Pt(),puwgt*trgwgt*(*xsec)*lumi*(*wgt)/(sumOfWeights.at(ifile)));


      if(*mmet < 250) continue;
      event_metCut.back()->Fill(*mmet,*xsec*lumi*(*wgt)/(sumOfWeights.at(ifile)));

      if(*nbjets != 0) continue;

      event_bVeto.back()->Fill(*mmet,*xsec*lumi*(*wgt)/(sumOfWeights.at(ifile)));

      if(chfrac->at(0) < 0.1) continue;
      if(nhfrac->at(0) > 0.8) continue;

      event_jetCleaning.back()->Fill(*mmet,*xsec*lumi*(*wgt)/(sumOfWeights.at(ifile)));

      if(*jmmdphi < 0.5) continue;

      event_jetMetDphi.back()->Fill(*mmet,*xsec*lumi*(*wgt)/(sumOfWeights.at(ifile)));

      if(category == Category::monojet){
	bool goodMonojet = true;
	if(boostedJetpt->size() > 0){
	  if(boostedJetpt->at(0) > 250 and
	     fabs(boostedJeteta->at(0)) < 2.4 and
	     boostedJettau2->at(0)/boostedJettau1->at(0) < 0.6 and
	     prunedJetm->at(0) > 65 and
	     prunedJetm->at(0) < 105) goodMonojet = false;
	}
	if(goodMonojet == false) continue;
	
	event_monoVVeto.back()->Fill(*mmet,*xsec*lumi*(*wgt)/(sumOfWeights.at(ifile)));
	
	if(*hltm90 == 0 or *hltm110 == 0 or *hltm120 == 0 or *hltmwm90 == 0 or *hltmwm120 == 0 or *hltmwm170 == 0 or *hltmwm300 == 0) continue;
	if(*fhbhe == 0 or *fhbiso == 0 or *feeb == 0 or *fetp == 0 or *fvtx == 0 or *fbadmu == 0 or *fbadch == 0 or *fcsc == 0) continue;
	
	event_expWeight.back()->Fill(*mmet,puwgt*trgwgt*(*xsec)*lumi*(*wgt)/(sumOfWeights.at(ifile)));                 
	final_templates.back()->Fill(*mmet,puwgt*trgwgt*(*xsec)*lumi*(*wgt)/(sumOfWeights.at(ifile)));
      }
      else if(category == Category::monoV){

	if(boostedJetpt->size() == 0) continue;
	if(boostedJetpt->at(0) < 250) continue;
	if(fabs(boostedJeteta->at(0)) > 2.4) continue;
	if(boostedJettau2->at(0)/boostedJettau1->at(0) > 0.6) continue;
	if(prunedJetm->at(0) < 65 or prunedJetm->at(0) > 105) continue;
	
	event_monoVVeto.back()->Fill(*mmet,*xsec*lumi*(*wgt)/(sumOfWeights.at(ifile)));

	if(*hltm90 == 0 or *hltm110 == 0 or *hltm120 == 0 or *hltmwm90 == 0 or *hltmwm120 == 0 or *hltmwm170 == 0 or *hltmwm300 == 0) continue;
	if(*fhbhe == 0 or *fhbiso == 0 or *feeb == 0 or *fetp == 0 or *fvtx == 0 or *fbadmu == 0 or *fbadch == 0 or *fcsc == 0) continue;
	
	event_expWeight.back()->Fill(*mmet,puwgt*trgwgt*(*xsec)*lumi*(*wgt)/(sumOfWeights.at(ifile)));                 
	final_templates.back()->Fill(*mmet,puwgt*trgwgt*(*xsec)*lumi*(*wgt)/(sumOfWeights.at(ifile)));

      }
    }
    ifile++;    
  }


  // calculate the efficiency
  vector<TH1F*> efficiency;
  TH1F* temp = (TH1F*) event_monoVVeto.front()->Clone("temp");

  for(int nhist = 0; nhist < event_expWeight.size(); nhist++){

    efficiency.push_back(new TH1F(TString(event_expWeight.at(nhist)->GetName()).ReplaceAll("event_expWeight","efficiency").Data(),"",11,0,12));
    efficiency.back()->Sumw2();

    temp->Divide(event_leptonVeto.at(nhist),event_inclusive.at(nhist),1,1,"B");
    efficiency.back()->SetBinContent(1,temp->GetBinContent(1));
    efficiency.back()->SetBinError(1,temp->GetBinError(1));
    efficiency.back()->GetXaxis()->SetBinLabel(1,"Lep veto");

    temp->Divide(event_tauVeto.at(nhist),event_inclusive.at(nhist),1,1,"B");
    efficiency.back()->SetBinContent(2,temp->GetBinContent(1));
    efficiency.back()->SetBinError(2,temp->GetBinError(1));
    efficiency.back()->GetXaxis()->SetBinLabel(2,"Tau veto");

    temp->Divide(event_photonVeto.at(nhist),event_inclusive.at(nhist),1,1,"B");
    efficiency.back()->SetBinContent(3,temp->GetBinContent(1));
    efficiency.back()->SetBinError(3,temp->GetBinError(1));
    efficiency.back()->GetXaxis()->SetBinLabel(3,"Photon veto");

    temp->Divide(event_caloMet.at(nhist),event_inclusive.at(nhist),1,1,"B");
    efficiency.back()->SetBinContent(4,temp->GetBinContent(1));
    efficiency.back()->SetBinError(4,temp->GetBinError(1));
    efficiency.back()->GetXaxis()->SetBinLabel(4,"Calo E_{T}^{miss}");

    temp->Divide(event_jetCut.at(nhist),event_inclusive.at(nhist),1,1,"B");
    efficiency.back()->SetBinContent(5,temp->GetBinContent(1));
    efficiency.back()->SetBinError(5,temp->GetBinError(1));
    efficiency.back()->GetXaxis()->SetBinLabel(5,"Jet Cut");

    temp->Divide(event_metCut.at(nhist),event_inclusive.at(nhist),1,1,"B");
    efficiency.back()->SetBinContent(6,temp->GetBinContent(1));
    efficiency.back()->SetBinError(6,temp->GetBinError(1));
    efficiency.back()->GetXaxis()->SetBinLabel(6,"E_{T}^{miss} cut");


    temp->Divide(event_bVeto.at(nhist),event_inclusive.at(nhist),1,1,"B");
    efficiency.back()->SetBinContent(7,temp->GetBinContent(1));
    efficiency.back()->SetBinError(7,temp->GetBinError(1));
    efficiency.back()->GetXaxis()->SetBinLabel(7,"B-jet veto");

    temp->Divide(event_jetCleaning.at(nhist),event_inclusive.at(nhist),1,1,"B");
    efficiency.back()->SetBinContent(8,temp->GetBinContent(1));
    efficiency.back()->SetBinError(8,temp->GetBinError(1));
    efficiency.back()->GetXaxis()->SetBinLabel(8,"Jet Cleaning");

    temp->Divide(event_jetMetDphi.at(nhist),event_inclusive.at(nhist),1,1,"B");
    efficiency.back()->SetBinContent(9,temp->GetBinContent(1));
    efficiency.back()->SetBinError(9,temp->GetBinError(1));
    efficiency.back()->GetXaxis()->SetBinLabel(9,"#Delta#phi_{jet,E_{T}^{miss}}");

    temp->Divide(event_monoVVeto.at(nhist),event_inclusive.at(nhist),1,1,"B");
    efficiency.back()->SetBinContent(10,temp->GetBinContent(1));
    efficiency.back()->SetBinError(10,temp->GetBinError(1));
    if(category == Category::monojet)
      efficiency.back()->GetXaxis()->SetBinLabel(10,"Mono-V veto");
    else if(category == Category::monoV)
      efficiency.back()->GetXaxis()->SetBinLabel(10,"Mono-V selection");
    
    temp->Divide(event_expWeight.at(nhist),event_inclusive.at(nhist),1,1,"B");
    efficiency.back()->SetBinContent(11,temp->GetBinContent(1));
    efficiency.back()->SetBinError(11,temp->GetBinError(1));
    efficiency.back()->GetXaxis()->SetBinLabel(11,"Exp. weights");

    efficiency.back()->GetXaxis()->LabelsOption("v");
  }

  /// Start plotting part
  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();
  canvas->SetBottomMargin(0.2);
  
  map<string,vector<TH1F*> > effciency_vs_ctau; // group all efficiency relative to the same mchi2
  map<string,vector<TH1F*> > effciency_vs_mchi2; // group all efficiency relative to the same ctau
  vector<string> label_mchi2;
  vector<string> label_ctau;

  for(int ieff = 0; ieff < efficiency.size();  ieff++){

    stringstream fileName(efficiency.at(ieff)->GetName());
    std::string segment;
    std::vector<std::string> seglist;

    while(std::getline(fileName, segment, '_')){
      seglist.push_back(segment);
    }

    if(std::find(label_mchi2.begin(),label_mchi2.end(),"M_{#chi2} = "+seglist.at(2)+" GeV") == label_mchi2.end())
      label_mchi2.push_back("M_{#chi2} = "+seglist.at(2)+" GeV");
    effciency_vs_ctau[seglist.at(2)].push_back(efficiency.at(ieff));
    if(std::find(label_ctau.begin(),label_ctau.end(),"c#tau = "+seglist.at(4)+" mm") == label_ctau.end())
      label_ctau.push_back("c#tau = "+seglist.at(4)+" mm");
    effciency_vs_mchi2[seglist.at(4)].push_back(efficiency.at(ieff));
    
  }
  
  // plot
  int ieff = 0;
  for(auto eff: effciency_vs_ctau)
    plotEfficiency(canvas,eff.second,label_ctau,"M_{#chi2} = "+eff.first+" GeV",outputDIR,lumi,"mchi_"+eff.first);
  // plot
  for(auto eff: effciency_vs_mchi2)
    plotEfficiency(canvas,eff.second,label_mchi2,"c#tau = "+eff.first+" mm",outputDIR,lumi,"ctau_"+eff.first);

  // Output file
  TFile* outputFile = new TFile((outputDIR+"/"+outputFileName+".root").c_str(),"RECREATE");
  outputFile->cd();  
  for(auto hist : final_templates) hist->Write();
  outputFile->mkdir("Observables");
  outputFile->cd("Observables");
  for(auto hist : mediator_pt){ hist->Scale(1./hist->Integral()); hist->Write();}
  for(auto hist : ratio_x2_pt) { hist->Scale(1./hist->Integral()); hist->Write();}
  for(auto hist : response_mediator_pt) { hist->Scale(1./hist->Integral()); hist->Write();}
  for(auto hist : ratio_x1_pt) { hist->Scale(1./hist->Integral()); hist->Write();}
  for(auto hist : ratio_x1x2_pt) { hist->Scale(1./hist->Integral()); hist->Write();}
  for(auto hist : response_x1_pt) { hist->Scale(1./hist->Integral()); hist->Write();}
  outputFile->cd();
  outputFile->Close();
}
