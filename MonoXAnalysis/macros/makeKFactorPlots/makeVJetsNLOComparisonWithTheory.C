#include "../CMS_lumi.h"

enum class Sample {zvv,wjet,zll,gam};

static float minX = 250;
static float maxX = 1400;

void makePlot(TH1* histo_cms, TH1* histo_theo, const string & outputDirectory, const Sample & sample, const string & postfix){

  TCanvas* canvas = new TCanvas(("canvas_"+postfix).c_str(),"",600,650);
  canvas->SetBottomMargin(0.3);

  TPad *pad = new TPad(("pad_"+postfix).c_str(),"pad",0,0.,1,0.9);
  pad->SetTopMargin(0.7);
  pad->SetFillColor(0);
  pad->SetGridy(1);
  pad->SetFillStyle(0);
  canvas->cd();

  histo_theo->GetXaxis()->SetTitleSize(0);
  histo_theo->GetXaxis()->SetLabelSize(0);
  histo_theo->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb/GeV]");
  histo_theo->GetYaxis()->SetTitleOffset(1.40);
  histo_theo->GetXaxis()->SetRangeUser(minX,maxX);

  histo_cms->SetLineColor(kRed);
  histo_cms->SetMarkerColor(kRed);
  histo_cms->SetMarkerStyle(20);
  histo_cms->SetMarkerSize(1);
  histo_cms->SetLineWidth(2);

  histo_theo->SetLineColor(kBlack);
  histo_theo->SetMarkerColor(kBlack);
  histo_theo->SetMarkerStyle(20);
  histo_theo->SetMarkerSize(1);
  histo_theo->SetLineWidth(2);
  histo_theo->GetYaxis()->SetRangeUser(min(histo_theo->GetBinContent(histo_theo->FindBin(maxX)),histo_cms->GetBinContent(histo_cms->FindBin(maxX)))*0.1,
				       max(histo_theo->GetBinContent(histo_theo->FindBin(minX)),histo_cms->GetBinContent(histo_cms->FindBin(minX)))*100);

  histo_theo->Draw("hist");
  TH1F* histo_theo_band = (TH1F*) histo_theo->Clone("histo_theo_band");
  histo_theo_band->SetFillColor(kGray);
  histo_theo_band->SetFillStyle(1001);
  histo_theo_band->Draw("E2same");
  histo_theo->Draw("hist same");
  histo_cms->Draw("hist same");
  CMS_lumi(canvas,"");

  TLegend leg (0.6,0.75,0.9,0.88);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  if(sample == Sample::zvv){
    leg.AddEntry(histo_theo_band,"Z+jets: theorists","FL");
    leg.AddEntry(histo_cms,"Z+jets: CMS","PL");
  }
  if(sample == Sample::wjet){
    leg.AddEntry(histo_theo_band,"W+jets: theorists","FL");
    leg.AddEntry(histo_cms,"W+jets: CMS","PL");
  }
  if(sample == Sample::zll){
    leg.AddEntry(histo_theo_band,"DY+jets: theorists","FL");
    leg.AddEntry(histo_cms,"DY+jets: CMS","PL");
  }
  if(sample == Sample::gam){
    leg.AddEntry(histo_theo_band,"#gamma+jets: theorists","FL");
    leg.AddEntry(histo_cms,"#gamma+jets: CMS","PL");
  }
  leg.Draw("same");
  canvas->SetLogy();

  pad->Draw();
  pad->cd();
  
  TH1F* ratio = (TH1F*) histo_theo->Clone("ratio");
  ratio->Divide(histo_cms);

  // set the error by hand to zero
  for(int iBin = 0; iBin <= ratio->GetNbinsX(); iBin++){
    ratio->SetBinError(iBin+1,0.);
  }

  ratio->SetLineColor(kRed);
  ratio->SetMarkerColor(kRed);
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerSize(1);
  ratio->SetLineWidth(2);
  ratio->GetXaxis()->SetLabelSize(0.04);
  ratio->GetXaxis()->SetTitleSize(0.04);

  TH1F* ratio_band = (TH1F*) histo_theo->Clone("ratio_band");
  ratio_band->Divide(histo_cms);
  for(int iBin = 0; iBin <= ratio_band->GetNbinsX(); iBin++)
    ratio_band->SetBinContent(iBin+1,1.);
  ratio_band->SetFillColor(kGray);
  ratio_band->SetFillStyle(1001);
    
  ratio->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
  ratio->GetXaxis()->SetTitleOffset(1.15);
  ratio->GetYaxis()->SetTitle("Theory/CMS");
  ratio->GetYaxis()->SetTitleOffset(1.30);
  ratio->GetYaxis()->SetNdivisions(504);
  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(2);
  ratio->GetYaxis()->SetRangeUser(0.7,1.3);
  ratio->GetXaxis()->SetRangeUser(minX,maxX);
  ratio->Draw("hist");
  ratio_band->SetMarkerSize(0);
  ratio_band->Draw("E2same");

  TF1* line = new TF1("line","1",minX,maxX);
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("L same");
  ratio->Draw("hist same");

  pad->RedrawAxis("same axis");
  
  canvas->cd();
  canvas->RedrawAxis("sameaxis");
  
  if(sample == Sample::zvv){
    canvas->SaveAs((outputDirectory+"/bosonpt_zvv_"+postfix+".png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/bosonpt_zvv_"+postfix+".pdf").c_str(),"pdf");
  }
  else if(sample == Sample::wjet){
    canvas->SaveAs((outputDirectory+"/bosonpt_wjet_"+postfix+".png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/bosonpt_wjet_"+postfix+".pdf").c_str(),"pdf");
  }
  else if(sample == Sample::zll){
    canvas->SaveAs((outputDirectory+"/bosonpt_zll_"+postfix+".png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/bosonpt_zll_"+postfix+".pdf").c_str(),"pdf");
  }
  else if(sample == Sample::gam){
    canvas->SaveAs((outputDirectory+"/bosonpt_gam_"+postfix+".png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/bosonpt_gam_"+postfix+".pdf").c_str(),"pdf");
  }

}

///////////////////
void makeVJetsNLOComparisonWithTheory(string inputCMSDirectory, string inputTheoryFile, string outputDirectory, Sample sample, bool isPhillTree = false, bool useFixedCone = false){

  TH1F* correction = NULL;
  if(not isPhillTree and not useFixedCone){
    TFile* input = TFile::Open("dynamic_cone_correction.root","READ");
    correction   = (TH1F*) input->Get("dynamic_over_fixed");
  }

  system(("mkdir -p "+outputDirectory).c_str());
  setTDRStyle();
  gROOT->SetBatch(kTRUE);

  float scale = 1;
  if(sample == Sample::wjet or sample == Sample::zll) scale = 0.5;
    
  //////////  
  TFile* file_theory = TFile::Open(inputTheoryFile.c_str(),"READ");
  TH1F* vjet_lo        = NULL; 
  TH1F* vjet_lo_kfact  = NULL; 
  TH1F* vjet_nlo_kfact = NULL;
  TH1F* vjet_nlo       = NULL;
  TH1F* vjet_d1K_NLO = NULL;
  TH1F* vjet_d2K_NLO = NULL;
  TH1F* vjet_d3K_NLO = NULL;
  TH1F* vjet_d1K_EWK = NULL;
  TH1F* vjet_d2K_EWK = NULL;
  TH1F* vjet_d3K_EWK = NULL;
  TH1F* vjet_d3K_NLO_mix = NULL;


  if(sample == Sample::zvv){
    vjet_lo = (TH1F*) file_theory->Get("vvj_pTV_LO");
    vjet_lo_kfact = (TH1F*) file_theory->Get("vvj_pTV_K_LO");
    vjet_nlo_kfact = (TH1F*) file_theory->Get("vvj_pTV_K_NLO");

    vjet_d1K_NLO = (TH1F*) file_theory->Get("vvj_pTV_d1K_NLO");
    vjet_d2K_NLO = (TH1F*) file_theory->Get("vvj_pTV_d2K_NLO");
    vjet_d3K_NLO = (TH1F*) file_theory->Get("vvj_pTV_d3K_NLO");
    vjet_d1K_EWK = (TH1F*) file_theory->Get("vvj_pTV_d1kappa_EW");
    vjet_d2K_EWK = (TH1F*) file_theory->Get("vvj_pTV_d2kappa_EW");
    vjet_d3K_EWK = (TH1F*) file_theory->Get("vvj_pTV_d3kappa_EW");
    vjet_d3K_NLO_mix = (TH1F*) file_theory->Get("vvj_pTV_dK_NLO_mix");
  }
  else if(sample == Sample::wjet){
    vjet_lo = (TH1F*) file_theory->Get("evj_pTV_LO");
    vjet_lo_kfact = (TH1F*) file_theory->Get("evj_pTV_K_LO");
    vjet_nlo_kfact = (TH1F*) file_theory->Get("evj_pTV_K_NLO");

    vjet_d1K_NLO = (TH1F*) file_theory->Get("evj_pTV_d1K_NLO");
    vjet_d2K_NLO = (TH1F*) file_theory->Get("evj_pTV_d2K_NLO");
    vjet_d3K_NLO = (TH1F*) file_theory->Get("evj_pTV_d3K_NLO");
    vjet_d1K_EWK = (TH1F*) file_theory->Get("evj_pTV_d1kappa_EW");
    vjet_d2K_EWK = (TH1F*) file_theory->Get("evj_pTV_d2kappa_EW");
    vjet_d3K_EWK = (TH1F*) file_theory->Get("evj_pTV_d3kappa_EW");
    vjet_d3K_NLO_mix = (TH1F*) file_theory->Get("evj_pTV_dK_NLO_mix");

  }
  else if(sample == Sample::zll){
    vjet_lo = (TH1F*) file_theory->Get("eej_pTV_LO");
    vjet_lo_kfact = (TH1F*) file_theory->Get("eej_pTV_K_LO");
    vjet_nlo_kfact = (TH1F*) file_theory->Get("eej_pTV_K_NLO");

    vjet_d1K_NLO = (TH1F*) file_theory->Get("eej_pTV_d1K_NLO");
    vjet_d2K_NLO = (TH1F*) file_theory->Get("eej_pTV_d2K_NLO");
    vjet_d3K_NLO = (TH1F*) file_theory->Get("eej_pTV_d3K_NLO");
    vjet_d1K_EWK = (TH1F*) file_theory->Get("eej_pTV_d1kappa_EW");
    vjet_d2K_EWK = (TH1F*) file_theory->Get("eej_pTV_d2kappa_EW");
    vjet_d3K_EWK = (TH1F*) file_theory->Get("eej_pTV_d3kappa_EW");
    vjet_d3K_NLO_mix = (TH1F*) file_theory->Get("eej_pTV_dK_NLO_mix");

  }
  else if(sample == Sample::gam){
    vjet_lo = (TH1F*) file_theory->Get("aj_pTV_LO");
    vjet_lo_kfact = (TH1F*) file_theory->Get("aj_pTV_K_LO");
    vjet_nlo_kfact = (TH1F*) file_theory->Get("aj_pTV_K_NLO");
    //vjet_nlo_kfact = (TH1F*) file_theory->Get("aj_pTV_K_NLO_fix");
    vjet_d1K_NLO = (TH1F*) file_theory->Get("aj_pTV_d1K_NLO");
    vjet_d2K_NLO = (TH1F*) file_theory->Get("aj_pTV_d2K_NLO");
    vjet_d3K_NLO = (TH1F*) file_theory->Get("aj_pTV_d3K_NLO");
    vjet_d1K_EWK = (TH1F*) file_theory->Get("aj_pTV_d1kappa_EW");
    vjet_d2K_EWK = (TH1F*) file_theory->Get("aj_pTV_d2kappa_EW");
    vjet_d3K_EWK = (TH1F*) file_theory->Get("aj_pTV_d3kappa_EW");
    vjet_d3K_NLO_mix = (TH1F*) file_theory->Get("aj_pTV_dK_NLO_mix");

  }

  vjet_nlo = (TH1F*) vjet_lo->Clone("vjet_nlo");
  vjet_nlo->Multiply(vjet_lo_kfact);
  vjet_nlo->Multiply(vjet_nlo_kfact);
  
  // set process errors
  for(int iBin = 0; iBin <= vjet_nlo->GetNbinsX(); iBin++){
    vjet_nlo->SetBinError(iBin+1,vjet_nlo->GetBinContent(iBin+1)*sqrt(pow(vjet_d1K_NLO->GetBinContent(iBin+1),2)+pow(vjet_d2K_NLO->GetBinContent(iBin+1),2)+
								      pow(vjet_d3K_NLO->GetBinContent(iBin+1),2)+pow(vjet_d1K_EWK->GetBinContent(iBin+1),2)+
								      pow(vjet_d2K_EWK->GetBinContent(iBin+1),2)+pow(vjet_d3K_EWK->GetBinContent(iBin+1),2)+
								      pow(vjet_d3K_NLO_mix->GetBinContent(iBin+1),2)));
  }


  //////
  vector<TTree*> tree_vjet_cms;
  vector<TFile*> file_vjet_cms;

  TH1F* vjet_nlo_cms_lhe   = (TH1F*) vjet_nlo->Clone("vjet_nlo_cms_lhe");
  TH1F* vjet_nlo_cms_ps    = (TH1F*) vjet_nlo->Clone("vjet_nlo_cms_ps");
  vjet_nlo_cms_lhe->Reset();
  vjet_nlo_cms_ps->Reset();

  TH1F* difference = new TH1F("difference","",100,-100,100);
  difference->Sumw2();

  system(("ls "+inputCMSDirectory+" | grep root > list.txt").c_str());
  ifstream infileQCD ("list.txt");
  string line;
  if(infileQCD.is_open()){
    while(!infileQCD.eof()){
      getline(infileQCD,line);
      if(not TString(line).Contains("root") or line == "") continue;
      file_vjet_cms.push_back(TFile::Open((inputCMSDirectory+"/"+line).c_str(),"READ"));
      if(not isPhillTree)
	tree_vjet_cms.push_back((TTree*) file_vjet_cms.back()->Get("gentree/tree"));
      else
	tree_vjet_cms.push_back((TTree*) file_vjet_cms.back()->Get("Events"));
	
    }
  }
  infileQCD.close();

  system("rm list.txt");
  
  // calculate sumwgt                                                                                                                                                                                 
  vector<double> sumwgt_vjet_cms;
  int ifile = 0;
  double eventsNoBoson = 0;

  for(auto tree : tree_vjet_cms){
    eventsNoBoson = 0;

    if(not isPhillTree){      
      TTreeReader reader (tree);
      TTreeReaderValue<float> wgt   (reader,"wgt");
      cout<<"Calculate sumwgt for QCD file "<<file_vjet_cms.at(ifile)->GetName()<<endl;
      double sumwgt = 0;
      while(reader.Next()){
	sumwgt += *wgt;
      }
      
      cout<<"Tree CMS with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<endl;
      sumwgt_vjet_cms.push_back(sumwgt);
      ifile++;
    }
    else{
      
      TTreeReader reader (tree);
      TTreeReaderArray<float> wgt   (reader,"evtweight.fEvtWeight");
      cout<<"Calculate sumwgt for QCD file "<<file_vjet_cms.at(ifile)->GetName()<<endl;
      double sumwgt = 0;
      while(reader.Next()){
	sumwgt += wgt[0];
      }      
      cout<<"Tree CMS with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<endl;
      sumwgt_vjet_cms.push_back(sumwgt);
      ifile++;
    }
  }

  /// start making selections                                                                                                                                                                        
  ifile = 0;
  for(auto tree: tree_vjet_cms){

    TTreeReader reader (tree);    
    if(not isPhillTree){

      TTreeReaderValue<float> xsec       (reader,"xsec");
      TTreeReaderValue<float> wgt        (reader,"wgt");
      TTreeReaderValue<int>   l1id       (reader,"l1id");
      TTreeReaderValue<int>   l2id       (reader,"l2id");
      // post-showered i.e. status 1
      TTreeReaderValue<float> wzmass_ps  (reader,"wzmass");
      TTreeReaderValue<float> wzeta_ps   (reader,"wzeta");
      TTreeReaderValue<float> wzphi_ps   (reader,"wzphi");
      TTreeReaderValue<float> wzpt_ps    (reader,"wzpt");
      TTreeReaderValue<int>   wzid_ps    (reader,"wzid");
      // pre-showered i.e. status 23
      TTreeReaderValue<float> wzmass_lhe (reader,"mvmass");
      TTreeReaderValue<float> wzeta_lhe  (reader,"mveta");
      TTreeReaderValue<float> wzphi_lhe  (reader,"mvphi");
      TTreeReaderValue<float> wzpt_lhe   (reader,"mvpt");
      TTreeReaderValue<int>   wzid_lhe   (reader,"mvid");
      
      cout<<"Loop on V+jets file "<<file_vjet_cms.at(ifile)->GetName()<<endl;      
      while(reader.Next()){     
	
	// filter away bad events with no matching                                                                                                                                                 
	if((sample == Sample::zvv or sample == Sample::zll) and fabs(*wzid_lhe) != 23){
	  eventsNoBoson++;
	  continue;
	}

	if(sample == Sample::wjet and fabs(*wzid_lhe) != 24){
	  eventsNoBoson++;
	  continue;
	}

	if(sample == Sample::gam and fabs(*wzid_lhe) != 22){
	  eventsNoBoson++;
	  continue;
	}

	// filter away bad events with no matching                                                                                                                                                 
	if((sample == Sample::zvv or sample == Sample::zll) and fabs(*wzid_ps) != 23){
	  eventsNoBoson++;
	  continue;
	}

	if(sample == Sample::wjet and fabs(*wzid_ps) != 24){
	  eventsNoBoson++;
	  continue;
	}

	if(sample == Sample::gam and fabs(*wzid_ps) != 22){
	  eventsNoBoson++;
	  continue;
	}
	
	// min-boson pt
	if(*wzpt_lhe < 30) continue;
	if(*wzpt_ps  < 30) continue;
	//remove tau events
	if(sample == Sample::wjet and fabs(*l1id) == 15) continue;
	// lower bound on the invariant mass
	if((sample == Sample::zvv or sample == Sample::zll) and *wzmass_lhe < 30) continue;
	if((sample == Sample::zvv or sample == Sample::zll) and *wzmass_ps  < 30) continue;
	
	//////////////////                                                                                                                                                                     
	vjet_nlo_cms_lhe->Fill(*wzpt_lhe,(*wgt)*(*xsec)*0.001*scale/sumwgt_vjet_cms.at(ifile)); // cross section in pb
	vjet_nlo_cms_ps->Fill (*wzpt_ps,(*wgt)*(*xsec)*0.001*scale/sumwgt_vjet_cms.at(ifile)); // cross section in pb
	difference->Fill(*wzpt_ps-*wzpt_lhe,(*wgt)*(*xsec)*0.001*scale/sumwgt_vjet_cms.at(ifile));
      }
      ifile++;
      cout<<"Event no Boson "<<eventsNoBoson<<" total events "<<tree->GetEntries()<<endl;
      eventsNoBoson = 0;
    }
    else{

      TTreeReaderArray<float> xsec  (reader,"xs.fXS");
      TTreeReaderArray<float> wgt   (reader,"evtweight.fEvtWeight");
      TTreeReaderArray<int>   l1id       (reader,"id_1.fId1");
      TTreeReaderArray<int>   l2id       (reader,"id_2.fId2");
      // post-shower boson
      TTreeReaderArray<int>   wzid_ps    (reader,"v_id.fVId");
      TTreeReaderArray<float> wzmass_ps  (reader,"v_m.fVM");
      TTreeReaderArray<float> wzpt_ps    (reader,"v_pt.fVPt");
      TTreeReaderArray<float> wzeta_ps   (reader,"v_y.fVEta");
      TTreeReaderArray<float> wzphi_ps   (reader,"v_phi.fVPhi");
      TTreeReaderArray<int>   wzstatus_ps  (reader,"v_status.fVStatus");
      //
      TTreeReaderArray<int>   wzid_lhe    (reader,"vme_id.fVMEId");
      TTreeReaderArray<float> wzmass_lhe  (reader,"vme_m.fVMEM");
      TTreeReaderArray<float> wzpt_lhe    (reader,"vme_pt.fVMEPt");
      TTreeReaderArray<float> wzeta_lhe   (reader,"vme_y.fVMEEta");
      TTreeReaderArray<float> wzphi_lhe   (reader,"vme_phi.fVMEPhi");
      TTreeReaderArray<int>   wzstatus_lhe  (reader,"vme_status.fVMEStatus");

      TTreeReaderArray<float> gamiso_04  (reader,"vme_iso04.fVMEIso");
      TTreeReaderArray<float> gamiso_dy  (reader,"vme_dyniso.fDVMEIso");

      cout<<"Loop on V+jets file "<<file_vjet_cms.at(ifile)->GetName()<<endl;      
      while(reader.Next()){     
	
	// filter away bad events with no matching                                                                                                                                                 
	if((sample == Sample::zvv or sample == Sample::zll) and fabs(wzid_lhe[0]) != 23){
	  eventsNoBoson++;
	  continue;
	}

	if(sample == Sample::wjet and fabs(wzid_lhe[0]) != 24){
	  eventsNoBoson++;
	  continue;
	}

	if(sample == Sample::gam and fabs(wzid_lhe[0]) != 22){
	  eventsNoBoson++;
	  continue;
	}

	// filter away bad events with no matching                                                                                                                                                 
	if((sample == Sample::zvv or sample == Sample::zll) and fabs(wzid_ps[0]) != 23){
	  eventsNoBoson++;
	  continue;
	}

	if(sample == Sample::wjet and fabs(wzid_ps[0]) != 24){
	  eventsNoBoson++;
	  continue;
	}

	if(sample == Sample::gam and fabs(wzid_ps[0]) != 22){
	  eventsNoBoson++;
	  continue;
	}

	// check the right status
	if(wzstatus_lhe[0] != 23 ) continue;
	      
	// min-boson pt
	if(wzpt_lhe[0] < 30) continue;
	if(wzpt_ps[0] < 30) continue;
	//remove tau events
	if(sample == Sample::wjet and (fabs(l1id[0]) == 15 or fabs(l2id[0]) == 15)) continue;
	//-----
	if((sample == Sample::zvv or sample == Sample::zll) and wzmass_lhe[0] < 30) continue;
	if((sample == Sample::zvv or sample == Sample::zll) and wzmass_ps[0] < 30) continue;
	//------
	if(useFixedCone){
	  if(sample == Sample::gam and gamiso_04[0] > 1.0*wzpt_lhe[0]) continue;
	}
	else{
	  if(sample == Sample::gam and gamiso_dy[0] > 0.1*wzpt_lhe[0]) continue;
	}

	//////////////////                                                                                                                                                                     
	vjet_nlo_cms_lhe->Fill(wzpt_lhe[0],(wgt[0])*(xsec[0])*scale/sumwgt_vjet_cms.at(ifile)); // cross section in pb
	vjet_nlo_cms_ps->Fill(wzpt_ps[0],(wgt[0])*(xsec[0])*scale/sumwgt_vjet_cms.at(ifile)); // cross section in pb	
	difference->Fill(wzpt_ps[0]-wzpt_lhe[0],(wgt[0])*(xsec[0])*scale/sumwgt_vjet_cms.at(ifile));
      }
      ifile++;
      cout<<"Event no Boson "<<eventsNoBoson<<" total events "<<tree->GetEntries()<<endl;
      eventsNoBoson = 0;
    }
  }


  vjet_nlo_cms_lhe->Scale(1,"width");
  vjet_nlo_cms_ps->Scale(1,"width");

  if(correction and sample == Sample::gam){
    for(int iBin = 0; iBin <= vjet_nlo_cms_lhe->GetNbinsX(); iBin++){
      if(correction->FindBin(vjet_nlo_cms_lhe->GetBinCenter(iBin+1)) != 0 and correction->FindBin(vjet_nlo_cms_lhe->GetBinCenter(iBin+1)) != vjet_nlo_cms_lhe->GetNbinsX()){
	if(vjet_nlo_cms_lhe->GetBinCenter(iBin+1) < 400){
	  vjet_nlo_cms_lhe->SetBinContent(iBin+1,vjet_nlo_cms_lhe->GetBinContent(iBin+1)*correction->GetBinContent(correction->FindBin(vjet_nlo_cms_lhe->GetBinCenter(iBin+1)))*0.95);
	  vjet_nlo_cms_ps->SetBinContent(iBin+1,vjet_nlo_cms_ps->GetBinContent(iBin+1)*correction->GetBinContent(correction->FindBin(vjet_nlo_cms_ps->GetBinCenter(iBin+1)))*0.95);
	}
	else{
	  vjet_nlo_cms_lhe->SetBinContent(iBin+1,vjet_nlo_cms_lhe->GetBinContent(iBin+1)*correction->GetBinContent(correction->FindBin(vjet_nlo_cms_lhe->GetBinCenter(iBin+1))));
	  vjet_nlo_cms_ps->SetBinContent(iBin+1,vjet_nlo_cms_ps->GetBinContent(iBin+1)*correction->GetBinContent(correction->FindBin(vjet_nlo_cms_ps->GetBinCenter(iBin+1))));
	}
      }
    }
  }

  string outputFileName;  
  if(sample == Sample::zvv)
    outputFileName = "bosonpt_zvv";
  else if(sample == Sample::wjet)
    outputFileName = "bosonpt_wjet";
  else if(sample == Sample::zll)
    outputFileName = "bosonpt_zll";
  else if(sample == Sample::gam)
    outputFileName = "bosonpt_gam";
  

  TFile* outputFile = new TFile((outputDirectory+"/"+outputFileName+".root").c_str(),"RECREATE");
  outputFile->cd();
  // make plot
  makePlot(vjet_nlo_cms_lhe,vjet_nlo,outputDirectory,sample,"lhe");
  makePlot(vjet_nlo_cms_ps,vjet_nlo,outputDirectory,sample,"ps");
  
  // difference lhe vs ps boson pt
  TCanvas* canvas = new TCanvas("canvas","",600,600);

  difference->GetXaxis()->SetTitle("(PS-LHE) Boson p_{T}");
  difference->GetYaxis()->SetTitle("Entries");
  difference->GetYaxis()->SetTitleOffset(1.30);
  
  difference->SetLineColor(kBlack);
  difference->SetMarkerColor(kBlack);
  difference->SetMarkerStyle(20);
  difference->SetMarkerSize(1);
  difference->SetLineWidth(2);
  difference->Draw("hist");

  CMS_lumi(canvas,"");
  canvas->SetLogy();

  canvas->SaveAs((outputDirectory+"/difference_boson_pt.png").c_str(),"png");
  canvas->SaveAs((outputDirectory+"/difference_boson_pt.pdf").c_str(),"pdf");

  vjet_nlo_cms_lhe->Write((outputFileName+"_cms_lhe").c_str());
  vjet_nlo_cms_ps->Write((outputFileName+"_cms_ps").c_str());
  vjet_nlo->Write((outputFileName+"_theo").c_str());

  outputFile->Close();

}
