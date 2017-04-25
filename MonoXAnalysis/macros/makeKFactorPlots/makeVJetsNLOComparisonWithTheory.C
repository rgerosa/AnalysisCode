#include "../CMS_lumi.h"

enum class Sample {zvv,wjet,zll,gam};

void makePlot(TH1* histo_cms, TH1* histo_theo, const string & outputDirectory, const Sample & sample, const string & postfix){

  TCanvas* canvas = new TCanvas(("canvas_"+postfix).c_str(),"",600,600);
  canvas->SetBottomMargin(0.3);
  TPad *pad = new TPad(("pad_"+postfix).c_str(),"pad",0,0.,1,0.9);
  pad->SetTopMargin(0.7);
  pad->SetRightMargin(0.06);
  pad->SetFillColor(0);
  pad->SetGridy(1);
  pad->SetFillStyle(0);
  canvas->cd();

  histo_cms->GetXaxis()->SetTitleSize(0);
  histo_cms->GetXaxis()->SetLabelSize(0);
  histo_cms->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb/GeV]");
  histo_cms->GetYaxis()->SetTitleOffset(1.30);
  histo_cms->GetXaxis()->SetRangeUser(200,1400);

  histo_theo->SetLineColor(kRed);
  histo_theo->SetMarkerColor(kRed);
  histo_theo->SetMarkerStyle(20);
  histo_theo->SetMarkerSize(1);
  histo_theo->SetLineWidth(2);

  histo_cms->SetLineColor(kBlack);
  histo_cms->SetMarkerColor(kBlack);
  histo_cms->SetMarkerStyle(20);
  histo_cms->SetMarkerSize(1);
  histo_cms->SetLineWidth(2);
  histo_cms->Draw("hist");

  TH1F* histo_cms_band = (TH1F*) histo_cms->Clone("histo_cms_band");
  histo_cms_band->SetFillColor(kBlack);
  histo_cms_band->SetFillStyle(3001);
  histo_cms_band->Draw("E2same");
  histo_cms->Draw("hist same");

  histo_theo->Draw("hist same");
  CMS_lumi(canvas,"");

  TLegend leg (0.7,0.7,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  if(sample == Sample::zvv){
    leg.AddEntry(histo_theo,"Z+jets: theorists","PL");
    leg.AddEntry(histo_cms,"Z+jets: CMS","PL");
  }
  if(sample == Sample::wjet){
    leg.AddEntry(histo_theo,"W+jets: theorists","PL");
    leg.AddEntry(histo_cms,"W+jets: CMS","PL");
  }
  if(sample == Sample::zll){
    leg.AddEntry(histo_theo,"DY+jets: theorists","PL");
    leg.AddEntry(histo_cms,"DY+jets: CMS","PL");
  }
  if(sample == Sample::gam){
    leg.AddEntry(histo_theo,"#gamma+jets: theorists","PL");
    leg.AddEntry(histo_cms,"#gamma+jets: CMS","PL");
  }
  leg.Draw("same");
  canvas->SetLogy();

  pad->Draw();
  pad->cd();
  
  TH1F* ratio = (TH1F*) histo_theo->Clone("ratio");
  ratio->Divide(histo_cms);
  ratio->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
  ratio->GetXaxis()->SetTitleOffset(1.15);
  ratio->GetYaxis()->SetTitle("Theory/CMS");
  ratio->GetYaxis()->SetTitleOffset(1.30);
  ratio->GetYaxis()->SetNdivisions(504);
  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(2);
  ratio->GetYaxis()->SetRangeUser(0.8,1.2);
  ratio->GetXaxis()->SetRangeUser(200,1400);
  ratio->Draw("EP");
  
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
void makeVJetsNLOComparisonWithTheory(string inputCMSDirectory, string inputTheoryFile, string outputDirectory, Sample sample, bool isPhillTree = false){

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
  if(sample == Sample::zvv){
    vjet_lo = (TH1F*) file_theory->Get("vvj_pTV_LO");
    vjet_lo_kfact = (TH1F*) file_theory->Get("vvj_pTV_K_LO");
    vjet_nlo_kfact = (TH1F*) file_theory->Get("vvj_pTV_K_NLO");
  }
  else if(sample == Sample::wjet){
    vjet_lo = (TH1F*) file_theory->Get("evj_pTV_LO");
    vjet_lo_kfact = (TH1F*) file_theory->Get("evj_pTV_K_LO");
    vjet_nlo_kfact = (TH1F*) file_theory->Get("evj_pTV_K_NLO");
  }
  else if(sample == Sample::zll){
    vjet_lo = (TH1F*) file_theory->Get("eej_pTV_LO");
    vjet_lo_kfact = (TH1F*) file_theory->Get("eej_pTV_K_LO");
    vjet_nlo_kfact = (TH1F*) file_theory->Get("eej_pTV_K_NLO");
  }
  else if(sample == Sample::gam){
    vjet_lo = (TH1F*) file_theory->Get("aj_pTV_LO");
    vjet_lo_kfact = (TH1F*) file_theory->Get("aj_pTV_K_LO");
    //vjet_nlo_kfact = (TH1F*) file_theory->Get("aj_pTV_K_NLO");
    vjet_nlo_kfact = (TH1F*) file_theory->Get("aj_pTV_K_NLO_fix");
  }

  vjet_nlo = (TH1F*) vjet_lo->Clone("vjet_nlo");
  vjet_nlo->Multiply(vjet_lo_kfact);
  vjet_nlo->Multiply(vjet_nlo_kfact);

  //set error to zero
  for(int iBin = 0; iBin < vjet_nlo->GetNbinsX()+1; iBin++)
    vjet_nlo->SetBinError(iBin+1,0.);

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
      
      cout<<"Tree CMS with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<" events no boson "<<eventsNoBoson/tree->GetEntries()<<endl;
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
      
      cout<<"Tree CMS with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<" events no boson "<<eventsNoBoson/tree->GetEntries()<<endl;
      sumwgt_vjet_cms.push_back(sumwgt);
      ifile++;
    }
  }

  /// start making selections                                                                                                                                                                        
  ifile = 0;
  for(auto tree: tree_vjet_cms){

    TTreeReader reader (tree);    
    if(not isPhillTree){

      TTreeReaderValue<float> xsec  (reader,"xsec");
      TTreeReaderValue<float> wgt   (reader,"wgt");
      TTreeReaderValue<int>   wzid_lhe  (reader,"mvid");
      TTreeReaderValue<int>   wzid_ps   (reader,"wzid");
      TTreeReaderValue<int>   l1id  (reader,"l1id");
      TTreeReaderValue<int>   l2id  (reader,"l2id");
      TTreeReaderValue<float> wzmass_lhe  (reader,"wzmass");
      TTreeReaderValue<float> wzmass_ps   (reader,"mvmass");
      TTreeReaderValue<float> wzeta_lhe  (reader,"wzeta");
      TTreeReaderValue<float> wzeta_ps   (reader,"mveta");
      TTreeReaderValue<float> wzphi_lhe  (reader,"wzphi");
      TTreeReaderValue<float> wzphi_ps   (reader,"mvphi");
      TTreeReaderValue<float> wzpt_ps    (reader,"wzpt");
      TTreeReaderValue<float> wzpt_lhe   (reader,"mvpt");
      
      cout<<"Loop on V+jets file "<<file_vjet_cms.at(ifile)->GetName()<<endl;      
      while(reader.Next()){     
	
	// filter away bad events with no matching                                                                                                                                                 
	if((sample == Sample::zvv or sample == Sample::zll) and fabs(*wzid_lhe) != 23){
	  eventsNoBoson++;
	  continue;
	}
	else if(sample == Sample::wjet and fabs(*wzid_lhe) != 24){
	  eventsNoBoson++;
	  continue;
	}
	else if(sample == Sample::gam and fabs(*wzid_lhe) != 22){
	  eventsNoBoson++;
	  continue;
	}
      
	// min-boson pt
	if(*wzpt_lhe < 30) continue;
	//remove tau events
	if(fabs(*l1id) == 15) continue;
	if((sample == Sample::zvv or sample == Sample::zll) and *wzmass_lhe < 30) continue;
	//////////////////                                                                                                                                                                     
	vjet_nlo_cms_lhe->Fill(*wzpt_lhe,(*wgt)*(*xsec)*0.001*scale/sumwgt_vjet_cms.at(ifile)); // cross section in pb
	vjet_nlo_cms_ps->Fill(*wzpt_ps,(*wgt)*(*xsec)*0.001*scale/sumwgt_vjet_cms.at(ifile)); // cross section in pb
	difference->Fill(*wzpt_ps-*wzpt_lhe,(*wgt)*(*xsec)*0.001*scale/sumwgt_vjet_cms.at(ifile));
      }
      ifile++;
    }
    else{

      TTreeReaderArray<float> xsec  (reader,"xs.fXS");
      TTreeReaderArray<float> wgt   (reader,"evtweight.fEvtWeight");
      TTreeReaderArray<int>   wzid_lhe  (reader,"vme_id.fVMEId");
      TTreeReaderArray<int>   wzid_ps   (reader,"v_id.fVId");
      TTreeReaderArray<int>   l1id  (reader,"id_1.fId1");
      TTreeReaderArray<int>   l2id  (reader,"id_2.fId2");
      TTreeReaderArray<float> wzmass_lhe (reader,"vme_m.fVMEM");
      TTreeReaderArray<float> wzmass_ps  (reader,"v_m.fVM");
      TTreeReaderArray<float> wzpt_lhe  (reader,"vme_pt.fVMEPt");
      TTreeReaderArray<float> wzpt_ps   (reader,"v_pt.fVPt");
      TTreeReaderArray<float> wzeta_lhe (reader,"vme_y.fVMEEta");
      TTreeReaderArray<float> wzeta_ps  (reader,"v_y.fVEta");
      TTreeReaderArray<float> wzphi_lhe (reader,"vme_phi.fVMEPhi");
      TTreeReaderArray<float> wzphi_ps  (reader,"v_phi.fVPhi");
      TTreeReaderArray<float> gamiso_04  (reader,"vme_iso04.fVMEIso");
      TTreeReaderArray<float> gamiso_dy  (reader,"vme_dyniso.fDVMEIso");

      cout<<"Loop on V+jets file "<<file_vjet_cms.at(ifile)->GetName()<<endl;      
      while(reader.Next()){     
	
	// filter away bad events with no matching                                                                                                                                                 
	if((sample == Sample::zvv or sample == Sample::zll) and fabs(wzid_lhe[0]) != 23){
	  eventsNoBoson++;
	  continue;
	}
	else if(sample == Sample::wjet and fabs(wzid_lhe[0]) != 24){
	  eventsNoBoson++;
	  continue;
	}
	else if(sample == Sample::gam and fabs(wzid_lhe[0]) != 22){
	  eventsNoBoson++;
	  continue;
	}
      
	// min-boson pt
	if(wzpt_lhe[0] < 30) continue;
	//remove tau events
	if(fabs(l1id[0]) == 15 or fabs(l2id[0]) == 15) continue;
	
	if((sample == Sample::zvv or sample == Sample::zll) and wzmass_lhe[0] < 30) continue;

	if(sample == Sample::gam and gamiso_04[0] > 0.1*wzpt_lhe[0]) continue;
	//if(sample == Sample::gam and gamiso_dy[0] > 0.1*wzpt_lhe[0]) continue;

	//////////////////                                                                                                                                                                     
	vjet_nlo_cms_ps->Fill(wzpt_ps[0],(wgt[0])*(xsec[0])*scale/sumwgt_vjet_cms.at(ifile)); // cross section in pb
	vjet_nlo_cms_lhe->Fill(wzpt_lhe[0],(wgt[0])*(xsec[0])*scale/sumwgt_vjet_cms.at(ifile)); // cross section in pb
	difference->Fill(wzpt_ps[0]-wzpt_lhe[0],(wgt[0])*(xsec[0])*scale/sumwgt_vjet_cms.at(ifile));
      }
      ifile++;
    }
  }

  vjet_nlo_cms_lhe->Scale(1,"width");
  vjet_nlo_cms_ps->Scale(1,"width");

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

}
