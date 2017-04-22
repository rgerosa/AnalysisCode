#include "../CMS_lumi.h"

enum class Sample {zvv,wjet,zll,gam};

///////////////////
void makeVJetsNLOComparisonWithTheory(string inputCMSDirectory, string inputTheoryFile, string outputDirectory, Sample sample, bool useLHELevelInfo = false, bool isPhillTree = false){

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

  TH1F* vjet_nlo_cms   = (TH1F*) vjet_nlo->Clone("vjet_nlo_cms");
  vjet_nlo_cms->Reset();

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
      TTreeReaderValue<float> xsec  (reader,"xsec");
      TTreeReaderValue<float> wgt   (reader,"wgt");
      string wzid_b = "wzid";
      if(useLHELevelInfo) wzid_b = "mvid";
      TTreeReaderValue<int>   wzid  (reader,wzid_b.c_str());
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
      TTreeReaderArray<float> xsec  (reader,"xs.fXS");
      TTreeReaderArray<float> wgt   (reader,"evtweight.fEvtWeight");
      TTreeReaderArray<int>   wzid  (reader,"v_id.fVId");

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
      string wzid_b = "wzid";
      if(useLHELevelInfo) wzid_b = "mvid";
      TTreeReaderValue<int>   wzid  (reader,wzid_b.c_str());
      TTreeReaderValue<int>   l1id  (reader,"l1id");
      TTreeReaderValue<int>   l2id  (reader,"l2id");
      string wzmass_b = "wzmass";
      if(useLHELevelInfo) wzmass_b = "mvmass";
      TTreeReaderValue<float> wzmass  (reader,wzmass_b.c_str());
      string wzpt_b = "wzpt";
      if(useLHELevelInfo) wzpt_b = "mvpt";
      TTreeReaderValue<float> wzpt    (reader,wzpt_b.c_str());
      string wzeta_b = "wzeta";
      if(useLHELevelInfo) wzeta_b = "mveta";
      TTreeReaderValue<float> wzeta   (reader,wzeta_b.c_str());
      string wzphi_b = "wzphi";
      if(useLHELevelInfo) wzphi_b = "mvphi";
      TTreeReaderValue<float> wzphi   (reader,wzphi_b.c_str());
      
      cout<<"Loop on V+jets file "<<file_vjet_cms.at(ifile)->GetName()<<endl;      
      while(reader.Next()){     
	
	// filter away bad events with no matching                                                                                                                                                 
	if((sample == Sample::zvv or sample == Sample::zll) and fabs(*wzid) != 23){
	  eventsNoBoson++;
	  continue;
	}
	else if(sample == Sample::wjet and fabs(*wzid) != 24){
	  eventsNoBoson++;
	  continue;
	}
	else if(sample == Sample::gam and fabs(*wzid) != 22){
	  eventsNoBoson++;
	  continue;
	}
      
	// min-boson pt
	if(*wzpt < 30) continue;
	//remove tau events
	if(fabs(*l1id) == 15) continue;
	if((sample == Sample::zvv or sample == Sample::zll) and *wzmass < 30) continue;
	//////////////////                                                                                                                                                                     
	vjet_nlo_cms->Fill(*wzpt,(*wgt)*(*xsec)*0.001*scale/sumwgt_vjet_cms.at(ifile)); // cross section in pb
      }
      ifile++;
    }
    else{

      TTreeReaderArray<float> xsec  (reader,"xs.fXS");
      TTreeReaderArray<float> wgt   (reader,"evtweight.fEvtWeight");
      TTreeReaderArray<int>   wzid  (reader,"v_id.fVId");
      TTreeReaderArray<int>   l1id  (reader,"id_1.fId1");
      TTreeReaderArray<int>   l2id  (reader,"id_2.fId2");
      TTreeReaderArray<float> wzmass  (reader,"v_m.fVM");
      TTreeReaderArray<float> wzpt    (reader,"v_pt.fVPt");
      TTreeReaderArray<float> wzeta   (reader,"v_y.fVEta");
      TTreeReaderArray<float> wzphi   (reader,"v_phi.fVPhi");

      cout<<"Loop on V+jets file "<<file_vjet_cms.at(ifile)->GetName()<<endl;      
      while(reader.Next()){     
	
	// filter away bad events with no matching                                                                                                                                                 
	if((sample == Sample::zvv or sample == Sample::zll) and fabs(wzid[0]) != 23){
	  eventsNoBoson++;
	  continue;
	}
	else if(sample == Sample::wjet and fabs(wzid[0]) != 24){
	  eventsNoBoson++;
	  continue;
	}
	else if(sample == Sample::gam and fabs(wzid[0]) != 22){
	  eventsNoBoson++;
	  continue;
	}
      
	// min-boson pt
	if(wzpt[0] < 30) continue;
	//remove tau events
	if(fabs(l1id[0]) == 15 or fabs(l2id[0]) == 15) continue;
	
	if((sample == Sample::zvv or sample == Sample::zll) and wzmass[0] < 30) continue;
	//////////////////                                                                                                                                                                     
	vjet_nlo_cms->Fill(wzpt[0],(wgt[0])*(xsec[0])*scale/sumwgt_vjet_cms.at(ifile)); // cross section in pb
      }
      ifile++;
    }
  }

  vjet_nlo_cms->Scale(1,"width");

  // make plot
  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->SetBottomMargin(0.3);
  TPad *pad = new TPad("pad","pad",0,0.,1,0.9);
  pad->SetTopMargin(0.7);
  pad->SetRightMargin(0.06);
  pad->SetFillColor(0);
  pad->SetGridy(1);
  pad->SetFillStyle(0);

  canvas->cd();

  vjet_nlo_cms->GetXaxis()->SetTitleSize(0);
  vjet_nlo_cms->GetXaxis()->SetLabelSize(0);
  vjet_nlo_cms->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb/GeV]");
  vjet_nlo_cms->GetYaxis()->SetTitleOffset(1.30);
  vjet_nlo_cms->GetXaxis()->SetRangeUser(200,1400);

  vjet_nlo->SetLineColor(kRed);
  vjet_nlo->SetMarkerColor(kRed);
  vjet_nlo->SetMarkerStyle(20);
  vjet_nlo->SetMarkerSize(1);
  vjet_nlo->SetLineWidth(2);

  vjet_nlo_cms->SetLineColor(kBlack);
  vjet_nlo_cms->SetMarkerColor(kBlack);
  vjet_nlo_cms->SetMarkerStyle(20);
  vjet_nlo_cms->SetMarkerSize(1);
  vjet_nlo_cms->SetLineWidth(2);
  vjet_nlo_cms->Draw("hist");

  TH1F* vjet_nlo_cms_band = (TH1F*) vjet_nlo_cms->Clone("vjet_nlo_cms_band");
  vjet_nlo_cms_band->SetFillColor(kBlack);
  vjet_nlo_cms_band->SetFillStyle(3001);
  vjet_nlo_cms_band->Draw("E2same");
  vjet_nlo_cms->Draw("hist same");

  vjet_nlo->Draw("hist same");

  CMS_lumi(canvas,"");

  TLegend leg (0.7,0.7,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  if(sample == Sample::zvv){
    leg.AddEntry(vjet_nlo,"Z+jets: theorists","PL");
    leg.AddEntry(vjet_nlo_cms,"Z+jets: CMS","PL");
  }
  if(sample == Sample::wjet){
    leg.AddEntry(vjet_nlo,"W+jets: theorists","PL");
    leg.AddEntry(vjet_nlo_cms,"W+jets: CMS","PL");
  }
  if(sample == Sample::zll){
    leg.AddEntry(vjet_nlo,"DY+jets: theorists","PL");
    leg.AddEntry(vjet_nlo_cms,"DY+jets: CMS","PL");
  }
  if(sample == Sample::gam){
    leg.AddEntry(vjet_nlo,"#gamma+jets: theorists","PL");
    leg.AddEntry(vjet_nlo_cms,"#gamma+jets: CMS","PL");
  }
  leg.Draw("same");
  canvas->SetLogy();

  pad->Draw();
  pad->cd();
  
  TH1F* ratio = (TH1F*) vjet_nlo->Clone("ratio");
  ratio->Divide(vjet_nlo_cms);
  ratio->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
  ratio->GetXaxis()->SetTitleOffset(1.15);
  ratio->GetYaxis()->SetTitle("Theory/CMS");
  ratio->GetYaxis()->SetTitleOffset(1.30);
  ratio->GetYaxis()->SetNdivisions(504);
  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(2);
  ratio->GetYaxis()->SetRangeUser(0.8,1.2);
  ratio->GetXaxis()->SetRangeUser(150,1400);
  ratio->Draw("hist");
  
  canvas->cd();
  canvas->RedrawAxis("sameaxis");
  
  if(sample == Sample::zvv){
    canvas->SaveAs((outputDirectory+"/bosonpt_zvv.png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/bosonpt_zvv.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::wjet){
    canvas->SaveAs((outputDirectory+"/bosonpt_wjet.png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/bosonpt_wjet.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::zll){
    canvas->SaveAs((outputDirectory+"/bosonpt_zll.png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/bosonpt_zll.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::gam){
    canvas->SaveAs((outputDirectory+"/bosonpt_gam.png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/bosonpt_gam.pdf").c_str(),"pdf");
  }
}
