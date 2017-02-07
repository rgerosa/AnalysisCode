#include "../CMS_lumi.h"

vector<float> metBin    = {200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1250};

void drawPlot(TCanvas* canvas, TH1F* centralValue_zjet, TH1F* maxValue_zjet, TH1F* minValue_zjet, string outputDIR, string label, string postfix, bool isLog = true){
  
  canvas->cd();

  centralValue_zjet->SetMarkerColor(kBlack);
  centralValue_zjet->SetMarkerStyle(20);
  centralValue_zjet->SetMarkerSize(1);

  if(not isLog)
    centralValue_zjet->GetYaxis()->SetRangeUser(-0.1,0.1);
  else
    centralValue_zjet->GetYaxis()->SetRangeUser(minValue_zjet->GetMinimum()*0.01,maxValue_zjet->GetMaximum()*100);

  
  if(isLog)
    centralValue_zjet->GetYaxis()->SetTitle("Events");  
  else
    centralValue_zjet->GetYaxis()->SetTitle("Z/W");  
  centralValue_zjet->GetYaxis()->SetTitleOffset(1.15);  
  centralValue_zjet->GetXaxis()->SetTitle("Boson p_{T} [GeV]");

  centralValue_zjet->Draw("P");
  maxValue_zjet->SetLineColor(kRed);
  maxValue_zjet->SetLineWidth(2);
  maxValue_zjet->Draw("hist same");
  minValue_zjet->SetLineColor(kBlue);
  minValue_zjet->SetLineWidth(2);
  minValue_zjet->Draw("hist same");

  TLatex* latex = new TLatex ();
  latex->SetNDC();
  latex->SetTextSize(0.6*gPad->GetTopMargin());
  latex->SetTextFont(62);
  latex->SetTextAlign(11);
  latex->DrawLatex(0.2,0.95,Form("CMS Simulation %s",label.c_str()));

  if(isLog) canvas->SetLogy();
  else canvas->SetLogy(0);

  canvas->SaveAs((outputDIR+"/"+postfix+"_qcdscale.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+postfix+"_qcdscale.pdf").c_str(),"pdf");

}


////
void makeQCDScaleUncertainty(string inputDIR_zjet, string inputDIR_wjet, string outputDIR, bool isZgamma = false){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  string postfix;
  if(isZgamma) postfix = "gam";
  else postfix = "wjet";

  float scale_wjet = 1;
  float scale_zjet = 3;

  system(("mkdir -p "+outputDIR).c_str());
  vector<TTree*> tree_zjet;
  vector<TTree*> tree_wjet;
  vector<TFile*> file_zjet;
  vector<TFile*> file_wjet;

  system(("ls "+inputDIR_zjet+" | grep root > list.txt").c_str());
  ifstream infile_zjet ("list.txt");
  string line;
  if(infile_zjet.is_open()){
    while(!infile_zjet.eof()){
      getline(infile_zjet,line);
      if(not TString(line).Contains("root") or line == "") continue;
      file_zjet.push_back(TFile::Open((inputDIR_zjet+"/"+line).c_str(),"READ"));
      tree_zjet.push_back((TTree*) file_zjet.back()->FindObjectAny("tree"));
    }
  }

  infile_zjet.close();

  system(("ls "+inputDIR_wjet+" | grep root > list.txt").c_str());
  ifstream infile_wjet ("list.txt");
  if(infile_wjet.is_open()){
    while(!infile_wjet.eof()){
      getline(infile_wjet,line);
      if(not TString(line).Contains("root") or line == "") continue;
      file_wjet.push_back(TFile::Open((inputDIR_wjet+"/"+line).c_str(),"READ"));
      tree_wjet.push_back((TTree*) file_wjet.back()->FindObjectAny("tree"));
    }
  }

  infile_wjet.close();
  system("rm list.txt");
  
  cout<<"Calculate sumwgt for Z-jets "<<endl;
  // calculate sumwgt                                                                                                                                                                            
  vector<double> sumwgt_zjet;
  int ifile = 0;
  for(auto tree : tree_zjet){
    TTreeReader reader (tree);
    TTreeReaderValue<float> xsec  (reader,"xsec");
    TTreeReaderValue<float> wgt   (reader,"wgt");
    TTreeReaderValue<int>   wzid  (reader,"wzid");

    cout<<"Calculate sumwgt for LO file "<<file_zjet.at(ifile)->GetName()<<endl;
    double sumwgt = 0;
    while(reader.Next()){
      // filter away bad events                                                                                                                                                                       
      if(fabs(*wzid) != 23) continue;
      sumwgt += *wgt;
    }
    cout<<"Tree LO with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<endl;
    sumwgt_zjet.push_back(sumwgt);
    ifile++;
  }

  cout<<"Calculate sumwgt for W-jets "<<endl;
  // calculate sumwgt                                                                                                                                                                            
  vector<double> sumwgt_wjet;
  ifile = 0;
  for(auto tree : tree_wjet){
    TTreeReader reader (tree);
    TTreeReaderValue<float> xsec  (reader,"xsec");
    TTreeReaderValue<float> wgt   (reader,"wgt");
    TTreeReaderValue<int>   wzid  (reader,"wzid");

    cout<<"Calculate sumwgt for LO file "<<file_wjet.at(ifile)->GetName()<<endl;
    double sumwgt = 0;
   
    while(reader.Next()){
      // filter away bad events                                                                                                                                                                       
      if( not isZgamma and fabs(*wzid) != 24) continue;
      else if(isZgamma and fabs(*wzid) != 22) continue;
      sumwgt += *wgt;
    }
    cout<<"Tree LO with entries "<<tree->GetEntries()<<" sumwgt "<<sumwgt<<endl;
    sumwgt_wjet.push_back(sumwgt);
    ifile++;
  }


  // loop on the events
  vector<TH1F*> bosonpt_zjet;
  vector<TH1F*> bosonpt_wjet;

  // Loop on Z-jet trees                                                                                                                                                                            
  ifile = 0;
  for(auto tree: tree_zjet){
    
    TTreeReader reader (tree);
    TTreeReaderValue<float> xsec  (reader,"xsec");
    TTreeReaderValue<float> wgt   (reader,"wgt");
    TTreeReaderValue<int>   wzid  (reader,"wzid");
    TTreeReaderValue<int>   l1id  (reader,"l1id");
    TTreeReaderValue<int>   l2id  (reader,"l2id");
    TTreeReaderValue<unsigned int>   njetsinc  (reader,"njetsinc");
    TTreeReaderValue<unsigned int>   njets  (reader,"njets");
    TTreeReaderValue<float> wzmass  (reader,"wzmass");
    TTreeReaderValue<float> wzpt    (reader,"wzpt");
    TTreeReaderValue<float> wzeta   (reader,"wzeta");
    TTreeReaderValue<float> wzphi   (reader,"wzphi");

    TTreeReaderValue<vector<float> > wgtqcd  (reader,"wgtqcd");

    TTreeReaderValue<float> l1pt   (reader,"l1pt");
    TTreeReaderValue<float> l1eta  (reader,"l1eta");
    TTreeReaderValue<float> l1phi  (reader,"l1phi");
    TTreeReaderValue<float> l2pt   (reader,"l2pt");
    TTreeReaderValue<float> l2eta  (reader,"l2eta");
    TTreeReaderValue<float> l2phi  (reader,"l2phi");

    TTreeReaderValue<vector<float> > jetpt  (reader,"jetpt");
    TTreeReaderValue<vector<float> > jeteta (reader,"jeteta");
    TTreeReaderValue<vector<float> > jetphi (reader,"jetphi");
    TTreeReaderValue<vector<float> > jetmass (reader,"jetmass");


    cout<<"Loop on Z-jet file "<<file_zjet.at(ifile)->GetName()<<endl;
    
    while(reader.Next()){
      if(ifile == 0 and bosonpt_zjet.size() == 0){
	for(size_t iqcd = 0; iqcd < wgtqcd->size(); iqcd++){
	  bosonpt_zjet.push_back(new TH1F(Form("bosonpt_zjet_iqcd_%d",int(iqcd)),"",metBin.size()-1,&metBin[0]));
	  bosonpt_zjet.back()->Sumw2();	  
	}	
      }

      // filter away bad events                                                                                                                                                                    
      if(fabs(*wzid) != 23) continue;
 
      vector<TLorentzVector> jets;
      for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
        TLorentzVector jet4V; jet4V.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetmass->at(ijet));
	float dphi1 = fabs(jetphi->at(ijet)-*l1phi);
	float dphi2 = fabs(jetphi->at(ijet)-*l2phi);
	if(dphi1 > TMath::Pi())
	  dphi1 = 2*TMath::Pi()-dphi1;
	if(dphi2 > TMath::Pi())
	  dphi2 = 2*TMath::Pi()-dphi2;
	if(sqrt(dphi1*dphi1+fabs(jeteta->at(ijet)-*l1eta)*fabs(jeteta->at(ijet)-*l1eta)) > 0.4 and
	   sqrt(dphi2*dphi2+fabs(jeteta->at(ijet)-*l2eta)*fabs(jeteta->at(ijet)-*l2eta)) > 0.4
	   ){ // check cleaning                                                                                                                                                                     
	  jets.push_back(jet4V);
	}
      }
      
      if(jets.size() < 1) continue;            
      if(*wzpt < 200) continue;
      if(jets.at(0).Pt() < 100) continue;
      if(fabs(jets.at(0).Eta()) > 2.5) continue;      

      for(size_t iqcd = 0; iqcd < wgtqcd->size(); iqcd++){
	bosonpt_zjet.at(iqcd)->Fill(*wzpt,*xsec*(*wgt)*scale_zjet*(wgtqcd->at(iqcd)/(*wgt))/sumwgt_zjet.at(ifile));
	
      }
    }
    ifile++;
  }


  // Loop on W-jet trees                                                                                                                                                                            
  ifile = 0;
  for(auto tree: tree_wjet){

    TTreeReader reader (tree);
    TTreeReaderValue<float> xsec  (reader,"xsec");
    TTreeReaderValue<float> wgt   (reader,"wgt");
    TTreeReaderValue<int>   wzid  (reader,"wzid");
    TTreeReaderValue<int>   l1id  (reader,"l1id");
    TTreeReaderValue<int>   l2id  (reader,"l2id");
    TTreeReaderValue<unsigned int>   njetsinc  (reader,"njetsinc");
    TTreeReaderValue<unsigned int>   njets  (reader,"njets");
    TTreeReaderValue<float> wzmass  (reader,"wzmass");
    TTreeReaderValue<float> wzpt    (reader,"wzpt");
    TTreeReaderValue<float> wzeta   (reader,"wzeta");
    TTreeReaderValue<float> wzphi   (reader,"wzphi");
    TTreeReaderValue<vector<float> > wgtqcd  (reader,"wgtqcd");

    TTreeReaderValue<float> l1pt   (reader,"l1pt");
    TTreeReaderValue<float> l1eta  (reader,"l1eta");
    TTreeReaderValue<float> l1phi  (reader,"l1phi");
    TTreeReaderValue<float> l2pt   (reader,"l2pt");
    TTreeReaderValue<float> l2eta  (reader,"l2eta");
    TTreeReaderValue<float> l2phi  (reader,"l2phi");

    TTreeReaderValue<vector<float> > jetpt  (reader,"jetpt");
    TTreeReaderValue<vector<float> > jeteta (reader,"jeteta");
    TTreeReaderValue<vector<float> > jetphi (reader,"jetphi");
    TTreeReaderValue<vector<float> > jetmass (reader,"jetmass");

    cout<<"Loop on W-jet file "<<file_wjet.at(ifile)->GetName()<<endl;
    
    int ientry = 0;
    while(reader.Next()){
      if(ifile == 0 and bosonpt_wjet.size() == 0){
	for(size_t iqcd = 0; iqcd < wgtqcd->size(); iqcd++){
	  bosonpt_wjet.push_back(new TH1F(Form("bosonpt_%s_iqcd_%d",postfix.c_str(),int(iqcd)),"",metBin.size()-1,&metBin[0]));
	  bosonpt_wjet.back()->Sumw2();	  
	}	
      }

      // filter away bad events                                                                                                                                                                    
      if(not isZgamma and fabs(*wzid) != 24) continue;
      else if(isZgamma and fabs(*wzid) != 22) continue;
 
      vector<TLorentzVector> jets;
      for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
        TLorentzVector jet4V; jet4V.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetmass->at(ijet));
	float dphi1 = fabs(jetphi->at(ijet)-*l1phi);
	float dphi2 = fabs(jetphi->at(ijet)-*l2phi);
	if(dphi1 > TMath::Pi())
	  dphi1 = 2*TMath::Pi()-dphi1;
	if(dphi2 > TMath::Pi())
            dphi2 = 2*TMath::Pi()-dphi2;
	if(sqrt(dphi1*dphi1+fabs(jeteta->at(ijet)-*l1eta)*fabs(jeteta->at(ijet)-*l1eta)) > 0.4 and
	   sqrt(dphi2*dphi2+fabs(jeteta->at(ijet)-*l2eta)*fabs(jeteta->at(ijet)-*l2eta)) > 0.4
	   ){ // check cleaning                                                                                                                                                                     
	  jets.push_back(jet4V);
	}
      }

      if(jets.size() < 1) continue;
      
      if(*wzpt < 200) continue;
      if(jets.at(0).Pt() < 100) continue;
      if(fabs(jets.at(0).Eta()) > 2.5) continue;      

      for(size_t iqcd = 0; iqcd < wgtqcd->size(); iqcd++){
	bosonpt_wjet.at(iqcd)->Fill(*wzpt,*xsec*(*wgt)*scale_wjet*(wgtqcd->at(iqcd)/(*wgt))/sumwgt_wjet.at(ifile));	
      }
    }
    ifile++;
  }

  // Plot results
  TFile* outputFile = new TFile((outputDIR+"/outputHistograms_qcd.root").c_str(),"RECREATE");
  outputFile->cd();
  for(auto ihist : bosonpt_zjet)
    ihist->Write();
  for(auto ihist : bosonpt_wjet)
    ihist->Write();

  //////////
  TCanvas* canvas =  new TCanvas("canvas","",600,625);
  canvas->cd();

  TH1F* maxValue_zjet = (TH1F*) bosonpt_zjet.at(0)->Clone("maxValue_zjet");
  maxValue_zjet->Reset();
  TH1F* maxValue_wjet = (TH1F*) bosonpt_wjet.at(0)->Clone(Form("maxValue_%s",postfix.c_str()));
  maxValue_zjet->Reset();
  TH1F* minValue_zjet = (TH1F*) bosonpt_zjet.at(0)->Clone("minValue_zjet");
  maxValue_zjet->Reset();
  TH1F* minValue_wjet = (TH1F*) bosonpt_wjet.at(0)->Clone(Form("maxValue_%s",postfix.c_str()));
  maxValue_zjet->Reset();
  
  for(size_t ihist = 0; ihist < bosonpt_zjet.size(); ihist++){
    for(int iBin = 0; iBin < bosonpt_zjet.at(ihist)->GetNbinsX()+1; iBin++){
      if(maxValue_zjet->GetBinContent(iBin+1) < bosonpt_zjet.at(ihist)->GetBinContent(iBin+1))
        maxValue_zjet->SetBinContent(iBin+1,bosonpt_zjet.at(ihist)->GetBinContent(iBin+1));
      if(minValue_zjet->GetBinContent(iBin+1) > bosonpt_zjet.at(ihist)->GetBinContent(iBin+1))
        minValue_zjet->SetBinContent(iBin+1,bosonpt_zjet.at(ihist)->GetBinContent(iBin+1));
    }
  }

  for(size_t ihist = 0; ihist < bosonpt_wjet.size(); ihist++){
    for(int iBin = 0; iBin < bosonpt_wjet.at(ihist)->GetNbinsX()+1; iBin++){
      if(maxValue_wjet->GetBinContent(iBin+1) < bosonpt_wjet.at(ihist)->GetBinContent(iBin+1))
        maxValue_wjet->SetBinContent(iBin+1,bosonpt_wjet.at(ihist)->GetBinContent(iBin+1));
      if(minValue_wjet->GetBinContent(iBin+1) > bosonpt_wjet.at(ihist)->GetBinContent(iBin+1))
        minValue_wjet->SetBinContent(iBin+1,bosonpt_wjet.at(ihist)->GetBinContent(iBin+1));
    }
  }
  
  vector<TH1F*> ratios;
  for(size_t ihisto = 0; ihisto < bosonpt_zjet.size(); ihisto++){
    if(isZgamma)
      ratios.push_back((TH1F*) bosonpt_zjet.at(ihisto)->Clone(Form("zog_ratio_bin_%d",int(ihisto))));
    else
      ratios.push_back((TH1F*) bosonpt_zjet.at(ihisto)->Clone(Form("zow_ratio_bin_%d",int(ihisto))));

    ratios.back()->Divide((TH1F*) bosonpt_wjet.at(ihisto));
  }

  // make plots
  drawPlot(canvas,bosonpt_zjet.at(0),maxValue_zjet,minValue_zjet,outputDIR,"Rate","zjets",true); // Z/W;
  if(not isZgamma)
    drawPlot(canvas,bosonpt_wjet.at(0),maxValue_wjet,minValue_wjet,outputDIR,"Rate","wjets",true); // Z/W;
  else
    drawPlot(canvas,bosonpt_wjet.at(0),maxValue_wjet,minValue_wjet,outputDIR,"Rate","gjets",true); // Z/W;

  maxValue_zjet->Divide(bosonpt_zjet.at(0));
  maxValue_wjet->Divide(bosonpt_wjet.at(0));
  minValue_zjet->Divide(bosonpt_zjet.at(0));
  minValue_wjet->Divide(bosonpt_wjet.at(0));

  TF1* maxValue_zjet_func = new TF1("maxValue_zjet_func","pol3",metBin.front(),metBin.back());
  TF1* minValue_zjet_func = new TF1("minValue_zjet_func","pol3",metBin.front(),metBin.back());
  TF1* maxValue_wjet_func = new TF1(Form("maxValue_%s_func",postfix.c_str()),"pol3",metBin.front(),metBin.back());
  TF1* minValue_wjet_func = new TF1(Form("minValue_%s_func",postfix.c_str()),"pol3",metBin.front(),metBin.back());

  maxValue_zjet->Fit(maxValue_zjet_func,"RME0");
  minValue_zjet->Fit(minValue_zjet_func,"RME0");
  maxValue_wjet->Fit(maxValue_wjet_func,"RME0");
  minValue_wjet->Fit(minValue_wjet_func,"RME0");

  maxValue_zjet_func->SetLineColor(kRed);
  maxValue_zjet_func->SetLineStyle(2);

  minValue_zjet_func->SetLineColor(kRed);
  minValue_zjet_func->SetLineStyle(2);

  maxValue_wjet_func->SetLineColor(kBlue);
  maxValue_wjet_func->SetLineStyle(2);

  minValue_wjet_func->SetLineColor(kBlue);
  minValue_wjet_func->SetLineStyle(2);

  maxValue_zjet_func->Write();
  minValue_zjet_func->Write();
  maxValue_wjet_func->Write();
  minValue_wjet_func->Write();

  TH1F* ratio_qcd_up = (TH1F*) ratios.at(1)->Clone("ratio_qcd_up");
  ratio_qcd_up->Reset();
  TH1F* ratio_qcd_dw = (TH1F*) ratios.at(1)->Clone("ratio_qcd_dw");
  ratio_qcd_dw->Reset();
  TH1F* ratio_qcd_sym_up = (TH1F*) ratios.at(1)->Clone("ratio_qcd_sym_up");
  ratio_qcd_sym_up->Reset();
  TH1F* ratio_qcd_sym_dw = (TH1F*) ratios.at(1)->Clone("ratio_qcd_sym_dw");
  ratio_qcd_sym_dw->Reset();

  for(int iBin = 0; iBin < ratio_qcd_up->GetNbinsX()+1; iBin++){
    ratio_qcd_up->SetBinContent(iBin+1,maxValue_zjet->GetBinContent(iBin+1)/maxValue_wjet->GetBinContent(iBin+1));
    ratio_qcd_dw->SetBinContent(iBin+1,minValue_zjet->GetBinContent(iBin+1)/minValue_wjet->GetBinContent(iBin+1));
  }

  for(int iBin = 0; iBin < ratio_qcd_up->GetNbinsX()+1; iBin++){
    ratio_qcd_sym_up->SetBinContent(iBin+1,(ratio_qcd_up->GetBinContent(iBin+1)-ratio_qcd_dw->GetBinContent(iBin+1))/2);
    ratio_qcd_sym_dw->SetBinContent(iBin+1,-(ratio_qcd_up->GetBinContent(iBin+1)-ratio_qcd_dw->GetBinContent(iBin+1))/2);
  }

  ratio_qcd_up->Write();
  ratio_qcd_dw->Write();
  ratio_qcd_sym_up->Write();
  ratio_qcd_sym_dw->Write();

  //
  canvas->SetLogy(0);
  ratio_qcd_sym_up->GetXaxis()->SetTitle("Boson p_T [GeV]");
  ratio_qcd_sym_up->GetYaxis()->SetTitle("QCD uncertainty");
  ratio_qcd_sym_up->SetLineColor(kRed);
  ratio_qcd_sym_up->SetLineWidth(2);
  ratio_qcd_sym_up->Draw("hist");
  if(not isZgamma)
    ratio_qcd_sym_up->GetYaxis()->SetRangeUser(-0.02,0.02);
  else
    ratio_qcd_sym_up->GetYaxis()->SetRangeUser(-0.10,0.10);
  ratio_qcd_sym_dw->GetXaxis()->SetTitle("Boson p_T [GeV]");
  ratio_qcd_sym_dw->GetYaxis()->SetTitle("QCD qcd uncertainty");
  ratio_qcd_sym_dw->SetLineColor(kBlue);
  ratio_qcd_sym_dw->SetLineWidth(2);
  ratio_qcd_sym_dw->Draw("hist same");
  canvas->SaveAs((outputDIR+"/qcd_uncertainty_final.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/qcd_uncertainty_final.png").c_str(),"png");

  
  // make QCD shape uncertainty:
  TF1* line_zjet = new TF1("line_zjet","[0]+[1]*x",metBin.front(),metBin.back());
  line_zjet->SetParameter(1,(minValue_zjet_func->Eval(metBin.back())-maxValue_zjet_func->Eval(metBin.front()))/(metBin.back()-metBin.front()));
  line_zjet->SetParameter(0,maxValue_zjet_func->Eval(metBin.front())-line_zjet->GetParameter(1)*metBin.front());

  TF1* line_wjet = new TF1("line_wjet","[0]+[1]*x",metBin.front(),metBin.back());
  line_wjet->SetParameter(1,(minValue_wjet_func->Eval(metBin.back())-maxValue_wjet_func->Eval(metBin.front()))/(metBin.back()-metBin.front()));
  line_wjet->SetParameter(0,maxValue_wjet_func->Eval(metBin.front())-line_wjet->GetParameter(1)*metBin.front());

  TF1* line_zjet_2 = new TF1("line_zjet_2","[0]+[1]*x",metBin.front(),metBin.back());
  line_zjet_2->SetParameter(1,(maxValue_zjet_func->Eval(metBin.back())-minValue_zjet_func->Eval(metBin.front()))/(metBin.back()-metBin.front()));
  line_zjet_2->SetParameter(0,minValue_zjet_func->Eval(metBin.front())-line_zjet_2->GetParameter(1)*metBin.front());

  TF1* line_wjet_2 = new TF1("line_wjet_2","[0]+[1]*x",metBin.front(),metBin.back());
  line_wjet_2->SetParameter(1,(maxValue_wjet_func->Eval(metBin.back())-minValue_wjet_func->Eval(metBin.front()))/(metBin.back()-metBin.front()));
  line_wjet_2->SetParameter(0,minValue_wjet_func->Eval(metBin.front())-line_wjet_2->GetParameter(1)*metBin.front());

  minValue_zjet->Write();
  maxValue_zjet->Write();
  minValue_wjet->Write();
  maxValue_wjet->Write();
  line_zjet->Write();
  line_wjet->Write();
  line_zjet_2->Write();
  line_wjet_2->Write();

  TH1F* ratio_shape_up = (TH1F*) ratios.at(1)->Clone("ratio_shape_up");
  ratio_shape_up->Reset();
  TH1F* ratio_shape_dw = (TH1F*) ratios.at(1)->Clone("ratio_shape_dw");
  ratio_shape_dw->Reset();
  TH1F* ratio_shape_sym_up = (TH1F*) ratios.at(1)->Clone("ratio_shape_sym_up");
  ratio_shape_sym_up->Reset();
  TH1F* ratio_shape_sym_dw = (TH1F*) ratios.at(1)->Clone("ratio_shape_sym_dw");
  ratio_shape_sym_dw->Reset();
  for(int iBin = 0; iBin < ratio_shape_up->GetNbinsX()+1; iBin++){
    ratio_shape_up->SetBinContent(iBin+1,line_zjet->Eval(ratio_shape_up->GetBinCenter(iBin+1))/line_wjet->Eval(ratio_shape_up->GetBinCenter(iBin+1)));
    ratio_shape_dw->SetBinContent(iBin+1,line_zjet_2->Eval(ratio_shape_dw->GetBinCenter(iBin+1))/line_wjet_2->Eval(ratio_shape_dw->GetBinCenter(iBin+1)));
    ratio_shape_sym_up->SetBinContent(iBin+1,(ratio_shape_up->GetBinContent(iBin+1)-ratio_shape_dw->GetBinContent(iBin+1))/2.);
    ratio_shape_sym_dw->SetBinContent(iBin+1,-((ratio_shape_up->GetBinContent(iBin+1)-ratio_shape_dw->GetBinContent(iBin+1))/2.));
  }
  ratio_shape_up->Write();
  ratio_shape_dw->Write();
  ratio_shape_sym_up->Write();
  ratio_shape_sym_dw->Write();

  // Last plot for the shape uncertainty
  canvas->cd();
  canvas->SetLogy(0);
  maxValue_zjet->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
  maxValue_zjet->GetYaxis()->SetTitle("Uncertainty");
  maxValue_zjet->SetLineColor(kRed);
  maxValue_zjet->SetMarkerStyle(20);
  maxValue_zjet->SetMarkerColor(kRed);
  maxValue_zjet->SetLineWidth(2);
  maxValue_zjet->Draw("hist");
  maxValue_zjet->Draw("Psame");
  minValue_zjet->SetLineColor(kRed);
  minValue_zjet->SetLineWidth(2);
  minValue_zjet->SetMarkerStyle(20);
  minValue_zjet->SetMarkerColor(kRed);
  minValue_zjet->Draw("hist same");
  minValue_zjet->Draw("P same");

  line_zjet->SetLineColor(kRed);
  line_zjet->SetLineWidth(2);
  line_zjet->SetLineStyle(7);
  line_zjet->Draw("L same");

  line_zjet_2->SetLineColor(kRed);
  line_zjet_2->SetLineWidth(2);
  line_zjet_2->SetLineStyle(7);
  line_zjet_2->Draw("L same");

  maxValue_wjet->SetLineColor(kBlue);
  maxValue_wjet->SetMarkerStyle(20);
  maxValue_wjet->SetMarkerColor(kBlue);
  maxValue_wjet->SetLineWidth(2);
  maxValue_wjet->Draw("hist same");
  maxValue_wjet->Draw("P same");
  minValue_wjet->SetLineColor(kBlue);
  minValue_wjet->SetLineWidth(2);
  minValue_wjet->SetMarkerStyle(20);
  minValue_wjet->SetMarkerColor(kBlue);
  minValue_wjet->Draw("hist same");
  minValue_wjet->Draw("P same");

  line_wjet->SetLineColor(kBlue);
  line_wjet->SetLineWidth(2);
  line_wjet->SetLineStyle(7);
  line_wjet->Draw("L same");

  line_wjet_2->SetLineColor(kBlue);
  line_wjet_2->SetLineWidth(2);
  line_wjet_2->SetLineStyle(7);
  line_wjet_2->Draw("L same");

  maxValue_zjet->GetYaxis()->SetRangeUser(min(minValue_zjet->GetMinimum(),minValue_wjet->GetMinimum())*0.9,max(maxValue_zjet->GetMaximum(),maxValue_wjet->GetMaximum())*1.3);
  
  TLatex* latex = new TLatex ();
  latex->SetNDC();
  latex->SetTextSize(0.6*gPad->GetTopMargin());
  latex->SetTextFont(62);
  latex->SetTextAlign(11);
  latex->DrawLatex(0.2,0.95,Form("CMS Simulation"));

  TLegend leg (0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(minValue_zjet,"Z+jets QCD-scale Unc.","LP");
  leg.AddEntry(line_zjet,"Z+jets QCD-scale Shape","L");
  if(isZgamma){
    leg.AddEntry(minValue_wjet,"#gamma+jets QCD-scale Unc.","LP");
    leg.AddEntry(line_wjet,"#gamma+jets QCD-scale Shape","L");
  }
  else{
    leg.AddEntry(minValue_wjet,"W+jets QCD-scale Unc.","LP");
    leg.AddEntry(line_wjet,"W+jets QCD-scale Shape","L");
  }
  leg.Draw("same");

  canvas->SaveAs((outputDIR+"/qcd_shape_uncertainty.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/qcd_shape_uncertainty.png").c_str(),"png");

  ////
  ratio_shape_sym_up->GetXaxis()->SetTitle("Boson p_T [GeV]");
  ratio_shape_sym_up->GetYaxis()->SetTitle("QCD shape uncertainty");
  ratio_shape_sym_up->SetLineColor(kRed);
  ratio_shape_sym_up->SetLineWidth(2);
  ratio_shape_sym_up->Draw("hist");
  if(not isZgamma)
    ratio_shape_sym_up->GetYaxis()->SetRangeUser(-0.02,0.02);
  else
    ratio_shape_sym_up->GetYaxis()->SetRangeUser(-0.10,0.10);
  ratio_shape_sym_dw->GetXaxis()->SetTitle("Boson p_T [GeV]");
  ratio_shape_sym_dw->GetYaxis()->SetTitle("QCD shape uncertainty");
  ratio_shape_sym_dw->SetLineColor(kBlue);
  ratio_shape_sym_dw->SetLineWidth(2);
  ratio_shape_sym_dw->Draw("hist same");
  canvas->SaveAs((outputDIR+"/qcd_shape_uncertainty_final.pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/qcd_shape_uncertainty_final.png").c_str(),"png");

    
  for(int iBin = 0; iBin < maxValue_zjet->GetNbinsX()+1; iBin++){
    maxValue_zjet->SetBinContent(iBin+1,maxValue_zjet->GetBinContent(iBin+1)-1);
    minValue_zjet->SetBinContent(iBin+1,minValue_zjet->GetBinContent(iBin+1)-1);
    maxValue_wjet->SetBinContent(iBin+1,maxValue_wjet->GetBinContent(iBin+1)-1);
    minValue_wjet->SetBinContent(iBin+1,minValue_wjet->GetBinContent(iBin+1)-1);
  }
  
  TH1F* maxValue_ratio = (TH1F*) ratios.at(0)->Clone("maxValue_ratio");
  TH1F* minValue_ratio = (TH1F*) ratios.at(0)->Clone("minValue_ratio");

  for(size_t ihist = 1; ihist < ratios.size(); ihist++){      
    for(int iBin = 0; iBin < ratios.at(0)->GetNbinsX()+1; iBin++){
      ratios.at(ihist)->SetBinContent(iBin+1,(ratios.at(ihist)->GetBinContent(iBin+1)/ratios.at(0)->GetBinContent(iBin+1))-1);
    }
  }

  for(int iBin = 0; iBin < ratios.at(0)->GetNbinsX()+1; iBin++){
    ratios.at(0)->SetBinContent(iBin+1,0.);
  }

  for(auto hist: ratios)
    hist->Write();
    
      
  for(size_t ihist = 1; ihist < ratios.size(); ihist++){
    for(int iBin = 0; iBin < ratios.at(ihist)->GetNbinsX()+1; iBin++){
      if(maxValue_ratio->GetBinContent(iBin+1) < ratios.at(ihist)->GetBinContent(iBin+1))
	maxValue_ratio->SetBinContent(iBin+1,ratios.at(ihist)->GetBinContent(iBin+1));
      if(minValue_ratio->GetBinContent(iBin+1) > ratios.at(ihist)->GetBinContent(iBin+1))
	minValue_ratio->SetBinContent(iBin+1,ratios.at(ihist)->GetBinContent(iBin+1));
    }
  }

  if(not isZgamma)
    drawPlot(canvas,ratios.at(0),maxValue_ratio,minValue_ratio,outputDIR,"Z/W","zoverw",false); // Z/W;
  else
    drawPlot(canvas,ratios.at(0),maxValue_ratio,minValue_ratio,outputDIR,"Z/#gamma","zoverg",false); // Z/W;

  outputFile->Close();

  
}
