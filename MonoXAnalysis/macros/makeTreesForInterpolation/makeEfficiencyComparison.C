#include "../CMS_lumi.h"

void printBin(TH1F* histo){

  cout<<"Bin content for histo: "<<histo->GetName()<<endl;
  for(int iBin = 0; iBin < histo->GetNbinsX()+1; iBin++)
    cout<<"Bin "<<iBin<<" content "<<histo->GetBinContent(iBin+1)<<" error "<<histo->GetBinError(iBin+1)<<endl;

  return;

}

void plot1D(TCanvas* canvas, TH1* histo, string name, string outputDIR, bool scalefactor = false){
  canvas->SetLogy(0);
  histo->GetXaxis()->SetTitle("Recoil [GeV]");
  if(not scalefactor){
    histo->GetZaxis()->SetTitle("a.u.");  
    canvas->SetLogy();
  }
  else
    histo->GetZaxis()->SetTitle("Scale Factor");

  if(scalefactor){
    for(int iBinX = 0; iBinX < histo->GetNbinsX()+1; iBinX++){
      for(int iBinY = 0; iBinY < histo->GetNbinsY()+1; iBinY++){
	if(histo->GetBinContent(iBinX+1,iBinY+1) > 5 or histo->GetBinContent(iBinX+1,iBinY+1) < 0.01)
	  histo->SetBinContent(iBinX+1,iBinY+1,0);
      }
    }
  }

  histo->SetMarkerStyle(20);
  histo->SetMarkerSize(1);
  if(scalefactor){
    histo->GetYaxis()->SetRangeUser(0.5,1.5);
    histo->Fit("pol0");
    gStyle->SetOptFit(1111);
    TF1* funz = (TF1*) histo->GetListOfFunctions()->At(0);
    funz->SetLineColor(kRed);
    funz->SetLineWidth(2);
  }
  histo->Draw("E1PFUNC");
  CMS_lumi(canvas,"2.77",true);
  canvas->SaveAs((outputDIR+"/"+name+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+name+".pdf").c_str(),"pdf");
  
}

void plot2D(TCanvas* canvas, TH2* histo, string name, string outputDIR, bool scalefactor = false){

  canvas->SetLogz(0);
  histo->GetYaxis()->SetTitle("mediator p_{T} [GeV]");
  histo->GetXaxis()->SetTitle("Recoil [GeV]");
  if(not scalefactor){
    histo->GetZaxis()->SetTitle("a.u.");
    canvas->SetLogz();
  }
  else
    histo->GetZaxis()->SetTitle("Scale Factor");

  if(scalefactor){
    for(int iBinX = 0; iBinX < histo->GetNbinsX()+1; iBinX++){
      for(int iBinY = 0; iBinY < histo->GetNbinsY()+1; iBinY++){
	if(histo->GetBinContent(iBinX+1,iBinY+1) > 5 or histo->GetBinContent(iBinX+1,iBinY+1) < 0.01)
	  histo->SetBinContent(iBinX+1,iBinY+1,0);
      }
    }
  }


  if(not scalefactor)
    histo->Draw("colz");
  else
    histo->Draw("colz texte");
  CMS_lumi(canvas,"2.77",true);
  canvas->SaveAs((outputDIR+"/"+name+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+name+".pdf").c_str(),"pdf");

}

void makeEfficiencyComparison(string namefile1, string namefile2, string outputDIR, bool reweight = true, bool useOnlyCommonPoints = false, float minMedMass = 0){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetPaintTextFormat(".2f");
  system(("mkdir -p "+outputDIR).c_str());

  TFile* file1 = TFile::Open(namefile1.c_str());
  TFile* file2 = TFile::Open(namefile2.c_str());

  TTree* tree1 = (TTree*) file1->Get("tree");
  TTree* tree2 = (TTree*) file2->Get("tree");

  // weight properly according to mediator pt:eta
  vector<float> recoilBin_2D = {200,250,350,450,650,1000};
  vector<float> ptBin_2D     = {100,200,300,1000};
  vector<float> recoilBin_1D = {200,300,450,650,850,1250};
  

  TH2F* mediator1_D_2D = new TH2F("mediator1_D_2D","",125,0,1250,125,0,1250);
  TH2F* mediator2_D_2D = new TH2F("mediator2_D_2D","",125,0,1250,125,0,1250);
  mediator1_D_2D->Sumw2();
  mediator2_D_2D->Sumw2();

  TH1F* mediator1_D_1D = new TH1F("mediator1_D_1D","",80,0,1250);
  TH1F* mediator2_D_1D = new TH1F("mediator2_D_1D","",80,0,1250);
  mediator1_D_1D->Sumw2();
  mediator2_D_1D->Sumw2();

  if(reweight){
    tree1->Draw("genMediatorPt:pfMetPt >> mediator1_D_2D",Form("(genWeight*weight)*(genMediatorPt > 0 && genMediatorMass > %f)",minMedMass),"goff");
    tree2->Draw("genMediatorPt:pfMetPt >> mediator2_D_2D",Form("(genWeight*weight)*(genMediatorPt > 0 && genMediatorMass > %f)",minMedMass),"goff");
    tree1->Draw("pfMetPt >> mediator1_D_1D",Form("(genWeight*weight)*(genMediatorPt > 0 && genMediatorMass > %f)",minMedMass),"goff");
    tree2->Draw("pfMetPt >> mediator2_D_1D",Form("(genWeight*weight)*(genMediatorPt > 0 && genMediatorMass > %f)",minMedMass),"goff");
  }
  mediator1_D_2D->Scale(1./mediator1_D_2D->Integral());
  mediator2_D_2D->Scale(2./mediator2_D_2D->Integral());
  mediator1_D_1D->Scale(1./mediator1_D_1D->Integral());
  mediator2_D_1D->Scale(2./mediator2_D_1D->Integral());

  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();
  if(reweight){
    gPad->SetRightMargin(0.2);
    plot2D(canvas,mediator1_D_2D,"mediator_pt_eta_1_2D",outputDIR);
    plot2D(canvas,mediator2_D_2D,"mediator_pt_eta_2_2D",outputDIR);
    gPad->SetRightMargin(0.06);
    plot1D(canvas,mediator1_D_1D,"mediator_pt_eta_1_1D",outputDIR);
    plot1D(canvas,mediator2_D_1D,"mediator_pt_eta_2_1D",outputDIR);
  }

  TH2F* mediator1_N_2D = new TH2F("mediator1_N_2D","",125,0,1250,125,0,1250);
  TH2F* mediator2_N_2D = new TH2F("mediator2_N_2D","",125,0,1250,125,0,1250);
  mediator1_N_2D->Sumw2();
  mediator2_N_2D->Sumw2();

  for(int iBinX = 0; iBinX < mediator1_N_2D->GetNbinsX()+1; iBinX++){
    for(int iBinY = 0; iBinY < mediator1_N_2D->GetNbinsY()+1; iBinY++){
      mediator1_N_2D->SetBinContent(iBinX+1,iBinY+1,1);
    }
  }
  mediator1_N_2D->Scale(1./mediator1_N_2D->Integral());
  mediator1_N_2D->Divide(mediator1_D_2D);
  
  for(int iBinX = 0; iBinX < mediator2_N_2D->GetNbinsX()+1; iBinX++){
    for(int iBinY = 0; iBinY < mediator2_N_2D->GetNbinsY()+1; iBinY++){
      mediator2_N_2D->SetBinContent(iBinX+1,iBinY+1,1);
    }
  }
  mediator2_N_2D->Scale(1./mediator2_N_2D->Integral());
  mediator2_N_2D->Divide(mediator2_D_2D);


  TH1F* mediator1_N_1D = new TH1F("mediator1_N_1D","",80,0,1250);
  TH1F* mediator2_N_1D = new TH1F("mediator2_N_1D","",80,0,1250);
  mediator1_N_1D->Sumw2();
  mediator2_N_1D->Sumw2();
  
  for(int iBinX = 0; iBinX < mediator1_N_1D->GetNbinsX()+1; iBinX++)
    mediator1_N_1D->SetBinContent(iBinX+1,1);
  
  mediator1_N_1D->Scale(1./mediator1_N_1D->Integral());
  mediator1_N_1D->Divide(mediator1_D_1D);
  
  for(int iBinX = 0; iBinX < mediator1_N_1D->GetNbinsX()+1; iBinX++)
    mediator2_N_1D->SetBinContent(iBinX+1,1);

  mediator2_N_1D->Scale(1./mediator2_N_1D->Integral());
  mediator2_N_1D->Divide(mediator2_D_1D);

  // use them as event weight
  TH2F* mediator1_w_2D = new TH2F("mediator1_w_2D","",125,0,1250,125,0,1250);
  TH2F* mediator2_w_2D = new TH2F("mediator2_w_2D","",125,0,1250,125,0,1250);
  mediator1_w_2D->Sumw2();
  mediator2_w_2D->Sumw2();

  TH1F* mediator1_w_1D = new TH1F("mediator1_w_1D","",80,0,1250);
  TH1F* mediator2_w_1D = new TH1F("mediator2_w_1D","",80,0,1250);
  mediator1_w_1D->Sumw2();
  mediator2_w_1D->Sumw2();

  TH2F* mediator1_weighted_2D = new TH2F("mediator1_weighted_2D","",recoilBin_2D.size()-1,&recoilBin_2D[0],ptBin_2D.size()-1,&ptBin_2D[0]);
  TH2F* mediator2_weighted_2D = new TH2F("mediator2_weighted_2D","",recoilBin_2D.size()-1,&recoilBin_2D[0],ptBin_2D.size()-1,&ptBin_2D[0]);
  mediator1_weighted_2D->Sumw2();
  mediator2_weighted_2D->Sumw2();

  TH1F* mediator1_weighted_1D = new TH1F("mediator1_weighted_1D","",recoilBin_1D.size()-1,&recoilBin_1D[0]);
  TH1F* mediator2_weighted_1D = new TH1F("mediator2_weighted_1D","",recoilBin_1D.size()-1,&recoilBin_1D[0]);
  mediator1_weighted_1D->Sumw2();
  mediator2_weighted_1D->Sumw2();

  TH2F* mediator1_weighted_id1 = new TH2F("mediator1_weighted_id1","",recoilBin_2D.size()-1,&recoilBin_2D[0],ptBin_2D.size()-1,&ptBin_2D[0]);
  TH2F* mediator2_weighted_id1 = new TH2F("mediator2_weighted_id1","",recoilBin_2D.size()-1,&recoilBin_2D[0],ptBin_2D.size()-1,&ptBin_2D[0]);
  mediator1_weighted_id1->Sumw2();
  mediator2_weighted_id1->Sumw2();

  TH1F* mediator1_weighted_1D_id1 = new TH1F("mediator1_weighted_1D_id1","",recoilBin_1D.size()-1,&recoilBin_1D[0]);
  TH1F* mediator2_weighted_1D_id1 = new TH1F("mediator2_weighted_1D_id1","",recoilBin_1D.size()-1,&recoilBin_1D[0]);
  mediator1_weighted_1D_id1->Sumw2();
  mediator2_weighted_1D_id1->Sumw2();

  TH2F* mediator1_weighted_id2 = new TH2F("mediator1_weighted_id2","",recoilBin_2D.size()-1,&recoilBin_2D[0],ptBin_2D.size()-1,&ptBin_2D[0]);
  TH2F* mediator2_weighted_id2 = new TH2F("mediator2_weighted_id2","",recoilBin_2D.size()-1,&recoilBin_2D[0],ptBin_2D.size()-1,&ptBin_2D[0]);
  mediator1_weighted_id2->Sumw2();
  mediator2_weighted_id2->Sumw2();

  TH1F* mediator1_weighted_1D_id2 = new TH1F("mediator1_weighted_1D_id2","",recoilBin_1D.size()-1,&recoilBin_1D[0]);
  TH1F* mediator2_weighted_1D_id2 = new TH1F("mediator2_weighted_1D_id2","",recoilBin_1D.size()-1,&recoilBin_1D[0]);
  mediator1_weighted_1D_id2->Sumw2();
  mediator2_weighted_1D_id2->Sumw2();

  TTreeReader reader1 (tree1);
  TTreeReaderValue<double> genMediatorPt1   (reader1,"genMediatorPt");
  TTreeReaderValue<double> recoilPt1        (reader1,"pfMetPt");
  TTreeReaderValue<double> genMediatorEta1  (reader1,"genMediatorEta");
  TTreeReaderValue<double> genMediatorMass1  (reader1,"genMediatorMass");
  TTreeReaderValue<int>    id1  (reader1,"id");
  TTreeReaderValue<double> weight1  (reader1,"weight");
  TTreeReaderValue<double> genweight1  (reader1,"genWeight");
 
  vector<double> massPoints;

  cout<<"start event loop on first set"<<endl;
  long int ientry = 0;
  while(reader1.Next()){
    ientry++;
    std::cout.flush();
    if(ientry %100000 == 0) cout<<"\r"<<"ientry == "<<double(ientry)/tree1->GetEntries()*100<<" %";

    if(*genMediatorMass1 <= minMedMass) continue;

    if(std::find(massPoints.begin(),massPoints.end(),*genMediatorMass1) == massPoints.end() and useOnlyCommonPoints)
      massPoints.push_back(*genMediatorMass1);

    if(reweight){
      mediator1_w_1D->Fill(*recoilPt1,(*weight1)*(*genweight1)*mediator1_N_1D->GetBinContent(mediator1_D_1D->FindBin(*recoilPt1)));
      mediator1_w_2D->Fill(*recoilPt1,*genMediatorPt1,(*weight1)*(*genweight1)*mediator1_N_2D->GetBinContent(mediator1_D_2D->FindBin(*recoilPt1,*genMediatorPt1)));
    }

    if(reweight){
      mediator1_weighted_1D->Fill(*recoilPt1,mediator1_N_1D->GetBinContent(mediator1_D_1D->FindBin(*recoilPt1))*(*weight1)*(*genweight1));
      mediator1_weighted_2D->Fill(*recoilPt1,*genMediatorPt1,mediator1_N_2D->GetBinContent(mediator1_D_2D->FindBin(*recoilPt1,*genMediatorPt1))*(*weight1)*(*genweight1));
      if(*id1 == 1){
	mediator1_weighted_1D_id1->Fill(*recoilPt1,mediator1_N_1D->GetBinContent(mediator1_D_1D->FindBin(*recoilPt1))*(*weight1)*(*genweight1));
	mediator1_weighted_id1->Fill(*recoilPt1,*genMediatorPt1,mediator1_N_2D->GetBinContent(mediator1_D_2D->FindBin(*recoilPt1,*genMediatorPt1))*(*weight1)*(*genweight1));
      }
      if(*id1 == 2){
	mediator1_weighted_1D_id2->Fill(*recoilPt1,mediator1_N_1D->GetBinContent(mediator1_D_1D->FindBin(*recoilPt1))*(*weight1)*(*genweight1));    
	mediator1_weighted_id2->Fill(*recoilPt1,*genMediatorPt1,mediator1_N_2D->GetBinContent(mediator1_D_2D->FindBin(*recoilPt1,*genMediatorPt1))*(*weight1)*(*genweight1));
      }
    }
    else{
      mediator1_weighted_1D->Fill(*recoilPt1,(*weight1)*(*genweight1));
      mediator1_weighted_2D->Fill(*recoilPt1,*genMediatorPt1,(*weight1)*(*genweight1));
      if(*id1 == 1){
	mediator1_weighted_1D_id1->Fill(*recoilPt1,(*weight1)*(*genweight1));
	mediator1_weighted_id1->Fill(*recoilPt1,*genMediatorPt1,(*weight1)*(*genweight1));
      }
      if(*id1 == 2){
	mediator1_weighted_1D_id2->Fill(*recoilPt1,(*weight1)*(*genweight1));    
	mediator1_weighted_id2->Fill(*recoilPt1,*genMediatorPt1,(*weight1)*(*genweight1));    
      }
      
    }    
  }

  TTreeReader reader2 (tree2);
  TTreeReaderValue<double> genMediatorPt2   (reader2,"genMediatorPt");
  TTreeReaderValue<double> recoilPt2   (reader2,"pfMetPt");
  TTreeReaderValue<double> genMediatorEta2  (reader2,"genMediatorEta");
  TTreeReaderValue<double> genMediatorMass2  (reader2,"genMediatorMass");
  TTreeReaderValue<int>    id2  (reader2,"id");
  TTreeReaderValue<double> weight2  (reader2,"weight");
  TTreeReaderValue<double> genweight2  (reader2,"genWeight");
    
  cout<<endl;
  cout<<"start event loop on the second set"<<endl;
  ientry = 0;
  while(reader2.Next()){
    ientry++;
    std::cout.flush();
    if(ientry %100000 == 0) cout<<"\r"<<"ientry == "<<double(ientry)/tree2->GetEntries()*100<<" %";;
    
    if(*genMediatorMass2 <= minMedMass) continue;

    if(std::find(massPoints.begin(),massPoints.end(),*genMediatorMass2) == massPoints.end() and useOnlyCommonPoints){
      continue;
    }

    if(reweight){
      mediator2_w_1D->Fill(*recoilPt2,mediator2_N_1D->GetBinContent(mediator2_D_1D->FindBin(*recoilPt2))*(*weight2)*(*genweight2));
      mediator2_w_2D->Fill(*recoilPt2,*genMediatorPt2,mediator2_N_2D->GetBinContent(mediator2_D_2D->FindBin(*recoilPt2,*genMediatorPt2))*(*weight2)*(*genweight2));
    }
    if(reweight){
      mediator2_weighted_1D->Fill(*recoilPt2,mediator2_N_1D->GetBinContent(mediator2_D_1D->FindBin(*recoilPt2))*(*weight2)*(*genweight2));
      mediator2_weighted_2D->Fill(*recoilPt2,*genMediatorPt2,mediator2_N_2D->GetBinContent(mediator2_D_2D->FindBin(*recoilPt2,*genMediatorPt2))*(*weight2)*(*genweight2));
    
      if(*id2 == 1){
	mediator2_weighted_1D_id1->Fill(*recoilPt2,mediator2_N_1D->GetBinContent(mediator2_D_1D->FindBin(*recoilPt2))*(*weight2)*(*genweight2));
	mediator2_weighted_id1->Fill(*recoilPt2,*genMediatorPt2,mediator2_N_2D->GetBinContent(mediator2_D_2D->FindBin(*recoilPt2,*genMediatorPt2))*(*weight2)*(*genweight2));
      }
      if(*id2 == 2){
	mediator2_weighted_1D_id2->Fill(*recoilPt2,mediator2_N_1D->GetBinContent(mediator2_D_1D->FindBin(*recoilPt2))*(*weight2)*(*genweight2));
	mediator2_weighted_id2->Fill(*recoilPt2,*genMediatorPt2,mediator2_N_2D->GetBinContent(mediator2_D_2D->FindBin(*recoilPt2,*genMediatorPt2))*(*weight2)*(*genweight2));
      }
    }
    else{
      mediator2_weighted_1D->Fill(*recoilPt2,(*weight2)*(*genweight2));
      mediator2_weighted_2D->Fill(*recoilPt2,*genMediatorPt2,(*weight2)*(*genweight2));
      if(*id2 == 1){
	mediator2_weighted_1D_id1->Fill(*recoilPt2,(*weight2)*(*genweight2));
	mediator2_weighted_id1->Fill(*recoilPt2,*genMediatorPt2,(*weight2)*(*genweight2));
      }
      if(*id2 == 2){
	mediator2_weighted_1D_id2->Fill(*recoilPt2,(*weight2)*(*genweight2));
	mediator2_weighted_id2->Fill(*recoilPt2,*genMediatorPt2,(*weight2)*(*genweight2));
      }
    }    
  }
  cout<<endl;

  mediator1_w_2D->Scale(1./mediator1_w_2D->Integral());
  mediator2_w_2D->Scale(1./mediator2_w_2D->Integral());
  mediator1_w_1D->Scale(1./mediator1_w_1D->Integral());
  mediator2_w_1D->Scale(1./mediator2_w_1D->Integral());
  
  if(reweight){
    gPad->SetRightMargin(0.2);
    plot2D(canvas,mediator1_w_2D,"recoil_medpt_weighted1_2D",outputDIR);
    plot2D(canvas,mediator2_w_2D,"recoil_medpt_weighted2_2D",outputDIR);
    gPad->SetRightMargin(0.06);
    plot1D(canvas,mediator1_w_1D,"recoil_weighted1_1D",outputDIR);
    plot1D(canvas,mediator2_w_1D,"recoil_weighted2_1D",outputDIR);
  }
  
  cout<<"make ids "<<endl;

  // make efficiency
  TH2* eff_med1_id1 = (TH2*) mediator1_weighted_id1->Clone("eff_med1_id1");
  eff_med1_id1->Reset();
  eff_med1_id1->Divide(mediator1_weighted_id1,mediator1_weighted_2D,1,1,"B");
  TH2* eff_med1_id2 = (TH2*) mediator1_weighted_id2->Clone("eff_med1_id2");
  eff_med1_id2->Reset();
  eff_med1_id2->Divide(mediator1_weighted_id2,mediator1_weighted_2D,1,1,"B");
  TH2* eff_med2_id1 = (TH2*) mediator2_weighted_id1->Clone("eff_med2_id1");
  eff_med2_id1->Reset();
  eff_med2_id1->Divide(mediator2_weighted_id1,mediator2_weighted_2D,1,1,"B");
  TH2* eff_med2_id2 = (TH2*) mediator2_weighted_id2->Clone("eff_med2_id2");
  eff_med2_id2->Reset();
  eff_med2_id2->Divide(mediator2_weighted_id2,mediator2_weighted_2D,1,1,"B");
  //make scale factors
  eff_med1_id1->Divide(eff_med2_id1);
  eff_med1_id2->Divide(eff_med2_id2);

  gPad->SetRightMargin(0.2);
  plot2D(canvas,eff_med1_id1,"recoil_medpt_SF1",outputDIR,true);
  plot2D(canvas,eff_med1_id2,"recoil_medpt_SF2",outputDIR,true);

  // make efficiency
  TH1* eff_med1_1D_id1 = (TH1*) mediator1_weighted_1D_id1->Clone("eff_med1_1D_id1");
  eff_med1_1D_id1->Reset();
  eff_med1_1D_id1->Divide(mediator1_weighted_1D_id1,mediator1_weighted_1D,1,1,"B");

  TH1* eff_med1_1D_id2 = (TH1*) mediator1_weighted_1D_id2->Clone("eff_med1_1D_id2");
  eff_med1_1D_id2->Reset();
  eff_med1_1D_id2->Divide(mediator1_weighted_1D_id2,mediator1_weighted_1D,1,1,"B");
  TH1* eff_med2_1D_id1 = (TH1*) mediator2_weighted_1D_id1->Clone("eff_med2_1D_id1");
  eff_med2_1D_id1->Reset();
  eff_med2_1D_id1->Divide(mediator2_weighted_1D_id1,mediator2_weighted_1D,1,1,"B");
  TH1* eff_med2_1D_id2 = (TH1*) mediator2_weighted_1D_id2->Clone("eff_med2_1D_id2");
  eff_med2_1D_id2->Reset();
  eff_med2_1D_id2->Divide(mediator2_weighted_1D_id2,mediator2_weighted_1D,1,1,"B");

  eff_med1_1D_id1->Divide(eff_med2_1D_id1);
  eff_med1_1D_id2->Divide(eff_med2_1D_id2);
  gPad->SetRightMargin(0.06);
  plot1D(canvas,eff_med1_1D_id1,"recoil_SF1",outputDIR,true);
  plot1D(canvas,eff_med1_1D_id2,"recoil_SF2",outputDIR,true);
  
}
