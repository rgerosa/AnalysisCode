#include "../CMS_lumi.h"

void plot2D(TCanvas* canvas, TH2* histo, string name, string outputDIR, bool text = false){

  histo->GetYaxis()->SetTitle("mediator p_{T} [GeV]");
  histo->GetXaxis()->SetTitle("mediator |#eta|");
  histo->GetZaxis()->SetTitle("a.u.");

  if(text){
    for(int iBinX = 0; iBinX < histo->GetNbinsX()+1; iBinX++){
      for(int iBinY = 0; iBinY < histo->GetNbinsY()+1; iBinY++){
	if(histo->GetBinContent(iBinX+1,iBinY+1) > 5 or histo->GetBinContent(iBinX+1,iBinY+1) < 0.01)
	  histo->SetBinContent(iBinX+1,iBinY+1,0);
      }
    }
  }


  if(not text)
    histo->Draw("colz");
  else
    histo->Draw("colz texte");
  CMS_lumi(canvas,"2.6",true);
  canvas->SaveAs((outputDIR+"/"+name+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/"+name+".pdf").c_str(),"pdf");

}

void makeEfficiencyComparison(string namefile1, string namefile2, string outputDIR, bool reweight = true){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  gStyle->SetPaintTextFormat(".2f");
  system(("mkdir -p "+outputDIR).c_str());

  TFile* file1 = TFile::Open(namefile1.c_str());
  TFile* file2 = TFile::Open(namefile2.c_str());

  TTree* tree1 = (TTree*) file1->Get("tree");
  TTree* tree2 = (TTree*) file2->Get("tree");

  // weight properly according to mediator pt:eta
  vector<float> etaBin = {0,1,2,4};
  vector<float> ptBin  = {100,130,160,200,240,280,320,400,480,560,700,1000};

  TH2F* mediator1_D = new TH2F("mediator1_D","",10,0,5,100,100,1000);
  TH2F* mediator2_D = new TH2F("mediator2_D","",10,0,5,100,100,1000);
  mediator1_D->Sumw2();
  mediator2_D->Sumw2();

  tree1->Draw("genMediatorPt:abs(genMediatorEta) >> mediator1_D","(genWeight*weight)*(genMediatorPt > 0)","goff");
  tree2->Draw("genMediatorPt:abs(genMediatorEta) >> mediator2_D","(genWeight*weight)*(genMediatorPt > 0)","goff");
  mediator1_D->Scale(1./mediator1_D->Integral());
  mediator2_D->Scale(2./mediator2_D->Integral());

  TCanvas* canvas = new TCanvas("canvas","",600,650);
  canvas->cd();
  gPad->SetRightMargin(0.2);
  plot2D(canvas,mediator1_D,"mediator_pt_eta_1",outputDIR);
  plot2D(canvas,mediator2_D,"mediator_pt_eta_2",outputDIR);

  TH2F* mediator1_N = new TH2F("mediator1_N","",10,0,5,100,100,1000);
  TH2F* mediator2_N = new TH2F("mediator2_N","",10,0,5,100,100,1000);
  mediator1_N->Sumw2();
  mediator2_N->Sumw2();

  for(int iBinX = 0; iBinX < mediator1_N->GetNbinsX()+1; iBinX++){
    for(int iBinY = 0; iBinY < mediator1_N->GetNbinsY()+1; iBinY++){
      mediator1_N->SetBinContent(iBinX+1,iBinY+1,1);
    }
  }
  mediator1_N->Scale(1./mediator1_N->Integral());
  mediator1_D->Divide(mediator1_N);
  
  for(int iBinX = 0; iBinX < mediator1_N->GetNbinsX()+1; iBinX++){
    for(int iBinY = 0; iBinY < mediator1_N->GetNbinsY()+1; iBinY++){
      mediator2_N->SetBinContent(iBinX+1,iBinY+1,1);
    }
  }
  mediator2_N->Scale(1./mediator2_N->Integral());
  mediator2_D->Divide(mediator2_N);

  // use them as event weight
  TH2F* mediator1_w = new TH2F("mediator1_w","",10,0,5,100,100,1000);
  TH2F* mediator2_w = new TH2F("mediator2_w","",10,0,5,100,100,1000);
  mediator1_w->Sumw2();
  mediator2_w->Sumw2();

  TH2F* mediator1_weighted = new TH2F("mediator1_weighted","",etaBin.size()-1,&etaBin[0],ptBin.size()-1,&ptBin[0]);
  TH2F* mediator2_weighted = new TH2F("mediator2_weighted","",etaBin.size()-1,&etaBin[0],ptBin.size()-1,&ptBin[0]);
  mediator1_weighted->Sumw2();
  mediator2_weighted->Sumw2();

  TH2F* mediator1_weighted_id1 = new TH2F("mediator1_weighted_id1","",etaBin.size()-1,&etaBin[0],ptBin.size()-1,&ptBin[0]);
  TH2F* mediator2_weighted_id1 = new TH2F("mediator2_weighted_id1","",etaBin.size()-1,&etaBin[0],ptBin.size()-1,&ptBin[0]);
  mediator1_weighted_id1->Sumw2();
  mediator2_weighted_id1->Sumw2();

  TH2F* mediator1_weighted_id2 = new TH2F("mediator1_weighted_id2","",etaBin.size()-1,&etaBin[0],ptBin.size()-1,&ptBin[0]);
  TH2F* mediator2_weighted_id2 = new TH2F("mediator2_weighted_id2","",etaBin.size()-1,&etaBin[0],ptBin.size()-1,&ptBin[0]);
  mediator1_weighted_id2->Sumw2();
  mediator2_weighted_id2->Sumw2();

  TTreeReader reader1 (tree1);
  TTreeReaderValue<double> genMediatorPt1   (reader1,"genMediatorPt");
  TTreeReaderValue<double> genMediatorEta1  (reader1,"genMediatorEta");
  TTreeReaderValue<int> id1  (reader1,"id");
  TTreeReaderValue<double> weight1  (reader1,"weight");
  TTreeReaderValue<double> genweight1  (reader1,"genWeight");
 
  TTreeReader reader2 (tree2);
  TTreeReaderValue<double> genMediatorPt2   (reader2,"genMediatorPt");
  TTreeReaderValue<double> genMediatorEta2  (reader2,"genMediatorEta");
  TTreeReaderValue<int> id2  (reader2,"id");
  TTreeReaderValue<double> weight2  (reader2,"weight");
  TTreeReaderValue<double> genweight2  (reader2,"genWeight");
  

  while(reader1.Next()){
    
    mediator1_w->Fill(fabs(*genMediatorEta1),*genMediatorPt1,1./mediator1_D->GetBinContent(mediator1_D->FindBin(fabs(*genMediatorEta1),*genMediatorPt1))
		    *(*weight1)*(*genweight1));

    if(reweight){
      mediator1_weighted->Fill(fabs(*genMediatorEta1),*genMediatorPt1,1./mediator1_D->GetBinContent(mediator1_D->FindBin(fabs(*genMediatorEta1),*genMediatorPt1))
			       *(*weight1)*(*genweight1));
      if(*id1 == 1)
	mediator1_weighted_id1->Fill(fabs(*genMediatorEta1),*genMediatorPt1,1./mediator1_D->GetBinContent(mediator1_D->FindBin(fabs(*genMediatorEta1),*genMediatorPt1))
				     *(*weight1)*(*genweight1));
      if(*id1 == 2)
	mediator1_weighted_id2->Fill(fabs(*genMediatorEta1),*genMediatorPt1,1./mediator1_D->GetBinContent(mediator1_D->FindBin(fabs(*genMediatorEta1),*genMediatorPt1))
				   *(*weight1)*(*genweight1));    
    }
    else{
      mediator1_weighted->Fill(fabs(*genMediatorEta1),*genMediatorPt1,(*weight1)*(*genweight1));
      if(*id1 == 1)
	mediator1_weighted_id1->Fill(fabs(*genMediatorEta1),*genMediatorPt1,(*weight1)*(*genweight1));
      if(*id1 == 2)
	mediator1_weighted_id2->Fill(fabs(*genMediatorEta1),*genMediatorPt1,(*weight1)*(*genweight1));    
      
    }    
  }

  while(reader2.Next()){
    mediator2_w->Fill(fabs(*genMediatorEta2),*genMediatorPt2,1./mediator2_D->GetBinContent(mediator2_D->FindBin(fabs(*genMediatorEta2),*genMediatorPt2))
		    *(*weight2)*(*genweight2));
    if(reweight){
      mediator2_weighted->Fill(fabs(*genMediatorEta2),*genMediatorPt2,1./mediator2_D->GetBinContent(mediator2_D->FindBin(fabs(*genMediatorEta2),*genMediatorPt2))
			       *(*weight2)*(*genweight2));
      if(*id2 == 1)
	mediator2_weighted_id1->Fill(fabs(*genMediatorEta2),*genMediatorPt2,1./mediator2_D->GetBinContent(mediator2_D->FindBin(fabs(*genMediatorEta2),*genMediatorPt2))
				     *(*weight2)*(*genweight2));
      if(*id2 == 2)
	mediator2_weighted_id2->Fill(fabs(*genMediatorEta2),*genMediatorPt2,1./mediator2_D->GetBinContent(mediator2_D->FindBin(fabs(*genMediatorEta2),*genMediatorPt2))
				     *(*weight2)*(*genweight2));
    }
    else{
      mediator2_weighted->Fill(fabs(*genMediatorEta2),*genMediatorPt2,(*weight2)*(*genweight2));
      if(*id2 == 1)
	mediator2_weighted_id1->Fill(fabs(*genMediatorEta2),*genMediatorPt2,(*weight2)*(*genweight2));
      if(*id2 == 2)
	mediator2_weighted_id2->Fill(fabs(*genMediatorEta2),*genMediatorPt2,(*weight2)*(*genweight2));
    }
  }

  mediator1_weighted_id1->Divide(mediator1_weighted);
  mediator1_weighted_id2->Divide(mediator1_weighted);
  mediator2_weighted_id1->Divide(mediator2_weighted);
  mediator2_weighted_id2->Divide(mediator2_weighted);

  mediator1_w->Scale(1./mediator1_w->Integral());
  mediator2_w->Scale(1./mediator2_w->Integral());
  plot2D(canvas,mediator1_w,"mediator_pt_eta_weighted1",outputDIR);
  plot2D(canvas,mediator2_w,"mediator_pt_eta_weighted2",outputDIR);

  mediator1_weighted_id1->Divide(mediator2_weighted_id1);
  mediator1_weighted_id2->Divide(mediator2_weighted_id2);

  plot2D(canvas,mediator1_weighted_id1,"mediator_pt_eta_SF1",outputDIR,true);
  plot2D(canvas,mediator1_weighted_id2,"mediator_pt_eta_SF2",outputDIR,true);

}
