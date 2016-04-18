#include "CMS_lumi.h"
#include "histoUtils.h"

void drawCorrelationSignalAndBackground(string inputFileName, 
					int category,
					vector<string> observables2D,
					vector<string> observablesLatex2D_X, 
					vector<string> observablesLatex2D_Y,
					string interaction  = "Vector",
					string mediatorMass = "1000", 
					string DMMass = "50",
					bool isHiggsInvisible = false,
					bool isLog = false,
					bool divideByBinWidth = false){

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  gStyle->SetOptStat(0);

  TFile* templateFile = TFile::Open(inputFileName.c_str());

  vector<TH2*> bkghist;
  vector<TH2*> monoJhist;
  vector<TH2*> monoWhist;
  vector<TH2*> monoZhist;

  for(auto obs : observables2D){    
    bkghist.push_back(roll2DHistograms((TH1*) templateFile->FindObjectAny(("zinvhist_"+obs).c_str()),obs,category,true));
    if(not isHiggsInvisible){
      monoJhist.push_back(roll2DHistograms((TH1*) templateFile->FindObjectAny(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+obs).c_str()),obs,category,true));
      monoWhist.push_back(roll2DHistograms((TH1*) templateFile->FindObjectAny(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+obs).c_str()),obs,category,true));
      monoZhist.push_back(roll2DHistograms((TH1*) templateFile->FindObjectAny(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+obs).c_str()),obs,category,true));
    }
    else{
      monoJhist.push_back(roll2DHistograms((TH1*) templateFile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+obs).c_str()),obs,category,true));
      monoWhist.push_back(roll2DHistograms((TH1*) templateFile->FindObjectAny(("wHhist_"+mediatorMass+"_"+obs).c_str()),obs,category,true));
      monoZhist.push_back(roll2DHistograms((TH1*) templateFile->FindObjectAny(("zHhist_"+mediatorMass+"_"+obs).c_str()),obs,category,true));
    }
  }

  TCanvas* canvas = new TCanvas("canvas","canvas",700,600);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetLeftMargin(0.11);
  canvas->SetRightMargin(0.15);
  canvas->SetTopMargin(0.06);
  canvas->SetBottomMargin(0.10);

  if(isLog)
    canvas->SetLogz();

  if(divideByBinWidth){
    for(auto hist : bkghist)
      hist->Scale(1.0,"width");
    for(auto hist : monoJhist)
      hist->Scale(1.0,"width");
    for(auto hist : monoWhist)
      hist->Scale(1.0,"width");
    for(auto hist : monoZhist)
      hist->Scale(1.0,"width");
  }
  
  int iObs = 0;
  for(auto obs : observables2D){    

    bin2D bins = selectBinning2D (obs,0);
    TH2* monojet = cloneHistoIncludingOverUnderFlow(monoJhist.at(iObs));
    monojet->Smooth();
    monoJhist.at(iObs)->Smooth();

    TGraph2D* graph2D_1 = new TGraph2D (monojet);
    graph2D_1->SetNpx(100);
    graph2D_1->SetNpy(100);
    for(int iBinX = 0; iBinX <= monojet->GetNbinsX(); iBinX++){
      for(int iBinY = 0; iBinY <= monojet->GetNbinsY(); iBinY++){
	graph2D_1->Interpolate(monojet->GetXaxis()->GetBinCenter(iBinX),monojet->GetYaxis()->GetBinCenter(iBinY));
      }
    }	

    TH2F* frame = new TH2F("","",int(bins.binX.size()-1),&bins.binX[0],int(bins.binY.size()-1),&bins.binY[0]);

    frame->GetXaxis()->SetTitle(observablesLatex2D_X.at(iObs).c_str());
    frame->GetYaxis()->SetTitle(observablesLatex2D_Y.at(iObs).c_str());
    frame->GetZaxis()->SetTitle("Events/GeV");
    frame->GetZaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetTitleSize(0.045);
          	      

    TH2* hist = graph2D_1->GetHistogram(); // make sure to fill all the bins
    
    TProfile* profileX = monoJhist.at(iObs)->ProfileX("_pfx");
    profileX->SetMarkerStyle(20);
    profileX->SetMarkerSize(0.8);
    profileX->SetMarkerColor(kBlack);
    profileX->SetLineColor(kBlack);
        
    frame->Draw();
    CMS_lumi(canvas,"2.30",true);    
    hist->Draw("colz same");
    profileX->Draw("EP same");
    
    canvas->RedrawAxis("sameaxis");
    canvas->SaveAs(("monoJ_"+obs+".pdf").c_str());
    canvas->SaveAs(("monoJ_"+obs+".png").c_str());

    frame->Draw();
    CMS_lumi(canvas,"2.30",true);
    monoWhist.at(iObs)->Smooth();

    TH2* monoW = cloneHistoIncludingOverUnderFlow(monoWhist.at(iObs));
    monoW->Smooth();
    monoWhist.at(iObs)->Smooth();

    TGraph2D* graph2D_2 = new TGraph2D (monoW);
    for(int iBinX = 0; iBinX < monoW->GetNbinsX()+1; iBinX++){
      for(int iBinY = 0; iBinY < monoW->GetNbinsY()+1; iBinY++){
	graph2D_2->Interpolate(monoW->GetXaxis()->GetBinCenter(iBinX+1),monoW->GetXaxis()->GetBinCenter(iBinY+1));
      }
    }

    graph2D_2->GetHistogram()->Draw("colz same");
    
    profileX = monoWhist.at(iObs)->ProfileX("_pfx");
    profileX->SetMarkerStyle(20);
    profileX->SetMarkerSize(0.8);
    profileX->SetMarkerColor(kBlack);
    profileX->SetLineColor(kBlack);
    
    profileX->Draw("EP same");
    canvas->RedrawAxis("sameaxis");
    
    canvas->SaveAs(("monoW_"+obs+".pdf").c_str());
    canvas->SaveAs(("monoW_"+obs+".png").c_str());

    frame->Draw();
    CMS_lumi(canvas,"2.30",true);

    TH2* monoZ = cloneHistoIncludingOverUnderFlow(monoZhist.at(iObs));
    monoZ->Smooth();
    monoZhist.at(iObs)->Smooth();

    TGraph2D* graph2D_3 = new TGraph2D (monoZ);

    for(int iBinX = 0; iBinX < monoZ->GetNbinsX()+1; iBinX++){
      for(int iBinY = 0; iBinY < monoZ->GetNbinsY()+1; iBinY++){
	graph2D_3->Interpolate(monoZ->GetXaxis()->GetBinCenter(iBinX+1),monoZ->GetXaxis()->GetBinCenter(iBinY+1));
      }
    }

    graph2D_3->GetHistogram()->Draw("colz same");

    profileX = monoZhist.at(iObs)->ProfileX("_pfx");
    profileX->SetMarkerStyle(20);
    profileX->SetMarkerSize(0.8);
    profileX->SetMarkerColor(kBlack);
    profileX->SetLineColor(kBlack);
    
    canvas->RedrawAxis("sameaxis");
    canvas->SaveAs(("monoZ_"+obs+".pdf").c_str());
    canvas->SaveAs(("monoZ_"+obs+".png").c_str());

    // plot MC   
    frame->Draw();
    CMS_lumi(canvas,"2.30",true);
    
    TH2* bkg = cloneHistoIncludingOverUnderFlow(bkghist.at(iObs));
    bkg->Smooth();
    graph2D_1  = new TGraph2D(bkg);
    for(int iBinX = 0; iBinX < bkg->GetNbinsX()+1; iBinX++){
      for(int iBinY = 0; iBinY < bkg->GetNbinsY()+1; iBinY++){
	graph2D_1->Interpolate(bkg->GetXaxis()->GetBinCenter(iBinX+1),bkg->GetXaxis()->GetBinCenter(iBinY+1));
      }
    }
    graph2D_1->GetHistogram()->Draw("colz same");
    profileX = bkghist.at(iObs)->ProfileX("_pfx",1,bkghist.at(iObs)->GetNbinsX());
    profileX->SetMarkerStyle(20);
    profileX->SetMarkerSize(0.8);
    profileX->SetMarkerColor(kBlack);
    profileX->SetLineColor(kBlack);
    profileX->Draw("EP same");
    canvas->RedrawAxis("sameaxis");

    canvas->SaveAs(("znn_"+obs+"_MC.pdf").c_str());
    canvas->SaveAs(("znn_"+obs+"_MC.png").c_str());
    iObs++;
  }
  

}

void drawCorrelationPlotsFrom2D(string inputFileName, string controlRegion, vector<string> observables2D, vector<string> observablesLatex2D_X, vector<string> observablesLatex2D_Y,string mediatorMass = "1000", string DMMass = "50",bool isLog = false){


  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  gStyle->SetOptStat(0);

  TFile* templateFile = TFile::Open(inputFileName.c_str());

  vector<TH2*> datahist;
  vector<TH2*> bkghist;
  vector<TH2*> monoJhist;
  vector<TH2*> monoWhist;
  vector<TH2*> monoZhist;
  
  for(auto obs : observables2D){
    if(controlRegion != "SR"){
      datahist.push_back((TH2*) templateFile->FindObjectAny(("datahist"+controlRegion+"_"+obs+"_2D").c_str()));

      if(controlRegion == "topmu" or controlRegion == "topel"){
	bkghist.push_back((TH2*) templateFile->FindObjectAny(("tbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
	bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("dbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
	bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("qbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
	bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("vlbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
	bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("vllbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
      }
      else if(controlRegion == "wmn" or controlRegion == "wen"){
	bkghist.push_back((TH2*) templateFile->FindObjectAny(("vlbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
        bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("dbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
        bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("qbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
        bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("tbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
        bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("vllbkghist"+controlRegion+"_"+obs+"_2D").c_str()));	
      }
      else if(controlRegion == "zmm" or controlRegion == "zee"){
	bkghist.push_back((TH2*) templateFile->FindObjectAny(("vllbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
        bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("dbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
        bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("qbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
        bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("tbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
        bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("vlbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
      }
      else if(controlRegion == "gam"){
	bkghist.push_back((TH2*) templateFile->FindObjectAny(("gbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
        bkghist.back()->Add((TH2*) templateFile->FindObjectAny(("qbkghist"+controlRegion+"_"+obs+"_2D").c_str()));
      }
    }
    if(controlRegion == "SR"){
      monoJhist.push_back((TH2*) templateFile->FindObjectAny(("monoJhist_"+mediatorMass+"_"+DMMass+"_"+obs+"_2D").c_str()));
      monoWhist.push_back((TH2*) templateFile->FindObjectAny(("monoWhist_"+mediatorMass+"_"+DMMass+"_"+obs+"_2D").c_str()));
      monoZhist.push_back((TH2*) templateFile->FindObjectAny(("monoZhist_"+mediatorMass+"_"+DMMass+"_"+obs+"_2D").c_str()));
    }
  }

  TCanvas* canvas = new TCanvas("canvas","canvas",700,600);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetLeftMargin(0.11);
  canvas->SetRightMargin(0.15);
  canvas->SetTopMargin(0.06);
  canvas->SetBottomMargin(0.10);

  if(isLog)
    canvas->SetLogz();

  for(auto hist : datahist)
    hist->Scale(1.0,"width");
  for(auto hist : bkghist)
    hist->Scale(1.0,"width");
  for(auto hist : monoJhist)
    hist->Scale(1.0,"width");
  for(auto hist : monoWhist)
    hist->Scale(1.0,"width");
  for(auto hist : monoZhist)
    hist->Scale(1.0,"width");

  int iObs = 0;
  for(auto obs : observables2D){    

    bin2D bins = selectBinning2D (obs,0);

    if(controlRegion == "SR"){ // just plot correlation matrxi + profile
      monoJhist.at(iObs)->Smooth();
      TGraph2D* graph2D_1 = new TGraph2D (monoJhist.at(iObs));
      graph2D_1->SetNpx(100);
      graph2D_1->SetNpy(100);
      for(int iBinX = 0; iBinX <= monoJhist.at(iObs)->GetNbinsX(); iBinX++){
	for(int iBinY = 0; iBinY <= monoJhist.at(iObs)->GetNbinsY(); iBinY++)
	  graph2D_1->Interpolate(monoJhist.at(iObs)->GetXaxis()->GetBinCenter(iBinX),monoJhist.at(iObs)->GetYaxis()->GetBinCenter(iBinY));
      }	

      TH2F* frame = new TH2F("","",int(bins.binX.size()-1),monoJhist.at(iObs)->GetXaxis()->GetBinCenter(1),monoJhist.at(iObs)->GetXaxis()->GetBinCenter(monoJhist.at(iObs)->GetNbinsX()),int(bins.binY.size()-1),monoJhist.at(iObs)->GetYaxis()->GetBinCenter(1),monoJhist.at(iObs)->GetYaxis()->GetBinCenter(monoJhist.at(iObs)->GetNbinsY()));

      frame->GetXaxis()->SetTitle(observablesLatex2D_X.at(iObs).c_str());
      frame->GetYaxis()->SetTitle(observablesLatex2D_Y.at(iObs).c_str());
      frame->GetZaxis()->SetTitle("Events/GeV");
      frame->GetZaxis()->SetTitleSize(0.045);
      frame->GetYaxis()->SetTitleSize(0.045);
      frame->GetXaxis()->SetTitleSize(0.045);
          	      

      TH2* hist = graph2D_1->GetHistogram(); // make sure to fill all the bins

      TProfile* profileX = monoJhist.at(iObs)->ProfileX("_pfx");
      profileX->SetMarkerStyle(20);
      profileX->SetMarkerSize(0.8);
      profileX->SetMarkerColor(kBlack);
      profileX->SetLineColor(kBlack);


      frame->Draw();
      CMS_lumi(canvas,"2.30",true);

      hist->Draw("colz same");
      profileX->Draw("EP same");

      TString corr = Form("corr. = %0.2f",monoJhist.at(iObs)->GetCorrelationFactor());
      TLatex correlation;
      correlation.SetNDC(); correlation.SetTextSize(0.035);
      correlation.DrawText(0.65,0.88,corr.Data());

      canvas->RedrawAxis("sameaxis");
      canvas->SaveAs(("monoJ_"+obs+".pdf").c_str());
      canvas->SaveAs(("monoJ_"+obs+".png").c_str());

      frame->Draw();
      CMS_lumi(canvas,"2.30",true);
      monoWhist.at(iObs)->Smooth();

      TGraph2D* graph2D_2 = new TGraph2D (monoWhist.at(iObs));
      for(int iBinX = 0; iBinX < monoWhist.at(iObs)->GetNbinsX()+1; iBinX++){
	for(int iBinY = 0; iBinY < monoWhist.at(iObs)->GetNbinsY()+1; iBinY++){
	  graph2D_2->Interpolate(monoWhist.at(iObs)->GetXaxis()->GetBinCenter(iBinX+1),monoWhist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1));
	}
      }

      graph2D_2->GetHistogram()->Draw("colz same");

      profileX = monoWhist.at(iObs)->ProfileX("_pfx");
      profileX->SetMarkerStyle(20);
      profileX->SetMarkerSize(0.8);
      profileX->SetMarkerColor(kBlack);
      profileX->SetLineColor(kBlack);

      profileX->Draw("EP same");
      corr = Form("corr = %0.2f",monoWhist.at(iObs)->GetCorrelationFactor());
      correlation.DrawText(0.7,0.7,corr.Data());
      canvas->RedrawAxis("sameaxis");

      canvas->SaveAs(("monoW_"+obs+".pdf").c_str());
      canvas->SaveAs(("monoW_"+obs+".png").c_str());

      frame->Draw();
      CMS_lumi(canvas,"2.30",true);
      profileX = monoZhist.at(iObs)->ProfileX("_pfx");
      profileX->SetMarkerStyle(20);
      profileX->SetMarkerSize(0.8);
      profileX->SetMarkerColor(kBlack);
      profileX->SetLineColor(kBlack);
      TGraph2D* graph2D_3 = new TGraph2D (monoZhist.at(iObs));
      for(int iBinX = 0; iBinX < monoZhist.at(iObs)->GetNbinsX()+1; iBinX++){
	for(int iBinY = 0; iBinY < monoZhist.at(iObs)->GetNbinsY()+1; iBinY++){
	  graph2D_3->Interpolate(monoZhist.at(iObs)->GetXaxis()->GetBinCenter(iBinX+1),monoZhist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1));
	}
      }

      graph2D_3->GetHistogram()->Draw("colz same");
      profileX->Draw("EP same");
      corr = Form("corr. = %0.2f",monoZhist.at(iObs)->GetCorrelationFactor());
      correlation.DrawText(0.7,0.7,corr.Data());
      canvas->RedrawAxis("sameaxis");

      canvas->SaveAs(("monoZ_"+obs+".pdf").c_str());
      canvas->SaveAs(("monoZ_"+obs+".png").c_str());
    }
    else{

      // plot MC
      TH2F* frame = new TH2F("","",int(bins.binX.size()-1),bkghist.at(iObs)->GetXaxis()->GetBinCenter(1),bkghist.at(iObs)->GetXaxis()->GetBinCenter(bkghist.at(iObs)->GetNbinsX()),int(bins.binY.size()-1),bkghist.at(iObs)->GetYaxis()->GetBinCenter(1),bkghist.at(iObs)->GetYaxis()->GetBinCenter(bkghist.at(iObs)->GetNbinsY()));

      frame->Draw();
      CMS_lumi(canvas,"2.30",true);

      TGraph2D* graph2D_1  = new TGraph2D(bkghist.at(iObs));
      for(int iBinX = 0; iBinX < bkghist.at(iObs)->GetNbinsX()+1; iBinX++){
        for(int iBinY = 0; iBinY < bkghist.at(iObs)->GetNbinsY()+1; iBinY++){
          graph2D_1->Interpolate(bkghist.at(iObs)->GetXaxis()->GetBinCenter(iBinX+1),bkghist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1));
        }
      }
      graph2D_1->GetHistogram()->Draw("colz same");
      TProfile* profileX = bkghist.at(iObs)->ProfileX("_pfx",1,bkghist.at(iObs)->GetNbinsX());
      profileX->SetMarkerStyle(20);
      profileX->SetMarkerSize(0.8);
      profileX->SetMarkerColor(kBlack);
      profileX->SetLineColor(kBlack);
      profileX->Draw("EP same");
      TString corr = Form("corr. = %0.2f",bkghist.at(iObs)->GetCorrelationFactor());
      TLatex* correlation = new TLatex(1.5,1.5,corr.Data());
      correlation->Draw("same");
      canvas->RedrawAxis("sameaxis");

      canvas->SaveAs((controlRegion+"_"+obs+"_MC.pdf").c_str());
      canvas->SaveAs((controlRegion+"_"+obs+"_MC.png").c_str());
 
      // plot data
      frame->Draw();
      CMS_lumi(canvas,"2.30",true);
      TGraph2D* graph2D_2 = new TGraph2D (datahist.at(iObs));
      for(int iBinX = 0; iBinX < datahist.at(iObs)->GetNbinsX()+1; iBinX++){
        for(int iBinY = 0; iBinY < datahist.at(iObs)->GetNbinsY()+1; iBinY++){
          graph2D_2->Interpolate(datahist.at(iObs)->GetXaxis()->GetBinCenter(iBinX+1),datahist.at(iObs)->GetXaxis()->GetBinCenter(iBinY+1));
        }
      }
      graph2D_2->GetHistogram()->Draw("colz same");
      profileX = datahist.at(iObs)->ProfileX("_pfx",1,datahist.at(iObs)->GetNbinsX());
      profileX->SetMarkerStyle(20);
      profileX->SetMarkerSize(0.8);
      profileX->SetMarkerColor(kBlack);
      profileX->SetLineColor(kBlack);
      profileX->Draw("EP same");
      corr = Form("#rho = %0.2f",datahist.at(iObs)->GetCorrelationFactor());
      correlation = new TLatex(0.7,0.7,corr.Data());
      correlation->Draw("same");
      canvas->RedrawAxis("sameaxis");

      canvas->SaveAs((controlRegion+"_"+obs+"_Data.pdf").c_str());
      canvas->SaveAs((controlRegion+"_"+obs+"_Data.png").c_str());

    }      
    iObs++;
  }
}
