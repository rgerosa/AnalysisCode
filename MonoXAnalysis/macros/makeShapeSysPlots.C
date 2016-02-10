#include "CMS_lumi.h"
#include "histoUtils.h"

void makeShapeSysPlots(string inputFileName, string controlRegion, string process, string observable, string observableLatex, int category,
		       string MediatorMass, string DMMass){

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);

  TFile* inputFile = TFile::Open(inputFileName.c_str());

  TH1* nominalHist = NULL;
  TH1* hist_bUp = NULL;
  TH1* hist_bDw = NULL;
  TH1* hist_jesUp = NULL;
  TH1* hist_jesDw = NULL;
  TH1* hist_jerUp = NULL;
  TH1* hist_jerDw = NULL;
  TH1* hist_uncUp = NULL;
  TH1* hist_uncDw = NULL;

  if(controlRegion == "ZM") controlRegion = "zmm";
  if(controlRegion == "ZE") controlRegion = "zee";
  if(controlRegion == "WM") controlRegion = "wmn";
  if(controlRegion == "WE") controlRegion = "wen";


  if(process == "Top" and controlRegion != "SR"){
    nominalHist = (TH1*) inputFile->Get(("tbkghist"+controlRegion+"_"+observable).c_str());
    hist_bUp = (TH1*) inputFile->Get(("tbkghist"+controlRegion+"_bUp_"+observable).c_str());
    hist_bDw = (TH1*) inputFile->Get(("tbkghist"+controlRegion+"_bDw_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->Get(("tbkghist"+controlRegion+"_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->Get(("tbkghist"+controlRegion+"_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->Get(("tbkghist"+controlRegion+"_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->Get(("tbkghist"+controlRegion+"_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->Get(("tbkghist"+controlRegion+"_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->Get(("tbkghist"+controlRegion+"_metUncDw_"+observable).c_str());
  }
  if(process == "Top" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->Get(("tbkghist_"+observable).c_str());
    hist_bUp = (TH1*) inputFile->Get(("tbkghist_bUp_"+observable).c_str());
    hist_bDw = (TH1*) inputFile->Get(("tbkghist_bDw_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->Get(("tbkghist_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->Get(("tbkghist_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->Get(("tbkghist_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->Get(("tbkghist_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->Get(("tbkghist_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->Get(("tbkghist_metUncDw_"+observable).c_str());
  }
  else if(process == "Dibosons" and controlRegion != "SR"){
    nominalHist = (TH1*) inputFile->Get(("dbkghist"+controlRegion+"_"+observable).c_str());
    hist_bUp = (TH1*) inputFile->Get(("dbkghist"+controlRegion+"_bUp_"+observable).c_str());
    hist_bDw = (TH1*) inputFile->Get(("dbkghist"+controlRegion+"_bDw_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->Get(("dbkghist"+controlRegion+"_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->Get(("dbkghist"+controlRegion+"_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->Get(("dbkghist"+controlRegion+"_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->Get(("dbkghist"+controlRegion+"_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->Get(("dbkghist"+controlRegion+"_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->Get(("dbkghist"+controlRegion+"_metUncDw_"+observable).c_str());
  }
  else if(process == "Dibosons" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->Get(("dbkghist_"+observable).c_str());
    hist_bUp = (TH1*) inputFile->Get(("dbkghist_bUp_"+observable).c_str());
    hist_bDw = (TH1*) inputFile->Get(("dbkghist_bDw_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->Get(("dbkghist_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->Get(("dbkghist_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->Get(("dbkghist_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->Get(("dbkghist_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->Get(("dbkghist_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->Get(("dbkghist_metUncDw_"+observable).c_str());
  }
  else if(process == "ZJets" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->Get(("zjethist_"+observable).c_str());
    hist_bUp = (TH1*) inputFile->Get(("zjethist_bUp_"+observable).c_str());
    hist_bDw = (TH1*) inputFile->Get(("zjethist_bDw_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->Get(("zjethist_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->Get(("zjethist_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->Get(("zjethist_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->Get(("zjethist_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->Get(("zjethist_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->Get(("zjethist_metUncDw_"+observable).c_str());
  }
  else if(process == "ZJets" and controlRegion != "SR"){
    nominalHist = (TH1*) inputFile->Get(("vllbkghist"+controlRegion+"_"+observable).c_str());
    hist_bUp = (TH1*) inputFile->Get(("vllbkghist"+controlRegion+"_bUp_"+observable).c_str());
    hist_bDw = (TH1*) inputFile->Get(("vllbkghist"+controlRegion+"_bDw_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->Get(("vllbkghist"+controlRegion+"_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->Get(("vllbkghist"+controlRegion+"_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->Get(("vllbkghist"+controlRegion+"_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->Get(("vllbkghist"+controlRegion+"_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->Get(("vllbkghist"+controlRegion+"_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->Get(("vllbkghist"+controlRegion+"_metUncDw_"+observable).c_str());
  }
  else if(process == "WJets" and controlRegion != "SR"){
    nominalHist = (TH1*) inputFile->Get(("vlbkghist"+controlRegion+"_"+observable).c_str());
    hist_bUp = (TH1*) inputFile->Get(("vlbkghist"+controlRegion+"_bUp_"+observable).c_str());
    hist_bDw = (TH1*) inputFile->Get(("vlbkghist"+controlRegion+"_bDw_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->Get(("vlbkghist"+controlRegion+"_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->Get(("vlbkghist"+controlRegion+"_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->Get(("vlbkghist"+controlRegion+"_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->Get(("vlbkghist"+controlRegion+"_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->Get(("vlbkghist"+controlRegion+"_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->Get(("vlbkghist"+controlRegion+"_metUncDw_"+observable).c_str());
  }
  else if(process == "MonoJ" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->Get(("monoJhist_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_bUp    = (TH1*) inputFile->Get(("monoJhist_bUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_bDw    = (TH1*) inputFile->Get(("monoJhist_bDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->Get(("monoJhist_metJetUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->Get(("monoJhist_metJetDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->Get(("monoJhist_metResUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->Get(("monoJhist_metResDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->Get(("monoJhist_metUncUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->Get(("monoJhist_metUncDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
  }
  else if(process == "MonoW" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->Get(("monoWhist_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_bUp    = (TH1*) inputFile->Get(("monoWhist_bUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_bDw    = (TH1*) inputFile->Get(("monoWhist_bDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->Get(("monoWhist_metJetUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->Get(("monoWhist_metJetDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->Get(("monoWhist_metResUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->Get(("monoWhist_metResDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->Get(("monoWhist_metUncUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->Get(("monoWhist_metUncDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
  }
  else if(process == "MonoZ" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->Get(("monoZhist_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_bUp    = (TH1*) inputFile->Get(("monoZhist_bUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_bDw    = (TH1*) inputFile->Get(("monoZhist_bDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->Get(("monoZhist_metJetUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->Get(("monoZhist_metJetDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->Get(("monoZhist_metResUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->Get(("monoZhist_metResDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->Get(("monoZhist_metUncUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->Get(("monoZhist_metUncDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
  }

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetLeftMargin(0.11);
  
  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetTickx();
  pad1->SetTicky();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.28);
  pad2->SetTickx();
  pad2->SetTicky();

  pad1->SetRightMargin(0.075);
  pad1->SetTopMargin(0.06);
  pad1->SetBottomMargin(0.0);
  pad1->Draw();
  pad1->cd();

  // set style
  nominalHist->SetLineColor(kBlack);
  nominalHist->SetMarkerColor(kBlack);
  nominalHist->SetMarkerStyle(20);
  nominalHist->SetMarkerSize(1.);

  nominalHist->Scale(1.,"width");
  hist_bUp->Scale(1.,"width");
  hist_bDw->Scale(1.,"width");
  hist_jesUp->Scale(1.,"width");
  hist_jesDw->Scale(1.,"width");
  hist_jerUp->Scale(1.,"width");
  hist_jerDw->Scale(1.,"width");
  hist_uncUp->Scale(1.,"width");
  hist_uncDw->Scale(1.,"width");

  if(category <= 1){

    if(process == "Top"){
      fixShapeUncertainty(nominalHist,hist_bUp,500.,1.06);
      fixShapeUncertainty(nominalHist,hist_bDw,500.,0.94);
      fixShapeUncertainty(nominalHist,hist_jesUp,500.,1.10);
      fixShapeUncertainty(nominalHist,hist_jesDw,500.,0.90);
      fixShapeUncertainty(nominalHist,hist_jerUp,500.,1.03);
      fixShapeUncertainty(nominalHist,hist_jerDw,500.,0.97);
      fixShapeUncertainty(nominalHist,hist_uncUp,500.,1.01);
      fixShapeUncertainty(nominalHist,hist_uncDw,500.,0.99);
    }
    else{
      fixShapeUncertainty(nominalHist,hist_bUp,500.,1.02);
      fixShapeUncertainty(nominalHist,hist_bDw,500.,0.98);
      fixShapeUncertainty(nominalHist,hist_jesUp,500.,1.06);
      fixShapeUncertainty(nominalHist,hist_jesDw,500.,0.94);
      fixShapeUncertainty(nominalHist,hist_jerUp,500.,1.02);
      fixShapeUncertainty(nominalHist,hist_jerDw,500.,0.98);
      fixShapeUncertainty(nominalHist,hist_uncUp,500.,1.01);
      fixShapeUncertainty(nominalHist,hist_uncDw,500.,0.99);
    }
  }

  hist_bUp->SetLineColor(kRed);
  hist_bDw->SetLineColor(kBlue);

  hist_jesUp->SetLineColor(kRed);
  hist_jesDw->SetLineColor(kBlue);

  hist_jerUp->SetLineColor(kRed);
  hist_jerDw->SetLineColor(kBlue);

  hist_uncUp->SetLineColor(kRed);
  hist_uncDw->SetLineColor(kBlue);

  TH1* frame =  pad1->DrawFrame(nominalHist->GetXaxis()->GetBinLowEdge(1), 1.5e-4, 
				nominalHist->GetXaxis()->GetBinLowEdge(nominalHist->GetNbinsX()+1), nominalHist->GetMaximum()*1000, "");
 
  frame->GetXaxis()->SetTitle(observableLatex.c_str());
  frame->GetYaxis()->SetTitle("Events / GeV");
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetLabelSize(0.);
  frame->GetXaxis()->SetLabelOffset(1.10);
  frame->GetXaxis()->SetTitleSize(0.);
  frame->GetYaxis()->SetTitleSize(0.050);

  frame ->Draw();
  CMS_lumi(pad1, 4, 0, true);
  nominalHist->Draw("P same");
  hist_bUp->Draw("hist same");
  hist_bDw->Draw("hist same");

  TLegend* leg = new TLegend(0.58, 0.62, 0.85, 0.85);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  leg->AddEntry(nominalHist,process.c_str(),"PL");
  leg->AddEntry(hist_bUp,(process+" bUp").c_str(),"L");
  leg->AddEntry(hist_bDw,(process+" bDw").c_str(),"L");
  leg->Draw("same");

  pad1->RedrawAxis("sameaxis");
  pad1->SetLogy();

  canvas->cd();
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->SetRightMargin(0.075);
  pad2->SetGridy();
  pad2->Draw();
  pad2->cd();
  
  TH1* frame2 =  pad2->DrawFrame(nominalHist->GetXaxis()->GetBinLowEdge(1), 0.5, nominalHist->GetXaxis()->GetBinLowEdge(nominalHist->GetNbinsX()+1), 1.5, "");
  frame2->GetXaxis()->SetLabelSize(0.10);
  frame2->GetXaxis()->SetLabelOffset(0.03);
  frame2->GetXaxis()->SetTitleSize(0.13);
  frame2->GetXaxis()->SetTitleOffset(1.05);
  frame2->GetYaxis()->SetLabelSize(0.08);
  frame2->GetYaxis()->SetTitleSize(0.10);
  frame2->GetXaxis()->SetTitle(observableLatex.c_str());
  frame2->GetYaxis()->SetNdivisions(504, false);
  frame2->GetYaxis()->SetTitle("Var./Nom.");
  frame2->GetYaxis()->SetTitleOffset(0.5);
  frame2->Draw();

  
  TH1* ratio_bUp = (TH1*) hist_bUp->Clone("ratio_bUp");
  TH1* ratio_bDw = (TH1*) hist_bDw->Clone("ratio_bDw");
  ratio_bUp->Divide(nominalHist);
  ratio_bDw->Divide(nominalHist);

  TH1* ratio_jesUp = (TH1*) hist_jesUp->Clone("ratio_jesUp");
  TH1* ratio_jesDw = (TH1*) hist_jesDw->Clone("ratio_jesDw");
  ratio_jesUp->Divide(nominalHist);
  ratio_jesDw->Divide(nominalHist);

  TH1* ratio_jerUp = (TH1*) hist_jerUp->Clone("ratio_jerUp");
  TH1* ratio_jerDw = (TH1*) hist_jerDw->Clone("ratio_jerDw");
  ratio_jerUp->Divide(nominalHist);
  ratio_jerDw->Divide(nominalHist);

  TH1* ratio_uncUp = (TH1*) hist_uncUp->Clone("ratio_uncUp");
  TH1* ratio_uncDw = (TH1*) hist_uncDw->Clone("ratio_uncDw");
  ratio_uncUp->Divide(nominalHist);
  ratio_uncDw->Divide(nominalHist);

  ratio_bUp->SetStats(kFALSE);
  ratio_bDw->SetStats(kFALSE);
  ratio_bUp->Draw("hist same");
  ratio_bUp->Draw("PE same");
  ratio_bDw->Draw("hist same");
  ratio_bDw->Draw("PE same");

  pad2->RedrawAxis("sameaxis");

  canvas->SaveAs((observable+"_"+controlRegion+"_"+process+"_btag.png").c_str());
  canvas->SaveAs((observable+"_"+controlRegion+"_"+process+"_btag.pdf").c_str());

  //
  leg->Clear();
  pad1->cd();
  frame->Draw();
  nominalHist->Draw("P same");
  hist_jesUp->Draw("hist same");
  hist_jesDw->Draw("hist same");

  leg->AddEntry(nominalHist,process.c_str(),"PL");
  leg->AddEntry(hist_bUp,(process+" jesUp").c_str(),"L");
  leg->AddEntry(hist_bDw,(process+" jesDw").c_str(),"L");
  leg->Draw("same");

  pad1->RedrawAxis("sameaxis");

  pad2->cd();
  frame2->Draw();

  ratio_jesUp->SetStats(kFALSE);
  ratio_jesDw->SetStats(kFALSE);
  ratio_jesUp->Draw("hist same");
  ratio_jesUp->Draw("PE same");
  ratio_jesDw->Draw("hist same");
  ratio_jesDw->Draw("PE same");

  pad2->RedrawAxis("sameaxis");

  canvas->SaveAs((observable+"_"+controlRegion+"_"+process+"_jes.png").c_str());
  canvas->SaveAs((observable+"_"+controlRegion+"_"+process+"_jes.pdf").c_str());
  
  ///
  leg->Clear();
  pad1->cd();
  frame->Draw();
  nominalHist->Draw("P same");
  hist_jerUp->Draw("hist same");
  hist_jerDw->Draw("hist same");

  leg->AddEntry(nominalHist,process.c_str(),"PL");
  leg->AddEntry(hist_bUp,(process+" jerUp").c_str(),"L");
  leg->AddEntry(hist_bDw,(process+" jerDw").c_str(),"L");
  leg->Draw("same");

  pad1->RedrawAxis("sameaxis");

  pad2->cd();
  frame2->Draw();

  ratio_jerUp->SetStats(kFALSE);
  ratio_jerDw->SetStats(kFALSE);
  ratio_jerUp->Draw("hist same");
  ratio_jerUp->Draw("PE same");
  ratio_jerDw->Draw("hist same");
  ratio_jerDw->Draw("PE same");

  pad2->RedrawAxis("sameaxis");

  canvas->SaveAs((observable+"_"+controlRegion+"_"+process+"_jer.png").c_str());
  canvas->SaveAs((observable+"_"+controlRegion+"_"+process+"_jer.pdf").c_str());
  
  ///
  leg->Clear();
  pad1->cd();
  frame->Draw();
  nominalHist->Draw("P same");
  hist_uncUp->Draw("hist same");
  hist_uncDw->Draw("hist same");

  leg->AddEntry(nominalHist,process.c_str(),"PL");
  leg->AddEntry(hist_bUp,(process+" uncUp").c_str(),"L");
  leg->AddEntry(hist_bDw,(process+" uncDw").c_str(),"L");
  leg->Draw("same");

  pad1->RedrawAxis("sameaxis");

  pad2->cd();
  frame2->Draw();

  ratio_uncUp->SetStats(kFALSE);
  ratio_uncDw->SetStats(kFALSE);
  ratio_uncUp->Draw("hist same");
  ratio_uncUp->Draw("PE same");
  ratio_uncDw->Draw("hist same");
  ratio_uncDw->Draw("PE same");

  pad2->RedrawAxis("sameaxis");

  canvas->SaveAs((observable+"_"+controlRegion+"_"+process+"_unc.png").c_str());
  canvas->SaveAs((observable+"_"+controlRegion+"_"+process+"_unc.pdf").c_str());
  

}
