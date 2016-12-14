#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"
#include "../makeTemplates/histoUtils2D.h"

void makeShapeSysPlots(string inputFileName,  // template file
		       Sample sample, // identify the control region
		       string process, // identify the process
		       string observable, 
		       string observableLatex,
		       string MediatorMass = "", 
		       string DMMass = ""){
  
  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);

  TFile* inputFile = TFile::Open(inputFileName.c_str());

  TH1* nominalHist = NULL;
  TH1* hist_jesUp = NULL;
  TH1* hist_jesDw = NULL;
  TH1* hist_jerUp = NULL;
  TH1* hist_jerDw = NULL;
  TH1* hist_uncUp = NULL;
  TH1* hist_uncDw = NULL;
  
  string controlRegion;
  if(sample == Sample::zmm) controlRegion = "zmm";
  else if(sample == Sample::zee) controlRegion = "zee";
  else if(sample == Sample::wmn) controlRegion = "wmn";
  else if(sample == Sample::wen) controlRegion = "wen";
  else if(sample == Sample::gam) controlRegion = "gam";
  else if(sample == Sample::sig) controlRegion = "SR";
  

  if(process == "Top" and controlRegion != "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("tbkghist"+controlRegion+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("tbkghist"+controlRegion+"_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("tbkghist"+controlRegion+"_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("tbkghist"+controlRegion+"_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("tbkghist"+controlRegion+"_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("tbkghist"+controlRegion+"_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("tbkghist"+controlRegion+"_metUncDw_"+observable).c_str());
  }
  if(process == "Top" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("tbkghist_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("tbkghist_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("tbkghist_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("tbkghist_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("tbkghist_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("tbkghist_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("tbkghist_metUncDw_"+observable).c_str());
  }
  else if(process == "Dibosons" and controlRegion != "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("dbkghist"+controlRegion+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("dbkghist"+controlRegion+"_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("dbkghist"+controlRegion+"_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("dbkghist"+controlRegion+"_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("dbkghist"+controlRegion+"_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("dbkghist"+controlRegion+"_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("dbkghist"+controlRegion+"_metUncDw_"+observable).c_str());
  }
  else if(process == "Dibosons" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("dbkghist_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("dbkghist_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("dbkghist_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("dbkghist_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("dbkghist_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("dbkghist_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("dbkghist_metUncDw_"+observable).c_str());
  }
  else if(process == "ZJets" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("zjethist_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("zjethist_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("zjethist_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("zjethist_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("zjethist_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("zjethist_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("zjethist_metUncDw_"+observable).c_str());
  }
  else if(process == "ZJets" and controlRegion != "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("vllbkghist"+controlRegion+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("vllbkghist"+controlRegion+"_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("vllbkghist"+controlRegion+"_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("vllbkghist"+controlRegion+"_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("vllbkghist"+controlRegion+"_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("vllbkghist"+controlRegion+"_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("vllbkghist"+controlRegion+"_metUncDw_"+observable).c_str());
  }
  else if(process == "WJets" and controlRegion != "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("vlbkghist"+controlRegion+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("vlbkghist"+controlRegion+"_metJetUp_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("vlbkghist"+controlRegion+"_metJetDw_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("vlbkghist"+controlRegion+"_metResUp_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("vlbkghist"+controlRegion+"_metResDw_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("vlbkghist"+controlRegion+"_metUncUp_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("vlbkghist"+controlRegion+"_metUncDw_"+observable).c_str());
  }
  else if(process == "MonoJ" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("monoJhist_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("monoJhist_metJetUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("monoJhist_metJetDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("monoJhist_metResUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("monoJhist_metResDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("monoJhist_metUncUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("monoJhist_metUncDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
  }
  else if(process == "MonoW" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("monoWhist_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("monoWhist_metJetUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("monoWhist_metJetDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("monoWhist_metResUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("monoWhist_metResDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("monoWhist_metUncUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("monoWhist_metUncDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
  }
  else if(process == "MonoZ" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("monoZhist_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("monoZhist_metJetUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("monoZhist_metJetDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("monoZhist_metResUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("monoZhist_metResDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("monoZhist_metUncUp_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("monoZhist_metUncDw_"+MediatorMass+"_"+DMMass+"_"+observable).c_str());
  }
  else if(process == "ggH" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("ggHhist_"+MediatorMass+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("ggHhist_metJetUp_"+MediatorMass+"_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("ggHhist_metJetDw_"+MediatorMass+"_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("ggHhist_metResUp_"+MediatorMass+"_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("ggHhist_metResDw_"+MediatorMass+"_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("ggHhist_metUncUp_"+MediatorMass+"_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("ggHhist_metUncDw_"+MediatorMass+"_"+observable).c_str());
  }
  else if(process == "vbfH" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("vbfHhist_"+MediatorMass+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("vbfHhist_metJetUp_"+MediatorMass+"_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("vbfHhist_metJetDw_"+MediatorMass+"_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("vbfHhist_metResUp_"+MediatorMass+"_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("vbfHhist_metResDw_"+MediatorMass+"_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("vbfHhist_metUncUp_"+MediatorMass+"_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("vbfHhist_metUncDw_"+MediatorMass+"_"+observable).c_str());
  }
  else if(process == "wH" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("wHhist_"+MediatorMass+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("wHhist_metJetUp_"+MediatorMass+"_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("wHhist_metJetDw_"+MediatorMass+"_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("wHhist_metResUp_"+MediatorMass+"_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("wHhist_metResDw_"+MediatorMass+"_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("wHhist_metUncUp_"+MediatorMass+"_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("wHhist_metUncDw_"+MediatorMass+"_"+observable).c_str());
  }
  else if(process == "zH" and controlRegion == "SR"){
    nominalHist = (TH1*) inputFile->FindObjectAny(("zHhist_"+MediatorMass+"_"+observable).c_str());
    hist_jesUp = (TH1*) inputFile->FindObjectAny(("zHhist_metJetUp_"+MediatorMass+"_"+observable).c_str());
    hist_jesDw = (TH1*) inputFile->FindObjectAny(("zHhist_metJetDw_"+MediatorMass+"_"+observable).c_str());
    hist_jerUp = (TH1*) inputFile->FindObjectAny(("zHhist_metResUp_"+MediatorMass+"_"+observable).c_str());
    hist_jerDw = (TH1*) inputFile->FindObjectAny(("zHhist_metResDw_"+MediatorMass+"_"+observable).c_str());
    hist_uncUp = (TH1*) inputFile->FindObjectAny(("zHhist_metUncUp_"+MediatorMass+"_"+observable).c_str());
    hist_uncDw = (TH1*) inputFile->FindObjectAny(("zHhist_metUncDw_"+MediatorMass+"_"+observable).c_str());
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
  hist_jesUp->Scale(1.,"width");
  hist_jesDw->Scale(1.,"width");
  hist_jerUp->Scale(1.,"width");
  hist_jerDw->Scale(1.,"width");
  hist_uncUp->Scale(1.,"width");
  hist_uncDw->Scale(1.,"width");

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
  CMS_lumi(pad1,"36.2",true);
  nominalHist->Draw("P same");

  TLegend* leg = new TLegend(0.58, 0.62, 0.85, 0.85);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

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

  //
  leg->Clear();
  pad1->cd();
  frame->Draw();
  nominalHist->Draw("P same");
  hist_jesUp->Draw("hist same");
  hist_jesDw->Draw("hist same");

  leg->AddEntry(nominalHist,process.c_str(),"PL");
  leg->AddEntry(hist_jesUp,(process+" jesUp").c_str(),"L");
  leg->AddEntry(hist_jesDw,(process+" jesDw").c_str(),"L");
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
  leg->AddEntry(hist_jerUp,(process+" jerUp").c_str(),"L");
  leg->AddEntry(hist_jerDw,(process+" jerDw").c_str(),"L");
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
  leg->AddEntry(hist_uncUp,(process+" uncUp").c_str(),"L");
  leg->AddEntry(hist_uncDw,(process+" uncDw").c_str(),"L");
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
