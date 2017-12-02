#include "../makeTemplates/histoUtils.h"

void makeModifyTemplateName(string inputFileName, string outputFileName, Sample sample){

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");

  outputFile->cd();
  if(sample == Sample::sig){
    outputFile->mkdir("shapes_prefit");
    outputFile->mkdir("shapes_prefit/ch1");
    outputFile->cd("shapes_prefit/ch1");

    TH1* hist = ((TH1F*) inputFile->Get("SR/zinvhist_mjj"));
    hist->Scale(1,"width");
    hist->Write("Znunu");

    hist = ((TH1F*) inputFile->Get("SR/wjethist_mjj"));
    hist->Scale(1,"width");
    hist->Write("WJets");

    hist = ((TH1F*) inputFile->Get("SR/ewkbkgzhist_mjj"));
    hist->Scale(1,"width");
    hist->Write("Znunu_EWK");

    hist = ((TH1F*) inputFile->Get("SR/ewkbkgwhist_mjj"));
    hist->Scale(1,"width");
    hist->Write("WJets_EWK");

    hist = ((TH1F*) inputFile->Get("SR/tbkghist_mjj"));
    hist->Scale(1,"width");
    hist->Write("Top");

    hist = ((TH1F*) inputFile->Get("SR/dbkghist_mjj"));
    hist->Scale(1,"width");
    hist->Write("Diboson");
    
  }
  else if(sample == Sample::zmm){
    outputFile->mkdir("shapes_prefit");
    outputFile->mkdir("shapes_prefit/ch2");
    outputFile->cd("shapes_prefit/ch2");

    TH1* hist = ((TH1F*) inputFile->Get("ZM/vllbkghistzmm_mjj"));
    hist->Scale(1,"width");
    hist->Write("Znunu");

    hist = ((TH1F*) inputFile->Get("ZM/ewkzbkghistzmm_mjj"));
    hist->Scale(1,"width");
    hist->Write("Znunu_EWK");

    hist = ((TH1F*) inputFile->Get("ZM/tbkghistzmm_mjj"));
    hist->Scale(1,"width");
    hist->Write("Top");

    hist = ((TH1F*) inputFile->Get("ZM/dbkghistzmm_mjj"));
    hist->Scale(1,"width");
    hist->Write("Diboson");

  }
  else if(sample == Sample::wmn){
    outputFile->mkdir("shapes_prefit");
    outputFile->mkdir("shapes_prefit/ch3");
    outputFile->cd("shapes_prefit/ch3");

    TH1* hist = ((TH1F*) inputFile->Get("WM/vlbkghistwmn_mjj"));
    hist->Scale(1,"width");
    hist->Write("WJets");

    hist = ((TH1F*) inputFile->Get("WM/ewkwbkghistwmn_mjj"));
    hist->Scale(1,"width");
    hist->Write("WJets_EWK");

    hist = ((TH1F*) inputFile->Get("WM/tbkghistwmn_mjj"));
    hist->Scale(1,"width");
    hist->Write("Top");

    hist = ((TH1F*) inputFile->Get("WM/dbkghistwmn_mjj"));
    hist->Scale(1,"width");
    hist->Write("Diboson");

    hist = ((TH1F*) inputFile->Get("WM/qbkghistwmn_mjj"));
    hist->Scale(1,"width");
    hist->Write("QCD");


  }
  else if(sample == Sample::wen){

    outputFile->mkdir("shapes_prefit");
    outputFile->mkdir("shapes_prefit/ch5");
    outputFile->cd("shapes_prefit/ch5");

    TH1* hist = ((TH1F*) inputFile->Get("WE/vlbkghistwen_mjj"));
    hist->Scale(1,"width");
    hist->Write("WJets");

    hist = ((TH1F*) inputFile->Get("WE/ewkwbkghistwen_mjj"));
    hist->Scale(1,"width");
    hist->Write("WJets_EWK");

    hist = ((TH1F*) inputFile->Get("WE/tbkghistwen_mjj"));
    hist->Scale(1,"width");
    hist->Write("Top");

    hist = ((TH1F*) inputFile->Get("WE/dbkghistwen_mjj"));
    hist->Scale(1,"width");
    hist->Write("Diboson");

    hist = ((TH1F*) inputFile->Get("WE/qbkghistwen_mjj"));
    hist->Scale(1,"width");
    hist->Write("QCD");

  }
  else if(sample == Sample::zee){

    outputFile->mkdir("shapes_prefit");
    outputFile->mkdir("shapes_prefit/ch4");
    outputFile->cd("shapes_prefit/ch4");

    TH1* hist = ((TH1F*) inputFile->Get("ZE/vllbkghistzee_mjj"));
    hist->Scale(1,"width");
    hist->Write("Znunu");

    hist = ((TH1F*) inputFile->Get("ZE/ewkzbkghistzee_mjj"));
    hist->Scale(1,"width");
    hist->Write("Znunu_EWK");

    hist = ((TH1F*) inputFile->Get("ZE/tbkghistzee_mjj"));
    hist->Scale(1,"width");
    hist->Write("Top");

    hist = ((TH1F*) inputFile->Get("ZE/dbkghistzee_mjj"));
    hist->Scale(1,"width");
    hist->Write("Diboson");

  }

  outputFile->Close();

}
