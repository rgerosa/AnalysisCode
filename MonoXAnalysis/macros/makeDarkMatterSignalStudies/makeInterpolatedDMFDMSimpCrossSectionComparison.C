#include "../CMS_lumi.h"

/////////////////////////
int mmed(double mh, int code){
  if (code == 800) return ((int)(mh-80000000000))/10000;
  if (code == 801) return ((int)(mh-80100000000))/10000;
  if (code == 805) return ((int)(mh-80500000000))/10000;
  if (code == 806) return ((int)(mh-80600000000))/10000;
  return -1;
}

////////////////////////
int mdm(double mh, int code){
  if (code == 800) return (mh-80000000000)  - ( ((Int_t)(mh-80000000000))/10000 )*10000;
  if (code == 801) return (mh-80100000000)  - ( ((Int_t)(mh-80100000000))/10000 )*10000;
  if (code == 805) return (mh-80500000000)  - ( ((Int_t)(mh-80500000000))/10000 )*10000;
  if (code == 806) return (mh-80600000000)  - ( ((Int_t)(mh-80600000000))/10000 )*10000;
  return -1;
}

///////////////////////
int code(double mh){
  return (int)(mh/100000000);
}

///////////////////////
static float scale_dmf    = 2.78;
static float scale_dmsimp = 1.;


///////////////////////
void fillMassPoint(TDirectory* dmf_dir, TDirectory* dmsimp_dir, TList* dmf_masspoint, TList* dmsimp_masspoint, vector<TH1F*> & histo_dmf, vector<TH1F*>& histo_dmsimp, const float & dmMass){

  if(dmf_masspoint == NULL) return;
  if(dmsimp_masspoint == NULL) return;

  ////////////////
  for(auto obj : *dmf_masspoint){
    TString name = Form("%s",obj->GetName());
    if(not name.Contains("signal_signal")) continue;
    name.ReplaceAll("signal_signal_","");
    double mh   = stod(name.Data());
    int c       = code(mh);
    int medmass = mmed(mh,c);
    int dmmass  = mdm(mh,c);

    if((c == 800 or c == 801) and medmass < 500) continue;
    if((c == 800 or c == 801) and medmass > 3000) continue;
    if(dmmass != dmMass) continue;
    
    histo_dmf.push_back((TH1F*) dmf_dir->Get(obj->GetName()));        
    histo_dmf.back()->Scale(scale_dmf);
  }

  //////////////// 
  for(auto obj : *dmsimp_masspoint){
    TString name = Form("%s",obj->GetName());
    if(not name.Contains("signal_signal")) continue;
    name.ReplaceAll("signal_signal_","");

    double mh   = stod(name.Data());
    int c       = code(mh);
    int medmass = mmed(mh,c);
    int dmmass  = mdm(mh,c);

    if((c == 800 or c == 801) and medmass < 500) continue;
    if((c == 800 or c == 801) and medmass > 3000) continue;
    if(dmmass != dmMass) continue;

    histo_dmsimp.push_back((TH1F*) dmsimp_dir->Get(obj->GetName()));        
    histo_dmsimp.back()->Scale(scale_dmsimp);
  }

  return;

}

/////////////
void fillGraph(TGraphErrors* ratio, const vector<TH1F*> & histo_dmf, const vector<TH1F*> & histo_dmsimp, int & ipoint){
  
  for(auto dmfhist : histo_dmf){
    TString name(dmfhist->GetName());
    TH1F* dmsimphist = NULL;
    for(size_t ihist = 0; ihist < histo_dmsimp.size(); ihist++){
      if(TString(histo_dmsimp.at(ihist)->GetName()) == name){// found common point
	dmsimphist = histo_dmsimp.at(ihist);
	break;
      }
    }
    if(dmsimphist != NULL){
      if(not name.Contains("signal_signal")) continue;
      name.ReplaceAll("signal_signal_","");
      double mh = stod(name.Data());
      int c     = code(mh);
      int medmass = mmed(mh,c);
      double dmsimp_error = 0;
      double dmsimp_intergral = dmsimphist->IntegralAndError(0,dmsimphist->GetNbinsX(),dmsimp_error);
      double dmf_error = 0;
      double dmf_intergral = dmfhist->IntegralAndError(0,dmfhist->GetNbinsX(),dmf_error);
      ratio->SetPoint(ipoint,medmass,dmsimp_intergral/dmf_intergral);
      ratio->SetPointError(ipoint,0,sqrt(dmsimp_error*dmsimp_error/(dmf_intergral*dmf_intergral)+dmf_error*dmf_error*(dmsimp_intergral*dmsimp_intergral)/(pow(dmf_intergral,4))));      
      ipoint++;
    }      
  }
  
  return;
}
  

///////////////////////
void makeInterpolatedDMFDMSimpCrossSectionComparison(int modelcode, string outputDIR, string postfix, float dmMass){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  //////////////////
  TFile* dmf_monojet_file = NULL;
  TFile* dmf_monow_file   = NULL;
  TFile* dmf_monoz_file   = NULL;

  if(modelcode == 800){ // vector model
    dmf_monojet_file = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMF/MonoJ_800_0.25_catmonojet_13TeV_v1.root","READ");
    dmf_monow_file   = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMF/MonoW_800_0.25_catmonojet_13TeV_v1.root","READ");
    dmf_monoz_file   = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMF/MonoZ_800_0.25_catmonojet_13TeV_v1.root","READ");
  }
  else if(modelcode == 801){
    dmf_monojet_file = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMF/MonoJ_801_0.25_catmonojet_13TeV_v1.root","READ");
    dmf_monow_file   = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMF/MonoW_801_0.25_catmonojet_13TeV_v1.root","READ");
    dmf_monoz_file   = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMF/MonoZ_801_0.25_catmonojet_13TeV_v1.root","READ");
    
  }
  else if(modelcode == 805){
    dmf_monojet_file = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMF/MonoJ_805_1.0_catmonojet_13TeV_v1.root","READ");
    //    dmf_monoz_file   = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMF/MonoZ_805_1.0_catmonojet_13TeV_v1.root","READ");
  }
  else if(modelcode == 806){
    dmf_monojet_file = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMF/MonoJ_806_1.0_catmonojet_13TeV_v1.root","READ");
  }

  TFile* dmsimp_monojet_file = NULL;
  TFile* dmsimp_monow_file   = NULL;
  TFile* dmsimp_monoz_file   = NULL;

  if(modelcode == 800){ // vector model
    dmsimp_monojet_file = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoJ_800_0.25_catmonojet_13TeV_v1.root","READ");
    dmsimp_monow_file   = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoW_800_0.25_catmonojet_13TeV_v1.root","READ");
    dmsimp_monoz_file   = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoZ_800_0.25_catmonojet_13TeV_v1.root","READ");
  }
  else if(modelcode == 801){
    dmsimp_monojet_file = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoJ_801_0.25_catmonojet_13TeV_v1.root","READ");
    dmsimp_monow_file   = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoW_801_0.25_catmonojet_13TeV_v1.root","READ");
    //    dmsimp_monoz_file   = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoZ_801_0.25_catmonojet_13TeV_v1.root","READ");    
  }
  else if(modelcode == 805){
    dmsimp_monojet_file = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoJ_805_1.0_catmonojet_13TeV_v1.root","READ");
    dmsimp_monoz_file   = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoZ_805_1.0_catmonojet_13TeV_v1.root","READ");
  }
  else if(modelcode == 806){
    dmsimp_monojet_file = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS_2016_Data/SignalTemplatesForLimit/Moriond_2016/DMSimp/MonoJ_806_1.0_catmonojet_13TeV_v1.root","READ");
  }
    
  
  ////////
  TDirectory* dmf_monojet_dir = NULL;
  TDirectory* dmf_monow_dir = NULL;
  TDirectory* dmf_monoz_dir = NULL;

  if(dmf_monojet_file)
    dmf_monojet_dir = (TDirectory*) dmf_monojet_file->Get("category_monojet");
  if(dmf_monow_file)
    dmf_monow_dir = (TDirectory*) dmf_monow_file->Get("category_monojet");
  if(dmf_monoz_file)
    dmf_monoz_dir = (TDirectory*) dmf_monoz_file->Get("category_monojet");

  TDirectory* dmsimp_monojet_dir = NULL;
  TDirectory* dmsimp_monow_dir = NULL;
  TDirectory* dmsimp_monoz_dir = NULL;

  if(dmsimp_monojet_file)
    dmsimp_monojet_dir = (TDirectory*) dmsimp_monojet_file->Get("category_monojet");
  if(dmsimp_monow_file)
    dmsimp_monow_dir = (TDirectory*) dmsimp_monow_file->Get("category_monojet");
  if(dmsimp_monoz_file)
    dmsimp_monoz_dir = (TDirectory*) dmsimp_monoz_file->Get("category_monojet");

  //////
  vector<TH1F*> histo_dmf_monojet;
  vector<TH1F*> histo_dmsimp_monojet;
  vector<TH1F*> histo_dmf_monow;
  vector<TH1F*> histo_dmsimp_monow;
  vector<TH1F*> histo_dmf_monoz;
  vector<TH1F*> histo_dmsimp_monoz;

  TList* dmf_monojet_masspoint = NULL;
  TList* dmf_monow_masspoint = NULL;
  TList* dmf_monoz_masspoint = NULL;
  
  if(dmf_monojet_dir) 
    dmf_monojet_masspoint = dmf_monojet_dir->GetListOfKeys();
  if(dmf_monow_dir) 
    dmf_monow_masspoint = dmf_monow_dir->GetListOfKeys();
  if(dmf_monoz_dir) 
    dmf_monoz_masspoint = dmf_monoz_dir->GetListOfKeys();

  TList* dmsimp_monojet_masspoint = NULL;
  TList* dmsimp_monow_masspoint = NULL;
  TList* dmsimp_monoz_masspoint = NULL;
  
  if(dmsimp_monojet_dir) 
    dmsimp_monojet_masspoint = dmsimp_monojet_dir->GetListOfKeys();
  if(dmsimp_monow_dir) 
    dmsimp_monow_masspoint = dmsimp_monow_dir->GetListOfKeys();
  if(dmsimp_monoz_dir) 
    dmsimp_monoz_masspoint = dmsimp_monoz_dir->GetListOfKeys();

  ///////
  fillMassPoint(dmf_monojet_dir,dmsimp_monojet_dir,dmf_monojet_masspoint,dmsimp_monojet_masspoint,histo_dmf_monojet,histo_dmsimp_monojet,dmMass);
  fillMassPoint(dmf_monow_dir,dmsimp_monow_dir,dmf_monow_masspoint,dmsimp_monow_masspoint,histo_dmf_monow,histo_dmsimp_monow,dmMass);
  fillMassPoint(dmf_monoz_dir,dmsimp_monoz_dir,dmf_monoz_masspoint,dmsimp_monoz_masspoint,histo_dmf_monoz,histo_dmsimp_monoz,dmMass);

  ////////////////////////////////
  TGraphErrors* ratio_monojet = new TGraphErrors();
  TGraphErrors* ratio_monow   = new TGraphErrors();
  TGraphErrors* ratio_monoz   = new TGraphErrors();

  int ipoint_monojet = 0;
  int ipoint_monow = 0;
  int ipoint_monoz = 0;

  fillGraph(ratio_monojet,histo_dmf_monojet,histo_dmsimp_monojet,ipoint_monojet);
  fillGraph(ratio_monow,histo_dmf_monow,histo_dmsimp_monow,ipoint_monow);
  fillGraph(ratio_monoz,histo_dmf_monoz,histo_dmsimp_monoz,ipoint_monoz);

  //////////////////////////// --> make final plot
  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();
    
  ratio_monojet->GetXaxis()->SetTitle("Mediator mass [GeV]");
  ratio_monojet->GetYaxis()->SetTitle("DMSimp/DMF");
  ratio_monojet->SetMarkerColor(kBlack);
  ratio_monojet->SetLineColor(kBlack);
  ratio_monojet->SetLineWidth(2);
  ratio_monojet->SetMarkerStyle(20);

  double min_monojet = TMath::MinElement(ratio_monojet->GetN(),ratio_monojet->GetY());
  double max_monojet = TMath::MaxElement(ratio_monojet->GetN(),ratio_monojet->GetY());
  double min_monow = 0;
  double max_monow = 0;
  if(ipoint_monow != 0){
    min_monow = TMath::MinElement(ratio_monow->GetN(),ratio_monow->GetY());
    max_monow = TMath::MaxElement(ratio_monow->GetN(),ratio_monow->GetY());
  }
  double min_monoz = 0;
  double max_monoz = 0;
  if(ipoint_monoz != 0){
    min_monoz = TMath::MinElement(ratio_monoz->GetN(),ratio_monoz->GetY());
    max_monoz = TMath::MaxElement(ratio_monoz->GetN(),ratio_monoz->GetY());    
  }

  if(ipoint_monow != 0 and ipoint_monoz != 0)
    ratio_monojet->GetYaxis()->SetRangeUser(min(min_monojet,min(min_monow,min_monoz))*0.75,max(max_monojet,max(max_monow,max_monoz))*1.25); 
  else if(ipoint_monow == 0 and ipoint_monoz != 0)
    ratio_monojet->GetYaxis()->SetRangeUser(min(min_monojet,min_monoz)*0.75,max(max_monojet,max_monoz)*1.25); 
  else if(ipoint_monoz == 0 and ipoint_monow != 0)
    ratio_monojet->GetYaxis()->SetRangeUser(min(min_monojet,min_monoz)*0.75,max(max_monojet,max_monoz)*1.25); 
  else if(ipoint_monoz == 0 and ipoint_monow != 0)
    ratio_monojet->GetYaxis()->SetRangeUser(min(min_monojet,min_monow)*0.75,max(max_monojet,max_monow)*1.25); 
  else if(ipoint_monoz == 0 and ipoint_monow == 0)
    ratio_monojet->GetYaxis()->SetRangeUser(min_monojet*0.75,max_monojet*1.25); 


  ratio_monojet->Draw("AEPLsame");

  if(ipoint_monow != 0){
    ratio_monow->SetMarkerColor(kBlue);
    ratio_monow->SetLineColor(kBlue);
    ratio_monow->SetLineWidth(2);
    ratio_monow->SetMarkerStyle(20);
    ratio_monow->Draw("EPLsame");
  }

  if(ipoint_monoz != 0){
    ratio_monoz->SetMarkerColor(kRed);
    ratio_monoz->SetLineColor(kRed);
    ratio_monoz->SetLineWidth(2);
    ratio_monoz->SetMarkerStyle(20);
    ratio_monoz->Draw("EPLsame");
  }

  TLegend leg (0.6,0.7,0.9,0.92);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(ratio_monojet,Form("mono-jet, m_{DM}=%d GeV",int(dmMass)),"PL");
  if(ipoint_monow != 0)
    leg.AddEntry(ratio_monow,Form("mono-W, m_{DM}=%d GeV",int(dmMass)),"PL");
  if(ipoint_monoz != 0)
    leg.AddEntry(ratio_monoz,Form("mono-Z, m_{DM}=%d GeV",int(dmMass)),"PL");
  leg.Draw("same");


  CMS_lumi(canvas,"35.9");
  canvas->SaveAs((outputDIR+"/dmsimp_vs_dmf_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/dmsimp_vs_dmf_"+postfix+".pdf").c_str(),"pdf");
}
