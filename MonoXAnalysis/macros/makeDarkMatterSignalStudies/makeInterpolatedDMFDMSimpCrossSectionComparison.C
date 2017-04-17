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
void makeInterpolatedDMFDMSimpCrossSectionComparison(string file_DMF, string file_DMSimp, string outputDIR, string postfix, float dmMass){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  //////////////////
  TFile* dmf_file = TFile::Open(file_DMF.c_str(),"READ");
  TFile* dmsimp_file = TFile::Open(file_DMSimp.c_str(),"READ");

  TDirectory* dmf_dir = (TDirectory*) dmf_file->Get("category_monojet");
  TDirectory* dmsimp_dir = (TDirectory*) dmsimp_file->Get("category_monojet");

  vector<TH1F*> histo_dmf;
  vector<TH1F*> histo_dmsimp;

  TList* dmf_masspoint = dmf_dir->GetListOfKeys();
  TList* dmsimp_masspoint = dmsimp_dir->GetListOfKeys();

  ////////////////
  for(auto obj : *dmf_masspoint){
    TString name = Form("%s",obj->GetName());
    if(not name.Contains("signal_signal")) continue;
    name.ReplaceAll("signal_signal_","");
    double mh = stod(name.Data());
    int c     = code(mh);
    int medmass = mmed(mh,c);
    int dmmass  = mdm(mh,c);

    if(dmmass != dmMass) continue;

    histo_dmf.push_back((TH1F*) dmf_dir->Get(obj->GetName()));        
    histo_dmf.back()->Scale(scale_dmf);
  }

  //////////////// 
  for(auto obj : *dmsimp_masspoint){
    TString name = Form("%s",obj->GetName());
    if(not name.Contains("signal_signal")) continue;
    name.ReplaceAll("signal_signal_","");

    double mh = stod(name.Data());
    int c     = code(mh);
    int medmass = mmed(mh,c);
    int dmmass  = mdm(mh,c);

    if(dmmass != dmMass) continue;
    histo_dmsimp.push_back((TH1F*) dmsimp_dir->Get(obj->GetName()));        
    histo_dmsimp.back()->Scale(scale_dmsimp);
  }

  ////////////////////////////////
  TGraph* ratio = new TGraph();
  int ipoint = 0;
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
      ratio->SetPoint(ipoint,medmass,dmsimphist->Integral()/dmfhist->Integral());
      cout<<dmfhist->GetName()<<" "<<dmfhist->Integral()<<" "<<dmsimphist->GetName()<<" "<<dmsimphist->Integral()<<endl;
      ipoint++;
    }      
  }
  
  ////////////////////////////
  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();
    
  ratio->GetXaxis()->SetTitle("Mediator mass [GeV]");
  ratio->GetYaxis()->SetTitle("DMSimp/DMF");
  ratio->SetMarkerColor(kBlack);
  ratio->SetLineColor(kBlack);
  ratio->SetLineWidth(2);
  ratio->SetMarkerStyle(20);
  ratio->GetYaxis()->SetRangeUser(TMath::MinElement(ratio->GetN(),ratio->GetY()),TMath::MaxElement(ratio->GetN(),ratio->GetY()));
  ratio->Draw("AEPLsame");

  TLegend leg (0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(ratio,"mono-jet, m_{DM}=1 GeV","PL");
  leg.Draw("same");


  CMS_lumi(canvas,"35.9");

  canvas->SaveAs((outputDIR+"/dmsimp_vs_dmf_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/dmsimp_vs_dmf_"+postfix+".pdf").c_str(),"pdf");
    
}
