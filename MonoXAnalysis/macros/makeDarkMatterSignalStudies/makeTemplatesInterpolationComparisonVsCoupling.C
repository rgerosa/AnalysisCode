#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void checkNegativeBin(TH1* histo){
  for(int iBin = 0; iBin < histo->GetNbinsX()+1; iBin++){
    if(histo->GetBinContent(iBin+1) <= 0)
      histo->SetBinContent(iBin+1,0.0001);
  }
}

int mmed(double mh, int code){

  if (code == 800) return ((int)(mh-80000000000))/10000;
  if (code == 801) return ((int)(mh-80100000000))/10000;
  if (code == 805) return ((int)(mh-80500000000))/10000;
  if (code == 806) return ((int)(mh-80600000000))/10000;
  return -1;
}

 
int mdm(double mh, int code){
  if (code == 800) return (mh-80000000000)  - ( ((Int_t)(mh-80000000000))/10000 )*10000;
  if (code == 801) return (mh-80100000000)  - ( ((Int_t)(mh-80100000000))/10000 )*10000;
  if (code == 805) return (mh-80500000000)  - ( ((Int_t)(mh-80500000000))/10000 )*10000;
  if (code == 806) return (mh-80600000000)  - ( ((Int_t)(mh-80600000000))/10000 )*10000;
  return -1;
}

int code(double mh){
  return (int)(mh/100000000);
}

///////////////////////////
void makeTemplatesInterpolationComparisonVsCoupling(string file_coupling_0p25,
						    string file_coupling_1p0,
						    string file_coupling_0p1,
						    string outputDIR,
						    Category category,
						    int    mDM = 1){

  system(("mkdir -p "+outputDIR).c_str());
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  
  TFile* inputFile_1 = TFile::Open(file_coupling_0p25.c_str(),"READ");
  TFile* inputFile_2 = TFile::Open(file_coupling_1p0.c_str(),"READ");
  TFile* inputFile_3 = TFile::Open(file_coupling_0p1.c_str(),"READ");

  string directory;
  if(category == Category::monojet)
    directory = "category_monojet";
  else if(category == Category::monoV)
    directory = "category_monov";

  TDirectory* dir_file_1 = (TDirectory*) inputFile_1->Get(directory.c_str());
  
  TIter next(dir_file_1->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    TString name = Form("%s",key->GetName());
    name.ReplaceAll("signal_signal_","");
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;   
    
    double masspoint = stod(string(name));
    int c =  code(masspoint);
    int medMass = mmed(masspoint,c);
    int dmMass  = mdm(masspoint,c);

    if(dmMass != mDM) continue;
    if(medMass >= 3000) continue;

    TH1* histo1 = (TH1*) key->ReadObj();    
    TH1* histo2 = (TH1*) inputFile_2->FindObjectAny(key->GetName());
    if(histo2 == 0 or histo2 == NULL){
      cerr<<"histo file 2 not found with name "<<key->GetName()<<" skip "<<endl;
      continue;
    }

    TH1* histo3 = (TH1*) inputFile_3->FindObjectAny(key->GetName());
    if(histo3 == 0 or histo3 == NULL){
      cerr<<"histo file 3 not found with name "<<key->GetName()<<" skip "<<endl;
      continue;
    }    

    checkNegativeBin(histo1);
    checkNegativeBin(histo2);
    checkNegativeBin(histo3);

    TCanvas* canvas = new TCanvas("canvas","",600,700);
    canvas->cd();
    canvas->SetBottomMargin(0.3);
    TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
    pad2->SetTopMargin(0.7);
    pad2->SetRightMargin(0.05);
    pad2->SetFillColor(0);
    pad2->SetFillStyle(0);

    canvas->cd();    
    histo1->Scale(1,"width");
    histo2->Scale(1,"width");
    histo3->Scale(1,"width");
    histo1->GetYaxis()->SetTitle("Events / GeV");
    histo1->GetYaxis()->SetTitleOffset(1.35);
    histo1->SetMarkerStyle(20);
    histo1->SetMarkerSize(1);
    histo1->SetMarkerColor(kBlack);
    histo1->SetLineColor(kBlack);
    histo1->SetLineWidth(2);
    histo1->GetXaxis()->SetLabelSize(0);
    histo1->GetXaxis()->SetTitleSize(0);
    histo1->Draw("hist");

    histo2->SetLineColor(kRed);
    histo2->SetLineWidth(2);
    histo2->Draw("hist same");

    histo3->SetLineColor(kBlue);
    histo3->SetLineWidth(2);
    histo3->Draw("hist same");

    histo1->GetYaxis()->SetRangeUser(min(histo1->GetMinimum(),min(histo2->GetMinimum(),histo3->GetMinimum()))*0.1,
				     max(histo1->GetMaximum(),min(histo2->GetMaximum(),histo3->GetMaximum()))*100);

    TLegend leg (0.6,0.60,0.92,0.92);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);
    leg.SetFillColor(0);
    leg.AddEntry((TObject*)(0),("m_{MED}="+to_string(medMass)+" GeV "+"m_{DM}="+to_string(dmMass)+" GeV").c_str(),"");
    leg.AddEntry(histo1,"g_{q} = 0.25, g_{DM} = 1","L");
    leg.AddEntry(histo2,"g_{q} = 1.0, g_{DM} = 1","L");
    leg.AddEntry(histo3,"g_{q} = 0.1, g_{DM} = 1","L");
    leg.Draw("same");
    
    CMS_lumi(canvas,"35.9");

    canvas->RedrawAxis("sameaxis");
    canvas->SetLogy();
    canvas->cd();

    pad2->Draw();
    pad2->cd();

    TH1* ratio1 = (TH1*) histo2->Clone("ratio1");
    TH1* ratio2 = (TH1*) histo3->Clone("ratio2");
    ratio1->Divide(histo1);
    ratio2->Divide(histo1);

    ratio1->SetLineColor(kRed);
    ratio2->SetLineColor(kBlue);

    ratio1->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
    ratio1->GetYaxis()->SetTitle("FullSIM/Interp.");
    ratio1->GetYaxis()->CenterTitle();
    ratio1->GetYaxis()->SetRangeUser(0.0,2);
    ratio1->GetYaxis()->SetNdivisions(504);
    ratio1->GetYaxis()->SetTitleOffset(1.5);
    ratio1->GetYaxis()->SetLabelSize(0.04);
    ratio1->GetYaxis()->SetTitleSize(0.04);
    ratio1->GetXaxis()->SetLabelSize(0.04);
    ratio1->GetXaxis()->SetTitleSize(0.05);
    ratio1->GetXaxis()->SetTitleOffset(1.1);
    ratio1->Draw("hist");
    ratio2->Draw("hist same");

    if(category == Category::monojet){
      canvas->SaveAs((outputDIR+"/comparison_monojet_"+string(name)+".png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/comparison_monojet_"+string(name)+".pdf").c_str(),"pdf");
    }
    else if(category == Category::monoV){
      canvas->SaveAs((outputDIR+"/comparison_monoV_"+string(name)+".png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/comparison_monoV_"+string(name)+".pdf").c_str(),"pdf");
    }

    if(canvas) delete canvas;
  }      
}
