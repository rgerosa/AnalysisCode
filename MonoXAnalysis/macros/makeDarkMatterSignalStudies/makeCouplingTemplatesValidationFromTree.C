#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void checkNegativeBin(TH1* histo){
  for(int iBin = 0; iBin < histo->GetNbinsX()+1; iBin++){
    if(histo->GetBinContent(iBin+1) < 0)
      histo->SetBinContent(iBin+1,0);
  }
}


static string luminosity = "35.9";

int mmed(double mh, int code){

  if (code == 800) return ((int)(mh-80000000000))/10000;
  if (code == 801) return ((int)(mh-80100000000))/10000;
  if (code == 805) return ((int)(mh-80500000000))/10000;
  if (code == 806) return ((int)(mh-80600000000))/10000;
  return -1;
}

//                                                                                                                                                                                                      
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



void makeCouplingTemplatesValidationFromTree(string inputTreeInterpolationName,  // tree used for the interpoaltion
					     string inputTemplateName,           // template file
					     Category category,                  // category of selections
					     string outputDirectory, 
					     string interpolationPoint,          //string with the mass point to extract the templated like 80010000001
					     string gq,
					     string gDM){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDirectory).c_str());
  
  /// find interpolated histogram
  TFile* inputFileTemplate   = TFile::Open(inputTemplateName.c_str());
  TH1F*  histo_interpolation  = (TH1F*) inputFileTemplate->FindObjectAny(("signal_signal_"+interpolationPoint).c_str());

  int c       = code(stod(interpolationPoint));
  int medMass = mmed(stod(interpolationPoint),c);
  int dmMass  = mdm(stod(interpolationPoint),c);
  
  // find file
  TFile* inputTreeInterpFile = TFile::Open(inputTreeInterpolationName.c_str());
  TTree* tree = (TTree*) inputTreeInterpFile->Get("tree");
  TH1F* histo_fullSIM = (TH1F*) histo_interpolation->Clone("histo_fullSIM");

  int pos = 0;
  vector<float> *gq_vec = NULL;
  vector<float> *gDM_vec = NULL;
  tree->SetBranchAddress("gSM",&gq_vec);
  tree->SetBranchAddress("gDM",&gDM_vec);
  tree->GetEntry(0);
  for(int size = 0; size < gq_vec->size(); size++){
    if((gq_vec->at(size) >= stod(gq)-0.001 and gq_vec->at(size) <= stod(gq)+0.001) and (gDM_vec->at(size) >= stod(gDM) - 0.001 and gDM_vec->at(size) <= stod(gDM)+0.001)){
      pos = size;
      break;
    }
  }

  histo_fullSIM->Reset();
  if(category == Category::monojet)
    tree->Draw("pfMetPt >> histo_fullSIM",("("+luminosity+"*weightPU*weightTurnOn*genWeight*xsec*(1./sumwgt)*couplingwgt["+to_string(pos)+"]/genWeight)*(id == 1 && genMediatorMass == "+to_string(medMass)+" && genX1Mass == "+to_string(dmMass)+")").c_str(),"goff");
  else if(category == Category::monoV)
    tree->Draw("pfMetPt >> histo_fullSIM",("("+luminosity+"*weightPU*weightTurnOn*genWeight*xsec*(1./sumwgt)*couplingwgt["+to_string(pos)+"]/genWeight)*(id == 2 && genMediatorMass == "+to_string(medMass)+" && genX1Mass == "+to_string(dmMass)+")").c_str(),"goff");
  
  // fix the overflow
  histo_fullSIM->SetBinContent(histo_fullSIM->GetNbinsX(),histo_fullSIM->GetBinContent(histo_fullSIM->GetNbinsX())+histo_fullSIM->GetBinContent(histo_fullSIM->GetNbinsX()+1));
  histo_fullSIM->SetBinError(histo_fullSIM->GetNbinsX(),sqrt(TMath::Power(histo_fullSIM->GetBinError(histo_fullSIM->GetNbinsX()),2)+
							     TMath::Power(histo_fullSIM->GetBinError(histo_fullSIM->GetNbinsX()+1),2)));
 


  checkNegativeBin(histo_fullSIM);
  checkNegativeBin(histo_interpolation);
  

  TCanvas* canvas = new TCanvas("canvas","",600,700);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  TPad* pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.05);
  pad2->SetFillColor(0);
  pad2->SetFillStyle(0);

  canvas->cd();
  histo_interpolation->Scale(1,"width");
  histo_fullSIM->Scale(1,"width");
  histo_interpolation->GetYaxis()->SetTitle("Events / GeV");
  histo_interpolation->GetYaxis()->SetTitleOffset(1.35);
  histo_interpolation->SetMarkerStyle(20);
  histo_interpolation->SetMarkerSize(1);
  histo_interpolation->SetMarkerColor(kBlack);
  histo_interpolation->SetLineColor(kBlack);
  histo_interpolation->GetXaxis()->SetLabelSize(0);
  histo_interpolation->GetXaxis()->SetTitleSize(0);
  histo_interpolation->Draw("EP");
 
  histo_fullSIM->SetLineColor(kRed);
  histo_fullSIM->SetLineWidth(2);

  TH1* histo_fullSIM_band = (TH1*) histo_fullSIM->Clone("histo_fullSIM_band");
  histo_fullSIM_band->SetFillColor(kGray);
  histo_fullSIM_band->Draw("E2 same");
  histo_fullSIM->Draw("hist same");
  histo_interpolation->Draw("EPsame");

  histo_interpolation->GetYaxis()->SetRangeUser(min(histo_interpolation->GetMinimum(),histo_fullSIM->GetMinimum())*0.01,
						max(histo_interpolation->GetMaximum(),histo_fullSIM->GetMaximum())*100);

  if(min(histo_interpolation->GetMinimum(),histo_fullSIM->GetMinimum()) == 0)
    histo_interpolation->GetYaxis()->SetRangeUser(0.0001,max(histo_interpolation->GetMaximum(),histo_fullSIM->GetMaximum())*100);

  TLegend leg (0.6,0.60,0.92,0.92);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.SetFillColor(0);
  leg.AddEntry((TObject*)(0),("m_{MED}="+to_string(medMass)+" GeV "+"m_{DM}="+to_string(dmMass)+" GeV").c_str(),"");
  leg.AddEntry((TObject*)(0),("gSM = "+gq+" gDM = "+gDM).c_str(),"");
  leg.AddEntry(histo_interpolation,"Interpolated template","EP");
  leg.AddEntry(histo_fullSIM_band,"Full SIM template","FL");
  leg.Draw("same");

  CMS_lumi(canvas,luminosity);

  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();
  canvas->cd();

  pad2->Draw();
  pad2->cd();

  TH1* ratio = (TH1*) histo_fullSIM->Clone("ratio");
  ratio->Divide(histo_interpolation);
  ratio->SetMarkerColor(kBlack);
  ratio->SetLineColor(kBlack);
  ratio->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
  ratio->GetYaxis()->SetTitle("FullSIM/Interp.");
  ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetRangeUser(0.3,1.7);
  ratio->GetYaxis()->SetNdivisions(504);
  ratio->GetYaxis()->SetTitleOffset(1.5);
  ratio->GetYaxis()->SetLabelSize(0.04);
  ratio->GetYaxis()->SetTitleSize(0.04);
  ratio->GetXaxis()->SetLabelSize(0.04);
  ratio->GetXaxis()->SetTitleSize(0.05);
  ratio->GetXaxis()->SetTitleOffset(1.1);
  ratio->Draw("EP");

  if(category == Category::monojet){
    canvas->SaveAs((outputDirectory+"/comparison_monojet_"+interpolationPoint+".png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/comparison_monojet_"+interpolationPoint+".pdf").c_str(),"pdf");
  }
  else if(category == Category::monoV){
    canvas->SaveAs((outputDirectory+"/comparison_monoV_"+interpolationPoint+".png").c_str(),"png");
    canvas->SaveAs((outputDirectory+"/comparison_monoV_"+interpolationPoint+".pdf").c_str(),"pdf");
  }
}
