#include "../CMS_lumi.h"

//////////////
vector<double> mediatorPtBinning   = {150,180,200,220,240,260,280,300,330,360,390,420,460,500,550,600,650,700,800,900,1000,1250};
vector<double> mediatorEtaBinning  = {0,1,2.5,5};
vector<double> skipMediator = {10.};

class massEfficiencyPoint {
public:
  massEfficiencyPoint();  
  massEfficiencyPoint(double mmed, double mdm, TH1* histo):
    mmed(mmed),
    mdm(mdm),
    histo(histo){}

  massEfficiencyPoint(double mmed, double mdm, TGraphAsymmErrors* eff):
    mmed(mmed),
    mdm(mdm),
    eff(eff){}

  bool operator == (const massEfficiencyPoint & b) {
    if(mmed == b.mmed and mdm == b.mdm) return true;
    else return false;
  }

  bool operator < (const massEfficiencyPoint & b) {
    if(mmed < b.mmed) return true;
    else if(mmed == b.mmed and mdm < b.mdm) return true;
    else return false;
  }
  
  double mmed;
  double mdm;
  TH1*   histo;
  TGraphAsymmErrors* eff;
  TF1* funz;
};

  
// use the small trees produced by makeTreesForInterpolation/makeSmallGenTree.C and used for gen-interpolation to make a spline
void makeSignalInterpolation (string inputFileName, string outputDIR, string category){
  
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());
  
  if(category != "monojet" and category != "monoV"){
    cerr<<"Problem with category definition --> should be monojet or monoV --> return "<<endl;
    return;
  }
    
  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TTree* inputTree  = (TTree*) inputFile->Get("tree");
  
  // for each mass point we need to store the correct efficiency after storing the related histograms for numerator and denominator
  TTreeReader reader(inputTree);
  TTreeReaderValue<int> id (reader,"id");
  TTreeReaderValue<double> genMediatorPt (reader,"genMediatorPt");
  TTreeReaderValue<double> genMediatorEta (reader,"genMediatorEta");
  TTreeReaderValue<double> genMediatorMass (reader,"genMediatorMass");
  TTreeReaderValue<double> genX1Pt (reader,"genX1Pt");
  TTreeReaderValue<double> genX1Eta (reader,"genX1Eta");
  TTreeReaderValue<double> genX1Mass (reader,"genX1Mass");
  TTreeReaderValue<double> genX2Pt (reader,"genX2Pt");
  TTreeReaderValue<double> genX2Eta (reader,"genX2Eta");
  TTreeReaderValue<double> genX2Mass (reader,"genX2Mass");
  TTreeReaderValue<double> pfMetPt (reader,"pfMetPt");
  TTreeReaderValue<double> weightPU (reader,"weightPU");
  TTreeReaderValue<double> weightTurnOn (reader,"weightTurnOn");
  TTreeReaderValue<double> genWeight (reader,"genWeight");

  // make map
  vector<massEfficiencyPoint> efficiencyNumerator;
  vector<massEfficiencyPoint> efficiencyDenominator;
  vector<massEfficiencyPoint> efficiency;

  TH1F* histtemp = NULL;
  while(reader.Next()){
    const massEfficiencyPoint temp(*genMediatorMass,*genX1Mass,histtemp);
    if(std::find(efficiencyNumerator.begin(),efficiencyNumerator.end(),temp) == efficiencyNumerator.end()){
      // use to skip low mediators
      if(std::find(skipMediator.begin(),skipMediator.end(),*genMediatorMass) != skipMediator.end()) continue;

      cout<<"Create Efficiency histogram for : mediator "<<temp.mmed<<" dm mass "<<temp.mdm<<endl;
      efficiencyNumerator.push_back(massEfficiencyPoint(*genMediatorMass,*genX1Mass,new TH1F(Form("numerator_%f_%f",*genMediatorMass,*genX1Mass),"",mediatorPtBinning.size()-1,&mediatorPtBinning[0])));
      efficiencyDenominator.push_back(massEfficiencyPoint(*genMediatorMass,*genX1Mass,new TH1F(Form("denominator_%f_%f",*genMediatorMass,*genX1Mass),"",mediatorPtBinning.size()-1,&mediatorPtBinning[0])));
      efficiency.push_back(massEfficiencyPoint(*genMediatorMass,*genX1Mass,new TGraphAsymmErrors()));
      efficiency.back().eff->SetName(Form("efficiency_%f_%f",efficiency.back().mmed,efficiency.back().mdm));
      efficiencyNumerator.back().histo->Sumw2();
      efficiencyDenominator.back().histo->Sumw2();
    }
  }

  std::sort(efficiencyNumerator.begin(),efficiencyNumerator.end());
  std::sort(efficiencyDenominator.begin(),efficiencyDenominator.end());
  std::sort(efficiency.begin(),efficiency.end());

  //////////////////////
  reader.SetEntry(0);
  while(reader.Next()){
    const massEfficiencyPoint temp(*genMediatorMass,*genX1Mass,histtemp);
    // use to skip low mediators
    if(std::find(skipMediator.begin(),skipMediator.end(),*genMediatorMass) != skipMediator.end()) continue;

    int pos = std::find(efficiencyNumerator.begin(),efficiencyNumerator.end(),temp)-efficiencyNumerator.begin();
    if(pos < efficiencyNumerator.size()){
      efficiencyDenominator.at(pos).histo->Fill(*genMediatorPt,(*weightPU)*(*genWeight));      
      if(category == "monojet" and *id == 1){
	efficiencyNumerator.at(pos).histo->Fill(*genMediatorPt,(*weightPU)*(*genWeight)*(*weightTurnOn));
      }
      if(category == "monoV" and *id == 2){
	efficiencyNumerator.at(pos).histo->Fill(*genMediatorPt,(*weightPU)*(*genWeight)*(*weightTurnOn));
      }
    } 
    else{
      cerr<<"Not able to find this mass point: mediator mass "<<*genMediatorMass<<" "<<" dmMass "<<*genX1Mass<<endl;
    }
  }
  
  /// make the ratio
  for(size_t ipos = 0; ipos < efficiency.size(); ipos++){
    cout<<"Efficiency map: med "<<efficiencyNumerator.at(ipos).mmed<<" dm "<<efficiencyNumerator.at(ipos).mdm<<" Integral denominator "<<efficiencyDenominator.at(ipos).histo->Integral()<<" Integral numerator "<<efficiencyNumerator.at(ipos).histo->Integral()<<endl;
    if(efficiencyNumerator.at(ipos).histo->Integral() == 0 and efficiencyDenominator.at(ipos).histo->Integral() == 0){ // remove the point
      efficiencyNumerator.erase(efficiencyNumerator.begin() + ipos);
      efficiencyDenominator.erase(efficiencyDenominator.begin() + ipos);
      efficiency.erase(efficiency.begin() + ipos);
      cout<<"erased "<<efficiency.size()<<endl;
      continue;
    }      
    efficiency.at(ipos).eff->BayesDivide(efficiencyNumerator.at(ipos).histo,efficiencyDenominator.at(ipos).histo);    
  }
  
  //Plot efficiency
  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();
  TLegend leg (0.7,0.25,0.9,0.5);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  vector<TGraph2DErrors*> efficiencyParam;
  int nPoint = 0;
  for(auto element : efficiency){
    element.eff->GetYaxis()->SetRangeUser(0,1);
    element.eff->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
    element.eff->GetXaxis()->SetTitleOffset(1.1);
    element.eff->GetYaxis()->SetTitle("Efficiency");
    element.eff->GetYaxis()->SetTitleOffset(1.1);
    element.eff->SetMarkerColor(kBlack);
    element.eff->SetMarkerSize(1);
    element.eff->SetMarkerStyle(20);
    element.histo = element.eff->GetHistogram();
    element.funz = new TF1(Form("funz_%f_%f",element.mmed,element.mmed),"[0]*(1+TMath::Erf((x-[1])/[2]))*TMath::Exp(-[3]*x)",element.histo->GetBinLowEdge(1),element.histo->GetBinLowEdge(element.histo->GetNbinsX()+1));
    element.funz->SetParameter(0,0.5);
    element.funz->SetParameter(1,200);
    element.funz->SetParameter(2,100);
    element.funz->SetParameter(3,-0.001);
    element.funz->SetLineColor(kBlue);
    element.funz->SetLineWidth(2);
    element.eff->Fit(element.funz,"QMR");
    element.histo->Draw("EP");
    element.eff->Draw("EPsame");

    if(efficiencyParam.size() == 0){
      efficiencyParam.push_back(new TGraph2DErrors());
      efficiencyParam.back()->SetName("efficiencyParam_amp");
      efficiencyParam.push_back(new TGraph2DErrors());
      efficiencyParam.back()->SetName("efficiencyParam_offset");
      efficiencyParam.push_back(new TGraph2DErrors());
      efficiencyParam.back()->SetName("efficiencyParam_width");
      efficiencyParam.push_back(new TGraph2DErrors());
      efficiencyParam.back()->SetName("efficiencyParam_exp");
    }

    // set point
    efficiencyParam.at(0)->SetPoint(nPoint,element.mmed,element.mdm,element.funz->GetParameter(0));
    efficiencyParam.at(1)->SetPoint(nPoint,element.mmed,element.mdm,element.funz->GetParameter(1));
    efficiencyParam.at(2)->SetPoint(nPoint,element.mmed,element.mdm,element.funz->GetParameter(2));
    efficiencyParam.at(3)->SetPoint(nPoint,element.mmed,element.mdm,element.funz->GetParameter(3));
    // set Error
    efficiencyParam.at(0)->SetPointError(nPoint,0,0,element.funz->GetParError(0));
    efficiencyParam.at(1)->SetPointError(nPoint,0,0,element.funz->GetParError(1));
    efficiencyParam.at(2)->SetPointError(nPoint,0,0,element.funz->GetParError(2));
    efficiencyParam.at(3)->SetPointError(nPoint,0,0,element.funz->GetParError(3));    
    nPoint++;
    cout<<"##############"<<endl;
    cout<<"Parameter [0] From the fit: mmed "<<element.mmed<<" mdm "<<element.mdm<<" param "<<element.funz->GetParameter(0)<<" error "<<element.funz->GetParError(0)<<endl;
    cout<<"Parameter [1] From the fit: mmed "<<element.mmed<<" mdm "<<element.mdm<<" param "<<element.funz->GetParameter(1)<<" error "<<element.funz->GetParError(1)<<endl;
    cout<<"Parameter [2] From the fit: mmed "<<element.mmed<<" mdm "<<element.mdm<<" param "<<element.funz->GetParameter(2)<<" error "<<element.funz->GetParError(2)<<endl;
    cout<<"Parameter [3] From the fit: mmed "<<element.mmed<<" mdm "<<element.mdm<<" param "<<element.funz->GetParameter(3)<<" error "<<element.funz->GetParError(3)<<endl;
    cout<<"##############"<<endl;

    CMS_lumi(canvas,"12.9");
    leg.Clear();
    leg.AddEntry((TObject*)0,Form("m_{MED} = %d, m_{DM} = %d",int(element.mmed),int(element.mdm)),"");
    canvas->SaveAs(Form("%s/efficiency_med_%d_dm_%d.png",outputDIR.c_str(),int(element.mmed),int(element.mdm)),"png");
    canvas->SaveAs(Form("%s/efficiency_med_%d_dm_%d.pdf",outputDIR.c_str(),int(element.mmed),int(element.mdm)),"pdf");
  }

  TCanvas *canvas2D = new TCanvas("canvas2D","",0,0,650,600);
  TH2* histoTemp = efficiencyParam.at(0)->GetHistogram();  
  histoTemp->GetXaxis()->SetTitle("m_{MED} [GeV]");
  histoTemp->GetYaxis()->SetTitle("m_{DM} [GeV]");
  histoTemp->GetZaxis()->SetTitle("Amplitude");
  histoTemp->GetXaxis()->SetTitleOffset(1.5);
  histoTemp->GetYaxis()->SetTitleOffset(1.5);
  histoTemp->GetZaxis()->SetTitleOffset(1.3);
  histoTemp->SetMarkerStyle(20);
  histoTemp->SetMarkerSize(1);
  histoTemp->SetMarkerColor(kBlack);
  histoTemp->GetZaxis()->SetRangeUser(0.8*TMath::MinElement(efficiencyParam.at(0)->GetN(),efficiencyParam.at(0)->GetZ()),1.2*TMath::MaxElement(efficiencyParam.at(0)->GetN(),efficiencyParam.at(0)->GetZ()));
  histoTemp->Reset();
  histoTemp->Draw("lego");
  efficiencyParam.at(0)->Draw("P0err same");  
  canvas2D->SaveAs((outputDIR+"/fitParameter_amplitude_2D.png").c_str(),"png");
  canvas2D->SaveAs((outputDIR+"/fitParameter_amplitude_2D.pdf").c_str(),"pdf");

  histoTemp = efficiencyParam.at(1)->GetHistogram();  
  histoTemp->GetXaxis()->SetTitle("m_{MED} [GeV]");
  histoTemp->GetYaxis()->SetTitle("m_{DM} [GeV]");
  histoTemp->GetZaxis()->SetTitle("Offset");
  histoTemp->GetXaxis()->SetTitleOffset(1.5);
  histoTemp->GetYaxis()->SetTitleOffset(1.5);
  histoTemp->GetZaxis()->SetTitleOffset(1.3);
  histoTemp->SetMarkerStyle(20);
  histoTemp->SetMarkerSize(1);
  histoTemp->SetMarkerColor(kBlack);
  histoTemp->GetZaxis()->SetRangeUser(0.8*TMath::MinElement(efficiencyParam.at(1)->GetN(),efficiencyParam.at(1)->GetZ()),1.2*TMath::MaxElement(efficiencyParam.at(1)->GetN(),efficiencyParam.at(1)->GetZ()));
  histoTemp->Reset();
  histoTemp->Draw("lego");
  efficiencyParam.at(1)->Draw("P0err same");  
  canvas2D->SaveAs((outputDIR+"/fitParameter_offset_2D.png").c_str(),"png");
  canvas2D->SaveAs((outputDIR+"/fitParameter_offset_2D.pdf").c_str(),"pdf");

  histoTemp = efficiencyParam.at(2)->GetHistogram();  
  histoTemp->GetXaxis()->SetTitle("m_{MED} [GeV]");
  histoTemp->GetYaxis()->SetTitle("m_{DM} [GeV]");
  histoTemp->GetZaxis()->SetTitle("Width");
  histoTemp->GetXaxis()->SetTitleOffset(1.5);
  histoTemp->GetYaxis()->SetTitleOffset(1.5);
  histoTemp->GetZaxis()->SetTitleOffset(1.3);
  histoTemp->SetMarkerStyle(20);
  histoTemp->SetMarkerSize(1);
  histoTemp->SetMarkerColor(kBlack);
  histoTemp->GetZaxis()->SetRangeUser(0.8*TMath::MinElement(efficiencyParam.at(2)->GetN(),efficiencyParam.at(2)->GetZ()),1.2*TMath::MaxElement(efficiencyParam.at(2)->GetN(),efficiencyParam.at(2)->GetZ()));
  histoTemp->Reset();
  histoTemp->Draw("lego");
  efficiencyParam.at(2)->Draw("P0err same");  
  canvas2D->SaveAs((outputDIR+"/fitParameter_width_2D.png").c_str(),"png");
  canvas2D->SaveAs((outputDIR+"/fitParameter_width_2D.pdf").c_str(),"pdf");

  histoTemp = efficiencyParam.at(3)->GetHistogram();  
  histoTemp->GetXaxis()->SetTitle("m_{MED} [GeV]");
  histoTemp->GetYaxis()->SetTitle("m_{DM} [GeV]");
  histoTemp->GetZaxis()->SetTitle("Exp");
  histoTemp->GetXaxis()->SetTitleOffset(1.5);
  histoTemp->GetYaxis()->SetTitleOffset(1.5);
  histoTemp->GetZaxis()->SetTitleOffset(1.3);
  histoTemp->SetMarkerStyle(20);
  histoTemp->SetMarkerSize(1);
  histoTemp->SetMarkerColor(kBlack);
  histoTemp->GetZaxis()->SetRangeUser(0.8*TMath::MinElement(efficiencyParam.at(3)->GetN(),efficiencyParam.at(3)->GetZ()),1.2*TMath::MaxElement(efficiencyParam.at(3)->GetN(),efficiencyParam.at(3)->GetZ()));
  histoTemp->Reset();
  histoTemp->Draw("lego");
  efficiencyParam.at(3)->Draw("P0err same");  
  canvas2D->SaveAs((outputDIR+"/fitParameter_exp_2D.png").c_str(),"png");
  canvas2D->SaveAs((outputDIR+"/fitParameter_exp_2D.pdf").c_str(),"pdf");
  
}
