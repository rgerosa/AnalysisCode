#include "../CMS_lumi.h"

// set of effective areas                                                                                                                                                                             
float getChargedHadronEA(float abseta){
  if(abseta < 1) return 0.0360;
  else if(abseta >  1 and abseta < 1.479) return 0.0377;
  else return 1;
}

float getNeutralHadronEA(float abseta){
  if(abseta < 1) return 0.0597;
  else if(abseta >  1 and abseta < 1.479) return 0.0807;
  else return 1;
}

float getPhotonEA(float abseta){
  if(abseta < 1) return 0.1210;
  else if(abseta >  1 and abseta < 1.479) return 0.1107;
  else return 1;
}


// photon id values --> class to handle it                                                                                                                                                           
class photonID {

 public:
 photonID(const float & HoE, const float & sigmaieie, const float & chadiso, const float & nhadiso0, const float & nhadiso1, const float & nhadiso2, const float & phiso0, const float & phiso1):
  HoE(HoE),
    sigmaieie(sigmaieie),
    chadiso(chadiso),
    nhadiso0(nhadiso0),
    nhadiso1(nhadiso1),
    nhadiso2(nhadiso2),
    phiso0(phiso0),
    phiso1(phiso1){
    };
  ~photonID(){};

  float HoE;
  float sigmaieie;
  float chadiso;
  float nhadiso0;
  float nhadiso1;
  float nhadiso2;
  float phiso0;
  float phiso1;
};

// class for bin fir                                                                                                                                                                                 
class fitPurity {
 public:
 fitPurity(const float & ptMin, const float & ptMax, TH1F* histo):
  ptMin(ptMin),
    ptMax(ptMax),
    phHisto(histo){
      ptMean = 0;
  }
    ~fitPurity(){};

  float ptMin;
  float ptMax;
  float ptMean;
  TH1F* phHisto;
};


////////////
void plotFitResult(TCanvas* canvas, TH1F* data, const RooAbsPdf & signal, const RooAbsPdf & background, const RooRealVar & x, 
		   const string & outputDIR, const int & ptMin, const int & ptMax, const string & postfix){

  //create histograms from pdfs
  TH1F* signal_hist = (TH1F*) signal.createHistogram(Form("signal%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax),x,RooFit::Binning(data->GetNbinsX(),data->GetXaxis()->GetBinLowEdge(1),data->GetXaxis()->GetBinLowEdge(data->GetNbinsX()+1)));
  signal_hist->Sumw2();
  TH1F* background_hist = (TH1F*) background.createHistogram(Form("background%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax),x,RooFit::Binning(data->GetNbinsX(),data->GetXaxis()->GetBinLowEdge(1),data->GetXaxis()->GetBinLowEdge(data->GetNbinsX()+1)));
  background_hist->Sumw2();

  TH1F* totalHist = (TH1F*) signal_hist->Clone(Form("totalHist%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax));
  totalHist->Add(background_hist);
  
  canvas->cd();
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);
  
  data->SetLineColor(kBlack);
  data->SetMarkerColor(kBlack);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(1);

  data->GetXaxis()->SetLabelSize(0);
  data->GetYaxis()->SetTitle("Events / GeV");
  data->GetYaxis()->SetRangeUser(data->GetMinimum()*0.1,data->GetMaximum()*100);
  data->Draw("EP");

  background_hist->SetLineColor(kGreen+1);
  background_hist->SetLineWidth(2);
  background_hist->Draw("hist same");

  signal_hist->SetLineColor(kBlue);
  signal_hist->SetLineWidth(2);
  signal_hist->Draw("hist same");

  totalHist->SetLineColor(kRed);
  totalHist->SetLineWidth(2);
  totalHist->Draw("hist same");

  TLegend leg (0.6,0.6,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(data,"Data","PE");
  leg.AddEntry(totalHist,"Total Fit","L");
  leg.AddEntry(signal_hist,"Signal Component","L");
  leg.AddEntry(background_hist,"Background Component","L");
  leg.Draw("same");
  canvas->RedrawAxis("sameaxis");

  pad2->Draw();
  pad2->cd();

  TH1* frame =  (TH1*) data->Clone(Form("frame_pt_%d_%d",ptMin,ptMax));
  frame->Reset("ICES");  
  frame->GetYaxis()->SetRangeUser(0.4,1.6);
  frame->GetYaxis()->SetNdivisions(5);
  frame->GetYaxis()->SetTitle("Data/Exp.");
  frame->GetXaxis()->SetTitle("Photon Isolation [GeV]");

  TH1F* ratio = (TH1F*) data->Clone(Form("ratio%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax));
  TH1F* denFit_noerr = (TH1F*) totalHist->Clone(Form("denFit_noerr%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax));
  TH1F* denFit = (TH1F*) totalHist->Clone(Form("denFit%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax));
  for(int iBin = 1; iBin < denFit->GetNbinsX()+1; iBin++)
    denFit_noerr->SetBinError(iBin,0.);

  ratio->Divide(denFit_noerr);
  denFit->Divide(denFit_noerr);
  denFit->SetFillColor(kGray);
  frame->Draw();
  ratio->Draw("EPsame");
  denFit->Draw("E2same");

  TF1* line = new TF1(Form("line%s_pt_%d_%d",postfix.c_str(),ptMin,ptMax),"1",ratio->GetXaxis()->GetBinLowEdge(1),ratio->GetXaxis()->GetBinLowEdge(ratio->GetNbinsX()+1));
  line->SetLineColor(kRed);
  line->SetLineWidth(2);
  line->Draw("Lsame");
  ratio->Draw("EPsame");
  pad2->RedrawAxis("sameaxis");
    

  CMS_lumi(canvas,"36.2");

  canvas->SetLogy();

  canvas->SaveAs((outputDIR+"/photonPurity"+postfix+"_pt_"+to_string(ptMin)+"_"+to_string(ptMax)+".png").c_str());
  canvas->SaveAs((outputDIR+"/photonPurity"+postfix+"_pt_"+to_string(ptMin)+"_"+to_string(ptMax)+".pdf").c_str());

}
