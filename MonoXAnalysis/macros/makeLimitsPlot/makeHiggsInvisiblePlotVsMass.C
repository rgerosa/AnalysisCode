#include "../CMS_lumi.h"

void makeHiggsInvisiblePlotVsMass(string inputFileName, string outputDIR, bool makeXSecLimit = true, bool isBlind = true){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  // for VBF production in case of XSEC limit vs Mass
  map<float,float> crossSectionMap;
  crossSectionMap[110] = 4.434E+00;
  crossSectionMap[125] = 3.925E+00;
  crossSectionMap[150] = 3.239E+00;
  crossSectionMap[200] = 2.282E+00;
  crossSectionMap[300] = 1.256E+00;
  crossSectionMap[400] = 7.580E-01;
  crossSectionMap[500] = 4.872E-01;
  crossSectionMap[600] = 3.274E-01;
  crossSectionMap[800] = 1.622E-01;
  crossSectionMap[1000] = 8.732E-02;
  
  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TTree* tree = (TTree*) inputFile->Get("limit");
  
  TTreeReader reader (tree);
  TTreeReaderValue<double> mh (reader,"mh");
  TTreeReaderValue<float> quantileExpected (reader,"quantileExpected");
  TTreeReaderValue<double> limit (reader,"limit");

  TGraph* graph_exp = new TGraph();
  TGraph* graph_obs = new TGraph();
  TGraph* graph_exp_1s_up = new TGraph();
  TGraph* graph_exp_1s_dw = new TGraph();
  TGraph* graph_exp_2s_up = new TGraph();
  TGraph* graph_exp_2s_dw = new TGraph();
  int ipoint_obs = 0;
  int ipoint_exp = 0;
  int ipoint_exp_1s_up = 0;
  int ipoint_exp_1s_dw = 0;
  int ipoint_exp_2s_up = 0;
  int ipoint_exp_2s_dw = 0;

  while(reader.Next()){

    if(*quantileExpected == -1 and not isBlind){
      if(not makeXSecLimit)
	graph_obs->SetPoint(ipoint_obs,*mh,*limit);
      else
	graph_obs->SetPoint(ipoint_obs,*mh,*limit*crossSectionMap[*mh]);

      ipoint_obs++;
    }

    else if(*quantileExpected == 0.5){
      if(not makeXSecLimit){
	graph_exp->SetPoint(ipoint_exp,*mh,*limit);
      }
      else
	graph_exp->SetPoint(ipoint_exp,*mh,*limit*crossSectionMap[*mh]);

      if(isBlind){
	if(not makeXSecLimit)
	  graph_obs->SetPoint(ipoint_exp,*mh,*limit);
	else
	  graph_obs->SetPoint(ipoint_obs,*mh,*limit*crossSectionMap[*mh]);
      }      
      ipoint_exp++;
    }
    else if(*quantileExpected >= 0.83 and *quantileExpected <= 0.85){
      if(not makeXSecLimit)
	graph_exp_1s_up->SetPoint(ipoint_exp_1s_up,*mh,*limit);
      else
	graph_exp_1s_up->SetPoint(ipoint_exp_1s_up,*mh,*limit*crossSectionMap[*mh]);
	
      ipoint_exp_1s_up++;
    }
    else if(*quantileExpected >= 0.9 and *quantileExpected <= 1){
      if(not makeXSecLimit)
	graph_exp_2s_up->SetPoint(ipoint_exp_2s_up,*mh,*limit);
      else
	graph_exp_2s_up->SetPoint(ipoint_exp_2s_up,*mh,*limit*crossSectionMap[*mh]);
      ipoint_exp_2s_up++;
    }
    else if(*quantileExpected >= 0.14 and *quantileExpected <= 0.17){
      if(not makeXSecLimit)
	graph_exp_1s_dw->SetPoint(ipoint_exp_1s_dw,*mh,*limit);
      else
	graph_exp_1s_dw->SetPoint(ipoint_exp_1s_dw,*mh,*limit*crossSectionMap[*mh]);
      ipoint_exp_1s_dw++;
    }
    else if(*quantileExpected >= 0 and *quantileExpected <= 0.09){
      if(not makeXSecLimit)
	graph_exp_2s_dw->SetPoint(ipoint_exp_2s_dw,*mh,*limit);
      else
	graph_exp_2s_dw->SetPoint(ipoint_exp_2s_dw,*mh,*limit*crossSectionMap[*mh]);
      ipoint_exp_2s_dw++;
    }    
  }

  // sort entries as a function of the mass
  graph_obs->Sort();
  graph_exp->Sort();
  graph_exp_1s_up->Sort();
  graph_exp_2s_up->Sort();
  graph_exp_1s_dw->Sort();
  graph_exp_2s_dw->Sort();

  // make theory line
  TGraph* theory_xsec = new TGraph();
  int ipoint = 0;
  for(auto item : crossSectionMap){
    theory_xsec->SetPoint(ipoint,item.first,item.second);
    ipoint++;
  }
  theory_xsec->Sort();

  //Make the brazilian band
  TGraphAsymmErrors* graph_1s = new TGraphAsymmErrors();
  TGraphAsymmErrors* graph_2s = new TGraphAsymmErrors();

  for(int ipoint = 0; ipoint < graph_exp->GetN(); ipoint++){
    double x,y;
    graph_exp->GetPoint(ipoint,x,y);
    graph_1s->SetPoint(ipoint,x,y);
    graph_2s->SetPoint(ipoint,x,y);
    double x_up,y_up;
    double x_dw,y_dw;
    graph_exp_1s_up->GetPoint(ipoint,x_up,y_up);
    graph_exp_1s_dw->GetPoint(ipoint,x_dw,y_dw);
    double y_temp, x_temp;
    graph_exp_1s_up->GetPoint(ipoint+1,x_temp,y_temp);    
    graph_1s->SetPointError(ipoint,0.,x_temp-x,fabs(y-y_dw),fabs(y-y_up));
    ///
    graph_exp_2s_up->GetPoint(ipoint,x_up,y_up);
    graph_exp_2s_dw->GetPoint(ipoint,x_dw,y_dw);
    graph_2s->SetPointError(ipoint,0.,x_temp-x,fabs(y-y_dw),fabs(y-y_up));
  }

  // make the plot
  TCanvas* canvas = new TCanvas("canvas","",600,600);
  canvas->cd();
  TH1* frame = canvas->DrawFrame(TMath::MinElement(graph_exp->GetN(),graph_exp->GetX())-50,TMath::MinElement(graph_exp->GetN(),graph_exp->GetY())*0.1,
				 TMath::MaxElement(graph_exp->GetN(),graph_exp->GetX())+50,TMath::MaxElement(graph_exp->GetN(),graph_exp->GetY())*100);

  if(not makeXSecLimit)
    frame->GetYaxis()->SetRangeUser(TMath::MinElement(graph_exp->GetN(),graph_exp->GetY())*0.5,TMath::MaxElement(graph_exp->GetN(),graph_exp->GetY())*1.5);

  frame->GetXaxis()->SetTitle("Higgs boson mass [GeV]");
  if(makeXSecLimit)
    frame->GetYaxis()->SetTitle("#sigma_{qq #rightarrow H} #times BR(H#rightarrow inv) [pb]");
  else
    frame->GetYaxis()->SetTitle("95% C.L. upper limit on #mu = #sigma x BR/#sigma_{SM}");  
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetXaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleSize(0.040);
  frame->GetXaxis()->SetTitleSize(0.040);
  frame->GetYaxis()->SetTitleOffset(1.25);
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->Draw();
  
  graph_1s->SetFillColor(kGreen+1);
  graph_2s->SetFillColor(kOrange);
  graph_2s->Draw("3same");
  graph_1s->Draw("3same");

  TF1* line = NULL;
  if(not makeXSecLimit){
    line = new TF1("line","1",TMath::MinElement(graph_exp->GetN(),graph_exp->GetX()),TMath::MaxElement(graph_exp->GetN(),graph_exp->GetX()));
    line->SetLineColor(kRed);
    line->SetLineWidth(2);
    line->Draw("Lsame");
  }
  else{
    theory_xsec->SetLineColor(kRed);
    theory_xsec->SetLineWidth(2);
    theory_xsec->Draw("Lsame");
  }

  graph_exp->SetLineColor(kBlack);
  graph_exp->SetLineWidth(2);
  graph_exp->SetLineStyle(2);
  graph_exp->Draw("PE0Lsame");

  graph_obs->SetLineColor(kBlack);
  graph_obs->SetLineWidth(2);
  graph_obs->SetMarkerStyle(20);
  graph_obs->SetMarkerSize(1);
  graph_obs->Draw("PE0Lsame");


  TLegend* leg = NULL;
  if(makeXSecLimit)
    leg = new TLegend(0.5,0.6,0.9,0.9);
  else
    leg = new TLegend(0.2,0.5,0.6,0.80);

  leg->SetBorderSize(0.);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(graph_obs,"Observed 95% CL","PL");
  leg->AddEntry(graph_exp,"Median expected 95% CL","PL");
  leg->AddEntry(graph_1s,"68% expected","F");
  leg->AddEntry(graph_2s,"95% expected","F");
  if(makeXSecLimit)
    leg->AddEntry(theory_xsec,"#sigma_{qq #rightarrow H}","L");
  else
    leg->AddEntry(line,"Signal strength = 1","L");
  leg->Draw("same");

  CMS_lumi(canvas,"35.9",false,true,false);
  
  canvas->RedrawAxis("sameaxis");
  
  if(makeXSecLimit)
    canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/limit_higgs_vs_mass.png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/limit_higgs_vs_mass.pdf").c_str(),"pdf");

}
