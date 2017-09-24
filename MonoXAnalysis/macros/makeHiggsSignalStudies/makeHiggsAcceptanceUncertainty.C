#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

static float luminosity = 35.9;

///////////////////
class theoryVariation {

public:
  theoryVariation(){
    histoCentral = NULL;
  };
  
  TH1F* histoCentral;
  vector<TH1F*> histoQCDScale;
  vector<TH1F*> histoPDF_1;
  vector<TH1F*> histoPDF_2;
  vector<TH1F*> histoPDF_3;

};

////////////////
TGraphAsymmErrors* makeGraphQCD(TH1F* histo, vector<TH1F*> variations){

  TGraphAsymmErrors* graph = new TGraphAsymmErrors();

  for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
    graph->SetPoint(iBin,histo->GetBinCenter(iBin+1),histo->GetBinContent(iBin+1));
    double maxValue = histo->GetBinContent(iBin+1);
    double minValue = histo->GetBinContent(iBin+1);
    for(auto hist: variations){
      if(hist->GetBinContent(iBin+1) < minValue)
	minValue = hist->GetBinContent(iBin+1);
      if(hist->GetBinContent(iBin+1) > maxValue)
	maxValue = hist->GetBinContent(iBin+1);
    }
    graph->SetPointError(iBin,histo->GetBinWidth(iBin+1)/2,histo->GetBinWidth(iBin+1)/2,(maxValue-minValue)/2.,(maxValue-minValue)/2.);
  }
  
  return graph;
}


////////////////
TGraphAsymmErrors* makeGraphPDF(TH1F* histo, vector<TH1F*> variations){
  
  TGraphAsymmErrors* graph = new TGraphAsymmErrors();
  TH1F* histogram = (TH1F*) histo->Clone("histogram");
  
  for(int iBin = 0; iBin < histo->GetNbinsX(); iBin++){
    graph->SetPoint(iBin,histo->GetBinCenter(iBin+1),histo->GetBinContent(iBin+1));
    vector<double> intervalAbove;
    vector<double> intervalBelow;
    histogram->Reset();
    for(auto hist: variations){
      if(hist->GetBinContent(iBin+1) > histo->GetBinContent(iBin+1)){
	intervalAbove.push_back(hist->GetBinContent(iBin+1)-histo->GetBinContent(iBin+1));
      }
      else{
	intervalBelow.push_back(hist->GetBinContent(iBin+1)-histo->GetBinContent(iBin+1));
      }
    }
    sort(intervalAbove.begin(),intervalAbove.end());
    sort(intervalBelow.begin(),intervalBelow.end());
    if(intervalBelow.size() != 0 and intervalAbove.size() != 0){     
      graph->SetPointError(iBin,histo->GetBinWidth(iBin+1)/2,histo->GetBinWidth(iBin+1)/2,
			   2*fabs(intervalBelow.at(int(0.66*intervalBelow.size())-1)),
			   2*fabs(intervalAbove.at(int(0.34*intervalAbove.size())+1)));    
    }
    else if(intervalBelow.size() == 0)
      graph->SetPointError(iBin,histo->GetBinWidth(iBin+1)/2,histo->GetBinWidth(iBin+1)/2,0.,fabs(intervalAbove.at(int(0.34*intervalAbove.size())+1)));    
    else if(intervalAbove.size() == 0)
      graph->SetPointError(iBin,histo->GetBinWidth(iBin+1)/2,histo->GetBinWidth(iBin+1)/2,fabs(intervalBelow.at(int(0.66*intervalBelow.size())-1)),0.);        
  }
    
  return graph;
}


//////////////////////////
void makeQCDPlot (TH1F* histo_inclusive,
		  TGraphAsymmErrors* graph_inclusive,
		  TH1F* histo_exclusive,
		  TGraphAsymmErrors* graph_exclusive, 
		  const string & outputDIR, const string & postfix){

  TCanvas* canvas = new TCanvas(Form("canvas_%s",postfix.c_str()),"",600,625);
  canvas->cd();
  canvas->SetBottomMargin(0.3);

  TPad *pad = new TPad(("pad_"+postfix).c_str(),"pad",0,0.,1,0.9);
  pad->SetTopMargin(0.7);
  pad->SetFillColor(0);
  pad->SetFillStyle(0);  
  canvas->cd();


  histo_inclusive->GetXaxis()->SetTitle("m_{jj} [GeV]");
  histo_inclusive->GetYaxis()->SetTitle("Events / GeV");  
  histo_inclusive->GetXaxis()->SetTitleOffset(1.20);
  histo_inclusive->GetYaxis()->SetTitleOffset(1.35);
  histo_inclusive->SetLineColor(kBlack);
  histo_inclusive->SetLineWidth(2);

  histo_inclusive->Draw("hist");
  CMS_lumi(canvas,Form("%.1f",luminosity));
  graph_inclusive->SetFillColor(kBlack);
  graph_inclusive->SetFillStyle(3001);
  graph_inclusive->SetMarkerSize(0);
  graph_inclusive->Draw("2same");
  histo_inclusive->Draw("hist same");


  histo_exclusive->Draw("hist same");
  histo_exclusive->SetLineColor(kRed);
  histo_exclusive->SetLineWidth(2);
  graph_exclusive->SetFillColor(kRed);
  graph_exclusive->SetFillStyle(3001);
  graph_exclusive->SetMarkerSize(0);
  graph_exclusive->Draw("2same");
  histo_exclusive->Draw("hist same");

  histo_inclusive->GetXaxis()->SetTitleSize(0);
  histo_inclusive->GetXaxis()->SetLabelSize(0);
  histo_inclusive->GetYaxis()->SetRangeUser(min(histo_inclusive->GetMinimum(),histo_exclusive->GetMinimum())*0.5,max(histo_inclusive->GetMaximum(),histo_exclusive->GetMaximum())*5.0);
  canvas->SetLogy();

  TLegend* leg = new TLegend (0.5,0.7,0.92,0.92);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->AddEntry(histo_inclusive,"Inclusive","L");
  leg->AddEntry(histo_exclusive,"After VBF selections","L");
  leg->Draw();


  pad->Draw();
  pad->cd();
  TH1* frame = (TH1*) histo_inclusive->Clone("frame");
  frame->Reset();

  /////////////////
  TGraphAsymmErrors* ratio = (TGraphAsymmErrors*) graph_exclusive->Clone("ratio");
  for(int iPoint = 0; iPoint < ratio->GetN(); iPoint++){
    double x,y1, y2;
    graph_exclusive->GetPoint(iPoint,x,y1);
    graph_inclusive->GetPoint(iPoint,x,y2);
    ratio->SetPoint(iPoint,x,1.);
    ratio->SetPointError(iPoint,ratio->GetErrorXlow(iPoint),ratio->GetErrorXhigh(iPoint),
			 sqrt(fabs((pow(graph_exclusive->GetErrorYlow(iPoint)/y1,2)-pow(graph_inclusive->GetErrorYlow(iPoint)/y2,2)))),
			 sqrt(fabs(pow(graph_exclusive->GetErrorYhigh(iPoint)/y1,2)-pow(graph_inclusive->GetErrorYhigh(iPoint)/y2,2))));
  }
  
  TF1* line = new TF1("line","1",frame->GetXaxis()->GetXmin(),frame->GetXaxis()->GetXmax());
  line->SetLineColor(kBlack);
  line->SetLineWidth(2);

  frame->GetXaxis()->SetTitle("m_{jj} [GeV]");
  frame->GetYaxis()->SetTitle("Uncertainty");
  frame->GetYaxis()->SetTitleOffset(1.40);
  frame->SetLineColor(kBlack);
  frame->SetMarkerStyle(20);
  frame->SetMarkerSize(1);
  frame->SetLineWidth(2);
  frame->GetXaxis()->SetLabelSize(0.04);
  frame->GetXaxis()->SetTitleSize(0.04);
  frame->GetYaxis()->SetRangeUser(0.95,1.05);
  frame->GetYaxis()->SetNdivisions(504);
  frame->Draw();
  ratio->Draw("2same");
  line->Draw("Lsame");
  
  
  canvas->SaveAs((outputDIR+"/envelope_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/envelope_"+postfix+".pdf").c_str(),"pdf");
}




//////////////////////////
void makePDFEigen (TH1F* histo_inclusive,
		  vector<TH1F*> histo_inclusive_PDF,
		  TH1F* histo_exclusive,
		  vector<TH1F*> histo_exclusive_PDF,
		  const string & outputDIR, 
		  const string & postfix){
  
  histo_inclusive->GetXaxis()->SetTitle("m_{jj} [GeV]");
  histo_inclusive->GetYaxis()->SetTitle("Events / GeV");  
  histo_inclusive->GetXaxis()->SetTitleOffset(1.20);
  histo_inclusive->GetYaxis()->SetTitleOffset(1.35);
  histo_inclusive->SetLineColor(kBlack);
  histo_inclusive->SetLineWidth(2);
  histo_exclusive->SetLineColor(kRed);
  histo_exclusive->SetLineWidth(2);

  TLegend*leg = new TLegend(0.5,0.6,0.9,0.9);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  histo_inclusive->GetYaxis()->SetRangeUser(min(histo_inclusive->GetMinimum(),histo_exclusive->GetMinimum())*0.5,max(histo_inclusive->GetMaximum(),histo_exclusive->GetMaximum())*5.0);
  
  for(int iPDF = 0; iPDF < histo_inclusive_PDF.size(); iPDF++){
    
    TCanvas* canvas = new TCanvas(Form("canvas_%s",postfix.c_str()),"",600,625);
    canvas->cd();
    canvas->SetBottomMargin(0.3);
    
    TPad *pad = new TPad(("pad_"+postfix).c_str(),"pad",0,0.,1,0.9);
    pad->SetTopMargin(0.7);
    pad->SetFillColor(0);
    pad->SetFillStyle(0);  
    canvas->cd();


    histo_inclusive->Draw("hist");    
    histo_inclusive_PDF.at(iPDF)->SetLineColor(kBlack);
    histo_inclusive_PDF.at(iPDF)->SetLineWidth(2);
    histo_inclusive_PDF.at(iPDF)->SetLineStyle(7);
    histo_inclusive_PDF.at(iPDF)->Draw("hist same");
    histo_exclusive->Draw("hist same");
    histo_exclusive_PDF.at(iPDF)->SetLineColor(kRed);
    histo_exclusive_PDF.at(iPDF)->SetLineWidth(2);
    histo_exclusive_PDF.at(iPDF)->SetLineStyle(7);
    histo_exclusive_PDF.at(iPDF)->Draw("hist same");
    leg->Clear();
    leg->AddEntry(histo_inclusive,"Inclusive","L");
    leg->AddEntry(histo_inclusive_PDF.at(iPDF),("NNPDF Eign = "+to_string(iPDF+1)).c_str(),"L");    
    leg->AddEntry(histo_exclusive,"After VBF selections","L");
    leg->AddEntry(histo_exclusive_PDF.at(iPDF),("NNPDF Eign = "+to_string(iPDF+1)).c_str(),"L");
    leg->Draw();
    CMS_lumi(canvas,Form("%.1f",luminosity));
    canvas->SetLogy();
    pad->Draw();
    pad->cd();

    TH1F* ratio_inclusive = (TH1F*) histo_inclusive_PDF.at(iPDF)->Clone("ratio_inclusive");
    ratio_inclusive->Divide(histo_inclusive);
    TH1F* ratio_exclusive = (TH1F*) histo_exclusive_PDF.at(iPDF)->Clone("ratio_exclusive");
    ratio_exclusive->Divide(histo_exclusive);
    ratio_inclusive->GetXaxis()->SetTitle("m_{jj} [GeV]");
    ratio_inclusive->GetYaxis()->SetTitle("Ratio");
    ratio_inclusive->GetYaxis()->SetTitleOffset(1.40);

    ratio_inclusive->GetYaxis()->SetNdivisions(504);
    ratio_inclusive->GetYaxis()->SetRangeUser(0.95,1.05);
    ratio_inclusive->Draw("hist");
    ratio_exclusive->Draw("hist same");
    
    canvas->SaveAs((outputDIR+"/"+postfix+"_"+to_string(iPDF)+".png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/"+postfix+"_"+to_string(iPDF)+".pdf").c_str(),"pdf");

    if(pad) delete pad;
    if(canvas) delete canvas;

  }

}


vector<float> binsPt  = {200.,225.,250,275.,300,350,400,450,500,600,700,800,1000,1250};
vector<float> binsMjj = {300.,600.,900,1200,1600.,2000.,2500,3000}; 

void  makeHiggsAcceptanceUncertainty(string inputFileName, Category category, string outputDIR, bool plotEigen = false){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TTree* tree = (TTree*) inputFile->Get("gentree/tree");

  float metCut, njetCut, jetpt1Cut, jeteta1Cut, jetpt2Cut, jeteta2Cut, detajjCut, dphijjCut, mjjCut; 

  if(category == Category::VBF or category == Category::VBFrelaxed){
    metCut      = 200;
    njetCut     = 2;
    jetpt1Cut   = 80;
    jetpt2Cut   = 40;
    jeteta1Cut  = 4.7;
    jeteta2Cut  = 4.7;
    dphijjCut   = 1.5;
    if(category == Category::VBFrelaxed){
      detajjCut = 1.0;
      mjjCut    = 300;
    }
    else if(category == Category::VBF){
      detajjCut = 4.0;
      mjjCut    = 1400;
    }
  }
  else{ 
    cerr<<"Code works only for VBF categories "<<endl;
  }

  // set branches
  TTreeReader reader (tree);
  TTreeReaderValue<float> xsec (reader,"xsec");
  TTreeReaderValue<float> wgt  (reader,"wgt");
  ///////////
  TTreeReaderValue<unsigned int>   njetsinc  (reader,"njetsinc");
  TTreeReaderValue<float> met       (reader,"met");
  //////////
  TTreeReaderValue<float> dmmass  (reader,"dmmass");
  TTreeReaderValue<float> dmpt    (reader,"dmpt");
  TTreeReaderValue<float> dmeta   (reader,"dmeta");
  TTreeReaderValue<float> dmphi   (reader,"dmphi");
  TTreeReaderValue<int>   dmid    (reader,"dmid");
  //////////
  TTreeReaderValue<vector<float> > jetpt  (reader,"jetpt");
  TTreeReaderValue<vector<float> > jeteta  (reader,"jeteta");
  TTreeReaderValue<vector<float> > jetphi  (reader,"jetphi");
  TTreeReaderValue<vector<float> > jetm  (reader,"jetmass");
  ///////////
  TTreeReaderValue<vector<float> > qcdscale_muR  (reader,"qcdscale_muR");
  TTreeReaderValue<vector<float> > qcdscale_muF  (reader,"qcdscale_muF");
  TTreeReaderValue<vector<float> > wgtqcd  (reader,"wgtqcd");
  TTreeReaderValue<vector<int> >   pdflhaid  (reader,"pdflhaid");
  TTreeReaderValue<vector<float> > wgtpdf  (reader,"wgtpdf");

  // calculate sum of weights
  double sumwgt = 0;
  while(reader.Next()){
    if(*dmid != 25){
      cerr<<"Higgs boson not found --> please check "<<endl;
      continue;
    }
    sumwgt += *wgt;
  }
  cout<<"Sum of weights "<<sumwgt<<endl;
  
  theoryVariation mjj_higgs_inclusive;
  theoryVariation mjj_higgs_exclusive;

  // Event loop
  reader.SetEntry(0);
  cout<<"Loop on events "<<endl;
  ///////////////////////////////////////////////
  while(reader.Next()){
    /////////////////
    if(*dmid != 25){
      cerr<<"Higgs boson not found --> please check "<<endl;
      continue;
    }

    // jet requriements
    if(*njetsinc < njetCut) continue;
    if(jetpt->size() < 2) continue;

    ///// CREATE histograms
    if(mjj_higgs_inclusive.histoQCDScale.size() == 0 and mjj_higgs_inclusive.histoPDF_1.size() == 0 and mjj_higgs_inclusive.histoPDF_2.size() == 0 and mjj_higgs_inclusive.histoPDF_3.size() == 0){

      /////////////////
      mjj_higgs_inclusive.histoCentral = new TH1F("mjj_higgs_inclusive_central","",binsMjj.size()-1,&binsMjj[0]);
      for(int iqcd = 0 ; iqcd < qcdscale_muR->size(); iqcd++){
	if(qcdscale_muR->at(iqcd) == qcdscale_muF->at(iqcd) and qcdscale_muR->at(iqcd) == 1) continue; // central value
	if(qcdscale_muR->at(iqcd) == 0.5 and qcdscale_muF->at(iqcd) == 2.0) continue; // central value
	if(qcdscale_muR->at(iqcd) == 2.0 and qcdscale_muF->at(iqcd) == 0.5) continue; // central value
	mjj_higgs_inclusive.histoQCDScale.push_back(new TH1F(Form("mjj_higgs_inclusive_muR_%1f_mf_%1f",qcdscale_muR->at(iqcd),qcdscale_muF->at(iqcd)),"",binsMjj.size()-1,&binsMjj[0]));
      }

      /////////////////
      for(int ipdf = 0; ipdf < pdflhaid->size(); ipdf++){
	if(pdflhaid->at(ipdf) >= 260000 and pdflhaid->at(ipdf) <= 267600)// NNPDF NLO 
	  mjj_higgs_inclusive.histoPDF_1.push_back(new TH1F(Form("mjj_higgs_inclusive_NNPDF_%d",pdflhaid->at(ipdf)),"",binsMjj.size()-1,&binsMjj[0]));
	else if(pdflhaid->at(ipdf) >= 11000 and pdflhaid->at(ipdf) < 11080 ) //CT10
	  mjj_higgs_inclusive.histoPDF_2.push_back(new TH1F(Form("mjj_higgs_inclusive_CT10_%d",pdflhaid->at(ipdf)),"",binsMjj.size()-1,&binsMjj[0]));
	else if(pdflhaid->at(ipdf) >= 25200 and pdflhaid->at(ipdf) <= 25300 ) // MMHT
	  mjj_higgs_inclusive.histoPDF_3.push_back(new TH1F(Form("mjj_higgs_inclusive_MMHT_%d",pdflhaid->at(ipdf)),"",binsMjj.size()-1,&binsMjj[0]));
      }	
    }
    
    // fill inclusive
    TLorentzVector jet1;
    TLorentzVector jet2;
    jet1.SetPtEtaPhiM(jetpt->at(0),jeteta->at(0),jetphi->at(0),jetm->at(0));
    jet2.SetPtEtaPhiM(jetpt->at(1),jeteta->at(1),jetphi->at(1),jetm->at(1));

    /////////////////
    float mjjval = (jet1+jet2).M();
    if(mjjval > binsMjj.back()) mjjval = binsMjj.back()-1; // overflow 
    mjj_higgs_inclusive.histoCentral->Fill(mjjval,*xsec*luminosity*(*wgt)/sumwgt);
    
    int ipos = 0;
    for(int iqcd = 0 ; iqcd < qcdscale_muR->size(); iqcd++){
      if(qcdscale_muR->at(iqcd) == qcdscale_muF->at(iqcd) and qcdscale_muR->at(iqcd) == 1) continue;
      if(qcdscale_muR->at(iqcd) == 0.5 and qcdscale_muF->at(iqcd) == 2.0) continue; // central value                                                                                                 
      if(qcdscale_muR->at(iqcd) == 2.0 and qcdscale_muF->at(iqcd) == 0.5) continue; // central value                                                                                                 
       mjj_higgs_inclusive.histoQCDScale.at(ipos)->Fill(mjjval,*xsec*luminosity*wgtqcd->at(ipos)/sumwgt);
      ipos++;
    }

    ipos = 0;
    for(int ipdf = 0 ; ipdf < pdflhaid->size(); ipdf++){
      if(pdflhaid->at(ipdf) >= 260000 and pdflhaid->at(ipdf) <= 267600){ //NNPDF                                                                                                                     
	if(ipos >= mjj_higgs_inclusive.histoPDF_1.size()) ipos = 0;
	mjj_higgs_inclusive.histoPDF_1.at(ipos)->Fill(mjjval,*xsec*luminosity*(wgtpdf->at(ipdf))/sumwgt);	
	ipos++;
      }
      else if(pdflhaid->at(ipdf) >= 11000 and pdflhaid->at(ipdf) < 11080 ){ // CT10                                                                                                                  
	if(ipos >= mjj_higgs_inclusive.histoPDF_2.size()) ipos = 0;
	mjj_higgs_inclusive.histoPDF_2.at(ipos)->Fill(mjjval,*xsec*luminosity*(wgtpdf->at(ipdf))/sumwgt);	
	ipos++;
      }
      else if(pdflhaid->at(ipdf) >= 25200 and pdflhaid->at(ipdf) <= 25300 ){ // MMHT                                                                                                                 
	if(ipos >= mjj_higgs_inclusive.histoPDF_2.size()) ipos = 0;
	mjj_higgs_inclusive.histoPDF_3.at(ipos)->Fill(mjjval,*xsec*luminosity*(wgtpdf->at(ipdf))/sumwgt);
        ipos++;	
      }
    }  


    /////////////////
    // dmpt cut
    /////////////////
    if(*dmpt < metCut) continue;    
    if(jetpt->at(0) < jetpt1Cut) continue;
    if(fabs(jeteta->at(0)) > jeteta1Cut) continue;

 
    /////////////////   
    if(category == Category::VBF or category == Category::VBFrelaxed){
      if(jetpt->at(1) < jetpt2Cut) continue;
      if(fabs(jeteta->at(1)) > jeteta2Cut) continue;
      if(jeteta->at(1)*jeteta->at(0) > 0) continue;
      if(fabs(jeteta->at(0)-jeteta->at(1)) < detajjCut) continue;
      if(fabs(jeteta->at(1)) > 3 and fabs(jeteta->at(0)) > 3) continue;
      if(fabs(jet1.DeltaPhi(jet2)) > dphijjCut) continue;
      if((jet1+jet2).M() < mjjCut) continue;      
    }

    ///// CREATE HISTOGRAMS
    if(mjj_higgs_exclusive.histoQCDScale.size() == 0 and mjj_higgs_exclusive.histoPDF_1.size() == 0 and mjj_higgs_exclusive.histoPDF_2.size() == 0 and mjj_higgs_exclusive.histoPDF_3.size() == 0){	
      mjj_higgs_exclusive.histoCentral = new TH1F("mjj_higgs_exclusive_central","",binsMjj.size()-1,&binsMjj[0]);
      /////////
      for(int iqcd = 0 ; iqcd < qcdscale_muR->size(); iqcd++){
	if(qcdscale_muR->at(iqcd) == qcdscale_muF->at(iqcd) and qcdscale_muR->at(iqcd) == 1) continue; // central value
	if(qcdscale_muR->at(iqcd) == 0.5 and qcdscale_muF->at(iqcd) == 2.0) continue; // central value                                                                                                
        if(qcdscale_muR->at(iqcd) == 2.0 and qcdscale_muF->at(iqcd) == 0.5) continue; // central value                                                                                                
 	mjj_higgs_exclusive.histoQCDScale.push_back(new TH1F(Form("mjj_higgs_exclusive_muR_%1f_mf_%1f",qcdscale_muR->at(iqcd),qcdscale_muF->at(iqcd)),"",binsMjj.size()-1,&binsMjj[0]));
      }
      
      ////////
      for(int ipdf = 0; ipdf < pdflhaid->size(); ipdf++){
	if(pdflhaid->at(ipdf) >= 260000 and pdflhaid->at(ipdf) <= 267600) //NNPDF
	  mjj_higgs_exclusive.histoPDF_1.push_back(new TH1F(Form("mjj_higgs_exclusive_NNPDF_%d",pdflhaid->at(ipdf)),"",binsMjj.size()-1,&binsMjj[0]));
	else if(pdflhaid->at(ipdf) >= 11000 and pdflhaid->at(ipdf) < 11080 ) // CT10
	  mjj_higgs_exclusive.histoPDF_2.push_back(new TH1F(Form("mjj_higgs_exclusive_CT10_%d",pdflhaid->at(ipdf)),"",binsMjj.size()-1,&binsMjj[0]));
	else if(pdflhaid->at(ipdf) >= 25200 and pdflhaid->at(ipdf) <= 25300 ) // MMHT
	  mjj_higgs_exclusive.histoPDF_3.push_back(new TH1F(Form("mjj_higgs_exclusive_MMHT_%d",pdflhaid->at(ipdf)),"",binsMjj.size()-1,&binsMjj[0]));
      }	
    }

    mjj_higgs_exclusive.histoCentral->Fill(mjjval,*xsec*luminosity*(*wgt)/sumwgt);

    ipos = 0;
    for(int iqcd = 0 ; iqcd < qcdscale_muR->size(); iqcd++){
      if(qcdscale_muR->at(iqcd) == qcdscale_muF->at(iqcd) and qcdscale_muR->at(iqcd) == 1) continue; // central value  
      if(qcdscale_muR->at(iqcd) == 0.5 and qcdscale_muF->at(iqcd) == 2.0) continue; // central value                                                                                                 
      if(qcdscale_muR->at(iqcd) == 2.0 and qcdscale_muF->at(iqcd) == 0.5) continue; // central value                                                                                                 
       mjj_higgs_exclusive.histoQCDScale.at(ipos)->Fill(mjjval,*xsec*luminosity*wgtqcd->at(ipos)/sumwgt);
      ipos++;
    }

    ipos = 0;
    for(int ipdf = 0 ; ipdf < pdflhaid->size(); ipdf++){
      if(pdflhaid->at(ipdf) >= 260000 and pdflhaid->at(ipdf) <= 267600){ //NNPDF                                                                                                                     
	if(ipos >= mjj_higgs_exclusive.histoPDF_1.size()) ipos = 0;
	mjj_higgs_exclusive.histoPDF_1.at(ipos)->Fill(mjjval,*xsec*luminosity*(wgtpdf->at(ipdf))/sumwgt);	
	ipos++;
      }
      else if(pdflhaid->at(ipdf) >= 11000 and pdflhaid->at(ipdf) < 11080 ){ // CT10                                                                                                                  
	if(ipos >= mjj_higgs_exclusive.histoPDF_2.size()) ipos = 0;
	mjj_higgs_exclusive.histoPDF_2.at(ipos)->Fill(mjjval,*xsec*luminosity*(wgtpdf->at(ipdf))/sumwgt);	
	ipos++;
      }
      else if(pdflhaid->at(ipdf) >= 25200 and pdflhaid->at(ipdf) <= 25300 ){ // MMHT                                                                                                                 
	if(ipos >= mjj_higgs_exclusive.histoPDF_3.size()) ipos = 0;
	mjj_higgs_exclusive.histoPDF_3.at(ipos)->Fill(mjjval,*xsec*luminosity*(wgtpdf->at(ipdf))/sumwgt);
        ipos++;	
      }
    }
  }  

  /// outpout
  TFile* outputFile = new TFile((outputDIR+"/outputHiggsDistributions.root").c_str(),"RECREATE");
  outputFile->cd();
  outputFile->cd();
  outputFile->mkdir("higgsMjj_inclusive");
  outputFile->cd("higgsMjj_inclusive");
  mjj_higgs_inclusive.histoCentral->Write();
  for(auto hist : mjj_higgs_inclusive.histoQCDScale)
    hist->Write();
  for(auto hist : mjj_higgs_inclusive.histoPDF_1)
    hist->Write();
  for(auto hist : mjj_higgs_inclusive.histoPDF_2)
    hist->Write();
  for(auto hist : mjj_higgs_inclusive.histoPDF_3)
    hist->Write();
  outputFile->cd();
  
  //////////////
  setTDRStyle();
  /*
  // Scale histograns
  mjj_higgs_inclusive.histoCentral->Scale(1./mjj_higgs_inclusive.histoCentral->Integral());
  mjj_higgs_exclusive.histoCentral->Scale(1./mjj_higgs_exclusive.histoCentral->Integral());
  for(auto hist : mjj_higgs_inclusive.histoQCDScale)
    hist->Scale(1./hist->Integral());
  for(auto hist : mjj_higgs_inclusive.histoPDF_1)
    hist->Scale(1./hist->Integral());
  for(auto hist : mjj_higgs_inclusive.histoPDF_2)
    hist->Scale(1./hist->Integral());
  for(auto hist : mjj_higgs_inclusive.histoPDF_3)
    hist->Scale(1./hist->Integral());
  for(auto hist : mjj_higgs_exclusive.histoQCDScale)
    hist->Scale(1./hist->Integral());
  for(auto hist : mjj_higgs_exclusive.histoPDF_1)
    hist->Scale(1./hist->Integral());
  for(auto hist : mjj_higgs_exclusive.histoPDF_2)
    hist->Scale(1./hist->Integral());
  for(auto hist : mjj_higgs_exclusive.histoPDF_3)
    hist->Scale(1./hist->Integral());
  */

  //// -- make the QCD scale variation plot  
  TGraphAsymmErrors* graphQCD_inclusive = makeGraphQCD(mjj_higgs_inclusive.histoCentral,mjj_higgs_inclusive.histoQCDScale);
  TGraphAsymmErrors* graphQCD_exclusive = makeGraphQCD(mjj_higgs_exclusive.histoCentral,mjj_higgs_exclusive.histoQCDScale);
  ////
  makeQCDPlot(mjj_higgs_inclusive.histoCentral,graphQCD_inclusive,mjj_higgs_exclusive.histoCentral,graphQCD_exclusive,outputDIR,"mjj_qcd_scale");
  ////
  if(plotEigen)
    makePDFEigen(mjj_higgs_inclusive.histoCentral,mjj_higgs_inclusive.histoPDF_1,mjj_higgs_exclusive.histoCentral,mjj_higgs_exclusive.histoPDF_1,outputDIR,"NNPDF_Eigen");
  /// 
  TGraphAsymmErrors* graphPDF_inclusive = makeGraphPDF(mjj_higgs_inclusive.histoCentral,mjj_higgs_inclusive.histoPDF_1);
  TGraphAsymmErrors* graphPDF_exclusive = makeGraphPDF(mjj_higgs_exclusive.histoCentral,mjj_higgs_exclusive.histoPDF_1); 
  makeQCDPlot(mjj_higgs_inclusive.histoCentral,graphPDF_inclusive,mjj_higgs_exclusive.histoCentral,graphPDF_exclusive,outputDIR,"mjj_NNPDF");
}
