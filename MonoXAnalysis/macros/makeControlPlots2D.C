#include "CMS_lumi.h"
#include "makehist.h"

void makeControlPlots2D(string templateFileName, 
			int    category, 
			string observable, 
			string observableLatex, 
			string controlRegion, 
			bool   blind, 
			bool   isLog,
			bool   plotResonant   = false,
			bool   isHiggsInvisible = false,
			bool   alongX = false,
			string interaction  = "Vector",
			string mediatorMass = "1000",
			string DMMass       = "50",
			int    signalScale     = 100) {

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  gStyle->SetOptStat(0);

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
  canvas->SetTickx();
  canvas->SetTicky();
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);

  setTDRStyle();

  TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
  pad1->SetTickx();
  pad1->SetTicky();

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.28);
  pad2->SetTickx();
  pad2->SetTicky();

  TFile* inputFile = new TFile(templateFileName.c_str());

  vector<TH1F*> datahist;
  vector<TH1F*> qcdhist ;
  vector<TH1F*> vllhist ;
  vector<TH1F*> vnnhist ;
  vector<TH1F*> vlhist  ;
  vector<TH1F*> dbhist  ;
  vector<TH1F*> tophist ;
  vector<TH1F*> tophist_matched  ;
  vector<TH1F*> tophist_unmatched;
  vector<TH1F*> gamhist ;

  vector<TH1F*> monoJhist ;
  vector<TH1F*> monoWhist ;
  vector<TH1F*> monoZhist ;

  vector<TH1F*> ggHhist ;
  vector<TH1F*> vbfHhist;
  vector<TH1F*> wHhist  ;
  vector<TH1F*> zHhist  ;
  
  // take the templates
  if(controlRegion == "gam"){  
    datahist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("datahistgam_"+observable).c_str()),observable,category,alongX);
    qcdhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("qbkghistgam_"+observable).c_str()),observable,category,alongX);
    gamhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("gbkghistgam_"+observable).c_str()),observable,category,alongX);
  }
  else if(controlRegion == "zmm"){
    datahist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("datahistzmm_"+observable).c_str()),observable,category,alongX);
    qcdhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("qbkghistzmm_"+observable).c_str()),observable,category,alongX);
    gamhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("gbkghistzmm_"+observable).c_str()),observable,category,alongX);
    tophist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("tbkghistzmm_"+observable).c_str()),observable,category,alongX);
    vlhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vlbkghistzmm_"+observable).c_str()),observable,category,alongX);
    vllhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vllbkghistzmm_"+observable).c_str()),observable,category,alongX);
    dbhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("dbkghistzmm_"+observable).c_str()),observable,category,alongX);
  }
  else if(controlRegion == "zee"){
    datahist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("datahistzee_"+observable).c_str()),observable,category,alongX);
    qcdhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("qbkghistzee_"+observable).c_str()),observable,category,alongX);
    gamhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("gbkghistzee_"+observable).c_str()),observable,category,alongX);
    tophist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("tbkghistzee_"+observable).c_str()),observable,category,alongX);
    vlhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vlbkghistzee_"+observable).c_str()),observable,category,alongX);
    vllhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vllbkghistzee_"+observable).c_str()),observable,category,alongX);
    dbhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("dbkghistzee_"+observable).c_str()),observable,category,alongX);
  }
  else if(controlRegion == "wmn"){
    datahist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("datahistwmn_"+observable).c_str()),observable,category,alongX);
    qcdhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("qbkghistwmn_"+observable).c_str()),observable,category,alongX);
    gamhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("gbkghistwmn_"+observable).c_str()),observable,category,alongX);
    tophist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("tbkghistwmn_"+observable).c_str()),observable,category,alongX);
    vlhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vlbkghistwmn_"+observable).c_str()),observable,category,alongX);
    vllhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vllbkghistwmn_"+observable).c_str()),observable,category,alongX);
    dbhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("dbkghistwmn_"+observable).c_str()),observable,category,alongX);
  }
  else if(controlRegion == "wen"){
    datahist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("datahistwen_"+observable).c_str()),observable,category,alongX);
    qcdhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("qbkghistwen_"+observable).c_str()),observable,category,alongX);
    gamhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("gbkghistwen_"+observable).c_str()),observable,category,alongX);
    tophist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("tbkghistwen_"+observable).c_str()),observable,category,alongX);
    vlhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vlbkghistwen_"+observable).c_str()),observable,category,alongX);
    vllhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vllbkghistwen_"+observable).c_str()),observable,category,alongX);
    dbhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("dbkghistwen_"+observable).c_str()),observable,category,alongX);
  }
  else if(controlRegion == "topmu" and plotResonant){    
    datahist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("datahisttopmu_"+observable).c_str()),observable,category,alongX);
    qcdhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("qbkghisttopmu_"+observable).c_str()),observable,category,alongX);
    tophist_matched   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("tbkghist_matchedtopmu_"+observable).c_str()),observable,category,alongX);
    tophist_unmatched = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("tbkghist_unmatchedtopmu_"+observable).c_str()),observable,category,alongX);
    vlhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vlbkghisttopmu_"+observable).c_str()),observable,category,alongX);
    vllhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vllbkghisttopmu_"+observable).c_str()),observable,category,alongX);
    dbhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("dbkghisttopmu_"+observable).c_str()),observable,category,alongX);
    gamhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("gbkghisttopmu_"+observable).c_str()),observable,category,alongX);
  }
  else if(controlRegion == "topmu" and not plotResonant){
    datahist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("datahisttopmu_"+observable).c_str()),observable,category,alongX);
    qcdhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("qbkghisttopmu_"+observable).c_str()),observable,category,alongX);
    tophist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("tbkghisttopmu_"+observable).c_str()),observable,category,alongX);
    vlhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vlbkghisttopmu_"+observable).c_str()),observable,category,alongX);
    vllhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vllbkghisttopmu_"+observable).c_str()),observable,category,alongX);
    dbhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("dbkghisttopmu_"+observable).c_str()),observable,category,alongX);
    gamhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("gbkghisttopmu_"+observable).c_str()),observable,category,alongX);
  }

  else if(controlRegion == "topel" and not plotResonant){
    datahist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("datahisttopel_"+observable).c_str()),observable,category,alongX);
    qcdhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("qbkghisttopel_"+observable).c_str()),observable,category, alongX);
    tophist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("tbkghisttopel_"+observable).c_str()),observable,category,alongX);
    vlhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vlbkghisttopel_"+observable).c_str()),observable,category,alongX);
    vllhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vllbkghisttopel_"+observable).c_str()),observable,category,alongX);
    dbhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("dbkghisttopel_"+observable).c_str()),observable,category,alongX);
    gamhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("gbkghisttopel_"+observable).c_str()),observable,category,alongX);
  }
  else if(controlRegion == "topel" and plotResonant){
    datahist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("datahisttopel_"+observable).c_str()),observable,category,alongX);
    qcdhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("qbkghisttopel_"+observable).c_str()),observable,category,alongX);
    tophist_matched   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("tbkghist_matchedtopel_"+observable).c_str()),observable,category,alongX);
    tophist_unmatched = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("tbkghist_unmatchedtopel_"+observable).c_str()),observable,category,alongX);
    vlhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vlbkghisttopel_"+observable).c_str()),observable,category,alongX);
    vllhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vllbkghisttopel_"+observable).c_str()),observable,category,alongX);
    dbhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("dbkghisttopel_"+observable).c_str()),observable,category,alongX);
    gamhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("gbkghisttopel_"+observable).c_str()),observable,category,alongX);
  }

  else if(controlRegion == "SR"){

    datahist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("datahist_"+observable).c_str()),observable,category,alongX);
    qcdhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("qbkghist_"+observable).c_str()),observable,category,alongX);
    tophist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("tbkghist_"+observable).c_str()),observable,category,alongX);
    vlhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("wjethist_"+observable).c_str()),observable,category,alongX);
    vllhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("zjethist_"+observable).c_str()),observable,category,alongX);
    dbhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("dbkghist_"+observable).c_str()),observable,category,alongX);
    gamhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("gbkghist_"+observable).c_str()),observable,category,alongX);
    vnnhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("zinvhist_"+observable).c_str()),observable,category,alongX);

    if(not isHiggsInvisible){
      monoJhist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()),observable,category,alongX);
      monoWhist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()),observable,category,alongX);
      monoZhist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()),observable,category,alongX);  
    }
    else{
      ggHhist  = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+observable).c_str()),observable,category,alongX);
      vbfHhist = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("vbfHhist_"+mediatorMass+"_"+observable).c_str()),observable,category,alongX);
      wHhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("wHhist_"+mediatorMass+"_"+observable).c_str()),observable,category,alongX);    
      zHhist   = transformUnrolledHistogram((TH1*)inputFile->FindObjectAny(("zHhist_"+mediatorMass+"_"+observable).c_str()),observable,category,alongX);    
    }
  }
  
  bin2D bins = selectBinning2D(observable,category);
  vector<double> bin ;
  if(alongX)
    bin = bins.binX ;
  else
    bin = bins.binY ;

  if(bin.size()-1 != datahist.size())
    cerr<<"Huston we have a problem with bin size ... "<<endl;

  pair<string,string> text = observableName(observable,alongX);

  for(size_t ihisto = 0; ihisto < datahist.size(); ihisto++){

    canvas->cd();  
    //SCALE BIN WIDTH
    if(TString(observableLatex).Contains("GeV")){
      
      if(controlRegion == "SR" and qcdhist.size() > ihisto and not TString(qcdhist.at(ihisto)->GetName()).Contains("qbkghistDD"))
	qcdhist.at(ihisto)->Scale(2.);

      if(datahist.size() > ihisto) datahist.at(ihisto)->Scale(1.0,"width");
      if(qcdhist.size() > ihisto)  qcdhist.at(ihisto)->Scale(1.0,"width");
      if(tophist.size() > ihisto)  tophist.at(ihisto)->Scale(1.0,"width");
      if(tophist_matched.size() > ihisto)   tophist_matched.at(ihisto)->Scale(1.0,"width");
      if(tophist_unmatched.size() > ihisto) tophist_unmatched.at(ihisto)->Scale(1.0,"width");
      if(vlhist.size() > ihisto)  vlhist.at(ihisto)->Scale(1.0,"width");
      if(vllhist.size() > ihisto) vllhist.at(ihisto)->Scale(1.0,"width");
      if(vnnhist.size() > ihisto) vnnhist.at(ihisto)->Scale(1.0,"width");
      if(dbhist.size() > ihisto)  dbhist.at(ihisto)->Scale(1.0,"width");
      if(gamhist.size() > ihisto) gamhist.at(ihisto)->Scale(1.0,"width");
      
      if(monoJhist.size() > ihisto) monoJhist.at(ihisto)->Scale(1.0,"width");
      
      if(monoWhist.size() > ihisto){
	monoWhist.at(ihisto)->Scale(1.0,"width");
	monoWhist.at(ihisto)->Scale(signalScale);
      }
      if(monoZhist.size() > ihisto){
	monoZhist.at(ihisto)->Scale(1.0,"width");
	monoZhist.at(ihisto)->Scale(signalScale);
      }
      
      if(ggHhist.size() > ihisto) ggHhist.at(ihisto)->Scale(1.0,"width");    
      if(vbfHhist.size() > ihisto) vbfHhist.at(ihisto)->Scale(1.0,"width");   
      if(wHhist.size() > ihisto) wHhist.at(ihisto)->Scale(1.0,"width");    
      if(zHhist.size() > ihisto) zHhist.at(ihisto)->Scale(1.0,"width");
    }
    else{
      
      if(controlRegion == "SR" and  qcdhist.size() > ihisto and not TString(qcdhist.at(ihisto)->GetName()).Contains("qbkghistDD"))
	qcdhist.at(ihisto)->Scale(2.);
      
      if(monoWhist.size() > ihisto) monoWhist.at(ihisto)->Scale(signalScale);    
      if(monoZhist.size() > ihisto) monoZhist.at(ihisto)->Scale(signalScale);
    }

    // BLIND OPTION
    if(blind and controlRegion == "SR" and datahist.size() > ihisto){

      for (int i = 0; i <= datahist.at(ihisto)->GetNbinsX()+1; i++) {

	double yield = 0.0;
	if(qcdhist.size() > ihisto)
	  yield += qcdhist.at(ihisto)->GetBinContent(i);
	if(gamhist.size() > ihisto)
	  yield += gamhist.at(ihisto)->GetBinContent(i);
	if(tophist.size() > ihisto)
	  yield += tophist.at(ihisto)->GetBinContent(i);
	if(dbhist.size() > ihisto)
	  yield += dbhist.at(ihisto)->GetBinContent(i);
	if(vllhist.size() > ihisto)
	  yield += vllhist.at(ihisto)->GetBinContent(i);
	if(vlhist.size() > ihisto)
	  yield += vlhist.at(ihisto)->GetBinContent(i);
	if(vlhist.size() > ihisto)
	  yield += vnnhist.at(ihisto)->GetBinContent(i);
	if(datahist.size() > ihisto){
	  datahist.at(ihisto)->SetBinContent(i, yield);
	  datahist.at(ihisto)->SetBinError(i, 0.);
	}
      }
    }

    // set colors
    if(datahist.size() > ihisto){
      datahist.at(ihisto)->SetLineColor(kBlack);
      datahist.at(ihisto)->SetMarkerColor(kBlack);
      datahist.at(ihisto)->SetMarkerStyle(20);
      datahist.at(ihisto)->SetMarkerSize(1.2);
    }

    if(vnnhist.size() > ihisto){
      vnnhist.at(ihisto)->SetFillColor(kGreen+1);
      vnnhist.at(ihisto)->SetLineColor(kBlack);
    }
    if(vllhist.size() > ihisto){
      vllhist.at(ihisto)->SetFillColor(kCyan);
      vllhist.at(ihisto)->SetLineColor(kBlack);
    }
    if(vlhist.size() > ihisto){
      vlhist.at(ihisto)->SetFillColor(kRed);
      vlhist.at(ihisto)->SetLineColor(kBlack);
    }
    if(tophist.size() > ihisto){
      tophist.at(ihisto)->SetFillColor(kBlue);
      tophist.at(ihisto)->SetLineColor(kBlack);
    }
    if(tophist_matched.size() > ihisto){
      tophist_matched.at(ihisto)->SetFillColor(kGreen+1);
      tophist_matched.at(ihisto)->SetLineColor(kBlack);
    }
    if(tophist_unmatched.size() > ihisto){
      tophist_unmatched.at(ihisto)->SetFillColor(kBlue);
      tophist_unmatched.at(ihisto)->SetLineColor(kBlack);
    }
    if(dbhist.size() > ihisto){
      dbhist.at(ihisto)->SetFillColor(kViolet);
      dbhist.at(ihisto)->SetLineColor(kBlack);
    }
    if(qcdhist.size() > ihisto) {
      qcdhist.at(ihisto)->SetFillColor(kGray+1);
      qcdhist.at(ihisto)->SetLineColor(kBlack);
    }
    if(gamhist.size() > ihisto){
      gamhist.at(ihisto)->SetFillColor(kOrange);
      gamhist.at(ihisto)->SetLineColor(kBlack);
    }

    if(monoJhist.size() > ihisto){
      monoJhist.at(ihisto)->SetFillColor(0);
      monoJhist.at(ihisto)->SetFillStyle(0);
      monoJhist.at(ihisto)->SetLineColor(kBlack);
      monoJhist.at(ihisto)->SetLineWidth(2);
    }

    if(monoWhist.size() > ihisto){
      monoWhist.at(ihisto)->SetFillColor(0);
      monoWhist.at(ihisto)->SetFillStyle(0);
      monoWhist.at(ihisto)->SetLineColor(kBlack);
      monoWhist.at(ihisto)->SetLineWidth(2);
      monoWhist.at(ihisto)->SetLineStyle(7);
    }

    if(monoZhist.size() > ihisto){
      monoZhist.at(ihisto)->SetFillColor(0);
      monoZhist.at(ihisto)->SetFillStyle(0);
      monoZhist.at(ihisto)->SetLineColor(kBlack);
      monoZhist.at(ihisto)->SetLineWidth(2);
      monoZhist.at(ihisto)->SetLineStyle(4);
    }

    if(ggHhist.size() > ihisto){
      ggHhist.at(ihisto)->SetFillColor(0);
      ggHhist.at(ihisto)->SetFillStyle(0);
      ggHhist.at(ihisto)->SetLineColor(kBlack);
      ggHhist.at(ihisto)->SetLineWidth(2);
    }

    if(vbfHhist.size() > ihisto){
      vbfHhist.at(ihisto)->SetFillColor(0);
      vbfHhist.at(ihisto)->SetFillStyle(0);
      vbfHhist.at(ihisto)->SetLineColor(kBlack);
      vbfHhist.at(ihisto)->SetLineWidth(2);
      vbfHhist.at(ihisto)->SetLineStyle(7);
    }

    if(wHhist.size() > ihisto){
      wHhist.at(ihisto)->SetFillColor(0);
      wHhist.at(ihisto)->SetFillStyle(0);
      wHhist.at(ihisto)->SetLineColor(kBlack);
      wHhist.at(ihisto)->SetLineWidth(2);
      wHhist.at(ihisto)->SetLineStyle(4);
    }

    if(zHhist.size() > ihisto){
      zHhist.at(ihisto)->SetFillColor(0);
      zHhist.at(ihisto)->SetFillStyle(0);
      zHhist.at(ihisto)->SetLineColor(kBlack);
      zHhist.at(ihisto)->SetLineWidth(2);
      zHhist.at(ihisto)->SetLineStyle(2);
    }
  
    THStack* stack = new THStack("stack", "stack");
    if(controlRegion == "gam"){
      stack->Add(qcdhist.at(ihisto));
      stack->Add(gamhist.at(ihisto));
    }
    else if(controlRegion == "zmm" or controlRegion == "zee"){
      stack->Add(qcdhist.at(ihisto));
      stack->Add(gamhist.at(ihisto));
      stack->Add(vlhist.at(ihisto));
      stack->Add(tophist.at(ihisto));
      stack->Add(dbhist.at(ihisto));
      stack->Add(vllhist.at(ihisto));
    }
    else if(controlRegion == "wmn" or controlRegion == "wen"){
      stack->Add(qcdhist.at(ihisto));
      stack->Add(gamhist.at(ihisto));
      stack->Add(vllhist.at(ihisto));
      stack->Add(tophist.at(ihisto));
      stack->Add(dbhist.at(ihisto));
      stack->Add(vlhist.at(ihisto));
    }
    else if((controlRegion == "topmu" or controlRegion == "topel") and not plotResonant){
      stack->Add(qcdhist.at(ihisto));
      stack->Add(gamhist.at(ihisto));
      stack->Add(vllhist.at(ihisto));
      stack->Add(dbhist.at(ihisto));
      stack->Add(vlhist.at(ihisto));
      stack->Add(tophist.at(ihisto));    
    }
    else if((controlRegion == "topmu" or controlRegion == "topel") and plotResonant){
      stack->Add(qcdhist.at(ihisto));
      stack->Add(gamhist.at(ihisto));
      stack->Add(vllhist.at(ihisto));
      stack->Add(dbhist.at(ihisto));
      stack->Add(vlhist.at(ihisto));
      stack->Add(tophist_unmatched.at(ihisto));    
      stack->Add(tophist_matched.at(ihisto));    
    }
    else if(controlRegion == "SR"){
      stack->Add(qcdhist.at(ihisto));
      stack->Add(gamhist.at(ihisto));
      stack->Add(dbhist.at(ihisto));
      stack->Add(tophist.at(ihisto));
      stack->Add(vllhist.at(ihisto));
      stack->Add(vlhist.at(ihisto));
      stack->Add(vnnhist.at(ihisto));
    }

    if(controlRegion == "SR" and observable == "met"){
    
      // write yields in a output in a text file 
      ofstream outputfile;
      outputfile.open(Form("preFitSR_bin_%d.txt",int(ihisto)));
      
      stringstream QCDRate;
      QCDRate << "Process: QCD";
      stringstream GJetsRate;
      GJetsRate << "Process: GJets";
      stringstream DiBosonRate;
      DiBosonRate << "Process: DiBoson";
      stringstream TopRate;
      TopRate << "Process: TopRate";
      stringstream ZJetsRate;
      ZJetsRate << "Process: ZJetsRate";
      stringstream WJetsRate;
      WJetsRate << "Process: WJetsRate";
      stringstream ZnunuRate;
      ZnunuRate << "Process: ZnunuRate";
      stringstream PreRate;
      PreRate << "Process: Pre-fit (total)";
      stringstream PostRate;
      PostRate << "Process: Post-fit (total)";
      stringstream DataRate;
      DataRate << "Process: Data";
      
      for(int iBin = 0; iBin < qcdhist.at(ihisto)->GetNbinsX(); iBin++){
	QCDRate << "   ";
	QCDRate << qcdhist.at(ihisto)->GetBinContent(iBin+1);
      }
      
      for(int iBin = 0; iBin < gamhist.at(ihisto)->GetNbinsX(); iBin++){
	GJetsRate << "   ";
	GJetsRate << gamhist.at(ihisto)->GetBinContent(iBin+1);
      }
      
      for(int iBin = 0; iBin < dbhist.at(ihisto)->GetNbinsX(); iBin++){
	DiBosonRate << "   ";
	DiBosonRate << dbhist.at(ihisto)->GetBinContent(iBin+1);
      }
      
      for(int iBin = 0; iBin < tophist.at(ihisto)->GetNbinsX(); iBin++){
	TopRate << "   ";
	TopRate << tophist.at(ihisto)->GetBinContent(iBin+1);
      }
    
      for(int iBin = 0; iBin < vllhist.at(ihisto)->GetNbinsX(); iBin++){
	ZJetsRate << "   ";
	ZJetsRate << vllhist.at(ihisto)->GetBinContent(iBin+1);
      }
      
      for(int iBin = 0; iBin < vlhist.at(ihisto)->GetNbinsX(); iBin++){
	WJetsRate << "   ";
	WJetsRate << vlhist.at(ihisto)->GetBinContent(iBin+1);
      }
      
      for(int iBin = 0; iBin < vnnhist.at(ihisto)->GetNbinsX(); iBin++){
	ZnunuRate << "   ";
	ZnunuRate << vnnhist.at(ihisto)->GetBinContent(iBin+1);
      }
      
      TH1* histoTotal = (TH1*) stack->GetStack()->At(stack->GetNhists()-1);
      
      for(int iBin = 0; iBin < histoTotal->GetNbinsX(); iBin++){
	PreRate << "   ";
	PreRate << histoTotal->GetBinContent(iBin+1);
      }
      
      for(int iBin = 0; iBin < datahist.at(ihisto)->GetNbinsX(); iBin++){
	DataRate << "   ";
	DataRate << datahist.at(ihisto)->GetBinContent(iBin+1);
      }

      outputfile<<"######################"<<endl;
      outputfile<<QCDRate.str()<<endl;
      outputfile<<"######################"<<endl;
      outputfile<<GJetsRate.str()<<endl;
      outputfile<<"######################"<<endl;
      outputfile<<DiBosonRate.str()<<endl;
      outputfile<<"######################"<<endl;
      outputfile<<TopRate.str()<<endl;
      outputfile<<"######################"<<endl;
      outputfile<<ZJetsRate.str()<<endl;
      outputfile<<"######################"<<endl;
      outputfile<<WJetsRate.str()<<endl;
      outputfile<<"######################"<<endl;
      outputfile<<ZnunuRate.str()<<endl;
      outputfile<<"######################"<<endl;
      outputfile<<PreRate.str()<<endl;
      outputfile<<"######################"<<endl;
      outputfile<<PostRate.str()<<endl;
      outputfile<<"######################"<<endl;
      outputfile<<DataRate.str()<<endl;
      outputfile<<"######################"<<endl;
      
      outputfile.close();
    }


    pad1->SetRightMargin(0.06);
    pad1->SetLeftMargin(0.12);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();
    pad1->cd();

    TH1* frame = (TH1*) datahist.at(ihisto)->Clone("frame");
    frame->Reset();

    // set Y-axis range                                                                                                                                                    
    if(category <= 1 and isLog)
      frame->GetYaxis()->SetRangeUser(5e-4,datahist.at(ihisto)->GetMaximum()*500);
    else if(category <= 1 and not isLog)
      frame->GetYaxis()->SetRangeUser(5e-4,datahist.at(ihisto)->GetMaximum()*1.5);
    else if(category > 1 and isLog)
      frame->GetYaxis()->SetRangeUser(5e-4,datahist.at(ihisto)->GetMaximum()*500);
    else
      frame->GetYaxis()->SetRangeUser(5e-4,datahist.at(ihisto)->GetMaximum()*2.5);
          
    frame->GetXaxis()->SetTitle(observableLatex.c_str());
    if(TString(observableLatex).Contains("GeV"))
      frame->GetYaxis()->SetTitle("Events / GeV");
    else
      frame->GetYaxis()->SetTitle("Events");

    frame->GetXaxis()->SetTitleSize(0);
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->Draw();

    if(controlRegion == "SR")
      CMS_lumi(pad1,"2.30",false,true);
    else
      CMS_lumi(pad1,"2.30");

    stack ->Draw("HIST SAME");
    datahist.at(ihisto)->Draw("PE SAME");

    if(controlRegion == "SR" and not isHiggsInvisible){
      monoJhist.at(ihisto)->Draw("hist same");
      monoWhist.at(ihisto)->Draw("hist same");
      monoZhist.at(ihisto)->Draw("hist same");
    }
    else if(controlRegion == "SR" and isHiggsInvisible){
      ggHhist.at(ihisto)->Draw("hist same");
      vbfHhist.at(ihisto)->Draw("hist same");
      wHhist.at(ihisto)->Draw("hist same");
      zHhist.at(ihisto)->Draw("hist same");
    }

    TLegend* leg = NULL;
    if(controlRegion == "gam")
      leg = new TLegend(0.62, 0.70, 0.85, 0.90);
    else if(controlRegion == "SR" and isLog)
      leg = new TLegend(0.42, 0.50, 0.88, 0.90);
    else if(controlRegion == "SR" and not isLog)
      leg = new TLegend(0.42, 0.50, 0.88, 0.90);
    else
      leg = new TLegend(0.62, 0.50, 0.85, 0.90);
    
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
      
    if(controlRegion == "gam"){
      leg->AddEntry(datahist.at(ihisto), "Data","PLE");
      leg->AddEntry(gamhist.at(ihisto), "#gamma+jets","F");
      leg->AddEntry(qcdhist.at(ihisto), "QCD","F");
    }
    
    else if(controlRegion == "zmm"){
      leg->AddEntry(datahist.at(ihisto),"Data","PLE");
      leg->AddEntry(vllhist.at(ihisto), "Z#rightarrow #mu#mu","F");
      leg->AddEntry(vlhist.at(ihisto),  "W #rightarrow #mu#nu","F");
      leg->AddEntry(tophist.at(ihisto), "Top","F");
      leg->AddEntry(dbhist.at(ihisto),  "Di-Boson","F");
      leg->AddEntry(gamhist.at(ihisto), "#gamma+jets","F");
      leg->AddEntry(qcdhist.at(ihisto), "QCD","F");
    }
    
    else if(controlRegion == "zee"){
      leg->AddEntry(datahist.at(ihisto),"Data","PLE");
      leg->AddEntry(vllhist.at(ihisto), "Z #rightarrow ee","F");
      leg->AddEntry(vlhist.at(ihisto),  "W #rightarrow e #nu","F");
      leg->AddEntry(tophist.at(ihisto), "Top","F");
      leg->AddEntry(dbhist.at(ihisto),  "Di-Boson","F");
      leg->AddEntry(gamhist.at(ihisto), "#gamma+jets","F");
      leg->AddEntry(qcdhist.at(ihisto), "QCD","F");
    }
    
    else if(controlRegion == "wmn"){
      leg->AddEntry(datahist.at(ihisto), "Data","PLE");
      leg->AddEntry(vlhist.at(ihisto),   "W #rightarrow #mu#nu","F");
      leg->AddEntry(vllhist.at(ihisto),  "Z #rightarrow #mu#mu","F");
      leg->AddEntry(tophist.at(ihisto),  "Top","F");
      leg->AddEntry(dbhist.at(ihisto),   "Di-Boson","F");
      leg->AddEntry(gamhist.at(ihisto),  "#gamma+jets","F");
      leg->AddEntry(qcdhist.at(ihisto),  "QCD","F");
    }
    
    else if(controlRegion == "wen"){
      leg->AddEntry(datahist.at(ihisto), "Data","PLE");
      leg->AddEntry(vlhist.at(ihisto), "W #rightarrow e#nu","F");
      leg->AddEntry(vllhist.at(ihisto),"Z #rightarrow ee","F");
      leg->AddEntry(tophist.at(ihisto),"Top","F");
      leg->AddEntry(dbhist.at(ihisto), "Di-Boson","F");
      leg->AddEntry(gamhist.at(ihisto),"#gamma+jets","F");
      leg->AddEntry(qcdhist.at(ihisto),"QCD","F");
    }
    
    else if(controlRegion == "topmu" and plotResonant){
      leg->AddEntry(datahist.at(ihisto), "Data","PLE");
      leg->AddEntry(tophist_matched.at(ihisto), "Top Resonant","F");
      leg->AddEntry(tophist_unmatched.at(ihisto), "Top non Resonant","F");
      leg->AddEntry(vlhist.at(ihisto), "W #rightarrow #mu#nu","F");
      leg->AddEntry(vllhist.at(ihisto),"Z #rightarrow #mu#mu","F");
      leg->AddEntry(dbhist.at(ihisto), "Di-Boson","F");
      leg->AddEntry(gamhist.at(ihisto),"#gamma+jets","F");
      leg->AddEntry(qcdhist.at(ihisto),"QCD","F");
      }
    
    else if(controlRegion == "topmu" and not plotResonant){
      leg->AddEntry(datahist.at(ihisto),"Data","PLE");
      leg->AddEntry(tophist.at(ihisto),"Top","F");
      leg->AddEntry(vlhist.at(ihisto), "W #rightarrow #mu#nu","F");
      leg->AddEntry(vllhist.at(ihisto),"Z #rightarrow #mu#mu","F");
      leg->AddEntry(dbhist.at(ihisto), "Di-Boson","F");
      leg->AddEntry(gamhist.at(ihisto),"#gamma+jets","F");
      leg->AddEntry(qcdhist.at(ihisto),"QCD","F");
    }
    
    else if(controlRegion == "topel" and not plotResonant){
      leg->AddEntry(datahist.at(ihisto),"Data","PLE");
      leg->AddEntry(tophist.at(ihisto), "Top","F");
      leg->AddEntry(vlhist.at(ihisto),  "W #rightarrow e#nu","F");
      leg->AddEntry(vllhist.at(ihisto), "Z #rightarrow e#mu","F");
      leg->AddEntry(dbhist.at(ihisto),  "Di-Boson","F");
      leg->AddEntry(gamhist.at(ihisto), "#gamma+jets","F");
      leg->AddEntry(qcdhist.at(ihisto), "QCD","F");
    }
    
    else if(controlRegion == "topel" and plotResonant){
      leg->AddEntry(datahist.at(ihisto), "Data","PLE");
      leg->AddEntry(tophist_matched.at(ihisto), "Top Resonant","F");
      leg->AddEntry(tophist_unmatched.at(ihisto), "Top non Resonant","F");
      leg->AddEntry(vlhist.at(ihisto), "W #rightarrow e#nu","F");
      leg->AddEntry(vllhist.at(ihisto),"Z #rightarrow e#mu","F");
      leg->AddEntry(dbhist.at(ihisto), "Di-Boson","F");
      leg->AddEntry(gamhist.at(ihisto),"#gamma+jets","F");
      leg->AddEntry(qcdhist.at(ihisto),"QCD","F");
    }
    
    else if(controlRegion == "SR"){
      leg->SetNColumns(2);
      leg->AddEntry(datahist.at(ihisto),"Data","PLE");
      leg->AddEntry(vnnhist.at(ihisto), "Z(#nu#nu)","F");
      leg->AddEntry(vlhist.at(ihisto),  "W(l#nu)", "F");
      leg->AddEntry(vllhist.at(ihisto), "Z(ll)", "F");
      leg->AddEntry(tophist.at(ihisto), "Top", "F");
      leg->AddEntry(dbhist.at(ihisto),  "Dibosons", "F");
      leg->AddEntry(gamhist.at(ihisto), "#gamma+jets","F");
      leg->AddEntry(qcdhist.at(ihisto), "QCD", "F");
      if( not isHiggsInvisible){
	TString mass = TString::Format("%.1f TeV",stof(mediatorMass)/1000); 
	leg->AddEntry(monoJhist.at(ihisto), ("Mono-J M_{Med} = "+string(mass)).c_str(), "L");
	leg->AddEntry(monoWhist.at(ihisto), ("Mono-W M_{Med} = "+string(mass)+" #times "+to_string(signalScale)).c_str(), "L");
	leg->AddEntry(monoZhist.at(ihisto), ("Mono-Z M_{Med} = "+string(mass)+" #times "+to_string(signalScale)).c_str(), "L");
      }
      else{
	leg->AddEntry(ggHhist.at(ihisto),"ggH(m_{H}=125 GeV)", "L");
	leg->AddEntry(vbfHhist.at(ihisto),"vbfH(m_{H}=125 GeV)", "L");
	leg->AddEntry(wHhist.at(ihisto),"wH(m_{H}=125 GeV)", "L");
	leg->AddEntry(zHhist.at(ihisto),"zH(m_{H}=125 GeV)", "L");
      }
    }  
    
    leg->Draw("SAME");
    
    pad1->RedrawAxis("sameaxis");
    if(isLog) pad1->SetLogy();
    
    // make data/MC ratio plot
    canvas->cd();
    pad2->SetTopMargin(0.08);
    pad2->SetRightMargin(0.06);
    pad2->SetLeftMargin(0.12);
    pad2->SetBottomMargin(0.35);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();
    
    TH1* frame2 = NULL;
    if(category <= 1)
      frame2 =  pad2->DrawFrame(datahist.at(ihisto)->GetBinLowEdge(1), 0.25, datahist.at(ihisto)->GetBinLowEdge(datahist.at(ihisto)->GetNbinsX()+1), 1.75, "");
    else if(category > 1)
      frame2 =  pad2->DrawFrame(datahist.at(ihisto)->GetBinLowEdge(1), 0.25, datahist.at(ihisto)->GetBinLowEdge(datahist.at(ihisto)->GetNbinsX()+1), 1.75, "");
    
    frame2->GetXaxis()->SetTitle(observableLatex.c_str());
    frame2->GetYaxis()->SetTitle("Data/Pred.");
    frame2->GetYaxis()->CenterTitle();
    frame2->GetXaxis()->SetLabelSize(0.11);
    frame2->GetYaxis()->SetLabelSize(0.10);
    frame2->GetXaxis()->SetTitleSize(0.135);
    frame2->GetYaxis()->SetTitleOffset(0.4);
    frame2->GetYaxis()->SetTitleSize(0.12);
    frame2->GetYaxis()->SetNdivisions(5);
    frame2->GetXaxis()->SetNdivisions(510);
    frame2->Draw();
    
    TH1* nhist = (TH1*) datahist.at(ihisto)->Clone("datahist_tot");
    TH1* unhist = (TH1*) datahist.at(ihisto)->Clone("unhist");
    TH1* dhist = (TH1*) stack->GetStack()->At(stack->GetNhists()-1)->Clone("mchist_tot");
    TH1* dhist_p = (TH1*) stack->GetStack()->At(stack->GetNhists()-1)->Clone("mchist_tot_p");
    
    nhist->SetStats(kFALSE);
    nhist->SetLineColor(kBlack);
    nhist->SetMarkerColor(kBlack);
    nhist->SetMarkerSize(1.2);
    
    // set to zero for plotting reasons of error bar and error band
    for (int i = 1; i <= dhist->GetNbinsX(); i++) dhist->SetBinError(i, 0);
    
    nhist->Divide(dhist);
    dhist_p->Divide(dhist);
      
    dhist_p->SetLineColor(0);
    dhist_p->SetMarkerColor(0);
    dhist_p->SetMarkerSize(0);
    dhist_p->SetFillColor(kGray);
    
    for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinContent(i, 1);
    for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinError(i, 0);
    unhist->SetMarkerSize(0);
    unhist->SetLineColor(kBlack);
    unhist->SetLineStyle(2);
    unhist->SetFillColor(0);
    
    nhist->Draw("PE1 SAME");
    dhist_p->Draw("E2 SAME");
    unhist->Draw("SAME");
    nhist->Draw("PE SAME");
    
    pad2->RedrawAxis("sameaxis");

    TLatex ttext;
    ttext.SetNDC();
    ttext.SetTextFont(42);
    ttext.SetTextAlign(31);
    ttext.SetTextSize(0.04);

    pad1->cd();
    if(ihisto < bin.size()-2)
      ttext.DrawLatex(0.45,0.75,Form("%d <= %s < %d ",int(bin.at(ihisto)),text.second.c_str(),int(bin.at(ihisto+1))));
    else
      ttext.DrawLatex(0.45,0.75,Form("%s >= %d ",text.second.c_str(),int(bin.at(ihisto))));
    
    canvas->SaveAs(Form("%s_%s_bin_%d.png",observable.c_str(),controlRegion.c_str(),int(ihisto)));
    canvas->SaveAs(Form("%s_%s_bin_%d.pdf",observable.c_str(),controlRegion.c_str(),int(ihisto)));
  }
}

