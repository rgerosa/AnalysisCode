#include "../CMS_lumi.h"
#include "../makeTemplates/makehist.h"

static float minTau2Tau1 = 0.1;
static bool saveTextYields = false;

void makeControlPlots(string templateFileName, 
		      Category category, 
		      string observable, 
		      string observableLatex, 
		      string controlRegion, 
		      bool blind, 
		      bool isLog,
		      bool plotResonant   = false,
		      bool isHiggsInvisible = false,
		      bool addSBPlots     = false,
		      bool addShapePlots  = false,
		      string interaction  = "Vector",
		      string mediatorMass = "1000",
		      string DMMass       = "50",
		      int signalScale     = 1) {

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(0);

  initializeBinning();

  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
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


  TFile* inputFile = new TFile(templateFileName.c_str());

  TH1* datahist = NULL;
  TH1* qcdhist  = NULL;
  TH1* vllhist  = NULL;
  TH1* ewkwhist  = NULL;
  TH1* ewkzhist  = NULL;
  TH1* vnnhist  = NULL;
  TH1* vlhist   = NULL;
  TH1* dbhist   = NULL;
  TH1* tophist  = NULL;
  TH1* tophist_matched   = NULL;
  TH1* tophist_unmatched = NULL;
  TH1* gamhist  = NULL;
  TH1* vghist   = NULL;

  TH1* monoJhist  = NULL;
  TH1* monoWhist  = NULL;
  TH1* monoZhist  = NULL;

  TH1* ggHhist  = NULL;
  TH1* vbfHhist = NULL;
  TH1* wHhist   = NULL;
  TH1* zHhist   = NULL;

  // take the templates
  if(controlRegion == "gam"){
    datahist = (TH1*)inputFile->FindObjectAny(("datahistgam_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghistgam_"+observable).c_str());
    gamhist  = (TH1*)inputFile->FindObjectAny(("gbkghistgam_"+observable).c_str());
    vghist   = (TH1*)inputFile->FindObjectAny(("vgbkghistgam_"+observable).c_str());
    vlhist   = (TH1*)inputFile->FindObjectAny(("vlbkghistgam_"+observable).c_str());
  }
  else if(controlRegion == "zmm"){
    datahist = (TH1*)inputFile->FindObjectAny(("datahistzmm_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghistzmm_"+observable).c_str());
    tophist  = (TH1*)inputFile->FindObjectAny(("tbkghistzmm_"+observable).c_str());
    vlhist   = (TH1*)inputFile->FindObjectAny(("vlbkghistzmm_"+observable).c_str());
    vllhist  = (TH1*)inputFile->FindObjectAny(("vllbkghistzmm_"+observable).c_str());
    dbhist   = (TH1*)inputFile->FindObjectAny(("dbkghistzmm_"+observable).c_str());  
    gamhist  = (TH1*)inputFile->FindObjectAny(("gbkghistzmm_"+observable).c_str());
    ewkwhist = (TH1*)inputFile->FindObjectAny(("ewkwbkghistzmm_"+observable).c_str());
    ewkzhist = (TH1*)inputFile->FindObjectAny(("ewkzbkghistzmm_"+observable).c_str());
  }
  else if(controlRegion == "zee"){
    datahist = (TH1*)inputFile->FindObjectAny(("datahistzee_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghistzee_"+observable).c_str());
    tophist  = (TH1*)inputFile->FindObjectAny(("tbkghistzee_"+observable).c_str());
    vlhist   = (TH1*)inputFile->FindObjectAny(("vlbkghistzee_"+observable).c_str());
    vllhist  = (TH1*)inputFile->FindObjectAny(("vllbkghistzee_"+observable).c_str());
    dbhist   = (TH1*)inputFile->FindObjectAny(("dbkghistzee_"+observable).c_str());  
    gamhist  = (TH1*)inputFile->FindObjectAny(("gbkghistzee_"+observable).c_str());
    ewkwhist = (TH1*)inputFile->FindObjectAny(("ewkwbkghistzee_"+observable).c_str());
    ewkzhist = (TH1*)inputFile->FindObjectAny(("ewkzbkghistzee_"+observable).c_str());
  }
  else if(controlRegion == "wmn"){
    datahist = (TH1*)inputFile->FindObjectAny(("datahistwmn_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghistwmn_"+observable).c_str());
    tophist  = (TH1*)inputFile->FindObjectAny(("tbkghistwmn_"+observable).c_str());
    vlhist   = (TH1*)inputFile->FindObjectAny(("vlbkghistwmn_"+observable).c_str());
    vllhist  = (TH1*)inputFile->FindObjectAny(("vllbkghistwmn_"+observable).c_str());
    dbhist   = (TH1*)inputFile->FindObjectAny(("dbkghistwmn_"+observable).c_str());  
    gamhist  = (TH1*)inputFile->FindObjectAny(("gbkghistwmn_"+observable).c_str());
    ewkwhist  = (TH1*)inputFile->FindObjectAny(("ewkwbkghistwmn_"+observable).c_str());
    ewkzhist  = (TH1*)inputFile->FindObjectAny(("ewkzbkghistwmn_"+observable).c_str());
  }
  else if(controlRegion == "wen"){
    datahist = (TH1*)inputFile->FindObjectAny(("datahistwen_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghistwen_"+observable).c_str());
    tophist  = (TH1*)inputFile->FindObjectAny(("tbkghistwen_"+observable).c_str());
    vlhist   = (TH1*)inputFile->FindObjectAny(("vlbkghistwen_"+observable).c_str());
    vllhist  = (TH1*)inputFile->FindObjectAny(("vllbkghistwen_"+observable).c_str());
    dbhist   = (TH1*)inputFile->FindObjectAny(("dbkghistwen_"+observable).c_str());  
    gamhist  = (TH1*)inputFile->FindObjectAny(("gbkghistwen_"+observable).c_str());
    ewkwhist  = (TH1*)inputFile->FindObjectAny(("ewkwbkghistwen_"+observable).c_str());
    ewkzhist  = (TH1*)inputFile->FindObjectAny(("ewkzbkghistwen_"+observable).c_str());
  }
  else if(controlRegion == "topmu" and plotResonant){
    datahist = (TH1*)inputFile->FindObjectAny(("datahisttopmu_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghisttopmu_"+observable).c_str());
    tophist_matched   = (TH1*)inputFile->FindObjectAny(("tbkghist_matchedtopmu_"+observable).c_str());
    tophist_unmatched = (TH1*)inputFile->FindObjectAny(("tbkghist_unmatchedtopmu_"+observable).c_str());
    vlhist   = (TH1*)inputFile->FindObjectAny(("vlbkghisttopmu_"+observable).c_str());
    vllhist  = (TH1*)inputFile->FindObjectAny(("vllbkghisttopmu_"+observable).c_str());
    dbhist   = (TH1*)inputFile->FindObjectAny(("dbkghisttopmu_"+observable).c_str());  
    gamhist  = (TH1*)inputFile->FindObjectAny(("gbkghisttopmu_"+observable).c_str());
  }
  else if(controlRegion == "topmu" and not plotResonant){
    datahist = (TH1*)inputFile->FindObjectAny(("datahisttopmu_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghisttopmu_"+observable).c_str());
    tophist  = (TH1*)inputFile->FindObjectAny(("tbkghisttopmu_"+observable).c_str());
    vlhist   = (TH1*)inputFile->FindObjectAny(("vlbkghisttopmu_"+observable).c_str());
    vllhist  = (TH1*)inputFile->FindObjectAny(("vllbkghisttopmu_"+observable).c_str());
    dbhist   = (TH1*)inputFile->FindObjectAny(("dbkghisttopmu_"+observable).c_str());  
    gamhist  = (TH1*)inputFile->FindObjectAny(("gbkghisttopmu_"+observable).c_str());
  }

  else if(controlRegion == "topel" and not plotResonant){
    datahist = (TH1*)inputFile->FindObjectAny(("datahisttopel_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghisttopel_"+observable).c_str());
    tophist  = (TH1*)inputFile->FindObjectAny(("tbkghisttopel_"+observable).c_str());
    vlhist   = (TH1*)inputFile->FindObjectAny(("vlbkghisttopel_"+observable).c_str());
    vllhist  = (TH1*)inputFile->FindObjectAny(("vllbkghisttopel_"+observable).c_str());
    dbhist   = (TH1*)inputFile->FindObjectAny(("dbkghisttopel_"+observable).c_str());  
    gamhist  = (TH1*)inputFile->FindObjectAny(("gbkghisttopel_"+observable).c_str());
  }

  else if(controlRegion == "topel" and plotResonant){
    datahist = (TH1*)inputFile->FindObjectAny(("datahisttopel_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghisttopel_"+observable).c_str());
    tophist_matched   = (TH1*)inputFile->FindObjectAny(("tbkghist_matchedtopel_"+observable).c_str());
    tophist_unmatched = (TH1*)inputFile->FindObjectAny(("tbkghist_unmatchedtopel_"+observable).c_str());
    vlhist   = (TH1*)inputFile->FindObjectAny(("vlbkghisttopel_"+observable).c_str());
    vllhist  = (TH1*)inputFile->FindObjectAny(("vllbkghisttopel_"+observable).c_str());
    dbhist   = (TH1*)inputFile->FindObjectAny(("dbkghisttopel_"+observable).c_str());  
    gamhist  = (TH1*)inputFile->FindObjectAny(("gbkghisttopel_"+observable).c_str());
  }
  else if(controlRegion == "qcd"){    
    datahist = (TH1*)inputFile->FindObjectAny(("datahistqcd_"+observable).c_str());
    qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghistqcd_"+observable).c_str());
    tophist  = (TH1*)inputFile->FindObjectAny(("tbkghistqcd_"+observable).c_str());
    dbhist   = (TH1*)inputFile->FindObjectAny(("dbkghistqcd_"+observable).c_str());
    vllhist  = (TH1*)inputFile->FindObjectAny(("vllbkghistqcd_"+observable).c_str());
    vlhist   = (TH1*)inputFile->FindObjectAny(("vlbkghistqcd_"+observable).c_str());
    gamhist  = (TH1*)inputFile->FindObjectAny(("gbkghistqcd_"+observable).c_str());
    ewkwhist  = (TH1*)inputFile->FindObjectAny(("ewkbkgwhistqcd_"+observable).c_str());
    ewkzhist  = (TH1*)inputFile->FindObjectAny(("ewkbkgzhistqcd_"+observable).c_str());
    vnnhist   = (TH1*)inputFile->FindObjectAny(("vnnbkghistqcd_"+observable).c_str());    
  }
  else if(controlRegion == "SR"){

    datahist = (TH1*)inputFile->FindObjectAny(("datahist_"+observable).c_str());

    if(category == Category::monojet and observable == "met")
      qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghistDD_"+observable).c_str());
    else if(category == Category::monoV and observable == "met")
      qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghistDD_"+observable).c_str());
    if(qcdhist == 0 or qcdhist == NULL)
      qcdhist  = (TH1*)inputFile->FindObjectAny(("qbkghist_"+observable).c_str());
    
    tophist  = (TH1*)inputFile->FindObjectAny(("tbkghist_"+observable).c_str());
    vlhist   = (TH1*)inputFile->FindObjectAny(("wjethist_"+observable).c_str());
    vllhist  = (TH1*)inputFile->FindObjectAny(("zjethist_"+observable).c_str());
    vnnhist  = (TH1*)inputFile->FindObjectAny(("zinvhist_"+observable).c_str());
    dbhist   = (TH1*)inputFile->FindObjectAny(("dbkghist_"+observable).c_str());  
    gamhist  = (TH1*)inputFile->FindObjectAny(("gbkghist_"+observable).c_str());
    ewkwhist  = (TH1*)inputFile->FindObjectAny(("ewkbkgwhist_"+observable).c_str());
    ewkzhist  = (TH1*)inputFile->FindObjectAny(("ewkbkgzhist_"+observable).c_str());
    
    if(not isHiggsInvisible){
      monoJhist = (TH1*)inputFile->FindObjectAny(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      monoWhist = (TH1*)inputFile->FindObjectAny(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());
      monoZhist = (TH1*)inputFile->FindObjectAny(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str());    
    }
    else{
      ggHhist  = (TH1*)inputFile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+observable).c_str());
      vbfHhist = (TH1*)inputFile->FindObjectAny(("vbfHhist_"+mediatorMass+"_"+observable).c_str());
      wHhist   = (TH1*)inputFile->FindObjectAny(("wHhist_"+mediatorMass+"_"+observable).c_str());    
      zHhist   = (TH1*)inputFile->FindObjectAny(("zHhist_"+mediatorMass+"_"+observable).c_str());    
    }
  }


  //  if((controlRegion == "topmu" or controlRegion == "topel") and (category == Category::monoV or category == Category::boosted)){
  //    tophist->Scale(0.92);
  //    vlhist->Scale(0.95);
  //  }

  if(saveTextYields){

    // write yields in a output in a text file 
    ofstream outputfile;
    outputfile.open(("preFit_"+controlRegion+".txt").c_str());

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
    stringstream DataRate;
    DataRate << "Process: Data";
    
    for(int iBin = 0; iBin < qcdhist->GetNbinsX(); iBin++){
      QCDRate << "   ";
      QCDRate << qcdhist->GetBinContent(iBin+1) << " \\pm "<<qcdhist->GetBinError(iBin+1);
    }
    
    for(int iBin = 0; iBin < gamhist->GetNbinsX(); iBin++){
      GJetsRate << "   ";
      GJetsRate << gamhist->GetBinContent(iBin+1)<< " \\pm "<<gamhist->GetBinError(iBin+1);
    }
    
    for(int iBin = 0; iBin < dbhist->GetNbinsX(); iBin++){
      DiBosonRate << "   ";
      DiBosonRate << dbhist->GetBinContent(iBin+1)<< " \\pm "<<dbhist->GetBinError(iBin+1);;
    }
    
    for(int iBin = 0; iBin < tophist->GetNbinsX(); iBin++){
      TopRate << "   ";
      TopRate << tophist->GetBinContent(iBin+1) << " \\pm "<<tophist->GetBinError(iBin+1);
    }
    
    for(int iBin = 0; iBin < vllhist->GetNbinsX(); iBin++){
      ZJetsRate << "   ";
      ZJetsRate << vllhist->GetBinContent(iBin+1)<< " \\pm "<<vllhist->GetBinError(iBin+1);
    }
    
    for(int iBin = 0; iBin < vlhist->GetNbinsX(); iBin++){
      WJetsRate << "   ";
      WJetsRate << vlhist->GetBinContent(iBin+1)<< " \\pm "<<vlhist->GetBinError(iBin+1);
    }
    
    for(int iBin = 0; iBin < vnnhist->GetNbinsX(); iBin++){
      ZnunuRate << "   ";
      ZnunuRate << vnnhist->GetBinContent(iBin+1) << " \\pm "<<vnnhist->GetBinError(iBin+1);
    }
    
    for(int iBin = 0; iBin < datahist->GetNbinsX(); iBin++){
      DataRate << "   ";
      DataRate << datahist->GetBinContent(iBin+1) << " \\pm "<<datahist->GetBinError(iBin+1);
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
    outputfile<<DataRate.str()<<endl;
    outputfile<<"######################"<<endl;

    outputfile.close();
  }


  cout<<"#########################"<<endl;
  cout<<"Total Yields for process "<<endl;
  cout<<"#########################"<<endl;
  if(qcdhist)
    cout<<"QCD Background        :"<<qcdhist->Integral()<<endl;
  if(dbhist)
    cout<<"Diboson Background    :"<<dbhist->Integral()<<endl;
  if(tophist)
    cout<<"Top Background        :"<<tophist->Integral()<<endl;
  if(gamhist)
    cout<<"#gamma+jet Background :"<<gamhist->Integral()<<endl;
  if(vllhist)
    cout<<"Z+jets Background     :"<<vllhist->Integral()<<endl;
  if(vlhist)
    cout<<"W+jets Background     :"<<vlhist->Integral()<<endl;
  if(ewkwhist)
    cout<<"EWK-W Background      :"<<ewkwhist->Integral()<<endl;
  if(ewkwhist)
    cout<<"EWK-Z Background      :"<<ewkzhist->Integral()<<endl;
  if(vnnhist)
    cout<<"Zvv   Background      :"<<vnnhist->Integral()<<endl;
  cout<<"-------------------------------------------"<<endl;
  if(datahist)
    cout<<"Data integral         :"<<datahist->Integral()<<endl;

  //SCALE BIN WIDTH
  if(TString(observableLatex).Contains("GeV")){

    if(controlRegion == "SR" and not TString(qcdhist->GetName()).Contains("qbkghistDD"))
      qcdhist->Scale(2.);
    
    if(datahist)
      datahist->Scale(1.0,"width");
    if(qcdhist)
      qcdhist->Scale(1.0,"width");
    if(tophist)
      tophist->Scale(1.0,"width");
    if(tophist_matched)
      tophist_matched->Scale(1.0,"width");
    if(tophist_unmatched)
      tophist_unmatched->Scale(1.0,"width");
    if(vlhist)
      vlhist->Scale(1.0,"width");
    if(vllhist)
      vllhist->Scale(1.0,"width");
    if(vnnhist)
      vnnhist->Scale(1.0,"width");
    if(dbhist)
      dbhist->Scale(1.0,"width");
    if(gamhist)
      gamhist->Scale(1.0,"width");
    if(ewkwhist)
      ewkwhist->Scale(1.0,"width");
    if(ewkzhist)
      ewkzhist->Scale(1.0,"width");
    if(vghist)
      vghist->Scale(1.0,"width");

    if(monoJhist){
      monoJhist->Scale(1.0,"width");
      monoJhist->Scale(signalScale);
    }
    if(monoWhist){
      monoWhist->Scale(1.0,"width");
      monoWhist->Scale(signalScale);
    }
    if(monoZhist){
      monoZhist->Scale(1.0,"width");
      monoZhist->Scale(signalScale);
    }

    if(ggHhist){
      ggHhist->Scale(1.0,"width");
      ggHhist->Scale(signalScale);
    }
    if(vbfHhist){
      vbfHhist->Scale(1.0,"width");
      vbfHhist->Scale(signalScale);
    }
    if(wHhist){
      wHhist->Scale(1.0,"width");
      wHhist->Scale(signalScale);
    }
    if(zHhist){
      zHhist->Scale(1.0,"width");
      zHhist->Scale(signalScale);
    }
  }
  else{
    if(controlRegion == "SR" and not TString(qcdhist->GetName()).Contains("qbkghistDD"))
      qcdhist->Scale(2.);

    if(ggHhist)
      ggHhist->Scale(signalScale);
    if(vbfHhist)
      vbfHhist->Scale(signalScale);
    if(wHhist)
      wHhist->Scale(signalScale);
    if(zHhist)
      zHhist->Scale(signalScale);
    if(monoJhist)
      monoJhist->Scale(signalScale);
    if(monoWhist)
      monoWhist->Scale(signalScale);
    if(monoZhist)
      monoZhist->Scale(signalScale);
  }

  // BLIND OPTION
  if(blind and controlRegion == "SR"){

    for (int i = 0; i <= datahist->GetNbinsX()+1; i++) {

      double yield = 0.0;
      yield += qcdhist->GetBinContent(i);
      yield += gamhist->GetBinContent(i);
      yield += tophist->GetBinContent(i);
      yield += dbhist->GetBinContent(i);
      if(category == Category::VBF){
	yield += ewkzhist->GetBinContent(i);
	yield += ewkwhist->GetBinContent(i);
      }
      yield += vllhist->GetBinContent(i);
      yield += vlhist->GetBinContent(i);
      yield += vnnhist->GetBinContent(i);
      datahist->SetBinContent(i,yield);
      datahist->SetBinError(i,0.);
    }
  }

  // set colors
  if(datahist){
    datahist->SetLineColor(kBlack);
    datahist->SetMarkerColor(kBlack);
    datahist->SetMarkerStyle(20);
    datahist->SetMarkerSize(1.2);
  }

  if(vnnhist) {
      vnnhist->SetFillColor(TColor::GetColor("#4D975D"));
      vnnhist->SetLineColor(kBlack);
  }
  if(vllhist){    
      vllhist->SetFillColor(TColor::GetColor("#9A9EAB"));
      vllhist->SetLineColor(kBlack);      
  }
  if(vlhist){
      vlhist->SetFillColor(TColor::GetColor("#FAAF08"));
      vlhist->SetLineColor(kBlack);
  }
  if(tophist){
      tophist->SetFillColor(TColor::GetColor("#CF3721"));
      tophist->SetLineColor(kBlack);
  }
  if(vghist){
    vghist->SetFillColor(kGreen+1);
    vghist->SetLineColor(kBlack);
  }
  if(tophist_matched){
    tophist_matched->SetFillColor(kGreen+1);
    tophist_matched->SetLineColor(kBlack);
  }
  if(tophist_unmatched){
    tophist_unmatched->SetFillColor(kBlue);
    tophist_unmatched->SetLineColor(kBlack);
  }
  if(dbhist){
    dbhist->SetFillColor(TColor::GetColor("#4897D8"));
  }
  if(qcdhist) {
    qcdhist->SetFillColor(TColor::GetColor("#F1F1F2"));
    qcdhist->SetLineColor(kBlack);
  }
  if(gamhist){
    if(controlRegion == "SR"){
      gamhist->SetFillColor(TColor::GetColor("#9A9EAB"));
      gamhist->SetLineColor(TColor::GetColor("#9A9EAB"));
    }
    else{
      gamhist->SetFillColor(TColor::GetColor("#db4dff"));
      gamhist->SetLineColor(kBlack);
    } 
  }

  if(ewkwhist){
    ewkwhist->SetFillColor(kCyan+1);
    ewkwhist->SetLineColor(kBlack);
  }

  if(ewkzhist){
    ewkzhist->SetLineColor(kCyan+1);
    ewkzhist->SetFillColor(kBlack);
  }

  if(monoJhist){
    monoJhist->SetFillColor(0);
    monoJhist->SetFillStyle(0);
    monoJhist->SetLineColor(kBlack);
    monoJhist->SetLineWidth(3);
    monoJhist->SetLineStyle(7);
  }

  if(monoWhist){
    monoWhist->SetFillColor(0);
    monoWhist->SetFillStyle(0);
    monoWhist->SetLineColor(kBlue);
    monoWhist->SetLineWidth(3);
  }

  if(monoZhist){
    monoZhist->SetFillColor(0);
    monoZhist->SetFillStyle(0);
    monoZhist->SetLineColor(TColor::GetColor("#A2C523"));
    monoZhist->SetLineWidth(3);
  }

  if(ggHhist){
    ggHhist->SetFillColor(0);
    ggHhist->SetFillStyle(0);
    ggHhist->SetLineColor(kBlack);
    ggHhist->SetLineWidth(3);
    ggHhist->SetLineStyle(7);
  }

  if(vbfHhist){
    vbfHhist->SetFillColor(0);
    vbfHhist->SetFillStyle(0);
    vbfHhist->SetLineColor(kBlue);
    vbfHhist->SetLineWidth(3);
  }

  if(wHhist){
    wHhist->SetFillColor(0);
    wHhist->SetFillStyle(0);
    wHhist->SetLineColor(TColor::GetColor("#A2C523"));
    wHhist->SetLineWidth(3);
  }

  if(zHhist){
    zHhist->SetFillColor(0);
    zHhist->SetFillStyle(0);
    zHhist->SetLineColor(TColor::GetColor("#A2C523"));
    zHhist->SetLineWidth(3);
  }

  if(wHhist and zHhist) // add them together in vH
    wHhist->Add(zHhist);

  THStack* stack = new THStack("stack", "stack");
  if(controlRegion == "gam"){
    stack->Add(vghist);
    stack->Add(vlhist);
    stack->Add(qcdhist);
    stack->Add(gamhist);
  }
  else if(controlRegion == "zmm" or controlRegion == "zee"){
    stack->Add(qcdhist);
    stack->Add(gamhist);
    stack->Add(vlhist);
    stack->Add(tophist);
    stack->Add(dbhist);
    if(category == Category::VBF){
      ewkwhist->Add(ewkzhist);
      stack->Add(ewkwhist);
    }
    stack->Add(vllhist);
  }
  else if(controlRegion == "wmn" or controlRegion == "wen"){
    stack->Add(qcdhist);
    stack->Add(gamhist);
    stack->Add(vllhist);
    stack->Add(tophist);
    stack->Add(dbhist);
    if(category == Category::VBF){
      ewkwhist->Add(ewkzhist);
      stack->Add(ewkwhist);
    }
    stack->Add(vlhist);
  }
  else if((controlRegion == "topmu" or controlRegion == "topel") and not plotResonant){
    stack->Add(qcdhist);
    stack->Add(gamhist);
    stack->Add(vllhist);
    if(category == Category::VBF){
      ewkwhist->Add(ewkzhist);
      stack->Add(ewkwhist);
    }
    stack->Add(dbhist);
    stack->Add(vlhist);
    stack->Add(tophist);    
  }
  else if((controlRegion == "topmu" or controlRegion == "topel") and plotResonant){
    stack->Add(qcdhist);
    stack->Add(gamhist);
    stack->Add(vllhist);
    stack->Add(dbhist);
    if(category == Category::VBF){
      ewkwhist->Add(ewkzhist);
      stack->Add(ewkwhist);
    }
    stack->Add(vlhist);
    stack->Add(tophist_unmatched);    
    stack->Add(tophist_matched);    
  }
  else if(controlRegion == "qcd"){
    if(not isnan(float(gamhist->Integral())))
      vllhist->Add(gamhist);
    if(not isnan(float(vllhist->Integral())))    
      stack->Add(vllhist);
    if(not isnan(float(dbhist->Integral())))    
      stack->Add(dbhist);
    if(not isnan(float(tophist->Integral())))    
      stack->Add(tophist);
    if(category == Category::VBF){
      ewkwhist->Add(ewkzhist);
      stack->Add(ewkwhist);
    }
    if(not isnan(float(vlhist->Integral())))    
      stack->Add(vlhist);
    if(not isnan(float(vnnhist->Integral())))    
      stack->Add(vnnhist);
    if(not isnan(float(qcdhist->Integral())))    
      stack->Add(qcdhist);
  }
  else if(controlRegion == "SR"){
    stack->Add(qcdhist);
    vllhist->Add(gamhist);
    stack->Add(vllhist);
    stack->Add(tophist);
    stack->Add(dbhist);
    if(category == Category::VBF){
      ewkwhist->Add(ewkzhist);
      stack->Add(ewkwhist);
    }
    stack->Add(vlhist);
    stack->Add(vnnhist);
  }


  TH1* frame = (TH1*) datahist->Clone("frame");
  frame->Reset();  
  vector<double> bins ;
  if(observable == "bosonPt")
    bins = selectBinning("bosonpt",category);
  else
    bins = selectBinning(observable,category);

  float xMin = bins.front();
  if(observable == "tau2tau1")
    xMin = minTau2Tau1;
  if(TString(observable).Contains("btag"))
    xMin = 0.01;
  float xMax = bins.back();

  // set Y-axis range
  if(category == Category::monojet and isLog)
    frame->GetYaxis()->SetRangeUser(1.5e-2,datahist->GetMaximum()*500);  
  else if(category == Category::inclusive and isLog)
    frame->GetYaxis()->SetRangeUser(1.5e-3,datahist->GetMaximum()*500);  
  else if(category == Category::monojet and not isLog)
    frame->GetYaxis()->SetRangeUser(1.5e-3,datahist->GetMaximum()*2.5);  
  else if(category == Category::inclusive and not isLog)
    frame->GetYaxis()->SetRangeUser(1.5e-3,datahist->GetMaximum()*2.5);  
  else if(category == Category::monoV and isLog)
    frame->GetYaxis()->SetRangeUser(1.5e-3,datahist->GetMaximum()*500);  
  else if(category == Category::boosted and isLog)
    frame->GetYaxis()->SetRangeUser(1.5e-3,datahist->GetMaximum()*500);  
  else if(category == Category::VBF and isLog)
    frame->GetYaxis()->SetRangeUser(1e-4,datahist->GetMaximum()*500);  
  else
    frame->GetYaxis()->SetRangeUser(1.5e-3,datahist->GetMaximum()*2.5);  
    
  frame->GetXaxis()->SetTitle(observableLatex.c_str());
  if(TString(observableLatex).Contains("GeV"))
    frame->GetYaxis()->SetTitle("Events / GeV");
  else
    frame->GetYaxis()->SetTitle("Events");

  frame->GetXaxis()->SetTitleSize(0);
  frame->GetXaxis()->SetLabelSize(0);
  frame->GetYaxis()->SetLabelSize(0.045);
  frame->GetYaxis()->SetTitleSize(0.055);
  frame->GetYaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetLabelSize(0.040);
  frame->GetYaxis()->SetTitleSize(0.050);
  if(category == Category::monojet)
    frame->GetXaxis()->SetNdivisions(510);
  else
    frame->GetXaxis()->SetNdivisions(504);
  
  frame->Draw();
  CMS_lumi(canvas,"36.4");
  
  stack ->Draw("HIST SAME");
  datahist->Draw("PE SAME");

  if(controlRegion == "SR" and not isHiggsInvisible){
    if(monoJhist)
      monoJhist->Draw("hist same");
    if(monoWhist)
      monoWhist->Draw("hist same");
    if(monoZhist)
    monoZhist->Draw("hist same");
  }
  else if(controlRegion == "SR" and isHiggsInvisible){
    if(ggHhist)
      ggHhist->Draw("hist same");
    if(vbfHhist)
      vbfHhist->Draw("hist same");
    if(wHhist)
      wHhist->Draw("hist same");
    if(zHhist)
      zHhist->Draw("hist same");
  }

  TLegend* leg = NULL;
  if(controlRegion == "gam")
    leg = new TLegend(0.60, 0.70, 0.92, 0.92);
  else if (observable == "chfrac" or observable == "nhfrac" or observable == "emfrac")
    leg = new TLegend(0.60, 0.35, 0.92, 0.65);  
  else if(controlRegion == "SR" and isLog)
    leg = new TLegend(0.60, 0.55, 0.92, 0.92);  
  else if(controlRegion == "SR" and not isLog)
    leg = new TLegend(0.60, 0.55, 0.92, 0.92);  
  else
    leg = new TLegend(0.60, 0.55, 0.92, 0.92);

  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  if(controlRegion == "gam"){
    leg->AddEntry(datahist, "Data","PLE");
    leg->AddEntry(gamhist, "#gamma+jets","F");
    leg->AddEntry(qcdhist, "QCD","F");
    leg->AddEntry(vlhist, "W+Jets","F");
    leg->AddEntry(vghist, "V#gamma","F");
  }

  else if(controlRegion == "zmm"){
    leg->AddEntry(datahist,"Data","PLE");
    leg->AddEntry(vllhist, "Z#rightarrow #mu#mu","F");
    leg->AddEntry(vlhist,  "W #rightarrow #mu#nu","F");
    if(category == Category::VBF)	  
      leg->AddEntry(ewkwhist,  "EWK W/Z + 2jet","F");
    leg->AddEntry(tophist, "Top","F");
    leg->AddEntry(dbhist,  "Di-Boson","F");
    leg->AddEntry(gamhist, "#gamma+jets","F");
    leg->AddEntry(qcdhist, "QCD","F");
  }

  else if(controlRegion == "zee"){
    leg->AddEntry(datahist,"Data","PLE");
    leg->AddEntry(vllhist, "Z #rightarrow ee","F");
    leg->AddEntry(vlhist,  "W #rightarrow e #nu","F");
    if(category == Category::VBF)	  
      leg->AddEntry(ewkwhist,  "EWK W/Z + 2jet","F");
    leg->AddEntry(tophist, "Top","F");
    leg->AddEntry(dbhist,  "Di-Boson","F");
    leg->AddEntry(gamhist, "#gamma+jets","F");
    leg->AddEntry(qcdhist, "QCD","F");
  }

  else if(controlRegion == "wmn"){
    leg->AddEntry(datahist, "Data","PLE");
    leg->AddEntry(vlhist,   "W #rightarrow #mu#nu","F");
    leg->AddEntry(vllhist,  "Z #rightarrow #mu#mu","F");
    if(category == Category::VBF)	  
      leg->AddEntry(ewkwhist,  "EWK W/Z + 2jet","F");
    leg->AddEntry(tophist,  "Top","F");
    leg->AddEntry(dbhist,   "Di-Boson","F");
    leg->AddEntry(gamhist,  "#gamma+jets","F");
    leg->AddEntry(qcdhist,  "QCD","F");
  }

  else if(controlRegion == "wen"){
    leg->AddEntry(datahist, "Data","PLE");
    leg->AddEntry(vlhist, "W #rightarrow e#nu","F");
    leg->AddEntry(vllhist,"Z #rightarrow ee","F");
    if(category == Category::VBF)	  
      leg->AddEntry(ewkwhist,  "EWK W/Z + 2jet","F");
    leg->AddEntry(tophist,"Top","F");
    leg->AddEntry(dbhist, "Di-Boson","F");
    leg->AddEntry(gamhist,"#gamma+jets","F");
    leg->AddEntry(qcdhist,"QCD","F");
  }

  else if(controlRegion == "topmu" and plotResonant){
    leg->AddEntry(datahist, "Data","PLE");
    leg->AddEntry(tophist_matched, "Top Resonant","F");
    leg->AddEntry(tophist_unmatched, "Top non Resonant","F");
    leg->AddEntry(vlhist, "W #rightarrow #mu#nu","F");
    leg->AddEntry(vllhist,"Z #rightarrow #mu#mu","F");
    leg->AddEntry(dbhist, "Di-Boson","F");
    leg->AddEntry(gamhist,"#gamma+jets","F");
    leg->AddEntry(qcdhist,"QCD","F");
  }

  else if(controlRegion == "topmu" and not plotResonant){
    leg->AddEntry(datahist,"Data","PLE");
    leg->AddEntry(tophist,"Top","F");
    leg->AddEntry(vlhist, "W #rightarrow #mu#nu","F");
    leg->AddEntry(vllhist,"Z #rightarrow #mu#mu","F");
    leg->AddEntry(dbhist, "Di-Boson","F");
    leg->AddEntry(gamhist,"#gamma+jets","F");
    leg->AddEntry(qcdhist,"QCD","F");
  }

  else if(controlRegion == "topel" and not plotResonant){
    leg->AddEntry(datahist,"Data","PLE");
    leg->AddEntry(tophist, "Top","F");
    leg->AddEntry(vlhist,  "W #rightarrow e#nu","F");
    leg->AddEntry(vllhist, "Z #rightarrow e#mu","F");
    leg->AddEntry(dbhist,  "Di-Boson","F");
    leg->AddEntry(gamhist, "#gamma+jets","F");
    leg->AddEntry(qcdhist, "QCD","F");
  }

  else if(controlRegion == "topel" and plotResonant){
    leg->AddEntry(datahist, "Data","PLE");
    leg->AddEntry(tophist_matched, "Top Resonant","F");
    leg->AddEntry(tophist_unmatched, "Top non Resonant","F");
    leg->AddEntry(vlhist, "W #rightarrow e#nu","F");
    leg->AddEntry(vllhist,"Z #rightarrow e#mu","F");
    leg->AddEntry(dbhist, "Di-Boson","F");
    leg->AddEntry(gamhist,"#gamma+jets","F");
    leg->AddEntry(qcdhist,"QCD","F");
  }

  else if(controlRegion == "qcd"){

    leg->AddEntry(datahist, "Data","PLE");
    leg->AddEntry(qcdhist,"QCD","F");
    leg->AddEntry(vnnhist,"Z #rightarrow #nu#nu","F");
    leg->AddEntry(vlhist,"W #rightarrow l#nu","F");
    if(category == Category::VBF)	  
      leg->AddEntry(ewkwhist,  "EWK W/Z + 2jet","F");
    leg->AddEntry(dbhist,  "WW/WZ/ZZ", "F");
    leg->AddEntry(tophist, "Top quark", "F");
    leg->AddEntry(vllhist, "Z #rightarrow ll, #gamma+jets","F");
  }

  else if(controlRegion == "SR"){
    leg->AddEntry(datahist,"Data","PLE");
    leg->AddEntry(vnnhist, "Z #rightarrow #nu#nu","F");
    leg->AddEntry(vlhist,  "W #rightarrow l#nu", "F");
    if(category == Category::VBF)	  
      leg->AddEntry(ewkwhist,  "EWK W/Z + 2jet","F");
    leg->AddEntry(dbhist,  "WW/WZ/ZZ", "F");
    leg->AddEntry(tophist, "Top quark", "F");
    leg->AddEntry(vllhist, "Z #rightarrow ll, #gamma+jets","F");
    leg->AddEntry(qcdhist, "QCD", "F");
    if( not isHiggsInvisible){
      TString mass = TString::Format("%.1f TeV",stof(mediatorMass)/1000); 
      if(monoJhist)
	leg->AddEntry(monoJhist, ("Mono-J M_{Med} = "+string(mass)).c_str(), "L");
      if(monoWhist)
	leg->AddEntry(monoWhist, ("Mono-W M_{Med} = "+string(mass)+" #times "+to_string(signalScale)).c_str(), "L");
      if(monoZhist)
	leg->AddEntry(monoZhist, ("Mono-Z M_{Med} = "+string(mass)+" #times "+to_string(signalScale)).c_str(), "L");
    }
    else{
      leg->AddEntry(ggHhist,"ggH(m_{H}=125 GeV)", "L");
      leg->AddEntry(vbfHhist,"vbfH(m_{H}=125 GeV)", "L");
      leg->AddEntry(wHhist,"vH(m_{H}=125 GeV)", "L");
    }
  }  

  leg->Draw("SAME");
  
  canvas->RedrawAxis("sameaxis");
  if(isLog)
    canvas->SetLogy();

  // make data/MC ratio plot
  canvas->cd();
  pad2->Draw();
  pad2->cd();

  TH1* frame2 = (TH1*) datahist->Clone("frame");
  frame2->Reset();
  if((category == Category::monojet or category == Category::inclusive) and controlRegion != "qcd")
    frame2->GetYaxis()->SetRangeUser(0.5,1.5);
  else if(category == Category::monoV and controlRegion != "qcd")
    frame2->GetYaxis()->SetRangeUser(0.25,1.75);
  else if(category == Category::twojet and controlRegion != "qcd")
    frame2->GetYaxis()->SetRangeUser(0.5,1.5);
  else if(category == Category::VBF and controlRegion != "qcd")
    frame2->GetYaxis()->SetRangeUser(0.5,1.5);
  else if((category == Category::boosted or category == Category::prunedMass or category == Category::tau2tau1) and controlRegion != "qcd")
    frame2->GetYaxis()->SetRangeUser(0.5,1.5);
  else if(controlRegion == "qcd")
    frame2->GetYaxis()->SetRangeUser(0.5,2.5);
  
  if(category == Category::monojet)
    frame2->GetXaxis()->SetNdivisions(510);
  else
    frame2->GetXaxis()->SetNdivisions(510);
  frame2->GetYaxis()->SetNdivisions(5);

  frame2->GetXaxis()->SetTitle(observableLatex.c_str());
  frame2->GetYaxis()->SetTitle("Data/Pred.");
  frame2->GetYaxis()->CenterTitle();
  frame2->GetYaxis()->SetTitleOffset(1.5);
  frame2->GetYaxis()->SetLabelSize(0.035);
  frame2->GetYaxis()->SetTitleSize(0.04);
  frame2->GetXaxis()->SetLabelSize(0.04);
  frame2->GetXaxis()->SetTitleSize(0.05);
  frame2->Draw();


  TH1* nhist = (TH1*) datahist->Clone("datahist_tot");
  TH1* unhist = (TH1*) datahist->Clone("unhist");
  TH1* dhist = (TH1*) stack->GetStack()->At(stack->GetNhists()-1)->Clone("mchist_tot");
  TH1* dhist_p = (TH1*) stack->GetStack()->At(stack->GetNhists()-1)->Clone("mchist_tot_p");

  nhist->SetStats(kFALSE);
  nhist->SetLineColor(kBlack);
  nhist->SetMarkerColor(kBlack);
  nhist->SetMarkerSize(0.8);

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

  canvas->SaveAs((observable+"_"+controlRegion+".png").c_str());
  canvas->SaveAs((observable+"_"+controlRegion+".pdf").c_str());

  if(addSBPlots){

    TH1* totalSignal = NULL;

    if(isHiggsInvisible){
      totalSignal = (TH1*) ggHhist->Clone("totalSignal");
      totalSignal->Add(vbfHhist);
      totalSignal->Add(wHhist);
      totalSignal->Add(zHhist);
    }
    else{
      totalSignal = (TH1*) monoJhist->Clone("monoJhist");
      totalSignal->Add(monoWhist);
      totalSignal->Add(monoZhist);
    }
    
    pad2->Draw();
    pad2->cd();
    frame2->GetYaxis()->SetTitle("(S+B)/B");

    TH1* SoverB_prefit = (TH1*) totalSignal->Clone("SoverB_prefit");
    SoverB_prefit->Add((TH1*) stack->GetStack()->At(stack->GetNhists()-1));
    SoverB_prefit->Divide((TH1*) stack->GetStack()->At(stack->GetNhists()-1));
    TH1* SoverB_prefit_d = (TH1*) SoverB_prefit->Clone("SoverB_prefit_d");
    for(int iBin = 0; iBin < SoverB_prefit_d->GetNbinsX(); iBin++)
      SoverB_prefit_d->SetBinContent(iBin+1,1);
    frame2->GetYaxis()->SetRangeUser(std::max(SoverB_prefit->GetMinimum(),SoverB_prefit_d->GetMinimum())*0.5,std::max(SoverB_prefit->GetMaximum(),SoverB_prefit_d->GetMaximum()) *1.2);
    frame2->Draw();
    SoverB_prefit->Draw("hist same");
    SoverB_prefit_d->SetLineColor(0);
    SoverB_prefit_d->SetMarkerColor(0);
    SoverB_prefit_d->SetMarkerSize(0);
    SoverB_prefit_d->SetFillColor(kGray);
    SoverB_prefit_d->SetFillStyle(1001);
    SoverB_prefit_d->Draw("E2 SAME");
    unhist->Draw("SAME");
    SoverB_prefit->Draw("hist same");
    pad2->RedrawAxis("sameaxis");
    canvas->SaveAs((observable+"_"+controlRegion+"_SoB.png").c_str());
    canvas->SaveAs((observable+"_"+controlRegion+"_SoB.pdf").c_str());   

  }

  if(addShapePlots and controlRegion == "SR"){

    canvas->cd();    
    TH1* totalBkg = (TH1*) stack->GetStack()->At(stack->GetNhists()-1);
    totalBkg->Scale(1./totalBkg->Integral());
    totalBkg->SetLineWidth(2);
    totalBkg->SetLineColor(kBlack);
    totalBkg->SetFillColor(kGray);
    totalBkg->SetFillStyle(3001);
    frame->GetYaxis()->SetTitle("A.U.");
    if(isLog)
      frame->GetYaxis()->SetRangeUser(totalBkg->GetMinimum()/10,totalBkg->GetMaximum()*30);
    else
      frame->GetYaxis()->SetRangeUser(0.,totalBkg->GetMaximum()*1.2);

    frame->Draw();
    CMS_lumi(canvas,"2.3");
    totalBkg->Draw("hist same");
    
    TLegend* leg2 = NULL;
    
    if(isHiggsInvisible){
      ggHhist->SetLineStyle(1);
      ggHhist->Scale(1./ggHhist->Integral());
      ggHhist->Draw("hist same");
      vbfHhist->Scale(1./vbfHhist->Integral());
      vbfHhist->Draw("hist same");
      wHhist->Scale(1./wHhist->Integral());
      wHhist->Draw("hist same");
      zHhist->Scale(1./wHhist->Integral());
      zHhist->Draw("hist same");

      leg2 = new TLegend(0.52, 0.65, 0.88, 0.90);
      leg2->SetFillColor(0);
      leg2->SetFillStyle(0);
      leg2->SetBorderSize(0);
      leg2->AddEntry(ggHhist,"ggH(m_{H}=125 GeV)", "L");
      leg2->AddEntry(vbfHhist,"vbfH(m_{H}=125 GeV)", "L");
      leg2->AddEntry(wHhist,"vH(m_{H}=125 GeV)", "L");
      leg2->AddEntry(totalBkg,"total background","FL");

    }
    else{
      monoJhist->SetLineStyle(1);
      monoJhist->Scale(1./monoJhist->Integral());
      monoJhist->Draw("hist same");
      monoWhist->Scale(1./monoWhist->Integral());
      monoWhist->Draw("hist same");
      monoZhist->Scale(1./monoZhist->Integral());
      monoZhist->Draw("hist same");

      leg2 = new TLegend(0.52, 0.65, 0.88, 0.90);
      leg2->SetFillColor(0);
      leg2->SetFillStyle(0);
      leg2->SetBorderSize(0);
      TString mass = TString::Format("%.1f TeV",stof(mediatorMass)/1000); 
      leg2->AddEntry(monoJhist, ("Mono-J M_{Med} = "+string(mass)).c_str(), "L");
      leg2->AddEntry(monoWhist, ("Mono-W M_{Med} = "+string(mass)+" #times "+to_string(signalScale)).c_str(), "L");
      leg2->AddEntry(monoZhist, ("Mono-Z M_{Med} = "+string(mass)+" #times "+to_string(signalScale)).c_str(), "L");
      leg2->AddEntry(totalBkg,"total background","FL");

    }
    
    leg2->Draw("same");
    
    TH1* totalSignal = NULL;
    if(isHiggsInvisible){
      totalSignal = (TH1*) ggHhist->Clone("totalSignal");
      totalSignal->Add(vbfHhist);
      totalSignal->Add(wHhist);
      totalSignal->Add(zHhist);
    }
    else{
      totalSignal = (TH1*) monoJhist->Clone("monoJhist");
      totalSignal->Add(monoWhist);
      totalSignal->Add(monoZhist);
    }

    TPad *pad3 = new TPad("pad3","pad3",0,0.,1,0.9);
    pad3->SetTopMargin(0.7);
    pad3->SetRightMargin(0.06);
    pad3->SetFillColor(0);
    pad3->SetGridy(1);
    pad3->SetFillStyle(0);
    pad3->Draw();
    pad3->cd();
    frame2->GetYaxis()->SetTitle("(S+B)/B");

    TH1* SoverB_prefit = (TH1*) totalSignal->Clone("SoverB_prefit");
    SoverB_prefit->Add(totalBkg);
    SoverB_prefit->Divide(totalBkg);
    TH1* SoverB_prefit_d = (TH1*) SoverB_prefit->Clone("SoverB_prefit_d");
    for(int iBin = 0; iBin < SoverB_prefit_d->GetNbinsX(); iBin++)
      SoverB_prefit_d->SetBinContent(iBin+1,1);
    frame2->GetYaxis()->SetRangeUser(std::max(SoverB_prefit->GetMinimum(),SoverB_prefit_d->GetMinimum())*0.5,std::max(SoverB_prefit->GetMaximum(),SoverB_prefit_d->GetMaximum())*1.2);
    frame2->Draw();
    SoverB_prefit->Draw("hist same");
    SoverB_prefit_d->SetLineColor(0);
    SoverB_prefit_d->SetMarkerColor(0);
    SoverB_prefit_d->SetMarkerSize(0);
    SoverB_prefit_d->SetFillColor(kGray);
    SoverB_prefit_d->SetFillStyle(1001);
    SoverB_prefit_d->Draw("E2 SAME");
    unhist->Draw("SAME");
    SoverB_prefit->Draw("hist same");
    pad3->RedrawAxis("sameaxis");
    
    canvas->SaveAs((observable+"_"+controlRegion+"_Shape.png").c_str());
    canvas->SaveAs((observable+"_"+controlRegion+"_Shape.pdf").c_str());

  }
}

