#include "CMS_lumi.h"
#include "makehist.h"

void prepostSig(string fitFilename, // from combine, mlfit
		string templateFileName, // input templates
		string observable, int category, 
		bool isHiggsInvisible, 
		int scaleSig = 1, 
		bool alongX  = false,
		bool blind   = true, 
		bool plotSBFit = false, 
		string interaction = "Vector", string mediatorMass = "2000", string DMMass = "10") {
  

  gROOT->SetBatch(kTRUE);

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

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.295);
  pad2->SetTickx();
  pad2->SetTicky();

  TFile* pfile = new TFile(fitFilename.c_str());
  TFile* dfile = new TFile(templateFileName.c_str());

  // in case of b-only fit just dispaly three possible signal on the stack
  vector<TH1F*> mjhist;
  vector<TH1F*> mwhist;
  vector<TH1F*> mzhist;
  vector<TH1F*> ggHhist;
  vector<TH1F*> vbfhist;
  vector<TH1F*> wHhist;
  vector<TH1F*> zHhist;
  vector<TH1F*> zhhist;
  
  if(! isHiggsInvisible){
    mjhist = transformUnrolledHistogram((TH1*) dfile->FindObjectAny(("monoJhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()),observable,category,alongX);
    mwhist = transformUnrolledHistogram((TH1*) dfile->FindObjectAny(("monoWhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()),observable,category,alongX);
    mzhist = transformUnrolledHistogram((TH1*) dfile->FindObjectAny(("monoZhist_"+interaction+"_"+mediatorMass+"_"+DMMass+"_"+observable).c_str()),observable,category,alongX);
    for(size_t ihist = 0; ihist < mjhist.size(); ihist++){
      mjhist.at(ihist)->Scale(1.0, "width");
      mwhist.at(ihist)->Scale(1.0, "width");
      mzhist.at(ihist)->Scale(1.0, "width");
    }
  }
  else{
    ggHhist  = transformUnrolledHistogram((TH1*) dfile->FindObjectAny(("ggHhist_"+mediatorMass+"_"+observable).c_str()),observable,category,alongX);
    vbfhist = transformUnrolledHistogram((TH1*) dfile->FindObjectAny(("vbfHhist_"+mediatorMass+"_"+observable).c_str()),observable,category,alongX);
    wHhist   = transformUnrolledHistogram((TH1*) dfile->FindObjectAny(("wHhist_"+mediatorMass+"_"+observable).c_str()),observable,category,alongX);
    zHhist   = transformUnrolledHistogram((TH1*) dfile->FindObjectAny(("zHhist_"+mediatorMass+"_"+observable).c_str()),observable,category,alongX);
    for(size_t ihist = 0; ihist < ggHhist.size(); ihist++){
      ggHhist.at(ihist)->Scale(1.0, "width");
      vbfhist.at(ihist)->Scale(1.0, "width");
      wHhist.at(ihist)->Scale(1.0, "width");
      zHhist.at(ihist)->Scale(1.0, "width");
    }
  }

  vector<TH1F*> znhist;
  vector<TH1F*> zlhist;
  vector<TH1F*> wlhist;
  vector<TH1F*> tthist;
  vector<TH1F*> dihist;
  vector<TH1F*> qchist;
  vector<TH1F*> gmhist;
  vector<TH1F*> tohist;
  vector<TH1F*> tphist;
  vector<TH1F*> sighist;

  if(!plotSBFit){
    znhist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_b/ch1/Znunu"),observable,category,alongX);    
    zlhist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_b/ch1/ZJets"),observable,category,alongX);    
    wlhist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_b/ch1/WJets"),observable,category,alongX);    
    tthist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_b/ch1/Top"),observable,category,alongX);    
    dihist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_b/ch1/Dibosons"),observable,category,alongX);
    qchist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_b/ch1/QCD"),observable,category,alongX);    
    gmhist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_b/ch1/GJets"),observable,category,alongX);    
    tohist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_b/ch1/total_background"),observable,category,alongX);    
    tphist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_prefit/ch1/total_background"),observable,category,alongX,"prefit");        
    for(size_t ihist = 0; ihist < znhist.size(); ihist++){
      znhist.at(ihist)->Scale(1.0,"width");
      zlhist.at(ihist)->Scale(1.0,"width");
      wlhist.at(ihist)->Scale(1.0,"width");
      tthist.at(ihist)->Scale(1.0,"width");
      dihist.at(ihist)->Scale(1.0,"width");
      qchist.at(ihist)->Scale(1.0,"width");
      gmhist.at(ihist)->Scale(1.0,"width");
      tohist.at(ihist)->Scale(1.0,"width");
      tphist.at(ihist)->Scale(1.0,"width");
    }
  }
  else{
    znhist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch1/Znunu"),observable,category,alongX);    
    zlhist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch1/ZJets"),observable,category,alongX);    
    wlhist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch1/WJets"),observable,category,alongX);    
    tthist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch1/Top"),observable,category,alongX);    
    dihist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch1/Dibosons"),observable,category,alongX);    
    qchist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch1/QCD"),observable,category,alongX);    
    gmhist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch1/GJets"),observable,category,alongX);    
    tohist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch1/total_background"),observable,category,alongX);    
    tphist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_prefit/ch1/total_background"),observable,category,alongX,"prefit");      
    sighist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch1/total_signal"),observable,category,alongX);
    for(size_t ihist = 0; ihist < znhist.size(); ihist++){
      znhist.at(ihist)->Scale(1.0,"width");
      zlhist.at(ihist)->Scale(1.0,"width");
      wlhist.at(ihist)->Scale(1.0,"width");
      tthist.at(ihist)->Scale(1.0,"width");
      dihist.at(ihist)->Scale(1.0,"width");
      qchist.at(ihist)->Scale(1.0,"width");
      gmhist.at(ihist)->Scale(1.0,"width");
      tohist.at(ihist)->Scale(1.0,"width");
      tphist.at(ihist)->Scale(1.0,"width");
      sighist.at(ihist)->Scale(1.0,"width");
    }
  }

  vector<TH1F*> dthist;
  if(!blind){    
    dthist = transformUnrolledHistogram((TH1*)dfile->FindObjectAny(("datahist_"+observable).c_str()),observable,category,alongX);
    for(size_t ihist = 0; ihist < dthist.size(); ihist++)
      dthist.at(ihist)->Scale(1.0,"width");
  }
  else{

    dthist = transformUnrolledHistogram((TH1*)dfile->FindObjectAny(("datahist_"+observable).c_str()),observable,category,alongX);

    for(size_t ihist = 0; ihist < dthist.size(); ihist++){
      dthist.at(ihist)->Scale(1.0,"width");
      for (int i = 0; i <= dthist.at(ihist)->GetNbinsX(); i++) {
	double yield = 0.0;
	yield += zlhist.at(ihist)->GetBinContent(i);
	yield += gmhist.at(ihist)->GetBinContent(i);
	yield += wlhist.at(ihist)->GetBinContent(i);
	yield += tthist.at(ihist)->GetBinContent(i);
	yield += dihist.at(ihist)->GetBinContent(i);
	yield += qchist.at(ihist)->GetBinContent(i);
	yield += znhist.at(ihist)->GetBinContent(i);
	dthist.at(ihist)->SetBinContent(i, yield);
	dthist.at(ihist)->SetBinError(i, 0.);
      }
    }
  }

  pair<string,string> text = observableName(observable,alongX);

  bin2D bins = selectBinning2D(observable,category);
  vector<double> bin ;
  if(alongX)
    bin = bins.binX ;
  else
    bin = bins.binY ;

  if(bin.size()-1 != dthist.size())
    cerr<<"Huston we have a problem with bin size ... "<<endl;

  for(size_t ihist = 0; ihist < dthist.size(); ihist++){
    // print in text file yields
    ofstream outputfile;
    outputfile.open(Form("prepostSR_bin_%d.txt",int(ihist)));
    
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
    stringstream PostRateUnc;
    PostRateUnc << "Process: Post-fit uncertainty (total)";
    stringstream DataRate;
    DataRate << "Process: Data";
    
    for(int iBin = 0; iBin < qchist.at(ihist)->GetNbinsX(); iBin++){
      QCDRate << "   ";
      QCDRate << qchist.at(ihist)->GetBinContent(iBin+1);
    }

    for(int iBin = 0; iBin < gmhist.at(ihist)->GetNbinsX(); iBin++){
      GJetsRate << "   ";
      GJetsRate << gmhist.at(ihist)->GetBinContent(iBin+1);
    }

    for(int iBin = 0; iBin < dihist.at(ihist)->GetNbinsX(); iBin++){
      DiBosonRate << "   ";
      DiBosonRate << dihist.at(ihist)->GetBinContent(iBin+1);
    }

    for(int iBin = 0; iBin < tthist.at(ihist)->GetNbinsX(); iBin++){
      TopRate << "   ";
      TopRate << tthist.at(ihist)->GetBinContent(iBin+1);
    }
    
    for(int iBin = 0; iBin < zlhist.at(ihist)->GetNbinsX(); iBin++){
      ZJetsRate << "   ";
      ZJetsRate << zlhist.at(ihist)->GetBinContent(iBin+1);
    }
    
    for(int iBin = 0; iBin < wlhist.at(ihist)->GetNbinsX(); iBin++){
      WJetsRate << "   ";
      WJetsRate << wlhist.at(ihist)->GetBinContent(iBin+1);
    }
    
    for(int iBin = 0; iBin < znhist.at(ihist)->GetNbinsX(); iBin++){
      ZnunuRate << "   ";
      ZnunuRate << znhist.at(ihist)->GetBinContent(iBin+1);
    }
    
    for(int iBin = 0; iBin < tphist.at(ihist)->GetNbinsX(); iBin++){
      PreRate << "   ";
      PreRate << tphist.at(ihist)->GetBinContent(iBin+1);
    }
    
    for(int iBin = 0; iBin < tohist.at(ihist)->GetNbinsX(); iBin++){
      PostRate << "   ";
      PostRate << tohist.at(ihist)->GetBinContent(iBin+1);
      PostRateUnc << "   ";
      PostRateUnc << tohist.at(ihist)->GetBinError(iBin+1);
    }
    
    for(int iBin = 0; iBin < dthist.at(ihist)->GetNbinsX(); iBin++){
      DataRate << "   ";
      DataRate << dthist.at(ihist)->GetBinContent(iBin+1);
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
    outputfile<<PostRateUnc.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<DataRate.str()<<endl;
    outputfile<<"######################"<<endl;

    outputfile.close();
    
    // pad 1
    pad1->SetRightMargin(0.06);
    pad1->SetLeftMargin(0.12);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();
    pad1->cd();
  
    //signal style  
    if(mjhist.size() > ihist){
      mjhist.at(ihist)->SetFillColor(0);
      mjhist.at(ihist)->SetFillStyle(0);
      mjhist.at(ihist)->SetLineColor(kBlack);
      mjhist.at(ihist)->SetLineWidth(3);
      mjhist.at(ihist)->Scale(scaleSig);
    }

    if(ggHhist.size() > ihist){
      ggHhist.at(ihist)->SetFillColor(0);
      ggHhist.at(ihist)->SetFillStyle(0);
      ggHhist.at(ihist)->SetLineColor(kBlack);
      ggHhist.at(ihist)->SetLineWidth(3);
      ggHhist.at(ihist)->Scale(scaleSig);
    }

    if(mwhist.size() > ihist){
      mwhist.at(ihist)->SetFillColor(0);
      mwhist.at(ihist)->SetFillStyle(0);
      mwhist.at(ihist)->SetLineColor(kViolet);
      mwhist.at(ihist)->SetLineWidth(3);
      mwhist.at(ihist)->SetLineStyle(7);
      mwhist.at(ihist)->Scale(scaleSig);
    }

    if(vbfhist.size() > ihist){
      vbfhist.at(ihist)->SetFillColor(0);
      vbfhist.at(ihist)->SetFillStyle(0);
      vbfhist.at(ihist)->SetLineColor(kViolet);
      vbfhist.at(ihist)->SetLineWidth(3);
      vbfhist.at(ihist)->SetLineStyle(7);
      vbfhist.at(ihist)->Scale(scaleSig);
    }

    if(mzhist.size() > ihist){
      mzhist.at(ihist)->SetFillColor(0);
      mzhist.at(ihist)->SetFillStyle(0);
      mzhist.at(ihist)->SetLineColor(kOrange+1);
      mzhist.at(ihist)->SetLineWidth(3);
      mzhist.at(ihist)->SetLineStyle(7);
      mzhist.at(ihist)->Scale(scaleSig);
    }

    if(wHhist.size() > ihist){
      wHhist.at(ihist)->SetFillColor(0);
      wHhist.at(ihist)->SetFillStyle(0);
      wHhist.at(ihist)->SetLineColor(kOrange+1);
      wHhist.at(ihist)->SetLineWidth(3);
      wHhist.at(ihist)->SetLineStyle(7);
      wHhist.at(ihist)->Scale(scaleSig);
    }

    if(zHhist.size() > ihist){
      zHhist.at(ihist)->SetFillColor(0);
      zHhist.at(ihist)->SetFillStyle(0);
      zHhist.at(ihist)->SetLineColor(kGreen);
      zHhist.at(ihist)->SetLineWidth(3);
      zHhist.at(ihist)->SetLineStyle(7);
      zHhist.at(ihist)->Scale(scaleSig);
    }

    znhist.at(ihist)->SetFillColor(kGreen+1);
    znhist.at(ihist)->SetLineColor(kBlack);
    
    gmhist.at(ihist)->SetFillColor(13);
    gmhist.at(ihist)->SetLineColor(13);
    
    zlhist.at(ihist)->SetFillColor(13);
    zlhist.at(ihist)->SetLineColor(kBlack);
    
    wlhist.at(ihist)->SetFillColor(kRed);
    wlhist.at(ihist)->SetLineColor(kBlack);

    tthist.at(ihist)->SetFillColor(kBlue);
    tthist.at(ihist)->SetLineColor(kBlack);
    
    dihist.at(ihist)->SetFillColor(13);
    dihist.at(ihist)->SetLineColor(13);
    
    qchist.at(ihist)->SetFillColor(kGray);
    qchist.at(ihist)->SetLineColor(kBlack);

    if(sighist.size() > ihist){
      sighist.at(ihist)->SetFillColor(kBlack);
      sighist.at(ihist)->SetLineColor(kBlack);
      sighist.at(ihist)->SetFillStyle(3001);
    }
    
    // make the stack
    THStack* stack = new THStack(Form("stack_bin_%d",int(ihist)), Form("stack_bin_%d",int(ihist)));
    stack->Add(gmhist.at(ihist));
    stack->Add(dihist.at(ihist));
    stack->Add(zlhist.at(ihist)); 
    stack->Add(tthist.at(ihist));
    stack->Add(qchist.at(ihist));
    stack->Add(wlhist.at(ihist));
    stack->Add(znhist.at(ihist));

    if(plotSBFit && sighist.size() > ihist)
      stack->Add(sighist.at(ihist));

    TH1* frame = (TH1*) dthist.at(ihist)->Clone(Form("frame_bin_%d",int(ihist)));
    frame->Reset();
    if(category <=1)
      frame->GetYaxis()->SetRangeUser(0.0005,1000000);
    else
      frame->GetYaxis()->SetRangeUser(0.002,9000);
    
    frame->GetXaxis()->SetTitleSize(0);
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetTitle("Events / GeV");
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame ->Draw();
    
    CMS_lumi(pad1,"2.30",false,true);
    
    stack ->Draw("HIST SAME");

    if(mjhist.size() > ihist && !plotSBFit)
      mjhist.at(ihist)->Draw("HIST SAME");
    if(mwhist.size() > ihist && !plotSBFit)
      mwhist.at(ihist)->Draw("HIST SAME");
    if(mzhist.size() > ihist && !plotSBFit)
      mzhist.at(ihist)->Draw("HIST SAME");    
    if(ggHhist.size() > ihist && !plotSBFit)
      ggHhist.at(ihist)->Draw("HIST SAME");
    if(vbfhist.size() > ihist && !plotSBFit)
      vbfhist.at(ihist)->Draw("HIST SAME");
    if(wHhist.size() > ihist && !plotSBFit)
      wHhist.at(ihist)->Draw("HIST SAME");
    if(zHhist.size() > ihist && !plotSBFit)
      zHhist.at(ihist)->Draw("HIST SAME");

    dthist.at(ihist)->SetMarkerSize(1.2);
    dthist.at(ihist)->SetMarkerStyle(20);
    dthist.at(ihist)->SetFillStyle(0);
    dthist.at(ihist)->SetFillColor(0);
    dthist.at(ihist)->SetLineColor(kBlack);
    dthist.at(ihist)->SetLineWidth(1);
    dthist.at(ihist)->SetMarkerColor(kBlack);
    dthist.at(ihist)->Draw("PE SAME");
    
    TLegend* leg = new TLegend(0.38, 0.38, 0.93, 0.92);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);  
    leg->SetNColumns(2);

    TLatex ttext;
    ttext.SetNDC();
    ttext.SetTextFont(42);
    ttext.SetTextAlign(31);
    ttext.SetTextSize(0.04);

    if(ihist < bin.size()-2)
      ttext.DrawLatex(0.35,0.75,Form("%.1f <= %s < %.1f ",bin.at(ihist),text.second.c_str(),bin.at(ihist+1)));
    else
      ttext.DrawLatex(0.35,0.75,Form("%s >= %.1f ",text.second.c_str(),bin.at(ihist)));

    
    canvas->cd();
    pad2->SetTopMargin(0.08);
    pad2->SetRightMargin(0.06);
    pad2->SetLeftMargin(0.12);
    pad2->SetBottomMargin(0.35);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    TH1* frame2 =  (TH1*) dthist.at(ihist)->Clone(Form("frame2_bin_%d",int(ihist)));
    frame2->Reset("ICES");
    
    if(category <=1)
      frame2->GetYaxis()->SetRangeUser(0.5,1.5);
    else
      frame2->GetYaxis()->SetRangeUser(0.25,1.75);
    
    frame2->GetXaxis()->SetTitle(text.first.c_str());
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
  
    // for post-fit pre-fit data/mc
    TH1* dphist = (TH1*)dthist.at(ihist)->Clone(Form("dahist_bin_%d",int(ihist)));
    TH1* dahist = (TH1*)dthist.at(ihist)->Clone(Form("dahist_bin_%d",int(ihist)));

    dphist->SetLineColor(kRed);
    dahist->SetLineColor(kBlue);
    dphist->SetMarkerColor(kRed);
    dphist->SetMarkerSize(1);
    dphist->SetMarkerStyle(20);
    dahist->SetMarkerColor(kBlue);
    dahist->SetMarkerSize(1);
    dahist->SetMarkerStyle(20);
    
    TH1* mphist = (TH1*)tphist.at(ihist)->Clone(Form("mphist_bin_%d",int(ihist)));
    TH1* mchist = (TH1*)zlhist.at(ihist)->Clone(Form("mchist_bin_%d",int(ihist)));
    TH1* unhist = (TH1*)zlhist.at(ihist)->Clone(Form("unhist_bin_%d",int(ihist)));
    mchist->Add(wlhist.at(ihist));
    mchist->Add(gmhist.at(ihist));
    mchist->Add(tthist.at(ihist));
    mchist->Add(dihist.at(ihist));
    mchist->Add(qchist.at(ihist));
    mchist->Add(znhist.at(ihist));
    if(sighist.size() > ihist && plotSBFit)
      mchist->Add(sighist.at(ihist));

    for (int i = 1; i <= mchist->GetNbinsX(); i++) mchist->SetBinError(i, 0);
    for (int i = 1; i <= mphist->GetNbinsX(); i++) mphist->SetBinError(i, 0);

    // ratio data/post-fit
    dahist->Divide(mchist);
    // ratio data/pre-fit
    dphist->Divide(mphist);
    
    // ratio post fit at 1 with uncertaitny
    tohist.at(ihist)->Divide(mchist);
    tohist.at(ihist)->SetLineColor(0);
    tohist.at(ihist)->SetMarkerColor(0);
    tohist.at(ihist)->SetMarkerSize(0);
    tohist.at(ihist)->SetFillColor(kGray);
    
    dahist->SetMarkerSize(1);
    dphist->SetMarkerSize(1);
    dahist->SetMarkerStyle(20);
    dphist->SetMarkerStyle(20);
    
    dahist->SetStats(kFALSE);
    dphist->SetStats(kFALSE);
    
    // line at 1
    for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinContent(i, 1);
    for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinError(i, 0);
    unhist->SetMarkerSize(0);
    unhist->SetLineColor(kBlack);
    unhist->SetLineStyle(2);
    unhist->SetLineWidth(2);
    unhist->SetFillColor(0);

    dahist->GetXaxis()->SetLabelOffset(999999);
    dahist->GetXaxis()->SetLabelSize(0);
    dahist->GetXaxis()->SetTitleOffset(999999);
    dahist->GetXaxis()->SetTitleSize(0);

    tohist.at(ihist)->Draw("E2 SAME");
    unhist->Draw("SAME");
    dahist->Draw("PE1 SAME");
    if(!blind)
      dphist->Draw("PE1 SAME");
    
    pad2->RedrawAxis("G sameaxis");

    canvas->cd();
    pad1->cd();
    
    leg->AddEntry(dthist.at(ihist), "Data", "PEL");
    if(sighist.size() > ihist && plotSBFit)
      leg->AddEntry(sighist.at(ihist), "Fitted Total Mono-X Signal", "F");

    leg->AddEntry(znhist.at(ihist), "Z #rightarrow #nu#nu", "F");
    leg->AddEntry(wlhist.at(ihist), "W #rightarrow l#nu", "F");
    leg->AddEntry(qchist.at(ihist), "QCD", "F");
    leg->AddEntry(tthist.at(ihist), "Top Quark", "F");
    leg->AddEntry(zlhist.at(ihist), "Others", "F");
    
    if(mjhist.size() > ihist && !plotSBFit)
      leg->AddEntry(mjhist.at(ihist), Form("Mono-J (V,2 TeV x%d)",scaleSig));
    
    if(mwhist.size() > ihist && !plotSBFit)
      leg->AddEntry(mwhist.at(ihist), Form("Mono-W (V,2 TeV x%d)",scaleSig));
    
    if(mzhist.size() > ihist && !plotSBFit)
      leg->AddEntry(mzhist.at(ihist), Form("Mono-Z (V,2 TeV x%d)",scaleSig));
    
    if(ggHhist.size() > ihist && !plotSBFit)
      leg->AddEntry(ggHhist.at(ihist), "gg #rightarrow H (m_{H} = 125 GeV)");
    
    if(vbfhist.size() > ihist && !plotSBFit)
      leg->AddEntry(vbfhist.at(ihist), "qq #rightarrow H (m_{H} = 125 GeV)");
    
    if(wHhist.size() > ihist && !plotSBFit)
      leg->AddEntry(wHhist.at(ihist), "qq #rightarrow WH (m_{H} = 125 GeV)");

    if(zHhist.size() > ihist && !plotSBFit)
      leg->AddEntry(zHhist.at(ihist), "qq #rightarrow ZH (m_{H} = 125 GeV)");

    leg->AddEntry(dahist,"Expected (post-fit)","PEL");
    leg->AddEntry(dphist,"Expected (pre-fit)","PEL");    
    leg->Draw("SAME");
  
    pad1->RedrawAxis("sameaxis");
    pad1->SetLogy();
    
    
    if(blind){
      canvas->SaveAs(Form("postfit_sig_blind_bin_%d.pdf",int(ihist)));
      canvas->SaveAs(Form("postfit_sig_blind_bin_%d.png",int(ihist)));
    }
    else{
      canvas->SaveAs(Form("postfit_sig_bin_%d.pdf",int(ihist)));
      canvas->SaveAs(Form("postfit_sig_bin_%d.png",int(ihist)));
    }   
  } 
}

