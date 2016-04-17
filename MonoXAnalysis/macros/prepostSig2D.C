#include "CMS_lumi.h"
#include "makehist.h"

void prepostSig2D(string fitFilename, // from combine, mlfit
		  string templateFileName, // input templates
		  string observable, 
		  int    category, 
		  bool   isHiggsInvisible, 
		  int    scaleSig = 1, 
		  bool   alongX  = false,
		  bool   blind   = true, 
		  bool   plotSBFit = false, 
		  string interaction = "Vector", 
		  string mediatorMass = "2000", 
		  string DMMass = "10") {
  

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

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

    TCanvas* canvas = new TCanvas(Form("canvas_%d",int(ihist)), "canvas", 600, 700);
    canvas->SetTickx(1);
    canvas->SetTicky(1);
    canvas->cd();
    canvas->SetBottomMargin(0.3);
    canvas->SetRightMargin(0.06);

    TPad *pad2 = new TPad(Form("pad2_%d",int(ihist)),"pad2",0,0.,1,0.9);
    pad2->SetTopMargin(0.7);
    pad2->SetRightMargin(0.06);
    pad2->SetFillColor(0);
    pad2->SetFillStyle(0);
    pad2->SetLineColor(0);
    pad2->SetGridy();

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

    //signal style                                                                                                                                                                  //signal style  
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
      mwhist.at(ihist)->SetLineColor(kBlue);
      mwhist.at(ihist)->SetLineWidth(3);
      mwhist.at(ihist)->Scale(scaleSig);
    }

    if(vbfhist.size() > ihist){
      vbfhist.at(ihist)->SetFillColor(0);
      vbfhist.at(ihist)->SetFillStyle(0);
      vbfhist.at(ihist)->SetLineColor(kViolet);
      vbfhist.at(ihist)->SetLineWidth(3);
      vbfhist.at(ihist)->Scale(scaleSig);
    }

    if(mzhist.size() > ihist){
      mzhist.at(ihist)->SetFillColor(0);
      mzhist.at(ihist)->SetFillStyle(0);
      mzhist.at(ihist)->SetLineColor(TColor::GetColor("#A2C523"));
      mzhist.at(ihist)->SetLineWidth(3);
      mzhist.at(ihist)->Scale(scaleSig);
    }

    if(wHhist.size() > ihist){
      wHhist.at(ihist)->SetFillColor(0);
      wHhist.at(ihist)->SetFillStyle(0);
      wHhist.at(ihist)->SetLineColor(TColor::GetColor("#A2C523"));
      wHhist.at(ihist)->SetLineWidth(3);
      wHhist.at(ihist)->Scale(scaleSig);
    }

    if(zHhist.size() > ihist){
      zHhist.at(ihist)->SetFillColor(0);
      zHhist.at(ihist)->SetFillStyle(0);
      zHhist.at(ihist)->SetLineColor(TColor::GetColor("#A2C523"));
      zHhist.at(ihist)->SetLineWidth(3);
      zHhist.at(ihist)->Scale(scaleSig);
    }

    znhist.at(ihist)->SetFillColor(TColor::GetColor("#258039"));
    znhist.at(ihist)->SetLineColor(kBlack);
    
    gmhist.at(ihist)->SetFillColor(TColor::GetColor("#9A9EAB"));
    gmhist.at(ihist)->SetLineColor(TColor::GetColor("#9A9EAB"));
    
    zlhist.at(ihist)->SetFillColor(TColor::GetColor("#9A9EAB"));
    zlhist.at(ihist)->SetLineColor(kBlack);
    
    wlhist.at(ihist)->SetFillColor(TColor::GetColor("#FAAF08"));
    wlhist.at(ihist)->SetLineColor(kBlack);

    tthist.at(ihist)->SetFillColor(TColor::GetColor("#CF3721"));
    tthist.at(ihist)->SetLineColor(kBlack);
    
    dihist.at(ihist)->SetFillColor(TColor::GetColor("#4897D8"));
    dihist.at(ihist)->SetLineColor(kBlack);
    
    qchist.at(ihist)->SetFillColor(TColor::GetColor("#F1F1F2"));
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
      frame->GetYaxis()->SetRangeUser(0.002,wlhist.at(ihist)->GetMaximum()*100);
    else
      frame->GetYaxis()->SetRangeUser(0.0005,wlhist.at(ihist)->GetMaximum()*200);

    frame->GetXaxis()->SetTitleSize(0);
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetTitle("Events / GeV");
    frame->GetYaxis()->SetTitleOffset(1.15);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetTitleSize(0.050);
    if(category <= 1)
      frame->GetXaxis()->SetNdivisions(510);
    else
      frame->GetXaxis()->SetNdivisions(504);
    frame ->Draw();

    CMS_lumi(canvas,"2.3");

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
    
    TLegend* leg = new TLegend(0.6, 0.55, 0.9, 0.9);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    TLatex ttext;
    ttext.SetNDC();
    ttext.SetTextFont(42);
    ttext.SetTextAlign(31);
    ttext.SetTextSize(0.04);

    if(ihist < bin.size()-2)
      ttext.DrawLatex(0.35,0.75,Form("%.1f <= %s < %.1f ",bin.at(ihist),text.second.c_str(),bin.at(ihist+1)));
    else
      ttext.DrawLatex(0.35,0.75,Form("%s >= %.1f ",text.second.c_str(),bin.at(ihist)));

    
    canvas->RedrawAxis("sameaxis");
    canvas->SetLogy();
    canvas->cd();

    pad2->Draw();
    pad2->cd();

    TH1* frame2 =  (TH1*) dthist.at(ihist)->Clone(Form("frame2_bin_%d",int(ihist)));
    frame2->Reset();
    if(category <=1)
      frame2->GetYaxis()->SetRangeUser(0.5,1.5);
    else
      frame2->GetYaxis()->SetRangeUser(-0.5,2);

    if(category <= 1)
      frame2->GetXaxis()->SetNdivisions(510);
    else
      frame2->GetXaxis()->SetNdivisions(210);
    frame2->GetXaxis()->SetTitle(text.first.c_str());
    frame2->GetYaxis()->SetTitle("Data/Pred.");
    frame2->GetYaxis()->CenterTitle();
    frame2->GetYaxis()->SetTitleOffset(1.5);
    frame2->GetYaxis()->SetLabelSize(0.03);
    frame2->GetYaxis()->SetTitleSize(0.04);
    frame2->GetXaxis()->SetLabelSize(0.04);
    frame2->GetXaxis()->SetTitleSize(0.05);
    frame2->GetYaxis()->SetNdivisions(5);
    frame2->Draw("AXIS");
    frame2->Draw("AXIG same");
  
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
    TH1F* htemp = (TH1F*) tohist.at(ihist)->Clone("htemp");
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

    TLegend* leg2 = new TLegend(0.14,0.14,0.65,0.21,NULL,"brNDC");
    leg2->SetFillColor(0);
    leg2->SetFillStyle(1);
    leg2->SetBorderSize(0);
    leg2->SetLineColor(0);
    leg2->SetNColumns(2);
    leg2->AddEntry(dahist,"Backgraound (Post-Fit)","PLE");
    leg2->AddEntry(dphist,"Backgraound (Pre-Fit)","PLE");
    leg2->Draw("same");


    pad2->RedrawAxis("G sameaxis");

    canvas->cd();
    
    leg->AddEntry(dthist.at(ihist), "Data", "PEL");
    if(sighist.size() > ihist && plotSBFit)
      leg->AddEntry(sighist.at(ihist), "Fitted Total Mono-X Signal", "F");

    leg->AddEntry(znhist.at(ihist), "Z #rightarrow #nu#nu", "F");
    leg->AddEntry(wlhist.at(ihist), "W #rightarrow l#nu", "F");
    leg->AddEntry(dihist.at(ihist), "WW/WZ/ZZ", "F");
    leg->AddEntry(tthist.at(ihist), "Top Quark", "F");
    leg->AddEntry(zlhist.at(ihist), "Z/#gamma #rightarrow ll, #gamma+jets", "F");
    leg->AddEntry(qchist.at(ihist), "QCD", "F");

    if(mjhist.size() > ihist && !plotSBFit)
      leg->AddEntry(mjhist.at(ihist), Form("Mono-J (V,2 TeV x%d)",scaleSig),"L");

    if(mwhist.size() > ihist && !plotSBFit)
      leg->AddEntry(mwhist.at(ihist), Form("Mono-W (V,2 TeV x%d)",scaleSig),"L");

    if(mzhist.size() > ihist && !plotSBFit)
      leg->AddEntry(mzhist.at(ihist), Form("Mono-Z (V,2 TeV x%d)",scaleSig),"L");

    if(ggHhist.size() > ihist && !plotSBFit)
      leg->AddEntry(ggHhist.at(ihist), "ggH m_{H} = 125 GeV","L");

    if(vbfhist.size() > ihist && !plotSBFit)
      leg->AddEntry(vbfhist.at(ihist), "qqH m_{H} = 125 GeV","L");
    if(wHhist.size() > ihist && !plotSBFit)
      leg->AddEntry(wHhist.at(ihist), "VH m_{H} = 125 GeV","L");


    leg->Draw("SAME");
    pad2->RedrawAxis("sameaxis");
    canvas->RedrawAxis("sameaxis");
    canvas->SetLogy();
    
    
    if(blind){
      canvas->SaveAs(Form("postfit_sig_blind_bin_%d.pdf",int(ihist)));
      canvas->SaveAs(Form("postfit_sig_blind_bin_%d.png",int(ihist)));
    }
    else{
      canvas->SaveAs(Form("postfit_sig_bin_%d.pdf",int(ihist)));
      canvas->SaveAs(Form("postfit_sig_bin_%d.png",int(ihist)));
    }   

    TH1* totalSignal = NULL;

    if(isHiggsInvisible and not plotSBFit){
      totalSignal = (TH1*) ggHhist.at(ihist)->Clone("totalSignal");
      totalSignal->Add(vbfhist.at(ihist));
      totalSignal->Add(wHhist.at(ihist));
      totalSignal->Add(zHhist.at(ihist));
    }
    else if(not isHiggsInvisible and not plotSBFit){
      totalSignal = (TH1*) mjhist.at(ihist)->Clone("totalSignal");
      totalSignal->Add(mwhist.at(ihist));
      totalSignal->Add(mzhist.at(ihist));
    }
    else if(plotSBFit)
      totalSignal = (TH1*) sighist.at(ihist)->Clone("totalSignal");

    canvas->cd();
    pad2->Draw();
    pad2->cd();
    if(not plotSBFit)
      frame2->GetYaxis()->SetTitle("(S+B)/B");
    else
      frame2->GetYaxis()->SetTitle("(S_{fit}+B)/B");

    TH1* SoverB_prefit = (TH1*) totalSignal->Clone("SoverB_prefit");
    TH1* SoverB_postfit = (TH1*) totalSignal->Clone("SoverB_postfit");
    SoverB_prefit->SetLineColor(kRed);
    SoverB_prefit->SetMarkerColor(kRed);
    SoverB_prefit->SetMarkerSize(1);
    SoverB_prefit->SetMarkerStyle(20);
    SoverB_postfit->SetLineColor(kBlue);
    SoverB_postfit->SetMarkerColor(kBlue);
    SoverB_postfit->SetMarkerSize(1);
    SoverB_postfit->SetMarkerStyle(20);

    SoverB_prefit->Add(tphist.at(ihist));
    SoverB_prefit->Divide(tphist.at(ihist));
    SoverB_postfit->Add(htemp);
    SoverB_postfit->Divide(htemp);

    frame2->GetYaxis()->SetRangeUser(0.5,SoverB_postfit->GetMaximum()*1.2);
    frame2->Draw();

    SoverB_postfit->Draw("hist same");
    TH1* SoverB_postfit_d = (TH1*) SoverB_postfit->Clone("SoverB_postfit_d");
    for(int iBin = 0; iBin < SoverB_postfit_d->GetNbinsX(); iBin++)
      SoverB_postfit_d->SetBinContent(iBin+1,1);
    SoverB_postfit_d->SetLineColor(0);

    SoverB_postfit->Draw("hist same");
    SoverB_postfit_d->SetMarkerColor(0);
    SoverB_postfit_d->SetMarkerSize(0);
    SoverB_postfit_d->SetFillColor(kGray);
    SoverB_postfit_d->SetFillStyle(1001);
    SoverB_postfit_d->Draw("E2 SAME");
    unhist->Draw("SAME");
    SoverB_prefit->Draw("hist same");
    SoverB_postfit->Draw("hist same");
    pad2->RedrawAxis("sameaxis");

    canvas->SaveAs(Form("postfit_sig_SoB_bin_%d.pdf",int(ihist)));
    canvas->SaveAs(Form("postfit_sig_SoB_bin_%d.png",int(ihist)));
   
  } 
}

