#include "../CMS_lumi.h"
#include "../makeTemplates/makehist.h"

void prepostZE2D(string fitFilename, string templateFileName, string observable, int category, bool alongX = false, bool plotSBFit = false) {

  gROOT->SetBatch(kTRUE); 
  setTDRStyle();
  
  TFile* pfile = new TFile(fitFilename.c_str());
  TFile* dfile = new TFile(templateFileName.c_str());

  vector<TH1F*> dthist; 
  vector<TH1F*> wlhist; 
  vector<TH1F*> tthist; 
  vector<TH1F*> dihist; 
  vector<TH1F*> pohist; 
  vector<TH1F*> prhist; 

  if(!plotSBFit){
    
    dthist = transformUnrolledHistogram((TH1*)dfile->FindObjectAny(("datahistzee_"+observable).c_str()),observable,category,alongX);
    wlhist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_b/ch5/WJets_ZE"),observable,category,alongX);
    tthist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_b/ch5/Top"),observable,category,alongX);
    dihist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_b/ch5/Dibosons"),observable,category,alongX);
    pohist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_b/ch5/total_background"),observable,category,alongX);
    prhist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_prefit/ch5/total_background"),observable,category,alongX,"prefit");
    
  }
  else{

    dthist = transformUnrolledHistogram((TH1*)dfile->FindObjectAny(("datahistzee_"+observable).c_str()),observable,category,alongX);
    wlhist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch5/WJets_ZE"),observable,category,alongX);
    tthist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch5/Top"),observable,category,alongX);
    dihist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch5/Dibosons"),observable,category,alongX);
    pohist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_fit_s/ch5/total_background"),observable,category,alongX);
    prhist = transformUnrolledHistogram((TH1*)pfile->Get("shapes_prefit/ch5/total_background"),observable,category,alongX,"prefit");

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


  for(size_t ihist = 0 ; ihist < dthist.size(); ihist++){

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
    pad2->SetGridy(1);
    pad2->SetFillStyle(0);


    dthist.at(ihist)->Scale(1.0,"width");
    wlhist.at(ihist)->Scale(1.0,"width");
    tthist.at(ihist)->Scale(1.0,"width");
    dihist.at(ihist)->Scale(1.0,"width");
    pohist.at(ihist)->Scale(1.0,"width");
    prhist.at(ihist)->Scale(1.0,"width");

    ofstream  outputfile;
    outputfile.open(Form("prepostZE_bin_%d.txt",int(ihist)));
    stringstream TopRate;
    TopRate << "Process: Top";
    stringstream VVRate;
    VVRate << "Process: DiBoson";
    stringstream WJetRate;
    WJetRate << "Process: WJet";
    stringstream PreRate;
    PreRate << "Process: Pre-fit (total)";
    stringstream PostRate;
    PostRate << "Process: Post-fit (total)";
    stringstream DataRate;
    DataRate << "Process: Data";
    
    for(int iBin = 0; iBin < tthist.at(ihist)->GetNbinsX(); iBin++){
      TopRate << "   ";
      TopRate << tthist.at(ihist)->GetBinContent(iBin);
    }

    for(int iBin = 0; iBin < dihist.at(ihist)->GetNbinsX(); iBin++){
      VVRate << "   ";
      VVRate << dihist.at(ihist)->GetBinContent(iBin);
    }
    
    for(int iBin = 0; iBin < wlhist.at(ihist)->GetNbinsX(); iBin++){
      WJetRate << "   ";
      WJetRate << wlhist.at(ihist)->GetBinContent(iBin);
    }
    
    for(int iBin = 0; iBin < prhist.at(ihist)->GetNbinsX(); iBin++){
      PreRate << "   ";
      PreRate << prhist.at(ihist)->GetBinContent(iBin);
    }
    
    for(int iBin = 0; iBin < pohist.at(ihist)->GetNbinsX(); iBin++){
      PostRate << "   ";
      PostRate << pohist.at(ihist)->GetBinContent(iBin);
    }  
    
    for(int iBin = 0; iBin < dthist.at(ihist)->GetNbinsX(); iBin++){
      DataRate << "   ";
      DataRate << dthist.at(ihist)->GetBinContent(iBin);
    }
    
    outputfile<<"######################"<<endl;
    outputfile<<TopRate.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<VVRate.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<WJetRate.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<PreRate.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<PostRate.str()<<endl;
    outputfile<<"######################"<<endl;
    outputfile<<DataRate.str()<<endl;
    
    outputfile.close();
    
    prhist.at(ihist)->SetLineColor(kRed);
    prhist.at(ihist)->SetLineWidth(2);
    pohist.at(ihist)->SetLineColor(kBlue);
    pohist.at(ihist)->SetLineWidth(2);
    prhist.at(ihist)->SetMarkerColor(kRed);
    pohist.at(ihist)->SetMarkerColor(kBlue);
    
    wlhist.at(ihist)->SetFillColor(kOrange+1);
    wlhist.at(ihist)->SetLineColor(kBlack);
    wlhist.at(ihist)->Add(tthist.at(ihist));
    wlhist.at(ihist)->Add(dihist.at(ihist));
    
    canvas->cd();
    TH1* frame = (TH1*) dthist.at(ihist)->Clone("frame");
    frame->Reset();
    if(category <=1)
      frame->GetYaxis()->SetRangeUser(0.001,wlhist.at(ihist)->GetMaximum()*100);
    else
      frame->GetYaxis()->SetRangeUser(0.0007,wlhist.at(ihist)->GetMaximum()*100);

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

    prhist.at(ihist)->Draw("HIST SAME");
    pohist.at(ihist)->Draw("HIST SAME");
    wlhist.at(ihist)->Draw("HIST SAME");
    
    dthist.at(ihist)->SetMarkerSize(1.2);
    dthist.at(ihist)->SetMarkerStyle(20);
    
    dthist.at(ihist)->Draw("EP SAME");
    
    TLegend* leg = new TLegend(0.5, 0.65, 0.90, 0.90);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->AddEntry(dthist.at(ihist), "Data","PEL");
    leg->AddEntry(pohist.at(ihist), "Post-fit (Z #rightarrow ee)","L");
    leg->AddEntry(prhist.at(ihist), "Pre-fit (Z #rightarrow ee)","L");
    leg->AddEntry(wlhist.at(ihist), "Other Backgrounds", "F");
    leg->Draw("SAME");
    
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
    TH1* frame2 =  (TH1*) dthist.at(ihist)->Clone("frame");
    frame2->Reset("ICES");

    if(category <=1)
      frame2->GetYaxis()->SetRangeUser(0.5,1.5);
    else
      frame2->GetYaxis()->SetRangeUser(0.25,1.75);

    if(category <= 1)
      frame2->GetXaxis()->SetNdivisions(510);
    else
      frame2->GetXaxis()->SetNdivisions(210);
    frame2->GetYaxis()->SetNdivisions(5);


    frame2->GetXaxis()->SetTitle(text.first.c_str());
    frame2->GetYaxis()->SetTitle("Data/Pred.");
    frame2->GetYaxis()->CenterTitle();
    frame2->GetYaxis()->SetTitleOffset(1.5);
    frame2->GetYaxis()->SetLabelSize(0.04);
    frame2->GetYaxis()->SetTitleSize(0.04);
    frame2->GetXaxis()->SetLabelSize(0.04);
    frame2->GetXaxis()->SetTitleSize(0.05);
    frame2->Draw();

    TH1* d1hist = (TH1*)dthist.at(ihist)->Clone(Form("d1hist_bin_%d",int(ihist)));
    TH1* d2hist = (TH1*)dthist.at(ihist)->Clone(Form("d2hist_bin_%d",int(ihist)));
    TH1* m1hist = (TH1*)prhist.at(ihist)->Clone(Form("m1hist_bin_%d",int(ihist)));
    TH1* m2hist = (TH1*)pohist.at(ihist)->Clone(Form("m2hist_bin_%d",int(ihist)));
    TH1* erhist = (TH1*)pohist.at(ihist)->Clone(Form("erhist_bin_%d",int(ihist)));
    
    d1hist->SetLineColor(kRed);
    d2hist->SetLineColor(kBlue);
    d1hist->SetMarkerColor(kRed);
    d1hist->SetMarkerSize(1);
    d1hist->SetMarkerStyle(20);
    d2hist->SetMarkerColor(kBlue);
    d2hist->SetMarkerSize(1);
    d2hist->SetMarkerStyle(20);
    
    for (int i = 1; i <= m1hist->GetNbinsX(); i++) m1hist->SetBinError(i, 0);
    for (int i = 1; i <= m2hist->GetNbinsX(); i++) m2hist->SetBinError(i, 0);
    
    d1hist->Divide(m1hist);
    d2hist->Divide(m2hist);
    erhist->Divide(m2hist);
    erhist->SetLineColor(0);
    erhist->SetMarkerColor(0);
    erhist->SetMarkerSize(0);
    erhist->SetFillColor(kGray);
    
    d1hist->SetMarkerSize(1.2);
    d2hist->SetMarkerSize(1.2);
    d1hist->SetStats(kFALSE);
  
    d1hist->GetXaxis()->SetLabelOffset(999999);
    d1hist->GetXaxis()->SetLabelSize(0);
    d1hist->GetXaxis()->SetTitleOffset(999999);
    d1hist->GetXaxis()->SetTitleSize(0);

    for(int iBin = 1; iBin <= erhist->GetNbinsX(); iBin++){
      if(erhist->GetBinError(iBin) > erhist->GetBinError(iBin+1) && iBin != erhist->GetNbinsX())
	erhist->SetBinError(iBin,erhist->GetBinError(iBin+1)*0.9);
    }
    
    d1hist->GetYaxis()->SetTitleOffset(0.3);
    d1hist->GetYaxis()->SetLabelSize(0.12);
    d1hist->GetYaxis()->SetTitleSize(0.15);
    d1hist->GetYaxis()->SetTitle("Data/Pred.");
  
    d1hist->Draw("PE1 SAME");    
    d2hist->Draw("PE1 SAME");
    erhist->Draw("E2 SAME");
    d1hist->Draw("PE SAME");
    d2hist->Draw("PE SAME");

    TH1* unhist = (TH1*)pohist.at(ihist)->Clone(Form("unhist_bin_%d",int(ihist)));
    
    for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinContent(i, 1);
    for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinError(i, 0);
    unhist->SetMarkerSize(0);
    unhist->SetLineColor(kBlack);
    unhist->SetLineStyle(2);
    unhist->SetLineWidth(2);
    unhist->SetFillColor(0);
    unhist->Draw("hist same");  
    pad2->RedrawAxis("G sameaxis");
  
    canvas->SaveAs(Form("prepostfit_zee_bin_%d.pdf",int(ihist)));
    canvas->SaveAs(Form("prepostfit_zee_bin_%d.png",int(ihist)));
  }
}

