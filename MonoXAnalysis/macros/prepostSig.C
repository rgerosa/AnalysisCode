void plotPostFitSig2(string fitFilename, string templateFileName, string observable, bool blind) {

    TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
    TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1.0);
    canvas->cd();
    TPad *pad2 = new TPad("pad2","pad2",0,0.1,1,0.3);

    TFile* pfile = new TFile(fitFilename.c_str());
    TFile* dfile = new TFile(templateFileName._str());

    TH1* dthist = NULL;
    if(not blind)
      dthist = (TH1*)dfile->Get(("datahist_"+observable).c_str());
    else
      dthist = (TH1*)pfile->Get("shapes_fit_b/ch1/total_background");

    TH1* mjhist = (TH1*)dfile->Get(("monoJhist_"+observable).c_str());
    TH1* mwhist = (TH1*)dfile->Get(("monoWhist_"+observable).c_str());
    TH1* mzhist = (TH1*)dfile->Get(("monoZhist_"+observable).c_str());

    /*
    TH1* znhist = (TH1*)pfile->Get("shapes_prefit/ch1/Znunu");    
    TH1* zlhist = (TH1*)pfile->Get("shapes_prefit/ch1/ZJets");    
    TH1* wlhist = (TH1*)pfile->Get("shapes_prefit/ch1/WJets");    
    TH1* tthist = (TH1*)pfile->Get("shapes_prefit/ch1/Top");    
    TH1* dihist = (TH1*)pfile->Get("shapes_prefit/ch1/Dibosons");    
    TH1* qchist = (TH1*)pfile->Get("shapes_prefit/ch1/QCD");    
    TH1* tohist = (TH1*)pfile->Get("shapes_prefit/ch1/total_background");    
    */

    TH1* znhist = (TH1*)pfile->Get("shapes_fit_b/ch1/Znunu");    
    TH1* zlhist = (TH1*)pfile->Get("shapes_fit_b/ch1/ZJets");    
    TH1* wlhist = (TH1*)pfile->Get("shapes_fit_b/ch1/WJets");    
    TH1* tthist = (TH1*)pfile->Get("shapes_fit_b/ch1/Top");    
    TH1* dihist = (TH1*)pfile->Get("shapes_fit_b/ch1/Dibosons");    
    TH1* qchist = (TH1*)pfile->Get("shapes_fit_b/ch1/QCD");    
    TH1* tohist = (TH1*)pfile->Get("shapes_fit_b/ch1/total_background");    
    TH1* tphist = (TH1*)pfile->Get("shapes_prefit/ch1/total_background");    

    for (int i = 0; i <= dthist->GetNbinsX(); i++) {
        double yield = 0.0;
        yield += zlhist->GetBinContent(i);
        yield += wlhist->GetBinContent(i);
        yield += tthist->GetBinContent(i);
        yield += dihist->GetBinContent(i);
        yield += qchist->GetBinContent(i);
        yield += znhist->GetBinContent(i);
        dthist->SetBinContent(i, yield);
        dthist->SetBinError(i, 0.);
    }

    mjhist->Scale(1.0, "width");
    mwhist->Scale(1.0, "width");
    mzhist->Scale(1.0, "width");

    TH1* srhist = (TH1*)tohist->Clone("srhist");
    for (int i = 0; i <= srhist->GetNbinsX(); i++) {
        srhist->SetBinContent(i, (srhist->GetBinContent(i)+sihist->GetBinContent(i))/srhist->GetBinContent(i));
        srhist->SetBinError(i, 0);
    }

    mjhist->SetLineColor(kBlack);
    mjhist->SetLineStyle(2);
    mjhist->SetLineWidth(4);
    mjhist->SetMarkerSize(0);

    mwhist->SetLineColor(kBlack);
    mwhist->SetLineStyle(3);
    mwhist->SetLineWidth(4);
    mwhist->SetMarkerSize(0);

    mzhist->SetLineColor(kBlack);
    mzhist->SetLineStyle(4);
    mzhist->SetLineWidth(4);
    mzhist->SetMarkerSize(0);

    srhist->SetLineStyle(2);
    srhist->SetLineWidth(2);
    srhist->SetMarkerSize(0);

    znhist->SetFillColor(kGreen+1);
    zlhist->SetFillColor(kCyan);
    wlhist->SetFillColor(kRed);
    tthist->SetFillColor(kViolet);
    dihist->SetFillColor(kBlue);
    qchist->SetFillColor(kYellow);

    THStack* stack = new THStack("stack", "stack");
    stack->Add(qchist);
    stack->Add(dihist);
    stack->Add(tthist);
    stack->Add(zlhist);
    stack->Add(wlhist);
    stack->Add(znhist);

    TH1* frame = canvas->DrawFrame(200., 1.5e-4, 1200., 80000.0, "");
    frame->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]    ");
    frame->GetYaxis()->SetTitle("Events / GeV");
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetLabelSize(0.9*frame->GetYaxis()->GetLabelSize());
    frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());

    pad1->SetRightMargin(0.075);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();
    pad1->cd();

    frame ->Draw();
    CMS_lumi(pad1, 4, 0);
    stack ->Draw("HIST SAME");
    dthist->Draw("P SAME");
    mjhist->Draw("HIST SAME");
    mwhist->Draw("HIST SAME");
    mzhist->Draw("HIST SAME");

    TH1* dahist = (TH1*)dthist->Clone("dahist");
    TH1* dphist = (TH1*)dthist->Clone("dahist");
    dphist->SetLineColor(kRed);
    dphist->SetMarkerColor(kRed);
    dphist->SetMarkerSize(0.7);
    dahist->SetLineColor(kBlue);
    dahist->SetMarkerColor(kBlue);
    dahist->SetMarkerSize(0.7);

    TLegend* leg = new TLegend(0.58, 0.42, 0.9, 0.92);
    leg->SetFillColor(0);
    leg->AddEntry(dthist, "Data");
    leg->AddEntry(mjhist, "Mono-J (V, 1TeV)");
    leg->AddEntry(mwhist, "Mono-W (V, 1TeV)");
    leg->AddEntry(mzhist, "Mono-Z (V, 1TeV)");
    leg->AddEntry(znhist, "Z(#nu#nu)", "F");
    leg->AddEntry(wlhist, "W(l#nu)", "F");
    leg->AddEntry(zlhist, "Z(ll)", "F");
    leg->AddEntry(tthist, "Top", "F");
    leg->AddEntry(dihist, "Dibosons", "F");
    leg->AddEntry(qchist, "QCD", "F");
    leg->Draw("SAME");

    pad1->RedrawAxis();
    pad1->SetLogy();

    canvas->cd();
    pad2->SetTopMargin(0.08);
    pad2->SetRightMargin(0.075);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    TH1* mphist = (TH1*)tphist->Clone("mphist");
    TH1* mchist = (TH1*)zlhist->Clone("mchist");
    TH1* unhist = (TH1*)zlhist->Clone("unhist");
    mchist->Add(wlhist);
    mchist->Add(tthist);
    mchist->Add(dihist);
    mchist->Add(qchist);
    mchist->Add(znhist);
    for (int i = 1; i <= mchist->GetNbinsX(); i++) mchist->SetBinError(i, 0);
    for (int i = 1; i <= mphist->GetNbinsX(); i++) mphist->SetBinError(i, 0);
    dahist->Divide(mchist);
    dphist->Divide(mphist);
    tohist->Divide(mchist);
    tohist->SetLineColor(0);
    tohist->SetMarkerColor(0);
    tohist->SetMarkerSize(0);
    tohist->SetFillColor(kGray);
    //tohist->SetFillStyle(3001);

    for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinContent(i, 1);
    for (int i = 1; i <= unhist->GetNbinsX(); i++) unhist->SetBinError(i, 0);
    unhist->SetMarkerSize(0);
    unhist->SetLineColor(kBlack);
    unhist->SetLineStyle(2);
    unhist->SetFillColor(0);

    dahist->GetXaxis()->SetLabelSize(3.0*dahist->GetXaxis()->GetLabelSize());
    dahist->GetYaxis()->SetLabelSize(3.0*dahist->GetYaxis()->GetLabelSize());
    dahist->GetYaxis()->SetRangeUser(0.5, 1.5);
    //dahist->GetYaxis()->SetRangeUser(0.6, 1.4);
    dahist->GetYaxis()->SetNdivisions(504, false);
    //dahist->SetLineWidth(2);
    dahist->Draw("P");
    //srhist->Draw("HIST SAME");
    dahist->Draw("P SAME");
    tohist->Draw("E2 SAME");
    unhist->Draw("SAME");
    dahist->Draw("P SAME");
    //dphist->Draw("P SAME");
    //unhist->Draw("SAME");
    dahist->GetYaxis()->SetTitleOffset(0.30);
    dahist->GetYaxis()->SetTitleSize(3.5*frame->GetYaxis()->GetTitleSize());
    dahist->GetYaxis()->SetTitle("Data/Pred.");

    pad1->cd();
    pad1->Draw();
    pad1->RedrawAxis();
    pad2->RedrawAxis();

    canvas->SaveAs("postfit_sig.pdf");
    canvas->SaveAs("postfit_sig.png");

}

