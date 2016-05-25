#ifndef DRAWPLOT_H
#define DRAWPLOT_H

void drawPlot(vector<TH1*> hists, double xmin, double ymin, double xmax, double ymax, const char* cname = "c", const char* xlabel = "", const char* ylabel = "", bool drawerrband = false) {

    if (hists.size() < 2) return;
    for (size_t i = 0; i < hists.size(); i++) {
        if (hists[i] == NULL) return;
    }

    THStack* stack = new THStack("stack", "stack");
    for (size_t i = 1; i < hists.size(); i++) {
        stack->Add(hists[i]);
    }
    
    //TLegend* leg = new TLegend(0.6, 0.55, 0.9, 0.9);
    TLegend* leg = new TLegend(0.6, 0.9-0.05*hists.size(), 0.9, 0.9);
    leg->SetFillColor(0);
    leg->AddEntry(hists[0], "Data");
    for (size_t i = 1; i < hists.size(); i++) {
        leg->AddEntry(hists[i], hists[i]->GetTitle(), "F");
    }
        
    TCanvas* canvas = new TCanvas(cname, cname, 600, 700);
    string spad1 = cname;
    string spad2 = cname;
    spad1 += "_pad1";
    spad2 += "_pad2";
    TPad *pad1 = new TPad(spad1.c_str(), spad1.c_str(), 0, 0.3, 1, 1.0);
    canvas->cd();
    TPad *pad2 = new TPad(spad2.c_str(), spad2.c_str(), 0, 0.1, 1, 0.3);

    TH1* frame = canvas->DrawFrame(xmin, ymin, xmax, ymax, "");
    frame->GetXaxis()->SetTitle(xlabel);
    frame->GetYaxis()->SetTitle(ylabel);
    frame->GetYaxis()->CenterTitle();
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetLabelSize(0.9*frame->GetYaxis()->GetLabelSize());
    frame->GetXaxis()->SetTitleSize(0.8*frame->GetXaxis()->GetTitleSize());

    pad1->SetRightMargin(0.075);
    pad1->SetTopMargin(0.06);
    pad1->SetBottomMargin(0.0);
    pad1->Draw();
    pad1->cd();
    frame->Draw();
    CMS_lumi(pad1, 4, 0);
    stack->Draw("HIST SAME");
    hists[0]->Draw("PE SAME");
    leg->Draw("SAME");

    pad1->RedrawAxis();
    pad1->SetLogy();
  
    canvas->cd();
    pad2->SetTopMargin(0.08);
    pad2->SetRightMargin(0.075);
    pad2->SetGridy();
    pad2->Draw();
    pad2->cd();

    TH1* dahist = (TH1*)hists[0]->Clone("dahist");
    TH1* mchist = (TH1*)hists[1]->Clone("mchist");
    for (size_t i = 1; i < hists.size(); i++) {
        mchist->Add(hists[i]);
    }
    if (drawerrband) {
        for (int i = 1; i <= mchist->GetNbinsX(); i++) mchist->SetBinError(i, 0.);
    }
    dahist->Divide(mchist);

    dahist->GetXaxis()->SetLabelSize(3.0*dahist->GetXaxis()->GetLabelSize());
    dahist->GetYaxis()->SetLabelSize(3.0*dahist->GetYaxis()->GetLabelSize());
    dahist->GetYaxis()->SetRangeUser(0., 2.0);
    dahist->GetYaxis()->SetNdivisions(504, false);
    dahist->SetMarkerSize(0);
    dahist->SetLineWidth(2);
    dahist->Draw("PE");
    dahist->GetYaxis()->SetTitleOffset(0.30);
    dahist->GetYaxis()->SetTitleSize(3.5*frame->GetYaxis()->GetTitleSize());
    dahist->GetYaxis()->SetTitle("Data/MC");

    if (drawerrband) {
        TH1* erhist = (TH1*)hists[1]->Clone("erhist");
        for (size_t i = 1; i < hists.size(); i++) {
            erhist->Add(hists[i]);
        }
        erhist->Divide(mchist);
        erhist->SetLineColor(0);
        erhist->SetMarkerColor(0);
        erhist->SetMarkerSize(0);
        erhist->SetFillStyle(3001);
        erhist->SetFillColor(kBlack);
        erhist->Draw("E2 SAME");
    }

    pad1->cd();
    pad1->Draw();
    pad1->RedrawAxis();

    canvas->Print((string(cname)+".pdf").c_str());
}

#endif
