#include "../CMS_lumi.h"
#include "../makeTemplates/histoUtils.h"

void drawUpperPlot(TPad* pad1, TH1* histo_1, TH1* histo_2, string label1, string label2){

  pad1->cd();
  TH1* frame = (TH1*) histo_1->Clone("htemp");
  frame->Reset("ICES");
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetLabelSize(0.);
  frame->GetXaxis()->SetLabelOffset(1.10);
  frame->GetXaxis()->SetTitleSize(0.);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleSize(0.042);
  frame->GetYaxis()->SetTitleOffset(1.35);
  frame->GetYaxis()->SetTitle("Events");
  if(histo_1->GetMinimum() != 0)
    frame->GetYaxis()->SetRangeUser(histo_1->GetMinimum()*0.5,histo_1->GetMaximum()*100);
  else 
    frame->GetYaxis()->SetRangeUser(0.001,histo_1->GetMaximum()*100);

  frame ->Draw();  

  CMS_lumi(pad1,"2.3");
  
  histo_1->SetLineColor(kBlack);
  histo_1->SetLineWidth(2);  
  TH1* histo_1_temp = (TH1*) histo_1->Clone("histo_1_temp");
  histo_1_temp->SetFillColor(kGray);

  histo_2->SetLineColor(kRed);
  histo_2->SetMarkerColor(kRed);
  histo_2->SetMarkerStyle(20);
  histo_2->SetMarkerSize(0.8);

  histo_1_temp->Draw("E2 SAME");
  histo_1->Draw("HIST SAME");
  histo_2->Draw("EP SAME");

  TLegend* leg = new TLegend(0.6,0.72,0.93,0.85);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(histo_1,label1.c_str(),"L");
  leg->AddEntry(histo_2,label2.c_str(),"PE");
  leg->AddEntry((TObject*)0,Form("Ratio norm = %f",histo_1->Integral()/histo_2->Integral()));
  leg->Draw("same");
  
  return ;
}

void drawDownPlot(TH1* histo_1, TH1* histo_2,string xAxisTitle){

  TPad *pad2 = new TPad("pad2","pad2",0,0.,1,0.9);
  pad2->SetTopMargin(0.7);
  pad2->SetRightMargin(0.06);
  pad2->SetFillColor(0);
  pad2->SetGridy(1);
  pad2->SetFillStyle(0);
  pad2->Draw();
  pad2->cd();

  TH1* frame = (TH1*) histo_1->Clone("htemp");
  frame->Reset("ICES");
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetLabelSize(0.035);
  frame->GetXaxis()->SetTitleSize(0.042);
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetLabelSize(0.035);
  frame->GetYaxis()->SetTitleSize(0.042);
  frame->GetYaxis()->SetTitle("Ratio");
  frame->GetXaxis()->SetTitle(xAxisTitle.c_str());
  frame->GetYaxis()->SetNdivisions(504);
  frame->GetYaxis()->SetRangeUser(0.75,1.25);
  frame->Draw();

  TH1* ratio = (TH1*) histo_1->Clone("ratio");
  TH1* histo_2_temp    = (TH1*) histo_2->Clone("histo_2_temp");
  TH1* histo_2_temp_v2 = (TH1*) histo_2->Clone("histo_2_temp_v2");
  for(int iBin = 0; iBin < histo_2_temp->GetNbinsX(); iBin++){    
    histo_2_temp->SetBinError(iBin+1,0.);
  }
  ratio->Divide(histo_2_temp);


  
  histo_2_temp_v2->Divide(histo_2_temp);
  histo_2_temp_v2->SetLineWidth(3);
  histo_2_temp_v2->SetLineColor(kRed);
  histo_2_temp_v2->SetFillColor(kGray);
  histo_2_temp_v2->SetMarkerSize(0);
  histo_2_temp_v2->Draw("E2 same");
  histo_2_temp->Divide(histo_2);
  histo_2_temp->SetLineWidth(2);
  histo_2_temp->SetLineColor(kRed);
  for(int iBin = 0; iBin < histo_2_temp->GetNbinsX(); iBin++) histo_2_temp->SetBinContent(iBin+1,1.);
  histo_2_temp->Draw("hist same");

  ratio->SetLineColor(kBlack);
  ratio->SetMarkerColor(kBlack);
  ratio->SetMarkerStyle(20);
  ratio->SetMarkerSize(1);
  ratio->SetLineWidth(1);

  ratio->Draw("PE1same");
  pad2->RedrawAxis("sameaxis");

  return ;
}

void makeTemplateComparison_withMIT( string templateFile_1, // template file 1
				     string templateFile_2, // template file 2
				     string controlRegion, // control region name
				     Category  category,  // category
				     string observable, 
				     bool   isHiggsInvisible,
				     string templateFile_1_Label, string templateFile_2_Label, string xAxisTitle){

  gROOT->SetBatch(kTRUE);
  gROOT->ForceStyle(kTRUE);
  setTDRStyle();
  gStyle->SetOptStat(0);

  TFile* template1 = TFile::Open(templateFile_1.c_str());
  TFile* template2 = TFile::Open(templateFile_2.c_str());

  TH1* data_1 = NULL;
  TH1* data_2 = NULL;

  string dir;
  if(category == Category::monojet)
    dir = "category_monojet";
  else if(category == Category::monoV)
    dir = "category_monov";

  string cat;
  if(controlRegion == "zmm")
    cat = "Zmm";
  else if(controlRegion == "zee")
    cat = "Zee";
  else if(controlRegion == "wen")
    cat = "Wen";
  else if(controlRegion == "wmn")
    cat = "Wmn";
  else if(controlRegion == "gam")
    cat = "gjets";
  else if(controlRegion == "SR")
    cat = "signal";

  if(controlRegion == "SR"){
    data_1 = (TH1*) template1->FindObjectAny(("datahist_"+observable).c_str());
    data_2 = (TH1*) template2->Get((dir+"/"+cat+"_data").c_str());
  }
  else{
    data_1 = (TH1*) template1->FindObjectAny(("datahist"+controlRegion+"_"+observable).c_str());
    data_2 = (TH1*) template2->Get((dir+"/"+cat+"_data").c_str());
  }

  TH1* zinv_1 = NULL; 
  TH1* zinv_2 = NULL;
  if(controlRegion == "SR"){
    zinv_1 = (TH1*) template1->FindObjectAny(("zinvhist_"+observable).c_str());
    zinv_2 = (TH1*) template2->Get((dir+"/"+cat+"_zjets").c_str());
  }

  TH1* wjet_1 = NULL;
  TH1* wjet_2 = NULL;
  if(controlRegion == "SR"){
    wjet_1 = (TH1*) template1->FindObjectAny(("wjethist_"+observable).c_str());
    wjet_2 = (TH1*) template2->Get((dir+"/"+cat+"_wjets").c_str());
  }
  else{
    wjet_1 = (TH1*) template1->FindObjectAny(("vlbkghist"+controlRegion+"_"+observable).c_str());
    wjet_2 = (TH1*) template2->Get((dir+"/"+cat+"_wjets").c_str());
  }

  TH1* zjet_1 = NULL;
  TH1* zjet_2 = NULL;
  if(controlRegion == "SR"){
    zjet_1 = (TH1*) template1->FindObjectAny(("zjethist_"+observable).c_str());
    zjet_2 = (TH1*) template2->Get((dir+"/"+cat+"_zll").c_str());
  }
  else{
    zjet_1 = (TH1*) template1->FindObjectAny(("vllbkghist"+controlRegion+"_"+observable).c_str());
    zjet_2 = (TH1*) template2->Get((dir+"/"+cat+"_zll").c_str());
  }

  TH1* diboson_1 = NULL;
  TH1* diboson_2 = NULL;
  if(controlRegion == "SR"){
    diboson_1 = (TH1*) template1->FindObjectAny(("dbkghist_"+observable).c_str());
    diboson_2 = (TH1*) template2->Get((dir+"/"+cat+"_diboson").c_str());
  }
  else{
    diboson_1 = (TH1*) template1->FindObjectAny(("dbkghist"+controlRegion+"_"+observable).c_str());
    diboson_2 = (TH1*) template2->Get((dir+"/"+cat+"_diboson").c_str());
  }

  TH1* qcd_1 = NULL;
  TH1* qcd_2 = NULL;
  if(controlRegion == "SR"){
    qcd_1 = (TH1*) template1->FindObjectAny(("qbkghist_"+observable).c_str());
    qcd_2 = (TH1*) template2->Get((dir+"/"+cat+"_qcd").c_str());
  }
  else{
    qcd_1 = (TH1*) template1->FindObjectAny(("qbkghist"+controlRegion+"_"+observable).c_str());
    qcd_2 = (TH1*) template2->Get((dir+"/"+cat+"_qcd").c_str());
  }

  TH1* gamma_1 = NULL;
  TH1* gamma_2 = NULL;
  if(controlRegion == "SR"){
    gamma_1 = (TH1*) template1->FindObjectAny(("gbkghist_"+observable).c_str());
    gamma_2 = (TH1*) template2->Get((dir+"/"+cat+"_gjets").c_str());
  }
  else{
    gamma_1 = (TH1*) template1->FindObjectAny(("gbkghist"+controlRegion+"_"+observable).c_str());
    gamma_2 = (TH1*) template2->Get((dir+"/"+cat+"_gjets").c_str());
  }

  TH1* top_1 = NULL;
  TH1* top_2 = NULL;
  if(controlRegion == "SR"){
    top_1 = (TH1*) template1->FindObjectAny(("tbkghist_"+observable).c_str());
    top_2 = (TH1*) template2->Get((dir+"/"+cat+"_top").c_str());
  }
  else{
    top_1 = (TH1*) template1->FindObjectAny(("tbkghist"+controlRegion+"_"+observable).c_str());
    top_2 = (TH1*) template2->Get((dir+"/"+cat+"_top").c_str());
  }

  TH1* ewk_1 = NULL;
  TH1* ewk_2 = NULL;
  if(controlRegion == "SR"){
    ewk_1 = (TH1*) template1->FindObjectAny(("ewkbkgzhist_"+observable).c_str());
    TH1* temp =  (TH1*) template1->FindObjectAny(("ewkbkgwhist_"+observable).c_str());
    if(temp)
      ewk_1->Add(temp);
    ewk_2 = (TH1*) template2->Get((dir+"/"+cat+"_ewk").c_str());
  }
  else{
    ewk_1 = (TH1*) template1->FindObjectAny(("ewkzbkghist"+controlRegion+"_"+observable).c_str());
    TH1* temp = (TH1*) template1->FindObjectAny(("ewkwbkghist"+controlRegion+"_"+observable).c_str());
    if(temp)
      ewk_1->Add(temp);
    ewk_2 = (TH1*) template2->Get((dir+"/"+cat+"_ewk").c_str());
  }

  // start comparison
  TCanvas* canvas = new TCanvas("canvas", "canvas", 600, 700);
  canvas->SetTickx(1);
  canvas->SetTicky(1);
  canvas->cd();
  canvas->SetBottomMargin(0.3);
  canvas->SetRightMargin(0.06);
  canvas->SetLogy();

  if(data_1 != NULL and data_2 != NULL and data_1 != 0 and data_2 != 0){
    drawUpperPlot(canvas,data_1,data_2,"Z #rightarrow #nu#nu "+ templateFile_1_Label,"Z #rightarrow #nu#nu "+ templateFile_2_Label);
    drawDownPlot(data_1,data_2,xAxisTitle);
    canvas->SaveAs(("Data_"+controlRegion+"_"+observable+".pdf").c_str(),"pdf");
    canvas->SaveAs(("Data_"+controlRegion+"_"+observable+".png").c_str(),"png");
  }

  if(zinv_1 != NULL and zinv_2 != NULL){
    drawUpperPlot(canvas,zinv_1,zinv_2,"Z #rightarrow #nu#nu "+ templateFile_1_Label,"Z #rightarrow #nu#nu "+ templateFile_2_Label);
    drawDownPlot(zinv_1,zinv_2,xAxisTitle);
    canvas->SaveAs(("Znunu_"+controlRegion+"_"+observable+".pdf").c_str(),"pdf");
    canvas->SaveAs(("Znunu_"+controlRegion+"_"+observable+".png").c_str(),"png");
  }

  if(wjet_1 != NULL and wjet_2 != NULL){
    drawUpperPlot(canvas,wjet_1,wjet_2,"W #rightarrow l#nu "+ templateFile_1_Label,"W #rightarrow l#nu "+ templateFile_2_Label);
    drawDownPlot(wjet_1,wjet_2,xAxisTitle);
    canvas->SaveAs(("Wjet_"+controlRegion+"_"+observable+".pdf").c_str(),"pdf");
    canvas->SaveAs(("Wjet_"+controlRegion+"_"+observable+".png").c_str(),"png");
  }
  
  if(zjet_1 != NULL and zjet_2 != NULL){
    drawUpperPlot(canvas,zjet_1,zjet_2,"Z #rightarrow ll "+ templateFile_1_Label,"Z #rightarrow ll "+ templateFile_2_Label);
    drawDownPlot(zjet_1,zjet_2,xAxisTitle);
    canvas->SaveAs(("Zjet_"+controlRegion+"_"+observable+".pdf").c_str(),"pdf");
    canvas->SaveAs(("Zjet_"+controlRegion+"_"+observable+".png").c_str(),"png");
  }

  if(qcd_1 != NULL and qcd_2 != NULL){
    drawUpperPlot(canvas,qcd_1,qcd_2,"QCD "+ templateFile_1_Label,"QCD "+ templateFile_2_Label);
    drawDownPlot(qcd_1,qcd_2,xAxisTitle);
    canvas->SaveAs(("QCD_"+controlRegion+"_"+observable+".pdf").c_str(),"pdf");
    canvas->SaveAs(("QDC_"+controlRegion+"_"+observable+".png").c_str(),"png");
  }

  if(diboson_1 != NULL and diboson_2 != NULL){
    drawUpperPlot(canvas,diboson_1,diboson_2,"DiBoson "+ templateFile_1_Label,"DiBoson "+ templateFile_2_Label);
    drawDownPlot(diboson_1,diboson_2,xAxisTitle);
    canvas->SaveAs(("DiBoson_"+controlRegion+"_"+observable+".pdf").c_str(),"pdf");
    canvas->SaveAs(("DiBoson_"+controlRegion+"_"+observable+".png").c_str(),"png");
  }


  if(gamma_1 != NULL and gamma_2 != NULL){
    drawUpperPlot(canvas,gamma_1,gamma_2,"#gamma+jet "+ templateFile_1_Label,"#gamma+jet "+ templateFile_2_Label);
    drawDownPlot(gamma_1,gamma_2,xAxisTitle);
    canvas->SaveAs(("GJets_"+controlRegion+"_"+observable+".pdf").c_str(),"pdf");
    canvas->SaveAs(("GJets_"+controlRegion+"_"+observable+".png").c_str(),"png");
  }

  if(top_1 != NULL and top_2 != NULL){
    drawUpperPlot(canvas,top_1,top_2,"top "+ templateFile_1_Label,"top "+ templateFile_2_Label);
    drawDownPlot(top_1,top_2,xAxisTitle);
    canvas->SaveAs(("Top_"+controlRegion+"_"+observable+".pdf").c_str(),"pdf");
    canvas->SaveAs(("Top_"+controlRegion+"_"+observable+".png").c_str(),"png");
  }

  if(ewk_1 != NULL and ewk_2 != NULL){
    drawUpperPlot(canvas,ewk_1,ewk_2,"W/Z ewk "+ templateFile_1_Label,"W/Z ewk "+ templateFile_2_Label);
    drawDownPlot(ewk_1,ewk_2,xAxisTitle);
    canvas->SaveAs(("EWK_"+controlRegion+"_"+observable+".pdf").c_str(),"pdf");
    canvas->SaveAs(("EWK_"+controlRegion+"_"+observable+".png").c_str(),"png");
  }
}
