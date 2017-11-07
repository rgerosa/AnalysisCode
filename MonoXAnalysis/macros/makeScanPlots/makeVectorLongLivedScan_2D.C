#include "../CMS_lumi.h"

int code(double mh){
  string mh_str = to_string(mh);
  return atoi(&mh_str[0]);
}

int mmed(double mh){
  string mh_str = to_string(mh);
  stringstream med_str;
  med_str << mh_str[1] << mh_str[2] << mh_str[3] << mh_str[4];
  return atoi((med_str.str()).c_str());
}

int mdm(double mh){
  string mh_str = to_string(mh);
  stringstream mdm_str;
  mdm_str << mh_str[5] << mh_str[6] << mh_str[7] << mh_str[8];
  return atoi(mdm_str.str().c_str());
}

int decaytime(double mh){
  string mh_str = to_string(mh);
  stringstream ctau_str;
  ctau_str << mh_str[9] << mh_str[10] << mh_str[11] << mh_str[12] << mh_str[13] << mh_str[14];
  return atoi(ctau_str.str().c_str());

}

TGraph* produceContour (const int & reduction){

  TObjArray *lContoursE = (TObjArray*) gROOT->GetListOfSpecials()->FindObject("contours");
  std::vector<double> lXE;
  std::vector<double> lYE;
  int lTotalContsE = lContoursE->GetSize();
  for(int i0 = 0; i0 < lTotalContsE; i0++){
    TList * pContLevel = (TList*)lContoursE->At(i0);
    TGraph *pCurv = (TGraph*)pContLevel->First();
    for(int i1 = 0; i1 < pContLevel->GetSize(); i1++){
      for(int i2  = 0; i2 < pCurv->GetN(); i2++) {
        if(i2%reduction != 0) continue; // reduce number of points                                                                                                                                     
        lXE.push_back(pCurv->GetX()[i2]);
        lYE.push_back(pCurv->GetY()[i2]);
      }
      pCurv->SetLineColor(kRed);
      pCurv = (TGraph*)pContLevel->After(pCurv);
    }
  }
  if(lXE.size() == 0) {
    lXE.push_back(0);
    lYE.push_back(0);
  }

  TGraph *lTotalE = new TGraph(lXE.size(),&lXE[0],&lYE[0]);
  return lTotalE;
}

static float nbinsX = 500;
static float nbinsY = 500;
static float minX = 10;
static float minY = 1;
static float maxX = 501;
static float maxY = 100000;
static float minZ = 0.5;
static float maxZ = 30;
static int reductionForContour = 20;

void makeVectorLongLivedScan_2D(string inputFileName, string outputDIR){

  int ncontours = 999;
  gStyle->SetNumberContours(ncontours);
  
  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TTree* tree = (TTree*) inputFile->Get("limit");
  
  TTreeReader reader(tree);
  TTreeReaderValue<double> mh (reader,"mh");
  TTreeReaderValue<float> quantile (reader,"quantileExpected");
  TTreeReaderValue<double> limit (reader,"limit");

  TGraph2D* grexp = new TGraph2D();
  TGraph2D* grexp_up   = new TGraph2D();
  TGraph2D* grexp_down = new TGraph2D();
  TGraph2D* grobs = new TGraph2D();

  int expcounter = 0;
  int obscounter = 0;
  int exp_down_counter = 0;
  int exp_up_counter = 0;

  int medmass = 0 ;
  float minmass = 1000;

  while(reader.Next()){
    
    medmass = mmed(*mh);
    int dmmass  = mdm(*mh);
    int ctau    = decaytime(*mh);

    if (*quantile == 0.5) { // expected limit                                                                                                                                                          
      grexp->SetPoint(expcounter, double(dmmass), double(ctau),*limit);
      expcounter++;
      if(dmmass < minmass) minmass =  dmmass;
    }
    
    if (*quantile < 0.17 && *quantile > 0.14 ) {
      grexp_down->SetPoint(exp_down_counter, double(dmmass), double(ctau),*limit);
      exp_down_counter++;
    }

    if (*quantile < 0.85 && *quantile > 0.83 ) {
      grexp_up->SetPoint(exp_up_counter, double(dmmass), double(ctau),*limit);
      exp_up_counter++;
    }

    if (*quantile == -1) { // observed                                                                                                                                                              
      grobs->SetPoint(obscounter, double(dmmass), double(ctau),*limit);
      obscounter++;
    }
    
  }

  TH2D* hexp       = new TH2D("hexp", "",      nbinsX, minX, maxX, nbinsY, minY, maxY);
  TH2D* hexp_up    = new TH2D("hexp_up", "",   nbinsX, minX, maxX, nbinsY, minY, maxY);
  TH2D* hexp_down  = new TH2D("hexp_down", "", nbinsX, minX, maxX, nbinsY, minY, maxY);
  TH2D* hobs = new TH2D("hobs", "", nbinsX, minX, maxX, nbinsY, minY, maxY);

  // make granularity                                                                                                                                                                                
  for (int i = 1; i   <= nbinsX; i++) {
    for (int j = 1; j <= nbinsY; j++) {
      hexp_up->SetBinContent(i,j,   grexp_up->Interpolate(hexp_up->GetXaxis()->GetBinCenter(i),hexp_up->GetYaxis()->GetBinCenter(j)));
      hexp_down->SetBinContent(i,j, grexp_down->Interpolate(hexp_down->GetXaxis()->GetBinCenter(i),hexp_down->GetYaxis()->GetBinCenter(j)));
      hexp->SetBinContent(i,j,grexp->Interpolate(hexp->GetXaxis()->GetBinCenter(i),hexp->GetYaxis()->GetBinCenter(j)));
      hobs->SetBinContent(i,j,grobs->Interpolate(hobs->GetXaxis()->GetBinCenter(i),hobs->GetYaxis()->GetBinCenter(j)));
    }
  }
  

  // fix mass points below min med mass generated                                                                                                                                                      
  for (int i = 1; i   <= nbinsX; i++) {
    for (int j = 1; j <= nbinsY; j++) {
      if(hexp->GetXaxis()->GetBinCenter(i) < minmass){
        hexp_up->SetBinContent(i,j,hexp_up->GetBinContent(hexp_up->FindBin(minmass,hexp_up->GetYaxis()->GetBinCenter(j))));
        hexp_down->SetBinContent(i,j,hexp_down->GetBinContent(hexp_down->FindBin(minmass,hexp_down->GetYaxis()->GetBinCenter(j))));
        hexp->SetBinContent(i,j,hexp->GetBinContent(hexp->FindBin(minmass,hexp->GetYaxis()->GetBinCenter(j))));
      }
      if(hexp->GetXaxis()->GetBinCenter(i) < minmass){
        hobs->SetBinContent(i,j,hobs->GetBinContent(hobs->FindBin(minmass,hobs->GetYaxis()->GetBinCenter(j))));
      }
    }
  }

  ////////////                                                                                                                                                                                       
  for(int i = 1; i <= nbinsX; i++){
    for(int j = 1; j <= nbinsY; j++){

      if(hexp -> GetBinContent(i,j) <= 0) hexp->SetBinContent(i,j,maxZ);
      if(hexp_down -> GetBinContent(i,j) <= 0) hexp_down->SetBinContent(i,j,maxZ);
      if(hexp_up -> GetBinContent(i,j) <= 0) hexp_up->SetBinContent(i,j,maxZ);
      if(hobs -> GetBinContent(i,j) <= 0) hobs->SetBinContent(i,j,maxZ);

      if(hexp -> GetBinContent(i,j) > maxZ) hexp->SetBinContent(i,j,maxZ);
      if(hexp_down -> GetBinContent(i,j) > maxZ) hexp_down->SetBinContent(i,j,maxZ);
      if(hexp_up -> GetBinContent(i,j) > maxZ) hexp_up->SetBinContent(i,j,maxZ);
      if(hobs -> GetBinContent(i,j) > maxZ) hobs->SetBinContent(i,j,maxZ);

      if(hexp -> GetBinContent(i,j) < minZ) hexp->SetBinContent(i,j,minZ);
      if(hexp_down -> GetBinContent(i,j) < minZ) hexp_down->SetBinContent(i,j,minZ);
      if(hexp_up -> GetBinContent(i,j) < minZ) hexp_up->SetBinContent(i,j,minZ);
      if(hobs -> GetBinContent(i,j) < minZ) hobs->SetBinContent(i,j,minZ);
    }
  }

  hexp->Smooth();
  hexp_down->Smooth();
  hexp_up->Smooth();
  hobs->Smooth();


  ////////////////                                                                                                                                                                                     
  TH2* hexp2 = (TH2*) hexp->Clone("hexp2");
  TH2* hexp2_up = (TH2*) hexp_up->Clone("hexp2_up");
  TH2* hexp2_down = (TH2*) hexp_down->Clone("hexp2_down");
  TH2* hobs2 = (TH2*) hobs->Clone("hobs2");

  //////////                                                                                                                                                                                           
  double contours[1]; contours[0]=1,
  hexp2->SetContour(1,contours);
  hexp2_up->SetContour(1,contours);
  hexp2_down->SetContour(1,contours);
  hobs2->SetContour(1,contours);

  // All the plotting and cosmetics                                                                                                                                                                    
  TCanvas* canvas = new TCanvas("canvas", "canvas",625,600);
  canvas->SetRightMargin(0.15);
  canvas->SetLeftMargin(0.13);
  canvas->SetLogz();

  TH1* frame = canvas->DrawFrame(minX,minY,maxX,maxY, "");
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle("m_{DM} [GeV]");
  frame->GetYaxis()->SetTitle("c#tau [mm]");
  frame->GetXaxis()->SetTitleOffset(1.15);
  frame->GetYaxis()->SetTitleOffset(1.20);
  frame->Draw();

  hobs->GetZaxis()->SetRangeUser(minZ,maxZ);
  hexp->GetZaxis()->SetRangeUser(minZ,maxZ);
  hexp_up->GetZaxis()->SetRangeUser(minZ,maxZ);
  hexp_down->GetZaxis()->SetRangeUser(minZ,maxZ);

  hexp2_up->GetZaxis()->SetLabelSize(0);
  hexp2_up->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp_up = produceContour(reductionForContour);

  hexp2_down->GetZaxis()->SetLabelSize(0);
  hexp2_down->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp_dw = produceContour(reductionForContour);

  hexp2->GetZaxis()->SetLabelSize(0);
  hexp2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_exp = produceContour(reductionForContour);

  hobs2->GetZaxis()->SetLabelSize(0);
  hobs2->Draw("contz list same");
  canvas->Update();
  TGraph* contour_obs = produceContour(reductionForContour);

  frame->Draw();
  hobs->Draw("COLZ SAME");

  contour_exp_up->SetLineColor(kBlack);
  contour_exp->SetLineColor(kBlack);
  contour_exp_dw->SetLineColor(kBlack);
  contour_exp_up->SetLineWidth(1);
  contour_exp->SetLineWidth(3);
  contour_exp_dw->SetLineWidth(1);
  contour_exp_up->SetLineStyle(7);
  contour_exp->SetLineStyle(7);
  contour_exp_dw->SetLineStyle(7);
  //contour_exp_up->Draw("Lsame");
  //contour_exp_dw->Draw("Lsame");
  contour_exp->Draw("Lsame");

  contour_obs->SetLineColor(kRed);
  contour_obs->SetLineWidth(3);
  //contour_obs->Draw("Lsame");

  CMS_lumi(canvas,"35.9",false,false,false,0,-0.09);
 
  TLegend *leg = new TLegend(0.50,0.20,0.85,0.50);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextFont(42);
  leg->AddEntry(contour_exp,"Median expected 95% CL","L");
  //leg->AddEntry(contour_exp_up,"68% expected","L");
  //leg->AddEntry(contour_obs,"Observed 95% CL","L");
  leg->Draw("SAME");

  TLatex * tex = new TLatex();
  tex->SetNDC();
  tex->SetTextFont(42);
  tex->SetLineWidth(2);
  tex->SetTextSize(0.030);
  tex->DrawLatex(0.175,0.80,("#bf{Vector med, Dirac DM, m_{med} = "+to_string(medmass)+" GeV}").c_str());

  TLatex *   tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetLineWidth(2);
  tex2->SetTextSize(0.042);
  tex2->SetTextAngle(90);
  tex2->DrawLatex(0.975,0.55,"Observed #sigma_{95% CL}/#sigma_{th}");

  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();

  canvas->SaveAs((outputDIR+"/scan_vector_mdm_ctau_mmed_"+to_string(medmass)+".pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_ctau_mmed_"+to_string(medmass)+".png").c_str());
}
