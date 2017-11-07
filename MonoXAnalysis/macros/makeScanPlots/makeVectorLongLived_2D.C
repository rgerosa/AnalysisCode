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


void makeVectorLongLived_2D(string inputFileName, string outputDIR){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TTree* tree = (TTree*) inputFile->Get("limit");
  
  TTreeReader reader(tree);
  TTreeReaderValue<double> mh (reader,"mh");
  TTreeReaderValue<float> quantile (reader,"quantileExpected");
  TTreeReaderValue<double> limit (reader,"limit");

  vector<double> dmMassPoints;
  vector<double> ctauPoints;
  int medmass = 0;
  while(reader.Next()){
    medmass = mmed(*mh);
    int dmmass  = mdm(*mh);
    int ctau    = decaytime(*mh);
    if(*quantile != 0.5) continue;

    if(std::find(dmMassPoints.begin(),dmMassPoints.end(),dmmass) == dmMassPoints.end())
      dmMassPoints.push_back(float(dmmass));
    if(std::find(ctauPoints.begin(),ctauPoints.end(),ctau) == ctauPoints.end())
       ctauPoints.push_back(float(ctau));
  }
  
  sort(dmMassPoints.begin(),dmMassPoints.end());
  sort(ctauPoints.begin(),ctauPoints.end());

  vector<float> dmMassBinning;
  dmMassBinning.push_back(0.);
  for(int mdm  = 0; mdm < dmMassPoints.size()-1; mdm++){
    dmMassBinning.push_back(dmMassPoints.at(mdm)+fabs(dmMassPoints.at(mdm)-dmMassPoints.at(mdm+1))/3);
    dmMassBinning.push_back(dmMassPoints.at(mdm+1)-fabs(dmMassPoints.at(mdm)-dmMassPoints.at(mdm+1))/3);
  }
  dmMassBinning.push_back(dmMassPoints.back()+(dmMassPoints.back()-dmMassPoints.at(dmMassPoints.size()-2)));


  vector<float> ctauBinning;
  ctauBinning.push_back(ctauPoints.front()/10);
  for(int mdm  = 0; mdm < ctauPoints.size()-1; mdm++){
    ctauBinning.push_back(ctauPoints.at(mdm)+fabs(ctauPoints.at(mdm)-ctauPoints.at(mdm+1))/3);
    ctauBinning.push_back(ctauPoints.at(mdm+1)-fabs(ctauPoints.at(mdm)-ctauPoints.at(mdm+1))/3);
  }
  ctauBinning.push_back(ctauPoints.back()+(ctauPoints.back()-ctauPoints.at(ctauPoints.size()-2)));

  TH2F* hexp = new TH2F("hexp","",dmMassBinning.size()-1,&dmMassBinning[0],ctauBinning.size()-1,&ctauBinning[0]);
  TH2F* hobs = new TH2F("hobs","",dmMassBinning.size()-1,&dmMassBinning[0],ctauBinning.size()-1,&ctauBinning[0]);

  int expcounter = 0;
  int obscounter = 0;

  reader.SetEntry(0);
  while(reader.Next()){
    
    medmass = mmed(*mh);
    int dmmass  = mdm(*mh);
    int ctau    = decaytime(*mh); 
    if (*quantile == 0.5) {
      hexp->SetBinContent(hexp->GetXaxis()->FindBin(dmmass),hexp->GetYaxis()->FindBin(ctau),*limit);
      expcounter++;
    }
    if (*quantile == -1) {
      hobs->SetBinContent(hobs->GetXaxis()->FindBin(dmmass),hobs->GetYaxis()->FindBin(ctau),*limit);
      obscounter++;
    }       
  }

  // All the plotting and cosmetics                                                                                                                                                                    
  TCanvas* canvas = new TCanvas("canvas", "canvas",625,600);
  canvas->SetRightMargin(0.15);
  canvas->SetLeftMargin(0.13);
  canvas->SetLogz();

  hexp->GetYaxis()->CenterTitle();
  hexp->GetXaxis()->SetTitle("m_{DM} [GeV]");
  hexp->GetYaxis()->SetTitle("c#tau [mm]");
  hexp->GetXaxis()->SetTitleOffset(1.15);
  hexp->GetYaxis()->SetTitleOffset(1.20);
  hexp->Draw("colz text");
  
  TLatex *   tex2 = new TLatex();
  tex2->SetNDC();
  tex2->SetTextFont(42);
  tex2->SetLineWidth(2);
  tex2->SetTextSize(0.042);
  tex2->SetTextAngle(90);
  tex2->DrawLatex(0.975,0.55,"Expected #sigma_{95% CL}/#sigma_{th}");


  CMS_lumi(canvas,"35.9",true,true,false,0,-0.09);
  
  canvas->RedrawAxis("sameaxis");
  canvas->SetLogy();

  canvas->SaveAs((outputDIR+"/scan_vector_mdm_ctau_mmed_"+to_string(medmass)+"_expected.pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_ctau_mmed_"+to_string(medmass)+"_expected.png").c_str());

  hobs->GetYaxis()->CenterTitle();
  hobs->GetXaxis()->SetTitle("m_{DM} [GeV]");
  hobs->GetYaxis()->SetTitle("c#tau [mm]");
  hobs->GetXaxis()->SetTitleOffset(1.15);
  hobs->GetYaxis()->SetTitleOffset(1.20);
  hobs->Draw("colz text");

  canvas->RedrawAxis("sameaxis");
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_ctau_mmed_"+to_string(medmass)+"_observed.pdf").c_str());
  canvas->SaveAs((outputDIR+"/scan_vector_mdm_ctau_mmed_"+to_string(medmass)+"_observed.png").c_str());
}
