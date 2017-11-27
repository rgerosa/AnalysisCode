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


void makeVectorLongLived_1D(string inputFileName, string outputDIR, bool alongCtau){

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

  vector<TGraph*> hexp;
  vector<TGraph*> hobs;
  vector<int> expcounter;
  vector<int> obscounter;

  if(alongCtau){
    for(int idm = 0; idm < dmMassPoints.size(); idm++){
      expcounter.push_back(0);
      obscounter.push_back(0);
      hexp.push_back(new TGraph());
      hexp.back()->SetName(Form("hexp_dm_%d",int(dmMassPoints.at(idm))));
      hobs.push_back(new TGraph());
      hobs.back()->SetName(Form("hobs_dm_%d",int(dmMassPoints.at(idm))));
    }
  }
  else{
    for(int idm = 0; idm < ctauPoints.size(); idm++){
      expcounter.push_back(0);
      obscounter.push_back(0);
      hexp.push_back(new TGraph());
      hexp.back()->SetName(Form("hexp_ctau_%d",int(ctauPoints.at(idm))));
      hobs.push_back(new TGraph());
      hobs.back()->SetName(Form("hobs_ctau_%d",int(ctauPoints.at(idm))));
    }
  }

  reader.SetEntry(0);

  double maxLimit = 0;
  double minLimit = 999;

  while(reader.Next()){
    
    medmass = mmed(*mh);
    int dmmass  = mdm(*mh);
    int ctau    = decaytime(*mh); 
    if (*quantile == 0.5) {
      int ipos = 0;
      for(auto hist: hexp){	
	if(TString(hist->GetName()).Contains(Form("dm_%d",dmmass))){
	  hist->SetPoint(expcounter.at(ipos),ctau,*limit);
	  if(*limit > maxLimit) maxLimit = *limit;
	  if(*limit < maxLimit) minLimit = *limit;
	  expcounter.at(ipos)++;
	}
	else if(TString(hist->GetName()).Contains(Form("ctau_%d",ctau))){
	  hist->SetPoint(expcounter.at(ipos),dmmass,*limit);
	  if(*limit > maxLimit) maxLimit = *limit;
	  if(*limit < maxLimit) minLimit = *limit;
	  expcounter.at(ipos)++;
	}
	ipos++;
      }
    }
    if (*quantile == -1) {
      int ipos = 0;
      for(auto hist: hobs){
	if(TString(hist->GetName()).Contains(Form("dm_%d",dmmass))){
	  if(*limit > maxLimit) maxLimit = *limit;
	  if(*limit < maxLimit) minLimit = *limit;
	  hist->SetPoint(obscounter.at(ipos),ctau,*limit);
	  obscounter.at(ipos)++;
	}
	else if(TString(hist->GetName()).Contains(Form("ctau_%d",ctau))){
	  if(*limit > maxLimit) maxLimit = *limit;
	  if(*limit < maxLimit) minLimit = *limit;
	  hist->SetPoint(obscounter.at(ipos),dmmass,*limit);
	  obscounter.at(ipos)++;
	}
	ipos++;
      }       
    }
  }
  
  // All the plotting and cosmetics                                                                                                                                                                    
  TCanvas* canvas = new TCanvas("canvas", "canvas",600,625);
  canvas->SetLogx();

  TLegend leg (0.7,0.7,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);

  int ihist = 0;
  int icolor = 1;

  for(auto hist : hexp){
    ihist++;    
    if(alongCtau)
      hist->GetXaxis()->SetTitle("c#tau [mm]");
    else
      hist->GetXaxis()->SetTitle("m_{X_{2}} [GeV]");

    hist->GetYaxis()->SetTitle("Expected #sigma_{95% CL}/#sigma_{th}");
    hist->GetXaxis()->SetTitleOffset(1.10);
    hist->GetYaxis()->SetTitleOffset(1.10);
    if(icolor == 3) icolor++;
    if(icolor == 5) icolor++;
    hist->SetLineColor(icolor);
    hist->SetLineWidth(2);    
    hist->SetMarkerStyle(19+icolor);    
    hist->SetMarkerSize(1.5);
    hist->SetMarkerColor(icolor);
    if(TString(hist->GetName()).Contains("dm"))
      leg.AddEntry(hist,Form("m_{X_{2}} = %d GeV",int(dmMassPoints.at(ihist-1))),"PL");
    else
      leg.AddEntry(hist,Form("c\tau = %d mm",int(ctauPoints.at(ihist-1))),"PL");
    if(ihist == 1){  
      hist->GetYaxis()->SetRangeUser(minLimit*0.1,maxLimit*10);
      hist->Draw("APL");      
    }
    else{
      hist->Draw("PLsame");
    }
    icolor++;
  }
  leg.Draw("same");
  CMS_lumi(canvas,"35.9");
  canvas->SetLogy();
  canvas->RedrawAxis("sameaxis");

  TLine* line = NULL;
  if(alongCtau)
    line = new TLine(ctauPoints.front(),1,ctauPoints.back(),1);
  else
    line = new TLine(dmMassPoints.front(),1,dmMassPoints.back(),1);
  line->SetLineStyle(7);
  line->Draw("Lsame");

  if(alongCtau){
    canvas->SaveAs((outputDIR+"/scan_vector_ctau_mmed_"+to_string(medmass)+"_expected.pdf").c_str());
    canvas->SaveAs((outputDIR+"/scan_vector_ctau_mmed_"+to_string(medmass)+"_expected.png").c_str());
  }
  else{
    canvas->SaveAs((outputDIR+"/scan_vector_mdm_mmed_"+to_string(medmass)+"_expected.pdf").c_str());
    canvas->SaveAs((outputDIR+"/scan_vector_mdm_mmed_"+to_string(medmass)+"_expected.png").c_str());
  }

}
