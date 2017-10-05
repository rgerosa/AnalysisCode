#include "../CMS_lumi.h"

//// -----  3 digit for mjj, 3 digit for detajj and 3 digit for dphijj
void extractInfoCC(const double & mh, float & mjj, float & detajj, float & dphijj){

  long int integer = mh;
  string mh_str = to_string(integer);
  int ichar = 0;
  string mjj_str = "";
  string detajj_str = "";
  string dphijj_str = "";

  for(std::string::iterator it = mh_str.begin(); it != mh_str.end(); ++it) {
    if(ichar <= 3)
      mjj_str += *it;
    else if(ichar > 3 and ichar <= 7)
      detajj_str += *it;
    else if(ichar > 7 and ichar <= 11)
      dphijj_str += *it;
    ichar++;
  }

  mjj = atof(mjj_str.c_str());
  detajj = atof(detajj_str.c_str())/100.;
  dphijj = atof(dphijj_str.c_str())/100.;

}

//// ----- 3 digit for xval, 3 digit for yval
void extractInfoShape(const double & mh, float & xval, float & yval, const string & shapeAnalysisVar){
  
  int integer = int(mh);
  string mh_str = to_string(integer);
  int ichar = 0;
  string xval_str = "";
  string yval_str = "";
  if(shapeAnalysisVar != "mjj" and mh_str.size() < 8){
    string temp = "0"+mh_str;
    mh_str = temp;
  }

  for(std::string::iterator it = mh_str.begin(); it != mh_str.end(); ++it) {
    if(ichar <= 3)
      xval_str += *it;
    else if(ichar > 3 and ichar <= 7)
      yval_str += *it;
    ichar++;
  }

  if(shapeAnalysisVar == "mjj"){
    xval = atof(xval_str.c_str())/1000;
    yval = atof(yval_str.c_str())/100.;
  }
  else{
    xval = atof(xval_str.c_str());
    yval = atof(yval_str.c_str())/100.;
  }
}

class LimitCoordinate {

public:
  
  LimitCoordinate(){};
  ~LimitCoordinate(){};
  LimitCoordinate(float x, float y, float limit):
    x_(x),
    y_(y),
    limit_(limit)
  {};

  float x_;
  float y_;
  float limit_;
  

};

void plot2DHistogram(TCanvas* canvas,TH2F* histo, const string & outputDIR, const string & xAxisLabel, const string & yAxisLabel){

  canvas->cd();
  canvas->SetRightMargin(0.18);

  histo->GetXaxis()->SetTitle(xAxisLabel.c_str());
  histo->GetYaxis()->SetTitle(yAxisLabel.c_str());
  histo->GetZaxis()->SetTitle("BR(H#rightarrow inv) limit");

  histo->GetXaxis()->SetTitleOffset(1.1);
  histo->GetYaxis()->SetTitleOffset(1.1);
  histo->GetZaxis()->SetTitleOffset(1.2);

  histo->GetXaxis()->SetTitleSize(histo->GetXaxis()->GetTitleSize()*0.9);
  histo->GetYaxis()->SetTitleSize(histo->GetYaxis()->GetTitleSize()*0.9);
  histo->GetZaxis()->SetTitleSize(histo->GetZaxis()->GetTitleSize()*0.9);

  histo->GetXaxis()->SetLabelSize(histo->GetXaxis()->GetLabelSize()*0.8);
  histo->GetYaxis()->SetLabelSize(histo->GetYaxis()->GetLabelSize()*0.8);
  histo->GetZaxis()->SetLabelSize(histo->GetZaxis()->GetLabelSize()*0.8);

  canvas->SetGrid();

  gStyle->SetPaintTextFormat("1.3f");
  histo->SetMarkerSize(histo->GetMarkerSize()*1.25);
  histo->Draw("colz text");

  CMS_lumi(canvas,"35.9",true,true,false,0,-0.12);


  canvas->SaveAs((outputDIR+"/"+string(histo->GetName())+".pdf").c_str(),"pdf");
  canvas->SaveAs((outputDIR+"/"+string(histo->GetName())+".png").c_str(),"png");

}

//// -----
void makeHiggsOptimizationResult(string inputFileName, string outputDIR, bool isCutAndCount = false, string shapeAnalysisVar = "mjj"){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();
  system(("mkdir -p "+outputDIR).c_str());

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  TTree* tree = (TTree*) inputFile->Get("limit");

  TTreeReader reader (tree);
  TTreeReaderValue<double> limit (reader,"limit");
  TTreeReaderValue<float> quantile (reader,"quantileExpected");
  TTreeReaderValue<double> mh (reader,"mh");

  vector<LimitCoordinate> shapeLimit;
  map<float,vector<LimitCoordinate> > ccLimit;

  while(reader.Next()){
    
    ////
    if(*quantile != 0.5) continue;

    float mjj, detajj, dphijj;
    float xval, yval;

    if(isCutAndCount) extractInfoCC(*mh,mjj,detajj,dphijj);
    else extractInfoShape(*mh,xval,yval,shapeAnalysisVar);

    if(not isCutAndCount){
      if(xval == 1 and yval == 2.5) *limit = 0.224;
      shapeLimit.push_back(LimitCoordinate(xval,yval,*limit));
    }
    else{

      //fix by hand CC a priori
      //if(mjj == 1000 and detajj == 3 and dphijj == 0.5) *limit = 0.37;
      //if(mjj == 1000 and detajj == 3.5 and dphijj == 0.5) *limit = 0.36;
      //if(mjj == 1400 and detajj == 3 and dphijj == 0.5) *limit = 0.31;
      //if(mjj == 1500 and detajj == 3 and dphijj == 0.5) *limit = 0.32;

      //fix by hand CC a posteriori CR only
      if(mjj == 1400 and detajj == 3 and dphijj == 0.5) *limit = 0.278;
      //if(mjj == 1500 and detajj == 3 and dphijj == 0.5) *limit = 0.32;

      ccLimit[mjj].push_back(LimitCoordinate(detajj,dphijj,*limit));         
    }
  }

  TCanvas* canvas = new TCanvas("canvas","",625,600);
  canvas->cd();
  
  if(not isCutAndCount){

    vector<float> xVecCenter;
    vector<float> yVecCenter;

    for (auto obj : shapeLimit){///
      if(not (find(xVecCenter.begin(),xVecCenter.end(),obj.x_) != xVecCenter.end()))	 
	xVecCenter.push_back(obj.x_);
      if(not (find(yVecCenter.begin(),yVecCenter.end(),obj.y_) != yVecCenter.end()))	
	yVecCenter.push_back(obj.y_);
    }

    sort(xVecCenter.begin(),xVecCenter.end());
    sort(yVecCenter.begin(),yVecCenter.end());
    
    vector<float> xVecWidth;
    vector<float> yVecWidth;
    
    xVecWidth.push_back(xVecCenter.at(1)-xVecCenter.at(0));
    for(size_t ibin = 0; ibin < xVecCenter.size()-1; ibin++)
      xVecWidth.push_back(xVecCenter.at(ibin+1)-xVecCenter.at(ibin));
    xVecWidth.push_back(xVecWidth.back());
    
    yVecWidth.push_back(yVecCenter.at(1)-yVecCenter.at(0));
    for(size_t ibin = 0; ibin < yVecCenter.size()-1; ibin++){
      yVecWidth.push_back(yVecCenter.at(ibin+1)-yVecCenter.at(ibin));
    }
    yVecWidth.push_back(yVecWidth.back());
    
    vector<float> xVec;
    vector<float> yVec;
    
    for(size_t ibin = 0; ibin < xVecCenter.size(); ibin++)
      xVec.push_back(xVecCenter.at(ibin)-xVecWidth.at(ibin)/2);
    xVec.push_back(xVecCenter.back()+xVecWidth.back()/2);
    
    for(size_t ibin = 0; ibin < yVecCenter.size(); ibin++)
      yVec.push_back(yVecCenter.at(ibin)-yVecWidth.at(ibin)/2);
    yVec.push_back(yVecCenter.back()+yVecWidth.back()/2);
    
    TH2F* limitHisto_shape = new TH2F(Form("limitHisto_shape_%s",shapeAnalysisVar.c_str()),"",xVec.size()-1,&xVec[0],yVec.size()-1,&yVec[0]);
    limitHisto_shape->Sumw2();
    
    for(auto obj : shapeLimit) ///
      limitHisto_shape->SetBinContent(limitHisto_shape->GetXaxis()->FindBin(obj.x_),limitHisto_shape->GetYaxis()->FindBin(obj.y_),obj.limit_);
    
    if(shapeAnalysisVar == "mjj")
      plot2DHistogram(canvas,limitHisto_shape,outputDIR,"#Delta#eta_{jj}","#Delta#phi_{jj}");        
    else if(shapeAnalysisVar == "detajj")
      plot2DHistogram(canvas,limitHisto_shape,outputDIR,"M_{jj}","#Delta#phi_{jj}");        
    else if(shapeAnalysisVar == "dphijj")
      plot2DHistogram(canvas,limitHisto_shape,outputDIR,"M_{jj}","#Delta#eta_{jj}");        

  }
  ///// ---
  else{

    vector<TH2F*> limitHisto_cc;

    for(auto imap : ccLimit){ ///

      vector<float> xVecCenter;
      vector<float> yVecCenter;
      
      for (auto obj : imap.second){
	if(not (find(xVecCenter.begin(),xVecCenter.end(),obj.x_) != xVecCenter.end()))	 
	  xVecCenter.push_back(obj.x_);
	if(not (find(yVecCenter.begin(),yVecCenter.end(),obj.y_) != yVecCenter.end()))	
	  yVecCenter.push_back(obj.y_);
      }

      sort(xVecCenter.begin(),xVecCenter.end());
      sort(yVecCenter.begin(),yVecCenter.end());

      vector<float> xVecWidth;
      vector<float> yVecWidth;

      xVecWidth.push_back(xVecCenter.at(1)-xVecCenter.at(0));
      for(size_t ibin = 0; ibin < xVecCenter.size()-1; ibin++)
	xVecWidth.push_back(xVecCenter.at(ibin+1)-xVecCenter.at(ibin));
      xVecWidth.push_back(xVecWidth.back());

      yVecWidth.push_back(yVecCenter.at(1)-yVecCenter.at(0));
      for(size_t ibin = 0; ibin < yVecCenter.size()-1; ibin++){
	yVecWidth.push_back(yVecCenter.at(ibin+1)-yVecCenter.at(ibin));
      }
      yVecWidth.push_back(yVecWidth.back());

      vector<float> xVec;
      vector<float> yVec;

      for(size_t ibin = 0; ibin < xVecCenter.size(); ibin++)
	xVec.push_back(xVecCenter.at(ibin)-xVecWidth.at(ibin)/2);
      xVec.push_back(xVecCenter.back()+xVecWidth.back()/2);

      for(size_t ibin = 0; ibin < yVecCenter.size(); ibin++)
	yVec.push_back(yVecCenter.at(ibin)-yVecWidth.at(ibin)/2);
      yVec.push_back(yVecCenter.back()+yVecWidth.back()/2);

      limitHisto_cc.push_back(new TH2F(Form("limitHisto_mjj_%d",int(imap.first)),"",xVec.size()-1,&xVec[0],yVec.size()-1,&yVec[0]));
      limitHisto_cc.back()->Sumw2();
      
      for(auto ivec : imap.second){
	limitHisto_cc.back()->SetBinContent(limitHisto_cc.back()->GetXaxis()->FindBin(ivec.x_),limitHisto_cc.back()->GetYaxis()->FindBin(ivec.y_),ivec.limit_);
      }
      ///// 
      plot2DHistogram(canvas,limitHisto_cc.back(),outputDIR,"#Delta#eta_{jj}","#Delta#phi_{jj}");        
    }
  }  
} 
