#include "../CMS_lumi.h"

void checkNegativeBin(TH1* histo){
  for(int iBin = 0; iBin < histo->GetNbinsX()+1; iBin++){
    if(histo->GetBinContent(iBin+1) <= 0)
      histo->SetBinContent(iBin+1,0.001);
  }
}


class massPoint {

public:

  massPoint(){};

  massPoint(float & medmass, float & dmmass){
    medmass_ =  medmass;
    dmmass_  =  dmmass;
  }

  bool operator < (const massPoint & a) const {
    if(medmass_ < a.medmass_) return true;
    else if(medmass_ > a.medmass_) return false;
    else if(medmass_ == a.medmass_){
      if(dmmass_ <= a.dmmass_) return true;
      else if(dmmass_ > a.dmmass_) return false;
    }
    
    return true;
  }

  bool operator == (const massPoint & a) const{
    if(medmass_ == a.medmass_ and dmmass_ == a.dmmass_) return true;
    else return false;
  }

  bool operator != (const massPoint & a) const{
    if(medmass_ != a.medmass_ or dmmass_ != a.dmmass_) return true;
    else return false;
  }
  
  float  medmass_;
  float  dmmass_;
  double sumwgt_;
  TH1F*  histo_nominal_;
  vector<TH1F*> histo_;
};

// binning in met
vector<float> bins = {250.,280.0,310.0,340.0,370.0,400.0,430.0,470.0,510.0,550.0,590.0,640.0, 690.0, 740.0, 790.0, 840.0, 900.0, 960.0, 1020.0, 1090.0, 1160.0, 1250.0, 1400.0};

vector<int> color  = { 1, 2, 4,862,858,870,867,864,856,852,800,804,807,797,793};
vector<int> marker = {20,21,24, 20, 21, 24, 20, 21, 24, 20, 21, 24, 20, 21, 24};

// main function
static float luminosity = 35.9;
void makeCouplingComparison(string inputFileName, string outputDIR, bool isSMM, bool isSpin1){

  gROOT->SetBatch(kTRUE);
  system(("mkdir -p "+outputDIR).c_str());
  setTDRStyle();

  TFile* inputFile = TFile::Open(inputFileName.c_str(),"READ");
  inputFile->cd();

  TTree* tree = (TTree*) inputFile->Get("tree");

  TTreeReader reader (tree);
  TTreeReaderValue<int> id (reader,"id");
  TTreeReaderValue<float> xsec (reader,"xsec");
  TTreeReaderValue<double> sumwgt (reader,"sumwgt");
  TTreeReaderValue<float> pfMetPt (reader,"pfMetPt");
  TTreeReaderValue<float> weightPU (reader,"weightPU");
  TTreeReaderValue<float> weightTurnOn (reader,"weightTurnOn");
  TTreeReaderValue<float> genWeight (reader,"genWeight");
  TTreeReaderValue<float> genMediatorMass (reader,"genMediatorMass");
  TTreeReaderValue<float> genX1Mass (reader,"genX1Mass");
  
  string gV_b = "gSM";
  string gDM_b = "gDM";
  if(isSMM){
    gV_b  = "theta";
    gDM_b = "yDM";
  }

  TTreeReaderValue<vector<float> > gSM (reader,gV_b.c_str());
  TTreeReaderValue<vector<float> > gDM (reader,gDM_b.c_str());
  TTreeReaderValue<vector<float> > couplingwgt (reader,"couplingwgt");

  vector<massPoint> massList;
  // first loop on the file dumping all mass points
  cout<<"Start loop on the file to get mass points "<<endl;
  while(reader.Next()){
    if(massList.size() == 0){
      massList.push_back(massPoint(*genMediatorMass,*genX1Mass));
      massList.back().sumwgt_ = 0;
      if(gSM->size() != 0 and gDM->size() !=0 and couplingwgt->size() !=0)
	massList.back().sumwgt_ += *genWeight;
    }
    else{
      if(massList.back().medmass_ == *genMediatorMass and massList.back().dmmass_ == *genX1Mass){
	if(gSM->size() != 0 and gDM->size() !=0 and couplingwgt->size() !=0)
	  massList.back().sumwgt_ += *genWeight;
	continue; // assumption is that mass points are in sequence
      }
      else{
	massList.push_back(massPoint(*genMediatorMass,*genX1Mass));
	massList.back().sumwgt_ = 0;
	if(gSM->size() != 0 and gDM->size() !=0 and couplingwgt->size() !=0)
	  massList.back().sumwgt_ += *genWeight;
      }
    }
  }
  // sorting mass points
  sort(massList.begin(),massList.end());

  cout<<"Mass point found: "<<endl;
  for(auto mass : massList)
    cout<<"Mediator mass "<<mass.medmass_<<" Dark matter mass "<<mass.dmmass_<<" weight "<<mass.sumwgt_<<endl;
  
  // loop on the even to fill the distributions
  cout<<"Filling the distributions "<<endl;
  long int firstEvent = 0;
  long int nTotal = tree->GetEntries();
  reader.SetEntry(0);
  massPoint currentMassPoint;
  
  while(reader.Next()){
    // skip bad mass points
    if(gSM->size() != gDM->size() or gSM->size() != couplingwgt->size() or gDM->size() != couplingwgt->size()) continue;
    if(gSM->size() == 0 or gDM->size() == 0 or couplingwgt->size() == 0) continue;
    
    massPoint a (*genMediatorMass,*genX1Mass);    
    if(firstEvent == 0){
      firstEvent ++ ;
      currentMassPoint = a;
    }
    else{
      if(currentMassPoint != a){
	cout<<"Finished last mass point: med mass "<<currentMassPoint.medmass_<<" dm mass "<<currentMassPoint.dmmass_<<endl;
	currentMassPoint  = a;
      }
    }
    
    vector<massPoint>::iterator itMass = find(massList.begin(),massList.end(),a);
    if(itMass != massList.end()){
      
      // Book the histogram
      if((*itMass).histo_.size() == 0){
	(*itMass).histo_nominal_ = new TH1F(Form("histo_medmass_%d_dmmass_%d_nominal",int(a.medmass_),int(a.dmmass_)),"",bins.size()-1,&bins[0]);
	(*itMass).histo_nominal_->Sumw2();
	for(int pos = 0; pos < gSM->size(); pos++){
	  if(gSM->at(pos) == 0 or gDM->at(pos) == 0) continue;
	  (*itMass).histo_.push_back(new TH1F(Form("histo_medmass_%d_dmmass_%d_gSM_%f_gDM_%f",int(a.medmass_),int(a.dmmass_),gSM->at(pos),gDM->at(pos)),"",bins.size()-1,&bins[0]));
	  (*itMass).histo_.back()->Sumw2();
	}
      }
      
      // Fill nominal histogram after monojet selections
      if(*id == 1){
	(*itMass).histo_nominal_->Fill(*pfMetPt,(*xsec)*luminosity*(*genWeight)*(*weightPU)*(*weightTurnOn)*1./((*itMass).sumwgt_));
	int realPos = 0;
	for(int pos = 0; pos < gSM->size(); pos++){
	  if(gSM->at(pos) == 0 or gDM->at(pos) == 0) continue;
	  if(realPos >= (*itMass).histo_.size()) continue;
	  (*itMass).histo_.at(realPos)->Fill(*pfMetPt,(*xsec)*luminosity*(*genWeight)*(*weightPU)*(*weightTurnOn)*(couplingwgt->at(pos)/(*genWeight))*1./((*itMass).sumwgt_));
	  realPos++;
	}
      }
    }
    else{
      cerr<<"Problem: mass point not found --> please check"<<endl;
    }    
  }

  ///////// make the plots and store histograms in a output file
  TFile* outputFile = new TFile((outputDIR+"/outputDistributions.root").c_str(),"RECREATE");
  outputFile->cd();
  for(auto mass : massList){
    TString dirName = Form("medmass_%d_dmmass_%d",int(mass.medmass_),int(mass.dmmass_));
    outputFile->mkdir(dirName);
    outputFile->cd(dirName);
    for(auto ihist : mass.histo_)
      ihist->Write();
    outputFile->cd();
  }

  if(not isSMM and isSpin1){

    // make plot for gSM = 0.25 vs gDM
    for(auto mass : massList){

      TCanvas* canvas = new TCanvas("canvas","canvas",600,625);
      canvas->cd();
      canvas->SetBottomMargin(0.3);
      TPad* pad = new TPad("pad","pad",0,0.,1,0.9);
      pad->SetTopMargin(0.7);
      pad->SetFillColor(0);
      pad->SetFillStyle(0);
      canvas->cd();
      
      TH1F* frame = new TH1F(Form("frame_gDM_med_%d_dm_%d",int(mass.medmass_),int(mass.dmmass_)),"",bins.size()-1,&bins[0]);
      TH1* frame2 =  (TH1*) frame->Clone("frame2");
      frame->GetYaxis()->SetTitle("Events");
      frame->GetXaxis()->SetTitleSize(0);
      frame->GetXaxis()->SetLabelSize(0);
      frame->Draw();
      CMS_lumi(canvas,"35.9");

      TLegend leg (0.6,0.7,0.99,0.99);
      leg.SetFillColor(0);
      leg.SetFillStyle(0);
      leg.SetBorderSize(0);
      leg.AddEntry((TObject*)(0),Form("m_{med} = %d, m_{dm} = %d",int(mass.medmass_),int(mass.dmmass_)),"");
      
      double minimum = 999999;
      double maximum = -1;
      unsigned int pos = 1;
      vector<TH1*> ratio;
      
      TH1* histo_nominal = mass.histo_nominal_;

      histo_nominal->SetLineWidth(2);
      histo_nominal->SetMarkerStyle(marker.at(0));
      histo_nominal->SetLineColor(color.at(0));
      histo_nominal->SetMarkerColor(color.at(0));
      histo_nominal->SetFillColor(kBlack);
      histo_nominal->SetFillStyle(3001);

      checkNegativeBin(histo_nominal);
      if(histo_nominal->GetMinimum() < minimum)
	minimum = histo_nominal->GetMinimum();
      if(histo_nominal->GetMaximum() > maximum)
	maximum = histo_nominal->GetMaximum();

      leg.AddEntry(histo_nominal,"nominal","PL");      
      histo_nominal->Draw("hist same");

      for(auto hist : mass.histo_){
	TString name = Form("%s",hist->GetName());
	if(not name.Contains("nominal") and not name.Contains("gSM_0.25")) continue;
	if(name.Contains("gDM_0.01")) continue;
	if(name.Contains("gDM_0.05")) continue;
	TH1* histo = hist;
	checkNegativeBin(histo);
	if(histo->GetMinimum() < minimum)
	  minimum = histo->GetMinimum();
	if(histo->GetMaximum() > maximum)
	  maximum = histo->GetMaximum();
	histo->SetLineWidth(2);
	histo->SetMarkerStyle(marker.at(pos));
	histo->SetLineColor(color.at(pos));
	histo->SetMarkerColor(color.at(pos));
	if(name.Contains("nominal")){
	  histo->SetFillColor(kBlack);
	  histo->SetFillStyle(3001);
	}
	canvas->cd();
	histo->Draw("hist same");
	histo->Draw("P same");
	pos++;
	name.ReplaceAll(Form("histo_medmass_%d_dmmass_%d_",int(mass.medmass_),int(mass.dmmass_)),"");
	name.ReplaceAll("_"," = ");
	name.ReplaceAll("0000","");
	name.ReplaceAll("000","");
	leg.AddEntry(histo,name,"PL");

	ratio.push_back((TH1*) histo->Clone("ratio"));
	ratio.back()->Divide(histo_nominal);
      }
      
      pad->Draw();
      pad->cd();      
      frame2->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
      frame2->GetYaxis()->SetTitle("Ratio");
      frame2->GetXaxis()->SetTitleOffset(1.1);
      frame2->GetYaxis()->SetTitleOffset(1.2);
      frame2->GetYaxis()->SetRangeUser(0,3);
      frame2->GetYaxis()->SetNdivisions(508);
      frame2->Draw();
      for(auto histo : ratio)
 	histo->Draw("hist same");	
      
      
      if(minimum != 0)
	frame->GetYaxis()->SetRangeUser(minimum*0.1,maximum*100);
      else
	frame->GetYaxis()->SetRangeUser(0.0001,maximum*100);
      
      leg.Draw("same");
      canvas->SetLogy();
            
      canvas->SaveAs((outputDIR+"/yields_vs_coupling_med_"+to_string(mass.medmass_)+"_dm_"+to_string(mass.dmmass_)+"_gSM0p25.png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/yields_vs_coupling_med_"+to_string(mass.medmass_)+"_dm_"+to_string(mass.dmmass_)+"_gSM0p25.pdf").c_str(),"pdf");
    }


    // make plot for gDM = 1.0 vs gSM
    for(auto mass : massList){

      TCanvas* canvas = new TCanvas("canvas","canvas",600,625);
      canvas->cd();
      canvas->SetBottomMargin(0.3);
      TPad* pad = new TPad("pad","pad",0,0.,1,0.9);
      pad->SetTopMargin(0.7);
      pad->SetFillColor(0);
      pad->SetFillStyle(0);
      canvas->cd();
      
      TH1F* frame = new TH1F(Form("frame_gSM_med_%d_dm_%d",int(mass.medmass_),int(mass.dmmass_)),"",bins.size()-1,&bins[0]);
      TH1* frame2 =  (TH1*) frame->Clone("frame2");
      frame->GetYaxis()->SetTitle("Events");
      frame->GetXaxis()->SetTitleSize(0);
      frame->GetXaxis()->SetLabelSize(0);
      frame->Draw();
      CMS_lumi(canvas,"35.9");

      TLegend leg (0.6,0.7,0.99,0.99);
      leg.SetFillColor(0);
      leg.SetFillStyle(0);
      leg.SetBorderSize(0);
      leg.AddEntry((TObject*)(0),Form("m_{med} = %d, m_{dm} = %d",int(mass.medmass_),int(mass.dmmass_)),"");
      
      double minimum = 999999;
      double maximum = -1;
      unsigned int pos = 1;
      vector<TH1*> ratio;
      
      TH1* histo_nominal = mass.histo_nominal_;

      histo_nominal->SetLineWidth(2);
      histo_nominal->SetMarkerStyle(marker.at(0));
      histo_nominal->SetLineColor(color.at(0));
      histo_nominal->SetMarkerColor(color.at(0));
      histo_nominal->SetFillColor(kBlack);
      histo_nominal->SetFillStyle(3001);

      checkNegativeBin(histo_nominal);
      if(histo_nominal->GetMinimum() < minimum)
	minimum = histo_nominal->GetMinimum();
      if(histo_nominal->GetMaximum() > maximum)
	maximum = histo_nominal->GetMaximum();

      leg.AddEntry(histo_nominal,"nominal","PL");      
      histo_nominal->Draw("hist same");

      for(auto hist : mass.histo_){
	TString name = Form("%s",hist->GetName());
	if(not name.Contains("nominal") and not name.Contains("gDM_1.0")) continue;
	TH1* histo = hist;
	checkNegativeBin(histo);
	if(histo->GetMinimum() < minimum)
	  minimum = histo->GetMinimum();
	if(histo->GetMaximum() > maximum)
	  maximum = histo->GetMaximum();
	histo->SetLineWidth(2);
	histo->SetMarkerStyle(marker.at(pos));
	histo->SetLineColor(color.at(pos));
	histo->SetMarkerColor(color.at(pos));
	if(name.Contains("nominal")){
	  histo->SetFillColor(kBlack);
	  histo->SetFillStyle(3001);
	}
	canvas->cd();
	histo->Draw("hist same");
	histo->Draw("P same");
	pos++;
	name.ReplaceAll(Form("histo_medmass_%d_dmmass_%d_",int(mass.medmass_),int(mass.dmmass_)),"");
	name.ReplaceAll("_"," = ");
	name.ReplaceAll("0000","");
	name.ReplaceAll("000","");
	leg.AddEntry(histo,name,"PL");

	ratio.push_back((TH1*) histo->Clone("ratio"));
	ratio.back()->Divide(histo_nominal);
      }
      
      pad->Draw();
      pad->cd();      
      frame2->GetXaxis()->SetTitle("E_{T}^{miss} [GeV]");
      frame2->GetYaxis()->SetTitle("Ratio");
      frame2->GetXaxis()->SetTitleOffset(1.1);
      frame2->GetYaxis()->SetTitleOffset(1.2);
      frame2->GetYaxis()->SetRangeUser(0,3);
      frame2->GetYaxis()->SetNdivisions(508);
      frame2->Draw();
      for(auto histo : ratio)
 	histo->Draw("hist same");	
      
      
      if(minimum != 0)
	frame->GetYaxis()->SetRangeUser(minimum*0.1,maximum*100);
      else
	frame->GetYaxis()->SetRangeUser(0.0001,maximum*100);
      
      leg.Draw("same");
      canvas->SetLogy();
            
      canvas->SaveAs((outputDIR+"/yields_vs_coupling_med_"+to_string(mass.medmass_)+"_dm_"+to_string(mass.dmmass_)+"_gDM1p0.png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/yields_vs_coupling_med_"+to_string(mass.medmass_)+"_dm_"+to_string(mass.dmmass_)+"_gDM1p0.pdf").c_str(),"pdf");
    }



  }
  outputFile->Close();
  
}
