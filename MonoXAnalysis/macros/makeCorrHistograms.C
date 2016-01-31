#include "makehist.h"

using namespace std;

void smoothEmptyBins(TH1* hist, int nsteps = 2){

  for(int iBin = 1 ; iBin <= hist->GetNbinsX(); iBin++){
    if(hist->GetBinContent(iBin) == 0){
      float average = 0.;
      for(int jBin = iBin -nsteps; jBin < iBin+nsteps; jBin++){
	if(jBin == iBin) continue;
	if(jBin > 0 and jBin <= hist->GetNbinsX()){
	  average += hist->GetBinContent(jBin);
	}
      }
      hist->SetBinContent(iBin,average/(nsteps*2)); 
      hist->SetBinError(iBin,hist->GetBinContent(iBin)*2); 
    }
  }
}


void smoothEmptyBins(TH2* hist, int nsteps = 1){

  for(int xBin = 1 ; xBin <= hist->GetNbinsX(); xBin++){
    for(int yBin = 1 ; yBin <= hist->GetNbinsY(); yBin++){
      if(hist->GetBinContent(xBin,yBin) == 0){
	float average = 0.;
	for(int jBin = xBin -nsteps; jBin < xBin+nsteps; jBin++){
	  for(int kBin = yBin -nsteps; kBin < yBin+nsteps; kBin++){
	    if(jBin == xBin and kBin == yBin) continue;
	    if((jBin > 0 and jBin <= hist->GetNbinsX()) and (kBin > 0 and kBin <= hist->GetNbinsY())){
	      average += hist->GetBinContent(jBin,kBin);
	    }
	  }
	}
	hist->SetBinContent(xBin,yBin,average/((nsteps*2+1)*(nsteps*2+1)-1)); 
	hist->SetBinError(xBin,yBin,hist->GetBinContent(xBin,yBin)*2); 
      }
    }
  }
}


// make histograms for Z->mumu to signal region correction                                                                                                                   
void makezmmcorhist( string  signalRegionFile,  string  zmumuFile,  string  kFactorFile, int category, vector<string> observables, double lumi, string outDir = "", string ext = "") {

  // open files                                                                                                                                                                
  TFile* nfile  = TFile::Open(signalRegionFile.c_str());
  TFile* dfile  = TFile::Open(zmumuFile.c_str());

  TTree* ntree = (TTree*) nfile->Get("tree/tree");
  TTree* dtree = (TTree*) dfile->Get("tree/tree");

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;

  vector<float> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    TH1F* nhist_temp = new TH1F(("nhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dhist_temp = new TH1F(("dhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
    dhist.push_back(dynamic_cast<TH1*>(dhist_temp));

  }

  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile(kFactorFile.c_str());
  TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
  TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
  TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");

  // Divide NLO/LO                                                                                                                                                             
  znlohist->Divide(zlohist);

  vector<TH1*> ehists;
  vector<TH1*> zhists;
  zhists.push_back(znlohist);
  zhists.push_back(zewkhist);

  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon, 
  makehist4(ntree, nhist, nhist_2D,  true, 0, category, false, 1.00, lumi, zhists, "", true, NULL);
  makehist4(dtree, dhist, dhist_2D,  true, 1, category, false, 1.00, lumi, zhists, "", true, NULL);

  string name = string("zmmcor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Divide(dhist.at(ihist));

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    nhist_2D.at(ihist)->Divide(dhist_2D.at(ihist));


  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    smoothEmptyBins(nhist.at(ihist),2);

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    smoothEmptyBins(nhist_2D.at(ihist),1);

  
  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    nhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist.at(ihist)->Write();
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    nhist_2D.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist_2D.at(ihist)->Write();
  }

  outfile.Close();
  nfile->Close();
  dfile->Close();

  cout << "Z(mumu)->Z(inv) transfer factor computed ..." << endl;
}


// make histograms for Z->ee to signal region correction                                                                                                                   
void makezeecorhist( string  signalRegionFile,  string  zeeFile,  string  kFactorFile, int category, vector<string> observables, double lumi, string outDir = "", string ext = "") {

  // open files                                                                                                                                                                
  TFile* nfile  = TFile::Open(signalRegionFile.c_str());
  TFile* dfile  = TFile::Open(zeeFile.c_str());

  TTree* ntree = (TTree*) nfile->Get("tree/tree");
  TTree* dtree = (TTree*) dfile->Get("tree/tree");

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;

  vector<float> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    TH1F* nhist_temp = new TH1F(("nhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dhist_temp = new TH1F(("dhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
    dhist.push_back(dynamic_cast<TH1*>(dhist_temp));

  }


  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile(kFactorFile.c_str());
  TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
  TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
  TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");

  // Divide NLO/LO                                                                                                                                                             
  znlohist->Divide(zlohist);

  vector<TH1*> ehists;
  vector<TH1*> zhists;
  zhists.push_back(znlohist);
  zhists.push_back(zewkhist);

  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon, 
  makehist4(ntree, nhist, nhist_2D,  true, 0, category, false, 1.00, lumi, zhists, "",true, NULL);
  makehist4(dtree, dhist, dhist_2D,  true, 3, category, false, 1.00, lumi, zhists, "",true, NULL);

  string name = string("zeecor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Divide(dhist.at(ihist));

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    nhist_2D.at(ihist)->Divide(dhist_2D.at(ihist));

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    smoothEmptyBins(nhist.at(ihist),2);

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    smoothEmptyBins(nhist_2D.at(ihist),1);
  
  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    nhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist.at(ihist)->Write();
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    nhist_2D.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist_2D.at(ihist)->Write();
  }

  outfile.Close();
  nfile->Close();
  dfile->Close();

  cout << "Z(ee)->Z(inv) transfer factor computed ..." << endl;
}



// make histograms for W->mnu to signal region correction                                                                                                                   
void makewmncorhist( string  signalRegionFile,  string  wmnFile,  string  kFactorFile, int category, vector<string> observables, double lumi, string outDir = "", string ext = "") {

  // open files                                                                                                                                                                
  TFile* nfile  = TFile::Open(signalRegionFile.c_str());
  TFile* dfile  = TFile::Open(wmnFile.c_str());

  TTree* ntree = (TTree*) nfile->Get("tree/tree");
  TTree* dtree = (TTree*) dfile->Get("tree/tree");

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;

  vector<float> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    TH1F* nhist_temp = new TH1F(("nhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dhist_temp = new TH1F(("dhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
    dhist.push_back(dynamic_cast<TH1*>(dhist_temp));

  }


  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile(kFactorFile.c_str());
  TH1* wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
  TH1*  wlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
  TH1* wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");

  wnlohist->Divide(wlohist);

  vector<TH1*> ehists;
  vector<TH1*> whists;
  whists.push_back(wnlohist);
  whists.push_back(wewkhist);

  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon, 
  makehist4(ntree, nhist, nhist_2D,  true, 0, category, false, 1.00, lumi, whists, "", true, NULL);
  makehist4(dtree, dhist, dhist_2D,  true, 2, category, false, 1.00, lumi, whists, "", true, NULL);

  string name = string("wmncor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Divide(dhist.at(ihist));

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    nhist_2D.at(ihist)->Divide(dhist_2D.at(ihist));

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    smoothEmptyBins(nhist.at(ihist),2);

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    smoothEmptyBins(nhist_2D.at(ihist),1);
  
  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    nhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist.at(ihist)->Write();
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    nhist_2D.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist_2D.at(ihist)->Write();
  }

  outfile.Close();
  nfile->Close();
  // dfile->Close();

  cout << "W(mnu)->W+Jets transfer factor computed ..." << endl;
}


// make histograms for W->enu to signal region correction                                                                                                                   
void makewencorhist( string  signalRegionFile,  string  wenFile,  string  kFactorFile, int category, vector<string> observables, double lumi, string outDir = "", string ext = "") {

  // open files                                                                                                                                                                
  TFile* nfile  = TFile::Open(signalRegionFile.c_str());
  TFile* dfile  = TFile::Open(wenFile.c_str());

  TTree* ntree = (TTree*) nfile->Get("tree/tree");
  TTree* dtree = (TTree*) dfile->Get("tree/tree");

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;

  vector<float> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    TH1F* nhist_temp = new TH1F(("nhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dhist_temp = new TH1F(("dhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
    dhist.push_back(dynamic_cast<TH1*>(dhist_temp));

  }

  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile(kFactorFile.c_str());
  TH1* wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
  TH1*  wlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
  TH1* wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");

  wnlohist->Divide(wlohist);

  vector<TH1*> ehists;
  vector<TH1*> whists;
  whists.push_back(wnlohist);
  whists.push_back(wewkhist);

  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon, 
  makehist4(ntree, nhist, nhist_2D,  true, 0, category, false, 1.00, lumi, whists, "", true, NULL);
  makehist4(dtree, dhist, dhist_2D,  true, 4, category, false, 1.00, lumi, whists, "", true, NULL);

  string name = string("wencor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Divide(dhist.at(ihist));

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    nhist_2D.at(ihist)->Divide(dhist_2D.at(ihist));

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    smoothEmptyBins(nhist.at(ihist),2);

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    smoothEmptyBins(nhist_2D.at(ihist),1);
  
  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    nhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist.at(ihist)->Write();
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    nhist_2D.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist_2D.at(ihist)->Write();
  }

  outfile.Close();
  nfile->Close();
  //dfile->Close();

  cout << "W(enu)->W+Jets transfer factor computed ..." << endl;
}


// make Z/W ratio
void  makezwjcorhist( string  znunuFile,  string  wlnuFile,  string  kFactorFile, int category, vector<string> observables, double lumi, string outDir = "",string ext = "",int kfact = 0) {

  // open files                                                                                                                                                                
  TFile* nfile  = TFile::Open(znunuFile.c_str());
  TFile* dfile  = TFile::Open(wlnuFile.c_str());

  TTree* ntree = (TTree*) nfile->Get("tree/tree");
  TTree* dtree = (TTree*) dfile->Get("tree/tree");

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;

  vector<float> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    TH1F* nhist_temp = new TH1F(("nhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dhist_temp = new TH1F(("dhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
    dhist.push_back(dynamic_cast<TH1*>(dhist_temp));

  }

  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile(kFactorFile.c_str());
  TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
  TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
  TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");
  TH1* zpdfhist = (TH1*)kffile.Get("znlo012/znlo012_pdfUp");

  TH1* wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
  TH1*  wlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
  TH1* wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");
  TH1* wpdfhist = (TH1*)kffile.Get("wnlo012/wnlo012_pdfUp");

  TH1* nomhist  = (TH1*)kffile.Get("znlo1_over_wnlo1/znlo1_over_wnlo1");
  TH1* re1hist  = (TH1*)kffile.Get("znlo1_over_wnlo1/znlo1_over_wnlo1_renCorrUp");
  TH1* re2hist  = (TH1*)kffile.Get("znlo1_over_wnlo1/znlo1_over_wnlo1_renAcorrUp");
  TH1* fa1hist  = (TH1*)kffile.Get("znlo1_over_wnlo1/znlo1_over_wnlo1_facCorrUp");
  TH1* fa2hist  = (TH1*)kffile.Get("znlo1_over_wnlo1/znlo1_over_wnlo1_facAcorrUp");


  // nlo correction for the PDF                                                                                                                                                
  zpdfhist->Divide(znlohist);
  wpdfhist->Divide(wnlohist);
  // Z/W NLO QCD re up / Z/W NLO QCD                                                                                                                                          
  re1hist->Divide(nomhist);
  // Z/W NLO QCD re EWK up / Z/W NLO QCD                                                                                                                                       
  re2hist->Divide(nomhist);
  // Z/W NLO QCD fac  up / Z/W NLO QCD                                                                                                                                         
  fa1hist->Divide(nomhist);
  // Z/W NLO QCD fac EWK up / Z/W NLO QCD                                                                                                                                       
  fa2hist->Divide(nomhist);

  // central value                                                                                                                                                              
  znlohist->Divide(zlohist);
  wnlohist->Divide(wlohist);

  vector<TH1*> zhists;
  vector<TH1*> whists;

  //kfact == 1 --> Znunu corrected for by NLO QCD, Wlnu by NLO QCD                                                                                                              
  if (kfact == 1) zhists.push_back(znlohist);
  if (kfact == 1) whists.push_back(wnlohist);
  //kfact == 2 --> Znunu corrected for by NLO QCD+EWK, Wlnu by NLO QCD+EWK                                                                                                      
  if (kfact == 2) {zhists.push_back(znlohist); zhists.push_back(zewkhist);}
  if (kfact == 2) {whists.push_back(wnlohist); whists.push_back(wewkhist);}
  //kfact == 3 --> Znunu and Wlnu by NLO QCD, ratio for ren scale up QCD                                                                                                        
  if (kfact == 3) {zhists.push_back(znlohist); zhists.push_back(re1hist) ;}
  if (kfact == 3) whists.push_back(wnlohist);
  //kfact == 4 --> Znunu and Wlnu by NLO QCD, ratio for fac scale up QCD                                                                                                        
  if (kfact == 4) {zhists.push_back(znlohist); zhists.push_back(fa1hist) ;}
  if (kfact == 4) whists.push_back(wnlohist);
  //kfact == 5 --> Znunu and Wlnu by NLO QCD, ratio for ren scale up EWK                                                                                                        
  if (kfact == 5) {zhists.push_back(znlohist); zhists.push_back(re2hist) ;}
  if (kfact == 5) whists.push_back(wnlohist);
  //kfact == 6 --> Znunu and Wlnu by NLO QCD, ratio for fac scale up EWK                                                                                                        
  if (kfact == 6) {zhists.push_back(znlohist); zhists.push_back(fa2hist) ;}
  if (kfact == 6) whists.push_back(wnlohist);
  //kfact == 7 --> Znunu corrected for by NLO NLO PDF, Wlnu by NLO                                                                                                              
  if (kfact == 7) {zhists.push_back(znlohist); zhists.push_back(zpdfhist);}
  if (kfact == 7) {whists.push_back(wnlohist); whists.push_back(wpdfhist);}

  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon, 
  makehist4(ntree, nhist, nhist_2D,  true, 0, category, false, 1.00, lumi, zhists, "", true, NULL);
  makehist4(dtree, dhist, dhist_2D,  true, 0, category, false, 1.00, lumi, whists, "", true, NULL);

  string name = string("zwjcor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Divide(dhist.at(ihist));

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    nhist_2D.at(ihist)->Divide(dhist_2D.at(ihist));

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    smoothEmptyBins(nhist.at(ihist),2);

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    smoothEmptyBins(nhist_2D.at(ihist),1);
  
  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    nhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist.at(ihist)->Write();
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    nhist_2D.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist_2D.at(ihist)->Write();
  }

  outfile.Close();
  nfile->Close();
  dfile->Close();

  cout << "W+Jets->Z+inv transfer factor computed ..." << endl;
}


// make Z/gamma ratio
void makegamcorhist( string  znunuFile,  string  photonFile,  string  kFactorFile,  string  fPfile, int category, vector<string> observables, double lumi, string outDir = "", string ext = "",int kfact = 0) {

  // open files                                                                                                                                                                
  TFile* nfile  = TFile::Open(znunuFile.c_str());
  TFile* dfile  = TFile::Open(photonFile.c_str());

  TTree* ntree = (TTree*) nfile->Get("tree/tree");
  TTree* dtree = (TTree*) dfile->Get("tree/tree");

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;

  vector<float> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    TH1F* nhist_temp = new TH1F(("nhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dhist_temp = new TH1F(("dhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
    dhist.push_back(dynamic_cast<TH1*>(dhist_temp));

  }

  // k-factors file from generator lebel: Z-boson pt at LO, NLO QCD and NLO QCD+EWK                                                                                         
  TFile kffile(kFactorFile.c_str());
  TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
  TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
  TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");
  TH1* zpdfhist = (TH1*)kffile.Get("znlo012/znlo012_pdfUp");

  TH1* anlohist = (TH1*)kffile.Get("anlo1/anlo1_nominal");
  TH1*  alohist = (TH1*)kffile.Get("alo/alo_nominal");
  TH1* aewkhist = (TH1*)kffile.Get("a_ewkcorr/a_ewkcorr");
  TH1* apdfhist = (TH1*)kffile.Get("anlo1/anlo1_pdfUp");

  TH1* nomhist  = (TH1*)kffile.Get("znlo1_over_anlo1/znlo1_over_anlo1");
  TH1* re1hist  = (TH1*)kffile.Get("znlo1_over_anlo1/znlo1_over_anlo1_renCorrUp");
  TH1* re2hist  = (TH1*)kffile.Get("znlo1_over_anlo1/znlo1_over_anlo1_renAcorrUp");
  TH1* fa1hist  = (TH1*)kffile.Get("znlo1_over_anlo1/znlo1_over_anlo1_facCorrUp");
  TH1* fa2hist  = (TH1*)kffile.Get("znlo1_over_anlo1/znlo1_over_anlo1_facAcorrUp");
  
  // ZNLO PDF UP / ZNLO                                                                                                                                                        
  zpdfhist->Divide(znlohist);

  // gamma NLO PDF UP / gamma NLO                                                                                                                                             
  apdfhist->Divide(anlohist);
  // Z/gam NLO re QCD Up / Z/gamma NLO                                                                                                                                       
  re1hist->Divide(nomhist);
  // Z/gam NLO re EWK Up / Z/gamma NLO                                                                                                                                        
  re2hist->Divide(nomhist);
  // Z/gam NLO fact QCD Up / Z/gamma NLO                                                                                                                                      
  fa1hist->Divide(nomhist);
  // Z/gam NLO fact EWK Up / Z/gamma NLO                                                                                                                                      
  fa2hist->Divide(nomhist);

  znlohist->Divide(zlohist);
  anlohist->Divide(alohist);

  TFile fpfile(fPfile.c_str());
  TH1* afpchist = (TH1*)fpfile.Get("FP_Down");

  vector<TH1*> zhists;
  vector<TH1*> ahists;
  vector<TH1*> ehists;
  // ZNLO QCD and Gamma NLO QCD                                                                                                                                                
  if (kfact == 1) zhists.push_back(znlohist);
  if (kfact == 1) ahists.push_back(anlohist);

  //ZNLO QCD+EWK and Gamma NLO QCD+EWK                                                                                                                                         
  if (kfact == 2) {zhists.push_back(znlohist); zhists.push_back(zewkhist);}
  if (kfact == 2) {ahists.push_back(anlohist); ahists.push_back(aewkhist);}

  // ZNLO QCD+Re up and Gamma NLO QCD                                                                                                                                          
  if (kfact == 3) {zhists.push_back(znlohist); zhists.push_back(re1hist) ;}
  if (kfact == 3) ahists.push_back(anlohist);

  // ZNLO QCD + fact Up and Gamma NLO QCD                                                                                                                                      
  if (kfact == 4) {zhists.push_back(znlohist); zhists.push_back(fa1hist) ;}
  if (kfact == 4) ahists.push_back(anlohist);

  // ZNLO QCD + re EWK up and Gamma NLO QCD                                                                                                                                   
  if (kfact == 5) {zhists.push_back(znlohist); zhists.push_back(re2hist) ;}
  if (kfact == 5) ahists.push_back(anlohist);

  // ZNLO QCD + fact EWK up and Gamma NLO QCD                                                                                                                                  
  if (kfact == 6) {zhists.push_back(znlohist); zhists.push_back(fa2hist) ;}
  if (kfact == 6) ahists.push_back(anlohist);

  // ZNLO QCD + PDF up and Gamma NLO QCD + PDF Up                                                                                                                               
  if (kfact == 7) {zhists.push_back(znlohist); zhists.push_back(zpdfhist);}
  if (kfact == 7) {ahists.push_back(anlohist); ahists.push_back(apdfhist);}

  // ZNLO QCD and Gamma NLO QCD + FP                                                                                                                                            
  if (kfact == 8) zhists.push_back(znlohist);
  if (kfact == 8) {ahists.push_back(anlohist); zhists.push_back(afpchist);}


  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 1 == di-muon, 
  makehist4(ntree, nhist, nhist_2D,  true, 0, category, false, 1.00, lumi, zhists, "", true, NULL);
  makehist4(dtree, dhist, dhist_2D,  true, 5, category, false, 1.00, lumi, ahists, "", true, NULL);

  string name = string("gamcor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Divide(dhist.at(ihist));

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    nhist_2D.at(ihist)->Divide(dhist_2D.at(ihist));

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    smoothEmptyBins(nhist.at(ihist),2);

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    smoothEmptyBins(nhist_2D.at(ihist),1);
  
  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    nhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist.at(ihist)->Write();
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    nhist_2D.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist_2D.at(ihist)->Write();
  }

  outfile.Close();
  nfile->Close();
  dfile->Close();

  cout << "Gamma+Jets->Z+inv transfer factor computed ..." << endl;
}

// correction for top
void maketopmucorhist( string  signalRegionFile,  string  topFile,  int category, vector<string> observables, double lumi, string outDir = "", string sys = "", string ext = ""){

  // open files                                                                                                                                                                
  TFile* nfile  = TFile::Open(signalRegionFile.c_str());
  TFile* dfile  = TFile::Open(topFile.c_str());

  TTree* ntree = (TTree*) nfile->Get("tree/tree");
  TTree* dtree = (TTree*) dfile->Get("tree/tree");

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;

  vector<float> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    TH1F* nhist_temp = new TH1F(("nhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dhist_temp = new TH1F(("dhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
    dhist.push_back(dynamic_cast<TH1*>(dhist_temp));

  }

  vector<TH1*> ehists;
  vector<TH1*> zhists;

  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 7 == b-tagged region, 
  makehist4(ntree, nhist, nhist_2D,  true, 0, category, false, 1.00, lumi, zhists, sys, true, NULL);
  makehist4(dtree, dhist, dhist_2D,  true, 7, category, false, 1.00, lumi, zhists, sys, true, NULL);

  string name = string("topmucor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Divide(dhist.at(ihist));

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    nhist_2D.at(ihist)->Divide(dhist_2D.at(ihist));

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    smoothEmptyBins(nhist.at(ihist),2);

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    smoothEmptyBins(nhist_2D.at(ihist),1);
  
  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    nhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist.at(ihist)->Write();
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    nhist_2D.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist_2D.at(ihist)->Write();
  }

  outfile.Close();
  nfile->Close();
  dfile->Close();

  cout << "Top(b-tag,mu)->Top(b-veto) transfer factor computed ..." << endl;
}


// correction for top
void maketopelcorhist( string  signalRegionFile,  string  topFile,  int category, vector<string> observables, double lumi, string outDir = "", string sys = "", string ext = "") {

  // open files                                                                                                                                                                
  TFile* nfile  = TFile::Open(signalRegionFile.c_str());
  TFile* dfile  = TFile::Open(topFile.c_str());

  TTree* ntree = (TTree*) nfile->Get("tree/tree");
  TTree* dtree = (TTree*) dfile->Get("tree/tree");

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;

  vector<float> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
      
    TH1F* nhist_temp = new TH1F(("nhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    TH1F* dhist_temp = new TH1F(("dhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
    nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
    dhist.push_back(dynamic_cast<TH1*>(dhist_temp));

  }

  vector<TH1*> ehists;
  vector<TH1*> zhists;

  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 7 == b-tagged region, 
  makehist4(ntree, nhist, nhist_2D,  true, 0, category, false, 1.00, lumi, zhists, sys, true, NULL);
  makehist4(dtree, dhist, dhist_2D,  true, 8, category, false, 1.00, lumi, zhists, sys, true, NULL);

  string name = string("topelcor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Divide(dhist.at(ihist));

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    nhist_2D.at(ihist)->Divide(dhist_2D.at(ihist));

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    smoothEmptyBins(nhist.at(ihist),2);

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    smoothEmptyBins(nhist_2D.at(ihist),1);
  
  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    nhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist.at(ihist)->Write();
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    nhist_2D.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist_2D.at(ihist)->Write();
  }

  outfile.Close();
  nfile->Close();
  dfile->Close();

  cout << "Top(b-tag,el)->Top(b-veto) transfer factor computed ..." << endl;
}


// correction for Z(nunu) or W+jets mass sidebaand
void makesidebandcorhist( string  signalRegionFile,  string  sidebandFile,  int category_num, int category_den, vector<string> observables, double lumi, string outDir = "", string ext = "") {

  // open files                                                                                                                                                                
  TFile* nfile  = TFile::Open(signalRegionFile.c_str());
  TFile* dfile  = TFile::Open(sidebandFile.c_str());

  TTree* ntree = (TTree*) nfile->Get("tree/tree");
  TTree* dtree = (TTree*) dfile->Get("tree/tree");

  // create histograms                                                                                                                                                         
  vector<TH1*> nhist;
  vector<TH1*> dhist;
  vector<TH2*> nhist_2D;
  vector<TH2*> dhist_2D;



  vector<float> bins;
  for(auto obs : observables){
    bins = selectBinning(obs,category_num);
    if(bins.empty())
      cout<<"No binning for this observable --> please define it"<<endl;
    
      TH1F* nhist_temp = new TH1F(("nhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dhist_temp = new TH1F(("dhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      nhist.push_back(dynamic_cast<TH1*>(nhist_temp));
      dhist.push_back(dynamic_cast<TH1*>(dhist_temp));
  }

  vector<TH1*> ehists;
  vector<TH1*> zhists;

  // loop over ntree and dtree events isMC=true, sample 0 == signal region, sample 7 == b-tagged region, 
  makehist4(ntree, nhist, nhist_2D,  true, 0, category_num, false, 1.00, lumi, zhists, "", true, NULL);
  makehist4(dtree, dhist, dhist_2D,  true, 0, category_den, false, 1.00, lumi, zhists, "", true, NULL);

  string name = string("sidebandcor")+ext;

  // divide the two                                                                                                                                                          
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    nhist.at(ihist)->Divide(dhist.at(ihist));

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    nhist_2D.at(ihist)->Divide(dhist_2D.at(ihist));

  //check for empty bins and apply smoothing
  for(size_t ihist = 0; ihist < nhist.size(); ihist++)
    smoothEmptyBins(nhist.at(ihist),2);

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++)
    smoothEmptyBins(nhist_2D.at(ihist),1);
  
  // create output file                                                                                                                                                        
  TFile outfile((outDir+"/"+name+".root").c_str(), "RECREATE");
  for(size_t ihist = 0; ihist < nhist.size(); ihist++){
    nhist.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist.at(ihist)->Write();
  }

  for(size_t ihist = 0; ihist < nhist_2D.size(); ihist++){
    nhist_2D.at(ihist)->SetName((name+"hist_"+observables.at(ihist)).c_str());
    nhist_2D.at(ihist)->Write();
  }

  outfile.Close();
  nfile->Close();
  dfile->Close();

  cout << "Sideband->Signal region transfer factor computed ..." << endl;
}
