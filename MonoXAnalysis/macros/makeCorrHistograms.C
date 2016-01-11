#include "makehist.h"

// make histograms for Z->mumu to signal region correction                                                                                                                   
void makezmmcorhist( std::string  signalRegionFile,  std::string  zmumuFile,  std::string  kFactorFile) {

  // open files                                                                                                                                                                
  TFile* nfile  = TFile::Open(signalRegionFile.c_str());
  TFile* dfile  = TFile::Open(zmumuFile.c_str());

  TTree* ntree = (TTree*) nfile->Get("tree");
  TTree* dtree = (TTree*) dfile->Get("tree");

  // create histograms                                                                                                                                                         
  TH1F nhist("nhist", "", nrbins, rbins);
  TH1F dhist("dhist", "", nrbins, rbins);

  // k-factors file from generator lebel                                                                                                                                       
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

  // loop over ntree and dtree events isMC=true, sample 0 == Znunu, sample 1 == Zmumu,                                                                                          
  makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
  makehist4(dtree, &dhist,  true, 1, 1.00, zhists, NULL);

  string name = string("zmmcor");
  // divide the two                                                                                                                                                          
  nhist.Divide(&dhist);

  // create output file                                                                                                                                                        
  TFile outfile((name+".root").c_str(), "RECREATE");
  nhist.SetName((name+"hist" ).c_str());
  nhist.Write();
  outfile.Close();

  nfile->Close();
  dfile->Close();

  cout << "Z(mumu)->Z(inv) transfer factor computed ..." << endl;
}


// make histograms for Z->ee to signal region correction                                                                                                                       
void makezeecorhist( std::string  signalRegionFile,  std::string  zeeFile,  std::string  kFactorFile) {
  
   TFile*  nfile = TFile::Open(signalRegionFile.c_str());
   TFile*  dfile = TFile::Open(zeeFile.c_str());

   TTree* ntree = (TTree*)nfile->Get("tree");
   TTree* dtree = (TTree*)dfile->Get("tree");

   TH1F nhist("nhist", "", nrbins, rbins);
   TH1F dhist("dhist", "", nrbins, rbins);

   TFile kffile(kFactorFile.c_str());
   TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
   TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
   TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");

   znlohist->Divide(zlohist);

   vector<TH1*> ehists;
   vector<TH1*> zhists;
   zhists.push_back(znlohist);
   zhists.push_back(zewkhist);

   makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
   // sample = 3 means Zee                                                                                                                                                     
   makehist4(dtree, &dhist,  true, 3, 1.00, zhists, NULL);

   string name = string("zeecor");

   nhist.Divide(&dhist);
   TFile outfile((name+".root").c_str(), "RECREATE");
   nhist.SetName((name+"hist" ).c_str());
   nhist.Write();
   outfile.Close();

   nfile->Close();
   dfile->Close();

   cout << "Z(ee)->Z(inv) transfer factor computed ..." << endl;
}

// make histograms for W->munu to signal region correction                                                                                                                     
void makewmncorhist( std::string  signalRegionFile,  std::string  wmunuFile,  std::string  kFactorFile) {

  TFile*  nfile = TFile::Open(signalRegionFile.c_str());
  TFile*  dfile = TFile::Open(wmunuFile.c_str());

  TTree* ntree = (TTree*)nfile->Get("tree");
  TTree* dtree = (TTree*)dfile->Get("tree");

  TH1F nhist("nhist", "", nrbins, rbins);
  TH1F dhist("dhist", "", nrbins, rbins);

  TFile kffile(kFactorFile.c_str());
  TH1* znlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
  TH1*  zlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
  TH1* zewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");

  znlohist->Divide(zlohist);

  vector<TH1*> ehists;
  vector<TH1*> zhists;
  zhists.push_back(znlohist);
  zhists.push_back(zewkhist);

  makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
  // sample = 2 means Wmunu                                                                                                                                                    
  makehist4(dtree, &dhist,  true, 2, 1.00, zhists, NULL);

  string name = string("wmncor");

  nhist.Divide(&dhist);
  TFile outfile((name+".root").c_str(), "RECREATE");
  nhist.SetName((name+"hist" ).c_str());
  nhist.Write();
  outfile.Close();

  nfile->Close();
  dfile->Close();

  cout << "W(munu)->W+jets transfer factor computed ..." << endl;
}

// make ratio between Z->nunu and W->munu in the signal region                                                                                                                 
void makewzmcorhist( std::string  znunuFile,  std::string  wmunuFile,  std::string  kFactorFile) {

  TFile*  nfile = TFile::Open(znunuFile.c_str());
  TFile*  dfile = TFile::Open(wmunuFile.c_str());

  TTree* ntree = (TTree*)nfile->Get("tree");
  TTree* dtree = (TTree*)dfile->Get("tree");

  TH1F nhist("nhist", "", nrbins, rbins);
  TH1F dhist("dhist", "", nrbins, rbins);

  TFile kffile(kFactorFile.c_str());
  TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
  TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
  TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");

  TH1* wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
  TH1*  wlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
  TH1* wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");

  znlohist->Divide(zlohist);
  wnlohist->Divide(wlohist);

  vector<TH1*> zhists;
  vector<TH1*> whists;
  zhists.push_back(znlohist);
  zhists.push_back(zewkhist);
  whists.push_back(wnlohist);
  whists.push_back(wewkhist);

  // sample 0 means Znunu                                                                                                                                                      
  makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
  // sample 2 means Wmunu                                                                                                                                                      
  makehist4(dtree, &dhist,  true, 2, 1.00, whists, NULL);

  string name = string("wzmcor");
  nhist.Divide(&dhist);
  TFile outfile((name+".root").c_str(), "RECREATE");
  nhist.SetName((name+"hist" ).c_str());
  nhist.Write();
  outfile.Close();

  nfile->Close();
  dfile->Close();

  cout << "W(munu)->Z(inv) transfer factor computed ..." << endl;
}


// W+jets to Znunu correction factor in the signal region                                                                                                                      
void makezwjcorhist( std::string  znunuFile,  std::string  wlnuFile,  std::string  kFactorFile,string ext="", int kfact=0) {

  TFile*  nfile =  TFile::Open(znunuFile.c_str());
  TFile*  dfile =  TFile::Open(wlnuFile.c_str());

  TTree* ntree = (TTree*)nfile->Get("tree");
  TTree* dtree = (TTree*)dfile->Get("tree");

  TH1F nhist("nhist", "", nrbins, rbins);
  TH1F dhist("dhist", "", nrbins, rbins);

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

  // sample = 0 means signal region cuts                                                                                                                                         
  makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
  makehist4(dtree, &dhist,  true, 0, 1.00, whists, NULL);

  string name = string("zwjcor")+ext;

  nhist.Divide(&dhist);
  TFile outfile((name+".root").c_str(), "RECREATE");
  nhist.SetName((name+"hist" ).c_str());
  nhist.Write();
  outfile.Close();

  nfile->Close();
  dfile->Close();
  kffile.Close();

  cout << "W+jets->Z(inv) transfer factor computed ..." << endl;
}


// Photon+jets / Znunu signal region                                                                                                                                            
void makegamcorhist( std::string  znunuFile,  std::string  photonFile,  std::string  kFactorFile,  std::string  fPfile, string ext="", int kfact=0) {

  TFile*  nfile = TFile::Open(znunuFile.c_str());
  TFile*  dfile = TFile::Open(photonFile.c_str());

  TTree* ntree = (TTree*)nfile->Get("tree");
  TTree* dtree = (TTree*)dfile->Get("tree");

  TH1F nhist("nhist", "", nrbins, rbins);
  TH1F dhist("dhist", "", nrbins, rbins);

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

  makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
  // sample 5 means gamma+jets                                                                                                                                                  
  makehist4(dtree, &dhist,  true, 5, 1.00, ahists, NULL);

  string name = string("gamcor")+ext;
  
  nhist.Divide(&dhist);
  TFile outfile((name+".root").c_str(), "RECREATE");
  nhist.SetName((name+"hist" ).c_str());
  nhist.Write();
  outfile.Close();

  nfile->Close();
  dfile->Close();
  kffile.Close();

  cout << "gamma+jets->Z(inv) transfer factor computed ..." << endl;
}



