#include "makehist.h"
#include "makeCorrHistograms.C"

using namespace std;

// Build templates for the signal region
void sigdatamchist(TFile* outfile, string kFactorFile, int category, 		   
		   vector<string> observables,
		   string interaction, string mediatorMass = "100", string DMMass = "1", 
		   double lumi = 2.11, bool blind = false) {
  
  // Files for Znunu, Wlnu, Zll, top, qcd , diboson, signal, data
  TFile* znfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root");
  TFile* wlfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root");
  TFile* zlfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/sigfilter/sig_tree_DYJetsToLL_M-50.root");
  TFile* ttfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root");
  TFile* qcdfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/QCD/sigfilter/sig_tree_QCD.root");
  TFile* dbfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DiBoson/sigfilter/sig_tree_DiBoson.root");

  TFile* monoJfile = NULL;
  TFile* monoWfile = NULL;
  TFile* monoZfile = NULL;

  if(interaction == "Vector"){
    monoJfile = TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DMV_Vector/sigfilter/sig_tree_DMV_NNPDF30_Vector_Mphi-"+mediatorMass+"_Mchi-"+DMMass+"_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str());
    monoWfile = TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Vector/sigfilter/sig_tree_VectorMonoW_Mphi-"+mediatorMass+"_Mchi-"+DMMass+"_gSM-1p0_gDM-1p0_13TeV-madgraph.root").c_str());
    monoZfile = TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoZ_Vector/sigfilter/sig_tree_VectorMonoZ_Mphi-"+mediatorMass+"_Mchi-"+DMMass+"_gSM-1p0_gDM-1p0_13TeV-madgraph.root").c_str());
  }
  else if(interaction == "Axial"){
    monoJfile = TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DMV_Axial/sigfilter/sig_tree_DMV_NNPDF30_Axial_Mphi-"+mediatorMass+"_Mchi-"+DMMass+"_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str());
    monoWfile = TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Axial/sigfilter/sig_tree_AxialMonoW_Mphi-"+mediatorMass+"_Mchi-"+DMMass+"_gSM-1p0_gDM-1p0_13TeV-madgraph.root").c_str());
    monoZfile = TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoZ_Axial/sigfilter/sig_tree_AxialMonoZ_Mphi-"+mediatorMass+"_Mchi-"+DMMass+"_gSM-1p0_gDM-1p0_13TeV-madgraph.root").c_str());
  }
  else if(interaction == "Scalar"){
    monoJfile = TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DMS_Scalar/sigfilter/sig_tree_DMS_NNPDF30_Scalar_Mphi-"+mediatorMass+"_Mchi-"+DMMass+"_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str());
    monoWfile = TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Scalar/sigfilter/sig_tree_DM_ScalarWH_Mphi-"+mediatorMass+"_Mchi-"+DMMass+"_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str());
  }
  else if(interaction == "Pseudoscalar"){
    monoJfile = TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DMS_Pseudoscalar/sigfilter/sig_tree_DMS_NNPDF30_Pseudoscalar_Mphi-"+mediatorMass+"_Mchi-"+DMMass+"_gSM-1p0_gDM-1p0_13TeV-powheg.root").c_str());
    monoWfile = TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoW_Pseudoscalar/sigfilter/sig_tree_DM_PseudoscalarWH_Mphi-"+mediatorMass+"_Mchi-"+DMMass+"_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str());
    monoZfile = TFile::Open(("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MonoZ_Pseudoscalar/sigfilter/sig_tree_DM_PseudoscalarZH_Mphi-"+mediatorMass+"_Mchi-"+DMMass+"_gSM-1p0_gDM-1p0_13TeV-JHUGen.root").c_str());
  }
  

  //data
  TFile* dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MET/sigfilter/sig_tree_crab_MET-Run2015.root");

  // make met histograms
  vector<TH1*> znhist;
  vector<TH1*> wlhist;
  vector<TH1*> zlhist;
  vector<TH1*> tthist;
  vector<TH1*> dihist;
  vector<TH1*> qcdhist;
  vector<TH1*> monoJhist;
  vector<TH1*> monoWhist;
  vector<TH1*> monoZhist;
  vector<TH1*> dthist;

  vector<float> bins;
  
  if(category <= 1){

    for(auto obs : observables){

      if(obs == "met")
        bins = bins_monoJ;
      else
        cout<<"No binning for this observable --> please define it"<<endl;

      TH1F* znhist_temp = new TH1F(("zinvhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* wlhist_temp = new TH1F(("wjethist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* zlhist_temp = new TH1F(("zjethist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_temp = new TH1F(("tbkghist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dihist_temp = new TH1F(("dbkghist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* qcdhist_temp = new TH1F(("qbkghist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* monoJhist_temp = new TH1F(("monoJhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* monoWhist_temp = new TH1F(("monoWhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* monoZhist_temp = new TH1F(("monoZhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dthist_temp = new TH1F(("datahist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      
      znhist.push_back(dynamic_cast<TH1*>(znhist_temp));
      wlhist.push_back(dynamic_cast<TH1*>(wlhist_temp));
      zlhist.push_back(dynamic_cast<TH1*>(zlhist_temp));
      tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
      qcdhist.push_back(dynamic_cast<TH1*>(qcdhist_temp));
      dihist.push_back(dynamic_cast<TH1*>(dihist_temp));
      monoJhist.push_back(dynamic_cast<TH1*>(monoJhist_temp));
      monoWhist.push_back(dynamic_cast<TH1*>(monoWhist_temp));
      monoZhist.push_back(dynamic_cast<TH1*>(monoZhist_temp));
      dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
    }
  }
  else{

    for(auto obs : observables){

      if(obs == "met")
        bins = bins_monoV;
      else
        cout<<"No binning for this observable --> please define it"<<endl;
      
      TH1F* znhist_temp = new TH1F(("zinvhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* wlhist_temp = new TH1F(("wjethist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* zlhist_temp = new TH1F(("zjethist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_temp = new TH1F(("tbkghist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dihist_temp = new TH1F(("dbkghist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* qcdhist_temp = new TH1F(("qbkghist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* monoJhist_temp = new TH1F(("monoJhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* monoWhist_temp = new TH1F(("monoWhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* monoZhist_temp = new TH1F(("monoZhist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dthist_temp = new TH1F(("datahist_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      
      znhist.push_back(dynamic_cast<TH1*>(znhist_temp));
      wlhist.push_back(dynamic_cast<TH1*>(wlhist_temp));
      zlhist.push_back(dynamic_cast<TH1*>(zlhist_temp));
      tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
      qcdhist.push_back(dynamic_cast<TH1*>(qcdhist_temp));
      dihist.push_back(dynamic_cast<TH1*>(dihist_temp));
      monoJhist.push_back(dynamic_cast<TH1*>(monoJhist_temp));
      monoWhist.push_back(dynamic_cast<TH1*>(monoWhist_temp));
      monoZhist.push_back(dynamic_cast<TH1*>(monoZhist_temp));
      dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
    }
  }

  vector<TH2*> znhist_2D;
  vector<TH2*> wlhist_2D;
  vector<TH2*> zlhist_2D;
  vector<TH2*> tthist_2D;
  vector<TH2*> dihist_2D;
  vector<TH2*> qcdhist_2D;
  vector<TH2*> monoJhist_2D;
  vector<TH2*> monoWhist_2D;
  vector<TH2*> monoZhist_2D;

  vector<TH2*> dthist_2D;


  // get trees
  TTree* zntree = (TTree*)znfile->Get("tree/tree");
  TTree* wltree = (TTree*)wlfile->Get("tree/tree");
  TTree* zltree = (TTree*)zlfile->Get("tree/tree");
  TTree* tttree = (TTree*)ttfile->Get("tree/tree");
  TTree* ditree = (TTree*)dbfile->Get("tree/tree");
  TTree* qcdtree = (TTree*)qcdfile->Get("tree/tree");
  TTree* monoJtree = NULL;
  TTree* monoWtree = NULL;
  TTree* monoZtree = NULL;

  if(monoJfile)
    monoJtree = (TTree*) monoJfile->Get("tree/tree");
  if(monoWfile)
    monoWtree = (TTree*) monoWfile->Get("tree/tree");
  if(monoZfile)
    monoZtree = (TTree*) monoZfile->Get("tree/tree");

  TTree* dttree = (TTree*)dtfile->Get("tree");

  // get k-factors NLO
  TFile kffile(kFactorFile.c_str());
  TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
  TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
  TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");
  TH1* wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
  TH1*  wlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
  TH1* wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");

  // NLO corrections for Z and W base mc
  znlohist->Divide(zlohist);
  wnlohist->Divide(wlohist);

  vector<TH1*> ehists;
  vector<TH1*> zhists;
  vector<TH1*> whists;
  // apply EWK and QCD corrections
  zhists.push_back(znlohist); zhists.push_back(zewkhist);
  whists.push_back(wnlohist); whists.push_back(wewkhist);

  bool isWJet = false;
  if(category == 2 or category == 3)
    isWJet = true;

  // make histograms for the signal region
  makehist4(zntree, znhist,  znhist_2D,  true, 0, category, false, 1.00, lumi, zhists, "", true, NULL);
  makehist4(wltree, wlhist,  wlhist_2D,  true, 0, category, false, 1.00, lumi, whists, "", true, NULL);
  makehist4(zltree, zlhist,  zlhist_2D,  true, 0, category, false, 1.00, lumi, zhists, "", true, NULL);
  makehist4(tttree, tthist,  tthist_2D,  true, 0, category, false, 1.00, lumi, ehists, "", true, NULL);
  makehist4(ditree, dihist,  dihist_2D,  true, 0, category, isWJet, 1.00, lumi, ehists, "", true, NULL);
  makehist4(qcdtree, qcdhist,  qcdhist_2D,  true, 0, category, false, 1.00, lumi, ehists, "", true, NULL);
  // signals
  if(monoJtree)
    makehist4(monoJtree, monoJhist,  monoJhist_2D,  true, 0, category, false, 1.00, lumi, ehists, "",  true, NULL);
  if(monoWtree)
    makehist4(monoWtree, monoWhist,  monoWhist_2D,  true, 0, category, isWJet, 1.00, lumi, ehists, "", true, NULL);
  if(monoZtree)
    makehist4(monoZtree, monoZhist,  monoZhist_2D,  true, 0, category, isWJet, 1.00, lumi, ehists, "", true, NULL);
  // data
  makehist4(dttree, dthist,  dthist_2D, false, 0, category, false, 1.00, lumi, ehists, "", true, NULL);


  if (blind) {
    for( size_t ihist = 0; ihist < dthist.size(); ihist++){
      for (int i = 1; i <= dthist.at(ihist)->GetNbinsX(); i++) {
	double binval = 0.0;
	binval += znhist.at(ihist)->GetBinContent(i); 
	binval += wlhist.at(ihist)->GetBinContent(i); 
	binval += zlhist.at(ihist)->GetBinContent(i); 
	binval += tthist.at(ihist)->GetBinContent(i); 
	binval += dihist.at(ihist)->GetBinContent(i); 
	binval += qcdhist.at(ihist)->GetBinContent(i); 
	dthist.at(ihist)->SetBinContent(i, int(binval));
      }
    }

    for( size_t ihist = 0; ihist < dthist_2D.size(); ihist++){
      for (int iX = 1; iX <= dthist_2D.at(ihist)->GetNbinsX(); iX++) {
	for (int iY = 1; iY <= dthist_2D.at(ihist)->GetNbinsY(); iY++) {
	  double binval = 0.0;
	  binval += znhist_2D.at(ihist)->GetBinContent(iX,iY); 
	  binval += wlhist_2D.at(ihist)->GetBinContent(iX,iY); 
	  binval += zlhist_2D.at(ihist)->GetBinContent(iX,iY); 
	  binval += tthist_2D.at(ihist)->GetBinContent(iX,iY); 
	  binval += dihist_2D.at(ihist)->GetBinContent(iX,iY); 
	  binval += qcdhist_2D.at(ihist)->GetBinContent(iX,iY); 
	  dthist_2D.at(ihist)->SetBinContent(iX,iY,int(binval));
	}
      }
    }
  }

  outfile->cd();
  // store histograms
  for(auto hist : znhist)
    hist->Write();
  for(auto hist : wlhist)
    hist->Write();
  for(auto hist : zlhist)
    hist->Write();
  for(auto hist : tthist)
    hist->Write();
  for(auto hist : dihist)
    hist->Write();
  for(auto hist : qcdhist)
    hist->Write();
  for(auto hist : monoJhist)
    hist->Write();
  for(auto hist : monoWhist)
    hist->Write();
  for(auto hist : monoZhist)
    hist->Write();
  for(auto hist : dthist)
    hist->Write();
  // store hist_2Dograms
  for(auto hist_2D : znhist_2D)
    hist_2D->Write();
  for(auto hist_2D : wlhist_2D)
    hist_2D->Write();
  for(auto hist_2D : zlhist_2D)
    hist_2D->Write();
  for(auto hist_2D : tthist_2D)
    hist_2D->Write();
  for(auto hist_2D : dihist_2D)
    hist_2D->Write();
  for(auto hist_2D : qcdhist_2D)
    hist_2D->Write();
  for(auto hist_2D : monoJhist_2D)
    hist_2D->Write();
  for(auto hist_2D : monoWhist_2D)
    hist_2D->Write();
  for(auto hist_2D : monoZhist_2D)
    hist_2D->Write();
  for(auto hist_2D : dthist_2D)
    hist_2D->Write();
  
  znfile->Close();
  wlfile->Close();
  zlfile->Close();
  ttfile->Close();
  dbfile->Close();
  qcdfile->Close();
  if(monoJfile)
    monoJfile->Close();
  if(monoWfile)
    monoWfile->Close();
  if(monoZfile)
    monoZfile->Close();
  dtfile->Close();  
  kffile.Close();

  cout << "Templates for the signal region computed ..." << endl;
}

// build templates for photon+jets control region
void gamdatamchist(TFile* outfile, string photonFile,
		   int category, vector<string> observables,
		   double lumi = 2.11) {


  TFile* dtfile =  TFile::Open(photonFile.c_str());

  vector<TH1*> dthist;
  vector<TH1*> qcdhist;

  vector<TH2*> dthist_2D;
  vector<TH2*> qcdhist_2D;

  vector<float> bins;

  if(category <=1){
    for(auto obs : observables){

      if(obs == "met")
        bins = bins_monoJ;
      else
        cout<<"No binning for this observable --> please define it"<<endl;

      TH1F* dthist_temp = new TH1F(("datahistgam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* qchist_temp = new TH1F(("qbkghistgam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
      qcdhist.push_back(dynamic_cast<TH1*>(qchist_temp));      
    }

  }
  else{
    for(auto obs : observables){

      if(obs == "met")
        bins = bins_monoV;
      else
        cout<<"No binning for this observable --> please define it"<<endl;

      TH1F* dthist_temp = new TH1F(("datahistgam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* qchist_temp = new TH1F(("qbkghistgam_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
      qcdhist.push_back(dynamic_cast<TH1*>(qchist_temp));      
    }
  }


  TTree* dttree = (TTree*)dtfile->Get("tree");
  vector<TH1*> ehists;
  

  makehist4(dttree, dthist, dthist_2D, false, 5, category, false, 1.00, lumi, ehists, "", true, NULL);
  // QCD template for photon + jets derived from photon purity
  makehist4(dttree, qcdhist, qcdhist_2D, false, 6, category, false, 1.00, lumi, ehists, "", true, NULL);
  
  outfile->cd();

  for(auto hist : dthist)
    hist->Write();
  for(auto hist : qcdhist)
    hist->Write();
  for(auto hist : dthist_2D)
    hist->Write();
  for(auto hist : qcdhist_2D)
    hist->Write();
  
  dtfile->Close();
  
  cout << "Templates for the gamma+jets control region computed ..." << endl;
}


//build templates for Zmumu, Zee, Wenu, Wmunu
void lepdatamchist(TFile* outfile, string kFactorFile, int sample, int category, vector<string> observables, double lumi=2.11) {

  if (sample != 1 && sample != 2 && sample != 3 && sample != 4) return;

  TFile* ttfile  = NULL;
  TFile* dbfile  = NULL;
  TFile* qcfile  = NULL;
  TFile* vlfile  = NULL;
  TFile* vllfile = NULL;
  TFile* dtfile  = NULL;

  string suffix;

  if(sample == 1){

    suffix = "zmm";

    vllfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/zmmfilter/zmm_tree_DYJetsToLL_M-50.root");
    vlfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/zmmfilter/zmm_tree_WJetsToLNu.root");
    qcfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/QCD/zmmfilter/zmm_tree_QCD.root");
    dbfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DiBoson/zmmfilter/zmm_tree_DiBoson.root");
    ttfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/zmmfilter/zmm_tree_Top.root");
    
    dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MET/zmmfilter/zmm_tree_crab_MET-Run2015.root");
  }
  else if(sample == 2){

    suffix = "wmn";

    vllfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/wmnfilter/wmn_tree_DYJetsToLL_M-50.root");
    vlfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/wmnfilter/wmn_tree_WJetsToLNu.root");
    qcfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/QCD/wmnfilter/wmn_tree_QCD.root");
    dbfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DiBoson/wmnfilter/wmn_tree_DiBoson.root");
    ttfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/wmnfilter/wmn_tree_Top.root");

    dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MET/wmnfilter/wmn_tree_crab_MET-Run2015.root");
  }
  else if(sample == 3){

    suffix = "zee";

    vllfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/zeefilter/zee_tree_DYJetsToLL_M-50.root");
    vlfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/zeefilter/zee_tree_WJetsToLNu.root");
    qcfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/QCD/zeefilter/zee_tree_QCD.root");
    dbfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DiBoson/zeefilter/zee_tree_DiBoson.root");
    ttfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/zeefilter/zee_tree_Top.root");

    dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/SingleElectron/zeefilter/zee_tree_crab_SingleEle-Run2015.root");
  }
  else if(sample == 4){

    suffix = "wen";

    vllfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/wenfilter/wen_tree_DYJetsToLL_M-50.root");
    vlfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/wenfilter/wen_tree_WJetsToLNu.root");
    qcfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/QCD/wenfilter/wen_tree_QCD.root");
    dbfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DiBoson/wenfilter/wen_tree_DiBoson.root");
    ttfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/wenfilter/wen_tree_Top.root");

    dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/SingleElectron/wenfilter/wen_tree_crab_SingleEle-Run2015.root");
  }

  vector<TH1*> dthist;
  vector<TH1*> tthist;
  vector<TH1*> qchist;
  vector<TH1*> dbhist;
  vector<TH1*> vlhist;
  vector<TH1*> vllhist;

  vector<float> bins;
  
  if(category <=1){
    
    for(auto obs : observables){      

      if(obs == "met")
        bins = bins_monoJ;
      else
        cout<<"No binning for this observable --> please define it"<<endl;

      TH1F* dthist_temp = new TH1F((string("datahist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_temp = new TH1F((string("tbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_temp = new TH1F((string("dbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* qchist_temp = new TH1F((string("qbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_temp = new TH1F((string("vlbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_temp = new TH1F((string("vllbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
      tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
      dbhist.push_back(dynamic_cast<TH1*>(dbhist_temp));
      qchist.push_back(dynamic_cast<TH1*>(qchist_temp));
      vlhist.push_back(dynamic_cast<TH1*>(vlhist_temp));
      vllhist.push_back(dynamic_cast<TH1*>(vllhist_temp));

    }
  }
  else{
    for(auto obs : observables){      

      if(obs == "met")
        bins = bins_monoV;
      else
        cout<<"No binning for this observable --> please define it"<<endl;

      TH1F* dthist_temp = new TH1F((string("datahist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_temp = new TH1F((string("tbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_temp = new TH1F((string("dbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* qchist_temp = new TH1F((string("qbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_temp = new TH1F((string("vlbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_temp = new TH1F((string("vllbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
      tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
      dbhist.push_back(dynamic_cast<TH1*>(dbhist_temp));
      qchist.push_back(dynamic_cast<TH1*>(qchist_temp));
      vlhist.push_back(dynamic_cast<TH1*>(vlhist_temp));
      vllhist.push_back(dynamic_cast<TH1*>(vllhist_temp));
    }
  }

  vector<TH2*> dthist_2D;
  vector<TH2*> tthist_2D;
  vector<TH2*> qchist_2D;
  vector<TH2*> dbhist_2D;
  vector<TH2*> vlhist_2D;
  vector<TH2*> vllhist_2D;

  TTree* dttree = (TTree*)dtfile->Get("tree");
  TTree* vltree = (TTree*)vlfile->Get("tree/tree");
  TTree* vlltree = (TTree*)vllfile->Get("tree/tree");
  TTree* tttree = (TTree*)ttfile->Get("tree/tree");
  TTree* dbtree = (TTree*)dbfile->Get("tree/tree");
  TTree* qctree = (TTree*)qcfile->Get("tree/tree");

  TFile kffile(kFactorFile.c_str());
  TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
  TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
  TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");
  TH1* wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
  TH1*  wlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
  TH1* wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");
  
  znlohist->Divide(zlohist);
  wnlohist->Divide(wlohist);

  vector<TH1*> ehists;
  vector<TH1*> vlhists;
  vector<TH1*> vllhists;
  // apply NLO QCD and EWK corrections for Zll and Wlnu
  vlhists.push_back(wnlohist); 
  vlhists.push_back(wewkhist);
  vllhists.push_back(znlohist); 
  vllhists.push_back(zewkhist);

  bool isWJet = false;
  if(category == 2 or category == 3)
    isWJet = true;

  makehist4(tttree, tthist,  tthist_2D,  true,  sample, category, false,  1.00, lumi, ehists, "", true, NULL);
  makehist4(dbtree, dbhist,  dbhist_2D,  true,  sample, category, isWJet, 1.00, lumi, ehists, "", true, NULL);
  makehist4(qctree, qchist,  qchist_2D,  true,  sample, category, false,  1.00, lumi, ehists, "", true, NULL);
  makehist4(vltree, vlhist,  vlhist_2D,  true,  sample, category, false,  1.00, lumi, vlhists, "", true, NULL);
  makehist4(vlltree,vllhist, vllhist_2D, true,  sample, category, false,  1.00, lumi, vllhists, "", true, NULL);
  makehist4(dttree, dthist, dthist_2D,   false, sample, category, false,  1.00, lumi, ehists, "", true, NULL);
  
  outfile->cd();
  for(auto hist :  dthist)
    hist->Write();
  for(auto hist :  tthist)
    hist->Write();
  for(auto hist :  dbhist)
    hist->Write();
  for(auto hist :  qchist)
    hist->Write();
  for(auto hist :  vlhist)
    hist->Write();
  for(auto hist :  vllhist)
    hist->Write();

  dtfile->Close();
  vlfile->Close();
  vllfile->Close();
  ttfile->Close();
  dbfile->Close();
  qcfile->Close();
  
  cout << "Templates for the lepton control region computed ..." << endl;
}




//build templates for Zmumu, Zee, Wenu, Wmunu
void topdatamchist(TFile* outfile, string kFactorFile, int sample, int category, vector<string> observables, double lumi=2.11) {

  if (sample != 7 && sample != 8) return;

  TFile* ttfile  = NULL;
  TFile* dbfile  = NULL;
  TFile* qcfile  = NULL;
  TFile* vlfile  = NULL;
  TFile* vllfile = NULL;
  TFile* dtfile  = NULL;

  string suffix;

  if(sample == 7)    
    suffix = "topmu";
  else if(sample == 8)
    suffix = "topel";
  
  vllfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/topfilter/top_tree_DYJetsToLL_M-50.root");
  vlfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/topfilter/top_tree_WJetsToLNu.root");
  qcfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/QCD/topfilter/top_tree_QCD.root");
  dbfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DiBoson/topfilter/top_tree_DiBoson.root");
  ttfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root");
    
  if(sample == 7)
    dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/MET/topfilter/top_tree_crab_MET-Run2015.root");
  else if(sample == 8)
    dtfile = TFile::Open("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/SingleElectron/topfilter/top_tree_crab_SingleEle-Run2015.root");

  vector<TH1*> dthist;
  vector<TH1*> tthist;
  vector<TH1*> qchist;
  vector<TH1*> dbhist;
  vector<TH1*> vlhist;
  vector<TH1*> vllhist;

  vector<float> bins;
  
  if(category <=1){
    
    for(auto obs : observables){      

      if(obs == "met")
        bins = bins_monoJ;
      else
        cout<<"No binning for this observable --> please define it"<<endl;

      TH1F* dthist_temp = new TH1F((string("datahist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_temp = new TH1F((string("tbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_temp = new TH1F((string("dbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* qchist_temp = new TH1F((string("qbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_temp = new TH1F((string("vlbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_temp = new TH1F((string("vllbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
      tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
      dbhist.push_back(dynamic_cast<TH1*>(dbhist_temp));
      qchist.push_back(dynamic_cast<TH1*>(qchist_temp));
      vlhist.push_back(dynamic_cast<TH1*>(vlhist_temp));
      vllhist.push_back(dynamic_cast<TH1*>(vllhist_temp));

    }
  }
  else{
    for(auto obs : observables){      

      if(obs == "met")
        bins = bins_monoV;
      else
        cout<<"No binning for this observable --> please define it"<<endl;

      TH1F* dthist_temp = new TH1F((string("datahist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* tthist_temp = new TH1F((string("tbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* dbhist_temp = new TH1F((string("dbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* qchist_temp = new TH1F((string("qbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vlhist_temp = new TH1F((string("vlbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);
      TH1F* vllhist_temp = new TH1F((string("vllbkghist")+suffix+"_"+obs).c_str(), "", int(bins.size()-1), &bins[0]);

      dthist.push_back(dynamic_cast<TH1*>(dthist_temp));
      tthist.push_back(dynamic_cast<TH1*>(tthist_temp));
      dbhist.push_back(dynamic_cast<TH1*>(dbhist_temp));
      qchist.push_back(dynamic_cast<TH1*>(qchist_temp));
      vlhist.push_back(dynamic_cast<TH1*>(vlhist_temp));
      vllhist.push_back(dynamic_cast<TH1*>(vllhist_temp));
    }
  }

  vector<TH2*> dthist_2D;
  vector<TH2*> tthist_2D;
  vector<TH2*> qchist_2D;
  vector<TH2*> dbhist_2D;
  vector<TH2*> vlhist_2D;
  vector<TH2*> vllhist_2D;

  TTree* dttree = (TTree*)dtfile->Get("tree");
  TTree* vltree = (TTree*)vlfile->Get("tree/tree");
  TTree* vlltree = (TTree*)vllfile->Get("tree/tree");
  TTree* tttree = (TTree*)ttfile->Get("tree/tree");
  TTree* dbtree = (TTree*)dbfile->Get("tree/tree");
  TTree* qctree = (TTree*)qcfile->Get("tree/tree");
  
  vector<TH1*> ehists;

  bool isWJet = false;
  if(category == 2 or category == 3)
    isWJet = true;

  makehist4(tttree,  tthist,  tthist_2D,  true,  sample, category, false,  1.00, lumi, ehists, "", true, NULL);
  makehist4(dbtree,  dbhist,  dbhist_2D,  true,  sample, category, isWJet, 1.00, lumi, ehists, "", true, NULL);
  makehist4(qctree,  qchist,  qchist_2D,  true,  sample, category, false,  1.00, lumi, ehists, "", true, NULL);
  makehist4(vltree,  vlhist,  vlhist_2D,  true,  sample, category, false,  1.00, lumi, ehists, "", true, NULL);
  makehist4(vlltree, vllhist, vllhist_2D, true,  sample, category, false,  1.00, lumi, ehists, "", true, NULL);
  makehist4(dttree,  dthist,  dthist_2D,  false, sample, category, false,  1.00, lumi, ehists, "", true, NULL);
  
  outfile->cd();
  for(auto hist :  dthist)
    hist->Write();
  for(auto hist :  tthist)
    hist->Write();
  for(auto hist :  dbhist)
    hist->Write();
  for(auto hist :  qchist)
    hist->Write();
  for(auto hist :  vlhist)
    hist->Write();
  for(auto hist :  vllhist)
    hist->Write();

  dtfile->Close();
  vlfile->Close();
  vllfile->Close();
  ttfile->Close();
  dbfile->Close();
  qcfile->Close();
  
  cout << "Templates for the top control region computed ..." << endl;
}


// Run the final analysis:
// 1) Store all corrections templates from input files (complient to combine)
// 2) Make data and expected yields templates for all the other processes

void hists(bool doCorrectionHistograms = false, int category = 0, double lumi = 2.11, string outDir = "", string ext ="") {

  string kfactorFile         = "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/kFactors/scalefactors_v4.root";
  vector<string> observables = {"met"};
  system(("mkdir -p "+outDir).c_str());

  if(doCorrectionHistograms){

    cout<<"make correction histogram for Zmm to Znn"<<endl;
    // make central values
    makezmmcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
    		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/zmmfilter/zmm_tree_DYJetsToLL_M-50.root",
    		   kfactorFile,
    		   category,	
		   observables,
		   lumi,
		   outDir,
    		   ext); 
    cout<<"make correction histogram for Zee to Znn"<<endl;
    makezeecorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/DYJets/zeefilter/zee_tree_DYJetsToLL_M-50.root",
		   kfactorFile,
                   category,
		   observables,
		   lumi,
		   outDir,
		   ext);

    cout<<"make correction histogram for Wmn to WJets"<<endl;
    makewmncorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/wmnfilter/wmn_tree_WJetsToLNu.root",
                   kfactorFile,
                   category,
		   observables,
		   lumi,
		   outDir,
		   ext);

    cout<<"make correction histogram for Wen to WJets"<<endl;
    makewencorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/wenfilter/wen_tree_WJetsToLNu.root",
                   kfactorFile,
                   category,
		   observables,
		   lumi,
		   outDir,
    		   ext);
   
    cout<<"make correction histogram for Gam+jets to Znn"<<endl;
    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root", 
		   kfactorFile,
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,
		   observables,
		   lumi,
		   outDir,
		   ext);

    cout<<"make Z/W ratio"<<endl;
    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   kfactorFile,
                   category,
		   observables,
		   lumi,
		   outDir,
		   ext);

    cout<<"make top mu ratio"<<endl;
    maketopmucorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root",
		     category,
		     observables,
		     lumi,
		     outDir,
		     ext);

    cout<<"make top el ratio"<<endl;
    maketopelcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root",
		     category,
		     observables,
		     lumi,
		     outDir,
		     ext);

    if(category == 2 or category == 3){

      cout<<"make sideband ratio for Zvv "<<std::endl;      
      makesidebandcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
			  "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
			  category,
			  category+2,
			  observables,
			  lumi,
			  outDir,
			  ext+"Z");

      cout<<"make sideband ratio for W+jets "<<std::endl;      
      makesidebandcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
			  "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
			  category,
			  category+2,
			  observables,
			  lumi,
			  outDir,
			  ext+"W");
    }
    
    // systematics
    cout<<"systematics on Z/gamma ratio "<<endl;
    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",
		   kfactorFile,
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,
		   observables,
		   lumi,
		   outDir,
		   "qcd"+ext,1);

    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",
		   kfactorFile,
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,
		   observables,
		   lumi,
		   outDir,
		   "ewk"+ext,2);
    
    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",
		   kfactorFile,
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,
		   observables,
		   lumi,
		   outDir,
		   "re1"+ext,3);

    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",
		   kfactorFile,
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,
		   observables,
		   lumi,
		   outDir,
		   "fa1"+ext,4);
    

    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",
		   kfactorFile,
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,
		   observables,
		   lumi,
		   outDir,
		   "re2"+ext,5);

    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",
		   kfactorFile,
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,
		   observables,
		   lumi,
		   outDir,
		   "fa2"+ext,6);

    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",
		   kfactorFile,
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,
		   observables,
		   lumi,
		   outDir,
		   "pdf"+ext,7);

    makegamcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/PhotonJets/gamfilter/gam_tree_GJets.root",
		   kfactorFile,
		   "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		   category,
		   observables,
		   lumi,
		   outDir,
		   "fpc"+ext,8);
    
    //
    cout<<"systematics on Z/W ratio "<<endl;
    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   kfactorFile,
                   category,
		   observables,
		   lumi,
		   outDir,
		   "qcd"+ext,1);

    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   kfactorFile,
                   category,
		   observables,
		   lumi,
		   outDir,
		   "ewk"+ext,2);

    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   kfactorFile,
                   category,
		   observables,
		   lumi,
		   outDir,
		   "re1"+ext,3);


    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   kfactorFile,
                   category,
		   observables,
		   lumi,
		   outDir,
		   "fa1"+ext,4);


    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   kfactorFile,
                   category,
		   observables,
		   lumi,
		   outDir,
		   "re2"+ext,5);


    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   kfactorFile,
                   category,
		   observables,
		   lumi,
		   outDir,
		   "fa2"+ext,6);

    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   kfactorFile,
                   category,
		   observables,
		   lumi,
		   outDir,
		   "fa2"+ext,6);

    makezwjcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/ZJets/sigfilter/sig_tree_ZJetsToNuNu.root",
		   "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/WJets/sigfilter/sig_tree_WJetsToLNu.root",
		   kfactorFile,
                   category,
		   observables,
		   lumi,
		   outDir,
		   "pdf"+ext,7);

    cout<<"systematics on top+mu ratio"<<endl;
    maketopmucorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root",
		     category,
		     observables,
		     lumi,
		     outDir,
		     "btagUp",
		     ext+"bUp");


    maketopmucorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root",
		     category,
		     observables,
		     lumi,
		     outDir,
		     "btagDown",
		     ext+"bDown");

    cout<<"systematics on top+el ratio"<<endl;
    maketopelcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root",
		     category,
		     observables,
		     lumi,
		     outDir,
		     "btagUp",		     
		     ext+"bUp");

    maketopelcorhist("/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/sigfilter/sig_tree_Top.root",
		     "/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/Top/topfilter/top_tree_Top.root",
		     category,
		     observables,
		     lumi,
		     outDir,
		     "btagDown",		     
  		     ext+"bDown");
  }
  
  // take correction files --> central value
  cout<<"Re-open file for correction histo"<<endl;
  TFile* zmmcorfile = TFile::Open((outDir+"/zmmcor"+ext+".root").c_str());
  TFile* zeecorfile = TFile::Open((outDir+"/zeecor"+ext+".root").c_str());
  TFile* wmncorfile = TFile::Open((outDir+"/wmncor"+ext+".root").c_str());
  TFile* wencorfile = TFile::Open((outDir+"/wencor"+ext+".root").c_str());
  TFile* zwjcorfile = TFile::Open((outDir+"/zwjcor"+ext+".root").c_str());
  TFile* gamcorfile = TFile::Open((outDir+"/gamcor"+ext+".root").c_str());
  TFile* topmucorfile = TFile::Open((outDir+"/topmucor"+ext+".root").c_str());
  TFile* topelcorfile = TFile::Open((outDir+"/topelcor"+ext+".root").c_str());

  TFile* sidebandfileZ = NULL;
  TFile* sidebandfileW = NULL;
  if(category == 2 or category == 3){
    sidebandfileZ = TFile::Open((outDir+"/sidebandcor"+ext+"Z.root").c_str());
    sidebandfileW = TFile::Open((outDir+"/sidebandcor"+ext+"W.root").c_str());
  }

  // QCD, EWK, factm re and footprint on Z/gamma
  TFile* gamcorqcdfile = TFile::Open((outDir+"/gamcorqcd"+ext+".root").c_str());
  TFile* gamcorewkfile = TFile::Open((outDir+"/gamcorewk"+ext+".root").c_str());
  TFile* gamcorre1file = TFile::Open((outDir+"/gamcorre1"+ext+".root").c_str());
  TFile* gamcorfa1file = TFile::Open((outDir+"/gamcorfa1"+ext+".root").c_str());
  TFile* gamcorre2file = TFile::Open((outDir+"/gamcorre2"+ext+".root").c_str());
  TFile* gamcorfa2file = TFile::Open((outDir+"/gamcorfa2"+ext+".root").c_str());
  TFile* gamcorpdffile = TFile::Open((outDir+"/gamcorpdf"+ext+".root").c_str());
  TFile* gamcorfpcfile = TFile::Open((outDir+"/gamcorfpc"+ext+".root").c_str());

  // QCD, EWK, factm re and footprint on Z/W
  TFile* zwjcorqcdfile = TFile::Open((outDir+"/zwjcorqcd"+ext+".root").c_str());
  TFile* zwjcorewkfile = TFile::Open((outDir+"/zwjcorewk"+ext+".root").c_str());
  TFile* zwjcorre1file = TFile::Open((outDir+"/zwjcorre1"+ext+".root").c_str());
  TFile* zwjcorfa1file = TFile::Open((outDir+"/zwjcorfa1"+ext+".root").c_str());
  TFile* zwjcorre2file = TFile::Open((outDir+"/zwjcorre2"+ext+".root").c_str());
  TFile* zwjcorfa2file = TFile::Open((outDir+"/zwjcorfa2"+ext+".root").c_str());
  TFile* zwjcorpdffile = TFile::Open((outDir+"/zwjcorpdf"+ext+".root").c_str());

  // Top btag up and down
  TFile* topmucorbupfile   = TFile::Open((outDir+"/topmucor"+ext+"bUp.root").c_str());
  TFile* topmucorbdownfile = TFile::Open((outDir+"/topmucor"+ext+"bDown.root").c_str());
  TFile* topelcorbupfile   = TFile::Open((outDir+"/topelcor"+ext+"bUp.root").c_str());
  TFile* topelcorbdownfile = TFile::Open((outDir+"/topelcor"+ext+"bDown.root").c_str());

  // get histograms  
  vector<TH1*> zmmcorhist;
  vector<TH1*> zeecorhist;
  vector<TH1*> wmncorhist;
  vector<TH1*> wencorhist;
  vector<TH1*> zwjcorhist;
  vector<TH1*> gamcorhist;
  vector<TH1*> topmucorhist;
  vector<TH1*> topelcorhist;

  vector<TH1*> gamcorewkhist;
  vector<TH1*> gamcorqcdhist;
  vector<TH1*> gamcorre1hist;
  vector<TH1*> gamcorfa1hist;
  vector<TH1*> gamcorre2hist;
  vector<TH1*> gamcorfa2hist;
  vector<TH1*> gamcorpdfhist;
  vector<TH1*> gamcorfpchist;

  vector<TH1*> zwjcorewkhist;
  vector<TH1*> zwjcorqcdhist;
  vector<TH1*> zwjcorre1hist;
  vector<TH1*> zwjcorre2hist;
  vector<TH1*> zwjcorfa1hist;
  vector<TH1*> zwjcorfa2hist;
  vector<TH1*> zwjcorpdfhist;

  vector<TH1F*> sidebandZhist;
  vector<TH1F*> sidebandWhist;

  vector<TH1*> topmucorbuphist;
  vector<TH1*> topmucorbdownhist;
  vector<TH1*> topelcorbuphist;
  vector<TH1*> topelcorbdownhist;

  for(auto obs : observables){

    cout<<"Get histograms for observable "<<obs<<endl;

    zmmcorhist .push_back( (TH1*)zmmcorfile->Get(("zmmcor"+ext+"hist_"+obs).c_str()));    
    zeecorhist .push_back( (TH1*)zeecorfile->Get(("zeecor"+ext+"hist_"+obs).c_str()));    
    wmncorhist .push_back( (TH1*)wmncorfile->Get(("wmncor"+ext+"hist_"+obs).c_str()));    
    wencorhist .push_back( (TH1*)wencorfile->Get(("wencor"+ext+"hist_"+obs).c_str()));    
    zwjcorhist .push_back( (TH1*)zwjcorfile->Get(("zwjcor"+ext+"hist_"+obs).c_str()));    
    gamcorhist .push_back( (TH1*)gamcorfile->Get(("gamcor"+ext+"hist_"+obs).c_str()));    
    topmucorhist .push_back( (TH1*)topmucorfile->Get(("topmucor"+ext+"hist_"+obs).c_str()));    
    topelcorhist .push_back( (TH1*)topelcorfile->Get(("topelcor"+ext+"hist_"+obs).c_str()));    

    if(category == 2 or category == 3){
      sidebandZhist.push_back((TH1F*) sidebandfileZ->Get(("sidebandcor"+ext+"Zhist_"+obs).c_str()));
      sidebandWhist.push_back((TH1F*) sidebandfileW->Get(("sidebandcor"+ext+"Whist_"+obs).c_str()));
    }
    

    // get histograms Z/gamma
    cout<<"Make Z/gamma sys histograms"<<endl;
    gamcorewkhist .push_back( (TH1*)gamcorewkfile->Get(("gamcor"+ext+"ewkhist_"+obs).c_str()));    
    gamcorqcdhist .push_back( (TH1*)gamcorqcdfile->Get(("gamcor"+ext+"qcdhist_"+obs).c_str()));    
    gamcorre1hist .push_back( (TH1*)gamcorre1file->Get(("gamcor"+ext+"re1hist_"+obs).c_str()));    
    gamcorfa1hist .push_back( (TH1*)gamcorfa1file->Get(("gamcor"+ext+"fa1hist_"+obs).c_str()));    
    gamcorre2hist .push_back( (TH1*)gamcorre2file->Get(("gamcor"+ext+"re2hist_"+obs).c_str()));    
    gamcorfa2hist .push_back( (TH1*)gamcorfa2file->Get(("gamcor"+ext+"fa2hist_"+obs).c_str()));    
    gamcorpdfhist .push_back( (TH1*)gamcorpdffile->Get(("gamcor"+ext+"pdfhist_"+obs).c_str()));    
    gamcorfpchist .push_back( (TH1*)gamcorfpcfile->Get(("gamcor"+ext+"fpchist_"+obs).c_str()));    
    
    // uncertainty histogram for combine
    TH1* gamuncewkhist = (TH1*)gamcorewkhist.back()->Clone(("gamuncewk"+ext+"hist_"+obs).c_str());    
    gamuncewkhist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncewkhist->GetNbinsX(); i++) 
      gamuncewkhist->SetBinContent(i, fabs(gamuncewkhist->GetBinContent(i)-1.0));
    gamuncewkhist->SetName("ZG_EWK");

    TH1* gamuncre1hist = (TH1*)gamcorre1hist.back()->Clone(("gamuncre1"+ext+"hist_"+obs).c_str());    
    gamuncre1hist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncre1hist->GetNbinsX(); i++) 
      gamuncre1hist->SetBinContent(i, fabs(gamuncre1hist->GetBinContent(i)-1.0));
    gamuncre1hist->SetName("ZG_RenScale1");
    
    TH1* gamuncfa1hist = (TH1*)gamcorfa1hist.back()->Clone(("gamuncfa1"+ext+"hist_"+obs).c_str());    
    gamuncfa1hist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncfa1hist->GetNbinsX(); i++) 
      gamuncfa1hist->SetBinContent(i, fabs(gamuncfa1hist->GetBinContent(i)-1.0));
    gamuncfa1hist->SetName("ZG_FactScale1");
    
    TH1* gamuncre2hist = (TH1*)gamcorre2hist.back()->Clone(("gamuncre2"+ext+"hist_"+obs).c_str());    
    gamuncre2hist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncre2hist->GetNbinsX(); i++) 
      gamuncre2hist->SetBinContent(i, fabs(gamuncre2hist->GetBinContent(i)-1.0));
    gamuncre2hist->SetName("ZG_RenScale2");
    
    TH1* gamuncfa2hist = (TH1*)gamcorfa2hist.back()->Clone(("gamuncfa2"+ext+"hist_"+obs).c_str());    
    gamuncfa2hist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncfa2hist->GetNbinsX(); i++) 
      gamuncfa2hist->SetBinContent(i, fabs(gamuncfa2hist->GetBinContent(i)-1.0));
    gamuncfa2hist->SetName("ZG_FactScale2");

    TH1* gamuncpdfhist = (TH1*)gamcorpdfhist.back()->Clone(("gamuncpdf"+ext+"hist_"+obs).c_str());    
    gamuncpdfhist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncpdfhist->GetNbinsX(); i++) 
      gamuncpdfhist->SetBinContent(i, fabs(gamuncpdfhist->GetBinContent(i)-1.0));
    gamuncpdfhist->SetName("ZG_PDF");
    
    TH1* gamuncfpchist = (TH1*)gamcorfpchist.back()->Clone(("gamuncfpc"+ext+"hist_"+obs).c_str());    
    gamuncfpchist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncfpchist->GetNbinsX(); i++) 
      gamuncfpchist->SetBinContent(i, fabs(gamuncfpchist->GetBinContent(i)-1.0));
    gamuncfpchist->SetName("ZG_Footprint");
    
    // Same thing for Z/W ratio
    cout<<"Make Z/W sys histograms"<<endl;
    zwjcorewkhist .push_back( (TH1*)zwjcorewkfile->Get(("zwjcorewk"+ext+"hist_"+obs).c_str()));    
    zwjcorqcdhist .push_back( (TH1*)zwjcorqcdfile->Get(("zwjcorqcd"+ext+"hist_"+obs).c_str()));    
    zwjcorre1hist .push_back( (TH1*)zwjcorre1file->Get(("zwjcorre1"+ext+"hist_"+obs).c_str()));    
    zwjcorfa1hist .push_back( (TH1*)zwjcorfa1file->Get(("zwjcorfa1"+ext+"hist_"+obs).c_str()));    
    zwjcorre2hist .push_back( (TH1*)zwjcorre2file->Get(("zwjcorre2"+ext+"hist_"+obs).c_str()));    
    zwjcorfa2hist .push_back( (TH1*)zwjcorfa2file->Get(("zwjcorfa2"+ext+"hist_"+obs).c_str()));    
    zwjcorpdfhist .push_back( (TH1*)zwjcorpdffile->Get(("zwjcorpdf"+ext+"hist_"+obs).c_str()));    

    TH1* zwjuncewkhist = (TH1*)zwjcorewkhist.back()->Clone(("zwjuncewk"+ext+"hist_"+obs).c_str());    
    zwjuncewkhist->Divide(zwjcorqcdhist.back());
    for (int i = 1; i <= zwjuncewkhist->GetNbinsX(); i++) 
      zwjuncewkhist->SetBinContent(i, fabs(zwjuncewkhist->GetBinContent(i)-1.0));
    zwjuncewkhist->SetName("ZW_EWK");
    
    TH1* zwjuncre1hist = (TH1*)zwjcorre1hist.back()->Clone(("zwjuncre1"+ext+"hist_"+obs).c_str());
    zwjuncre1hist->Divide(zwjcorqcdhist.back());
    for (int i = 1; i <= zwjuncre1hist->GetNbinsX(); i++) 
      zwjuncre1hist->SetBinContent(i, fabs(zwjuncre1hist->GetBinContent(i)-1.0));
    zwjuncre1hist->SetName("ZW_RenScale1");

    TH1* zwjuncfa1hist = (TH1*)zwjcorfa1hist.back()->Clone(("zwjuncfa1"+ext+"hist_"+obs).c_str());
    zwjuncfa1hist->Divide(zwjcorqcdhist.back());
    for (int i = 1; i <= zwjuncfa1hist->GetNbinsX(); i++) 
      zwjuncfa1hist->SetBinContent(i, fabs(zwjuncfa1hist->GetBinContent(i)-1.0));
    zwjuncfa1hist->SetName("ZW_FactScale1");
    
    TH1* zwjuncre2hist = (TH1*)zwjcorre2hist.back()->Clone(("zwjuncre2"+ext+"hist_"+obs).c_str());
    zwjuncre2hist->Divide(zwjcorqcdhist.back());
    for (int i = 1; i <= zwjuncre2hist->GetNbinsX(); i++) 
    zwjuncre2hist->SetBinContent(i, fabs(zwjuncre2hist->GetBinContent(i)-1.0));
    zwjuncre2hist->SetName("ZW_RenScale2");
    
    TH1* zwjuncfa2hist = (TH1*)zwjcorfa2hist.back()->Clone(("zwjuncfa2"+ext+"hist_"+obs).c_str());
    zwjuncfa2hist->Divide(zwjcorqcdhist.back());
    for (int i = 1; i <= zwjuncfa2hist->GetNbinsX(); i++) 
      zwjuncfa2hist->SetBinContent(i, fabs(zwjuncfa2hist->GetBinContent(i)-1.0));
    zwjuncfa2hist->SetName("ZW_FactScale2");
    
    TH1* zwjuncpdfhist = (TH1*)zwjcorpdfhist.back()->Clone(("zwjuncpdf"+ext+"hist_"+obs).c_str());
    zwjuncpdfhist->Divide(zwjcorqcdhist.back());
    for (int i = 1; i <= zwjuncpdfhist->GetNbinsX(); i++) 
      zwjuncpdfhist->SetBinContent(i, fabs(zwjuncpdfhist->GetBinContent(i)-1.0));
    zwjuncpdfhist->SetName("ZW_PDF");

    // make b-tagging top
    cout<<"Make top sys histograms"<<endl;
    topmucorbuphist .push_back( (TH1*) topmucorbupfile->Get(("topmucor"+ext+"bUphist_"+obs).c_str()));    
    topmucorbdownhist .push_back( (TH1*) topmucorbdownfile->Get(("topmucor"+ext+"bDownhist_"+obs).c_str()));    
    topelcorbuphist .push_back( (TH1*) topelcorbupfile->Get(("topelcor"+ext+"bUphist_"+obs).c_str()));    
    topelcorbdownhist .push_back( (TH1*) topelcorbdownfile->Get(("topelcor"+ext+"bDownhist_"+obs).c_str()));    

    // make symmetrization
    TH1* topmucorbuphist_tmp = (TH1*) topmucorbuphist.back()->Clone(("topmucorbup_tmp"+ext+"hist_"+obs).c_str());    
    topmucorbuphist_tmp->Divide(topmucorhist.back());
    TH1* topmucorbdownhist_tmp = (TH1*) topmucorbdownhist.back()->Clone(("topmucorbdown_tmp"+ext+"hist_"+obs).c_str());    
    topmucorbdownhist_tmp->Divide(topmucorhist.back());
    
    TH1* topmucoruncbhist = (TH1*) topmucorhist.back()->Clone(("topmucoruncbhist"+ext+"hist_"+obs).c_str());    
    for (int i = 1; i <= topmucoruncbhist->GetNbinsX(); i++) {      
      topmucoruncbhist->SetBinContent(i, fabs(fabs(topmucorbuphist_tmp->GetBinContent(i)+topmucorbdownhist_tmp->GetBinContent(i))/2-1.0));
    }
    topmucoruncbhist->SetName("TOP_MU_B");

    // make symmetrization
    TH1* topelcorbuphist_tmp = (TH1*) topelcorbuphist.back()->Clone(("topelcorbup_tmp"+ext+"hist_"+obs).c_str());    
    topelcorbuphist_tmp->Divide(topelcorhist.back());
    TH1* topelcorbdownhist_tmp = (TH1*) topelcorbdownhist.back()->Clone(("topelcorbdown_tmp"+ext+"hist_"+obs).c_str());    
    topelcorbdownhist_tmp->Divide(topelcorhist.back());
    
    TH1* topelcoruncbhist = (TH1*) topelcorhist.back()->Clone(("topelcoruncbhist"+ext+"hist_"+obs).c_str());    
    for (int i = 1; i <= topelcoruncbhist->GetNbinsX(); i++) {      
      topelcoruncbhist->SetBinContent(i, fabs(fabs(topelcorbuphist_tmp->GetBinContent(i)+topelcorbdownhist_tmp->GetBinContent(i))/2-1.0));
    }
    topelcoruncbhist->SetName("TOP_EL_B");

    // output file
    TFile outfile((outDir+"/templates"+ext+"_"+obs+".root").c_str(), "RECREATE");
    outfile.cd();

    cout<<"Save transfer factor"<<endl;
    zmmcorhist.back()->Write();
    zeecorhist.back()->Write();
    wmncorhist.back()->Write();
    wencorhist.back()->Write();
    zwjcorhist.back()->Write();
    gamcorhist.back()->Write();
    topmucorhist.back()->Write();
    topelcorhist.back()->Write();

    if(category == 2 or category == 3){
      sidebandZhist.back()->Write();
      sidebandWhist.back()->Write();
    }
    
    gamcorqcdhist.back()->Write();
    gamcorewkhist.back()->Write();
    gamcorre1hist.back()->Write();
    gamcorfa1hist.back()->Write();
    gamcorre2hist.back()->Write();
    gamcorfa2hist.back()->Write();
    gamcorpdfhist.back()->Write();
    gamcorfpchist.back()->Write();

    gamuncewkhist->Write();
    gamuncre1hist->Write();
    gamuncfa1hist->Write();
    gamuncre2hist->Write();
    gamuncfa2hist->Write();
    gamuncpdfhist->Write();
    gamuncfpchist->Write();
    
    zwjcorqcdhist.back()->Write();
    zwjcorewkhist.back()->Write();
    zwjcorre1hist.back()->Write();
    zwjcorfa1hist.back()->Write();
    zwjcorre2hist.back()->Write();
    zwjcorfa2hist.back()->Write();
    zwjcorpdfhist.back()->Write();

    zwjuncewkhist->Write();
    zwjuncre1hist->Write();
    zwjuncfa1hist->Write();
    zwjuncre2hist->Write();
    zwjuncfa2hist->Write();
    zwjuncpdfhist->Write();

    topmucoruncbhist->Write();
    topelcoruncbhist->Write();

    // signal region templates
    cout<<"start signal region data"<<endl;
    sigdatamchist(&outfile,kfactorFile,category,observables,"Vector","1000","100",lumi,false);
    // gamma + jets
    cout<<"start gamma+jets region data"<<endl;
    gamdatamchist(&outfile,"/home/rgerosa/MONOJET_ANALYSIS/Production-24-1-2016/SinglePhoton/gamfilter/gam_tree_crab_SinglePhoton-Run2015.root",
		  category,observables,
		  lumi);

    // lepton control regions
    cout<<"start zmumu region data"<<endl;
    lepdatamchist(&outfile,kfactorFile,1,category,observables,lumi); 
    cout<<"start wmunu region data"<<endl;
    lepdatamchist(&outfile,kfactorFile,2,category,observables,lumi); 
    cout<<"start zee region data"<<endl;
    lepdatamchist(&outfile,kfactorFile,3,category,observables,lumi); 
    cout<<"start wenu region data"<<endl;
    lepdatamchist(&outfile,kfactorFile,4,category,observables,lumi);     
    // top control regions
    cout<<"start top+mu region data"<<endl;
    topdatamchist(&outfile,kfactorFile,7,category,observables,lumi);
    cout<<"start Top+el region data"<<endl;
    topdatamchist(&outfile,kfactorFile,8,category,observables,lumi);

    outfile.Close();
  }
}
