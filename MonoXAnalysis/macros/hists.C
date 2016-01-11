#include "makehist.h"

// Build templates for the signal region
void sigdatamchist(TFile* outfile, bool blind = false, std::string kFactorFile) {

  // Files for Znunu, Wlnu, Zll, top, qcd , diboson, signal, data
  TFile* znfile = TFile::Open("/Users/avartak/CMS/MonoX/FinalTrees/znn100toinf/sigtree.root");
  TFile* wlfile = TFile::Open("/Users/avartak/CMS/MonoX/FinalTrees/wln100toinf/sigtree.root");
  TFile* zlfile = TFile::Open("/Users/avartak/CMS/MonoX/FinalTrees/zll100toinf/sigtree.root");
  TFile* ttfile = TFile::Open("/Users/avartak/CMS/MonoX/FinalTrees/top/sigtree.root");
  TFile* difile = TFile::Open("/Users/avartak/CMS/MonoX/FinalTrees/qcd/sigtree.root");
  TFile* qcfile = TFile::Open("/Users/avartak/CMS/MonoX/FinalTrees/dibosons/sigtree.root");
  TFile* sifile = TFile::Open("/Users/avartak/CMS/MonoX/FinalTrees/psM200m1/sigtree.root");
  //data
  TFile* dtfile = TFile::Open("/Users/avartak/CMS/MonoX/FinalTrees/met/sigtree.root");

  // make met histograms
  TH1F znhist("zinvhist", "", nbins, bins);
  TH1F wlhist("wjethist", "", nbins, bins);
  TH1F zlhist("zjethist", "", nbins, bins);
  TH1F tthist("tbkghist", "", nbins, bins);
  TH1F dihist("dbkghist", "", nbins, bins);
  TH1F qchist("qbkghist", "", nbins, bins);
  TH1F sihist("sig1hist", "", nbins, bins);
  TH1F dthist("datahist", "", nbins, bins);
  
  // get trees
  TTree* zntree = (TTree*)znfile->Get("tree");
  TTree* wltree = (TTree*)wlfile->Get("tree");
  TTree* zltree = (TTree*)zlfile->Get("tree");
  TTree* tttree = (TTree*)ttfile->Get("tree");
  TTree* ditree = (TTree*)difile->Get("tree");
  TTree* qctree = (TTree*)qcfile->Get("tree");
  TTree* sitree = (TTree*)sifile->Get("tree");
  TTree* dttree = (TTree*)dtfile->Get("tree");

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

  // make histograms for the signal region
  makehist4(zntree, &znhist,  true, 0, 1.00, zhists, NULL);
  makehist4(wltree, &wlhist,  true, 0, 1.00, whists, NULL);
  makehist4(zltree, &zlhist,  true, 0, 1.00, zhists, NULL);
  makehist4(tttree, &tthist,  true, 0, 1.00, ehists, NULL);
  makehist4(ditree, &dihist,  true, 0, 1.00, ehists, NULL);
  makehist4(qctree, &qchist,  true, 0, 1.00, ehists, NULL);
  makehist4(sitree, &sihist,  true, 0, 1.00, ehists, NULL);
  makehist4(dttree, &dthist, false, 0, 1.00, ehists, NULL);

  if (blind) {
    for (int i = 1; i <= dthist.GetNbinsX(); i++) {
      double binval = 0.0;
            binval += znhist.GetBinContent(i); 
            binval += wlhist.GetBinContent(i); 
            binval += zlhist.GetBinContent(i); 
            binval += tthist.GetBinContent(i); 
            binval += dihist.GetBinContent(i); 
            binval += qchist.GetBinContent(i); 
            dthist.SetBinContent(i, int(binval));
    }
    }
  
  outfile->cd();
  znhist.Write();
  wlhist.Write();
  zlhist.Write();
  tthist.Write();
  dihist.Write();
  qchist.Write();
  sihist.Write();
  dthist.Write();
  
  znfile->Close();
  wlfile->Close();
  zlfile->Close();
  ttfile->Close();
  difile->Close();
  qcfile->Close();
  sifile->Close();
  dtfile->Close();
  
  kffile.Close();

  cout << "Templates for the signal region computed ..." << endl;
}

// build templates for photon+jets control region
void gamdatamchist(TFile* outfile, std::string photonFile) {

  TFile dtfile =  TFile::Open(photonFile.c_str());

  TH1F dthist("datahistgam", "", nbins, bins);
  TH1F qchist("qbkghistgam", "", nbins, bins);
  
  TTree* dttree = (TTree*)dtfile->Get("tree");

  vector<TH1*> ehists;
  
  makehist4(dttree, &dthist, false, 5, 1.00, ehists, NULL);
  // QCD template for photon + jets derived from photon purity
  makehist4(dttree, &qchist, false, 6, 1.00, ehists, NULL);

  outfile->cd();
  dthist.Write();
  qchist.Write();

  dtfile->Close();
  
  cout << "Templates for the gamma+jets control region computed ..." << endl;
}

//build templates for Zmumu, Zee, Wenu, Wmunu
void lepdatamchist(TFile* outfile, int sample, std::string kFactorFile) {

  if (sample != 1 && sample != 2 && sample != 3 && sample != 4) return;
  
  string filename;
  if      (sample == 1)   filename = "zmmtree.root";
  else if (sample == 2)   filename = "wmntree.root";
  else if (sample == 3)   filename = "zeetree.root";
  else if (sample == 4)   filename = "wentree.root";

  string vlfilename;
  if      (sample == 1) vlfilename = "wln100toinf/";
  else if (sample == 2) vlfilename = "zll100toinf/";
  else if (sample == 3) vlfilename = "wln100toinf/";
  else if (sample == 4) vlfilename = "zll100toinf/";
  
  string dtfilename;
  if      (sample == 1) dtfilename = "met/";
  else if (sample == 2) dtfilename = "met/";
  else if (sample == 3) dtfilename = "singleel/";
  else if (sample == 4) dtfilename = "singleel/";
  
  string suffix;
  if      (sample == 1) suffix = "zmm";
  else if (sample == 2) suffix = "wmn";
  else if (sample == 3) suffix = "zee";
  else if (sample == 4) suffix = "wen";

  TFile* ttfile = TFile::Open((string("/Users/avartak/CMS/MonoX/FinalTrees/top/")          + filename).c_str());
  TFile* difile = TFile::Open((string("/Users/avartak/CMS/MonoX/FinalTrees/dibosons/")     + filename).c_str());
  TFile* qcfile = TFile::Open((string("/Users/avartak/CMS/MonoX/FinalTrees/qcd/")          + filename).c_str());
  TFile* vlfile = TFile::Open((string("/Users/avartak/CMS/MonoX/FinalTrees/") + vlfilename + filename).c_str());
  TFile* dtfile = TFile::Open((string("/Users/avartak/CMS/MonoX/FinalTrees/") + dtfilename + filename).c_str());

  TH1F dthist((string("datahist")+suffix).c_str(), "", nbins, bins);
  TH1F vlhist((string("lbkghist")+suffix).c_str(), "", nbins, bins);
  TH1F tthist((string("tbkghist")+suffix).c_str(), "", nbins, bins);
  TH1F dihist((string("dbkghist")+suffix).c_str(), "", nbins, bins);
  TH1F qchist((string("qbkghist")+suffix).c_str(), "", nbins, bins);

  TTree* dttree = (TTree*)dtfile->Get("tree");
  TTree* vltree = (TTree*)vlfile->Get("tree");
  TTree* tttree = (TTree*)ttfile->Get("tree");
  TTree* ditree = (TTree*)difile->Get("tree");
  TTree* qctree = (TTree*)qcfile->Get("tree");
  
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
  vector<TH1*> vhists;
  // apply NLO QCD and EWK corrections for Zll and Wlnu
  if (vlfilename == "wln100toinf/") {vhists.push_back(wnlohist); vhists.push_back(wewkhist);}
  if (vlfilename == "zll100toinf/") {vhists.push_back(znlohist); vhists.push_back(zewkhist);}

  makehist4(dttree, &dthist, false, sample, 1.00, ehists, NULL);
  makehist4(vltree, &vlhist,  true, sample, 1.00, ehists, NULL);
  makehist4(tttree, &tthist,  true, sample, 1.00, ehists, NULL);
  makehist4(ditree, &dihist,  true, sample, 1.00, ehists, NULL);
  makehist4(qctree, &qchist,  true, sample, 1.00, ehists, NULL);

  outfile->cd();
  dthist.Write();
  tthist.Write();
  dihist.Write();
  qchist.Write();
  vlhist.Write();
  
  dtfile->Close();
  vlfile->Close();
  ttfile->Close();
  difile->Close();
  qcfile->Close();
  
  cout << "Templates for the lepton control region computed ..." << endl;
}

// Run the final analysis:
// 1) Store all corrections templates from input files (complient to combine)
// 2) Make data and expected yields templates for all the other processes

void hists() {

  // take correction files --> central value
  TFile* zmmcorfile = TFile::Open("zmmcor.root");
  TFile* zeecorfile = TFile::Open("zeecor.root");
  TFile* wmncorfile = TFile::Open("wmncor.root");
  TFile* wencorfile = TFile::Open("wencor.root");
  TFile* zwjcorfile = TFile::Open("zwjcor.root");
  TFile* gamcorfile = TFile::Open("gamcor.root");

  // QCD, EWK, factm re and footprint on Z/gamma
  TFile* gamcorqcdfile = TFile::Open("gamcorqcd.root");
  TFile* gamcorewkfile = TFile::Open("gamcorewk.root");
  TFile* gamcorre1file = TFile::Open("gamcorre1.root");
  TFile* gamcorfa1file = TFile::Open("gamcorfa1.root");
  TFile* gamcorre2file = TFile::Open("gamcorre2.root");
  TFile* gamcorfa2file = TFile::Open("gamcorfa2.root");
  TFile* gamcorpdffile = TFile::Open("gamcorpdf.root");
  TFile* gamcorfpcfile = TFile::Open("gamcorfpc.root");
  
  // QCD, EWK, factm re and footprint on Z/W
  TFile* zwjcorqcdfile = TFile::Open("zwjcorqcd.root");
  TFile* zwjcorewkfile = TFile::Open("zwjcorewk.root");
  TFile* zwjcorre1file = TFile::Open("zwjcorre1.root");
  TFile* zwjcorfa1file = TFile::Open("zwjcorfa1.root");
  TFile* zwjcorre2file = TFile::Open("zwjcorre2.root");
  TFile* zwjcorfa2file = TFile::Open("zwjcorfa2.root");
  TFile* zwjcorpdffile = TFile::Open("zwjcorpdf.root");

  // get histograms  
  TH1* zmmcorhist = (TH1*)zmmcorfile->Get("zmmcorhist");    
  TH1* zeecorhist = (TH1*)zeecorfile->Get("zeecorhist");    
  TH1* wmncorhist = (TH1*)wmncorfile->Get("wmncorhist");    
  TH1* wencorhist = (TH1*)wencorfile->Get("wencorhist");    
  TH1* zwjcorhist = (TH1*)zwjcorfile->Get("zwjcorhist");    
  TH1* gamcorhist = (TH1*)gamcorfile->Get("gamcorhist");    
  
  // get histograms Z/gamma
  TH1* gamcorewkhist = (TH1*)gamcorewkfile->Get("gamcorewkhist");    
  TH1* gamcorqcdhist = (TH1*)gamcorqcdfile->Get("gamcorqcdhist");    
  TH1* gamcorre1hist = (TH1*)gamcorre1file->Get("gamcorre1hist");    
  TH1* gamcorfa1hist = (TH1*)gamcorfa1file->Get("gamcorfa1hist");    
  TH1* gamcorre2hist = (TH1*)gamcorre2file->Get("gamcorre2hist");    
  TH1* gamcorfa2hist = (TH1*)gamcorfa2file->Get("gamcorfa2hist");    
  TH1* gamcorpdfhist = (TH1*)gamcorpdffile->Get("gamcorpdfhist");    
  TH1* gamcorfpchist = (TH1*)gamcorfpcfile->Get("gamcorfpchist");    
  
  // uncertainty histogram for combine
  TH1* gamuncewkhist = (TH1*)gamcorewkhist->Clone("gamuncewkhist");    
  gamuncewkhist->Divide(gamcorqcdhist);
  for (int i = 1; i <= gamuncewkhist->GetNbinsX(); i++) 
    gamuncewkhist->SetBinContent(i, fabs(gamuncewkhist->GetBinContent(i)-1.0));
  gamuncewkhist->SetName("ZG_EWK");

  TH1* gamuncre1hist = (TH1*)gamcorre1hist->Clone("gamuncre1hist");    
  gamuncre1hist->Divide(gamcorqcdhist);
  for (int i = 1; i <= gamuncre1hist->GetNbinsX(); i++) 
    gamuncre1hist->SetBinContent(i, fabs(gamuncre1hist->GetBinContent(i)-1.0));
  gamuncre1hist->SetName("ZG_RenScale1");

  TH1* gamuncfa1hist = (TH1*)gamcorfa1hist->Clone("gamuncfa1hist");    
  gamuncfa1hist->Divide(gamcorqcdhist);
  for (int i = 1; i <= gamuncfa1hist->GetNbinsX(); i++) 
    gamuncfa1hist->SetBinContent(i, fabs(gamuncfa1hist->GetBinContent(i)-1.0));
  gamuncfa1hist->SetName("ZG_FactScale1");

  TH1* gamuncre2hist = (TH1*)gamcorre2hist->Clone("gamuncre2hist");    
  gamuncre2hist->Divide(gamcorqcdhist);
  for (int i = 1; i <= gamuncre2hist->GetNbinsX(); i++) 
    gamuncre2hist->SetBinContent(i, fabs(gamuncre2hist->GetBinContent(i)-1.0));
  gamuncre2hist->SetName("ZG_RenScale2");

  TH1* gamuncfa2hist = (TH1*)gamcorfa2hist->Clone("gamuncfa2hist");    
  gamuncfa2hist->Divide(gamcorqcdhist);
  for (int i = 1; i <= gamuncfa2hist->GetNbinsX(); i++) 
    gamuncfa2hist->SetBinContent(i, fabs(gamuncfa2hist->GetBinContent(i)-1.0));
  gamuncfa2hist->SetName("ZG_FactScale2");

  TH1* gamuncpdfhist = (TH1*)gamcorpdfhist->Clone("gamuncpdfhist");    
  gamuncpdfhist->Divide(gamcorqcdhist);
  for (int i = 1; i <= gamuncpdfhist->GetNbinsX(); i++) 
    gamuncpdfhist->SetBinContent(i, fabs(gamuncpdfhist->GetBinContent(i)-1.0));
  gamuncpdfhist->SetName("ZG_PDF");

  TH1* gamuncfpchist = (TH1*)gamcorfpchist->Clone("gamuncfpchist");    
  gamuncfpchist->Divide(gamcorqcdhist);
  for (int i = 1; i <= gamuncfpchist->GetNbinsX(); i++) 
    gamuncfpchist->SetBinContent(i, fabs(gamuncfpchist->GetBinContent(i)-1.0));
  gamuncfpchist->SetName("ZG_Footprint");

  // Same thing for Z/W ratio
  TH1* zwjcorewkhist = (TH1*)zwjcorewkfile->Get("zwjcorewkhist");    
  TH1* zwjcorqcdhist = (TH1*)zwjcorqcdfile->Get("zwjcorqcdhist");    
  TH1* zwjcorre1hist = (TH1*)zwjcorre1file->Get("zwjcorre1hist");    
  TH1* zwjcorfa1hist = (TH1*)zwjcorfa1file->Get("zwjcorfa1hist");    
  TH1* zwjcorre2hist = (TH1*)zwjcorre2file->Get("zwjcorre2hist");    
  TH1* zwjcorfa2hist = (TH1*)zwjcorfa2file->Get("zwjcorfa2hist");    
  TH1* zwjcorpdfhist = (TH1*)zwjcorpdffile->Get("zwjcorpdfhist");    
  
  TH1* zwjuncewkhist = (TH1*)zwjcorewkhist->Clone("zwjuncewkhist");    
  zwjuncewkhist->Divide(zwjcorqcdhist);
  for (int i = 1; i <= zwjuncewkhist->GetNbinsX(); i++) 
    zwjuncewkhist->SetBinContent(i, fabs(zwjuncewkhist->GetBinContent(i)-1.0));
  zwjuncewkhist->SetName("ZW_EWK");

  TH1* zwjuncre1hist = (TH1*)zwjcorre1hist->Clone("zwjuncre1hist");
  zwjuncre1hist->Divide(zwjcorqcdhist);
  for (int i = 1; i <= zwjuncre1hist->GetNbinsX(); i++) 
    zwjuncre1hist->SetBinContent(i, fabs(zwjuncre1hist->GetBinContent(i)-1.0));
  zwjuncre1hist->SetName("ZW_RenScale1");

  TH1* zwjuncfa1hist = (TH1*)zwjcorfa1hist->Clone("zwjuncfa1hist");
  zwjuncfa1hist->Divide(zwjcorqcdhist);
  for (int i = 1; i <= zwjuncfa1hist->GetNbinsX(); i++) 
    zwjuncfa1hist->SetBinContent(i, fabs(zwjuncfa1hist->GetBinContent(i)-1.0));
  zwjuncfa1hist->SetName("ZW_FactScale1");

  TH1* zwjuncre2hist = (TH1*)zwjcorre2hist->Clone("zwjuncre2hist");
  zwjuncre2hist->Divide(zwjcorqcdhist);
  for (int i = 1; i <= zwjuncre2hist->GetNbinsX(); i++) 
    zwjuncre2hist->SetBinContent(i, fabs(zwjuncre2hist->GetBinContent(i)-1.0));
  zwjuncre2hist->SetName("ZW_RenScale2");

  TH1* zwjuncfa2hist = (TH1*)zwjcorfa2hist->Clone("zwjuncfa2hist");
  zwjuncfa2hist->Divide(zwjcorqcdhist);
  for (int i = 1; i <= zwjuncfa2hist->GetNbinsX(); i++) 
    zwjuncfa2hist->SetBinContent(i, fabs(zwjuncfa2hist->GetBinContent(i)-1.0));
  zwjuncfa2hist->SetName("ZW_FactScale2");
  
  TH1* zwjuncpdfhist = (TH1*)zwjcorpdfhist->Clone("zwjuncpdfhist");
  zwjuncpdfhist->Divide(zwjcorqcdhist);
  for (int i = 1; i <= zwjuncpdfhist->GetNbinsX(); i++) 
    zwjuncpdfhist->SetBinContent(i, fabs(zwjuncpdfhist->GetBinContent(i)-1.0));
  zwjuncpdfhist->SetName("ZW_PDF");
    
  TFile outfile("templates.root", "RECREATE");

  zmmcorhist->Write();
  zeecorhist->Write();
  wmncorhist->Write();
  wencorhist->Write();
  zwjcorhist->Write();
  gamcorhist->Write();

  gamcorqcdhist->Write();
  gamcorewkhist->Write();
  gamcorre1hist->Write();
  gamcorfa1hist->Write();
  gamcorre2hist->Write();
  gamcorfa2hist->Write();
  gamcorpdfhist->Write();
  gamcorfpchist->Write();
  gamuncewkhist->Write();
  gamuncre1hist->Write();
  gamuncfa1hist->Write();
  gamuncre2hist->Write();
  gamuncfa2hist->Write();
  gamuncpdfhist->Write();
  gamuncfpchist->Write();

  zwjcorqcdhist->Write();
  zwjcorewkhist->Write();
  zwjcorre1hist->Write();
  zwjcorfa1hist->Write();
  zwjcorre2hist->Write();
  zwjcorfa2hist->Write();
  zwjcorpdfhist->Write();
  zwjuncewkhist->Write();
  zwjuncre1hist->Write();
  zwjuncfa1hist->Write();
  zwjuncre2hist->Write();
  zwjuncfa2hist->Write();
  zwjuncpdfhist->Write();
  
  sigdatamchist(&outfile);
  gamdatamchist(&outfile);
  lepdatamchist(&outfile, 1);
  lepdatamchist(&outfile, 2);
  lepdatamchist(&outfile, 3);
  lepdatamchist(&outfile, 4);
  
  outfile.Close();
}
