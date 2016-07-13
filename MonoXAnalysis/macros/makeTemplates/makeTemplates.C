#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <sstream>

#include "makehist.h"
#include "makeCorrHistograms.C"
#include "makeDataHistograms.C"
#include "makeSignalHistograms.C"

using namespace std;

///////////////////////////// from TFs files to the template one
void fillAndSaveCorrHistograms(const vector<string> & observables, TFile & outputFile, const string & outDir, const bool & addTop, const string & ext, 
			       const bool & addHistoForCutAndCount = false){

  // take correction files --> central value
  cout<<"Re-open file for correction histo"<<endl;
  TFile* zmmcorfile = TFile::Open((outDir+"/zmmcor"+ext+".root").c_str());
  TFile* zeecorfile = TFile::Open((outDir+"/zeecor"+ext+".root").c_str());
  TFile* wmncorfile = TFile::Open((outDir+"/wmncor"+ext+".root").c_str());
  TFile* wencorfile = TFile::Open((outDir+"/wencor"+ext+".root").c_str());
  TFile* zwjcorfile = TFile::Open((outDir+"/zwjcor"+ext+".root").c_str());
  TFile* gamcorfile = TFile::Open((outDir+"/gamcor"+ext+".root").c_str());
  TFile* wgamcorfile = TFile::Open((outDir+"/wgamcor"+ext+".root").c_str());
  TFile* topmucorfile = NULL;
  TFile* topelcorfile = NULL;
  if(addTop){
    topmucorfile = TFile::Open((outDir+"/topmucor"+ext+".root").c_str());
    topelcorfile = TFile::Open((outDir+"/topelcor"+ext+".root").c_str());
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

  // QCD, EWK, factm re and footprint on W/gamma
  TFile* wgamcorqcdfile = TFile::Open((outDir+"/wgamcorqcd"+ext+".root").c_str());
  TFile* wgamcorewkfile = TFile::Open((outDir+"/wgamcorewk"+ext+".root").c_str());
  TFile* wgamcorre1file = TFile::Open((outDir+"/wgamcorre1"+ext+".root").c_str());
  TFile* wgamcorfa1file = TFile::Open((outDir+"/wgamcorfa1"+ext+".root").c_str());
  TFile* wgamcorre2file = TFile::Open((outDir+"/wgamcorre2"+ext+".root").c_str());
  TFile* wgamcorfa2file = TFile::Open((outDir+"/wgamcorfa2"+ext+".root").c_str());
  TFile* wgamcorpdffile = TFile::Open((outDir+"/wgamcorpdf"+ext+".root").c_str());
  TFile* wgamcorfpcfile = TFile::Open((outDir+"/wgamcorfpc"+ext+".root").c_str());
  
  // QCD, EWK, factm re and footprint on Z/W
  TFile* zwjcorqcdfile = TFile::Open((outDir+"/zwjcorqcd"+ext+".root").c_str());
  TFile* zwjcorewkfile = TFile::Open((outDir+"/zwjcorewk"+ext+".root").c_str());
  TFile* zwjcorre1file = TFile::Open((outDir+"/zwjcorre1"+ext+".root").c_str());
  TFile* zwjcorfa1file = TFile::Open((outDir+"/zwjcorfa1"+ext+".root").c_str());
  TFile* zwjcorre2file = TFile::Open((outDir+"/zwjcorre2"+ext+".root").c_str());
  TFile* zwjcorfa2file = TFile::Open((outDir+"/zwjcorfa2"+ext+".root").c_str());
  TFile* zwjcorpdffile = TFile::Open((outDir+"/zwjcorpdf"+ext+".root").c_str());
  
  // Top btag up and down
  TFile* topmucorbupfile   = NULL;
  TFile* topmucorbdownfile = NULL;
  TFile* topelcorbupfile   = NULL;
  TFile* topelcorbdownfile = NULL;
  if(addTop){
    topmucorbupfile   = TFile::Open((outDir+"/topmucor"+ext+"bUp.root").c_str());
    topmucorbdownfile = TFile::Open((outDir+"/topmucor"+ext+"bDown.root").c_str());
    topelcorbupfile   = TFile::Open((outDir+"/topelcor"+ext+"bUp.root").c_str());
    topelcorbdownfile = TFile::Open((outDir+"/topelcor"+ext+"bDown.root").c_str());
  }
  // get histograms  
  vector<TH1*> zmmcorhist;
  vector<TH1*> zmmcorhist_num;
  vector<TH1*> zmmcorhist_den;

  vector<TH1*> zeecorhist;
  vector<TH1*> zeecorhist_num;
  vector<TH1*> zeecorhist_den;

  vector<TH1*> wmncorhist;
  vector<TH1*> wmncorhist_num;
  vector<TH1*> wmncorhist_den;

  vector<TH1*> wencorhist;
  vector<TH1*> wencorhist_num;
  vector<TH1*> wencorhist_den;

  vector<TH1*> zwjcorhist;
  vector<TH1*> zwjcorhist_num;
  vector<TH1*> zwjcorhist_den;

  vector<TH1*> gamcorhist;
  vector<TH1*> gamcorhist_num;
  vector<TH1*> gamcorhist_den;

  vector<TH1*> wgamcorhist;
  vector<TH1*> wgamcorhist_num;
  vector<TH1*> wgamcorhist_den;

  vector<TH1*> topmucorhist;
  vector<TH1*> topmucorhist_num;
  vector<TH1*> topmucorhist_den;

  vector<TH1*> topelcorhist;
  vector<TH1*> topelcorhist_num;
  vector<TH1*> topelcorhist_den;
  
  vector<TH1*> gamcorewkhist;
  vector<TH1*> gamcorewkhist_num;
  vector<TH1*> gamcorewkhist_den;

  vector<TH1*> gamcorqcdhist;
  vector<TH1*> gamcorqcdhist_num;
  vector<TH1*> gamcorqcdhist_den;

  vector<TH1*> gamcorre1hist;
  vector<TH1*> gamcorre1hist_num;
  vector<TH1*> gamcorre1hist_den;

  vector<TH1*> gamcorfa1hist;
  vector<TH1*> gamcorfa1hist_num;
  vector<TH1*> gamcorfa1hist_den;

  vector<TH1*> gamcorre2hist;
  vector<TH1*> gamcorre2hist_num;
  vector<TH1*> gamcorre2hist_den;

  vector<TH1*> gamcorfa2hist;
  vector<TH1*> gamcorfa2hist_num;
  vector<TH1*> gamcorfa2hist_den;

  vector<TH1*> gamcorpdfhist;
  vector<TH1*> gamcorpdfhist_num;
  vector<TH1*> gamcorpdfhist_den;

  vector<TH1*> gamcorfpchist;
  vector<TH1*> gamcorfpchist_num;
  vector<TH1*> gamcorfpchist_den;
  
  vector<TH1*> zwjcorewkhist;
  vector<TH1*> zwjcorewkhist_num;
  vector<TH1*> zwjcorewkhist_den;

  vector<TH1*> zwjcorqcdhist;
  vector<TH1*> zwjcorqcdhist_num;
  vector<TH1*> zwjcorqcdhist_den;

  vector<TH1*> zwjcorre1hist;
  vector<TH1*> zwjcorre1hist_num;
  vector<TH1*> zwjcorre1hist_den;

  vector<TH1*> zwjcorre2hist;
  vector<TH1*> zwjcorre2hist_num;
  vector<TH1*> zwjcorre2hist_den;

  vector<TH1*> zwjcorfa1hist;
  vector<TH1*> zwjcorfa1hist_num;
  vector<TH1*> zwjcorfa1hist_den;

  vector<TH1*> zwjcorfa2hist;
  vector<TH1*> zwjcorfa2hist_num;
  vector<TH1*> zwjcorfa2hist_den;

  vector<TH1*> zwjcorpdfhist;
  vector<TH1*> zwjcorpdfhist_num;
  vector<TH1*> zwjcorpdfhist_den;

  vector<TH1*> wgamcorewkhist;
  vector<TH1*> wgamcorewkhist_num;
  vector<TH1*> wgamcorewkhist_den;

  vector<TH1*> wgamcorqcdhist;
  vector<TH1*> wgamcorqcdhist_num;
  vector<TH1*> wgamcorqcdhist_den;

  vector<TH1*> wgamcorre1hist;
  vector<TH1*> wgamcorre1hist_num;
  vector<TH1*> wgamcorre1hist_den;

  vector<TH1*> wgamcorfa1hist;
  vector<TH1*> wgamcorfa1hist_num;
  vector<TH1*> wgamcorfa1hist_den;

  vector<TH1*> wgamcorre2hist;
  vector<TH1*> wgamcorre2hist_num;
  vector<TH1*> wgamcorre2hist_den;

  vector<TH1*> wgamcorfa2hist;
  vector<TH1*> wgamcorfa2hist_num;
  vector<TH1*> wgamcorfa2hist_den;

  vector<TH1*> wgamcorpdfhist;
  vector<TH1*> wgamcorpdfhist_num;
  vector<TH1*> wgamcorpdfhist_den;

  vector<TH1*> wgamcorfpchist;
  vector<TH1*> wgamcorfpchist_num;
  vector<TH1*> wgamcorfpchist_den;
  
  vector<TH1*> topmucorbuphist;
  vector<TH1*> topmucorbuphist_num;
  vector<TH1*> topmucorbuphist_den;

  vector<TH1*> topmucorbdownhist;
  vector<TH1*> topmucorbdownhist_num;
  vector<TH1*> topmucorbdownhist_den;

  vector<TH1*> topelcorbuphist;
  vector<TH1*> topelcorbuphist_num;
  vector<TH1*> topelcorbuphist_den;

  vector<TH1*> topelcorbdownhist;
  vector<TH1*> topelcorbdownhist_num;
  vector<TH1*> topelcorbdownhist_den;
  
  // output file
  for(auto obs : observables){
    
    cout<<"Get histograms for observable "<<obs<<endl;
    
    zmmcorhist.push_back( (TH1*)zmmcorfile->FindObjectAny(("zmmcor"+ext+"hist_"+obs).c_str()));    
    zmmcorhist_num.push_back( (TH1*)zmmcorfile->FindObjectAny(("nhist_zmm"+ext+"_"+obs).c_str()));    
    zmmcorhist_den.push_back( (TH1*)zmmcorfile->FindObjectAny(("dhist_zmm"+ext+"_"+obs).c_str()));    

    zeecorhist.push_back( (TH1*)zeecorfile->FindObjectAny(("zeecor"+ext+"hist_"+obs).c_str()));    
    zeecorhist_num.push_back( (TH1*)zeecorfile->FindObjectAny(("nhist_zee"+ext+"_"+obs).c_str()));    
    zeecorhist_den.push_back( (TH1*)zeecorfile->FindObjectAny(("dhist_zee"+ext+"_"+obs).c_str()));    

    wmncorhist.push_back( (TH1*)wmncorfile->FindObjectAny(("wmncor"+ext+"hist_"+obs).c_str()));    
    wmncorhist_num.push_back( (TH1*)wmncorfile->FindObjectAny(("nhist_wmn"+ext+"_"+obs).c_str()));    
    wmncorhist_den.push_back( (TH1*)wmncorfile->FindObjectAny(("dhist_wmn"+ext+"_"+obs).c_str()));    

    wencorhist.push_back( (TH1*)wencorfile->FindObjectAny(("wencor"+ext+"hist_"+obs).c_str()));    
    wencorhist_num.push_back( (TH1*)wencorfile->FindObjectAny(("nhist_wen"+ext+"_"+obs).c_str()));    
    wencorhist_den.push_back( (TH1*)wencorfile->FindObjectAny(("dhist_wen"+ext+"_"+obs).c_str()));    

    zwjcorhist.push_back( (TH1*)zwjcorfile->FindObjectAny(("zwjcor"+ext+"hist_"+obs).c_str()));    
    zwjcorhist_num.push_back( (TH1*)zwjcorfile->FindObjectAny(("nhist_zwj"+ext+"_"+obs).c_str()));    
    zwjcorhist_den.push_back( (TH1*)zwjcorfile->FindObjectAny(("dhist_zwj"+ext+"_"+obs).c_str()));    

    gamcorhist.push_back( (TH1*)gamcorfile->FindObjectAny(("gamcor"+ext+"hist_"+obs).c_str()));  
    gamcorhist_num.push_back( (TH1*)gamcorfile->FindObjectAny(("nhist_gam"+ext+"_"+obs).c_str()));    
    gamcorhist_den.push_back( (TH1*)gamcorfile->FindObjectAny(("dhist_gam"+ext+"_"+obs).c_str()));    

    wgamcorhist.push_back( (TH1*) wgamcorfile->FindObjectAny(("wgamcor"+ext+"hist_"+obs).c_str()));  
    wgamcorhist_num.push_back( (TH1*) wgamcorfile->FindObjectAny(("nhist_wgam"+ext+"_"+obs).c_str()));    
    wgamcorhist_den.push_back( (TH1*) wgamcorfile->FindObjectAny(("dhist_wgam"+ext+"_"+obs).c_str()));    

    if(addTop){
      topmucorhist.push_back( (TH1*)topmucorfile->FindObjectAny(("topmucor"+ext+"hist_"+obs).c_str()));    
      topmucorhist_num.push_back( (TH1*)topmucorfile->FindObjectAny(("nhist_topmu"+ext+"_"+obs).c_str()));    
      topmucorhist_den.push_back( (TH1*)topmucorfile->FindObjectAny(("dhist_topmu"+ext+"_"+obs).c_str()));    

      topelcorhist.push_back( (TH1*)topelcorfile->FindObjectAny(("topelcor"+ext+"hist_"+obs).c_str()));    
      topelcorhist_num.push_back( (TH1*)topelcorfile->FindObjectAny(("nhist_topel"+ext+"_"+obs).c_str()));    
      topelcorhist_den.push_back( (TH1*)topelcorfile->FindObjectAny(("dhist_topel"+ext+"_"+obs).c_str()));    
    }

    // get histograms Z/gamma
    cout<<"Make Z/gamma sys histograms"<<endl;
    gamcorewkhist.push_back( (TH1*)gamcorewkfile->FindObjectAny(("gamcor"+ext+"ewkhist_"+obs).c_str()));    
    gamcorewkhist_num.push_back( (TH1*)gamcorewkfile->FindObjectAny(("nhist_gam_ewk"+ext+"_"+obs).c_str()));    
    gamcorewkhist_den.push_back( (TH1*)gamcorewkfile->FindObjectAny(("dhist_gam_ewk"+ext+"_"+obs).c_str()));    

    gamcorqcdhist.push_back( (TH1*)gamcorqcdfile->FindObjectAny(("gamcor"+ext+"qcdhist_"+obs).c_str()));    
    gamcorqcdhist_num.push_back( (TH1*)gamcorqcdfile->FindObjectAny(("nhist_gam_qcd"+ext+"_"+obs).c_str()));    
    gamcorqcdhist_den.push_back( (TH1*)gamcorqcdfile->FindObjectAny(("dhist_gam_qcd"+ext+"_"+obs).c_str()));    

    gamcorre1hist.push_back( (TH1*)gamcorre1file->FindObjectAny(("gamcor"+ext+"re1hist_"+obs).c_str()));    
    gamcorre1hist_num.push_back( (TH1*)gamcorre1file->FindObjectAny(("nhist_gam_re1"+ext+"_"+obs).c_str()));    
    gamcorre1hist_den.push_back( (TH1*)gamcorre1file->FindObjectAny(("dhist_gam_re1"+ext+"_"+obs).c_str()));    

    gamcorfa1hist.push_back( (TH1*)gamcorfa1file->FindObjectAny(("gamcor"+ext+"fa1hist_"+obs).c_str()));    
    gamcorfa1hist_num.push_back( (TH1*)gamcorfa1file->FindObjectAny(("nhist_gam_fa1"+ext+"_"+obs).c_str()));    
    gamcorfa1hist_den.push_back( (TH1*)gamcorfa1file->FindObjectAny(("dhist_gam_fa1"+ext+"_"+obs).c_str()));    

    gamcorre2hist.push_back( (TH1*)gamcorre2file->FindObjectAny(("gamcor"+ext+"re2hist_"+obs).c_str()));    
    gamcorre2hist_num.push_back( (TH1*)gamcorre2file->FindObjectAny(("nhist_gam_re2"+ext+"_"+obs).c_str()));    
    gamcorre2hist_den.push_back( (TH1*)gamcorre2file->FindObjectAny(("dhist_gam_re2"+ext+"_"+obs).c_str()));    

    gamcorfa2hist.push_back( (TH1*)gamcorfa2file->FindObjectAny(("gamcor"+ext+"fa2hist_"+obs).c_str()));    
    gamcorfa2hist_num.push_back( (TH1*)gamcorfa2file->FindObjectAny(("nhist_gam_fa2"+ext+"_"+obs).c_str()));    
    gamcorfa2hist_den.push_back( (TH1*)gamcorfa2file->FindObjectAny(("dhist_gam_fa2"+ext+"_"+obs).c_str()));    

    gamcorpdfhist.push_back( (TH1*)gamcorpdffile->FindObjectAny(("gamcor"+ext+"pdfhist_"+obs).c_str()));    
    gamcorpdfhist_num.push_back( (TH1*)gamcorpdffile->FindObjectAny(("nhist_gam_pdf"+ext+"_"+obs).c_str()));    
    gamcorpdfhist_den.push_back( (TH1*)gamcorpdffile->FindObjectAny(("dhist_gam_pdf"+ext+"_"+obs).c_str()));    

    gamcorfpchist.push_back( (TH1*)gamcorfpcfile->FindObjectAny(("gamcor"+ext+"fpchist_"+obs).c_str()));    
    gamcorfpchist_num.push_back( (TH1*)gamcorfpcfile->FindObjectAny(("nhist_gam_fpc"+ext+"_"+obs).c_str()));    
    gamcorfpchist_den.push_back( (TH1*)gamcorfpcfile->FindObjectAny(("dhist_gam_fpc"+ext+"_"+obs).c_str()));    

    // uncertainty histogram for combine
    TH1* gamuncewkhist = (TH1*)gamcorewkhist.back()->Clone(("gamuncewk"+ext+"hist_"+obs).c_str());    
    gamuncewkhist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncewkhist->GetNbinsX(); i++) 
      gamuncewkhist->SetBinContent(i, fabs(gamuncewkhist->GetBinContent(i)-1.0));
    gamuncewkhist->SetName(("ZG_EWK_"+obs).c_str());

    ///
    TH1* gamuncre1hist = (TH1*)gamcorre1hist.back()->Clone(("gamuncre1"+ext+"hist_"+obs).c_str());    
    gamuncre1hist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncre1hist->GetNbinsX(); i++) 
      gamuncre1hist->SetBinContent(i, fabs(gamuncre1hist->GetBinContent(i)-1.0));
    gamuncre1hist->SetName(("ZG_RenScale1_"+obs).c_str());

    ////
    TH1* gamuncfa1hist = (TH1*)gamcorfa1hist.back()->Clone(("gamuncfa1"+ext+"hist_"+obs).c_str());    
    gamuncfa1hist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncfa1hist->GetNbinsX(); i++) 
      gamuncfa1hist->SetBinContent(i, fabs(gamuncfa1hist->GetBinContent(i)-1.0));
    gamuncfa1hist->SetName(("ZG_FactScale1_"+obs).c_str());
    
    ////
    TH1* gamuncre2hist = (TH1*)gamcorre2hist.back()->Clone(("gamuncre2"+ext+"hist_"+obs).c_str());    
    gamuncre2hist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncre2hist->GetNbinsX(); i++) 
      gamuncre2hist->SetBinContent(i, fabs(gamuncre2hist->GetBinContent(i)-1.0));
    gamuncre2hist->SetName(("ZG_RenScale2_"+obs).c_str());

    ///
    TH1* gamuncfa2hist = (TH1*)gamcorfa2hist.back()->Clone(("gamuncfa2"+ext+"hist_"+obs).c_str());    
    gamuncfa2hist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncfa2hist->GetNbinsX(); i++) 
      gamuncfa2hist->SetBinContent(i, fabs(gamuncfa2hist->GetBinContent(i)-1.0));
    gamuncfa2hist->SetName(("ZG_FactScale2_"+obs).c_str());

    ///
    TH1* gamuncpdfhist = (TH1*)gamcorpdfhist.back()->Clone(("gamuncpdf"+ext+"hist_"+obs).c_str());    
    gamuncpdfhist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncpdfhist->GetNbinsX(); i++) 
      gamuncpdfhist->SetBinContent(i, fabs(gamuncpdfhist->GetBinContent(i)-1.0));
    gamuncpdfhist->SetName(("ZG_PDF_"+obs).c_str());

    ///    
    TH1* gamuncfpchist = (TH1*)gamcorfpchist.back()->Clone(("gamuncfpc"+ext+"hist_"+obs).c_str());    
    gamuncfpchist->Divide(gamcorqcdhist.back());
    for (int i = 1; i <= gamuncfpchist->GetNbinsX(); i++) 
      gamuncfpchist->SetBinContent(i, fabs(gamuncfpchist->GetBinContent(i)-1.0));
    gamuncfpchist->SetName(("ZG_Footprint_"+obs).c_str());

    // Same thing for Z/W ratio
    cout<<"Make Z/W sys histograms"<<endl;
    zwjcorewkhist.push_back( (TH1*)zwjcorewkfile->FindObjectAny(("zwjcorewk"+ext+"hist_"+obs).c_str()));    
    zwjcorewkhist_num.push_back( (TH1*)zwjcorewkfile->FindObjectAny(("nhist_zwj_ewk"+ext+"_"+obs).c_str()));    
    zwjcorewkhist_den.push_back( (TH1*)zwjcorewkfile->FindObjectAny(("dhist_zwj_ewk"+ext+"_"+obs).c_str()));    

    zwjcorqcdhist.push_back( (TH1*)zwjcorqcdfile->FindObjectAny(("zwjcorqcd"+ext+"hist_"+obs).c_str()));    
    zwjcorqcdhist_num.push_back( (TH1*)zwjcorqcdfile->FindObjectAny(("nhist_zwj_qcd"+ext+"_"+obs).c_str()));    
    zwjcorqcdhist_den.push_back( (TH1*)zwjcorqcdfile->FindObjectAny(("dhist_zwj_qcd"+ext+"_"+obs).c_str()));    

    zwjcorre1hist.push_back( (TH1*)zwjcorre1file->FindObjectAny(("zwjcorre1"+ext+"hist_"+obs).c_str()));    
    zwjcorre1hist_num.push_back( (TH1*)zwjcorre1file->FindObjectAny(("nhist_zwj_re1"+ext+"_"+obs).c_str()));    
    zwjcorre1hist_den.push_back( (TH1*)zwjcorre1file->FindObjectAny(("dhist_zwj_re1"+ext+"_"+obs).c_str()));    

    zwjcorfa1hist.push_back( (TH1*)zwjcorfa1file->FindObjectAny(("zwjcorfa1"+ext+"hist_"+obs).c_str()));    
    zwjcorfa1hist_num.push_back( (TH1*)zwjcorfa1file->FindObjectAny(("nhist_zwj_fa1"+ext+"_"+obs).c_str()));    
    zwjcorfa1hist_den.push_back( (TH1*)zwjcorfa1file->FindObjectAny(("dhist_zwj_fa1"+ext+"_"+obs).c_str()));    

    zwjcorre2hist.push_back( (TH1*)zwjcorre2file->FindObjectAny(("zwjcorre2"+ext+"hist_"+obs).c_str()));    
    zwjcorre2hist_num.push_back( (TH1*)zwjcorre2file->FindObjectAny(("nhist_zwj_re2"+ext+"_"+obs).c_str()));    
    zwjcorre2hist_den.push_back( (TH1*)zwjcorre2file->FindObjectAny(("dhist_zwj_re2"+ext+"_"+obs).c_str()));    

    zwjcorfa2hist.push_back( (TH1*)zwjcorfa2file->FindObjectAny(("zwjcorfa2"+ext+"hist_"+obs).c_str()));    
    zwjcorfa2hist_num.push_back( (TH1*)zwjcorfa2file->FindObjectAny(("nhist_zwj_fa2"+ext+"_"+obs).c_str()));    
    zwjcorfa2hist_den.push_back( (TH1*)zwjcorfa2file->FindObjectAny(("dhist_zwj_fa2"+ext+"_"+obs).c_str()));    

    zwjcorpdfhist.push_back( (TH1*)zwjcorpdffile->FindObjectAny(("zwjcorpdf"+ext+"hist_"+obs).c_str()));    
    zwjcorpdfhist_num.push_back( (TH1*)zwjcorpdffile->FindObjectAny(("nhist_zwj_pdf"+ext+"_"+obs).c_str()));    
    zwjcorpdfhist_den.push_back( (TH1*)zwjcorpdffile->FindObjectAny(("dhist_zwj_pdf"+ext+"_"+obs).c_str()));    
      
    TH1* zwjuncewkhist = (TH1*)zwjcorewkhist.back()->Clone(("zwjuncewk"+ext+"hist_"+obs).c_str());    
    zwjuncewkhist->Divide(zwjcorqcdhist.back());
    for (int i = 1; i <= zwjuncewkhist->GetNbinsX(); i++) 
      zwjuncewkhist->SetBinContent(i, fabs(zwjuncewkhist->GetBinContent(i)-1.0));
    zwjuncewkhist->SetName(("ZW_EWK_"+obs).c_str());
    
    ////
    TH1* zwjuncre1hist = (TH1*)zwjcorre1hist.back()->Clone(("zwjuncre1"+ext+"hist_"+obs).c_str());
    zwjuncre1hist->Divide(zwjcorqcdhist.back());
    for (int i = 1; i <= zwjuncre1hist->GetNbinsX(); i++) 
      zwjuncre1hist->SetBinContent(i, fabs(zwjuncre1hist->GetBinContent(i)-1.0));
    zwjuncre1hist->SetName(("ZW_RenScale1_"+obs).c_str());

    ///
    TH1* zwjuncfa1hist = (TH1*)zwjcorfa1hist.back()->Clone(("zwjuncfa1"+ext+"hist_"+obs).c_str());
    zwjuncfa1hist->Divide(zwjcorqcdhist.back());
    for (int i = 1; i <= zwjuncfa1hist->GetNbinsX(); i++) 
      zwjuncfa1hist->SetBinContent(i, fabs(zwjuncfa1hist->GetBinContent(i)-1.0));
    zwjuncfa1hist->SetName(("ZW_FactScale1_"+obs).c_str());

    TH1* zwjuncre2hist = (TH1*)zwjcorre2hist.back()->Clone(("zwjuncre2"+ext+"hist_"+obs).c_str());
    zwjuncre2hist->Divide(zwjcorqcdhist.back());
    for (int i = 1; i <= zwjuncre2hist->GetNbinsX(); i++) 
      zwjuncre2hist->SetBinContent(i, fabs(zwjuncre2hist->GetBinContent(i)-1.0));
    zwjuncre2hist->SetName(("ZW_RenScale2_"+obs).c_str());
    
    TH1* zwjuncfa2hist = (TH1*)zwjcorfa2hist.back()->Clone(("zwjuncfa2"+ext+"hist_"+obs).c_str());
    zwjuncfa2hist->Divide(zwjcorqcdhist.back());
    for (int i = 1; i <= zwjuncfa2hist->GetNbinsX(); i++) 
      zwjuncfa2hist->SetBinContent(i, fabs(zwjuncfa2hist->GetBinContent(i)-1.0));
    zwjuncfa2hist->SetName(("ZW_FactScale2_"+obs).c_str());
    
    TH1* zwjuncpdfhist = (TH1*)zwjcorpdfhist.back()->Clone(("zwjuncpdf"+ext+"hist_"+obs).c_str());
    zwjuncpdfhist->Divide(zwjcorqcdhist.back());
    for (int i = 1; i <= zwjuncpdfhist->GetNbinsX(); i++) 
      zwjuncpdfhist->SetBinContent(i, fabs(zwjuncpdfhist->GetBinContent(i)-1.0));
    zwjuncpdfhist->SetName(("ZW_PDF_"+obs).c_str());

    /////
    cout<<"Make W/gamma sys histograms"<<endl;
    wgamcorewkhist.push_back( (TH1*) wgamcorewkfile->FindObjectAny(("wgamcor"+ext+"ewkhist_"+obs).c_str()));    
    wgamcorewkhist_num.push_back( (TH1*) wgamcorewkfile->FindObjectAny(("nhist_wgam_ewk"+ext+"_"+obs).c_str()));    
    wgamcorewkhist_den.push_back( (TH1*) wgamcorewkfile->FindObjectAny(("dhist_wgam_ewk"+ext+"_"+obs).c_str()));    

    wgamcorqcdhist.push_back( (TH1*) wgamcorqcdfile->FindObjectAny(("wgamcor"+ext+"qcdhist_"+obs).c_str()));    
    wgamcorqcdhist_num.push_back( (TH1*) wgamcorqcdfile->FindObjectAny(("nhist_wgam_qcd"+ext+"_"+obs).c_str()));    
    wgamcorqcdhist_den.push_back( (TH1*) wgamcorqcdfile->FindObjectAny(("dhist_wgam_qcd"+ext+"_"+obs).c_str()));    

    wgamcorre1hist.push_back( (TH1*) wgamcorre1file->FindObjectAny(("wgamcor"+ext+"re1hist_"+obs).c_str()));    
    wgamcorre1hist_num.push_back( (TH1*) wgamcorre1file->FindObjectAny(("nhist_wgam_re1"+ext+"_"+obs).c_str()));    
    wgamcorre1hist_den.push_back( (TH1*) wgamcorre1file->FindObjectAny(("dhist_wgam_re1"+ext+"_"+obs).c_str()));    

    wgamcorfa1hist.push_back( (TH1*) wgamcorfa1file->FindObjectAny(("wgamcor"+ext+"fa1hist_"+obs).c_str()));    
    wgamcorfa1hist_num.push_back( (TH1*) wgamcorfa1file->FindObjectAny(("nhist_wgam_fa1"+ext+"_"+obs).c_str()));    
    wgamcorfa1hist_den.push_back( (TH1*) wgamcorfa1file->FindObjectAny(("dhist_wgam_fa1"+ext+"_"+obs).c_str()));    

    wgamcorre2hist.push_back( (TH1*) wgamcorre2file->FindObjectAny(("wgamcor"+ext+"re2hist_"+obs).c_str()));    
    wgamcorre2hist_num.push_back( (TH1*) wgamcorre2file->FindObjectAny(("nhist_wgam_re2"+ext+"_"+obs).c_str()));    
    wgamcorre2hist_den.push_back( (TH1*) wgamcorre2file->FindObjectAny(("dhist_wgam_re2"+ext+"_"+obs).c_str()));    

    wgamcorfa2hist.push_back( (TH1*) wgamcorfa2file->FindObjectAny(("wgamcor"+ext+"fa2hist_"+obs).c_str()));    
    wgamcorfa2hist_num.push_back( (TH1*) wgamcorfa2file->FindObjectAny(("nhist_wgam_fa2"+ext+"_"+obs).c_str()));    
    wgamcorfa2hist_den.push_back( (TH1*) wgamcorfa2file->FindObjectAny(("dhist_wgam_fa2"+ext+"_"+obs).c_str()));    

    wgamcorpdfhist.push_back( (TH1*) wgamcorpdffile->FindObjectAny(("wgamcor"+ext+"pdfhist_"+obs).c_str()));    
    wgamcorpdfhist_num.push_back( (TH1*) wgamcorpdffile->FindObjectAny(("nhist_wgam_pdf"+ext+"_"+obs).c_str()));    
    wgamcorpdfhist_den.push_back( (TH1*) wgamcorpdffile->FindObjectAny(("dhist_wgam_pdf"+ext+"_"+obs).c_str()));    

    wgamcorfpchist.push_back( (TH1*) wgamcorfpcfile->FindObjectAny(("wgamcor"+ext+"fpchist_"+obs).c_str()));    
    wgamcorfpchist_num.push_back( (TH1*) wgamcorfpcfile->FindObjectAny(("nhist_wgam_fpc"+ext+"_"+obs).c_str()));    
    wgamcorfpchist_den.push_back( (TH1*) wgamcorfpcfile->FindObjectAny(("dhist_wgam_fpc"+ext+"_"+obs).c_str()));    

    // uncertainty histogram for combine
    TH1* wgamuncewkhist = (TH1*) wgamcorewkhist.back()->Clone(("wgamuncewk"+ext+"hist_"+obs).c_str());    
    wgamuncewkhist->Divide(wgamcorqcdhist.back());
    for (int i = 1; i <= wgamuncewkhist->GetNbinsX(); i++) 
      wgamuncewkhist->SetBinContent(i, 0.5*fabs(wgamuncewkhist->GetBinContent(i)-1.0));
    wgamuncewkhist->SetName(("WG_EWK_"+obs).c_str());

    ///
    TH1* wgamuncre1hist = (TH1*) wgamcorre1hist.back()->Clone(("wgamuncre1"+ext+"hist_"+obs).c_str());    
    wgamuncre1hist->Divide(wgamcorqcdhist.back());
    for (int i = 1; i <= wgamuncre1hist->GetNbinsX(); i++) 
      wgamuncre1hist->SetBinContent(i, 0.5*fabs(wgamuncre1hist->GetBinContent(i)-1.0));
    wgamuncre1hist->SetName(("WG_RenScale1_"+obs).c_str());

    ////
    TH1* wgamuncfa1hist = (TH1*) wgamcorfa1hist.back()->Clone(("wgamuncfa1"+ext+"hist_"+obs).c_str());    
    wgamuncfa1hist->Divide(wgamcorqcdhist.back());
    for (int i = 1; i <= wgamuncfa1hist->GetNbinsX(); i++) 
      wgamuncfa1hist->SetBinContent(i, 0.5*fabs(wgamuncfa1hist->GetBinContent(i)-1.0));
    wgamuncfa1hist->SetName(("WG_FactScale1_"+obs).c_str());
    
    ////
    TH1* wgamuncre2hist = (TH1*) wgamcorre2hist.back()->Clone(("wgamuncre2"+ext+"hist_"+obs).c_str());    
    wgamuncre2hist->Divide(wgamcorqcdhist.back());
    for (int i = 1; i <= wgamuncre2hist->GetNbinsX(); i++) 
      wgamuncre2hist->SetBinContent(i, 0.5*fabs(wgamuncre2hist->GetBinContent(i)-1.0));
    wgamuncre2hist->SetName(("WG_RenScale2_"+obs).c_str());

    ///
    TH1* wgamuncfa2hist = (TH1*) wgamcorfa2hist.back()->Clone(("wgamuncfa2"+ext+"hist_"+obs).c_str());    
    wgamuncfa2hist->Divide(wgamcorqcdhist.back());
    for (int i = 1; i <= wgamuncfa2hist->GetNbinsX(); i++) 
      wgamuncfa2hist->SetBinContent(i, 0.5*fabs(wgamuncfa2hist->GetBinContent(i)-1.0));
    wgamuncfa2hist->SetName(("WG_FactScale2_"+obs).c_str());

    ///
    TH1* wgamuncpdfhist = (TH1*) wgamcorpdfhist.back()->Clone(("wgamuncpdf"+ext+"hist_"+obs).c_str());    
    wgamuncpdfhist->Divide(wgamcorqcdhist.back());
    for (int i = 1; i <= wgamuncpdfhist->GetNbinsX(); i++) 
      wgamuncpdfhist->SetBinContent(i, 0.5*fabs(wgamuncpdfhist->GetBinContent(i)-1.0));
    wgamuncpdfhist->SetName(("WG_PDF_"+obs).c_str());

    ///    
    TH1* wgamuncfpchist = (TH1*) wgamcorfpchist.back()->Clone(("wgamuncfpc"+ext+"hist_"+obs).c_str());    
    wgamuncfpchist->Divide(wgamcorqcdhist.back());
    for (int i = 1; i <= wgamuncfpchist->GetNbinsX(); i++) 
      wgamuncfpchist->SetBinContent(i, 0.5*fabs(wgamuncfpchist->GetBinContent(i)-1.0));
    wgamuncfpchist->SetName(("WG_Footprint_"+obs).c_str());
    
    // make b-tagging top
    TH1* topmucoruncbhist = NULL;
    TH1* topelcoruncbhist = NULL;

    if(addTop){
      cout<<"Make top sys histograms"<<endl;
      topmucorbuphist.push_back( (TH1*) topmucorbupfile->FindObjectAny(("topmucor"+ext+"bUphist_"+obs).c_str()));    
      topmucorbuphist_num.push_back( (TH1*) topmucorbupfile->FindObjectAny(("nhist_topmu_bup"+ext+"_"+obs).c_str()));    
      topmucorbuphist_den.push_back( (TH1*) topmucorbupfile->FindObjectAny(("dhist_topmu_bup"+ext+"_"+obs).c_str()));    
      
      topmucorbdownhist.push_back( (TH1*) topmucorbdownfile->FindObjectAny(("topmucor"+ext+"bDownhist_"+obs).c_str()));    
      topmucorbdownhist_num.push_back( (TH1*) topmucorbdownfile->FindObjectAny(("nhist_topmu_bdown"+ext+"_"+obs).c_str()));    
      topmucorbdownhist_den.push_back( (TH1*) topmucorbdownfile->FindObjectAny(("dhist_topmu_bdown"+ext+"_"+obs).c_str()));    

      topelcorbuphist.push_back( (TH1*) topelcorbupfile->FindObjectAny(("topelcor"+ext+"bUphist_"+obs).c_str()));    
      topelcorbuphist_num.push_back( (TH1*) topelcorbupfile->FindObjectAny(("nhist_topel_bup"+ext+"_"+obs).c_str()));    
      topelcorbuphist_den.push_back( (TH1*) topelcorbupfile->FindObjectAny(("dhist_topel_bup"+ext+"_"+obs).c_str()));    
      
      topelcorbdownhist.push_back( (TH1*) topelcorbdownfile->FindObjectAny(("topelcor"+ext+"bDownhist_"+obs).c_str()));    
      topelcorbdownhist_num.push_back( (TH1*) topelcorbdownfile->FindObjectAny(("nhist_topel_bdown"+ext+"_"+obs).c_str()));    
      topelcorbdownhist_den.push_back( (TH1*) topelcorbdownfile->FindObjectAny(("dhist_topel_bdown"+ext+"_"+obs).c_str()));    

      // make symmetrization
      TH1* topmucorbuphist_tmp = (TH1*) topmucorbuphist.back()->Clone(("topmucorbup_tmp"+ext+"hist_"+obs).c_str());    
      topmucorbuphist_tmp->Divide(topmucorhist.back());
      TH1* topmucorbdownhist_tmp = (TH1*) topmucorbdownhist.back()->Clone(("topmucorbdown_tmp"+ext+"hist_"+obs).c_str());    
      topmucorbdownhist_tmp->Divide(topmucorhist.back());
      
      topmucoruncbhist = (TH1*) topmucorhist.back()->Clone(("topmucoruncbhist"+ext+"hist_"+obs).c_str());    
      for (int i = 1; i <= topmucoruncbhist->GetNbinsX(); i++) 
	topmucoruncbhist->SetBinContent(i, fabs(fabs(topmucorbuphist_tmp->GetBinContent(i)+topmucorbdownhist_tmp->GetBinContent(i))/2-1.0));
      
      topmucoruncbhist->SetName(("TOP_MU_B_"+obs).c_str());
      
      // make symmetrization
      TH1* topelcorbuphist_tmp = (TH1*) topelcorbuphist.back()->Clone(("topelcorbup_tmp"+ext+"hist_"+obs).c_str());    
      topelcorbuphist_tmp->Divide(topelcorhist.back());
      TH1* topelcorbdownhist_tmp = (TH1*) topelcorbdownhist.back()->Clone(("topelcorbdown_tmp"+ext+"hist_"+obs).c_str());    
      topelcorbdownhist_tmp->Divide(topelcorhist.back());
      
      topelcoruncbhist = (TH1*) topelcorhist.back()->Clone(("topelcoruncbhist"+ext+"hist_"+obs).c_str());    
      for (int i = 1; i <= topelcoruncbhist->GetNbinsX(); i++) 
	topelcoruncbhist->SetBinContent(i, fabs(fabs(topelcorbuphist_tmp->GetBinContent(i)+topelcorbdownhist_tmp->GetBinContent(i))/2-1.0));
      
      topelcoruncbhist->SetName(("TOP_EL_B_"+obs).c_str());
    }

    outputFile.cd();
    
    cout<<"Save transfer factor"<<endl;
    if(not outputFile.GetDirectory("TF_ZM"))
      outputFile.mkdir("TF_ZM");
    outputFile.cd("TF_ZM");
    zmmcorhist.back()->Write();
    if(addHistoForCutAndCount){
      zmmcorhist_num.back()->Write();
      zmmcorhist_den.back()->Write();
    }

    outputFile.cd();
    if(not outputFile.GetDirectory("TF_ZE"))
      outputFile.mkdir("TF_ZE");
    outputFile.cd("TF_ZE");
    zeecorhist.back()->Write();
    if(addHistoForCutAndCount){
      zeecorhist_num.back()->Write();
      zeecorhist_den.back()->Write();
    }

    outputFile.cd();
    if(not outputFile.GetDirectory("TF_WM"))
      outputFile.mkdir("TF_WM");
    outputFile.cd("TF_WM");
    wmncorhist.back()->Write();
    if(addHistoForCutAndCount){
      wmncorhist_num.back()->Write();
      wmncorhist_den.back()->Write();
    }

    outputFile.cd();
    if(not outputFile.GetDirectory("TF_WE"))
      outputFile.mkdir("TF_WE");
    outputFile.cd("TF_WE");
    wencorhist.back()->Write();
    if(addHistoForCutAndCount){
      wencorhist_num.back()->Write();
      wencorhist_den.back()->Write();
    }

    outputFile.cd();
    if(not outputFile.GetDirectory("TF_WZ"))
      outputFile.mkdir("TF_WZ");
    outputFile.cd("TF_WZ");
    zwjcorhist.back()->Write();
    zwjcorqcdhist.back()->Write();
    zwjcorewkhist.back()->Write();
    zwjcorre1hist.back()->Write();
    zwjcorfa1hist.back()->Write();
    zwjcorre2hist.back()->Write();
    zwjcorfa2hist.back()->Write();
    zwjcorpdfhist.back()->Write();    

    if(addHistoForCutAndCount){
      zwjcorhist_num.back()->Write();
      zwjcorhist_den.back()->Write();
      zwjcorqcdhist_num.back()->Write();
      zwjcorqcdhist_den.back()->Write();
      zwjcorewkhist_num.back()->Write();
      zwjcorewkhist_den.back()->Write();
      zwjcorre1hist_num.back()->Write();
      zwjcorre1hist_den.back()->Write();
      zwjcorfa1hist_num.back()->Write();
      zwjcorfa1hist_den.back()->Write();
      zwjcorre2hist_num.back()->Write();
      zwjcorre2hist_den.back()->Write();
      zwjcorfa2hist_num.back()->Write();
      zwjcorfa2hist_den.back()->Write();
      zwjcorpdfhist_num.back()->Write();    
      zwjcorpdfhist_den.back()->Write();    
    }

    zwjuncewkhist->Write();
    zwjuncre1hist->Write();
    zwjuncfa1hist->Write();
    zwjuncre2hist->Write();
    zwjuncfa2hist->Write();
    zwjuncpdfhist->Write();

    outputFile.cd();
    if(not outputFile.GetDirectory("TF_GJ"))
      outputFile.mkdir("TF_GJ");
    outputFile.cd("TF_GJ");
    gamcorhist.back()->Write();
    gamcorqcdhist.back()->Write();
    gamcorewkhist.back()->Write();
    gamcorre1hist.back()->Write();
    gamcorfa1hist.back()->Write();
    gamcorre2hist.back()->Write();
    gamcorfa2hist.back()->Write();
    gamcorpdfhist.back()->Write();
    gamcorfpchist.back()->Write();

    if(addHistoForCutAndCount){
      gamcorqcdhist_num.back()->Write();
      gamcorqcdhist_den.back()->Write();
      gamcorewkhist_num.back()->Write();
      gamcorewkhist_den.back()->Write();
      gamcorre1hist_num.back()->Write();
      gamcorre1hist_den.back()->Write();
      gamcorfa1hist_num.back()->Write();
      gamcorfa1hist_den.back()->Write();
      gamcorre2hist_num.back()->Write();
      gamcorre2hist_den.back()->Write();
      gamcorfa2hist_num.back()->Write();
      gamcorfa2hist_den.back()->Write();
      gamcorpdfhist_num.back()->Write();
      gamcorpdfhist_den.back()->Write();
      gamcorfpchist_num.back()->Write();
      gamcorfpchist_den.back()->Write();
    }

    gamuncewkhist->Write();
    gamuncre1hist->Write();
    gamuncfa1hist->Write();
    gamuncre2hist->Write();
    gamuncfa2hist->Write();
    gamuncpdfhist->Write();
    gamuncfpchist->Write();

    outputFile.cd();
    if(not outputFile.GetDirectory("TF_WG"))
      outputFile.mkdir("TF_WG");
    outputFile.cd("TF_WG");
    wgamcorhist.back()->Write();
    wgamcorqcdhist.back()->Write();
    wgamcorewkhist.back()->Write();
    wgamcorre1hist.back()->Write();
    wgamcorfa1hist.back()->Write();
    wgamcorre2hist.back()->Write();
    wgamcorfa2hist.back()->Write();
    wgamcorpdfhist.back()->Write();
    wgamcorfpchist.back()->Write();

    if(addHistoForCutAndCount){
      wgamcorqcdhist_num.back()->Write();
      wgamcorqcdhist_den.back()->Write();
      wgamcorewkhist_num.back()->Write();
      wgamcorewkhist_den.back()->Write();
      wgamcorre1hist_num.back()->Write();
      wgamcorre1hist_den.back()->Write();
      wgamcorfa1hist_num.back()->Write();
      wgamcorfa1hist_den.back()->Write();
      wgamcorre2hist_num.back()->Write();
      wgamcorre2hist_den.back()->Write();
      wgamcorfa2hist_num.back()->Write();
      wgamcorfa2hist_den.back()->Write();
      wgamcorpdfhist_num.back()->Write();
      wgamcorpdfhist_den.back()->Write();
      wgamcorfpchist_num.back()->Write();
      wgamcorfpchist_den.back()->Write();
    }

    wgamuncewkhist->Write();
    wgamuncre1hist->Write();
    wgamuncfa1hist->Write();
    wgamuncre2hist->Write();
    wgamuncfa2hist->Write();
    wgamuncpdfhist->Write();
    wgamuncfpchist->Write();
    
    outputFile.cd();
    if(addTop){
      if(not outputFile.GetDirectory("TF_TM"))
	outputFile.mkdir("TF_TM");
      outputFile.cd("TF_TM");
      topmucorhist.back()->Write();
      topmucoruncbhist->Write();

      if(addHistoForCutAndCount){
	topmucorhist_num.back()->Write();
	topmucorhist_den.back()->Write();
	topmucorbuphist_num.back()->Write();
	topmucorbuphist_den.back()->Write();
	topmucorbdownhist_num.back()->Write();
	topmucorbdownhist_num.back()->Write();
      }

      outputFile.cd();
      if(not outputFile.GetDirectory("TF_TE"))
	outputFile.mkdir("TF_TE");
      outputFile.cd("TF_TE");
      topelcorhist.back()->Write();          
      topelcoruncbhist->Write();      

      if(addHistoForCutAndCount){
	topelcorhist_num.back()->Write();
	topelcorhist_den.back()->Write();
	topelcorbuphist_num.back()->Write();
	topelcorbuphist_den.back()->Write();
	topelcorbdownhist_num.back()->Write();
	topelcorbdownhist_num.back()->Write();
      }
      
      outputFile.cd();
    }
  }
}


// Run the final analysis:
// 1) Store all corrections templates from input files (complient to combine)
// 2) Make data and expected yields templates for all the other processes
void makeTemplates(bool doCorrectionHistograms   = false,  // calculate transfer factors and sys
		   bool skipCorrectionHistograms = false,  // skip to open and dump transfer factors
		   Category category             = Category::monojet,  // 0 = inclusive mono-j, 1 = exclsuive mono-j, 2 V-tag HP ..
		   double lumi                   = 2.30, // 
		   string outDir                 = "", // output dir for template file
		   string templateSuffix         = "",  // suffix for the output file
		   vector<string> observables    = {"met"}, // 1D histo
		   vector<string> observables_2D = {},  // 2D histo
		   bool doShapeSystematics    = false, // run all the met, b-tag shape variations
		   bool makeResonantSelection = false, // split top in resonant and non resonant
		   int  typeOfDMSignal        = 0,     // 0 means both mono-j and mono-V, 1 is mono-j, 2 is mono-V
		   bool runHiggsInvisible     = false, // run Higgs invisible analysis
		   bool runOnlySignal         = false, // produce a file with only signal templates
		   bool runOnlyBackground     = false, // produce a file with only background templates
		   bool useNLOSamples         = false,
		   bool applyPostFitWeights   = false,
		   bool addTop                = false,
		   bool addHistoForCutAndCount= false) {

  system(("mkdir -p "+outDir).c_str());

  // to initialize the binning map
  initializeBinning();

  // find all possible mass pont to use in the analysis for each Model: Vector, Axial, Scalar and Pseudoscalar .. if onlyMonoJetSignal is true just use all the available mono-j signals in the directories
  vector<signalSample> signalMassPoint;
  if(not runHiggsInvisible and not runOnlyBackground){
    findAllPossibleMassPoints(signalMassPoint,"Vector",typeOfDMSignal);  
    findAllPossibleMassPoints(signalMassPoint,"Axial",typeOfDMSignal);
    findAllPossibleMassPoints(signalMassPoint,"Scalar",typeOfDMSignal);
    findAllPossibleMassPoints(signalMassPoint,"Pseudoscalar",typeOfDMSignal);
  }

  if(doCorrectionHistograms){    
    // make central values
    if(category == Category::VBF){

      cout<<"make correction histogram for Zmm to Znn"<<endl;
      makezmmcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",		   
		     baseInputTreePath+"/DYJets/zmmfilter/",baseInputTreePath+"/ZJetsToLLEWK/zmmfilter/",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false); 
      
      cout<<"make correction histogram for Zee to Znn"<<endl;
      makezeecorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/DYJets/zeefilter/",baseInputTreePath+"/ZJetsToLLEWK/zeefilter/",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false); 
      
      cout<<"make correction histogram for Wmn to WJets"<<endl;
      makewmncorhist(baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		     baseInputTreePath+"/WJetsNLO/wmnfilter/",baseInputTreePath+"/WJetsEWK/wmnfilter/",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 
      
      cout<<"make correction histogram for Wen to WJets"<<endl;
      makewencorhist(baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		     baseInputTreePath+"/WJetsNLO/wenfilter/",baseInputTreePath+"/WJetsEWK/wenfilter/",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 
      
      cout<<"make correction histogram for Gam+jets to Znn"<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/PhotonJets/gamfilter/", 		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false);
      
      
      cout<<"systematics on Z/gamma ratio --> NLO QCD "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"qcd",1);
      
      cout<<"systematics on Z/gamma ratio --> NLO EWK "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"ewk",2);
      
      cout<<"systematics on Z/gamma ratio --> RE 1 "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"re1",3);
      
      cout<<"systematics on Z/gamma ratio --> FA 1 "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"fa1",4);
      
      
      cout<<"systematics on Z/gamma ratio --> RE 2 "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"re2",5);
      
      cout<<"systematics on Z/gamma ratio --> FA 2 "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"fa2",6);
      
      cout<<"systematics on Z/gamma ratio --> PDF "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"pdf",7);
      
      cout<<"systematics on Z/gamma ratio --> FP "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"fpc",8);
      
      cout<<"make Z/W ratio"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 
      
      cout<<"systematics on Z/W ratio --> NLO QCD"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"qcd",1);
      
      cout<<"systematics on Z/W ratio --> NLO EWK"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"ewk",2);
      
      cout<<"systematics on Z/W ratio --> RE 1"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",		   
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"re1",3);
      
      
      cout<<"systematics on Z/W ratio --> FA 1"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"fa1",4);
      
      
      cout<<"systematics on Z/W ratio --> RE 2"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"re2",5);
      
      
      cout<<"systematics on Z/W ratio --> FA 2"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"fa2",6);
      
      cout<<"systematics on Z/W ratio --> PDF"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/",baseInputTreePath+"/ZJetsToNuNuEWK/sigfilter/",
		     baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"pdf",7);

      
      cout<<"make W/gamma ratio "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true);

      cout<<"systematics W/gamma ratio --> NLO QCD "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"qcd",1);
      
      cout<<"systematics W/gamma ratio --> NLO QCD+EWK "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"ewk",2);

      cout<<"systematics W/gamma ratio --> RE 1 "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"re1",3);
      
      cout<<"systematics W/gamma ratio --> FA 1 "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"fa1",4);
		    
      cout<<"systematics W/gamma ratio --> RE 2 "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"re2",5);
      
      cout<<"systematics W/gamma ratio --> FA 2 "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"fa2",6);
      
      cout<<"systematics W/gamma ratio --> PDF "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"pdf",7);
      
      cout<<"systematics on W/gamma ratio --> FP "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/",baseInputTreePath+"/WJetsEWK/sigfilter/",
		      baseInputTreePath+"/PhotonJets/gamfilter/",		   
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"fpc",8);

    }
    else{
      
      
      cout<<"make correction histogram for Zmm to Znn"<<endl;
      makezmmcorhist(baseInputTreePath+"/ZJets/sigfilter/","",		   
		   baseInputTreePath+"/DYJets/zmmfilter/","",
		   category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false); 
    
      cout<<"make correction histogram for Zee to Znn"<<endl;
      makezeecorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/DYJets/zeefilter/","",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false); 
      
      cout<<"make correction histogram for Wmn to WJets"<<endl;
      makewmncorhist(baseInputTreePath+"/WJetsNLO/sigfilter/","",
		     baseInputTreePath+"/WJetsNLO/wmnfilter/","",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 
      
      cout<<"make correction histogram for Wen to WJets"<<endl;
      makewencorhist(baseInputTreePath+"/WJetsNLO/sigfilter/","",
		     baseInputTreePath+"/WJetsNLO/wenfilter/","",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 
      
      cout<<"make correction histogram for Gam+jets to Znn"<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/PhotonJets/gamfilter/", 		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false);
    
      
      cout<<"systematics on Z/gamma ratio --> NLO QCD "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"qcd",1);
      
      cout<<"systematics on Z/gamma ratio --> NLO EWK "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"ewk",2);
      
      cout<<"systematics on Z/gamma ratio --> RE 1 "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"re1",3);
      
      cout<<"systematics on Z/gamma ratio --> FA 1 "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"fa1",4);
    
      
      cout<<"systematics on Z/gamma ratio --> RE 2 "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"re2",5);
      
      cout<<"systematics on Z/gamma ratio --> FA 2 "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"fa2",6);
      
      cout<<"systematics on Z/gamma ratio --> PDF "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"pdf",7);

      cout<<"systematics on Z/gamma ratio --> FP "<<endl;
      makegamcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/PhotonJets/gamfilter/",		   
		     "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,false,"fpc",8);
      
      cout<<"make Z/W ratio"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/WJetsNLO/sigfilter/","",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true); 
      
      cout<<"systematics on Z/W ratio --> NLO QCD"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/WJetsNLO/sigfilter/","",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"qcd",1);
      
      cout<<"systematics on Z/W ratio --> NLO EWK"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/WJetsNLO/sigfilter/","",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"ewk",2);
      
      cout<<"systematics on Z/W ratio --> RE 1"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/WJetsNLO/sigfilter/","",		   
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"re1",3);
      
      
      cout<<"systematics on Z/W ratio --> FA 1"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/WJetsNLO/sigfilter/","",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"fa1",4);
      
      
      cout<<"systematics on Z/W ratio --> RE 2"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/WJetsNLO/sigfilter/","",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"re2",5);
      
      
      cout<<"systematics on Z/W ratio --> FA 2"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		   baseInputTreePath+"/WJetsNLO/sigfilter/","",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"fa2",6);
      
      cout<<"systematics on Z/W ratio --> PDF"<<endl;
      makezwjcorhist(baseInputTreePath+"/ZJets/sigfilter/","",
		     baseInputTreePath+"/WJetsNLO/sigfilter/","",
		     category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"pdf",7);
      
      cout<<"make W/gamma ratio "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/","",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true);
      
      cout<<"systematics W/gamma ratio --> NLO QCD "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/","",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"qcd",1);
      
      cout<<"systematics W/gamma ratio --> NLO QCD+EWK "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/","",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"ewk",2);
      
      cout<<"systematics W/gamma ratio --> RE 1 "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/","",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"re1",3);
      
      cout<<"systematics W/gamma ratio --> FA 1 "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/","",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"fa1",4);
		    
      cout<<"systematics W/gamma ratio --> RE 2 "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/","",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"re2",5);
      
      cout<<"systematics W/gamma ratio --> FA 2 "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/","",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"fa2",6);
      
      cout<<"systematics W/gamma ratio --> PDF "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/","",
		      baseInputTreePath+"/PhotonJets/gamfilter/",
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"pdf",7);

      cout<<"systematics on W/gamma ratio --> FP "<<endl;
      makewgamcorhist(baseInputTreePath+"/WJetsNLO/sigfilter/","",
		      baseInputTreePath+"/PhotonJets/gamfilter/",		   
		      "$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/photonSF/FP_v2.root",
		      category,observables,observables_2D,lumi,outDir,"",runHiggsInvisible,true,"fpc",8);
    }
    
    if(addTop){
      cout<<"make TOP+MU ratio"<<endl;
      maketopmucorhist(baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       category,observables,observables_2D,lumi,
		       baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       outDir,"",runHiggsInvisible);
      
      cout<<"systematics on TOP+MU ratio --> bUp"<<endl;
      maketopmucorhist(baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       category,observables,observables_2D,lumi,
		       baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       outDir,"btagUp",runHiggsInvisible,"bUp");
      
      
      cout<<"systematics on TOP+MU ratio --> bDw"<<endl;
      maketopmucorhist(baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       category,observables,observables_2D,lumi,
		       baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       outDir,"btagDown",runHiggsInvisible,"bDown");
      
      cout<<"make TOP+EL ratio"<<endl;
      maketopelcorhist(baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       category,observables,observables_2D,lumi,
		       baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       outDir,"",runHiggsInvisible);
      
      
      cout<<"systematics on TOP+EL ratio --> bUp"<<endl;
      maketopelcorhist(baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       category,observables,observables_2D,lumi,
		       baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       outDir,"btagUp",runHiggsInvisible,"bUp");
      
      cout<<"systematics on TOP+EL ratio --> bDw"<<endl;
      maketopelcorhist(baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       category,observables,observables_2D,lumi,
		       baseInputTreePath+"/Top/sigfilter/",
		       baseInputTreePath+"/Top/topfilter/",
		       outDir,"btagDown",runHiggsInvisible,"bDown");
      
      
    }
  }

  TFile outfile((outDir+"/templates_"+templateSuffix+".root").c_str(), "RECREATE");
  
  if(not skipCorrectionHistograms){
    fillAndSaveCorrHistograms(observables,outfile,outDir,addTop,"",addHistoForCutAndCount);
    fillAndSaveCorrHistograms(observables_2D,outfile,outDir,addTop,"",addHistoForCutAndCount);
  }
  
  // signal region templates
  cout<<"start signal region shapes for signal"<<endl;
  if(not runHiggsInvisible and not runOnlyBackground){
    signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Vector",lumi,doShapeSystematics);
    signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Axial",lumi,doShapeSystematics);
    signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Scalar",lumi,doShapeSystematics);
    signalmchist(&outfile,category,observables,observables_2D,signalMassPoint,"Pseudoscalar",lumi,doShapeSystematics);
  }
  else if(runHiggsInvisible and not runOnlyBackground){
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"110",{5.507E+04,4.434E+03,2.194E+03,1.309E+03},1);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"125",{4.414E+04,3.782E+03,1.373E+03,7.612E+02,1.227E+02},0);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"150",{3.210E+04,3.239E+03,8.154E+02,5.279E+02},1);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"200",{1.812E+04,2.282E+03,3.023E+02,2.054E+02},1);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"300",{9.823E+03,1.256E+03,6.724E+01,4.132E+01},1);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"400",{9.516E+04,7.580E+02,2.163E+01,1.273E+01},1);
    signalHiggshist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,"500",{4.538E+03,4.872E+02,8.621E+00,5.256E+00},1);
  }

  if(not runOnlySignal){

    cout<<"start signal region data"<<endl;
    sigdatamchist(&outfile,category,observables,observables_2D,lumi,doShapeSystematics,false,false,runHiggsInvisible,false,applyPostFitWeights);
    // gamma + jets
    cout<<"start gamma+jets region data"<<endl;
    gamdatamchist(&outfile,category,observables,observables_2D,lumi,runHiggsInvisible,true,applyPostFitWeights);
    // lepton control regions
    cout<<"start zmumu region data"<<endl;
    lepdatamchist(&outfile,Sample::zmm,category,observables,observables_2D,lumi,doShapeSystematics,runHiggsInvisible,true,false,false,applyPostFitWeights); 
    cout<<"start wmunu region data"<<endl;
    lepdatamchist(&outfile,Sample::wmn,category,observables,observables_2D,lumi,doShapeSystematics,runHiggsInvisible,true,false,false,applyPostFitWeights); 
    cout<<"start zee region data"<<endl;
    lepdatamchist(&outfile,Sample::zee,category,observables,observables_2D,lumi,doShapeSystematics,runHiggsInvisible,true,false,true,applyPostFitWeights); 
    cout<<"start wenu region data"<<endl;
    lepdatamchist(&outfile,Sample::wen,category,observables,observables_2D,lumi,doShapeSystematics,runHiggsInvisible,true,false,true,applyPostFitWeights);     
    // top control regions
    if(addTop){
      cout<<"start top+mu region data"<<endl;
      topdatamchist(&outfile,Sample::topmu,category,observables,observables_2D,lumi,makeResonantSelection,doShapeSystematics,runHiggsInvisible,true,applyPostFitWeights);
      cout<<"start Top+el region data"<<endl;
      topdatamchist(&outfile,Sample::topel,category,observables,observables_2D,lumi,makeResonantSelection,doShapeSystematics,runHiggsInvisible,true,applyPostFitWeights);
    }

    //add qcd data templates
    TFile* qcdfile_data = TFile::Open("$CMSSW_BASE/src/AnalysisCode/MonoXAnalysis/data/QCD/templates_2016.root");
    if(qcdfile_data and (category == Category::monojet or category == Category::monoV)){
      cout<<"Take templates QCD from data"<<endl;
      vector<double> met_bins = selectBinning("met",category);
      TH1F*  qcd_nominal    = new TH1F("qbkghistDD_met","",int(met_bins.size()-1),&met_bins[0]);
      TH1F*  qcd_nominal_up = new TH1F("qbkghistDD_shapeUp_met","",int(met_bins.size()-1),&met_bins[0]);
      TH1F*  qcd_nominal_dw = new TH1F("qbkghistDD_shapeDw_met","",int(met_bins.size()-1),&met_bins[0]);
      
      TH1F* temp = NULL;
      if(category ==  Category::monojet)
	temp = (TH1F*) qcdfile_data->FindObjectAny("hQCD_MonoJ_nominal");
      else if(category ==  Category::monoV)
	temp = (TH1F*) qcdfile_data->FindObjectAny("hQCD_MonoV_nominal");
      
      for(int iBinX = 0; iBinX < qcd_nominal->GetNbinsX(); iBinX++)   
	qcd_nominal->SetBinContent(iBinX+1,temp->GetBinContent(temp->FindBin(qcd_nominal->GetBinCenter(iBinX+1))));
      
      if(category == Category::monojet)
	temp = (TH1F*) qcdfile_data->FindObjectAny("hQCD_MonoJ_FitUp");
      else if(category ==  Category::monoV)
	temp = (TH1F*) qcdfile_data->FindObjectAny("hQCD_MonoV_FitUp");
      
      for(int iBinX = 0; iBinX < qcd_nominal->GetNbinsX(); iBinX++)   
	qcd_nominal_up->SetBinContent(iBinX+1,temp->GetBinContent(temp->FindBin(qcd_nominal->GetBinCenter(iBinX+1))));
      
      if(category == Category::monojet)
	temp = (TH1F*) qcdfile_data->FindObjectAny("hQCD_MonoJ_FitDown");
      else if(category ==  Category::monoV)
	temp = (TH1F*) qcdfile_data->FindObjectAny("hQCD_MonoV_FitDown");
      
      for(int iBinX = 0; iBinX < qcd_nominal->GetNbinsX(); iBinX++)   
	qcd_nominal_dw->SetBinContent(iBinX+1,temp->GetBinContent(temp->FindBin(qcd_nominal->GetBinCenter(iBinX+1))));
      
      outfile.cd();
      outfile.cd("SR");
      qcd_nominal->Write();
      qcd_nominal_up->Write();
      qcd_nominal_dw->Write();      
      outfile.cd();
    }
  }  
  outfile.Close();
}

