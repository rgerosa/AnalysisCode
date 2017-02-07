///////////////////////////// from TFs files to the template one                                                                                                                                      
void fillAndSaveCorrQCDHistograms(const vector<string> & observables, // observables to be considered
				  TFile & outputFile, // output file
				  const string & outDir,  // output directory
				  const Category & category,
				  const bool & addZgamma,
				  const bool & addWgamma, 
				  const bool & addTop, 
				  const string & ext,
				  const bool & addHistoForCutAndCount = false,
				  const bool & useNewTheoryUncertainty = false){
  
  cout<<"Re-open file for correction histo"<<endl;
  TFile* zmmcorfile = TFile::Open((outDir+"/zmmcor"+ext+".root").c_str());
  TFile* zeecorfile = TFile::Open((outDir+"/zeecor"+ext+".root").c_str());
  TFile* wmncorfile = TFile::Open((outDir+"/wmncor"+ext+".root").c_str());
  TFile* wencorfile = TFile::Open((outDir+"/wencor"+ext+".root").c_str());
  TFile* zwjcorfile = TFile::Open((outDir+"/zwjcor"+ext+".root").c_str());
  TFile* gamcorfile  = NULL; TFile::Open((outDir+"/gamcor"+ext+".root").c_str());
  TFile* wgamcorfile = NULL;
  if(category != Category::VBF or (category == Category::VBF and addZgamma)){
    gamcorfile = TFile::Open((outDir+"/gamcor"+ext+".root").c_str());
    if(addWgamma)
      wgamcorfile = TFile::Open((outDir+"/wgamcor"+ext+".root").c_str());
  }
  TFile* topmucorfile = NULL;
  TFile* topelcorfile = NULL;
  if(addTop){
    topmucorfile = TFile::Open((outDir+"/topmucor"+ext+".root").c_str());
    topelcorfile = TFile::Open((outDir+"/topelcor"+ext+".root").c_str());
  }

  // QCD, EWK, factm re and footprint on Z/gamma                                                                                                                                                    
  TFile* gamcorqcdfile = NULL;
  TFile* gamcorewkfile = NULL;
  TFile* gamcorre1file = NULL;
  TFile* gamcorfa1file = NULL;
  TFile* gamcorre2file = NULL;
  TFile* gamcorfa2file = NULL;
  TFile* gamcorpdffile = NULL;
  TFile* gamcorfpcfile = NULL;

  if(category != Category::VBF or (category == Category::VBF and addZgamma)){
    gamcorqcdfile = TFile::Open((outDir+"/gamcorqcd"+ext+".root").c_str());
    gamcorewkfile = TFile::Open((outDir+"/gamcorewk"+ext+".root").c_str());
    gamcorre1file = TFile::Open((outDir+"/gamcorre1"+ext+".root").c_str());
    gamcorfa1file = TFile::Open((outDir+"/gamcorfa1"+ext+".root").c_str());
    gamcorre2file = TFile::Open((outDir+"/gamcorre2"+ext+".root").c_str());
    gamcorfa2file = TFile::Open((outDir+"/gamcorfa2"+ext+".root").c_str());
    gamcorpdffile = TFile::Open((outDir+"/gamcorpdf"+ext+".root").c_str());
    gamcorfpcfile = TFile::Open((outDir+"/gamcorfpc"+ext+".root").c_str());
  }
  // QCD, EWK, factm re and footprint on W/gamma                                                                                                                                                      
  TFile* wgamcorqcdfile = NULL;
  TFile* wgamcorewkfile = NULL;
  TFile* wgamcorre1file = NULL;
  TFile* wgamcorfa1file = NULL;
  TFile* wgamcorre2file = NULL;
  TFile* wgamcorfa2file = NULL;
  TFile* wgamcorpdffile = NULL;
  TFile* wgamcorfpcfile = NULL;

  if(addWgamma){
    wgamcorqcdfile = TFile::Open((outDir+"/wgamcorqcd"+ext+".root").c_str());
    wgamcorewkfile = TFile::Open((outDir+"/wgamcorewk"+ext+".root").c_str());
    wgamcorre1file = TFile::Open((outDir+"/wgamcorre1"+ext+".root").c_str());
    wgamcorfa1file = TFile::Open((outDir+"/wgamcorfa1"+ext+".root").c_str());
    wgamcorre2file = TFile::Open((outDir+"/wgamcorre2"+ext+".root").c_str());
    wgamcorfa2file = TFile::Open((outDir+"/wgamcorfa2"+ext+".root").c_str());
    wgamcorpdffile = TFile::Open((outDir+"/wgamcorpdf"+ext+".root").c_str());
    wgamcorfpcfile = TFile::Open((outDir+"/wgamcorfpc"+ext+".root").c_str());
  }

  // QCD, EWK, factm re and footprint on Z/W                                                                                                                                                          
  TFile* zwjcorqcdfile = TFile::Open((outDir+"/zwjcorqcd"+ext+".root").c_str());
  TFile* zwjcorewkfile = TFile::Open((outDir+"/zwjcorewk"+ext+".root").c_str());
  TFile* zwjcorre1file = NULL;
  TFile* zwjcorfa1file = NULL;
  TFile* zwjcorre2file = NULL;
  TFile* zwjcorfa2file = NULL;
  TFile* zwjcorpdffile = NULL;
  TFile* zwjcorqcdscaleupfile = NULL;
  TFile* zwjcorqcdscaledwfile = NULL;
  TFile* zwjcornloewkupfile = NULL;
  TFile* zwjcornloewkdwfile = NULL;
  TFile* zwjcorsudewkupfile = NULL;
  TFile* zwjcorsudewkdwfile = NULL;
  TFile* zwjcorewkqcdupfile = NULL;
  TFile* zwjcorewkqcddwfile = NULL;
  
  if(not useNewTheoryUncertainty){
    zwjcorre1file = TFile::Open((outDir+"/zwjcorre1"+ext+".root").c_str());
    zwjcorfa1file = TFile::Open((outDir+"/zwjcorfa1"+ext+".root").c_str());
    zwjcorre2file = TFile::Open((outDir+"/zwjcorre2"+ext+".root").c_str());
    zwjcorfa2file = TFile::Open((outDir+"/zwjcorfa2"+ext+".root").c_str());
    zwjcorpdffile = TFile::Open((outDir+"/zwjcorpdf"+ext+".root").c_str());
  }
  else{
    zwjcorqcdscaleupfile = TFile::Open((outDir+"/zwjcorqcd_scaleup"+ext+".root").c_str());
    zwjcorqcdscaledwfile = TFile::Open((outDir+"/zwjcorqcd_scaledw"+ext+".root").c_str());
    zwjcornloewkupfile = TFile::Open((outDir+"/zwjcornloewk_up"+ext+".root").c_str());
    zwjcornloewkdwfile = TFile::Open((outDir+"/zwjcornloewk_dw"+ext+".root").c_str());
    zwjcorsudewkupfile = TFile::Open((outDir+"/zwjcorsudewk_up"+ext+".root").c_str());
    zwjcorsudewkdwfile = TFile::Open((outDir+"/zwjcorsudewk_dw"+ext+".root").c_str());
    zwjcorewkqcdupfile = TFile::Open((outDir+"/zwjcorewkqcd_up"+ext+".root").c_str());
    zwjcorewkqcddwfile = TFile::Open((outDir+"/zwjcorewkqcd_dw"+ext+".root").c_str());
  }

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
  vector<TH1*> zmmcorhist, zeecorhist, wmncorhist, wencorhist, zwjcorhist, gamcorhist, wgamcorhist;
  vector<TH1*> zmmcorhist_num, zeecorhist_num, wmncorhist_num, wencorhist_num, gamcorhist_num, wgamcorhist_num, zwjcorhist_num;
  vector<TH1*> zmmcorhist_den, zeecorhist_den, wmncorhist_den, wencorhist_den, gamcorhist_den, wgamcorhist_den, zwjcorhist_den;
  vector<TH1*> topmucorhist, topelcorhist;
  vector<TH1*> topmucorhist_num, topelcorhist_num;
  vector<TH1*> topmucorhist_den, topelcorhist_den;

  vector<TH1*> gamcorewkhist, gamcorqcdhist, gamcorre1hist, gamcorfa1hist, gamcorre2hist, gamcorfa2hist, gamcorpdfhist, gamcorfpchist;
  vector<TH1*> gamcorewkhist_num, gamcorqcdhist_num, gamcorre1hist_num, gamcorfa1hist_num, gamcorre2hist_num, gamcorpdfhist_num, gamcorfa2hist_num, gamcorfpchist_num;
  vector<TH1*> gamcorewkhist_den, gamcorqcdhist_den, gamcorre1hist_den, gamcorfa1hist_den, gamcorre2hist_den, gamcorpdfhist_den, gamcorfa2hist_den, gamcorfpchist_den;

  vector<TH1*> zwjcorewkhist, zwjcorqcdhist, zwjcorre1hist, zwjcorre2hist, zwjcorfa1hist, zwjcorfa2hist, zwjcorpdfhist;
  vector<TH1*> zwjcorewkhist_num, zwjcorqcdhist_num, zwjcorre1hist_num, zwjcorre2hist_num, zwjcorfa1hist_num, zwjcorfa2hist_num, zwjcorpdfhist_num;
  vector<TH1*> zwjcorewkhist_den, zwjcorqcdhist_den, zwjcorre1hist_den, zwjcorre2hist_den, zwjcorfa1hist_den, zwjcorfa2hist_den, zwjcorpdfhist_den;

  vector<TH1*> zwjcorqcdscaleuphist, zwjcorqcdscaledwhist, zwjcornloewkuphist, zwjcornloewkdwhist, zwjcorsudewkuphist, zwjcorsudewkdwhist, zwjcorewkqcduphist, zwjcorewkqcddwhist;
  vector<TH1*> zwjcorqcdscaleuphist_num, zwjcorqcdscaledwhist_num, zwjcornloewkuphist_num, zwjcornloewkdwhist_num, zwjcorsudewkuphist_num, zwjcorsudewkdwhist_num;
  vector<TH1*> zwjcorewkqcduphist_num, zwjcorewkqcddwhist_num;
  vector<TH1*> zwjcorqcdscaleuphist_den, zwjcorqcdscaledwhist_den, zwjcornloewkuphist_den, zwjcornloewkdwhist_den, zwjcorsudewkuphist_den, zwjcorsudewkdwhist_den;
  vector<TH1*> zwjcorewkqcduphist_den, zwjcorewkqcddwhist_den;

  vector<TH1*> wgamcorewkhist, wgamcorqcdhist, wgamcorre1hist, wgamcorfa1hist, wgamcorre2hist, wgamcorfa2hist, wgamcorpdfhist, wgamcorfpchist;
  vector<TH1*> wgamcorewkhist_num, wgamcorqcdhist_num, wgamcorre1hist_num, wgamcorfa1hist_num, wgamcorre2hist_num, wgamcorfa2hist_num, wgamcorpdfhist_num, wgamcorfpchist_num;
  vector<TH1*> wgamcorewkhist_den, wgamcorqcdhist_den, wgamcorre1hist_den, wgamcorfa1hist_den, wgamcorre2hist_den, wgamcorfa2hist_den, wgamcorpdfhist_den, wgamcorfpchist_den;

  vector<TH1*> topmucorbuphist, topmucorbdownhist, topelcorbuphist, topelcorbdownhist;
  vector<TH1*> topmucorbuphist_num, topmucorbdownhist_num, topelcorbuphist_num, topelcorbdownhist_num;
  vector<TH1*> topmucorbuphist_den, topmucorbdownhist_den, topelcorbuphist_den, topelcorbdownhist_den;

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

    if(category != Category::VBF or (category == Category::VBF and addZgamma)){
      gamcorhist.push_back( (TH1*)gamcorfile->FindObjectAny(("gamcor"+ext+"hist_"+obs).c_str()));
      gamcorhist_num.push_back( (TH1*)gamcorfile->FindObjectAny(("nhist_gam"+ext+"_"+obs).c_str()));
      gamcorhist_den.push_back( (TH1*)gamcorfile->FindObjectAny(("dhist_gam"+ext+"_"+obs).c_str()));
    }

    if(addWgamma){
      wgamcorhist.push_back( (TH1*) wgamcorfile->FindObjectAny(("wgamcor"+ext+"hist_"+obs).c_str()));
      wgamcorhist_num.push_back( (TH1*) wgamcorfile->FindObjectAny(("nhist_wgam"+ext+"_"+obs).c_str()));
      wgamcorhist_den.push_back( (TH1*) wgamcorfile->FindObjectAny(("dhist_wgam"+ext+"_"+obs).c_str()));
    }

    if(addTop){
      topmucorhist.push_back( (TH1*)topmucorfile->FindObjectAny(("topmucor"+ext+"hist_"+obs).c_str()));
      topmucorhist_num.push_back( (TH1*)topmucorfile->FindObjectAny(("nhist_topmu"+ext+"_"+obs).c_str()));
      topmucorhist_den.push_back( (TH1*)topmucorfile->FindObjectAny(("dhist_topmu"+ext+"_"+obs).c_str()));

      topelcorhist.push_back( (TH1*)topelcorfile->FindObjectAny(("topelcor"+ext+"hist_"+obs).c_str()));
      topelcorhist_num.push_back( (TH1*)topelcorfile->FindObjectAny(("nhist_topel"+ext+"_"+obs).c_str()));
      topelcorhist_den.push_back( (TH1*)topelcorfile->FindObjectAny(("dhist_topel"+ext+"_"+obs).c_str()));
    }

    // get histograms Z/gamma                 
    TH1* gamuncewkhist = NULL;
    TH1* gamuncre1hist = NULL;
    TH1* gamuncfa1hist = NULL;
    TH1* gamuncre2hist = NULL;
    TH1* gamuncfa2hist = NULL;
    TH1* gamuncpdfhist = NULL;
    TH1* gamuncfpchist = NULL;

    if(category != Category::VBF or (category == Category::VBF and addZgamma)){

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
      gamuncewkhist = (TH1*)gamcorewkhist.back()->Clone(("gamuncewk"+ext+"hist_"+obs).c_str());
      gamuncewkhist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncewkhist->GetNbinsX(); i++)
	gamuncewkhist->SetBinContent(i, fabs(gamuncewkhist->GetBinContent(i)-1.0));
      gamuncewkhist->SetName(("ZG_EWK_"+obs).c_str());
      
      ///                                                                                                                                                                                             
      gamuncre1hist = (TH1*)gamcorre1hist.back()->Clone(("gamuncre1"+ext+"hist_"+obs).c_str());
      gamuncre1hist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncre1hist->GetNbinsX(); i++)
	gamuncre1hist->SetBinContent(i, fabs(gamuncre1hist->GetBinContent(i)-1.0));
      gamuncre1hist->SetName(("ZG_RenScale1_"+obs).c_str());

      ////                                                                                                                                                                                           
      gamuncfa1hist = (TH1*)gamcorfa1hist.back()->Clone(("gamuncfa1"+ext+"hist_"+obs).c_str());
      gamuncfa1hist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncfa1hist->GetNbinsX(); i++)
	gamuncfa1hist->SetBinContent(i, fabs(gamuncfa1hist->GetBinContent(i)-1.0));
      gamuncfa1hist->SetName(("ZG_FactScale1_"+obs).c_str());

      ////                                                                                                                                                                                            
      gamuncre2hist = (TH1*)gamcorre2hist.back()->Clone(("gamuncre2"+ext+"hist_"+obs).c_str());
      gamuncre2hist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncre2hist->GetNbinsX(); i++)
	gamuncre2hist->SetBinContent(i, fabs(gamuncre2hist->GetBinContent(i)-1.0));
      gamuncre2hist->SetName(("ZG_RenScale2_"+obs).c_str());
      
      ///                                                                                                                                                                                           
      gamuncfa2hist = (TH1*)gamcorfa2hist.back()->Clone(("gamuncfa2"+ext+"hist_"+obs).c_str());
      gamuncfa2hist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncfa2hist->GetNbinsX(); i++)
	gamuncfa2hist->SetBinContent(i, fabs(gamuncfa2hist->GetBinContent(i)-1.0));
      gamuncfa2hist->SetName(("ZG_FactScale2_"+obs).c_str());

      ///                                                                                                                                                                                             
      gamuncpdfhist = (TH1*)gamcorpdfhist.back()->Clone(("gamuncpdf"+ext+"hist_"+obs).c_str());
      gamuncpdfhist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncpdfhist->GetNbinsX(); i++)
	gamuncpdfhist->SetBinContent(i, fabs(gamuncpdfhist->GetBinContent(i)-1.0));
      gamuncpdfhist->SetName(("ZG_PDF_"+obs).c_str());
      
      ///                                                                                                                                                                                            
      gamuncfpchist = (TH1*)gamcorfpchist.back()->Clone(("gamuncfpc"+ext+"hist_"+obs).c_str());
      gamuncfpchist->Divide(gamcorqcdhist.back());
      for (int i = 1; i <= gamuncfpchist->GetNbinsX(); i++)
      gamuncfpchist->SetBinContent(i, fabs(gamuncfpchist->GetBinContent(i)-1.0));
      gamuncfpchist->SetName(("ZG_Footprint_"+obs).c_str());
    }

    // Same thing for Z/W ratio                                                                                                                                                                     
    cout<<"Make Z/W sys histograms"<<endl;
    zwjcorewkhist.push_back( (TH1*)zwjcorewkfile->FindObjectAny(("zwjcorewk"+ext+"hist_"+obs).c_str()));
    zwjcorewkhist_num.push_back( (TH1*)zwjcorewkfile->FindObjectAny(("nhist_zwj_ewk"+ext+"_"+obs).c_str()));
    zwjcorewkhist_den.push_back( (TH1*)zwjcorewkfile->FindObjectAny(("dhist_zwj_ewk"+ext+"_"+obs).c_str()));

    zwjcorqcdhist.push_back( (TH1*)zwjcorqcdfile->FindObjectAny(("zwjcorqcd"+ext+"hist_"+obs).c_str()));
    zwjcorqcdhist_num.push_back( (TH1*)zwjcorqcdfile->FindObjectAny(("nhist_zwj_qcd"+ext+"_"+obs).c_str()));
    zwjcorqcdhist_den.push_back( (TH1*)zwjcorqcdfile->FindObjectAny(("dhist_zwj_qcd"+ext+"_"+obs).c_str()));

    TH1* zwjuncewkhist  = NULL;
    TH1* zwjuncre1hist  = NULL;
    TH1* zwjuncfa1hist  = NULL;
    TH1* zwjuncre2hist  = NULL;
    TH1* zwjuncfa2hist  = NULL;
    TH1* zwjuncpdfhist  = NULL;
    TH1* zwjuncqcdscalehist = NULL;
    TH1* zwjuncnloewkhist = NULL;
    TH1* zwjuncsudewkhist = NULL;
    TH1* zwjuncewkqcdhist = NULL;

    if(not useNewTheoryUncertainty){
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
      
      zwjuncewkhist = (TH1*)zwjcorewkhist.back()->Clone(("zwjuncewk"+ext+"hist_"+obs).c_str());
      zwjuncewkhist->Divide(zwjcorqcdhist.back());
      for (int i = 1; i <= zwjuncewkhist->GetNbinsX(); i++)
	zwjuncewkhist->SetBinContent(i, fabs(zwjuncewkhist->GetBinContent(i)-1.0));
      zwjuncewkhist->SetName(("ZW_EWK_"+obs).c_str());

      ////                                                                                                                                                                                            
      zwjuncre1hist = (TH1*)zwjcorre1hist.back()->Clone(("zwjuncre1"+ext+"hist_"+obs).c_str());
      zwjuncre1hist->Divide(zwjcorqcdhist.back());
      for (int i = 1; i <= zwjuncre1hist->GetNbinsX(); i++)
	zwjuncre1hist->SetBinContent(i, fabs(zwjuncre1hist->GetBinContent(i)-1.0));
      zwjuncre1hist->SetName(("ZW_RenScale1_"+obs).c_str());
      
      ///                                                                                                                                                                                            
      zwjuncfa1hist = (TH1*)zwjcorfa1hist.back()->Clone(("zwjuncfa1"+ext+"hist_"+obs).c_str());
      zwjuncfa1hist->Divide(zwjcorqcdhist.back());
      for (int i = 1; i <= zwjuncfa1hist->GetNbinsX(); i++)
	zwjuncfa1hist->SetBinContent(i, fabs(zwjuncfa1hist->GetBinContent(i)-1.0));
      zwjuncfa1hist->SetName(("ZW_FactScale1_"+obs).c_str());

      zwjuncre2hist = (TH1*)zwjcorre2hist.back()->Clone(("zwjuncre2"+ext+"hist_"+obs).c_str());
      zwjuncre2hist->Divide(zwjcorqcdhist.back());
      for (int i = 1; i <= zwjuncre2hist->GetNbinsX(); i++)
	zwjuncre2hist->SetBinContent(i, fabs(zwjuncre2hist->GetBinContent(i)-1.0));
      zwjuncre2hist->SetName(("ZW_RenScale2_"+obs).c_str());
      
      zwjuncfa2hist = (TH1*)zwjcorfa2hist.back()->Clone(("zwjuncfa2"+ext+"hist_"+obs).c_str());
      zwjuncfa2hist->Divide(zwjcorqcdhist.back());
      for (int i = 1; i <= zwjuncfa2hist->GetNbinsX(); i++)
	zwjuncfa2hist->SetBinContent(i, fabs(zwjuncfa2hist->GetBinContent(i)-1.0));
      zwjuncfa2hist->SetName(("ZW_FactScale2_"+obs).c_str());

      zwjuncpdfhist = (TH1*)zwjcorpdfhist.back()->Clone(("zwjuncpdf"+ext+"hist_"+obs).c_str());
      zwjuncpdfhist->Divide(zwjcorqcdhist.back());
      for (int i = 1; i <= zwjuncpdfhist->GetNbinsX(); i++)
	zwjuncpdfhist->SetBinContent(i, fabs(zwjuncpdfhist->GetBinContent(i)-1.0));
      zwjuncpdfhist->SetName(("ZW_PDF_"+obs).c_str());      
    }
    else{ 

      // new theory uncertainties
      zwjcorqcdscaleuphist.push_back( (TH1*)zwjcorqcdscaleupfile->FindObjectAny(("zwjcorqcd_scaleup"+ext+"hist_"+obs).c_str()));
      zwjcorqcdscaleuphist_num.push_back( (TH1*)zwjcorqcdscaleupfile->FindObjectAny(("nhist_zwjqcd_scaleup"+ext+"_"+obs).c_str()));
      zwjcorqcdscaleuphist_den.push_back( (TH1*)zwjcorqcdscaleupfile->FindObjectAny(("dhist_zwjqcd_scaleup"+ext+"_"+obs).c_str()));

      zwjcorqcdscaledwhist.push_back( (TH1*)zwjcorqcdscaledwfile->FindObjectAny(("zwjcorqcd_scaledw"+ext+"hist_"+obs).c_str()));
      zwjcorqcdscaledwhist_num.push_back( (TH1*)zwjcorqcdscaledwfile->FindObjectAny(("nhist_zwjqcd_scaledw"+ext+"_"+obs).c_str()));
      zwjcorqcdscaledwhist_den.push_back( (TH1*)zwjcorqcdscaledwfile->FindObjectAny(("dhist_zwjqcd_scaledw"+ext+"_"+obs).c_str()));

      zwjcornloewkuphist.push_back( (TH1*)zwjcornloewkupfile->FindObjectAny(("zwjcornloewk_up"+ext+"hist_"+obs).c_str()));
      zwjcornloewkuphist_num.push_back( (TH1*)zwjcornloewkupfile->FindObjectAny(("nhist_zwjnloewk_up"+ext+"_"+obs).c_str()));
      zwjcornloewkuphist_den.push_back( (TH1*)zwjcornloewkupfile->FindObjectAny(("dhist_zwjnloewk_up"+ext+"_"+obs).c_str()));

      zwjcornloewkdwhist.push_back( (TH1*)zwjcornloewkdwfile->FindObjectAny(("zwjcornloewk_dw"+ext+"hist_"+obs).c_str()));
      zwjcornloewkdwhist_num.push_back( (TH1*)zwjcornloewkdwfile->FindObjectAny(("nhist_zwjnloewk_dw"+ext+"_"+obs).c_str()));
      zwjcornloewkdwhist_den.push_back( (TH1*)zwjcornloewkdwfile->FindObjectAny(("dhist_zwjnloewk_dw"+ext+"_"+obs).c_str()));

      zwjcorsudewkuphist.push_back( (TH1*)zwjcorsudewkupfile->FindObjectAny(("zwjcorsudewk_up"+ext+"hist_"+obs).c_str()));
      zwjcorsudewkuphist_num.push_back( (TH1*)zwjcorsudewkupfile->FindObjectAny(("nhist_zwjsudewk_up"+ext+"_"+obs).c_str()));
      zwjcorsudewkuphist_den.push_back( (TH1*)zwjcorsudewkupfile->FindObjectAny(("dhist_zwjsudewk_up"+ext+"_"+obs).c_str()));

      zwjcorsudewkdwhist.push_back( (TH1*)zwjcorsudewkdwfile->FindObjectAny(("zwjcorsudewk_dw"+ext+"hist_"+obs).c_str()));
      zwjcorsudewkdwhist_num.push_back( (TH1*)zwjcorsudewkdwfile->FindObjectAny(("nhist_zwjsudewk_dw"+ext+"_"+obs).c_str()));
      zwjcorsudewkdwhist_den.push_back( (TH1*)zwjcorsudewkdwfile->FindObjectAny(("dhist_zwjsudewk_dw"+ext+"_"+obs).c_str()));

      zwjcorewkqcduphist.push_back( (TH1*)zwjcorewkqcdupfile->FindObjectAny(("zwjcorewkqcd_up"+ext+"hist_"+obs).c_str()));
      zwjcorewkqcduphist_num.push_back( (TH1*)zwjcorewkqcdupfile->FindObjectAny(("nhist_zwjewkqcd_up"+ext+"_"+obs).c_str()));
      zwjcorewkqcduphist_den.push_back( (TH1*)zwjcorewkqcdupfile->FindObjectAny(("dhist_zwjewkqcd_up"+ext+"_"+obs).c_str()));

      zwjcorewkqcddwhist.push_back( (TH1*)zwjcorewkqcddwfile->FindObjectAny(("zwjcorewkqcd_dw"+ext+"hist_"+obs).c_str()));
      zwjcorewkqcddwhist_num.push_back( (TH1*)zwjcorewkqcddwfile->FindObjectAny(("nhist_zwjewkqcd_dw"+ext+"_"+obs).c_str()));
      zwjcorewkqcddwhist_den.push_back( (TH1*)zwjcorewkqcddwfile->FindObjectAny(("dhist_zwjewkqcd_dw"+ext+"_"+obs).c_str()));

      // QCD scale --> to be symmetrized
      zwjuncqcdscalehist = (TH1*) zwjcorqcdscaleuphist.back()->Clone(("zwjuncqcdscale"+ext+"hist_"+obs).c_str());
      zwjuncqcdscalehist->Reset("ICES");
      for (int i = 0; i <= zwjuncqcdscalehist->GetNbinsX()+1; i++)
        zwjuncqcdscalehist->SetBinContent(i,fabs(zwjcorqcdscaleuphist.back()->GetBinContent(i)-zwjcorqcdscaledwhist.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
      zwjuncqcdscalehist->SetName(("ZW_QCDScale_"+obs).c_str());
      //temp fix
      zwjuncqcdscalehist->SetBinContent(1,zwjuncqcdscalehist->GetBinContent(2)+0.003);

      // NLO EWK --> to be symmetrized
      zwjuncnloewkhist = (TH1*) zwjcornloewkuphist.back()->Clone(("zwjuncnloewk"+ext+"hist_"+obs).c_str());
      zwjuncnloewkhist->Reset("ICES");
      for (int i = 0; i <= zwjuncnloewkhist->GetNbinsX()+1; i++)
        zwjuncnloewkhist->SetBinContent(i,fabs(zwjcornloewkuphist.back()->GetBinContent(i)-zwjcornloewkdwhist.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
      zwjuncnloewkhist->SetName(("ZW_NLOEWK_"+obs).c_str());
      
      // NLO EWK Sud --> to be symmetrized
      zwjuncsudewkhist = (TH1*) zwjcorsudewkuphist.back()->Clone(("zwjuncsudewk"+ext+"hist_"+obs).c_str());
      zwjuncsudewkhist->Reset("ICES");
      for (int i = 0; i <= zwjuncsudewkhist->GetNbinsX()+1; i++)
        zwjuncsudewkhist->SetBinContent(i,fabs(zwjcorsudewkuphist.back()->GetBinContent(i)-zwjcorsudewkdwhist.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
      zwjuncsudewkhist->SetName(("ZW_EWKSudakov_"+obs).c_str());

      // NLO QCD+EWK --> to be symmetrized
      zwjuncewkqcdhist = (TH1*) zwjcorewkqcduphist.back()->Clone(("zwjuncewkqcd"+ext+"hist_"+obs).c_str());
      zwjuncewkqcdhist->Reset("ICES");
      for (int i = 0; i <= zwjuncewkqcdhist->GetNbinsX()+1; i++)
        zwjuncewkqcdhist->SetBinContent(i,fabs(zwjcorewkqcduphist.back()->GetBinContent(i)-zwjcorewkqcddwhist.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
      zwjuncewkqcdhist->SetName(("ZW_QCDEWK_"+obs).c_str());


    }


    /////                                                                                                                                                                                             
    TH1* wgamuncewkhist = NULL;
    TH1* wgamuncre1hist = NULL;
    TH1* wgamuncfa1hist = NULL;
    TH1* wgamuncre2hist = NULL;
    TH1* wgamuncfa2hist = NULL;
    TH1* wgamuncpdfhist = NULL;
    TH1* wgamuncfpchist = NULL;


    if(addWgamma and category != Category::VBF){

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
      wgamuncewkhist = (TH1*) wgamcorewkhist.back()->Clone(("wgamuncewk"+ext+"hist_"+obs).c_str());
      wgamuncewkhist->Divide(wgamcorqcdhist.back());
      for (int i = 1; i <= wgamuncewkhist->GetNbinsX(); i++)
        wgamuncewkhist->SetBinContent(i, 0.5*fabs(wgamuncewkhist->GetBinContent(i)-1.0));
      wgamuncewkhist->SetName(("WG_EWK_"+obs).c_str());

      ///                                                                                                                                                                                             
      wgamuncre1hist = (TH1*) wgamcorre1hist.back()->Clone(("wgamuncre1"+ext+"hist_"+obs).c_str());
      wgamuncre1hist->Divide(wgamcorqcdhist.back());
      for (int i = 1; i <= wgamuncre1hist->GetNbinsX(); i++)
        wgamuncre1hist->SetBinContent(i, 0.5*fabs(wgamuncre1hist->GetBinContent(i)-1.0));
      wgamuncre1hist->SetName(("WG_RenScale1_"+obs).c_str());

      ////                                                                                                                                                                                            
      wgamuncfa1hist = (TH1*) wgamcorfa1hist.back()->Clone(("wgamuncfa1"+ext+"hist_"+obs).c_str());
      wgamuncfa1hist->Divide(wgamcorqcdhist.back());
      for (int i = 1; i <= wgamuncfa1hist->GetNbinsX(); i++)
	wgamuncfa1hist->SetBinContent(i, 0.5*fabs(wgamuncfa1hist->GetBinContent(i)-1.0));
      wgamuncfa1hist->SetName(("WG_FactScale1_"+obs).c_str());

      ////                                                                                                                                                                                            
      wgamuncre2hist = (TH1*) wgamcorre2hist.back()->Clone(("wgamuncre2"+ext+"hist_"+obs).c_str());
      wgamuncre2hist->Divide(wgamcorqcdhist.back());
      for (int i = 1; i <= wgamuncre2hist->GetNbinsX(); i++)
        wgamuncre2hist->SetBinContent(i, 0.5*fabs(wgamuncre2hist->GetBinContent(i)-1.0));
      wgamuncre2hist->SetName(("WG_RenScale2_"+obs).c_str());

      ///                                                                                                                                                                                             
      wgamuncfa2hist = (TH1*) wgamcorfa2hist.back()->Clone(("wgamuncfa2"+ext+"hist_"+obs).c_str());
      wgamuncfa2hist->Divide(wgamcorqcdhist.back());
      for (int i = 1; i <= wgamuncfa2hist->GetNbinsX(); i++)
        wgamuncfa2hist->SetBinContent(i, 0.5*fabs(wgamuncfa2hist->GetBinContent(i)-1.0));
      wgamuncfa2hist->SetName(("WG_FactScale2_"+obs).c_str());

      ///                                                                                                                                                                                             
      wgamuncpdfhist = (TH1*) wgamcorpdfhist.back()->Clone(("wgamuncpdf"+ext+"hist_"+obs).c_str());
      wgamuncpdfhist->Divide(wgamcorqcdhist.back());
      for (int i = 1; i <= wgamuncpdfhist->GetNbinsX(); i++)
        wgamuncpdfhist->SetBinContent(i, 0.5*fabs(wgamuncpdfhist->GetBinContent(i)-1.0));
      wgamuncpdfhist->SetName(("WG_PDF_"+obs).c_str());

      ///                                                                                                                                                                                             
      wgamuncfpchist = (TH1*) wgamcorfpchist.back()->Clone(("wgamuncfpc"+ext+"hist_"+obs).c_str());
      wgamuncfpchist->Divide(wgamcorqcdhist.back());
      for (int i = 1; i <= wgamuncfpchist->GetNbinsX(); i++)
        wgamuncfpchist->SetBinContent(i, 0.5*fabs(wgamuncfpchist->GetBinContent(i)-1.0));
      wgamuncfpchist->SetName(("WG_Footprint_"+obs).c_str());
    }

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

    if(not useNewTheoryUncertainty){
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
    }
    else{

      zwjcorqcdscaleuphist.back()->Write();
      zwjcorqcdscaledwhist.back()->Write();
      zwjcornloewkuphist.back()->Write();
      zwjcornloewkdwhist.back()->Write();
      zwjcorsudewkuphist.back()->Write();
      zwjcorsudewkdwhist.back()->Write();
      zwjcorewkqcduphist.back()->Write();
      zwjcorewkqcddwhist.back()->Write();

      if(addHistoForCutAndCount){
	zwjcorqcdscaleuphist_num.back()->Write();
	zwjcorqcdscaledwhist_num.back()->Write();
	zwjcornloewkuphist_num.back()->Write();
	zwjcornloewkdwhist_num.back()->Write();
	zwjcorsudewkuphist_num.back()->Write();
	zwjcorsudewkdwhist_num.back()->Write();
	zwjcorewkqcduphist_num.back()->Write();
	zwjcorewkqcddwhist_num.back()->Write();

	zwjcorqcdscaleuphist_den.back()->Write();
	zwjcorqcdscaledwhist_den.back()->Write();
	zwjcornloewkuphist_den.back()->Write();
	zwjcornloewkdwhist_den.back()->Write();
	zwjcorsudewkuphist_den.back()->Write();
	zwjcorsudewkdwhist_den.back()->Write();
	zwjcorewkqcduphist_den.back()->Write();
	zwjcorewkqcddwhist_den.back()->Write();
      }
      
      zwjuncqcdscalehist->Write();
      zwjuncnloewkhist->Write();
      zwjuncsudewkhist->Write();
      zwjuncewkqcdhist->Write();      
    }
  
    outputFile.cd();
    if(category != Category::VBF or (category == Category::VBF and addZgamma)){
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
    }

    outputFile.cd();
    if(addWgamma and category != Category::VBF){
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
    }

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

  // clear vectors
  zmmcorhist.clear(); zeecorhist.clear(); wmncorhist.clear(); wencorhist.clear(); zwjcorhist.clear(); gamcorhist.clear(); wgamcorhist.clear();
  topmucorhist.clear(); topelcorhist.clear();
  gamcorewkhist.clear(); gamcorqcdhist.clear(); gamcorre1hist.clear(); gamcorfa1hist.clear(); gamcorre2hist.clear(); gamcorfa2hist.clear(); gamcorpdfhist.clear(); gamcorfpchist.clear();
  zwjcorewkhist.clear(); zwjcorqcdhist.clear(); zwjcorre1hist.clear(); zwjcorre2hist.clear(); zwjcorfa1hist.clear(); zwjcorfa2hist.clear(); zwjcorpdfhist.clear();
  zwjcorqcdscaleuphist.clear(); zwjcorqcdscaledwhist.clear(); zwjcornloewkuphist.clear(); zwjcornloewkdwhist.clear(); zwjcorsudewkuphist.clear(); 
  zwjcorsudewkdwhist.clear(); zwjcorewkqcduphist.clear(); zwjcorewkqcddwhist.clear();
  wgamcorewkhist.clear(); wgamcorqcdhist.clear(); wgamcorre1hist.clear(); wgamcorfa1hist.clear(); wgamcorre2hist.clear(); wgamcorfa2hist.clear(); wgamcorpdfhist.clear(); wgamcorfpchist.clear();
  topmucorbuphist.clear(); topmucorbdownhist.clear(); topelcorbuphist.clear(); topelcorbdownhist.clear();

  zmmcorhist_num.clear(); zeecorhist_num.clear(); wmncorhist_num.clear(); wencorhist_num.clear(); gamcorhist_num.clear(); wgamcorhist_num.clear(); zwjcorhist_num.clear();
  zmmcorhist_den.clear(); zeecorhist_den.clear(); wmncorhist_den.clear(); wencorhist_den.clear(); gamcorhist_den.clear(); wgamcorhist_den.clear(); zwjcorhist_den.clear();
  topmucorhist_num.clear(); topelcorhist_num.clear();
  topmucorhist_den.clear(); topelcorhist_den.clear();
  
  gamcorewkhist_num.clear(); gamcorqcdhist_num.clear(); gamcorre1hist_num.clear(); gamcorfa1hist_num.clear(); gamcorre2hist_num.clear(); gamcorpdfhist_num.clear(); 
  gamcorfa2hist_num.clear(); gamcorfpchist_num.clear();
  gamcorewkhist_den.clear(); gamcorqcdhist_den.clear(); gamcorre1hist_den.clear(); gamcorfa1hist_den.clear(); gamcorre2hist_den.clear(); gamcorpdfhist_den.clear(); 
  gamcorfa2hist_den.clear(); gamcorfpchist_den.clear();
  
  zwjcorewkhist_num.clear(); zwjcorqcdhist_num.clear(); zwjcorre1hist_num.clear(); zwjcorre2hist_num.clear(); zwjcorfa1hist_num.clear(); zwjcorfa2hist_num.clear(); zwjcorpdfhist_num.clear();
  zwjcorewkhist_den.clear(); zwjcorqcdhist_den.clear(); zwjcorre1hist_den.clear(); zwjcorre2hist_den.clear(); zwjcorfa1hist_den.clear(); zwjcorfa2hist_den.clear(); zwjcorpdfhist_den.clear();
  
  zwjcorqcdscaleuphist_num.clear(); zwjcorqcdscaledwhist_num.clear(); zwjcornloewkuphist_num.clear(); zwjcornloewkdwhist_num.clear(); zwjcorsudewkuphist_num.clear(); 
  zwjcorsudewkdwhist_num.clear(); zwjcorewkqcduphist_num.clear(); zwjcorewkqcddwhist_num.clear();
  zwjcorqcdscaleuphist_den.clear(); zwjcorqcdscaledwhist_den.clear(); zwjcornloewkuphist_den.clear(); zwjcornloewkdwhist_den.clear(); zwjcorsudewkuphist_den.clear(); 
  zwjcorsudewkdwhist_den.clear(); zwjcorewkqcduphist_den.clear(); zwjcorewkqcddwhist_den.clear();
  
  wgamcorewkhist_num.clear(); wgamcorqcdhist_num.clear(); wgamcorre1hist_num.clear(); wgamcorfa1hist_num.clear(); wgamcorre2hist_num.clear(); wgamcorfa2hist_num.clear(); 
  wgamcorpdfhist_num.clear(); wgamcorfpchist_num.clear();
  wgamcorewkhist_den.clear(); wgamcorqcdhist_den.clear(); wgamcorre1hist_den.clear(); wgamcorfa1hist_den.clear(); wgamcorre2hist_den.clear(); wgamcorfa2hist_den.clear();
  wgamcorpdfhist_den.clear(); wgamcorfpchist_den.clear();

  topmucorbuphist_num.clear(); topmucorbdownhist_num.clear(); topelcorbuphist_num.clear(); topelcorbdownhist_num.clear();
  topmucorbuphist_den.clear(); topmucorbdownhist_den.clear(); topelcorbuphist_den.clear(); topelcorbdownhist_den.clear();
  
  /// closing files
  if(zmmcorfile) zmmcorfile->Close();
  if(zeecorfile) zeecorfile->Close();
  if(wmncorfile) wmncorfile->Close();
  if(wencorfile) wencorfile->Close();
  if(gamcorfile) gamcorfile->Close();
  if(wgamcorfile) wgamcorfile->Close();
  if(topmucorfile) topmucorfile->Close();
  if(topelcorfile) topelcorfile->Close();
  if(gamcorqcdfile) gamcorqcdfile->Close();
  if(gamcorewkfile) gamcorewkfile->Close();
  if(gamcorre1file) gamcorre1file->Close();
  if(gamcorre2file) gamcorre2file->Close();
  if(gamcorfa2file) gamcorfa2file->Close();
  if(gamcorpdffile) gamcorpdffile->Close();
  if(gamcorfpcfile) gamcorfpcfile->Close();
  if(wgamcorqcdfile) wgamcorqcdfile->Close();
  if(wgamcorewkfile) wgamcorewkfile->Close();
  if(wgamcorre1file) wgamcorre1file->Close();
  if(wgamcorfa1file) wgamcorfa1file->Close();
  if(wgamcorre2file) wgamcorre2file->Close();
  if(wgamcorfa2file) wgamcorfa2file->Close();
  if(wgamcorpdffile) wgamcorpdffile->Close();
  if(wgamcorfpcfile) wgamcorfpcfile->Close();
  if(zwjcorqcdfile) zwjcorqcdfile->Close();
  if(zwjcorewkfile) zwjcorewkfile->Close();
  if(zwjcorre1file) zwjcorre1file->Close();
  if(zwjcorfa1file) zwjcorfa1file->Close();
  if(zwjcorre2file) zwjcorre2file->Close();
  if(zwjcorfa2file) zwjcorfa2file->Close();
  if(zwjcorpdffile) zwjcorpdffile->Close();
  if(topmucorbupfile) topmucorbupfile->Close();
  if(topmucorbdownfile) topmucorbdownfile->Close();
  if(topelcorbupfile) topelcorbupfile->Close();
  if(topelcorbdownfile) topelcorbdownfile->Close();
  if(zwjcorqcdscaleupfile) zwjcorqcdscaleupfile->Close();
  if(zwjcorqcdscaledwfile) zwjcorqcdscaledwfile->Close();
  if(zwjcornloewkupfile) zwjcornloewkupfile->Close();
  if(zwjcornloewkdwfile) zwjcornloewkdwfile->Close();
  if(zwjcorsudewkupfile) zwjcorsudewkupfile->Close();
  if(zwjcorsudewkdwfile) zwjcorsudewkdwfile->Close();
  if(zwjcorewkqcdupfile) zwjcorewkqcdupfile->Close();
  if(zwjcorewkqcddwfile) zwjcorewkqcddwfile->Close();
}


///////////////////////////// from TFs files to the template one                                                                                                                                      
void fillAndSaveCorrEWKHistograms(const vector<string> & observables, // observables to be considered
				  TFile & outputFile, // output file
				  const string & outDir,  // output directory
				  const Category & category,
				  const bool & addWgamma, 
				  const bool & addTop, 
				  const string & ext,
				  const bool & addHistoForCutAndCount = false){

  cout<<"Re-open file for correction histo"<<endl;
  TFile* zmmcorfile = TFile::Open((outDir+"/zewkmmcor"+ext+".root").c_str());
  TFile* zeecorfile = TFile::Open((outDir+"/zewkeecor"+ext+".root").c_str());
  TFile* wmncorfile = TFile::Open((outDir+"/wewkmncor"+ext+".root").c_str());
  TFile* wencorfile = TFile::Open((outDir+"/wewkencor"+ext+".root").c_str());
  TFile* zwjcorfile = TFile::Open((outDir+"/zwjewkcor"+ext+".root").c_str());
  TFile* gamcorfile = NULL;
  TFile* wgamcorfile = NULL;

  if(addWgamma){
    gamcorfile = TFile::Open((outDir+"/gamewkcor"+ext+".root").c_str());
    wgamcorfile = TFile::Open((outDir+"/wgamewkcor"+ext+".root").c_str());
  }

  // get histograms                                                                                                                                                                                   
  vector<TH1*> zmmcorhist, zeecorhist, wmncorhist, wencorhist, zwjcorhist, gamcorhist, wgamcorhist;
  vector<TH1*> zmmcorhist_num, zeecorhist_num, wmncorhist_num, wencorhist_num, zwjcorhist_num, gamcorhist_num, wgamcorhist_num;
  vector<TH1*> zmmcorhist_den, zeecorhist_den, wmncorhist_den, wencorhist_den, zwjcorhist_den, gamcorhist_den, wgamcorhist_den;

  // output file                                                                                                                                                                                      
  for(auto obs : observables){

    cout<<"Get histograms for observable "<<obs<<endl;
    zmmcorhist.push_back( (TH1*)zmmcorfile->FindObjectAny(("zewkmmcor"+ext+"hist_"+obs).c_str()));
    zmmcorhist_num.push_back( (TH1*)zmmcorfile->FindObjectAny(("nhist_ewk_zmm"+ext+"_"+obs).c_str()));
    zmmcorhist_den.push_back( (TH1*)zmmcorfile->FindObjectAny(("dhist_ewk_zmm"+ext+"_"+obs).c_str()));

    zeecorhist.push_back( (TH1*)zeecorfile->FindObjectAny(("zewkeecor"+ext+"hist_"+obs).c_str()));
    zeecorhist_num.push_back( (TH1*)zeecorfile->FindObjectAny(("nhist_ewk_zee"+ext+"_"+obs).c_str()));
    zeecorhist_den.push_back( (TH1*)zeecorfile->FindObjectAny(("dhist_ewk_zee"+ext+"_"+obs).c_str()));

    wmncorhist.push_back( (TH1*)wmncorfile->FindObjectAny(("wewkmncor"+ext+"hist_"+obs).c_str()));
    wmncorhist_num.push_back( (TH1*)wmncorfile->FindObjectAny(("nhist_ewk_wmn"+ext+"_"+obs).c_str()));
    wmncorhist_den.push_back( (TH1*)wmncorfile->FindObjectAny(("dhist_ewk_wmn"+ext+"_"+obs).c_str()));

    wencorhist.push_back( (TH1*)wencorfile->FindObjectAny(("wewkencor"+ext+"hist_"+obs).c_str()));
    wencorhist_num.push_back( (TH1*)wencorfile->FindObjectAny(("nhist_ewk_wen"+ext+"_"+obs).c_str()));
    wencorhist_den.push_back( (TH1*)wencorfile->FindObjectAny(("dhist_ewk_wen"+ext+"_"+obs).c_str()));

    zwjcorhist.push_back( (TH1*)zwjcorfile->FindObjectAny(("zwjewkcor"+ext+"hist_"+obs).c_str()));
    zwjcorhist_num.push_back( (TH1*)zwjcorfile->FindObjectAny(("nhist_ewk_zwj"+ext+"_"+obs).c_str()));
    zwjcorhist_den.push_back( (TH1*)zwjcorfile->FindObjectAny(("dhist_ewk_zwj"+ext+"_"+obs).c_str()));

    if(addWgamma){

      gamcorhist.push_back( (TH1*)gamcorfile->FindObjectAny(("gamewkcor"+ext+"hist_"+obs).c_str()));
      gamcorhist_num.push_back( (TH1*)gamcorfile->FindObjectAny(("nhist_ewk_gam"+ext+"_"+obs).c_str()));
      gamcorhist_den.push_back( (TH1*)gamcorfile->FindObjectAny(("dhist_ewk_gam"+ext+"_"+obs).c_str()));
      
      wgamcorhist.push_back( (TH1*) wgamcorfile->FindObjectAny(("wgamewkcor"+ext+"hist_"+obs).c_str()));
      wgamcorhist_num.push_back( (TH1*) wgamcorfile->FindObjectAny(("nhist_ewk_wgam"+ext+"_"+obs).c_str()));
      wgamcorhist_den.push_back( (TH1*) wgamcorfile->FindObjectAny(("dhist_ewk_wgam"+ext+"_"+obs).c_str()));
    }
    
    
    outputFile.cd();
    cout<<"Save transfer factor"<<endl;
    if(not outputFile.GetDirectory("TF_ZM_EWK"))
      outputFile.mkdir("TF_ZM_EWK");
    outputFile.cd("TF_ZM_EWK");
    zmmcorhist.back()->Write();
    if(addHistoForCutAndCount){
      zmmcorhist_num.back()->Write();
      zmmcorhist_den.back()->Write();
    }
    cout<<"Zee"<<endl;
    outputFile.cd();
    if(not outputFile.GetDirectory("TF_ZE_EWK"))
      outputFile.mkdir("TF_ZE_EWK");
    outputFile.cd("TF_ZE_EWK");
    zeecorhist.back()->Write();
    if(addHistoForCutAndCount){
      zeecorhist_num.back()->Write();
      zeecorhist_den.back()->Write();
    }
    cout<<"Wmn"<<endl;

    outputFile.cd();
    if(not outputFile.GetDirectory("TF_WM_EWK"))
      outputFile.mkdir("TF_WM_EWK");
    outputFile.cd("TF_WM_EWK");
    wmncorhist.back()->Write();
    if(addHistoForCutAndCount){
      wmncorhist_num.back()->Write();
      wmncorhist_den.back()->Write();
    }

    cout<<"Wen"<<endl;
    outputFile.cd();
    if(not outputFile.GetDirectory("TF_WE_EWK"))
      outputFile.mkdir("TF_WE_EWK");
    outputFile.cd("TF_WE_EWK");
    wencorhist.back()->Write();
    if(addHistoForCutAndCount){
      wencorhist_num.back()->Write();
      wencorhist_den.back()->Write();
    }

    cout<<"Z/W"<<endl;
    outputFile.cd();
    if(not outputFile.GetDirectory("TF_WZ_EWK"))
      outputFile.mkdir("TF_WZ_EWK");
    outputFile.cd("TF_WZ_EWK");
    zwjcorhist.back()->Write();

    if(addHistoForCutAndCount){
      zwjcorhist_num.back()->Write();
      zwjcorhist_den.back()->Write();
    }

    outputFile.cd();

    if(addWgamma){
      if(not outputFile.GetDirectory("TF_GJ_EWK"))
	outputFile.mkdir("TF_GJ_EWK");
      outputFile.cd("TF_GJ_EWK");
      if(gamcorhist.size() !=0 and gamcorhist.back() != NULL)
	gamcorhist.back()->Write();
      
      if(addHistoForCutAndCount){
	if(gamcorhist_num.back())
	  gamcorhist_num.back()->Write();
	if(gamcorhist_den.back())
	  gamcorhist_den.back()->Write();
      }
      
      cout<<"WGamma+jets"<<endl;
      outputFile.cd();
      if(not outputFile.GetDirectory("TF_WG_EWK"))
        outputFile.mkdir("TF_WG_EWK");
      outputFile.cd("TF_WG_EWK");
      if(wgamcorhist.size() != 0 and wgamcorhist.back() != NULL)
	wgamcorhist.back()->Write();
      
      if(addHistoForCutAndCount){
	if(wgamcorhist_num.back())
	  wgamcorhist_num.back()->Write();
	if(wgamcorhist_den.back())
	  wgamcorhist_den.back()->Write();
      }
    }
    outputFile.cd();
  }

  // clear vectors
  zmmcorhist.clear(); zeecorhist.clear(); wmncorhist.clear(); wencorhist.clear(); zwjcorhist.clear(); gamcorhist.clear(); wgamcorhist.clear();
  zmmcorhist_num.clear(); zeecorhist_num.clear(); wmncorhist_num.clear(); wencorhist_num.clear(); zwjcorhist_num.clear(); gamcorhist_num.clear(); wgamcorhist_num.clear();
  zmmcorhist_den.clear(); zeecorhist_den.clear(); wmncorhist_den.clear(); wencorhist_den.clear(); zwjcorhist_den.clear(); gamcorhist_den.clear(); wgamcorhist_den.clear();


  // Close files
  if(zmmcorfile) zmmcorfile->Close();
  if(zeecorfile) zeecorfile->Close();
  if(wmncorfile) wmncorfile->Close();
  if(wencorfile) wencorfile->Close();
  if(zwjcorfile) zwjcorfile->Close();
  if(gamcorfile) gamcorfile->Close();
  if(wgamcorfile) wgamcorfile->Close();
  
}
