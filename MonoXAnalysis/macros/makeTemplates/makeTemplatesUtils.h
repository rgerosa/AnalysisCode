
bool hasSameSign (TH1* histo1){
  bool sign = true;
  bool firstsign = false;
  for(int iBin = 0; iBin < histo1->GetNbinsX()+1; iBin++){
    if(iBin == 0 and histo1->GetBinContent(iBin+1) > 0){
      firstsign = true;
      continue;
    }
    else if(iBin == 0 and histo1->GetBinContent(iBin+1) < 0){
      firstsign = false;
      continue;
    }    
    else if(iBin != 0 and firstsign and histo1->GetBinContent(iBin+1) < 0){
      sign = false;
    }
    else if(iBin != 0 and not firstsign and histo1->GetBinContent(iBin+1) > 0){
      sign = false;
    }
  }
  return sign;    
}

///////////////////////////// from TFs files to the template one                                                                                                                                      
void fillAndSaveCorrQCDHistograms(const vector<string> & observables, // observables to be considered
				  TFile & outputFile, // output file
				  const string & outDir,  // output directory
				  const Category & category,
				  const bool & addZWratio,
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
  TFile* gamcorfile  = NULL; 
  TFile* wgamcorfile = NULL;

  if(addZgamma){
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
  TFile* gamcorqcdscaleupfile = NULL;
  TFile* gamcorqcdscaledwfile = NULL;
  TFile* gamcorqcdshapeupfile = NULL;
  TFile* gamcorqcdshapedwfile = NULL;
  TFile* gamcorqcdprocupfile = NULL;
  TFile* gamcorqcdprocdwfile = NULL;
  TFile* gamcornnloewkupfile = NULL;
  TFile* gamcornnloewkdwfile = NULL;
  TFile* gamcorsudakovupfile_1 = NULL;
  TFile* gamcorsudakovdwfile_1 = NULL;
  TFile* gamcorsudakovupfile_2 = NULL;
  TFile* gamcorsudakovdwfile_2 = NULL;
  TFile* gamcornnlomissupfile_1 = NULL;
  TFile* gamcornnlomissdwfile_1 = NULL;
  TFile* gamcornnlomissupfile_2 = NULL;
  TFile* gamcornnlomissdwfile_2 = NULL;
  TFile* gamcormixupfile = NULL;
  TFile* gamcormixdwfile = NULL;

  if(addZgamma){
    if(not useNewTheoryUncertainty){
      gamcorqcdfile = TFile::Open((outDir+"/gamcorqcd"+ext+".root").c_str());
      gamcorewkfile = TFile::Open((outDir+"/gamcorewk"+ext+".root").c_str());
      gamcorre1file = TFile::Open((outDir+"/gamcorre1"+ext+".root").c_str());
      gamcorfa1file = TFile::Open((outDir+"/gamcorfa1"+ext+".root").c_str());
      gamcorre2file = TFile::Open((outDir+"/gamcorre2"+ext+".root").c_str());
      gamcorfa2file = TFile::Open((outDir+"/gamcorfa2"+ext+".root").c_str());
      gamcorpdffile = TFile::Open((outDir+"/gamcorpdf"+ext+".root").c_str());
      gamcorfpcfile = TFile::Open((outDir+"/gamcorfpc"+ext+".root").c_str());
    }
    else{
      gamcorewkfile = TFile::Open((outDir+"/gamcorewk"+ext+".root").c_str());
      gamcorqcdfile = TFile::Open((outDir+"/gamcorqcd"+ext+".root").c_str());
      gamcorpdffile = TFile::Open((outDir+"/gamcorpdf"+ext+".root").c_str());
      gamcorqcdscaleupfile = TFile::Open((outDir+"/gamcorqcdscale_up"+ext+".root").c_str());
      gamcorqcdscaledwfile = TFile::Open((outDir+"/gamcorqcdscale_dw"+ext+".root").c_str());
      gamcorqcdshapeupfile = TFile::Open((outDir+"/gamcorqcdshape_up"+ext+".root").c_str());
      gamcorqcdshapedwfile = TFile::Open((outDir+"/gamcorqcdshape_dw"+ext+".root").c_str());
      gamcorqcdprocupfile = TFile::Open((outDir+"/gamcorqcdproc_up"+ext+".root").c_str());
      gamcorqcdprocdwfile = TFile::Open((outDir+"/gamcorqcdproc_dw"+ext+".root").c_str());
      gamcornnloewkupfile = TFile::Open((outDir+"/gamcornnloewk_up"+ext+".root").c_str());
      gamcornnloewkdwfile = TFile::Open((outDir+"/gamcornnloewk_dw"+ext+".root").c_str());
      gamcorsudakovupfile_1 = TFile::Open((outDir+"/gamcorsudakov_up_1"+ext+".root").c_str());
      gamcorsudakovdwfile_1 = TFile::Open((outDir+"/gamcorsudakov_dw_1"+ext+".root").c_str());
      gamcorsudakovupfile_2 = TFile::Open((outDir+"/gamcorsudakov_up_2"+ext+".root").c_str());
      gamcorsudakovdwfile_2 = TFile::Open((outDir+"/gamcorsudakov_dw_2"+ext+".root").c_str());
      gamcornnlomissupfile_1 = TFile::Open((outDir+"/gamcornnlomiss_up_1"+ext+".root").c_str());
      gamcornnlomissdwfile_1 = TFile::Open((outDir+"/gamcornnlomiss_dw_1"+ext+".root").c_str());
      gamcornnlomissupfile_2 = TFile::Open((outDir+"/gamcornnlomiss_up_2"+ext+".root").c_str());
      gamcornnlomissdwfile_2 = TFile::Open((outDir+"/gamcornnlomiss_dw_2"+ext+".root").c_str());
      gamcormixupfile = TFile::Open((outDir+"/gamcormix_up"+ext+".root").c_str());
      gamcormixdwfile = TFile::Open((outDir+"/gamcormix_dw"+ext+".root").c_str());
    }
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
  TFile* wgamcorqcdscaleupfile = NULL;
  TFile* wgamcorqcdscaledwfile = NULL;
  TFile* wgamcorqcdshapeupfile = NULL;
  TFile* wgamcorqcdshapedwfile = NULL;
  TFile* wgamcorqcdprocupfile = NULL;
  TFile* wgamcorqcdprocdwfile = NULL;
  TFile* wgamcornnloewkupfile = NULL;
  TFile* wgamcornnloewkdwfile = NULL;
  TFile* wgamcorsudakovupfile_1 = NULL;
  TFile* wgamcorsudakovdwfile_1 = NULL;
  TFile* wgamcorsudakovupfile_2 = NULL;
  TFile* wgamcorsudakovdwfile_2 = NULL;
  TFile* wgamcornnlomissupfile_1 = NULL;
  TFile* wgamcornnlomissdwfile_1 = NULL;
  TFile* wgamcornnlomissupfile_2 = NULL;
  TFile* wgamcornnlomissdwfile_2 = NULL;
  TFile* wgamcormixupfile = NULL;
  TFile* wgamcormixdwfile = NULL;

  if(addWgamma){
    if(not useNewTheoryUncertainty){
      wgamcorqcdfile = TFile::Open((outDir+"/wgamcorqcd"+ext+".root").c_str());
      wgamcorewkfile = TFile::Open((outDir+"/wgamcorewk"+ext+".root").c_str());
      wgamcorre1file = TFile::Open((outDir+"/wgamcorre1"+ext+".root").c_str());
      wgamcorfa1file = TFile::Open((outDir+"/wgamcorfa1"+ext+".root").c_str());
      wgamcorre2file = TFile::Open((outDir+"/wgamcorre2"+ext+".root").c_str());
      wgamcorfa2file = TFile::Open((outDir+"/wgamcorfa2"+ext+".root").c_str());
      wgamcorpdffile = TFile::Open((outDir+"/wgamcorpdf"+ext+".root").c_str());
      wgamcorfpcfile = TFile::Open((outDir+"/wgamcorfpc"+ext+".root").c_str());
    }
    else{
      wgamcorewkfile = TFile::Open((outDir+"/wgamcorewk"+ext+".root").c_str());
      wgamcorqcdfile = TFile::Open((outDir+"/wgamcorqcd"+ext+".root").c_str());
      wgamcorpdffile = TFile::Open((outDir+"/wgamcorpdf"+ext+".root").c_str());
      wgamcorqcdscaleupfile = TFile::Open((outDir+"/wgamcorqcdscale_up"+ext+".root").c_str());
      wgamcorqcdscaledwfile = TFile::Open((outDir+"/wgamcorqcdscale_dw"+ext+".root").c_str());
      wgamcorqcdshapeupfile = TFile::Open((outDir+"/wgamcorqcdshape_up"+ext+".root").c_str());
      wgamcorqcdshapedwfile = TFile::Open((outDir+"/wgamcorqcdshape_dw"+ext+".root").c_str());
      wgamcorqcdprocupfile = TFile::Open((outDir+"/wgamcorqcdproc_up"+ext+".root").c_str());
      wgamcorqcdprocdwfile = TFile::Open((outDir+"/wgamcorqcdproc_dw"+ext+".root").c_str());
      wgamcornnloewkupfile = TFile::Open((outDir+"/wgamcornnloewk_up"+ext+".root").c_str());
      wgamcornnloewkdwfile = TFile::Open((outDir+"/wgamcornnloewk_dw"+ext+".root").c_str());
      wgamcorsudakovupfile_1 = TFile::Open((outDir+"/wgamcorsudakov_up_1"+ext+".root").c_str());
      wgamcorsudakovdwfile_1 = TFile::Open((outDir+"/wgamcorsudakov_dw_1"+ext+".root").c_str());
      wgamcorsudakovupfile_2 = TFile::Open((outDir+"/wgamcorsudakov_up_2"+ext+".root").c_str());
      wgamcorsudakovdwfile_2 = TFile::Open((outDir+"/wgamcorsudakov_dw_2"+ext+".root").c_str());
      wgamcornnlomissupfile_1 = TFile::Open((outDir+"/wgamcornnlomiss_up_1"+ext+".root").c_str());
      wgamcornnlomissdwfile_1 = TFile::Open((outDir+"/wgamcornnlomiss_dw_1"+ext+".root").c_str());
      wgamcornnlomissupfile_2 = TFile::Open((outDir+"/wgamcornnlomiss_up_2"+ext+".root").c_str());
      wgamcornnlomissdwfile_2 = TFile::Open((outDir+"/wgamcornnlomiss_dw_2"+ext+".root").c_str());
      wgamcormixupfile = TFile::Open((outDir+"/wgamcormix_up"+ext+".root").c_str());
      wgamcormixdwfile = TFile::Open((outDir+"/wgamcormix_dw"+ext+".root").c_str());
    }
  }


  // QCD, EWK, factm re and footprint on Z/W                                                                                                                                                          
  TFile* zwjcorqcdfile = NULL;
  TFile* zwjcorewkfile = NULL;
  TFile* zwjcorre1file = NULL;
  TFile* zwjcorfa1file = NULL;
  TFile* zwjcorre2file = NULL;
  TFile* zwjcorfa2file = NULL;
  TFile* zwjcorpdffile = NULL;

  TFile* zwjcorqcdscaleupfile = NULL;
  TFile* zwjcorqcdscaledwfile = NULL;
  TFile* zwjcorqcdshapeupfile = NULL;
  TFile* zwjcorqcdshapedwfile = NULL;
  TFile* zwjcorqcdprocupfile = NULL;
  TFile* zwjcorqcdprocdwfile = NULL;

  TFile* zwjcornnloewkupfile = NULL;
  TFile* zwjcornnloewkdwfile = NULL;
  TFile* zwjcornnlomissupfile_1 = NULL;
  TFile* zwjcornnlomissdwfile_1 = NULL;
  TFile* zwjcornnlomissupfile_2 = NULL;
  TFile* zwjcornnlomissdwfile_2 = NULL;
  TFile* zwjcorsudakovupfile_1 = NULL;
  TFile* zwjcorsudakovdwfile_1 = NULL;
  TFile* zwjcorsudakovupfile_2 = NULL;
  TFile* zwjcorsudakovdwfile_2 = NULL;

  TFile* zwjcormixupfile = NULL;
  TFile* zwjcormixdwfile = NULL;

  if(addZWratio){
    if(not useNewTheoryUncertainty){
      zwjcorqcdfile = TFile::Open((outDir+"/zwjcorqcd"+ext+".root").c_str());
      zwjcorewkfile = TFile::Open((outDir+"/zwjcorewk"+ext+".root").c_str());
      zwjcorre1file = TFile::Open((outDir+"/zwjcorre1"+ext+".root").c_str());
      zwjcorfa1file = TFile::Open((outDir+"/zwjcorfa1"+ext+".root").c_str());
      zwjcorre2file = TFile::Open((outDir+"/zwjcorre2"+ext+".root").c_str());
      zwjcorfa2file = TFile::Open((outDir+"/zwjcorfa2"+ext+".root").c_str());
      zwjcorpdffile = TFile::Open((outDir+"/zwjcorpdf"+ext+".root").c_str());
    }
    else{
      zwjcorqcdfile = TFile::Open((outDir+"/zwjcorqcd"+ext+".root").c_str());
      zwjcorewkfile = TFile::Open((outDir+"/zwjcorewk"+ext+".root").c_str());
      
      zwjcorqcdscaleupfile = TFile::Open((outDir+"/zwjcorqcd_scaleup"+ext+".root").c_str());
      zwjcorqcdscaledwfile = TFile::Open((outDir+"/zwjcorqcd_scaledw"+ext+".root").c_str());
      
      zwjcorqcdshapeupfile = TFile::Open((outDir+"/zwjcorqcdshape_up"+ext+".root").c_str());
      zwjcorqcdshapedwfile = TFile::Open((outDir+"/zwjcorqcdshape_dw"+ext+".root").c_str());
      
      zwjcorqcdprocupfile  = TFile::Open((outDir+"/zwjcorqcdproc_up"+ext+".root").c_str());
      zwjcorqcdprocdwfile  = TFile::Open((outDir+"/zwjcorqcdproc_dw"+ext+".root").c_str());
      
      zwjcornnloewkupfile = TFile::Open((outDir+"/zwjcornnloewk_up"+ext+".root").c_str());
      zwjcornnloewkdwfile = TFile::Open((outDir+"/zwjcornnloewk_dw"+ext+".root").c_str());
      
      zwjcornnlomissupfile_1 = TFile::Open((outDir+"/zwjcornnlomiss_up_1"+ext+".root").c_str());
      zwjcornnlomissupfile_2 = TFile::Open((outDir+"/zwjcornnlomiss_up_2"+ext+".root").c_str());
      zwjcornnlomissdwfile_1 = TFile::Open((outDir+"/zwjcornnlomiss_dw_1"+ext+".root").c_str());
      zwjcornnlomissdwfile_2 = TFile::Open((outDir+"/zwjcornnlomiss_dw_2"+ext+".root").c_str());
      
      zwjcorsudakovupfile_1 = TFile::Open((outDir+"/zwjcorsudakov_up_1"+ext+".root").c_str());
      zwjcorsudakovupfile_2 = TFile::Open((outDir+"/zwjcorsudakov_up_2"+ext+".root").c_str());
      zwjcorsudakovdwfile_1 = TFile::Open((outDir+"/zwjcorsudakov_dw_1"+ext+".root").c_str());
      zwjcorsudakovdwfile_2 = TFile::Open((outDir+"/zwjcorsudakov_dw_2"+ext+".root").c_str());
      
      zwjcormixupfile = TFile::Open((outDir+"/zwjcormix_up"+ext+".root").c_str());
      zwjcormixdwfile = TFile::Open((outDir+"/zwjcormix_dw"+ext+".root").c_str());
      zwjcorpdffile   = TFile::Open((outDir+"/zwjcorpdf"+ext+".root").c_str());
    }
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
  vector<TH1*> zmmcorhist, zeecorhist, wmncorhist, wencorhist, zwjcorhist, gamcorhist, wgamcorhist, topmucorhist, topelcorhist;

  ////
  vector<TH1*> zmmcorhist_num, zeecorhist_num, wmncorhist_num, wencorhist_num, gamcorhist_num, wgamcorhist_num, zwjcorhist_num;
  vector<TH1*> zmmcorhist_den, zeecorhist_den, wmncorhist_den, wencorhist_den, gamcorhist_den, wgamcorhist_den, zwjcorhist_den;
  vector<TH1*> topmucorhist_num, topelcorhist_num, topmucorhist_den, topelcorhist_den;

  /// Z/Gamma ratio
  vector<TH1*> gamcorewkhist, gamcorqcdhist, gamcorre1hist, gamcorfa1hist, gamcorre2hist, gamcorfa2hist, gamcorpdfhist, gamcorfpchist;
  vector<TH1*> gamcorewkhist_num, gamcorqcdhist_num, gamcorre1hist_num, gamcorfa1hist_num, gamcorre2hist_num, gamcorpdfhist_num, gamcorfa2hist_num, gamcorfpchist_num;
  vector<TH1*> gamcorewkhist_den, gamcorqcdhist_den, gamcorre1hist_den, gamcorfa1hist_den, gamcorre2hist_den, gamcorpdfhist_den, gamcorfa2hist_den, gamcorfpchist_den;

  vector<TH1*> gamcorqcdscaleuphist, gamcorqcdscaledwhist, gamcorqcdshapeuphist, gamcorqcdshapedwhist, gamcorqcdprocuphist, gamcorqcdprocdwhist;
  vector<TH1*> gamcornnloewkuphist, gamcornnloewkdwhist, gamcornnlomissuphist_1, gamcornnlomissuphist_2, gamcornnlomissdwhist_1, gamcornnlomissdwhist_2;
  vector<TH1*> gamcorsudakovuphist_1, gamcorsudakovdwhist_1, gamcorsudakovuphist_2, gamcorsudakovdwhist_2;
  vector<TH1*> gamcormixuphist, gamcormixdwhist;

  vector<TH1*> gamcorqcdscaleuphist_num, gamcorqcdscaledwhist_num, gamcorqcdshapeuphist_num,   gamcorqcdshapedwhist_num,   gamcorqcdprocuphist_num,    gamcorqcdprocdwhist_num;
  vector<TH1*> gamcornnloewkuphist_num,  gamcornnloewkdwhist_num,  gamcornnlomissuphist_1_num, gamcornnlomissuphist_2_num, gamcornnlomissdwhist_1_num, gamcornnlomissdwhist_2_num;
  vector<TH1*> gamcorsudakovuphist_1_num, gamcorsudakovdwhist_1_num, gamcorsudakovuphist_2_num,   gamcorsudakovdwhist_2_num;
  vector<TH1*> gamcormixuphist_num, gamcormixdwhist_num;

  vector<TH1*> gamcorqcdscaleuphist_den, gamcorqcdscaledwhist_den, gamcorqcdshapeuphist_den,   gamcorqcdshapedwhist_den,   gamcorqcdprocuphist_den,    gamcorqcdprocdwhist_den;
  vector<TH1*> gamcornnloewkuphist_den,  gamcornnloewkdwhist_den,  gamcornnlomissuphist_1_den, gamcornnlomissuphist_2_den, gamcornnlomissdwhist_1_den, gamcornnlomissdwhist_2_den;
  vector<TH1*> gamcorsudakovuphist_1_den, gamcorsudakovdwhist_1_den, gamcorsudakovuphist_2_den,   gamcorsudakovdwhist_2_den;
  vector<TH1*> gamcormixuphist_den, gamcormixdwhist_den;

  /// Z/W ratio
  vector<TH1*> zwjcorewkhist, zwjcorqcdhist, zwjcorre1hist, zwjcorre2hist, zwjcorfa1hist, zwjcorfa2hist, zwjcorpdfhist;
  vector<TH1*> zwjcorewkhist_num, zwjcorqcdhist_num, zwjcorre1hist_num, zwjcorre2hist_num, zwjcorfa1hist_num, zwjcorfa2hist_num, zwjcorpdfhist_num;
  vector<TH1*> zwjcorewkhist_den, zwjcorqcdhist_den, zwjcorre1hist_den, zwjcorre2hist_den, zwjcorfa1hist_den, zwjcorfa2hist_den, zwjcorpdfhist_den;

  vector<TH1*> zwjcorqcdscaleuphist, zwjcorqcdscaledwhist, zwjcorqcdshapeuphist, zwjcorqcdshapedwhist, zwjcorqcdprocuphist, zwjcorqcdprocdwhist;
  vector<TH1*> zwjcornnloewkuphist, zwjcornnloewkdwhist, zwjcornnlomissuphist_1, zwjcornnlomissuphist_2, zwjcornnlomissdwhist_1, zwjcornnlomissdwhist_2;
  vector<TH1*> zwjcorsudakovuphist_1, zwjcorsudakovdwhist_1, zwjcorsudakovuphist_2, zwjcorsudakovdwhist_2;
  vector<TH1*> zwjcormixuphist, zwjcormixdwhist;

  vector<TH1*> zwjcorqcdscaleuphist_num, zwjcorqcdscaledwhist_num, zwjcorqcdshapeuphist_num,   zwjcorqcdshapedwhist_num,   zwjcorqcdprocuphist_num,    zwjcorqcdprocdwhist_num;
  vector<TH1*> zwjcornnloewkuphist_num,  zwjcornnloewkdwhist_num,  zwjcornnlomissuphist_1_num, zwjcornnlomissuphist_2_num, zwjcornnlomissdwhist_1_num, zwjcornnlomissdwhist_2_num;
  vector<TH1*> zwjcorsudakovuphist_1_num, zwjcorsudakovdwhist_1_num, zwjcorsudakovuphist_2_num,   zwjcorsudakovdwhist_2_num;
  vector<TH1*> zwjcormixuphist_num, zwjcormixdwhist_num;

  vector<TH1*> zwjcorqcdscaleuphist_den, zwjcorqcdscaledwhist_den, zwjcorqcdshapeuphist_den,   zwjcorqcdshapedwhist_den,   zwjcorqcdprocuphist_den,    zwjcorqcdprocdwhist_den;
  vector<TH1*> zwjcornnloewkuphist_den,  zwjcornnloewkdwhist_den,  zwjcornnlomissuphist_1_den, zwjcornnlomissuphist_2_den, zwjcornnlomissdwhist_1_den, zwjcornnlomissdwhist_2_den;
  vector<TH1*> zwjcorsudakovuphist_1_den, zwjcorsudakovdwhist_1_den, zwjcorsudakovuphist_2_den,   zwjcorsudakovdwhist_2_den;
  vector<TH1*> zwjcormixuphist_den, zwjcormixdwhist_den;
  
  ///W/wgamma
  vector<TH1*> wgamcorewkhist, wgamcorqcdhist, wgamcorre1hist, wgamcorfa1hist, wgamcorre2hist, wgamcorfa2hist, wgamcorpdfhist, wgamcorfpchist;
  vector<TH1*> wgamcorewkhist_num, wgamcorqcdhist_num, wgamcorre1hist_num, wgamcorfa1hist_num, wgamcorre2hist_num, wgamcorpdfhist_num, wgamcorfa2hist_num, wgamcorfpchist_num;
  vector<TH1*> wgamcorewkhist_den, wgamcorqcdhist_den, wgamcorre1hist_den, wgamcorfa1hist_den, wgamcorre2hist_den, wgamcorpdfhist_den, wgamcorfa2hist_den, wgamcorfpchist_den;

  vector<TH1*> wgamcorqcdscaleuphist, wgamcorqcdscaledwhist, wgamcorqcdshapeuphist, wgamcorqcdshapedwhist, wgamcorqcdprocuphist, wgamcorqcdprocdwhist;
  vector<TH1*> wgamcornnloewkuphist, wgamcornnloewkdwhist, wgamcornnlomissuphist_1, wgamcornnlomissuphist_2, wgamcornnlomissdwhist_1, wgamcornnlomissdwhist_2;
  vector<TH1*> wgamcorsudakovuphist_1, wgamcorsudakovdwhist_1, wgamcorsudakovuphist_2, wgamcorsudakovdwhist_2;
  vector<TH1*> wgamcormixuphist, wgamcormixdwhist;

  vector<TH1*> wgamcorqcdscaleuphist_num, wgamcorqcdscaledwhist_num, wgamcorqcdshapeuphist_num,   wgamcorqcdshapedwhist_num,   wgamcorqcdprocuphist_num,    wgamcorqcdprocdwhist_num;
  vector<TH1*> wgamcornnloewkuphist_num,  wgamcornnloewkdwhist_num,  wgamcornnlomissuphist_1_num, wgamcornnlomissuphist_2_num, wgamcornnlomissdwhist_1_num, wgamcornnlomissdwhist_2_num;
  vector<TH1*> wgamcorsudakovuphist_1_num, wgamcorsudakovdwhist_1_num, wgamcorsudakovuphist_2_num,   wgamcorsudakovdwhist_2_num;
  vector<TH1*> wgamcormixuphist_num, wgamcormixdwhist_num;

  vector<TH1*> wgamcorqcdscaleuphist_den, wgamcorqcdscaledwhist_den, wgamcorqcdshapeuphist_den,   wgamcorqcdshapedwhist_den,   wgamcorqcdprocuphist_den,    wgamcorqcdprocdwhist_den;
  vector<TH1*> wgamcornnloewkuphist_den,  wgamcornnloewkdwhist_den,  wgamcornnlomissuphist_1_den, wgamcornnlomissuphist_2_den, wgamcornnlomissdwhist_1_den, wgamcornnlomissdwhist_2_den;
  vector<TH1*> wgamcorsudakovuphist_1_den, wgamcorsudakovdwhist_1_den, wgamcorsudakovuphist_2_den,   wgamcorsudakovdwhist_2_den;
  vector<TH1*> wgamcormixuphist_den, wgamcormixdwhist_den;

  ////
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

    if(addZgamma){
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
    TH1* gamuncqcdscalehist = NULL;
    TH1* gamuncqcdshapehist = NULL;
    TH1* gamuncqcdprochist = NULL;
    TH1* gamuncnnloewkhist = NULL;
    TH1* gamuncsudakovhist_1 = NULL;
    TH1* gamuncsudakovhist_2 = NULL;
    TH1* gamuncnnlomisshist_1 = NULL;
    TH1* gamuncnnlomisshist_2 = NULL;
    TH1* gamuncmixhist = NULL;

    //////
    if(addZgamma){

      if(not useNewTheoryUncertainty){
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
	
	// ZG ewk uncertainty
	gamuncewkhist = (TH1*)gamcorewkhist.back()->Clone(("gamuncewk"+ext+"hist_"+obs).c_str());
	gamuncewkhist->Divide(gamcorqcdhist.back());
	for (int i = 1; i <= gamuncewkhist->GetNbinsX(); i++)
	  gamuncewkhist->SetBinContent(i, fabs(gamuncewkhist->GetBinContent(i)-1.0));
	gamuncewkhist->SetName(("ZG_EWK_"+obs).c_str());
	
	/// ZG renormalization scale                                                                                                                                                           
	gamuncre1hist = (TH1*)gamcorre1hist.back()->Clone(("gamuncre1"+ext+"hist_"+obs).c_str());
	gamuncre1hist->Divide(gamcorqcdhist.back());
	for (int i = 1; i <= gamuncre1hist->GetNbinsX(); i++)
	  gamuncre1hist->SetBinContent(i, fabs(gamuncre1hist->GetBinContent(i)-1.0));
	gamuncre1hist->SetName(("ZG_RenScale1_"+obs).c_str());
	
	/// ZG renormalization scale                                                                                                                                                           
	gamuncfa1hist = (TH1*)gamcorfa1hist.back()->Clone(("gamuncfa1"+ext+"hist_"+obs).c_str());
	gamuncfa1hist->Divide(gamcorqcdhist.back());
	for (int i = 1; i <= gamuncfa1hist->GetNbinsX(); i++)
	  gamuncfa1hist->SetBinContent(i, fabs(gamuncfa1hist->GetBinContent(i)-1.0));
	gamuncfa1hist->SetName(("ZG_FactScale1_"+obs).c_str());
	
	/// ZG renormalization scale                                                                                                                                                           
	gamuncre2hist = (TH1*)gamcorre2hist.back()->Clone(("gamuncre2"+ext+"hist_"+obs).c_str());
	gamuncre2hist->Divide(gamcorqcdhist.back());
	for (int i = 1; i <= gamuncre2hist->GetNbinsX(); i++)
	  gamuncre2hist->SetBinContent(i, fabs(gamuncre2hist->GetBinContent(i)-1.0));
	gamuncre2hist->SetName(("ZG_RenScale2_"+obs).c_str());
      
	/// ZG factorization scale                                                                                                                                                                 
	gamuncfa2hist = (TH1*)gamcorfa2hist.back()->Clone(("gamuncfa2"+ext+"hist_"+obs).c_str());
	gamuncfa2hist->Divide(gamcorqcdhist.back());
	for (int i = 1; i <= gamuncfa2hist->GetNbinsX(); i++)
	  gamuncfa2hist->SetBinContent(i, fabs(gamuncfa2hist->GetBinContent(i)-1.0));
	gamuncfa2hist->SetName(("ZG_FactScale2_"+obs).c_str());
	
	/// ZG pdfs                                                                                                                                                                 
	gamuncpdfhist = (TH1*)gamcorpdfhist.back()->Clone(("gamuncpdf"+ext+"hist_"+obs).c_str());
	gamuncpdfhist->Divide(gamcorqcdhist.back());
	for (int i = 1; i <= gamuncpdfhist->GetNbinsX(); i++)
	  gamuncpdfhist->SetBinContent(i, fabs(gamuncpdfhist->GetBinContent(i)-1.0));
	gamuncpdfhist->SetName(("ZG_PDF_"+obs).c_str());
      
	/// ZG footprint                                                                                                                                                                            
	gamuncfpchist = (TH1*)gamcorfpchist.back()->Clone(("gamuncfpc"+ext+"hist_"+obs).c_str());
	gamuncfpchist->Divide(gamcorqcdhist.back());
	for (int i = 1; i <= gamuncfpchist->GetNbinsX(); i++)
	  gamuncfpchist->SetBinContent(i, fabs(gamuncfpchist->GetBinContent(i)-1.0));
	gamuncfpchist->SetName(("ZG_Footprint_"+obs).c_str());
      }
      else{

	cout<<"Make Z/gamma sys histograms"<<endl;    
	gamcorewkhist.push_back( (TH1*)gamcorewkfile->FindObjectAny(("gamcor"+ext+"ewkhist_"+obs).c_str()));
	gamcorewkhist_num.push_back( (TH1*)gamcorewkfile->FindObjectAny(("nhist_gam_ewk"+ext+"_"+obs).c_str()));
	gamcorewkhist_den.push_back( (TH1*)gamcorewkfile->FindObjectAny(("dhist_gam_ewk"+ext+"_"+obs).c_str()));

	gamcorqcdhist.push_back( (TH1*)gamcorqcdfile->FindObjectAny(("gamcor"+ext+"qcdhist_"+obs).c_str()));
	gamcorqcdhist_num.push_back( (TH1*)gamcorqcdfile->FindObjectAny(("nhist_gam_qcd"+ext+"_"+obs).c_str()));
	gamcorqcdhist_den.push_back( (TH1*)gamcorqcdfile->FindObjectAny(("dhist_gam_qcd"+ext+"_"+obs).c_str()));
	
	gamcorpdfhist.push_back( (TH1*)gamcorpdffile->FindObjectAny(("gamcor"+ext+"pdfhist_"+obs).c_str()));
        gamcorpdfhist_num.push_back( (TH1*)gamcorpdffile->FindObjectAny(("nhist_gam_pdf"+ext+"_"+obs).c_str()));
        gamcorpdfhist_den.push_back( (TH1*)gamcorpdffile->FindObjectAny(("dhist_gam_pdf"+ext+"_"+obs).c_str()));

	gamcorqcdscaleuphist.push_back( (TH1*)gamcorqcdscaleupfile->FindObjectAny(("gamcorqcdscale_up"+ext+"hist_"+obs).c_str()));
	gamcorqcdscaleuphist_num.push_back( (TH1*)gamcorqcdscaleupfile->FindObjectAny(("nhist_gam_qcdscale_up"+ext+"_"+obs).c_str()));
	gamcorqcdscaleuphist_den.push_back( (TH1*)gamcorqcdscaleupfile->FindObjectAny(("dhist_gam_qcdscale_up"+ext+"_"+obs).c_str()));

	gamcorqcdscaledwhist.push_back( (TH1*)gamcorqcdscaledwfile->FindObjectAny(("gamcorqcdscale_dw"+ext+"hist_"+obs).c_str()));
	gamcorqcdscaledwhist_num.push_back( (TH1*)gamcorqcdscaledwfile->FindObjectAny(("nhist_gam_qcdscale_dw"+ext+"_"+obs).c_str()));
	gamcorqcdscaledwhist_den.push_back( (TH1*)gamcorqcdscaledwfile->FindObjectAny(("dhist_gam_qcdscale_dw"+ext+"_"+obs).c_str()));

	gamcorqcdshapeuphist.push_back( (TH1*)gamcorqcdshapeupfile->FindObjectAny(("gamcorqcdshape_up"+ext+"hist_"+obs).c_str()));
	gamcorqcdshapeuphist_num.push_back( (TH1*)gamcorqcdshapeupfile->FindObjectAny(("nhist_gam_qcdshape_up"+ext+"_"+obs).c_str()));
	gamcorqcdshapeuphist_den.push_back( (TH1*)gamcorqcdshapeupfile->FindObjectAny(("dhist_gam_qcdshape_up"+ext+"_"+obs).c_str()));

	gamcorqcdshapedwhist.push_back( (TH1*)gamcorqcdshapedwfile->FindObjectAny(("gamcorqcdshape_dw"+ext+"hist_"+obs).c_str()));
	gamcorqcdshapedwhist_num.push_back( (TH1*)gamcorqcdshapedwfile->FindObjectAny(("nhist_gam_qcdshape_dw"+ext+"_"+obs).c_str()));
	gamcorqcdshapedwhist_den.push_back( (TH1*)gamcorqcdshapedwfile->FindObjectAny(("dhist_gam_qcdshape_dw"+ext+"_"+obs).c_str()));

	gamcorqcdprocuphist.push_back( (TH1*)gamcorqcdprocupfile->FindObjectAny(("gamcorqcdproc_up"+ext+"hist_"+obs).c_str()));
	gamcorqcdprocuphist_num.push_back( (TH1*)gamcorqcdprocupfile->FindObjectAny(("nhist_gam_qcdproc_up"+ext+"_"+obs).c_str()));
	gamcorqcdprocuphist_den.push_back( (TH1*)gamcorqcdprocupfile->FindObjectAny(("dhist_gam_qcdproc_up"+ext+"_"+obs).c_str()));

	gamcorqcdprocdwhist.push_back( (TH1*)gamcorqcdprocdwfile->FindObjectAny(("gamcorqcdproc_dw"+ext+"hist_"+obs).c_str()));
	gamcorqcdprocdwhist_num.push_back( (TH1*)gamcorqcdprocdwfile->FindObjectAny(("nhist_gam_qcdproc_dw"+ext+"_"+obs).c_str()));
	gamcorqcdprocdwhist_den.push_back( (TH1*)gamcorqcdprocdwfile->FindObjectAny(("dhist_gam_qcdproc_dw"+ext+"_"+obs).c_str()));

	gamcornnloewkuphist.push_back( (TH1*)gamcornnloewkupfile->FindObjectAny(("gamcornnloewk_up"+ext+"hist_"+obs).c_str()));
	gamcornnloewkuphist_num.push_back( (TH1*)gamcornnloewkupfile->FindObjectAny(("nhist_gam_nnloewk_up"+ext+"_"+obs).c_str()));
	gamcornnloewkuphist_den.push_back( (TH1*)gamcornnloewkupfile->FindObjectAny(("dhist_gam_nnloewk_up"+ext+"_"+obs).c_str()));

	gamcornnloewkdwhist.push_back( (TH1*)gamcornnloewkdwfile->FindObjectAny(("gamcornnloewk_dw"+ext+"hist_"+obs).c_str()));
	gamcornnloewkdwhist_num.push_back( (TH1*)gamcornnloewkdwfile->FindObjectAny(("nhist_gam_nnloewk_dw"+ext+"_"+obs).c_str()));
	gamcornnloewkdwhist_den.push_back( (TH1*)gamcornnloewkdwfile->FindObjectAny(("dhist_gam_nnloewk_dw"+ext+"_"+obs).c_str()));

	gamcorsudakovuphist_1.push_back( (TH1*)gamcorsudakovupfile_1->FindObjectAny(("gamcorsudakov_up_1"+ext+"hist_"+obs).c_str()));
	gamcorsudakovuphist_1_num.push_back( (TH1*)gamcorsudakovupfile_1->FindObjectAny(("nhist_gam_sudakov_up_1"+ext+"_"+obs).c_str()));
	gamcorsudakovuphist_1_den.push_back( (TH1*)gamcorsudakovupfile_1->FindObjectAny(("dhist_gam_sudakov_up_1"+ext+"_"+obs).c_str()));
	
	gamcorsudakovdwhist_1.push_back( (TH1*)gamcorsudakovdwfile_1->FindObjectAny(("gamcorsudakov_dw_1"+ext+"hist_"+obs).c_str()));
	gamcorsudakovdwhist_1_num.push_back( (TH1*)gamcorsudakovdwfile_1->FindObjectAny(("nhist_gam_sudakov_dw_1"+ext+"_"+obs).c_str()));
	gamcorsudakovdwhist_1_den.push_back( (TH1*)gamcorsudakovdwfile_1->FindObjectAny(("dhist_gam_sudakov_dw_1"+ext+"_"+obs).c_str()));
	
	gamcorsudakovuphist_2.push_back( (TH1*)gamcorsudakovupfile_2->FindObjectAny(("gamcorsudakov_up_2"+ext+"hist_"+obs).c_str()));
	gamcorsudakovuphist_2_num.push_back( (TH1*)gamcorsudakovupfile_2->FindObjectAny(("nhist_gam_sudakov_up_2"+ext+"_"+obs).c_str()));
	gamcorsudakovuphist_2_den.push_back( (TH1*)gamcorsudakovupfile_2->FindObjectAny(("dhist_gam_sudakov_up_2"+ext+"_"+obs).c_str()));

	gamcorsudakovdwhist_2.push_back( (TH1*)gamcorsudakovdwfile_2->FindObjectAny(("gamcorsudakov_dw_2"+ext+"hist_"+obs).c_str()));
	gamcorsudakovdwhist_2_num.push_back( (TH1*)gamcorsudakovdwfile_2->FindObjectAny(("nhist_gam_sudakov_dw_2"+ext+"_"+obs).c_str()));
	gamcorsudakovdwhist_2_den.push_back( (TH1*)gamcorsudakovdwfile_2->FindObjectAny(("dhist_gam_sudakov_dw_2"+ext+"_"+obs).c_str()));
	
	gamcornnlomissuphist_1.push_back( (TH1*)gamcornnlomissupfile_1->FindObjectAny(("gamcornnlomiss_up_1"+ext+"hist_"+obs).c_str()));
	gamcornnlomissuphist_1_num.push_back( (TH1*)gamcornnlomissupfile_1->FindObjectAny(("nhist_gam_nnlomiss_up_1"+ext+"_"+obs).c_str()));
	gamcornnlomissuphist_1_den.push_back( (TH1*)gamcornnlomissupfile_1->FindObjectAny(("dhist_gam_nnlomiss_up_1"+ext+"_"+obs).c_str()));
	
	gamcornnlomissdwhist_1.push_back( (TH1*)gamcornnlomissdwfile_1->FindObjectAny(("gamcornnlomiss_dw_1"+ext+"hist_"+obs).c_str()));
	gamcornnlomissdwhist_1_num.push_back( (TH1*)gamcornnlomissdwfile_1->FindObjectAny(("nhist_gam_nnlomiss_dw_1"+ext+"_"+obs).c_str()));
	gamcornnlomissdwhist_1_den.push_back( (TH1*)gamcornnlomissdwfile_1->FindObjectAny(("dhist_gam_nnlomiss_dw_1"+ext+"_"+obs).c_str()));
	
	gamcornnlomissuphist_2.push_back( (TH1*)gamcornnlomissupfile_2->FindObjectAny(("gamcornnlomiss_up_2"+ext+"hist_"+obs).c_str()));
	gamcornnlomissuphist_2_num.push_back( (TH1*)gamcornnlomissupfile_2->FindObjectAny(("nhist_gam_nnlomiss_up_2"+ext+"_"+obs).c_str()));
	gamcornnlomissuphist_2_den.push_back( (TH1*)gamcornnlomissupfile_2->FindObjectAny(("dhist_gam_nnlomiss_up_2"+ext+"_"+obs).c_str()));
	
	gamcornnlomissdwhist_2.push_back( (TH1*)gamcornnlomissdwfile_2->FindObjectAny(("gamcornnlomiss_dw_2"+ext+"hist_"+obs).c_str()));
	gamcornnlomissdwhist_2_num.push_back( (TH1*)gamcornnlomissdwfile_2->FindObjectAny(("nhist_gam_nnlomiss_dw_2"+ext+"_"+obs).c_str()));
	gamcornnlomissdwhist_2_den.push_back( (TH1*)gamcornnlomissdwfile_2->FindObjectAny(("dhist_gam_nnlomiss_dw_2"+ext+"_"+obs).c_str()));
	
	gamcormixuphist.push_back( (TH1*)gamcormixupfile->FindObjectAny(("gamcormix_up"+ext+"hist_"+obs).c_str()));
	gamcormixuphist_num.push_back( (TH1*)gamcormixupfile->FindObjectAny(("nhist_gam_mix_up"+ext+"_"+obs).c_str()));
	gamcormixuphist_den.push_back( (TH1*)gamcormixupfile->FindObjectAny(("dhist_gam_mix_up"+ext+"_"+obs).c_str()));

	gamcormixdwhist.push_back( (TH1*)gamcormixdwfile->FindObjectAny(("gamcormix_dw"+ext+"hist_"+obs).c_str()));
	gamcormixdwhist_num.push_back( (TH1*)gamcormixdwfile->FindObjectAny(("nhist_gam_mix_dw"+ext+"_"+obs).c_str()));
	gamcormixdwhist_den.push_back( (TH1*)gamcormixdwfile->FindObjectAny(("dhist_gam_mix_dw"+ext+"_"+obs).c_str()));

	// QCD scale --> to be symmetrized
	gamuncqcdscalehist = (TH1*) gamcorqcdscaleuphist.back()->Clone(("gamuncqcdscale"+ext+"hist_"+obs).c_str());	
	gamuncqcdscalehist->Reset("ICES");
	// fill first
	for (int i = 0; i <= gamuncqcdscalehist->GetNbinsX()+1; i++)
	  gamuncqcdscalehist->SetBinContent(i,(gamcorqcdscaleuphist.back()->GetBinContent(i)-gamcorqcdscaledwhist.back()->GetBinContent(i))/(2*gamcorewkhist.back()->GetBinContent(i)));
	// smooth
	gamuncqcdscalehist->Smooth(1,"R");
	// evaluate sign
	bool sign = hasSameSign(gamuncqcdscalehist);
	if(sign){
	  for (int i = 0; i <= gamuncqcdscalehist->GetNbinsX()+1; i++)
	    gamuncqcdscalehist->SetBinContent(i,fabs(gamuncqcdscalehist->GetBinContent(i)));
	}
	gamuncqcdscalehist->SetName(("ZG_QCDScale_"+obs).c_str());
	
	// QCD Shape --> to be symmetrized
	gamuncqcdshapehist = (TH1*) gamcorqcdshapeuphist.back()->Clone(("gamuncqcdshape"+ext+"hist_"+obs).c_str());
	gamuncqcdshapehist->Reset("ICES");
	for (int i = 0; i <= gamuncqcdshapehist->GetNbinsX()+1; i++)
	    gamuncqcdshapehist->SetBinContent(i,(gamcorqcdshapeuphist.back()->GetBinContent(i)-gamcorqcdshapedwhist.back()->GetBinContent(i))/(2*gamcorewkhist.back()->GetBinContent(i)));
	gamuncqcdshapehist->Smooth(1,"R");
	sign = hasSameSign(gamuncqcdshapehist);
	if(sign){
	  for (int i = 0; i <= gamuncqcdshapehist->GetNbinsX()+1; i++)
	    gamuncqcdshapehist->SetBinContent(i,fabs(gamuncqcdshapehist->GetBinContent(i)));
	}	
	gamuncqcdshapehist->SetName(("ZG_QCDShape_"+obs).c_str());

	// QCD Process --> to be symmetrized
	gamuncqcdprochist = (TH1*) gamcorqcdprocuphist.back()->Clone(("gamuncqcdproc"+ext+"hist_"+obs).c_str());
	gamuncqcdprochist->Reset("ICES");
	for (int i = 0; i <= gamuncqcdprochist->GetNbinsX()+1; i++)
	  gamuncqcdprochist->SetBinContent(i,(gamcorqcdprocuphist.back()->GetBinContent(i)-gamcorqcdprocdwhist.back()->GetBinContent(i))/(2*gamcorewkhist.back()->GetBinContent(i)));
	gamuncqcdprochist->Smooth(1,"R");	
	sign = hasSameSign(gamuncqcdprochist);
	if(sign){
	  for (int i = 0; i <= gamuncqcdprochist->GetNbinsX()+1; i++)
	    gamuncqcdshapehist->SetBinContent(i,fabs(gamuncqcdprochist->GetBinContent(i)));
	}
	gamuncqcdprochist->SetName(("ZG_QCDProcess_"+obs).c_str());

	// NNLO EWK --> to be symmetrized
	gamuncnnloewkhist = (TH1*) gamcornnloewkuphist.back()->Clone(("gamuncnnloewk"+ext+"hist_"+obs).c_str());
	gamuncnnloewkhist->Reset("ICES");
	for (int i = 0; i <= gamuncnnloewkhist->GetNbinsX()+1; i++){
	  gamuncnnloewkhist->SetBinContent(i,fabs(gamcornnloewkuphist.back()->GetBinContent(i)-gamcornnloewkdwhist.back()->GetBinContent(i))/(2*gamcorewkhist.back()->GetBinContent(i)));
	}	
	gamuncnnloewkhist->Smooth(1,"R");
	gamuncnnloewkhist->SetName(("ZG_NNLOEWK_"+obs).c_str());
	
	// NLO EWK Sud --> to be symmetrized
	gamuncsudakovhist_1 = (TH1*) gamcorsudakovuphist_1.back()->Clone(("gamuncsudakov1"+ext+"hist_"+obs).c_str());
	gamuncsudakovhist_1->Reset("ICES");
	for (int i = 0; i <= gamuncsudakovhist_1->GetNbinsX()+1; i++){
	  gamuncsudakovhist_1->SetBinContent(i,fabs(gamcorsudakovuphist_1.back()->GetBinContent(i)-gamcorsudakovdwhist_1.back()->GetBinContent(i))/(2*gamcorewkhist.back()->GetBinContent(i)));
	}
	gamuncsudakovhist_1->Smooth(1,"R");
	gamuncsudakovhist_1->SetName(("ZG_Sudakov1_"+obs).c_str());

	// NLO EWK Sud --> to be symmetrized
	gamuncsudakovhist_2 = (TH1*) gamcorsudakovuphist_2.back()->Clone(("gamuncsudakov2"+ext+"hist_"+obs).c_str());
	gamuncsudakovhist_2->Reset("ICES");
	for (int i = 0; i <= gamuncsudakovhist_2->GetNbinsX()+1; i++){
	  gamuncsudakovhist_2->SetBinContent(i,fabs(gamcorsudakovuphist_2.back()->GetBinContent(i)-gamcorsudakovdwhist_2.back()->GetBinContent(i))/(2*gamcorewkhist.back()->GetBinContent(i)));
	}
	gamuncsudakovhist_2->Smooth(1,"R");
	gamuncsudakovhist_2->SetName(("ZG_Sudakov2_"+obs).c_str());

	// NNLO Miss --> to be symmetrized
	gamuncnnlomisshist_1 = (TH1*) gamcornnlomissuphist_1.back()->Clone(("gamuncnnlomiss1"+ext+"hist_"+obs).c_str());
	gamuncnnlomisshist_1->Reset("ICES");
	for (int i = 0; i <= gamuncnnlomisshist_1->GetNbinsX()+1; i++){
	  gamuncnnlomisshist_1->SetBinContent(i,fabs(gamcornnlomissuphist_1.back()->GetBinContent(i)-gamcornnlomissdwhist_1.back()->GetBinContent(i))/(2*gamcorewkhist.back()->GetBinContent(i)));
	}
	gamuncnnlomisshist_1->Smooth(1,"R");
	gamuncnnlomisshist_1->SetName(("ZG_NNLOMiss1_"+obs).c_str());
	
	// NNLO Miss --> to be symmetrized
	gamuncnnlomisshist_2 = (TH1*) gamcornnlomissuphist_2.back()->Clone(("gamuncnnlomiss2"+ext+"hist_"+obs).c_str());
	gamuncnnlomisshist_2->Reset("ICES");
	for (int i = 0; i <= gamuncnnlomisshist_2->GetNbinsX()+1; i++){
	  gamuncnnlomisshist_2->SetBinContent(i,fabs(gamcornnlomissuphist_2.back()->GetBinContent(i)-gamcornnlomissdwhist_2.back()->GetBinContent(i))/(2*gamcorewkhist.back()->GetBinContent(i)));
	}
	gamuncnnlomisshist_2->Smooth(1,"R");
	gamuncnnlomisshist_2->SetName(("ZG_NNLOMiss2_"+obs).c_str());

	// NLO QCD+EWK --> to be symmetrized
	gamuncmixhist = (TH1*) gamcormixuphist.back()->Clone(("gamuncmix"+ext+"hist_"+obs).c_str());
	gamuncmixhist->Reset("ICES");
	for (int i = 0; i <= gamuncmixhist->GetNbinsX()+1; i++){
	  gamuncmixhist->SetBinContent(i,fabs(gamcormixuphist.back()->GetBinContent(i)-gamcormixdwhist.back()->GetBinContent(i))/(2*gamcorewkhist.back()->GetBinContent(i)));
	}
	gamuncmixhist->Smooth(1,"R");
	gamuncmixhist->SetName(("ZG_MIX_"+obs).c_str());

	// NLO PDF
	gamuncpdfhist = (TH1*)gamcorpdfhist.back()->Clone(("gamuncpdf"+ext+"hist_"+obs).c_str());
	gamuncpdfhist->Divide(gamcorqcdhist.back());
	for (int i = 1; i <= gamuncpdfhist->GetNbinsX(); i++)
	  gamuncpdfhist->SetBinContent(i,fabs(gamuncpdfhist->GetBinContent(i)-1.0));
	gamuncpdfhist->SetName(("ZG_PDF_"+obs).c_str());
      }
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
    TH1* zwjuncqcdshapehist = NULL;
    TH1* zwjuncqcdprochist = NULL;
    TH1* zwjuncnnloewkhist = NULL;
    TH1* zwjuncnnlomisshist_1 = NULL;
    TH1* zwjuncnnlomisshist_2 = NULL;
    TH1* zwjuncsudakovhist_1 = NULL;
    TH1* zwjuncsudakovhist_2 = NULL;
    TH1* zwjuncmixhist = NULL;

    if(addZWratio){
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
	zwjcorqcdscaleuphist_num.push_back( (TH1*)zwjcorqcdscaleupfile->FindObjectAny(("nhist_zwj_qcd_scaleup"+ext+"_"+obs).c_str()));
	zwjcorqcdscaleuphist_den.push_back( (TH1*)zwjcorqcdscaleupfile->FindObjectAny(("dhist_zwj_qcd_scaleup"+ext+"_"+obs).c_str()));

	zwjcorqcdscaledwhist.push_back( (TH1*)zwjcorqcdscaledwfile->FindObjectAny(("zwjcorqcd_scaledw"+ext+"hist_"+obs).c_str()));
	zwjcorqcdscaledwhist_num.push_back( (TH1*)zwjcorqcdscaledwfile->FindObjectAny(("nhist_zwj_qcd_scaledw"+ext+"_"+obs).c_str()));
	zwjcorqcdscaledwhist_den.push_back( (TH1*)zwjcorqcdscaledwfile->FindObjectAny(("dhist_zwj_qcd_scaledw"+ext+"_"+obs).c_str()));
	
	zwjcorqcdshapeuphist.push_back( (TH1*)zwjcorqcdshapeupfile->FindObjectAny(("zwjcorqcdshape_up"+ext+"hist_"+obs).c_str()));
	zwjcorqcdshapeuphist_num.push_back( (TH1*)zwjcorqcdshapeupfile->FindObjectAny(("nhist_zwj_qcdshape_up"+ext+"_"+obs).c_str()));
	zwjcorqcdshapeuphist_den.push_back( (TH1*)zwjcorqcdshapeupfile->FindObjectAny(("dhist_zwj_qcdshape_up"+ext+"_"+obs).c_str()));
	
	zwjcorqcdshapedwhist.push_back( (TH1*)zwjcorqcdshapedwfile->FindObjectAny(("zwjcorqcdshape_dw"+ext+"hist_"+obs).c_str()));
	zwjcorqcdshapedwhist_num.push_back( (TH1*)zwjcorqcdshapedwfile->FindObjectAny(("nhist_zwj_qcdshape_dw"+ext+"_"+obs).c_str()));
	zwjcorqcdshapedwhist_den.push_back( (TH1*)zwjcorqcdshapedwfile->FindObjectAny(("dhist_zwj_qcdshape_dw"+ext+"_"+obs).c_str()));
	
	zwjcorqcdprocuphist.push_back( (TH1*)zwjcorqcdprocupfile->FindObjectAny(("zwjcorqcdproc_up"+ext+"hist_"+obs).c_str()));
	zwjcorqcdprocuphist_num.push_back( (TH1*)zwjcorqcdprocupfile->FindObjectAny(("nhist_zwj_qcdproc_up"+ext+"_"+obs).c_str()));
	zwjcorqcdprocuphist_den.push_back( (TH1*)zwjcorqcdprocupfile->FindObjectAny(("dhist_zwj_qcdproc_up"+ext+"_"+obs).c_str()));
	
	zwjcorqcdprocdwhist.push_back( (TH1*)zwjcorqcdprocdwfile->FindObjectAny(("zwjcorqcdproc_dw"+ext+"hist_"+obs).c_str()));
	zwjcorqcdprocdwhist_num.push_back( (TH1*)zwjcorqcdprocdwfile->FindObjectAny(("nhist_zwj_qcdproc_dw"+ext+"_"+obs).c_str()));
	zwjcorqcdprocdwhist_den.push_back( (TH1*)zwjcorqcdprocdwfile->FindObjectAny(("dhist_zwj_qcdproc_dw"+ext+"_"+obs).c_str()));
	
	zwjcornnloewkuphist.push_back( (TH1*)zwjcornnloewkupfile->FindObjectAny(("zwjcornnloewk_up"+ext+"hist_"+obs).c_str()));
	zwjcornnloewkuphist_num.push_back( (TH1*)zwjcornnloewkupfile->FindObjectAny(("nhist_zwj_nnloewk_up"+ext+"_"+obs).c_str()));
	zwjcornnloewkuphist_den.push_back( (TH1*)zwjcornnloewkupfile->FindObjectAny(("dhist_zwj_nnloewk_up"+ext+"_"+obs).c_str()));
	
	zwjcornnloewkdwhist.push_back( (TH1*)zwjcornnloewkdwfile->FindObjectAny(("zwjcornnloewk_dw"+ext+"hist_"+obs).c_str()));
	zwjcornnloewkdwhist_num.push_back( (TH1*)zwjcornnloewkdwfile->FindObjectAny(("nhist_zwj_nnloewk_dw"+ext+"_"+obs).c_str()));
	zwjcornnloewkdwhist_den.push_back( (TH1*)zwjcornnloewkdwfile->FindObjectAny(("dhist_zwj_nnloewk_dw"+ext+"_"+obs).c_str()));
	
	zwjcorsudakovuphist_1.push_back( (TH1*)zwjcorsudakovupfile_1->FindObjectAny(("zwjcorsudakov_up_1"+ext+"hist_"+obs).c_str()));
	zwjcorsudakovuphist_1_num.push_back( (TH1*)zwjcorsudakovupfile_1->FindObjectAny(("nhist_zwj_sudakov_up_1"+ext+"_"+obs).c_str()));
	zwjcorsudakovuphist_1_den.push_back( (TH1*)zwjcorsudakovupfile_1->FindObjectAny(("dhist_zwj_sudakov_up_1"+ext+"_"+obs).c_str()));

	zwjcorsudakovdwhist_1.push_back( (TH1*)zwjcorsudakovdwfile_1->FindObjectAny(("zwjcorsudakov_dw_1"+ext+"hist_"+obs).c_str()));
	zwjcorsudakovdwhist_1_num.push_back( (TH1*)zwjcorsudakovdwfile_1->FindObjectAny(("nhist_zwj_sudakov_dw_1"+ext+"_"+obs).c_str()));
	zwjcorsudakovdwhist_1_den.push_back( (TH1*)zwjcorsudakovdwfile_1->FindObjectAny(("dhist_zwj_sudakov_dw_1"+ext+"_"+obs).c_str()));
	
	zwjcorsudakovuphist_2.push_back( (TH1*)zwjcorsudakovupfile_2->FindObjectAny(("zwjcorsudakov_up_2"+ext+"hist_"+obs).c_str()));
	zwjcorsudakovuphist_2_num.push_back( (TH1*)zwjcorsudakovupfile_2->FindObjectAny(("nhist_zwj_sudakov_up_2"+ext+"_"+obs).c_str()));
	zwjcorsudakovuphist_2_den.push_back( (TH1*)zwjcorsudakovupfile_2->FindObjectAny(("dhist_zwj_sudakov_up_2"+ext+"_"+obs).c_str()));
	
	zwjcorsudakovdwhist_2.push_back( (TH1*)zwjcorsudakovdwfile_2->FindObjectAny(("zwjcorsudakov_dw_2"+ext+"hist_"+obs).c_str()));
	zwjcorsudakovdwhist_2_num.push_back( (TH1*)zwjcorsudakovdwfile_2->FindObjectAny(("nhist_zwj_sudakov_dw_2"+ext+"_"+obs).c_str()));
	zwjcorsudakovdwhist_2_den.push_back( (TH1*)zwjcorsudakovdwfile_2->FindObjectAny(("dhist_zwj_sudakov_dw_2"+ext+"_"+obs).c_str()));
	
	zwjcornnlomissuphist_1.push_back( (TH1*)zwjcornnlomissupfile_1->FindObjectAny(("zwjcornnlomiss_up_1"+ext+"hist_"+obs).c_str()));
	zwjcornnlomissuphist_1_num.push_back( (TH1*)zwjcornnlomissupfile_1->FindObjectAny(("nhist_zwj_nnlomiss_up_1"+ext+"_"+obs).c_str()));
	zwjcornnlomissuphist_1_den.push_back( (TH1*)zwjcornnlomissupfile_1->FindObjectAny(("dhist_zwj_nnlomiss_up_1"+ext+"_"+obs).c_str()));
	
	zwjcornnlomissdwhist_1.push_back( (TH1*)zwjcornnlomissdwfile_1->FindObjectAny(("zwjcornnlomiss_dw_1"+ext+"hist_"+obs).c_str()));
	zwjcornnlomissdwhist_1_num.push_back( (TH1*)zwjcornnlomissdwfile_1->FindObjectAny(("nhist_zwj_nnlomiss_dw_1"+ext+"_"+obs).c_str()));
	zwjcornnlomissdwhist_1_den.push_back( (TH1*)zwjcornnlomissdwfile_1->FindObjectAny(("dhist_zwj_nnlomiss_dw_1"+ext+"_"+obs).c_str()));
	
	zwjcornnlomissuphist_2.push_back( (TH1*)zwjcornnlomissupfile_2->FindObjectAny(("zwjcornnlomiss_up_2"+ext+"hist_"+obs).c_str()));
	zwjcornnlomissuphist_2_num.push_back( (TH1*)zwjcornnlomissupfile_2->FindObjectAny(("nhist_zwj_nnlomiss_up_2"+ext+"_"+obs).c_str()));
	zwjcornnlomissuphist_2_den.push_back( (TH1*)zwjcornnlomissupfile_2->FindObjectAny(("dhist_zwj_nnlomiss_up_2"+ext+"_"+obs).c_str()));
	
	zwjcornnlomissdwhist_2.push_back( (TH1*)zwjcornnlomissdwfile_2->FindObjectAny(("zwjcornnlomiss_dw_2"+ext+"hist_"+obs).c_str()));
	zwjcornnlomissdwhist_2_num.push_back( (TH1*)zwjcornnlomissdwfile_2->FindObjectAny(("nhist_zwj_nnlomiss_dw_2"+ext+"_"+obs).c_str()));
	zwjcornnlomissdwhist_2_den.push_back( (TH1*)zwjcornnlomissdwfile_2->FindObjectAny(("dhist_zwj_nnlomiss_dw_2"+ext+"_"+obs).c_str()));
	
	zwjcormixuphist.push_back( (TH1*)zwjcormixupfile->FindObjectAny(("zwjcormix_up"+ext+"hist_"+obs).c_str()));
	zwjcormixuphist_num.push_back( (TH1*)zwjcormixupfile->FindObjectAny(("nhist_zwj_mix_up"+ext+"_"+obs).c_str()));
	zwjcormixuphist_den.push_back( (TH1*)zwjcormixupfile->FindObjectAny(("dhist_zwj_mix_up"+ext+"_"+obs).c_str()));
	
	zwjcormixdwhist.push_back( (TH1*)zwjcormixdwfile->FindObjectAny(("zwjcormix_dw"+ext+"hist_"+obs).c_str()));
	zwjcormixdwhist_num.push_back( (TH1*)zwjcormixdwfile->FindObjectAny(("nhist_zwj_mix_dw"+ext+"_"+obs).c_str()));
	zwjcormixdwhist_den.push_back( (TH1*)zwjcormixdwfile->FindObjectAny(("dhist_zwj_mix_dw"+ext+"_"+obs).c_str()));
	
	zwjcorpdfhist.push_back( (TH1*)zwjcorpdffile->FindObjectAny(("zwjcorpdf"+ext+"hist_"+obs).c_str()));
	zwjcorpdfhist_num.push_back( (TH1*)zwjcorpdffile->FindObjectAny(("nhist_zwj_pdf"+ext+"_"+obs).c_str()));
	zwjcorpdfhist_den.push_back( (TH1*)zwjcorpdffile->FindObjectAny(("dhist_zwj_pdf"+ext+"_"+obs).c_str()));
	
	// QCD scale --> to be symmetrized
	zwjuncqcdscalehist = (TH1*) zwjcorqcdscaleuphist.back()->Clone(("zwjuncqcdscale"+ext+"hist_"+obs).c_str());
	zwjuncqcdscalehist->Reset("ICES");
	for (int i = 0; i <= zwjuncqcdscalehist->GetNbinsX()+1; i++)
	  zwjuncqcdscalehist->SetBinContent(i,(zwjcorqcdscaleuphist.back()->GetBinContent(i)-zwjcorqcdscaledwhist.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
	zwjuncqcdscalehist->Smooth(1,"R");
	bool sign = hasSameSign(zwjuncqcdscalehist);
	if(sign){
	  for (int i = 0; i <= zwjuncqcdscalehist->GetNbinsX()+1; i++){
	    zwjuncqcdscalehist->SetBinContent(i,fabs(zwjuncqcdscalehist->GetBinContent(i)));
	  }
	}
	zwjuncqcdscalehist->SetName(("ZW_QCDScale_"+obs).c_str());
	
	// QCD Shape --> to be symmetrized
	zwjuncqcdshapehist = (TH1*) zwjcorqcdshapeuphist.back()->Clone(("zwjuncqcdshape"+ext+"hist_"+obs).c_str());
	zwjuncqcdshapehist->Reset("ICES");
	for (int i = 0; i <= zwjuncqcdshapehist->GetNbinsX()+1; i++)
	  zwjuncqcdshapehist->SetBinContent(i,(zwjcorqcdshapeuphist.back()->GetBinContent(i)-zwjcorqcdshapedwhist.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
	zwjuncqcdshapehist->Smooth(1,"R");	
	sign = hasSameSign(zwjuncqcdshapehist);
	if(sign){
	  for (int i = 0; i <= zwjuncqcdshapehist->GetNbinsX()+1; i++)
	    zwjuncqcdshapehist->SetBinContent(i,fabs(zwjuncqcdshapehist->GetBinContent(i)));
	}
	zwjuncqcdshapehist->SetName(("ZW_QCDShape_"+obs).c_str());
	
	// QCD Process --> to be symmetrized
	zwjuncqcdprochist = (TH1*) zwjcorqcdprocuphist.back()->Clone(("zwjuncqcdproc"+ext+"hist_"+obs).c_str());
	zwjuncqcdprochist->Reset("ICES");
	for (int i = 0; i <= zwjuncqcdprochist->GetNbinsX()+1; i++)
	  zwjuncqcdprochist->SetBinContent(i,(zwjcorqcdprocuphist.back()->GetBinContent(i)-zwjcorqcdprocdwhist.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
	zwjuncqcdprochist->Smooth(1,"R");	
	sign = hasSameSign(zwjuncqcdprochist);
	if(sign){
	  for (int i = 0; i <= zwjuncqcdprochist->GetNbinsX()+1; i++)
	    zwjuncqcdprochist->SetBinContent(i,fabs(zwjuncqcdprochist->GetBinContent(i)));
	}
	zwjuncqcdprochist->SetName(("ZW_QCDProcess_"+obs).c_str());
	
	// NNLO EWK --> to be symmetrized
	zwjuncnnloewkhist = (TH1*) zwjcornnloewkuphist.back()->Clone(("zwjuncnnloewk"+ext+"hist_"+obs).c_str());
	zwjuncnnloewkhist->Reset("ICES");
	for (int i = 0; i <= zwjuncnnloewkhist->GetNbinsX()+1; i++){
	  zwjuncnnloewkhist->SetBinContent(i,fabs(zwjcornnloewkuphist.back()->GetBinContent(i)-zwjcornnloewkdwhist.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
	}
	zwjuncnnloewkhist->Smooth(1,"R");
	zwjuncnnloewkhist->SetName(("ZW_NNLOEWK_"+obs).c_str());
	
	// NLO EWK Sud --> to be symmetrized
	zwjuncsudakovhist_1 = (TH1*) zwjcorsudakovuphist_1.back()->Clone(("zwjuncsudakov1"+ext+"hist_"+obs).c_str());
	zwjuncsudakovhist_1->Reset("ICES");
	for (int i = 0; i <= zwjuncsudakovhist_1->GetNbinsX()+1; i++){
	  zwjuncsudakovhist_1->SetBinContent(i,fabs(zwjcorsudakovuphist_1.back()->GetBinContent(i)-zwjcorsudakovdwhist_1.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
	}
	zwjuncsudakovhist_1->Smooth(1,"R");
	zwjuncsudakovhist_1->SetName(("ZW_Sudakov1_"+obs).c_str());
	
	// NLO EWK Sud --> to be symmetrized
	zwjuncsudakovhist_2 = (TH1*) zwjcorsudakovuphist_2.back()->Clone(("zwjuncsudakov2"+ext+"hist_"+obs).c_str());
	zwjuncsudakovhist_2->Reset("ICES");
	for (int i = 0; i <= zwjuncsudakovhist_2->GetNbinsX()+1; i++){
	  zwjuncsudakovhist_2->SetBinContent(i,fabs(zwjcorsudakovuphist_2.back()->GetBinContent(i)-zwjcorsudakovdwhist_2.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
	}
	zwjuncsudakovhist_2->Smooth(1,"R");
	zwjuncsudakovhist_2->SetName(("ZW_Sudakov2_"+obs).c_str());
	
	// NNLO Miss --> to be symmetrized
	zwjuncnnlomisshist_1 = (TH1*) zwjcornnlomissuphist_1.back()->Clone(("zwjuncnnlomiss1"+ext+"hist_"+obs).c_str());
	zwjuncnnlomisshist_1->Reset("ICES");
	for (int i = 0; i <= zwjuncnnlomisshist_1->GetNbinsX()+1; i++){
	  zwjuncnnlomisshist_1->SetBinContent(i,fabs(zwjcornnlomissuphist_1.back()->GetBinContent(i)-zwjcornnlomissdwhist_1.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
	}
	zwjuncnnlomisshist_1->Smooth(1,"R");
	zwjuncnnlomisshist_1->SetName(("ZW_NNLOMiss1_"+obs).c_str());
	
	// NNLO Miss --> to be symmetrized
	zwjuncnnlomisshist_2 = (TH1*) zwjcornnlomissuphist_2.back()->Clone(("zwjuncnnlomiss2"+ext+"hist_"+obs).c_str());
	zwjuncnnlomisshist_2->Reset("ICES");
	for (int i = 0; i <= zwjuncnnlomisshist_2->GetNbinsX()+1; i++){
	  zwjuncnnlomisshist_2->SetBinContent(i,fabs(zwjcornnlomissuphist_2.back()->GetBinContent(i)-zwjcornnlomissdwhist_2.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
	}
	zwjuncnnlomisshist_2->Smooth(1,"R");
	zwjuncnnlomisshist_2->SetName(("ZW_NNLOMiss2_"+obs).c_str());
	
	// NLO QCD+EWK --> to be symmetrized
	zwjuncmixhist = (TH1*) zwjcormixuphist.back()->Clone(("zwjuncmix"+ext+"hist_"+obs).c_str());
	zwjuncmixhist->Reset("ICES");
	for (int i = 0; i <= zwjuncmixhist->GetNbinsX()+1; i++){
	  zwjuncmixhist->SetBinContent(i,fabs(zwjcormixuphist.back()->GetBinContent(i)-zwjcormixdwhist.back()->GetBinContent(i))/(2*zwjcorewkhist.back()->GetBinContent(i)));
	}
	zwjuncmixhist->Smooth(1,"R");
	zwjuncmixhist->SetName(("ZW_MIX_"+obs).c_str());
      
	// NLO PDF
	zwjuncpdfhist = (TH1*)zwjcorpdfhist.back()->Clone(("zwjuncpdf"+ext+"hist_"+obs).c_str());
	zwjuncpdfhist->Divide(zwjcorqcdhist.back());
	for (int i = 1; i <= zwjuncpdfhist->GetNbinsX(); i++)
	  zwjuncpdfhist->SetBinContent(i,fabs(zwjuncpdfhist->GetBinContent(i)-1.0));
	zwjuncpdfhist->SetName(("ZW_PDF_"+obs).c_str());
      }
    }

    /////////////

    // get histograms W/gamma               
    TH1* wgamuncewkhist = NULL;
    TH1* wgamuncre1hist = NULL;
    TH1* wgamuncfa1hist = NULL;
    TH1* wgamuncre2hist = NULL;
    TH1* wgamuncfa2hist = NULL;
    TH1* wgamuncpdfhist = NULL;
    TH1* wgamuncfpchist = NULL;
    TH1* wgamuncqcdscalehist = NULL;
    TH1* wgamuncqcdshapehist = NULL;
    TH1* wgamuncqcdprochist = NULL;
    TH1* wgamuncnnloewkhist = NULL;
    TH1* wgamuncsudakovhist_1 = NULL;
    TH1* wgamuncsudakovhist_2 = NULL;
    TH1* wgamuncnnlomisshist_1 = NULL;
    TH1* wgamuncnnlomisshist_2 = NULL;
    TH1* wgamuncmixhist = NULL;

    //////
    if(addWgamma){

      if(not useNewTheoryUncertainty){
	cout<<"Make W/gamma sys histograms"<<endl;    
	wgamcorewkhist.push_back( (TH1*)wgamcorewkfile->FindObjectAny(("wgamcor"+ext+"ewkhist_"+obs).c_str()));
	wgamcorewkhist_num.push_back( (TH1*)wgamcorewkfile->FindObjectAny(("nhist_wgam_ewk"+ext+"_"+obs).c_str()));
	wgamcorewkhist_den.push_back( (TH1*)wgamcorewkfile->FindObjectAny(("dhist_wgam_ewk"+ext+"_"+obs).c_str()));
	
	wgamcorqcdhist.push_back( (TH1*)wgamcorqcdfile->FindObjectAny(("wgamcor"+ext+"qcdhist_"+obs).c_str()));
	wgamcorqcdhist_num.push_back( (TH1*)wgamcorqcdfile->FindObjectAny(("nhist_wgam_qcd"+ext+"_"+obs).c_str()));
	wgamcorqcdhist_den.push_back( (TH1*)wgamcorqcdfile->FindObjectAny(("dhist_wgam_qcd"+ext+"_"+obs).c_str()));
	
	wgamcorre1hist.push_back( (TH1*)wgamcorre1file->FindObjectAny(("wgamcor"+ext+"re1hist_"+obs).c_str()));
	wgamcorre1hist_num.push_back( (TH1*)wgamcorre1file->FindObjectAny(("nhist_wgam_re1"+ext+"_"+obs).c_str()));
	wgamcorre1hist_den.push_back( (TH1*)wgamcorre1file->FindObjectAny(("dhist_wgam_re1"+ext+"_"+obs).c_str()));
	
	wgamcorfa1hist.push_back( (TH1*)wgamcorfa1file->FindObjectAny(("wgamcor"+ext+"fa1hist_"+obs).c_str()));
	wgamcorfa1hist_num.push_back( (TH1*)wgamcorfa1file->FindObjectAny(("nhist_wgam_fa1"+ext+"_"+obs).c_str()));
	wgamcorfa1hist_den.push_back( (TH1*)wgamcorfa1file->FindObjectAny(("dhist_wgam_fa1"+ext+"_"+obs).c_str()));
	
	wgamcorre2hist.push_back( (TH1*)wgamcorre2file->FindObjectAny(("wgamcor"+ext+"re2hist_"+obs).c_str()));
	wgamcorre2hist_num.push_back( (TH1*)wgamcorre2file->FindObjectAny(("nhist_wgam_re2"+ext+"_"+obs).c_str()));
	wgamcorre2hist_den.push_back( (TH1*)wgamcorre2file->FindObjectAny(("dhist_wgam_re2"+ext+"_"+obs).c_str()));
	
	wgamcorfa2hist.push_back( (TH1*)wgamcorfa2file->FindObjectAny(("wgamcor"+ext+"fa2hist_"+obs).c_str()));
	wgamcorfa2hist_num.push_back( (TH1*)wgamcorfa2file->FindObjectAny(("nhist_wgam_fa2"+ext+"_"+obs).c_str()));
	wgamcorfa2hist_den.push_back( (TH1*)wgamcorfa2file->FindObjectAny(("dhist_wgam_fa2"+ext+"_"+obs).c_str()));
	
	wgamcorpdfhist.push_back( (TH1*)wgamcorpdffile->FindObjectAny(("wgamcor"+ext+"pdfhist_"+obs).c_str()));
	wgamcorpdfhist_num.push_back( (TH1*)wgamcorpdffile->FindObjectAny(("nhist_wgam_pdf"+ext+"_"+obs).c_str()));
	wgamcorpdfhist_den.push_back( (TH1*)wgamcorpdffile->FindObjectAny(("dhist_wgam_pdf"+ext+"_"+obs).c_str()));
	
	wgamcorfpchist.push_back( (TH1*)wgamcorfpcfile->FindObjectAny(("wgamcor"+ext+"fpchist_"+obs).c_str()));
	wgamcorfpchist_num.push_back( (TH1*)wgamcorfpcfile->FindObjectAny(("nhist_wgam_fpc"+ext+"_"+obs).c_str()));
	wgamcorfpchist_den.push_back( (TH1*)wgamcorfpcfile->FindObjectAny(("dhist_wgam_fpc"+ext+"_"+obs).c_str()));
	
	// ZG ewk uncertainty
	wgamuncewkhist = (TH1*)wgamcorewkhist.back()->Clone(("wgamuncewk"+ext+"hist_"+obs).c_str());
	wgamuncewkhist->Divide(wgamcorqcdhist.back());
	for (int i = 1; i <= wgamuncewkhist->GetNbinsX(); i++)
	  wgamuncewkhist->SetBinContent(i, fabs(wgamuncewkhist->GetBinContent(i)-1.0));
	wgamuncewkhist->SetName(("WG_EWK_"+obs).c_str());
	
	/// ZG renormalization scale                                                                                                                                                           
	wgamuncre1hist = (TH1*)wgamcorre1hist.back()->Clone(("wgamuncre1"+ext+"hist_"+obs).c_str());
	wgamuncre1hist->Divide(wgamcorqcdhist.back());
	for (int i = 1; i <= wgamuncre1hist->GetNbinsX(); i++)
	  wgamuncre1hist->SetBinContent(i, fabs(wgamuncre1hist->GetBinContent(i)-1.0));
	wgamuncre1hist->SetName(("WG_RenScale1_"+obs).c_str());
	
	/// ZG renormalization scale                                                                                                                                                           
	wgamuncfa1hist = (TH1*)wgamcorfa1hist.back()->Clone(("wgamuncfa1"+ext+"hist_"+obs).c_str());
	wgamuncfa1hist->Divide(wgamcorqcdhist.back());
	for (int i = 1; i <= wgamuncfa1hist->GetNbinsX(); i++)
	  wgamuncfa1hist->SetBinContent(i, fabs(wgamuncfa1hist->GetBinContent(i)-1.0));
	wgamuncfa1hist->SetName(("WG_FactScale1_"+obs).c_str());
	
	/// ZG renormalization scale                                                                                                                                                           
	wgamuncre2hist = (TH1*)wgamcorre2hist.back()->Clone(("wgamuncre2"+ext+"hist_"+obs).c_str());
	wgamuncre2hist->Divide(wgamcorqcdhist.back());
	for (int i = 1; i <= wgamuncre2hist->GetNbinsX(); i++)
	  wgamuncre2hist->SetBinContent(i, fabs(wgamuncre2hist->GetBinContent(i)-1.0));
	wgamuncre2hist->SetName(("WG_RenScale2_"+obs).c_str());
      
	/// ZG factorization scale                                                                                                                                                                 
	wgamuncfa2hist = (TH1*)wgamcorfa2hist.back()->Clone(("wgamuncfa2"+ext+"hist_"+obs).c_str());
	wgamuncfa2hist->Divide(wgamcorqcdhist.back());
	for (int i = 1; i <= wgamuncfa2hist->GetNbinsX(); i++)
	  wgamuncfa2hist->SetBinContent(i, fabs(wgamuncfa2hist->GetBinContent(i)-1.0));
	wgamuncfa2hist->SetName(("WG_FactScale2_"+obs).c_str());
	
	/// ZG pdfs                                                                                                                                                                 
	wgamuncpdfhist = (TH1*)wgamcorpdfhist.back()->Clone(("wgamuncpdf"+ext+"hist_"+obs).c_str());
	wgamuncpdfhist->Divide(wgamcorqcdhist.back());
	for (int i = 1; i <= wgamuncpdfhist->GetNbinsX(); i++)
	  wgamuncpdfhist->SetBinContent(i, fabs(wgamuncpdfhist->GetBinContent(i)-1.0));
	wgamuncpdfhist->SetName(("WG_PDF_"+obs).c_str());
      
	/// ZG footprint                                                                                                                                                                            
	wgamuncfpchist = (TH1*)wgamcorfpchist.back()->Clone(("wgamuncfpc"+ext+"hist_"+obs).c_str());
	wgamuncfpchist->Divide(wgamcorqcdhist.back());
	for (int i = 1; i <= wgamuncfpchist->GetNbinsX(); i++)
	  wgamuncfpchist->SetBinContent(i, fabs(wgamuncfpchist->GetBinContent(i)-1.0));
	wgamuncfpchist->SetName(("WG_Footprint_"+obs).c_str());
      }
      else{

	cout<<"Make W/gamma sys histograms"<<endl;    
	wgamcorewkhist.push_back( (TH1*)wgamcorewkfile->FindObjectAny(("wgamcor"+ext+"ewkhist_"+obs).c_str()));
	wgamcorewkhist_num.push_back( (TH1*)wgamcorewkfile->FindObjectAny(("nhist_wgam_ewk"+ext+"_"+obs).c_str()));
	wgamcorewkhist_den.push_back( (TH1*)wgamcorewkfile->FindObjectAny(("dhist_wgam_ewk"+ext+"_"+obs).c_str()));

	wgamcorqcdhist.push_back( (TH1*)wgamcorqcdfile->FindObjectAny(("wgamcor"+ext+"qcdhist_"+obs).c_str()));
	wgamcorqcdhist_num.push_back( (TH1*)wgamcorqcdfile->FindObjectAny(("nhist_wgam_qcd"+ext+"_"+obs).c_str()));
	wgamcorqcdhist_den.push_back( (TH1*)wgamcorqcdfile->FindObjectAny(("dhist_wgam_qcd"+ext+"_"+obs).c_str()));
	
	wgamcorpdfhist.push_back( (TH1*)wgamcorpdffile->FindObjectAny(("wgamcor"+ext+"pdfhist_"+obs).c_str()));
        wgamcorpdfhist_num.push_back( (TH1*)wgamcorpdffile->FindObjectAny(("nhist_wgam_pdf"+ext+"_"+obs).c_str()));
        wgamcorpdfhist_den.push_back( (TH1*)wgamcorpdffile->FindObjectAny(("dhist_wgam_pdf"+ext+"_"+obs).c_str()));

	wgamcorqcdscaleuphist.push_back( (TH1*)wgamcorqcdscaleupfile->FindObjectAny(("wgamcorqcdscale_up"+ext+"hist_"+obs).c_str()));
	wgamcorqcdscaleuphist_num.push_back( (TH1*)wgamcorqcdscaleupfile->FindObjectAny(("nhist_wgam_qcdscale_up"+ext+"_"+obs).c_str()));
	wgamcorqcdscaleuphist_den.push_back( (TH1*)wgamcorqcdscaleupfile->FindObjectAny(("dhist_wgam_qcdscale_up"+ext+"_"+obs).c_str()));

	wgamcorqcdscaledwhist.push_back( (TH1*)wgamcorqcdscaledwfile->FindObjectAny(("wgamcorqcdscale_dw"+ext+"hist_"+obs).c_str()));
	wgamcorqcdscaledwhist_num.push_back( (TH1*)wgamcorqcdscaledwfile->FindObjectAny(("nhist_wgam_qcdscale_dw"+ext+"_"+obs).c_str()));
	wgamcorqcdscaledwhist_den.push_back( (TH1*)wgamcorqcdscaledwfile->FindObjectAny(("dhist_wgam_qcdscale_dw"+ext+"_"+obs).c_str()));

	wgamcorqcdshapeuphist.push_back( (TH1*)wgamcorqcdshapeupfile->FindObjectAny(("wgamcorqcdshape_up"+ext+"hist_"+obs).c_str()));
	wgamcorqcdshapeuphist_num.push_back( (TH1*)wgamcorqcdshapeupfile->FindObjectAny(("nhist_wgam_qcdshape_up"+ext+"_"+obs).c_str()));
	wgamcorqcdshapeuphist_den.push_back( (TH1*)wgamcorqcdshapeupfile->FindObjectAny(("dhist_wgam_qcdshape_up"+ext+"_"+obs).c_str()));

	wgamcorqcdshapedwhist.push_back( (TH1*)wgamcorqcdshapedwfile->FindObjectAny(("wgamcorqcdshape_dw"+ext+"hist_"+obs).c_str()));
	wgamcorqcdshapedwhist_num.push_back( (TH1*)wgamcorqcdshapedwfile->FindObjectAny(("nhist_wgam_qcdshape_dw"+ext+"_"+obs).c_str()));
	wgamcorqcdshapedwhist_den.push_back( (TH1*)wgamcorqcdshapedwfile->FindObjectAny(("dhist_wgam_qcdshape_dw"+ext+"_"+obs).c_str()));

	wgamcorqcdprocuphist.push_back( (TH1*)wgamcorqcdprocupfile->FindObjectAny(("wgamcorqcdproc_up"+ext+"hist_"+obs).c_str()));
	wgamcorqcdprocuphist_num.push_back( (TH1*)wgamcorqcdprocupfile->FindObjectAny(("nhist_wgam_qcdproc_up"+ext+"_"+obs).c_str()));
	wgamcorqcdprocuphist_den.push_back( (TH1*)wgamcorqcdprocupfile->FindObjectAny(("dhist_wgam_qcdproc_up"+ext+"_"+obs).c_str()));

	wgamcorqcdprocdwhist.push_back( (TH1*)wgamcorqcdprocdwfile->FindObjectAny(("wgamcorqcdproc_dw"+ext+"hist_"+obs).c_str()));
	wgamcorqcdprocdwhist_num.push_back( (TH1*)wgamcorqcdprocdwfile->FindObjectAny(("nhist_wgam_qcdproc_dw"+ext+"_"+obs).c_str()));
	wgamcorqcdprocdwhist_den.push_back( (TH1*)wgamcorqcdprocdwfile->FindObjectAny(("dhist_wgam_qcdproc_dw"+ext+"_"+obs).c_str()));

	wgamcornnloewkuphist.push_back( (TH1*)wgamcornnloewkupfile->FindObjectAny(("wgamcornnloewk_up"+ext+"hist_"+obs).c_str()));
	wgamcornnloewkuphist_num.push_back( (TH1*)wgamcornnloewkupfile->FindObjectAny(("nhist_wgam_nnloewk_up"+ext+"_"+obs).c_str()));
	wgamcornnloewkuphist_den.push_back( (TH1*)wgamcornnloewkupfile->FindObjectAny(("dhist_wgam_nnloewk_up"+ext+"_"+obs).c_str()));

	wgamcornnloewkdwhist.push_back( (TH1*)wgamcornnloewkdwfile->FindObjectAny(("wgamcornnloewk_dw"+ext+"hist_"+obs).c_str()));
	wgamcornnloewkdwhist_num.push_back( (TH1*)wgamcornnloewkdwfile->FindObjectAny(("nhist_wgam_nnloewk_dw"+ext+"_"+obs).c_str()));
	wgamcornnloewkdwhist_den.push_back( (TH1*)wgamcornnloewkdwfile->FindObjectAny(("dhist_wgam_nnloewk_dw"+ext+"_"+obs).c_str()));

	wgamcorsudakovuphist_1.push_back( (TH1*)wgamcorsudakovupfile_1->FindObjectAny(("wgamcorsudakov_up_1"+ext+"hist_"+obs).c_str()));
	wgamcorsudakovuphist_1_num.push_back( (TH1*)wgamcorsudakovupfile_1->FindObjectAny(("nhist_wgam_sudakov_up_1"+ext+"_"+obs).c_str()));
	wgamcorsudakovuphist_1_den.push_back( (TH1*)wgamcorsudakovupfile_1->FindObjectAny(("dhist_wgam_sudakov_up_1"+ext+"_"+obs).c_str()));
	
	wgamcorsudakovdwhist_1.push_back( (TH1*)wgamcorsudakovdwfile_1->FindObjectAny(("wgamcorsudakov_dw_1"+ext+"hist_"+obs).c_str()));
	wgamcorsudakovdwhist_1_num.push_back( (TH1*)wgamcorsudakovdwfile_1->FindObjectAny(("nhist_wgam_sudakov_dw_1"+ext+"_"+obs).c_str()));
	wgamcorsudakovdwhist_1_den.push_back( (TH1*)wgamcorsudakovdwfile_1->FindObjectAny(("dhist_wgam_sudakov_dw_1"+ext+"_"+obs).c_str()));
	
	wgamcorsudakovuphist_2.push_back( (TH1*)wgamcorsudakovupfile_2->FindObjectAny(("wgamcorsudakov_up_2"+ext+"hist_"+obs).c_str()));
	wgamcorsudakovuphist_2_num.push_back( (TH1*)wgamcorsudakovupfile_2->FindObjectAny(("nhist_wgam_sudakov_up_2"+ext+"_"+obs).c_str()));
	wgamcorsudakovuphist_2_den.push_back( (TH1*)wgamcorsudakovupfile_2->FindObjectAny(("dhist_wgam_sudakov_up_2"+ext+"_"+obs).c_str()));

	wgamcorsudakovdwhist_2.push_back( (TH1*)wgamcorsudakovdwfile_2->FindObjectAny(("wgamcorsudakov_dw_2"+ext+"hist_"+obs).c_str()));
	wgamcorsudakovdwhist_2_num.push_back( (TH1*)wgamcorsudakovdwfile_2->FindObjectAny(("nhist_wgam_sudakov_dw_2"+ext+"_"+obs).c_str()));
	wgamcorsudakovdwhist_2_den.push_back( (TH1*)wgamcorsudakovdwfile_2->FindObjectAny(("dhist_wgam_sudakov_dw_2"+ext+"_"+obs).c_str()));
	
	wgamcornnlomissuphist_1.push_back( (TH1*)wgamcornnlomissupfile_1->FindObjectAny(("wgamcornnlomiss_up_1"+ext+"hist_"+obs).c_str()));
	wgamcornnlomissuphist_1_num.push_back( (TH1*)wgamcornnlomissupfile_1->FindObjectAny(("nhist_wgam_nnlomiss_up_1"+ext+"_"+obs).c_str()));
	wgamcornnlomissuphist_1_den.push_back( (TH1*)wgamcornnlomissupfile_1->FindObjectAny(("dhist_wgam_nnlomiss_up_1"+ext+"_"+obs).c_str()));
	
	wgamcornnlomissdwhist_1.push_back( (TH1*)wgamcornnlomissdwfile_1->FindObjectAny(("wgamcornnlomiss_dw_1"+ext+"hist_"+obs).c_str()));
	wgamcornnlomissdwhist_1_num.push_back( (TH1*)wgamcornnlomissdwfile_1->FindObjectAny(("nhist_wgam_nnlomiss_dw_1"+ext+"_"+obs).c_str()));
	wgamcornnlomissdwhist_1_den.push_back( (TH1*)wgamcornnlomissdwfile_1->FindObjectAny(("dhist_wgam_nnlomiss_dw_1"+ext+"_"+obs).c_str()));
	
	wgamcornnlomissuphist_2.push_back( (TH1*)wgamcornnlomissupfile_2->FindObjectAny(("wgamcornnlomiss_up_2"+ext+"hist_"+obs).c_str()));
	wgamcornnlomissuphist_2_num.push_back( (TH1*)wgamcornnlomissupfile_2->FindObjectAny(("nhist_wgam_nnlomiss_up_2"+ext+"_"+obs).c_str()));
	wgamcornnlomissuphist_2_den.push_back( (TH1*)wgamcornnlomissupfile_2->FindObjectAny(("dhist_wgam_nnlomiss_up_2"+ext+"_"+obs).c_str()));
	
	wgamcornnlomissdwhist_2.push_back( (TH1*)wgamcornnlomissdwfile_2->FindObjectAny(("wgamcornnlomiss_dw_2"+ext+"hist_"+obs).c_str()));
	wgamcornnlomissdwhist_2_num.push_back( (TH1*)wgamcornnlomissdwfile_2->FindObjectAny(("nhist_wgam_nnlomiss_dw_2"+ext+"_"+obs).c_str()));
	wgamcornnlomissdwhist_2_den.push_back( (TH1*)wgamcornnlomissdwfile_2->FindObjectAny(("dhist_wgam_nnlomiss_dw_2"+ext+"_"+obs).c_str()));
	
	wgamcormixuphist.push_back( (TH1*)wgamcormixupfile->FindObjectAny(("wgamcormix_up"+ext+"hist_"+obs).c_str()));
	wgamcormixuphist_num.push_back( (TH1*)wgamcormixupfile->FindObjectAny(("nhist_wgam_mix_up"+ext+"_"+obs).c_str()));
	wgamcormixuphist_den.push_back( (TH1*)wgamcormixupfile->FindObjectAny(("dhist_wgam_mix_up"+ext+"_"+obs).c_str()));

	wgamcormixdwhist.push_back( (TH1*)wgamcormixdwfile->FindObjectAny(("wgamcormix_dw"+ext+"hist_"+obs).c_str()));
	wgamcormixdwhist_num.push_back( (TH1*)wgamcormixdwfile->FindObjectAny(("nhist_wgam_mix_dw"+ext+"_"+obs).c_str()));
	wgamcormixdwhist_den.push_back( (TH1*)wgamcormixdwfile->FindObjectAny(("dhist_wgam_mix_dw"+ext+"_"+obs).c_str()));

	// QCD scale --> to be symmetrized
	wgamuncqcdscalehist = (TH1*) wgamcorqcdscaleuphist.back()->Clone(("wgamuncqcdscale"+ext+"hist_"+obs).c_str());	
	wgamuncqcdscalehist->Reset("ICES");
	for (int i = 0; i <= wgamuncqcdscalehist->GetNbinsX()+1; i++)
	  wgamuncqcdscalehist->SetBinContent(i,(wgamcorqcdscaleuphist.back()->GetBinContent(i)-wgamcorqcdscaledwhist.back()->GetBinContent(i))/(2*wgamcorewkhist.back()->GetBinContent(i)));	
	wgamuncqcdscalehist->Smooth(1,"R");
	bool sign = hasSameSign(wgamuncqcdscalehist);
	if(sign){
	  for (int i = 0; i <= wgamuncqcdscalehist->GetNbinsX()+1; i++)
	    wgamuncqcdscalehist->SetBinContent(i,fabs(wgamuncqcdscalehist->GetBinContent(i)));
	}
	wgamuncqcdscalehist->SetName(("WG_QCDScale_"+obs).c_str());
	
	// QCD Shape --> to be symmetrized
	wgamuncqcdshapehist = (TH1*) wgamcorqcdshapeuphist.back()->Clone(("wgamuncqcdshape"+ext+"hist_"+obs).c_str());
	wgamuncqcdshapehist->Reset("ICES");
	for (int i = 0; i <= wgamuncqcdshapehist->GetNbinsX()+1; i++)
	    wgamuncqcdshapehist->SetBinContent(i,(wgamcorqcdshapeuphist.back()->GetBinContent(i)-wgamcorqcdshapedwhist.back()->GetBinContent(i))/(2*wgamcorewkhist.back()->GetBinContent(i)));
	wgamuncqcdshapehist->Smooth(1,"R");	
	sign = hasSameSign(wgamuncqcdshapehist);
	if(sign){
	  for (int i = 0; i <= wgamuncqcdshapehist->GetNbinsX()+1; i++)
	    wgamuncqcdshapehist->SetBinContent(i,fabs(wgamuncqcdshapehist->GetBinContent(i)));
	}
	wgamuncqcdshapehist->SetName(("WG_QCDShape_"+obs).c_str());

	// QCD Process --> to be symmetrized
	wgamuncqcdprochist = (TH1*) wgamcorqcdprocuphist.back()->Clone(("wgamuncqcdproc"+ext+"hist_"+obs).c_str());
	wgamuncqcdprochist->Reset("ICES");
	for (int i = 0; i <= wgamuncqcdprochist->GetNbinsX()+1; i++)
	  wgamuncqcdprochist->SetBinContent(i,(wgamcorqcdprocuphist.back()->GetBinContent(i)-wgamcorqcdprocdwhist.back()->GetBinContent(i))/(2*wgamcorewkhist.back()->GetBinContent(i)));
	wgamuncqcdprochist->Smooth(1,"R");
	sign = hasSameSign(wgamuncqcdprochist);
	if(sign){
	  for (int i = 0; i <= wgamuncqcdprochist->GetNbinsX()+1; i++)
	    wgamuncqcdprochist->SetBinContent(i,fabs(wgamuncqcdprochist->GetBinContent(i)));
	}
	wgamuncqcdprochist->SetName(("WG_QCDProcess_"+obs).c_str());
	
	// NNLO EWK --> to be symmetrized
	wgamuncnnloewkhist = (TH1*) wgamcornnloewkuphist.back()->Clone(("wgamuncnnloewk"+ext+"hist_"+obs).c_str());
	wgamuncnnloewkhist->Reset("ICES");
	for (int i = 0; i <= wgamuncnnloewkhist->GetNbinsX()+1; i++)
	    wgamuncnnloewkhist->SetBinContent(i,fabs(wgamcornnloewkuphist.back()->GetBinContent(i)-wgamcornnloewkdwhist.back()->GetBinContent(i))/(2*wgamcorewkhist.back()->GetBinContent(i)));
	wgamuncnnloewkhist->Smooth(1,"R");
	wgamuncnnloewkhist->SetName(("WG_NNLOEWK_"+obs).c_str());
	
	// NLO EWK Sud --> to be symmetrized
	wgamuncsudakovhist_1 = (TH1*) wgamcorsudakovuphist_1.back()->Clone(("wgamuncsudakov1"+ext+"hist_"+obs).c_str());
	wgamuncsudakovhist_1->Reset("ICES");
	for (int i = 0; i <= wgamuncsudakovhist_1->GetNbinsX()+1; i++)
	  wgamuncsudakovhist_1->SetBinContent(i,fabs(wgamcorsudakovuphist_1.back()->GetBinContent(i)-wgamcorsudakovdwhist_1.back()->GetBinContent(i))/(2*wgamcorewkhist.back()->GetBinContent(i)));
	wgamuncsudakovhist_1->Smooth(1,"R");
	wgamuncsudakovhist_1->SetName(("WG_Sudakov1_"+obs).c_str());

	// NLO EWK Sud --> to be symmetrized
	wgamuncsudakovhist_2 = (TH1*) wgamcorsudakovuphist_2.back()->Clone(("wgamuncsudakov2"+ext+"hist_"+obs).c_str());
	wgamuncsudakovhist_2->Reset("ICES");
	for (int i = 0; i <= wgamuncsudakovhist_2->GetNbinsX()+1; i++)
	  wgamuncsudakovhist_2->SetBinContent(i,fabs(wgamcorsudakovuphist_2.back()->GetBinContent(i)-wgamcorsudakovdwhist_2.back()->GetBinContent(i))/(2*wgamcorewkhist.back()->GetBinContent(i)));
	wgamuncsudakovhist_2->Smooth(1,"R");
      	wgamuncsudakovhist_2->SetName(("WG_Sudakov2_"+obs).c_str());

	// NNLO Miss --> to be symmetrized
	wgamuncnnlomisshist_1 = (TH1*) wgamcornnlomissuphist_1.back()->Clone(("wgamuncnnlomiss1"+ext+"hist_"+obs).c_str());
	wgamuncnnlomisshist_1->Reset("ICES");
	for (int i = 0; i <= wgamuncnnlomisshist_1->GetNbinsX()+1; i++)
	  wgamuncnnlomisshist_1->SetBinContent(i,fabs(wgamcornnlomissuphist_1.back()->GetBinContent(i)-wgamcornnlomissdwhist_1.back()->GetBinContent(i))/(2*wgamcorewkhist.back()->GetBinContent(i)));
	wgamuncnnlomisshist_1->Smooth(1,"R");
	wgamuncnnlomisshist_1->SetName(("WG_NNLOMiss1_"+obs).c_str());
	
	// NNLO Miss --> to be symmetrized
	wgamuncnnlomisshist_2 = (TH1*) wgamcornnlomissuphist_2.back()->Clone(("wgamuncnnlomiss2"+ext+"hist_"+obs).c_str());
	wgamuncnnlomisshist_2->Reset("ICES");
	for (int i = 0; i <= wgamuncnnlomisshist_2->GetNbinsX()+1; i++)
	  wgamuncnnlomisshist_2->SetBinContent(i,fabs(wgamcornnlomissuphist_2.back()->GetBinContent(i)-wgamcornnlomissdwhist_2.back()->GetBinContent(i))/(2*wgamcorewkhist.back()->GetBinContent(i)));
	wgamuncnnlomisshist_2->Smooth(1,"R");
	wgamuncnnlomisshist_2->SetName(("WG_NNLOMiss2_"+obs).c_str());

	// NLO QCD+EWK --> to be symmetrized
	wgamuncmixhist = (TH1*) wgamcormixuphist.back()->Clone(("wgamuncmix"+ext+"hist_"+obs).c_str());
	wgamuncmixhist->Reset("ICES");
	for (int i = 0; i <= wgamuncmixhist->GetNbinsX()+1; i++)
	  wgamuncmixhist->SetBinContent(i,fabs(wgamcormixuphist.back()->GetBinContent(i)-wgamcormixdwhist.back()->GetBinContent(i))/(2*wgamcorewkhist.back()->GetBinContent(i)));
	wgamuncmixhist->Smooth(1,"R");
	wgamuncmixhist->SetName(("WG_MIX_"+obs).c_str());
	
	// NLO PDF
	wgamuncpdfhist = (TH1*)wgamcorpdfhist.back()->Clone(("wgamuncpdf"+ext+"hist_"+obs).c_str());
	wgamuncpdfhist->Divide(wgamcorqcdhist.back());
	for (int i = 1; i <= wgamuncpdfhist->GetNbinsX(); i++)
	  wgamuncpdfhist->SetBinContent(i,fabs(wgamuncpdfhist->GetBinContent(i)-1.0));
	wgamuncpdfhist->SetName(("WG_PDF_"+obs).c_str());
      }
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
      zwjcorqcdshapeuphist.back()->Write();
      zwjcorqcdshapedwhist.back()->Write();
      zwjcorqcdprocuphist.back()->Write();
      zwjcorqcdprocdwhist.back()->Write();
      zwjcornnloewkuphist.back()->Write();
      zwjcornnloewkdwhist.back()->Write();
      zwjcorsudakovuphist_1.back()->Write();
      zwjcorsudakovdwhist_1.back()->Write();
      zwjcorsudakovuphist_2.back()->Write();
      zwjcorsudakovdwhist_2.back()->Write();
      zwjcornnlomissuphist_1.back()->Write();
      zwjcornnlomissdwhist_1.back()->Write();
      zwjcornnlomissuphist_2.back()->Write();
      zwjcornnlomissdwhist_2.back()->Write();
      zwjcormixuphist.back()->Write();
      zwjcormixdwhist.back()->Write();
      zwjcorpdfhist.back()->Write();

      if(addHistoForCutAndCount){
	zwjcorqcdscaleuphist_num.back()->Write();
	zwjcorqcdscaledwhist_num.back()->Write();
	zwjcorqcdshapeuphist_num.back()->Write();
	zwjcorqcdshapedwhist_num.back()->Write();
	zwjcorqcdprocuphist_num.back()->Write();
	zwjcorqcdprocdwhist_num.back()->Write();
	zwjcornnloewkuphist_num.back()->Write();
	zwjcornnloewkdwhist_num.back()->Write();
	zwjcorsudakovuphist_1_num.back()->Write();
	zwjcorsudakovdwhist_1_num.back()->Write();
	zwjcorsudakovuphist_2_num.back()->Write();
	zwjcorsudakovdwhist_2_num.back()->Write();
	zwjcornnlomissuphist_1_num.back()->Write();
	zwjcornnlomissdwhist_1_num.back()->Write();
	zwjcornnlomissuphist_2_num.back()->Write();
	zwjcornnlomissdwhist_2_num.back()->Write();
	zwjcormixuphist_num.back()->Write();
	zwjcormixdwhist_num.back()->Write();

	zwjcorqcdscaleuphist_den.back()->Write();
	zwjcorqcdscaledwhist_den.back()->Write();
	zwjcorqcdshapeuphist_den.back()->Write();
	zwjcorqcdshapedwhist_den.back()->Write();
	zwjcorqcdprocuphist_den.back()->Write();
	zwjcorqcdprocdwhist_den.back()->Write();
	zwjcornnloewkuphist_den.back()->Write();
	zwjcornnloewkdwhist_den.back()->Write();
	zwjcorsudakovuphist_1_den.back()->Write();
	zwjcorsudakovdwhist_1_den.back()->Write();
	zwjcorsudakovuphist_2_den.back()->Write();
	zwjcorsudakovdwhist_2_den.back()->Write();
	zwjcornnlomissuphist_1_den.back()->Write();
	zwjcornnlomissdwhist_1_den.back()->Write();
	zwjcornnlomissuphist_2_den.back()->Write();
	zwjcornnlomissdwhist_2_den.back()->Write();
	zwjcormixuphist_den.back()->Write();
	zwjcormixdwhist_den.back()->Write();
      }
      
      zwjuncqcdscalehist->Write();
      zwjuncqcdshapehist->Write();
      zwjuncqcdprochist->Write();
      zwjuncnnloewkhist->Write();
      zwjuncsudakovhist_1->Write();
      zwjuncsudakovhist_2->Write();
      zwjuncnnlomisshist_1->Write();
      zwjuncnnlomisshist_2->Write();
      zwjuncmixhist->Write();      
      zwjuncpdfhist->Write();
    }
  
    outputFile.cd();
    if(addZgamma){
      if(not outputFile.GetDirectory("TF_GJ"))
	outputFile.mkdir("TF_GJ");
      outputFile.cd("TF_GJ");
      gamcorhist.back()->Write();
      gamcorqcdhist.back()->Write();
      gamcorewkhist.back()->Write();

      if(not useNewTheoryUncertainty){
	
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

      else{
	
	gamcorqcdscaleuphist.back()->Write();
	gamcorqcdscaledwhist.back()->Write();
	gamcorqcdshapeuphist.back()->Write();
	gamcorqcdshapedwhist.back()->Write();
	gamcorqcdprocuphist.back()->Write();
	gamcorqcdprocdwhist.back()->Write();
	gamcornnloewkuphist.back()->Write();
	gamcornnloewkdwhist.back()->Write();
	gamcorsudakovuphist_1.back()->Write();
	gamcorsudakovdwhist_1.back()->Write();
	gamcorsudakovuphist_2.back()->Write();
	gamcorsudakovdwhist_2.back()->Write();
	gamcornnlomissuphist_1.back()->Write();
	gamcornnlomissdwhist_1.back()->Write();
	gamcornnlomissuphist_2.back()->Write();
	gamcornnlomissdwhist_2.back()->Write();
	gamcormixuphist.back()->Write();
	gamcormixdwhist.back()->Write();
	gamcorpdfhist.back()->Write();
	
	if(addHistoForCutAndCount){
	  gamcorqcdscaleuphist_num.back()->Write();
	  gamcorqcdscaledwhist_num.back()->Write();
	  gamcorqcdshapeuphist_num.back()->Write();
	  gamcorqcdshapedwhist_num.back()->Write();
	  gamcorqcdprocuphist_num.back()->Write();
	  gamcorqcdprocdwhist_num.back()->Write();
	  gamcornnloewkuphist_num.back()->Write();
	  gamcornnloewkdwhist_num.back()->Write();
	  gamcorsudakovuphist_1_num.back()->Write();
	  gamcorsudakovdwhist_1_num.back()->Write();
	  gamcorsudakovuphist_2_num.back()->Write();
	  gamcorsudakovdwhist_2_num.back()->Write();
	  gamcornnlomissuphist_1_num.back()->Write();
	  gamcornnlomissdwhist_1_num.back()->Write();
	  gamcornnlomissuphist_2_num.back()->Write();
	  gamcornnlomissdwhist_2_num.back()->Write();
	  gamcormixuphist_num.back()->Write();
	  gamcormixdwhist_num.back()->Write();
	  
	  gamcorqcdscaleuphist_den.back()->Write();
	  gamcorqcdscaledwhist_den.back()->Write();
	  gamcorqcdshapeuphist_den.back()->Write();
	  gamcorqcdshapedwhist_den.back()->Write();
	  gamcorqcdprocuphist_den.back()->Write();
	  gamcorqcdprocdwhist_den.back()->Write();
	  gamcornnloewkuphist_den.back()->Write();
	  gamcornnloewkdwhist_den.back()->Write();
	  gamcorsudakovuphist_1_den.back()->Write();
	  gamcorsudakovdwhist_1_den.back()->Write();
	  gamcorsudakovuphist_2_den.back()->Write();
	  gamcorsudakovdwhist_2_den.back()->Write();
	  gamcornnlomissuphist_1_den.back()->Write();
	  gamcornnlomissdwhist_1_den.back()->Write();
	  gamcornnlomissuphist_2_den.back()->Write();
	  gamcornnlomissdwhist_2_den.back()->Write();
	  gamcormixuphist_den.back()->Write();
	  gamcormixdwhist_den.back()->Write();
	}
	
	gamuncqcdscalehist->Write();
	gamuncqcdshapehist->Write();
	gamuncqcdprochist->Write();
	gamuncnnloewkhist->Write();
	gamuncsudakovhist_1->Write();
	gamuncsudakovhist_2->Write();
	gamuncnnlomisshist_1->Write();
	gamuncnnlomisshist_2->Write();
	gamuncmixhist->Write();      
	gamuncpdfhist->Write();
      }
    }

    outputFile.cd();
    if(addWgamma){
      if(not outputFile.GetDirectory("TF_WG"))
	outputFile.mkdir("TF_WG");
      outputFile.cd("TF_WG");
      wgamcorhist.back()->Write();
      wgamcorqcdhist.back()->Write();
      wgamcorewkhist.back()->Write();

      if(not useNewTheoryUncertainty){
	
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

      else{
	
	wgamcorqcdscaleuphist.back()->Write();
	wgamcorqcdscaledwhist.back()->Write();
	wgamcorqcdshapeuphist.back()->Write();
	wgamcorqcdshapedwhist.back()->Write();
	wgamcorqcdprocuphist.back()->Write();
	wgamcorqcdprocdwhist.back()->Write();
	wgamcornnloewkuphist.back()->Write();
	wgamcornnloewkdwhist.back()->Write();
	wgamcorsudakovuphist_1.back()->Write();
	wgamcorsudakovdwhist_1.back()->Write();
	wgamcorsudakovuphist_2.back()->Write();
	wgamcorsudakovdwhist_2.back()->Write();
	wgamcornnlomissuphist_1.back()->Write();
	wgamcornnlomissdwhist_1.back()->Write();
	wgamcornnlomissuphist_2.back()->Write();
	wgamcornnlomissdwhist_2.back()->Write();
	wgamcormixuphist.back()->Write();
	wgamcormixdwhist.back()->Write();
	wgamcorpdfhist.back()->Write();
	
	if(addHistoForCutAndCount){
	  wgamcorqcdscaleuphist_num.back()->Write();
	  wgamcorqcdscaledwhist_num.back()->Write();
	  wgamcorqcdshapeuphist_num.back()->Write();
	  wgamcorqcdshapedwhist_num.back()->Write();
	  wgamcorqcdprocuphist_num.back()->Write();
	  wgamcorqcdprocdwhist_num.back()->Write();
	  wgamcornnloewkuphist_num.back()->Write();
	  wgamcornnloewkdwhist_num.back()->Write();
	  wgamcorsudakovuphist_1_num.back()->Write();
	  wgamcorsudakovdwhist_1_num.back()->Write();
	  wgamcorsudakovuphist_2_num.back()->Write();
	  wgamcorsudakovdwhist_2_num.back()->Write();
	  wgamcornnlomissuphist_1_num.back()->Write();
	  wgamcornnlomissdwhist_1_num.back()->Write();
	  wgamcornnlomissuphist_2_num.back()->Write();
	  wgamcornnlomissdwhist_2_num.back()->Write();
	  wgamcormixuphist_num.back()->Write();
	  wgamcormixdwhist_num.back()->Write();
	  
	  wgamcorqcdscaleuphist_den.back()->Write();
	  wgamcorqcdscaledwhist_den.back()->Write();
	  wgamcorqcdshapeuphist_den.back()->Write();
	  wgamcorqcdshapedwhist_den.back()->Write();
	  wgamcorqcdprocuphist_den.back()->Write();
	  wgamcorqcdprocdwhist_den.back()->Write();
	  wgamcornnloewkuphist_den.back()->Write();
	  wgamcornnloewkdwhist_den.back()->Write();
	  wgamcorsudakovuphist_1_den.back()->Write();
	  wgamcorsudakovdwhist_1_den.back()->Write();
	  wgamcorsudakovuphist_2_den.back()->Write();
	  wgamcorsudakovdwhist_2_den.back()->Write();
	  wgamcornnlomissuphist_1_den.back()->Write();
	  wgamcornnlomissdwhist_1_den.back()->Write();
	  wgamcornnlomissuphist_2_den.back()->Write();
	  wgamcornnlomissdwhist_2_den.back()->Write();
	  wgamcormixuphist_den.back()->Write();
	  wgamcormixdwhist_den.back()->Write();
	}
	
	wgamuncqcdscalehist->Write();
	wgamuncqcdshapehist->Write();
	wgamuncqcdprochist->Write();
	wgamuncnnloewkhist->Write();
	wgamuncsudakovhist_1->Write();
	wgamuncsudakovhist_2->Write();
	wgamuncnnlomisshist_1->Write();
	wgamuncnnlomisshist_2->Write();
	wgamuncmixhist->Write();      
	wgamuncpdfhist->Write();
      }
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
  zmmcorhist.clear(); zeecorhist.clear(); wmncorhist.clear(); wencorhist.clear(); zwjcorhist.clear(); wgamcorhist.clear(); wgamcorhist.clear(); topmucorhist.clear(); topelcorhist.clear();
  topmucorbuphist.clear(); topmucorbdownhist.clear(); topelcorbuphist.clear(); topelcorbdownhist.clear();

  gamcorewkhist.clear(); gamcorqcdhist.clear(); gamcorre1hist.clear(); gamcorfa1hist.clear(); gamcorre2hist.clear(); gamcorfa2hist.clear(); gamcorpdfhist.clear(); gamcorfpchist.clear();
  gamcorqcdscaleuphist.clear(); gamcorqcdscaledwhist.clear(); gamcorqcdshapeuphist.clear(); gamcorqcdshapedwhist.clear(); gamcorqcdprocuphist.clear(); gamcorqcdprocdwhist.clear();
  gamcornnloewkuphist.clear(); gamcornnloewkdwhist.clear();   gamcorsudakovuphist_1.clear(); gamcorsudakovuphist_2.clear(); gamcorsudakovdwhist_1.clear(); gamcorsudakovdwhist_2.clear(); 
  gamcormixuphist.clear(); gamcormixdwhist.clear(); 

  zwjcorewkhist.clear(); zwjcorqcdhist.clear(); zwjcorre1hist.clear(); zwjcorre2hist.clear(); zwjcorfa1hist.clear(); zwjcorfa2hist.clear(); zwjcorpdfhist.clear();
  zwjcorqcdscaleuphist.clear(); zwjcorqcdscaledwhist.clear(); zwjcorqcdshapeuphist.clear(); zwjcorqcdshapedwhist.clear(); zwjcorqcdprocuphist.clear(); zwjcorqcdprocdwhist.clear();
  zwjcornnloewkuphist.clear(); zwjcornnloewkdwhist.clear();   zwjcorsudakovuphist_1.clear(); zwjcorsudakovuphist_2.clear(); zwjcorsudakovdwhist_1.clear(); zwjcorsudakovdwhist_2.clear(); 
  zwjcormixuphist.clear(); zwjcormixdwhist.clear(); 

  wgamcorewkhist.clear(); wgamcorqcdhist.clear(); wgamcorre1hist.clear(); wgamcorfa1hist.clear(); wgamcorre2hist.clear(); wgamcorfa2hist.clear(); wgamcorpdfhist.clear(); wgamcorfpchist.clear();
  wgamcorqcdscaleuphist.clear(); wgamcorqcdscaledwhist.clear(); wgamcorqcdshapeuphist.clear(); wgamcorqcdshapedwhist.clear(); wgamcorqcdprocuphist.clear(); wgamcorqcdprocdwhist.clear();
  wgamcornnloewkuphist.clear(); wgamcornnloewkdwhist.clear();   wgamcorsudakovuphist_1.clear(); wgamcorsudakovuphist_2.clear(); wgamcorsudakovdwhist_1.clear(); wgamcorsudakovdwhist_2.clear(); 
  wgamcormixuphist.clear(); wgamcormixdwhist.clear(); 

  zmmcorhist_num.clear(); zeecorhist_num.clear(); wmncorhist_num.clear(); wencorhist_num.clear(); gamcorhist_num.clear(); wgamcorhist_num.clear(); zwjcorhist_num.clear();
  zmmcorhist_den.clear(); zeecorhist_den.clear(); wmncorhist_den.clear(); wencorhist_den.clear(); gamcorhist_den.clear(); wgamcorhist_den.clear(); zwjcorhist_den.clear();
  topmucorhist_num.clear(); topelcorhist_num.clear();
  topmucorhist_den.clear(); topelcorhist_den.clear();
  
  gamcorewkhist_num.clear(); gamcorqcdhist_num.clear(); gamcorre1hist_num.clear(); gamcorfa1hist_num.clear(); gamcorre2hist_num.clear(); gamcorpdfhist_num.clear(); 
  gamcorfa2hist_num.clear(); gamcorfpchist_num.clear();
  gamcorewkhist_den.clear(); gamcorqcdhist_den.clear(); gamcorre1hist_den.clear(); gamcorfa1hist_den.clear(); gamcorre2hist_den.clear(); gamcorpdfhist_den.clear(); 
  gamcorfa2hist_den.clear(); gamcorfpchist_den.clear();
  gamcorqcdscaleuphist_num.clear(); gamcorqcdscaledwhist_num.clear(); gamcorqcdshapeuphist_num.clear(); gamcorqcdshapedwhist_num.clear(); gamcorqcdprocuphist_num.clear(); 
  gamcorqcdprocdwhist_num.clear();
  gamcornnloewkuphist_num.clear(); gamcornnloewkdwhist_num.clear();   gamcorsudakovuphist_1_num.clear(); gamcorsudakovuphist_2_num.clear(); gamcorsudakovdwhist_1_num.clear(); 
  gamcorsudakovdwhist_2_num.clear(); 
  gamcormixuphist_num.clear(); gamcormixdwhist_num.clear(); 
  gamcorqcdscaleuphist_den.clear(); gamcorqcdscaledwhist_den.clear(); gamcorqcdshapeuphist_den.clear(); gamcorqcdshapedwhist_den.clear(); gamcorqcdprocuphist_den.clear(); 
  gamcorqcdprocdwhist_den.clear();
  gamcornnloewkuphist_den.clear(); gamcornnloewkdwhist_den.clear();   gamcorsudakovuphist_1_den.clear(); gamcorsudakovuphist_2_den.clear(); gamcorsudakovdwhist_1_den.clear(); 
  gamcorsudakovdwhist_2_den.clear(); 
  gamcormixuphist_den.clear(); gamcormixdwhist_den.clear(); 

  zwjcorewkhist_num.clear(); zwjcorqcdhist_num.clear(); zwjcorre1hist_num.clear(); zwjcorfa1hist_num.clear(); zwjcorre2hist_num.clear(); zwjcorpdfhist_num.clear(); zwjcorfa2hist_num.clear(); 
  zwjcorewkhist_den.clear(); zwjcorqcdhist_den.clear(); zwjcorre1hist_den.clear(); zwjcorfa1hist_den.clear(); zwjcorre2hist_den.clear(); zwjcorpdfhist_den.clear(); zwjcorfa2hist_den.clear();
  zwjcorqcdscaleuphist_num.clear(); zwjcorqcdscaledwhist_num.clear(); zwjcorqcdshapeuphist_num.clear(); zwjcorqcdshapedwhist_num.clear(); zwjcorqcdprocuphist_num.clear(); 
  zwjcorqcdprocdwhist_num.clear();
  zwjcornnloewkuphist_num.clear(); zwjcornnloewkdwhist_num.clear();   zwjcorsudakovuphist_1_num.clear(); zwjcorsudakovuphist_2_num.clear(); zwjcorsudakovdwhist_1_num.clear(); 
  zwjcorsudakovdwhist_2_num.clear(); 
  zwjcormixuphist_num.clear(); zwjcormixdwhist_num.clear(); 
  zwjcorqcdscaleuphist_den.clear(); zwjcorqcdscaledwhist_den.clear(); zwjcorqcdshapeuphist_den.clear(); zwjcorqcdshapedwhist_den.clear(); zwjcorqcdprocuphist_den.clear(); 
  zwjcorqcdprocdwhist_den.clear();
  zwjcornnloewkuphist_den.clear(); zwjcornnloewkdwhist_den.clear();   zwjcorsudakovuphist_1_den.clear(); zwjcorsudakovuphist_2_den.clear(); zwjcorsudakovdwhist_1_den.clear(); 
  zwjcorsudakovdwhist_2_den.clear(); 
  zwjcormixuphist_den.clear(); zwjcormixdwhist_den.clear(); 

  wgamcorewkhist_num.clear(); wgamcorqcdhist_num.clear(); wgamcorre1hist_num.clear(); wgamcorfa1hist_num.clear(); wgamcorre2hist_num.clear(); wgamcorpdfhist_num.clear(); 
  wgamcorfa2hist_num.clear(); wgamcorfpchist_num.clear();
  wgamcorewkhist_den.clear(); wgamcorqcdhist_den.clear(); wgamcorre1hist_den.clear(); wgamcorfa1hist_den.clear(); wgamcorre2hist_den.clear(); wgamcorpdfhist_den.clear(); 
  wgamcorfa2hist_den.clear(); wgamcorfpchist_den.clear();
  wgamcorqcdscaleuphist_num.clear(); wgamcorqcdscaledwhist_num.clear(); wgamcorqcdshapeuphist_num.clear(); wgamcorqcdshapedwhist_num.clear(); wgamcorqcdprocuphist_num.clear(); 
  wgamcorqcdprocdwhist_num.clear();
  wgamcornnloewkuphist_num.clear(); wgamcornnloewkdwhist_num.clear();   wgamcorsudakovuphist_1_num.clear(); wgamcorsudakovuphist_2_num.clear(); wgamcorsudakovdwhist_1_num.clear(); 
  wgamcorsudakovdwhist_2_num.clear(); 
  wgamcormixuphist_num.clear(); wgamcormixdwhist_num.clear(); 
  wgamcorqcdscaleuphist_den.clear(); wgamcorqcdscaledwhist_den.clear(); wgamcorqcdshapeuphist_den.clear(); wgamcorqcdshapedwhist_den.clear(); wgamcorqcdprocuphist_den.clear(); 
  wgamcorqcdprocdwhist_den.clear();
  wgamcornnloewkuphist_den.clear(); wgamcornnloewkdwhist_den.clear();   wgamcorsudakovuphist_1_den.clear(); wgamcorsudakovuphist_2_den.clear(); wgamcorsudakovdwhist_1_den.clear(); 
  wgamcorsudakovdwhist_2_den.clear(); 
  wgamcormixuphist_den.clear(); wgamcormixdwhist_den.clear(); 

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
  if(gamcorqcdscaleupfile) gamcorqcdscaleupfile->Close();
  if(gamcorqcdscaledwfile) gamcorqcdscaledwfile->Close();
  if(gamcorqcdprocupfile) gamcorqcdprocupfile->Close();
  if(gamcorqcdprocdwfile) gamcorqcdprocdwfile->Close();
  if(gamcornnloewkupfile) gamcornnloewkupfile->Close();
  if(gamcornnloewkdwfile) gamcornnloewkdwfile->Close();
  if(gamcorsudakovupfile_1) gamcorsudakovupfile_1->Close();
  if(gamcorsudakovdwfile_1) gamcorsudakovdwfile_1->Close();
  if(gamcorsudakovupfile_2) gamcorsudakovupfile_2->Close();
  if(gamcorsudakovdwfile_2) gamcorsudakovdwfile_2->Close();
  if(gamcornnlomissupfile_1) gamcornnlomissupfile_1->Close();
  if(gamcornnlomissdwfile_1) gamcornnlomissdwfile_1->Close();
  if(gamcornnlomissupfile_2) gamcornnlomissupfile_2->Close();
  if(gamcornnlomissdwfile_2) gamcornnlomissdwfile_2->Close();
  if(gamcormixupfile) gamcormixupfile->Close();
  if(gamcormixdwfile) gamcormixdwfile->Close();

  if(wgamcorqcdfile) wgamcorqcdfile->Close();
  if(wgamcorewkfile) wgamcorewkfile->Close();
  if(wgamcorre1file) wgamcorre1file->Close();
  if(wgamcorre2file) wgamcorre2file->Close();
  if(wgamcorfa2file) wgamcorfa2file->Close();
  if(wgamcorpdffile) wgamcorpdffile->Close();
  if(wgamcorfpcfile) wgamcorfpcfile->Close();
  if(wgamcorqcdscaleupfile) wgamcorqcdscaleupfile->Close();
  if(wgamcorqcdscaledwfile) wgamcorqcdscaledwfile->Close();
  if(wgamcorqcdprocupfile) wgamcorqcdprocupfile->Close();
  if(wgamcorqcdprocdwfile) wgamcorqcdprocdwfile->Close();
  if(wgamcornnloewkupfile) wgamcornnloewkupfile->Close();
  if(wgamcornnloewkdwfile) wgamcornnloewkdwfile->Close();
  if(wgamcorsudakovupfile_1) wgamcorsudakovupfile_1->Close();
  if(wgamcorsudakovdwfile_1) wgamcorsudakovdwfile_1->Close();
  if(wgamcorsudakovupfile_2) wgamcorsudakovupfile_2->Close();
  if(wgamcorsudakovdwfile_2) wgamcorsudakovdwfile_2->Close();
  if(wgamcornnlomissupfile_1) wgamcornnlomissupfile_1->Close();
  if(wgamcornnlomissdwfile_1) wgamcornnlomissdwfile_1->Close();
  if(wgamcornnlomissupfile_2) wgamcornnlomissupfile_2->Close();
  if(wgamcornnlomissdwfile_2) wgamcornnlomissdwfile_2->Close();
  if(wgamcormixupfile) wgamcormixupfile->Close();
  if(wgamcormixdwfile) wgamcormixdwfile->Close();

  if(zwjcorqcdfile) zwjcorqcdfile->Close();
  if(zwjcorewkfile) zwjcorewkfile->Close();
  if(zwjcorre1file) zwjcorre1file->Close();
  if(zwjcorfa1file) zwjcorfa1file->Close();
  if(zwjcorre2file) zwjcorre2file->Close();
  if(zwjcorfa2file) zwjcorfa2file->Close();
  if(zwjcorpdffile) zwjcorpdffile->Close();
  if(zwjcorqcdscaleupfile) zwjcorqcdscaleupfile->Close();
  if(zwjcorqcdscaledwfile) zwjcorqcdscaledwfile->Close();
  if(zwjcorqcdprocupfile) zwjcorqcdprocupfile->Close();
  if(zwjcorqcdprocdwfile) zwjcorqcdprocdwfile->Close();
  if(zwjcornnloewkupfile) zwjcornnloewkupfile->Close();
  if(zwjcornnloewkdwfile) zwjcornnloewkdwfile->Close();
  if(zwjcorsudakovupfile_1) zwjcorsudakovupfile_1->Close();
  if(zwjcorsudakovdwfile_1) zwjcorsudakovdwfile_1->Close();
  if(zwjcorsudakovupfile_2) zwjcorsudakovupfile_2->Close();
  if(zwjcorsudakovdwfile_2) zwjcorsudakovdwfile_2->Close();
  if(zwjcornnlomissupfile_1) zwjcornnlomissupfile_1->Close();
  if(zwjcornnlomissdwfile_1) zwjcornnlomissdwfile_1->Close();
  if(zwjcornnlomissupfile_2) zwjcornnlomissupfile_2->Close();
  if(zwjcornnlomissdwfile_2) zwjcornnlomissdwfile_2->Close();
  if(zwjcormixupfile) zwjcormixupfile->Close();
  if(zwjcormixdwfile) zwjcormixdwfile->Close();

  if(topmucorbupfile) topmucorbupfile->Close();
  if(topmucorbdownfile) topmucorbdownfile->Close();
  if(topelcorbupfile) topelcorbupfile->Close();
  if(topelcorbdownfile) topelcorbdownfile->Close();
}


///////////////////////////// from TFs files to the template one                                                                                                                                      
void fillAndSaveCorrEWKHistograms(const vector<string> & observables, // observables to be considered
				  TFile & outputFile, // output file
				  const string & outDir,  // output directory
				  const Category & category,
				  const bool & addZWratio, 
				  const bool & addTop, 
				  const string & ext,
				  const bool & addHistoForCutAndCount = false){

  cout<<"Re-open file for correction histo"<<endl;
  TFile* zmmcorfile = TFile::Open((outDir+"/zewkmmcor"+ext+".root").c_str());
  TFile* zeecorfile = TFile::Open((outDir+"/zewkeecor"+ext+".root").c_str());
  TFile* wmncorfile = TFile::Open((outDir+"/wewkmncor"+ext+".root").c_str());
  TFile* wencorfile = TFile::Open((outDir+"/wewkencor"+ext+".root").c_str());
  TFile* zwjcorfile = TFile::Open((outDir+"/zwjewkcor"+ext+".root").c_str());

  // get histograms                                                                                                                                                                                   
  vector<TH1*> zmmcorhist, zeecorhist, wmncorhist, wencorhist, zwjcorhist;
  vector<TH1*> zmmcorhist_num, zeecorhist_num, wmncorhist_num, wencorhist_num, zwjcorhist_num;
  vector<TH1*> zmmcorhist_den, zeecorhist_den, wmncorhist_den, wencorhist_den, zwjcorhist_den;

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

    if(addZWratio){
      zwjcorhist.push_back( (TH1*)zwjcorfile->FindObjectAny(("zwjewkcor"+ext+"hist_"+obs).c_str()));
      zwjcorhist_num.push_back( (TH1*)zwjcorfile->FindObjectAny(("nhist_ewk_zwj"+ext+"_"+obs).c_str()));
      zwjcorhist_den.push_back( (TH1*)zwjcorfile->FindObjectAny(("dhist_ewk_zwj"+ext+"_"+obs).c_str()));
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

    if(addZWratio){
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
    }
    outputFile.cd();
  }

  // clear vectors
  zmmcorhist.clear(); zeecorhist.clear(); wmncorhist.clear(); wencorhist.clear(); zwjcorhist.clear(); 
  zmmcorhist_num.clear(); zeecorhist_num.clear(); wmncorhist_num.clear(); wencorhist_num.clear(); zwjcorhist_num.clear(); 
  zmmcorhist_den.clear(); zeecorhist_den.clear(); wmncorhist_den.clear(); wencorhist_den.clear(); zwjcorhist_den.clear(); 


  // Close files
  if(zmmcorfile) zmmcorfile->Close();
  if(zeecorfile) zeecorfile->Close();
  if(wmncorfile) wmncorfile->Close();
  if(wencorfile) wencorfile->Close();
  if(zwjcorfile) zwjcorfile->Close();
  
}
