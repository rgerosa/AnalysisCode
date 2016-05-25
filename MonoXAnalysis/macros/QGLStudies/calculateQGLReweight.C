void calcQGLWeight(TFile* inputFile, string category, vector<TH1*> & gam, vector<TH1*>& wln, vector<TH1*>& zll, vector<TH1*>& top){

  TH1* templateQGL_data_wmn = (TH1*) inputFile->FindObjectAny(("datahistwmn_QGL_AK8"+category).c_str());
  TH1* templateQGL_data_wen = (TH1*) inputFile->FindObjectAny(("datahistwen_QGL_AK8"+category).c_str());
  TH1* templateQGL_data_zee = (TH1*) inputFile->FindObjectAny(("datahistzee_QGL_AK8"+category).c_str());
  TH1* templateQGL_data_zmm = (TH1*) inputFile->FindObjectAny(("datahistzmm_QGL_AK8"+category).c_str());
  TH1* templateQGL_data_gam = (TH1*) inputFile->FindObjectAny(("datahistgam_QGL_AK8"+category).c_str());
  TH1* templateQGL_data_topmu = (TH1*) inputFile->FindObjectAny(("datahisttopmu_QGL_AK8"+category).c_str());
  TH1* templateQGL_data_topel = (TH1*) inputFile->FindObjectAny(("datahisttopel_QGL_AK8"+category).c_str());

  TH1* templateQGL_bkg_wmn   = (TH1*) inputFile->FindObjectAny(("vlbkghistwmn_QGL_AK8"+category).c_str());
  TH1* templateQGL_bkg_wen   = (TH1*) inputFile->FindObjectAny(("vlbkghistwen_QGL_AK8"+category).c_str());
  TH1* templateQGL_bkg_zee   = (TH1*) inputFile->FindObjectAny(("vllbkghistzee_QGL_AK8"+category).c_str());
  TH1* templateQGL_bkg_zmm   = (TH1*) inputFile->FindObjectAny(("vllbkghistzmm_QGL_AK8"+category).c_str());
  TH1* templateQGL_bkg_gam   = (TH1*) inputFile->FindObjectAny(("gbkghistgam_QGL_AK8"+category).c_str());
  TH1* templateQGL_bkg_topmu = (TH1*) inputFile->FindObjectAny(("tbkghisttopmu_QGL_AK8"+category).c_str());
  TH1* templateQGL_bkg_topel = (TH1*) inputFile->FindObjectAny(("tbkghisttopel_QGL_AK8"+category).c_str());

  
  TH1* templateQGL_data_wln = (TH1*) templateQGL_data_wmn->Clone(("templateQGL_data_wln"+category).c_str());
  templateQGL_data_wln->Add(templateQGL_data_wen);

  TH1* templateQGL_data_zll = (TH1*) templateQGL_data_zmm->Clone(("templateQGL_data_zll"+category).c_str());
  templateQGL_data_zll->Add(templateQGL_data_zee);

  TH1* templateQGL_data_top = (TH1*) templateQGL_data_topmu->Clone(("templateQGL_data_top"+category).c_str());
  templateQGL_data_top->Add(templateQGL_data_topel);

  TH1* templateQGL_bkg_wln = (TH1*) templateQGL_bkg_wmn->Clone(("templateQGL_bkg_wln"+category).c_str());
  templateQGL_bkg_wln->Add(templateQGL_bkg_wen);

  TH1* templateQGL_bkg_zll = (TH1*) templateQGL_bkg_zmm->Clone(("templateQGL_bkg_zll"+category).c_str());
  templateQGL_bkg_zll->Add(templateQGL_bkg_zee);

  TH1* templateQGL_bkg_top = (TH1*) templateQGL_bkg_topmu->Clone(("templateQGL_bkg_top"+category).c_str());
  templateQGL_bkg_top->Add(templateQGL_bkg_topel);

  // normalize to get shapes
  templateQGL_data_wln->Scale(1./templateQGL_data_wln->Integral());
  templateQGL_data_zll->Scale(1./templateQGL_data_zll->Integral());
  templateQGL_data_gam->Scale(1./templateQGL_data_gam->Integral());
  templateQGL_data_top->Scale(1./templateQGL_data_top->Integral());

  templateQGL_bkg_wln->Scale(1./templateQGL_bkg_wln->Integral());
  templateQGL_bkg_zll->Scale(1./templateQGL_bkg_zll->Integral());
  templateQGL_bkg_gam->Scale(1./templateQGL_bkg_gam->Integral());
  templateQGL_bkg_top->Scale(1./templateQGL_bkg_top->Integral());

  
  TH1* QGL_ratio_wln = (TH1*) templateQGL_data_wln->Clone(("QGL_weight_wln"+category).c_str());
  TH1* QGL_ratio_wln_tmp = (TH1*) templateQGL_bkg_wln->Clone(("QGL_weight_wln_tmp"+category).c_str());
  QGL_ratio_wln->Divide(templateQGL_bkg_wln);
  QGL_ratio_wln_tmp->Divide(templateQGL_bkg_wln);

  TH1* QGL_ratio_wln_up = (TH1*) QGL_ratio_wln->Clone(("QGL_weight_wln_up"+category).c_str());
  TH1* QGL_ratio_wln_dw = (TH1*) QGL_ratio_wln->Clone(("QGL_weight_wln_dw"+category).c_str());
  for(int iBin = 1; iBin <= QGL_ratio_wln_up->GetNbinsX(); iBin++){
    QGL_ratio_wln_up->SetBinContent(iBin,QGL_ratio_wln_up->GetBinContent(iBin)+
				    sqrt(pow(QGL_ratio_wln_up->GetBinError(iBin),2)+pow(QGL_ratio_wln_tmp->GetBinError(iBin),2)));
    QGL_ratio_wln_dw->SetBinContent(iBin,QGL_ratio_wln_dw->GetBinContent(iBin)-
				    sqrt(pow(QGL_ratio_wln_dw->GetBinError(iBin),2)+pow(QGL_ratio_wln_tmp->GetBinError(iBin),2)));
  }

  wln.push_back(QGL_ratio_wln);
  wln.push_back(QGL_ratio_wln_up);
  wln.push_back(QGL_ratio_wln_dw);

  QGL_ratio_wln->Write();

  TH1* QGL_ratio_zll = (TH1*) templateQGL_data_zll->Clone(("QGL_weight_zll"+category).c_str());
  TH1* QGL_ratio_zll_tmp = (TH1*) templateQGL_bkg_zll->Clone(("QGL_weight_zll_tmp"+category).c_str());
  QGL_ratio_zll->Divide(templateQGL_bkg_zll);
  QGL_ratio_zll_tmp->Divide(templateQGL_bkg_zll);

  TH1* QGL_ratio_zll_up = (TH1*) QGL_ratio_zll->Clone(("QGL_weight_zll_up"+category).c_str());
  TH1* QGL_ratio_zll_dw = (TH1*) QGL_ratio_zll->Clone(("QGL_weight_zll_dw"+category).c_str());
  for(int iBin = 1; iBin <= QGL_ratio_zll_up->GetNbinsX(); iBin++){
    QGL_ratio_zll_up->SetBinContent(iBin,QGL_ratio_zll_up->GetBinContent(iBin)+
				    sqrt(pow(QGL_ratio_zll_up->GetBinError(iBin),2)+pow(QGL_ratio_zll_tmp->GetBinError(iBin),2)));
    QGL_ratio_zll_dw->SetBinContent(iBin,QGL_ratio_zll_dw->GetBinContent(iBin)-
				    sqrt(pow(QGL_ratio_zll_dw->GetBinError(iBin),2)+pow(QGL_ratio_zll_tmp->GetBinError(iBin),2)));
  }

  zll.push_back(QGL_ratio_zll);
  zll.push_back(QGL_ratio_zll_up);
  zll.push_back(QGL_ratio_zll_dw);

  QGL_ratio_zll->Write();

  TH1* QGL_ratio_top = (TH1*) templateQGL_data_top->Clone(("QGL_weight_top"+category).c_str());
  TH1* QGL_ratio_top_tmp = (TH1*) templateQGL_bkg_top->Clone(("QGL_weight_top_tmp"+category).c_str());
  QGL_ratio_top->Divide(templateQGL_bkg_top);
  QGL_ratio_top_tmp->Divide(templateQGL_bkg_top);

  TH1* QGL_ratio_top_up = (TH1*) QGL_ratio_top->Clone(("QGL_weight_top_up"+category).c_str());
  TH1* QGL_ratio_top_dw = (TH1*) QGL_ratio_top->Clone(("QGL_weight_top_dw"+category).c_str());
  for(int iBin = 1; iBin <= QGL_ratio_top_up->GetNbinsX(); iBin++){
    QGL_ratio_top_up->SetBinContent(iBin,QGL_ratio_top_up->GetBinContent(iBin)+
				    sqrt(pow(QGL_ratio_top_up->GetBinError(iBin),2)+pow(QGL_ratio_top_tmp->GetBinError(iBin),2)));
    QGL_ratio_top_dw->SetBinContent(iBin,QGL_ratio_top_dw->GetBinContent(iBin)-
				    sqrt(pow(QGL_ratio_top_dw->GetBinError(iBin),2)+pow(QGL_ratio_top_tmp->GetBinError(iBin),2)));
  }

  top.push_back(QGL_ratio_top);
  top.push_back(QGL_ratio_top_up);
  top.push_back(QGL_ratio_top_dw);

  QGL_ratio_top->Write();


  TH1* QGL_ratio_gam = (TH1*) templateQGL_data_gam->Clone(("QGL_weight_gam"+category).c_str());
  TH1* QGL_ratio_gam_tmp = (TH1*) templateQGL_bkg_gam->Clone(("QGL_weight_gam_tmp"+category).c_str());
  QGL_ratio_gam->Divide(templateQGL_bkg_gam);
  QGL_ratio_gam_tmp->Divide(templateQGL_bkg_gam);

  TH1* QGL_ratio_gam_up = (TH1*) QGL_ratio_gam->Clone(("QGL_weight_gam_up"+category).c_str());
  TH1* QGL_ratio_gam_dw = (TH1*) QGL_ratio_gam->Clone(("QGL_weight_gam_dw"+category).c_str());
  for(int iBin = 1; iBin <= QGL_ratio_gam_up->GetNbinsX(); iBin++){
    QGL_ratio_gam_up->SetBinContent(iBin,QGL_ratio_gam_up->GetBinContent(iBin)+
				    sqrt(pow(QGL_ratio_gam_up->GetBinError(iBin),2)+pow(QGL_ratio_gam_tmp->GetBinError(iBin),2)));
    QGL_ratio_gam_dw->SetBinContent(iBin,QGL_ratio_gam_dw->GetBinContent(iBin)-
				    sqrt(pow(QGL_ratio_gam_dw->GetBinError(iBin),2)+pow(QGL_ratio_gam_tmp->GetBinError(iBin),2)));
  }

  gam.push_back(QGL_ratio_gam);
  gam.push_back(QGL_ratio_gam_up);
  gam.push_back(QGL_ratio_gam_dw);

  QGL_ratio_gam->Write();

}

void fillHisto(TH2* histo, TH1* histo_1 = NULL, TH1* histo_2 = NULL, TH1* histo_3 = NULL){

  for(int iBinX = 0; iBinX < histo->GetNbinsX(); iBinX++){
    for(int iBinY = 0; iBinY < histo->GetNbinsY(); iBinY++){
      if(iBinX == 0){
	histo->SetBinContent(iBinX+1,iBinY+1,histo_1->GetBinContent(iBinY+1));	
      }
      else if(iBinX == 1){
	histo->SetBinContent(iBinX+1,iBinY+1,histo_2->GetBinContent(iBinY+1));
      }
      else if(iBinX == 2){
	histo->SetBinContent(iBinX+1,iBinY+1,histo_3->GetBinContent(iBinY+1));
      }
    }
  }
}

void calculateQGLReweight(string inputFileName, string outputFileName){

  TFile* inputFile = TFile::Open(inputFileName.c_str());
  TFile* outputFile = new TFile(outputFileName.c_str(),"RECREATE");

  vector<TH1*> gam_1;
  vector<TH1*> top_1;
  vector<TH1*> zll_1;
  vector<TH1*> wln_1;

  calcQGLWeight(inputFile,"",gam_1,wln_1,zll_1,top_1);
  /*
  vector<TH1*> gam_2;
  vector<TH1*> top_2;
  vector<TH1*> zll_2;
  vector<TH1*> wln_2;

  calcQGLWeight(inputFile,"_2",gam_2,wln_2,zll_2,top_2);

  vector<TH1*> gam_3;
  vector<TH1*> top_3;
  vector<TH1*> zll_3;
  vector<TH1*> wln_3;

  calcQGLWeight(inputFile,"_3",gam_3,wln_3,zll_3,top_3);
  */  
  vector<double> binX = {100.,2000.};
  double* array       = const_cast<double*>(zll_1.at(0)->GetXaxis()->GetXbins()->GetArray());

  TH2* QGL_weight_Z    = dynamic_cast<TH2*>(new TH2D("QGL_weight_Z","",binX.size()-1,&binX[0],zll_1.at(0)->GetNbinsX(),array));
  //  fillHisto(QGL_weight_Z,zll_1.at(0),zll_2.at(0),zll_3.at(0));
  fillHisto(QGL_weight_Z,zll_1.at(0));
  QGL_weight_Z->Write();

  TH2* QGL_weight_Z_up = dynamic_cast<TH2*>(new TH2D("QGL_weight_Z_up","",binX.size()-1,&binX[0],zll_1.at(0)->GetNbinsX(),array));
  //  fillHisto(QGL_weight_Z_up,zll_1.at(1),zll_2.at(1),zll_3.at(1));
  fillHisto(QGL_weight_Z_up,zll_1.at(1));
  QGL_weight_Z_up->Write();

  TH2* QGL_weight_Z_dw = dynamic_cast<TH2*>(new TH2D("QGL_weight_Z_dw","",binX.size()-1,&binX[0],zll_1.at(0)->GetNbinsX(),array));
  fillHisto(QGL_weight_Z_dw,zll_1.at(2));
  //  fillHisto(QGL_weight_Z_dw,zll_1.at(2),zll_2.at(2),zll_3.at(2));
  QGL_weight_Z_dw->Write();

  TH2* QGL_weight_W = dynamic_cast<TH2*>(new TH2D("QGL_weight_W","",binX.size()-1,&binX[0],wln_1.at(0)->GetNbinsX(),array));
  fillHisto(QGL_weight_W,wln_1.at(0));
  //  fillHisto(QGL_weight_W,wln_1.at(0),wln_2.at(0),wln_3.at(0));
  QGL_weight_W->Write();

  TH2* QGL_weight_W_up = dynamic_cast<TH2*>(new TH2D("QGL_weight_W_up","",binX.size()-1,&binX[0],wln_1.at(0)->GetNbinsX(),array));
  //  fillHisto(QGL_weight_W_up,wln_1.at(1),wln_2.at(1),wln_3.at(1));
  fillHisto(QGL_weight_W_up,wln_1.at(1));
  QGL_weight_W_up->Write();
  TH2* QGL_weight_W_dw = dynamic_cast<TH2*>(new TH2D("QGL_weight_W_dw","",binX.size()-1,&binX[0],wln_1.at(0)->GetNbinsX(),array));
  fillHisto(QGL_weight_W_dw,wln_1.at(2));
  //  fillHisto(QGL_weight_W_dw,wln_1.at(2),wln_2.at(2),wln_3.at(2));
  QGL_weight_W_dw->Write();

  TH2* QGL_weight_T = dynamic_cast<TH2*>(new TH2D("QGL_weight_T","",binX.size()-1,&binX[0],top_1.at(0)->GetNbinsX(),array));
  fillHisto(QGL_weight_T,top_1.at(0));
  //  fillHisto(QGL_weight_T,top_1.at(0),top_2.at(0),top_3.at(0));
  QGL_weight_T->Write();

  TH2* QGL_weight_T_up = dynamic_cast<TH2*>(new TH2D("QGL_weight_T_up","",binX.size()-1,&binX[0],top_1.at(0)->GetNbinsX(),array));
  fillHisto(QGL_weight_T_up,top_1.at(1));
  //  fillHisto(QGL_weight_T_up,top_1.at(1),top_2.at(1),top_3.at(1));
  QGL_weight_T_up->Write();

  TH2* QGL_weight_T_dw = dynamic_cast<TH2*>(new TH2D("QGL_weight_T_dw","",binX.size()-1,&binX[0],top_1.at(0)->GetNbinsX(),array));
  fillHisto(QGL_weight_T_dw,top_1.at(2));
  //  fillHisto(QGL_weight_T_dw,top_1.at(2),top_2.at(2),top_3.at(2)); 
  QGL_weight_T_dw->Write();

  TH2* QGL_weight_G = dynamic_cast<TH2*>(new TH2D("QGL_weight_G","",binX.size()-1,&binX[0],gam_1.at(0)->GetNbinsX(),array));
  fillHisto(QGL_weight_G,gam_1.at(0));
  //  fillHisto(QGL_weight_G,gam_1.at(0),gam_2.at(0),gam_3.at(0));
  QGL_weight_G->Write();

  TH2* QGL_weight_G_up = dynamic_cast<TH2*>(new TH2D("QGL_weight_G_up","",binX.size()-1,&binX[0],gam_1.at(0)->GetNbinsX(),array));
  fillHisto(QGL_weight_G_up,gam_1.at(1));
  //  fillHisto(QGL_weight_G_up,gam_1.at(1),gam_2.at(1),gam_3.at(1));
  QGL_weight_G_up->Write();
  TH2* QGL_weight_G_dw = dynamic_cast<TH2*>(new TH2D("QGL_weight_G_dw","",binX.size()-1,&binX[0],gam_1.at(0)->GetNbinsX(),array));
  fillHisto(QGL_weight_G_dw,gam_1.at(2));
  //  fillHisto(QGL_weight_G_dw,gam_1.at(2),gam_2.at(2),gam_3.at(2));
  QGL_weight_G_dw->Write();

  outputFile->Close();
}
