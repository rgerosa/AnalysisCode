#include "../CMS_lumi.h"

vector<float> metBin    = {200.,250.,300,350,400,450,500,550,600,650,700,800,900,1000,1125,1250,1400.};
vector<float> mjjBin    = {200,400,600,900,1200,1500,2000,2750,3500,5000};
static bool symmetrize = true;
enum class Sample {zjet,wjet,gam};
enum class Category {monojet,VBF};

float zjet_ren_threshold_max = 1.09;
float zjet_ren_threshold_min = 1.06;
float wjet_ren_threshold_min = 1.06;
float zjet_fac_threshold_min = 1.04;
float wjet_fac_threshold_min = 1.04;

////
void makeProcessUncertaintyFromGen(string inputDIR,  string outputDIR, Sample sample, Category category){

  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  TRandom3* rand = new TRandom3();

  string postfix;
  if(sample == Sample::gam) postfix = "gam";
  else if(sample == Sample::wjet) postfix = "wjet";
  else if(sample == Sample::zjet) postfix = "zjet";

  float scale = 1;
  if(sample == Sample::zjet)
    scale = 3;

  system(("mkdir -p "+outputDIR).c_str());
  vector<TTree*> treeList;
  vector<TFile*> fileList;

  system(("ls "+inputDIR+" | grep root > list.txt").c_str());
  ifstream infileList ("list.txt");
  string line;
  if(infileList.is_open()){
    while(!infileList.eof()){
      getline(infileList,line);
      if(not TString(line).Contains("root") or line == "") continue;
      fileList.push_back(TFile::Open((inputDIR+"/"+line).c_str(),"READ"));
      treeList.push_back((TTree*) fileList.back()->Get("gentree/tree"));
    }
  }

  infileList.close();
  system("rm list.txt");

  // calculate sumwgt                                                                                                                                                                            
  vector<double> sumwgt;
  int ifile = 0;
  for(auto tree : treeList){
    TTreeReader reader (tree);
    TTreeReaderValue<float> xsec  (reader,"xsec");
    TTreeReaderValue<float> wgt   (reader,"wgt");
    TTreeReaderValue<int>   wzid  (reader,"wzid");

    cout<<"Calculate sumwgt for LO file "<<fileList.at(ifile)->GetName()<<endl;
    double sum = 0;
    while(reader.Next()){
      // filter away bad events                                                                                                                                                                       
      if(sample == Sample::zjet and fabs(*wzid) != 23) continue;
      if(sample == Sample::wjet and fabs(*wzid) != 24) continue;
      if(sample == Sample::gam  and fabs(*wzid) != 21) continue;
      sum += *wgt;
    }
    cout<<"Tree NLO with entries "<<tree->GetEntries()<<" sumwgt "<<sum<<endl;
    sumwgt.push_back(sum);
    ifile++;
  }
  

  
  // loop on the events
  vector<TH1F*> bosonpt_qcd_ren;
  vector<TH1F*> mjj_qcd_ren;
  vector<TH1F*> bosonpt_qcd_fac;
  vector<TH1F*> mjj_qcd_fac;
  vector<TH1F*> bosonpt_pdf;
  vector<TH1F*> mjj_pdf;

  // Loop on trees                                                                                                                                                                            
  ifile = 0;
  for(auto tree: treeList){
    
    TTreeReader reader (tree);
    TTreeReaderValue<float> xsec  (reader,"xsec");
    TTreeReaderValue<float> wgt   (reader,"wgt");
    TTreeReaderValue<int>   wzid  (reader,"wzid");
    TTreeReaderValue<unsigned int>   njetsinc  (reader,"njetsinc");
    TTreeReaderValue<unsigned int>   njets  (reader,"njets");
    TTreeReaderValue<float> wzmass  (reader,"wzmass");
    TTreeReaderValue<float> wzpt    (reader,"wzpt");
    TTreeReaderValue<float> wzpt_lhe    (reader,"mvpt");
    TTreeReaderValue<float> wzeta   (reader,"wzeta");
    TTreeReaderValue<float> wzphi   (reader,"wzphi");
    TTreeReaderValue<vector<float> > wgtqcd  (reader,"wgtqcd");
    TTreeReaderValue<vector<float> > wgtpdf  (reader,"wgtpdf");
    TTreeReaderValue<int> l1id   (reader,"l1id");
    TTreeReaderValue<float> l1pt   (reader,"l1pt");
    TTreeReaderValue<float> l1eta  (reader,"l1eta");
    TTreeReaderValue<float> l1phi  (reader,"l1phi");
    TTreeReaderValue<int> l2id   (reader,"l2id");
    TTreeReaderValue<float> l2pt   (reader,"l2pt");
    TTreeReaderValue<float> l2eta  (reader,"l2eta");
    TTreeReaderValue<float> l2phi  (reader,"l2phi");
    TTreeReaderValue<vector<float> > jetpt  (reader,"jetpt");
    TTreeReaderValue<vector<float> > jeteta (reader,"jeteta");
    TTreeReaderValue<vector<float> > jetphi (reader,"jetphi");
    TTreeReaderValue<vector<float> > jetmass (reader,"jetmass");


    cout<<"Loop on file "<<fileList.at(ifile)->GetName()<<endl;
    while(reader.Next()){
      if(ifile == 0 and bosonpt_qcd_ren.size() == 0 and mjj_qcd_ren.size() == 0){
	for(size_t iqcd = 0; iqcd < wgtqcd->size(); iqcd++){
	  if(iqcd == 5 or iqcd == 7) continue; // skip extreme scale variations
	  if(iqcd == 3 or iqcd == 6 or iqcd == 0){
	    bosonpt_qcd_ren.push_back(new TH1F(Form("bosonpt_iqcd_ren_%d",int(iqcd)),"",metBin.size()-1,&metBin[0]));
	    bosonpt_qcd_ren.back()->Sumw2();	  
	    mjj_qcd_ren.push_back(new TH1F(Form("mjj_iqcd_ren_%d",int(iqcd)),"",mjjBin.size()-1,&mjjBin[0]));
	    mjj_qcd_ren.back()->Sumw2();	  
	  }
	  if(iqcd == 0 or iqcd == 1 or iqcd == 2){
	    bosonpt_qcd_fac.push_back(new TH1F(Form("bosonpt_iqcd_fac_%d",int(iqcd)),"",metBin.size()-1,&metBin[0]));
	    bosonpt_qcd_fac.back()->Sumw2();	  
	    mjj_qcd_fac.push_back(new TH1F(Form("mjj_iqcd_fac_%d",int(iqcd)),"",mjjBin.size()-1,&mjjBin[0]));
	    mjj_qcd_fac.back()->Sumw2();	  
	  }
	}	
      }
      if(ifile == 0 and bosonpt_pdf.size() == 0 and mjj_pdf.size() == 0){
	for(size_t ipdf = 0; ipdf < wgtpdf->size(); ipdf++){
	  bosonpt_pdf.push_back(new TH1F(Form("bosonpt_ipdf_%d",int(ipdf)),"",metBin.size()-1,&metBin[0]));
	  bosonpt_pdf.back()->Sumw2();	  
	  mjj_pdf.push_back(new TH1F(Form("mjj_ipdf_%d",int(ipdf)),"",mjjBin.size()-1,&mjjBin[0]));
	  mjj_pdf.back()->Sumw2();	  
	}	
      }

      // filter away bad events                                                                                                                                                                    
      if(sample == Sample::zjet and fabs(*wzid) != 23) continue;
      else if(sample == Sample::wjet and fabs(*wzid) != 24) continue;
      else if(sample == Sample::gam  and fabs(*wzid) != 21) continue;

      if(sample == Sample::wjet and fabs(*l1id) == 15 and fabs(*l2id) == 16) continue;
      if(sample == Sample::wjet and fabs(*l2id) == 15 and fabs(*l1id) == 16) continue;
      
      if(*wzpt < 200) continue;
      
      vector<TLorentzVector> jets;
      for(size_t ijet = 0; ijet < jetpt->size(); ijet++){
        TLorentzVector jet4V; jet4V.SetPtEtaPhiM(jetpt->at(ijet),jeteta->at(ijet),jetphi->at(ijet),jetmass->at(ijet));
	float dphi1 = fabs(jetphi->at(ijet)-*l1phi);
	float dphi2 = fabs(jetphi->at(ijet)-*l2phi);
	if(dphi1 > TMath::Pi())
	  dphi1 = 2*TMath::Pi()-dphi1;
	if(dphi2 > TMath::Pi())
	  dphi2 = 2*TMath::Pi()-dphi2;
	if(sqrt(dphi1*dphi1+fabs(jeteta->at(ijet)-*l1eta)*fabs(jeteta->at(ijet)-*l1eta)) > 0.4 and
	   sqrt(dphi2*dphi2+fabs(jeteta->at(ijet)-*l2eta)*fabs(jeteta->at(ijet)-*l2eta)) > 0.4
	   ){ // check cleaning                                                                                                                                                                     
	  jets.push_back(jet4V);
	}
      }
      
      if(category == Category::monojet){
	if(jets.size() < 1) continue;            
	if(jets.at(0).Pt() < 100) continue;
	if(fabs(jets.at(0).Eta()) > 2.5) continue;      
      }
      else if(category == Category::VBF){
	if(jets.size() < 2) continue;
	if(jets.at(0).Pt() < 80) continue;
	if(jets.at(1).Pt() < 40) continue;
	if(fabs(jets.at(0).Eta()) > 4.7) continue;
	if(fabs(jets.at(1).Eta()) > 4.7) continue;
	if(jets.at(0).Eta()*jets.at(1).Eta() > 0) continue;
	if(fabs(jets.at(0).Eta()-jets.at(1).Eta()) < 1) continue;
	if((jets.at(0)+jets.at(1)).M() < 200) continue;
	if(fabs(jets.at(0).DeltaPhi(jets.at(1))) > 1.5) continue;
      }      

      if(sample == Sample::wjet and *wzpt_lhe >= 400 and fabs(*wgt) > 100) continue;
      if(sample == Sample::wjet and *wzpt_lhe >= 600 and fabs(*wgt) > 10)  continue;
      if(sample == Sample::zjet and *wzpt_lhe >= 400 and fabs(*wgt) > 10)  continue;
      if(sample == Sample::zjet and *wzpt_lhe >= 650 and fabs(*wgt) > 1)   continue;
      // QCD 
      int ipos_ren = 0;
      int ipos_fac = 0;

      float bosonptval = min(*wzpt,metBin.back()-1);
      float mjjval = 0;
      if(category == Category::VBF)
	mjjval = min((jets.at(0)+jets.at(1)).M(),double(mjjBin.back()-1));


      for(size_t iqcd = 0; iqcd < wgtqcd->size(); iqcd++){
	
	if(iqcd == 3 or iqcd == 6 or iqcd == 0){
	  bosonpt_qcd_ren.at(ipos_ren)->Fill(bosonptval,*xsec*(*wgt)*scale*(wgtqcd->at(iqcd)/(*wgt))/sumwgt.at(ifile));
	  if(category == Category::VBF and bosonptval > 200)
	    mjj_qcd_ren.at(ipos_ren)->Fill(mjjval,*xsec*(*wgt)*scale*(wgtqcd->at(iqcd)/(*wgt))/sumwgt.at(ifile));		  
	  ipos_ren++;
	}
	if(iqcd == 0 or iqcd == 1 or iqcd == 2){
	  bosonpt_qcd_fac.at(ipos_fac)->Fill(bosonptval,*xsec*(*wgt)*scale*(wgtqcd->at(iqcd)/(*wgt))/sumwgt.at(ifile));
	  if(category == Category::VBF and bosonptval > 200)
	    mjj_qcd_fac.at(ipos_fac)->Fill(mjjval,*xsec*(*wgt)*scale*(wgtqcd->at(iqcd)/(*wgt))/sumwgt.at(ifile));		  	
	  ipos_fac++;
	}
      }
      // PDF 
      for(size_t ipdf = 0; ipdf < wgtpdf->size(); ipdf++){	  
	if(category == Category::VBF and bosonptval > 200)
	  mjj_pdf.at(ipdf)->Fill(mjjval,*xsec*(*wgt)*scale*(wgtpdf->at(ipdf)/(*wgt))/sumwgt.at(ifile));		
	bosonpt_pdf.at(ipdf)->Fill(bosonptval,*xsec*(*wgt)*scale*(wgtpdf->at(ipdf)/(*wgt))/sumwgt.at(ifile));
      }            
    } 
    ifile++;
  }

  // Build the renormalization scale uncertainty
  TH1F* bosonpt_qcd_ren_uncup = (TH1F*) bosonpt_qcd_ren.at(0)->Clone("bosonpt_qcd_ren_uncup");
  bosonpt_qcd_ren_uncup->Reset();
  TH1F* bosonpt_qcd_ren_uncdw = (TH1F*) bosonpt_qcd_ren.at(0)->Clone("bosonpt_qcd_ren_uncdw");
  bosonpt_qcd_ren_uncdw->Reset();
  TH1F* mjj_qcd_ren_uncup = (TH1F*) mjj_qcd_ren.at(0)->Clone("mjj_qcd_ren_uncup");
  mjj_qcd_ren_uncup->Reset();
  TH1F* mjj_qcd_ren_uncdw = (TH1F*) mjj_qcd_ren.at(0)->Clone("mjj_qcd_ren_uncdw");
  mjj_qcd_ren_uncdw->Reset();

  for(int iBin = 0; iBin < bosonpt_qcd_ren_uncup->GetNbinsX(); iBin++){
    float maxValue = -100;
    float minValue = 999999;
    for(int ihisto = 0; ihisto < bosonpt_qcd_ren.size(); ihisto++){
      if(bosonpt_qcd_ren.at(ihisto)->GetBinContent(iBin+1) > maxValue) maxValue = bosonpt_qcd_ren.at(ihisto)->GetBinContent(iBin+1);
      if(bosonpt_qcd_ren.at(ihisto)->GetBinContent(iBin+1) < minValue) minValue = bosonpt_qcd_ren.at(ihisto)->GetBinContent(iBin+1);
    }
    if(not symmetrize){
      bosonpt_qcd_ren_uncup->SetBinContent(iBin+1,maxValue/bosonpt_qcd_ren.at(0)->GetBinContent(iBin+1));
      bosonpt_qcd_ren_uncdw->SetBinContent(iBin+1,minValue/bosonpt_qcd_ren.at(0)->GetBinContent(iBin+1));    
    }
    else{
      bosonpt_qcd_ren_uncup->SetBinContent(iBin+1,1+(fabs(maxValue/bosonpt_qcd_ren.at(0)->GetBinContent(iBin+1)-1)+fabs(minValue/bosonpt_qcd_ren.at(0)->GetBinContent(iBin+1)-1))/2);
      bosonpt_qcd_ren_uncdw->SetBinContent(iBin+1,1-(fabs(maxValue/bosonpt_qcd_ren.at(0)->GetBinContent(iBin+1)-1)+fabs(minValue/bosonpt_qcd_ren.at(0)->GetBinContent(iBin+1)-1))/2);
    }
  }

  bosonpt_qcd_ren_uncup->Smooth();
  bosonpt_qcd_ren_uncdw->Smooth();

  if(category == Category::VBF){
    for(int iBin = 0; iBin < mjj_qcd_ren_uncup->GetNbinsX(); iBin++){
      float maxValue = -100;
      float minValue = 999999;
      for(int ihisto = 0; ihisto < mjj_qcd_ren.size(); ihisto++){
	if(mjj_qcd_ren.at(ihisto)->GetBinContent(iBin+1) > maxValue) maxValue = mjj_qcd_ren.at(ihisto)->GetBinContent(iBin+1);
	if(mjj_qcd_ren.at(ihisto)->GetBinContent(iBin+1) < minValue) minValue = mjj_qcd_ren.at(ihisto)->GetBinContent(iBin+1);
      }
      if(not symmetrize){
	mjj_qcd_ren_uncup->SetBinContent(iBin+1,maxValue/mjj_qcd_ren.at(0)->GetBinContent(iBin+1));
	mjj_qcd_ren_uncdw->SetBinContent(iBin+1,minValue/mjj_qcd_ren.at(0)->GetBinContent(iBin+1));    
      }
      else{
	mjj_qcd_ren_uncup->SetBinContent(iBin+1,1+(fabs(maxValue/mjj_qcd_ren.at(0)->GetBinContent(iBin+1)-1)+fabs(minValue/mjj_qcd_ren.at(0)->GetBinContent(iBin+1)-1))/2);
	mjj_qcd_ren_uncdw->SetBinContent(iBin+1,1-(fabs(maxValue/mjj_qcd_ren.at(0)->GetBinContent(iBin+1)-1)+fabs(minValue/mjj_qcd_ren.at(0)->GetBinContent(iBin+1)-1))/2);
      }

      float number = rand->Uniform(0,0.01);
      if(sample == Sample::zjet and mjj_qcd_ren_uncup->GetBinContent(iBin+1) < zjet_ren_threshold_min) mjj_qcd_ren_uncup->SetBinContent(iBin+1,zjet_ren_threshold_min+number);
      if(sample == Sample::zjet and mjj_qcd_ren_uncup->GetBinContent(iBin+1) > zjet_ren_threshold_max) mjj_qcd_ren_uncup->SetBinContent(iBin+1,zjet_ren_threshold_max+number);
      if(sample == Sample::wjet and mjj_qcd_ren_uncup->GetBinContent(iBin+1) < wjet_ren_threshold_min) mjj_qcd_ren_uncup->SetBinContent(iBin+1,wjet_ren_threshold_min+number);
      if(sample == Sample::zjet and mjj_qcd_ren_uncdw->GetBinContent(iBin+1) > 1-(zjet_ren_threshold_min-1)) mjj_qcd_ren_uncdw->SetBinContent(iBin+1,1-(zjet_ren_threshold_min-1)-number);
      if(sample == Sample::zjet and mjj_qcd_ren_uncdw->GetBinContent(iBin+1) < 1-(zjet_ren_threshold_max-1)) mjj_qcd_ren_uncdw->SetBinContent(iBin+1,1-(zjet_ren_threshold_max-1)-number);
      if(sample == Sample::wjet and mjj_qcd_ren_uncdw->GetBinContent(iBin+1) > 1-(wjet_ren_threshold_min-1)) mjj_qcd_ren_uncdw->SetBinContent(iBin+1,1-(wjet_ren_threshold_min-1)-number);
    }

    mjj_qcd_ren_uncup->Smooth();
    mjj_qcd_ren_uncdw->Smooth();    
  }
  
  // Build the factorization scale uncertainty
  TH1F* bosonpt_qcd_fac_uncup = (TH1F*) bosonpt_qcd_fac.at(0)->Clone("bosonpt_qcd_fac_uncup");
  bosonpt_qcd_fac_uncup->Reset();
  TH1F* bosonpt_qcd_fac_uncdw = (TH1F*) bosonpt_qcd_fac.at(0)->Clone("bosonpt_qcd_fac_uncdw");
  bosonpt_qcd_fac_uncdw->Reset();
  TH1F* mjj_qcd_fac_uncup = (TH1F*) mjj_qcd_fac.at(0)->Clone("mjj_qcd_fac_uncup");
  mjj_qcd_fac_uncup->Reset();
  TH1F* mjj_qcd_fac_uncdw = (TH1F*) mjj_qcd_fac.at(0)->Clone("mjj_qcd_fac_uncdw");
  mjj_qcd_fac_uncdw->Reset();

  for(int iBin = 0; iBin < bosonpt_qcd_fac_uncup->GetNbinsX(); iBin++){
    float maxValue = -100;
    float minValue = 999999;
    for(int ihisto = 0; ihisto < bosonpt_qcd_fac.size(); ihisto++){
      if(bosonpt_qcd_fac.at(ihisto)->GetBinContent(iBin+1) > maxValue) maxValue = bosonpt_qcd_fac.at(ihisto)->GetBinContent(iBin+1);
      if(bosonpt_qcd_fac.at(ihisto)->GetBinContent(iBin+1) < minValue) minValue = bosonpt_qcd_fac.at(ihisto)->GetBinContent(iBin+1);
    }

    if(not symmetrize){
      bosonpt_qcd_fac_uncup->SetBinContent(iBin+1,maxValue/bosonpt_qcd_fac.at(0)->GetBinContent(iBin+1));
      bosonpt_qcd_fac_uncdw->SetBinContent(iBin+1,minValue/bosonpt_qcd_fac.at(0)->GetBinContent(iBin+1));    
    }
    else{
      bosonpt_qcd_fac_uncup->SetBinContent(iBin+1,1+(fabs(maxValue/bosonpt_qcd_fac.at(0)->GetBinContent(iBin+1)-1)+fabs(minValue/bosonpt_qcd_fac.at(0)->GetBinContent(iBin+1)-1))/2);
      bosonpt_qcd_fac_uncdw->SetBinContent(iBin+1,1-(fabs(maxValue/bosonpt_qcd_fac.at(0)->GetBinContent(iBin+1)-1)+fabs(minValue/bosonpt_qcd_fac.at(0)->GetBinContent(iBin+1)-1))/2);
    }
  }

  bosonpt_qcd_fac_uncup->Smooth();
  bosonpt_qcd_fac_uncdw->Smooth();

  if(category == Category::VBF){
    for(int iBin = 0; iBin < mjj_qcd_fac_uncup->GetNbinsX(); iBin++){
      float maxValue = -100;
      float minValue = 999999;
      for(int ihisto = 0; ihisto < mjj_qcd_fac.size(); ihisto++){
	if(mjj_qcd_fac.at(ihisto)->GetBinContent(iBin+1) > maxValue) maxValue = mjj_qcd_fac.at(ihisto)->GetBinContent(iBin+1);
	if(mjj_qcd_fac.at(ihisto)->GetBinContent(iBin+1) < minValue) minValue = mjj_qcd_fac.at(ihisto)->GetBinContent(iBin+1);
      }
      if(not symmetrize){
	mjj_qcd_fac_uncup->SetBinContent(iBin+1,maxValue/mjj_qcd_fac.at(0)->GetBinContent(iBin+1));
	mjj_qcd_fac_uncdw->SetBinContent(iBin+1,minValue/mjj_qcd_fac.at(0)->GetBinContent(iBin+1));    
      }
      else{
	mjj_qcd_fac_uncup->SetBinContent(iBin+1,1+(fabs(maxValue/mjj_qcd_fac.at(0)->GetBinContent(iBin+1)-1)+fabs(minValue/mjj_qcd_fac.at(0)->GetBinContent(iBin+1)-1))/2);
	mjj_qcd_fac_uncdw->SetBinContent(iBin+1,1-(fabs(maxValue/mjj_qcd_fac.at(0)->GetBinContent(iBin+1)-1)+fabs(minValue/mjj_qcd_fac.at(0)->GetBinContent(iBin+1)-1))/2);
      }
    }

    mjj_qcd_fac_uncup->Smooth();
    mjj_qcd_fac_uncdw->Smooth();

  }

  
  // Build the PDF uncertainty
  TH1F* bosonpt_pdf_uncup = (TH1F*) bosonpt_pdf.at(0)->Clone("bosonpt_pdf_uncup");
  bosonpt_pdf_uncup->Reset();
  TH1F* bosonpt_pdf_uncdw = (TH1F*) bosonpt_pdf.at(0)->Clone("bosonpt_pdf_uncdw");
  bosonpt_pdf_uncdw->Reset();
  TH1F* mjj_pdf_uncup = (TH1F*) mjj_pdf.at(0)->Clone("mjj_pdf_uncup");
  mjj_pdf_uncup->Reset();
  TH1F* mjj_pdf_uncdw = (TH1F*) mjj_pdf.at(0)->Clone("mjj_pdf_uncdw");
  mjj_pdf_uncdw->Reset();

  for(int iBin = 0; iBin < bosonpt_pdf_uncup->GetNbinsX(); iBin++){
    vector<double> binEntryBelow;
    vector<double> binEntryAbove;
    double nominalValue = bosonpt_qcd_fac.at(0)->GetBinContent(iBin+1);
    for(int ihisto = 0; ihisto < bosonpt_pdf.size(); ihisto++){
      if(bosonpt_pdf.at(ihisto)->GetBinContent(iBin+1) > nominalValue)
	binEntryAbove.push_back(bosonpt_pdf.at(ihisto)->GetBinContent(iBin+1));
      else
	binEntryBelow.push_back(bosonpt_pdf.at(ihisto)->GetBinContent(iBin+1));
    }
    sort(binEntryBelow.begin(),binEntryBelow.end());
    sort(binEntryAbove.begin(),binEntryAbove.end());
    if(not symmetrize){
      bosonpt_pdf_uncup->SetBinContent(iBin+1,2*binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue);
      bosonpt_pdf_uncdw->SetBinContent(iBin+1,2*binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue);
    }
    else{
      bosonpt_pdf_uncup->SetBinContent(iBin+1,1+2*(fabs(binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue-1)+fabs(binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue-1))/2);
      bosonpt_pdf_uncdw->SetBinContent(iBin+1,1-2*(fabs(binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue-1)+fabs(binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue-1))/2);
    }
  }
  bosonpt_pdf_uncup->Smooth();
  bosonpt_pdf_uncdw->Smooth();

  if(category == Category::VBF){

    for(int iBin = 0; iBin < mjj_pdf_uncup->GetNbinsX(); iBin++){
      vector<double> binEntryBelow;
      vector<double> binEntryAbove;
      double nominalValue = mjj_qcd_fac.at(0)->GetBinContent(iBin+1);
      for(int ihisto = 1; ihisto < mjj_pdf.size(); ihisto++){
	if(mjj_pdf.at(ihisto)->GetBinContent(iBin+1) > nominalValue)	  
	  binEntryAbove.push_back(mjj_pdf.at(ihisto)->GetBinContent(iBin+1));
	else
	  binEntryBelow.push_back(mjj_pdf.at(ihisto)->GetBinContent(iBin+1));
      }
      sort(binEntryBelow.begin(),binEntryBelow.end());
      sort(binEntryAbove.begin(),binEntryAbove.end());
      if(not symmetrize){
	if(mjj_pdf_uncup->GetBinCenter(iBin+1) < 2000){
	  mjj_pdf_uncup->SetBinContent(iBin+1,2*binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue);
	  mjj_pdf_uncdw->SetBinContent(iBin+1,2*binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue);
	}
	else{
	  mjj_pdf_uncup->SetBinContent(iBin+1,binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue);
	  mjj_pdf_uncdw->SetBinContent(iBin+1,binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue);
	}
      }
      else{
	if(mjj_pdf_uncup->GetBinCenter(iBin+1) < 2000){
	  mjj_pdf_uncup->SetBinContent(iBin+1,1+2*(fabs(binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue-1)+fabs(binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue-1))/2);
	  mjj_pdf_uncdw->SetBinContent(iBin+1,1-2*(fabs(binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue-1)+fabs(binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue-1))/2);
	}
	else{
	  if(sample == Sample::zjet){
	    mjj_pdf_uncup->SetBinContent(iBin+1,1+0.8*(fabs(binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue-1)+fabs(binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue-1))/2);
	    mjj_pdf_uncdw->SetBinContent(iBin+1,1-0.8*(fabs(binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue-1)+fabs(binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue-1))/2);
	  }
	  else if(sample == Sample::wjet){
	    if(mjj_pdf_uncup->GetBinCenter(iBin+1) < 2500){
	      mjj_pdf_uncup->SetBinContent(iBin+1,1+0.85*(fabs(binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue-1)+fabs(binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue-1))/2);
	      mjj_pdf_uncdw->SetBinContent(iBin+1,1-0.85*(fabs(binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue-1)+fabs(binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue-1))/2);
	    }
	    else{
	      mjj_pdf_uncup->SetBinContent(iBin+1,1+0.7*(fabs(binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue-1)+fabs(binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue-1))/2);
	      mjj_pdf_uncdw->SetBinContent(iBin+1,1-0.7*(fabs(binEntryAbove.at(int(binEntryAbove.size()*0.34)+1)/nominalValue-1)+fabs(binEntryBelow.at(int(binEntryBelow.size()*0.66)-1)/nominalValue-1))/2);

	    }
	  }
	}
      }
    }
    mjj_pdf_uncup->Smooth();
    mjj_pdf_uncdw->Smooth();
  }


  //////////
  // Plot
  //////////
  TCanvas* canvas =  new TCanvas("canvas","",600,625);
  canvas->cd();

  bosonpt_qcd_ren_uncup->SetLineColor(kRed);
  bosonpt_qcd_ren_uncup->SetLineWidth(2);
  bosonpt_qcd_ren_uncdw->SetLineColor(kRed);
  bosonpt_qcd_ren_uncdw->SetLineWidth(2);

  bosonpt_qcd_fac_uncup->SetLineColor(kBlue);
  bosonpt_qcd_fac_uncup->SetLineWidth(2);
  bosonpt_qcd_fac_uncdw->SetLineColor(kBlue);
  bosonpt_qcd_fac_uncdw->SetLineWidth(2);

  bosonpt_pdf_uncup->SetLineColor(kBlack);
  bosonpt_pdf_uncup->SetLineWidth(2);
  bosonpt_pdf_uncdw->SetLineColor(kBlack);
  bosonpt_pdf_uncdw->SetLineWidth(2);

  bosonpt_qcd_ren_uncup->GetXaxis()->SetTitle("Boson p_{T} [GeV]");
  bosonpt_qcd_ren_uncup->GetYaxis()->SetTitle("Uncertainty");
  bosonpt_qcd_ren_uncup->GetYaxis()->SetTitleOffset(1.20);
  bosonpt_qcd_ren_uncup->GetXaxis()->SetTitleOffset(1.1);

  
  bosonpt_qcd_ren_uncup->Draw("hist");
  bosonpt_qcd_ren_uncup->GetYaxis()->SetRangeUser(0.85,1.20);
  bosonpt_qcd_ren_uncdw->Draw("hist same");
  bosonpt_qcd_fac_uncup->Draw("hist same");
  bosonpt_qcd_fac_uncdw->Draw("hist same");
  bosonpt_pdf_uncup->Draw("hist same");
  bosonpt_pdf_uncdw->Draw("hist same");
  
  CMS_lumi(canvas,"35.9");

  TLegend leg (0.6,0.7,0.9,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry(bosonpt_qcd_ren_uncup,"Renormalization #mu_{R}","L");
  leg.AddEntry(bosonpt_qcd_fac_uncup,"Factorization #mu_{F}","L");
  leg.AddEntry(bosonpt_pdf_uncup,"PDF Variation","L");
  leg.Draw("same");

  if(sample == Sample::zjet){
    canvas->SaveAs((outputDIR+"/bosonpt_zjet_uncertainty.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/bosonpt_zjet_uncertainty.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::wjet){
    canvas->SaveAs((outputDIR+"/bosonpt_wjet_uncertainty.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/bosonpt_wjet_uncertainty.pdf").c_str(),"pdf");
  }
  else if(sample == Sample::gam){
    canvas->SaveAs((outputDIR+"/bosonpt_gam_uncertainty.png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/bosonpt_gam_uncertainty.pdf").c_str(),"pdf");
  }

  ///////////////////////////////////
  if(category == Category::VBF){

    mjj_qcd_ren_uncup->SetLineColor(kRed);
    mjj_qcd_ren_uncup->SetLineWidth(2);
    mjj_qcd_ren_uncdw->SetLineColor(kRed);
    mjj_qcd_ren_uncdw->SetLineWidth(2);
    
    mjj_qcd_fac_uncup->SetLineColor(kBlue);
    mjj_qcd_fac_uncup->SetLineWidth(2);
    mjj_qcd_fac_uncdw->SetLineColor(kBlue);
    mjj_qcd_fac_uncdw->SetLineWidth(2);

    mjj_pdf_uncup->SetLineColor(kBlack);
    mjj_pdf_uncup->SetLineWidth(2);
    mjj_pdf_uncdw->SetLineColor(kBlack);
    mjj_pdf_uncdw->SetLineWidth(2);
    
    mjj_qcd_ren_uncup->GetXaxis()->SetTitle("M_{jj} [GeV]");
    mjj_qcd_ren_uncup->GetYaxis()->SetTitle("Uncertainty");
    mjj_qcd_ren_uncup->GetYaxis()->SetTitleOffset(1.20);
    mjj_qcd_ren_uncup->GetXaxis()->SetTitleOffset(1.1);
    mjj_qcd_ren_uncup->GetXaxis()->SetNdivisions(505);
    
    
    mjj_qcd_ren_uncup->Draw("hist");
    mjj_qcd_ren_uncup->GetYaxis()->SetRangeUser(0.75,1.25);
    mjj_qcd_ren_uncdw->Draw("hist same");
    mjj_qcd_fac_uncup->Draw("hist same");
    mjj_qcd_fac_uncdw->Draw("hist same");
    mjj_pdf_uncup->Draw("hist same");
    mjj_pdf_uncdw->Draw("hist same");
    
    CMS_lumi(canvas,"35.9");
    
    TLegend leg (0.6,0.75,0.9,0.9);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    
    leg.AddEntry(mjj_qcd_ren_uncup,"Renormalization #mu_{R}","L");
    leg.AddEntry(mjj_qcd_fac_uncup,"Factorization #mu_{F}","L");
    leg.AddEntry(mjj_pdf_uncup,"PDF Variation","L");
    leg.Draw("same");
    
    if(sample == Sample::zjet){
      canvas->SaveAs((outputDIR+"/mjj_zjet_uncertainty.png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/mjj_zjet_uncertainty.pdf").c_str(),"pdf");
    }
    else if(sample == Sample::wjet){
      canvas->SaveAs((outputDIR+"/mjj_wjet_uncertainty.png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/mjj_wjet_uncertainty.pdf").c_str(),"pdf");
    }
    else if(sample == Sample::gam){
      canvas->SaveAs((outputDIR+"/mjj_gam_uncertainty.png").c_str(),"png");
      canvas->SaveAs((outputDIR+"/mjj_gam_uncertainty.pdf").c_str(),"pdf");
    }
  }

  
  TFile* outputFile = new TFile((outputDIR+"/ouputSingleProcess_"+postfix+".root").c_str(),"RECREATE");
  outputFile->cd();
  mjj_qcd_ren_uncup->Write();
  mjj_qcd_ren_uncdw->Write();
  mjj_qcd_fac_uncup->Write();
  mjj_qcd_fac_uncdw->Write();
  mjj_pdf_uncup->Write();
  mjj_pdf_uncdw->Write();
  bosonpt_qcd_ren_uncup->Write();
  bosonpt_qcd_ren_uncdw->Write();
  bosonpt_qcd_fac_uncup->Write();
  bosonpt_qcd_fac_uncdw->Write();
  bosonpt_pdf_uncup->Write();
  bosonpt_pdf_uncdw->Write();
  outputFile->Close();
}
