#include "../CMS_lumi.h"

static float minBosonPt = 200;

class massPoint {
public:
  massPoint(){};
  ~massPoint(){};
  massPoint(int medMass, int dmMass){
    medMass_ = medMass;
    dmMass_  = dmMass;
  }
  
  int medMass_;
  int dmMass_;
  vector<float> gdm;
  vector<float> gtheta;
  vector<TH1F*> bosonpt; // to embed coupling variations
};

/// for coupling variations
vector<float> gDMValue    = {1.0,1.0,1.0,0.4,0.4,0.4,0.1,0.1,0.1};
vector<float> gThetaValue = {0.4,0.2,0.1,0.4,0.2,0.1,0.4,0.2,0.1};

vector<int>   color_S    = {1,1,1,2,2,2,4,4,4};
vector<int>   style_S    = {20,24,25,20,24,25,20,24,25};

vector<float> metBin    = {200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1250};

void makeSMMDistributions(string inputDIR,
			  string outputDIR,
			  string postfix,
			  vector<int> medMass,
			  vector<int> dmMass){


  float lumi = 36.6;
  
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputDIR).c_str());

  vector<TTree*> chain_smm;
  vector<TFile*> file_smm;

  if(medMass.size() != dmMass.size()){
    cerr<<"DM mediator mass point size is different form DM mass one --> please check"<<endl;
    return;
  }
  
  
  if(inputDIR != ""){
    for(size_t imass = 0; imass < medMass.size(); imass++){      
      system(("ls "+inputDIR+" | grep Mphi-"+to_string(medMass.at(imass))+"_ | grep Mchi-"+to_string(dmMass.at(imass))+"_gSM | grep root > list.temp ").c_str());
      ifstream file ("list.temp");
      if(file.is_open()){
	string line;
	while(!file.eof()){
	  getline(file,line);
	  if(line == "") continue;
	  file_smm.push_back(TFile::Open((inputDIR+"/"+line).c_str()));
	  if(file_smm.back())
	    chain_smm.push_back((TTree*) file_smm.back()->Get("gentree/tree"));
	}
      }
    }
  }

  system("rm list.temp");

  cout<<"Sum of weights smm "<<endl;
  vector<double> sumwgt_smm;
  for(auto tree : chain_smm){

    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");
    TTreeReaderValue<float> xsec (reader,"xsec");
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");
    TTreeReaderValue<float> dmpt (reader,"dmpt");
    TTreeReaderValue<vector<float> > gDMV (reader,"gDMV");
    TTreeReaderValue<vector<float> >gTheta (reader,"gTheta");

    sumwgt_smm.push_back(0);
    while(reader.Next()){
      if(gDMV->size() == 0 or gTheta->size() == 0) continue;
      sumwgt_smm.back() += *wgt;
    }
  }


  // each loop is 1 mass point
  cout<<"Loop on SMM chain "<<endl; 
  vector<massPoint> bosonPtSpectrum;
  int ifile = 0;
  for(auto tree : chain_smm){
    ifile++;
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");
    TTreeReaderValue<float> xsec (reader,"xsec");
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");
    TTreeReaderValue<float> dmpt (reader,"dmpt");
    TTreeReaderValue<vector<float> > couplingwgt (reader,"couplingwgt");
    TTreeReaderValue<vector<float> > gDMV (reader,"gDMV");
    TTreeReaderValue<vector<float> > gV   (reader,"gV");
    TTreeReaderValue<vector<float> >gTheta (reader,"gTheta");

    long int iEvent = 0;
    while(reader.Next()){

      if(gDMV->size() == 0 or gTheta->size() == 0) continue;
      iEvent++;
            
      if(iEvent == 1)
	bosonPtSpectrum.push_back(massPoint(int(*samplemedM),int(*sampledmM)));
      
      if(bosonPtSpectrum.back().bosonpt.size() == 0){
	if(gDMV->size() != gTheta->size()){
	  cerr<<"problem with vector coupling dimension --> to be checked "<<endl;
	  return;
	}
	
	for(size_t iwgt = 0; iwgt < gDMV->size(); iwgt++){
	  if(gTheta->at(iwgt) == 0 and gDMV->at(iwgt) == 0) break;

	  bosonPtSpectrum.back().bosonpt.push_back(new TH1F(Form("bosonpt_med_%d_dm_%d_gtheta_%.3f_gdm_%.3f",int(*samplemedM),int(*sampledmM),gTheta->at(iwgt),gDMV->at(iwgt)),"",metBin.size()-1,&metBin[0]));
	  bosonPtSpectrum.back().bosonpt.back()->Sumw2();
	  bosonPtSpectrum.back().gtheta.push_back(gTheta->at(iwgt));
	  bosonPtSpectrum.back().gdm.push_back(gDMV->at(iwgt));	    
	}

	bosonPtSpectrum.back().bosonpt.push_back(new TH1F(Form("bosonpt_med_%d_dm_%d_nominal",int(*samplemedM),int(*sampledmM)),"",metBin.size()-1,&metBin[0]));
	bosonPtSpectrum.back().bosonpt.back()->Sumw2();
	bosonPtSpectrum.back().gtheta.push_back(0.4);
	bosonPtSpectrum.back().gdm.push_back(1.0);	    
      }	
      
      
      // appply boson ot cut
      if(*dmpt < minBosonPt) continue;
      
      double wgt_central = *wgt;
      for(size_t iwgt = 0; iwgt < gDMV->size(); iwgt++){
	if(gDMV->at(iwgt) == float(1.0) and gTheta->at(iwgt) == float(0.4))
	  wgt_central = couplingwgt->at(iwgt);	
      }
      bosonPtSpectrum.back().bosonpt.back()->Fill(*dmpt,*xsec*(*wgt)*lumi/sumwgt_smm.at(ifile-1));
      for(size_t iwgt = 0; iwgt < gDMV->size(); iwgt++){
	if(gTheta->at(iwgt) == 0 and gDMV->at(iwgt) == 0) break;
	bosonPtSpectrum.back().bosonpt.at(iwgt)->Fill(*dmpt,*xsec*(*wgt)*lumi*couplingwgt->at(iwgt)/(wgt_central)*1./sumwgt_smm.at(ifile-1));
      }
    }
  }
  
  TFile* output = new TFile((outputDIR+"/templates_mass_coupling_"+postfix+".root").c_str(),"RECREATE");
  output->cd();
  
  for(auto element : bosonPtSpectrum){
    output->mkdir(Form("med_%d_dm_%d",element.medMass_,element.dmMass_));
    output->cd(Form("med_%d_dm_%d",element.medMass_,element.dmMass_));
    for(auto hist : element.bosonpt){
      hist->Write();
    }
    output->cd();
  }

  // make interesting plots
  TCanvas* canvas = new TCanvas("canvas","",600,625);
  canvas->cd();

  TH1F* frame = new TH1F("frame","",metBin.size()-1,&metBin[0]);
  frame->GetXaxis()->SetTitle("Mediator p_{T} [GeV]");
  frame->GetYaxis()->SetTitle("Events");                                                                                                                                                             
  frame->Draw();                                                                                                                                                                                     
  CMS_lumi(canvas,"36.6");                                                                                                                                                                           

  TLegend leg (0.5,0.6,0.92,0.9);
  leg.SetFillColor(0);
  leg.SetFillStyle(0);
  leg.SetBorderSize(0);
  leg.AddEntry((TObject*)(0),"sin^{2}#theta = 0.4, y_{DM} = 1.0","");
  leg.AddEntry((TObject*)(0),"","");
  int icolor = 1;
  double min = 100000000;
  double max = 0;

  for(auto element : bosonPtSpectrum){
    if(element.bosonpt.back()->GetMinimum() < min)
      min = element.bosonpt.back()->GetMinimum();
    if(element.bosonpt.back()->GetMaximum() > max)
      max = element.bosonpt.back()->GetMaximum();
    if(icolor == 3 or icolor == 5) icolor++;
    element.bosonpt.back()->SetLineColor(icolor);
    element.bosonpt.back()->SetLineWidth(2);
    element.bosonpt.back()->SetMarkerColor(icolor);
    if(icolor == 1)
      element.bosonpt.back()->SetMarkerStyle(20);
    else
      element.bosonpt.back()->SetMarkerStyle(24);

    leg.AddEntry(element.bosonpt.back(),Form("m_{med} = %d GeV, m_{dm} = %d GeV",element.medMass_,element.dmMass_),"L");
    element.bosonpt.back()->Draw("hist same");
    element.bosonpt.back()->Draw("EP same");    
    icolor++;
  }

  frame->GetYaxis()->SetRangeUser(min*0.5,max*100);
  
  leg.Draw("same");
  canvas->SetLogy();
  canvas->SaveAs((outputDIR+"/yields_mmed_mdm_"+postfix+".png").c_str(),"png");
  canvas->SaveAs((outputDIR+"/yields_mmed_mdm_"+postfix+".pdf").c_str(),"pdf");

  // coupling plot
  
  canvas->cd();
  for(auto element : bosonPtSpectrum){
    bool isfound = false; 
    icolor = 0;
    
    TH1F* frame = new TH1F(Form("frame_med_%d_dm_%d",element.medMass_,element.dmMass_),"",metBin.size()-1,&metBin[0]);
    frame->GetXaxis()->SetTitle("Mediator p_{T} [GeV]");
    frame->GetYaxis()->SetTitle("Events");
    frame->Draw();
    CMS_lumi(canvas,"36.6");
    
    TLegend leg (0.5,0.6,0.92,0.92);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry((TObject*)(0),Form("m_{med} = %d GeV, m_{dm} = %d GeV",element.medMass_,element.dmMass_),"");
    leg.AddEntry((TObject*)(0),"","");
    
    double minimum = 999999;
    double maximum = -1;
    unsigned int pos = 0;
        
    for(size_t ig = 0; ig < gDMValue.size(); ig++){
      if(std::find(element.gdm.begin(),element.gdm.end(),gDMValue.at(ig)) != element.gdm.end()){
	if(std::find(element.gtheta.begin(),element.gtheta.end(),gThetaValue.at(ig)) != element.gtheta.end()){
	  for(size_t val = 0; val < element.gtheta.size(); val++){
	    if(element.gtheta.at(val) == gThetaValue.at(ig) and element.gdm.at(val) == gDMValue.at(ig))
	      pos = val;
	  }
	  
	  TH1* histo = element.bosonpt.at(pos);	  
	  if(histo->GetMinimum() < minimum)
	    minimum = histo->GetMinimum();
	  if(histo->GetMaximum() > maximum)
	    maximum = histo->GetMaximum();
	  
	  histo->SetLineWidth(2);
	  histo->SetMarkerStyle(style_S.at(icolor));
	  histo->SetLineColor(color_S.at(icolor));
	  histo->SetMarkerColor(color_S.at(icolor));
	  histo->Draw("hist same");
	  histo->Draw("P same");
	  cout<<"histo name "<<histo->GetName()<<" integral "<<histo->Integral()<<endl;
	  leg.AddEntry(histo,Form("sin^{2}#theta_{SM} = %.2f, y_{DM} = %.2f",gThetaValue.at(ig),gDMValue.at(ig)),"PL");
	  icolor++;	  	  
	}	
      }
    }

    if(minimum != 0)
      frame->GetYaxis()->SetRangeUser(minimum*0.1,maximum*100);
    else
      frame->GetYaxis()->SetRangeUser(0.0001,maximum*100);
    
    leg.Draw("same");      
    canvas->SetLogy();
    canvas->SaveAs((outputDIR+"/yields_vs_coupling_med_"+to_string(element.medMass_)+"_dm_"+to_string(element.dmMass_)+"_"+postfix+".png").c_str(),"png");
    canvas->SaveAs((outputDIR+"/yields_vs_coupling_med_"+to_string(element.medMass_)+"_dm_"+to_string(element.dmMass_)+"_"+postfix+".pdf").c_str(),"pdf");
    
  }
  
  output->Close();
}
