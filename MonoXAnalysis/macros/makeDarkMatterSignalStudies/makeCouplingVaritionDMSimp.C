#include "../CMS_lumi.h"

static float minBosonPt = 200;

enum class Sample {vector,axial,scalar,pseudoscalar};

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
  vector<float> gsm;
  vector<float> gdm;
  vector<TH1F*> bosonpt; // to embed coupling variations
};


vector<float> gSMValue_V = {0.25,0.25,0.25,0.1,0.1,0.1,0.50,0.50,0.50};
vector<float> gDMValue_V = {1.00,0.50,2.00,1.0,0.5,2.0,1.00,0.50,2.00};

vector<int>   color_V    = {1,1,1,2,2,2,4,4,4};
vector<int>   style_V    = {20,24,25,20,24,25,20,24,25};

vector<float> gSMValue_S = {1.0,0.5,2.0,3.5};
vector<float> gDMValue_S = {1.0,1.0,1.0,1.0};

vector<int>   color_S    = {1,2,4,6};
vector<int>   style_S    = {20,20,20,20};


vector<float> metBin    = {200.,230.,260,290,320,350,390,430,470,510,550,590,640,690,740,790,840,900,960,1020,1090,1160,1250};

void makeCouplingVaritionDMSimp(string inputDIR,
				string outputDIR,
				Sample sample,
				string postfix,
				vector<int> medMass,
				vector<int> dmMass){


  float lumi = 36.6;
  
  gROOT->SetBatch(kTRUE);
  setTDRStyle();

  system(("mkdir -p "+outputDIR).c_str());

  vector<TTree*> chain_monojet;
  vector<TFile*> file_monojet;

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
	  file_monojet.push_back(TFile::Open((inputDIR+"/"+line).c_str()));
	  if(file_monojet.back())
	    chain_monojet.push_back((TTree*) file_monojet.back()->Get("gentree/tree"));
	}
      }
    }
  }


  cout<<"Sum of weights monojet "<<endl;
  vector<double> sumwgt_monojet;
  for(auto tree : chain_monojet){

    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");
    TTreeReaderValue<float> xsec (reader,"xsec");
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");
    TTreeReaderValue<float> dmpt (reader,"dmpt");
    TTreeReaderValue<vector<float> > gDMV (reader,"gDMV");
    TTreeReaderValue<vector<float> > gDMA (reader,"gDMA");
    TTreeReaderValue<vector<float> > gV (reader,"gV");
    TTreeReaderValue<vector<float> > gA (reader,"gA");
    TTreeReaderValue<vector<float> >gTheta (reader,"gTheta");

    sumwgt_monojet.push_back(0);
    while(reader.Next()){

      if(sample == Sample::vector and (gDMV->size() == 0 or gV->size() == 0)) continue;
      if(sample == Sample::scalar and (gDMV->size() == 0 or gV->size() == 0)) continue;
      if(sample == Sample::axial and (gDMV->size() == 0 or gV->size() == 0)) continue;
      if(sample == Sample::pseudoscalar and (gDMA->size() == 0 or gA->size() == 0)) continue;

      sumwgt_monojet.back() += *wgt;
    }
  }


  system("rm list.temp");

  // each loop is 1 mass point
  cout<<"Loop on monojet chain "<<endl; 
  vector<massPoint> bosonPtSpectrum;
  int ifile = 0;
  for(auto tree : chain_monojet){
    ifile++;
    TTreeReader reader (tree);
    TTreeReaderValue<float> wgt (reader,"wgt");
    TTreeReaderValue<float> xsec (reader,"xsec");
    TTreeReaderValue<float> sampledmM (reader,"sampledmM");
    TTreeReaderValue<float> samplemedM (reader,"samplemedM");
    TTreeReaderValue<float> dmpt (reader,"dmpt");
    TTreeReaderValue<vector<float> > couplingwgt (reader,"couplingwgt");
    TTreeReaderValue<vector<float> > gDMV (reader,"gDMV");
    TTreeReaderValue<vector<float> > gDMA (reader,"gDMA");
    TTreeReaderValue<vector<float> > gV (reader,"gV");
    TTreeReaderValue<vector<float> > gA (reader,"gA");
    TTreeReaderValue<vector<float> >gTheta (reader,"gTheta");

    long int iEvent = 0;
    while(reader.Next()){
      iEvent++;
      
      if(iEvent == 1)
	bosonPtSpectrum.push_back(massPoint(int(*samplemedM),int(*sampledmM)));


      if(sample == Sample::vector and (gDMV->size() == 0 or gV->size() == 0)) continue;
      if(sample == Sample::scalar and (gDMV->size() == 0 or gV->size() == 0)) continue;
      if(sample == Sample::axial and (gDMV->size() == 0 or gV->size() == 0)) continue;
      if(sample == Sample::pseudoscalar and (gDMA->size() == 0 or gA->size() == 0)) continue;
      
      if(bosonPtSpectrum.back().bosonpt.size() == 0){
	if(sample == Sample::vector or sample == Sample::scalar){
	  if(gDMV->size() != gV->size()){
	    cerr<<"problem with vector coupling dimension --> to be checked "<<endl;
	    return;
	  }
	  for(size_t iwgt = 0; iwgt < gDMV->size(); iwgt++){
	    if(gV->at(iwgt) == 0 and gDMV->at(iwgt) == 0) break;
	    bosonPtSpectrum.back().bosonpt.push_back(new TH1F(Form("bosonpt_med_%d_dm_%d_gsm_%.2f_gdm_%.2f",int(*samplemedM),int(*sampledmM),gV->at(iwgt),gDMV->at(iwgt)),"",metBin.size()-1,&metBin[0]));
	    bosonPtSpectrum.back().bosonpt.back()->Sumw2();
	    bosonPtSpectrum.back().gsm.push_back(gV->at(iwgt));
	    bosonPtSpectrum.back().gdm.push_back(gDMV->at(iwgt));	    
	  }
	  
	  if(sample == Sample::vector){
	    bosonPtSpectrum.back().bosonpt.push_back(new TH1F(Form("bosonpt_med_%d_dm_%d_nominal",int(*samplemedM),int(*sampledmM)),"",metBin.size()-1,&metBin[0]));
	    bosonPtSpectrum.back().bosonpt.back()->Sumw2();
	    bosonPtSpectrum.back().gsm.push_back(0.25);
	    bosonPtSpectrum.back().gdm.push_back(1.0);	    
	  }
	  if(sample == Sample::scalar){
	    bosonPtSpectrum.back().bosonpt.push_back(new TH1F(Form("bosonpt_med_%d_dm_%d_nominal",int(*samplemedM),int(*sampledmM)),"",metBin.size()-1,&metBin[0]));
            bosonPtSpectrum.back().bosonpt.back()->Sumw2();
            bosonPtSpectrum.back().gsm.push_back(1.0);
            bosonPtSpectrum.back().gdm.push_back(1.0);
	  }
	}
	else if(sample == Sample::axial or sample == Sample::pseudoscalar){
	  if(gDMA->size() != gA->size()){
	    cerr<<"problem with vector coupling dimension --> to be checked "<<endl;
	    return;
	  }
	  for(size_t iwgt = 0; iwgt < gDMA->size(); iwgt++){
	    if(gA->at(iwgt) == 0 and gDMA->at(iwgt) == 0) break;
	    bosonPtSpectrum.back().bosonpt.push_back(new TH1F(Form("bosonpt_med_%d_dm_%d_gsm_%.2f_gdm_%.2f",int(*samplemedM),int(*sampledmM),gA->at(iwgt),gDMA->at(iwgt)),"",metBin.size()-1,&metBin[0]));
	    bosonPtSpectrum.back().bosonpt.back()->Sumw2();
	    bosonPtSpectrum.back().gsm.push_back(gA->at(iwgt));
	    bosonPtSpectrum.back().gdm.push_back(gDMA->at(iwgt));
	  }
	  if(sample == Sample::axial){
	    bosonPtSpectrum.back().bosonpt.push_back(new TH1F(Form("bosonpt_med_%d_dm_%d_nominal",int(*samplemedM),int(*sampledmM)),"",metBin.size()-1,&metBin[0]));
	    bosonPtSpectrum.back().bosonpt.back()->Sumw2();
	    bosonPtSpectrum.back().gsm.push_back(0.25);
	    bosonPtSpectrum.back().gdm.push_back(1.0);
	  }
	  if(sample == Sample::pseudoscalar){
	    bosonPtSpectrum.back().bosonpt.push_back(new TH1F(Form("bosonpt_med_%d_dm_%d_nominal",int(*samplemedM),int(*sampledmM)),"",metBin.size()-1,&metBin[0]));
	    bosonPtSpectrum.back().bosonpt.back()->Sumw2();
	    bosonPtSpectrum.back().gsm.push_back(1.0);
	    bosonPtSpectrum.back().gdm.push_back(1.0);
	  }
	}	
      }
      
      // appply boson ot cut
      if(*dmpt < minBosonPt) continue;
      
      double wgt_central = *wgt;
      if(sample == Sample::scalar or sample == Sample::pseudoscalar){ // find weight in a more robust way
	for(size_t iwgt = 0; iwgt < gDMV->size(); iwgt++){
	  if(sample == Sample::scalar and gV->at(iwgt) == float(1.0) and gDMV->at(iwgt) == float(1.0))
	    wgt_central = couplingwgt->at(iwgt);
	  else if(sample == Sample::pseudoscalar and gA->at(iwgt) == float(1.0) and gDMA->at(iwgt) == float(1.0))
	    wgt_central = couplingwgt->at(iwgt);
	}	  
      }

      if(sample == Sample::vector or sample == Sample::scalar){
	bosonPtSpectrum.back().bosonpt.back()->Fill(*dmpt,*xsec*(*wgt)*lumi/sumwgt_monojet.at(ifile-1));
	for(size_t iwgt = 0; iwgt < gDMV->size(); iwgt++){
	  if(gV->at(iwgt) == 0 and gDMV->at(iwgt) == 0) break;
	  bosonPtSpectrum.back().bosonpt.at(iwgt)->Fill(*dmpt,*xsec*(*wgt)*lumi*couplingwgt->at(iwgt)/(wgt_central)*1./sumwgt_monojet.at(ifile-1));
	}
      }
      if(sample == Sample::axial or sample == Sample::pseudoscalar){
	bosonPtSpectrum.back().bosonpt.back()->Fill(*dmpt,*xsec*(*wgt)*lumi/sumwgt_monojet.at(ifile-1));
	for(size_t iwgt = 0; iwgt < gDMA->size(); iwgt++){
	  if(gA->at(iwgt) == 0 and gDMA->at(iwgt) == 0) break;
	  bosonPtSpectrum.back().bosonpt.at(iwgt)->Fill(*dmpt,*xsec*(*wgt)*lumi*couplingwgt->at(iwgt)/(wgt_central)*1./sumwgt_monojet.at(ifile-1));
	}
      }      
    }
  }

  TFile* output = new TFile((outputDIR+"/templates_coupling_"+postfix+".root").c_str(),"RECREATE");
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

  for(auto element : bosonPtSpectrum){
    bool isfound = false; 
    int icolor = 0;

    TH1F* frame = new TH1F(Form("frame_med_%d_dm_%d",element.medMass_,element.dmMass_),"",metBin.size()-1,&metBin[0]);
    frame->GetXaxis()->SetTitle("Mediator p_{T} [GeV]");
    frame->GetYaxis()->SetTitle("Events");
    frame->Draw();
    CMS_lumi(canvas,"36.6");

    TLegend leg (0.5,0.6,0.92,0.92);
    leg.SetFillColor(0);
    leg.SetFillStyle(0);
    leg.SetBorderSize(0);
    leg.AddEntry((TObject*)(0),Form("m_{med} = %d, m_{dm} = %d",element.medMass_,element.dmMass_),"");

    double minimum = 999999;
    double maximum = -1;
    unsigned int pos = 0;

    vector<float> gSM;
    vector<float> gDM;

    if(sample == Sample::vector or sample == Sample::axial){
      gSM = gSMValue_V;
      gDM = gDMValue_V;
    }
    else if(sample == Sample::scalar or sample == Sample::pseudoscalar){
      gSM = gSMValue_S;
      gDM = gDMValue_S;
    }

    for(size_t ig = 0; ig < gSM.size(); ig++){
      if(std::find(element.gsm.begin(),element.gsm.end(),gSM.at(ig)) != element.gsm.end()){
	if(std::find(element.gdm.begin(),element.gdm.end(),gDM.at(ig)) != element.gdm.end()){
	  for(size_t val = 0; val < element.gsm.size(); val++){
	    if(element.gsm.at(val) == gSM.at(ig) and element.gdm.at(val) == gDM.at(ig))
	      pos = val;
	  }
	  
	  TH1* histo = element.bosonpt.at(pos);	  
	  if(histo->GetMinimum() < minimum)
	    minimum = histo->GetMinimum();
	  if(histo->GetMaximum() > maximum)
	    maximum = histo->GetMaximum();
	  
	  if(sample == Sample::vector or sample == Sample::axial){
	    histo->SetLineWidth(2);
	    histo->SetMarkerStyle(style_V.at(icolor));
	    histo->SetLineColor(color_V.at(icolor));
	    histo->SetMarkerColor(color_V.at(icolor));
	    histo->Draw("hist same");
	    histo->Draw("P same");
	    cout<<"histo name "<<histo->GetName()<<" integral "<<histo->Integral()<<endl;
	  }
	  else if(sample == Sample::scalar or sample == Sample::pseudoscalar){
	    histo->SetLineWidth(2);
	    histo->SetMarkerStyle(style_S.at(icolor));
	    histo->SetLineColor(color_S.at(icolor));
	    histo->SetMarkerColor(color_S.at(icolor));
	    histo->Draw("hist same");
	    histo->Draw("P same");
	    cout<<"histo name "<<histo->GetName()<<" integral "<<histo->Integral()<<endl;
	  }
	  leg.AddEntry(histo,Form("g_{SM} = %.2f, g_{DM} = %.2f",gSM.at(ig),gDM.at(ig)),"PL");
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
