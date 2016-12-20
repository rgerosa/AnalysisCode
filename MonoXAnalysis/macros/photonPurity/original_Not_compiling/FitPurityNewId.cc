gSystem->Load("libRoofit"); 
using namespace RooFit;
using namespace RooStats ;
//#include "./CMS_lumi.C"
static const Int_t NCAT = 4;  
TString mainCut;
TString mainCutQCD;
TString sieieCut;
TString sieieCutSideBand;
double PHIsoAbsCut;
TString PHIsoRelCut;
TString RescaledVar;
TString RescaledRNDVar;
TString Tregion;
TString Twp;
TString Tptcut;
Float_t minFit= 0.;
Float_t maxFit = 12;
//Int_t nbins =30;//40;//80
Int_t nbins =24;//40;//80
//Int_t nbins =40;//40;//80

//TString treePath = "/u2/mciprian/TREES_PHOTON_PURITY/merged";
void AddData(RooWorkspace*);
void runfits(RooWorkspace* w,std::string region, std::string phid, std::string ptcut);
RooFitResult*  ModelFit(RooWorkspace*);

double deltaPhi(double phi1,double phi2)
{
  double result = phi1 - phi2;
  while (result > TMath::Pi()) result -= 2*TMath::Pi();
  while (result <= -TMath::Pi()) result += 2*TMath::Pi();
  return result;
}

RooArgSet* defineVariables() {
  RooRealVar* ph_eta_sig_livia = new RooRealVar("pPurityheta", "eta(g)",-5,5,"");
  RooRealVar* ph_phi_sig_livia = new RooRealVar("phPurityphi", "pt(g)",-5,5,"");
  RooRealVar* ph_pt_sig_livia = new RooRealVar("phPuritypt", "pt(g)",0,1000,""); 
  RooRealVar* ph_elveto_sig_livia = new RooRealVar("phPurityElectronVeto", "veto",0,1000,""); 
  RooRealVar* ph_iso_sig_livia  = new RooRealVar("phPurityPHiso", "PFIsophoton",minFit,maxFit,"GeV");//PH iso
  RooRealVar* ph_rndiso04_sig_livia  = new RooRealVar("phPurityRND04PHiso", "RNDPFIsophoton",minFit,maxFit,"GeV");//RND PH iso
  RooRealVar* ph_rndiso08_sig_livia  = new RooRealVar("phPurityRND08PHiso", "RNDPFIsophoton",minFit,maxFit,"GeV");//RND PH iso
  RooRealVar* ph_EA_sig_livia  = new RooRealVar("phPurityEA", "EAs",0,10,"GeV");//RND PH iso
  RooRealVar* rho  = new RooRealVar("rho", "rho",0, 100,"GeV");//RND PH iso
  RooRealVar* ph_chiso_sig_livia  = new RooRealVar("phPurityCHiso", "PFCHIsophoton",minFit,maxFit,"GeV");// CH iso
  RooRealVar* ph_nhiso_sig_livia  = new RooRealVar("phPurityNHiso", "PFNHIsophoton",minFit,maxFit,"GeV");// NH iso
  RooRealVar* ph_hoe_sig_livia  = new RooRealVar("phPurityhoe", "HoEphoton",minFit,maxFit,"GeV");//HoE
  RooRealVar* ph_sieie_sig_livia  = new RooRealVar("phPuritysieie", "sieie",minFit,maxFit,"");//Sieie

  RooRealVar* weight = new RooRealVar("weight","weighting",0,1000,"");

  RooRealVar* nphotons = new RooRealVar("nphotons","nphotons",0,10,"");
  RooRealVar* nphotonsPurity = new RooRealVar("nphotonsPurity","nphotonsPurity",0,10,"");
  RooRealVar* hltphoton165 = new RooRealVar("hltphoton165","hltphoton165",0,1,"");
  RooRealVar* hltphoton175 = new RooRealVar("hltphoton175","hltphoton175",0,1,"");
  RooRealVar* signaljetpt = new RooRealVar("signaljetpt","signaljetpt",0,2000.,"");
  RooRealVar* signaljeteta = new RooRealVar("signaljeteta","signaljeteta",-5.,5.,"");
  RooRealVar* njets = new RooRealVar("njets","njets",0,20,"");
  /* RooRealVar* nvtx = new RooRealVar("nvtx", "nvtx", 0, 50, "");
  */
  RooRealVar* wzid   = new RooRealVar("wzid", "wzid", 0, 1000, "");

  RooArgSet* ntplVars = new RooArgSet(*nphotons, *hltphoton165, *hltphoton175, *ph_elveto_sig_livia, *rho, *ph_EA_sig_livia);
  RooArgSet* ntplVarsSpectators = new RooArgSet(*signaljetpt,*signaljeteta,*njets, *weight,*ph_nhiso_sig_livia,*ph_hoe_sig_livia, *ph_rndiso08_sig_livia);
  RooArgSet* ntplVarsSpectatorsPurity = new RooArgSet( *ph_eta_sig_livia, *ph_phi_sig_livia, *ph_pt_sig_livia, *ph_iso_sig_livia, *ph_rndiso04_sig_livia,*ph_chiso_sig_livia, *ph_sieie_sig_livia, *nphotonsPurity);
  RooArgSet* ntplVarsGen = new RooArgSet(*wzid);
  ntplVars->add(*ntplVarsSpectators);
  ntplVars->add(*ntplVarsSpectatorsPurity);
  ntplVars->add(*ntplVarsGen);
  return ntplVars;
}

RooArgSet* defineVariables_match() {
  RooRealVar* ph_eta_sig_livia = new RooRealVar("pPurityheta", "eta(g)",-5,5,"");
  RooRealVar* ph_phi_sig_livia = new RooRealVar("phPurityphi", "pt(g)",-5,5,"");
  RooRealVar* ph_pt_sig_livia = new RooRealVar("phPuritypt", "pt(g)",0,1000,""); 
  RooRealVar* ph_iso_sig_livia  = new RooRealVar("phPurityPHiso", "PFIsophoton",minFit,maxFit,"GeV");
  RooRealVar* ph_elveto_sig_livia = new RooRealVar("phPurityElectronVeto", "veto",0,1000,""); 
 
  RooRealVar* ph_rndiso04_sig_livia  = new RooRealVar("phPurityRND04PHiso", "RNDPFIsophoton",minFit,maxFit,"GeV");//RND PH iso
  RooRealVar* ph_rndiso08_sig_livia  = new RooRealVar("phPurityRND08PHiso", "RNDPFIsophoton",minFit,maxFit,"GeV");//RND PH iso
  RooRealVar* ph_EA_sig_livia  = new RooRealVar("phPurityEA", "EAs",0,10,"GeV");//RND PH iso
  RooRealVar* rho  = new RooRealVar("rho", "rho",0, 100,"GeV");//RND PH iso
  RooRealVar* ph_chiso_sig_livia  = new RooRealVar("phPurityCHiso", "PFCHIsophoton",minFit,maxFit,"GeV");
  RooRealVar* ph_nhiso_sig_livia  = new RooRealVar("phPurityNHiso", "PFNHIsophoton",minFit,maxFit,"GeV");
  RooRealVar* ph_hoe_sig_livia  = new RooRealVar("phPurityhoe", "HoEphoton",minFit,maxFit,"GeV");

  RooRealVar* ph_sieie_sig_livia  = new RooRealVar("phPuritysieie", "sieie",minFit,maxFit,"");

  RooRealVar* weight = new RooRealVar("weight","weighting",0,1000,"");

  RooRealVar* nphotons = new RooRealVar("nphotons","nphotons",0,10,"");
  RooRealVar* nphotonsPurity = new RooRealVar("nphotonsPurity","nphotonsPurity",0,10,"");
  RooRealVar* hltphoton165 = new RooRealVar("hltphoton165","hltphoton165",0,1,"");
  RooRealVar* hltphoton175 = new RooRealVar("hltphoton175","hltphoton175",0,1,"");
  RooRealVar* signaljetpt = new RooRealVar("signaljetpt","signaljetpt",0,2000.,"");
  RooRealVar* signaljeteta = new RooRealVar("signaljeteta","signaljeteta",-5.,5.,"");
  RooRealVar* njets = new RooRealVar("njets","njets",0,20,"");
  /* RooRealVar* nvtx = new RooRealVar("nvtx", "nvtx", 0, 50, "");
  */
  RooRealVar* wzid   = new RooRealVar("wzid", "wzid", 0, 1000, "");
  RooRealVar* ismatch   = new RooRealVar("ismatch", "ismatch", 0, 1000, "");

  RooArgSet* ntplVars = new RooArgSet(*nphotons, *hltphoton165, *hltphoton175, *ph_elveto_sig_livia, *rho, *ph_EA_sig_livia);
  RooArgSet* ntplVarsSpectators = new RooArgSet(*signaljetpt,*signaljeteta,*njets, *weight,*ph_nhiso_sig_livia,*ph_hoe_sig_livia, *ph_rndiso08_sig_livia);
  RooArgSet* ntplVarsSpectatorsPurity = new RooArgSet( *ph_eta_sig_livia, *ph_phi_sig_livia, *ph_pt_sig_livia, *ph_iso_sig_livia, *ph_rndiso04_sig_livia,*ph_chiso_sig_livia, *ph_sieie_sig_livia, *nphotonsPurity);
  RooArgSet* ntplVarsGen = new RooArgSet(*wzid, *ismatch);
  ntplVars->add(*ntplVarsSpectators);
  ntplVars->add(*ntplVarsSpectatorsPurity);
  ntplVars->add(*ntplVarsGen);
  return ntplVars;
}


void runAllfits(std::string phid){
  TString card_name("./InitWS.rs");

  HLFactory hlfptcut1("HLFactoryptcut1", card_name, false);
  HLFactory hlfptcut2("HLFactoryptcut2", card_name, false);
  HLFactory hlfptcut3("HLFactoryptcut3", card_name, false);
  HLFactory hlfptcut4("HLFactoryptcut4", card_name, false);
  HLFactory hlfptcut5("HLFactoryptcut5", card_name, false);
  RooWorkspace* wptcut1 = hlfptcut1.GetWs(); 
  RooWorkspace* wptcut2 = hlfptcut2.GetWs(); 
  RooWorkspace* wptcut3 = hlfptcut3.GetWs(); 
  RooWorkspace* wptcut4 = hlfptcut4.GetWs(); 
  RooWorkspace* wptcut5 = hlfptcut5.GetWs(); 

  HLFactory hlfEBNEWptcut1("HLFactoryEBNEWptcut1", card_name, false);
  HLFactory hlfEBNEWptcut2("HLFactoryEBNEWptcut2", card_name, false);
  HLFactory hlfEBNEWptcut3("HLFactoryEBNEWptcut3", card_name, false);
  HLFactory hlfEBNEWptcut4("HLFactoryEBNEWptcut4", card_name, false);
  HLFactory hlfEBNEWptcut5("HLFactoryEBNEWptcut5", card_name, false);
  RooWorkspace* wEBNEWptcut1 = hlfEBNEWptcut1.GetWs(); 
  RooWorkspace* wEBNEWptcut2 = hlfEBNEWptcut2.GetWs(); 
  RooWorkspace* wEBNEWptcut3 = hlfEBNEWptcut3.GetWs(); 
  RooWorkspace* wEBNEWptcut4 = hlfEBNEWptcut4.GetWs(); 
  RooWorkspace* wEBNEWptcut5 = hlfEBNEWptcut5.GetWs(); 


   //Add Variables
  RooArgSet* ntplVars = defineVariables();  
  TChain* Tree  = new TChain("tree/tree"); 
  //Tree->Add("../rootfiles/Single_Photon_Run2B_50ns_0825_Loose_LooseRC_0000.root");////
  //Tree->Add("../rootfiles/Single_Photon_Run2D_25ns_2301_0000.root");////
  //  Tree->Add("SinglePhoton-Run2016B-PromptReco-v2_0000.root");////
  //Tree->Add("SinglePhoton-Run2016B-PromptReco-v2_0001.root");////
  //Tree->Add("SinglePhoton-Run2016C-PromptReco-v2_0000.root");////
  // Tree->Add("datatree.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016B_23Sep2016_v2_161109_090401.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016C-23Sep2016-v1_161109_090505.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016D-23Sep2016-v1_161109_090340.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016E-23Sep2016-v1_161109_090258.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016F-23Sep2016-v1_161116_161001.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016G-23Sep2016-v1_161109_090319.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016H-PromptReco-v2_161109_090444.root");////
  Tree->SetTitle("Tree");
  Tree->SetName("Tree");
  wptcut1->import(*ntplVars);
  wptcut2->import(*ntplVars);
  wptcut3->import(*ntplVars);
  wptcut4->import(*ntplVars);
  wptcut5->import(*ntplVars);
  
  wEBNEWptcut1->import(*ntplVars);
  wEBNEWptcut2->import(*ntplVars);
  wEBNEWptcut3->import(*ntplVars);
  wEBNEWptcut4->import(*ntplVars);
  wEBNEWptcut5->import(*ntplVars);

  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "Entering runfits() ..." << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;


  //EB medium
  runfits(wptcut1,"EB",phid,"ptcut1");
  RooRealVar* purityptcut1 = wptcut1->var(("purityEB"+phid+"ptcut1").c_str());
  RooRealVar* puritySystptcut1 = wptcut1->var(("puritySystEB"+phid+"ptcut1").c_str());

  // cout << "First fit done. Now exiting" << endl;
  // exit(0);
    
  runfits(wptcut2,"EB",phid,"ptcut2");
  RooRealVar* purityptcut2 = wptcut2->var(("purityEB"+phid+"ptcut2").c_str());
  RooRealVar* puritySystptcut2 = wptcut2->var(("puritySystEB"+phid+"ptcut2").c_str());
  
  runfits(wptcut3,"EB",phid,"ptcut3");
  RooRealVar* purityptcut3 = wptcut3->var(("purityEB"+phid+"ptcut3").c_str());
  RooRealVar* puritySystptcut3 = wptcut3->var(("puritySystEB"+phid+"ptcut3").c_str());
 
  runfits(wptcut4,"EB",phid,"ptcut4");
  RooRealVar* purityptcut4 = wptcut4->var(("purityEB"+phid+"ptcut4").c_str());
  RooRealVar* puritySystptcut4 = wptcut4->var(("puritySystEB"+phid+"ptcut4").c_str());

  runfits(wptcut5,"EB",phid,"ptcut5");
  RooRealVar* purityptcut5 = wptcut5->var(("purityEB"+phid+"ptcut5").c_str());
  RooRealVar* puritySystptcut5 = wptcut5->var(("puritySystEB"+phid+"ptcut5").c_str());

  //print summary table MEDIUM
  std::cout<<"\\begin{table}[tb]"<<std::endl;
  std::cout<<"\\begin{tabular}{||c|c||}"<<std::endl;
  std::cout<<"$\\p_T$ -- region & Purity  \\\\"<<std::endl;
  std::cout<<"\\hline"<<std::endl;
  // std::cout<<" 175-250 GeV & "<<*purityptcut1->format(2, "EXP") <<"(stat) \\pm "<<"0.0285(syst)    \\\\"<<std::endl;  
  //std::cout<<" 250-300 GeV & "<<*purityptcut2->format(2, "EXP") <<"(stat) \\pm "<<"0.0217(syst)   \\\\"<<std::endl;   
  //std::cout<<" 300-Inf GeV & "<<*purityptcut3->format(2, "EXP") <<"(stat) \\pm "<<"0.0156(syst)   \\\\"<<std::endl;   
  //std::cout<<" 350-Inf GeV & "<<*purityptcut4->format(2, "EXP") <<"(stat) \\pm "<<"0.0156(syst)   \\\\"<<std::endl;   
  std::cout<<" 175-250 GeV & "<<*(purityptcut1->format(2, "EXP")) <<"(stat) \\pm "<<*(puritySystptcut1->format(3, "EXP")) <<" (syst)    \\\\"<<std::endl;  
  std::cout<<" 250-300 GeV & "<<*(purityptcut2->format(2, "EXP")) <<"(stat) \\pm "<<*(puritySystptcut2->format(3, "EXP")) <<" (syst)   \\\\"<<std::endl;   
  std::cout<<" 300-400 GeV & "<<*(purityptcut3->format(2, "EXP")) <<"(stat) \\pm "<<*(puritySystptcut3->format(3, "EXP")) <<" (syst)   \\\\"<<std::endl;   
  std::cout<<" 400-600 GeV & "<<*(purityptcut4->format(2, "EXP")) <<"(stat) \\pm "<<*(puritySystptcut4->format(3, "EXP")) <<" (syst)   \\\\"<<std::endl;   
  std::cout<<" 600-Inf GeV & "<<*(purityptcut5->format(2, "EXP")) <<"(stat) \\pm "<<*(puritySystptcut5->format(3, "EXP")) <<" (syst)   \\\\"<<std::endl;   
 
  std::cout<<"\\hline"<<std::endl;
  std::cout<<"\\end{tabular}"<<std::endl;
  std::cout<<"\\end{table}"<<std::endl;


  TCanvas* c = new TCanvas("c","c",1);


  double x[5] = {225,275, 350, 500, 800};  // center of bin (Inf is 1000)
  double x_err[5] = {25., 25, 50, 100, 200};//10., 10., 10., 20, 30, 335  // value so that the +,- error band covers the whole bin
  double y[5] ={purityptcut1->getVal(),purityptcut2->getVal(),purityptcut3->getVal(), purityptcut4->getVal(),purityptcut5->getVal()};
  double y_err[5] = {purityptcut1->getError(),purityptcut2->getError(),purityptcut3->getError(), purityptcut4->getError(),purityptcut5->getError()};
  //double y_errSyst[5] = {0.0286,0.021,0.0156};
  double y_errSyst[5] = {puritySystptcut1->getVal(),puritySystptcut2->getVal(),puritySystptcut3->getVal(), puritySystptcut4->getVal(),puritySystptcut5->getVal()};
 
 
  //double x[4] = {225,275, 325.,675};
  //double x_err[4] = {25., 25, 25., 350};//10., 10., 10., 20, 30, 335
  //double y[4] ={purityptcut1->getVal(),purityptcut2->getVal(),purityptcut3->getVal(),purityptcut4->getVal()};
  //double y_err[4] = {purityptcut1->getError(),purityptcut2->getError(),purityptcut3->getError(),purityptcut4->getError()};
  //double y_errSyst[4] = {0.0286,0.021,0.0156,0.0156};
 
  TGraphErrors* EB = new TGraphErrors(5, x, y, x_err, y_err);
  TGraphErrors* EBSyst = new TGraphErrors(5, x, y, x_err, y_errSyst);
  EB->SetMarkerColor(kAzure+7);
  EB->SetLineColor(kAzure+7);
  EBSyst->SetMarkerColor(kAzure+7);
  EBSyst->SetLineColor(kAzure+7);

  TH1F* h = new TH1F("h", "", 3, 175, 1000);
  h ->Draw("hist"); 
  h ->GetYaxis()->SetRangeUser(0.3, 1.);
  h ->GetYaxis()->SetTitle("Photon Purity");
  h ->GetXaxis()->SetTitle("p^{#gamma}_{T} [GeV]");
  EB->Draw("PESAME");
  EBSyst->Draw("L[]same");
  TLegend *leg = new TLegend(0.43834677,0.200339,0.7945968,0.3958475, "","brNDC");
  leg->SetTextSize(0.045);
  leg->SetTextFont(42);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->AddEntry(EBSyst,"Barrel = MEDIUM ID ","LPE");
 
  leg->Draw("same");
  // CMS_lumi(c,false,0 );
  c->SaveAs(("purity_vs_Pt_"+phid+".png").c_str());
  c->SaveAs(("purity_vs_Pt_"+phid+".pdf").c_str());

  
 
}





void runfitsPlot(std::string region, std::string phid) {
  TString mainCut;
  TString sieieCut;
  TString sieieCutSideBand;
  TString ptCut;
  Tptcut = "all";
   if(region=="EB"){
    Tregion="EB";
     
    if(phid=="NEW"){
      Twp = "NEW";
      mainCut="(hltphoton165 > 0 || hltphoton175 > 0) && abs(phPurityeta)>0. && abs(phPurityeta)<1.4442 && phPurityElectronVeto && nphotonsPurity>=1 && phPuritypt>175 && phPurityCHiso<5.";//abs(phPurityeta)<1.4442 && 
      mainCutQCD="(hltphoton165 > 0 || hltphoton175 > 0) && abs(phPurityeta)>0. && abs(phPurityeta)<1.4442 && phPurityElectronVeto && nphotonsPurity>=1 && phPuritypt>175 && phPurityCHiso<5.";
      sieieCut=" && phPuritysieie<0.0105";
      sieieCutSideBand=" && phPuritysieie>0.0105 && phPuritysieie<0.014";
      PHIsoAbsCut= 2.75; 
      PHIsoRelCut = "0.0045";  
      RescaledVar =TString::Format("(2.5-rho*phPurityEA+phPurityPHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");  
      RescaledRNDVar =TString::Format("(2.5-rho*phPurityEA+phPurityRND04PHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");  
    }else if(phid=="MEDIUM"){
      Twp = "MEDIUM";
      mainCut="(hltphoton165 > 0 || hltphoton175 > 0) && abs(pPurityheta)<1.4442 && phPurityElectronVeto && nphotonsPurity>=1 && phPuritypt>175 && phPurityCHiso<1.37  &&  phPurityNHiso< (1.06+0.014*phPuritypt + 0.000019*phPuritypt*phPuritypt)";
      mainCutQCD="(hltphoton165 > 0 || hltphoton175 > 0) && abs(pPurityheta)<1.4442 && nphotonsPurity>=1 && phPuritypt>175 && phPurityCHiso<1.37  &&  phPurityNHiso<(1.06+0.014*phPuritypt + 0.000019*phPuritypt*phPuritypt )";
      sieieCut=" && phPuritysieie<0.0102";
      sieieCutSideBand=" && phPuritysieie>0.0102 && phPuritysieie<0.014";
      PHIsoAbsCut= 0.28;   
      PHIsoRelCut = "0.0053"; 
      RescaledVar =TString::Format("(-rho*phPurityEAEGamma+phPurityPHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");  
      RescaledRNDVar =TString::Format("(-rho*phPurityEAEGamma+phPurityPHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");  
      //RescaledRNDVar =TString::Format("(-rho*phPurityEAEGamma+phPurityRND08PHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");  
      //RescaledRNDVar =TString::Format("(-rho*phPurityEAEGamma+phPurityRND04PHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");  
    }else if(phid=="MEDIUM50ns"){
      Twp = "MEDIUM50ns";
      mainCut="(hltphoton165 > 0 || hltphoton175 > 0) && abs(phPurityeta)<1.4442 && phPurityElectronVeto && nphotonsPurity>=1 && phPuritypt>175 && phPurityCHiso<1.31  &&  phPurityNHiso< 0.6+exp(0.0044*phPuritypt+0.5809)";
      mainCutQCD="(hltphoton165 > 0 || hltphoton175 > 0) && abs(phPurityeta)<1.4442 && nphotonsPurity>=1 && phPuritypt>175 && phPurityCHiso<1.31  &&  phPurityNHiso< 0.6+exp(0.0044*phPuritypt+0.5809)";
      sieieCut=" && phPuritysieie<0.010";
      sieieCutSideBand=" && phPuritysieie>0.010 && phPuritysieie<0.014";
      PHIsoAbsCut= 1.33;   
      PHIsoRelCut = "0.0043"; 
      RescaledVar =TString::Format("(-rho*phPurityEAEGamma+phPurityPHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");  
      RescaledRNDVar =TString::Format("(-rho*phPurityEAEGamma+phPurityPHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");  
    }
   
  }
  
 
  mainCut;
  std::cout<<mainCut<<std::endl;
  TString card_name("./InitWS.rs");
  HLFactory hlf("HLFactory", card_name, false);
  RooWorkspace* w = hlf.GetWs();  

  //Add Variables
  RooArgSet* ntplVars = defineVariables();  
  w->import(*ntplVars);
 

  bool makeplots=1;
  if(makeplots){
    // PlotSignalTemplates(w);
    PlotBkgTemplates(w);
    //   PlotSignalTemplateMCs(w);   // to fix
    //PlotBkgTemplateMCs(w);   
  }


}



void runfits(RooWorkspace* w,std::string region, std::string phid, std::string ptcut) {
  TString mainCut;
  TString mainCutQCD;
  TString sieieCut;
  TString sieieCutSideBand;
  TString ptCut;
  Tptcut = ptcut;
  if(ptcut=="ptcut1") ptCut = " && phPuritypt<250";
  if(ptcut=="ptcut2") ptCut = " && phPuritypt>250 && phPuritypt<300";
  if(ptcut=="ptcut3") ptCut = " && phPuritypt>=300 && phPuritypt<400";
  if(ptcut=="ptcut4") ptCut = " && phPuritypt>400 && phPuritypt<600";
  if(ptcut=="ptcut5") ptCut = " && phPuritypt>600 && phPuritypt<=1000";
  if(ptcut=="ptcut6") ptCut = " && phPuritypt>200";


  if(region=="EB"){
    Tregion="EB";
     
    if(phid=="NEW"){
      Twp = "NEW";
      mainCut="(hltphoton165 > 0 || hltphoton175 > 0) && abs(phPurityeta)>0. && abs(phPurityeta)<1.4442 && phPurityElectronVeto && nphotonsPurity>=1 && phPuritypt>175 && phPurityCHiso<5.";//abs(phPurityeta)<1.4442 && 
      mainCutQCD="(hltphoton165 > 0 || hltphoton175 > 0) && abs(phPurityeta)>0. && abs(phPurityeta)<1.4442 && nphotonsPurity>=1 && phPuritypt>175 && phPurityCHiso<5.";
      sieieCut=" && phPuritysieie<0.0105";
      sieieCutSideBand=" && phPuritysieie>0.0105 && phPuritysieie<0.014";
      PHIsoAbsCut= 2.75; 
      PHIsoRelCut = "0.0045";  
      RescaledVar =TString::Format("(2.5-rho*phPurityEA+phPurityPHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");  
      RescaledRNDVar =TString::Format("(2.5-rho*phPurityEA+phPurityRND04PHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");  
    }else if(phid=="MEDIUM"){
      Twp = "MEDIUM";
      mainCut="combinejetpt[0]>100 && abs(combinejeteta[0])<2.5 && (hltphoton165 > 0 || hltphoton175 > 0) && abs(pPurityheta)<1.4442 && phPurityElectronVeto && nphotonsPurity>=1 && phPuritypt>175 && phPurityCHiso<1.37  &&  phPurityNHiso< (1.06+0.014*phPuritypt + 0.000019*phPuritypt*phPuritypt)";
      mainCutQCD="combinejetpt[0]>100 && abs(combinejeteta[0])<2.5 && (hltphoton165 > 0 || hltphoton175 > 0) && abs(pPurityheta)<1.4442 && nphotonsPurity>=1 && phPuritypt>175 && phPurityCHiso<1.37  &&  phPurityNHiso<(1.06+0.014*phPuritypt + 0.000019*phPuritypt*phPuritypt )";
      sieieCut=" && phPuritysieie<0.0102";
      sieieCutSideBand=" && phPuritysieie>0.0102 && phPuritysieie<0.014";
      PHIsoAbsCut= 0.28;   
      PHIsoRelCut = "0.0053"; 
      RescaledVar =TString::Format("(-rho*phPurityEAEGamma+phPurityPHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");  
      RescaledRNDVar =TString::Format("(-rho*phPurityEAEGamma+phPurityPHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");  
    }else if(phid=="MEDIUM50ns"){
      Twp = "MEDIUM50ns";
      mainCut="(hltphoton165 > 0 || hltphoton175 > 0) && abs(phPurityeta)<1.4442 && phPurityElectronVeto && nphotonsPurity==1 && phPuritypt>175 && phPurityCHiso<1.31  &&  phPurityNHiso< (0.6+exp(0.0044*phPuritypt+0.5809))";
      mainCutQCD="(hltphoton165 > 0 || hltphoton175 > 0) && abs(phPurityeta)<1.4442 && nphotonsPurity==1 && phPuritypt>175 && phPurityCHiso<1.31  &&  phPurityNHiso< (0.6+exp(0.0044*phPuritypt+0.5809))";
      sieieCut=" && phPuritysieie<0.010";
      sieieCutSideBand=" && phPuritysieie>0.010 && phPuritysieie<0.014";
      PHIsoAbsCut= 1.33;   
      PHIsoRelCut = "0.0043"; 
      RescaledVar =TString::Format("(-rho*phPurityEAEGamma+phPurityPHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)"); 
      RescaledRNDVar =TString::Format("(-rho*phPurityEAEGamma+phPurityRND04PHiso-")+PHIsoRelCut+TString::Format("*phPuritypt)");   
    }
   
  }
  
  mainCut += ptCut;
  mainCutQCD += ptCut;
  std::cout<<mainCut<<std::endl;
 
  TChain* Tree  = new TChain("tree/tree"); 
  //Tree->Add("../rootfiles/Single_Photon_Run2D_25ns_2301_0000.root");////
  //  Tree->Add("SinglePhoton_data_20062016.root");////
  // Tree->Add("datatree.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016B_23Sep2016_v2_161109_090401.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016C-23Sep2016-v1_161109_090505.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016D-23Sep2016-v1_161109_090340.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016E-23Sep2016-v1_161109_090258.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016F-23Sep2016-v1_161116_161001.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016G-23Sep2016-v1_161109_090319.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016H-PromptReco-v2_161109_090444.root");////


  //  Tree->Add("SinglePhoton-Run2016B-PromptReco-v2_0000.root");////
  //Tree->Add("SinglePhoton-Run2016B-PromptReco-v2_0001.root");////
  //Tree->Add("SinglePhoton-Run2016C-PromptReco-v2_0000.root");////

  Tree->SetTitle("Tree");
  Tree->SetName("Tree");
  w->import(*Tree);
  bool makefit=1;
  AddData(w, Tree);  
  cout << endl; cout << "Now Add Data" << endl;
  

  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "Going to fit ..." << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;


  if(makefit){  
    //AddSignalTemplateRND08(w, Tree);  
    AddSignalTemplateRND04(w, Tree);  
    AddBkgTemplate(w, Tree);    
    ModelFit(w);
    // >>>>>>> For now no systematics
    //ComputeAllSystUnc(w, "EB", "MEDIUM");
    //PlotAllSystUnc(w, "EB", "MEDIUM");
 
  }
 
  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "Fit done ..." << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;

  std::cout<<w->var(TString::Format("purity")+Tregion+Twp+Tptcut)<<std::endl;

}



 

// Add Data Set
void AddData(RooWorkspace* w, TTree* tree) {  

  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "Entered AddData() ..." << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;

  TH1F* h_data = new TH1F("h_data", "h_data", nbins, minFit, maxFit );
  h_data->Sumw2();
  tree->Draw("phPurityPHiso>>h_data",mainCut+sieieCut+TString::Format(" && njets<=3"),"goff");//nJets<=3
  RooDataHist* ReducedDataHist = new RooDataHist("DataHist","DataHist",*w->var("phPurityPHiso"), Import(*h_data) );
  w->import(*ReducedDataHist);

  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "Exiting AddData() ..." << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;
  cout << "============================" << endl;

  
  }



// Add Signal Templates
void AddSignalTemplate(RooWorkspace* w, TTree* tree) {  
  TH1F* h_sig = new TH1F("h_sig", "h_sig", nbins, minFit, maxFit );
  h_sig->Sumw2();
  tree->Draw("phPurityPHiso>>h_sig",mainCut+sieieCut+TString::Format(" && njets<=3"),"goff");
  RooDataHist* ReducedDataHist = new RooDataHist("SignalDataHist","SignalDataHist",*w->var("phPurityPHiso"), Import(*h_sig) );
  RooHistPdf* ReducedHistPdf = new RooHistPdf("SignalDataHistPdf","SignalDataHistPdf",*w->var("phPurityPHiso"),*ReducedDataHist);
  w->import(*ReducedDataHist);  
  w->import(*ReducedHistPdf); 
 

  }

void AddSignalTemplateRND04(RooWorkspace* w, TTree* tree) {  
  TH1F* h_sigRND = new TH1F("h_sigRND", "h_sigRND", nbins, minFit, maxFit );
  h_sigRND->Sumw2();  
  tree->Draw("phPurityRND04PHiso>>h_sigRND",mainCut+sieieCut+TString::Format(""),"goff");
  RooDataHist* ReducedDataHist = new RooDataHist("SignalDataHistRND","SignalDataHistRND",*w->var("phPurityRND04PHiso"), Import(*h_sigRND) );
  RooHistPdf* ReducedHistPdf = new RooHistPdf("SignalDataHistPdfRND","SignalDataHistPdfRND",*w->var("phPurityRND04PHiso"),*ReducedDataHist);
  w->import(*ReducedDataHist);  
  w->import(*ReducedHistPdf);
   

  TH1F* h_sigRNDShape = new TH1F("h_sigRNDShape", "h_sigRNDShape", nbins, minFit, maxFit );
  h_sigRNDShape->Sumw2();
  tree->Draw("phPurityRND04PHiso>>h_sigRNDShape",mainCut+sieieCut+TString::Format(""),"goff");
  RooDataHist* ReducedDataHistRNDShape = new RooDataHist("SignalDataHistRNDShape","SignalDataHistRNDShape",*w->var("phPurityPHiso"), Import(*h_sigRNDShape) ); 
  RooHistPdf* ReducedHistPdfRNDShape = new RooHistPdf("SignalDataHistPdfRNDShape","SignalDataHistPdfRNDShape",*w->var("phPurityPHiso"),*ReducedDataHistRNDShape);
  w->import(*ReducedDataHistRNDShape);  
  w->import(*ReducedHistPdfRNDShape);  
  std::cout<<"========================> "<<h_sigRNDShape->Integral()<<"  "<<ReducedDataHistRNDShape->sumEntries()<<"   "<<ReducedDataHistRNDShape->numEntries()<<std::endl;

 
  }
void AddSignalTemplateRND08(RooWorkspace* w, TTree* tree) {  
  TH1F* h_sigRND = new TH1F("h_sigRND", "h_sigRND", nbins, minFit, maxFit );
  h_sigRND->Sumw2();  
  tree->Draw("phPurityRND08PHiso>>h_sigRND",mainCut+sieieCut+TString::Format(""),"goff");
  RooDataHist* ReducedDataHist = new RooDataHist("SignalDataHistRND","SignalDataHistRND",*w->var("phPurityRND08PHiso"), Import(*h_sigRND) );
  RooHistPdf* ReducedHistPdf = new RooHistPdf("SignalDataHistPdfRND","SignalDataHistPdfRND",*w->var("phPurityRND08PHiso"),*ReducedDataHist);
  w->import(*ReducedDataHist);  
  w->import(*ReducedHistPdf);
   

  TH1F* h_sigRNDShape = new TH1F("h_sigRNDShape", "h_sigRNDShape", nbins, minFit, maxFit );
  h_sigRNDShape->Sumw2();
  tree->Draw("phPurityRND08PHiso>>h_sigRNDShape",mainCut+sieieCut+TString::Format(""),"goff");
  RooDataHist* ReducedDataHistRNDShape = new RooDataHist("SignalDataHistRNDShape","SignalDataHistRNDShape",*w->var("phPurityPHiso"), Import(*h_sigRNDShape) ); 
  RooHistPdf* ReducedHistPdfRNDShape = new RooHistPdf("SignalDataHistPdfRNDShape","SignalDataHistPdfRNDShape",*w->var("phPurityPHiso"),*ReducedDataHistRNDShape);
  w->import(*ReducedDataHistRNDShape);  
  w->import(*ReducedHistPdfRNDShape);  
  std::cout<<"========================> "<<h_sigRNDShape->Integral()<<"  "<<ReducedDataHistRNDShape->sumEntries()<<"   "<<ReducedDataHistRNDShape->numEntries()<<std::endl;

  
}


// Add Signal Template
void AddBkgTemplate(RooWorkspace* w, TTree* tree) {   
  TH1F* h_bkg = new TH1F("h_bkg", "h_bkg", nbins, minFit, maxFit );
  h_bkg->Sumw2();
  tree->Draw("phPurityPHiso>>h_bkg",mainCut+sieieCutSideBand+TString::Format(""),"goff");
  RooDataHist* ReducedDataHist = new RooDataHist("BkgDataHist","BkgDataHist",*w->var("phPurityPHiso"), Import(*h_bkg) );
  RooHistPdf* ReducedHistPdf = new RooHistPdf("BkgDataHistPdf","BkgDataHistPdf",*w->var("phPurityPHiso"),*ReducedDataHist);
  w->import(*ReducedDataHist);  
  w->import(*ReducedHistPdf);  


  }

void AddSignalMCRNDRC08(RooWorkspace* w){

  RooArgSet* ntplVars_match = defineVariables_match();  
  TChain* Tree  = new TChain("tree/tree"); 
  //  Tree->Add("../rootfiles/GJets_HTs_RC08.root");////
  //Tree->Add("../rootfiles/GJets_0611_05Oct2015_Loose_0000_MEDIUMIDok.root");////
  //Tree->Add("GJet_MC_20062016.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-40To100.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-100To200.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-200To400.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-400To600.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-600ToInf.root");////
  Tree->SetTitle("TreeMCRC08");
  Tree->SetName("TreeMCRC08");
  // RooDataSet SignalMC("MCRC08","SignalMCRC08",Tree,*ntplVars_match,mainCutQCD+sieieCut+TString::Format(""),"1");
 

  TH1F* h_sigMCRNDRC08 = new TH1F("h_sigMCRNDRC08", "h_sigMCRNDRC08", nbins, minFit, maxFit );
  h_sigMCRNDRC08->Sumw2();  
  Tree->Draw("phPurityRND08PHiso>>h_sigMCRNDRC08",mainCutQCD+TString::Format("*xsec*wgt/wgtsum"),"goff");
  //Tree->Draw("phPurityRND08PHiso>>h_sigMCRNDRC08",mainCutQCD+TString::Format("*1"));
 
  TH1F* h_sigMCRNDShapeRC08 = new TH1F("h_sigMCRNDShapeRC08", "h_sigMCRNDShapeRC08", nbins, minFit, maxFit );
  h_sigMCRNDShapeRC08->Sumw2();  
  Tree->Draw("phPurityRND08PHiso>>h_sigMCRNDShapeRC08",TString::Format("(")+mainCutQCD+TString::Format("*xsec*wgt/wgtsum)"),"goff");
  //Tree->Draw("phPurityRND08PHiso>>h_sigMCRNDShapeRC08",TString::Format("(")+mainCutQCD+TString::Format("*1)"));
  RooDataHist* ReducedDataHistShapeRC08 = new RooDataHist("SignalMCHistRNDShapeRC08","SignalMCHistRNDShapeRC08",*w->var("phPurityPHiso"), Import(*h_sigMCRNDShapeRC08) );
  RooHistPdf* ReducedHistPdfShapeRC08 = new RooHistPdf("SignalMCHistPdfRNDShapeRC08","SignalMCHistPdfRNDShapeRC08",*w->var("phPurityPHiso"),*ReducedDataHistShapeRC08);
  w->import(*ReducedDataHistShapeRC08);  
  w->import(*ReducedHistPdfShapeRC08);

  std::cout<<"flag"<<std::endl;
  ///add MC matching
  TH1F* h_sigMCRNDmatch = new TH1F("h_sigMCRNDmatch", "h_sigMCRNDmatch", nbins, minFit, maxFit );
  h_sigMCRNDmatch->Sumw2();  
  Tree->Draw("phPurityRND08PHiso>>h_sigMCRNDmatch",TString::Format("(ismatch && ")+mainCutQCD+TString::Format(")*xsec*wgt/wgtsum"),"goff");
  //Tree->Draw("phPurityRND08PHiso>>h_sigMCRNDmatch",TString::Format("(ismatch && ")+mainCutQCD+TString::Format(")*1"));
  TH1F* h_sigMCRNDShapeMatch = new TH1F("h_sigMCRNDShapeMatch", "h_sigMCRNDShapeMatch", nbins, minFit, maxFit );
  h_sigMCRNDShapeMatch->Sumw2();  
  //  Tree->Draw("phPurityRND08PHiso>>h_sigMCRNDShapeMatch",TString::Format("(ismatch && ")+mainCutQCD+TString::Format(")*xsec*wgt/wgtsum"));
  Tree->Draw("phPurityRND08PHiso>>h_sigMCRNDShapeMatch",TString::Format("(sqrt(pow((wzeta-pPurityheta),2)+pow(deltaPhi(wzphi,phPurityphi),2))<0.3 && ")+mainCutQCD+TString::Format(")*xsec*wgt/wgtsum"),"goff");
  //Tree->Draw("phPurityRND08PHiso>>h_sigMCRNDShapeMatch",TString::Format("(ismatch && ")+mainCutQCD+TString::Format(")*1"));
  RooDataHist* ReducedDataHistShapeMatch = new RooDataHist("SignalMCHistRNDShapeMatch","SignalMCHistRNDShapeMatch",*w->var("phPurityPHiso"), Import(*h_sigMCRNDShapeMatch) );
  RooHistPdf* ReducedHistPdfShapeMatch = new RooHistPdf("SignalMCHistPdfRNDShapeMatch","SignalMCHistPdfRNDShapeMatch",*w->var("phPurityPHiso"),*ReducedDataHistShapeMatch);
  w->import(*ReducedDataHistShapeMatch);  
  w->import(*ReducedHistPdfShapeMatch);


    std::cout<<"flag"<<std::endl;

}



void AddSignalMCRND(RooWorkspace* w){
 std::cout<<"flag"<<std::endl;
  RooArgSet* ntplVars_match = defineVariables_match();  
  TChain* Tree  = new TChain("tree/tree"); 
  //Tree->Add("../rootfiles/GJets_0611_05Oct2015_Loose_0000_MEDIUMIDok.root");////
  //Tree->Add("GJet_MC_20062016.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-40To100.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-100To200.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-200To400.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-400To600.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-600ToInf.root");////


  Tree->SetTitle("TreeMC");
  Tree->SetName("TreeMC");
   std::cout<<"flag"<<std::endl;
  TH1F* h_sigMCRND = new TH1F("h_sigMCRND", "h_sigMCRND", nbins, minFit, maxFit );
  h_sigMCRND->Sumw2();  
  Tree->Draw("phPurityRND04PHiso>>h_sigMCRND",mainCutQCD+TString::Format("*xsec*wgt/wgtsum"),"goff");
  //Tree->Draw("phPurityRND04PHiso>>h_sigMCRND",mainCutQCD+TString::Format("*1"));
  std::cout<<"flag22"<<std::endl;

  TH1F* h_sigMCRNDShape = new TH1F("h_sigMCRNDShape", "h_sigMCRNDShape", nbins, minFit, maxFit );
  h_sigMCRNDShape->Sumw2();  
  Tree->Draw("phPurityRND04PHiso>>h_sigMCRNDShape",mainCutQCD+TString::Format("*xsec*wgt/wgtsum"),"goff");
  //Tree->Draw("phPurityRND04PHiso>>h_sigMCRNDShape",mainCutQCD+TString::Format("*1"));
  RooDataHist* ReducedDataHistShapeForPlot = new RooDataHist("SignalMCHistRNDShapeForPlot","SignalMCHistRNDShapeForPlot",RooArgList(*w->var("phPurityPHiso"),*w->var("phPuritypt")), Import(*h_sigMCRNDShape) );
  RooDataHist* ReducedDataHistShape = new RooDataHist("SignalMCHistRNDShape","SignalMCHistRNDShape",*w->var("phPurityPHiso"), Import(*h_sigMCRNDShape) );
  RooHistPdf* ReducedHistPdfShape = new RooHistPdf("SignalMCHistPdfRNDShape","SignalMCHistPdfRNDShape",*w->var("phPurityPHiso"),*ReducedDataHistShape);
  w->import(*ReducedDataHistShapeForPlot);  
  w->import(*ReducedDataHistShape);  
  w->import(*ReducedHistPdfShape);
 std::cout<<"flag"<<std::endl;

  ///add MC matching
  TH1F* h_sigMCRNDmatch = new TH1F("h_sigMCRNDmatch", "h_sigMCRNDmatch", nbins, minFit, maxFit );
  h_sigMCRNDmatch->Sumw2();  
  Tree->Draw("phPurityRND04PHiso>>h_sigMCRNDmatch",TString::Format("(ismatch && ")+mainCutQCD+TString::Format(")*xsec*wgt/wgtsum"),"goff");
  //Tree->Draw("phPurityRND04PHiso>>h_sigMCRNDmatch",TString::Format("(ismatch && ")+mainCutQCD+TString::Format(")*1"));

  TH1F* h_sigMCRNDShapeMatch = new TH1F("h_sigMCRNDShapeMatch", "h_sigMCRNDShapeMatch", nbins, minFit, maxFit );
  h_sigMCRNDShapeMatch->Sumw2();  
  Tree->Draw("phPurityRND04PHiso>>h_sigMCRNDShapeMatch",TString::Format("(sqrt(pow((wzeta-pPurityheta),2)+pow(deltaPhi(wzphi,phPurityphi),2))<0.3 && ")+mainCutQCD+TString::Format(")*xsec*wgt/wgtsum"),"goff");
  //Tree->Draw("phPurityRND04PHiso>>h_sigMCRNDShapeMatch",TString::Format("(ismatch && ")+mainCutQCD+TString::Format(")*1"));
  RooDataHist* ReducedDataHistShapeMatch = new RooDataHist("SignalMCHistRNDShapeMatch","SignalMCHistRNDShapeMatch",*w->var("phPurityPHiso"), Import(*h_sigMCRNDShapeMatch) );
  RooHistPdf* ReducedHistPdfShapeMatch = new RooHistPdf("SignalMCHistPdfRNDShapeMatch","SignalMCHistPdfRNDShapeMatch",*w->var("phPurityPHiso"),*ReducedDataHistShapeMatch);
  w->import(*ReducedDataHistShapeMatch);  
  w->import(*ReducedHistPdfShapeMatch);

 
  

}

void AddSignalMC(RooWorkspace* w){
  std::cout<<"flag"<<std::endl;
  RooArgSet* ntplVars_match = defineVariables_match();  
  TChain* Tree  = new TChain("tree/tree"); 
  //Tree->Add("../rootfiles/GJets_0611_05Oct2015_Loose_0000_MEDIUMIDok.root");////
  //  Tree->Add("MCmatching/GJet_Pt-15To6000_0000_wMatching_Loose.root/tree");
  //  Tree->Add("GJet_MC_20062016.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-40To100.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-100To200.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-200To400.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-400To600.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_GJets_HT-600ToInf.root");////
  Tree->SetTitle("TreeMC");
  Tree->SetName("TreeMC");
  std::cout<<"flag"<<std::endl;
  TH1F* h_sigMC = new TH1F("h_sigMC", "h_sigMC", nbins, minFit, maxFit );
  h_sigMC->Sumw2();  
  Tree->Draw("phPurityPHiso>>h_sigMC",mainCutQCD+TString::Format("*xsec*wgt/wgtsum"),"goff");// && njets==1
  //Tree->Draw("phPurityPHiso>>h_sigMC",mainCutQCD+TString::Format("*1"));// && njets==1
  RooDataHist* ReducedDataHist = new RooDataHist("SignalMCHist","SignalMCHist",*w->var("phPurityPHiso"), Import(*h_sigMC) );
  RooHistPdf* ReducedHistPdf = new RooHistPdf("SignalMCHistPdf","SignalMCHistPdf",*w->var("phPurityPHiso"),*ReducedDataHist);
  w->import(*ReducedDataHist);  
  w->import(*ReducedHistPdf);

  std::cout<<"flag"<<std::endl;
  ///add MC matching
  TH1F* h_sigMCmatch = new TH1F("h_sigMCmatch", "h_sigMCmatch", nbins, minFit, maxFit );
  h_sigMCmatch->Sumw2();  
  Tree->Draw("phPurityPHiso>>h_sigMCmatch",TString::Format("(sqrt(pow((wzeta-pPurityheta),2)+pow(deltaPhi(wzphi,phPurityphi),2))<0.3 && ")+mainCutQCD+TString::Format(")*xsec*wgt/wgtsum"),"goff");// && njets==1
  //Tree->Draw("phPurityPHiso>>h_sigMCmatch",TString::Format("(ismatch && ")+mainCutQCD+TString::Format(")*1"));// && njets==1
  RooDataHist* ReducedDataHistMatch = new RooDataHist("SignalMCHistmatch","SignalMCHistmatch",*w->var("phPurityPHiso"), Import(*h_sigMCmatch) );
  RooHistPdf* ReducedHistPdfMatch = new RooHistPdf("SignalMCHistPdfmatch","SignalMCHistPdfmatch",*w->var("phPurityPHiso"),*ReducedDataHistMatch);
  w->import(*ReducedDataHistMatch);  
  w->import(*ReducedHistPdfMatch);


    std::cout<<"flag"<<std::endl;

}

void AddBkgMC(RooWorkspace* w){

  RooArgSet* ntplVars_match = defineVariables_match();  
  TChain* Tree  = new TChain("tree/tree"); 
  //Tree->Add("../rootfiles/QCD_EMEnriched_2810_0000.root");////
  //  Tree->Add("MCmatching/GJet_Pt-15To6000_0000_wMatching_Loose.root/tree");
  //Tree->Add("QCD_MC_20062016.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_QCD_Pt-120to170_EMEnriched.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_QCD_Pt-170to300_EMEnriched.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/tree_crab_QCD_Pt-300toInf_EMEnriched.root");////
  Tree->SetTitle("TreeMC");
  Tree->SetName("TreeMC");
  //  RooDataSet BkgMC("MCbkg","BkgMC",Tree,*ntplVars_match,mainCutQCD+sieieCutSideBand+TString::Format(""),"1");
  // w->import(BkgMC);
 
  TH1F* h_bkgMCNOsieie = new TH1F("h_bkgMCNOsieie", "h_bkgMCNOsieie", nbins, minFit, maxFit );
  h_bkgMCNOsieie->Sumw2();  
  Tree->Draw("phPurityPHiso>>h_bkgMCNOsieie",TString::Format("(")+mainCutQCD+sieieCutSideBand+TString::Format(")*xsec*wgt/wgtsum"),"goff");
  //Tree->Draw("phPurityPHiso>>h_bkgMCNOsieie",TString::Format("(")+mainCutQCD+sieieCutSideBand+TString::Format(")*1"));
  RooDataHist* ReducedDataHist = new RooDataHist("BkgMCHistNOsieie","BkgMCHistNOsieie",*w->var("phPurityPHiso"), Import(*h_bkgMCNOsieie) );
  RooHistPdf* ReducedHistPdf = new RooHistPdf("BkgMCHistPdfNOsieie","BkgMCHistPdfNOsieie",*w->var("phPurityPHiso"),*ReducedDataHist);
  w->import(*ReducedDataHist);  
  w->import(*ReducedHistPdf);
 
  ///inverting MC matching
  TH1F* h_bkgMCNOmatch = new TH1F("h_bkgMCNOmatch", "h_bkgMCNOmatch", nbins, minFit, maxFit );
  h_bkgMCNOmatch->Sumw2();  
  Tree->Draw("phPurityPHiso>>h_bkgMCNOmatch",TString::Format("(sqrt(pow((wzeta-pPurityheta),2)+pow(deltaPhi(wzphi,phPurityphi),2))>2 && ")+mainCutQCD+TString::Format(")*xsec*wgt/wgtsum"),"goff");
  //Tree->Draw("phPurityPHiso>>h_bkgMCNOmatch",TString::Format("(ismatch==0 && ")+mainCutQCD+TString::Format(")*1"));
  RooDataHist* ReducedDataHistMatch = new RooDataHist("BkgMCHistNOmatch","BkgMCHistNOmatch",*w->var("phPurityPHiso"), Import(*h_bkgMCNOmatch) );
  RooHistPdf* ReducedHistPdfMatch = new RooHistPdf("BkgMCHistPdfNOmatch","BkgMCHistPdfNOmatch",*w->var("phPurityPHiso"),*ReducedDataHistMatch);
  w->import(*ReducedDataHistMatch);  
  w->import(*ReducedHistPdfMatch);


  

}

void ReduceSignalMCRND(RooWorkspace* w, TString ptrange, TString suffix){

  RooDataSet* ReducedMC;  
  ReducedMC = (RooDataSet*) w->data("SignalMCHistRNDShapeForPlot")->reduce(*w->var("phPurityRNDPHiso"),ptrange);//// 
  RooDataHist* ReducedMCHist = (RooDataHist*)ReducedMC->binnedClone(TString::Format("SignalMCHist")+suffix,TString::Format("SignalMCHist")+suffix);
  RooHistPdf* ReducedHistPdf = new RooHistPdf(TString::Format("SignalMCHistPdf")+suffix,TString::Format("SignalMCHistPdf")+suffix,RooArgList(*w->var("phPurityPHiso"),*w->var("phPuritypt")),*ReducedMCHist);
  w->import(*ReducedMCHist);  
  w->import(*ReducedHistPdf);  
  cout << endl;
  cout << TString::Format("SignalMC")+suffix << endl;
  ReducedMCHist->Print("v");
  cout << "---- nX:  " <<ReducedMCHist->sumEntries() << endl;   
  cout << endl; 
}

void ReduceBkgMC(RooWorkspace* w, TString ptrange, TString suffix){

  RooDataSet* ReducedMC;  
  ReducedMC = (RooDataSet*) w->data("SignalMCHistRNDShape")->reduce(*w->var("phPurityRNDPHiso"),ptrange);//// 
  RooDataHist* ReducedMCHist = (RooDataHist*)ReducedMC->binnedClone(TString::Format("BkgMCHist")+suffix,TString::Format("BkgMCHist")+suffix);
  RooHistPdf* ReducedHistPdf = new RooHistPdf(TString::Format("BkgMCHistPdf")+suffix,TString::Format("BkgMCHistPdf")+suffix,*w->var("phPurityRND04PHiso"),*ReducedMCHist);
  w->import(*ReducedMCHist);  
  w->import(*ReducedHistPdf);  
  cout << endl;
  cout << TString::Format("BkgMC")+suffix << endl;
  ReducedMCHist->Print("v");
  cout << "---- nX:  " <<ReducedMCHist->sumEntries() << endl;   
  cout << endl; 
}



void PlotSignalTemplateMCs(RooWorkspace* w){
  AddSignalMCRND(w);
  ReduceSignalMCRND(w, "phPuritypt>100 && phPuritypt<=200", "_100to200");  
  ReduceSignalMCRND(w, "phPuritypt>200 && phPuritypt<=250", "_200to250");  
  ReduceSignalMCRND(w, "phPuritypt>250", "_250toInf");  
 

  RooHistPdf* MC100to300 =(RooHistPdf*) w->pdf("SignalMCHistPdf_100to200");
  RooHistPdf* MC300to600 =(RooHistPdf*) w->pdf("SignalMCHistPdf_200to250");
  RooHistPdf* MC600toInf =(RooHistPdf*) w->pdf("SignalMCHistPdf_250toInf");
 

  TCanvas* ctmp = new TCanvas("ctmp","PH Iso Signal Shape",1);
  RooPlot* plot;
  RooPlot* plotrnd;
  RooRealVar* ph_iso = w->var("phPurityRND04PHiso");  
  ph_iso->setUnit("GeV");
  ph_iso->setRange("fit range", minFit, maxFit); 
  plot = ph_iso->frame(minFit, maxFit,40);
  plot->GetYaxis()->SetTitle("a.u.");
 

  plot->GetXaxis()->SetTitle("PF Photon Iso [GeV]");
  MC100to300->plotOn(plot,LineColor(kAzure));
  MC300to600->plotOn(plot,LineColor(kAzure+4));
  MC600toInf->plotOn(plot,LineColor(kAzure+9));
 

  plot->Draw();
  plot->GetYaxis()->SetRangeUser(0.0001, 0.5);

  TLegend *legdata = new TLegend(0.4834677,0.680339,0.7945968,0.8958475, "","brNDC");
  legdata->SetHeader("MC Signal - RC");

  legdata->AddEntry(plot->getObject(0),"p_{T} in [100-200] GeV","L");
  legdata->AddEntry(plot->getObject(1),"p_{T} in ]200-250] GeV","L");
  legdata->AddEntry(plot->getObject(2),"p_{T} in ]250-Inf] GeV","L");
  
  legdata->SetTextSize(0.035);
  legdata->SetTextFont(42);
  legdata->SetBorderSize(0);
  legdata->SetFillStyle(0);
  legdata->Draw("same");
  //CMS_lumi( ctmp,true,0 );
  ctmp->cd();
  ctmp->SetLogy(0);
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_SignalShape_MC_ptBins.pdf"));
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_SignalShape_MC_ptBins.png"));
  ctmp->SetLogy(1);
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_SignalShape_LOG_MC_ptBins.pdf"));
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_SignalShape_LOG_MC_ptBins.png"));


}



void PlotBkgTemplateMCs(RooWorkspace* w){
  AddBkgMC(w);
  ReduceBkgMC(w, "phPuritypt>100 && phPuritypt<=200", "_100to100");  
  ReduceBkgMC(w, "phPuritypt>200 && phPuritypt<=250", "_100to250");  
  ReduceBkgMC(w, "phPuritypt>250", "_250toInf");  
 

  RooHistPdf* MC100to300 =(RooHistPdf*) w->pdf("BkgMCHistPdf_100to200");
  RooHistPdf* MC300to600 =(RooHistPdf*) w->pdf("BkgMCHistPdf_200to250");
  RooHistPdf* MC600toInf =(RooHistPdf*) w->pdf("BkgMCHistPdf_250toInf");
 

  TCanvas* ctmp = new TCanvas("ctmp","PH Iso Signal Shape",1);
  RooPlot* plot;
  RooPlot* plotrnd;
  RooRealVar* ph_iso = w->var("phPurityRND04PHiso");  
  ph_iso->setUnit("GeV");
  ph_iso->setRange("fit range", minFit, maxFit); 
  plot = ph_iso->frame(minFit, maxFit,40);
  plot->GetYaxis()->SetTitle("a.u.");
 

  plot->GetXaxis()->SetTitle("PF Photon Iso [GeV]");
  MC100to300->plotOn(plot,LineColor(kAzure));
  MC300to600->plotOn(plot,LineColor(kAzure+4));
  MC600toInf->plotOn(plot,LineColor(kAzure+9));
 

  plot->Draw();
  plot->GetYaxis()->SetRangeUser(0.0001, 0.5);

  TLegend *legdata = new TLegend(0.4834677,0.680339,0.7945968,0.8958475, "","brNDC");
  legdata->SetHeader("MC Signal - RC");

  legdata->AddEntry(plot->getObject(0),"p_{T} in [100-200] GeV","L");
  legdata->AddEntry(plot->getObject(1),"p_{T} in ]200-250] GeV","L");
  legdata->AddEntry(plot->getObject(2),"p_{T} in ]250-Inf] GeV","L");
  
  legdata->SetTextSize(0.035);
  legdata->SetTextFont(42);
  legdata->SetBorderSize(0);
  legdata->SetFillStyle(0);
  legdata->Draw("same");
  //CMS_lumi( ctmp,true,0 );
  ctmp->cd();
  ctmp->SetLogy(0);
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_BkgShape_MC_ptBins.pdf"));
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_BkgShape_MC_ptBins.png"));
  ctmp->SetLogy(1);
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_BkgShape_LOG_MC_ptBins.pdf"));
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_BkgShape_LOG_MC_ptBins.png"));


}



void PlotSignalTemplates(RooWorkspace* w){
  TChain* Tree  = new TChain("tree/tree"); 
  //  Tree->Add("SinglePhoton_data_20062016.root");////
  // Tree->Add("datatree.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016B_23Sep2016_v2_161109_090401.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016C-23Sep2016-v1_161109_090505.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016D-23Sep2016-v1_161109_090340.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016E-23Sep2016-v1_161109_090258.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016F-23Sep2016-v1_161116_161001.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016G-23Sep2016-v1_161109_090319.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016H-PromptReco-v2_161109_090444.root");////


  //  Tree->Add("SinglePhoton-Run2016B-PromptReco-v2_0000.root");////
  //  Tree->Add("SinglePhoton-Run2016B-PromptReco-v2_0001.root");////
  //Tree->Add("SinglePhoton-Run2016C-PromptReco-v2_0000.root");////

  Tree->SetTitle("Tree");
  Tree->SetName("Tree");
  
  //AddSignalTemplate( w, Tree);
  AddSignalTemplateRND04(w, Tree);

  std::cout<<"flag"<<std::endl;
  AddSignalMCRND(w);
  AddSignalMCRNDRC08(w);
  // AddSignalMCRNDRC04(w);
  AddSignalMC(w);

  //RooHistPdf* DataPdf =(RooHistPdf*) w->pdf("SignalDataHistPdf");
  RooHistPdf* DataPdfRND =(RooHistPdf*) w->pdf("SignalDataHistPdfRNDShape");
  RooHistPdf* MCPdfRND =(RooHistPdf*) w->pdf("SignalMCHistPdfRNDShape");
  RooHistPdf* MCPdfRNDRC08 =(RooHistPdf*) w->pdf("SignalMCHistPdfRNDShapeRC08");
  // RooHistPdf* MCPdfRNDRC04 =(RooHistPdf*) w->pdf("SignalMCHistPdfRNDShapeRC04");
  RooHistPdf* MCPdfmatch =(RooHistPdf*) w->pdf("SignalMCHistPdfmatch");

  TCanvas* ctmp = new TCanvas("ctmp","Signal Shape",1);
  RooPlot* plot;
  RooPlot* plotrnd;
  RooRealVar* ph_iso = w->var("phPurityPHiso");  
  ph_iso->setUnit("GeV");
  ph_iso->setRange("fit range", minFit, maxFit); 
  plot = ph_iso->frame(minFit, maxFit,40);
 
  //DataPdf->plotOn(plot,LineColor(kOrange+8));
  DataPdfRND->plotOn(plot,LineColor(kSpring+9));
  MCPdfRND->plotOn(plot,LineColor(kAzure+9));
  MCPdfRNDRC08->plotOn(plot,LineColor(kBlue-6));
  // MCPdfRNDRC04->plotOn(plot,LineColor(kBlue-6));
  MCPdfmatch->plotOn(plot,LineColor(kMagenta));
 
  plot->GetXaxis()->SetTitle("PF Photon Iso [GeV]");
  plot->GetYaxis()->SetTitle("a. u.");
  plot->Draw();
  
  TLegend *legdata = new TLegend(0.4734677,0.640339,0.7945968,0.8958475, "","brNDC");
  // legdata->AddEntry(plot->getObject(0),"Data #sigma_{#eta #eta}< 0.0103","L");
  legdata->AddEntry(plot->getObject(0),"Data RC #Delta R > 0.4 ","L");
  legdata->AddEntry(plot->getObject(1),"MC RC #Delta R > 0.4","L");
  legdata->AddEntry(plot->getObject(2),"MC RC #Delta R > 0.8","L");
  legdata->AddEntry(plot->getObject(3),"MC Prompt","L");
  legdata->SetTextSize(0.035);
  legdata->SetTextFont(42);
  legdata->SetBorderSize(0);
  legdata->SetFillStyle(0);
  legdata->Draw("same");
  //  CMS_lumi( ctmp,false,0 );
  ctmp->cd();
  ctmp->SetLogy(0);
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_SignalShape_MC_RC0804.pdf"));
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_SignalShape_MC_RC0804.png"));
  ctmp->SetLogy(1);
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_SignalShape_LOG_MC_RC0804.pdf"));
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_SignalShape_LOG_MC_RC0804.png"));

}


void PlotBkgTemplates(RooWorkspace* w){
  TChain* Tree  = new TChain("tree/tree"); 
  //  Tree->Add("SinglePhoton_data_20062016.root");////
  //  Tree->Add("SinglePhoton-Run2016B-PromptReco-v2_0000.root");////
  //Tree->Add("SinglePhoton-Run2016B-PromptReco-v2_0001.root");////
  //Tree->Add("SinglePhoton-Run2016C-PromptReco-v2_0000.root");////

  //Tree->Add("datatree.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016B_23Sep2016_v2_161109_090401.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016C-23Sep2016-v1_161109_090505.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016D-23Sep2016-v1_161109_090340.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016E-23Sep2016-v1_161109_090258.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016F-23Sep2016-v1_161116_161001.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016G-23Sep2016-v1_161109_090319.root");////
  Tree->Add("/u2/mciprian/TREES_PHOTON_PURITY/merged/SinglePhoton_Run2016H-PromptReco-v2_161109_090444.root");////

  Tree->SetTitle("Tree");
  Tree->SetName("Tree");
  AddBkgTemplate(w, Tree);
  AddBkgMC(w);
  RooHistPdf* DataPdf =(RooHistPdf*) w->pdf("BkgDataHistPdf");
  RooHistPdf* MCPdfNOmatch =(RooHistPdf*) w->pdf("BkgMCHistPdfNOmatch");
  RooHistPdf* MCPdfNOsieie =(RooHistPdf*) w->pdf("BkgMCHistPdfNOsieie");

  TCanvas* ctmp = new TCanvas("ctmp","Signal Shape",1);
  RooPlot* plot;
  RooPlot* plotrnd;
  RooRealVar* ph_iso = w->var("phPurityPHiso");  
  ph_iso->setUnit("GeV");
  ph_iso->setRange("fit range", minFit, maxFit); 
  plot = ph_iso->frame(minFit, maxFit,10);
 
  // DataPdf->plotOn(plot,LineColor(kOrange+8));
  // DataPdf->plotOn(plot,LineColor(kSpring+9));
  
  MCPdfNOmatch->plotOn(plot,LineColor(kMagenta));
  MCPdfNOsieie->plotOn(plot,LineColor(kBlue-6));
 
  plot->GetXaxis()->SetTitle("PF Photon Iso [GeV]");
  plot->GetYaxis()->SetTitle("a. u.");
  plot->Draw();
  
  TLegend *legdata = new TLegend(0.3734677,0.340339,0.6945968,0.6958475, "","brNDC");
  // legdata->AddEntry(plot->getObject(0),"Data #sigma_{#eta #eta} sidebands","LP");
  // legdata->AddEntry(plot->getObject(1),"Data Signal RC","L");
  legdata->AddEntry(plot->getObject(0),"MC Non-Prompt","L");
  //legdata->AddEntry(plot->getObject(2),"MC Signal RC - matched","L");
  legdata->AddEntry(plot->getObject(1),"MC #sigma_{#eta #eta} sidebands","L");
  legdata->SetTextSize(0.035);
  legdata->SetTextFont(42);
  legdata->SetBorderSize(0);
  legdata->SetFillStyle(0);
  legdata->Draw("same");
  //CMS_lumi( ctmp,false,0 );
  ctmp->cd();
  ctmp->SetLogy(0);
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_BackgroundShape_DATAMC.pdf"));
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_BackgroundShape_DATAMC.png"));
  ctmp->SetLogy(1);
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_BackgroundShape_LOG_DATAMC.pdf"));
  ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/PFIsoPhoton_BackgroundShape_LOG_DATAMC.png"));

}


RooRealVar* ModelFit(RooWorkspace* w) {
  
  RooFitResult* fitresult;
  RooPlot* plot;
  
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);
 
  // Fit data with background pdf for data limit
  RooRealVar* ph_iso = w->var("phPurityPHiso");  
  ph_iso->setUnit("GeV");
  ph_iso->setRange("fit range", minFit, maxFit); 
  RooDataHist* data;
  data = (RooDataHist*) w->data("DataHist");
  std::cout<<"--------------> "<<data->sumEntries()<<std::endl;

  RooDataHist* dataSignal;
  dataSignal = (RooDataHist*) w->data("SignalDataHistRNDShape");
  RooDataHist* dataBkg;
  dataBkg = (RooDataHist*) w->data("BkgDataHist");
  double initFrac = 0.99;
    
  RooPlot* p = ph_iso->frame();
  data->plotOn(p);
  dataSignal->plotOn(p, RooFit::LineColor(kRed));
  dataBkg->plotOn(p, RooFit::LineColor(kBlue));
  TCanvas* cc = new TCanvas("cc", "cc",1);
  cc->cd();
  p->Draw();
  cc->SetLogy();
  cc->SaveAs("test.png");
  RooRealVar bkgYield("bkgYield","bkgYield",0.01*data->sumEntries(),0.,data->sumEntries());
  RooRealVar sigYield("sigYield","sigYield",0.99*data->sumEntries(),0., data->sumEntries());
  RooRealVar* frac = new RooRealVar("frac","frac",initFrac,0.8, 1.);

  RooHistPdf* fitFuncSig =(RooHistPdf*) w->pdf("SignalDataHistPdfRNDShape");
  RooHistPdf* fitFuncBkg =(RooHistPdf*) w->pdf("BkgDataHistPdf");
  RooAddPdf* fitFunc = new RooAddPdf("fitFunc","fitFunc",*fitFuncSig,*fitFuncBkg,*frac);
  fitresult = fitFunc->fitTo(*data, Range(minFit,maxFit),RooFit::FitOptions("MHTER"), SumW2Error(kTRUE), Save(kTRUE));   //

  fitresult->Print("V");
  w->import(*fitFunc);

  //scale considering PH Iso cut
  TTree* tree = (TTree*) w->obj("Tree");
  TH1F* h_dataRescaled = new TH1F("h_dataRescaled", "h_dataRescaled", 5*nbins, minFit-5, maxFit+5 );
  h_dataRescaled->Sumw2();
  tree->Draw(RescaledVar+TString::Format(">>h_dataRescaled"),mainCut+sieieCutSideBand+TString::Format(" && njets<=3"),"goff");// && njets==1
  std::cout<< h_dataRescaled->Integral()<<"   "<<data->sumEntries()<<std::endl;
  h_dataRescaled->Scale(data->sumEntries()/h_dataRescaled->Integral());
  
  TH1F* h_sigRescaled = new TH1F("h_sigRescaled", "h_sigRescaled", 5*nbins, minFit-5, maxFit+5 );
  h_sigRescaled->Sumw2();  
  tree->Draw(RescaledRNDVar+TString::Format(">>h_sigRescaled"),mainCut+sieieCut+TString::Format(""),"goff");
  std::cout<< h_sigRescaled->Integral()<<"   "<<dataSignal->sumEntries()<<std::endl;
  h_sigRescaled->Scale(data->sumEntries()*frac->getVal()/h_sigRescaled->Integral());
  
  TH1F* h_bkgRescaled = new TH1F("h_bkgRescaled", "h_bkgRescaled", 5*nbins, minFit-5, maxFit+5 );
  h_bkgRescaled->Sumw2();
  tree->Draw(RescaledVar+TString::Format(">>h_bkgRescaled"),mainCut+sieieCutSideBand+TString::Format(""),"goff");
  std::cout<< h_bkgRescaled->Integral()<<"   "<<dataBkg->sumEntries()<<std::endl;
  h_bkgRescaled->Scale(data->sumEntries()*(1-frac->getVal())/h_bkgRescaled->Integral());
 
 
  double dataIntErr;
  double sigIntErr;
  double bkgIntErr;
  double dataIntRedErr;
  double sigIntRedErr;
  double bkgIntRedErr;

  double dataInt = h_dataRescaled->IntegralAndError(h_dataRescaled->GetXaxis()->FindBin(-6.),h_dataRescaled->GetNbinsX(), dataIntErr);
  double sigInt = h_sigRescaled->IntegralAndError(h_dataRescaled->GetXaxis()->FindBin(-6.),h_sigRescaled->GetNbinsX(), sigIntErr);
  double bkgInt = h_bkgRescaled->IntegralAndError(h_dataRescaled->GetXaxis()->FindBin(-6.),h_bkgRescaled->GetNbinsX(),bkgIntErr);
  double dataIntRed = h_dataRescaled->IntegralAndError(h_dataRescaled->GetXaxis()->FindBin(-6.),h_dataRescaled->GetXaxis()->FindBin(PHIsoAbsCut), dataIntRedErr);
  double sigIntRed = h_sigRescaled->IntegralAndError(h_dataRescaled->GetXaxis()->FindBin(-6.),h_sigRescaled->GetXaxis()->FindBin(PHIsoAbsCut), sigIntRedErr);
  double bkgIntRed = h_bkgRescaled->IntegralAndError(h_dataRescaled->GetXaxis()->FindBin(-6.),h_bkgRescaled->GetXaxis()->FindBin(PHIsoAbsCut), bkgIntRedErr);
  
  double purity = sigInt/(sigInt+bkgInt);
  double purityErr= sqrt(sigIntErr*sigIntErr*(bkgInt*bkgInt/pow((sigInt+bkgInt),4))+ bkgIntErr*bkgIntErr*(sigInt*sigInt/pow((sigInt+bkgInt),4)));
  double purityRed = sigIntRed/(sigIntRed+bkgIntRed);
  // double purityRedDATA = sigIntRed/dataIntRed;
  double purityRedErr= sqrt(sigIntRedErr*sigIntRedErr*(bkgIntRed*bkgIntRed/pow((sigIntRed+bkgIntRed),4))+ bkgIntRedErr*bkgIntRedErr*(sigIntRed*sigIntRed/pow((sigIntRed+bkgIntRed),4)));
  std::cout<<"data: "<<dataInt<<" "<<dataIntErr<<" "<<dataIntRed<<" "<<dataIntRedErr<<std::endl;
  std::cout<<"sig: "<<sigInt<<" "<<sigIntErr<<" "<<sigIntRed<<" "<<sigIntRedErr<<std::endl;
  std::cout<<"bkg: "<<bkgInt<<" "<<bkgIntErr<<" "<<bkgIntRed<<" "<<bkgIntRedErr<<std::endl;
  std::cout<<"========================> "<<purity<<" " << purityErr<<std::endl;
  std::cout<<"========================> "<<purityRed<<" "<<purityRedErr<<std::endl;
  // std::cout<<"========================> "<<purityRedDATA<<" "<<purityRedErr<<std::endl;
  RooRealVar* RooPurity = new RooRealVar(TString::Format("purity")+Tregion+Twp+Tptcut,TString::Format("purity")+Tregion+Twp+Tptcut , purityRed);
  RooRealVar* RooNoCutPurity = new RooRealVar(TString::Format("NoCutPurity")+Tregion+Twp+Tptcut,TString::Format("NoCutpurity")+Tregion+Twp+Tptcut , purity);
  RooPurity->setError(purityRedErr);
  RooNoCutPurity->setError(purityErr);
  w->import(*RooPurity); 
  w->import(*RooNoCutPurity); 

 
  //************************************************
  // Plot  fit results
  TCanvas* ctmp = new TCanvas("ctmp","Photon Iso Fit ",1);
  Int_t nbinsx(10);
 

  h_dataRescaled->Draw("pe");
  h_sigRescaled->Add(h_bkgRescaled);
  h_sigRescaled->SetLineColor(kBlue);
  h_bkgRescaled->SetLineColor(kRed);
  h_sigRescaled->Draw("histsame");
  h_bkgRescaled->Draw("histsame");

  std::cout<< h_dataRescaled->Integral()<<" "<< h_sigRescaled->Integral()<<" "<<h_bkgRescaled->Integral()<<std::endl;
  //ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PFIsoPhoton_fit_REDUCED.pdf"));
  //ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PFIsoPhoton_fit_REDUCED.png"));
  ctmp->SetLogy(1);
  h_dataRescaled->GetYaxis()->SetRangeUser(0.01,h_data->GetMaximum()*100);
  //ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PFIsoPhoton_fit_LOG_REDUCED.pdf"));
  //ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PFIsoPhoton_fit_LOG_REDUCED.png"));
 

  ctmp->SetLogy(0);
  TH1F* h_data;// = new TH1F("h_data","h_data", nbinsx,minFit, maxFit);
  TH1F* h_pdf;
  TH1F* h_sig;
  TH1F* h_bkg;
  h_data = (TH1F*) data->createHistogram("phPurityPHiso");
  std::cout<<"--------------> "<<h_data->Integral()<<std::endl;
  h_pdf =  (TH1F*) fitFunc->createHistogram("phPurityPHiso");
  h_sig =  (TH1F*) fitFuncSig->createHistogram("phPurityPHiso");
  h_bkg =  (TH1F*) fitFuncBkg->createHistogram("phPurityPHiso");
 
  h_pdf ->Scale(h_data->Integral()/h_pdf->Integral());
  h_sig ->Scale(frac->getVal()*h_data->Integral()/h_sig->Integral());
  h_bkg ->Scale((1-frac->getVal())*h_data->Integral()/h_bkg->Integral());
  std::cout<<h_sig->GetEntries()<<" " <<h_sig->Integral()<<std::endl;
  /*h_data->Rebin(2);
  h_pdf->Rebin(2);
  h_sig->Rebin(2);
  h_bkg->Rebin(2);*/
  h_data->SetMarkerSize(0.7);
  h_pdf->SetLineColor(kBlue);
  h_sig->SetLineColor(8);
  h_bkg->SetLineColor(kRed);
  h_pdf->SetLineWidth(2);
  h_sig->SetLineWidth(2);
  h_bkg->SetLineWidth(2);

  ctmp->cd();
  
  //-------pad 1-------//
 
  TPad * pad1 = new TPad("pad1", "pad1",0.01, 0.16, 1., 1.);  
  pad1->SetRightMargin(0.1);  
  pad1->SetLogy();
  pad1->Draw();
  pad1->cd();
  
  h_data->GetXaxis()->SetTitle("PFIso - Photon [GeV]");
  h_data->GetYaxis()->SetTitle("Events");
  // h_data->GetYaxis()->SetRangeUser(0.01,220);
  // plot->Draw(); 
  
  h_data->Draw("pe");
  h_pdf->Draw("histsame"); 
  h_sig->Draw("histsame");
  h_bkg->Draw("histsame");
 

  TLegend *legdata = new TLegend(0.5334677,0.680339,0.8245968,0.8958475, "","brNDC");
  legdata->AddEntry(h_data,"Data","LP");
  legdata->AddEntry(h_sig,"Signal Component","L");
  legdata->AddEntry(h_bkg,"Background Component","L");
  legdata->AddEntry(h_pdf,"Fit Model","L");
  legdata->SetTextSize(0.035);
  legdata->SetTextFont(42);
  legdata->SetBorderSize(0);
  legdata->SetFillStyle(0);
  legdata->Draw("same");
  
  int iPos=11 ;
  //CMS_lumi( pad1,false,iPos );
    
    

    ctmp->cd();
    
    //-------pad 2------//
    TPad * pad2 = new TPad("pad2", "pad2",0.01,0.02,1.,0.23);

    pad2->SetLeftMargin(0.1789709);
    pad2->SetRightMargin(0.01565995);
    pad2->SetGrid();
    pad2->SetBottomMargin(0.4);
    pad2->SetRightMargin(0.1);
    pad2->Draw();
   
    pad2->cd();
    h_data->Sumw2();
    h_pdf->Sumw2();
    TH1F* h1_ratio1 = (TH1F*) h_data->Clone();
    TH1F* h1_ratio2 = (TH1F*)h_pdf->Clone();
    h1_ratio2->SetFillColor(kBlack);
    h1_ratio2->SetFillStyle(3003);
    h1_ratio2->SetMarkerSize(0);

    Double_t xmax =maxFit;
    Double_t xmin = minFit;
    TLine* line = new TLine(xmin,1.,xmax,1.);
  
    std::cout<<h1_ratio1->GetXaxis()->GetXmax()<< "  " << h1_ratio1->GetXaxis()->GetXmin()<<std::endl;

    h1_ratio1->SetStats(0);
    h1_ratio1->Divide(h_pdf);
    h1_ratio2->Divide(h_pdf);
    h1_ratio1->SetMarkerStyle(20);
    h1_ratio1->SetMarkerSize(0.5);
    //  h1_ratio1->GetXaxis()->SetTitle(xAxis.c_str());
    h1_ratio1->GetYaxis()->SetRangeUser(0., 2); //
    // h1_ratio1->GetXaxis()->SetRangeUser(plot->GetMinimum(),plot->GetMaximum()); //
    std::cout<<h1_ratio1->GetXaxis()->GetXmax()<< "  " << h1_ratio1->GetXaxis()->GetXmin()<<std::endl; 
    h1_ratio1->GetYaxis()->SetNdivisions(2,false);
    h1_ratio1->GetYaxis()->SetTitle("Data/Fit");
    h1_ratio1->GetYaxis()->SetTitleFont(42);
    h1_ratio1->GetXaxis()->SetTitleFont(42);
    h1_ratio1->GetXaxis()->SetTitle("PFIso - Photon [GeV]");
    
    h1_ratio1->GetXaxis()->SetTitleSize(0.18);
    h1_ratio1->GetXaxis()->SetLabelSize(0.16);
    h1_ratio1->GetXaxis()->SetLabelOffset(5.16);
    h1_ratio1->GetXaxis()->SetTitleOffset(0.6);

    h1_ratio1->GetYaxis()->SetLabelSize(0.16);
    h1_ratio1->GetYaxis()->SetTitleSize(0.18);
    h1_ratio1->GetYaxis()->SetTitleOffset(0.35);
    
    
  for(int j = 0;j<=h1_ratio1->GetNbinsX();j++){
      
      if(h_pdf->GetBinContent(j))  h1_ratio1->SetBinError(j,sqrt(pow(h_data->GetBinError(j)/h_pdf->GetBinContent(j), 2)+ pow(h_data->GetBinContent(j)*h_pdf->GetBinError(j)/(h_pdf->GetBinContent(j)*h_pdf->GetBinContent(j)),2)));
      else h1_ratio1->SetBinError(j,0.);
    }
    h1_ratio1->Draw("PEX0");
  
for(int j = 0;j<=h1_ratio2->GetNbinsX();j++){
      if(h_pdf->GetBinContent(j)) h1_ratio2->SetBinError(j,h_pdf->GetBinError(j)/h_pdf->GetBinContent(j));
      else h1_ratio2->SetBinError(j,0.);
    }
//   h1_ratio2->Draw("E2same");

 line->SetLineWidth(1.);  h_data->GetYaxis()->SetRangeUser(0.01,h_data->GetMaximum()*2.);
 line->SetLineColor(kRed);
 line->Draw("same");

 ctmp->cd();
 pad1->SetLogy(0);
 ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PFIsoPhoton_fit.pdf"));
 ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PFIsoPhoton_fit.png"));
 pad1->SetLogy(1);
 h_data->GetYaxis()->SetRangeUser(0.01,h_data->GetMaximum()*100);
 ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PFIsoPhoton_fit_LOG.pdf"));
 ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PFIsoPhoton_fit_LOG.png"));
   
 return RooPurity;
}


void PlotAllSystUnc(RooWorkspace* w,std::string region, std::string id){ //
  //ComputeAllSystUnc(w, region,id);

  TFile* fin1=  TFile::Open(TString::Format("fout_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_MCPdfRND_MCPdfNOmatch.root"));
  TFile* fin2=  TFile::Open(TString::Format("fout_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_MCPdfmatch_MCPdfNOsieie.root"));
  TFile* fin3=  TFile::Open(TString::Format("fout_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_MCPdfmatch_MCPdfNOmatch.root"));
  TFile* fin4=  TFile::Open(TString::Format("fout_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_MCPdfRND_MCPdfNOsieie.root"));
  TCanvas *c1 = (TCanvas*) fin1->Get("ctmp");	
  TCanvas *c2 = (TCanvas*) fin2->Get("ctmp");	
  TCanvas *c3 = (TCanvas*) fin3->Get("ctmp");	
  TCanvas *c4 = (TCanvas*) fin4->Get("ctmp");	
  TIter next1(c1->GetListOfPrimitives());
  TIter next2(c2->GetListOfPrimitives());
  TIter next3(c3->GetListOfPrimitives());
  TIter next4(c4->GetListOfPrimitives());

  TH1F* h_MCPdfRND_MCPdfNOmatch =(TH1F*)next1();//(TH1F*)c1->GetPrimitive("h_purity");
  TH1F* h_MCPdfmatch_MCPdfNOsieie =(TH1F*)next2(); //(TH1F*)c2->GetPrimitive("h_purity");
  TH1F* h_MCPdfmatch_MCPdfNOmatch =(TH1F*)next3();//(TH1F*)c3->GetPrimitive("h_purity");
  TH1F* h_MCPdfRND_MCPdfNOsieie =(TH1F*)next4();//(TH1F*)c3->GetPrimitive("h_purity");

  TCanvas* c = new TCanvas("c","Signal Shape",1);
  c->cd(); 
  int iPos=11 ;
  //CMS_lumi(c,false,iPos );
  TF1 *f1 = new TF1("f1", "gaus", -0.5, 0.5);
  h_MCPdfRND_MCPdfNOmatch->GetXaxis()->SetTitle("#Delta f");

  h_MCPdfRND_MCPdfNOmatch->Fit("f1");
  h_MCPdfRND_MCPdfNOmatch->Draw("hist");
  double delta1= h_MCPdfRND_MCPdfNOmatch->GetMean();
  double mean1=f1->GetParameter(1);
  double RMS1=f1->GetParameter(2);
  c->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PurityToy_MCPdfRND_MCPdfNOmatch.png"));
  c->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PurityToy_MCPdfRND_MCPdfNOmatch.pdf"));

  h_MCPdfmatch_MCPdfNOsieie->GetXaxis()->SetTitle("#Delta f");
 
  h_MCPdfmatch_MCPdfNOsieie->Fit("f1");
  h_MCPdfmatch_MCPdfNOsieie->Draw("hist");
  double delta2= h_MCPdfmatch_MCPdfNOsieie->GetMean();
  double mean2=f1->GetParameter(1);
  double RMS2=f1->GetParameter(2);
  c->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PurityToy_MCPdfmatch_MCPdfNOsieie.png"));
  c->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PurityToy_MCPdfmatch_MCPdfNOsieie.pdf"));
  
  h_MCPdfmatch_MCPdfNOmatch->GetXaxis()->SetTitle("#Delta f");
 
  h_MCPdfmatch_MCPdfNOmatch->Fit("f1");
  h_MCPdfmatch_MCPdfNOmatch->Draw("hist");
  double delta3= h_MCPdfmatch_MCPdfNOmatch->GetMean();

  double mean3=f1->GetParameter(1);
  double RMS3=f1->GetParameter(2);
  c->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PurityToy_MCPdfmatch_MCPdfNOmatch.png"));
  c->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PurityToy_MCPdfmatch_MCPdfNOmatch.pdf"));


  h_MCPdfRND_MCPdfNOsieie->GetXaxis()->SetTitle("#Delta f");
 
  h_MCPdfRND_MCPdfNOsieie->Fit("f1");
  h_MCPdfRND_MCPdfNOsieie->Draw("hist");
  double delta4= h_MCPdfRND_MCPdfNOsieie->GetMean();

  double mean4=f1->GetParameter(1);
  double RMS4=f1->GetParameter(2);
  c->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PurityToy_MCPdfRND_MCPdfNOsieie.png"));
  c->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_PurityToy_MCPdfRND_MCPdfNOsieie.pdf"));

 
  // std::cout<<mean1<<"  "<<mean2<<" "<<mean3<<std::endl;
  double deltaMax=fabs(mean1);
  //if(fabs(mean2)>deltaMax)deltaMax=fabs(mean2);
  if(fabs(mean3)>deltaMax)deltaMax=fabs(mean3);
  //if(fabs(mean4)>deltaMax)deltaMax=fabs(mean4);
  
  std::cout<<"MCPdfRND_MCPdfNOmatch: "<<mean1<< " +/- "<<RMS1<<std::endl;
  std::cout<<"MCPdfmatch_MCPdfNOsieie: "<<mean2<< " +/- "<<RMS2<<std::endl;
  std::cout<<"MCPdfRNDh_MCPdfNOsieie: "<<mean4<< " +/- "<<RMS4<<std::endl;
  std::cout<<"MCPdfmatch_MCPdfNOmatch: "<<mean3<< " +/- "<<RMS3<<std::endl;
  std::cout<<deltaMax<<std::endl;

  RooRealVar* puritySyst = new RooRealVar(TString::Format("puritySyst")+Tregion+Twp+Tptcut,TString::Format("puritySyst")+Tregion+Twp+Tptcut, deltaMax);
  w->import(*puritySyst); 

}
void ComputeAllSystUnc(RooWorkspace* w, std::string region, std::string id){
  AddBkgMC(w);
  std::cout<<"flag"<<stsd::endl;
  AddSignalMC(w);
  AddSignalMCRND(w);
  std::cout<<"flag"<<stsd::endl;
  
  std::cout<<"flag"<<stsd::endl;
  RooHistPdf* MCPdfRND =(RooHistPdf*) w->pdf("SignalMCHistPdfRNDShape");
  RooHistPdf* MCPdfmatch =(RooHistPdf*) w->pdf("SignalMCHistPdfmatch");
  RooHistPdf* MCPdfNOmatch =(RooHistPdf*) w->pdf("BkgMCHistPdfNOmatch");
  RooHistPdf* MCPdfNOsieie =(RooHistPdf*) w->pdf("BkgMCHistPdfNOsieie");

  TString Tsig1="MCPdfRND";
  TString Tbkg1="MCPdfNOmatch";
  ComputeSystUnc(w,  region, id,MCPdfRND,MCPdfNOsieie,MCPdfRND,MCPdfNOmatch, Tsig1, Tbkg1);
  TString Tsig2="MCPdfmatch";
  TString Tbkg2="MCPdfNOsieie";
  ComputeSystUnc(w,  region,  id, MCPdfRND,MCPdfNOsieie,MCPdfmatch,MCPdfNOsieie, Tsig2, Tbkg2);
  TString Tsig3="MCPdfmatch";
  TString Tbkg3="MCPdfNOmatch";
  ComputeSystUnc(w,  region,  id, MCPdfRND,MCPdfNOsieie,MCPdfmatch,MCPdfNOmatch, Tsig3, Tbkg3);
  TString Tsig4="MCPdfRND";
  TString Tbkg4="MCPdfNOsieie";
  ComputeSystUnc(w,  region,  id,MCPdfmatch,MCPdfNOmatch, MCPdfRND,MCPdfNOsieie, Tsig4, Tbkg4);
  
}



void ComputeSystUnc(RooWorkspace*w, std::string region, std::string id, RooHistPdf* genSig, RooHistPdf* genBkg, RooHistPdf* fitSig, RooHistPdf* fitBkg, TString Tsig, TString Tbkg){
  RooFitResult* res;
  TCanvas* ctmp = new TCanvas("ctmp","Photon Iso Fit ",1);
  ctmp->cd();
  RooRealVar* purity = w->var(TString::Format("NoCutPurity")+Tregion+Twp+Tptcut);
  TLatex *text = new TLatex();
  text->SetNDC();
  text->SetTextSize(0.04);

  // Fit data with background pdf for data limit
  RooRealVar* ph_iso = w->var("phPurityPHiso");  
  ph_iso->setUnit("GeV");
  ph_iso->setRange("fit range", minFit, maxFit); 
 
  RooDataHist* data;
  data = (RooDataHist*) w->data("DataHist");
  std::cout<<"--------------> "<<data->sumEntries()<<std::endl;
  double totEvents = data->sumEntries();
  RooPlot* plot= ph_iso->frame();

  //generate: MCPdfRND+MCPdfNosieie fit: MCPdfRND+MCPdfNOmatch
  double frac_base=purity->getVal();
  RooRealVar* rooPurity = new RooRealVar("rooPurity","rooPurity",purity->getVal());
  RooRealVar* frac = new RooRealVar("frac","frac",frac_base, 0.8, 1);
  RooAddPdf* templateFunc = new RooAddPdf("templateFunc","templateFunc",*genSig,*genBkg,*rooPurity);
  RooAddPdf* fitFunc = new RooAddPdf("fitFunc","fitFunc",*fitSig,*fitBkg,*frac);
  TH1F*  h_data = (TH1F*) templateFunc->createHistogram("phPurityPHiso");
  h_data->Scale(totEvents/h_data->Integral());
  RooDataSet* genDataUnBin;
  RooDataHist* genData;
  RooRandom::randomGenerator()->SetSeed(0);
  double systPurity;
  TH1F* h_purity= new TH1F(" h_purity", "h_purity", 200, -0.2, 0.2);
  for (int iToy=0; iToy<1000; iToy++){
  //for (int iToy=0; iToy<2; iToy++){
    systPurity=999.;
    frac->setVal(purity->getVal());
    cout << "----------------------------------------------" << endl;
    cout << "------------- TOY: " << iToy << "------------" << endl;
    cout << "----------------------------------------------" << endl;
    std::cout<<totEvents<<std::endl;
    genData = (RooDataHist*) templateFunc->generate(*ph_iso,totEvents, kFALSE,kFALSE);
    //genData = new RooDataHist("gen","gen",*ph_iso,*genDataUnBin);
  
    res = fitFunc->fitTo(*genData, Range(minFit,maxFit),SumW2Error(kTRUE), Save(kTRUE));   //
    res->Print("v");
    bool makePlot = true;
    if(makePlot && (iToy==0||iToy==10 ||iToy==20||iToy==100 ||iToy==500)){
     
      TH1F* h_pdf;
      TH1F* h_sig;
      TH1F* h_bkg;

      // std::cout<<"--------------> "<<genDataUnBin->sumEntries()<<std::endl;
      h_pdf =  (TH1F*) fitFunc->createHistogram("phPurityPHiso");
      h_pdf->Scale(h_data->Integral()/h_pdf->Integral()); 
      h_data->SetMarkerSize(0.7);
      h_pdf->SetLineColor(kBlue);
      h_pdf->SetLineWidth(2);
      h_data->GetXaxis()->SetTitle("PFIso - Photon [GeV]");
      h_data->GetYaxis()->SetTitle("Events");
      //   h_data->GetYaxis()->SetRangeUser(0.01,220);
      //
      h_data->Draw("pe");
      h_pdf->Draw("histsame"); 
      TLegend *legdata = new TLegend(0.5334677,0.680339,0.8245968,0.8958475, "","brNDC");
      legdata->AddEntry(h_data,"Pseudo-Data","LP");
      legdata->AddEntry(h_pdf,"Fit Model","L");
      legdata->SetTextSize(0.035);
      legdata->SetTextFont(42);
      legdata->SetBorderSize(0);
      legdata->SetFillStyle(0);
      legdata->Draw("same");
  
      int iPos=11 ;
      //     CMS_lumi(ctmp,false,iPos );
      ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/")+Tsig+TString::Format("_")+Tbkg+TString::Format("_PFIsoPhoton_fit_Toy%d.pdf", iToy));
      ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/")+Tsig+TString::Format("_")+Tbkg+TString::Format("_PFIsoPhoton_fit_Toy%d.png", iToy));
      ctmp->SetLogy(1);
      h_data->GetYaxis()->SetRangeUser(0.01,h_data->GetMaximum()*100);
      ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/")+Tsig+TString::Format("_")+Tbkg+TString::Format("_PFIsoPhoton_fit_LOG_Toy%d.pdf", iToy));
      ctmp->SaveAs(TString::Format("plots_")+Tregion+TString::Format("_")+Twp+TString::Format("/")+Tsig+TString::Format("_")+Tbkg+TString::Format("_PFIsoPhoton_fit_LOG_Toy%d.png", iToy));

    }


    systPurity=((RooRealVar *)(res->floatParsFinal().find("frac")))->getVal();
    // double systPurityErr=((RooRealVar *)(res->floatParsFinal().find("frac")))->getError();
    std::cout<<systPurity<<std::endl;
    if(systPurity>0.5 && systPurity<1.2)h_purity->Fill(systPurity-purity->getVal());   

  }
 
  TFile* fout =  new TFile(TString::Format("fout_")+Tregion+TString::Format("_")+Twp+TString::Format("_")+Tptcut+TString::Format("_")+Tsig+TString::Format("_")+Tbkg+TString::Format(".root"), "RECREATE");
  TCanvas* ctmp = new TCanvas("ctmp","Signal Shape",1);
  ctmp->cd();
  h_purity->Draw();
  fout->cd();
  ctmp->Write();
 
  fout->Close();

}
