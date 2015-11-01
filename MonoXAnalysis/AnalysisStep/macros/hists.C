#include "makehist.h"

int nbins  = 7;
float bins[]  = {200., 250., 300., 350., 400., 500., 600., 1000.};

int nrbins = 7;
float rbins[] = {200., 250., 300., 350., 400., 500., 600., 1000.};

void makezmmcorhist(string ext="", TH1* khist=NULL) {
    TFile  nfile("/Users/avartak/CMS/MonoX/Trees/znn100toinf/sigtree.root");     
    TFile  dfile("/Users/avartak/CMS/MonoX/Trees/zll100toinf/zmmtree.root");  

    TTree* ntree = (TTree*)nfile.Get("tree/tree");  
    TTree* dtree = (TTree*)dfile.Get("tree/tree");  

    TH1F nhist("nhist", "", nrbins, rbins);
    TH1F dhist("dhist", "", nrbins, rbins);

    TFile* sffile = new TFile("leptonIDsfs.root");
    TH2*  sflhist = (TH2*)sffile->Get("muon_loose_SF");
    TH2*  sfthist = (TH2*)sffile->Get("muon_tight_SF");

    makehist(ntree, &nhist,  true, 0, 1.00, sflhist, sfthist, NULL, khist);
    makehist(dtree, &dhist,  true, 1, 1.00, sflhist, sfthist, NULL, khist);

    string name = string("zmmcor")+ext;    

    nhist.Divide(&dhist);
    TFile outfile((name+".root").c_str(), "RECREATE");
    nhist.SetName((name+"hist" ).c_str());
    nhist.Write();
    outfile.Close();
}

void makezeecorhist(string ext="", TH1* khist=NULL) {
    TFile  nfile("/Users/avartak/CMS/MonoX/Trees/znn100toinf/sigtree.root");
    TFile  dfile("/Users/avartak/CMS/MonoX/Trees/zll100toinf/zeetree.root");

    TTree* ntree = (TTree*)nfile.Get("tree/tree");
    TTree* dtree = (TTree*)dfile.Get("tree/tree");

    TH1F nhist("nhist", "", nrbins, rbins);
    TH1F dhist("dhist", "", nrbins, rbins);

    TFile* sffile = new TFile("leptonIDsfs.root");
    TH2*  sflhist = (TH2*)sffile->Get("electron_veto_SF");
    TH2*  sfthist = (TH2*)sffile->Get("electron_tight_SF");

    makehist(ntree, &nhist,  true, 0, 1.00, sflhist, sfthist, NULL, khist);
    makehist(dtree, &dhist,  true, 3, 1.00, sflhist, sfthist, NULL, khist);

    string name = string("zeecor")+ext;    

    nhist.Divide(&dhist);
    TFile outfile((name+".root").c_str(), "RECREATE");
    nhist.SetName((name+"hist" ).c_str());
    nhist.Write();
    outfile.Close();
}

void makewmncorhist(string ext="", TH1* khist=NULL) {
    TFile  nfile("/Users/avartak/CMS/MonoX/Trees/wln100toinf/sigtree.root");
    TFile  dfile("/Users/avartak/CMS/MonoX/Trees/wln100toinf/wmntree.root");

    TTree* ntree = (TTree*)nfile.Get("tree/tree");
    TTree* dtree = (TTree*)dfile.Get("tree/tree");

    TH1F nhist("nhist", "", nrbins, rbins);
    TH1F dhist("dhist", "", nrbins, rbins);

    TFile* sffile = new TFile("leptonIDsfs.root");
    TH2*  sflhist = (TH2*)sffile->Get("muon_loose_SF");
    TH2*  sfthist = (TH2*)sffile->Get("muon_tight_SF");

    makehist(ntree, &nhist,  true, 0, 1.00, sflhist, sfthist, NULL, khist);
    makehist(dtree, &dhist,  true, 2, 1.00, sflhist, sfthist, NULL, khist);

    string name = string("wmncor")+ext;    

    nhist.Divide(&dhist);
    TFile outfile((name+".root").c_str(), "RECREATE");
    nhist.SetName((name+"hist" ).c_str());
    nhist.Write();
    outfile.Close();
}

void makewzmcorhist(string ext="", TH1* zkhist=NULL, TH1* wkhist=NULL) {
    TFile  nfile("/Users/avartak/CMS/MonoX/Trees/znn100toinf/sigtree.root");
    TFile  dfile("/Users/avartak/CMS/MonoX/Trees/wln100toinf/wmntree.root");

    TTree* ntree = (TTree*)nfile.Get("tree/tree");
    TTree* dtree = (TTree*)dfile.Get("tree/tree");

    TH1F nhist("nhist", "", nrbins, rbins);
    TH1F dhist("dhist", "", nrbins, rbins);

    TFile* sffile = new TFile("leptonIDsfs.root");
    TH2*  sflhist = (TH2*)sffile->Get("muon_loose_SF");
    TH2*  sfthist = (TH2*)sffile->Get("muon_tight_SF");

    makehist(ntree, &nhist,  true, 0, 1.00, sflhist, sfthist, NULL, zkhist);
    makehist(dtree, &dhist,  true, 2, 1.00, sflhist, sfthist, NULL, wkhist);

    string name = string("wzmcor")+ext;    

    nhist.Divide(&dhist);
    TFile outfile((name+".root").c_str(), "RECREATE");
    nhist.SetName((name+"hist" ).c_str());
    nhist.Write();
    outfile.Close();
}

void makegamcorhist(string ext="", TH1* zkhist=NULL, TH1* gkhist=NULL) {
    TFile  nfile("/Users/avartak/CMS/MonoX/Trees/znn100toinf/sigtree.root");
    TFile  dfile("/Users/avartak/CMS/MonoX/Trees/gam100toinf/gamtree.root");

    TTree* ntree = (TTree*)nfile.Get("tree/tree");
    TTree* dtree = (TTree*)dfile.Get("tree/tree");

    TH1F nhist("nhist", "", nrbins, rbins);
    TH1F dhist("dhist", "", nrbins, rbins);

    makehist(ntree, &nhist,  true, 0, 1.00, NULL, NULL, NULL, zkhist);
    makehist(dtree, &dhist,  true, 5, 1.00, NULL, NULL, NULL, gkhist);

    string name = string("gamcor")+ext;    

    nhist.Divide(&dhist);
    TFile outfile((name+".root").c_str(), "RECREATE");
    nhist.SetName((name+"hist" ).c_str());
    nhist.Write();
    outfile.Close();
}

void sigdatamchist(TFile* outfile, bool blind=false) {

    TFile znfile("/Users/avartak/CMS/MonoX/Trees/znn100toinf/sigtree.root");
    TFile wlfile("/Users/avartak/CMS/MonoX/Trees/wln100toinf/sigtree.root");
    TFile zlfile("/Users/avartak/CMS/MonoX/Trees/zll100toinf/sigtree.root");
    TFile ttfile("/Users/avartak/CMS/MonoX/Trees/top/sigtree.root");
    TFile difile("/Users/avartak/CMS/MonoX/Trees/qcd/sigtree.root");
    TFile qcfile("/Users/avartak/CMS/MonoX/Trees/dibosons/sigtree.root");
    TFile sifile("/Users/avartak/CMS/MonoX/Trees/avM10000m1/sigtree.root");
    TFile dtfile("/Users/avartak/CMS/MonoX/Trees/metD/sigtree.root");

    TH1F znhist("zinvhist", "", nbins, bins);
    TH1F wlhist("wjethist", "", nbins, bins);
    TH1F zlhist("zjethist", "", nbins, bins);
    TH1F tthist("tbkghist", "", nbins, bins);
    TH1F dihist("dbkghist", "", nbins, bins);
    TH1F qchist("qbkghist", "", nbins, bins);
    TH1F sihist("sig1hist", "", nbins, bins);
    TH1F dthist("datahist", "", nbins, bins);

    TTree* zntree = (TTree*)znfile.Get("tree/tree");
    TTree* wltree = (TTree*)wlfile.Get("tree/tree");
    TTree* zltree = (TTree*)zlfile.Get("tree/tree");
    TTree* tttree = (TTree*)ttfile.Get("tree/tree");
    TTree* ditree = (TTree*)difile.Get("tree/tree");
    TTree* qctree = (TTree*)qcfile.Get("tree/tree");
    TTree* sitree = (TTree*)sifile.Get("tree/tree");
    TTree* dttree = (TTree*)dtfile.Get("tree/tree");

    makehist(zntree, &znhist,  true, 0, 1.23, NULL, NULL, NULL, NULL);
    makehist(wltree, &wlhist,  true, 0, 1.21, NULL, NULL, NULL, NULL);
    makehist(zltree, &zlhist,  true, 0, 1.23, NULL, NULL, NULL, NULL);
    makehist(tttree, &tthist,  true, 0, 1.00, NULL, NULL, NULL, NULL);
    makehist(ditree, &dihist,  true, 0, 1.00, NULL, NULL, NULL, NULL);
    makehist(qctree, &qchist,  true, 0, 1.00, NULL, NULL, NULL, NULL);
    makehist(sitree, &sihist,  true, 0, 1.00, NULL, NULL, NULL, NULL);
    makehist(dttree, &dthist, false, 0, 1.00, NULL, NULL, NULL, NULL);

    sihist.Scale(1e5);

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
}

void gamdatamchist(TFile* outfile) {

    TFile dtfile("/Users/avartak/CMS/MonoX/Trees/singlephD/gamtree.root");

    TH1F dthist("datahistgam", "", nbins, bins);
    TH1F qchist("qbkghistgam", "", nbins, bins);

    TTree* dttree = (TTree*)dtfile.Get("tree/tree");

    makehist(dttree, &dthist, false, 5, 1.00, NULL, NULL, NULL, NULL);

    for (int i = 1; i <= dthist.GetNbinsX(); i++) qchist.SetBinContent(i, dthist.GetBinContent(i)*0.05);

    outfile->cd();
    dthist.Write();
    qchist.Write();

}

void lepdatamchist(TFile* outfile, int sample) {

    if (sample != 1 && sample != 2 && sample != 3) return;

    string filename;
    if      (sample == 1)   filename = "zmmtree.root";
    else if (sample == 3)   filename = "zeetree.root";
    else if (sample == 2)   filename = "wmntree.root";

    string vlfilename;
    if      (sample == 1) vlfilename = "wln100toinf/";
    else if (sample == 3) vlfilename = "wln100toinf/";
    else if (sample == 2) vlfilename = "zll100toinf/";

    string dtfilename;
    if      (sample == 1) dtfilename = "metD/";
    else if (sample == 3) dtfilename = "singleelD/";
    else if (sample == 2) dtfilename = "metD/";

    string suffix;
    if      (sample == 1) suffix = "zmm";
    else if (sample == 3) suffix = "zee";
    else if (sample == 2) suffix = "wmn";

    TFile ttfile((string("/Users/avartak/CMS/MonoX/Trees/top/")          + filename).c_str());
    TFile difile((string("/Users/avartak/CMS/MonoX/Trees/dibosons/")     + filename).c_str());
    TFile qcfile((string("/Users/avartak/CMS/MonoX/Trees/qcd/")          + filename).c_str());
    TFile vlfile((string("/Users/avartak/CMS/MonoX/Trees/") + vlfilename + filename).c_str());
    TFile dtfile((string("/Users/avartak/CMS/MonoX/Trees/") + dtfilename + filename).c_str());

    TH1F dthist((string("datahist")+suffix).c_str(), "", nbins, bins);
    TH1F vlhist((string("lbkghist")+suffix).c_str(), "", nbins, bins);
    TH1F tthist((string("tbkghist")+suffix).c_str(), "", nbins, bins);
    TH1F dihist((string("dbkghist")+suffix).c_str(), "", nbins, bins);
    TH1F qchist((string("qbkghist")+suffix).c_str(), "", nbins, bins);

    TTree* dttree = (TTree*)dtfile.Get("tree/tree");
    TTree* vltree = (TTree*)vlfile.Get("tree/tree");
    TTree* tttree = (TTree*)ttfile.Get("tree/tree");
    TTree* ditree = (TTree*)difile.Get("tree/tree");
    TTree* qctree = (TTree*)qcfile.Get("tree/tree");

    makehist(dttree, &dthist, false, sample, 1.00, NULL, NULL, NULL, NULL);
    makehist(vltree, &vlhist,  true, sample, 1.23, NULL, NULL, NULL, NULL);
    makehist(tttree, &tthist,  true, sample, 1.00, NULL, NULL, NULL, NULL);
    makehist(ditree, &dihist,  true, sample, 1.00, NULL, NULL, NULL, NULL);
    makehist(qctree, &qchist,  true, sample, 1.00, NULL, NULL, NULL, NULL);

    outfile->cd();
    dthist.Write();
    tthist.Write();
    dihist.Write();
    qchist.Write();
    vlhist.Write();

}

void hists() {

    TFile* zkfile = new TFile("zkfactor.root");
    TFile* gkfile = new TFile("gammakfactor.root");

    TH1F*  zkhist = (TH1F*)zkfile->Get("zkfactor"); 
    TH1F*  gkhist = (TH1F*)zkfile->Get("gammakfactor"); 

    makezmmcorhist();    
    makezeecorhist();
    makewmncorhist();    
    makewzmcorhist();
    makegamcorhist();
    makegamcorhist("nlo", zkhist, gkhist);

    TFile* zmmcorfile = new TFile("zmmcor.root");
    TFile* zeecorfile = new TFile("zeecor.root");
    TFile* wmncorfile = new TFile("wmncor.root");
    TFile* wzmcorfile = new TFile("wzmcor.root");
    TFile* gamcorfile = new TFile("gamcor.root");

    TFile* gamcornlofile = new TFile("gamcornlo.root");

    TH1* zmmcorhist = (TH1*)zmmcorfile->Get("zmmcorhist");    
    TH1* zeecorhist = (TH1*)zeecorfile->Get("zeecorhist");    
    TH1* wmncorhist = (TH1*)wmncorfile->Get("wmncorhist");    
    TH1* wzmcorhist = (TH1*)wzmcorfile->Get("wzmcorhist");    
    TH1* gamcorhist = (TH1*)gamcorfile->Get("gamcorhist");    

    TH1* gamcornlohist = (TH1*)gamcornlofile->Get("gamcornlohist");    
    gamcornlohist->Divide(gamcorhist);
    for (int i = 1; i <= gamcornlohist->GetNbinsX(); i++) gamcornlohist->SetBinContent(i, gamcornlohist->GetBinContent(i)-1.0);
    gamcornlohist->SetName("Theory");

    TFile outfile("templates.root", "RECREATE");

    zmmcorhist->Write();
    zeecorhist->Write();
    wmncorhist->Write();
    wzmcorhist->Write();
    gamcorhist->Write();
    gamcornlohist->Write();

    sigdatamchist(&outfile);
    gamdatamchist(&outfile);
    lepdatamchist(&outfile, 1);
    lepdatamchist(&outfile, 2);
    lepdatamchist(&outfile, 3);

    outfile.Close();

}
