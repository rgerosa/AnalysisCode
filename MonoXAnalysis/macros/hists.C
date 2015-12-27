#include "makehist4.h"

int nbins  = 7;
float bins[]  = {200., 250., 300., 350., 400., 500., 600., 1000.};

int nrbins = 7;
float rbins[] = {200., 250., 300., 350., 400., 500., 600., 1000.};

void makezmmcorhist() {
    TFile  nfile("/Users/avartak/CMS/MonoX/FinalTrees/znn100toinf/sigtree.root");     
    TFile  dfile("/Users/avartak/CMS/MonoX/FinalTrees/zll100toinf/zmmtree.root");  

    TTree* ntree = (TTree*)nfile.Get("tree");  
    TTree* dtree = (TTree*)dfile.Get("tree");  

    TH1F nhist("nhist", "", nrbins, rbins);
    TH1F dhist("dhist", "", nrbins, rbins);

    TFile kffile("scalefactors_v4.root");
    TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
    TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
    TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");

    znlohist->Divide(zlohist);

    vector<TH1*> ehists;
    vector<TH1*> zhists;
    zhists.push_back(znlohist);
    zhists.push_back(zewkhist);

    makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
    makehist4(dtree, &dhist,  true, 1, 1.00, zhists, NULL);

    string name = string("zmmcor");    

    nhist.Divide(&dhist);
    TFile outfile((name+".root").c_str(), "RECREATE");
    nhist.SetName((name+"hist" ).c_str());
    nhist.Write();
    outfile.Close();

    nfile.Close();
    dfile.Close();

    cout << "Z(mumu)->Z(inv) transfer factor computed ..." << endl;
}

void makezeecorhist() {
    TFile  nfile("/Users/avartak/CMS/MonoX/FinalTrees/znn100toinf/sigtree.root");
    TFile  dfile("/Users/avartak/CMS/MonoX/FinalTrees/zll100toinf/zeetree.root");

    TTree* ntree = (TTree*)nfile.Get("tree");
    TTree* dtree = (TTree*)dfile.Get("tree");

    TH1F nhist("nhist", "", nrbins, rbins);
    TH1F dhist("dhist", "", nrbins, rbins);

    TFile kffile("scalefactors_v4.root");
    TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
    TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
    TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");

    znlohist->Divide(zlohist);

    vector<TH1*> ehists;
    vector<TH1*> zhists;
    zhists.push_back(znlohist);
    zhists.push_back(zewkhist);

    makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
    makehist4(dtree, &dhist,  true, 3, 1.00, zhists, NULL);

    string name = string("zeecor");    

    nhist.Divide(&dhist);
    TFile outfile((name+".root").c_str(), "RECREATE");
    nhist.SetName((name+"hist" ).c_str());
    nhist.Write();
    outfile.Close();

    nfile.Close();
    dfile.Close();

    cout << "Z(ee)->Z(inv) transfer factor computed ..." << endl;
}

void makewmncorhist() {
    TFile  nfile("/Users/avartak/CMS/MonoX/FinalTrees/wln100toinf/sigtree.root");
    TFile  dfile("/Users/avartak/CMS/MonoX/FinalTrees/wln100toinf/wmntree.root");

    TTree* ntree = (TTree*)nfile.Get("tree");
    TTree* dtree = (TTree*)dfile.Get("tree");

    TH1F nhist("nhist", "", nrbins, rbins);
    TH1F dhist("dhist", "", nrbins, rbins);

    TFile kffile("scalefactors_v4.root");
    TH1* znlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
    TH1*  zlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
    TH1* zewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");

    znlohist->Divide(zlohist);

    vector<TH1*> ehists;
    vector<TH1*> zhists;
    zhists.push_back(znlohist);
    zhists.push_back(zewkhist);


    makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
    makehist4(dtree, &dhist,  true, 2, 1.00, zhists, NULL);

    string name = string("wmncor");    

    nhist.Divide(&dhist);
    TFile outfile((name+".root").c_str(), "RECREATE");
    nhist.SetName((name+"hist" ).c_str());
    nhist.Write();
    outfile.Close();

    nfile.Close();
    dfile.Close();

    cout << "W(munu)->W+jets transfer factor computed ..." << endl;
}

void makewencorhist(string ext="", int pdf=0) {
    TFile  nfile("/Users/avartak/CMS/MonoX/FinalTrees/wln100toinf/sigtree.root");
    TFile  dfile("/Users/avartak/CMS/MonoX/FinalTrees/wln100toinf/wentree.root");

    TTree* ntree = (TTree*)nfile.Get("tree");
    TTree* dtree = (TTree*)dfile.Get("tree");

    TH1F nhist("nhist", "", nrbins, rbins);
    TH1F dhist("dhist", "", nrbins, rbins);

    TFile kffile("scalefactors_v4.root");
    TH1* znlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
    TH1*  zlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
    TH1* zewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");

    znlohist->Divide(zlohist);

    vector<TH1*> ehists;
    vector<TH1*> zhists;
    zhists.push_back(znlohist);
    zhists.push_back(zewkhist);

    makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
    makehist4(dtree, &dhist,  true, 4, 1.00, zhists, NULL);

    string name = string("wencor")+ext;

    nhist.Divide(&dhist);
    TFile outfile((name+".root").c_str(), "RECREATE");
    nhist.SetName((name+"hist" ).c_str());
    nhist.Write();
    outfile.Close();

    nfile.Close();
    dfile.Close();

    cout << "W(enu)->W+jets transfer factor computed ..." << endl;
}

void makewzmcorhist() {
    TFile  nfile("/Users/avartak/CMS/MonoX/FinalTrees/znn100toinf/sigtree.root");
    TFile  dfile("/Users/avartak/CMS/MonoX/FinalTrees/wln100toinf/wmntree.root");

    TTree* ntree = (TTree*)nfile.Get("tree");
    TTree* dtree = (TTree*)dfile.Get("tree");

    TH1F nhist("nhist", "", nrbins, rbins);
    TH1F dhist("dhist", "", nrbins, rbins);

    TFile kffile("scalefactors_v4.root");
    TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
    //TH1* znlohist = (TH1*)kffile.Get("znlo1/znlo1_nominal");
    TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
    TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");

    TH1* wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
    //TH1* wnlohist = (TH1*)kffile.Get("wnlo1/wnlo1_nominal");
    TH1*  wlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
    TH1* wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");

    znlohist->Divide(zlohist);
    wnlohist->Divide(wlohist);

    vector<TH1*> zhists;
    vector<TH1*> whists;
    zhists.push_back(znlohist);
    zhists.push_back(zewkhist);
    whists.push_back(wnlohist);
    whists.push_back(wewkhist);

    makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
    makehist4(dtree, &dhist,  true, 2, 1.00, whists, NULL);

    string name = string("wzmcor");    

    nhist.Divide(&dhist);
    TFile outfile((name+".root").c_str(), "RECREATE");
    nhist.SetName((name+"hist" ).c_str());
    nhist.Write();
    outfile.Close();

    nfile.Close();
    dfile.Close();

    cout << "W(munu)->Z(inv) transfer factor computed ..." << endl;
}

void makezwjcorhist(string ext="", int kfact=0) {
    TFile  nfile("/Users/avartak/CMS/MonoX/FinalTrees/znn100toinf/sigtree.root");
    TFile  dfile("/Users/avartak/CMS/MonoX/FinalTrees/wln100toinf/sigtree.root");

    TTree* ntree = (TTree*)nfile.Get("tree");
    TTree* dtree = (TTree*)dfile.Get("tree");

    TH1F nhist("nhist", "", nrbins, rbins);
    TH1F dhist("dhist", "", nrbins, rbins);

    TFile kffile("scalefactors_v4.root");
    TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
    TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
    TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");
    TH1* zpdfhist = (TH1*)kffile.Get("znlo012/znlo012_pdfUp");

    TH1* wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
    TH1*  wlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
    TH1* wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");
    TH1* wpdfhist = (TH1*)kffile.Get("wnlo012/wnlo012_pdfUp");

    TH1* nomhist  = (TH1*)kffile.Get("znlo1_over_wnlo1/znlo1_over_wnlo1");
    TH1* re1hist  = (TH1*)kffile.Get("znlo1_over_wnlo1/znlo1_over_wnlo1_renCorrUp");
    TH1* re2hist  = (TH1*)kffile.Get("znlo1_over_wnlo1/znlo1_over_wnlo1_renAcorrUp");
    TH1* fa1hist  = (TH1*)kffile.Get("znlo1_over_wnlo1/znlo1_over_wnlo1_facCorrUp");
    TH1* fa2hist  = (TH1*)kffile.Get("znlo1_over_wnlo1/znlo1_over_wnlo1_facAcorrUp");

    zpdfhist->Divide(znlohist);
    wpdfhist->Divide(wnlohist);

    re1hist->Divide(nomhist);
    re2hist->Divide(nomhist);
    fa1hist->Divide(nomhist);
    fa2hist->Divide(nomhist);

    znlohist->Divide(zlohist);
    wnlohist->Divide(wlohist);

    vector<TH1*> zhists;
    vector<TH1*> whists;
    if (kfact == 1) zhists.push_back(znlohist);
    if (kfact == 2) {zhists.push_back(znlohist); zhists.push_back(zewkhist);}
    if (kfact == 3) {zhists.push_back(znlohist); zhists.push_back(re1hist) ;}
    if (kfact == 4) {zhists.push_back(znlohist); zhists.push_back(fa1hist) ;}
    if (kfact == 5) {zhists.push_back(znlohist); zhists.push_back(re2hist) ;}
    if (kfact == 6) {zhists.push_back(znlohist); zhists.push_back(fa2hist) ;}
    if (kfact == 7) {zhists.push_back(znlohist); zhists.push_back(zpdfhist);}

    if (kfact == 1) whists.push_back(wnlohist);
    if (kfact == 2) {whists.push_back(wnlohist); whists.push_back(wewkhist);}
    if (kfact == 3) whists.push_back(wnlohist);
    if (kfact == 4) whists.push_back(wnlohist);
    if (kfact == 5) whists.push_back(wnlohist);
    if (kfact == 6) whists.push_back(wnlohist);
    if (kfact == 7) {whists.push_back(wnlohist); whists.push_back(wpdfhist);}


    makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
    makehist4(dtree, &dhist,  true, 0, 1.00, whists, NULL);

    string name = string("zwjcor")+ext;

    nhist.Divide(&dhist);
    TFile outfile((name+".root").c_str(), "RECREATE");
    nhist.SetName((name+"hist" ).c_str());
    nhist.Write();
    outfile.Close();

    nfile.Close();
    dfile.Close();
    kffile.Close();

    cout << "W+jets->Z(inv) transfer factor computed ..." << endl;
}

void makegamcorhist(string ext="", int kfact=0) {
    TFile  nfile("/Users/avartak/CMS/MonoX/FinalTrees/znn100toinf/sigtree.root");
    TFile  dfile("/Users/avartak/CMS/MonoX/FinalTrees/gam100toinf/gamtree.root");
    //TFile  dfile("/Users/avartak/CMS/MonoX/FinalTrees/gam100toinf/loosegamtree.root");

    TTree* ntree = (TTree*)nfile.Get("tree");
    TTree* dtree = (TTree*)dfile.Get("tree");

    TH1F nhist("nhist", "", nrbins, rbins);
    TH1F dhist("dhist", "", nrbins, rbins);

    TFile kffile("scalefactors_v4.root");
    TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
    TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
    TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");
    TH1* zpdfhist = (TH1*)kffile.Get("znlo012/znlo012_pdfUp");

    TH1* anlohist = (TH1*)kffile.Get("anlo1/anlo1_nominal");
    TH1*  alohist = (TH1*)kffile.Get("alo/alo_nominal");
    TH1* aewkhist = (TH1*)kffile.Get("a_ewkcorr/a_ewkcorr");
    TH1* apdfhist = (TH1*)kffile.Get("anlo1/anlo1_pdfUp");

    TH1* nomhist  = (TH1*)kffile.Get("znlo1_over_anlo1/znlo1_over_anlo1");
    TH1* re1hist  = (TH1*)kffile.Get("znlo1_over_anlo1/znlo1_over_anlo1_renCorrUp");
    TH1* re2hist  = (TH1*)kffile.Get("znlo1_over_anlo1/znlo1_over_anlo1_renAcorrUp");
    TH1* fa1hist  = (TH1*)kffile.Get("znlo1_over_anlo1/znlo1_over_anlo1_facCorrUp");
    TH1* fa2hist  = (TH1*)kffile.Get("znlo1_over_anlo1/znlo1_over_anlo1_facAcorrUp");

    zpdfhist->Divide(znlohist);
    apdfhist->Divide(anlohist);

    re1hist->Divide(nomhist);
    re2hist->Divide(nomhist);
    fa1hist->Divide(nomhist);
    fa2hist->Divide(nomhist);

    znlohist->Divide(zlohist);
    anlohist->Divide(alohist);

    TFile fpfile("FP_v2.root");
    TH1* afpchist = (TH1*)fpfile.Get("FP_Down");

    vector<TH1*> zhists;
    vector<TH1*> ahists;
    vector<TH1*> ehists;
    if (kfact == 1) zhists.push_back(znlohist);
    if (kfact == 2) {zhists.push_back(znlohist); zhists.push_back(zewkhist);}
    if (kfact == 3) {zhists.push_back(znlohist); zhists.push_back(re1hist) ;}
    if (kfact == 4) {zhists.push_back(znlohist); zhists.push_back(fa1hist) ;}
    if (kfact == 5) {zhists.push_back(znlohist); zhists.push_back(re2hist) ;}
    if (kfact == 6) {zhists.push_back(znlohist); zhists.push_back(fa2hist) ;}
    if (kfact == 7) {zhists.push_back(znlohist); zhists.push_back(zpdfhist);}
    if (kfact == 8) zhists.push_back(znlohist);

    if (kfact == 1) ahists.push_back(anlohist);
    if (kfact == 2) {ahists.push_back(anlohist); ahists.push_back(aewkhist);}
    if (kfact == 3) ahists.push_back(anlohist);
    if (kfact == 4) ahists.push_back(anlohist);
    if (kfact == 5) ahists.push_back(anlohist);
    if (kfact == 6) ahists.push_back(anlohist);
    if (kfact == 7) {ahists.push_back(anlohist); ahists.push_back(apdfhist);}
    if (kfact == 8) {ahists.push_back(anlohist); zhists.push_back(afpchist);}

    makehist4(ntree, &nhist,  true, 0, 1.00, zhists, NULL);
    makehist4(dtree, &dhist,  true, 5, 1.00, ahists, NULL);

    string name = string("gamcor")+ext;    

    nhist.Divide(&dhist);
    TFile outfile((name+".root").c_str(), "RECREATE");
    nhist.SetName((name+"hist" ).c_str());
    nhist.Write();
    outfile.Close();

    nfile.Close();
    dfile.Close();
    kffile.Close();

    cout << "gamma+jets->Z(inv) transfer factor computed ..." << endl;
}

void sigdatamchist(TFile* outfile, bool blind=false) {

    TFile znfile("/Users/avartak/CMS/MonoX/FinalTrees/znn100toinf/sigtree.root");
    TFile wlfile("/Users/avartak/CMS/MonoX/FinalTrees/wln100toinf/sigtree.root");
    TFile zlfile("/Users/avartak/CMS/MonoX/FinalTrees/zll100toinf/sigtree.root");
    TFile ttfile("/Users/avartak/CMS/MonoX/FinalTrees/top/sigtree.root");
    TFile difile("/Users/avartak/CMS/MonoX/FinalTrees/qcd/sigtree.root");
    TFile qcfile("/Users/avartak/CMS/MonoX/FinalTrees/dibosons/sigtree.root");
    TFile sifile("/Users/avartak/CMS/MonoX/FinalTrees/psM200m1/sigtree.root");
    TFile dtfile("/Users/avartak/CMS/MonoX/FinalTrees/met/sigtree.root");

    TH1F znhist("zinvhist", "", nbins, bins);
    TH1F wlhist("wjethist", "", nbins, bins);
    TH1F zlhist("zjethist", "", nbins, bins);
    TH1F tthist("tbkghist", "", nbins, bins);
    TH1F dihist("dbkghist", "", nbins, bins);
    TH1F qchist("qbkghist", "", nbins, bins);
    TH1F sihist("sig1hist", "", nbins, bins);
    TH1F dthist("datahist", "", nbins, bins);

    TTree* zntree = (TTree*)znfile.Get("tree");
    TTree* wltree = (TTree*)wlfile.Get("tree");
    TTree* zltree = (TTree*)zlfile.Get("tree");
    TTree* tttree = (TTree*)ttfile.Get("tree");
    TTree* ditree = (TTree*)difile.Get("tree");
    TTree* qctree = (TTree*)qcfile.Get("tree");
    TTree* sitree = (TTree*)sifile.Get("tree");
    TTree* dttree = (TTree*)dtfile.Get("tree");

    TFile kffile("scalefactors_v4.root");
    TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
    TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
    TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");
    TH1* wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
    TH1*  wlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
    TH1* wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");

    znlohist->Divide(zlohist);
    wnlohist->Divide(wlohist);

    vector<TH1*> ehists;
    vector<TH1*> zhists;
    vector<TH1*> whists;
    zhists.push_back(znlohist); zhists.push_back(zewkhist);
    whists.push_back(wnlohist); whists.push_back(wewkhist);


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

    znfile.Close();
    wlfile.Close();
    zlfile.Close();
    ttfile.Close();
    difile.Close();
    qcfile.Close();
    sifile.Close();
    dtfile.Close();

    kffile.Close();

    cout << "Templates for the signal region computed ..." << endl;
}

void gamdatamchist(TFile* outfile) {

    TFile dtfile("/Users/avartak/CMS/MonoX/FinalTrees/singleph/gamtree.root");

    TH1F dthist("datahistgam", "", nbins, bins);
    TH1F qchist("qbkghistgam", "", nbins, bins);

    TTree* dttree = (TTree*)dtfile.Get("tree");

    vector<TH1*> ehists;

    makehist4(dttree, &dthist, false, 5, 1.00, ehists, NULL);
    makehist4(dttree, &qchist, false, 6, 1.00, ehists, NULL);

    //for (int i = 1; i <= dthist.GetNbinsX(); i++) qchist.SetBinContent(i, dthist.GetBinContent(i)*0.05);

    outfile->cd();
    dthist.Write();
    qchist.Write();

    dtfile.Close();

    cout << "Templates for the gamma+jets control region computed ..." << endl;
}

void lepdatamchist(TFile* outfile, int sample) {

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

    TFile ttfile((string("/Users/avartak/CMS/MonoX/FinalTrees/top/")          + filename).c_str());
    TFile difile((string("/Users/avartak/CMS/MonoX/FinalTrees/dibosons/")     + filename).c_str());
    TFile qcfile((string("/Users/avartak/CMS/MonoX/FinalTrees/qcd/")          + filename).c_str());
    TFile vlfile((string("/Users/avartak/CMS/MonoX/FinalTrees/") + vlfilename + filename).c_str());
    TFile dtfile((string("/Users/avartak/CMS/MonoX/FinalTrees/") + dtfilename + filename).c_str());

    TH1F dthist((string("datahist")+suffix).c_str(), "", nbins, bins);
    TH1F vlhist((string("lbkghist")+suffix).c_str(), "", nbins, bins);
    TH1F tthist((string("tbkghist")+suffix).c_str(), "", nbins, bins);
    TH1F dihist((string("dbkghist")+suffix).c_str(), "", nbins, bins);
    TH1F qchist((string("qbkghist")+suffix).c_str(), "", nbins, bins);

    TTree* dttree = (TTree*)dtfile.Get("tree");
    TTree* vltree = (TTree*)vlfile.Get("tree");
    TTree* tttree = (TTree*)ttfile.Get("tree");
    TTree* ditree = (TTree*)difile.Get("tree");
    TTree* qctree = (TTree*)qcfile.Get("tree");

    TFile kffile("scalefactors_v4.root");
    TH1* znlohist = (TH1*)kffile.Get("znlo012/znlo012_nominal");
    //TH1* znlohist = (TH1*)kffile.Get("znlo1/znlo1_nominal");
    TH1*  zlohist = (TH1*)kffile.Get("zlo/zlo_nominal");
    TH1* zewkhist = (TH1*)kffile.Get("z_ewkcorr/z_ewkcorr");
    TH1* wnlohist = (TH1*)kffile.Get("wnlo012/wnlo012_nominal");
    //TH1* wnlohist = (TH1*)kffile.Get("wnlo1/wnlo1_nominal");
    TH1*  wlohist = (TH1*)kffile.Get("wlo/wlo_nominal");
    TH1* wewkhist = (TH1*)kffile.Get("w_ewkcorr/w_ewkcorr");

    znlohist->Divide(zlohist);
    wnlohist->Divide(wlohist);

    vector<TH1*> ehists;
    vector<TH1*> vhists;
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

    dtfile.Close();
    vlfile.Close();
    ttfile.Close();
    difile.Close();
    qcfile.Close();

    cout << "Templates for the lepton control region computed ..." << endl;
}

void hists() {
    /*
    makezmmcorhist(); 
    makezeecorhist();
    makewmncorhist();    
    makewencorhist();    
    makegamcorhist();
    makezwjcorhist();

    makegamcorhist("qcd", 1);
    makegamcorhist("ewk", 2);
    makegamcorhist("re1", 3);
    makegamcorhist("fa1", 4);
    makegamcorhist("re2", 5);
    makegamcorhist("fa2", 6);
    makegamcorhist("pdf", 7);
    makegamcorhist("fpc", 8);

    makezwjcorhist("qcd", 1);
    makezwjcorhist("ewk", 2);
    makezwjcorhist("re1", 3);
    makezwjcorhist("fa1", 4);
    makezwjcorhist("re2", 5);
    makezwjcorhist("fa2", 6);
    makezwjcorhist("pdf", 7);
    */

    TFile* zmmcorfile = new TFile("zmmcor.root");
    TFile* zeecorfile = new TFile("zeecor.root");
    TFile* wmncorfile = new TFile("wmncor.root");
    TFile* wencorfile = new TFile("wencor.root");
    TFile* zwjcorfile = new TFile("zwjcor.root");
    TFile* gamcorfile = new TFile("gamcor.root");

    TFile* gamcorqcdfile = new TFile("gamcorqcd.root");
    TFile* gamcorewkfile = new TFile("gamcorewk.root");
    TFile* gamcorre1file = new TFile("gamcorre1.root");
    TFile* gamcorfa1file = new TFile("gamcorfa1.root");
    TFile* gamcorre2file = new TFile("gamcorre2.root");
    TFile* gamcorfa2file = new TFile("gamcorfa2.root");
    TFile* gamcorpdffile = new TFile("gamcorpdf.root");
    TFile* gamcorfpcfile = new TFile("gamcorfpc.root");

    TFile* zwjcorqcdfile = new TFile("zwjcorqcd.root");
    TFile* zwjcorewkfile = new TFile("zwjcorewk.root");
    TFile* zwjcorre1file = new TFile("zwjcorre1.root");
    TFile* zwjcorfa1file = new TFile("zwjcorfa1.root");
    TFile* zwjcorre2file = new TFile("zwjcorre2.root");
    TFile* zwjcorfa2file = new TFile("zwjcorfa2.root");
    TFile* zwjcorpdffile = new TFile("zwjcorpdf.root");

    TH1* zmmcorhist = (TH1*)zmmcorfile->Get("zmmcorhist");    
    TH1* zeecorhist = (TH1*)zeecorfile->Get("zeecorhist");    
    TH1* wmncorhist = (TH1*)wmncorfile->Get("wmncorhist");    
    TH1* wencorhist = (TH1*)wencorfile->Get("wencorhist");    
    TH1* zwjcorhist = (TH1*)zwjcorfile->Get("zwjcorhist");    
    TH1* gamcorhist = (TH1*)gamcorfile->Get("gamcorhist");    

    TH1* gamcorewkhist = (TH1*)gamcorewkfile->Get("gamcorewkhist");    
    TH1* gamcorqcdhist = (TH1*)gamcorqcdfile->Get("gamcorqcdhist");    
    TH1* gamcorre1hist = (TH1*)gamcorre1file->Get("gamcorre1hist");    
    TH1* gamcorfa1hist = (TH1*)gamcorfa1file->Get("gamcorfa1hist");    
    TH1* gamcorre2hist = (TH1*)gamcorre2file->Get("gamcorre2hist");    
    TH1* gamcorfa2hist = (TH1*)gamcorfa2file->Get("gamcorfa2hist");    
    TH1* gamcorpdfhist = (TH1*)gamcorpdffile->Get("gamcorpdfhist");    
    TH1* gamcorfpchist = (TH1*)gamcorfpcfile->Get("gamcorfpchist");    

    TH1* gamuncewkhist = (TH1*)gamcorewkhist->Clone("gamuncewkhist");    
    gamuncewkhist->Divide(gamcorqcdhist);
    for (int i = 1; i <= gamuncewkhist->GetNbinsX(); i++) gamuncewkhist->SetBinContent(i, fabs(gamuncewkhist->GetBinContent(i)-1.0));
    gamuncewkhist->SetName("ZG_EWK");

    TH1* gamuncre1hist = (TH1*)gamcorre1hist->Clone("gamuncre1hist");    
    gamuncre1hist->Divide(gamcorqcdhist);
    for (int i = 1; i <= gamuncre1hist->GetNbinsX(); i++) gamuncre1hist->SetBinContent(i, fabs(gamuncre1hist->GetBinContent(i)-1.0));
    gamuncre1hist->SetName("ZG_RenScale1");

    TH1* gamuncfa1hist = (TH1*)gamcorfa1hist->Clone("gamuncfa1hist");    
    gamuncfa1hist->Divide(gamcorqcdhist);
    for (int i = 1; i <= gamuncfa1hist->GetNbinsX(); i++) gamuncfa1hist->SetBinContent(i, fabs(gamuncfa1hist->GetBinContent(i)-1.0));
    gamuncfa1hist->SetName("ZG_FactScale1");

    TH1* gamuncre2hist = (TH1*)gamcorre2hist->Clone("gamuncre2hist");    
    gamuncre2hist->Divide(gamcorqcdhist);
    for (int i = 1; i <= gamuncre2hist->GetNbinsX(); i++) gamuncre2hist->SetBinContent(i, fabs(gamuncre2hist->GetBinContent(i)-1.0));
    gamuncre2hist->SetName("ZG_RenScale2");

    TH1* gamuncfa2hist = (TH1*)gamcorfa2hist->Clone("gamuncfa2hist");    
    gamuncfa2hist->Divide(gamcorqcdhist);
    for (int i = 1; i <= gamuncfa2hist->GetNbinsX(); i++) gamuncfa2hist->SetBinContent(i, fabs(gamuncfa2hist->GetBinContent(i)-1.0));
    gamuncfa2hist->SetName("ZG_FactScale2");

    TH1* gamuncpdfhist = (TH1*)gamcorpdfhist->Clone("gamuncpdfhist");    
    gamuncpdfhist->Divide(gamcorqcdhist);
    for (int i = 1; i <= gamuncpdfhist->GetNbinsX(); i++) gamuncpdfhist->SetBinContent(i, fabs(gamuncpdfhist->GetBinContent(i)-1.0));
    gamuncpdfhist->SetName("ZG_PDF");

    TH1* gamuncfpchist = (TH1*)gamcorfpchist->Clone("gamuncfpchist");    
    gamuncfpchist->Divide(gamcorqcdhist);
    for (int i = 1; i <= gamuncfpchist->GetNbinsX(); i++) gamuncfpchist->SetBinContent(i, fabs(gamuncfpchist->GetBinContent(i)-1.0));
    gamuncfpchist->SetName("ZG_Footprint");

    TH1* zwjcorewkhist = (TH1*)zwjcorewkfile->Get("zwjcorewkhist");    
    TH1* zwjcorqcdhist = (TH1*)zwjcorqcdfile->Get("zwjcorqcdhist");    
    TH1* zwjcorre1hist = (TH1*)zwjcorre1file->Get("zwjcorre1hist");    
    TH1* zwjcorfa1hist = (TH1*)zwjcorfa1file->Get("zwjcorfa1hist");    
    TH1* zwjcorre2hist = (TH1*)zwjcorre2file->Get("zwjcorre2hist");    
    TH1* zwjcorfa2hist = (TH1*)zwjcorfa2file->Get("zwjcorfa2hist");    
    TH1* zwjcorpdfhist = (TH1*)zwjcorpdffile->Get("zwjcorpdfhist");    

    TH1* zwjuncewkhist = (TH1*)zwjcorewkhist->Clone("zwjuncewkhist");    
    zwjuncewkhist->Divide(zwjcorqcdhist);
    for (int i = 1; i <= zwjuncewkhist->GetNbinsX(); i++) zwjuncewkhist->SetBinContent(i, fabs(zwjuncewkhist->GetBinContent(i)-1.0));
    zwjuncewkhist->SetName("ZW_EWK");

    TH1* zwjuncre1hist = (TH1*)zwjcorre1hist->Clone("zwjuncre1hist");
    zwjuncre1hist->Divide(zwjcorqcdhist);
    for (int i = 1; i <= zwjuncre1hist->GetNbinsX(); i++) zwjuncre1hist->SetBinContent(i, fabs(zwjuncre1hist->GetBinContent(i)-1.0));
    zwjuncre1hist->SetName("ZW_RenScale1");

    TH1* zwjuncfa1hist = (TH1*)zwjcorfa1hist->Clone("zwjuncfa1hist");
    zwjuncfa1hist->Divide(zwjcorqcdhist);
    for (int i = 1; i <= zwjuncfa1hist->GetNbinsX(); i++) zwjuncfa1hist->SetBinContent(i, fabs(zwjuncfa1hist->GetBinContent(i)-1.0));
    zwjuncfa1hist->SetName("ZW_FactScale1");

    TH1* zwjuncre2hist = (TH1*)zwjcorre2hist->Clone("zwjuncre2hist");
    zwjuncre2hist->Divide(zwjcorqcdhist);
    for (int i = 1; i <= zwjuncre2hist->GetNbinsX(); i++) zwjuncre2hist->SetBinContent(i, fabs(zwjuncre2hist->GetBinContent(i)-1.0));
    zwjuncre2hist->SetName("ZW_RenScale2");

    TH1* zwjuncfa2hist = (TH1*)zwjcorfa2hist->Clone("zwjuncfa2hist");
    zwjuncfa2hist->Divide(zwjcorqcdhist);
    for (int i = 1; i <= zwjuncfa2hist->GetNbinsX(); i++) zwjuncfa2hist->SetBinContent(i, fabs(zwjuncfa2hist->GetBinContent(i)-1.0));
    zwjuncfa2hist->SetName("ZW_FactScale2");

    TH1* zwjuncpdfhist = (TH1*)zwjcorpdfhist->Clone("zwjuncpdfhist");
    zwjuncpdfhist->Divide(zwjcorqcdhist);
    for (int i = 1; i <= zwjuncpdfhist->GetNbinsX(); i++) zwjuncpdfhist->SetBinContent(i, fabs(zwjuncpdfhist->GetBinContent(i)-1.0));
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
