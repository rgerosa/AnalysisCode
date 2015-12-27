#ifndef MAKEHIST4_H
#define MAKEHIST4_H

#include <vector>

void makehist4(TTree* tree, TH1* hist, bool isMC, int sample, double scale, vector<TH1*> khists, TH1* rhist=NULL) {
    double lumi = 2.11;

    TFile* pufile = new TFile("purwt.root");
    TH1*   puhist = (TH1*)pufile->Get("puhist");

    TFile* sffile = new TFile("leptonIDsfs.root");
    TH2*  msflhist = (TH2*)sffile->Get("muon_loose_SF");
    TH2*  msfthist = (TH2*)sffile->Get("muon_tight_SF");

    TH2*  esflhist = (TH2*)sffile->Get("electron_veto_SF");
    TH2*  esfthist = (TH2*)sffile->Get("electron_tight_SF");

    TFile* psffile = new TFile("PhotonSFandEffandPurity_Lumi2p1fb_2211.root");
    TH2*  psfhist = (TH2*)psffile->Get("PhotonSF");

    TFile* purfile = new TFile("PhotonSFandEffandPurity_Lumi2p1fb_2211.root");
    TH2*  purhist = (TH2*)psffile->Get("PhotonPurity");

    TFile* trefile = new TFile("leptonTrigsfs.root");
    TH2*  trehist = (TH2*)trefile->Get("hltel27_SF");

    TFile* trmfile = new TFile("mettrigSF.root");
    TH1*  trmhist = (TH1*)trefile->Get("mettrigSF");

    hist->Sumw2();

    TBranch  *brun        = tree->GetBranch("run");
    TBranch  *bnvtx       = tree->GetBranch("nvtx");
    TBranch  *bxsec       = tree->GetBranch("xsec");
    TBranch  *bwgt        = tree->GetBranch("wgt");
    TBranch  *bwgtsum     = tree->GetBranch("wgtsum"); 

    TBranch  *bhltm       = tree->GetBranch("hltmet90");
    TBranch  *bhlte       = tree->GetBranch("hltsingleel");
    TBranch  *bhltp       = tree->GetBranch("hltphoton165");
    TBranch  *bhltp2      = tree->GetBranch("hltphoton175");

    TBranch  *bfhbhe      = tree->GetBranch("flaghbheloose");
    TBranch  *bfhbiso     = tree->GetBranch("flaghbheiso");
    TBranch  *bfcsc;
    TBranch  *bfeeb;
    if (isMC) {
        bfcsc             = tree->GetBranch("flagcsctight");
        bfeeb             = tree->GetBranch("flageebadsc");
    }
    else {
        bfcsc             = tree->GetBranch("flagcscnew");
        bfeeb             = tree->GetBranch("flageescnew");
    }

    TBranch  *bnjets      = tree->GetBranch("njets");
    TBranch  *bnbjets     = tree->GetBranch("nbjetslowpt");
    TBranch  *bjetpt      = tree->GetBranch("signaljetpt");
    TBranch  *bj1pt       = tree->GetBranch("leadingjetpt");
    TBranch  *bchfrac     = tree->GetBranch("signaljetCHfrac");
    TBranch  *bnhfrac     = tree->GetBranch("signaljetNHfrac");
    TBranch  *bemfrac     = tree->GetBranch("signaljetEMfrac");

    TBranch  *bmet        = tree->GetBranch("t1pfmet");
    TBranch  *bmetphi     = tree->GetBranch("t1pfmetphi");
    TBranch  *bmmet       = tree->GetBranch("t1mumet");
    TBranch  *bmmetphi    = tree->GetBranch("t1mumetphi");
    TBranch  *bemet       = tree->GetBranch("t1elmet");
    TBranch  *bemetphi    = tree->GetBranch("t1elmetphi");
    TBranch  *bpmet       = tree->GetBranch("t1phmet");
    TBranch  *bpmetphi    = tree->GetBranch("t1phmetphi");
    TBranch  *bjmmdphi    = tree->GetBranch("incjetmetdphimin4");
    TBranch  *bjemdphi    = tree->GetBranch("incjetelmetdphimin4");
    TBranch  *bjpmdphi    = tree->GetBranch("incjetphmetdphimin4");

    TBranch  *bmu1pid     = tree->GetBranch("mu1pid");
    TBranch  *bmu2pid     = tree->GetBranch("mu2pid");
    TBranch  *bmu1id      = tree->GetBranch("mu1id");
    TBranch  *bmu2id      = tree->GetBranch("mu2id");
    TBranch  *bmu1pt      = tree->GetBranch("mu1pt");
    TBranch  *bmu2pt      = tree->GetBranch("mu2pt");
    TBranch  *bmu1eta     = tree->GetBranch("mu1eta");
    TBranch  *bmu2eta     = tree->GetBranch("mu2eta");
    TBranch  *bel1pid     = tree->GetBranch("el1pid");
    TBranch  *bel2pid     = tree->GetBranch("el2pid");
    TBranch  *bel1id      = tree->GetBranch("el1id");
    TBranch  *bel2id      = tree->GetBranch("el2id");
    TBranch  *bel1pt      = tree->GetBranch("el1pt");
    TBranch  *bel2pt      = tree->GetBranch("el2pt");
    TBranch  *bel1eta     = tree->GetBranch("el1eta");
    TBranch  *bel2eta     = tree->GetBranch("el2eta");
    TBranch  *bphid       = tree->GetBranch("phidm");
    TBranch  *bphpt       = tree->GetBranch("phpt");
    TBranch  *bpheta      = tree->GetBranch("pheta");
    TBranch  *bphphi      = tree->GetBranch("phphi");

    TBranch  *bwzpt       = tree->GetBranch("wzpt");
    TBranch  *bwzeta      = tree->GetBranch("wzeta");
    TBranch  *bzmass      = tree->GetBranch("zmass");
    TBranch  *bzmmpt      = tree->GetBranch("zpt");
    TBranch  *bzeept      = tree->GetBranch("zeept");
    TBranch  *bzmmeta     = tree->GetBranch("zeta");
    TBranch  *bzeeeta     = tree->GetBranch("zeeeta");

    UInt_t   run          = 0;
    UInt_t   nvtx         = 0;
    Double_t xsec         = 0.0;
    Double_t wgt          = 0.0;
    Double_t wgtsum       = 0.0;

    UChar_t  hltm         = 0;
    UChar_t  hlte         = 0;
    UChar_t  hltp         = 0;
    UChar_t  hltp2        = 0;

    UChar_t  fhbhe        = 0;
    UChar_t  fhbiso       = 0;
    UChar_t  fcsc         = 0;
    UChar_t  feeb         = 0;

    UInt_t   njets        = 0;
    UInt_t   nbjets       = 0;
    Double_t jetpt        = 0.0;
    Double_t j1pt         = 0.0;
    Double_t chfrac       = 0.0;
    Double_t nhfrac       = 0.0;
    Double_t emfrac       = 0.0;

    Double_t pfmet        = 0.0;
    Double_t pfmetphi     = 0.0;
    Double_t mmet         = 0.0;
    Double_t mmetphi      = 0.0;
    Double_t emet         = 0.0;
    Double_t emetphi      = 0.0;
    Double_t pmet         = 0.0;
    Double_t pmetphi      = 0.0;
    Double_t jmmdphi      = 0.0;
    Double_t jemdphi      = 0.0;
    Double_t jpmdphi      = 0.0;

    Int_t    mu1pid       = 0;
    Int_t    mu2pid       = 0;
    Int_t    mu1id        = 0;
    Int_t    mu2id        = 0;
    Double_t mu1pt        = 0.0;
    Double_t mu2pt        = 0.0;
    Double_t mu1eta       = 0.0;
    Double_t mu2eta       = 0.0;
    Int_t    el1pid       = 0;
    Int_t    el2pid       = 0;
    Int_t    el1id        = 0;
    Int_t    el2id        = 0;
    Double_t el1pt        = 0.0;
    Double_t el2pt        = 0.0;
    Double_t el1eta       = 0.0;
    Double_t el2eta       = 0.0;
    Int_t    phid         = 0;
    Double_t phpt         = 0.0;
    Double_t pheta        = 0.0;
    Double_t phphi        = 0.0;

    Double_t wzpt         = 0.0;
    Double_t wzeta        = 0.0;
    Double_t zmass        = 0.0;
    Double_t zmmpt        = 0.0;
    Double_t zeept        = 0.0;
    Double_t zmmeta       = 0.0;
    Double_t zeeeta       = 0.0;

    brun                  ->SetAddress(&run);
    bnvtx                 ->SetAddress(&nvtx);
    bxsec                 ->SetAddress(&xsec);
    bwgt                  ->SetAddress(&wgt);
    bwgtsum               ->SetAddress(&wgtsum);

    bhltm                 ->SetAddress(&hltm);
    bhlte                 ->SetAddress(&hlte);
    bhltp                 ->SetAddress(&hltp);
    bhltp2                ->SetAddress(&hltp2);

    bfhbhe                ->SetAddress(&fhbhe); 
    bfhbiso               ->SetAddress(&fhbiso);
    bfcsc                 ->SetAddress(&fcsc);
    bfeeb                 ->SetAddress(&feeb);

    bnjets                ->SetAddress(&njets);
    bnbjets               ->SetAddress(&nbjets);
    bjetpt                ->SetAddress(&jetpt);
    bj1pt                 ->SetAddress(&j1pt);
    bchfrac               ->SetAddress(&chfrac);
    bnhfrac               ->SetAddress(&nhfrac);
    bemfrac               ->SetAddress(&emfrac);

    bmet                  ->SetAddress(&pfmet);
    bmetphi               ->SetAddress(&pfmetphi);
    bmmet                 ->SetAddress(&mmet);
    bmmetphi              ->SetAddress(&mmetphi);
    bemet                 ->SetAddress(&emet);
    bemetphi              ->SetAddress(&emetphi);
    bpmet                 ->SetAddress(&pmet);
    bpmetphi              ->SetAddress(&pmetphi);
    bjmmdphi              ->SetAddress(&jmmdphi);
    bjemdphi              ->SetAddress(&jemdphi);
    bjpmdphi              ->SetAddress(&jpmdphi);

    bmu1pid               ->SetAddress(&mu1pid);
    bmu2pid               ->SetAddress(&mu2pid);
    bmu1id                ->SetAddress(&mu1id);
    bmu2id                ->SetAddress(&mu2id);
    bmu1pt                ->SetAddress(&mu1pt);
    bmu2pt                ->SetAddress(&mu2pt);
    bmu1eta               ->SetAddress(&mu1eta);
    bmu2eta               ->SetAddress(&mu2eta);
    bel1pid               ->SetAddress(&el1pid);
    bel2pid               ->SetAddress(&el2pid);
    bel1id                ->SetAddress(&el1id);
    bel2id                ->SetAddress(&el2id);
    bel1pt                ->SetAddress(&el1pt);
    bel2pt                ->SetAddress(&el2pt);
    bel1eta               ->SetAddress(&el1eta);
    bel2eta               ->SetAddress(&el2eta);
    bphid                 ->SetAddress(&phid);
    bphpt                 ->SetAddress(&phpt);
    bpheta                ->SetAddress(&pheta);
    bphphi                ->SetAddress(&phphi);

    bwzpt                 ->SetAddress(&wzpt);
    bwzeta                ->SetAddress(&wzeta);
    bzmass                ->SetAddress(&zmass);
    bzmmpt                ->SetAddress(&zmmpt);
    bzeept                ->SetAddress(&zeept);
    bzmmeta               ->SetAddress(&zmmeta);
    bzeeeta               ->SetAddress(&zeeeta);

    for (Long64_t i = 0; i < tree->GetEntries(); i++) {
        brun              ->GetEvent(i);
        bnvtx             ->GetEvent(i);
        bxsec             ->GetEvent(i);
        bwgt              ->GetEvent(i);
        bwgtsum           ->GetEvent(i);

        bhltm             ->GetEvent(i);
        bhlte             ->GetEvent(i);
        bhltp             ->GetEvent(i);
        bhltp2            ->GetEvent(i);

        bfhbhe            ->GetEvent(i);
        bfhbiso           ->GetEvent(i);
        bfcsc             ->GetEvent(i);
        bfeeb             ->GetEvent(i);

        bnjets            ->GetEvent(i);
        bnbjets           ->GetEvent(i);
        bjetpt            ->GetEvent(i);
        bj1pt             ->GetEvent(i);
        bchfrac           ->GetEvent(i);
        bnhfrac           ->GetEvent(i);
        bemfrac           ->GetEvent(i);

        bmet              ->GetEvent(i);
        bmetphi           ->GetEvent(i);
        bmmet             ->GetEvent(i);
        bmmetphi          ->GetEvent(i);
        bemet             ->GetEvent(i);
        bemetphi          ->GetEvent(i);
        bpmet             ->GetEvent(i);
        bpmetphi          ->GetEvent(i);
        bjmmdphi          ->GetEvent(i);
        bjemdphi          ->GetEvent(i);
        bjpmdphi          ->GetEvent(i);

        bmu1pid           ->GetEvent(i);
        bmu2pid           ->GetEvent(i);
        bmu1id            ->GetEvent(i);
        bmu2id            ->GetEvent(i);
        bmu1pt            ->GetEvent(i);
        bmu2pt            ->GetEvent(i);
        bmu1eta           ->GetEvent(i);
        bmu2eta           ->GetEvent(i);

        bel1pid           ->GetEvent(i);
        bel2pid           ->GetEvent(i);
        bel1id            ->GetEvent(i);
        bel2id            ->GetEvent(i);
        bel1pt            ->GetEvent(i);
        bel2pt            ->GetEvent(i);
        bel1eta           ->GetEvent(i);
        bel2eta           ->GetEvent(i);

        bphid             ->GetEvent(i);
        bphpt             ->GetEvent(i);
        bpheta            ->GetEvent(i);
        bphphi            ->GetEvent(i);

        bwzpt             ->GetEvent(i);
        bwzeta            ->GetEvent(i);
        bzmass            ->GetEvent(i);
        bzmmpt            ->GetEvent(i);
        bzeept            ->GetEvent(i);
        bzmmeta           ->GetEvent(i);
        bzeeeta           ->GetEvent(i);
        
        Double_t hlt = 0.0;
        if      (sample == 0 || sample == 1 || sample == 2) hlt = hltm;
        else if (sample == 3 || sample == 4)                hlt = hlte;
        else if (sample == 5 || sample == 6)                hlt = hltp;

        if ((sample == 5 || sample == 6) && hltp2 > 0)      hlt = hltp2;

        Double_t jmdphi = 0.0;
        if      (sample == 0 || sample == 1 || sample == 2) jmdphi = fabs(jmmdphi);
        else if (sample == 3 || sample == 4)                jmdphi = fabs(jemdphi);
        else if (sample == 5 || sample == 6)                jmdphi = fabs(jpmdphi);

        Double_t met = 0.0;
        if      (sample == 0 || sample == 1 || sample == 2) met = mmet;
        else if (sample == 3 || sample == 4)                met = emet;
        else if (sample == 5 || sample == 6)                met = pmet;

        Double_t zpt = 0.0;
        if      (sample == 1) zpt = zmmpt;
        else if (sample == 3) zpt = zeept;
        else if (sample == 5) zpt = phpt;

        Double_t puwgt = 0.;
        if (nvtx <= 35) puwgt = puhist->GetBinContent(nvtx);
        //puwgt = 1.0;

        Int_t    id1   = 0;
        Int_t    id2   = 0;
        Double_t pt1   = 0.0;
        Double_t pt2   = 0.0;
        Double_t eta1  = 0.0;
        Double_t eta2  = 0.0;

        if (sample == 1 || sample == 2) {
            id1  = mu1id;
            id2  = mu2id;
            pt1  = mu1pt;
            pt2  = mu2pt;
            eta1 = fabs(mu1eta);
            eta2 = fabs(mu2eta);
        }
        else if (sample == 3 || sample == 4) {
            id1  = el1id;
            id2  = el2id;
            pt1  = el1pt;
            pt2  = el2pt;
            eta1 = fabs(el1eta);
            eta2 = fabs(el2eta);
        }
        else if (sample == 5 || sample == 6) {
            id1  = 1.0;
            id2  = 1.0;
            pt1  = phpt;
            eta1 = fabs(pheta);
        }

        if (pt1 >= 1000.) pt1 = 999.0;
        if (pt2 >= 1000.) pt2 = 999.0;

        TH2* sflhist = NULL;
        TH2* sfthist = NULL;

        if (sample == 1 || sample == 2) {
            sflhist = msflhist;
            sfthist = msfthist;
        }
        if (sample == 3 || sample == 4) {
            sflhist = esflhist;
            sfthist = esfthist;
        }

        Double_t sfwgt = 1.0;
        if (isMC && sflhist && sfthist) {
            if (pt1 > 0.) {
                if (id1 == 1) sfwgt *= sfthist->GetBinContent(sfthist->FindBin(pt1, eta1)); 
                else          sfwgt *= sflhist->GetBinContent(sflhist->FindBin(pt1, eta1)); 
            }
            if (pt2 > 0.) {
                if (id2 == 1) sfwgt *= sfthist->GetBinContent(sfthist->FindBin(pt2, eta2)); 
                else          sfwgt *= sflhist->GetBinContent(sflhist->FindBin(pt2, eta2)); 
            }
        }
        if (isMC && trehist && sample == 4) {
            if (pt1 > 0. && id1 == 1) {
                sfwgt *= trehist->GetBinContent(trehist->FindBin(pt1, eta1));
            }
        }
        if (isMC && psfhist && (sample == 5 || sample == 6)) {
            if (pt1 > 0. && id1 == 1) {
                sfwgt *= psfhist->GetBinContent(psfhist->FindBin(pt1, eta1));
            }
        }
        if (!isMC && purhist && sample == 6) {
            if (pt1 > 175. && id1 == 1) {
                sfwgt *= (1.0 - purhist->GetBinContent(purhist->FindBin(pt1, eta1)));
            }
        }
        if (isMC && trmhist && (sample == 0 || sample == 1 || sample == 2)) {
            sfwgt *= trmhist->GetBinContent(trmhist->FindBin(met));
        }

        Double_t kvar = wzpt;
        if (wzpt < 100. ) wzpt = 100.;
        if (wzpt > 1000.) wzpt = 999.;

        Double_t rwgt = 1.0;
        if (rhist) rwgt = rhist->GetBinContent(rhist->FindBin(wzpt));

        Double_t kwgt = 1.0;
        for (unsigned i = 0; i < khists.size(); i++) {
            if (isMC && khists[i]) kwgt *= khists[i]->GetBinContent(khists[i]->FindBin(wzpt));
        }

        if (hlt  == 0) continue;
        if (fhbhe == 0 || fhbiso == 0 || fcsc == 0 || feeb == 0) continue;
        if (chfrac < 0.1) continue;
        if (nhfrac > 0.8) continue;
        if (njets  < 1) continue;
        if (nbjets > 0) continue;
        if (jetpt  < 100.) continue;
        if (jetpt  < j1pt) continue;
        if (jmdphi < 0.5) continue;
        if (sample == 1 && mu1pid == mu2pid) continue;
        if (sample == 3 && el1pid == el2pid) continue;
        if ((sample == 5 || sample == 6) && phpt < 175.) continue;
        if ((sample == 5 || sample == 6) && fabs(pheta) > 1.4442) continue;
        if (sample == 4 && pfmet < 50.) continue;
        if (met < 200.) continue;

        double fillvar = met;
        if (fillvar >= hist->GetBinLowEdge(hist->GetNbinsX())+hist->GetBinWidth(hist->GetNbinsX())) fillvar = hist->GetXaxis()->GetBinCenter(hist->GetNbinsX());

        double evtwgt = 1.0;
        if (isMC) evtwgt = xsec*scale*lumi*wgt*puwgt*sfwgt*rwgt*kwgt/wgtsum;
        if (!isMC && sample == 6) evtwgt = sfwgt;

        hist->Fill(fillvar, evtwgt);
    }

    pufile  ->Close();
    sffile  ->Close();
    psffile ->Close();
    trefile ->Close();

}

#endif
