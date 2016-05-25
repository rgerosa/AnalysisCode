#ifndef DATASETS_H
#define DATASETS_H

void fillDataSetNames(vector<string>& files, string sample) {
    if (sample == "znn") {
        files.push_back("ZJets");
        files.push_back("ZJetsToNuNu_HT-100To200_13TeV-madgraph");
        files.push_back("ZJetsToNuNu_HT-200To400_13TeV-madgraph");
        files.push_back("ZJetsToNuNu_HT-400To600_13TeV-madgraph");
        files.push_back("ZJetsToNuNu_HT-600To800_13TeV-madgraph");
        files.push_back("ZJetsToNuNu_HT-800To1200_13TeV-madgraph");
        files.push_back("ZJetsToNuNu_HT-1200To2500_13TeV-madgraph");
        files.push_back("ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph");
    }

    if (sample == "wln") {
        files.push_back("WJets");
        files.push_back("WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("WJetsToLNu_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("WJetsToLNu_HT-600To800_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("WJetsToLNu_HT-800To1200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("WJetsToLNu_HT-1200To2500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("WJetsToLNu_HT-2500ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
    }

    if (sample == "zll") {
        files.push_back("DYJets");
        files.push_back("DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
    }

    if (sample == "top") {
        files.push_back("Top");
        files.push_back("TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8");
        files.push_back("ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1");
        files.push_back("ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1");
        files.push_back("ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1");
        files.push_back("ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1");
        files.push_back("ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1");
    }

    if (sample == "dib") {
        files.push_back("DiBoson");
        files.push_back("WWTo2L2Nu_13TeV-powheg");
        files.push_back("WWTo4Q_13TeV-powheg");
        files.push_back("WWToLNuQQ_13TeV-powheg");
        files.push_back("WZ_TuneCUETP8M1_13TeV-pythia8");
        files.push_back("ZZ_TuneCUETP8M1_13TeV-pythia8");
    }

    if (sample == "qcd") {
        files.push_back("QCD");
        files.push_back("QCD_HT100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("QCD_HT200to300_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("QCD_HT300to500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("QCD_HT500to700_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("QCD_HT700to1000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("QCD_HT1000to1500_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("QCD_HT1500to2000_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("QCD_HT2000toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
    }

    if (sample == "gam") {
        files.push_back("PhotonJets");
        files.push_back("GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
        files.push_back("GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8");
    }

}


#endif 
