{
    gROOT->ProcessLine(".L PDFs/RooErfExpPdf.cc+");
    gROOT->ProcessLine(".L makeTnPTemplates.C");
    gROOT->ProcessLine("makeTnPTemplates\(\"/home/rgerosa/MONOJET_ANALYSIS_2016_Data/TagAndProbe/TagAndProbe_MC/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8\",\"muon\",\"TemplatesNominalReco\",true)");

}

