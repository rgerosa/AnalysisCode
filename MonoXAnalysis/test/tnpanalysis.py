import FWCore.ParameterSet.Config as cms

## CMSSW command line parameter parser                                                                                                                                          
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('python')

## data or MC options                                                                                                                                                          
options.register (
    'isMC',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'flag to indicate data or MC');

options.register(
    'typeID',"tight",VarParsing.multiplicity.singleton,VarParsing.varType.string,
    'string to indicate the lepton type to be considered');

options.register(
    'leptonPID',11,VarParsing.multiplicity.singleton,VarParsing.varType.int,
    '11 means electrons, 13 means muons');

options.register(
    'floatShapeParameters',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,
    'float shape parameters in the fit');

options.register(
    'binsForMassPlots',25,VarParsing.multiplicity.singleton,VarParsing.varType.int,
    'number of bins for the output plot');

options.register(
    'backgroundType',0,VarParsing.multiplicity.singleton,VarParsing.varType.int,
    'Available options at the moment: 0 means exponential, 1 means chebychev polynomial');

## parsing command line arguments                                                                                                                                            
options.parseArguments()

process = cms.Process("TNP")
process.source = cms.Source("EmptySource")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = ['cout', 'cerr']
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

#### name for the output file
if options.isMC : 
    OutputFilePrefix = "efficiency-mc"
else :
    OutputFilePrefix = "efficiency-data"

PDFName = "pdfSignalPlusBackground"

### Binnning in eta and pt for the efficiency fit
EfficiencyBins = cms.PSet(
    pt     = cms.vdouble( 10., 20., 30., 40., 50., 70., 100.),
    abseta = cms.vdouble( 1.2, 2.4 )
    )
### Defining the observable as well as the PDF names for each bin (total PDF for fitting pass and fail samples)-> One per pt
EfficiencyBinningSpecification = cms.PSet(
    UnbinnedVariables = cms.vstring("mass"),
    BinnedVariables = cms.PSet(
        EfficiencyBins,
    ),
    BinToPDFmap = cms.vstring([
        "pdfSignalPlusBackground", 
        "*pt_bin0*", "pdfSignalPlusBackgroundb0",
        "*pt_bin1*", "pdfSignalPlusBackgroundb1",
        "*pt_bin2*", "pdfSignalPlusBackgroundb2",
        "*pt_bin3*", "pdfSignalPlusBackgroundb3",
        "*pt_bin4*", "pdfSignalPlusBackgroundb4",
        "*pt_bin5*", "pdfSignalPlusBackgroundb5",
    ])  
)

IDModule = cms.PSet(
    Id = cms.PSet(
        EfficiencyBinningSpecification,
        EfficiencyCategoryAndState = cms.vstring(options.typeID,"pass"),
    )
)

leptonTree = "";
if options.leptonPID == 11:
    leptonTree = "electrontnptree";
elif options.leptonPID == 13:    
    leptonTree = "muontnptree";
else:
    sys.exit('leptind pdgId not recognized --> not possible to run tag and probe --> exit');

process.muonIdTnP = cms.EDAnalyzer("TagProbeFitTreeAnalyzer",
                                   InputFileNames = cms.vstring(options.inputFiles),
                                   InputDirectoryName = cms.string(leptonTree),
                                   InputTreeName = cms.string("fitter_tree"),
                                   OutputFileName = cms.string(OutputFilePrefix+".root"),
                                   NumCPU = cms.uint32(1),
                                   SaveWorkspace = cms.bool(True),
                                   floatShapeParameters = cms.bool(options.floatShapeParameters),
                                   binsForMassPlots = cms.uint32(options.binsForMassPlots),
                                   Variables = cms.PSet(
        mass   = cms.vstring("Tag-Probe Mass", "65.0", "115.0", "GeV/c^{2}"),
        pt     = cms.vstring("Probe p_{T}", "0", "1000", "GeV/c"),
        abseta = cms.vstring("Probe #eta",  "0", "2.4", ""),                
        ),
                                   
                                   Categories = cms.PSet(
        mcTrue  = cms.vstring("MC Truth", "dummy[true=1,false=0]"),
        tightid = cms.vstring("Tight ID", "dummy[pass=1,fail=0]"), 
        ),
                                   PDFs = cms.PSet(
        pdfSignalPlusBackground = cms.vstring(
            "BreitWigner::signalPass(mass, mP[91, 89, 93], sP[3.0, 0.0, 7.0])", ## simple Breit Wigner for pass sample
            "BreitWigner::signalFail(mass, mF[91, 89, 93], sF[3.0, 0.0, 7.0])", ## simple Breit Wigner for fail sample
            "Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])", ## exp 1 parameter
            "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])", ## exp 1 parameter
            "efficiency[0.8,0,1]"
            ),
        pdfSignalPlusBackgroundb0 = cms.vstring(
            "#import templates0.root:w:signalPassMC", ## need a file with a workspace with signalPass RooDataHist
            "#import templates0.root:w:signalFailMC", ## need a file with a workspace with signalFail RooDataHist
            "Gaussian::signalPassSmear(mass, mP[0., -1, 1.], sP[0.01, 0.0, 0.5])", ## Gaussian Convolution
            "Gaussian::signalFailSmear(mass, mF[0., -1, 1.], sF[0.01, 0.0, 0.5])",
            "FCONV::signalPass(mass, signalPassMC, signalPassSmear)",
            "FCONV::signalFail(mass, signalFailMC, signalFailSmear)",
        ),
        pdfSignalPlusBackgroundb1 = cms.vstring(
            "#import templates1.root:w:signalPassMC",
            "#import templates1.root:w:signalFailMC",
            "Gaussian::signalPassSmear(mass, mP[0., -1, 1.], sP[0.01, 0.0, 0.5])",
            "Gaussian::signalFailSmear(mass, mF[0., -1, 1.], sF[0.01, 0.0, 0.5])",
            "FCONV::signalPass(mass, signalPassMC, signalPassSmear)",
            "FCONV::signalFail(mass, signalFailMC, signalFailSmear)",
        ),
        pdfSignalPlusBackgroundb2 = cms.vstring(
            "#import templates2.root:w:signalPassMC",
            "#import templates2.root:w:signalFailMC",
            "Gaussian::signalPassSmear(mass, mP[0., -1, 1.], sP[0.01, 0.0, 0.5])",
            "Gaussian::signalFailSmear(mass, mF[0., -1, 1.], sF[0.01, 0.0, 0.5])",
            "FCONV::signalPass(mass, signalPassMC, signalPassSmear)",
            "FCONV::signalFail(mass, signalFailMC, signalFailSmear)",
        ),
        pdfSignalPlusBackgroundb3 = cms.vstring(
            "#import templates3.root:w:signalPassMC",
            "#import templates3.root:w:signalFailMC",
            "Gaussian::signalPassSmear(mass, mP[0., -1, 1.], sP[0.01, 0.0, 0.5])",
            "Gaussian::signalFailSmear(mass, mF[0., -1, 1.], sF[0.01, 0.0, 0.5])",
            "FCONV::signalPass(mass, signalPassMC, signalPassSmear)",
            "FCONV::signalFail(mass, signalFailMC, signalFailSmear)",
        ),
        pdfSignalPlusBackgroundb4 = cms.vstring(
            "#import templates4.root:w:signalPassMC",
            "#import templates4.root:w:signalFailMC",
            "Gaussian::signalPassSmear(mass, mP[0., -1, 1.], sP[0.01, 0.0, 0.5])",
            "Gaussian::signalFailSmear(mass, mF[0., -1, 1.], sF[0.01, 0.0, 0.5])",
            "FCONV::signalPass(mass, signalPassMC, signalPassSmear)",
            "FCONV::signalFail(mass, signalFailMC, signalFailSmear)",
        ),
        pdfSignalPlusBackgroundb5 = cms.vstring(
            "#import templates5.root:w:signalPassMC",
            "#import templates5.root:w:signalFailMC",
            "Gaussian::signalPassSmear(mass, mP[0., -1, 1.], sP[0.01, 0.0, 0.5])",
            "Gaussian::signalFailSmear(mass, mF[0., -1, 1.], sF[0.01, 0.0, 0.5])",
            "FCONV::signalPass(mass, signalPassMC, signalPassSmear)",
            "FCONV::signalFail(mass, signalFailMC, signalFailSmear)",
        ),
    ),

    Efficiencies = cms.PSet(
        IDModule
    )
)

if options.backgroundType == 0:

    process.muonIdTnP.PDFs.pdfSignalPlusBackgroundb0 += cms.vstring("Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
                                                                    "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
                                                                    "efficiency[0.8,0,1]")
    process.muonIdTnP.PDFs.pdfSignalPlusBackgroundb1 += cms.vstring("Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
                                                                    "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
                                                                    "efficiency[0.8,0,1]")
    process.muonIdTnP.PDFs.pdfSignalPlusBackgroundb2 += cms.vstring("Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
                                                                    "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
                                                                    "efficiency[0.8,0,1]")
    process.muonIdTnP.PDFs.pdfSignalPlusBackgroundb3 += cms.vstring("Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
                                                                    "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
                                                                    "efficiency[0.8,0,1]")
    process.muonIdTnP.PDFs.pdfSignalPlusBackgroundb4 += cms.vstring("Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
                                                                    "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
                                                                    "efficiency[0.8,0,1]")
    process.muonIdTnP.PDFs.pdfSignalPlusBackgroundb5 += cms.vstring("Exponential::backgroundPass(mass, cP[0.0, -1.0, 0.0])",
                                                                    "Exponential::backgroundFail(mass, cF[0.0, -1.0, 0.0])",
                                                                    "efficiency[0.8,0,1]")

elif options.backgroundType == 1:
    process.muonIdTnP.PDFs.pdfSignalPlusBackgroundb0 += cms.vstring("Chebychev::backgroundPass(mass, {c1P[0.0, -1.0, 0.0], c2P[0.0, -1.0, 1.0]})",
                                                                    "Chebychev::backgroundFail(mass, {c1F[0.0, -1.0, 0.0], c2F[0.0, -1.0, 1.0]})");
    process.muonIdTnP.PDFs.pdfSignalPlusBackgroundb1 += cms.vstring("Chebychev::backgroundPass(mass, {c1P[0.0, -1.0, 0.0], c2P[0.0, -1.0, 1.0]})",
                                                                    "Chebychev::backgroundFail(mass, {c1F[0.0, -1.0, 0.0], c2F[0.0, -1.0, 1.0]})");
    process.muonIdTnP.PDFs.pdfSignalPlusBackgroundb2 += cms.vstring("Chebychev::backgroundPass(mass, {c1P[0.0, -1.0, 0.0], c2P[0.0, -1.0, 1.0]})",
                                                                    "Chebychev::backgroundFail(mass, {c1F[0.0, -1.0, 0.0], c2F[0.0, -1.0, 1.0]})");
    process.muonIdTnP.PDFs.pdfSignalPlusBackgroundb3 += cms.vstring("Chebychev::backgroundPass(mass, {c1P[0.0, -1.0, 0.0], c2P[0.0, -1.0, 1.0]})",
                                                                    "Chebychev::backgroundFail(mass, {c1F[0.0, -1.0, 0.0], c2F[0.0, -1.0, 1.0]})");
    process.muonIdTnP.PDFs.pdfSignalPlusBackgroundb4 += cms.vstring("Chebychev::backgroundPass(mass, {c1P[0.0, -1.0, 0.0], c2P[0.0, -1.0, 1.0]})",
                                                                    "Chebychev::backgroundFail(mass, {c1F[0.0, -1.0, 0.0], c2F[0.0, -1.0, 1.0]})");
    process.muonIdTnP.PDFs.pdfSignalPlusBackgroundb5 += cms.vstring("Chebychev::backgroundPass(mass, {c1P[0.0, -1.0, 0.0], c2P[0.0, -1.0, 1.0]})",
                                                                    "Chebychev::backgroundFail(mass, {c1F[0.0, -1.0, 0.0], c2F[0.0, -1.0, 1.0]})");

else:
    sys.exit('options not recognized --> exit');

process.fit = cms.Path(process.muonIdTnP)
