#include "TSystem.h"
#include "TROOT.h"

#include <string.h>
#include <cstdlib> //as stdlib.h    
#include <cstdio>

using namespace std;

void runMakePhotonPurity(const string& inputDirectory = "", // directory with dataFiles (file names must match with RunEra: all matched files, and only them, are used)    
                         const string& lumi = "36.2", // luminosity                                                                      
                         const string& outputDIR = "",
			 const string& addSystematics = "false",
                         const string& inputDirectorySignalMC = "",   // here there must be only MC signal root files                                     
                         const string& inputDirectoryBackgroundMC = "",  // here there must be only MC background root files              
                         const string& makeFitBasedOnlyOnTemplates = "false"
                         ){


  string command = "makePhotonPurityFit(\"" + inputDirectory + "\"," + lumi + ",\"" + outputDIR + "\"," + addSystematics + ",\"" + inputDirectorySignalMC + "\",\"" + inputDirectoryBackgroundMC + "\"," + makeFitBasedOnlyOnTemplates + ")";

  gROOT->ProcessLine(".L makePhotonPurityPDFs.cc+");
  gSystem->Load("makePhotonPurityPDFs_cc.so");
  gROOT->ProcessLine(".L makePhotonPurityFit.C");
  gROOT->ProcessLine(command.c_str());


}
