#include "TSystem.h"
#include "TROOT.h"

#include <string.h>
#include <cstdlib> //as stdlib.h    
#include <cstdio>

using namespace std;

// usage :
// root -l -q -b 'runMakePhotonPurity.C++("argument1","argument2",...)' &> logfile.log  

// &> logfile.log is useful if you want to save output in a log for debugging
// You call this macro with the same arguments as makePhotonPurityFit.C, but using only strings (to build command in an easier way)
// e.g., if you pass a number as the lumi or a bool as addSystematics, you write these as if they were strings


void runMakePhotonPurity(const string& inputDirectory = "", // directory with dataFiles (file names must match with RunEra: all matched files, and only them, are used)    
                         const string& lumi = "36.2", // luminosity                                                                      
                         const string& outputDIR = "",
			 const string& addSystematics = "false",
                         const string& inputDirectorySignalMC = "",   // here there must be only MC signal root files                                     
                         const string& inputDirectoryBackgroundMC = "",  // here there must be only MC background root files              
                         const string& makeFitBasedOnlyOnTemplates = "false",
			 const string& uniformIsoBinning = "true"
                         ){


  string command = "makePhotonPurityFit(\"" + inputDirectory + "\"," + lumi + ",\"" + outputDIR + "\"," + addSystematics + ",\"" + inputDirectorySignalMC + "\",\"" + inputDirectoryBackgroundMC + "\"," + makeFitBasedOnlyOnTemplates + "," + uniformIsoBinning + ")";

  gROOT->ProcessLine(".L makePhotonPurityPDFs.cc+");
  gSystem->Load("makePhotonPurityPDFs_cc.so");
  gROOT->ProcessLine(".L makePhotonPurityFit.C");
  gROOT->ProcessLine(command.c_str());


}
