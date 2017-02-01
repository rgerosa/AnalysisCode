import os
import FWCore.ParameterSet.Config as cms
from CondCore.DBCommon.CondDBSetup_cfi import *
from CondCore.CondDB.CondDB_cfi import*

## setup JEC on PAT jets from miniAOD
def JECConfiguration(process,usePrivateSQlite,JECEra,isMC,applyL2L3Residuals,isCrab):

	## if true look to a local SQLite file instead of GT entry
	if usePrivateSQlite:
		era_file    = JECEra;
		era_payload = JECEra;
		if isMC : 
   		 	era_file += "_MC";
			era_payload += "_MC";
  		else :
			era_file += "_DATA";
			if JECEra == "Spring16_25nsV10": 
				era_payload = "Spring16_25nsV10All_DATA";
			elif JECEra == "Spring16_23Sep2016V1":
				era_payload = "Spring16_23Sep2016AllV1_DATA";
			elif JECEra == "Summer16_23Sep2016V3":
				era_payload = "Summer16_23Sep2016AllV3_DATA";
			else:
				era_payload += "_DATA";
			
  		dBFile = os.path.expandvars(era_file+".db")


		## connect to local SQLite file
		process.jec = cms.ESSource("PoolDBESSource",
					   CondDBSetup,
					   connect = cms.string("sqlite_file:../../data/JEC/"+dBFile),
					   toGet =  cms.VPSet(
				## AK4PF corrections
				cms.PSet(
					record = cms.string("JetCorrectionsRecord"),
					tag = cms.string("JetCorrectorParametersCollection_"+era_payload+"_AK4PF"),
					label= cms.untracked.string("AK4PF")
					),
				## AK4PFchs corrections
				cms.PSet(
					record = cms.string("JetCorrectionsRecord"),
					tag = cms.string("JetCorrectorParametersCollection_"+era_payload+"_AK4PFchs"),
					label= cms.untracked.string("AK4PFchs")
					),
				## AK4PFPUPPI corrections
				cms.PSet(
					record = cms.string("JetCorrectionsRecord"),
					tag = cms.string("JetCorrectorParametersCollection_"+era_payload+"_AK4PFPuppi"),
					label= cms.untracked.string("AK4PFPuppi")
					),		
				## AK8PFchs corrections            
				cms.PSet(
					record = cms.string("JetCorrectionsRecord"),
					tag = cms.string("JetCorrectorParametersCollection_"+era_payload+"_AK8PFchs"),
					label= cms.untracked.string("AK8PFchs")
					),	
				## AK8PFPuppi corrections            
				cms.PSet(
					record = cms.string("JetCorrectionsRecord"),
					tag = cms.string("JetCorrectorParametersCollection_"+era_payload+"_AK8PFPuppi"),
					label= cms.untracked.string("AK8PFPuppi")
					),	
				)
					   )

		if isCrab:
			process.jec.connect = cms.string("sqlite_file:src/AnalysisCode/MonoXAnalysis/data/JEC/"+dBFile)

     	## give preference wrt to GT
	if usePrivateSQlite == True:
	     	process.es_prefer_jec = cms.ESPrefer("PoolDBESSource",'jec')
		

	# JEC levels when redoing jets and MET
	if isMC:
		process.JECLevels = cms.PSet(labels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute'))
	else :
		if not applyL2L3Residuals : 
			process.JECLevels = cms.PSet(labels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute'))
		else : 
			process.JECLevels = cms.PSet(labels = cms.vstring('L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'))
 
	
