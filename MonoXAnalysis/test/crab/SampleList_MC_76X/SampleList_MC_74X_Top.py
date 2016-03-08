## TTbar: LO Madgraph inclusive, NLO amcatnlo, powheg                                                                                                                         
#samples['TTJets_madgraph']  = ['/TTJets_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',                                 #                               ['useLHEWeights=True','addQCDPDFWeights=True','isSignalSample=True','crossSection=831.76']]                                                      
samples['TTJets_amcatnlo']  = ['/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v3/MINIAODSIM',
                               ['useLHEWeights=True','addQCDPDFWeights=True','isSignalSample=True','crossSection=831.76']]
samples['TTJets_powheg']    = ['/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
                               ['useLHEWeights=True','addQCDPDFWeights=True','isSignalSample=True','crossSection=831.76']]
## Single-top: t-channel                                                                                                                                                     
samples['STop_t-channel_powheg']     = ['/ST_t-channel_top_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
                                        ['useLHEWeights=False','addQCDPDFWeights=False','isSignalSample=True','crossSection=44.0748']]
samples['SAntiTop_t-channel_powheg'] = ['/ST_t-channel_antitop_4f_leptonDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
                                        ['useLHEWeights=False','addQCDPDFWeights=False','isSignalSample=True','crossSection=26.2278']]
## Single-top: t-channel                                                                                                                                                       
samples['STop_s-channel_amcatnlo'] = ['/ST_s-channel_4f_leptonDecays_13TeV-amcatnlo-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
                                     ['useLHEWeights=True','addQCDPDFWeights=True','isSignalSample=True','crossSection=3.34']]
## Single-top: tW channel                                                                                                                                                      
samples['STop_tW-channel_powheg']     = ['/ST_tW_top_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v2/MINIAODSIM',
                                         ['useLHEWeights=True','addQCDPDFWeights=True','isSignalSample=True','crossSection=35.6']]
samples['SAntiTop_tW-channel_powheg'] = ['/ST_tW_antitop_5f_inclusiveDecays_13TeV-powheg-pythia8_TuneCUETP8M1/RunIISpring15MiniAODv2-74X_mcRun2_asymptotic_v2-v1/MINIAODSIM',
                                         ['useLHEWeights=True','addQCDPDFWeights=True','isSignalSample=True','crossSection=35.6']]
