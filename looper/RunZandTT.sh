#root -b -q 'processNtuple.C("tt532","/tas/gzevi/files/babies/RelValTTbar/CMSSW_5_3_6-PU_START53_V14-v1/GEN-SIM-RECO/ntuple*.root")' > tt5.log
#root -b -q 'processNtuple.C("Zee532","/tas/gzevi/files/babies/RelValZEE/CMSSW_5_3_6-PU_START53_V14-v1/GEN-SIM-RECO/ntuple*.root")' > Zee5.log
#root -b -q 'processNtuple.C("tt700pre12","/tas/gzevi/files/babies/RelValTTbar/CMSSW_7_0_0_pre12-PU_START70_V5/GEN-SIM-RECO/ntuple*.root")' > tt7.log 
#root -b -q 'processNtuple.C("Zee700pre12","/tas/gzevi/files/babies/RelValZEE/CMSSW_7_0_0_pre12-PU_START70_V5/GEN-SIM-RECO/ntuple*.root")' > Zee7.log

#root -b -q 'processNtuple.C("tt700pre12patch","/tas/gzevi/files/babies/RelValTTbar/CMSSW_7_0_0_pre12-PU_START70_V5/GEN-SIM-RECO/TestLindsey/ntuple*.root")' > tt7b.log 
#root -b -q 'processNtuple.C("Zee700pre12patch","/tas/gzevi/files/babies/RelValZEE/CMSSW_7_0_0_pre12-PU_START70_V5/GEN-SIM-RECO/TestLindsey/ntuple*.root")' > Zee7b.log

#root -b -q 'processNtuple.C("ZeePlusTT700","/tas/gzevi/files/babies/ttPlusZ700/ntuple*.root", 0, 1)' 
#root -b -q 'processNtuple.C("ZeePlusTT700patch","/tas/gzevi/files/babies/ttPlusZpatch/ntuple*.root", 0, 1)' 

#root -b -q 'processNtuple.C("ZeePlusTT536_10to20_BasicLoop","/tas/gzevi/files/babies/ttPlusZ536/ntuple*.root", 0, 0)' 
#root -b -q 'processNtuple.C("ZeePlusTT536_10to20_AODarea","/tas/gzevi/files/babies/ttPlusZ536/ntuple*.root", 0, 1)'
#root -b -q 'processNtuple.C("ZeePlusTT536_10to20_DeltaBetaSimple","/tas/gzevi/files/babies/ttPlusZ536/ntuple*.root", 0, 2)' 
#root -b -q 'processNtuple.C("ZeePlusTT536_10to20_DeltaBetaWeights","/tas/gzevi/files/babies/ttPlusZ536/ntuple*.root", 0, 3)' 
#root -b -q 'processNtuple.C("ZeePlusTT536_10to20_DeltaBetaWeights01","/tas/gzevi/files/babies/ttPlusZ536/ntuple*.root", 0, 4)' 

#/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root
#/tas/gzevi/files/miniAOD/DYJetsToLL_M-50_13TeV-madgraph-pythia8_PAT_CMS3.root
#/tas/gzevi/files/miniAOD/TT_Tune4C_13TeV-pythia8-tauola-Flat20to50_POSTLS170_V5-v1_Big_CMS3.root 
#root -b -q 'processNtuple.C("TTCMS3_noCorr_noMap","/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 5)' 
#root -b -q 'processNtuple.C("TTCMS3_noCorr_map"  ,"/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 15)' 

#root -b -q 'processNtuple.C("TTCMS3_noIsoLoop"  ,"/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 0)' 


#root -b -q 'processNtuple.C("TTCMS3_40bx50_2basicDeltaBeta_ROC2_withHF_LARGEdZ"  ,"/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 2)' 
#root -b -q 'processNtuple.C("TTCMS3_40bx50_3DeltaBetaDR"  ,"/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 3)' 
#root -b -q 'processNtuple.C("TTCMS3_40bx50_5uncorr"  ,"/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 5)' 
#root -b -q 'processNtuple.C("TTCMS3_40bx50_6DeltaBetaLogPtDR"  ,"/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 6)' 

#root -b -q 'processNtuple.C("tt532_2basicDeltaBeta_ROC_withHF","/tas/gzevi/files/babies/RelValTTbar/CMSSW_5_3_6-PU_START53_V14-v1/GEN-SIM-RECO/ntuple*.root", 0, 2)' 
##root -b -q 'processNtuple.C("tt532_1area","/tas/gzevi/files/babies/RelValTTbar/CMSSW_5_3_6-PU_START53_V14-v1/GEN-SIM-RECO/ntuple*.root", 0, 1)' 

#root -b -q 'processNtuple.C("TTCMS3_25bx2050_2basicDeltaBeta"  ,"/tas/gzevi/files/miniAOD/TT_Tune4C_13TeV-pythia8-tauola-Flat20to50_POSTLS170_V5-v1_Big_CMS3.root ", 0, 2)' 

#root -b -q 'processNtuple.C("TTCMS3_20bx25_2basicDeltaBeta_stdIDstdISOROC_withHF"  ,"/afs/cern.ch/user/g/gzevi/work/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_PU20bx25_POSTLS170_V5-v1.root", 0, 2)'
#With PUiso: /tas/gzevi/files/miniAOD/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_PU20bx25_POSTLS170_V5-v1_PUiso.root
#root -b -q 'processNtuple.C("TTCMS3_20bx25_2basicDeltaBeta_ROC2_withHF"  ,"/tas/gzevi/files/miniAOD/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola_PU20bx25_POSTLS170_V5-v1_PUiso.root", 0, 2)'
#root -b -q 'processNtuple.C("dataDoubleElectronD_CMS3_test","/tas/gzevi/files/miniAOD/DoubleEleMiniAOD.root",0,2)'


# OOT PU REDUCTION
#3D timing (Scenario2) /tas/gzevi/CMSSW_7_1_3/src/CMS2/NtupleMaker/test/reco_nochi2_10_miniAOD_CMS3_Full.root
#root -b -q 'processNtuple.C("TTCMS3_40bx25_2basicDeltaBeta_test"  ,"/tas/gzevi/CMSSW_7_1_3/src/CMS2/NtupleMaker/test/merged/reco_PF71X_ONLYTIMECUT_miniAOD_CMS3.root", 0, 2)' 
root -b -q 'processNtuple.C("TTCMS3_40bx25_2basicDeltaBeta_OOTPUScenario1"  ,"/tas/gzevi/CMSSW_7_1_3/src/CMS2/NtupleMaker/test/merged/reco_PF71X_ONLYTIMECUT_miniAOD_CMS3*.root", 0, 2)' 
root -b -q 'processNtuple.C("TTCMS3_40bx25_2basicDeltaBeta_OOTPUScenario2"  ,"/tas/gzevi/CMSSW_7_1_3/src/CMS2/NtupleMaker/test/merged/reco_PF71X_NOCHI2_miniAOD_CMS3*.root", 0, 2)' 
root -b -q 'processNtuple.C("TTCMS3_40bx25_2basicDeltaBeta_OOTPUScenario0"  ,"/tas/gzevi/CMSSW_7_1_3/src/CMS2/NtupleMaker/test/merged/reco_PF71X_REALNOTIME_miniAOD_CMS3*.root", 0, 2)' 
