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


#root -b -q 'processNtuple.C("TTCMS3_gt20_2basicDeltaBeta"  ,"/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 2)' 
#root -b -q 'processNtuple.C("TTCMS3_gt20_3DeltaBetaDR"  ,"/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 3)' 
#root -b -q 'processNtuple.C("TTCMS3_gt20_5uncorr"  ,"/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 5)' 
#root -b -q 'processNtuple.C("TTCMS3_gt20_6DeltaBetaLogPtDR"  ,"/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 6)' 

root -b -q 'processNtuple.C("tt532_2basicDeltaBeta","/tas/gzevi/files/babies/RelValTTbar/CMSSW_5_3_6-PU_START53_V14-v1/GEN-SIM-RECO/ntuple*.root", 0, 2)' 
#root -b -q 'processNtuple.C("TTCMS3_2basicDeltaBeta"  ,"/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 2)' 
#root -b -q 'processNtuple.C("TTCMS3_25ns_2basicDeltaBeta"  ,"/tas/gzevi/files/miniAOD/TT_Tune4C_13TeV-pythia8-tauola-Flat20to50_POSTLS170_V5-v1_Big_CMS3.root ", 0, 2)' 





#root -b -q 'processNtuple.C("TTCMS3_12basicDeltaBetaMap"  ,"/tas/gzevi/miniAOD/CMSSW_7_0_4/src/CMS2/NtupleMaker/test/TT_Tune4C_13TeV-pythia8-tauola-PU_S14_PAT_big_CMS3.root", 0, 12)' 