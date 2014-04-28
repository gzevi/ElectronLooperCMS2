#root -b -q 'processNtuple.C("tt532","/tas/gzevi/files/babies/RelValTTbar/CMSSW_5_3_6-PU_START53_V14-v1/GEN-SIM-RECO/ntuple*.root")' > tt5.log
#root -b -q 'processNtuple.C("Zee532","/tas/gzevi/files/babies/RelValZEE/CMSSW_5_3_6-PU_START53_V14-v1/GEN-SIM-RECO/ntuple*.root")' > Zee5.log
#root -b -q 'processNtuple.C("tt700pre12","/tas/gzevi/files/babies/RelValTTbar/CMSSW_7_0_0_pre12-PU_START70_V5/GEN-SIM-RECO/ntuple*.root")' > tt7.log 
#root -b -q 'processNtuple.C("tt700pre12patch","/tas/gzevi/files/babies/RelValTTbar/CMSSW_7_0_0_pre12-PU_START70_V5/GEN-SIM-RECO/TestLindsey/ntuple*.root")' > tt7b.log 
#root -b -q 'processNtuple.C("Zee700pre12","/tas/gzevi/files/babies/RelValZEE/CMSSW_7_0_0_pre12-PU_START70_V5/GEN-SIM-RECO/ntuple*.root")' > Zee7.log
#root -b -q 'processNtuple.C("Zee700pre12patch","/tas/gzevi/files/babies/RelValZEE/CMSSW_7_0_0_pre12-PU_START70_V5/GEN-SIM-RECO/TestLindsey/ntuple*.root")' > Zee7b.log

root -b -q 'processNtuple.C("ZeePlusTT700","/tas/gzevi/files/babies/ttPlusZ700/ntuple*.root")' 
root -b -q 'processNtuple.C("ZeePlusTT700patch","/tas/gzevi/files/babies/ttPlusZpatch/ntuple*.root")' 

#root -b -q 'processNtuple.C("ZeePlusTT536_10to20_BasicLoop","/tas/gzevi/files/babies/ttPlusZ536/ntuple*.root", 0, 0)' 
#root -b -q 'processNtuple.C("ZeePlusTT536_10to20_AODarea","/tas/gzevi/files/babies/ttPlusZ536/ntuple*.root", 0, 1)'
#root -b -q 'processNtuple.C("ZeePlusTT536_10to20_DeltaBetaSimple","/tas/gzevi/files/babies/ttPlusZ536/ntuple*.root", 0, 2)' 
#root -b -q 'processNtuple.C("ZeePlusTT536_10to20_DeltaBetaWeights","/tas/gzevi/files/babies/ttPlusZ536/ntuple*.root", 0, 3)' 
#root -b -q 'processNtuple.C("ZeePlusTT536_10to20_DeltaBetaWeights01","/tas/gzevi/files/babies/ttPlusZ536/ntuple*.root", 0, 4)' 

