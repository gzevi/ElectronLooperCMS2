#!/bin/bash

#
# All MT2 related datasets available on hadoop
#

#TAG="V00-00-00"
#TAG="V00-00-01" ## iso with pT > 0.5 GeV 
#TAG="V00-00-02" ## iso with pT > 0 GeV, new plots
#TAG="V00-00-03" ## iso with pT > 0 GeV, rho plots, rho weight for QCD
#TAG="V00-00-04" ## iso with pT > 0 GeV, rho plots, rho weight for QCD, another weight based on SumPt, genJet isolation
#TAG="V00-00-05" ## iso with pT > 0 GeV, rho plots, rho weight for QCD, another weight based on SumPt, genJet isolation, 20 < pT < 25 
TAG="V00-00-06" ## added HT with pT > 2 GeV particles, GenJetHT (qcd only has been run)

#./writeConfig.sh /hadoop/cms/store/user/fgolf/CMS2_V05-03-23/ttbar_pythia6_8tev_fgolf-ttbar_pythia6_8tev-5b7f5e26b2c20e06a23122b93fee9986/  ${TAG}_ttpythia
./writeConfigSep.sh /hadoop/cms/store/user/fgolf/CMS2_V05-03-23/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19_ext1-v1/ ${TAG}_ttsemilep1
./writeConfigSep.sh /hadoop/cms/store/user/fgolf/CMS2_V05-03-23/TTJets_SemiLeptMGDecays_8TeV-madgraph-tauola_Summer12_DR53X-PU_S10_START53_V19_ext2-v1/ ${TAG}_ttsemilep2
./writeConfigSep.sh /hadoop/cms/store/user/fgolf/CMS2_V05-03-23/QCD_Pt_20_MuEnrichedPt_15_TuneZ2star_8TeV_pythia6_Summer12_DR53X-PU_S10_START53_V7A-v3/ ${TAG}_qcdmuenriched
#./writeConfig.sh /hadoop/cms/store/user/gzevi/CMS2_V05-03-32X/BJets_HT-100To250_8TeV-madgraph_Summer12_DR53X-PU_S10_START53_V7A-v1/ ${TAG}_bjets

# --- write submit script ---
mkdir -p configs_${TAG}

mv condor_${TAG}*.cmd configs_${TAG}
echo "#!/bin/bash" > submitAll.sh
echo "voms-proxy-init -voms cms -valid 240:00" >> submitAll.sh
for file in configs_${TAG}/*.cmd
do 
    echo "condor_submit ${file}" >> submitAll.sh
done
chmod +x submitAll.sh
echo "[writeAllConfig] wrote submit script submitAll.sh"
