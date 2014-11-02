#!/bin/bash

function run () {
    echo root -b -q mergeHadoopFiles.C\(\"${HADOOPDIR}/${TAG}_$1/\",\"${OUTPUTDIR}/$1.root\"\)
    root -b -q mergeHadoopFiles.C\(\"${HADOOPDIR}/${TAG}_$1/\",\"${OUTPUTDIR}/$1.root\"\) >& log_merge_${TAG}_$1.txt &
}

function roothadd () {
    rm ${OUTPUTDIR}/${TAG}/$1.root
    mkdir ${OUTPUTDIR}/${TAG}
    echo "hadd -k ${OUTPUTDIR}/$1.root ${HADOOPDIR}/${TAG}_$1/*.root"
    hadd -k ${OUTPUTDIR}/${TAG}/$1.root ${HADOOPDIR}/${TAG}_$1/*.root
    # now remake the efficiency plot
    echo root -b -q makeEffMergedFiles.C\(\"${OUTPUTDIR}/${TAG}/$1.root\"\)
    root -b -q makeEffMergedFiles.C\(\"${OUTPUTDIR}/${TAG}/$1.root\"\)
}


TAG=V00-00-04

HADOOPDIR=/hadoop/cms/store/user/${USER}/genFakeHistos/
OUTPUTDIR=/home/users/gzevi/FakeLepton/genCMS2/ElectronLooperCMS2/looper/batchsubmit/merged/

mkdir -p $OUTPUTDIR
chmod -R a+wrx $OUTPUTDIR


roothadd ttpythia
#roothadd ttsemilep1
#roothadd ttsemilep2
roothadd ttsemilep*
mv merged/${TAG}/ttsemilep\*.root merged/${TAG}/ttsemilep.root
roothadd qcdmuenriched
roothadd bjets




