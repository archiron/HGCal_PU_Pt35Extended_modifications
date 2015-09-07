#!/bin/sh
# This file is called ./HGCal.sh

if [ "$1" == "" ] 
then
	echo "no argument"
	exit
fi

LOG_SOURCE=/afs/cern.ch/work/a/archiron/private/CMSSW_6_2_0_SLHC25_integration/src/

cd $LOG_SOURCE
eval `scramv1 runtime -sh`
cd -
echo "valeur : $1"
cmsRun /afs/cern.ch/work/a/archiron/private/CMSSW_6_2_0_SLHC25_integration/src/RecoEgamma/Examples/test/ElectronsPt35-2023HGCALUpg_SLHC25_step1_cfg.py $1


