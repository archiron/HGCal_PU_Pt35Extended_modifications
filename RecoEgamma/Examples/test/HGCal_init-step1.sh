#!/bin/sh
# This file is called ./HGCal_step1.sh

LOG_SOURCE=/afs/cern.ch/work/a/archiron/private/CMSSW_6_2_0_SLHC25_integration/src/

cd $LOG_SOURCE
eval `scramv1 runtime -sh`
cd -
for (( i=0; i<10; i++ ))
do
   bsub -q 8nh HGCal-step1.sh $i
done


