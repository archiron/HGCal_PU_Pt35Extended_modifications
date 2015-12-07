#!/bin/sh
# This file is called ./HGCal_step2.sh

LOG_SOURCE=/afs/cern.ch/work/a/archiron/private/CMSSW_6_2_0_SLHC25_integration/src/

cd $LOG_SOURCE
eval `scramv1 runtime -sh`
cd -
for (( i=10; i<20; i++ ))
do
   bsub -q 1nd HGCal-step2.sh $i
#   ./HGCal-step2.sh $i
done


