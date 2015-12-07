#!/usr/bin/env python

import os,sys
import optparse
import commands
import time
import random

usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,    dest='queue'              , help='batch queue'                                        , default='1nh')
parser.add_option('-n', '--njobs'      ,    dest='njobs'              , help='number of jobs'                                     , default=1,  type=int)
parser.add_option('-o', '--out'        ,    dest='output'             , help='output directory'                                   , default='/data_CMS2/cms/SLHC22_DYToEE_PU140')
parser.add_option('-c', '--cfg'        ,    dest='cfg'                , help='cfg file'                                           , default='DYToEE_PU140_RECO_step3_cfg.py')
parser.add_option('-r', '--rep'        ,    dest='customReplacements' , help='sed replacements for cfg  key1:val1,key2:val2,...'  , default=None)
(opt, args) = parser.parse_args()


#prepare output
cmsswBase=os.environ['CMSSW_BASE']
jobsDir='/data_CMS2/cms/FARM%s'%(time.time())
os.system('mkdir -p %s'%jobsDir)

def replfunc(match):
    return repldict[match.group(0)]


#loop over the required number of jobs
for n in xrange(0,opt.njobs):

    jobSeed=n+1
    # here increment the jobSeed from a fixed number
    #jobSeed=n+74
    
    
    #sed the cfg template 
    inCfg = open(opt.cfg).read()
    outCfg = open('%s/cmssw_%d_cfg.py'%(jobsDir,jobSeed), 'w')
    replacements = {'XXX_SEED_XXX':str(jobSeed),'XXX_SKIP_XXX':str(jobSeed*50 - 50),'XXX_RND_XXX':str(random.randint(1,100000000))}
    if opt.customReplacements is not None:
        for rep in opt.customReplacements.split(','):
            repKeys=rep.split(':')
            replacements[repKeys[0]]=repKeys[1]
    
    for i in replacements.keys():
        inCfg = inCfg.replace(i, replacements[i])
    outCfg.write(inCfg)
    outCfg.close()
    
    
    #create a wrapper for standalone cmssw job
    scriptFile = open('%s/runJob_%d.sh'%(jobsDir,jobSeed), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('export X509_USER_PROXY=/home/llr/cms/salerno/.t3/proxy.cert\n')
    scriptFile.write('source /cvmfs/cms.cern.ch/cmsset_default.sh\n')
    scriptFile.write('export SCRAM_ARCH=slc6_amd64_gcc472\n')
    scriptFile.write('cd /home/llr/cms/salerno/CMSSW/HGCAL/CMSSW_6_2_0_SLHC22/src\n')
    scriptFile.write('eval `scram r -sh`\n')
    scriptFile.write('cd %s\n'%jobsDir)
    scriptFile.write('cmsRun cmssw_%d_cfg.py\n'%jobSeed)
    scriptFile.write('cp DYToEE_PU140_RECO_%d.root %s\n'%(jobSeed,opt.output))
    scriptFile.write('rm DYToEE_PU140_RECO_%d.root\n'%jobSeed)    
    scriptFile.write('echo "All done for job %d" \n'%jobSeed)
    scriptFile.close()
    os.system('chmod u+rwx %s/runJob_%d.sh'%(jobsDir,jobSeed))

    #submit it to the batch 
    os.system("/opt/exp_soft/cms/t3/t3submit -k eo -q cms \'%s/runJob_%d.sh\'"%(jobsDir,jobSeed))
    
