#!/bin/bash
echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
source /cvmfs/cms.cern.ch/cmsset_default.sh  ## if a tcsh script, use .csh instead of .sh

tar -xf sueps.tar 
cd SUEP_Analysis
source /cvmfs/sft.cern.ch/lcg/views/LCG_97rc4/x86_64-centos7-gcc9-opt/setup.sh

./doHistos $1
