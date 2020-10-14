#!/bin/bash

#dir=/uscms/home/mreid/nobackup/CRAB/signal_region_analysis/CMSSW_10_2_18/src/iDMSkimmer/washAOD
dir=/uscms/home/mreid/nobackup/sueps/CMSSW_10_6_13/src/SUEP_Analysis/plots
for d in $(find ${dir} -type d)
do
  #Do something, the directory is accessible with $d:
  echo $d
  python3 $dir/../make_html_listing.py $d
done
cp -r $dir/* /publicweb/m/mreid/SUEPs/

