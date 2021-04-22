#!/bin/bash

#dir=/uscms/home/mreid/nobackup/sueps/CMSSW_10_6_13/src/SUEP_Analysis/jet_clustering/qcd_output
#dir=/uscms/home/mreid/nobackup/sueps/CMSSW_10_6_13/src/SUEP_Analysis/jet_clustering/ROC
#dir=/uscms/home/mreid/nobackup/sueps/CMSSW_10_6_13/src/SUEP_Analysis/jet_clustering/DIST
#dir=/uscms/home/mreid/nobackup/sueps/CMSSW_10_6_13/src/SUEP_Analysis/jet_clustering/COMBO
#dir=/uscms/home/mreid/nobackup/sueps/CMSSW_10_6_13/src/SUEP_Analysis/jet_clustering/DIST2d
dir=/uscms/home/mreid/nobackup/sueps/CMSSW_10_6_13/src/SUEP_Analysis/jet_clustering/Plots
for d in $(find ${dir} -type d)
do
  #Do something, the directory is accessible with $d:
  echo $d
  python3 $dir/../../make_html_listing.py $d
done
#outdir=/publicweb/m/mreid/SUEPs/darkPhoHad
#outdir=/publicweb/m/mreid/SUEPs/QCD
#outdir=/publicweb/m/mreid/SUEPs/COMBO
#outdir=/publicweb/m/mreid/SUEPs/ROC
#outdir=/publicweb/m/mreid/SUEPs/DIST
#outdir=/publicweb/m/mreid/SUEPs/DIST2d
#outdir=/publicweb/m/mreid/SUEPs/JetAlgoComparisons_baseline_pt/fullbox
#outdir=/publicweb/m/mreid/SUEPs/JetAlgoComparisons_baseline_multi/pure_containment
#outdir=/publicweb/m/mreid/SUEPs/mass_test
#outdir=/publicweb/m/mreid/SUEPs/testdir
#outdir=/publicweb/m/mreid/SUEPs/JetAlgoComparisons_baseline_ptv2
#outdir=/publicweb/m/mreid/SUEPs/JetAlgoComparisons_baseline_multiv2/scalar
outdir=/publicweb/m/mreid/SUEPs/nsubjettiness
rm -r $outdir
mkdir $outdir
cp -r $dir/* $outdir

