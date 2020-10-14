#!/bin/bash

export SUEPS_DIR=`pwd`

rm $SUEPS_DIR/sueps.tar
cd $SUEPS_DIR/../
tar -cf sueps.tar $SUEPS_DIR
mv $SUEPS_DIR/../sueps.tar $SUEPS_DIR/
cd $SUEPS_DIR

condor_submit jobs.jdl
