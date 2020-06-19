# Define path variables
export SUEP_ANALYSIS=`pwd`
if [ -z $FASTJET_PATH ]
then
  export FASTJET_PATH=~/fastjet-install
fi

# Compile doHistos
g++ -I $SUEP_ANALYSIS -Wall $(root-config --cflags --libs) $($FASTJET_PATH/bin/fastjet-config \
    --cxxflags --libs --plugins) -o doHistos Root/doHistos.C
