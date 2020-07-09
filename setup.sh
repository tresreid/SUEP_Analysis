# Export project base directory
export SUEP_BASE=$PWD

# Setup environment
source /cvmfs/sft.cern.ch/lcg/views/LCG_97rc4/x86_64-centos7-gcc9-opt/setup.sh

# Create necessary directories
if [ ! -d "$SUEP_BASE/output" ]
then
  mkdir $SUEP_BASE/output
fi
if [ ! -d "$SUEP_BASE/plots" ]
then
  mkdir $SUEP_BASE/plots
fi

# Install fastjet locally
if [ ! -d "$SUEP_BASE/fastjet-install"]
then
  curl -O http://fastjet.fr/repo/fastjet-3.3.4.tar.gz
  tar zxvf fastjet-3.3.4.tar.gz
  cd $SUEP_BASE/fastjet-3.3.4
  ./configure --prefix=$SUEP_BASE/fastjet-install
  make -j 8
  make check
  make install
  cd $SUEP_BASE
fi
