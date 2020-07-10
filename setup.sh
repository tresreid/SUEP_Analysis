# Export project base directory
export SUEP_BASE=$PWD

# Setup environment
source /cvmfs/sft.cern.ch/lcg/views/LCG_97rc4/x86_64-centos7-gcc9-opt/setup.sh

# Create necessary directories
if [ ! -d "$SUEP_BASE/output" ]
then
  mkdir $SUEP_BASE/output
  echo "Created $SUEP_BASE/output directory."
fi
if [ ! -d "$SUEP_BASE/plots" ]
then
  mkdir $SUEP_BASE/plots
  echo "Created $SUEP_BASE/plots directory."
fi

# Install fastjet locally
function FASTJET {
  curl -O http://fastjet.fr/repo/fastjet-3.3.4.tar.gz
  tar zxvf fastjet-3.3.4.tar.gz
  cd $SUEP_BASE/fastjet-3.3.4
  ./configure --prefix=$SUEP_BASE/fastjet-install
  make -j 8
  make check
  make install
  cd $SUEP_BASE
  exit 1
}

while getopts "f" OPT; do
  case $OPT in
    f)
      FASTJET
      ;;
  esac
done
