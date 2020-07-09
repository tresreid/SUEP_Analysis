g++ -I $SUEP_BASE -Wno-deprecated $(root-config --cflags --libs) $($SUEP_BASE/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins) -o doHistos Root/doHistos.C
