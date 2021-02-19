#g++ -I $PWD $(root-config --cflags --libs) $($PWD/../../../../fastjet-install/bin/fastjet-config --cxxflags --libs --plugins) -fopenmp -O3 -lm -lgomp -o test jetAlgoCompv2.cc  -g
g++ -I $PWD $(root-config --cflags --libs) $($PWD/../../../../fastjet-install/bin/fastjet-config --cxxflags --libs --plugins) -fopenmp -O3 -lm -lgomp -o test clusterTest.cc  -g
