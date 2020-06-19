g++ -I $(pwd) -Wall $(root-config --cflags --libs) $(~/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins) -o doHistos Root/doHistos.C
