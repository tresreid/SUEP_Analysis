#!/bin/bash

./test 0 0
./test 1 0
./test 2 0

#./test 3 0
#./test 4 0
#./test 5 0
#./test 6 0
#./test 7 0
#./test 8 0
#
#./test 3 1
#./test 4 1
#./test 5 1
#./test 6 1
#./test 7 1
#./test 8 1

cat data/qcd_300_v1.txt >> data/qcd_300_v0.txt
cat data/qcd_500_v1.txt >> data/qcd_500_v0.txt
cat data/qcd_700_v1.txt >> data/qcd_700_v0.txt
cat data/qcd_1000_v1.txt >> data/qcd_1000_v0.txt
cat data/qcd_1500_v1.txt >> data/qcd_1500_v0.txt
cat data/qcd_2000_v1.txt >> data/qcd_2000_v0.txt

