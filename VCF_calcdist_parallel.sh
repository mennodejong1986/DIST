#!/bin/bash

echo "Starting runs..."

./VCF_calcdist.sh 1 4 &
./VCF_calcdist.sh 5 8 &
./VCF_calcdist.sh 9 12 &
./VCF_calcdist.sh 13 16 &
./VCF_calcdist.sh 17 20 &
./VCF_calcdist.sh 21 24 &
./VCF_calcdist.sh 25 28 &
./VCF_calcdist.sh 29 32 &
./VCF_calcdist.sh 33 36 &
./VCF_calcdist.sh 37 40 &
./VCF_calcdist.sh 41 44 &
./VCF_calcdist.sh 45 48 &
./VCF_calcdist.sh 49 52 &
./VCF_calcdist.sh 53 56 &
./VCF_calcdist.sh 57 60 &
./VCF_calcdist.sh 61 64 &
./VCF_calcdist.sh 65 68 &
./VCF_calcdist.sh 69 72 &
./VCF_calcdist.sh 73 76 &
./VCF_calcdist.sh 77 80 &
./VCF_calcdist.sh 81 84 &
./VCF_calcdist.sh 85 88 &
./VCF_calcdist.sh 89 92 &
./VCF_calcdist.sh 93 96 &
./VCF_calcdist.sh 97 100 &
./VCF_calcdist.sh 101 104 &
./VCF_calcdist.sh 105 108 &
./VCF_calcdist.sh 109 112 &
./VCF_calcdist.sh 113 116 &
./VCF_calcdist.sh 117 120 &
./VCF_calcdist.sh 121 124 &
./VCF_calcdist.sh 125 128 &
./VCF_calcdist.sh 129 132 &
./VCF_calcdist.sh 133 136 &
./VCF_calcdist.sh 137 140 &
./VCF_calcdist.sh 141 144 &
./VCF_calcdist.sh 145 148 &
./VCF_calcdist.sh 149 152 &
./VCF_calcdist.sh 153 153 &

wait
echo "All analyses finished."

echo "Combining..."
cat vcfdist*txt | grep -v 'name1' | sed '1i ind1\tind2\tname1\tname2\tnsites\tnmiss\tpmiss\td_snps\tploidy\tchrom\tstartbp' > allvcfdist.txt

echo "All data has been combined in file 'allvcfdist.txt'."
