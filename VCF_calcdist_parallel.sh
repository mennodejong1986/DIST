echo "Starting runs..."

./VCF_calcdist.sh 1 5 &
./VCF_calcdist.sh 6 10 &
./VCF_calcdist.sh 11 15 &
./VCF_calcdist.sh 16 20 &
./VCF_calcdist.sh 21 25 &
./VCF_calcdist.sh 26 30 &
./VCF_calcdist.sh 31 35 &
./VCF_calcdist.sh 36 40 &
./VCF_calcdist.sh 41 45 &
./VCF_calcdist.sh 46 50 &
./VCF_calcdist.sh 51 55 &
./VCF_calcdist.sh 56 60 &
./VCF_calcdist.sh 61 65 &
./VCF_calcdist.sh 66 70 &
./VCF_calcdist.sh 71 75 &
./VCF_calcdist.sh 76 80 &
./VCF_calcdist.sh 81 85 &
./VCF_calcdist.sh 86 90 &
./VCF_calcdist.sh 91 95 &
./VCF_calcdist.sh 96 100 &
./VCF_calcdist.sh 101 105 &
./VCF_calcdist.sh 106 110 &
./VCF_calcdist.sh 111 115 &
./VCF_calcdist.sh 116 120 &
./VCF_calcdist.sh 121 123 &

wait
echo "All analyses finished."