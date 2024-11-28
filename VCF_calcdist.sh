#!/bin/bash
# A shell script to calculate pairwise sequence dissimilarity between any two samples in a vcf file.

# INPUT FILE
# As input, use a vcf-file with an unbiased representation of monomorphic and polymorphic sites. So, not a SNP dataset!
# This input file should not be too large, otherwise the computation takes ages, and intermediate files will be huge.
# For mammalian genomes of 3Gb, I typically use a randomly thinned dataset which contains one datapoint every 100bp, generated using the following command:
# ${VCFTOOLS} --thin 100 --gzvcf allsites.globalfilter.vcf.gz -c --recode --recode-INFO-all | gzip > allsites.globalfilter.thinned.100.vcf.gz &
# Also, make sure indels have been removed from the file:
# bcftools view --threads 20 --exclude-types indels allsites.globalfilter.thinned.100.vcf.gz -O z -o allsites.globalfilter.thinned.100.noindels.vcf.gz

# USAGE:
# Execute by typing, for example:
# dos2unix VCF_calcdist.sh
# chmod +x VCF_calcdist.sh
# ./VCF_calcdist.sh 1 10
# This will calculate all pairwise scores for individuals 1 to 10.
# For example, if the dataset contains in total 100 individuals, then the command above will calculate the scores for 10*100 = 1000 pairwise comparisons.

# If you forget to add the two numbers after VCF_calcdist.sh, you will receive the error:
# line 76: ((: i = : syntax error: operand expected (error token is "= ")

# Run the command in two steps: 
# 1. convert the data (convertdata=TRUE,run_loop=FALSE), which will generate the files 'myinput.samples.txt', 'myinput.samples.txt' and 'myinput.pos.txt' 
# 2. afterwards run the actual analysis (convertdata=FALSE,run_loop=TRUE)

# PARALLEL USAGE
# To speed up the calculations, run simulatenously for subsets as such:
# ./VCF_calcdist.sh 1 10 &		
# ./VCF_calcdist.sh 11 20 &
# IMPORTANT!! If running simulatenously, the flag convertdata should be set to FALSE, and the output files of the data conversion step should already be present.

# The script automatically avoids double calculations (i.e, i vs j, and j vs i) by only performing calculations if ( i < j && (i + j)%%2!=0) or ( i > j && (i + j)%%2==0)
# For instance, for individuals 1 and 5, the sum is even (6), and hence the script performs calculations for i=5 and j=1, but not for i=1 and j=5.
# In contrast, for individuals 1 and 6, the sum is odd (7), and hence the script performs calculations for i=1 and j=6, but not for i=6 and j=1.     


##############################################################
# CONTROL PANEL
start1=$1       		# Specify this number on the command line.
end1=$2         		# Specify this number on the command line.

start2=1
end2=270        		# Set this number to the total of individuals in the input vcf-file.

convertdata=TRUE		# This step will create the intermediate files 'myinput.samples.txt', 'myinput.samples.txt' and 'myinput.pos.txt'
replacebar=FALSE		# If data is phased (meaning alleles are separated by a vertical bar rather than forward slash), set this is TRUE when running convertdata step

run_loop=FALSE
haploiddata=FALSE		# Set to TRUE in case input vcf-file contains haploid data (e.g. MT-DNA or Y-chromosome)	
haplodiploiddata=FALSE	# Set to TRUE in case input vcf-file contains haplodiploid data (e.g. X-chromosome)

MYVCF=allsites.thinned.100.vcf.gz									# specify here name to input vcf-file (which contains an unbiased representation of monomorphic and polymorphic sites)
BCFTOOLS=/home/mdejong/software/bcftools/bcftools-1.20/bcftools		# specify here path to bcftools executable
#BCFTOOLS=/opt/software/bcftools-1.20/bin/bcftools
##############################################################


if [[ "$convertdata" == TRUE ]]
	then
	echo "Converting vcf to geno..."
	$BCFTOOLS query --list-samples $MYVCF > myinput.samples.txt
	$BCFTOOLS query -f '[\t%GT]\n' $MYVCF | sed 's/^[ \t]*//' > myinput.geno.txt
	zgrep -v '#' $MYVCF | head -1 | cut -f1-2 > myinput.pos.txt
	#
	if [[ "$replacebar" == TRUE ]]
		then
		sed -i 's:|:/:g' myinput.geno.txt
	fi
	else
	echo "Assuming input files 'myinput.geno.txt' and 'myinput.samples.txt' are already present."
fi

if [[ "$run_loop" == TRUE ]]
	then
	mychrom=$(cut -f1 myinput.pos.txt)
	mypos=$(cut -f2 myinput.pos.txt)
	nsites=$(wc -l myinput.geno.txt | cut -f1 -d ' ')
	
	echo "Starting pairwise comparisons..."
	echo "ind1 ind2 name1 name2 nsites nmiss pmiss d ploidy chrom startbp" > vcfdist.${start1}_${end1}.txt
	for (( i = $start1; i <= $end1; i++ ))
	do
		echo "$i"
		ind1=$(awk -v myline="$i" 'NR==myline' myinput.samples.txt)
		cut -f${i} myinput.geno.txt > myind1.${i}.txt
		
		for (( j = $start2 ; j <= $end2; j++ ))
			do
			
			# run analysis for this pair or skip?
			docalc=TRUE
			k=$(( $i + $j ))
            if [[ $i -lt $j && $(($k % 2)) = 0 ]]; then docalc=FALSE; fi
            if [[ $i -gt $j && $(($k % 2)) != 0 ]]; then docalc=FALSE; fi
			
			if [[ "$docalc" == TRUE ]]
				then
				# select data:
				# echo "$i $j"
				ind2=$(awk -v myline="$j" 'NR==myline' myinput.samples.txt)
				cut -f${j} myinput.geno.txt > myind2.${i}_${j}.txt
				
				paste myind1.${i}.txt myind2.${i}_${j}.txt | sed 's|\/|\t|g' > mypair.${i}_${j}.txt
				grep -v '\.' mypair.${i}_${j}.txt > mypair.${i}_${j}.nomissing.txt
				nmiss=$(grep '\.' mypair.${i}_${j}.txt | wc -l)
				pmiss=$(awk "BEGIN {print $nmiss/$nsites}")
				
				if [[ "$haploiddata" == TRUE ]]
					then
					# haploid data:
					awk -v OFS="\t" '{ if($1==$2) $3=0; else $3=1; print $3; }' mypair.${i}_${j}.nomissing.txt > mypair.${i}_${j}.sum.txt
					dist=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' mypair.${i}_${j}.sum.txt)
					ploidy="HA_HA"
					else
					if [[ "$haplodiploiddata" == TRUE ]]
						then
						# haplodiploid data:
						indbool1=$(head -1 myind1.${i}.txt | grep '/' | wc -l)
						indbool2=$(head -1 myind2.${i}_${j}.txt | grep '/' | wc -l)
						
						if [[ "$indbool1" == 0 && "$indbool2" == 0 ]]
							then
							echo "Individual $i and individual $j are both haploid."
							awk -v OFS="\t" '{ if($1==$2) $3=0; else $3=1; print $3; }' mypair.${i}_${j}.nomissing.txt > mypair.${i}_${j}.sum.txt
							dist=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' mypair.${i}_${j}.sum.txt)
							ploidy="HA_HA"
						fi	
						if [[ "$indbool1" == 1 && "$indbool2" == 0  ]]
							then
							echo "Individual $i is diploid and $j is haploid."
							awk -v OFS="\t" '{ if($1==$3) $4=0; else $4=1; if($2==$3) $5=0; else $5=1; $6=($4+$5)/2; print $6; }' mypair.${i}_${j}.nomissing.txt > mypair.${i}_${j}.sum.txt
							dist=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' mypair.${i}_${j}.sum.txt)	
							ploidy="DI_HA"
						fi	
						if [[ "$indbool1" == 0 && "$indbool2" == 1  ]]
							then
							echo "Individual $i is haploid and $j is diploid."
							awk -v OFS="\t" '{ if($1==$2) $4=0; else $4=1; if($1==$3) $5=0; else $5=1; $6=($4+$5)/2; print $6; }' mypair.${i}_${j}.nomissing.txt > mypair.${i}_${j}.sum.txt
							dist=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' mypair.${i}_${j}.sum.txt)
							ploidy="HA_DI"
						fi	
						if [[ "$indbool1" == 1 && "$indbool2" == 1  ]]
							then
							echo "Individual $i and $j are both diploid."
							awk -v OFS="\t" '{ if($1==$3) $5=0; else $5=1; if($1==$4) $6=0; else $6=1; if($2==$3) $7=0; else $7=1; if($2==$4) $8=0; else $8=1; $9=($5+$6+$7+$8)/4; print $9; }' mypair.${i}_${j}.nomissing.txt > mypair.${i}_${j}.sum.txt
							dist=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' mypair.${i}_${j}.sum.txt)	
							ploidy="DI_DI"
						fi	
						else
						# diploid data:
						awk -v OFS="\t" '{ if($1==$3) $5=0; else $5=1; if($1==$4) $6=0; else $6=1; if($2==$3) $7=0; else $7=1; if($2==$4) $8=0; else $8=1; $9=($5+$6+$7+$8)/4; print $9; }' mypair.${i}_${j}.nomissing.txt > mypair.${i}_${j}.sum.txt
						dist=$(awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' mypair.${i}_${j}.sum.txt)
						ploidy="DI_DI"	
					fi
				fi
			echo "$i $j $ind1 $ind2 $nsites $nmiss $pmiss $dist $ploidy $mychrom $mypos" >> vcfdist.${start1}_${end1}.txt
			rm myind2.${i}_${j}.txt mypair.${i}_${j}.txt mypair.${i}_${j}.nomissing.txt mypair.${i}_${j}.sum.txt
			fi
			done
		rm myind1.${i}.txt
	done
	sed -i 's/ /\t/g' vcfdist.${start1}_${end1}.txt
	#
	echo "Correcting within-individual comparisons..."
	awk -v OFS="\t" '{ if($3==$4 && $9=="DI_DI") $8=$8+$8; else $8=$8; print $0; }' vcfdist.${start1}_${end1}.txt > vcfdist.${start1}_${end1}.tmp.txt
	mv vcfdist.${start1}_${end1}.tmp.txt vcfdist.${start1}_${end1}.txt
	#
	echo "Finished analyses :-)."
	echo "Pairwise sequence dissimilarity estimates stored in file 'vcfdist.txt'."
	echo "These estimates can be used for distance-based phylogenetic analyses in SambaR." 
	else
	echo "Not running analyses because the flag 'run_loop' is set to FALSE."
fi
