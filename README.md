# Calculate pairwise distances

This script allows to calculate sequence dissimilarity estimates between all pairs of individuals within a VCF-file.
The input-file can be a VCF-file containing variable sites only. 
However, in order to scale the obtained estimates, it is important to know the total number of sites from which the variable sites have been extracted. 
The genome-wide distance can be obtained using the formula: d = d_snps*n _snps/n_sites. 

For instance: say that you have a gVCF-file with a total number of sites of 2Gb (monomorphic and polymorphic, but without indels). 
Say furthermore that after selecting variable sites, you obtained a VCF-file containing 100Mb SNPs (biallelic AND multiallelic!).
The smaller VCF-file with SNP data is to be used as input for the VCF_calcdist script.
Say that for two individuals, A and B, the obtained SNP distance (d_snps) is 0.2. This will be the value in the output file of the VCF_calcdist.sh.
The genome-wide distance is: d = 0.2x(10^7)/(2x10^9) = 0.001 = 0.1%.


# Usage
Execute by typing, for example:
dos2unix VCF_calcdist.sh
chmod +x VCF_calcdist.sh
./VCF_calcdist.sh 1 10
This will calculate all pairwise scores for individuals 1 to 10.

Run the script in two steps: 
1. convert the data (convertdata=TRUE,run_loop=FALSE), which will generate the files 'myinput.samples.txt', 'myinput.samples.txt' and 'myinput.pos.txt' 
2. afterwards run the actual analysis (convertdata=FALSE,run_loop=TRUE)

For the second step, if you forget to add the two numbers after VCF_calcdist.sh, you will receive the error:
line 76: ((: i = : syntax error: operand expected (error token is "= ")

# Parallel usage
To speed up the calculations, run simulatenously for subsets as such:
./VCF_calcdist.sh 1 10 &		
./VCF_calcdist.sh 11 20 &
etc.
IMPORTANT!! If running simulatenously, the flag convertdata should be set to FALSE, and the output files of the data conversion step should already be present.

The script automatically avoids double calculations (i.e, i vs j, and j vs i) by only performing calculations if ( i < j && (i + j)%%2!=0) or ( i > j && (i + j)%%2==0)
For instance, for individuals 1 and 5, the sum is even (6), and hence the script performs calculations for i=5 and j=1, but not for i=1 and j=5.
In contrast, for individuals 1 and 6, the sum is odd (7), and hence the script performs calculations for i=1 and j=6, but not for i=6 and j=1.     

# Use distance estimates for population-genetic analyses

The obtained distance-estimates can be used for distance-based population-genetic analyses, such as tree reconstruction, pcoa-analyses and heterozygosity and Hudson Fst estimation.
To do, make first a dummy dataset in SambaR:

mydf	          <- read.table("popfile.txt",header=TRUE)			# specify here name of tab-separated file, which should contain two columns: name and pop; names should correspond to names in allvcfdist.txt file.

mymat			      <- matrix(sample(c(0,1),nrow(mydf)*100,replace=TRUE),nrow=nrow(mydf),ncol=100)

rownames(mymat)	<- mydf$name

colnames(mymat)	<- paste("snp",c(1:100),sep="_")	

mygl			      <- as.genlight(mymat)

genlight2sambar(genlight_object="mygl",do_confirm=TRUE)

mydf$pop		    <- NULL

mydf$popcol		  <- NULL

inds			      <- merge(inds,mydf,by="name")

filterdata(min_spacing=0,min_mac=1,dohefilter=FALSE)

Next, run the analyses:

add2inds2(myfile="allvcfdist.txt",miss_filter=inds2$pmiss<=0.8,ntotalsites=2000000000)

calcdxy()
