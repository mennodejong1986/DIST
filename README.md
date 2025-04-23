# DIST: Distance-based Inference of Species Trees

DIST is a method to infer species trees from distance matrices containing genome-wide sequence dissimilarity estimates.
The workflow of DIST contains three steps:

1. DIST measures genetic distance for each pair of individuals in terms of mean pairwise sequence dissimilarity, E(p), preferably from an unbiased selection of monomorphic and polymorphic sites. These distances are optionally depicted in a concatenation tree, traditionally known as 'tree of individuals'.

2. DIST uses basic equations from coalescent theory to analytically infer the species tree which best predicts, in terms of population split times measured in coalescent units, these observed genetic distances.

3. During the optional third step, gene flow edges are added to resolve cases of non-additivity. 

This github-repository contains a Unix-script to perform the first step: calculate E(p)-estimates for all pairs of individuals within an input VCF-file. 
The input-file can be a VCF-file containing variable sites only. 
However, in order to scale the obtained estimates, it is important to know the total number of sites from which the variable sites have been extracted. 
The genome-wide distance can be obtained using the formula: d = d_snps*n _snps/n_sites. 

For instance: say that you have a gVCF-file with a total number of sites of 2Gb (monomorphic and polymorphic, but without indels). 
Say furthermore that after selecting variable sites, you obtained a VCF-file containing 100Mb SNPs (biallelic AND multiallelic!).
The smaller VCF-file with SNP data is to be used as input for the VCF_calcdist script.
Say that for two individuals, A and B, the obtained SNP distance (d_snps) is 0.2. This will be the value in the main output file, called 'allvcfdist.txt'.
The genome-wide distance is: d = 0.2x(10^7)/(2x10^9) = 0.001 = 0.1%.

# Usage
Execute by typing, for example:

*dos2unix VCF_calcdist.sh*

*chmod +x VCF_calcdist.sh*

*./VCF_calcdist.sh 1 10*

This will calculate all pairwise scores for individuals 1 to 10.

Run the script in two steps: 
1. convert the data (convertdata=TRUE,run_loop=FALSE), which will generate the files 'myinput.samples.txt', 'myinput.samples.txt' and 'myinput.pos.txt' 
2. afterwards run the actual analysis (convertdata=FALSE,run_loop=TRUE)

For the second step, if you forget to add the two numbers after VCF_calcdist.sh, you will receive the error:
line 76: ((: i = : syntax error: operand expected (error token is "= ")

# Parallel usage
To speed up the calculations, run simulatenously for subsets as such:

*./VCF_calcdist.sh 1 10 &*		

*./VCF_calcdist.sh 11 20 &*

etc.

IMPORTANT!! If running simulatenously, the flag convertdata should be set to FALSE, and the output files of the data conversion step should already be present.

The script automatically avoids double calculations (i.e, i vs j, and j vs i) by only performing calculations if ( i < j && (i + j)%%2!=0) or ( i > j && (i + j)%%2==0)
For instance, for individuals 1 and 5, the sum is even (6), and hence the script performs calculations for i=5 and j=1, but not for i=1 and j=5.
In contrast, for individuals 1 and 6, the sum is odd (7), and hence the script performs calculations for i=1 and j=6, but not for i=6 and j=1.     

As a rough estimation of computation time: a gVCF-file of 100 individuals and 50Mb variable sites can be processed within 8 hours, if the workload is divided over 25 parallel runs (i.e., 4 individuals per run).     
The combined output will be stored in the output file 'allvcfdist.txt'.

# Species tree inference in SambaR

The file 'allvcfdist.txt' can be as input for SambaR for species tree inference.
To do so, first create a dummy genlight object in R:

*source("https://github.com/mennodejong1986/SambaR/raw/master/SAMBAR_v1.xx.txt")* # replace xx with the correct version

*getpackages()*

*gldummy(popfile="popfile.txt")*   # popfile.txt should be a tab-separated file, which should contain two columns: name and pop; all names should correspond to the names in the file 'allvcfdist.txt'.


Convert this genlight object into SambaR objects:

*genlight2sambar(genlight_object="mygl",do_confirm=TRUE,popvector=as.character(inddf$pop),pop_order=NULL,colourvector=NULL)*

*filterdata(min_spacing=0,min_mac=1,dohefilter=FALSE)*


Next, add to SambaR's inds2 dataframe the information from the file 'allvcfdist.txt':

*add2inds2(myfile="allvcfdist.txt",miss_filter=inds2$pmiss<=0.8,ntotalsites=2000000000)*    # replace ntotalsites with the the correct number. See the explanation above.


Next, run the analyses:

*runDIST()*
