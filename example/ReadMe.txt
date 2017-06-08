=== Example data in this folder ===

The files in this example folder include genotype data from the Human Genome Diversity Panel (Li et al. 2008 Science). 
In total, genotypes at 9,608 SNP markers on chromosome 22 for 938 individuals from 53 worldwide populations are included. 
We also include sequencing data for 6 HapMap samples, piled up to the 9,608 SNP markers in the HGDP example data.


The HGDP_700_chr22.geno file contains genotype data for 700 individuals randomly selected from the HGDP sample.

The HGDP_700_chr22.site file provides detailed information of the 9608 SNP markers in the HGDP_700_chr22.geno file.

The HGDP_238_chr22.geno file contains genotype data for the remaining 238 HGDP individuals.

The HGDP_238_chr22.site file provides detailed information of the 9608 SNP markers in the HGDP_238_chr22.geno file.

The HGDP_238_chr22.RefPC file contains coordinates of the first 8 PCs for the 238 HGDP samples, obtained by performing PCA on genotypes at 9608 loci.

The HapMap_6_chr22.seq file contains sequencing reads for 6 HapMap samples, including 1 CEU trio and 1 YRI trio.


=== Run the program ===

To run the TRACE program on the example data, type:  

./trace -p ./example/trace.conf

which will run the TRACE program with parameter values taken from the parameter file "./example/trace.conf".

--

To run the LASER program on the example data, type:

./laser -p ./example/laser.conf

which will run the LASER program with parameter values taken from the parameter file "./example/laser.conf".


=== Questions ===

If you have questions, please read the User Manuals first.

If the manual does not answer your question, contact Chaolong Wang at chaolong@umich.edu.
