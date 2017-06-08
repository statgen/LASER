# Create bed file from site file:

awk '{if(NR>1){print $1, $2-1, $2, $3;}}' ../example/HGDP_238_chr22.site > HGDP_238_chr22.bed

# Generate .pileup file:

samtools mpileup -q 30 -Q 20 -f PATH/hs37d5.fa -l HGDP_238_chr22.bed exampleBAM/NA12878.chrom22.recal.bam > NA12878.chrom22.pileup

# Convert .pileup file to LASER .seq format:

python pileup2seq.py  -f PATH/hs37d5.fa -m ../example/HGDP_238_chr22.site -o test NA12878.chrom22.pileup 

## If running pileup2seq.py without the -f option, the program will proceed without checking the reference alleles. 
## It is highly recommended to use the -f option so that possible strand flipping issue in the ancestry reference panel can be detected.
