./laser -g ../example/HGDP_238_chr22.geno -c ../example/HGDP_238_chr22.RefPC.coord -s ../example/HapMap_6_chr22.seq -o test_laser

./laser -g ../example/HGDP_238_chr22.geno -k 20 -pca 1 -o test_pca

./trace -g ../example/HGDP_238_chr22.geno -c ../example/HGDP_238_chr22.RefPC.coord -s ../example/HGDP_700_chr22.geno -x 1 -y 700 -o test_trace

