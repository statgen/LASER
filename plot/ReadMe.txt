This folder include example R codes to generate figures for results based on the HGDP reference panel. 

This script "hgdp_cols_symb.r" provides the list of colors and symbols for 53 HGDP populations. 
This set of colors and symbols has been consistently used in a series of publications based on the HGDP data:

	- The color scheme for the HGDP samples was first used in Rosenberg et al Science (2002) 298: 2381-2385.
	- The list of symbols was first used for 29 HGDP populations in Jakobsson et al Nature (2008) 451:998-1003.
	- The list of symbols was augmented to include 53 HGDP populations in Wang et al (2010) Stat Appl Genet Mol Biol 9:13.
	- The color scheme and symbols have been consistently used in many papers, especially those from Noah Rosenberg's lab.
	- Results in the LASER paper (Wang et al. 2014, Nature Genetics) were based on the same colors and symbols.

To generate your own figure:

	1. Modify the first two lines of the script "plotHGDP.r" to specify your coordinate files.
	2. Run source("plotHGDP.r") in the R environment. (See www.r-project.org for information about R.)
	3. A figure file named "Results_on_HGDP.pdf" will be generated in the current directory.

========
Chaolong Wang (chaolong@umich.edu)
May 15, 2014

