/*  
	LASER: Locating Ancestry from SEquence Reads
    Copyright (C) 2013-2016  Chaolong Wang

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <map>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics_double.h>
#include <cmath>
#include <algorithm>
#include <time.h>
#include <stdlib.h>
#include <omp.h>
#include <cblas.h>
#define  ARMA_DONT_USE_WRAPPER    
#include "armadillo"
using namespace arma;
using namespace std;

const string ARG_PARAM_FILE = "-p";
const string ARG_GENO_FILE = "-g";
const string ARG_COORD_FILE = "-c";
const string ARG_SEQ_FILE = "-s";
const string ARG_OUT_PREFIX = "-o";
const string ARG_MIN_LOCI = "-l";
const string ARG_MAX_LOCI = "-L";
const string ARG_DIM = "-k";
const string ARG_DIM_HIGH = "-K";
const string ARG_SEQ_ERR = "-e";

const string ARG_FIRST_IND = "-x";
const string ARG_LAST_IND = "-y";
const string ARG_ALPHA = "-a";
const string ARG_THRESHOLD = "-t";
const string ARG_REPS = "-r";
const string ARG_OUTPUT_REPS = "-R";
const string ARG_CHECK_COVERAGE = "-cov";
const string ARG_CHECK_FORMAT = "-fmt";
const string ARG_PCA_MODE = "-pca";
const string ARG_EXCLUDE_LIST = "-ex";
const string ARG_TRIM_PROP = "-M";
const string ARG_MIN_COVERAGE = "-minc";
const string ARG_MAX_COVERAGE = "-maxc";
const string ARG_RANDOM_SEED= "-seed";
const string ARG_PROCRUSTES_SCALE = "-rho";
const string ARG_REF_SIZE = "-N";
const string ARG_KNN_ZSCORE= "-knn";
const string ARG_NUM_THREADS = "-nt";

const string default_str = "---this-is-a-default-string---";
const int default_int = -999999998;
const double default_double = -9.99999999;
char tmpstr[2000];

string LOG_FILE;
string TMP_LOG_FILE = "Temporary.";

string PARAM_FILE = default_str;    // parameter file name
string GENO_FILE = default_str;     // Reference genotype datafile name
string COORD_FILE = default_str;    // Reference coordinates datafile name
string SEQ_FILE = default_str;      // Sequence datafile name
string OUT_PREFIX = default_str;    // Prefix for output files

double SEQ_ERR = default_double;    // Sequencing error rate per read;
									// By default set to -1, which will use the estimated value from .seq file

int DIM = default_int;              // Number of PCs to match;
int DIM_HIGH = default_int;         // Number of PCs in from sample-specific PCA;
int MIN_LOCI = default_int;         // Minimum loci that have at least one read;
int MAX_LOCI = default_int;         // Maximum loci for each study sample to include in the analysis;

int REF_SIZE = default_int;         // Number of individuals randomly selected from GENO_FILE as anchor points;

int FIRST_IND = default_int;	    // First individual in the list sample to be tested; 
int LAST_IND = default_int;	        // Last individual in the list of sample to be tested;
double ALPHA = default_double;      // Significance level of the TW statistic when setting DIM_HIGH
double THRESHOLD = default_double;  // Convergence criterion for the projection Procrustes analysis
int REPS = default_int;             // Number of replicates to run for each sample;
int OUTPUT_REPS = default_int;      // 0: Only output the mean and sd; 
                                    // 1: Output results from all replicates;
int CHECK_COVERAGE = default_int;   // 0: do not check coverage, proceed to major computation;
                                    // 1: check coverage first, proceed to major computation;
                                    // 2: check coverage only.
int CHECK_FORMAT = default_int;     // 0: Do not check format of the input data files, and proceed to major computation;
                                    // 1: Check format of all files and stop;  10: Proceed after checking all files
                                    // 2: Check format of GENO_FILE and stop;  20: Proceed after checking GENO_FILE
                                    // 3: Check format of SEQ_FILE and stop;   30: Proceed after checking SEQ_FILE
                                    // 4: Check format of COORD_FILE and stop; 40: Proceed after checking COORD_FILE
int PCA_MODE = default_int;         // 0: LASER;
                                    // 1: PCA on genotypes;
									// 2: PCA on genotypes (slow but memory efficient);
									// 3: PCA based on SVD and will output SNP loadings;
double TRIM_PROP = default_double;  // Proportion of random loci to be removed from analyses (same set across all samples)
string EXCLUDE_LIST = default_str;  // File name of a list of SNPs to exclude
double MIN_COVERAGE = default_double;  // Excluding sites with mean coverage smaller than MIN_COVERAGE; 
double MAX_COVERAGE = default_double;  // Excluding sites with mean coverage greater than MAX_COVERAGE;
int PROCRUSTES_SCALE = default_int;  // 0: Fit the scaling parameter to maximize similarity   
									 // 1: Fix the scaling to match the variance between X and Y
int RANDOM_SEED = default_int;        // Random seed used in the program
int KNN_ZSCORE = default_int;       // Number of nearest neigbors used to calculate the Z score for each study individual.
int NUM_THREADS = default_int;        // Number of CPU cores for multi-threading parallel analysis 
									 
// The following parameters will be determined from the input data files					 
int REF_INDS = default_int;         // Number of reference individuals
int SEQ_INDS = default_int;         // Number of sequence samples

int LOCI_G = default_int;      // Number of loci in GENO_FILE;
int LOCI_S = default_int;      // Number of loci in SEQ_FILE;
int LOCI = default_int;        // Number of shared loci

int NUM_PCS = default_int;          // Number of PCs in the COORD_FILE;
int GENO_NON_DATA_ROWS = 0;    // Number of non-data rows in the GENO_FILE;
int GENO_NON_DATA_COLS = 2;    // Number of non-data columns in the GENO_FILE;
int COORD_NON_DATA_ROWS = 1;   // Number of non-data rows in the COORD_FILE;
int COORD_NON_DATA_COLS = 2;   // Number of non-data columns in the COORD_FILE;
int SEQ_NON_DATA_ROWS = 0;     // Number of non-data rows in the SEQ_FILE;
int SEQ_NON_DATA_COLS = 2;     // Number of non-data columns in the SEQ_FILE;

bool AUTO_MODE = false; // If the program will determine DIM_HIGH automatically;
int MAX_ITER = 10000;    // Maximum iterations for the projection Procrustes analysis
double TW = default_double;    // Threshold to determine significant Tracy-Widom statistic

string GENO_SITE_FILE = default_str;      // Sitefile of the reference data
string SEQ_SITE_FILE = default_str;      // Sitefile of the sequence data

//=======================================================================================
bool is_int(string str);
bool is_numeric(string str);
bool parse_cmd_line(int argc, char* argv[], map<string,string> &args, map<string,int> &argi, map<string,double> &argd);
int read_paramfile(string filename);
int create_paramfile(string filename);
int check_parameters();
void print_configuration();

int pca_geno(Mat<char> &G, int nPCs, mat &PC, rowvec &PCvar);
int pca(fmat G, int nPCs, mat &PC, rowvec &PCvar, mat &M);
int pca_svd(fmat G, int nPCs, mat &PC, rowvec &PCvar, fmat &Gm, fmat &Gsd, fmat &W);
int normalize(fmat &G, fmat &Gm, fmat &Gsd);
int procrustes(mat &X, mat &Y, mat &Xnew, double &t, double &rho, mat &A, rowvec &b, int ps);
double pprocrustes(mat &X, mat &Y, mat &Xnew, double &t, double &rho, mat &A, rowvec &b, int iter, double eps, int ps);
int simuseq(Mat<char> &G, urowvec &C, uvec &Loc, double e, fmat &S, gsl_rng *rng);
int simuseq2(Mat<char> &G, urowvec &C, uvec &Loc, frowvec &Q, fmat &S, gsl_rng *rng);

int check_coverage(int output, int first_ind, int last_ind, uvec cmnS, urowvec &ExLoci, int &Ls, int &Lg);
int check_format_geno(string filename, int inds, int loci);
int check_format_seq(string filename, int inds, int loci);
int check_format_coord(string filename, int inds, int npcs);
bool get_table_dim(int &nrow, int &ncol, string filename, char separator);

ofstream foutLog;
//=========================================================================================================
int main(int argc, char* argv[]){
	int i=0;
	int j=0;
	int k=0;
	int tmp=0;
	
	string str;
	ifstream fin;
	string outfile;
	ofstream fout;
	string outfile2;
	ofstream fout2;
	string outfile3;
	ofstream fout3;

	wall_clock timer;
	double runningtime;
	timer.tic();   // Program starting time
	
	time_t rawtime;
  	struct tm * timeinfo;
 	time ( &rawtime );
  	timeinfo = localtime ( &rawtime );

	stringstream ss;
	ss << getpid();
	string strpid = ss.str();
	TMP_LOG_FILE.append(strpid);
	TMP_LOG_FILE.append(".log");
	
	foutLog.open(TMP_LOG_FILE.c_str());
	if(foutLog.fail()){
		cerr << "Error: cannot create a temporary log file." << endl;
		return 1;
	}

	cout << endl;
	cout << "====================================================================" <<endl;
	cout << "====        LASER: Locating Ancestry from SEquence Reads        ====" <<endl; 
	cout << "====          Version 2.04, Last updated on Jan/11/2017         ====" <<endl;	
	cout << "====          (C) 2013-2017 Chaolong Wang, GNU GPL v3.0         ====" <<endl;
	cout << "====================================================================" <<endl;
  	cout << "Started at: " << asctime (timeinfo) << endl;

	foutLog << "====================================================================" <<endl;
	foutLog << "====        LASER: Locating Ancestry from SEquence Reads        ====" <<endl; 
	foutLog << "====          Version 2.04, Last updated on Jan/11/2017         ====" <<endl;
	foutLog << "====          (C) 2013-2017 Chaolong Wang, GNU GPL v3.0         ====" <<endl;		
	foutLog << "====================================================================" <<endl;
  	foutLog << "Started at: " << asctime (timeinfo) << endl;
		
	// ################ Read in command line ##########################
	map<string,string> args;
	map<string,int> argi;
	map<string,double> argd;	
	bool cmd_flag = parse_cmd_line(argc, argv, args, argi, argd);
	if(args[ARG_PARAM_FILE].compare(default_str)!=0){PARAM_FILE = args[ARG_PARAM_FILE];}
	if(args[ARG_GENO_FILE].compare(default_str)!=0){GENO_FILE = args[ARG_GENO_FILE];}
	if(args[ARG_COORD_FILE].compare(default_str)!=0){COORD_FILE = args[ARG_COORD_FILE];}
	if(args[ARG_SEQ_FILE].compare(default_str)!=0){SEQ_FILE = args[ARG_SEQ_FILE];}
	if(args[ARG_OUT_PREFIX].compare(default_str)!=0){OUT_PREFIX = args[ARG_OUT_PREFIX];}
	if(argi[ARG_MIN_LOCI]!=default_int){MIN_LOCI = argi[ARG_MIN_LOCI];}
	if(argi[ARG_MAX_LOCI]!=default_int){MAX_LOCI = argi[ARG_MAX_LOCI];}
	if(argi[ARG_DIM]!=default_int){DIM = argi[ARG_DIM];}
	if(argi[ARG_DIM_HIGH]!=default_int){DIM_HIGH = argi[ARG_DIM_HIGH];}	
	if(argd[ARG_SEQ_ERR]!=default_double){SEQ_ERR = argd[ARG_SEQ_ERR];}
	if(argi[ARG_REF_SIZE]!=default_int){REF_SIZE = argi[ARG_REF_SIZE];}
	if(argi[ARG_FIRST_IND]!=default_int){FIRST_IND = argi[ARG_FIRST_IND];}
	if(argi[ARG_LAST_IND]!=default_int){LAST_IND = argi[ARG_LAST_IND];}
	if(argd[ARG_ALPHA]!=default_double){ALPHA = argd[ARG_ALPHA];}
	if(argd[ARG_THRESHOLD]!=default_double){THRESHOLD = argd[ARG_THRESHOLD];}		
	if(argi[ARG_REPS]!=default_int){REPS = argi[ARG_REPS];}
	if(argi[ARG_OUTPUT_REPS]!=default_int){OUTPUT_REPS = argi[ARG_OUTPUT_REPS];}
	if(argi[ARG_CHECK_COVERAGE]!=default_int){CHECK_COVERAGE = argi[ARG_CHECK_COVERAGE];}
	if(argi[ARG_CHECK_FORMAT]!=default_int){CHECK_FORMAT = argi[ARG_CHECK_FORMAT];}
	if(argi[ARG_PCA_MODE]!=default_int){PCA_MODE = argi[ARG_PCA_MODE];}	
	if(argd[ARG_MIN_COVERAGE]!=default_double){MIN_COVERAGE = argd[ARG_MIN_COVERAGE];}		
	if(argd[ARG_MAX_COVERAGE]!=default_double){MAX_COVERAGE = argd[ARG_MAX_COVERAGE];}
	if(args[ARG_EXCLUDE_LIST]!=default_str){EXCLUDE_LIST = args[ARG_EXCLUDE_LIST];}		
	if(argd[ARG_TRIM_PROP]!=default_double){TRIM_PROP = argd[ARG_TRIM_PROP];}
	if(argi[ARG_PROCRUSTES_SCALE]!=default_int){PROCRUSTES_SCALE = argi[ARG_PROCRUSTES_SCALE];}	
	if(argi[ARG_RANDOM_SEED]!=default_int){RANDOM_SEED = argi[ARG_RANDOM_SEED];}
	if(argi[ARG_KNN_ZSCORE]!=default_int){KNN_ZSCORE = argi[ARG_KNN_ZSCORE];}
	if(argi[ARG_NUM_THREADS]!=default_int){NUM_THREADS = argi[ARG_NUM_THREADS];}
	//##################     Read in and check parameter values #######################
	if(PARAM_FILE.compare(default_str)==0){ PARAM_FILE = "laser.conf"; }
	int flag = read_paramfile(PARAM_FILE);
	if(flag==0 || cmd_flag==0){
		foutLog.close();
		if(OUT_PREFIX.compare(default_str)==0){
			LOG_FILE = "laser.log";
		}else{
			LOG_FILE = OUT_PREFIX;
			LOG_FILE.append(".log");
		}
		sprintf(tmpstr, "%s%s%s%s", "mv ", TMP_LOG_FILE.c_str()," ", LOG_FILE.c_str());
		int sys_msg = system(tmpstr);
		return 1;
	}
	//################  Set default values to some parameters #######################
	if(SEQ_ERR==default_double){ SEQ_ERR = -1; }
	if(MIN_LOCI==default_int){ MIN_LOCI = 100; }
	if(MAX_LOCI==default_int){ MAX_LOCI = 1000000; }
	if(DIM==default_int){ DIM = 2; }
	if(REPS == default_int){ REPS = 1; }
	if(OUTPUT_REPS == default_int){ OUTPUT_REPS = 0; }
	if(CHECK_COVERAGE==default_int){ CHECK_COVERAGE = 0; }
	if(CHECK_FORMAT==default_int){ CHECK_FORMAT = 10; }
	if(PCA_MODE==default_int){ PCA_MODE = 0; }
	if(DIM_HIGH == default_int){ DIM_HIGH = 20; }  // 0 means DIM_HIGH will be automatically set by the program
	if(ALPHA == default_double){ ALPHA = 0.1; }	
	if(THRESHOLD == default_double) { THRESHOLD = 0.000001; }
	if(MIN_COVERAGE==default_double){ MIN_COVERAGE = 0; }	
	if(MAX_COVERAGE==default_double){ MAX_COVERAGE = -1; }
	if(TRIM_PROP==default_double){ TRIM_PROP = 0; }
	if(PROCRUSTES_SCALE==default_int){ PROCRUSTES_SCALE = 0; }
	if(RANDOM_SEED==default_int){ RANDOM_SEED = 0; }
	if(KNN_ZSCORE==default_int){ KNN_ZSCORE = 10; }
	if(NUM_THREADS==default_int){ NUM_THREADS = 8; }
	//###############################################################################
	if(OUT_PREFIX.compare(default_str)==0){ OUT_PREFIX = "laser"; }
	foutLog.close();
	LOG_FILE = OUT_PREFIX;
	LOG_FILE.append(".log");
	sprintf(tmpstr, "%s%s%s%s", "mv ", TMP_LOG_FILE.c_str()," ", LOG_FILE.c_str());
	int sys_msg = system(tmpstr);
	foutLog.open(LOG_FILE.c_str(), ios::app);
	if(foutLog.fail()){
		cerr << "Error: cannot create the log file." << endl;	
		return 1;
	}
	//############## Get values for REF_INDS, LOCI, SEQ_INDS, NUM_PCS ################
	int nrow = 0;
	int ncol = 0;
	flag = 1;
	map<string,int> idxS;
	map<string,int> idxG;
	map<string,string> alleleS;
	map<string,string> alleleG;	
	vector<string> cmnsnp;
	uvec cmnS;
	uvec cmnG;
	map<string,int> exSNP;	

	int Lex = 0;
	int LOCI_trim = 0;
	int LOCI_ex = 0;
	int unmatchSite = 0;
	
	gsl_rng *rng;
	rng = gsl_rng_alloc(gsl_rng_taus);
	// long seed = time(NULL)*getpid();
	gsl_rng_set(rng, RANDOM_SEED);
	
	if(SEQ_FILE.compare(default_str) != 0 && flag == 1){
		if(!get_table_dim(nrow, ncol, SEQ_FILE, '\t')){
			cerr << "Error: cannot open the SEQ_FILE '" << SEQ_FILE << "'." << endl;    
			foutLog << "Error: cannot open the SEQ_FILE '" << SEQ_FILE << "'." << endl; 
			foutLog.close();	
			return 1;
		}	
		SEQ_INDS = nrow - SEQ_NON_DATA_ROWS;
		int tmpLOCI = ncol - SEQ_NON_DATA_COLS;
		cout << SEQ_INDS << " individuals in the SEQ_FILE." << endl;
		foutLog << SEQ_INDS << " individuals in the SEQ_FILE." << endl;
		if(SEQ_INDS < 0){
			cerr << "Error: Invalid number of rows in the SEQ_FILE '" << SEQ_FILE << "'." << endl;
			foutLog << "Error: Invalid number of rows in the SEQ_FILE '" << SEQ_FILE << "'." << endl;
			flag = 0;
		}
		if(tmpLOCI < 0){
			cerr << "Error: Invalid number of columns in the SEQ_FILE '" << SEQ_FILE << "'." << endl;
			foutLog << "Error: Invalid number of columns in the SEQ_FILE '" << SEQ_FILE << "'." << endl;
			flag = 0;
		}		
		SEQ_SITE_FILE = SEQ_FILE;
		SEQ_SITE_FILE.replace(SEQ_SITE_FILE.length()-4, 5, ".site");
		fin.open(SEQ_SITE_FILE.c_str());
		if(fin.fail()){
			cerr << "Error: cannot open the file '" << SEQ_SITE_FILE << "'." << endl;    
			foutLog << "Error: cannot open the file '" << SEQ_SITE_FILE << "'." << endl;   
			foutLog.close();
			return 1;
		}else{
			getline(fin, str);
			LOCI_S = 0;
			while(!fin.eof()){
				getline(fin, str);
				if(str.length()>0 && str!=" "){
					j=0;
					int tabpos[5];
					for(i=0; i<str.length(); i++){
						if(str[i]=='\t' && j<5){
							tabpos[j] = i;
							j++;
						}
					}
					str[tabpos[0]] = ':';
					str[tabpos[3]] = ',';
					string allele = str.substr(tabpos[2]+1,i-tabpos[2]-1);
					str.resize(tabpos[1]);
					idxS[str] = LOCI_S;
					alleleS[str] = allele;
					LOCI_S++;
				}
			}
			fin.close();
		}
		cout << LOCI_S << " loci in the SEQ_FILE." << endl; 
		foutLog << LOCI_S << " loci in the SEQ_FILE." << endl; 		
		if(tmpLOCI < 0 || tmpLOCI != LOCI_S){
			cerr << "Error: Number of loci doesn't match in '" << SEQ_SITE_FILE << "' and '" << SEQ_FILE << "'." << endl;
			foutLog << "Error: Number of loci doesn't match in '" << SEQ_SITE_FILE << "' and '" << SEQ_FILE << "'." << endl;
			flag = 0;
		}
	}		

	if(GENO_FILE.compare(default_str) != 0 && flag == 1){
		if(!get_table_dim(nrow, ncol, GENO_FILE, '\t')){
			cerr << "Error: cannot open the file '" << GENO_FILE << "'." << endl;    
			foutLog << "Error: cannot open the file '" << GENO_FILE << "'." << endl; 
			foutLog.close();
			return 1;
		}	
		REF_INDS = nrow - GENO_NON_DATA_ROWS;
		int tmpLOCI = ncol - GENO_NON_DATA_COLS;
		cout << REF_INDS << " individuals in the GENO_FILE." << endl; 
		foutLog << REF_INDS << " individuals in the GENO_FILE." << endl; 		 		
		if(REF_INDS < 0){
			cerr << "Error: Invalid number of rows in '" << GENO_FILE << "'." << endl;
			foutLog << "Error: Invalid number of rows in '" << GENO_FILE << "'." << endl;
			flag = 0;
		}
		if(tmpLOCI < 0){
			cerr << "Error: Invalid number of columns in the GENO_FILE '" << GENO_FILE << "'." << endl;
			foutLog << "Error: Invalid number of columns in the GENO_FILE '" << GENO_FILE << "'." << endl;
			flag = 0;
		}			
		GENO_SITE_FILE = GENO_FILE;
		GENO_SITE_FILE.replace(GENO_SITE_FILE.length()-5, 5, ".site");
		fin.open(GENO_SITE_FILE.c_str());
		if(fin.fail()){
			cerr << "Error: cannot open the file '" << GENO_SITE_FILE << "'." << endl;    
			foutLog << "Error: cannot open the file '" << GENO_SITE_FILE << "'." << endl;   
			foutLog.close();
			return 1;
		}else{
			getline(fin, str);
			LOCI_G = 0;
			LOCI = 0;
			while(!fin.eof()){
				getline(fin, str);
				if(str.length()>0 && str!=" "){
					j=0;
					int tabpos[5];
					for(i=0; i<str.length(); i++){
						if(str[i]=='\t' && j<5){
							tabpos[j] = i;
							j++;
						}
					}
					str[tabpos[0]] = ':';
					str[tabpos[3]] = ',';
					string allele = str.substr(tabpos[2]+1,i-tabpos[2]-1);
					string snpID = str.substr(tabpos[1]+1,tabpos[2]-tabpos[1]-1);
					str.resize(tabpos[1]);
					idxG[str] = LOCI_G;
					alleleG[str] = allele;
					LOCI_G++;
					
					if(SEQ_FILE.compare(default_str)!=0 && PCA_MODE==0){					
						if(idxS.count(str)>0){
							if(alleleS[str].compare(alleleG[str]) != 0){
								cerr << "Warning: Two datasets have different alleles at locus [" << str << "]: ";
								cerr << "[" << alleleG[str] << "] vs [" << alleleS[str] << "]." << endl;
								foutLog << "Warning: Two datasets have different alleles at locus [" << str << "]: ";
								foutLog << "[" << alleleG[str] << "] vs [" << alleleS[str] << "]." << endl;
								unmatchSite++;
							}else{
								cmnsnp.push_back(str);
								LOCI++;
							}
						}
					}else{
						cmnsnp.push_back(str);
						LOCI++;
					}
				}	
			}
			fin.close();
		}
		cout << LOCI_G << " loci in the GENO_FILE." << endl; 
		foutLog << LOCI_G << " loci in the GENO_FILE." << endl;
		
		if(tmpLOCI < 0 || tmpLOCI != LOCI_G){
			cerr << "Error: Number of loci doesn't match in '" << GENO_SITE_FILE << "' and '" << GENO_FILE << "'." << endl;
			foutLog << "Error: Number of loci doesn't match in '" << GENO_SITE_FILE << "' and '" << GENO_FILE << "'." << endl;	
			flag = 0;
		}else if(LOCI>0){
			if(SEQ_FILE.compare(default_str)!=0 && PCA_MODE==0){
				cmnS.set_size(LOCI);
				cmnG.set_size(LOCI);
				for(i=0; i<LOCI; i++){
					cmnG(i)=idxG[cmnsnp[i]];
					cmnS(i)=idxS[cmnsnp[i]];
				}
			}else{
				cmnG.set_size(LOCI);
				for(i=0; i<LOCI; i++){
					cmnG(i)=idxG[cmnsnp[i]];
				}				
			}
		}
		cmnsnp.clear();
		idxG.clear();
		idxS.clear();
		alleleG.clear();
		alleleS.clear();
	}else{
		cerr << "Error: GENO_FILE (-g) is not specified." << endl;
		foutLog << "Error: GENO_FILE (-g) is not specified." << endl;
		foutLog.close();
		gsl_rng_free(rng);
		return 1;
	}
	
	if(COORD_FILE.compare(default_str) != 0 && GENO_FILE.compare(default_str) != 0 && flag == 1){
		if(!get_table_dim(nrow, ncol, COORD_FILE, '\t')){
			cerr << "Error: cannot open the COORD_FILE '" << COORD_FILE << "'." << endl;    
			foutLog << "Error: cannot open the COORD_FILE '" << COORD_FILE << "'." << endl; 
			foutLog.close();
			return 1;
		}	
		int tmpINDS = nrow - COORD_NON_DATA_ROWS;
		NUM_PCS = ncol - COORD_NON_DATA_COLS;
		cout << tmpINDS << " individuals in the COORD_FILE." << endl;
		cout << NUM_PCS << " PCs in the COORD_FILE." << endl;
		foutLog << tmpINDS << " individuals in the COORD_FILE." << endl;
		foutLog << NUM_PCS << " PCs in the COORD_FILE." << endl;
		
		if(tmpINDS < 0){
			cerr << "Error: Invalid number of rows in the COORD_FILE " << COORD_FILE << "." << endl;
			foutLog << "Error: Invalid number of rows in the COORD_FILE " << COORD_FILE << "." << endl;
			flag = 0;
		}else if(tmpINDS != REF_INDS && REF_INDS >= 0){
			cout << tmpINDS << "\t" << REF_INDS << endl;
			cerr << "Error: Number of individuals in the COORD_FILE is not the same as in the GENO_FILE." << endl;
			foutLog << "Error: Number of individuals in the COORD_FILE is not the same as in the GENO_FILE." << endl;
			flag = 0;
		}
		if(NUM_PCS < 0){
			cerr << "Error: Invalid number of columns in the COORD_FILE " << COORD_FILE << "." << endl;
			foutLog << "Error: Invalid number of columns in the COORD_FILE " << COORD_FILE << "." << endl;
			flag = 0;
		}	
	}
	if(flag == 0){
		foutLog.close();
		return 1;
	}
	
	//################  Set default values to some parameters #######################
	if(REF_SIZE==default_int){ REF_SIZE = REF_INDS; }
	if(FIRST_IND==default_int){ FIRST_IND = 1; }
	if(LAST_IND==default_int){ LAST_IND = SEQ_INDS; }	

	// ################### Check Parameters #############################
	flag = check_parameters();
	if(MIN_LOCI <= DIM && PCA_MODE==0){
		cerr << "Warning: DIM>=MIN_LOCI is found; DIM=" << DIM << ", MIN_LOCI=" << MIN_LOCI << "." << endl; 
		cerr << "Reset MIN_LOCI to DIM+1: MIN_LOCI=" << DIM+1 << "." << endl;
		foutLog << "Warning: DIM>=MIN_LOCI is found; DIM=" << DIM << ", MIN_LOCI=" << MIN_LOCI << "." << endl; 
		foutLog << "Reset MIN_LOCI to DIM+1: MIN_LOCI=" << DIM+1 << "." << endl;
		MIN_LOCI = DIM+1;
	}
	if(REF_SIZE > REF_INDS){
		cerr << "Warning: REF_SIZE>REF_INDS is found; REF_SIZE=" << REF_SIZE << ", REF_INDS=" << REF_INDS << "." << endl; 
		cerr << "Reset REF_SIZE to REF_INDS: REF_SIZE=" << REF_INDS << "." << endl;
		foutLog << "Warning: REF_SIZE>REF_INDS is found; REF_SIZE=" << REF_SIZE << ", REF_INDS=" << REF_INDS << "." << endl; 
		foutLog << "Reset REF_SIZE to REF_INDS: REF_SIZE=" << REF_INDS << "." << endl;
		REF_SIZE = REF_INDS;
	}	
	if(REPS==1){
		OUTPUT_REPS = 0;
	}
	if(flag==0){
		foutLog.close();
		return 1;
	}else{
		print_configuration();
	}
	// Setting significance cutoff for the Tracy-Widom statistic
	if(ALPHA==0.2){
		TW = -0.1653;
	}else if(ALPHA==0.15){
		TW = 0.1038;
	}else if(ALPHA==0.1){
		TW = 0.4501;
	}else if(ALPHA==0.05){
		TW = 0.9793;
	}else if(ALPHA==0.01){
		TW = 2.0233;
	}else if(ALPHA==0.005){
		TW = 2.4221;
	}else if(ALPHA==0.001){
		TW = 3.2712;
	}
	// Set number of threads 
	openblas_set_num_threads(NUM_THREADS);
	
	// #####################  Check data format ############################
	if(CHECK_FORMAT != 0){
		int flag1 = 1;
		int flag2 = 1;
		int flag3 = 1;
 		time ( &rawtime );
  		timeinfo = localtime ( &rawtime );
		cout << endl << asctime (timeinfo);
		cout << "Checking data format ..." << endl;
		foutLog << endl << asctime (timeinfo);
		foutLog << "Checking data format ..." << endl;
		if(CHECK_FORMAT==1 || CHECK_FORMAT==2 || CHECK_FORMAT==10 || CHECK_FORMAT==20){
			if(GENO_FILE.compare(default_str)!=0){
				flag1 = check_format_geno(GENO_FILE, REF_INDS, LOCI_G);
				if(flag1==1){
					cout << "GENO_FILE: OK." << endl;
					foutLog << "GENO_FILE: OK." << endl;				
				}
			}else{
					cout << "GENO_FILE: not specified." << endl;
					foutLog << "GENO_FILE: not specified." << endl;
			}
		}
		if(CHECK_FORMAT==1 || CHECK_FORMAT==3 || CHECK_FORMAT==10 || CHECK_FORMAT==30){
			if(SEQ_FILE.compare(default_str)!=0){
				flag2 = check_format_seq(SEQ_FILE, SEQ_INDS, LOCI_S);
				if(flag2==1){
					cout << "SEQ_FILE: OK." << endl;
					foutLog << "SEQ_FILE: OK." << endl;
				}
			}else{
					cout << "SEQ_FILE: not specified." << endl;
					foutLog << "SEQ_FILE: not specified." << endl;
			}
		}
		if(CHECK_FORMAT==1 || CHECK_FORMAT==4 || CHECK_FORMAT==10 || CHECK_FORMAT==40){
			if(COORD_FILE.compare(default_str)!=0){
				flag3 = check_format_coord(COORD_FILE, REF_INDS, NUM_PCS);
				if(flag3==1){
					cout << "COORD_FILE: OK." << endl;
					foutLog << "COORD_FILE: OK." << endl;
				}
			}else{
					cout << "COORD_FILE: not specified." << endl;
					foutLog << "COORD_FILE: not specified." << endl;
			}
		}
		if(flag1==0 || flag2==0 || flag3==0){
			foutLog.close();
			gsl_rng_free(rng);	
			return 1;
		}
		if(CHECK_FORMAT < 5){
			time ( &rawtime );
 		 	timeinfo = localtime ( &rawtime );			
			runningtime = timer.toc();
			cout << endl << "Finished at: " << asctime (timeinfo);
			cout << "Total wall clock time: " << runningtime << " seconds." << endl;
			cout << "====================================================================" <<endl;
			foutLog << endl << "Finished at: " << asctime (timeinfo);
			foutLog << "Total wall clock time: " << runningtime << " seconds." << endl;
			foutLog << "====================================================================" <<endl;
			foutLog.close(); 
			gsl_rng_free(rng);	
			return 0;
		}	
	}
	
	// === get the index of the reference individuals ===
	uvec Refset(REF_SIZE);
	int Index[REF_INDS], subset[REF_SIZE];
	for(i=0; i<REF_INDS; i++){
		Index[i] = i;
	}
	if(REF_SIZE >0){
		if(REF_SIZE < REF_INDS){
			gsl_ran_choose(rng, subset, REF_SIZE, Index, REF_INDS, sizeof (int));
			for(i=0; i<REF_SIZE; i++){
				Refset(i) = subset[i];
			}
			Refset = sort(Refset);
			time ( &rawtime );
			timeinfo = localtime ( &rawtime );
			cout << endl << asctime (timeinfo);
			cout << "Randomly select " << REF_SIZE << " reference individuals (-N)." << endl;
			foutLog << "Randomly select " << REF_SIZE << " reference individuals (-N)." << endl;
		}else{
			for(i=0; i<REF_INDS; i++){
				Refset(i) = i;
			}
		}
	}
	
	// ###################### Exclude SNPs in the excluding list #################
	urowvec ExLoci = zeros<urowvec>(LOCI);	

	if(EXCLUDE_LIST.compare(default_str) != 0){
		fin.open(EXCLUDE_LIST.c_str());
		if(fin.fail()){
			cerr << "Error: cannot open the file '" << EXCLUDE_LIST << "'." << endl;    
			foutLog << "Error: cannot open the file '" << EXCLUDE_LIST << "'." << endl;   
			foutLog.close();
			gsl_rng_free(rng);	
			return 1;
		}else{
			map<string,int> exSNP;
			while(!fin.eof()){
				fin >> str;
				if(str.length()>0 && str!=" "){
					exSNP[str] = 1;
				}				
			}
			fin.close();
					
			fin.open(GENO_SITE_FILE.c_str());
			if(fin.fail()){
				cerr << "Error: cannot open the file '" << GENO_SITE_FILE << "'." << endl;    
				foutLog << "Error: cannot open the file '" << GENO_SITE_FILE << "'." << endl;   
				foutLog.close();
				gsl_rng_free(rng);	
				return 1;
			}else if(exSNP.size()>0){
				getline(fin, str);
				string snpID;
				k = 0;
				for(j=0; j<LOCI_G; j++){
					fin >> str >> str >> snpID;
					getline(fin, str);
					if(j==cmnG(k)){
						if(exSNP.count(snpID)>0 && ExLoci(k) == 0){
							ExLoci(k) = 1;
							LOCI_ex++;
						}
						k++;
					}
					if(k==LOCI){ break;}
				}
			}
			fin.close();
			Lex += LOCI_ex;
		}
	}
	if(TRIM_PROP > 0){
		for(int j=0; j<LOCI; j++){
			if(gsl_rng_uniform(rng)<TRIM_PROP && ExLoci(j) == 0){
				ExLoci(j) = 1;
				LOCI_trim++;
			}
		}
		Lex += LOCI_trim;
	}
		
	// #####################  Check coverage ############################
	int Ls = 0;
	int Lg = 0;	
	if((CHECK_COVERAGE != 0 || MAX_COVERAGE>0 || MIN_COVERAGE>0) && PCA_MODE == 0){
		time ( &rawtime );
  		timeinfo = localtime ( &rawtime );
  		cout << endl << asctime (timeinfo);
		cout << "Checking coverage in the sequence data ..." << endl;
  		foutLog << endl << asctime (timeinfo);
		foutLog << "Checking coverage in the sequence data ..." << endl;

		int flag = check_coverage(CHECK_COVERAGE, FIRST_IND, LAST_IND, cmnS, ExLoci, Ls, Lg);
		if(flag == 0){	
			foutLog.close();
			gsl_rng_free(rng);	
			return 1;
		}else if(CHECK_COVERAGE == 2){	
			time ( &rawtime );
 		 	timeinfo = localtime ( &rawtime );
			runningtime = timer.toc();
			cout << endl << "Finished at: " << asctime (timeinfo);
			cout << "Total wall clock time: " << runningtime << " seconds." << endl;
			cout << "====================================================================" <<endl;
			foutLog << endl << "Finished at: " << asctime (timeinfo);
			foutLog << "Total wall clock time: " << runningtime << " seconds." << endl;
			foutLog << "====================================================================" <<endl;
			foutLog.close();
			gsl_rng_free(rng);	
			return 0;			
		}else{
			Lex += Lg;
			Lex += Ls;
		}
	}

	// ##################################################################
	
	int LOCI_in = LOCI-Lex;	
	if(LOCI_in<LOCI_G || SEQ_FILE.compare(default_str)!=0){
		time ( &rawtime );
		timeinfo = localtime ( &rawtime );
		cout << endl << asctime (timeinfo);	
		if(SEQ_FILE.compare(default_str) != 0){
			cout << "Identify " << LOCI+unmatchSite << " loci shared by SEQ_FILE and GENO_FILE." << endl;
			foutLog << "Identify " << LOCI+unmatchSite << " loci shared by SEQ_FILE and GENO_FILE." << endl;
		}
		if(unmatchSite>0){
			cout << "Exclude "<< unmatchSite << " loci that have different alleles in two datasets." << endl;
			foutLog << "Exclude "<< unmatchSite << " loci that have different alleles in two datasets." << endl;
		}
		if(LOCI_ex>0){
			cout << "Exclude " << LOCI_ex << " loci given by the EXCLUDE_LIST (-ex)." << endl;
			foutLog << "Exclude " << LOCI_ex << " loci given by the EXCLUDE_LIST (-ex)." << endl;
		}
		if(TRIM_PROP>0){
			cout << "Exclude " << LOCI_trim << " random loci by the TRIM_PROP (-M) option." << endl;
			foutLog << "Exclude " << LOCI_trim << " random loci by the TRIM_PROP (-M) option." << endl;
		}		
		if(MAX_COVERAGE>0 && PCA_MODE==0){
			cout << "Exclude " << Lg << " loci with mean coverage >" << MAX_COVERAGE << "X (-maxc)." << endl;
			foutLog << "Exclude " << Lg << " loci with mean coverage >" << MAX_COVERAGE << "X (-maxc)." << endl;	
		}
		if(MIN_COVERAGE>0 && PCA_MODE==0){
			cout << "Exclude " << Ls << " loci with mean coverage <" << MIN_COVERAGE << "X (-minc)." << endl;
			foutLog << "Exclude " << Ls << " loci with mean coverage <" << MIN_COVERAGE << "X (-minc)." << endl;
		}
		cout << "The analysis will base on the remaining " << LOCI_in << " loci." << endl;		
		foutLog << "The analysis will base on the remaining " << LOCI_in << " loci." << endl;
		
		if(LOCI_in==0){
			cout << "Error: No data for the analysis. Program exit!" << endl;
			foutLog << "Error: No data for the analysis. Program exit!" << endl;
			foutLog.close();
			gsl_rng_free(rng);
			return 1;
		}
	}
	// ###################################################################
	// == Variables to be saved for reused ==
	string *RefInfo1 = new string [REF_SIZE];
	string *RefInfo2 = new string [REF_SIZE];
	Mat<char> RefG(REF_SIZE, LOCI_in);
	mat refPC = zeros<mat>(REF_SIZE,DIM);
	//========================= Read reference data ==========================
	fin.open(GENO_FILE.c_str());
	if(fin.fail()){
		cerr << "Error: cannot find the GENO_FILE '" << GENO_FILE << "'." << endl;    
		foutLog << "Error: cannot find the GENO_FILE '" << GENO_FILE << "'." << endl;   
		foutLog.close();
		gsl_rng_free(rng);	
		return 1;
	}
 	time ( &rawtime );
  	timeinfo = localtime ( &rawtime );
  	cout << endl << asctime (timeinfo);	
	cout << "Reading reference genotypes ..." << endl;
  	foutLog << endl << asctime (timeinfo);
	foutLog << "Reading reference genotypes ..." << endl;
	for(i=0; i<GENO_NON_DATA_ROWS; i++){
		getline(fin, str);          // Read non-data rows
	}
	int ii = 0;
	for(i=0; i<REF_INDS; i++){
		if(i==Refset[ii]){
			fin >> RefInfo1[ii] >> RefInfo2[ii];
			for(j=2; j<GENO_NON_DATA_COLS; j++){
				fin >> str;
			}
			frowvec tmpG = zeros<frowvec>(LOCI_G);
			for(j=0; j<LOCI_G; j++){
				fin >> tmpG(j);    // Read genotype data
			}
			k = 0;
			for(j=0; j<LOCI; j++){
				if(ExLoci(j)==0){
					RefG(ii,k) = tmpG(cmnG(j));
					k++;
				}
			}
			getline(fin, str); // get the rest of the line
			ii++;
		}else{
			getline(fin, str);
		}
		if(ii==REF_SIZE){
			break;
		}
	}
		
	if(!fin.good()){
		fin.close();
		cerr << "Error: ifstream error occurs when reading the GENO_FILE." << endl;
		cerr << "Run 'laser -fmt 2' to check the GENO_FILE '" << GENO_FILE << "'." << endl;
		foutLog << "Error: ifstream error occurs when reading the GENO_FILE." << endl;
		foutLog << "Run 'laser -fmt 2' to check the GENO_FILE '" << GENO_FILE << "'." << endl;
		foutLog.close();
		gsl_rng_free(rng);	
		return 1;
	}		
	fin.close();
	//========================= Get reference coordinates  ==========================
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	if(COORD_FILE.compare(default_str)!=0 && PCA_MODE==0){		
		fin.open(COORD_FILE.c_str());
		if(fin.fail()){
			cerr << "Error: cannot find the COORD_FILE '" << COORD_FILE << "'." << endl;
			foutLog << "Error: cannot find the COORD_FILE '" << COORD_FILE << "'." << endl;  
			foutLog.close();
			gsl_rng_free(rng);	
			return 1;
		}else{
			cout << endl << asctime (timeinfo);
			cout << "Reading reference PCA coordinates ..." << endl;
			foutLog << endl << asctime (timeinfo);
			foutLog << "Reading reference PCA coordinates ..." << endl;
			for(i=0; i<COORD_NON_DATA_ROWS; i++){
				getline(fin, str);          // Read non-data rows
			}
			int ii = 0;
			for(i=0; i<REF_INDS; i++){
				if(i==Refset[ii]){
					string popstr;
					string indstr;
					fin >> popstr >> indstr;
					if(popstr.compare(RefInfo1[ii])!=0 || indstr.compare(RefInfo2[ii])!=0){
						cerr << "Error: ID of individual " << i+1 << " in the COORD_FILE differs from that in the GENO_FILE." << endl;
						foutLog << "Error: ID of individual " << i+1 << " in the COORD_FILE differs from that in the GENO_FILE." << endl;
						fin.close();
						foutLog.close();
						gsl_rng_free(rng);	
						return 1;
					}
					for(j=2; j<COORD_NON_DATA_COLS; j++){
						fin >> str;
					}
					for(j=0; j<DIM; j++){
						fin >> refPC(ii,j);    // Read reference coordinates
					}
					getline(fin, str);  // Read the rest of the line
					ii++;
				}else{
					getline(fin, str);
				}
				if(ii==REF_SIZE){
					break;
				}
			}
		}
		if(!fin.good()){
			fin.close();
			cerr << "Error: ifstream error occurs when reading the COORD_FILE." << endl;
			cerr << "Run 'laser -fmt 4' to check the COORD_FILE '" << COORD_FILE << "'." << endl;
			foutLog << "Error: ifstream error occurs when reading the COORD_FILE." << endl;
			foutLog << "Run 'laser -fmt 4' to check the COORD_FILE '" << COORD_FILE << "'." << endl;
			foutLog.close();
			gsl_rng_free(rng);	
			return 1;
		}
		fin.close();
	}else{                
		cout << endl << asctime (timeinfo);
		cout << "Performing PCA on reference genotypes ..." << endl;
		foutLog << endl << asctime (timeinfo);
		foutLog << "Performing PCA on reference genotypes ..." << endl;	
		rowvec PCvar = zeros<rowvec>(DIM);
		fmat Gm;
		fmat Gsd;
		fmat W;
		if(PCA_MODE==1){
			mat GRM(REF_SIZE, REF_SIZE);
			pca(conv_to<fmat>::from(RefG), DIM, refPC, PCvar, GRM);   // Perform PCA based on EVD
			GRM = GRM/LOCI_in;
			outfile = OUT_PREFIX;
			outfile.append(".RefPC.grm");
			fout.open(outfile.c_str());
			if(fout.fail()){
				cerr << "Error: cannot create a file named " << outfile << "." << endl;
				foutLog << "Error: cannot create a file named " << outfile << "." << endl;  
				foutLog.close();
				gsl_rng_free(rng);	
				return 1;
			}		
			for(i=0; i<REF_SIZE; i++){
				for(j=0; j<REF_SIZE-1; j++){
					fout << GRM(i,j) << "\t";
				}
				fout << GRM(i,REF_SIZE-1) << endl;	
			}		
			fout.close();
			cout << "Genetic relationship matrix is output to '" << outfile << "'." << endl;
			foutLog << "Genetic relationship matrix is output to '" << outfile << "'." << endl;				
			
			//======== Calculating Z scores for reference individuals =========
			outfile = OUT_PREFIX;
			outfile.append(".RefPC.zscore");
			fout.open(outfile.c_str());
			if(fout.fail()){
				cerr << "Error: cannot create a file named " << outfile << "." << endl;
				foutLog << "Error: cannot create a file named " << outfile << "." << endl;  
				foutLog.close();
				gsl_rng_free(rng);	
				return 1;
			}				
			fout << "popID\tindivID\tZ\tW" << endl;
			for(i=0; i<REF_SIZE; i++){
				vec d1 = zeros<vec>(REF_SIZE);
				for(j=0; j<REF_SIZE; j++){
					if(j!=i){
						rowvec v = refPC.row(i)-refPC.row(j);
						for(k=0; k<DIM; k++) d1(j) += v(k)*v(k);
					}else{
						d1(j) = -9;
					}
				}
				uvec idx = sort_index(d1);
				vec Mk = zeros<vec>(KNN_ZSCORE);
				for(j=0; j<KNN_ZSCORE; j++) Mk(j) = GRM(idx(j+1),idx(j+1));
				idx = sort_index(Mk);
				double IQR = Mk(idx(round(KNN_ZSCORE*0.75)))-Mk(idx(round(KNN_ZSCORE*0.25)));
				double Z = (GRM(i,i)-mean(Mk))/stddev(Mk);
				double W = (GRM(i,i)-median(Mk))/IQR;
				fout << RefInfo1[i] << "\t" << RefInfo2[i] << "\t" << Z << "\t" << W << endl;
			}
			fout.close();
			cout << "Z scores of reference individuals are output to '" << outfile << "'." << endl;
			foutLog << "Z scores of reference individuals are output to '" << outfile << "'." << endl;
			//=====================================================================
			
		}else if(PCA_MODE==2){
			pca_geno(RefG, DIM, refPC, PCvar);   // Perform PCA	(slow but memory efficient)
		}else if(PCA_MODE==3){
			pca_svd(conv_to<fmat>::from(RefG), DIM, refPC, PCvar, Gm, Gsd, W); // Perform PCA based on SVD
		}else{
			mat GRM(REF_SIZE, REF_SIZE);
			pca(conv_to<fmat>::from(RefG), DIM, refPC, PCvar, GRM);   // Perform PCA based on EVD
		}
		//==================== Output reference PCs ==========================
		outfile = OUT_PREFIX;
		outfile.append(".RefPC.coord");
		fout.open(outfile.c_str());
		if(fout.fail()){
			cerr << "Error: cannot create a file named " << outfile << "." << endl;
			foutLog << "Error: cannot create a file named " << outfile << "." << endl;  
			foutLog.close();
			gsl_rng_free(rng);	
			return 1;
		}		
		fout << "popID\tindivID\t";
		for(j=0; j<DIM-1; j++){ fout << "PC" << j+1 << "\t"; }
		fout << "PC" << DIM << endl;
		for(i=0; i<REF_SIZE; i++){
			fout << RefInfo1[i] << "\t" << RefInfo2[i] << "\t";
			for(j=0; j<DIM-1; j++){
				fout << refPC(i,j) << "\t";
			}
			fout << refPC(i,DIM-1) << endl;	
		}		
		fout.close();
		cout << "Reference PCA coordinates are output to '" << outfile << "'." << endl;
		foutLog << "Reference PCA coordinates are output to '" << outfile << "'." << endl;
		//==================================================================
		outfile = OUT_PREFIX;
		outfile.append(".RefPC.var");
		fout.open(outfile.c_str());
		if(fout.fail()){
			cerr << "Error: cannot create a file named " << outfile << "." << endl;
			foutLog << "Error: cannot create a file named " << outfile << "." << endl;  
			foutLog.close();
			gsl_rng_free(rng);	
			return 1;
		}
		fout << "PC" << "\t" << "Variance(%)" << endl;
		for(j=0; j<DIM; j++){ 
			fout << j+1 << "\t" << PCvar(j) << endl;
		}	
		fout.close();
		PCvar.clear();
		cout << "Variances explained by PCs are output to '" << outfile << "'." << endl;
		foutLog << "Variances explained by PCs are output to '" << outfile << "'." << endl;
		//===================================================================
		if(PCA_MODE==3){
			fin.open(GENO_SITE_FILE.c_str());
			if(fin.fail()){
				cerr << "Error: cannot find the SITE_FILE '" << GENO_SITE_FILE << "'." << endl;
				foutLog << "Error: cannot find the SITE_FILE '" << GENO_SITE_FILE << "'." << endl;
				foutLog.close();
				gsl_rng_free(rng);	
				return 1;
			}
			getline(fin, str);
			
			outfile = OUT_PREFIX;
			outfile.append(".RefPC.load");
			fout.open(outfile.c_str());
			if(fout.fail()){
				cerr << "Error: cannot create a file named " << outfile << "." << endl;
				foutLog << "Error: cannot create a file named " << outfile << "." << endl;  
				foutLog.close();
				gsl_rng_free(rng);	
				return 1;
			}		
			fout << "ID\tMean\tSd\t";
			for(k=0; k<DIM-1; k++){ fout << "PC" << k+1 << "\t"; }
			fout << "PC" << DIM << endl;
			i=0;
			for(j=0; j<LOCI; j++){
				fin >> str >> str >> str;
				if(ExLoci(j)==0){
					fout << str << "\t" << Gm(i) << "\t" << Gsd(i) << "\t";
					for(k=0; k<DIM-1; k++){ fout << W(i,k) << "\t"; }; 
					fout << W(i,DIM-1) << endl;
					i++;
				}
				getline(fin, str);
			}		
			fout.close();
			cout << "Reference PCA loadings are output to '" << outfile << "'." << endl;
			foutLog << "Reference PCA loadings are output to '" << outfile << "'." << endl;
			fin.close();
		}	
		//===================================================================
		if(PCA_MODE > 0){                                 // If performing PCA only
 	 		delete [] RefInfo1;
			delete [] RefInfo2;
			time ( &rawtime );
 		 	timeinfo = localtime ( &rawtime );
			runningtime = timer.toc();
			cout << endl << "Finished at: " << asctime (timeinfo);
			cout << "Total wall clock time: " << runningtime << " seconds." << endl;
			cout << "====================================================================" <<endl;
			foutLog << endl << "Finished at: " << asctime (timeinfo);
			foutLog << "Total wall clock time: " << runningtime << " seconds." << endl;
			foutLog << "====================================================================" <<endl;
			foutLog.close();
			gsl_rng_free(rng);	
			return 0;	
		}
	}
	delete [] RefInfo1;
	delete [] RefInfo2;
					
	//========================= Read sequence data ==========================
		
	fin.open(SEQ_FILE.c_str());
	if(fin.fail()){
		cerr << "Error: cannot find the SEQ_FILE '"<< SEQ_FILE <<"'." << endl;
		foutLog << "Error: cannot find the SEQ_FILE '"<< SEQ_FILE <<"'." << endl;
		foutLog.close();
		gsl_rng_free(rng);		
		return 1;
	}
	//==== Open output file ====
	outfile = OUT_PREFIX;
	outfile.append(".SeqPC.coord");
	fout.open(outfile.c_str());
	if(fout.fail()){
		cerr << "Error: cannot create a file named " << outfile << "." << endl;
		foutLog << "Error: cannot create a file named " << outfile << "." << endl;
		foutLog.close();
		gsl_rng_free(rng);		
		return 1;
	}
	fout << "popID\t" << "indivID\t" << "L1\t" << "Ci\t" << "K\t" << "t\t" << "Z\t";
	for(j=0; j<DIM-1; j++){ fout << "PC" << j+1 << "\t"; }
	fout << "PC" << DIM << endl;

	if(REPS>1){
		outfile2 = outfile;
		outfile2.append(".sd");
		fout2.open(outfile2.c_str());
		if(fout2.fail()){
			cerr << "Error: cannot create a file named " << outfile2 << "." << endl;
			foutLog << "Error: cannot create a file named " << outfile2 << "." << endl;
			foutLog.close(); 
			fout.close();
			gsl_rng_free(rng);		
			return 1;
		}	
		fout2 << "popID\t" << "indivID\t" << "t.sd\t" << "Z.sd\t";
		for(j=0; j<DIM-1; j++){
			fout2 << "PC" << j+1 << ".sd\t";
		}
		fout2 << "PC" << DIM << ".sd" << endl;
		if(OUTPUT_REPS==1){
			outfile3 = outfile;
			outfile3.append(".reps");
			fout3.open(outfile3.c_str());
			if(fout3.fail()){
				cerr << "Error: cannot create a file named " << outfile3 << "." << endl;
				foutLog << "Error: cannot create a file named " << outfile3 << "." << endl;
				foutLog.close(); 
				fout.close();
				fout2.close();
				gsl_rng_free(rng);		
				return 1;
			}	
			fout3 << "popID\t" << "indivID\t" << "L1\t" << "Ci\t" << "K\t" << "t\t" << "Z\t";
			for(j=0; j<DIM-1; j++){
				fout3 << "PC" << j+1 << "\t";
			}
			fout3 << "PC" << DIM << endl;
		}
	}
	//==========================
 	time ( &rawtime );
  	timeinfo = localtime ( &rawtime );
  	cout << endl << asctime (timeinfo);
	cout << "Analyzing sequence samples ..." << endl;
  	foutLog << endl << asctime (timeinfo);	
	foutLog << "Analyzing sequence samples ..." << endl;
	
	if(DIM_HIGH == 0){
		AUTO_MODE = true;
	}
	
	for(i=0; i<SEQ_NON_DATA_ROWS; i++){
		getline(fin, str);          // Read non-data rows
	}

//	#pragma omp parallel for private(i)
	for(i=1; i<=LAST_IND; i++){
		if(i<FIRST_IND){
			getline(fin, str);
		}else{
			string SeqInfo1;
			string SeqInfo2;
			urowvec tmpC(LOCI_S);    // Coverage of one sample
			frowvec tmpS(LOCI_S);    // Sequence read of one sample
			frowvec tmpQ(LOCI_S);    // Base quality of one sample
						
			fin >> SeqInfo1 >> SeqInfo2;	
			for(j=2; j<SEQ_NON_DATA_COLS; j++){
				fin >> str;
			}
			for(j=0; j<LOCI_S; j++){				
				fin >> tmpC(j) >> tmpS(j) >> tmpQ(j);					
 				if(tmpC(j)<0 || tmpS(j)<0 || tmpS(j)>tmpC(j) || tmpQ(j)<0){
					if(!fin.good()){
						fin.close();
						cerr << "Error: ifstream error occurs when reading the SEQ_FILE." << endl;
						cerr << "Run 'laser -fmt 3' to check the SEQ_FILE '" << SEQ_FILE << "'." << endl;
						foutLog << "Error: ifstream error occurs when reading the SEQ_FILE." << endl;
						foutLog << "Run 'laser -fmt 3' to check the SEQ_FILE '" << SEQ_FILE << "'." << endl;
						foutLog.close();
						gsl_rng_free(rng);	
						return 1;
					}else{
						fin.close();
						cerr << "Error: invalid value at locus "<< j+1 << " of individual " << i << " in the SEQ_FILE." << endl;
						foutLog << "Error: invalid value at locus "<< j+1 << " of individual " << i << " in the SEQ_FILE." << endl;
						foutLog.close();
						gsl_rng_free(rng);		
						return 1;
					}
				} 
			}
			
			urowvec C(LOCI_in);    // Coverage of one sample
			frowvec S(LOCI_in);    // Sequence read of one sample
			frowvec Q(LOCI_in);    // Base quality of one sample
			uvec LOC(LOCI_in);    // Indices of covered loci
			double meanC = 0;
			int Lcov = 0;     // Number of loci with nonzero coverage
			k = 0;
			for(j=0; j<LOCI; j++){				
				if(ExLoci(j)==0){
					C(k) = tmpC(cmnS(j));
					C(k) = tmpC(cmnS(j));
					S(k) = tmpS(cmnS(j));
					Q(k) = tmpQ(cmnS(j));					
					if(C(k)>0){
						LOC(Lcov) = k;
						Lcov++;
						meanC += C(k);
					}
					k++;
				}
			}
			tmpC.clear();
			tmpS.clear();
			tmpQ.clear();
			LOC.resize(Lcov);
			meanC = meanC/LOCI_in;
			
			uvec Loc;
			int Linc;    // Number of loci to include in the computation of PCA
			if(Lcov>MAX_LOCI){  // Randomly excluding loci if Lcov >  MAX_LOCI
				Linc = MAX_LOCI;
				vec v = randu<vec>(Lcov);
				uvec idx = sort_index(v);			
				idx.resize(MAX_LOCI);
				idx = sort(idx, "ascend");
				Loc = LOC(idx);
				cout << "Randomly select " << Linc << " out of " << Lcov << " covered loci for " << SeqInfo1 << ":" << SeqInfo2 << "." << endl;
				foutLog << "Randomly select " << Linc << " out of " << Lcov << " covered loci for " << SeqInfo1 << ":" << SeqInfo2 << "." << endl;
			}else{
				Linc = Lcov;
				Loc = LOC;
			}
			LOC.clear();
			urowvec Cc(Linc);
			frowvec Sc(Linc);
			frowvec Qc(Linc);
			for(j=0; j<Linc; j++){
				Cc(j) = C(Loc(j));
				Sc(j) = S(Loc(j));
				Qc(j) = Q(Loc(j));
			}
			C.clear();
			S.clear();
			Q.clear();
			
			if(Linc >= MIN_LOCI){
				double t_m1 = 0;
				double t_m2 = 0;
				double Z_m1 = 0;
				double Z_m2 = 0;
				double dim_high = 0;
				rowvec rotPC_m1 = zeros<rowvec>(DIM);
				rowvec rotPC_m2 = zeros<rowvec>(DIM);
				for(int rep=0; rep<REPS; rep++){
					//=================== Simulate sequence reads ======================
					fmat SS(REF_SIZE, Linc);
					if(SEQ_ERR != -1){
						simuseq(RefG, Cc, Loc, SEQ_ERR, SS, rng);
					}else{
						simuseq2(RefG, Cc, Loc, Qc, SS, rng);
					}					  
					SS.insert_rows(REF_SIZE, Sc);
					//=================== Perform PCA =================================	
					fmat SSm;
					fmat SSsd;
					normalize(SS, SSm, SSsd);	
					mat M = conv_to<mat>::from(SS*SS.t());					
					vec eigval;
					mat eigvec;
                    eig_sym(eigval, eigvec, M, "dc");	// use "divide & conquer" algorithm
					//M.clear();
					SS.clear();									
					if(AUTO_MODE){						
						// ####    Calculate Tracy-Widom Statistics and determine DIM_HIGH  ####
						// Calculation of TW statistic follows Patterson et al 2006 PLoS Genetics 
						DIM_HIGH = 0;
						double eigsum = 0;
						double eig2sum = 0;
						double eigsum2 = 0;
						for(j=0; j<REF_SIZE; j++){     // The length of eigval is REF_SIZE+1;
							eigsum += eigval(j+1);
							eig2sum += pow(eigval(j+1), 2);
						}
						for(j=0; j<REF_SIZE; j++){
							int m = REF_SIZE-j;
							if(j>0){
								eigsum -= eigval(m+1);
								eig2sum -= pow(eigval(m+1),2);
							}
							eigsum2 = eigsum*eigsum;
							double n = (m+1)*eigsum2/((m-1)*eig2sum-eigsum2);
							double nsqrt = sqrt(n-1);
							double msqrt = sqrt(m);
							double mu = pow(nsqrt+msqrt, 2)/n;
							double sigma = (nsqrt+msqrt)/n*pow(1/nsqrt+1/msqrt, 1.0/3);
							double x = (m*eigval(m)/eigsum-mu)/sigma;  // Tracy-Widom statistic
							if(x>TW){           // TW is the threshold for the Tracy-Widom statisic
								DIM_HIGH++;
							}else{
								break;
							}
						}
						if(DIM_HIGH<DIM){
							DIM_HIGH = DIM;
							cout << "Warning: DIM is greater than the number of significant PCs for study sample " << i << "." << endl;
							foutLog << "Warning: DIM is greater than the number of significant PCs for study sample " << i << "." << endl;
						}
					}
					mat simuPC(REF_SIZE, DIM_HIGH);
					rowvec PC_one = zeros<rowvec>(DIM_HIGH);					
					for(j=0; j<DIM_HIGH; j++){
						for(k=0; k<REF_SIZE; k++){
							simuPC(k,j) = eigvec(k, REF_SIZE-j)*sqrt(eigval(REF_SIZE-j));
						}
						PC_one(j) = eigvec(REF_SIZE, REF_SIZE-j)*sqrt(eigval(REF_SIZE-j));
					}				
					
					//=================  Procrustes Analysis =======================
					mat simuPC_rot(REF_SIZE, DIM_HIGH);
					double t;
					double rho;
					mat A(DIM_HIGH, DIM_HIGH);
					rowvec b(DIM_HIGH);
					double epsilon = pprocrustes(simuPC, refPC, simuPC_rot, t, rho, A, b, MAX_ITER, THRESHOLD, PROCRUSTES_SCALE);
					if(epsilon>THRESHOLD){
						cout << "Warning: Projection Procrustes analysis doesn't converge in " << MAX_ITER << " iterations for " << SeqInfo2 <<", epsilon=" << epsilon << "." << endl;
						foutLog << "Warning: Projection Procrustes analysis doesn't converge in " << MAX_ITER << " iterations for " << SeqInfo2 <<", epsilon=" << epsilon << "." << endl;
					}	
					simuPC.clear();
					simuPC_rot.clear();
					rowvec rotPC_one = rho*PC_one*A+b;
					if(DIM_HIGH > DIM){
						rotPC_one.shed_cols(DIM, DIM_HIGH-1);
					}
					
					//== Calculating Z score to indicate if an individual's ancestry is represented in the reference ==
					vec d1 = zeros<vec>(REF_SIZE);
					for(j=0; j<REF_SIZE; j++){
						rowvec v = rotPC_one-refPC.row(j);
						for(k=0; k<DIM; k++) d1(j) += v(k)*v(k);
					}
					uvec idx = sort_index(d1);
					vec Mk = zeros<vec>(KNN_ZSCORE);
					for(j=0; j<KNN_ZSCORE; j++) Mk(j) = M(idx(j),idx(j));
					double Z = (M(REF_SIZE,REF_SIZE)-mean(Mk))/stddev(Mk);
					
					//================= Output Procrustes Results for one repated run ===================	
					if(REPS == 1){
						fout << SeqInfo1 << "\t" << SeqInfo2 << "\t" << Lcov << "\t" << meanC << "\t" << DIM_HIGH << "\t" << t << "\t" << Z << "\t"; 
						for(j=0; j<DIM-1; j++){ fout << rotPC_one(j) << "\t"; }
						fout << rotPC_one(DIM-1) << endl;
					}else if(REPS > 1){
						if(OUTPUT_REPS == 1){
							fout3 << SeqInfo1 << "\t" << SeqInfo2 << "\t" << Lcov << "\t" << meanC << "\t" << DIM_HIGH << "\t" << t << "\t" << Z << "\t"; 
							for(j=0; j<DIM-1; j++){ fout3 << rotPC_one(j) << "\t"; }
							fout3 << rotPC_one(DIM-1) << endl;
						}	
						t_m1 += t;
						t_m2 += pow(t,2);
						Z_m1 += Z;
						Z_m2 += pow(Z,2);
						rotPC_m1 = rotPC_m1 + rotPC_one;
						rotPC_m2 = rotPC_m2 + rotPC_one%rotPC_one;
						dim_high = dim_high + DIM_HIGH;
					}
				}
				Cc.clear();	
				Sc.clear();
				Qc.clear();
				Loc.clear();		
				//================= Output Procrustes Results ===================
				if(REPS>1){
					// calculate mean and sd
					double t_mean = t_m1/REPS;
					double t_sd = sqrt((t_m2-REPS*pow(t_mean,2))/(REPS-1));
					double Z_mean = Z_m1/REPS;
					double Z_sd = sqrt((Z_m2-REPS*pow(Z_mean,2))/(REPS-1));
					dim_high = dim_high/REPS;
					rowvec rotPC_mean = rotPC_m1/REPS;
					rowvec rotPC_sd = sqrt((rotPC_m2-REPS*rotPC_mean%rotPC_mean)/(REPS-1));	
					// output mean values of results
					fout << SeqInfo1 << "\t" << SeqInfo2 << "\t" << Lcov << "\t" << meanC << "\t" << dim_high << "\t" << t_mean << "\t" << Z_mean << "\t"; 
					for(j=0; j<DIM-1; j++){
						fout << rotPC_mean(j) << "\t";
					}
					fout << rotPC_mean(DIM-1) << endl;
					// output sd values of results
					fout2 << SeqInfo1 << "\t" << SeqInfo2 << "\t" << t_sd << "\t" << Z_sd << "\t"; 
					for(j=0; j<DIM-1; j++){
						fout2 << rotPC_sd(j) << "\t";
					}
					fout2 << rotPC_sd(DIM-1) << endl;
					// declare vectors
					rotPC_mean.clear();
					rotPC_sd.clear();
				}
				rotPC_m1.clear();
				rotPC_m2.clear();
			}else{
				//Too few number of loci covered. Skip computation and output "NA".
				cout << "Warning: skipping sample "<< SeqInfo2 << " (# covered loci < MIN_LOCI)." << endl;
				foutLog << "Warning: skipping sample "<< SeqInfo2 << " (# covered loci < MIN_LOCI)." << endl;
				fout << SeqInfo1 << "\t" << SeqInfo2 << "\t" << Lcov << "\t" << meanC << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t"; 
				for(j=0; j<DIM-1; j++){
					fout << "NA" << "\t";
				}
				fout << "NA" << endl;
				if(REPS>1){
					fout2 << SeqInfo1 << "\t" << SeqInfo2 << "\t" << "NA" << "\t"; 
					for(j=0; j<DIM-1; j++){ fout2 << "NA" << "\t"; }
					fout2 << "NA" << endl;
					if(OUTPUT_REPS==1){
						for(int rep=0; rep<REPS; rep++){
							fout3 << SeqInfo1 << "\t" << SeqInfo2 << "\t" << Lcov << "\t" << meanC << "\t" << "NA" << "\t" << "NA" << "\t" << "NA" << "\t"; 
							for(j=0; j<DIM-1; j++){ fout3 << "NA" << "\t"; }
							fout3 << "NA" << endl;
						}
					}
				}
			}		
			if(i%50==0){	
				cout << "Progress: finish analysis of individual " << i << "." << endl;
				foutLog << "Progress: finish analysis of individual " << i << "." << endl;			
			}
		}
	}
	if(!fin.good()){
		fin.close();
		fout.close();
		if(REPS>1){
			fout2.close();
			if(OUTPUT_REPS==1){
				fout3.close();
			}
		}	
		cerr << "Error: ifstream error occurs when reading the SEQ_FILE." << endl;
		cerr << "Run 'laser -fmt 3' to check the SEQ_FILE '" << SEQ_FILE << "'." << endl;
		foutLog << "Error: ifstream error occurs when reading the SEQ_FILE." << endl;
		foutLog << "Run 'laser -fmt 3' to check the SEQ_FILE '" << SEQ_FILE << "'." << endl;
		foutLog.close();
		gsl_rng_free(rng);	
		return 1;
	}
	fin.close();
	fout.close();
	if(REPS>1){
		fout2.close();
		cout << "Results for the sequence samples are output to:" << endl; 
		cout << "'" << outfile << "' (mean across " << REPS << " repeated runs)" << endl;
		cout << "'" << outfile2 << "' (standard deviation across " << REPS << " repeated runs)" << endl;
		foutLog << "Results for the sequence samples are output to:" << endl; 
		foutLog << "'" << outfile << "' (mean across " << REPS << " repeated runs)" << endl;
		foutLog << "'" << outfile2 << "' (standard deviation across " << REPS << " repeated runs)" << endl;
		if(OUTPUT_REPS==1){
			fout3.close();	
			cout << "'" << outfile3 << "' (results from all " << REPS << " repeated runs)" << endl;	
			foutLog << "'" << outfile3 << "' (results from all " << REPS << " repeated runs)" << endl;	
		}
	}else{
		cout << "Results for the sequence samples are output to '" << outfile << "'." << endl;
		foutLog << "Results for the sequence samples are output to '" << outfile << "'." << endl;
	}

	//#########################################################################	
	time ( &rawtime );
 	timeinfo = localtime ( &rawtime );
	runningtime = timer.toc();
	cout << endl << "Finished at: " << asctime (timeinfo);
	cout << "Total wall clock time: " << runningtime << " seconds." << endl;
	cout << "====================================================================" <<endl;
	foutLog << endl << "Finished at: " << asctime (timeinfo);
	foutLog << "Total wall clock time: " << runningtime << " seconds." << endl;
	foutLog << "====================================================================" <<endl;
	foutLog.close();
	gsl_rng_free(rng);		
	return 0;
}
//##########################################################################################################
bool parse_cmd_line(int argc, char* argv[], map<string,string> &args, map<string,int> &argi, map<string,double> &argd){
	bool flag=1;
	//Populate with default values
	args[ARG_PARAM_FILE] = default_str;
	args[ARG_GENO_FILE] = default_str;
	args[ARG_SEQ_FILE] = default_str;
	args[ARG_COORD_FILE] = default_str;
	args[ARG_EXCLUDE_LIST] = default_str;
	args[ARG_OUT_PREFIX] = default_str;	
 	argi[ARG_DIM] = default_int;
	argi[ARG_DIM_HIGH] = default_int;	
	argi[ARG_MIN_LOCI] = default_int;
	argi[ARG_MAX_LOCI] = default_int;
 	argd[ARG_SEQ_ERR] = default_double;
	argi[ARG_REF_SIZE] = default_int;
	argi[ARG_FIRST_IND] = default_int;	
	argi[ARG_LAST_IND] = default_int;
	argi[ARG_REPS] = default_int;
	argi[ARG_OUTPUT_REPS] = default_int;
	argd[ARG_TRIM_PROP] = default_double;
	argd[ARG_MIN_COVERAGE] = default_double;	
	argd[ARG_MAX_COVERAGE] = default_double;
	argi[ARG_CHECK_COVERAGE] = default_int;
	argi[ARG_CHECK_FORMAT] = default_int;	
	argi[ARG_PCA_MODE] = default_int;
	argd[ARG_ALPHA] = default_double;	
	argd[ARG_THRESHOLD] = default_double;
	argi[ARG_PROCRUSTES_SCALE] = default_int;
	argi[ARG_RANDOM_SEED] = default_int;
	argi[ARG_NUM_THREADS] = default_int;
	argi[ARG_KNN_ZSCORE] = default_int;
	
	for(int i = 1; i < argc-1; i++){
		if(args.count(argv[i]) > 0){
	  		args[argv[i]] = argv[i+1];
			i++;
		}else if(argi.count(argv[i]) > 0){
			if(is_int(argv[i+1])){
	  			argi[argv[i]] = atoi(argv[i+1]);
				i++;
			}else{
				cerr <<"Error: "<<"invalid value for "<<argv[i]<<"." << endl;
				foutLog <<"Error: "<<"invalid value for "<<argv[i]<<"." << endl;
				flag=0;
			}
		}else if(argd.count(argv[i]) > 0){
			if(is_numeric(argv[i+1])){
	  			argd[argv[i]] = atof(argv[i+1]);
				i++;
			}else{
				cerr <<"Error: "<<"invalid value for "<<argv[i]<<"." << endl;
				foutLog <<"Error: "<<"invalid value for "<<argv[i]<<"." << endl;
				flag=0;
			}
		}else{
			cerr << "Error: " << argv[i] << " is not recognized as a valid argument." << endl;
			foutLog << "Error: " << argv[i] << " is not recognized as a valid argument." << endl;
			flag=0;
		}
	}
	return flag;
}
//############## Check if a string is an integer #######################
bool is_int(string str){
	bool flag=1;
	for(int i=0; i<str.length(); i++){
		if( str[i] < '0' || str[i] > '9' && !(i==0 && str[i]=='-')){
			flag=0;         // not an integer number
		}
	}
	return flag;
}
//############## Check if a string is a number #######################
bool is_numeric(string str){
	bool flag=1;
	bool dot_flag=0;

	for(int i=0; i<str.length(); i++){
		if( str[i] < '0' || str[i] > '9' ){
			if(str[i] == '.' && dot_flag==0){
				dot_flag = 1;
			}else if(!(i==0 && str[i]=='-')){
				flag = 0;
			}
		}
	}
	return flag;
}
//######################### Standard Procrustes Analysis ##########################
int procrustes(mat &X, mat &Y, mat &Xnew, double &t, double &rho, mat &A, rowvec &b, int ps){
	int NUM = X.n_rows;
	//======================= Center to mean =======================
	mat Xm = mean(X);
	mat Ym = mean(Y);
	mat Xc = X-repmat(Xm, NUM, 1);
	mat Yc = Y-repmat(Ym, NUM, 1);
	//======================  SVD =====================
	mat C = Yc.t()*Xc;
	mat U;
	vec s;
	mat V;
	bool bflag = svd(U, s, V, C, "dc");	// use "divide & conquer" algorithm
	//bool bflag = svd(U, s, V, C);
	if(!bflag){
		cout << "Error: singular value decomposition in procrustes() fails." << endl;
		return 0;
	}
	//===================== Transformation ===================
	double trXX = trace(Xc.t()*Xc);
	double trYY = trace(Yc.t()*Yc);
	double trS = sum(s);
	A = V*U.t(); 
	if(ps==1){     // Orthogonal Procrustes analysis, match variance between X and Y
		rho = sqrt(trYY/trXX);
	}else{ 
		rho = trS/trXX;
	}
	b = Ym-rho*Xm*A;
	//============= New coordinates and similarity score ========
	Xnew = rho*X*A+repmat(b, NUM, 1);	
	mat Z = Y-Xnew;
	double d = trace(Z.t()*Z);
	double D = d/trYY;
	t = sqrt(1-D);
	return 1;
}
//######################### Projection Procrustes Analysis ##########################
double pprocrustes(mat &X, mat &Y, mat &Xnew, double &t, double &rho, mat &A, rowvec &b, int iter, double eps, int ps){
	double epsilon = 0;
	int NUM = X.n_rows;
	int DimX = X.n_cols;
	int DimY = Y.n_cols;
	int i = 0;
	if(DimX<DimY){
		cout << "Error: dimension of Y cannot be higher than dimension of X." <<endl;
		return 0;
	}else if(DimX==DimY){
		procrustes(X, Y, Xnew, t, rho, A, b, ps);
		return 0;
	}else{
		mat Z = zeros<mat>(NUM, DimX-DimY);
		for(i=0; i<iter; i++){
			mat W = Y;
			W.insert_cols(DimY, Z);
			double tt;
			procrustes(X, W, Xnew, tt, rho, A, b, ps);
			mat Znew = Xnew;
			Znew.shed_cols(0,(DimY-1));
			mat Zm = mean(Znew);
			mat Zc = Znew-repmat(Zm, NUM, 1);
			mat Zd = Znew-Z;
			epsilon = trace(Zd.t()*Zd)/trace(Zc.t()*Zc);
			if(epsilon<=eps){
				break;
			}else{
				Z = Znew;
			}
		}	
		mat Xnew2;
		mat A2;
		rowvec b2;
		double rho2;
		mat X2 = Xnew;
		X2.shed_cols(DimY, (DimX-1));
		procrustes(X2, Y, Xnew2, t, rho2, A2, b2, ps);
		return epsilon; 
	}
}
//#########################     Normalization      ##########################
int normalize(fmat &G, fmat &Gm, fmat &Gsd){
	int i=0;
	int j=0;
	int N = G.n_rows;
	int L = G.n_cols;
	Gm = mean(G,0);
	for(j=0; j<L; j++){
		fvec Gj = G.col(j);
		uvec mis = find(Gj==-9);  //find missing elements
		int M = mis.n_elem;
		if(M>0){		
			Gm(j) = (Gm(j)*N+9*M)/(N-M);
			for(i=0; i<M; i++){
				G(mis(i),j) = Gm(j);
			}
		}
	}
	Gsd = stddev(G,0);
 	for(int j=0; j<L; j++){
		if(Gsd(j)==0){    // Monophmorphic sites are set to 0
			G.col(j) = zeros<fvec>(N);
		}else{
			G.col(j) = (G.col(j)-Gm(j))/Gsd(j);
		}
	}
	return 1;
}
//#########################     PCA      ##########################
int pca(fmat G, int nPCs, mat &PC, rowvec &PCvar, mat &M){
	int i=0;
	int j=0;
	int N = G.n_rows;
	int L = G.n_cols;
	// Use eigen decomposition to get PCA results
	fmat Gm(L,1);
	fmat Gsd(L,1);
	//time_t rawtime;
  	//struct tm * timeinfo;
	//time ( &rawtime );
 	//timeinfo = localtime ( &rawtime );
	//cout << "start normalization at: " << asctime (timeinfo);
	normalize(G, Gm, Gsd);
	//time ( &rawtime );
 	//timeinfo = localtime ( &rawtime );
	//cout << "end normalization at: " << asctime (timeinfo);
	//======================
	M = conv_to<mat>::from(G*G.t());
	Gm.clear();
	Gsd.clear();
	//time ( &rawtime );
 	//timeinfo = localtime ( &rawtime );
	//cout << "get GRM at: " << asctime (timeinfo);
	
	//======== testing multi-threading =========
	/*int Lb = 10000; // Every 10000 SNPs as a block
	int Nb = int((L-0.5)/Lb)+1; // Number of blocks
	fmat Mb = zeros<fmat>(N,N);
	#pragma omp parallel for reduction(+:Mb)
	for(i=0; i<Nb; i++){
		cout << i << endl;
		int a = i*Lb;
		int b = a+Lb-1;
		if(b>L-1) b = L-1;
		fmat Gb = G.cols(a,b);  
		fmat Gm;
		fmat Gsd;
		normalize(Gb, Gm, Gsd);
		Mb += Gb*Gb.t();
	}
	M = conv_to<mat>::from(Mb);
	Mb.clear();*/
	//==========================================
	vec eigval;
	mat eigvec;
	bool bflag = eig_sym(eigval, eigvec, M, "dc");	// use "divide & conquer" algorithm	
	if(!bflag){
		eigval.clear();
		eigvec.clear();
		cout << "Error: eigen decomposition in pca() fails." << endl;
		return 0;
	}
	double eigsum = sum(eigval);
	vec propvar = eigval/eigsum*100;	
	PCvar = zeros<rowvec>(nPCs);
	for(j=0; j<nPCs; j++){
		if(eigval(N-1-j)<0){
			PCvar(j) = 0;
			for(i=0; i<N; i++){
				PC(i,j) = 0;
			}
		}else{
			PCvar(j) = propvar(N-1-j);
			for(i=0; i<N; i++){
				PC(i,j) = eigvec(i, N-1-j)*sqrt(eigval(N-1-j));
			}	
		}
	}
	return 1;
}

//#########################     PCA based on SVD     ##########################
int pca_svd(fmat G, int nPCs, mat &PC, rowvec &PCvar, fmat &Gm, fmat &Gsd, fmat &W){
	int i=0;
	int j=0;
	int N = G.n_rows;
	int L = G.n_cols;
	// Use eigen decomposition to get PCA results
	normalize(G, Gm, Gsd);
	mat eigvec;
	vec S;
	mat V;
	bool bflag = svd_econ(eigvec, S, V, conv_to<mat>::from(G));	// use "divide & conquer" algorithm	
	if(!bflag){
		cout << "Error: singular value decomposition in pca_svd() fails." << endl;
		return 0;
	}
	W = conv_to<fmat>::from(V.cols(0,nPCs-1));
	V.clear();
	vec eigval = S%S;
	double eigsum = sum(eigval);
	vec propvar = eigval/eigsum*100;	
	PCvar = zeros<rowvec>(nPCs);
	for(j=0; j<nPCs; j++){
		if(eigval(j)<0){
			PCvar(j) = 0;
			for(i=0; i<N; i++){
				PC(i,j) = 0;
			}
		}else{
			PCvar(j) = propvar(j);
			for(i=0; i<N; i++){
				PC(i,j) = eigvec(i, j)*S(j);
			}	
		}
	}
	return 1;	
}
//################### Simulate sequence reads from genotypes ##########################
int simuseq(Mat<char> &G, urowvec &C, uvec &Loc, double e, fmat &S, gsl_rng *rng){
	// This function simulates sequence reads S from genotypes G
	// Loc is a list of loci to be simulated
	// Coverage C has to be greater than 0;
	// e is the sequencing error rate per read;
	int N = G.n_rows;
	int L = C.n_cols;
	double P[3];
	P[0] = e;
	P[1] = 0.5;
	P[2] = 1-e;
	for(int i=0; i<N; i++){
		for(int j=0; j<L; j++){
			if(G(i,Loc(j))==-9){
				S(i,j) = -9;
			}else{
				if(P[G(i,Loc(j))]>0 && P[G(i,Loc(j))]<1){
					S(i,j) = gsl_ran_binomial(rng, P[G(i,Loc(j))], C(j));
				}else if(P[G(i,Loc(j))]==0){
					S(i,j) = 0;
				}else{
					S(i,j) = float(C(j));
				}
			}
		}
	}
	return 1;
}

int simuseq2(Mat<char> &G, urowvec &C, uvec &Loc, frowvec &Q, fmat &S, gsl_rng *rng){
	// This function simulates sequence reads S from genotypes G
	// Loc is a list of loci to be simulated
	// Coverage C has to be greater than 0;
	// E is the sequencing error rate per read for for each locus;
	int N = G.n_rows;
	int L = C.n_cols;
	double P[3];
	P[1] = 0.5;
	for(int j=0; j<L; j++){
		P[0] = pow(0.1, Q(j)/10);
		P[2] = 1-P[0];
		for(int i=0; i<N; i++){
			if(G(i,Loc(j))==-9){
				S(i,j) = -9;
			}else{
				if(P[G(i,Loc(j))]>0 && P[G(i,Loc(j))]<1){
					S(i,j) = gsl_ran_binomial(rng, P[G(i,Loc(j))], C(j));
				}else if(P[G(i,Loc(j))]==0){
					S(i,j) = 0;
				}else{
					S(i,j) = float(C(j));
				}
			}
		}
	}
	return 1;
}

//################### Check the average coverage per sample and per locus ####################
int check_coverage(int output, int first_ind, int last_ind, uvec cmnS, urowvec &ExLoci, int &Ls, int &Lg){	
	int i=0;
	int j=0;
	int k=0;
	string str;
	ifstream fin;
	ofstream fout;
	fin.open(SEQ_FILE.c_str());
	if(fin.fail()){
		cerr << "Error: cannot find the SEQ_FILE '" << SEQ_FILE << "'." << endl;
		foutLog << "Error: cannot find the SEQ_FILE '" << SEQ_FILE << "'." << endl;         
		return 0;
	}
	for(i=0; i<SEQ_NON_DATA_ROWS; i++){
		getline(fin, str);          // Read non-data rows
	}

	string outfile = OUT_PREFIX;
	outfile.append(".ind.cov");	
	if(output>0){
		fout.open(outfile.c_str());
		if(fout.fail()){
			cerr << "Error: cannot create a file named " << outfile << "." << endl;
			return 0;
		}
		fout << "popID" << "\t" << "indivID" << "\t" << "L1" << "\t"  << "Ci" << endl;
	}

	int L = LOCI-sum(ExLoci);
	uvec idx(L);
	uvec idx2(L);
	vec C_loc = zeros<vec>(L);	// average coverage per marker
	vec Ncov = zeros<vec>(L);   // number of samples with non-zero coverage
	
	i = 0;
	k = 0;
	for(j=0; j<LOCI_S; j++){
		if(j==cmnS(k)){
			if(ExLoci(k)==0){
				idx(i) = j;
				idx2(i) = k;
				i++;
				if(i==L) break;
			}
			k++;
		}
	}
	
	for(i=1; i<=last_ind; i++)
	{
		if(i<first_ind){
			getline(fin, str);
		}else{
			string SeqInfo1;
			string SeqInfo2;
			vec C(LOCI_S);  // Coverage of one sample
			double S;       // Sequence read
			double Q;       // Base quality
			int Lcov=0;      // number of markers with non-zero coverage
			double C_ind=0;  // average coverage per sample		
			fin >> SeqInfo1 >> SeqInfo2;
			for(j=2; j<SEQ_NON_DATA_COLS; j++){
				fin >> str;
			}
			for(j=0; j<LOCI_S; j++){
				fin >> C(j) >> S >> Q;
			}
			for(j=0; j<L; j++){
				if(C(idx(j))>0){
					Lcov++;
					C_ind += C(idx(j));
					Ncov(j)++;
					C_loc(j) += C(idx(j));
				}				
			}
			C_ind = C_ind/L;
			if(output>0){
				fout << SeqInfo1 << "\t" << SeqInfo2 << "\t" << Lcov << "\t" << C_ind << endl;
			}		
		}
	}
	fin.close();
	if(output>0){
		fout.close();
		cout << "Results of the mean coverage per individual are output to '" << outfile << "'." << endl;
		foutLog << "Results of the mean coverage per individual are output to '" << outfile << "'." << endl;
	}
	C_loc = C_loc/(last_ind-first_ind+1);

	outfile = OUT_PREFIX;
	outfile.append(".loc.cov");
	if(output>0){
		fout.open(outfile.c_str());
		if(fout.fail()){
			cerr << "Error: cannot create a file named " << outfile << "." << endl;
			return 0;
		}
		fout << "ID" << "\t" << "N1" << "\t"  << "Cl" << endl;
	}
	fin.open(SEQ_SITE_FILE.c_str());
	if(fin.fail()){
		cerr << "Error: cannot find the SITE_FILE '" << SEQ_SITE_FILE << "'." << endl;
		foutLog << "Error: cannot find the SITE_FILE '" << SEQ_SITE_FILE << "'." << endl;         
		return 0;
	}
	getline(fin, str);
	k = 0;
	for(j=0; j<LOCI_S; j++){
		if(j==idx(k)){
			if(output>0){
				fin >> str >> str >> str;
				fout << str << "\t" << Ncov(k) << "\t" << C_loc(k) << endl;
			}
			if(MAX_COVERAGE>0 && C_loc(k)>MAX_COVERAGE){
				ExLoci(idx2(k)) = 1;
				Lg++;
			}else if(MIN_COVERAGE>0 && C_loc(k)<MIN_COVERAGE){
				ExLoci(idx2(k)) = 1;
				Ls++;
			}
			k++;
		}
		getline(fin, str);
		if(k==L){ break; }
	}
	if(output>0){
		fout.close();
		cout << "Results of the mean coverage per locus are output to '" << outfile << "'." << endl;
		foutLog << "Results of the mean coverage per locus are output to '" << outfile << "'." << endl;
	}
	Ncov.clear();
	C_loc.clear();

	return 1;
}
//################# Function to check the GENO_FILE format  ##################
int check_format_geno(string filename, int inds, int loci){
	string str;
	int nrow = 0;
	int ncol = 0;
	ifstream fin;
	fin.open(filename.c_str());
	if(fin.fail()){
		cerr << "Error: cannot find the GENO_FILE '" << filename << "'." << endl;    
		foutLog << "Error: cannot find the GENO_FILE '" << filename << "'." << endl;   
		return 0;
	}
	//==========================================================		
	while(nrow < GENO_NON_DATA_ROWS){
		getline(fin, str);     // Read in non-data rows
		nrow+=1;
	}
	while(!fin.eof()){
		getline(fin, str);
		if(str.length()>0 && str!=" "){
			nrow+=1;
			ncol=0;
			bool tab=true;    //Previous character is a tab
			for(int i=0; i<str.length(); i++){
				bool missing=false;     
				if(str[i]!='\t' && i==0){        //Read in the first element
					ncol+=1;
					tab=false;
				}else if(str[i]!='\t' && i>0 && tab){
					ncol+=1;
					tab=false;
					if(ncol>GENO_NON_DATA_COLS && (str[i]!='0' && str[i]!='1'&& str[i]!='2')){
						if(i<(str.length()-2)){
							if(str[i]=='-'&&str[i+1]=='9'&&str[i+2]=='\t'){
								missing = true;
							}
						}else if(i==(str.length()-2)){
							if(str[i]=='-'&&str[i+1]=='9'){
								missing = true;
							}
						}
						if(missing == false){
							cerr<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<<") in the GENO_FILE."<<endl;
							foutLog<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<<") in the GENO_FILE."<<endl;
							fin.close();
							return 0;
						}
					}
				}else if(str[i]!='\t' && i>0 && !tab){
					if(ncol>GENO_NON_DATA_COLS){
						if(str[i-1]!='-' || str[i]!='9'){
							cerr<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<<") in the GENO_FILE."<<endl;
							foutLog<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<<") in the GENO_FILE."<<endl;
							fin.close();
							return 0;
						}
					}
				}else if(str[i]=='\t'){
					tab=true;
				}
			}
			if(ncol!=(loci+GENO_NON_DATA_COLS)){
				cerr << "Error: incorrect number of loci in row " << nrow << " in the GENO_FILE." <<endl;
				foutLog << "Error: incorrect number of loci in row "<< nrow << " in the GENO_FILE." <<endl;
				fin.close();
				return 0;
			}
		}
	}
	if(nrow!=(inds+GENO_NON_DATA_ROWS)){
		cerr << "Error: incorrect number of individuals in the GENO_FILE." << endl;
		foutLog << "Error: incorrect number of individuals in the GENO_FILE." << endl;
		fin.close();
		return 0;
	}
	fin.close();
	return 1;
}
//################# Function to check the SEQ_FILE format  ##################
int check_format_seq(string filename, int inds, int loci){
	string str;
	int nrow = 0;
	int ncol = 0;
	ifstream fin;
	fin.open(filename.c_str());
	if(fin.fail()){
		cerr << "Error: cannot find the SEQ_FILE '" << filename << "'." << endl;    
		foutLog << "Error: cannot find the SEQ_FILE '" << filename << "'." << endl;   
		return 0;
	}
	//==========================================================	
	while(nrow < SEQ_NON_DATA_ROWS){
		getline(fin, str);     // Read in non-data rows
		nrow+=1;
	}
	while(!fin.eof()){
		getline(fin, str);
		if(str.length()>0 && str!=" "){
			nrow+=1;
			ncol=0;
			bool tab=true;   //Previous character is a tab
			int space=0;    // Number of spaces found after the last tab 
			for(int i=0; i<str.length(); i++){    
				if(str[i]!='\t' && i==0){        //Read in the first element
					ncol+=1;
					tab=false;
				}else if(str[i]!='\t' && i>0 && tab){
					ncol+=1;
					tab=false;
					if(ncol>SEQ_NON_DATA_COLS && (str[i]<'0' || str[i]>'9')){
						cerr<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<<") in the SEQ_FILE."<<endl;
						foutLog<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<<") in the SEQ_FILE."<<endl;
						fin.close();
						return 0;
					}
				}else if(str[i]!='\t' && i>0 && !tab){
					if(ncol>SEQ_NON_DATA_COLS && (str[i]<'0' || str[i]>'9')){
						if(str[i]==' ' && space<2){
							space+=1;
						}else{
							cerr<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<<") in the SEQ_FILE."<<endl;
							foutLog<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<<") in the SEQ_FILE."<<endl;
							fin.close();
							return 0;
						}
					}
				}else if(str[i]=='\t'){
					if(ncol>SEQ_NON_DATA_COLS && (space<2 || str[i-1]==' ')){
						cerr<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<<") in the SEQ_FILE."<<endl;
						foutLog<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<<") in the SEQ_FILE."<<endl;
						cerr << "Note: the format of SEQ_FILE has been changed in LASER 2.0 to include base quality scores (see manual for details)." << endl;
						foutLog << "Note: the format of SEQ_FILE has been changed in LASER 2.0 to include base quality scores (see manual for details)." << endl;						
						fin.close();
						return 0;
					}else{
						tab=true;
						space=0;
					}
				}
			}
			if(ncol!=(loci+SEQ_NON_DATA_COLS)){
				cerr << "Error: incorrect number of loci in row " << nrow << " in the SEQ_FILE." <<endl;
				foutLog << "Error: incorrect number of loci in row " << nrow << " in the SEQ_FILE." <<endl;
				fin.close();
				return 0;
			}
		}
	}
	if(nrow!=(inds+SEQ_NON_DATA_ROWS)){
		cerr << "Error: incorrect number of individuals in the SEQ_FILE." << endl;
		foutLog << "Error: incorrect number of individuals in the SEQ_FILE." << endl;
		fin.close();
		return 0;
	}
	fin.close();
	return 1;
}
//################# Function to check the COORD_FILE format  ##################
int check_format_coord(string filename, int inds, int npcs){
	string str;
	int nrow = 0;
	int ncol = 0;
	ifstream fin;
	fin.open(filename.c_str());
	if(fin.fail()){
		cerr << "Error: cannot find the file '" << COORD_FILE << "'." << endl;    
		foutLog << "Error: cannot find the file '" << COORD_FILE << "'." << endl;   
		return 0;
	}
	//==========================================================	
	while(nrow < COORD_NON_DATA_ROWS){
		getline(fin, str);     // Read in non-data rows
		nrow+=1;
	}
	while(!fin.eof()){
		getline(fin, str);
		if(str.length()>0 && str!=" "){
			nrow+=1;
			ncol=0;
			bool tab=true;   //Previous character is a tab
			bool dot=false;    // A dot has been found after the last space 
			bool dash=false;
			bool exp=false;
			for(int i=0; i<str.length(); i++){    
				if(str[i]!='\t' && i==0){        //Read in the first element
					ncol+=1;
					tab=false;
				}else if(str[i]!='\t' && i>0 && tab){
					ncol+=1;
					tab=false;
					if(ncol>COORD_NON_DATA_COLS && (str[i]<'0' || str[i]>'9') && str[i]!='-'){
						cerr<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<< ") in the file '" << filename << "'."<<endl;
						foutLog<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<< ") in the file '" << filename << "'."<<endl;
						fin.close();
						return 0;
					}
				}else if(str[i]!='\t' && i>0 && !tab){
					if(ncol>COORD_NON_DATA_COLS && (str[i]<'0' || str[i]>'9')){
						if(str[i]=='.' && dot==false){
							dot=true;
						}else if(str[i]=='-' && (str[i-1]=='e' || str[i-1]=='E') && dash==false){
							dash=true;
						}else if(str[i]=='e' || str[i]=='E' && exp==false){
							exp=true;
						}else{
							cerr<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<< ") in the file '" << filename << "'."<<endl;
							foutLog<<"Error: invalid value in (row "<<nrow<<", column "<<ncol<< ") in the file '" << filename << "'."<<endl;
							fin.close();
							return 0;
						}
					}
				}else if(str[i]=='\t'){
					tab=true;
					dot=false;
					dash=false;
					exp=false;
				}
			}
			if(ncol!=(npcs+COORD_NON_DATA_COLS)){
				cerr << "Error: incorrect number of PCs in row " << nrow << ") in the file '" << filename << "'."<<endl;
				foutLog << "Error: incorrect number of PCs in row " << nrow << ") in the file '" << filename << "'."<<endl;
				fin.close();
				return 0;
			}
		}
	}
	if(nrow!=(inds+COORD_NON_DATA_ROWS)){
		cerr << "Error: incorrect number of individuals in the file '" << filename << "'." << endl;
		foutLog << "Error: incorrect number of individuals in the file '" << filename << "'." << endl;
		fin.close();
		return 0;
	}
	fin.close();
	return 1;
}
//################# Function to create an empty paramfile  ##################
int create_paramfile(string filename){
	ofstream fout;
	fout.open(filename.c_str());
	if(fout.fail()){
		cerr << "Error: cannot create a file named '" << filename << "'." << endl;
		return 0;
	}
	fout << "# This is a parameter file for LASER v2.04." <<endl;
	fout << "# The entire line after a '#' will be ignored." <<endl;
	fout << "\n" << "###----Main Parameters----###" <<endl;
	fout << endl << "GENO_FILE          # File name of the reference genotype data (include path if in a different directory)" <<endl;
	fout << endl << "SEQ_FILE           # File name of the study sequence data (include path if in a different directory)" <<endl;
	fout << endl << "COORD_FILE         # File name of the reference coordinates (include path if in a different directory)" <<endl;
	fout << endl << "OUT_PREFIX         # Prefix of output files (include path if output to a different directory, default \"laser\")" <<endl;
	fout << endl << "DIM                # Number of PCs to compute (must be a positive integer; default 2)" <<endl;
	fout << endl << "DIM_HIGH           # Number of informative PCs for projection (must be a positive integer >= DIM; default 20)" <<endl;	
	fout << endl << "MIN_LOCI           # Minimum number of covered loci in a sample (must be a positive integer; default 100)" <<endl;

	fout << "\n\n" << "###----Advanced Parameters----###" <<endl;
	fout << endl << "MAX_LOCI           # Maximum number of covered loci in a sample to include (must be a positive integer > MIN_LOCI; default 1000000)" <<endl;
	fout << endl << "SEQ_ERR            # Sequencing error rate per base (must be a number between 0 and 1, or -1; default -1)" <<endl;
	fout <<	        "                   # -1: Use the locus-specific error rates provided in the SEQ_FILE (Phred scale)" <<endl; 
	fout <<	        "                   # Otherwise: Use the specified error rate for all loci and individuals" <<endl;	
	fout << endl << "ALPHA              # Significance level to determine informative PCs (must be a number between 0 and 1; default 0.1)" <<endl;	
	fout <<         "                   # This parameter is effective only if DIM_HIGH is undefined or set to 0." <<endl;
	fout << endl << "THRESHOLD          # Convergence criterion of the projection Procrustes analysis (must be a positive number; default 0.000001)" <<endl;
	fout << endl << "FIRST_IND          # Index of the first sample to analyze (must be a positive integer; default 1)" <<endl;
	fout << endl << "LAST_IND           # Index of the last sample to analyze (must be a positive integer; default [last sample in the SEQ_FILE])" <<endl;
	fout << endl << "REPS               # Number of repeated runs in analyzing each sample (must be a positive integer; default 1)" <<endl;
	fout << endl << "OUTPUT_REPS        # Output results from each repeated run (must be 0 or 1; default 0)" <<endl;
	fout <<	        "                   # 0: Only output mean and standardard deviation across repeated runs" <<endl; 
	fout <<         "                   # 1: Also output results from each repeated run" <<endl;
	fout << endl << "CHECK_COVERAGE     # Check the sequencing coverage (must be 0, 1, or 2; default 0)" <<endl; 
	fout <<	        "                   # 0: Do not check the coverage, and proceed to major computation" <<endl; 
	fout <<         "                   # 1: Check the coverage and proceed to major computation" <<endl;
	fout <<         "                   # 2: Check the coverage and stop" <<endl;
	fout << endl << "CHECK_FORMAT       # Check the format of input files (must be 0, 1, 2, 3, 4, 10, 20, 30, or 40; default 10)" <<endl; 
	fout <<	        "                   # 0: Do not check the format of input files, and proceed to major computation" <<endl; 
	fout <<         "                   # 1: Check the format of all files and stop;     10: Proceed after checking all files" <<endl;
	fout <<         "                   # 2: Check the format of GENO_FILE and stop;     20: Proceed after checking GENO_FILE" <<endl;
	fout <<         "                   # 3: Check the format of SEQ_FILE and stop;      30: Proceed after checking SEQ_FILE" <<endl;
	fout <<         "                   # 4: Check the format of COORD_FILE and stop;    40: Proceed after checking COORD_FILE" <<endl;
	fout << endl << "PCA_MODE           # Switch to the PCA mode (must be 0, 1, 2, or 3; default 0)" <<endl; 
	fout <<	        "                   # 0: Perform LASER to estimate ancestry from sequencing data" <<endl; 
	fout <<         "                   # 1: Perform PCA on the reference genotypes (based on EVD) and output genetic relationship matrix" <<endl;
	fout <<         "                   # 2: Perform PCA on the reference genotypes (memory efficient but slow algorithm)" <<endl;
	fout <<         "                   # 3: Perform PCA on the reference genotypes (based on SVD) and output SNP weights/loadings" <<endl;
	fout << endl << "REF_SIZE           # Number of individuals randomly selected as the reference (must be a positive integer; default [sample size in the GENO_FILE])" <<endl;	
	fout << endl << "TRIM_PROP          # Proportion of shared loci to be trimmed off for all samples (must be a number between 0 and 1; default 0)" <<endl;
	fout << endl << "EXCLUDE_LIST       # File name of a list of SNPs to exclude from the analysis (include path if in a different directory)" <<endl;	
	fout << endl << "MIN_COVERAGE       # Minimum mean coverage for a locus to be included in the analysis (must be a non-negative number; default 0)" <<endl; 	
	fout << endl << "MAX_COVERAGE       # Maximum mean coverage for a locus to be included in the analysis (must be a positive number or -1; default -1)" <<endl; 
	fout <<	        "                   # -1: Include all loci in the analysis without removal based on mean coverage" <<endl;
	fout << endl << "PROCRUSTES_SCALE   # Methods to calculate the scaling parameter in Procrustes analysis (must be 0 or 1; default 0)" <<endl; 
	fout <<	        "                   # 0: Calculate the scaling parameter to maximize the Procrustes similarity" <<endl; 
	fout <<         "                   # 1: Fix the scaling parameter to match the variance of two sets of coordinates in Procrustes analysis" <<endl;
	fout << endl << "KNN_ZSCORE         # Number of nearest neigbors used to calculate the Z score for each study individual (must be an integer >2; default 10)" <<endl;
	fout << endl << "RANDOM_SEED        # Seed for the random number generator in the program (must be a non-negative integer; default 0)" <<endl; 
	fout << endl << "NUM_THREADS        # Number of CPU cores for multi-threading parallel analysis (must be a positive integer; default 8)" <<endl; 


 	fout << "\n\n" << "###----Command line arguments----###" <<endl <<endl;
	fout << "# -p     parameterfile (this file)" << endl;
	fout << "# -g     GENO_FILE" <<endl;
	fout << "# -s     SEQ_FILE" <<endl;
	fout << "# -c     COORD_FILE" <<endl;
	fout << "# -o     OUT_PREFIX" <<endl;
	fout << "# -k     DIM" <<endl;
	fout << "# -K     DIM_HIGH" <<endl;	
	fout << "# -l     MIN_LOCI" <<endl;
	fout << "# -L     MAX_LOCI" <<endl;
	fout << "# -e     SEQ_ERR" <<endl;
	fout << "# -a     ALPHA" << endl;
	fout << "# -t     THRESHOLD" <<endl;	
	fout << "# -x     FIRST_IND" <<endl;
	fout << "# -y     LAST_IND" <<endl;
	fout << "# -r     REPS" <<endl;
	fout << "# -R     OUTPUT_REPS" <<endl;
	fout << "# -cov   CHECK_COVERAGE" <<endl;
	fout << "# -fmt   CHECK_FORMAT" <<endl;	
	fout << "# -pca   PCA_MODE" <<endl;
	fout << "# -N     REF_SIZE" <<endl;
	fout << "# -M     TRIM_PROP" <<endl;
	fout << "# -ex    EXCLUDE_LIST" <<endl;	
	fout << "# -minc  MAX_COVERAGE" << endl;	
	fout << "# -maxc  MAX_COVERAGE" << endl;
	fout << "# -rho   PROCRUSTES_SCALE" << endl;
	fout << "# -knn   KNN_ZSCORE" << endl;
	fout << "# -seed  RANDOM_SEED" << endl;
	fout << "# -nt    NUM_THREADS" << endl;

	fout << "\n" << "###----end of file----###";
	fout.close();
	cout << "An empty template parameter file named '"<< filename << "' has been created." << endl;
	foutLog << "An empty template parameter file named '"<< filename << "' has been created." << endl;  
	return 1;
}		
//################# Function to read and check the paramfile  ##################
int read_paramfile(string filename){	
	int flag = 1;
	ifstream fin;
	fin.open(filename.c_str());
	if(fin.fail()){
		cerr << "Warning: cannot find the PARAM_FILE '" << filename << "'." << endl;
		foutLog << "Warning: cannot find the PARAM_FILE '" << filename << "'." << endl;
		create_paramfile(filename);
		return flag;
	}
	string str;
	while(!fin.eof()){
		fin>>str;
		if(str[0]=='#'){
			getline(fin, str);
		}else if(str.compare("GENO_FILE")==0){
			fin>>str;
			if(str[0]!='#'){
				if(GENO_FILE == default_str){
					GENO_FILE = str;
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("COORD_FILE")==0){
			fin>>str;
			if(str[0]!='#'){
				if(COORD_FILE == default_str){
					COORD_FILE = str;
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("SEQ_FILE")==0){
			fin>>str;
			if(str[0]!='#'){
				if(SEQ_FILE == default_str){
					SEQ_FILE = str;
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("EXCLUDE_LIST")==0){
			fin>>str;
			if(str[0]!='#'){
				if(EXCLUDE_LIST == default_str){
					EXCLUDE_LIST = str;
				}
			}else{
				getline(fin, str);
			}			
		}else if(str.compare("OUT_PREFIX")==0){
			fin>>str;
			if(str[0]!='#'){
				if(OUT_PREFIX == default_str){
					OUT_PREFIX = str;
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("MIN_LOCI")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_int(str) && atoi(str.c_str())>0){
					if(MIN_LOCI == default_int){
						MIN_LOCI = atoi(str.c_str());
					}
				}else{
					if(MIN_LOCI != default_int){
						cerr<< "Warning: MIN_LOCI in the parameter file is not a positive integer." <<endl;
						foutLog<< "Warning: MIN_LOCI in the parameter file is not a positive integer." <<endl;
					}else{
						MIN_LOCI = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("MAX_LOCI")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_int(str) && atoi(str.c_str())>0){
					if(MAX_LOCI == default_int){
						MAX_LOCI = atoi(str.c_str());
					}
				}else{
					if(MAX_LOCI != default_int){
						cerr<< "Warning: MAX_LOCI in the parameter file is not a positive integer." <<endl;
						foutLog<< "Warning: MAX_LOCI in the parameter file is not a positive integer." <<endl;
					}else{
						MAX_LOCI = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}			
		}else if(str.compare("DIM")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_int(str) && atoi(str.c_str())>0){
					if(DIM == default_int){
						DIM = atoi(str.c_str());
					}
				}else{
					if(DIM != default_int){
						cerr<< "Warning: DIM in the parameter file is not a positive integer." <<endl;
						foutLog<< "Warning: DIM in the parameter file is not a positive integer." <<endl;
					}else{
						DIM = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("DIM_HIGH")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_int(str) && atoi(str.c_str())>=0){
					if(DIM_HIGH == default_int){
						DIM_HIGH = atoi(str.c_str());
					}
				}else{
					if(DIM_HIGH != default_int){
						cerr<< "Warning: DIM_HIGH in the parameter file is not a positive integer." <<endl;
						foutLog<< "Warning: DIM_HIGH in the parameter file is not a positive integer." <<endl;
					}else{
						DIM_HIGH = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}				
		}else if(str.compare("SEQ_ERR")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_numeric(str) && (atof(str.c_str())>=0 && atof(str.c_str())<=1) || atof(str.c_str())==-1){
					if(SEQ_ERR == default_double){
						SEQ_ERR = atof(str.c_str());
					}
				}else{
					if(SEQ_ERR != default_double){
						cerr<< "Warning: SEQ_ERR in the parameter file is not between 0 and 1." <<endl;
						foutLog<< "Warning: SEQ_ERR in the parameter file is not between 0 and 1." <<endl;
					}else{
						SEQ_ERR  = default_double-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("ALPHA")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_numeric(str) && atof(str.c_str())>=0 && atof(str.c_str())<=1){
					if(ALPHA == default_double){
						ALPHA = atof(str.c_str());
					}
				}else{
					if(ALPHA != default_double){
						cerr<< "Warning: ALPHA in the parameter file is not between 0 and 1." <<endl;
						foutLog<< "Warning: ALPHA in the parameter file is not between 0 and 1." <<endl;
					}else{
						ALPHA  = default_double-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("THRESHOLD")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_numeric(str) && atof(str.c_str())>=0){
					if(THRESHOLD == default_double){
						THRESHOLD = atof(str.c_str());
					}
				}else{
					if(THRESHOLD != default_double){
						cerr<< "Warning: THRESHOLD in the parameter file is not a positive number." <<endl;
						foutLog<< "Warning: THRESHOLD in the parameter file is not a positive number." <<endl;
					}else{
						THRESHOLD  = default_double-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("REF_SIZE")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_int(str) && atoi(str.c_str())>=0){
					if(REF_SIZE == default_int){
						REF_SIZE = atoi(str.c_str());
					}
				}else{
					if(REF_SIZE != default_int){
						cerr<< "Warning: REF_SIZE in the parameter file is not a positive integer." <<endl;
						foutLog<< "Warning: REF_SIZE in the parameter file is not a positive integer." <<endl;
					}else{
						REF_SIZE = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}			
		}else if(str.compare("FIRST_IND")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_int(str) && atoi(str.c_str())>0){
					if(FIRST_IND == default_int){
						FIRST_IND = atoi(str.c_str());
					}
				}else{
					if(FIRST_IND != default_int){
						cerr<< "Warning: FIRST_IND in the parameter file is not a positive integer." <<endl;
						foutLog<< "Warning: FIRST_IND in the parameter file is not a positive integer." <<endl;
					}else{
						 FIRST_IND = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("LAST_IND")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_int(str) && atoi(str.c_str())>0){
					if(LAST_IND == default_int){
						LAST_IND = atoi(str.c_str());
					}
				}else{
					if(LAST_IND != default_int){
						cerr<< "Warning: LAST_IND in the parameter file is not a positive integer." <<endl;
						foutLog<< "Warning: LAST_IND in the parameter file is not a positive integer." <<endl;
					}else{
						 LAST_IND = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("REPS")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_int(str) && atoi(str.c_str())>0){
					if(REPS == default_int){
						REPS = atoi(str.c_str());
					}
				}else{
					if(REPS != default_int){
						cerr<< "Warning: REPS in the parameter file is not a positive integer." <<endl;
						foutLog<< "Warning: REPS in the parameter file is not a positive integer." <<endl;
					}else{
						 REPS = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("OUTPUT_REPS")==0){
			fin>>str;
			if(str[0]!='#'){
				if(str=="0" || str=="1"){
					if(OUTPUT_REPS == default_int){
						OUTPUT_REPS = atoi(str.c_str());
					}
				}else{
					if(OUTPUT_REPS != default_int){
						cerr<< "Warning: OUTPUT_REPS in the parameter file is not 0 or 1." <<endl;
						foutLog<< "Warning: OUTPUT_REPS in the parameter file is not 0 or 1." <<endl;
					}else{
						OUTPUT_REPS  = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("CHECK_FORMAT")==0){
			fin>>str;
			if(str[0]!='#'){
				if(str=="0" || str=="1" || str=="2" || str=="3" || str=="4" || str=="10" || str=="20" || str=="30" || str=="40"){
					if(CHECK_FORMAT == default_int){
						CHECK_FORMAT = atoi(str.c_str());
					}
				}else{
					if(CHECK_FORMAT != default_int){
						cerr<< "Warning: CHECK_FORMAT in the parameter file is not 0, 1, 2, 3, 4, 10, 20, 30, or 40." <<endl;
						foutLog<< "Warning: CHECK_FORMAT in the parameter file is not 0, 1, 2, 3, 4, 10, 20, 30, or 40." <<endl;
					}else{
						CHECK_FORMAT  = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("TRIM_PROP")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_numeric(str) && atof(str.c_str())>=0 && atof(str.c_str())<=1){
					if(TRIM_PROP == default_double){
						TRIM_PROP = atof(str.c_str());
					}
				}else{
					if(TRIM_PROP != default_double){
						cerr<< "Warning: TRIM_PROP in the parameter file is not between 0 and 1." <<endl;
						foutLog<< "Warning: TRIM_PROP in the parameter file is not between 0 and 1." <<endl;
					}else{
						TRIM_PROP  = default_double-1;
					}
				}
			}else{
				getline(fin, str);
			}			
		}else if(str.compare("MIN_COVERAGE")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_numeric(str)){
					if(MIN_COVERAGE == default_double){
						MIN_COVERAGE = atof(str.c_str());
					}
				}else{
					if(MIN_COVERAGE != default_double){
						cerr<< "Warning: MIN_COVERAGE in the parameter file is not a non-negative number." <<endl;
						foutLog<< "Warning: MIN_COVERAGE in the parameter file is not a non-negative number." <<endl;
					}else{
						MIN_COVERAGE  = default_double-1;
					}
				}
			}else{
				getline(fin, str);
			}			
		}else if(str.compare("MAX_COVERAGE")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_numeric(str) && (atof(str.c_str())>0 || atof(str.c_str())==-1)){
					if(MAX_COVERAGE == default_double){
						MAX_COVERAGE = atof(str.c_str());
					}
				}else{
					if(MAX_COVERAGE != default_double){
						cerr<< "Warning: MAX_COVERAGE in the parameter file is not a positive number or -1." <<endl;
						foutLog<< "Warning: MAX_COVERAGE in the parameter file is not a positive number or -1." <<endl;
					}else{
						MAX_COVERAGE  = default_double-1;
					}
				}
			}else{
				getline(fin, str);
			}			
		}else if(str.compare("CHECK_COVERAGE")==0){
			fin>>str;
			if(str[0]!='#'){
				if(str=="0" || str=="1" || str=="2"){
					if(CHECK_COVERAGE == default_int){
						CHECK_COVERAGE = atoi(str.c_str());
					}
				}else{
					if(CHECK_COVERAGE != default_int){
						cerr<< "Warning: CHECK_COVERAGE in the parameter file is not 0, 1, or 2." <<endl;
						foutLog<< "Warning: CHECK_COVERAGE in the parameter file is not 0, 1, or 2." <<endl;
					}else{
						CHECK_COVERAGE  = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("PCA_MODE")==0){
			fin>>str;
			if(str[0]!='#'){
				if(str=="0" || str=="1"  || str=="2"  || str=="3"){
					if(PCA_MODE == default_int){
						PCA_MODE = atoi(str.c_str());
					}
				}else{
					if(PCA_MODE != default_int){
						cerr<< "Warning: PCA_MODE in the parameter file is not 0, 1, 2, or 3." <<endl;
						foutLog<< "Warning: PCA_MODE in the parameter file is not 0, 1, 2, or 3." <<endl;
					}else{
						PCA_MODE  = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("PROCRUSTES_SCALE")==0){
			fin>>str;
			if(str[0]!='#'){
				if(str=="0" || str=="1"){
					if(PROCRUSTES_SCALE == default_int){
						PROCRUSTES_SCALE = atoi(str.c_str());
					}
				}else{
					if(PROCRUSTES_SCALE != default_int){
						cerr<< "Warning: PROCRUSTES_SCALE in the parameter file is not 0 or 1." <<endl;
						foutLog<< "Warning: PROCRUSTES_SCALE in the parameter file is not 0 or 1." <<endl;
					}else{
						PROCRUSTES_SCALE  = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("RANDOM_SEED")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_int(str) && atoi(str.c_str())>=0){
					if(RANDOM_SEED == default_int){
						RANDOM_SEED = atoi(str.c_str());
					}
				}else{
					if(RANDOM_SEED != default_int){
						cerr<< "Warning: RANDOM_SEED in the parameter file is not a non-negative integer." <<endl;
						foutLog<< "Warning: RANDOM_SEED in the parameter file is not a non-negative integer." <<endl;
					}else{
						RANDOM_SEED = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}
		}else if(str.compare("KNN_ZSCORE")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_int(str) && atoi(str.c_str())>2){
					if(KNN_ZSCORE == default_int){
						KNN_ZSCORE = atoi(str.c_str());
					}
				}else{
					if(KNN_ZSCORE != default_int){
						cerr<< "Warning: KNN_ZSCORE in the parameter file is not an integer >2." <<endl;
						foutLog<< "Warning: KNN_ZSCORE in the parameter file is not an integer >2." <<endl;
					}else{
						 KNN_ZSCORE = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}			
		}else if(str.compare("NUM_THREADS")==0){
			fin>>str;
			if(str[0]!='#'){
				if(is_int(str) && atoi(str.c_str())>0){
					if(NUM_THREADS == default_int){
						NUM_THREADS = atoi(str.c_str());
					}
				}else{
					if(NUM_THREADS != default_int){
						cerr<< "Warning: NUM_THREADS in the parameter file is not a positive integer." <<endl;
						foutLog<< "Warning: NUM_THREADS in the parameter file is not a positive integer." <<endl;
					}else{
						NUM_THREADS = default_int-1;
					}
				}
			}else{
				getline(fin, str);
			}	
		}		
	}
	fin.close();
	return flag;
}
//################# Function to print parameters in execution  ##################
void print_configuration(){
	cout <<endl << "Parameter values used in execution:" <<endl;
	cout << "-------------------------------------------------" << endl;
	if(PCA_MODE == 0){
		cout << "SEQ_FILE (-s)" << "\t" << SEQ_FILE << endl;
		cout << "GENO_FILE (-g)" << "\t" << GENO_FILE <<endl;
		if(COORD_FILE.compare(default_str)!=0){
			cout << "COORD_FILE (-c)" << "\t" << COORD_FILE << endl;
		}
		cout << "OUT_PREFIX (-o)" << "\t" << OUT_PREFIX << endl;		
		cout << "DIM (-k)" << "\t" << DIM << endl;
		if(DIM_HIGH != 0){
			cout << "DIM_HIGH (-K)" << "\t" << DIM_HIGH << endl;
		}else{
			cout << "ALPHA (-a)" << "\t" << ALPHA << endl;
		}
		cout << "THRESHOLD (-t)" << "\t" << THRESHOLD << endl;			
		cout << "MIN_LOCI (-l)" << "\t" << MIN_LOCI << endl;
		cout << "MAX_LOCI (-L)" << "\t" << MAX_LOCI << endl;
		cout << "SEQ_ERR (-e)" << "\t" << SEQ_ERR << endl;
		cout << "FIRST_IND (-x)" << "\t" << FIRST_IND << endl;
		cout << "LAST_IND (-y)" << "\t" << LAST_IND << endl;
		cout << "REPS (-r)" << "\t" << REPS << endl;	
		cout << "OUTPUT_REPS (-R)" << "\t" << OUTPUT_REPS << endl;			
		cout << "CHECK_COVERAGE (-cov)" << "\t" << CHECK_COVERAGE <<endl;
		cout << "CHECK_FORMAT (-fmt)" << "\t" << CHECK_FORMAT << endl; 	 
		cout << "PCA_MODE (-pca)" << "\t" << PCA_MODE << endl;
		if(REF_SIZE != REF_INDS){
			cout << "REF_SIZE (-N)" << "\t" << REF_SIZE << endl;
		}
		if(TRIM_PROP>0){
			cout << "TRIM_PROP (-M)" << "\t" << TRIM_PROP << endl;
		}
		if(EXCLUDE_LIST.compare(default_str)!=0){
			cout << "EXCLUDE_LIST (-ex)" << "\t" << EXCLUDE_LIST << endl;
		}
		if(MIN_COVERAGE>0){
			cout << "MIN_COVERAGE (-minc)" << "\t" << MIN_COVERAGE <<endl;
		}
		if(MAX_COVERAGE>0){
			cout << "MAX_COVERAGE (-maxc)" << "\t" << MAX_COVERAGE <<endl;
		}
		if(PROCRUSTES_SCALE>0){
			cout << "PROCRUSTES_SCALE (-rho)" << "\t" << PROCRUSTES_SCALE << endl;
		}
		cout << "KNN_ZSCORE (-knn)" << "\t" << KNN_ZSCORE << endl;
		cout << "RANDOM_SEED (-seed)" << "\t" << RANDOM_SEED << endl;
		cout << "NUM_THREADS (-nt)" << "\t" << NUM_THREADS << endl;
	}else{
		cout << "GENO_FILE (-g)" << "\t" << GENO_FILE <<endl;
		cout << "DIM (-k)" << "\t" << DIM << endl;
		cout << "OUT_PREFIX (-o)" << "\t" << OUT_PREFIX << endl;
		cout << "CHECK_FORMAT (-fmt)" << "\t" << CHECK_FORMAT << endl; 
		cout << "PCA_MODE (-pca)" << "\t" << PCA_MODE << endl;
		if(REF_SIZE != REF_INDS){
			cout << "REF_SIZE (-N)" << "\t" << REF_SIZE << endl;
		}		
		if(TRIM_PROP>0){
			cout << "TRIM_PROP (-M)" << "\t" << TRIM_PROP << endl;
		}
		if(EXCLUDE_LIST.compare(default_str)!=0){
			cout << "EXCLUDE_LIST (-ex)" << "\t" << EXCLUDE_LIST << endl;
		}
		cout << "RANDOM_SEED (-seed)" << "\t" << RANDOM_SEED << endl;
		cout << "NUM_THREADS (-nt)" << "\t" << NUM_THREADS << endl;
	}
	cout << "-------------------------------------------------" << endl; 

	foutLog <<endl << "Parameter values used in execution:" <<endl;
	foutLog << "-------------------------------------------------" << endl;
	if(PCA_MODE == 0){
		foutLog << "GENO_FILE (-g)" << "\t" << GENO_FILE <<endl;
		foutLog << "SEQ_FILE (-s)" << "\t" << SEQ_FILE << endl;
		if(COORD_FILE.compare(default_str)!=0){
			foutLog << "COORD_FILE (-c)" << "\t" << COORD_FILE << endl;
		}
		foutLog << "OUT_PREFIX (-o)" << "\t" << OUT_PREFIX << endl;
		foutLog << "DIM (-k)" << "\t" << DIM << endl;
		if(DIM_HIGH != 0){
			foutLog << "DIM_HIGH (-K)" << "\t" << DIM_HIGH << endl;
		}else{
			foutLog << "ALPHA (-a)" << "\t" << ALPHA << endl;
		}
		foutLog << "THRESHOLD (-t)" << "\t" << THRESHOLD << endl;			
		foutLog << "MIN_LOCI (-l)" << "\t" << MIN_LOCI << endl;
		foutLog << "MAX_LOCI (-L)" << "\t" << MAX_LOCI << endl;
		foutLog << "SEQ_ERR (-e)" << "\t" << SEQ_ERR << endl;
		foutLog << "FIRST_IND (-x)" << "\t" << FIRST_IND << endl;
		foutLog << "LAST_IND (-y)" << "\t" << LAST_IND << endl;
		foutLog << "REPS (-r)" << "\t" << REPS << endl;	
		foutLog << "OUTPUT_REPS (-R)" << "\t" << OUTPUT_REPS << endl;							
		foutLog << "CHECK_COVERAGE (-cov)" << "\t" << CHECK_COVERAGE <<endl;
		foutLog << "CHECK_FORMAT (-fmt)" << "\t" << CHECK_FORMAT << endl; 	 
		foutLog << "PCA_MODE (-pca)" << "\t" << PCA_MODE << endl;
		if(REF_SIZE != REF_INDS){
			foutLog << "REF_SIZE (-N)" << "\t" << REF_SIZE << endl;
		}		
		if(TRIM_PROP>0){
			foutLog << "TRIM_PROP (-M)" << "\t" << TRIM_PROP << endl;
		}
		if(EXCLUDE_LIST.compare(default_str)!=0){
			foutLog << "EXCLUDE_LIST (-ex)" << "\t" << EXCLUDE_LIST << endl;
		}		
		if(MIN_COVERAGE>0){
			foutLog << "MIN_COVERAGE (-minc)" << "\t" << MIN_COVERAGE <<endl;
		}
		if(MAX_COVERAGE>0){
			foutLog << "MAX_COVERAGE (-maxc)" << "\t" << MAX_COVERAGE <<endl;
		}
		if(PROCRUSTES_SCALE>0){
			foutLog << "PROCRUSTES_SCALE (-rho)" << "\t" << PROCRUSTES_SCALE << endl;
		}
		foutLog << "KNN_ZSCORE (-knn)" << "\t" << KNN_ZSCORE << endl;
		foutLog << "RANDOM_SEED (-seed)" << "\t" << RANDOM_SEED << endl;
		foutLog << "NUM_THREADS (-nt)" << "\t" << NUM_THREADS << endl;
	}else{
		foutLog << "GENO_FILE (-g)" << "\t" << GENO_FILE <<endl;
		foutLog << "DIM (-k)" << "\t" << DIM << endl;
		foutLog << "OUT_PREFIX (-o)" << "\t" << OUT_PREFIX << endl;
		foutLog << "CHECK_FORMAT (-fmt)" << "\t" << CHECK_FORMAT << endl; 
		foutLog << "PCA_MODE (-pca)" << "\t" << PCA_MODE << endl;
		if(REF_SIZE != REF_INDS){
			foutLog << "REF_SIZE (-N)" << "\t" << REF_SIZE << endl;
		}		
		if(TRIM_PROP>0){
			foutLog << "TRIM_PROP (-M)" << "\t" << TRIM_PROP << endl;
		}
		if(EXCLUDE_LIST.compare(default_str)!=0){
			foutLog << "EXCLUDE_LIST (-ex)" << "\t" << EXCLUDE_LIST << endl;
		}
		foutLog << "RANDOM_SEED (-seed)" << "\t" << RANDOM_SEED << endl;
		foutLog << "NUM_THREADS (-nt)" << "\t" << NUM_THREADS << endl;		
	}
	foutLog << "-------------------------------------------------" << endl; 
}
//################# Function to check parameter values  ##################
int check_parameters(){
	int flag = 1;
	if(GENO_FILE.compare(default_str)==0){
		cerr << "Error: GENO_FILE (-g) is not specified." << endl;
		foutLog << "Error: GENO_FILE (-g) is not specified." << endl;
		flag = 0;
	}
	if(SEQ_FILE.compare(default_str)==0 && PCA_MODE==0){
		cerr << "Error: SEQ_FILE (-s) is not specified." << endl;
		foutLog << "Error: SEQ_FILE (-s) is not specified." << endl;
		flag = 0;
	}
	if(DIM==default_int){
		cerr << "Error: DIM (-k) is not specified." << endl;
		foutLog << "Error: DIM (-k) is not specified." << endl;
		flag = 0;
	}else if(DIM<1){ 
		cerr << "Error: invalid value for DIM (-k)." << endl;
		foutLog << "Error: invalid value for DIM (-k)." << endl;
		flag = 0;
	}else if(REF_SIZE!=default_int && DIM>=REF_SIZE){
		cerr << "Error: invalid value for DIM (-k)." << endl;
		cerr << "DIM must be smaller than REF_SIZE;" << endl;
		foutLog << "Error: invalid value for DIM (-k)." << endl;
		foutLog << "DIM must be smaller than REF_SIZE;" << endl;
		flag = 0;		
	}else if(LOCI!=default_int && DIM>=LOCI){
		cerr << "Error: invalid value for DIM (-k)." << endl;
		cerr << "DIM must be smaller than the number of loci in the GENO_FILE;" << endl;
		foutLog << "Error: invalid value for DIM (-k)." << endl;
		foutLog << "DIM must be smaller than the number of loci in the GENO_FILE;" << endl;
		flag = 0;
	}else if(NUM_PCS!=default_int && DIM>NUM_PCS && PCA_MODE==0){
		cerr << "Error: invalid value for DIM (-k)." << endl;
		cerr << "DIM cannot be greater than the number of PCs in the COORD_FILE;" << endl;
		foutLog << "Error: invalid value for DIM (-k)." << endl;
		foutLog << "DIM cannot be greater than the number of PCs in the COORD_FILE;" << endl;
		flag = 0;
	}
	if(DIM_HIGH==default_int && PCA_MODE==0){
		cerr << "Error: DIM_HIGH (-K) is not specified." << endl;
		foutLog << "Error: DIM_HIGH (-K) is not specified." << endl;
		flag = 0;
	}else if(DIM_HIGH<DIM && DIM_HIGH!=0 && PCA_MODE==0){ 
		cerr << "Error: invalid value for DIM_HIGH (-K)." << endl;
		cerr << "DIM_HIGH cannot be smaller than DIM." << endl;
		foutLog << "Error: invalid value for DIM_HIGH (-K)." << endl;
		foutLog << "DIM_HIGH cannot be smaller than DIM." << endl;
		flag = 0;
	}else if(REF_SIZE!=default_int && DIM_HIGH>=REF_SIZE && PCA_MODE==0){
		cerr << "Error: invalid value for DIM_HIGH (-K)." << endl;
		cerr << "DIM_HIGH must be smaller than REF_SIZE." << endl;
		foutLog << "Error: invalid value for DIM_HIGH (-K)." << endl;
		foutLog << "DIM_HIGH must be smaller than REF_SIZE." << endl;
		flag = 0;
	}else if(LOCI!=default_int && DIM_HIGH>=LOCI && PCA_MODE==0){
		cerr << "Error: invalid value for DIM_HIGH (-K)." << endl;
		cerr << "DIM_HIGH must be smaller than the total number of shared loci." << endl;
		foutLog << "Error: invalid value for DIM_HIGH (-K)." << endl;
		foutLog << "DIM_HIGH must be smaller than the total number of shared loci." << endl;
		flag = 0;		
	}		
	if(MIN_LOCI < 1){
		cerr << "Error: invalid value for MIN_LOCI (-l)." << endl;
		foutLog << "Error: invalid value for MIN_LOCI (-l)." << endl;
		flag = 0;
	}else if(MIN_LOCI>LOCI  && LOCI!=default_int){
		cerr << "Error: invalid value for MIN_LOCI (-l)." << endl;
		cerr << "MIN_LOCI cannot be greater than the total number of loci." << endl;
		foutLog << "Error: invalid value for MIN_LOCI (-l)." << endl;
		foutLog << "MIN_LOCI cannot be greater than the total number of loci." << endl;
		flag = 0;
	}
	if(MAX_LOCI < MIN_LOCI){
		cerr << "Error: invalid value for MAX_LOCI (-L)." << endl;
		cerr << "MAX_LOCI must be greater than MIN_LOCI." << endl;
		foutLog << "Error: invalid value for MAX_LOCI (-L)." << endl;
		foutLog << "MAX_LOCI must be greater than MIN_LOCI." << endl;
		flag = 0;
	}	
	if((SEQ_ERR < 0 || SEQ_ERR > 1) && SEQ_ERR != -1){
		cerr << "Error: invalid value for SEQ_ERR (-e)." << endl;
		foutLog << "Error: invalid value for SEQ_ERR (-e)." << endl;
		flag = 0;
	}
	if(REF_SIZE < 0){
		cerr << "Error: invalid value for REF_SIZE (-N)." << endl;
		foutLog << "Error: invalid value for REF_SIZE (-N)." << endl;
		flag = 0;	
	}else if(REF_SIZE<=DIM_HIGH){
		cerr << "Error: invalid value for REF_SIZE (-N)." << endl;
		cerr << "REF_SIZE must be greater than DIM_HIGH." << endl;
		foutLog << "Error: invalid value for REF_SIZE (-N)." << endl;
		foutLog << "REF_SIZE must be greater than DIM_HIGH." << endl;
		flag = 0;
	}else if(REF_SIZE>REF_INDS){
		cerr << "Error: invalid value for REF_SIZE (-N)." << endl;
		cerr << "REF_SIZE cannot be greater than the number of individuals in the GENO_FILE." << endl;
		foutLog << "Error: invalid value for REF_SIZE (-N)." << endl;
		foutLog << "REF_SIZE cannot be greater than the number of individuals in the GENO_FILE." << endl;
		flag = 0;	
	}
	if(FIRST_IND < 0){
		cerr << "Error: invalid value for FIRST_IND (-x)." << endl;
		foutLog << "Error: invalid value for FIRST_IND (-x)." << endl;
		flag = 0;
	}else if(FIRST_IND>SEQ_INDS && SEQ_INDS!=default_int){
		cerr << "Error: invalid value for FIRST_IND (-x)." << endl;
		cerr << "FIRST_IND cannot be greater than the number of individuals in the SEQ_FILE." << endl;
		foutLog << "Error: invalid value for FIRST_IND (-x)." << endl;
		foutLog << "FIRST_IND cannot be greater than the number of individuals in the SEQ_FILE." << endl;
		flag = 0;
	}
	if(LAST_IND<0 && LAST_IND!=default_int){
		cerr << "Error: invalid value for LAST_IND (-y)." << endl;
		foutLog << "Error: invalid value for LAST_IND (-y)." << endl;
		flag = 0;	
	}else if(LAST_IND<FIRST_IND && FIRST_IND!=default_int && LAST_IND!=default_int){
		cerr << "Error: invalid value for LAST_IND (-y)." << endl;
		cerr << "LAST_IND cannot be smaller than FIRST_IND." << endl;
		foutLog << "Error: invalid value for LAST_IND (-y)." << endl;
		foutLog << "LAST_IND cannot be smaller than FIRST_IND." << endl;
		flag = 0;
	}else if(LAST_IND>SEQ_INDS && SEQ_INDS!=default_int){
		LAST_IND = SEQ_INDS;
	}
	if(REPS < 1){
		cerr << "Error: invalid value for REPS (-r)." << endl;
		foutLog << "Error: invalid value for REPS (-r)." << endl;
		flag = 0;
	}
	if(OUTPUT_REPS!=0 && OUTPUT_REPS!=1){
		cerr << "Error: invalid value for OUTPUT_REPS (-R)." << endl;
		foutLog << "Error: invalid value for OUTPUT_REPS (-R)." << endl;
		flag = 0;
	}
	if(CHECK_COVERAGE!=0 && CHECK_COVERAGE!=1 && CHECK_COVERAGE!=2){
		cerr << "Error: invalid value for CHECK_COVERAGE (-cov)." << endl;
		foutLog << "Error: invalid value for CHECK_COVERAGE (-cov)." << endl;
		flag = 0;
	}
	if(CHECK_FORMAT!=0 && CHECK_FORMAT!=1 && CHECK_FORMAT!=2 && CHECK_FORMAT!=3 && CHECK_FORMAT!=4){
	 	if(CHECK_FORMAT!=10 && CHECK_FORMAT!=20 && CHECK_FORMAT!=30 && CHECK_FORMAT!=40){
			cerr << "Error: invalid value for CHECK_FORMAT (-fmt)." << endl;
			foutLog << "Error: invalid value for CHECK_FORMAT (-fmt)." << endl;
			flag = 0;
		}
	}
	if(PCA_MODE!=0 && PCA_MODE!=1 && PCA_MODE!=2 && PCA_MODE!=3){
		cerr << "Error: invalid value for PCA_MODE (-pca)." << endl;
		foutLog << "Error: invalid value for PCA_MODE (-pca)." << endl;
		flag = 0;
	}	
	if(TRIM_PROP < 0 || TRIM_PROP >= 1){
		cerr << "Error: invalid value for TRIM_PROP (-M)." << endl;
		foutLog << "Error: invalid value for TRIM_PROP (-M)." << endl;
		foutLog << "Error: invalid value for TRIM_PROP (-M)." << endl;
		flag = 0;
	}	
	if(MIN_COVERAGE<0){
		cerr << "Error: invalid value for MIN_COVERAGE (-minc)." << endl;
		foutLog << "Error: invalid value for MIN_COVERAGE (-minc)." << endl;
		flag = 0;
	}	
	if(MAX_COVERAGE<=0 && MAX_COVERAGE!=-1){
		cerr << "Error: invalid value for MAX_COVERAGE (-maxc)." << endl;
		foutLog << "Error: invalid value for MAX_COVERAGE (-maxc)." << endl;
		flag = 0;
	}else if(MAX_COVERAGE<MIN_COVERAGE && MAX_COVERAGE!=-1){
		cerr << "Error: invalid value for MAX_COVERAGE (-maxc)." << endl;
		cerr << "MAX_COVERAGE cannot be smaller than MIN_COVERAGE." << endl;		
		foutLog << "Error: invalid value for MAX_COVERAGE (-maxc)." << endl;
		foutLog << "MAX_COVERAGE cannot be smaller than MIN_COVERAGE." << endl;
		flag = 0;
	}	
	if(ALPHA!=0.2 && ALPHA!=0.15 && ALPHA!=0.1 && ALPHA!=0.05 && ALPHA!=0.01 && ALPHA!=0.005 && ALPHA!=0.001){
		cerr << "Error: invalid value for ALPHA (-a)." << endl;
		cerr << "Current version only allows ALPHA to be 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, or 0.2." << endl;
		foutLog << "Error: invalid value for ALPHA (-a)." << endl;
		foutLog << "Current version only allows ALPHA to be 0.001, 0.005, 0.01, 0.05, 0.1, 0.15, or 0.2." << endl;
		flag = 0;
	}
	if(THRESHOLD<=0){
		cerr << "Error: invalid value for THRESHOLD (-t)." << endl;
		foutLog << "Error: invalid value for THRESHOLD (-t)." << endl;
		flag = 0;
	}
	if(PROCRUSTES_SCALE!=0 && PROCRUSTES_SCALE!=1){
		cerr << "Error: invalid value for PROCRUSTES_SCALE (-rho)." << endl;
		foutLog << "Error: invalid value for PROCRUSTES_SCALE (-rho)." << endl;
		flag = 0;
	}
	if(KNN_ZSCORE < 3){
		cerr << "Error: invalid value for KNN_ZSCORE (-knn)." << endl;
		cerr << "KNN_ZSCORE need to be an integer greater than 2." << endl;
		foutLog << "Error: invalid value for KNN_ZSCORE (-knn)." << endl;
		foutLog << "KNN_ZSCORE need to be an integer greater than 2." << endl;
		flag = 0;	
	}else if(KNN_ZSCORE>REF_SIZE){
		cerr << "Error: invalid value for KNN_ZSCORE (-knn)." << endl;
		cerr << "KNN_ZSCORE cannot be greater than REF_SIZE." << endl;
		foutLog << "Error: invalid value for KNN_ZSCORE (-knn)." << endl;
		foutLog << "KNN_ZSCORE cannot be greater than REF_SIZE." << endl;
		flag = 0;	
	}		
	if(RANDOM_SEED < 0){
		cerr << "Error: invalid value for RANDOM_SEED (-seed)." << endl;
		foutLog << "Error: invalid value for RANDOM_SEED (-seed)." << endl;
		flag = 0;
	}
	if(NUM_THREADS < 1){
		cerr << "Error: invalid value for NUM_THREADS (-nt)." << endl;
		foutLog << "Error: invalid value for NUM_THREADS (-nt)." << endl;
		flag = 0;
	}	
	//============================================================================
	return flag;
}
//################# Function to calculate input table file dimension  ##################
bool get_table_dim(int &nrow, int &ncol, string filename, char separator){
	ifstream fin;
	string str;
	nrow = 0;
	ncol = 0;
	fin.open(filename.c_str());
	if(fin.fail()){
		return false;
	}
	while(!fin.eof()){
		getline(fin, str);
		if(str.length()>0 && str!=" "){
			nrow+=1;
			if(ncol==0){
				bool is_sep=true;    //Previous character is a separator
				for(int i=0; i<str.length(); i++){
					if(str[i]!=separator && i==0){        //Read in the first element
						ncol+=1;
						is_sep=false;
					}else if(str[i]!=separator && i>0 && is_sep){
						ncol+=1;
						is_sep=false;
					}else if(str[i]==separator){
						is_sep=true;
					}	
				}
			}
		}
	}
	return true;
}
//################# Function to perform PCA on genotype matrix (coded as 0, 1, 2, -9)  ##################
int pca_geno(Mat<char> &G, int nPCs, mat &PC, rowvec &PCvar){
	int i=0;
	int j=0;
	int k=0;
	int N = G.n_rows;
	int L = G.n_cols;
	mat M = zeros<mat>(N, N);
	// Normalization and calculate M
	for(j=0; j<L; j++){
		double X1=0;
		double X2=0;
		double NM=0;
		for(i=0; i<N; i++){
			if(G(i,j) != -9){
				X1 += G(i,j);
				X2 += G(i,j)*G(i,j);
				NM++;
			}
		}			
		double colM=X1/NM;
		double colVAR=(X2-NM*pow(colM,2))/(NM-1);
		mat table = zeros<mat>(3, 3);
		if(colVAR != 0){
			for(int g1=0; g1<3; g1++){
				for(int g2=0; g2<3; g2++){
					table(g1,g2) = (g1-colM)*(g2-colM)/colVAR;
				}
			}
			for(i=0; i<N; i++){
				for(k=0; k<=i; k++){
					if(G(i,j)!=-9 && G(k,j)!=-9){
						M(i,k) += table(G(i,j),G(k,j));
					}
				}
			}
		}
	}
	for(i=0; i<N; i++){
		for(k=0; k<i; k++){
			M(k,i) = M(i,k);
		}
	}	
	// Use eigen decomposition to get PCA results	
	vec eigval;
	mat eigvec;
	eig_sym(eigval, eigvec, M, "dc");	// use "divide & conquer" algorithm
	M.clear();
	double eigsum = sum(eigval);
	vec propvar = eigval/eigsum*100;	
	PCvar = zeros<rowvec>(nPCs);
	for(j=0; j<nPCs; j++){
		if(eigval(N-1-j)<0){
			PCvar(j) = 0;
			for(i=0; i<N; i++){
				PC(i,j) = 0;
			}
		}else{
			PCvar(j) = propvar(N-1-j);
			for(i=0; i<N; i++){
				PC(i,j) = eigvec(i, N-1-j)*sqrt(eigval(N-1-j));
			}	
		}
	}
	return 1;
}
