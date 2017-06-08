########
# This script was written by Chaolong Wang (chaolong@umich.edu) in 2009
# This script include a function to assign colors and symbols for plotting spatial maps of 53 HGDP populations.
# The color scheme for the HGDP samples was first used in Rosenberg et al Science (2002) 298: 2381-2385.
# The list of symbols was first used for 29 HGDP populations in Jakobsson et al Nature (2008) 451:998-1003.
# The list of symbols was augmented to include 53 HGDP populations in Wang et al (2010) Stat Appl Genet Mol Biol 9:13.
# The color scheme and symbols have been consistently used in many papers on the HGDP data. Recent examples include
# Wang et al (2012) PLoS Genetics 8:e1002886 
# Wang et al (2014) Nature Genetics 46:409-415
########

source("hgdp_cols_symb.r");

assign_cols_symb <- function(popID){
	N_ind <- length(popID);
	cols <-	rep("darkgrey", N_ind);
	symb <- rep(1, N_ind);
	for (i in 1:N_ind){
		if(is.na(popID[i])==FALSE){
			if(popID[i]==430 || as.character(popID[i])=="BantuSouthAfrica"){  	#########  HGDP #########
				symb[i] <- BantuSouthAfrica;
				cols[i] <- Africa;
			}else if(popID[i]==441 || as.character(popID[i])=="BantuKenya"){
				symb[i] <- BantuKenya;
				cols[i] <- Africa;
			}else if(popID[i]==464 || as.character(popID[i])=="Mandenka"){
				symb[i] <- Mandenka;
				cols[i] <- Africa;
			}else if(popID[i]==465 || as.character(popID[i])=="Yoruba"){
				symb[i] <- Yoruba;
				cols[i] <- Africa;
			}else if(popID[i]==488 || as.character(popID[i])=="BiakaPygmy"){
				symb[i] <- BiakaPygmy;
				cols[i] <- Africa;
			}else if(popID[i]==489 || as.character(popID[i])=="MbutiPygmy"){
				symb[i] <- MbutiPygmy;
				cols[i] <- Africa;
			}else if(popID[i]==494 || as.character(popID[i])=="San"){
				symb[i] <- San;
				cols[i] <- Africa;
			}else if(popID[i]==20 || as.character(popID[i])=="Orcadian"){
				symb[i] <- Orcadian;
				cols[i] <- Europe;
			}else if(popID[i]==21 || as.character(popID[i])=="Adygei"){
				symb[i] <- Adygei;
				cols[i] <- Europe;
			}else if(popID[i]==22 || as.character(popID[i])=="Russian"){
				symb[i] <- Russian;
				cols[i] <- Europe;
			}else if(popID[i]==24 || as.character(popID[i])=="Basque"){
				symb[i] <- Basque;
				cols[i] <- Europe;
			}else if(popID[i]==25 || as.character(popID[i])=="French"){
				symb[i] <- French;
				cols[i] <- Europe;
			}else if(popID[i]==27 || as.character(popID[i])=="Italian"){
				symb[i] <- Italian;
				cols[i] <- Europe;
			}else if(popID[i]==28 || as.character(popID[i])=="Sardinian"){
				symb[i] <- Sardinian;
				cols[i] <- Europe;
			}else if(popID[i]==29 || as.character(popID[i])=="Tuscan"){
				symb[i] <- Tuscan;
				cols[i] <- Europe;
			}else if(popID[i]==34 || as.character(popID[i])=="Mozabite"){
				symb[i] <- Mozabite;
				cols[i] <- MiddleEast;
			}else if(popID[i]==36 || as.character(popID[i])=="Bedouin"){
				symb[i] <- Bedouin;
				cols[i] <- MiddleEast;
			}else if(popID[i]==37 || as.character(popID[i])=="Druze"){
				symb[i] <- Druze;
				cols[i] <- MiddleEast;
			}else if(popID[i]==38 || as.character(popID[i])=="Palestinian"){
				symb[i] <- Palestinian;
				cols[i] <- MiddleEast;
			}else if(popID[i]==50 || as.character(popID[i])=="Balochi"){
				symb[i] <- Balochi;
				cols[i] <- CentralSouthAsia;
			}else if(popID[i]==51 || as.character(popID[i])=="Brahui"){
				symb[i] <- Brahui;
				cols[i] <- CentralSouthAsia;
			}else if(popID[i]==52 || as.character(popID[i])=="Burusho"){
				symb[i] <- Burusho;
				cols[i] <- CentralSouthAsia;
			}else if(popID[i]==54 || as.character(popID[i])=="Hazara"){
				symb[i] <- Hazara;
				cols[i] <- CentralSouthAsia;
			}else if(popID[i]==56 || as.character(popID[i])=="Kalash"){
				symb[i] <- Kalash;
				cols[i] <- CentralSouthAsia;
			}else if(popID[i]==57 || as.character(popID[i])=="Makrani"){
				symb[i] <- Makrani;
				cols[i] <- CentralSouthAsia;
			}else if(popID[i]==58 || as.character(popID[i])=="Pathan"){
				symb[i] <- Pathan;
				cols[i] <- CentralSouthAsia;
			}else if(popID[i]==59 || as.character(popID[i])=="Sindhi"){
				symb[i] <- Sindhi;
				cols[i] <- CentralSouthAsia;
			}else if(popID[i]==601 || as.character(popID[i])=="Han"){
				symb[i] <- Han;
				cols[i] <- EastAsia;
			}else if(popID[i]==602 || as.character(popID[i])=="Han-NChina"){
				symb[i] <- HanNChina;
				cols[i] <- EastAsia;
			}else if(popID[i]==606 || as.character(popID[i])=="Dai"){
				symb[i] <- Dai;
				cols[i] <- EastAsia;
			}else if(popID[i]==607 || as.character(popID[i])=="Daur"){
				symb[i] <- Daur;
				cols[i] <- EastAsia;
			}else if(popID[i]==608 || as.character(popID[i])=="Hezhen"){
				symb[i] <- Hezhen;
				cols[i] <- EastAsia;
			}else if(popID[i]==611 || as.character(popID[i])=="Lahu"){
				symb[i] <- Lahu;
				cols[i] <- EastAsia;
			}else if(popID[i]==612 || as.character(popID[i])=="Miao"){
				symb[i] <- Miao;
				cols[i] <- EastAsia;
			}else if(popID[i]==613 || as.character(popID[i])=="Oroqen"){
				symb[i] <- Oroqen;
				cols[i] <- EastAsia;
			}else if(popID[i]==615 || as.character(popID[i])=="She"){
				symb[i] <- She;
				cols[i] <- EastAsia;
			}else if(popID[i]==616 || as.character(popID[i])=="Tujia"){
				symb[i] <- Tujia;
				cols[i] <- EastAsia;
			}else if(popID[i]==617 || as.character(popID[i])=="Tu"){
				symb[i] <- Tu;
				cols[i] <- EastAsia;
			}else if(popID[i]==618 || as.character(popID[i])=="Xibo"){
				symb[i] <- Xibo;
				cols[i] <- EastAsia;
			}else if(popID[i]==619 || as.character(popID[i])=="Yi"){
				symb[i] <- Yi;
				cols[i] <- EastAsia;
			}else if(popID[i]==622 || as.character(popID[i])=="Mongola"){
				symb[i] <- Mongola;
				cols[i] <- EastAsia;
			}else if(popID[i]==625 || as.character(popID[i])=="Naxi"){
				symb[i] <- Naxi;
				cols[i] <- EastAsia;
			}else if(popID[i]==629 || as.character(popID[i])=="Uygur"){
				symb[i] <- Uygur;
				cols[i] <- CentralSouthAsia;
			}else if(popID[i]==677 || as.character(popID[i])=="Cambodian"){
				symb[i] <- Cambodian;
				cols[i] <- EastAsia;
			}else if(popID[i]==684 || as.character(popID[i])=="Japanese"){
				symb[i] <- Japanese;
				cols[i] <- EastAsia;
			}else if(popID[i]==699 || as.character(popID[i])=="Yakut"){
				symb[i] <- Yakut;
				cols[i] <- EastAsia;
			}else if(popID[i]==71 || as.character(popID[i])=="Melanesian"){
				symb[i] <- Melanesian;
				cols[i] <- Oceania;
			}else if(popID[i]==75 || as.character(popID[i])=="Papuan"){
				symb[i] <- Papuan;
				cols[i] <- Oceania;
			}else if(popID[i]==81 || as.character(popID[i])=="Colombian"){
				symb[i] <- Colombian;
				cols[i] <- America;
			}else if(popID[i]==82 || as.character(popID[i])=="Karitiana"){
				symb[i] <- Karitiana;
				cols[i] <- America;
			}else if(popID[i]==83 || as.character(popID[i])=="Surui"){
				symb[i] <- Surui;
				cols[i] <- America;
			}else if(popID[i]==86 || as.character(popID[i])=="Maya"){
				symb[i] <- Maya;
				cols[i] <- America;
			}else if(popID[i]==87 || as.character(popID[i])=="Pima"){
				symb[i] <- Pima;
				cols[i] <- America;
			}else if(popID[i]==900 || as.character(popID[i])=="ASW"){   #########  HapMap #########
				symb[i] <- ASW;
				cols[i] <- Africa;
			}else if(popID[i]==901 || as.character(popID[i])=="CEU"){
				symb[i] <- CEU;
				cols[i] <- Europe;
			}else if(popID[i]==902 || as.character(popID[i])=="CHB"){
				symb[i] <- CHB;
				cols[i] <- EastAsia;
			}else if(popID[i]==903 || as.character(popID[i])=="CHD"){
				symb[i] <- CHD;
				cols[i] <- EastAsia;
			}else if(popID[i]==904 || as.character(popID[i])=="GIH"){
				symb[i] <- GIH;
				cols[i] <- CentralSouthAsia;
			}else if(popID[i]==905 || as.character(popID[i])=="JPT"){
				symb[i] <- JPT;
				cols[i] <- EastAsia;
			}else if(popID[i]==906 || as.character(popID[i])=="LWK"){
				symb[i] <- LWK;
				cols[i] <- Africa;
			}else if(popID[i]==907 || as.character(popID[i])=="MXL" || as.character(popID[i])=="MEX"){
				symb[i] <- MEX;
				cols[i] <- admix;
			}else if(popID[i]==908 || as.character(popID[i])=="MKK"){
				symb[i] <- MKK;
				cols[i] <- Africa;
			}else if(popID[i]==909 || as.character(popID[i])=="TSI"){
				symb[i] <- TSI;
				cols[i] <- Europe;
			}else if(popID[i]==910 || as.character(popID[i])=="YRI"){
				symb[i] <- YRI;
				cols[i] <- Africa;
			}else if(as.character(popID[i])=="AL"){       #########  POPRES #########
				cols[i] <- AL;
			}else if(as.character(popID[i])=="AT"){
				cols[i] <- AT;
			}else if(as.character(popID[i])=="BA"){
				cols[i] <- BA;
			}else if(as.character(popID[i])=="BE"){
				cols[i] <- BE;
			}else if(as.character(popID[i])=="BG"){
				cols[i] <- BG;
			}else if(as.character(popID[i])=="CH"){
				cols[i] <- CH;
			}else if(as.character(popID[i])=="CH-F"){
				cols[i] <- CH;
			}else if(as.character(popID[i])=="CH-G"){
				cols[i] <- CH;
			}else if(as.character(popID[i])=="CH-I"){
				cols[i] <- CH;
			}else if(as.character(popID[i])=="CY"){
				cols[i] <- CY;
			}else if(as.character(popID[i])=="CZ"){
				cols[i] <- CZ;
			}else if(as.character(popID[i])=="DE"){
				cols[i] <- DE;
			}else if(as.character(popID[i])=="DK"){
				cols[i] <- DK;
			}else if(as.character(popID[i])=="ES"){
				cols[i] <- ES;
			}else if(as.character(popID[i])=="FI"){
				cols[i] <- FI;
			}else if(as.character(popID[i])=="FR"){
				cols[i] <- FR;
			}else if(as.character(popID[i])=="GB"){
				cols[i] <- GB;
			}else if(as.character(popID[i])=="GR"){
				cols[i] <- GR;
			}else if(as.character(popID[i])=="HR"){
				cols[i] <- HR;
			}else if(as.character(popID[i])=="HU"){
				cols[i] <- HU;
			}else if(as.character(popID[i])=="IE"){
				cols[i] <- IE;
			}else if(as.character(popID[i])=="IT"){
				cols[i] <- IT;
			}else if(as.character(popID[i])=="KS"){
				cols[i] <- KS;
			}else if(as.character(popID[i])=="LV"){
				cols[i] <- LV;
			}else if(as.character(popID[i])=="MK"){
				cols[i] <- MK;
			}else if(as.character(popID[i])=="NL"){
				cols[i] <- NL;
			}else if(as.character(popID[i])=="NO"){
				cols[i] <- NO;
			}else if(as.character(popID[i])=="PL"){
				cols[i] <- PL;
			}else if(as.character(popID[i])=="PT"){
				cols[i] <- PT;
			}else if(as.character(popID[i])=="RO"){
				cols[i] <- RO;
			}else if(as.character(popID[i])=="RS"){
				cols[i] <- RS;
			}else if(as.character(popID[i])=="RU"){
				cols[i] <- RU;
			}else if(as.character(popID[i])=="Sct"){
				cols[i] <- Sct;
			}else if(as.character(popID[i])=="SE"){
				cols[i] <- SE;
			}else if(as.character(popID[i])=="SI"){
				cols[i] <- SI;
			}else if(as.character(popID[i])=="SK"){
				cols[i] <- SK;
			}else if(as.character(popID[i])=="TR"){
				cols[i] <- TR;
			}else if(as.character(popID[i])=="UA"){
				cols[i] <- UA;
			}else if(as.character(popID[i])=="YG"){
				cols[i] <- YG;
			}else{
				print(paste("Warning: ", popID[i], " is not in the population list:", i));
			}
		}else{
			print(paste("Error: NA found in the popID:", i));
		}
	}
	return(list(cols=cols, symb=symb));
}
