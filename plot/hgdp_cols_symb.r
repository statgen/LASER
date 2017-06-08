########
# This script was written by Chaolong Wang (chaolong@umich.edu) in 2009.
# This script provides a list of colors and symbols for 53 HGDP populations (and 11 HapMap populations).
# The color scheme for the HGDP samples was first used in Rosenberg et al Science (2002) 298: 2381-2385.
# The list of symbols was first used for 29 HGDP populations in Jakobsson et al Nature (2008) 451:998-1003.
# The list of symbols was augmented to include 53 HGDP populations in Wang et al (2010) Stat Appl Genet Mol Biol 9:13.
# The color scheme and symbols have been consistently used in many papers on the HGDP data. Recent examples include
# Wang et al (2012) PLoS Genetics 8:e1002886 
# Wang et al (2014) Nature Genetics 46:409-415
########

myorange <- rgb(1,.4,0)
myblue <- rgb(0,.6,.9)
myyellow <- rgb(220/256,239/256,6/256)
myred <- rgb(250/256,0,38/256)
mypink <- rgb(1,.6,.9)
mygreen <- rgb(.2,.6,.2)
mypurple <- rgb(.5,0,.5)

myoldgold <- rgb(207/255,181/255,59/255);   # For admixtured populations (ASW, MEX)

# population names
popNames=list(
  "430"="BantuSouthAfrica",
  "441"="BantuKenya",
  "464"="Mandenka",
  "465"="Yoruba",
  "488"="BiakaPygmy",
  "489"="MbutiPygmy",
  "494"="San",
  "20"="Orcadian",
  "21"="Adygei",
  "22"="Russian",
  "24"="Basque",
  "25"="French",
  "27"="Italian",
  "28"="Sardinian",
  "29"="Tuscan",
  "34"="Mozabite",
  "36"="Bedouin",
  "37"="Druze",
  "38"="Palestinian",
  "50"="Balochi",
  "51"="Brahui",
  "52"="Burusho",
  "54"="Hazara",
  "56"="Kalash",
  "57"="Makrani",
  "58"="Pathan",
  "59"="Sindhi",
  "601"="Han",
  "602"="HanNChina",
  "606"="Dai",
  "607"="Daur",
  "608"="Hezhen",
  "611"="Lahu",
  "612"="Miao",
  "613"="Oroqen",
  "615"="She",
  "616"="Tujia",
  "617"="Tu",
  "618"="Xibo",
  "619"="Yi",
  "622"="Mongola",
  "625"="Naxi",
  "629"="Uygur",
  "677"="Cambodian",
  "684"="Japanese",
  "699"="Yakut",
  "71"="Melanesian",
  "75"="Papuan",
  "81"="Colombian",
  "82"="Karitiana",
  "83"="Surui",
  "86"="Maya",
  "87"="Pima",
  "900"="ASW",
  "901"="CEU",
  "902"="CHB",
  "903"="CHD",
  "904"="GIH",
  "905"="JPT",
  "906"="LWK",
  "907"="MEX",
  "908"="MKK",
  "909"="TSI",
  "910"="YRI");

# population regions
popRegions = list(
  BantuSouthAfrica="Africa",
  BantuKenya="Africa",
  Mandenka="Africa",
  Yoruba="Africa",
  BiakaPygmy="Africa",
  MbutiPygmy="Africa",
  San="Africa",
  Orcadian="Europe",
  Adygei="Europe",
  Russian="Europe",
  Basque="Europe",
  French="Europe",
  Italian="Europe",
  Sardinian="Europe",
  Tuscan="Europe",
  Mozabite="MiddleEast",
  Bedouin="MiddleEast",
  Druze="MiddleEast",
  Palestinian="MiddleEast",
  Balochi="CentralSouthAsia",
  Brahui="CentralSouthAsia",
  Burusho="CentralSouthAsia",
  Hazara="CentralSouthAsia",
  Kalash="CentralSouthAsia",
  Makrani="CentralSouthAsia",
  Pathan="CentralSouthAsia",
  Sindhi="CentralSouthAsia",
  Han="EastAsia",
  HanNChina="EastAsia",
  Dai="EastAsia",
  Daur="EastAsia",
  Hezhen="EastAsia",
  Lahu="EastAsia",
  Miao="EastAsia",
  Oroqen="EastAsia",
  She="EastAsia",
  Tujia="EastAsia",
  Tu="EastAsia",
  Xibo="EastAsia",
  Yi="EastAsia",
  Mongola="EastAsia",
  Naxi="EastAsia",
  Uygur="CentralSouthAsia",
  Cambodian="EastAsia",
  Japanese="EastAsia",
  Yakut="EastAsia",
  Melanesian="Oceania",
  Papuan="Oceania",
  Colombian="America",
  Karitiana="America",
  Surui="America",
  Maya="America",
  Pima="America",
  ASW="admix1",  # Admixed population
  CEU="Europe",
  CHB="EastAsia",
  CHD="EastAsia",
  GIH="CentralSouthAsia",
  JPT="EastAsia",
  LWK="Africa",
  MEX="admix2",    # Admixed population
  MKK="Africa",
  TSI="Europe",
  YRI="Africa");

# color references

  Africa <- myorange;
  Europe <- myblue;
  MiddleEast <- myyellow;
  CentralSouthAsia <- myred;
  EastAsia <- mypink;
  Oceania <- mygreen;
  America <- mypurple;
  admix <- myoldgold;

# plot symbols

# Africa
  BantuSouthAfrica <- 2;
  BantuKenya <- 17;
  Mandenka <- 5;
  Yoruba <- 1;
  BiakaPygmy <- 18;
  MbutiPygmy <- 6;
  San <- 0;
  

# Europe  
  Orcadian <- 1;
  Adygei <- 17;
  Russian <- 2;
  Basque <- 23;
  French <- 3;
  Italian <- 4;
  Sardinian <- 18;
  Tuscan <- 6;
  
  
# MiddleEast  
  Mozabite <- 1;
  Bedouin <- 15;
  Druze <- 19;
  Palestinian <- 0;
  

# CSAsia  
  Balochi <- 4;
  Brahui <- 1;
  Burusho <- 15;
  Hazara <- 2;
  Kalash <- 18;
  Makrani <- 5;
  Pathan <- 6;
  Sindhi <- 16;
  Uygur <- 3;
  

# EAsia  
  Han <- 16;
  HanNChina <- 4;
  Dai <- 3;
  Daur <- 15;
  Hezhen <- 17;
  Lahu <- 5;
  Miao <- 18;
  Oroqen <- 13;
  She <- 7;
  Tujia <- 8;
  Tu <- 12;
  Xibo <- 14;
  Yi <- 0;
  Mongola <- 2;
  Naxi <- 10;
  Cambodian <- 1;
  Japanese <- 9;
  Yakut <- 6;  

# Oceania  
  Melanesian <- 1;
  Papuan <- 0;

# America  
  Colombian <- 5;
  Karitiana <- 1;
  Surui <- 0;
  Maya <- 2;
  Pima <- 6;
    
# HapMap
  ASW <- 65; #"A";  
  CEU <- 67; #"C"; 
  CHB <- 66; #"B"; 
  CHD <- 68; #"D"; 
  GIH <- 71; #"G"; 
  JPT <- 74; #"J"; 
  LWK <- 76; #"L"; 
  MEX <- 77; #"M"; 
  MKK <- 75; #"K"; 
  TSI <- 84; #"T"; 
  YRI <- 89; #"Y"; 