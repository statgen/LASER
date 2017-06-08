args <- commandArgs(trailingOnly = TRUE)

load_prefix <- args[1];
geno_prefix <- args[2];
out_prefix <- args[3];

geno <- read.table(paste(geno_prefix, ".geno", sep=""), header=F, sep="\t", na.strings="-9");
proj <- read.table(paste(load_prefix, ".RefPC.load", sep=""), header=T, sep="\t");

info <- geno[,1:2];
colnames(info) <- c("popID", "indivID");
geno <- geno[,-c(1:2)];

rmlist <- which(proj$Sd==0);
if(length(rmlist)>0){
	geno <- geno[,-rmlist];
	proj <- proj[-rmlist];
}
X <- apply(geno, 1, function(x) (x-proj$Mean)/proj$Sd)
X[which(is.na(X))] <- 0;

W <- as.matrix(proj[,-c(1:3)]);
Y <- t(X)%*%W;

outtable <- cbind(info, Y);
outfile <- paste(out_prefix, ".proj.coord", sep="")
write.table(outtable, file=outfile, sep="\t", row.names=F, col.names=T, quote=F)

