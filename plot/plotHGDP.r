mydata <- read.table("../example/test.SeqPC.coord", header=TRUE, sep="\t");         ## Output coordinate file of LASER or TRACE for the study samples
Ref <- read.table("../example/HGDP_238_chr22.RefPC.coord", header=TRUE, sep="\t");   ## Reference coordinate file used in LASER or TRACE

source("assign_cols_symb.r")
ref <- assign_cols_symb(Ref$popID);

#############################################################################################
leg.region.txt <- c("Africa", "Europe", "Middle East", "C/S Asia", "East Asia", "Oceania", "America");
leg.region.col <- c(Africa, Europe, MiddleEast, CentralSouthAsia, EastAsia, Oceania, America);
#############################################################################################

pdf("Results_on_HGDP.pdf", width=7, height=7.9)
par(mfrow=c(2,2));

rangePC1 <- c(min(c(Ref$PC1,mydata$PC1)), max(c(Ref$PC1,mydata$PC1)));
rangePC2 <- c(min(c(Ref$PC2,mydata$PC2)), max(c(Ref$PC2,mydata$PC2)));
rangePC3 <- c(min(c(Ref$PC3,mydata$PC3)), max(c(Ref$PC3,mydata$PC3)));
rangePC4 <- c(min(c(Ref$PC4,mydata$PC4)), max(c(Ref$PC4,mydata$PC4)));

plot(Ref$PC1,  Ref$PC2, xlab="PC1", ylab="PC2", xlim=rangePC1, ylim=rangePC2, main="Reference map", type="n", asp=1)
points(Ref$PC1, Ref$PC2, pch=ref$symb, col=ref$cols, cex=0.7);

plot(Ref$PC3,  Ref$PC4, xlab="PC3", ylab="PC4", xlim=rangePC3, ylim=rangePC4, main="Reference map", type="n", asp=1)
points(Ref$PC3, Ref$PC4, pch=ref$symb, col=ref$cols, cex=0.7);
legend("topleft", leg.region.txt, fill=leg.region.col, ncol=1, cex=0.8);

plot(Ref$PC1,  Ref$PC2, xlab="PC1", ylab="PC2", xlim=rangePC1, ylim=rangePC2, main="Study samples", type="n", asp=1)
points(Ref$PC1, Ref$PC2, pch=ref$symb, col="grey", cex=0.7);
points(mydata$PC1, mydata$PC2, pch=4, col="black", cex=0.7);

plot(Ref$PC3,  Ref$PC4, xlab="PC3", ylab="PC4", xlim=rangePC3, ylim=rangePC4, main="Study samples", type="n", asp=1)
points(Ref$PC3, Ref$PC4, pch=ref$symb, col="grey", cex=0.7);
points(mydata$PC3, mydata$PC4, pch=4, col="black", cex=0.7);

dev.off();
