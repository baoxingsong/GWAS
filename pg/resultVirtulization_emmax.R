library(compiler)
enableJIT(3)
library(ggplot2)
library("Cairo")

snp <- read.table("snp_emmax_abn.ps")
colnames(snp) <- c("id", "beta", "pvalue")
snpmarkers <- read.table("snp.map")
colnames(snpmarkers) <- c("chr", "id", "zero", "position")
snp<- merge(x=snp, y=snpmarkers, by="id")
snp <- snp[(order(snp$pvalue)),]
snp$chr <- paste("Chr", snp$chr, sep="")

indel <- read.table("indel3_song_emmax_abn.ps")
indel = indel[,-4]
colnames(indel) <- c("id", "beta", "pvalue")
indelmarkers <- read.table("../../indel/indel.map")
colnames(indelmarkers) <- c("chr", "id", "zero", "position")
indel <- merge(x=indel, y=indelmarkers, by="id")
indel <- indel[(order(indel$pvalue)),]
indel$chr <- paste("Chr", indel$chr, sep="")

orf <- read.table("orfn_emmax_abn.ps")
colnames(orf) <- c("id", "beta", "pvalue")
orfmarkers <- read.table("orfn.map")
colnames(orfmarkers) <- c("chr", "id", "zero", "position")
orf<- merge(x=orf, y=orfmarkers, by="id")
orf <- orf[(order(orf$pvalue)),]
orf$chr <- paste("Chr", orf$chr, sep="")


genotype = rbind(snp, indel)
genotype = rbind(genotype, orf)
genotype <- genotype[(order(genotype$pvalue)),]


ymax<-max(-log10(genotype$pvalue))
if (ymax<8) ymax<-8

nMrk <- as.numeric(nrow(genotype))
sign <- 0.05/nMrk

permutation_threshold = 6.85974791698633

changetoM <- function ( position ){
	position=position/1000000;
	paste(position, "M", sep="")
}

qqPlot <- function (pval, truncate = FALSE, ylim = NULL,
    ci = TRUE, ...) {
	chisqs = qchisq(1-pval, 1)
    lambda <- median(chisqs)/qchisq(0.5, 1)
    pval <- -log10(sort(pval))
    n <- length(pval)
    x <- ppoints(n)
    a <- 1:n
    ind <- 1:n
    upper <- qbeta(0.025, a[ind], rev(a)[ind])
    lower <- qbeta(0.975, a[ind], rev(a)[ind])
    polyArea <- data.frame(x=-log10(c(x[ind], rev(x[ind]))), y=-log10(c(upper, rev(lower))))

    points <- data.frame(observed=pval, expected=-log10(x))
    ymax = max(pval)

    ggplot()+geom_polygon(data=polyArea, aes(x=x, y=y), fill="gray")+geom_point(data=points, aes(x=expected, y=observed))+geom_abline(intercept = 0, slope=1, colour="red")+
    annotate("text", x = 1, y = ymax-0.5, label = paste("λ=", round(lambda, 3), sep=""), size=10)+labs(x="-log10(Expected P)", y="-log10(Observed P)", fill="", title="") +
    theme_grey(base_size = 36) + 
    theme(axis.line = element_line(colour = "black"), panel.background = element_blank(),panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid") )
}


manhattanPlot <- function(genotype, ymax, sign, addBT=TRUE){
	nMrkT <- as.numeric(nrow(genotype))
    signT <- 0.05/nMrkT
    y_hline = -log10(sign)
    if (ymax < y_hline+0.2){
            ymax = y_hline + 0.2
    }
	if( addBT ){
    	ggplot(data=genotype, aes(x=position, y=-log10(pvalue), group=chr, color=chr))+geom_point(data=genotype, aes(x=position, y=-log10(pvalue), group=chr, color=chr))+
    		guides(colour=FALSE)+facet_grid(~chr, scales="free_x", space="free_x")+scale_x_continuous(labels=changetoM) +
	    	labs(x="", y="-log10(pvalue)", title="") + ylim(0, ymax)+ geom_hline(aes(yintercept=-log10(sign)), color="red")+geom_hline(aes(yintercept=permutation_threshold), color="black", linetype=2)+
	    	geom_hline(aes(yintercept=-log10(signT)), color="green", linetype=2)+
	    	theme_bw() +theme_grey(base_size = 36) + theme(axis.line = element_blank(),
	        panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank(),
	        panel.border = element_blank(),
	        panel.background = element_blank(),
	        axis.text.y = element_text( colour = "black"),
	    	axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
	}else{
    	ggplot(data=genotype, aes(x=position, y=-log10(pvalue), group=chr, color=chr))+geom_point(data=genotype, aes(x=position, y=-log10(pvalue), group=chr, color=chr))+
    		guides(colour=FALSE)+facet_grid(~chr, scales="free_x", space="free_x")+scale_x_continuous(labels=changetoM) +
	    	labs(x="", y="-log10(pvalue)", title="") + ylim(0, ymax)+ geom_hline(aes(yintercept=-log10(sign)), color="red")+geom_hline(aes(yintercept=permutation_threshold), color="black", linetype=2)+
	    	theme_bw() +theme_grey(base_size = 36) + theme(axis.line = element_blank(),
	        panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank(),
	        panel.border = element_blank(),
	        panel.background = element_blank(),
	        axis.text.y = element_text( colour = "black"),
	    	axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
	}
}


genotypeT = snp
prefix = "snp"
plot = manhattanPlot(genotypeT, ymax, sign)
CairoPNG(paste(prefix, "_mahm.png", sep="") , width=2000, height=450)
plot
dev.off()

CairoPDF(paste(prefix, "_mahm.pdf", sep=""), width=30, height=6.9)
plot
dev.off()

qqp <- qqPlot(genotypeT$pvalue)
png(paste(prefix, "_qq.png", sep=""), width=500,height=460)
qqp
dev.off()

pdf(paste(prefix, "_qq.pdf", sep=""), width=7.5,height=6.9)
qqp
dev.off()

load("../../../tair10")
top2 <- snp[which( (-log10(snp$pvalue))>permutation_threshold ),]

sink(paste(prefix, "_summary.txt",sep=""))
for (i in  1:nrow(top2) ) {
	print (top2[i,])
	print (tair10[which(tair10$V1 == top2[i,]$chr  & tair10$V4 <= top2[i,]$position & tair10$V5 >= top2[i,]$position), ]$V9, max.levels=0)
	print ("")
	print ("")
}
sink()

genotypeT = indel
prefix = "indel"
plot = manhattanPlot(genotypeT, ymax, sign)
CairoPNG(paste(prefix, "_mahm.png", sep="") , width=2000, height=450)
plot
dev.off()

CairoPDF(paste(prefix, "_mahm.pdf", sep=""), width=30, height=6.9)
plot
dev.off()

qqp <- qqPlot(genotypeT$pvalue)
png(paste(prefix, "_qq.png", sep=""), width=500,height=460)
qqp
dev.off()

pdf(paste(prefix, "_qq.pdf", sep=""), width=7.5,height=6.9)
qqp
dev.off()

load("../../../tair10")
top2 <- indel[which( (-log10(indel$pvalue))>permutation_threshold ),]

sink(paste(prefix, "_summary.txt",sep=""))
for (i in  1:nrow(top2) ) {
    print (top2[i,])
    print (tair10[which(tair10$V1 == top2[i,]$chr  & tair10$V4 <= top2[i,]$position & tair10$V5 >= top2[i,]$position), ]$V9, max.levels=0)
    print ("")
    print ("")
}
sink()


genotypeT = orf
prefix = "orf"
plot = manhattanPlot(genotypeT, ymax, sign)
CairoPNG(paste(prefix, "_mahm.png", sep="") , width=2000, height=450)
plot
dev.off()

CairoPDF(paste(prefix, "_mahm.pdf", sep=""), width=30, height=6.9)
plot
dev.off()

qqp <- qqPlot(genotypeT$pvalue)
png(paste(prefix, "_qq.png", sep=""), width=500,height=460)
qqp
dev.off()

pdf(paste(prefix, "_qq.pdf", sep=""), width=7.5,height=6.9)
qqp
dev.off()

load("../../../tair10")
top2 <- orf[which( (-log10(orf$pvalue))>permutation_threshold ),]

sink(paste(prefix, "_summary.txt",sep=""))
for (i in  1:nrow(top2) ) {
    print (top2[i,])
    print (tair10[which(tair10$V1 == top2[i,]$chr  & tair10$V4 <= top2[i,]$position & tair10$V5 >= top2[i,]$position), ]$V9, max.levels=0)
    print ("")
    print ("")
}
sink()
all = genotype
genotypeT = genotype
prefix = "all"
plot = manhattanPlot(genotypeT, ymax, sign, addBT=FALSE)
CairoPNG(paste(prefix, "_mahm.png", sep="") , width=2000, height=450)
plot
dev.off()

CairoPDF(paste(prefix, "_mahm.pdf", sep=""), width=30, height=6.9)
plot
dev.off()

qqp <- qqPlot(genotypeT$pvalue)
png(paste(prefix, "_qq.png", sep=""), width=500,height=460)
qqp
dev.off()

pdf(paste(prefix, "_qq.pdf", sep=""), width=7.5,height=6.9)
qqp
dev.off()

load("../../../tair10")
top2 <- all[which( (-log10(all$pvalue))>permutation_threshold ),]

sink(paste(prefix, "_summary.txt",sep=""))
for (i in  1:nrow(top2) ) {
    print (top2[i,])
    print (tair10[which(tair10$V1 == top2[i,]$chr  & tair10$V4 <= top2[i,]$position & tair10$V5 >= top2[i,]$position), ]$V9, max.levels=0)
    print ("")
    print ("")
}
sink()
