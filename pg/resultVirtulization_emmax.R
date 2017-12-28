library(compiler)
enableJIT(3)
library(ggplot2)
library("Cairo")
args <- commandArgs(trailingOnly = TRUE)

print ("usage: Rscript resultVirtulization_emmax.R emmax_result map_file output_prefix")

print(args)

snp <- read.table(args[1]) # the output of emmax, end with .ps
colnames(snp) <- c("id", "beta", "pvalue")
markers <- read.table(args[2]) # .map file
colnames(markers) <- c("chr", "id", "zero", "position")

snp<- merge(x=snp, y=markers, by="id")

snp <- snp[(order(snp$pvalue)),]
nMrk <- nrow(snp)
nMrk <- as.numeric(nMrk)
sign <- 0.05/nMrk
#print(paste("sign", sign))
changetoM <- function ( position ){
	position=position/1000000;
	paste(position, "M", sep="")
}


qqPlot <- function (pval, truncate = FALSE, ylim = NULL,
    ci = TRUE, ...) 
{

    lambda <- median(-2*log(pval)) / 1.39
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
#    ggplot()+geom_polygon(data=polyArea, aes(x=x, y=y), fill="gray")+geom_point(data=points, aes(x=expected, y=observed))+geom_abline(intercept = 0, slope=1, colour="red")+annotate("text", x = 1, y = ymax-0.5, label = paste(expression(lambda), sep=""), size=10)+labs(x="-log10(Expected P)", y="-log10(Observed P)", fill="", title="") +theme_grey(base_size = 36) + theme(axis.line = element_line(colour = "black"), panel.background = element_blank(),panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid") )

}

snp$chr <- paste("Chr", snp$chr, sep="")
ymax<-max(-log10(snp$pvalue))
if (ymax<8) ymax<-8
plot <-ggplot(data=snp, aes(x=position, y=-log10(pvalue), group=chr, color=chr))+geom_point(data=snp, aes(x=position, y=-log10(pvalue), group=chr, color=chr))+
	guides(colour=FALSE)+facet_grid(~chr, scales="free_x", space="free_x")+scale_x_continuous(labels=changetoM) +
	labs(x="", y="-log10(pvalue)", title="") + ylim(0, ymax)+ geom_hline(aes(yintercept=-log10(sign)), color="red")+
	theme_bw() +theme_grey(base_size = 36) + theme(axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text( colour = "black"),
	axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))


CairoPNG(paste(args[3], "_mahm.png", sep=""), width=2000, height=450)
plot
dev.off()

CairoPDF(paste(args[3], "_mahm.pdf", sep=""), width=30, height=6.9)
plot
dev.off()

qqp <- qqPlot(snp$pvalue)
png(paste(args[3], "_qq.png", sep=""), width=500,height=460)
qqp
dev.off()

pdf(paste(args[3], "_qq.pdf", sep=""), width=7.5,height=6.9)
qqp
dev.off()


#load("../../../tair10")
tair10 =  read.table("TAIR10_GFF3_genes_no_UM.gff") #  annotation in gff/gtf format
top2 <- snp[which( (-log10(snp$pvalue))>-log10(sign) ),]

sink(paste(args[3], "_summary.txt",sep=""))
for (i in  1:nrow(top2) ) {
	print (top2[i,])
	print (tair10[which(tair10$V1 == top2[i,]$chr  & tair10$V4 <= top2[i,]$position & tair10$V5 >= top2[i,]$position), ]$V9, max.levels=0)
	print ("")
	print ("")
}
sink()

