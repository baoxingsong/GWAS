library(compiler)
enableJIT(3)
library(ggplot2)
library("Cairo")

args = commandArgs(trailingOnly=TRUE)

genotype <- read.table(args[1], header=T)

#genotype <- read.table("WC_S.result.assoc.txt", header=T)


colnames(genotype) <- c("chr", "id", "position",  colnames(genotype)[-c(1,2,3)])
genotype$chr <- paste("Chr", genotype$chr, sep="")
genotype <- genotype[(order(genotype$p_wald)),]
genotype$chr = factor(genotype$chr, levels=c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10"))
ymax<-max(-log10(genotype$p_wald))
if (ymax<8) ymax<-8

nMrk <- as.numeric(nrow(genotype))
sign <- 0.05/nMrk

permutation_threshold = 6.4

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
    theme(axis.line = element_blank(), panel.background = element_blank(),panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid") )
}


manhattanPlot <- function(genotype, ymax, sign, addBT=TRUE){
	nMrkT <- as.numeric(nrow(genotype))
    signT <- 0.05/nMrkT
    y_hline = -log10(sign)
    if (ymax < y_hline+0.2){
            ymax = y_hline + 0.2
    }
	if( addBT ){
    	ggplot(data=genotype, aes(x=position, y=-log10(p_wald), group=chr, color=chr))+geom_point(data=genotype, aes(x=position, y=-log10(p_wald), group=chr, color=chr))+
    		guides(colour=FALSE)+facet_grid(~chr, scales="free_x", space="free_x")+scale_x_continuous(labels=changetoM) +
	    	labs(x="", y="-log10(pvalue)", title="") +geom_hline(aes(yintercept=permutation_threshold), color="black", linetype=2)+
	    	theme_bw() +theme_grey(base_size = 36) + theme(axis.line = element_blank(),
	        panel.grid.major = element_blank(),
	        panel.grid.minor = element_blank(),
	        panel.border = element_blank(),
	        panel.background = element_blank(),
	        axis.text.y = element_text( colour = "black"),
	    	axis.text.x = element_text(angle=300, hjust=0, vjust=1, colour = "black"))
	}else{
    	ggplot(data=genotype, aes(x=position, y=-log10(p_wald), group=chr, color=chr))+geom_point(data=genotype, aes(x=position, y=-log10(p_wald), group=chr, color=chr))+
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



prefix = args[1]
plot = manhattanPlot(genotype, ymax, sign)
CairoPNG(paste(prefix, "_mahm.png", sep="") , width=2000, height=450)
print(plot)
dev.off()


qqp <- qqPlot(genotype$p_wald)
png(paste(prefix, "_qq.png", sep=""), width=500,height=460)
print(qqp)
dev.off()


