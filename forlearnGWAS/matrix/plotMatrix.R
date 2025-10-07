system("perl individualList.pl ./snp_matrix.fam > fullIndividualList");

require(reshape2)
library(ggplot2)
library(Cairo)

data <- read.table("plink.mibs")
list <- read.table("fullIndividualList")
data <- as.matrix(data)
names <- list$V3

#names[which((c(1:length(names)) %%2) ==0)] = paste(names[which((c(1:length(names)) %%2) ==0)], "      ", sep="")

colnames(data) <- names
rownames(data) <- names

list$V3 <- sort(list$V3)
names <- list$V3
oddNames <- names[which((c(1:length(names)) %%2) ==0)]

dat <- melt(data)
dat$Var1 <- as.character(dat$Var1)
dat$Var1[which(dat$Var1 %in% oddNames)] = paste(dat$Var1[which(dat$Var1 %in% oddNames)],"         ", sep="")
dat$Var2 <- as.character(dat$Var2)
dat$Var2[which(dat$Var2 %in% oddNames)] = paste(dat$Var2[which(dat$Var2 %in% oddNames)],"         ", sep="")


p <- ggplot(data =  dat, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value), colour = "white") + 
	scale_fill_gradient(low = "green", high = "red")+theme_bw()+labs(x="", y="") +theme_grey(base_size = 36)+ theme(axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
	axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))

CairoPNG(file="matrix.png", width=3200, height=3200)
p
dev.off()

CairoPDF(file="matrix.pdf", width=96, height=96)
p
dev.off()




