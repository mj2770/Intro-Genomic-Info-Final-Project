
install.packages("gplots")
installed.packages("heatmap.plus")
installed.packages("RColorBrewer")
installed.packages("data.table")

library("data.table")
library("gplots")
library("heatmap.plus")
library("RColorBrewer")


#---------------------------------------------------------------------
#Read the results
SampleA <- read.csv("D:/E4060/SampleA.csv")
SampleB <- read.csv("D:/E4060/SampleB.csv")
SampleC <- read.csv("D:/E4060/SampleC.csv")
SampleD <- read.csv("D:/E4060/SampleD.csv")
SampleE <- read.csv("D:/E4060/SampleE.csv")

#Remove unclassified taxonomy.
SampleA <- SampleA[!(SampleA$taxon=="unclassified"),]
SampleB <- SampleB[!(SampleB$taxon=="unclassified"),]
SampleC <- SampleC[!(SampleC$taxon=="unclassified"),]
SampleD <- SampleD[!(SampleD$taxon=="unclassified"),]
SampleE <- SampleE[!(SampleE$taxon=="unclassified"),]

#-----------------------------------------------------------------------------------------
#Define a function to explore data in different taxonomy level
#1 is kindom, 2 is phylum, 3 is class, 4 is order, 5 is family and 6 is genus level
taxon <- function(INP1,INP2,INP3,INP4,INP5,taxonlevel){
  S1 <-INP1[(INP1$taxlevel==taxonlevel),]
  S2 <-INP2[(INP2$taxlevel==taxonlevel),]
  S3 <-INP3[(INP3$taxlevel==taxonlevel),]
  S4 <-INP4[(INP4$taxlevel==taxonlevel),]
  S5 <-INP5[(INP5$taxlevel==taxonlevel),]
  S_All <- merge(S1,S2,by="taxon",all=TRUE)
  S_All <- merge(S_All,S3,by="taxon",all=TRUE)
  S_All <- merge(S_All,S4,by="taxon",all=TRUE)
  S_All <- merge(S_All,S5,by="taxon",all=TRUE)
  S_All <- S_All[,c(1,5,9,13,17,21)]
  S_All[is.na(S_All)] <- 0
  names(S_All)[2:6]<-c('SampleA','SampleB','SampleC','SampleD','SampleE')
  return(S_All)
  }

kindom <- taxon(SampleA,SampleB,SampleC,SampleD,SampleE,taxonlevel=1)
phylum <- taxon(SampleA,SampleB,SampleC,SampleD,SampleE,taxonlevel=2)
class <- taxon(SampleA,SampleB,SampleC,SampleD,SampleE,taxonlevel=3)
order <- taxon(SampleA,SampleB,SampleC,SampleD,SampleE,taxonlevel=4)
family <- taxon(SampleA,SampleB,SampleC,SampleD,SampleE,taxonlevel=5)
genus <- taxon(SampleA,SampleB,SampleC,SampleD,SampleE,taxonlevel=6)


#Because there are sequencings not belong to 16s rRNA, we need to recalculate the relative abundance
RelativeAbundance <- function(INPFILE){
  INPFILE$RA1 <- INPFILE$SampleA/sum(INPFILE$SampleA)
  INPFILE$RA2 <- INPFILE$SampleB/sum(INPFILE$SampleB)
  INPFILE$RA3 <- INPFILE$SampleC/sum(INPFILE$SampleC)
  INPFILE$RA4 <- INPFILE$SampleD/sum(INPFILE$SampleD)
  INPFILE$RA5 <- INPFILE$SampleE/sum(INPFILE$SampleE)
  names(INPFILE)[7:11]<-c('SampleA','SampleB','SampleC','SampleD','SampleE')
  rownames(INPFILE) <-INPFILE$taxon
  OUTPUT <- INPFILE[,-c(1:6)]
  OUTPUT <- as.matrix(OUTPUT)
  return(OUTPUT)
}

kindom_ <- RelativeAbundance(kindom)
phylum_ <- RelativeAbundance(phylum)
class_ <- RelativeAbundance(class)
order_ <- RelativeAbundance(order)
family_ <- RelativeAbundance(family)
genus_ <- RelativeAbundance(genus)

#----------------------------------------------------------------------------------------------------------
#Barplot in different taxonomy level

dev.off()
barplot(kindom_, space = NULL, width = 130, xlim = c(-10,1000),
        ylab = "Relative Abundance",col = brewer.pal(n = length(rownames(kindom_)), name = 'RdBu'))
legend ("right",c(rownames(kindom_)),fill = brewer.pal(n = length(rownames(kindom_)), name = 'RdBu'))


dev.off()
barplot(phylum_, space = NULL, width = 120, xlim = c(-10,1000),
        ylab = "Relative Abundance",col = brewer.pal(n = length(rownames(phylum_)), name = 'RdBu'))
legend ("right",c(rownames(phylum_)),fill = brewer.pal(n = length(rownames(phylum_)), name = 'RdBu'))

dev.off()
barplot(class_, space = NULL, width = 120, xlim = c(-10,1000),
        ylab = "Relative Abundance",col = brewer.pal(n = length(rownames(class_)), name = 'RdBu'))
legend ("right",c(rownames(class_)),fill = brewer.pal(n = length(rownames(class_)), name = 'RdBu'))

dev.off()
barplot(order_, space = NULL, width = 125, xlim = c(-10,1000),
        ylab = "Relative Abundance",col = brewer.pal(n = length(rownames(order_)), name = 'RdBu'))
legend ("right",c(rownames(order_)),fill = brewer.pal(n = length(rownames(order_)), name = 'RdBu'), cex = 0.75)

dev.off()
barplot(family_, space = NULL, width = 140, xlim = c(-10,1000),
        ylab = "Relative Abundance",col = brewer.pal(n = length(rownames(family_)), name = 'RdBu'))
legend ("right",c(rownames(family_)),fill = brewer.pal(n = length(rownames(family_)), name = 'RdBu'), cex = 0.4)

#----------------------------------------------------------------------------------------------------------
#Make heatmap in different taxonomy level

dev.off()
colors<-seq(0,1,length.out=21)
heatmap.2(class_, trace="none",density="none", col=heat.colors(20),margins = c(7,10), breaks=colors)

dev.off()
colors<-seq(0,1,length.out=21)
heatmap.2(order_, trace="none",density="none", col=heat.colors(20),margins = c(7,10), breaks=colors)


dev.off()
colors<-seq(0,1,length.out=21)
heatmap.2(family_, trace="none",density="none", col=heat.colors(20),margins = c(7,10), breaks=colors)

