#

#figure5

#

# (c) Martha Anita Demba (MRC at LSHTM)

# Date: 2022-07-12

#

# Purpose:

# Calculate LD

# Calculate and generate LD heatmap plots using the ldHeatmap Package in R.

# ==========================================================

ldheatmapSNP = function(dir, file){

  library(LDheatmap)

  require(snpStats)

  library(data.table)

  setwd(dir)

  dat = as.data.frame(fread(file, header = T))

  row.names(dat) =dat$V1

  dat1=dat[,-1]

  

  dat1[dat1=="./."] = NA

  dat2=as.data.frame(t(dat1))

  dat3 <- Filter(function(x)!all(is.na(x)), dat2)

  for (i in 1:ncol(dat3)) {

    dat3[,i]=as.genotype(dat3[,i])

  }

  

  snpName = gsub('(.*)_(.*)_\\w+', '\\1', colnames(dat3))

  snpName = gsub(',(.*)\\w+', '', snpName)

  colnames(dat3)=snpName

  rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")

  pdf("LdheatmapSnp.pdf", width = 10)

  MyHeatmap <- LDheatmap(dat3, LDmeasure="r",

                         title="Pairwise LD in r^2", add.map=F, name="myLDgrob",

                         add.key=TRUE, color=rgb.palette(18), SNP.name=colnames(dat3))

  dev.off()

  

}

ldheatmapSNP("/Users/marthaanitademba/Documents/amplicon_eniyou/bwaSingleEndVcfsFreebayes", "Output.g5mac3dp3.recodeSnp.Cleaned.txt")

