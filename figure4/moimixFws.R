#

#figure4

#

# (c) Martha Anita Demba (MRC at LSHTM)

# Date: 2022-07-05

#

# Purpose:

# Calculate fws

# Generate boxplot for fws results

# ==========================================================

calculateFWS = function(dir, pop, vcfFile){

  #load libraries

  library(SeqArray)

  library(moimix)

  library(ggplot2)

  library(ggpubr)

  library(ggsignif)

  

  setwd(dir) #set directory

  #read and order file containing information about the populations

  pop = read.table(pop, header = T)

  pop = pop[order(pop$Location),]

  region = c("WAF", "CAF", "EAF", "ASIA")

  #order by region

  pop2=NULL

  for (j in region) 

    pop2 = rbind.data.frame(pop2, pop[which(pop$Region == j),])

  

  gds_file=gsub("vcf", "gds", vcfFile)

  seqVCF2GDS(vcfFile, gds_file) #convert vcf to gds format

  isolates <- seqOpen(gds_file) #read gds file

  fws_all <- getFws(isolates) #calculate fws

  names(fws_all)=gsub("_mergedTrimmed.fastq.gz", "", names(fws_all)) 

  filtered_fws_all = fws_all[which(names(fws_all) %in% pop$Barcode)] #get loci that are present in both pop and fws_all

  #hist(filtered_fws_all) #histogram plot

  #order fws results tp match population data

  filtered_fws_all1 = filtered_fws_all[match(pop2$Barcode, names(filtered_fws_all))] 

  pop2$fws = filtered_fws_all1 #add fes results to population data

  pop2$Location <- factor(pop2$Location, levels = unique(pop2$Location))

  pop2$Region <- factor(pop2$Region, levels = unique(pop2$Region))

  #plot boxplot using ggplot2

  p = ggplot(pop2, aes(Location, fws, fill=factor(Region))) + geom_boxplot(position=position_dodge(1))

  

  #save plot in pdf

  pdf(gsub(".vcf", "Fws.pdf", vcfFile), width = 8)

  p+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),

          panel.background = element_blank(), axis.line = element_line(colour = "black"))+

    ggtitle("Fws per country") + ylim(0, 1) +

    guides(fill=guide_legend(title="Data type")) +

    # scale_fill_manual(values=c("#999999", "#c76b2e"))+

    theme(text = element_text(size = 10, face="bold"),

          axis.text.x = element_text(angle=45, hjust=1)) + 

    stat_compare_means(label = "p.signif", method = "t.test",label.y = 1.01, paired = T)

  dev.off()  

}

calculateFWS("/Users/marthaanitademba/Documents/amplicon_eniyou/bwaSingleEndVcfsFreebayes", "../popDataPM.txt", "Output.g5mac3dp3.recodeSnp.vcf")

