#
#figure1
#
# (c) Martha Anita Demba (MRC at LSHTM)
# Date: 2022-07-05
#
# Purpose:
# Calculate genetic distances
# Generate heatmap for genetic distance results
# ==========================================================


genetic.dist=function(data, method){
  library(poppr)
  if(method=="bruvo")
    hh=as.matrix(bruvo.dist(data))
  else if (method=="provesti")
    hh=as.matrix(provesti.dist(data))
  else if (method=="nei")
    hh=as.matrix(nei.dist(data, warning = T))
  else{  
  message("Error: Please use \"bruvo\",  \"provesti\" or \"nei\".")
    stop()
  }
}


complexPheatmap= function(hdist, pop2, type){
  library(ComplexHeatmap)
  ha = HeatmapAnnotation(
    Region = pop2$Region,
    Country = pop2$Location, 
    col = list(Region = c("WAF"="#1f9463","CAF"="#00ff00","EAF"="#ff0000","ASIA"="#ffccff"),
               Country = c("Nigeria"="#9003FE","Mali"="#FEFD03","Guinea"="#3258F2","Burkina_Faso"="#FA3D15","Ghana"="#14F6E3","Cameroon"="#6FDF31","Tanzania"="#900C3F","Asia"="#000000")
    ), gp = gpar(col = "black", lwd = 0.5), show_legend = FALSE
  )
  rha= rowAnnotation(
    Region = pop2$Region,
    Country = pop2$Location, 
    col = list(Region = c("WAF"="#1f9463","CAF"="#00ff00","EAF"="#ff0000","ASIA"="#ffccff"),
               Country = c("Nigeria"="#9003FE","Mali"="#FEFD03","Guinea"="#3258F2","Burkina_Faso"="#FA3D15","Ghana"="#14F6E3","Cameroon"="#6FDF31","Tanzania"="#900C3F","Asia"="#000000")
    ), gp = gpar(col = "black", lwd = 0.5)
  )
  #pdf("Genetic_distance_SNP.pdf", width = 8)
  colnames(hdist)=NULL
  rownames(hdist)=NULL
  Heatmap(hdist, col = c( "white","blue","yellow", "red"), rect_gp = gpar(col = "black", lwd = 0.3),
          column_title = paste0( type, "Genetic distance"), name = "Distance", top_annotation = ha, left_annotation = rha)
  #dev.off()
  
}

Pheatmap=function(hdist, pop2, type){
  library(pheatmap)
  newmeta=pop2[c(4,3)]
  newmeta[,1] = factor(newmeta[,1])
  newmeta[,2] = factor(newmeta[,2])
  rownames(newmeta)=colnames(hdist)
  Location <- c(Nigeria="#9003FE",Mali="#FEFD03",Guinea="#3258F2",Burkina_Faso="#FA3D15",Ghana="#14F6E3",Cameroon="#6FDF31",Tanzania="#900C3F",Asia="#000000")
  Region <- c(WAF="#1f9463",CAF="#00ff00",EAF="#ff0000",ASIA="#ffccff")
  annotation_colors <- list(Region=Region, Location=Location)
  
  pheatmap::pheatmap(mat = hdist, 
           color = colorRampPalette(c( "white","blue","yellow", "red"))(100), 
           show_colnames     =  FALSE,
           show_rownames     =  FALSE,
           annotation_col    = newmeta,
           annotation_row    = newmeta,
           annotation_colors = annotation_colors,
           drop_levels       = TRUE,
           fontsize          = 14,
           main              = paste0(type, " Genetic distance"), cex=1
  )
}


GeneticDist = function(dir, file, popData, type)
{
  library(data.table)
  library(adegenet)
  setwd(dir)
  x = as.data.frame(fread(file), header=T) #read Genotype file
  pop = read.table(popData, header = T) #read pop data       "../popDataPM.txt"
  
  #Order populations
  country=c("Nigeria", "Mali", "Guinea", "Burkina_Faso", "Ghana", "Cameroon", "Tanzania", "Asia")
  pop2=NULL
  for (j in country) 
    pop2 = rbind.data.frame(pop2, pop[which(pop$Location == j),])
  
  pop2$Sample_ID=gsub("/", "_", pop2$Sample_ID)
  # #set annotation colours
  # Location <- c(Nigeria="#9003FE",Mali="#FEFD03",Guinea="#3258F2",Burkina_Faso="#FA3D15",Ghana="#14F6E3",Cameroon="#6FDF31",Tanzania="#900C3F",Asia="#000000")
  # Region <- c(WAF="#1f9463",CAF="#00ff00",EAF="#ff0000",ASIA="#ffccff")
  # #Create population data for heatmap
  # newmeta=pop2[c(4,3)]
  # newmeta[,1] = factor(newmeta[,1])
  # newmeta[,2] = factor(newmeta[,2])
  # rownames(newmeta)=pop2$Barcode
  
  if(type == "msat")
  {
    x=as.data.frame(t(x))
    info = x[1:3,]
    xx=as.data.frame(t(x[-1:-2,-1:-3]))[-6:-9]
    #colnames(xx) = as.vector(x[3:7, 3])
    rownames(xx) = as.vector(x[1,4:ncol(x)])
    pop3 = pop2[which(pop2$Sample_ID %in% rownames(xx)),]
    xx1=xx[match(pop3$Sample_ID, rownames(xx)),] #reorder genotype data to match population data
    p2 = df2genind(xx1, ploidy = 1, pop=as.factor(pop3$Location)) # convert dataframe to genind
    method="bruvo" #preferred method either "provesti" or "bruvo"
    hdist = genetic.dist(p2, method) #calculate genetic distance. choose brovo dist or provesti for microsatellites e.g hdist = genetic.dist(p2, "provesti")
    pdf(paste0("msatGeneticDistance_", method, ".pdf"), width = 10, height = 8)
        #Plot heatmap with complex heatmap
        #print(complexPheatmap(hdist, pop2, type))
        #plot heatmap with pheatmap
        rownames(pop3)=colnames(hdist)
        print(Pheatmap(hdist, pop3, type))
        dev.off()
    
  }
  else if (type == "snp")
  {
    rownames(x) = gsub(" ", "_", x$V1)
    x=x[,-1]
    colnames(x) = x[1,]#pop$Sample_ID#samples$V1
    x = x[-1,]
    df2 = as.data.frame(t(x))
    #df2 = df2[which(rownames(df2) %in% pop$Barcode),]
    pop3 = pop2[which(pop2$Barcode %in% rownames(df2)),]
    #df3 = df3[which(rownames(x) %in% pop1$Barcode),]
    xx1=df2[match(pop3$Barcode, rownames(df2)),] #reorder genotype data to match population data
    p2 = df2genind(xx1, ploidy = 1, pop=as.factor(pop3$Location)) # convert dataframe to genind
    method="nei" #preferred method either "provesti" or "nei"
    hdist = genetic.dist(p2, method) #calculate genetic distance. choose brovo dist or provesti for microsatellites e.g hdist = genetic.dist(p2, "provesti")
    pdf(paste0("snpGeneticDistance_", method, ".pdf"), width = 10, height = 8)
        #Plot heatmap with complex heatmap
        #print(complexPheatmap(hdist, pop2, type))
        #plot heatmap with pheatmap
        print(Pheatmap(hdist, pop3, type))
        dev.off()
  }
  else
    stop()
}
GeneticDist("/Users/marthaanitademba/Documents/amplicon_eniyou/bwaSingleEndVcfsFreebayes", "Output.g5mac3dp3.recodeSnp.Cleaned.txt", "../popDataPM.txt", "snp") 
GeneticDist("/Users/marthaanitademba/Documents/amplicon_eniyou/bwaSingleEndVcfsFreebayes", "../Pm_msat_data.csv", "../popDataPM.txt", "msat") #msatData="../Pm_msat_data.csv"
#msatData="../Pm_msat_data.csv"
