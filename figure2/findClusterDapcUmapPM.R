#

# figure2

#

# (c) Martha Anita Demba (MRC at LSHTM)

# Date: 2022-07-07

#

# Purpose:

# find cluster, Discriminant analysis of principal component (DAPC) plus Uniform Manifold approximation and Projection (UMAP)

# Generate scatterplot for results.

# ==========================================================

dapc_umap = function(dir, popData, file, samplesFile)

{

  library(adegenet)

  library(ggplot2)

  library(data.table)

  library(splitstackshape)

  library(ggplot2)

  library(umap)

  

  setwd(dir) #set directory

  #read input files

  h = as.data.frame(fread(file, header = F))

  pop = read.table(popData, header = T)

  samples=fread(samplesFile, header = F)

  

  h1 = h[,which(samples$V1 %in% pop$Barcode)] #get common samples found in samples and population meta data

  hh = as.data.frame(t(h1)) #transpose data

  colnames(hh) = gsub(" ", "_", h$V1)

  #hh1 = cSplit(hh, colnames(hh), "/") # split by "/

  # hhgenind = df2genind(hh1, ploidy=1, ncode=1)

  hhgenind = df2genind(hh, ploidy=1, ncode=1) #convert dataframe to genind

  hhgenind$pop = as.factor(pop$Location)

  

  mycol = c("#000000", "#FA3D15", "#6FDF31", "#14F6E3", "#3258F2", "#FEFD03", "#9003FE", "#900C3F") #set colours for different populations in plot

  #find cluster groups

  grp <- find.clusters(hhgenind, max.n.clust=40, n.pca = 4, n.clust = 5) 

  #dapc with results from find cluster

  dapc1 <- dapc(hhgenind, grp$grp, n.pca = 4, n.da = 3) 

  #scatter(dapc1)

  

  #umap analysis on dapc results

  umap1 = umap(as.matrix(tab1), init="spca")

  points=as.data.frame(umap1$layout)

  colnames(points) = paste0("dim", 1:2)

  points$pop=as.factor(pop$Location)

  

  #plot results

  pdf("snpPmDapc_Umap1.pdf", width = 10)

  q=ggplot(points, aes(dim1, dim2)) + 

    geom_point(colour="Black", shape=21, size = 5, 

               aes(fill = factor(pop))) + 

    scale_fill_manual("Population",values=c("#000000", "#eb1f10", "#30ff08", "#05f0fc", "#0509fc", "#FEFD03", "#9003FE", "#900C3F"))+

    theme(text = element_text(size = 15, face = "bold"),  panel.grid.minor = element_blank(), 

          panel.background = element_blank(), axis.line = element_line(colour = "black")) #+theme_classic()

  print(q+scale_y_continuous(breaks=seq(min(round(points$dim2))-1, max(round(points$dim2))+1, by = 1))+scale_x_continuous(breaks=seq(min(round(points$dim1))-1, max(round(points$dim1))+1, by = 1)))

  dev.off()  

}

dapc_umap("/Users/marthaanitademba/Documents/amplicon_eniyou/bwaSingleEndVcfsFreebayes", "../popDataPM.txt", "Output.g5mac3dp3.recodeSnp.txt", "samples.txt")

