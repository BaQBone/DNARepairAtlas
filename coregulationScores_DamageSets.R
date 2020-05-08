library(ggplot2)
library(dplyr)
library(data.table)
library(WGCNA)
library(igraph)
library(Hmisc)
library(tidyverse)
library(ggpubr)
library(RcmdrMisc)
library(WGCNA)
library(brainwaver)
library(intergraph)
library(PerseusR)

set.seed(42)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("./ZScored")

dpc<-read.table("DPC_Imp.txt", header=T, sep="\t", check.names=F)
dsb<-read.table("DSB_Imp.txt", header=T, sep="\t", check.names=F)
icl<-read.table("ICL_Imp.txt", header=T, sep="\t", check.names=F)
fc<-read.table("FC_Imp.txt", header=T, sep="\t", check.names=F)
rt<-read.table("RT_Imp.txt", header=T, sep="\t", check.names=F)


doWGCNA<-function(df, name, adj = "power", m = 0.91, a = 37){
  
  options(warn=-1)
  
  setwd("C:\\Users\\BaQBone\\Desktop\\filterOptimization\\ZScored\\rawSigSum")
  
  #Filtering for significance
  filtered<-df
  filtered$Sig<-NULL
  filtered[is.na(filtered)] <- 0 #IMPUTATION OF MNAR DATA:replace NA with 0 to imprve correlation performance       MNAR=Missing not at random
  rownames(filtered)<-filtered$Name
  filtered$Name<-NULL
  filtered<-t(filtered)
  
  #Calculating correlation matrix
  flat_cor_mat <- function(cor_r, cor_p){
    library(tidyr)
    library(tibble)
    cor_r <- rownames_to_column(as.data.frame(cor_r), var = "row")
    cor_r <- gather(cor_r, column, cor, -1)
    cor_p <- rownames_to_column(as.data.frame(cor_p), var = "row")
    cor_p <- gather(cor_p, column, p, -1)
    cor_p_matrix <- left_join(cor_r, cor_p, by = c("row", "column"))
    cor_p_matrix    
  }
  input<-as.matrix(filtered)
  corr.matrix<-rcorr(input, type="pearson")
  corr.vals<-as.matrix(corr.matrix$r)
  corr.p<-corr.matrix$P
  rm(corr.matrix)
  flattened<-flat_cor_mat(corr.vals,corr.p)
  out<-filter(flattened, p <= .01)
  out$p<-NULL
  names(out)[names(out) == "row"] <- "Source" 
  names(out)[names(out) == "column"] <- "Target"
  names(out)[names(out) == "cor"] <- "Weight"
  out <- dcast(out, Source ~ Target, value.var = "Weight")
  corr_data <- as.matrix( data.frame(out, row.names = "Source"))
  
  #Calculating the adjacency matrix using the sigmoid adjacency function
  adj_mat <- sigmoidAdjacencyFunction(corr_data, mu = m, alpha =a )
  
  #Get the Topological Overlap Matrix
  tom_sim <- TOMsimilarity( adj_mat, TOMDenom = "mean" )
  colnames(tom_sim) <- colnames(corr_data)
  row.names(tom_sim) <- colnames(corr_data)
  
  ## Test if network is scale-free
  connectivity <- colSums(tom_sim, na.rm = TRUE) - 1
  png(paste(name,"_scaleFreeness.png", sep = ""))
  scaleFreePlot(connectivity)
  dev.off()
  
  #Turning similarity matrices into long-format, remove duplicates and merge into final table
  tc_sim_dt <- as.data.table( melt( corr_data ))
  tc_sim_dt <- tc_sim_dt[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tc_sim = value ) ]
  tc_sim_dt <- tc_sim_dt[ Protein_1 > Protein_2 ]                  # Removes duplicate pairs (incl. self-comparisons)
  
  adj_mat_dt <- as.data.table( melt( adj_mat ))
  adj_mat_dt <- adj_mat_dt[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tc_pow = value ) ]
  adj_mat_dt <- adj_mat_dt[ Protein_1 > Protein_2 ]                # Removes duplicate pairs (incl. self-comparisons)
  
  tc_tom_dt <- as.data.table( melt( tom_sim ))
  tc_tom_dt <- tc_tom_dt[, .( Protein_1 = as.character(Var1), Protein_2 = as.character(Var2), tc_tom = value ) ]
  tc_tom_dt <- tc_tom_dt[ Protein_1 > Protein_2 ]                  # Removes duplicate pairs (incl. self-comparisons)
  
  tc_dt <- merge( tc_sim_dt, adj_mat_dt, by = c("Protein_1", "Protein_2"))
  tc_dt <- merge( tc_dt, tc_tom_dt,      by = c("Protein_1", "Protein_2"))
  
  names(tc_dt)[names(tc_dt) == "Protein_1"] <- "Source"
  names(tc_dt)[names(tc_dt) == "Protein_2"] <- "Target"
  names(tc_dt)[names(tc_dt) == "tc_sim"] <- "SIM"
  names(tc_dt)[names(tc_dt) == "tc_pow"] <- "ADJ"
  names(tc_dt)[names(tc_dt) == "tc_tom"] <- "TOM"
  
  tc_dt<-tc_dt[!(is.na(SIM))]
  tc_dt<-tc_dt[!(SIM <= 0)]
  
  message("Done.")
  setwd("./WGCNA/")
  
  tc_dt$SIM<-NULL
  tc_dt$ADJ<-NULL
  
  fwrite(tc_dt, paste(name,"_coregulationScores.txt", sep = ""), sep = "\t")
  
  message("Cleaning up...")
  options(warn=0)
  message("Resetting Working Directory")
  setwd("C:\\Users\\BaQBone\\Desktop\\filterOptimization\\ZScored\\rawSigSum")
  message("...Done.")
  
  return(tc_dt)
  
}

all<-doWGCNA(all, name = "allSets")
dpc<-doWGCNA(dpc, name = "DPC")
dsb<-doWGCNA(dsb, name = "DSB")
fc<-doWGCNA(fc, name = "FC")
icl<-doWGCNA(icl, name = "ICL")
rt<-doWGCNA(rt, name = "RT")