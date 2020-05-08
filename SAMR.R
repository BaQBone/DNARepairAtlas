library(samr)
library(dplyr)
library(reshape)
#library(reshape2)
library(ggplot2)
library(vsn)
library(limma)
library(smoothmest)
library(Biobase)
library(MSnbase)
library(matrixStats)
library(gplots)
library(fdrtool)
library(purrr)
library(LSD)
library(labeling)
library(munspace)
library(taRifx)
library(plotly)

#library(org.Sc.sgd.db)
library(stringr)
#install.packages("dplyr")
#install.packages("reshape")
#library(reshape2)
install.packages("ggplot2")
BiocManager::install("vsn")
BiocManager::install("limma")
BiocManager::install("smoothmest")
BiocManager::install("Biobase")
BiocManager::install("MSnbase")
BiocManager::install("matrixStats")
BiocManager::install("gplots")
BiocManager::install("fdrtool")
BiocManager::install("purrr")
BiocManager::install("LSD")
BiocManager::install("org.Sc.sgd.db")
BiocManager::install("stringr")
BiocManager::install("labeling")
BiocManager::install("taRifx")


##############################################################################################
########################################## FUNCTIONS #########################################
##############################################################################################

#Function Imputation: replace missing values by normal distribution
#-----------------------------------------------------------------

Imputation <- function(x, shift=1.8, width=0.3){
  
  temp <- c()
  inf.values <- is.na(x)
  m <- mean(x[!inf.values])
  s <- sd(x[!inf.values])
  imputed <- rnorm(length(grep(T, is.na(x))), mean = (m-shift*s), sd = (s*width))
  x[inf.values] <- imputed
  if (length(temp) == 0){
    temp <- x
  }else{
    temp <- cbind(temp, x)
  }
  return(temp)
}

##############################################################################################


# Data Import
#------------------
setwd("R:\\MS-Data\\DNA REPAIR ATLAS\\ATLAS INPUT") 


# Import of the Metadata
#---------------------------

conditions <- read.table("CHROMASS - Experimental Conditions V5.txt",
                         header = T,
                         sep = "\t",
                         fill = T,
                         stringsAsFactors = F,
                         comment.char = "")

name.rep <- paste(conditions$Experiment.Series,conditions$Treatment,conditions$Time,sep=".")

# Import of the LFQ data
#-------------------------

LFQ <- read.table("Repair_Atlas_LFQ_Intensities_flipped_V7.txt", header = T, sep = "\t", fill = T, stringsAsFactors = F, comment.char = "")


# Transpose LFQ and map column labels from MetaData
#---------------------------------------------------------

LFQ <- as.data.frame(t(LFQ))

# Generate a lookup table to retrieve Experimental condition from MetaData with FILE_ID
columnHeader <- conditions$Group
names(columnHeader) <- paste("id_",conditions$FILE_ID,sep="")   


# Map column name onto LFQ (Remember: last row of LFQ contains the FILE_ID)
names(LFQ) <- columnHeader[paste("id_",LFQ[nrow(LFQ),],sep="")]

# remove last row from LFQ containing the FILE_ID
LFQ<-LFQ[1:nrow(LFQ)-1,]     

# Fix some column names:
colnames(LFQ)[grep("EXP08.NCA.30_01", colnames(LFQ))] <- "EXP07.NCA.30_03"
colnames(LFQ)[grep("EXP08.NCA.60_01", colnames(LFQ))] <- "EXP07.NCA.60_01"
colnames(LFQ)[grep("EXP08.NCA.60_02", colnames(LFQ))] <- "EXP07.NCA.60_02"
colnames(LFQ)[grep("EXP08.NCA.60_03", colnames(LFQ))] <- "EXP07.NCA.60_03"



groups <- unique(sub("_[0-9]{2}","",colnames(LFQ)))
write.table(groups,
            "RepairAtlas name of replicates.txt", sep="\t",
            col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA") 

write.table(LFQ,
            "LFQ_temp.txt", sep="\t",
            col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA") 




# Generate tables to append the reults of the T-Tests and imputed LFQ values
#-----------------------------------------------------------------------------
results <- LFQ[,1:2]
results[,1] <- 1:nrow(LFQ)
results[,2] <- rownames(LFQ)
colnames(results) <- c("prot.id", "prot.name")
rownames(results) <- results$prot.id

lfq.imp  <- results




# Generate matrix for all the T-Test comparisons
#-------------------------------------------------------
#tTest <- matrix(c("LFQ.intensity.CAP","LFQ.intensity.CTR", "TRUE",
#                 "LFQ.intensity.LIN","LFQ.intensity.CTR", "TRUE",
#                 "LFQ.intensity.PHS","LFQ.intensity.CTR", "TRUE"),
#               nrow  = 3,
#               ncol  = 3,
#               byrow = T)
#colnames(tTest) <- c("right_side","left_side","impute_left" )




tTest <- read.table("comparisons.test.txt",
                         header = T,
                         sep = "\t",
                         fill = T,
                         stringsAsFactors = F,
                         comment.char = "")



setwd("R:\\MS-Data\\DNA REPAIR ATLAS\\ATLAS INPUT\\TTEST ANALYSIS\\RESULTS") 




for (i in 50:nrow(tTest)){

       Int <- cbind(results[,1:2],LFQ[,c(grep(tTest[i,1],colnames(LFQ)),
                                   grep(tTest[i,2],colnames(LFQ)))]) 


       # Filter for at least three valid values in the right replicate
       repRight.valid <- as.matrix(Int[,c(grep(tTest[i,2],colnames(Int)))])
       repRight.valid[which(repRight.valid != 0)] <- 1
       valid <- rowSums(repRight.valid)
       Int <- subset(Int, valid > 2)
       Int$valid <- NULL

      
       # Impute values into the replicates of the "left side" of the volcano plot
       if (tTest[i,3]){Int[,grep(tTest[i,1],names(Int))] <- Imputation(Int[,grep(tTest[i,1],names(Int))])}
       lfq.imp <- merge(lfq.imp,Int, by.x = "prot.id", by.y = "prot.id", all.x = T )
       
 
       # Define variables for the T-Tests
       desc <- ""
       min.fold <- 1.5
       max.fdr  <- 0.05
       plots <- T
       
       x<-as.matrix(Int[,c(grep(tTest[i,1],names(Int)),
                           grep(tTest[i,2],names(Int)))])
       
       
       rep <- unique(sub("(.*)_[0-9]{2}","\\1",colnames(x)))
       r1  <- grep(rep[1], colnames(x))
       r2  <- grep(rep[2], colnames(x))
       
       test_desc <- paste(tTest[i,2],"vs",tTest[i,1],sep="") 
       test_desc
       input = list(
          x = x,
          y = c(rep(1, length(grep(tTest[i,1],colnames(Int)))), rep(2, length(grep(tTest[i,2],colnames(Int))))),
          geneid = Int$prot.name,
          genenames = Int$prot$id,
          logged2 = T)

        samr.obj <- samr(input,
                         resp.type = "Two class unpaired",#,center.arrays=TRUE)  if arrays should be normalized
                         nperms = 1000)
   
        output <- Int[,1:2]
        output$DIFF      <- apply(x,1,function(x){mean(x[r2]) - mean(x[r1])})
        output$PVAL      <- unlist(apply(x, 1, function(x) {-log10(t.test(
                                                            x = x[r1],
                                                            y = x[r2],
                                                            var.equal = T)$p.value)}))

        colnames(output)[3] <-  paste(test_desc, "_DIFF", sep="")
        colnames(output)[4] <-  paste(test_desc, "_PVAL", sep="")  
        
        write.table(output,
                    "output.txt", sep="\t",
                    col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA")        
        
        qplot(output[,3],output[,4])
      
        delta.table <- samr.compute.delta.table(samr.obj, min.foldchange = min.fold, dels = c((1:200)/100))
        s0 <- round(samr.obj$s0,2)
        write.table(delta.table,
                  paste("Delta_table_", test_desc, "_","MFC_",min.fold,"_S0=",s0,".txt", sep=""), sep="\t",
                  col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA")
      
    
        if (min(na.omit(delta.table[,5])) <= max.fdr){
            
            # Retrieve the highest Delta so that fdr < max.fdr
            j=1
            while(na.omit(delta.table[j,5]) > max.fdr) {j<-j+1}
            delta<-delta.table[j,1]

            # compute Significant genes 
            siggenes     <-samr.compute.siggenes.table(samr.obj, del = delta, input, delta.table, min.foldchange = min.fold)  
            siggenes_all <-samr.compute.siggenes.table(samr.obj, del = delta, input, delta.table, min.foldchange = min.fold, all.genes=T)
            all_genes    <-rbind(siggenes_all$genes.up,siggenes_all$genes.lo) 

            
            # Prepare Protein ID lookup table
            # ----------------------------------            
            prot.id <- Int$prot.id
            names(prot.id) <- paste("g",1:nrow(Int), sep="")

                    
            # Printing out tables  

            
            
            if (siggenes$ngenes.up > 0){
            genes_up<-siggenes$genes.up
            genes_up<-data.frame(genes_up)
            genes_up$UP <- 1
            colnames(genes_up)[ncol(genes_up)] <- paste(test_desc, "_UP", sep="")
            genes_up$prot.id <- prot.id[paste(genes_up$Gene.ID, sep=",")]

            write.table(genes_up,
                        paste("Genes_up_",test_desc,".txt", sep=""), sep="\t",
                        col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA")}
            
            
            if (siggenes$ngenes.lo > 0){
            genes_lo<-siggenes$genes.lo
            genes_lo<-data.frame(genes_lo)
            genes_lo$DOWN <- 1
            colnames(genes_lo)[ncol(genes_lo)] <- paste(test_desc, "_DOWN", sep="")
            genes_lo$prot.id <- prot.id[paste(genes_lo$Gene.ID, sep=",")]           

            write.table(data.frame(siggenes$genes.lo),
                        paste("Genes_lo_",test_desc,".txt", sep=""), sep="\t",
                        col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA")}
        
            write.table(data.frame(all_genes),
                        paste("Genes_all_",test_desc,".txt", sep=""), sep="\t",
                        col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA") 
        
            # Plot
            if (plots == T){
            file.name <- paste("SAM_Plots_", test_desc, "_","MFC_",min.fold,".pdf", sep="")
            pdf(file.name)
            samr.plot(samr.obj, del=delta, min.foldchange=min.fold)
            title(main=paste("SAM PLOT: ",
                             test_desc,
                             "  (FDR=",
                             max.fdr,
                             "; s0=",
                             round(samr.obj$s0, digits = 2),
                             ")",sep="" ))
             delta.frame <- data.frame(delta.table)
             plot(delta.frame$delta, delta.frame$median.FDR, cex = 0.2, xlab="Delta", ylab = "Median FDR")
             title(main= "FDR PLOT", sub= paste("UP:",siggenes[[4]],", DOWN:",siggenes[[5]],", False Positives=",round(delta.table[i,2]),sep=""))
             points(delta.frame$delta[i],delta.frame$median.FDR[i], col="red", cex=0.6)
             text(delta.frame$delta[i],delta.frame$median.FDR[i], paste("Delta:",delta.table[i,1],", Fold Change:",min.fold,", FDR:",max.fdr,sep=""), col="red",pos=4, cex=0.75)
             dev.off()}
        
            
            if (siggenes$ngenes.up > 0 && siggenes$ngenes.lo > 0){
                output <- merge(output, genes_up[,c(4,9:10)], by.x = "prot.id", by.y= "prot.id", all.x = T)
                output <- merge(output, genes_lo[,c(4,9:10)], by.x = "prot.id", by.y= "prot.id", all.x = T)              

                colnames(output)[5] <- paste(test_desc, "_SCORE_UP", sep="")
                colnames(output)[7] <- paste(test_desc, "_SCORE_DOWN", sep="")

                output <- remove.factors(output)                             # requires library taRifix
                output[,5:8][is.na(output[,5:8])] <- 0
                output[,3]<-round(as.numeric(output[,3]), digits=3) 
                output[,4]<-round(as.numeric(output[,4]), digits=3) 
                output[,5]<-round(as.numeric(output[,5]), digits=3)  
                output[,7]<-round(as.numeric(output[,7]), digits=3)}                
                              
              
            if (siggenes$ngenes.up > 0 && siggenes$ngenes.lo == 0){
                output <- merge(output, genes_up[,c(4,9:10)], by.x = "prot.id", by.y= "prot.id", all.x = T)

                output$SCORE_DOWN <- 0
                output$DOWN <- 0
                
                colnames(output)[5] <- paste(test_desc, "_SCORE_UP", sep="")
                colnames(output)[7] <- paste(test_desc, "_SCORE_DOWN", sep="")   
                colnames(output)[8] <- paste(test_desc, "_DOWN", sep="")

                output <- remove.factors(output)                             # requires library taRifix
                output[,5:8][is.na(output[,5:8])]<-0}                
               
                
            if (siggenes$ngenes.up == 0 && siggenes$ngenes.lo > 0){
                output$SCORE_UP <- 0
                output$UP <- 0
                output <- merge(output, genes_lo[,c(4,9:10)], by.x = "prot.id", by.y= "prot.id", all.x = T)

                colnames(output)[5] <- paste(test_desc, "_SCORE_UP", sep="")
                colnames(output)[6] <- paste(test_desc, "_UP", sep="")
                colnames(output)[8] <- paste(test_desc, "_DOWN", sep="")
                
                output <- remove.factors(output)                             # requires library taRifix
                output[,5:8][is.na(output[,5:8])]<-0}               

            
            if (siggenes$ngenes.up == 0 && siggenes$ngenes.lo == 0){
                output$SCORE_UP <- 0
                output$UP <- 0
                output$SCORE_DOWN <- 0
                output$DOWN <- 0
                
                colnames(output)[5] <- paste(test_desc, "_SCORE_UP", sep="")               
                colnames(output)[6] <- paste(test_desc, "_UP", sep="")
                colnames(output)[7] <- paste(test_desc, "_SCORE_DOWN", sep="")
                colnames(output)[8] <- paste(test_desc, "_DOWN", sep="")
               
               output <- remove.factors(output)                             # requires library taRifix
               output[,5:8][is.na(output[,5:8])]<-0}               
        }  

        if (min(na.omit(delta.table[,5]))>max.fdr){
        output$SCORE_UP <- 0  
        output$UP <- 0    
        output$SCORE_DOWN <- 0  
        output$DOWN <- 0           
        colnames(output)[5] <- paste(test_desc, "_SCORE_UP", sep="")               
        colnames(output)[6] <- paste(test_desc, "_UP", sep="")
        colnames(output)[7] <- paste(test_desc, "_SCORE_DOWN", sep="")
        colnames(output)[8] <- paste(test_desc, "_DOWN", sep="")
        output <- remove.factors(output)                                  # requires library taRifix
        output[,5:8][is.na(output[,5:8])]<-0
        }   
   

        output$prot.name <- NULL
        results <- merge(results, output, by.x = "prot.id", by.y = "prot.id", all.x = T)}
        

results[,3:ncol(results)][is.na(results[,3:ncol(results)])] <- 0

write.table(results,
            "Repair_Atlas_Results.txt", sep="\t",
            col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA")


# Volcano
res<- results[,c(1:4,6)]
res<- filter(res, res[,4] != 0)
res$SIG <- as.factor(res$SIG)
colnames(res)<-c("id","name","DIFF","PVAL","SIG")
p <- ggplot(data=res, aes(x=DIFF, y= PVAL, text=name, color=SIG )) +
  geom_point(alpha=0.3, size=1.5) +
  xlab("Log2 fold change") + 
  ylab("-Log10 p-value") +
  ggtitle("Psoralen +/- Geminin, 45 Min")
ggplotly(p)




results.brief <- results[,1:2]

prot.sig <- results[,c(1:2,seq(from =6, to = nrow(tTest)*6 , by =6))]
prot.sig$COUNT_UP <- apply(prot.sig[,3:ncol(prot.sig)],1,sum)
results.brief$COUNT_UP  <- prot.sig$COUNT_UP

# Number of proteins that score significant in at least one condition: (optimized S0)
nrow(subset(prot.sig, prot.sig$COUNT_UP >0))

# Calculate the protduct of all p-values that are associated with upregulated proteins
prot.pval <- results[,seq(from =6, to = nrow(tTest)*6 , by =6)-2]
prot.pval <- 10^-prot.pval

prot.diff <- results[,seq(from =6, to = nrow(tTest)*6 , by =6)-3]
prot.pval[prot.diff < 0] <- 1
prot.pval$prod <- -log10(apply(prot.pval,1,prod))
results.brief$PVAL_PROD <- prot.pval$prod

# Add tables with individual counts, pval-products and scores 

prot.pval <- results[,seq(from =6, to = nrow(tTest)*6 , by =6)-2]
prot.score <-results[,seq(from =6, to = nrow(tTest)*6 , by =6)-1]

results.brief <- cbind(results.brief, prot.sig)
results.brief <- cbind(results.brief, prot.pval)
results.brief <- cbind(results.brief, prot.score)

write.table(results.brief,
            "Repair_Atlas_Results.brief.txt", sep="\t",
            col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA") 











       
        
        cols <- sapply(results, is.logical)
results[,cols] <- lapply(results[,cols], as.numeric)



write.table(lfq.imp,
            "Repair_Atlas_Normalized_LFQ.txt", sep="\t",
            col.names=TRUE, row.names=FALSE, quote=FALSE, na="NA")


               


     


#################
###### APPEND function
#################


Append.df<- function(df,data=mydata){
  temp <- as.data.frame(data)
  temp$id <- rownames(temp)
  temp <- merge(df, temp, by.x = "id", by.y = "id", all.x = T)
}






