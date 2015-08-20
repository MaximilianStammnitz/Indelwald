############ INDELWALD - HYBRID INDEL CALLING ##############
####################### Version 1.0 ########################

## STRINGENT FILTERING FUNCTIONS
## Last Update - 14/08/2015 ##
## mrs72 / Maximilian Stammnitz ##

# 1. Extracting the Contig from GL
From_GL <- function(List, GLname){
  tmp <- which(List[,2]==GLname)
  as.character(List[tmp,1])
}


# 2. Extracting the GL from Contig
From_Contig <- function(List, Contig){
  tmp <- which(List[,1]==Contig)
  as.character(List[tmp,2])
}


# 3. Compare the Raw Calls of Pindel and Platypus
raw.calls <- function(pindel.vcf, platypus.vcf, 
                      chromosome, sample){
  
  # Required Packages
  require(VennDiagram)
  
  # Output folder
  setwd(paste0(main.path,"/Toy/Output"))
          
  
  # a. Merge Contig and Positions for obtaining Unique Identifiers
  pindel.vcf[,1] <- paste(pindel.vcf[,1], pindel.vcf[,2], sep="_")
  platypus.vcf[,1] <- paste(platypus.vcf[,1], platypus.vcf[,2], sep="_")
  raw.union <- length(union(pindel.vcf[,1], platypus.vcf[,1]))
  
  
  # b. Total number of RAW calls
  full.pindel <- length(pindel.vcf[,"Contig"]) 
  full.platypus <- length(platypus.vcf[,"Contig"])
  intersection <- intersect(pindel.vcf[,"Contig"],platypus.vcf[,"Contig"])
  intersection <- length(intersection)
    
  
  # c. Draw Venn between Pindel and Platypus Calls
  png(paste0("Intersection_", sample, "_chr", chromosome, ".png"), width=400,height=350)
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 5, 0, 0))
  draw.pairwise.venn(area1=full.pindel,
                     area2=full.platypus,
                     cross.area=intersection,
                     category = c(paste0("Pindel"),
                                  paste0("Platypus")),
                     fontfamily = rep("sans", 3), cat.fontfamily = "sans",
                     fill = c("blue", "orange"), euler.d=TRUE, scaled=TRUE, 
                     ind = T, cex=1.3, cat.cex=0.9, cat.dist=0.03)
  dev.off()
    
  
  # d. Display level of concordance
  concordance <- round(intersection/raw.union, 3)*100
  cat("\n Sample:\t", sample, "\n Concordance:\t", concordance, "%")
  
  
  # e. Restore File Path
  setwd(main.path)
  
}


# 4. Compare the Raw Sizes of Pindel and Platypus Calls
raw.sizes <- function(pindel.vcf, platypus.vcf,
                      chromosome, sample){

  
  # Output folder
  setwd(paste0(main.path,"/Toy/Output"))
  
  
  ## a. Length of all INDELs
  
  # Pindel
  pindel.ref.lenghts <- nchar(as.character(pindel.vcf[,3]))
  pindel.alt.lenghts <- nchar(as.character(pindel.vcf[,4]))
  pindel.alt.InsDels <- pindel.alt.lenghts-pindel.ref.lenghts
  
  # Platypus
  platypus.ref.lenghts <- nchar(as.character(platypus.vcf[,3]))
  platypus.alt.lenghts <- nchar(as.character(platypus.vcf[,4]))
  platypus.alt.InsDels <- platypus.alt.lenghts-platypus.ref.lenghts
  
  
  ## b. Plot Sensitivity-Curves

  pdf(paste0("Size_Plot_", sample, "_chr", chromosome, ".pdf"), width=12, height=4)
  par(mfrow=c(1,3),oma = c(0, 0, 3, 0))
  hist.pindel.sizes <- hist(pindel.alt.InsDels[which(pindel.alt.InsDels> -100 & pindel.alt.InsDels< 100,T)], 
                            plot=F, breaks=seq(f= -100, t=100, by=2))
  hist.platypus.sizes <- hist(platypus.alt.InsDels[which(platypus.alt.InsDels> -100 & platypus.alt.InsDels< 100,T)], 
                              plot=F, breaks=seq(f= -100, t=100, by=2))
  
  ## Add small Indel distributions
  options(warn=-1)
  plot(hist.pindel.sizes$counts~hist.pindel.sizes$mids, log="y", type="l",
       lwd=2, pch=16, col="blue", lty=1, main=NA, ylab="Frequency", xlab="Size [Bp]",
       ylim=c(1,100000), cex=0.7, cex.lab=1.3, cex.main=1.8, xlim=c(-100,100))
  lines(hist.platypus.sizes$counts~hist.platypus.sizes$mids, type="l",
        lwd=2, pch=16, cex=0.7, col="orange", lty=1)
  legend("topleft", legend=c("Platypus", "Pindel"), col=c("orange", "blue"), lty=c(1,1), lwd=c(2,2), cex=1.5)
  
  hist.pindel.sizes <- hist(pindel.alt.InsDels, 
                            plot=F, breaks=seq(f= round(min(pindel.alt.InsDels, -9900)-100, -2),
                                               t= round(max(pindel.alt.InsDels, 500)+100, -2), by=100))
  hist.platypus.sizes <- hist(platypus.alt.InsDels, 
                              plot=F, breaks=seq(f= round(min(platypus.alt.InsDels, -9900)-100, -2),
                                                 t = round(max(platypus.alt.InsDels, 500)+100, -2), by=100))
  
  ## Add large Deletion distributions
  plot(hist.pindel.sizes$counts~hist.pindel.sizes$mids, log="y", type="l",
       lwd=2, pch=16, col="blue", lty=1, main=NA, xlab="Size [Bp]",
       ylim=c(1,100000), cex=0.7, cex.lab=1.3, cex.main=1.8, ylab="Frequency",
       xlim=c(round(min(pindel.alt.InsDels, -9900)-100, -2), round(max(platypus.alt.InsDels, 500)+100, -2)))
  lines(hist.platypus.sizes$counts~hist.platypus.sizes$mids, type="l",
        lwd=2, pch=16, cex=0.7, col="orange", lty=1)
  options(warn=0)
  legend("topleft", legend=c("Platypus", "Pindel"), col=c("orange", "blue"), lty=c(1,1), lwd=c(2,2), cex=1.5)
  
  ## Add Size Boxplot
  boxplot(list(platypus.alt.InsDels, pindel.alt.InsDels),
          outline = T, names=c("Platypus", "cgpPindel"), outcol=c(rgb(1,0.5,0,0.05),rgb(0,0,1,0.05)),
          pch=rep(16,2), cex=1.2, ylab="Size [Bp]", xlab="Variant Caller", cex.lab=1.3)
  
  ## Add Title
  mtext(paste0("Indel Sizes: ", sample, " - Chromosome ", chromosome), 
        side = 3, outer = T, cex = 1.8)
  dev.off()
  
  
  # c. Restore File Path
  setwd(main.path)
  
}


# 5. Taking the Overlap of two different VCF Inputs
overlap.vcfs <- function(pindel.vcf, 
                         platypus.vcf){
  
  
  # a. Merge Contig & Position for memory-efficient Filtering
  pindel.vcf[,1] <- paste(pindel.vcf[,1], pindel.vcf[,2], sep="_")
  platypus.vcf[,1] <- paste(platypus.vcf[,1], platypus.vcf[,2], sep="_")
  

  # b. Remove multiple Calls on single Positions
  doubles.pindel <- duplicated(pindel.vcf[,"Contig"])
  if (any(doubles.pindel==T)){
    pindel.vcf <- pindel.vcf[-which(doubles.pindel==T),]
  }
  doubles.platypus <- duplicated(platypus.vcf[,"Contig"])
  if (any(doubles.platypus==T)){
    platypus.vcf <- platypus.vcf[-which(doubles.platypus==T),]    
  }
  
  
  # c. Build Overlap Sets
  both <- c()
  both <- intersect(pindel.vcf[,"Contig"], platypus.vcf[,"Contig"])
  pindel.vcf <- pindel.vcf[match(both,pindel.vcf[,"Contig"]),]
  platypus.vcf <- platypus.vcf[match(both,platypus.vcf[,"Contig"]),]
   
  
  # d. Change "Contig"-Column back to original state
  pindel.tmp <- matrix(NA, nrow = length(pindel.vcf[,"Contig"]), ncol = 4)
  platypus.tmp <- matrix(NA, nrow = length(platypus.vcf[,"Contig"]), ncol = 4)
  pindel.tmp <- str_split_fixed(pindel.vcf[,"Contig"],"_", 4)
  pindel.vcf[,1] <- paste(pindel.tmp[,1],pindel.tmp[,2],pindel.tmp[,3],sep="_")
  platypus.tmp <- str_split_fixed(platypus.vcf[,"Contig"],"_", 4)
  platypus.vcf[,1] <- paste(platypus.tmp[,1],platypus.tmp[,2],platypus.tmp[,3],sep="_")

  
  # e. Output
  return(list("Pindel" = pindel.vcf,
              "Platypus" = platypus.vcf))

}


# 6. Filtering out Calls with bad Pindel Quality
pindel.filter <- function(overlap, threshold){
  
  
  # a. Higher than the Threshold
  pindel.qualities <- overlap$"Pindel"[,"Quality"]
  pass <- pindel.qualities>threshold
  
  
  # b. Take Low-Q calls off the Pindel and Playpus Set
  if(length(which(pass==F))!=0){
    overlap$"Pindel" <- overlap$"Pindel"[-which(pass==F),]
    overlap$"Platypus" <- overlap$"Platypus"[-which(pass==F),]
  }
  
  # c. Output
  return(overlap)
}


# 7. Filtering out Calls with bad Platypus Quality
platypus.filter <- function(overlap, threshold){
  
  
  # a. Hard Q-Threshold  
  platypus.qualities <- overlap$"Platypus"[,"Quality"]
  high.quality.pass <- which(platypus.qualities>threshold)
  
  
  # b. Flag-Threshold: label "PASS" and 75 % Hard Q-Threshold
  soft.threshold <- 0.75*threshold
  platypus.filter <- as.character(overlap$"Platypus"[,"Filter"])
  filter.pass <- grep("PASS", platypus.filter)
  low.quality.pass.bottom <- which(platypus.qualities>soft.threshold)
  low.quality.pass.top <- which(platypus.qualities<threshold)
  low.quality.pass <- intersect(low.quality.pass.bottom, low.quality.pass.top)
  low.quality.pass <- intersect(filter.pass, low.quality.pass)

  
  # c. Concatenate the two sets
  pass <- c(low.quality.pass, high.quality.pass)
  pass <- sort(pass)
  
  
  # d. Take Low-Q calls off the Pindel and Playpus Set
  if(length(which(pass==T))!=0){
    overlap$"Pindel" <- overlap$"Pindel"[pass,]
    overlap$"Platypus" <- overlap$"Platypus"[pass,]
  }

  
  # e. Output
  return(overlap)
}


# 8. Summarise calls
indel.summary <- function(contigs.coordinates, 
                          platypus.vcf, 
                          pindel.vcf){
  
  # a. Pre-allocate a Matrix with two Columns, filled with the actual Positions
  hits <- matrix(NA, ncol = 7, nrow = length(platypus.vcf[,1]))
  colnames(hits) <- c("Contig", "Contig Position", "Chromosome Position", "Type", "Length", "Platypus_Calls", "Pindel_Calls")
  
  
  # b. Take Data from original Output
  hits[,"Contig"] <- as.character(platypus.vcf[,"Contig"])
  hits[,"Contig Position"] <- platypus.vcf[,"Position"]
  hits[,"Length"] <- nchar(as.character(platypus.vcf[,"Pattern"]))-nchar(as.character(platypus.vcf[,3]))
  hits[,"Platypus_Calls"] <- as.character(platypus.vcf[,"Tumour"])
  split <- matrix(NA, nrow = length(hits[,"Platypus_Calls"]), ncol = 6)
  split <- str_split_fixed(hits[,"Platypus_Calls"], ":", n=6)
  hits[,"Platypus_Calls"] <- paste(split[,6],split[,5], sep="/")
  hits[,"Pindel_Calls"] <- as.character(pindel.vcf[,"Tumour"])
  rm(split)
  split <- matrix(NA, nrow = length(hits[,"Pindel_Calls"]), ncol = 10)
  split <- str_split_fixed(hits[,"Pindel_Calls"], ":", n=10)
  hits[,"Pindel_Calls"] <- paste(as.numeric(split[,2])+
                                 as.numeric(split[,3]), 
                                 as.numeric(split[,8])+
                                 as.numeric(split[,9]), sep="/")

  
  # c. Exchange chromosomal for genomic Position
  
  # Take "000000-contig" into account
  if(length(grep("000000000",unique(hits[,"Contig"])[1]))!=0){
    starters <- which(hits[,"Contig"]==unique(hits[,"Contig"])[1])
    match.contigs.lengths <- match(hits[-starters,"Contig"],contig.sizes[,"Contig"])
    
    # Update Starters
    hits[starters,"Chromosome Position"] <- hits[starters,"Contig Position"]
    hits[-starters,"Chromosome Position"] <- as.numeric(hits[-starters,"Contig Position"])+
      contig.sizes[match.contigs.lengths-1,"Chromosome_Position"]
 
  }else{
    match.contigs.lengths <- match(hits[,"Contig"],contigs.coordinates[,"Contig"])
    hits[,"Chromosome Position"] <- as.numeric(hits[,"Contig Position"])+
      contigs.coordinates[match.contigs.lengths-1,"Chromosome_Position"]
  }
   
  
  # d. Define INDEL-type
  
  call.ref.lenghts <- nchar(as.character(platypus.vcf[,3]))
  call.alt.lenghts <- nchar(as.character(platypus.vcf[,4]))
  call.alt.InsDels <- call.alt.lenghts-call.ref.lenghts
  call.insertions <- which(call.alt.InsDels>0)
  call.deletions <- which(call.alt.InsDels<0)
  hits[call.insertions,"Type"] <- "I"
  hits[call.deletions,"Type"] <- "D"
  
  # e. Output
  return(hits)
}


# 9. Density Plot on Chromosome of interest
indel.density <- function(hits,contigs.coordinates, 
                          genes, sample, chromosome){
  
  
  # Output folder
  setwd(paste0(main.path,"/Toy/Output"))
  
  
  # a. Preparation
  full.chr.length <- contigs.coordinates[last(grep(paste0("Chr",chromosome), 
                                                   contigs.coordinates[,1])),"Chromosome_Position"]
  
  
  # b. Create Gene Annotations
  genes.coords <- cbind(rep(genes[,"Gene.Start"],2), rep(genes[,"Gene.End"],2))
  genes.coords <- apply(genes.coords, 1, function(x) (x[1]+x[2])/2)
  
  
  # c. Plotting
  pdf(paste0("Density_", sample, "_chr", chromosome, ".pdf"), width=15,height=8)
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 5, 0, 0)) 
  
  # Empty Chromosome String
  empty <- c(1, full.chr.length)
  plot(rep(1,2)~empty, type="l", lwd=7, ylab="Genes", xlab="Position [Bp]", cex.lab=1.7,
       main=paste0("Indel Density: ", sample, " - Chromosome ", chromosome), cex.main=1.8, 
       ylim=c(0,2), yaxt="n", col="white")
  
  ## Add genes and corresponding Annotations
  points(x=genes.coords, y=rep(1,length(genes.coords)), col=rgb(0.1,0.1,0.1,0.05), pch=3,
         cex=7)
  
  # All Insertions
  ins.pos <- as.integer(hits[which(hits[,"Type"]=="I"),"Chromosome Position"])
  ins.height <- rep(1.3, length(ins.pos))
  points(x=ins.pos, y=ins.height, cex=2, pch=18, col=rgb(1,0,0,0.01))
  
  # All Deletions: increased transparency - to facilitate relative comparisons
  # with insertions (see alpha.adjust)
  dels.pos <- as.integer(hits[which(hits[,"Type"]=="D"),"Chromosome Position"])
  dels.height <- rep(0.7, length(dels.pos))
  alpha.adjust <- length(ins.height)/length(dels.height)
  points(x=dels.pos, y=dels.height, cex=2, pch=18, col=rgb(0,0,1,0.01*alpha.adjust))
  
  # Legend
  legend("topleft", legend=c("Insertions", "Deletions"), col=c(rgb(1,0,0,0.8), rgb(0,0,1,0.8)), 
         pch=18, lwd=2, cex=2, pt.cex=2, bty="n")
  dev.off()
  
  
  # d. Restore File Path
  setwd(main.path)
  
}


# 10. Rainfall Plot on Chromosome
rainfall <- function(hits, sample, chromosome, 
                     contigs.coordinates = contig.lengths,
                     split){  
  
  
  # Output folder
  setwd(paste0(main.path,"/Toy/Output"))
  
  
  # a. Take all Distances
  names <- colnames(hits)
  hits <- cbind(hits, c(0,diff(as.integer(hits[,"Chromosome Position"]))))
  colnames(hits) <- c(names, "Indel Distances")
  hits <- hits[,c("Chromosome Position", "Indel Distances", "Type")]
  full.chr.length <- contigs.coordinates[last(grep(paste0("Chr",chromosome), 
                                                   contigs.coordinates[,1])),"Chromosome_Position"]
  
  
  # b. Rainfall Plots
  pdf(paste0("Rainfall_", sample, "_chr", chromosome, ".pdf"), width=18,height=12)
  par(mfrow=c(2,1))
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 5, 0, 0)) 
  
  # Shared Rainfall Plots
  
  ## Add Shared, relative distance
  plot(x=which(hits[,"Type"]=="D")[-1], 
       y=as.numeric(hits[which(hits[,"Type"]=="D"),2])[-1], ylab="Intermutation Distances", 
       xlab="Indel Number", pch=16, cex=0.6, log="y", 
       main=paste0("Relative-Distance Rainfall: ", sample, " - Chromosome ", chromosome), 
       cex.main=1.8, cex.lab=1.3,yaxt="n", col=rgb(0,0,1,0.5), 
       xlim=c(0,length(hits[,1])), ylim=c(10,1000000))
  points(x=which(hits[,"Type"]=="I")[-1], 
         y=as.numeric(hits[which(hits[,"Type"]=="I"),2])[-1],
         pch=16, cex=0.6, col=rgb(1,0,0,0.5))
  
  # Legend and X-axis
  axis(2, at=c(100, 10000, 1000000), labels=c("100", "10000", "1000000"), cex.axis=0.8)
  legend("bottomleft", legend=c("Insertions", "Deletions"), col=c(rgb(1,0,0,0.8), rgb(0,0,1,0.8)), 
       pch=16, cex=1.2, pt.cex=2, bty="o", bg="white")
  
  ## Add Shared, chromosomal distance
  plot(x=as.numeric(hits[which(hits[,"Type"]=="D"),1])[-1], 
     y=as.numeric(hits[which(hits[,"Type"]=="D"),2])[-1], ylab="Intermutation Distances", 
     xlab="Chromosomal Distance", pch=16, cex=0.6, log="y", 
     main=paste0("Chromosomal-Distance Rainfall: ", sample, " - Chromosome ", chromosome), 
     cex.main=1.8, cex.lab=1.3,yaxt="n", col=rgb(0,0,1,0.5), 
     xlim=c(0,full.chr.length), ylim=c(10,1000000))
  
  points(x=as.numeric(hits[which(hits[,"Type"]=="I"),1])[-1], 
       y=as.numeric(hits[which(hits[,"Type"]=="I"),2])[-1],
       pch=16, cex=0.6, col=rgb(1,0,0,0.5))
  
  # Legend and X-axis 
  axis(2, at=c(100, 10000, 1000000), labels=c("100", "10000", "1000000"), cex.axis=0.8)
  legend("bottomleft", legend=c("Insertions", "Deletions"), col=c(rgb(1,0,0,0.8), rgb(0,0,1,0.8)), 
         pch=16, cex=1.2, pt.cex=2, bty="o", bg="white")
  
  if(split == "N"){
    dev.off()
    
  }else if(split == "Y"){
    dev.off()
    
    # Add Single Rainfall-Plots
    pdf(paste0("Rainfall_split_", sample, "_chr", chromosome, ".pdf"), width=18,height=12)
    par(mfrow=c(2,1))
    mar.default <- c(5,4,4,2) + 0.1
    par(mar = mar.default + c(0, 5, 0, 0)) 
    
    # Insertions, relative distance
    plot(x=which(hits[,"Type"]=="I")[-1], 
         y=as.numeric(hits[which(hits[,"Type"]=="I"),2])[-1], 
         ylab="Intermutation Distances", xlab="Indel Number",
         pch=16, cex=0.6, log="y", 
         main=paste0("Insertions: Relative-Distance Rainfall: ", sample, " - Chromosome ", chromosome), 
         cex.main=1.8, cex.lab=1.3,yaxt="n", col=rgb(1,0,0,0.5), 
         xlim=c(0,length(hits[,1])), ylim=c(10,1000000))
    
    # Legend and X-axis
    axis(2, at=c(100, 10000, 1000000), labels=c("100", "10000", "1000000"), cex.axis=0.8)
    legend("bottomleft", legend=c("Insertions"), col=rgb(1,0,0,0.8), 
           pch=16, cex=1.2, pt.cex=2, bty="o", bg="white")
    
    # Insertions, genomic distance
    plot(x=as.numeric(hits[which(hits[,"Type"]=="I"),1])[-1], 
         y=as.numeric(hits[which(hits[,"Type"]=="I"),2])[-1], ylab="Intermutation Distances", 
         xlab="Chromosomal Distance", pch=16, cex=0.6, log="y", 
         main=paste0("Insertions: Chromosomal-Distance Rainfall: ", sample, " - Chromosome ", chromosome), 
         cex.main=1.8, cex.lab=1.3,yaxt="n", col=rgb(1,0,0,0.5), 
         xlim=c(0,full.chr.length), ylim=c(10,1000000))
    
    # Legend and X-axis 
    axis(2, at=c(100, 10000, 1000000), labels=c("100", "10000", "1000000"), cex.axis=0.8)
    legend("bottomleft", legend=c("Insertions"), col=rgb(1,0,0,0.8), 
           pch=16, cex=1.2, pt.cex=2, bty="o", bg="white")
    
    # Deletions, relative distance
    plot(x=which(hits[,"Type"]=="D")[-1], 
         y=as.numeric(hits[which(hits[,"Type"]=="D"),2])[-1], 
         ylab="Intermutation Distances", xlab="Indel Number", pch=16, cex=0.6, log="y", 
         main=paste0("Deletions: Relative-Distance Rainfall: ", sample, " - Chromosome ", chromosome), 
         cex.main=1.8, cex.lab=1.3,yaxt="n", col=rgb(0,0,1,0.5), 
         xlim=c(0,length(hits[,1])), ylim=c(10,1000000))
    
    # Legend and X-axis
    axis(2, at=c(100, 10000, 1000000), labels=c("100", "10000", "1000000"), cex.axis=0.8)
    legend("bottomleft", legend=c("Deletions"), col=rgb(0,0,1,0.8), 
           pch=16, cex=1.2, pt.cex=2, bty="o", bg="white")
    
    # Deletions, genomic distance
    plot(x=as.numeric(hits[which(hits[,"Type"]=="D"),1])[-1], 
         y=as.numeric(hits[which(hits[,"Type"]=="D"),2])[-1], 
         ylab="Intermutation Distances", xlab="Chromosomal Distance", pch=16, cex=0.6, log="y", 
         main=paste0("Deletions: Chromosomal-Distance Rainfall: ", sample, " - Chromosome ", chromosome), 
         cex.main=1.8, cex.lab=1.3,yaxt="n", col=rgb(0,0,1,0.5), 
         xlim=c(0,full.chr.length), ylim=c(10,1000000))
    
    # Legend and X-axis 
    axis(2, at=c(100, 10000, 1000000), labels=c("100", "10000", "1000000"), cex.axis=0.8)
    legend("bottomleft", legend=c("Deletions"), col=rgb(0,0,1,0.8), 
           pch=16, cex=1.2, pt.cex=2, bty="o", bg="white")
    
    dev.off()
  }
  
  
  # d. Restore File Path
  setwd(main.path)
  
}


# 11. Rainfall: by Indeltype
rainfall.types <- function(hits, sample, chromosome){ 
  
  
  # Output folder
  setwd(paste0(main.path,"/Toy/Output"))
  
  
  # a. Take all Distances
  indel.distances <- diff(as.integer(hits[,"Chromosome Position"]))
  names(indel.distances) <- hits[c(-1),"Type"]
  indel.lengths <- as.integer(hits[c(-1),"Length"])
  names(indel.lengths) <- hits[c(-1),"Type"]
  
  
  # b. Plot: Deletion-Types
  pdf(paste0("Rainfall_Types_", sample, "_chr", chromosome, ".pdf"), width=18,height=12)
  mar.default <- c(5,4,4,2) + 0.1
  par(mfrow=c(2,1))
  par(mar = mar.default + c(0, 5, 0, 0)) 
  
  plot(x=1:length(indel.distances), y=indel.distances, ylab="Intermutation Distances", 
       xlab="Relative Genomic Position", pch=16, cex=0.6, log="y", 
       main=paste0("Size-Type Deletion Rainfall: ",  sample, " - Chromosome ", chromosome), 
       cex.main=1.8, cex.lab=1.3,
       yaxt="n", col=rgb(0,0,1,0.0), xlim=c(0,length(indel.distances)), ylim=c(10,1000000))
  
  # Deletionsize: >= -1 bp, < -3bp
  check <- -3 < indel.lengths & indel.lengths <= -1
  points(x=which(check==T), y=indel.distances[which(check==T)], pch=16,
         col="red", cex=0.7)
  
  # Deletionsize: > -5 bp, <= -3bp
  check <- -5 <= indel.lengths & indel.lengths <= -3
  points(x=which(check==T), y=indel.distances[which(check==T)], pch=16,
         col="orange", cex=0.7)
  
  # Deletionsize: > -10 bp, <= -5bp
  check <- -10 <= indel.lengths & indel.lengths < -5
  points(x=which(check==T), y=indel.distances[which(check==T)], pch=16,
         col="yellow", cex=0.7)
  
  # Deletionsize: > -20 bp, <= -10bp
  check <- -20 <= indel.lengths & indel.lengths < -10
  points(x=which(check==T), y=indel.distances[which(check==T)], pch=16,
         col="green", cex=0.7)
  
  # Deletionsize: < -20 bp
  check <- indel.lengths < -20
  points(x=which(check==T), y=indel.distances[which(check==T)], pch=16,
         col="blue", cex=0.7)  
  axis(2, at=c(100, 10000, 1000000), labels=c("100", "10000", "1000000"), cex.axis=0.8)
  
  # Legend
  legend("bottomleft", legend=c("1 to 2 bp", "3 to 5 bp","6 to 10 bp","11 to 20 bp","> 20 bp"), 
         col=c("red", "orange", "yellow", "green", "blue"), pch=16, cex=1.2, pt.cex=2, bty="o", bg="white")
  
  
  # c. Plot: Insertion-Types
  
  plot(1:length(indel.distances), y=indel.distances, ylab="Intermutation Distances", 
       xlab="Relative Genomic Position", pch=16, cex=0.6, log="y", 
       main=paste0("Size-Type Insertion Rainfall: ", sample, " - Chromosome ", chromosome), 
       cex.main=1.8, cex.lab=1.3,
       yaxt="n", col=rgb(0,0,1,0.0), xlim=c(0,length(indel.distances)), ylim=c(10,1000000))
  
  # Insertionsize: = 1 bp, < 3bp
  
  check <- 1 <= indel.lengths & indel.lengths < 3
  points(x=which(check==T), y=indel.distances[which(check==T)], pch=16,
         col="red", cex=0.7)
  
  # Insertionsize: < 5 bp, >= 3bp
  
  check <- 3 <= indel.lengths & indel.lengths <= 5
  points(x=which(check==T), y=indel.distances[which(check==T)], pch=16,
         col="orange", cex=0.7)
  
  # Insertionsize: < 10 bp, >= 5bp
  
  check <- 5 < indel.lengths & indel.lengths <= 10
  points(x=which(check==T), y=indel.distances[which(check==T)], pch=16,
         col="yellow", cex=0.7)
  
  # Insertionsize: < 20 bp, >= 10bp
  
  check <- 10 < indel.lengths & indel.lengths <= 20
  points(x=which(check==T), y=indel.distances[which(check==T)], pch=16,
         col="green", cex=0.7)
  
  # Insertionsize: > 20 bp
  
  check <- 20 < indel.lengths 
  points(x=which(check==T), y=indel.distances[which(check==T)], pch=16,
         col="blue", cex=0.7)  
  axis(2, at=c(100, 10000, 1000000), labels=c("100", "10000", "1000000"), cex.axis=0.8)
  
  # Legend
  legend("bottomleft", legend=c("1 to 2 bp", "3 to 5 bp","6 to 10 bp","11 to 20 bp","> 20 bp"), 
         col=c("red", "orange", "yellow", "green", "blue"), pch=16, cex=1.2, pt.cex=2, bty="o", bg="white")
  dev.off()
  
  
  # d. Restore File Path
  setwd(main.path)
  
}


# 12. Bi-allelic Frequency Plot
indel.baf <- function(hits, contigs.coordinates,
                      sample, chromosome, pindel = "N"){
  
  
  # Output folder
  setwd(paste0(main.path,"/Toy/Output"))
  
  
  # a. Preparation
  full.chr.length <- contigs.coordinates[last(grep(paste0("Chr",chromosome), contigs.coordinates[,1])),
                                    "Chromosome_Position"]
  
  
  # b. Calculate the Ratios (Platypus)
  
  hits.pindel <- hits.platypus <- hits
  
  # First Remove ambiguous calls from BAF analysis
  platypus.call <- str_split_fixed(hits.platypus[,"Platypus_Calls"], "/", n=2)[,1]
  hits.platypus <- hits.platypus[-grep(",", platypus.call),]
  platypus.call <- str_split_fixed(hits.platypus[,"Platypus_Calls"], "/", n=2)[,1]
  platypus.cov <- str_split_fixed(hits.platypus[,"Platypus_Calls"], "/", n=2)[,2] 
  platypus.call <- as.numeric(platypus.call)
  platypus.cov <- as.numeric(platypus.cov)
  platypus.baf <- platypus.call/platypus.cov
  hits.platypus[,"Platypus_Calls"] <- platypus.baf
  
  
  # c. Calculate the Ratios (Pindel)
  
  pindel.call <- str_split_fixed(hits.pindel[,"Pindel_Calls"], "/", n=2)[,1]
  pindel.cov <- str_split_fixed(hits.pindel[,"Pindel_Calls"], "/", n=2)[,2] 
  pindel.call <- as.numeric(pindel.call)
  pindel.cov <- as.numeric(pindel.cov)
  pindel.baf <- pindel.call/pindel.cov
  hits.pindel[,"Pindel_Calls"] <- pindel.baf
  
  
  # d. Plot BAFs
  if(pindel == "N"){
    
    pdf(paste0("BAF_", sample, "_chr", chromosome, ".pdf"), width=18,height=7)
    mar.default <- c(5,4,4,2) + 0.1
    par(mar = mar.default + c(0, 5, 0, 0)) 
    
    ## Platypus
    
    # Empty Chromosome String
    empty <- c(1, full.chr.length)
    plot(rep(1,2)~empty, type="l", lwd=7, ylab="Call:Coverage", xlab="Position [Bp]", cex.lab=1.3,
         main=paste0("Platypus-called BAF - ", sample, " - on Chr. ", chromosome), cex.main=1.8, ylim=c(0,1),
         col="white", yaxt="n")
    axis(2, at=c(0, 1/3, 0.5, 2/3, 1), labels=c("0", "1/3", "1/2", "2/3", "1"), cex.axis=0.8)
    
    # All Insertions
    ins.pos <- as.integer(hits.platypus[which(hits.platypus[,"Type"]=="I"),"Chromosome Position"])
    ins.height <- as.numeric(hits.platypus[which(hits.platypus[,"Type"]=="I"),"Platypus_Calls"]) 
    points(x=ins.pos, y=ins.height, cex=0.8, pch=16, col=rgb(1,0,0,0.5))
    
    # All Deletions  
    dels.pos <- as.integer(hits.platypus[which(hits.platypus[,"Type"]=="D"),"Chromosome Position"])
    dels.height <- as.numeric(hits.platypus[which(hits.platypus[,"Type"]=="D"),"Platypus_Calls"]) 
    alpha.adjust <- length(ins.height)/length(dels.height)
    points(x=dels.pos, y=dels.height, cex=0.8, pch=16, col=rgb(0,0,1,0.5))
    
    # Legend
    legend("bottomleft", legend=c("Insertions", "Deletions"), col=c(rgb(1,0,0,0.8), rgb(0,0,1,0.8)), 
           pch=16, cex=1.5, pt.cex=2, bty="o", bg="white") 
    
    dev.off()
    
  }else if(pindel == "Y"){
    
    pdf(paste0("BAF_", sample, "_chr", chromosome, ".pdf"), width=18,height=14)
    par(mfrow=c(2,1))
    mar.default <- c(5,4,4,2) + 0.1
    par(mar = mar.default + c(0, 5, 0, 0)) 
    
    ## Platypus
    
    # Empty Chromosome String
    empty <- c(1, full.chr.length)
    plot(rep(1,2)~empty, type="l", lwd=7, ylab="Call:Coverage", xlab="Position [Bp]", cex.lab=1.3,
         main=paste0("Platypus-called BAF - ", sample, " - on Chr. ", chromosome), cex.main=1.8, ylim=c(0,1),
         col="white", yaxt="n")
    axis(2, at=c(0, 1/3, 0.5, 2/3, 1), labels=c("0", "1/3", "1/2", "2/3", "1"), cex.axis=0.8)
    
    # All Insertions
    ins.pos <- as.integer(hits.platypus[which(hits.platypus[,"Type"]=="I"),"Chromosome Position"])
    ins.height <- as.numeric(hits.platypus[which(hits.platypus[,"Type"]=="I"),"Platypus_Calls"]) 
    points(x=ins.pos, y=ins.height, cex=0.8, pch=16, col=rgb(1,0,0,0.5))
    
    # All Deletions  
    dels.pos <- as.integer(hits.platypus[which(hits.platypus[,"Type"]=="D"),"Chromosome Position"])
    dels.height <- as.numeric(hits.platypus[which(hits.platypus[,"Type"]=="D"),"Platypus_Calls"]) 
    alpha.adjust <- length(ins.height)/length(dels.height)
    points(x=dels.pos, y=dels.height, cex=0.8, pch=16, col=rgb(0,0,1,0.5))
    
    # Legend
    legend("bottomleft", legend=c("Insertions", "Deletions"), col=c(rgb(1,0,0,0.8), rgb(0,0,1,0.8)), 
           pch=16, cex=1.5, pt.cex=2, bty="o", bg="white") 
    
    
    ## Pindel
    
    # Empty Chromosome String
    empty <- c(1, full.chr.length)
    plot(rep(1,2)~empty, type="l", lwd=7, ylab="Call:Coverage", xlab="Position [Bp]", cex.lab=1.3,
         main=paste0("Pindel-called BAF - ", sample, " - on Chr. ", chromosome), cex.main=1.8, ylim=c(0,1), yaxt="n",
         col="white", yaxt="n")
    axis(2, at=c(0, 1/3, 0.5, 2/3, 1), labels=c("0", "1/3", "1/2", "2/3", "1"), cex.axis=0.8)
    
    # All Insertions
    ins.pos <- as.integer(hits.pindel[which(hits.pindel[,"Type"]=="I"),"Chromosome Position"])
    ins.height <- as.numeric(hits.pindel[which(hits.pindel[,"Type"]=="I"),"Pindel_Calls"]) 
    points(x=ins.pos, y=ins.height, cex=0.8, pch=16, col=rgb(1,0,0,0.5))
    
    # All Deletions  
    dels.pos <- as.integer(hits.pindel[which(hits.pindel[,"Type"]=="D"),"Chromosome Position"])
    dels.height <- as.numeric(hits.pindel[which(hits.pindel[,"Type"]=="D"),"Pindel_Calls"]) 
    alpha.adjust <- length(ins.height)/length(dels.height)
    points(x=dels.pos, y=dels.height, cex=0.8, pch=16, col=rgb(0,0,1,0.5))
    
    # Legend
    legend("bottomleft", legend=c("Insertions", "Deletions"), col=c(rgb(1,0,0,0.8), rgb(0,0,1,0.8)), 
           pch=16, cex=1.5, pt.cex=2, bty="o", bg="white")
    
    dev.off()
  }
  
  
  # e. Restore File Path
  setwd(main.path)
}


# 13. Gene for Indel Presence in Genes
transcript.hits <- function(hits, chromosome){
  
  
  # a. Check only against genes from particular Chromosome 
  devil.transcripts.chr <- devil.transcripts[grep(paste0("Chr", chromosome), devil.transcripts[,1]),]
  
  
  # b. Take the hit-list Position
  names <- colnames(hits)
  hits <- cbind(hits, rep(NA, length(hits[,1])), rep(NA, length(hits[,1])), rep(NA, length(hits[,1])), 
                rep(NA, length(hits[,1])), rep(NA, length(hits[,1])), rep(NA, length(hits[,1])))
  colnames(hits) <- c(names, "Coding Type", "RNA Type", "Exon/Intron Number", 
                      "Strand", "Gene ID", "Transcript ID")
  
  
  # c. Define the Ranges
  Hit.Ranges <- GRanges(seqnames = Rle(rep(paste0("chr", chromosome), length(hits[,1]))),
                        ranges = IRanges(start = as.integer(hits[,"Chromosome Position"]), 
                                         end = as.integer(hits[,"Chromosome Position"])))
  Target.Ranges <- GRanges(seqnames = Rle(rep(paste0("chr", chromosome), length(devil.transcripts.chr[,1]))),
                           ranges = IRanges(start = devil.transcripts.chr[,"Start"], 
                                            end = devil.transcripts.chr[,"End"]))
  
  
  # d. Find Overlaps
  Overlaps <- findOverlaps(Hit.Ranges, Target.Ranges)
  Overlaps <- as.matrix(Overlaps)
  colnames(Overlaps) <- c("Hitsample", "Targetlocation")
  devil.transcripts.chr <- as.matrix(devil.transcripts.chr)

  
  # e. Extract the corresponding Values 
  if(length(as.integer(Overlaps[,1]))>0){
    for (i in 1:length(as.integer(Overlaps[,1]))){
      hits[as.integer(Overlaps[i,"Hitsample"]),"Coding Type"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Coding.Type")]
      hits[as.integer(Overlaps[i,"Hitsample"]),"RNA Type"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("RNA.Type")]
      hits[as.integer(Overlaps[i,"Hitsample"]),"Exon/Intron Number"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Exon.Intron.Number")]
      hits[as.integer(Overlaps[i,"Hitsample"]),"Strand"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Strand")]
      hits[as.integer(Overlaps[i,"Hitsample"]),"Gene ID"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Gene.ID")]
      hits[as.integer(Overlaps[i,"Hitsample"]),"Transcript ID"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Transcript.ID")]    
    }
  }

  
  # f. Shortening the table to variants that were really called in Transcripts
  hits <- hits[-which(is.na(hits[,"Coding Type"])==T),,drop=F]
  
  
  # g. Output
  return(hits)
  
}


# 14. Translate particular Gene Hits with the Reference
hits.translate <- function(hits){
  
  
  # a. Add Colums to hits
  names <- colnames(hits)
  hits.translated <- cbind(hits, rep(NA, length(hits[,1])), rep(NA, length(hits[,1])))
  colnames(hits.translated) <- c(names, "Gene Translation", "Transcript Translation")
  rownames(devil.genes.trans) <- devil.genes.trans[,1]
  rownames(devil.transcripts.trans) <- devil.transcripts.trans[,1]
  
  
  # b. Add External Gene and Transcriptname, also in case of NA
  gene.matches <- match(hits[,"Gene ID"],devil.genes.trans[,"ID"])
  hits.translated[,"Gene Translation"] <- as.character(devil.genes.trans[gene.matches,"External.Name"])
  transcript.matches <- match(hits[,"Transcript ID"],devil.transcripts.trans[,"ID"])
  hits.translated[,"Transcript Translation"] <- as.character(devil.transcripts.trans[transcript.matches,"External.Name"])
  
  
  # c. Merge first and second column, to represent IGV/JBROWSE address
  hits.translated[,1] <- paste(hits.translated[,1], hits.translated[,2], sep=":")
  
  
  # d. Output
  return(hits.translated)
  
}


# 15. Final Results Statistics
final.results.plots <- function(pindel.vcf,
                                platypus.vcf,
                                summary.list,
                                sample,
                                chromosome){
  
  
  # Output folder
  setwd(paste0(main.path,"/Toy/Output"))
  
  
  # a. Retrieve Raw Numbers
  results <- matrix(0, nrow = 1, ncol = 7)
  colnames(results) <- c("Union", "Pindel", "Platypus", "Overlap",
                         "Pindel-Q", "Platypus-Q", "Transcripts")
  
  pindel.vcf[,1] <- paste(pindel.vcf[,1], pindel.vcf[,2], sep="_")
  platypus.vcf[,1] <- paste(platypus.vcf[,1], platypus.vcf[,2], sep="_")  
  results[1] <- length(union(pindel.vcf[,1], platypus.vcf[,1]))
  results[2] <- length(pindel.vcf[,1]) 
  results[3] <- length(platypus.vcf[,1])
  results[4] <- length(summary.list$'overlap'[[1]][,1])
  results[5] <- length(summary.list$'pindel-Q'[[1]][,1])
  results[6] <- length(summary.list$'platypus-Q'[[1]][,1])
  results[7] <- length(summary.list$'in-transcript'[,1])
  
  
  # b. Draw the Bar-Plots
  pdf(paste0("Stats_Plots_", sample, "_", chromosome, ".pdf"), width=15, height=12)
  stats.plot <- barplot(results, beside = T,
                        col = c(rgb(1,0,0,1), rgb(1,0,0,0.85), rgb(1,0,0,0.7),
                                rgb(1,0,0,0.55), rgb(1,0,0,0.4), rgb(1,0,0,0.25),
                                rgb(1,0,0,0.1)),
                        main=paste0("Filtering Results: ", sample, " - Chromosome ", chromosome),
                        xlab=NA, ylab="Number of Calls", cex.lab=1.5, cex.main=2,
                        cex.names=1.5)
  
  # Add Values
  text(x=stats.plot, y=results+3000, labels=as.character(results), xpd=TRUE, cex=1.2, srt=90)
  dev.off()
  
}

# 16. Helper-Files
main.path <- getwd()
setwd(paste0(main.path,"/Helpers"))
contig.lengths <- read.table("Devil_Reference_7.1_Contigs.txt", header=T)
devil.genes <- read.table("Devil_Genes_7.1_chr5.txt", header=T)
devil.transcripts <- read.table("Devil_Transcript_7.1.txt", header=T, sep=";")
devil.genes.trans <- read.table("Devil_Genes_7.1_ENSEMBL_translated.txt", header=T, sep=";")
devil.transcripts.trans <- read.table("Devil_Transcript_7.1_ENSEMBL_translated.txt", header=T, sep=";")
<<<<<<< HEAD
<<<<<<< HEAD
#all.contigs <- read.table("Contig_Translate.txt")
#cosmic.genes <- read.csv("COSMIC_Drivers.csv")
=======
>>>>>>> master
=======
>>>>>>> master
setwd(main.path)