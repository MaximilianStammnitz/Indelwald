############ INDELWALD - HYBRID INDEL CALLING ##############
####################### Version 1.0 ########################

## CROSS-REFERENCING FUNCTIONS
## Last Update - 15/08/2015 ##
## mrs72 / Maximilian Stammnitz ##


# 1. Prepare an integrated Pindel & Platypus - Cross Reference
cross.ref.prepare <- function(filtered_high,
                               platypus,
                               pindel){
  
  
  # a. Merge Contig & Contig Position
  cat("\n Merging Contigs and Contig Positions.")
  
  platypus_filtered <- platypus
  pindel_filtered <- pindel
  filtered_high.firsts <- lapply(filtered_high, function(x) paste(x[,1],x[,2], sep=":"))
  platypus_filtered.firsts <- lapply(platypus_filtered, function(x) paste(x[,1],x[,2], sep=":"))
  pindel_filtered.firsts <- lapply(pindel_filtered, function(x) paste(x[,1],x[,2], sep=":"))
  
  # can be done with lapply.
  for (i in 1:length(filtered_high)){
    filtered_high[[i]][,1] <- filtered_high.firsts[[i]]
    platypus_filtered[[i]][,1] <- platypus_filtered.firsts[[i]]
    pindel_filtered[[i]][,1] <- pindel_filtered.firsts[[i]]
  }
  
  
  # b. Remove Raw Contig Positions & Pindel Calls
  cat("\n Remove Raw Contig Positions & Pindel Calls.")
  
  filtered_high <- lapply(filtered_high, function(x) x[,-c(2,7)])
  platypus_filtered <- lapply(platypus_filtered, function(x) x[,-2])
  pindel_filtered <- lapply(pindel_filtered, function(x) x[,-2])
  
  
  # c. Normalise Double Calls
  cat("\n Normalise Double Calls of stringent Hits.")
  
  ambi_calls.high <- lapply(filtered_high, function(x) grep("[,]", x[,"Platypus_Calls"]))
  ambi_calls.high.first <- mapply(function(x,y) str_split_fixed(x[y,"Platypus_Calls"],"/", 2)[,1], 
                                  x = filtered_high, y = ambi_calls.high)  
  ambi_calls.high.first <- lapply(ambi_calls.high.first, function(x) str_split_fixed(x,",", 2)[,1])
  ambi_calls.high.first <- lapply(ambi_calls.high.first, function(x) as.integer(x))
  ambi_calls.high.second <- mapply(function(x,y) str_split_fixed(x[y,"Platypus_Calls"],"/", 2)[,2], 
                                   x = filtered_high, y = ambi_calls.high)  
  ambi_calls.high.second <- lapply(ambi_calls.high.second, function(x) str_split_fixed(x,",", 2)[,1])
  ambi_calls.high.second <- lapply(ambi_calls.high.second, function(x) as.integer(x))  
  filtered_high.ambis <- mapply(function(x,y) paste(x,y,sep="/"),
                                x = ambi_calls.high.first, y = ambi_calls.high.second)
  
  for (i in 1:length(filtered_high)){
    filtered_high[[i]][ambi_calls.high[[i]],"Platypus_Calls"] <- filtered_high.ambis[[i]]
  }
  
  # d. Output
  cat("\n Done.")
  
  
  return(list("Stringent" = filtered_high,
              "Platypus" = platypus_filtered, 
              "Pindel" = pindel_filtered))
}


# 2. Build an integrated Pindel & Platypus - Cross Reference
build.cross.ref <- function(input){
  
  
  # a. Prepair full List
  cat("\n Pre-allocate list.")
  cross_reference <- vector(mode="list", length=length(names(input[[1]])))
  setnames <- names(input[[1]])
  names(cross_reference) <- setnames
  
  
  # b. Run functions over each set
  cat("\n Set-up Individual Tables. \n")
  for (j in 1:length(setnames)){
    
    cat("\n Preparing ", setnames[j], " ...")
    
    # Pre-allocate and name individual Table
    cross_reference[[j]] <- matrix(0, ncol = length(input[[1]])+4, nrow = length(input$"Stringent"[[j]][,1]))
    colnames(cross_reference[[j]]) <- c("Contig", "Chromosome Position",  "Type",  "Length", setnames)
    cross_reference[[j]][,1:4] <- input$"Stringent"[[j]][,1:4]
    
    # Retrieve calls from particular set of interest, which was filtered stringently
    set.pos <- which(colnames(cross_reference[[j]])==setnames[j])
    cross_reference[[j]][,set.pos] <- input$"Stringent"[[j]][,"Platypus_Calls"] ## Stringent
    cross_reference[[j]][,1] <- paste(cross_reference[[j]][,1], cross_reference[[j]][,3],
                                      cross_reference[[j]][,4], sep=";")
    
    # Extend table with all other sets, which were filtered non-stringently
    set.pos <- which(setnames==setnames[j])
    for (i in c(1:length(setnames))[-set.pos]){
      
      # a. Make a positional overlap between the set at hand and both other non-stringent sets   
      tmp.prep.set.pindel <- paste(prepared.sets$"Pindel"[[i]][,1], prepared.sets$"Pindel"[[i]][,3],
                                   prepared.sets$"Pindel"[[i]][,4], sep=";")
      int.sets.pindel <- intersect(cross_reference[[j]][,1], tmp.prep.set.pindel)
      
      tmp.prep.set.platypus <- paste(prepared.sets$"Platypus"[[i]][,1], prepared.sets$"Platypus"[[i]][,3],
                                     prepared.sets$"Platypus"[[i]][,4], sep=";")
      int.sets.platypus <- intersect(cross_reference[[j]][,1], tmp.prep.set.platypus)
      
      
      # b. Which hits match the intersection
      match.new.pindel <- match(int.sets.pindel, tmp.prep.set.pindel)
      match.old.pindel <- match(int.sets.pindel, cross_reference[[j]][,1])
      match.new.platypus <- match(int.sets.platypus, tmp.prep.set.platypus)
      match.old.platypus <- match(int.sets.platypus, cross_reference[[j]][,1])
      
      
      # c. Update Table: enter Pindel calls and potentiall overwrite them with Platypus Calls
      cross_reference[[j]][match.old.pindel,i+4] <- prepared.sets$"Pindel"[[i]][match.new.pindel,"Pindel Calls"]
      cross_reference[[j]][match.old.platypus,i+4] <- prepared.sets$"Platypus"[[i]][match.new.platypus,"Platypus Calls"]
    }
    
  }
  
  
  # d. Cleaning first lines
  cat("\n Cleaning first lines.")
  cross_reference <- lapply(cross_reference, function(x) {x[,1] <- str_split_fixed(x[,1],";", 3)[,1]
                                                          return(x)})
  
  
  # e. Output
  return(cross_reference)
}


# 3. Obtain integreated Transcript Hits
transcript.hits.cross <- function(hits,
                                  devil.transcripts, 
                                  devil.genes.translation,
                                  devil.transcripts.translation){
  
  
  # a. Define the Transcript Ranges for each chromosome   
  cat("\n Define the Transcript Ranges for each chromosome.")
  
  chromosome <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chrx", "ChrU")
  chr.in.transcripts <- vector(mode="list", length=length(chromosome))
  names(chr.in.transcripts) <- chromosome
  chr.in.transcripts <- lapply(chromosome, function(x) {grep(x, devil.transcripts[,1])})
  
  
  for (i in 1:length(hits)){
    
    # b. Take the hit-list Positions
    cat(paste0("\n Take the hit-list Position: ", names(hits)[i], "."))
    
    names <- colnames(hits[[i]])
    hits[[i]] <- cbind(hits[[i]], rep(NA, length(hits[[i]][,1])), rep(NA, length(hits[[i]][,1])), rep(NA, length(hits[[i]][,1])),
                       rep(NA, length(hits[[i]][,1])), rep(NA, length(hits[[i]][,1])), rep(NA, length(hits[[i]][,1])))
    colnames(hits[[i]]) <- c(names, "Coding Type", "RNA Type", "Exon/Intron Number", "Strand", "Gene ID", "Transcript ID")
    
    
    # c. Define the Hit Ranges for each chromosome  
    cat(paste0("\n Define the Hit Ranges for each chromosome: ", names(hits)[i], "."))
    
    chr.in.hits <- vector(mode="list", length=length(chromosome))
    names(chr.in.hits) <- chromosome
    chr.in.hits <- lapply(chromosome, function(x) {grep(x, hits[[i]][,1])})
    
    
    # d. Build GRanges Objects
    cat(paste0("\n Build GRanges Objects: ", names(hits)[i], "."))
    
    Hit.Ranges <- GRanges(seqnames = Rle(c(rep(chromosome[1], length(chr.in.hits[[1]])),
                                           rep(chromosome[2], length(chr.in.hits[[2]])),
                                           rep(chromosome[3], length(chr.in.hits[[3]])),
                                           rep(chromosome[4], length(chr.in.hits[[4]])),
                                           rep(chromosome[5], length(chr.in.hits[[5]])),
                                           rep(chromosome[6], length(chr.in.hits[[6]])),
                                           rep(chromosome[7], length(chr.in.hits[[7]])),
                                           rep(chromosome[8], length(chr.in.hits[[8]])))),
                          ranges = IRanges(start = as.integer(hits[[i]][,"Chromosome Position"]),
                                           end = as.integer(hits[[i]][,"Chromosome Position"])))
    
    Target.Ranges <- GRanges(seqnames = Rle(c(rep(chromosome[1], length(chr.in.transcripts[[1]])),
                                              rep(chromosome[2], length(chr.in.transcripts[[2]])),
                                              rep(chromosome[3], length(chr.in.transcripts[[3]])),
                                              rep(chromosome[4], length(chr.in.transcripts[[4]])),
                                              rep(chromosome[5], length(chr.in.transcripts[[5]])),
                                              rep(chromosome[6], length(chr.in.transcripts[[6]])),
                                              rep(chromosome[7], length(chr.in.transcripts[[7]])),
                                              rep(chromosome[8], length(chr.in.transcripts[[8]])))),
                             ranges = IRanges(start = devil.transcripts[,"Start"],
                                              end = devil.transcripts[,"End"]))
    
    
    # e. Find Overlaps     
    cat(paste0("\n Find Overlaps: ", names(hits)[i], "."))
    
    Overlaps <- findOverlaps(Hit.Ranges, Target.Ranges)
    Overlaps <- as.matrix(Overlaps)
    colnames(Overlaps) <- c("Hitsample", "Targetlocation")
    
    
    # f. Extract the corresponding Values
    cat(paste0("\n Extract the corresponding Values: ", names(hits)[i], "."))
    
    if(length(as.integer(Overlaps[,1]))>0){
      want.extract <- c("Coding Type", "RNA Type", "Exon/Intron Number", "Strand", "Gene ID", "Transcript ID")
      will.extract <- c("Coding.Type", "RNA.Type", "Exon.Intron.Number", "Strand", "Gene.ID", "Transcript.ID")
      devil.transcripts.tmp <- as.matrix(devil.transcripts)
      hits[[i]][Overlaps[,1],want.extract] <- devil.transcripts.tmp[Overlaps[,2],will.extract]
    }
    
    
    # g. Shortening the table to variants that were really called
    cat(paste0("\n Shortening the table: ", names(hits)[i], "."))
    
    hits[[i]] <- hits[[i]][-which(is.na(hits[[i]][,"Coding Type"])==T),,drop=F]
    
    
    # h. Translating Gene Names to HGNC
    cat(paste0("\n Translating Gene Names to HGNC: ", names(hits)[i], "."))
    
    names <- colnames(hits[[i]])
    hits[[i]] <- cbind(hits[[i]], rep(NA, length(hits[[i]][,1])), rep(NA, length(hits[[i]][,1])))
    colnames(hits[[i]]) <- c(names, "Gene Translation", "Transcript Translation")
    rownames(devil.genes.trans) <- devil.genes.trans[,1]
    rownames(devil.transcripts.trans) <- devil.transcripts.trans[,1]
    
    gene.matches <- match(hits[[i]][,"Gene ID"],devil.genes.trans[,"ID"])
    hits[[i]][,"Gene Translation"] <- as.character(devil.genes.trans[gene.matches,"External.Name"])
    transcript.matches <- match(hits[[i]][,"Transcript ID"],devil.transcripts.trans[,"ID"])
    hits[[i]][,"Transcript Translation"] <- as.character(devil.transcripts.trans[transcript.matches,"External.Name"])
    
    
    cat("\n Done. \n")
    
  }
  
  # h. Output
  return(hits)
}


# 4. Build a Sophisticated Venn for DFT2, based on Cross-References
cross.DFTD.venn <- function(hits = cross.transcripts, 
                            threshold = 4, 
                            type,
                            contig.lengths = contig.size,
                            titulo){
  
  
  # Output folder
  setwd(paste0(main.path,"/Toy/Output"))
  
  
  # a. Extract the Exons ##
  hits <- lapply(hits, function(x) x[which(x[,"RNA Type"]=="exon"),])
  hits <- lapply(hits, function(x) {x[,1] <- paste(x[,1], x[,2], x[,3], x[,4], x[,16], sep=";")
                                    return(x)})
  
  # b. Build-up a binary copy of the Mastertable
  all.positions <- sapply(hits, function(x) as.character(x[,1]))
  all.positions <- Reduce(union, list(all.positions[[1]],
                                      all.positions[[2]],
                                      all.positions[[3]],
                                      all.positions[[4]],
                                      all.positions[[5]],
                                      all.positions[[6]]))
  all.positions <- sort(all.positions) 
  bin.venn <- matrix(0, ncol = 3+length(names(hits)), nrow = length(all.positions))
  colnames(bin.venn) <- c("Contig", names(hits),
                          "Gene ID", "Gene Translation")
  bin.venn[,1] <- all.positions
  
  
  # c. Extract strict calls
  sets <- names(hits)  
  for (i in 1:length(sets)){
    hit.set <- which(names(hits)==sets[i])
    assign(paste0("n.", sets[i]), match(hits[[hit.set]][,"Contig"],bin.venn[,"Contig"]))
    tmp <- paste0("n.", sets[i])
    bin.venn[get(tmp),sets[i]] <- hits[[hit.set]][,sets[i]]
  }
  bin.stringent.only <- bin.venn
  
  
  # d. Add non-stringent calls, in case the corresponding bin value is still 0
  for (i in 1:length(sets)){
    hit.set <- which(names(hits)==sets[i])
    tmp <- get(paste0("n.", sets[i]))
    bin.venn[tmp,sets[-i]] <- hits[[hit.set]][,sets[-i]]
  }
  
  
  # e. Remove fully empty lines
  #for(i in 1:length(sets)){
  #  assign(paste0("empty.",sets[i]),which(bin.venn[,sets[i]]==0))
  #}
  #empty.full <- Reduce(intersect, list(empty.Host1, 
  #                                     empty.Host2, 
  #                                     empty.Host3, 
  #                                     empty.Tumour1, 
  #                                     empty.Tumour2, 
  #                                     empty.Tumour3))
  
  
  # g. Retrieve Gene Name/ID
  for (i in 1:length(hits)){
    match.pos <- match(hits[[i]][,1],all.positions)
    bin.venn[match.pos,c(length(hits)+2:3)] <- bin.stringent.only[match.pos,c(length(hits)+2:3)] <- hits[[i]][,c("Gene ID", "Gene Translation")]
  }
  summary.venn <- bin.venn
  #bin.venn <- bin.venn[-empty.full,]
  #bin.stringent.only <- bin.stringent.only[-empty.full,]
  
  
  # g. Make a binary Table, based on Threshold and Duplicate-Hit Specification
  if(type=="all"){
    
    names <- colnames(summary.venn)
    bin.venn <- cbind(summary.venn, rep(0, length(summary.venn[,1])))
    colnames(bin.venn) <- c(names, "Germline")
    germ <- vector(mode="list", length=3)
    names(germ) <- c("Host1", "Host2", "Host3")
    for (i in 1:length(names(germ))){
      germ[[i]] <- which(as.integer(str_split_fixed(bin.venn[,c(1+i)],"/",2)[,1])>=threshold)
      bin.venn[germ[[i]],"Germline"] <- 1
    }
    
    bin.venn <- bin.venn[,-c(2:4)]
    bin.venn <- bin.venn[,c(1,7,2:6)]
    
    pdf(paste0("DFTD_Exon_Overlap_all.pdf"), width = 10, height = 8)
    
    
  }else if (type=="cleared"){
    
    names <- colnames(summary.venn)
    bin.venn <- cbind(summary.venn, rep(0, length(summary.venn[,1])))
    colnames(bin.venn) <- c(names, "Germline")
    germ <- vector(mode="list", length=3)
    names(germ) <- c("Host1", "Host2", "Host3")
    for (i in 1:length(names(germ))){
      germ[[i]] <- which(as.integer(str_split_fixed(bin.venn[,c(1+i)],"/",2)[,1])>=threshold)
      bin.venn[germ[[i]],"Germline"] <- 1
    }
    
    ## Take duplicated reads, remove them (also from the Summary)
    multi.pos <- which(duplicated(str_split_fixed(bin.venn[,1], ";", 4)[,1])==T)
    bin.venn <- bin.venn[-sort(c(multi.pos, multi.pos-1)),]           # here
    bin.stringent.only <- bin.stringent.only[-sort(c(multi.pos, multi.pos-1)),]   # here
    summary.venn <- summary.venn[-sort(c(multi.pos, multi.pos-1)),]   # here
    bin.venn <- bin.venn[,-c(2:4)]
    bin.venn <- bin.venn[,c(1,7,2:6)]
    
    pdf(paste0("DFTD_Exon_Overlap_cleared.pdf"), width = 10, height = 8)
    
    
  }else if (type=="unified"){
    
    ## Take duplicated reads, merge them (contigname change, calls from all; overwritten)
    multi.pos <- which(duplicated(str_split_fixed(bin.venn[,1], ";", 4)[,1])==T)
    multi.pos.contig <- str_split_fixed(bin.venn[sort(unique(c(multi.pos, multi.pos-1))),1], ";", 5)[,1:2]
    bin.venn[sort(unique(c(multi.pos, multi.pos-1))),1] <- paste(multi.pos.contig[,1],multi.pos.contig[,2],sep=";")
    multi.pos.names <- unique(bin.venn[multi.pos,1])
    
    # For-loop: how many copies are there per Duplicate? - merge up these values
    multi.pos.id <- vector(mode="list",length=length(multi.pos.names))  
    for (i in 1:length(multi.pos.names)){
      multi.pos.id[[i]] <- grep(multi.pos.names[i], bin.venn[,1])
      bin.venn[multi.pos.id[[i]][1],2] <- paste(sum(as.integer(str_split_fixed(bin.venn[multi.pos.id[[i]],2],"/",2)[,1])),
                                                sum(as.integer(str_split_fixed(bin.venn[multi.pos.id[[i]],2],"/",2)[,2]), na.rm=T),
                                                sep="/")
      bin.venn[multi.pos.id[[i]][1],3] <- paste(sum(as.integer(str_split_fixed(bin.venn[multi.pos.id[[i]],3],"/",2)[,1])),
                                                sum(as.integer(str_split_fixed(bin.venn[multi.pos.id[[i]],3],"/",2)[,2]), na.rm=T),
                                                sep="/")
      bin.venn[multi.pos.id[[i]][1],4] <- paste(sum(as.integer(str_split_fixed(bin.venn[multi.pos.id[[i]],4],"/",2)[,1])),
                                                sum(as.integer(str_split_fixed(bin.venn[multi.pos.id[[i]],4],"/",2)[,2]), na.rm=T),
                                                sep="/") 
      bin.venn[multi.pos.id[[i]][1],5] <- paste(sum(as.integer(str_split_fixed(bin.venn[multi.pos.id[[i]],5],"/",2)[,1])),
                                                sum(as.integer(str_split_fixed(bin.venn[multi.pos.id[[i]],5],"/",2)[,2]), na.rm=T),
                                                sep="/") 
      bin.venn[multi.pos.id[[i]][1],6] <- paste(sum(as.integer(str_split_fixed(bin.venn[multi.pos.id[[i]],6],"/",2)[,1])),
                                                sum(as.integer(str_split_fixed(bin.venn[multi.pos.id[[i]],6],"/",2)[,2]), na.rm=T),
                                                sep="/") 
      bin.venn[multi.pos.id[[i]][1],7] <- paste(sum(as.integer(str_split_fixed(bin.venn[multi.pos.id[[i]],7],"/",2)[,1])),
                                                sum(as.integer(str_split_fixed(bin.venn[multi.pos.id[[i]],7],"/",2)[,2]), na.rm=T),
                                                sep="/")
    }
    summary.venn <- bin.venn
    
    # Mark as germline
    names <- colnames(summary.venn)
    bin.venn <- cbind(summary.venn, rep(0, length(summary.venn[,1])))
    colnames(bin.venn) <- c(names, "Germline")
    germ <- vector(mode="list", length=3)
    names(germ) <- c("Host1", "Host2", "Host3")
    for (i in 1:length(names(germ))){
      germ[[i]] <- which(as.integer(str_split_fixed(bin.venn[,c(1+i)],"/",2)[,1])>=threshold)
      bin.venn[germ[[i]],"Germline"] <- 1
    }
    
    # After merging, remove the duplicaes (also from the Summary)
    bin.venn <- bin.venn[-multi.pos,]
    bin.stringent.only <-bin.stringent.only[-multi.pos,]
    summary.venn <- summary.venn[-multi.pos,]
    bin.venn <- bin.venn[,-c(2:4)]
    bin.venn <- bin.venn[,c(1,7,2:6)]
    
    pdf(paste0("DFTD_Exon_Overlap_unified.pdf"), width = 10, height = 8)
  }  
  
  ## Binarise tumour data  
  bin.venn[which(as.integer(str_split_fixed(bin.venn[,"Tumour1"],"/",2)[,1])<=threshold),"Tumour1"] <- 0
  bin.venn[which(as.integer(str_split_fixed(bin.venn[,"Tumour1"],"/",2)[,1])>threshold),"Tumour1"] <- 1
  bin.venn[which(as.integer(str_split_fixed(bin.venn[,"Tumour2"],"/",2)[,1])<=threshold),"Tumour2"] <- 0
  bin.venn[which(as.integer(str_split_fixed(bin.venn[,"Tumour2"],"/",2)[,1])>threshold),"Tumour2"] <- 1
  bin.venn[which(as.integer(str_split_fixed(bin.venn[,"Tumour3"],"/",2)[,1])<=threshold),"Tumour3"] <- 0
  bin.venn[which(as.integer(str_split_fixed(bin.venn[,"Tumour3"],"/",2)[,1])>threshold),"Tumour3"] <- 1    
  
  
  # h. Calculate the overlap areas
  
  n1 <- hit.Tumour2 <- which(bin.venn[,"Tumour2"]==1) # n1
  n2 <- hit.germ <- which(bin.venn[,"Germline"]==1) # n2
  n3 <- hit.Tumour3 <- which(bin.venn[,"Tumour3"]==1) # n3
  n4 <- hit.Tumour1 <- which(bin.venn[,"Tumour1"]==1) # n4
  
  n12 <- intersect(n1,n2)
  n13 <- intersect(n1,n3)
  n14 <- intersect(n1,n4)
  n23 <- intersect(n2,n3)
  n24 <- intersect(n2,n4)
  n34 <- intersect(n3,n4)
  n123 <- intersect(intersect(n1,n2),n3)
  n124 <- intersect(intersect(n1,n2),n4)
  n134 <- intersect(intersect(n1,n3),n4)
  n234 <- intersect(intersect(n2,n3),n4)
  n1234 <- Reduce(intersect, list(n1,n2,n3,n4))
  
  # Amounts
  hit.Tumour2.n <- length(n1)
  hit.germ.n <- length(n2)
  hit.Tumour3.n <- length(n3)
  hit.Tumour1.n <- length(n4)
  n12.n <- length(n12)
  n13.n <- length(n13)
  n14.n <- length(n14)
  n23.n <- length(n23)
  n24.n <- length(n24)
  n34.n <- length(n34)
  n123.n <- length(n123)
  n124.n <- length(n124)
  n134.n <- length(n134)
  n234.n <- length(n234)
  n1234.n <- length(n1234)
  
  
  # i. Draw Venn
  
  draw.quad.venn(area1 = hit.Tumour2.n,
                 area2 = hit.germ.n,
                 area3 = hit.Tumour3.n,
                 area4 = hit.Tumour1.n,
                 n12 = n12.n,
                 n13 = n13.n,
                 n14 = n14.n,
                 n23 = n23.n,
                 n24 = n24.n,
                 n34 = n34.n,
                 n123 = n123.n,
                 n124 = n124.n,
                 n134 = n134.n,
                 n234 = n234.n,
                 n1234 = n1234.n,
                 category = c("Tumour 2", "Germline", "Tumour 3", "Tumour 1"),
                 lwd = c(4,4,4,2),
                 lty = c(3,3,3,1),
                 col = rep("black", 4),
                 fill = c(rgb(1,0,0,0.5), rgb(0,0,1,0.7), rgb(1,0,0,0.5), rgb(1,0,0,0.5)),
                 alpha = c(0.7,0.7,0.7,0.7),
                 cat.cex = rep(1.8,4),
                 cat.fontfamily = rep("serif", 4),
                 cex = rep(1.5,15),
                 fontfamily = rep("serif", 15),
                 ind=T)
    dev.off()
  
  
  # j. Hierarchical Clustering
  require(cluster)
  
  ## On Binary Exon Indel Tables
  colnames(bin.venn)[3:5] <- c("Tumour 1", "Tumour 2", "Tumour 3")
  cluster.n <- as.matrix(bin.venn[,2:5])
  class(cluster.n) <- "numeric"
  cluster.samples <- agnes(t(cluster.n))
  heat.cols <- colorRampPalette(c("azure", "darkgreen"))(n = 2)
  
  pdf(paste0("DFTD_Hierarchical_Clustering.pdf"), width=10, height=7)
  plot(cluster.samples, main = "DFTD-Clustering by Sample",
       xlab = "Samples", which.plots=2)
  options(warn=-1)
  heatmap(t(cluster.n), Colv=F, Rowv=F, cexCol=0.01, scale='none',
          xlab="Indels", ylab="Samples", col=heat.cols, na.rm=T,
          cexRow = 1)
  options(warn=0)
  dev.off()
  
  
  # k. Final Sets
  n1 <- hit.Tumour2 <- which(bin.venn[,"Tumour 2"]==1) # n1
  n2 <- hit.germ <- which(bin.venn[,"Germline"]==1) # n2
  n3 <- hit.Tumour3 <- which(bin.venn[,"Tumour 3"]==1) # n3
  n4 <- hit.Tumour1 <- which(bin.venn[,"Tumour 1"]==1) # n4
  
  # Tumour 1 
  unique.Tumour1 <- hit.Tumour1[which(!(hit.Tumour1 %in% Reduce(union, list(hit.Tumour2,hit.Tumour3,hit.germ)))==T)]
  som.shared.all <- Reduce(intersect, list(hit.Tumour1,hit.Tumour2,hit.Tumour3))
  som.shared.all <- som.shared.all[-which(bin.venn[som.shared.all,"Germline"]==1)]
  som.shared.Tumour1.Tumour2 <- intersect(hit.Tumour1,hit.Tumour2)
  som.shared.Tumour1.Tumour2 <- som.shared.Tumour1.Tumour2[-which((som.shared.Tumour1.Tumour2 %in% Reduce(union, list(hit.Tumour3, hit.germ)))==T)]
  som.shared.Tumour1.Tumour3 <- intersect(hit.Tumour1,hit.Tumour3)
  som.shared.Tumour1.Tumour3 <- som.shared.Tumour1.Tumour3[-which((som.shared.Tumour1.Tumour3 %in% Reduce(union, list(hit.Tumour2, hit.germ)))==T)]
  
  # Tumour2
  unique.Tumour2 <- hit.Tumour2[which(!(hit.Tumour2 %in% Reduce(union, list(hit.Tumour1,hit.Tumour3,hit.germ)))==T)]
  som.shared.Tumour2.Tumour3 <- intersect(hit.Tumour2,hit.Tumour3)
  som.shared.Tumour2.Tumour3 <- som.shared.Tumour2.Tumour3[-which((som.shared.Tumour2.Tumour3 %in% Reduce(union, list(hit.Tumour1, hit.germ)))==T)]
  
  # Tumour3
  unique.Tumour3 <- hit.Tumour3[which(!(hit.Tumour3 %in% Reduce(union, list(hit.Tumour1,hit.Tumour2,hit.germ)))==T)]
  
  # Germline
  germ.all <- which(bin.venn[,"Germline"]==1)
  
  # Add Set Index to Binary-Table and Summary
  names <- colnames(bin.venn)
  bin.venn <- cbind(bin.venn,rep(NA,length(bin.venn[,1])))
  colnames(bin.venn) <- c(names, "Type")
  bin.venn[unique.Tumour1,"Type"] <- "Unique Tumour1"
  bin.venn[unique.Tumour2,"Type"] <- "Unique Tumour2"
  bin.venn[unique.Tumour3,"Type"] <- "Unique Tumour3"
  bin.venn[germ.all,"Type"] <- "Germline"
  bin.venn[som.shared.all,"Type"] <- "Shared Som."
  bin.venn[som.shared.Tumour1.Tumour2,"Type"] <- "Shared Tumour1-Tumour2"
  bin.venn[som.shared.Tumour1.Tumour3,"Type"] <- "Shared Tumour1-Tumour3"
  bin.venn[som.shared.Tumour2.Tumour3,"Type"] <- "Shared Tumour2-Tumour3"
  bin.stringent.only <- bin.stringent.only[-which(is.na(bin.venn[,"Type"])==T),]
  bin.venn <- bin.venn[-which(is.na(bin.venn[,"Type"])==T),]
  
  # Add Set Index to Summary
  names <- colnames(summary.venn)
  summary.venn <- cbind(summary.venn,rep(NA,length(summary.venn[,1])))
  colnames(summary.venn) <- c(names, "Type")
  summary.venn[unique.Tumour1,"Type"] <- "Unique Tumour1"
  summary.venn[unique.Tumour2,"Type"] <- "Unique Tumour2"
  summary.venn[unique.Tumour3,"Type"] <- "Unique Tumour3"
  summary.venn[germ.all,"Type"] <- "Germline"
  summary.venn[som.shared.all,"Type"] <- "Shared Som."
  summary.venn[som.shared.Tumour1.Tumour2,"Type"] <- "Shared Tumour1-Tumour2"
  summary.venn[som.shared.Tumour1.Tumour3,"Type"] <- "Shared Tumour1-Tumour3"
  summary.venn[som.shared.Tumour2.Tumour3,"Type"] <- "Shared Tumour2-Tumour3"
  
  
  # l. Check the Calls vs. COSMIC database of Cancer Drivers
  genes <- unique(sort(c(summary.venn[which(summary.venn[,"Type"]=="Shared Som."),"Gene Translation"],
                         summary.venn[which(summary.venn[,"Type"]=="Shared Tumour1-Tumour2"),"Gene Translation"],
                         summary.venn[which(summary.venn[,"Type"]=="Shared Tumour1-Tumour3"),"Gene Translation"],
                         summary.venn[which(summary.venn[,"Type"]=="Shared Tumour2-Tumour3"),"Gene Translation"],
                         summary.venn[which(summary.venn[,"Type"]=="unique.Tumour1"),"Gene Translation"],
                         summary.venn[which(summary.venn[,"Type"]=="unique.Tumour2"),"Gene Translation"],
                         summary.venn[which(summary.venn[,"Type"]=="unique.Tumour3"),"Gene Translation"])))
  cosmic.check <- c()
  cosmic.check <- sapply(genes, function(x) which(cosmic.genes[,"Gene.Symbol"]==x))
  cosmic.pos <- sapply(cosmic.check, function(x) length(x))
  cat("\n COSMIC Driver Gene Hits: ", names(which(cosmic.pos!=0)))
  cosmic.lines <- c()
  for (i in 1:length(which(cosmic.pos!=0))){
    cosmic.lines <- c(cosmic.lines, as.integer(unlist(cosmic.check[which(cosmic.pos!=0)][i])))
  }
  cosmic.genes <- cosmic.genes[cosmic.lines,]
  cosmic.names <- names(which(cosmic.pos!=0))
  
  # Set path
  setwd(main.path)
}


# 5. Platypus Raw Output Transcript Summary
platypus.transcript.summary <- function(hits,
                                        chromosome,
                                        contig.lengths,
                                        devil.transcripts){
  
  
  # a. Build a Raw Matrix
  
  hits.summary <- matrix(NA, ncol = 8, nrow = length(hits[,1]))
  colnames(hits.summary) <- c("Contig", "Contig Position", "Chromosome Position",
                              "Type", "Length", "Platypus Calls", "Platypus Quality", 
                              "Platypus Filter")
  
  
  # b. Fill all available cells
  
  hits.summary[,1] <- hits[,1]
  hits.summary[,2] <- hits[,2]
  hits.summary[,7] <- hits[,5]
  hits.summary[,8] <- hits[,6]
  hits.summary[,5] <- nchar(hits[,4])-nchar(hits[,3])
  hits.summary[which(hits.summary[,5]<0) ,4] <- "D" 
  hits.summary[which(hits.summary[,5]>0) ,4] <- "I" 
  
  
  # c. Adding Platypus calls
  
  tmp.1 <- str_split_fixed(hits[,"Tumour"], ":", 6)[,6]
  tmp.2 <- str_split_fixed(hits[,"Tumour"], ":", 6)[,5]
  hits.summary[,6] <- paste(tmp.1, tmp.2, sep="/")
  tmp.1 <- str_split_fixed(hits.summary[grep(",", hits.summary[,6]),6], "/", 2)[,1]
  tmp.2 <- str_split_fixed(hits.summary[grep(",", hits.summary[,6]),6], "/", 2)[,2]
  tmp.1 <- str_split_fixed(tmp.1, ",", 2)[,1]
  tmp.2 <- str_split_fixed(tmp.2, ",", 2)[,1]
  hits.summary[grep(",", hits.summary[,6]),6] <- paste(tmp.1, tmp.2, sep="/")
  
  
  # d. Adding Chromosome Position
  
  unique.hit.contigs <- unique(hits.summary[,1])
  
  if(length(grep("000000000",unique.hit.contigs[1]))!=0){
    starters <- which(hits.summary[,"Contig"]==unique.hit.contigs[1])
    match.contigs.lengths <- match(hits.summary[-starters,"Contig"],contig.lengths[,"Contig"])
    
    # Update Starters
    hits.summary[starters,"Chromosome Position"] <- hits.summary[starters,"Contig Position"]
    
    # Update non-Starters
    hits.summary[-starters,"Chromosome Position"] <- as.numeric(hits.summary[-starters,"Contig Position"])+
      contig.lengths[match.contigs.lengths-1,"Chromosome_Position"]
    
  }else{
    match.contigs.lengths <- match(hits.summary[,"Contig"],contig.lengths[,"Contig"])
    
    # Update all positions
    hits.summary[,"Chromosome Position"] <- as.numeric(hits.summary[,"Contig Position"])+
      contig.lengths[match.contigs.lengths-1,"Chromosome_Position"]
  }
  
  
  # e. Checking for Transcript Hits
  
  devil.transcripts.chr <- devil.transcripts[grep(paste0("Chr", chromosome), devil.transcripts[,1]),]
  names <- colnames(hits.summary)
  hits.summary <- cbind(hits.summary, rep(NA, length(hits.summary[,1])), rep(NA, length(hits.summary[,1])), 
                        rep(NA, length(hits.summary[,1])), rep(NA, length(hits.summary[,1])), 
                        rep(NA, length(hits.summary[,1])), rep(NA, length(hits.summary[,1])))
  colnames(hits.summary) <- c(names, "Coding Type", "RNA Type", "Exon/Intron Number", "Strand", "Gene ID", "Transcript ID")
  Hit.Ranges <- GRanges(seqnames = Rle(rep(paste0("chr", chromosome), length(hits.summary[,1]))),
                        ranges = IRanges(start = as.integer(hits.summary[,"Chromosome Position"]), 
                                         end = as.integer(hits.summary[,"Chromosome Position"])))
  
  Target.Ranges <- GRanges(seqnames = Rle(rep(paste0("chr", chromosome), length(devil.transcripts.chr[,1]))),
                           ranges = IRanges(start = devil.transcripts.chr[,"Start"], 
                                            end = devil.transcripts.chr[,"End"]))
  Overlaps <- findOverlaps(Hit.Ranges, Target.Ranges)
  Overlaps <- as.matrix(Overlaps)
  colnames(Overlaps) <- c("Hitsample", "Targetlocation")
  devil.transcripts.chr <- as.matrix(devil.transcripts.chr)
  #options(warn=-1)
  if(length(as.integer(Overlaps[,1]))>0){
    for (i in 1:length(as.integer(Overlaps[,1]))){
      hits.summary[as.integer(Overlaps[i,"Hitsample"]),"Coding Type"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Coding.Type")]
      hits.summary[as.integer(Overlaps[i,"Hitsample"]),"RNA Type"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("RNA.Type")]
      hits.summary[as.integer(Overlaps[i,"Hitsample"]),"Exon/Intron Number"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Exon.Intron.Number")]
      hits.summary[as.integer(Overlaps[i,"Hitsample"]),"Strand"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Strand")]
      hits.summary[as.integer(Overlaps[i,"Hitsample"]),"Gene ID"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Gene.ID")]
      hits.summary[as.integer(Overlaps[i,"Hitsample"]),"Transcript ID"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Transcript.ID")]    
    }
    #options(warn=0)  
  }
  hits.summary <- hits.summary[-which(is.na(hits.summary[,"Coding Type"])==T),,drop=F]
  
  # e. Output
  
  return(hits.summary)
  
}


# 6. Pindel Raw Output Transcript Summary
pindel.transcript.summary <- function(hits,
                                      chromosome ,
                                      contig.lengths,
                                      devil.transcripts){
  
  
  # a. Build a Raw Matrix
  
  hits.summary <- matrix(NA, ncol = 8, nrow = length(hits[,1]))
  colnames(hits.summary) <- c("Contig", "Contig Position", "Chromosome Position",
                              "Type", "Length", "Pindel Calls", "Pindel Quality", "Pindel Filter")
  
  
  # b. Fill all available cells
  
  hits.summary[,1] <- hits[,1]
  hits.summary[,2] <- hits[,2]
  hits.summary[,7] <- hits[,5]
  hits.summary[,8] <- hits[,6]
  hits.summary[,5] <- nchar(hits[,4])-nchar(hits[,3])
  hits.summary[which(hits.summary[,5]<0) ,4] <- "D" 
  hits.summary[which(hits.summary[,5]>0) ,4] <- "I" 
  
  
  # c. Adding Pindel calls
  
  split <- matrix(NA, ncol = 10, nrow = length(hits[,"Tumour"]))
  split <- str_split_fixed(hits[,"Tumour"], ":", n=10)
  hits.summary[,"Pindel Calls"] <- paste(as.numeric(split[,2])+as.numeric(split[,3]),
                                         as.numeric(split[,8])+as.numeric(split[,9]), sep="/")
  
  
  # d. Adding Chromosome Position
  
  unique.hit.contigs <- unique(hits.summary[,1])
  
  if(length(grep("000000000",unique.hit.contigs[1]))!=0){
    starters <- which(hits.summary[,"Contig"]==unique.hit.contigs[1])
    match.contigs.lengths <- match(hits.summary[-starters,"Contig"],contig.lengths[,"Contig"])
    
    # Update Starters
    hits.summary[starters,"Chromosome Position"] <- hits.summary[starters,"Contig Position"]
    
    # Update non-Starters
    hits.summary[-starters,"Chromosome Position"] <- as.numeric(hits.summary[-starters,"Contig Position"])+
      contig.lengths[match.contigs.lengths-1,"Chromosome_Position"]
    
  }else{
    match.contigs.lengths <- match(hits.summary[,"Contig"],contig.lengths[,"Contig"])
    
    # Update all positions
    hits.summary[,"Chromosome Position"] <- as.numeric(hits.summary[,"Contig Position"])+
      contig.lengths[match.contigs.lengths-1,"Chromosome_Position"]
  }
  
  
  # e. Checking for Transcript Hits
  
  devil.transcripts.chr <- devil.transcripts[grep(paste0("Chr", chromosome), devil.transcripts[,1]),]
  names <- colnames(hits.summary)
  hits.summary <- cbind(hits.summary, rep(NA, length(hits.summary[,1])), rep(NA, length(hits.summary[,1])), 
                        rep(NA, length(hits.summary[,1])), rep(NA, length(hits.summary[,1])), 
                        rep(NA, length(hits.summary[,1])), rep(NA, length(hits.summary[,1])))
  colnames(hits.summary) <- c(names, "Coding Type", "RNA Type", "Exon/Intron Number", "Strand", "Gene ID", "Transcript ID")
  Hit.Ranges <- GRanges(seqnames = Rle(rep(paste0("chr", chromosome), length(hits.summary[,1]))),
                        ranges = IRanges(start = as.integer(hits.summary[,"Chromosome Position"]), 
                                         end = as.integer(hits.summary[,"Chromosome Position"])))
  
  Target.Ranges <- GRanges(seqnames = Rle(rep(paste0("chr", chromosome), length(devil.transcripts.chr[,1]))),
                           ranges = IRanges(start = devil.transcripts.chr[,"Start"], 
                                            end = devil.transcripts.chr[,"End"]))
  Overlaps <- findOverlaps(Hit.Ranges, Target.Ranges)
  Overlaps <- as.matrix(Overlaps)
  colnames(Overlaps) <- c("Hitsample", "Targetlocation")
  devil.transcripts.chr <- as.matrix(devil.transcripts.chr)
  #options(warn=-1)
  if(length(as.integer(Overlaps[,1]))>0){
    for (i in 1:length(as.integer(Overlaps[,1]))){
      hits.summary[as.integer(Overlaps[i,"Hitsample"]),"Coding Type"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Coding.Type")]
      hits.summary[as.integer(Overlaps[i,"Hitsample"]),"RNA Type"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("RNA.Type")]
      hits.summary[as.integer(Overlaps[i,"Hitsample"]),"Exon/Intron Number"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Exon.Intron.Number")]
      hits.summary[as.integer(Overlaps[i,"Hitsample"]),"Strand"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Strand")]
      hits.summary[as.integer(Overlaps[i,"Hitsample"]),"Gene ID"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Gene.ID")]
      hits.summary[as.integer(Overlaps[i,"Hitsample"]),"Transcript ID"] <- devil.transcripts.chr[as.integer(Overlaps[i,"Targetlocation"]),c("Transcript.ID")]    
    }
    #options(warn=0)  
  }
  hits.summary <- hits.summary[-which(is.na(hits.summary[,"Coding Type"])==T),,drop=F]
  
  # e. Output
  
  return(hits.summary)
  
}

# 7. Helper-Files
main.path <- getwd()
setwd(paste0(main.path,"/Helpers"))
contig.lengths <- read.table("Devil_Reference_7.1_Contigs.txt", header=T)
devil.genes <- read.table("Devil_Genes_7.1_chr5.txt", header=T)
devil.transcripts <- read.table("Devil_Transcript_7.1.txt", header=T, sep=";")
devil.genes.trans <- read.table("Devil_Genes_7.1_ENSEMBL_translated.txt", header=T, sep=";")
devil.transcripts.trans <- read.table("Devil_Transcript_7.1_ENSEMBL_translated.txt", header=T, sep=";")
cosmic.genes <- read.csv("COSMIC_Drivers.csv")
setwd(main.path)