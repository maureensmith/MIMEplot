readCountData_1d <- function(countDir, ref.df, nucl.df) {
  L=nrow(ref.df)
  
  countDir1d=paste0(countDir, "/1d")
  coverage1d_all.df=data.frame()
  
  if(dir.exists(countDir1d)) {
    files1d = list.files(countDir1d)
    index=0
    for(f in files1d) {
      #filename without extension
      f_name=tools::file_path_sans_ext(f)
      
      # read in 1d file
      file1d=paste0(countDir1d,"/", f)
      if(file.exists(file1d)) {
        index=index+1
        
        counts1dTable =  read.table(file1d, header = T, sep = "\t")
        
        pos = counts1dTable$pos1
        #for some reason counts prints an additional line 
        validPos = which(pos <= L)
        pos= pos[validPos]
        
        wtNucls = ref.df$asInt[pos]
        
        counts_1d= as.matrix(setNames(counts1dTable[validPos,-1],NULL))
        dimnames(counts_1d)<-NULL
        mutations_1d = do.call(rbind, lapply(1:length(wtNucls), function(x) counts_1d[x,-wtNucls[x]]))
        wtCounts_1d= do.call(rbind, lapply(1:length(wtNucls), function(x) counts_1d[x,wtNucls[x]]))
        
        mut_1 =sapply(ref.df$asInt, function(x) nucl.df$asInt[-x][1])
        mut_2 = sapply(ref.df$asInt, function(x) nucl.df$asInt[-x][2])
        mut_3 = sapply(ref.df$asInt, function(x) nucl.df$asInt[-x][3])
        mut_counts.df=data.frame(wt=rep(ref.df$asChar,3), 
                                 counts=c(mutations_1d[,1], mutations_1d[,2],mutations_1d[,3]), 
                                 mutNo=c(rep("mut1", L), rep("mut2", L),rep("mut3", L)), 
                                 #maxMut=rep(max_mut, 3),
                                 mut=nucl.df$asChar[c(mut_1, mut_2, mut_3)])
        
        #for each row, the mutation with the maximum counts
        max_mut=sapply(seq(L), function(x) which.max(mutations_1d[x,]))
        max_mut_nucl=sapply(seq(L), function(x) nucl.df$asChar[nucl.df$asInt[-ref.df$asInt[x]][ which.max(mutations_1d[x,])]])
        
        #coverage
        coverage_1d = as.vector(rowSums(counts_1d))
        
        #mutation rate
        mutFreq1d = mutations_1d/coverage_1d
        
        
        #df_mutRate = data.frame(position=pos, mutFrequency=apply(mutFreq1d, 1, median), data="1d")
        #df_mutRate= df_mutRate[!is.na(df_mutRate$mutFrequency),]
        
        #total mutation rate per position
        totMutRate=rowSums(mutations_1d)/coverage_1d
        
        # Shannon entropy per position
        mutFreqs_all <- counts_1d/rowSums(counts_1d)
        entropy_perPos <- apply(mutFreqs_all,1, computeShannonEntropy)
        
        # collect all data
        coverage1d_all.df = rbind(coverage1d_all.df, data.frame(index=index, 
                                                                sample=f_name, 
                                                                position=pos,
                                                                wt = wtNucls,
                                                                coverage=coverage_1d, 
                                                                perc_of_maxCov=coverage_1d/max(coverage_1d), 
                                                                wtCount = wtCounts_1d,
                                                                mutCount1 = mutations_1d[,1],
                                                                mutCount2 = mutations_1d[,2],
                                                                mutCount3 = mutations_1d[,3],
                                                                mutRate=totMutRate, 
                                                                mutRate1 = mutFreq1d[,1],
                                                                mutRate2 = mutFreq1d[,2],
                                                                mutRate3 = mutFreq1d[,3],
                                                                entropy=entropy_perPos))
      } #end if file exists
    } # end for 1d files
  } # end dir exsists
  return(coverage1d_all.df)
}# end function


readCountData_2d <- function(countDir, ref.df, nucl.df) {
  L=nrow(ref.df)
  
  countDir2d=paste0(countDir, "/2d")
  coverage2d_all.df=data.frame()
  
  
  if(dir.exists(countDir2d)) {
    files2d = list.files(countDir2d)
    index=0
    for(f in files2d) {
      #filename without extension
      f_name=tools::file_path_sans_ext(f)
      # read in 2d file
      file2d=paste0(countDir2d,"/", f)
      if(file.exists(file2d)) {
        index=index+1
        counts2dTable = read.table(file2d, header = T, sep = "\t")
  
        counts_2d = as.matrix(setNames(counts2dTable[,-(1:2)],NULL))
        dimnames(counts_2d)<-NULL
        #coverage
        coverage_2d_all = as.vector(rowSums(counts_2d))
        #poistions
        pos = unique(counts2dTable$pos1)
        
        coverage2d.df= data.frame()
        #mutFreq2d.df=data.frame()
        for(p in pos) {
          quantiles=quantile(coverage_2d_all[counts2dTable$pos1==p | counts2dTable$pos2==p], probs = c(0.25, 0.5, 0.75))
          maxCov=max(coverage_2d_all[counts2dTable$pos1==p | counts2dTable$pos2==p])
          coverage2d.df= rbind(coverage2d.df, data.frame(position=p, coverage_max=maxCov, coverage_median=quantiles[2], coverage_25=quantiles[1], coverage_75=quantiles[3])) 
          # # # mutation rate
          # mutations_2d=counts_2d[counts2dTable$pos1==p, wtNucls[p]]
          # mutations_2d = do.call(rbind, lapply(1:length(wtNucls), function(x) counts_1d[x,-wtNucls[x]]))
          # quantiles=quantile(coverage_2d_all[counts2dTable$pos1==p | counts2dTable$pos2==p], probs = c(0.05, 0.5, 0.95))
          # coverage2d.df= rbind(coverage2d.df, data.frame(position=p, coverage_median=quantiles[2], coverage_5=quantiles[1], coverage_95=quantiles[3]))
        }
        coverage2d_all.df = rbind(coverage2d_all.df, data.frame(index=index, 
                                                                sample=f_name, 
                                                                position=pos, 
                                                                maxCov=coverage2d.df$coverage_max, 
                                                                covMedian=coverage2d.df$coverage_median, 
                                                                cov25=coverage2d.df$coverage_25, 
                                                                cov75=coverage2d.df$coverage_75))
      } # end if file exists
    }# end for 2d files
  }# end if dir exists
  return(coverage2d_all.df)
} # end function

getNucleotide <- function() {
  require("scales")
  # for illustration of nucleotides
  nucl.df = data.frame(asChar=c("A","C","G","T"), asInt=seq(4), color=hue_pal()(4))
  return(nucl.df)
}

getRefSeq <- function(referenceFile, nucl.df) {
  require(seqinr)
  # read and prepare reference information
  refSeqFasta = read.fasta(referenceFile)
  refSeqChar = toupper(as.vector(refSeqFasta[[1]]))
  refSeq = as.vector(vapply(refSeqChar, function(x) which(nucl.df$asChar == x), numeric(1)))
  refSeq.df <- data.frame(asChar=refSeqChar, asInt=refSeq)
  return(refSeq.df)
}

computeShannonEntropy <- function(mutFreqs){
  m <- mutFreqs[mutFreqs>0]
  return(-sum(m*log(m)))
}



