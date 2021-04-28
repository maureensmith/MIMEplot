library(ggplot2)
library(ggthemes)

getLongPlot <- function(data.table, refSeq.df) {
  nucls=c("A", "C", "G", "T")
  nucls_idx=seq(4)
  
  plot.table <- data.frame()
  for (i in seq(nrow(data.table))) {
    pos= data.table$pos1[i]
    wt = data.table$wt.base[i]
    mut_idx = nucls_idx[-wt]
    for(m in mut_idx) {
      plot.table <- rbind(plot.table, data.frame(pos=pos, wt=wt, wt_nucl=nucls[wt], mut=m, mut_nucl=nucls[m], 
                                                 kd=data.table[i, 4+(m-1)*7], 
                                                 kd5=data.table[i, 9+(m-1)*7],
                                                 kd95=data.table[i, 10+(m-1)*7],
                                                 pvalue=data.table[i, 5+(m-1)*7]))
    }
  }
  
  plot.table <- na.omit(plot.table)
  pvalue_pos=rep(NA, nrow(plot.table)) 
  pvalue_pos[plot.table$pvalue< 0.05 & log2(plot.table$kd) > 0] = 1
  pvalue_pos[plot.table$pvalue< 0.05 & log2(plot.table$kd) < 0] = -1
  
  min_kd = min(log2(plot.table$kd5), na.rm = T)
  range_kd= abs(max(log2(plot.table$kd95), na.rm = T) - min_kd )
  range_pos= max(plot.table$pos, na.rm = T) - min(plot.table$pos, na.rm = T) 
  
  size_shape=5 #range_kd*1.5
  #size_shape=max(4,range_kd*1.2)
  
  
  p <- ggplot()+
    geom_hline(yintercept = 0, size=1)+
    #geom_point(aes(x=pos, y=log2(kd), col=mut_nucl))+
    geom_pointrange(data=plot.table, aes(x=pos, y=log2(kd), ymin = log2(kd5), ymax =log2(kd95), col=mut_nucl), position = position_dodge(width = 0.5), show.legend = F)+
    geom_point(aes(x=seq(nrow(refSeq.df)), y=min_kd-0.05*range_kd, col=refSeq.df$asChar), shape=refSeq.df$asChar, size=size_shape, show.legend = F)+
    geom_point(data=plot.table, aes(x=pos, y=(min_kd-(0.05*range_kd)+pvalue_pos*range_kd*0.03), col=mut_nucl), shape="*", size=size_shape, position = position_dodge(width = 0.5), show.legend = F)+
    theme_calc()+ scale_colour_calc()+
    ylab("rel. Kd (log2)") +
    xlab("Position") +    
    scale_x_continuous(breaks = pretty(plot.table$pos, n= range_pos/10))+
    scale_y_continuous(breaks = breaks_pretty(n =range_kd))+
    theme(axis.text = element_text(size=20),
          axis.title = element_text(size=20,face="bold"))
  
  return(p)
}

getLongPlot_severalSamples <- function(data.table, refSeq.df) {
  nucls=c("A", "C", "G", "T")
  nucls_idx=seq(4)
  
  plot.table <- data.frame()
  for (i in seq(nrow(data.table))) {
    pos= data.table$pos1[i]
    wt = data.table$wt.base[i]
    mut_idx = nucls_idx[-wt]
    sample=data.table$comparison[i]
    for(m in mut_idx) {
      plot.table <- rbind(plot.table, data.frame(pos=pos, wt=wt, wt_nucl=nucls[wt], mut=m, mut_nucl=nucls[m], 
                                                 kd=data.table[i, 4+(m-1)*7], 
                                                 kd5=data.table[i, 9+(m-1)*7],
                                                 kd95=data.table[i, 10+(m-1)*7],
                                                 pvalue=data.table[i, 5+(m-1)*7],
                                                 sample=sample))
    }
  }
  
  pvalue_pos=rep(NA, nrow(plot.table)) 
  pvalue_pos[plot.table$pvalue< 0.05 & log2(plot.table$kd) > 0] = 1
  pvalue_pos[plot.table$pvalue< 0.05 & log2(plot.table$kd) < 0] = -1
  
  min_kd = min(log2(plot.table$kd5), na.rm = T)
  range_kd= abs(max(log2(plot.table$kd95), na.rm = T) - min_kd )
  range_pos= max(plot.table$pos, na.rm = T) - min(plot.table$pos, na.rm = T) 
  
  #size_shape=4 #range_kd*1.5
  size_shape=range_kd*1.5
  
  
  p <- ggplot()+
    geom_hline(yintercept = 0, size=1)+
    #geom_point(aes(x=pos, y=log2(kd), col=mut_nucl))+
    geom_pointrange(data=plot.table,aes(x=pos, y=log2(kd), ymin = log2(kd5), ymax =log2(kd95), 
                        shape=mut_nucl, col=sample), position = position_dodge(width = 0.5), show.legend = T)+
    #geom_point(aes(x=pos, y=min_kd-0.05*range_kd), shape=plot.table$wt_nucl, size=size_shape, show.legend = F)+
    geom_point(aes(x=seq(nrow(refSeq.df)), y=min_kd-0.05*range_kd), shape=refSeq.df$asChar, size=size_shape, show.legend = F)+
    
    # geom_point(aes(x=pos, y=(min_kd-0.05*range_kd)+pvalue_pos*size_shape/50, shape=mut_nucl), 
    #            size=size_shape, position = position_dodge(width = 0.5), show.legend = F)+
    theme_calc()+ #scale_colour_calc()+
    ylab("rel. Kd (log2)") +
    xlab("Position") +    
    scale_x_continuous(breaks = pretty(plot.table$pos, n= range_pos/10))+
    #scale_y_continuous(breaks = pretty(plot.table$pos, n= range_kd))+
    theme(axis.text = element_text(size=20),
          axis.title = element_text(size=20,face="bold"))
  
  return(p)
}

getPlot <- function(data.table, withRegion=F) {
  library(scales)
  #define important regions where signal is expected
  region = rep("other", nrow(data.table))
  #region_polyA = region_SL2 = rep(NaN,length(positions))
  region[1:51] = "TAR"
  region[58:104] = "PolyA"
  region[119:238] = "PBS"
  region[243:277] = "SL1"
  region[282:300] = "SL2"
  region[312:325] = "SL3"
  
  regionStart=c(1,58,119,243,282,312)
  regionEnd=c(51,104,238,277,300,325)
  
  regionOrder = c("TAR","PolyA","PBS", "SL1", "SL2", "SL3", "other")
  defColors <- hue_pal(l=40)(length(regionOrder)-1)
  colBlack <- "#000000"
  
  
  maxMaxKd = max(log2(data.table$smoothedMedian), na.rm = T) + 1
  ### Max Kd
  p <- ggplot(data=data.table) +
    geom_hline(aes(yintercept=0), size=1)
  
  if(withRegion) {
   p<-  p+ geom_segment(aes(x=which(region=="TAR")[1], xend=tail(which(region=="TAR"),n=1), y=maxMaxKd, yend=maxMaxKd), size=13, col=defColors[1]) +
      geom_segment(aes(x=which(region=="PolyA")[1], xend=tail(which(region=="PolyA"),n=1), y=maxMaxKd, yend=maxMaxKd), size=13, col=defColors[2]) +
      geom_segment(aes(x=which(region=="PBS")[1], xend=tail(which(region=="PBS"),n=1), y=maxMaxKd, yend=maxMaxKd), size=13, col=defColors[3]) +
      geom_segment(aes(x=which(region=="SL1")[1], xend=tail(which(region=="SL1"),n=1), y=maxMaxKd, yend=maxMaxKd), size=13, col=defColors[4]) +
      geom_segment(aes(x=which(region=="SL2")[1], xend=tail(which(region=="SL2"),n=1), y=maxMaxKd, yend=maxMaxKd), size=13, col=defColors[5]) +
      geom_segment(aes(x=which(region=="SL3")[1], xend=tail(which(region=="SL3"),n=1), y=maxMaxKd, yend=maxMaxKd), size=13, col=defColors[6]) +
      geom_text(aes(label="TAR", x=median(which(region=="TAR")), y= maxMaxKd-0.5), col=defColors[1])+
      geom_text(aes(label="PolyA", x=median(which(region=="PolyA")), y= maxMaxKd-0.5), col=defColors[2])+
      geom_text(aes(label="PBS", x=median(which(region=="PBS")), y= maxMaxKd-0.5), col=defColors[3])+
      geom_text(aes(label="SL1", x=median(which(region=="SL1")), y= maxMaxKd-0.5), col=defColors[4])+
      geom_text(aes(label="SL2", x=median(which(region=="SL2")), y= maxMaxKd-0.5), col=defColors[5])+
      geom_text(aes(label="SL3", x=median(which(region=="SL3")), y= maxMaxKd-0.5), col=defColors[6])
  }
  p<-p+geom_ribbon(aes(x=pos1, ymin=log2(smoothed5), ymax=log2(smoothed95), fill=comparison, group=comparison), alpha=0.3)+
    geom_line(aes(x=pos1, y=log2(smoothedMedian), col=comparison), size=1, alpha=0.8)+
    ylab("rel. Kd (log2)") +
    xlab("Position") +
    theme_bw()+
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=14,face="bold"))
  return(p)
}
