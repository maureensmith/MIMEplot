# clear the environment variables
rm(list = ls())

# install needed packages
dynamic_require <- function(package){
  if(eval(parse(text=paste("require(",package,")"))))
    return(TRUE)
  
  install.packages(package)
  return (eval(parse(text=paste("require(",package,")"))))
}

#"ggformula"
for(p in c("ggplot2", "seqinr", "stringr" , "ggthemes", "RColorBrewer")) {
  dynamic_require(p)
}


#Read in arguments and normalize paths
args = commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  cat("\nCall the script with at least 2 arguments: MIMEAnToDirectory referenceFile suffix\n
The MIMEAnToDirectory contains the one or several MIMEAnTo result directories with output files.\n
The reference file contains the reference sequence in a fasta file.\n
The plots are saved in the resultDirectory.\n
An optional fourth argument suffix can be given, if the MIMAnTo files were saved with a suffix. \n\n")

  #terminate without saving workspace
  quit("no")
}

cat(c("Arguments: ", args, "\n"), sep = "\n")

# set the absolute paths
data_dir <- normalizePath(args[1])
referenceFile <- normalizePath(args[2])
suffix=""
if(length(args) > 2)
  suffix=args[3]

cat(c("Arguments: ", data_dir, referenceFile, suffix, "\n"), sep = "\n")


# set working directory to call other R Scripts
getScriptPath <- function(){
  cmd.args <- commandArgs()
  m <- regexpr("(?<=^--file=).+", cmd.args, perl=TRUE)
  script.dir <- dirname(regmatches(cmd.args, m))
  if(length(script.dir) == 0) stop("can't determine script dir: please call the script with Rscript")
  if(length(script.dir) > 1) stop("can't determine script dir: more than one '--file' argument detected")
  return(script.dir)
}

this.dir <- getSrcDirectory(function(x) {x})

if (rstudioapi::isAvailable()) {
  this.dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
}else {
  this.dir <- getScriptPath()
}
setwd(this.dir)

source("mimeantoPlotRoutines.R")
source("readCountData.R")

kd.table_all= data.frame()
max_kd.table_all = data.frame()

nucl.df = getNucleotide() 
refSeq.df <- getRefSeq(referenceFile, nucl.df)

evaluations=list.files(data_dir)
evaluations <- evaluations[dir.exists(paste0(data_dir,"/",evaluations))]

for(eval_dir in evaluations) {
  results_dir <- paste0(data_dir,"/",eval_dir)
  kd_file <- paste0(results_dir,"/PositionWiseKdEstimates",suffix,".csv")
  max_kd_file <- paste0(results_dir,"/PositionWiseMaxKd",suffix,".csv")
  if(file.exists(kd_file) & file.exists(max_kd_file)) {
    kd.table <- read.csv(kd_file, header=T, sep="\t")
    kd.table_all <- rbind(kd.table_all, data.frame(kd.table, comparison=eval_dir))
    
    max_kd.table <- read.csv(max_kd_file, header=T, sep="\t")
    
    width=5
    smoothed_median = filter(max_kd.table$median.Kd, rep(1,width)/width, sides=2)
    smoothed_5 = filter(max_kd.table$X5..percentil, rep(1,width)/width, sides=2)
    smoothed_95= filter(max_kd.table$X95..percentil, rep(1,width)/width, sides=2)
    max_kd_smooth.table <- data.frame(max_kd.table, 
                                  smoothedMedian=as.vector(smoothed_median),
                                  smoothed5=as.vector(smoothed_5),
                                  smoothed95=as.vector(smoothed_95),
                                  comparison=eval_dir)
    max_kd.table_all <- rbind(max_kd.table_all, max_kd_smooth.table)
    #############Plot max Kd
    #kd_sub.table <- max_kd.table_all[str_detect(max_kd.table_all$comparison, paste0(eval_dir,"$")),]
    p<-getPlot(max_kd_smooth.table)
    outputFile = paste0(results_dir, "/maxKd_",eval_dir,suffix,".pdf")
    ggsave(filename = outputFile, plot=p, device = "pdf", width=40, height=10, units = "cm")
    
    
    ############ Plot long plot with all Kds
    position_range= max(kd.table$pos1, na.rm = T) - min(kd.table$pos1, na.rm = T)
    p<-getLongPlot(kd.table, refSeq.df)
    outputFile = paste0(results_dir, "/positionwise_Kds_",eval_dir,suffix,".pdf")
    ggsave(filename = outputFile, plot=p, device = "pdf", width=40*position_range/100, height=10, limitsize = FALSE)

  }
}

# Plot all smoothed Kds in one Plot

data.table <- max_kd.table_all

p<- getPlot(data.table)
outputFile = paste0(data_dir, "/maxKd_all.pdf")
ggsave(filename = outputFile, plot=p, device = "pdf", width=40, height=10, units = "cm")
