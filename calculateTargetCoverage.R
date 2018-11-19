# Christopher Medway
# Rscript to reformat coverage and gaps metrics

seqid    <- commandArgs(trailingOnly = T)[1]
sampleid <- commandArgs(trailingOnly = T)[2]

message(sampleid)

# generated using TSO ROI bed files
coverageIn <- read.table(paste0(seqid,"_",sampleid,"_customClinicalCoverageTargetMetrics.txt"), stringsAsFactors = F)

# extract gene name from custom coverage file
gene <- stringr::str_extract(string = coverageIn$V4, pattern = "^\\w+")

# loop over list of genes and reformat coverage data
coverageInList <- split(coverageIn, gene)
out <- lapply(names(coverageInList), function(g) {
  
  message(g)
  
  df <- coverageInList[[g]]

  vec <- c(
    "V1" = df$V1[1],
    "V2" = min(df$V2),
    "V3" = max(df$V3),
    "V4" = paste0(g,"_allTargets"),
    "V5" = 1,
    "V6" = sum(df$V6),
    "V7" = sum(df$V7),
    "V8" = sum(df$V6) / sum(df$V7) 
    )

  return(rbind(vec, df))  
})

# collapse list object and write
out <- do.call(rbind,out)
out$V8 <- round(as.numeric(out$V8) * 100, digits = 2)
write.table(out, paste0("./",seqid,"_",sampleid,"_customClinicalCoverageTargetCoverage.txt"), row.names = F, sep = "\t", quote = F, col.names=F)

# reformat gaps file
gaps <- read.table(paste0("./",seqid,"_",sampleid,"_customGaps.bed"), stringsAsFactors=F)
gaps[,5:8] <- ""
gaps$V9 <- unlist(lapply(stringr::str_split(gaps$V4, "\\."), function(x) x[1]))
names(gaps) <- c("CHR","START","STOP","TARGET","","","","","GENE")
write.table(gaps, paste0("./",seqid,"_",sampleid,"_customGaps.bed"), row.names = F, sep = "\t", quote = F)
