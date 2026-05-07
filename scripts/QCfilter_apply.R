# Apply QC flag cutoffs to filter contigs from a region BED file.
# Use QCfilter_explore.R to choose cutoffs before running this script.
# Usage: Rscript QCfilter_apply.R <input_bed> <filter_tsv> <filter_type> [FLAG=cutoff ...]
# Example: Rscript QCfilter_apply.R input.bed qc.tsv centromere MISJOIN=4 Col=2
# 2/19/2026: For KaryoScope, recommend to use Err=0, COLLAPSE=0, COLLAPSE_VAR=0
QC_filter <- function(input_path,filter_path,filter_type,
                      Col=Inf,Dup=Inf, Err=Inf, NNN=Inf,COLLAPSE=Inf,
                      COLLAPSE_OTHER=Inf,COLLAPSE_VAR=Inf,MISJOIN=Inf){
  # filter_type=  'centromere', 'subtelo_q', 'subtelo_p'
  # I/O
  in_data <- read.delim(input_path, header = F)
  colnames(in_data) <- c("contig","start","end","feature","chrom")
  out_bed <- sub("\\.bed$", ".pass.bed", basename(input_path))                              
  out_report <- sub("\\.bed$", ".pass.report.txt", basename(input_path))
  filter_tb <- read.delim(filter_path)
  
  # Main
  chrs<- c(paste0("chr",1:22),"chrX","chrY")
  chrs <- factor(chrs, levels = c(paste0("chr", 1:22), "chrX", "chrY"))
  filter_tb <- subset(filter_tb,
                      filter_tb$region == filter_type &
                      filter_tb$X_n_Col <= Col & 
                      filter_tb$X_n_Dup <= Dup &
                      filter_tb$X_n_Err <= Err &
                      filter_tb$X_n_NNN <= NNN &
                      filter_tb$X_n_COLLAPSE <= COLLAPSE &
                      filter_tb$X_n_COLLAPSE_OTHER <= COLLAPSE_OTHER &
                      filter_tb$X_n_COLLAPSE_VAR <= COLLAPSE_VAR &
                      filter_tb$X_n_MISJOIN <= MISJOIN
                      )
  
  for (i in 1:length(chrs)){
    chr <- chrs[i]
    sub_filter <- subset(filter_tb, filter_tb$chromosome == chr)
    sub_in_data <- subset(in_data, in_data$chrom == chr & in_data$contig %in% sub_filter$contig)
    
    if(i == 1){
      out_df <- sub_in_data
    } else {
      out_df <- rbind(out_df, sub_in_data)
    }
  }
  
  write.table(out_df, out_bed, sep = "\t", quote = F, row.names = F, col.names = F)
  
  sink(out_report)
  cat(paste0("Input: ",input_path,"\n"))
  cat(paste0("Output: ",out_bed,"\n"))
  cat(paste0("filter file:", filter_path,"\n"))
  cat("Filter settings\n")
  cat(paste0("Collinear (X_n_Col): ",Col,"\n"))
  cat(paste0("Duplicated (X_n_Dup): ",Dup,"\n"))
  cat(paste0("Erroneous (X_n_Err): ",Err,"\n"))
  cat(paste0("NNN regions (X_n_NNN): ",NNN,"\n"))
  cat(paste0("Collapsed (X_n_COLLAPSE): ",COLLAPSE,"\n"))
  cat(paste0("Collapsed other (X_n_COLLAPSE_OTHER): ",COLLAPSE_OTHER,"\n"))
  cat(paste0("Collapsed with variants (X_n_COLLAPSE_VAR): ",COLLAPSE_VAR,"\n"))
  cat(paste0("Misjoined (X_n_MISJOIN): ",MISJOIN,"\n"))
  sink()
}

# --- CLI ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  cat("Usage: Rscript QCfilter_region.bed.R <input_bed> <filter_tsv> <filter_type> [Col=Inf] [Dup=Inf] [Err=Inf] [NNN=Inf] [COLLAPSE=Inf] [COLLAPSE_OTHER=Inf] [COLLAPSE_VAR=Inf] [MISJOIN=Inf]\n")
  cat("  filter_type: 'centromere', 'subtelo_q', or 'subtelo_p'\n")
  cat("  Only specify flags you want to change, e.g.: MISJOIN=4 Col=2\n")
  quit(status = 1)
}

# Parse named args (e.g. "MISJOIN=4") into a list
cutoffs <- list()
for (a in args[-(1:3)]) {
  kv <- strsplit(a, "=")[[1]]
  if (length(kv) == 2) cutoffs[[kv[1]]] <- as.numeric(kv[2])
}

QC_filter(
  input_path      = args[1],
  filter_path     = args[2],
  filter_type     = args[3],
  Col             = ifelse(is.null(cutoffs$Col), Inf, cutoffs$Col),
  Dup             = ifelse(is.null(cutoffs$Dup), Inf, cutoffs$Dup),
  Err             = ifelse(is.null(cutoffs$Err), Inf, cutoffs$Err),
  NNN             = ifelse(is.null(cutoffs$NNN), Inf, cutoffs$NNN),
  COLLAPSE        = ifelse(is.null(cutoffs$COLLAPSE), Inf, cutoffs$COLLAPSE),
  COLLAPSE_OTHER  = ifelse(is.null(cutoffs$COLLAPSE_OTHER), Inf, cutoffs$COLLAPSE_OTHER),
  COLLAPSE_VAR    = ifelse(is.null(cutoffs$COLLAPSE_VAR), Inf, cutoffs$COLLAPSE_VAR),
  MISJOIN         = ifelse(is.null(cutoffs$MISJOIN), Inf, cutoffs$MISJOIN)
)

# For Dev
# QC_filter(input_path = "../centromere_region_beds/pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.bed",
#           filter_path = "../aggregate_qc_v4.manual_curation.tsv",
#           filter_type = "centromere")
