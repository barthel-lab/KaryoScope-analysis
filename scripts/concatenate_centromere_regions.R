# Concatenate per-chromosome centromere region BED files into a single
# all-chromosome BED file for downstream KaryoScope analysis.

inDir <- "../centromere_region_beds/"
chrs<- c(paste0("chr",1:22),"chrX","chrY")
chrs <- factor(chrs, levels = c(paste0("chr", 1:22), "chrX", "chrY"))

for(i in 1:24){
  chr <- chrs[i]
  bed <- read.delim(paste0(inDir,"pangenome.",chr,".centromere.KS_human_CHM13.presmoothed.region.bed"), header = F)
  bed$V5 <- as.character(chr)
  if(i == 1){
    cat_bed = bed
  } else {cat_bed <- rbind(cat_bed, bed)}
}

write.table(cat_bed, "pangenome.ALLchr.centromere.KS_human_CHM13.presmoothed.region.bed", sep = "\t", quote = F, row.names = F, col.names = F)
