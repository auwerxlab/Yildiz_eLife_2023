# Parameters
n_clust <- 16
sample_dir <- "./Data/RAW/F20FTSEUHT0958_MOUipaE_12sample_lowyield_ece/Filter_SOAPnuke/Clean/"
fasta_ref <- "./Resources/Reference/Mus_musculus.GRCm38.dna.primary_assembly.fa"
gtf_ref <- "./Resources/GTF/Mus_musculus.GRCm38.102.gtf"
star_genome <- "./Data/STARgenome"

# -----------------------
# Read mapping using STAR
# -----------------------
# raw file names
fq_files <- list.files(sample_dir, full.names = F, recursive = T, pattern = "\\.fq\\.gz$")
fq_files <- paste0(sample_dir, fq_files)
# load metadata from xls sheet
metadata <- data.frame(read_xlsx(path = "Data/Metadata/20200917_RNAseq sample info.xlsx", 
                                 sheet = "Sayfa1", n_max = 13))
metadata$Raw.file <- c("1", "2", "3", "6", "7", "10", "11", "12", "15", "16", "17", "18")
# Add batch information
metadata$Batch <- as.integer(gsub(".*_(\\d*)$", "\\1", metadata$SampleID))

# save R Object
save(metadata, file = "./Data/R_Data/metadata.RData")

# check if all samples are accounted for
read_Pairs <- lapply(metadata$Raw.file, function(x){
  rp <- paste0(x, c("_1","_2"),".fq.gz")
  setNames(sort(fq_files[grepl(paste0("Clean/",x, "/"), fq_files)]), rp)
})
names(read_Pairs) <- metadata$Raw.file
if(any(sapply(read_Pairs, length) != 2)){
  stop("Some samples are missing")
}

# Track Progress
# check which samples are already done
System("mkdir ./Data/fq_temp", Log_name = "./Log/Star_Mapping.log", Append = FALSE)
System("mkdir ./Data/Star_alignments_temp", Log_name = "./Log/Star_Mapping.log", Append = TRUE)
done <- list.files("./Data/Star_alignments/", pattern = "Aligned.sortedByCoord.out.bam")
done_samples <- gsub("refAlign\\.", "", gsub(".Aligned.sortedByCoord.out.bam", "", done))

# Map remaining samples
for( s in metadata$Raw.file[!metadata$Raw.file %in% done_samples]){
  tmp_Sample <- s
  tmp_Pair <- read_Pairs[[s]]
  # make a local copy of fq
  system(paste0("cp ", paste(tmp_Pair, collapse = " "), " ./Data/fq_temp"))
  local_fq <- paste0("./Data/fq_temp/", names(tmp_Pair))
  align_Command <- paste0('STAR --runThreadN ', n_clust, 
                          ' --genomeDir ', star_genome,
                          ' --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --readFilesIn ', 
                          paste(local_fq, collapse = " "),
                          ' --outFileNamePrefix ./Data/Star_alignments_temp/refAlign.',s,'.')
  System(align_Command, Log_name = "./Log/Star_Mapping.log", Append = TRUE)
  star_Output <- list.files("./Data/Star_alignments_temp", pattern = s, full.names = T)  
  star_Output_Copy_Command <- paste0("cp ", star_Output, " ./Data/Star_alignments/")
  sapply(star_Output_Copy_Command, system)
  star_Output_Remove_Command <- paste0("rm -rf ", star_Output)
  sapply(star_Output_Remove_Command, system)
  system(paste0("rm -rf ", paste(local_fq, collapse = " ")))
}

# Check Log for All Jobs
Mapping_Log_data <- Log_data("STAR --runThreadN 16 --genomeDir .* --readFilesIn .*fq_temp/(.*)_[12].fq.gz .* --out.*", 
                             Log_name = "./Log/Star_Mapping.log")
mean(QC_Log_data$Exit.Status)

# ---------------
# Index bam files
# ---------------
allBams <- list.files("./Data/Star_alignments/", pattern = "sortedByCoord.out.bam$", full.names = T)
parallel::mclapply(mc.cores = 8, X = allBams, indexBam)

# -----------
# Read Counts
# -----------
read_count <- list.files(path = "./Data/Star_alignments/")
read_count <- read_count[grepl("ReadsPerGene.out.tab",read_count)]
list_count <- lapply(read_count, function(x){
  count_sample <- read.delim(file = paste0("./Data/Star_alignments/",x),header = FALSE)
  count_sample <- count_sample[5:dim(count_sample)[1],c(1,2)]
  colnames(count_sample) <- c("gene_id","count")
  count_sample$SampleID <- gsub("refAlign.(\\d*).ReadsPerGene.out.tab","\\1",x)
  return(count_sample)
})
count_star <- do.call(rbind, list_count)
count_star <- maditr::dcast(count_star, gene_id ~ SampleID, value.var = "count")
count_mat <- data.matrix(count_star[,-1])
rownames(count_mat) <- count_star$gene_id

# save count matrix
save(count_mat, file = "./Data/R_Data/count_mat.RData")
