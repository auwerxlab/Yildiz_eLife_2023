# Parameters
n_clust <- 16
sample_dir <- "./Data/RAW/F20FTSEUHT0958_MOUipaE_12sample_lowyield_ece/Filter_SOAPnuke/Clean/"

# Read Fastq File names
fq_files <- list.files(sample_dir, full.names = F, recursive = T, pattern = "\\.fq.gz$")
fq_files <- paste0(sample_dir, fq_files)
names(fq_files) <- gsub(".fastq", "", gsub("(.)+\\/", "", fq_files))

# --------------------------
# Perform QC for each sample
# --------------------------
# Create Directory for Quality Control
System("mkdir ./Data/QC", Log_name = "./Log/QC.log", Append = FALSE)
# Track Progress (Needed bz of Large Sample Size)
done_dirs <- list.dirs("./Data/QC/", recursive = F)
html <- list.files(done_dirs, pattern = "html")
done_file <- gsub("_fastqc.html", "", html)
length(done_file)
length(fq_files)
remaining <- setdiff(names(fq_files), done_file)
length(remaining)
remaining_files <- fq_files[remaining]

mclapply(mc.cores = n_clust, X = remaining_files, function(f){
  s <- gsub("(.)+\\/", "", f)
  print(paste0("Performing QC for sequencing reads of sample: ", s))
  fqc_dir <- paste0("./Data/QC/", s)
  dir.create(fqc_dir)
  fq_command <- paste0("./Tools/FastQC/fastqc ", f," -o ./Data/QC/", s)
  System(fq_command, Log_name = "./Log/QC.log", Append = TRUE)
  print(paste0("Fast QC complete for ", s))
  return(s)
})

QC_Log_data <- Log_data("./Tools/FastQC/fastqc ./Data/RAW/.*.fq.gz -o", 
                        Log_name = "./Log/QC.log")
mean(QC_Log_data$Exit.Status)

# --------------------------
# MultiQC Report
# --------------------------
System("docker run -t -v `pwd`:`pwd` -w `pwd` -u $(id -u ${USER}):$(id -g ${USER}) ewels/multiqc ./Data --title 'MultiQC Report'", 
       Log_name = "./Log/MultiQC.log", Append = FALSE)




