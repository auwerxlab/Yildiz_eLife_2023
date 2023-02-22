# Parameters
n_clust <- 16
sample_dir <- "./Data/RAW/F20FTSEUHT0958_MOUipaE_12sample_lowyield_ece/Filter_SOAPnuke/Clean/"
fasta_ref <- "./Resources/Reference/Mus_musculus.GRCm38.dna.primary_assembly.fa"
gtf_ref <- "./Resources/GTF/Mus_musculus.GRCm38.102.gtf"
bed_ref <- "./Resources/BED/Mus_musculus.GRCm38.102.bed"
star_genome <- "./Data/STARgenome"

# GRCm38 Latest Reference / Release 102
# https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/
# rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz .
# rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz .

# .toplevel - Includes haplotype information (not sure how aligners deal with this)
# .primary_assembly - Single reference base per position

# Retrieve GTF file
# rsync -av rsync://ftp.ensembl.org/ensembl/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz

# ---------------------------------
# Build Reference Genome using STAR
# ---------------------------------
command <- paste0("STAR --runThreadN ", n_clust, 
                  " --runMode genomeGenerate --genomeDir ", star_genome,
                  " --genomeFastaFiles ", fasta_ref, 
                  " --sjdbGTFfile ", gtf_ref, 
                  " --sjdbOverhang 99 --limitGenomeGenerateRAM 35924399488")
System(command, Log_name = "./Log/Genome_Preperation.log", Append = FALSE)

# ------------------
# Convert GTF to Bed
# ------------------
Command_gtf2bed <- paste0(BedOps, "/gtf2bed < ",
                          gtf_ref, " > ", 
                          bed_ref)
# Remove Chr from bed file
Command_GENCODE_bed_mod <- paste0("sed 's/^chr//' ./Resources/BED/mm10_GENCODE_vm25.bed > ./Resources/BED/mm10_GENCODE_vm25_MOD.bed")

# -----------------------
# Convert GTF to R-Object
# -----------------------
ens <- ensemblGenome()
read.gtf(ens, "./Resources/GTF/Mus_musculus.GRCm38.102.gtf")
Ens.gtf <- getGenePositions(ens)

tableFeatures(ens)
table(Ens.gtf$seqid)

# save R Object
save(Ens.gtf, file = "./Data/R_Data/Ens_gtf.RData")

# ------------
# Collapse GTF
# ------------
# using https://github.com/broadinstitute/gtex-pipeline/tree/master/gene_model
# conda activate RNA-seq
# python ./Tools/GTEX_Pipeline/collapse_annotation.py 
# ./Resources/GTF/Mus_musculus.GRCm38.102.gtf Mus_musculus.GRCm38.102.Collapsed.gtf






