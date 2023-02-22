NJOBS=16
RAW_FOLDER="./Data/RAW/F20FTSEUHT0958_MOUipaE_12sample_lowyield_ece/Filter_SOAPnuke/Clean/"
FASTA_REF="./Resources/Reference/Mus_musculus.GRCm38.dna.primary_assembly.fa"
GTF_REF="./Resources/GTF/Mus_musculus.GRCm38.102.gtf"
BED_REF="./Resources/BED/mm10_GENCODE_vm25_MOD.bed"

# Create QC Read Mapping folder
mkdir ./Data/QC_Mapping

# Retireve sample IDs
SAMPLEID="$(ls $RAW_FOLDER)"

# QualiMap
# http://qualimap.conesalab.org/doc_html/analysis.html
# -----
for SAMPLE in $SAMPLEID
do
echo "Analyzing sample $SAMPLE"
qualimap rnaseq --java-mem-size=4G \
-bam ./Data/Star_alignments/refAlign.${SAMPLE}.Aligned.sortedByCoord.out.bam \
-gtf ./Resources/GTF/Mus_musculus.GRCm38.102.gtf \
-outdir ./Data/QC_Mapping/QualiMap/${SAMPLE}
done

# RNA-SeQC
# https://github.com/getzlab/rnaseqc
# -----
for SAMPLE in $SAMPLEID
do
echo "Analyzing sample $SAMPLE"
./Tools/rnaseqc.v2.4.2.linux/rnaseqc ./Resources/GTF/Mus_musculus.GRCm38.102.Collapsed.gtf \
./Data/Star_alignments/refAlign.${SAMPLE}.Aligned.sortedByCoord.out.bam \
./Data/QC_Mapping/RNA-SeQC/ -s refAlign.${SAMPLE}
done

# ----------------------------
# Read mapping Quality Control
# ----------------------------
# RSeQC 
# http://rseqc.sourceforge.net/#tin-py
# -----
mkdir ./Data/QC_Mapping/RSeQC

# Read Distribution
mkdir ./Data/QC_Mapping/RSeQC/read_distribution

for SAMPLE in $SAMPLEID
do
echo "Analyzing sample $SAMPLE"
read_distribution.py -r ${BED_REF} -i ./Data/Star_alignments/refAlign.${SAMPLE}.Aligned.sortedByCoord.out.bam > \
./Data/QC_Mapping/RSeQC/read_distribution/${SAMPLE}.RD.log &
done

# Bam Stat
mkdir ./Data/QC_Mapping/RSeQC/bam_stat
for SAMPLE in $SAMPLEID
do
echo "Analyzing sample $SAMPLE"
bam_stat.py -i ./Data/Star_alignments/refAlign.${SAMPLE}.Aligned.sortedByCoord.out.bam > \
./Data/QC_Mapping/RSeQC/bam_stat/refAlign.${SAMPLE} &
done

# Junction Annotation
mkdir ./Data/QC_Mapping/RSeQC/junction_annotation
for SAMPLE in $SAMPLEID
do
echo "Analyzing sample $SAMPLE"
junction_annotation.py -i ./Data/Star_alignments/refAlign.${SAMPLE}.Aligned.sortedByCoord.out.bam -r ${BED_REF} -o \
./Data/QC_Mapping/RSeQC/junction_annotation/${SAMPLE}.JA 2> ./Data/QC_Mapping/RSeQC/junction_annotation/${SAMPLE}.JA &
done

# Gene Body Coverage
mkdir ./Data/QC_Mapping/RSeQC/geneBody_coverage
for SAMPLE in $SAMPLEID
do
echo "Analyzing sample $SAMPLE"
geneBody_coverage.py -r ${BED_REF} -i ./Data/Star_alignments/refAlign.${SAMPLE}.Aligned.sortedByCoord.out.bam -o \
./Data/QC_Mapping/RSeQC/geneBody_coverage/${SAMPLE}.geneBodyCoverage.txt &
done

# Inner Distance Frequency
mkdir ./Data/QC_Mapping/RSeQC/inner_distance
for SAMPLE in $SAMPLEID
do
echo "Analyzing sample $SAMPLE"
inner_distance.py -r ${BED_REF} -i ./Data/Star_alignments/refAlign.${SAMPLE}.Aligned.sortedByCoord.out.bam -o \
./Data/QC_Mapping/RSeQC/inner_distance/${SAMPLE}.inner_distance_freq.txt &
done

# Infer Experiment
mkdir ./Data/QC_Mapping/RSeQC/infer_experiment
for SAMPLE in $SAMPLEID
do
echo "Analyzing sample $SAMPLE"
infer_experiment.py -r ${BED_REF} -i ./Data/Star_alignments/refAlign.${SAMPLE}.Aligned.sortedByCoord.out.bam > \
./Data/QC_Mapping/RSeQC/infer_experiment/${SAMPLE}.infer_experiment.txt &
done

# Infer Experiment
mkdir ./Data/QC_Mapping/RSeQC/tin
for SAMPLE in $SAMPLEID
do
echo "Analyzing sample $SAMPLE"
tin.py -r ${BED_REF} -i ./Data/Star_alignments/refAlign.${SAMPLE}.Aligned.sortedByCoord.out.bam > \
./Data/QC_Mapping/RSeQC/tin/${SAMPLE}.summary.txt &
done

