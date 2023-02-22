# #####################################################################################################################################
# 1. Heatmap for the hepatocyte markers, progenitor-associated markers, and BEC markers
# #####################################################################################################################################
compare.genes <- TRUE

hepatocyte.markers <- c("CYP1A1", "ASGR1", "ALB", "CYP3A11", "SERPINA7") # SERPINA7 orthologue of TBG
BEC.markers <- c("HNF1B", "SOX9", "SPP1", "PROM1", "EPCAM", "KRT19", "CD44")

all.markers <- c(hepatocyte.markers, BEC.markers)

# Load expression data
load("./Data/R_Data/Exp.RData", verbose = TRUE)
#exp$counts <- exp$counts[edgeR::filterByExpr(exp, design = exp$design.Treatment),]
exp.tmm.cpm <- log2(edgeR::cpm(exp) + 1)

# Retrieve all E2f genes
genes.info <- exp$gtf[toupper(exp$gtf$gene_name) %in% c(hepatocyte.markers, BEC.markers), 
                      c('gene_name', 'gene_id')]

genes.counts <- data.frame(exp.tmm.cpm[genes.info$gene_id,])
rownames(genes.counts) <- genes.info[match(rownames(genes.counts), genes.info$gene_id), 'gene_name']
samples.order <- colnames(genes.counts)[order(colnames(genes.counts))]

write.csv(genes.counts, "./Reports/markers_counts.csv", quote = FALSE, row.names = FALSE)

# Calculate Z-score
if(compare.genes){
  genes.counts <- data.frame(apply(genes.counts, 2, Zscore))
} else {
  genes.counts <- data.frame(t(apply(genes.counts, 1, Zscore)))
}

genes.counts$gene <- tolower(rownames(genes.counts))
genes.counts <- reshape2::melt(genes.counts, id = "gene")
genes.counts$cond <- gsub("(.*)_.*", "\\1", genes.counts$variable)
genes.counts$variable <- factor(genes.counts$variable, levels = samples.order)
genes.counts$gene <- factor(genes.counts$gene, levels = rev(tolower(all.markers[all.markers %in% toupper(genes.counts$gene)])))

g <- ggplot(genes.counts, aes(x = variable, y = gene, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",
                       high="red", space ="Lab" )

if(compare.genes){
  ggsave(g, filename = "./Plots/Rebuttal/In_vivo_markers_heatmap_genes.pdf", width = 10, height = 6)
} else {
  ggsave(g, filename = "./Plots/Rebuttal/In_vivo_markers_heatmap_treatment.pdf", width = 10, height = 6)
}


# #####################################################################################################################################
# 2. DR markers fold change box plot
# #####################################################################################################################################
DR.markers <- c("SOX9", "EPCAM", "PROM1", "CXCR4", "NCAM1", "TNFRSF12A", "TACSTD2",
                "KRT23", "KRT7", "KRT19") #"NCAM1" Prom1

# Load expression data
load("./Data/R_Data/Exp.RData", verbose = TRUE)
#exp$counts <- exp$counts[edgeR::filterByExpr(exp, design = exp$design.Treatment),]
exp.tmm.cpm <- log2(edgeR::cpm(exp) + 1)

# Retrieve all E2f genes
genes.info <- exp$gtf[toupper(exp$gtf$gene_name) %in% DR.markers, c('gene_name', 'gene_id')]

genes.counts <- data.frame(exp.tmm.cpm[genes.info$gene_id,])
rownames(genes.counts) <- genes.info[match(rownames(genes.counts), genes.info$gene_id), 'gene_name']
genes.counts$gene <- tolower(rownames(genes.counts))
genes.counts <- reshape2::melt(genes.counts, id = "gene")
genes.counts$cond <- gsub("(.*)_.*", "\\1", genes.counts$variable)

write.table(genes.counts, file = "./Reports/S2_table.tsv")

g <- ggplot(genes.counts, aes(x = gene, y = value, fill = cond)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', position = "dodge") +
  facet_wrap(~ gene, scales = "free", nrow = 1) +
  labs(y = "cpm", x = "Gene")+
  theme_classic() +
  ggpubr::stat_compare_means(method = "t.test")
ggsave(g, filename = "./Plots/Rebuttal/In_vivo_DC_markers_boxplot.pdf", width = 16, height = 3.5)


# #####################################################################################################################################
# 3. Epithelial marker changes
# #####################################################################################################################################
# GO:0043236 (laminin binding)
# GO:0005540 (hyaluronic acid binding)
# GO:0045296 (cadherin binding)
# GO:0044331 (cell-cell adhesion mediated by cadherin)
# GO:0033631 (cell-cell adhesion mediated by integrin)
# GO:0016339 (calcium-dependent cell-cell adhesion via plasma membrane cell adhesion molecules)
# GO:0016338 (calcium-independent cell-cell adhesion via plasma membrane cell-adhesion molecules)
# GO:0090136 (epithelial cell adhesion)

go.selected <- c("GOMF_LAMININ_BINDING",
"GOMF_HYALURONIC_ACID_BINDING",
"GOMF_CADHERIN_BINDING",
"GOBP_CELL_CELL_ADHESION_MEDIATED_BY_CADHERIN",
"GOBP_CELL_CELL_ADHESION_MEDIATED_BY_INTEGRIN",
"GOBP_CALCIUM_DEPENDENT_CELL_CELL_ADHESION_VIA_PLASMA_MEMBRANE_CELL_ADHESION_MOLECULES",
"GOBP_CALCIUM_INDEPENDENT_CELL_CELL_ADHESION_VIA_PLASMA_MEMBRANE_CELL_ADHESION_MOLECULES",
"GOBP_EPITHELIAL_CELL_CELL_ADHESION")

GO_gmt <- clusterProfiler::read.gmt("./Resources/GeneSets/c5.go.v7.5.1.symbols.gmt")

# Load the expression data
load("./Data/R_Data/Exp_DESeq2.RData", verbose = T)

go.selected.genes <- mapply(geneset = go.selected, FUN = function(geneset){
  GO_gmt[GO_gmt$term == geneset,]
}, SIMPLIFY = FALSE)

# Create dataframe with expression of genes
gene.list <- unique(unlist(lapply(go.selected.genes, function(x){unlist(x$gene)})))
gene.list <- gene.list[gene.list %in% toupper(DESeq2.Results.Table$gene_name)]
gene.list.exp <- DESeq2.Results.Table[match(gene.list, toupper(DESeq2.Results.Table$gene_name)),]
gene.list.exp <- gene.list.exp[abs(gene.list.exp$log2FoldChange) > 0.5,]
gene.list.exp <- gene.list.exp[!is.na(gene.list.exp$log2FoldChange),]z

df <- data.frame(gene = gene.list.exp$gene_name)

for(geneset in go.selected){
  df[, geneset] <- gene.list.exp$log2FoldChange
  # Set to NA genes not in geneset
  df[!toupper(df$gene) %in% go.selected.genes[[geneset]]$gene, geneset] <- 0
}

rownames(df) <- df$gene
df <- t(df[, -1])

pdf("./Plots/Rebuttal/In_vivo_Epithelial_marker_genesets_Heatmap.pdf", 
    height = 5, width = 15, useDingbats=FALSE)
Heatmap(df, row_names_gp = gpar(fontsize = 9),
        clustering_method_columns = "median")
dev.off()


# #####################################################################################################################################
# 4. qPCR Z-score plot
# #####################################################################################################################################
qPCR.data <- data.frame(readxl::read_xlsx("./Reports/Rebuttal_Data_from_Ece/CD-HFD org HEP&amp;CHOL expression.xlsx", skip = 1))
rownames(qPCR.data) <- qPCR.data$Replicate..

qPCR.data <- qPCR.data[, -1]

# Calculate Z-score
genes.counts <- data.frame(apply(qPCR.data, 2, Zscore))
# Invert the z-score
genes.counts <- -1 * genes.counts

genes.counts$gene <- tolower(rownames(genes.counts))
genes.counts <- reshape2::melt(genes.counts, id = "gene")
genes.counts$cond <- gsub("(.*)_.*", "\\1", genes.counts$variable)
genes.counts$variable <- factor(genes.counts$variable, levels = c("CD1", "CD2", "CD3", "CD4", "CD5",
                                                                  "HFD1", "HFD2", "HFD3", "HFD4", "HFD5"))
genes.counts$gene <- factor(genes.counts$gene, 
                            levels = tolower(rownames(qPCR.data)))

g <- ggplot(genes.counts, aes(x = variable, y = gene, fill = value)) + 
  geom_tile() +
  scale_fill_gradient2(midpoint=0, low="blue", mid="white",
                       high="red", space ="Lab" )

if(compare.genes){
  ggsave(g, filename = "./Plots/Rebuttal/In_vivo_markers_heatmap_genes.pdf", width = 10, height = 6)
} else {
  ggsave(g, filename = "./Plots/Rebuttal/In_vivo_markers_heatmap_treatment.pdf", width = 10, height = 6)
}











































