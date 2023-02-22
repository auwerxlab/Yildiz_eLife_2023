# #####################################################################################################################################
# 1. GSEA Plot Top / Down 15 Changing Gene Sets
# #####################################################################################################################################
# Load GSEA results to plot
load("./Data/R_Data/Exp_DESeq2_GSEA.RData", verbose = T)

gsea.go <- DE.gsea$GO@result
gsea.go$Type <- gsub("GO_(.*)_MM_.*", "\\1", DE.gsea$GO@result$ID)

# Proper GSEA Plot
result <- gsea.go[gsea.go$Type %in% c("MF", "BP"),]
result$core_enrichmentNgenes <- sapply(strsplit(unlist(result$core_enrichment), 
                                                split = "/", fixed = T), length)
result$core_enrichmentgenes <- unlist(sapply(strsplit(unlist(result$core_enrichment), 
                                                      split = "/", fixed = T), 
                                             function(x){paste(sort(x), collapse = " ")}))
result$signficant <- result$qvalues <= 0.01
result$label <- gsub("\\(.*\\)", "", gsub("_", " ", result$Description))
result$gene_ratio <- result$core_enrichmentNgenes / result$setSize

# sort by abs value of NES and take top 15 values
tmp_data <- head(result[order(result$qvalues),], 15)

GSEA.plot <-
  ggplot(tmp_data, aes(x = NES, y = forcats::fct_reorder(label, NES, .desc = FALSE))) + 
  geom_point(aes(size = abs(gene_ratio), col =  qvalues )) +
  scale_colour_gradient("q-value", limits=c(0, 0.2), low="chartreuse3", high = "brown1") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("NES")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )

ggsave(plot = GSEA.plot, 
       paste0("Plots/Alternative_Figures/", "GSEA_Top_15_qValue_MF_BP.pdf"), 
       height = 6, width = 8)


# #####################################################################################################################################
# 2. Specific Gene Set Heatmap Plots
# #####################################################################################################################################
GO_gmt <- read.gmt("./Resources/GeneSets/c5.go.v7.5.1.symbols.gmt")
# Load the expression data
load("./Data/R_Data/Exp_DESeq2.RData", verbose = T)
# Get list of genes in the 2 genesets GOBP_LIPID_IMPORT_INTO_CELL / GOBP_FATTY_ACID_METABOLIC_PROCESS
GOBP_LIPID_IMPORT_INTO_CELL <- GO_gmt[GO_gmt$term == "GOBP_LIPID_IMPORT_INTO_CELL", ]
GOBP_FATTY_ACID_METABOLIC_PROCESS <- GO_gmt[GO_gmt$term == "GOBP_FATTY_ACID_METABOLIC_PROCESS",]

# Create dataframe with expression of genes
gene.list <- unique(c(GOBP_LIPID_IMPORT_INTO_CELL[, 2], GOBP_FATTY_ACID_METABOLIC_PROCESS[, 2]))
gene.list <- gene.list[gene.list %in% toupper(DESeq2.Results.Table$gene_name)]
gene.list.exp <- DESeq2.Results.Table[match(gene.list, toupper(DESeq2.Results.Table$gene_name)),]
gene.list.exp <- gene.list.exp[abs(gene.list.exp$log2FoldChange) > 0.5,]
gene.list.exp <- gene.list.exp[!is.na(gene.list.exp$log2FoldChange),]
df <- 
data.frame(gene = gene.list.exp$gene_name,
           Lipid_import = gene.list.exp$log2FoldChange,
           Fatty_acid = gene.list.exp$log2FoldChange)
# Set to NA genes not in geneset
df[!toupper(df$gene) %in% GOBP_LIPID_IMPORT_INTO_CELL$gene, 2] <- 0
df[!toupper(df$gene) %in% GOBP_FATTY_ACID_METABOLIC_PROCESS$gene, 3] <- 0

pdf("./Plots/Alternative_Figures/Lipid_GO_Heatmap.pdf", 
    height = 10, width = 3, useDingbats=FALSE)
Heatmap(df[, c(2, 3)], row_labels = df$gene, row_names_gp = gpar(fontsize = 9))
dev.off()


# #####################################################################################################################################
# 3. Cnet Plots for specific genesets of interest
# #####################################################################################################################################
# Load GSEA results to plot
load("./Data/R_Data/Exp_DESeq2_GSEA.RData", verbose = T)
gsea.go <- DE.gsea$GO@result
gsea.go$Type <- gsub("GO_(.*)_MM_.*", "\\1", DE.gsea$GO@result$ID)

# Up Cnet
tmp <- DE.gsea$GO
tmp@result <- tmp@result[match(c("GO_BP_MM_CELL_DIVISION", "GO_BP_MM_DNA_REPAIR", "GO_BP_MM_DNA_REPLICATION"), tmp@result$ID),]

g <- clusterProfiler::cnetplot(tmp, foldChange = tmp@geneList, showCategory = 3, circular = TRUE, colorEdge = TRUE)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Cnet_UP.pdf"), 
       height = 12, width = 16,  useDingbats=FALSE)

# Down Cnet
# Up Cnet
tmp <- DE.gsea$GO
tmp@result <- tmp@result[match(c("GO_BP_MM_EXTRACELLULAR_MATRIX_ORGANIZATION", 
                                 "GO_BP_MM_CELL_ADHESION", 
                                 "GO_MF_MM_GLUTATHIONE_TRANSFERASE_ACTIVITY"), tmp@result$ID),]

g <- clusterProfiler::cnetplot(tmp, foldChange = tmp@geneList, showCategory = 3, circular = TRUE, colorEdge = TRUE)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Cnet_Down.pdf"), 
       height = 12, width = 16,  useDingbats=FALSE)


# #####################################################################################################################################
# 4. Joint Enrichment Plot for TF Enrichment HFD / DDC
# #####################################################################################################################################
# Load analysis results
DDC_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/DDC_Enrichment_Overlap_Analysis.RDS")

DDC.TF.Up <- DDC_Enrichment_Overlap_Analysis$DDC_up$TF
DESeq2.TF.Up <- DDC_Enrichment_Overlap_Analysis$DESeq2_up$TF

DDC.TF.Up$comp <- "DDC_Up"
DESeq2.TF.Up$comp <- "DESeq2_Up"

df_merged <- rbind(DDC.TF.Up, DESeq2.TF.Up)
df_merged$signficant <- df_merged$p.adjust <= 0.05
df_merged$Description <- tolower(gsub("TFACTS_MM_","", df_merged$Description))
df_merged$Gene_Ratio <- sapply(df_merged$GeneRatio, function(x) round(eval(parse(text=x)), digits = 2))
df_merged <- df_merged[rownames(df_merged)[!rownames(df_merged) %in% "TF_MM_WINGENDER_E2F-1"],]

# plot only significant ones
df_merged <- df_merged[df_merged$signficant,]

g <- ggplot(df_merged, aes(x = Gene_Ratio, y = forcats::fct_reorder(Description, Gene_Ratio, .desc = FALSE))) + 
  geom_point(aes(size = Count, col =  p.adjust )) +
  scale_colour_gradient("p.adjust", limits=c(0, 0.05), low="red", high = "blue") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("Gene Ratio")) +
  facet_wrap(~comp) +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_DDC_TF_DESeq2_Up_Merged.pdf"), 
       height = 4, width = 6)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_DDC_TF_DESeq2_Up_Merged.svg"), 
       height = 4, width = 6)

# DOWN TF
# Load analysis results
DDC_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/DDC_Enrichment_Overlap_Analysis.RDS")

DDC.TF.Down <- DDC_Enrichment_Overlap_Analysis$DDC_down$TF
DESeq2.TF.Down <- DDC_Enrichment_Overlap_Analysis$DESeq2_down$TF

DDC.TF.Down$comp <- "DDC_Down"
DESeq2.TF.Down$comp <- "DESeq2_Down"

df_merged <- rbind(DDC.TF.Down, DESeq2.TF.Down)
df_merged$signficant <- df_merged$p.adjust <= 0.05
df_merged$Description <- tolower(gsub("TFACTS_MM_","", df_merged$Description))
df_merged$Gene_Ratio <- sapply(df_merged$GeneRatio, function(x) round(eval(parse(text=x)), digits = 2))

# plot only significant ones
df_merged <- df_merged[df_merged$signficant,]

g <- ggplot(df_merged, aes(x = Gene_Ratio, y = forcats::fct_reorder(Description, Gene_Ratio, .desc = FALSE))) + 
  geom_point(aes(size = Count, col =  p.adjust )) +
  scale_colour_gradient("p.adjust", limits=c(0, 0.05), low="red", high = "blue") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("Gene Ratio")) +
  facet_wrap(~comp) +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_DDC_DESeq2_TF_Down_Merged.pdf"), 
       height = 4, width = 6)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_DDC_DESeq2_TF_Down_Merged.svg"), 
       height = 4, width = 6)


# #####################################################################################################################################
# 5. Joint Enrichment Plot for GO Enrichment HFD / DDC
# #####################################################################################################################################
# Load analysis results
DDC_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/DDC_Enrichment_Overlap_Analysis.RDS")

DDC.GO.Up <- DDC_Enrichment_Overlap_Analysis$DDC_up$GO
DESeq2.GO.Up <- DDC_Enrichment_Overlap_Analysis$DESeq2_up$GO

DDC.GO.Up$comp <- "DDC_Up"
DESeq2.GO.Up$comp <- "DESeq2_Up"

df_merged <- rbind(DDC.GO.Up[1:20,], DESeq2.GO.Up[1:20,])
df_merged$signficant <- df_merged$p.adjust <= 0.05
df_merged$Description <- tolower(gsub("GO_.._MM_","", df_merged$Description))
df_merged$Gene_Ratio <- sapply(df_merged$GeneRatio, function(x) round(eval(parse(text=x)), digits = 2))

# plot only significant ones
df_merged <- df_merged[df_merged$signficant,]

g <- ggplot(df_merged, aes(x = Gene_Ratio, y = forcats::fct_reorder(Description, Gene_Ratio, .desc = FALSE))) + 
  geom_point(aes(size = Count, col =  p.adjust )) +
  scale_colour_gradient("p.adjust", limits=c(0, 0.05), low="red", high = "blue") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("Gene Ratio")) +
  facet_wrap(~comp) +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_DDC_GO_DESeq2_Up_Merged.pdf"), 
       height = 4, width = 8)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_DDC_GO_DESeq2_Up_Merged.svg"), 
       height = 4, width = 8)

# DOWN
# Load analysis results
DDC_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/DDC_Enrichment_Overlap_Analysis.RDS")

DDC.GO.Down <- DDC_Enrichment_Overlap_Analysis$DDC_down$GO
DESeq2.GO.Down <- DDC_Enrichment_Overlap_Analysis$DESeq2_down$GO

DDC.GO.Down$comp <- "DDC_Down"
DESeq2.GO.Down$comp <- "DESeq2_Down"

df_merged <- rbind(DDC.GO.Down[1:20,], DESeq2.GO.Down[1:20,])
df_merged$signficant <- df_merged$p.adjust <= 0.05
df_merged$Description <- tolower(gsub("GO_.._MM_","", df_merged$Description))
df_merged$Gene_Ratio <- sapply(df_merged$GeneRatio, function(x) round(eval(parse(text=x)), digits = 2))

# plot only significant ones
df_merged <- df_merged[df_merged$signficant,]

g <- ggplot(df_merged, aes(x = Gene_Ratio, y = forcats::fct_reorder(Description, Gene_Ratio, .desc = FALSE))) + 
  geom_point(aes(size = Count, col =  p.adjust )) +
  scale_colour_gradient("p.adjust", limits=c(0, 0.05), low="red", high = "blue") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("Gene Ratio")) +
  facet_wrap(~comp) +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_DDC_GO_DESeq2_Down_Merged.pdf"), 
       height = 4, width = 10)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_DDC_GO_DESeq2_Down_Merged.svg"), 
       height = 4, width = 10)


# #####################################################################################################################################
# 6. TF enrichment heatmap 
# #####################################################################################################################################
TF_gmt <- tail(read.gmt("./Resources/GeneSets/MousePath_TF_gmt.gmt"), -1)
# Load the expression data
load("./Data/R_Data/Exp_DESeq2.RData", verbose = T)
# Get list of genes in the 2 genesets GOBP_LIPID_IMPORT_INTO_CELL / GOBP_FATTY_ACID_METABOLIC_PROCESS
TFACTS_MM_E2F1 <- TF_gmt[TF_gmt$term == "TFACTS_MM_E2F1", ]
TFACTS_MM_E2F4 <- TF_gmt[TF_gmt$term == "TFACTS_MM_E2F4", ]
TFACTS_MM_E2F2 <- TF_gmt[TF_gmt$term == "TFACTS_MM_E2F2", ]
TFACTS_MM_E2F3 <- TF_gmt[TF_gmt$term == "TFACTS_MM_E2F3", ]

# Create dataframe with expression of genes
gene.list <- unique(c(TFACTS_MM_E2F1[, 2], 
                      TFACTS_MM_E2F4[, 2],
                      TFACTS_MM_E2F2[, 2],
                      TFACTS_MM_E2F3[, 2]))
gene.list <- gene.list[gene.list %in% toupper(DESeq2.Results.Table$gene_name)]
gene.list.exp <- DESeq2.Results.Table[match(gene.list, toupper(DESeq2.Results.Table$gene_name)),]
gene.list.exp <- gene.list.exp[abs(gene.list.exp$log2FoldChange) > 0.5,]
gene.list.exp <- gene.list.exp[!is.na(gene.list.exp$log2FoldChange),]
df <- 
  data.frame(gene = gene.list.exp$gene_name,
             E2F1 = gene.list.exp$log2FoldChange,
             E2F4 = gene.list.exp$log2FoldChange,
             E2F2 = gene.list.exp$log2FoldChange,
             E2F3 = gene.list.exp$log2FoldChange)
# Set to NA genes not in geneset
df[!toupper(df$gene) %in% TFACTS_MM_E2F1$gene, 2] <- 0
df[!toupper(df$gene) %in% TFACTS_MM_E2F4$gene, 3] <- 0
df[!toupper(df$gene) %in% TFACTS_MM_E2F2$gene, 4] <- 0
df[!toupper(df$gene) %in% TFACTS_MM_E2F3$gene, 5] <- 0

pdf("./Plots/Alternative_Figures/E2F_Heatmap_2.pdf", 
    height = 6, width = 5, useDingbats=FALSE)
Heatmap(df[, c(2, 3, 4, 5)], row_labels = df$gene, row_names_gp = gpar(fontsize = 9))
dev.off()


# #####################################################################################################################################
# 7. Chosen Fat Genes Heatmap
# #####################################################################################################################################
# Load the expression data
load("./Data/R_Data/Exp_DESeq2.RData", verbose = T)
# Chosen Lipid genes to plot
Lipid.Genes <- c("cd36", "acot1", "pdk4", "hmgcs2", "angptl4", "aldh1a1", "acot2", "cpt1a", "crot", "acot4", "scd1")

# Create dataframe with expression of genes
gene.list <- unique(Lipid.Genes)
gene.list <- gene.list[gene.list %in% tolower(DESeq2.Results.Table$gene_name)]
gene.list.exp <- DESeq2.Results.Table[match(gene.list, tolower(DESeq2.Results.Table$gene_name)),]
df <- 
  data.frame(gene = gene.list.exp$gene_name,
             Lipid.Genes = gene.list.exp$log2FoldChange,
             p.value = gene.list.exp$padj,
             y = "Lipid.Genes")
df$label <- ""
df$label[df$p.value < 0.05] <- "*"
# Order by decreasing value
df <- df[with(df, order(-Lipid.Genes)),]
df$gene <- factor(df$gene, levels = df$gene)

g <- ggplot(df, aes(y, gene, fill= Lipid.Genes)) + 
  geom_tile() +
  geom_text(aes(label= label),
            position = position_dodge(width=0),  size=5) +
  scale_fill_gradient2(low = "blue", mid = "white",high = "red") +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 12, angle = 90),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.spacing.x=unit(0, "lines"),
        panel.spacing.y=unit(0, "lines"),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 5),
        strip.text.x = element_text(size = 7, 
                                    face = "bold", color = "black"),
        legend.position="left") +
  coord_fixed(ratio = 1/4) + 
  xlab("") +
  scale_y_discrete(position = "right", expand = expand_scale(mult = c(0, 0))) +
  scale_x_discrete(expand = expand_scale(mult = c(0, 0)))

pdf("./Plots/Alternative_Figures/Selected_Lipid_Genes_Heatmap.pdf", 
    height = 3, width = 2.5, useDingbats=FALSE)
g
dev.off()


# #####################################################################################################################################
# 8. Boxplot of ncam1 gene in CDvsHFD
# #####################################################################################################################################
# Load the expression data
load("./Data/R_Data/Exp.RData", verbose = T)

# CPM & Log transform counts
data.Expr <- log2(edgeR::cpm(exp) + 1)
# Get only expression of the gene of interest
data.Expr <- data.Expr['ENSMUSG00000039542',]
# Prepare dataframe for plotting
Expr.df.ncam1 <-
data.frame(sample = names(data.Expr),
           `log(cpm + 1)` = data.Expr,
           condition = gsub('(.*)_.*', '\\1', names(data.Expr)))
g <- 
ggplot(Expr.df.ncam1, aes(x=condition, y=`log.cpm...1.`, fill = condition)) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="Log2(cpm + 1)", x = "Condition")+
  theme_classic() +
  ggpubr::stat_compare_means(method = "t.test")

ggsave(plot = g, paste0("./Plots/Alternative_Figures/Ncam1_Boxplot.pdf"), 
       height = 4, width = 4)


# #####################################################################################################################################
# 9. Joint Enrichment Plot for TF Enrichment HFD / DDC / T0_Organoids
# #####################################################################################################################################
# Load analysis results
DDC_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/DDC_Enrichment_Overlap_Analysis.RDS")
OrgT0_vs_InVivo_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/OrgT0_vs_InVivo_Enrichment_Overlap_Analysis.RDS")

DDC.GO.Up <- DDC_Enrichment_Overlap_Analysis$DDC_up$TF
DESeq2.GO.Up <- DDC_Enrichment_Overlap_Analysis$DESeq2_up$TF
OrgT0_vs_InVivo <- OrgT0_vs_InVivo_Enrichment_Overlap_Analysis$OrgT0.vs.InVivo_up$TF

DDC.GO.Up$comp <- "DDC_Up"
DESeq2.GO.Up$comp <- "DESeq2_Up"
OrgT0_vs_InVivo$comp <- "OrgT0_vs_InVivo"

df_merged <- rbind(DDC.GO.Up[1:10,], DESeq2.GO.Up[1:10,], OrgT0_vs_InVivo[1:10,])
df_merged <- df_merged[!is.na(df_merged$Count),]
df_merged$signficant <- df_merged$p.adjust <= 0.05
df_merged$Description <- tolower(gsub("TFACTS_MM_|TF_MM_","", df_merged$Description))
df_merged$Gene_Ratio <- sapply(df_merged$GeneRatio, function(x) round(eval(parse(text=x)), digits = 2))

# plot only significant ones
df_merged <- df_merged[df_merged$signficant,]

g <- ggplot(df_merged, aes(x = Gene_Ratio, y = forcats::fct_reorder(Description, Gene_Ratio, .desc = FALSE))) + 
  geom_point(aes(size = Count, col =  p.adjust )) +
  scale_colour_gradient("p.adjust", limits=c(0, 0.05), low="red", high = "blue") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("Gene Ratio")) +
  facet_wrap(~comp) +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_TF_ALL_Merged.pdf"), 
       height = 4, width = 8)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_TF_ALL_Merged.svg"), 
       height = 4, width = 8)


# #####################################################################################################################################
# 10. Joint Enrichment Plot for GO Enrichment HFD / DDC / T0_Organoids
# #####################################################################################################################################
# Load analysis results
DDC_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/DDC_Enrichment_Overlap_Analysis.RDS")
OrgT0_vs_InVivo_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/OrgT0_vs_InVivo_Enrichment_Overlap_Analysis.RDS")

DDC.GO.Down <- DDC_Enrichment_Overlap_Analysis$DDC_down$GO
DESeq2.GO.Down <- DDC_Enrichment_Overlap_Analysis$DESeq2_down$GO
OrgT0_vs_InVivo.Down <- OrgT0_vs_InVivo_Enrichment_Overlap_Analysis$OrgT0.vs.InVivo_down$GO

DDC.GO.Down$comp <- "DDC_Down"
DESeq2.GO.Down$comp <- "DESeq2_Down"
OrgT0_vs_InVivo.Down$comp <- "OrgT0_vs_InVivo_Down"

df_merged <- rbind(DDC.GO.Down[1:10,], DESeq2.GO.Down[1:10,], OrgT0_vs_InVivo.Down[1:10,])
df_merged <- df_merged[!is.na(df_merged$Count),]
df_merged$signficant <- df_merged$p.adjust <= 0.05
df_merged$Description <- tolower(gsub("GO_.._MM_","", df_merged$Description))
df_merged$Gene_Ratio <- sapply(df_merged$GeneRatio, function(x) round(eval(parse(text=x)), digits = 2))

# plot only significant ones
df_merged <- df_merged[df_merged$signficant,]

g <- ggplot(df_merged, aes(x = Gene_Ratio, y = forcats::fct_reorder(Description, Gene_Ratio, .desc = FALSE))) + 
  geom_point(aes(size = Count, col =  p.adjust )) +
  scale_colour_gradient("p.adjust", limits=c(0, 0.05), low="red", high = "blue") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("Gene Ratio")) +
  facet_wrap(~comp) +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_GO_ALL_Merged.pdf"), 
       height = 4, width = 8)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_GO_ALL_Merged.svg"), 
       height = 4, width = 8)

# ----------------------------------
# Alternative only Chosen - Keywords
# ----------------------------------
# Load analysis results
DDC_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/DDC_Enrichment_Overlap_Analysis.RDS")
OrgT0_vs_InVivo_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/OrgT0_vs_InVivo_Enrichment_Overlap_Analysis.RDS")

DDC.GO.Down <- DDC_Enrichment_Overlap_Analysis$DDC_down$GO
DESeq2.GO.Down <- DDC_Enrichment_Overlap_Analysis$DESeq2_down$GO
OrgT0_vs_InVivo.Down <- OrgT0_vs_InVivo_Enrichment_Overlap_Analysis$OrgT0.vs.InVivo_down$GO

DDC.GO.Down$comp <- "DDC_Down"
DESeq2.GO.Down$comp <- "DESeq2_Down"
OrgT0_vs_InVivo.Down$comp <- "OrgT0_vs_InVivo_Down"

pattern <- "adhesion|extracellular_matrix|collagen|cell_surface"
DDC.GO.Down <- DDC.GO.Down[grepl(pattern, tolower(DDC.GO.Down$Description)),]
DESeq2.GO.Down <- DESeq2.GO.Down[grepl(pattern, tolower(DESeq2.GO.Down$Description)),]
OrgT0_vs_InVivo.Down <- OrgT0_vs_InVivo.Down[grepl(pattern, tolower(OrgT0_vs_InVivo.Down$Description)),]

df_merged <- rbind(DDC.GO.Down, DESeq2.GO.Down, OrgT0_vs_InVivo.Down)
df_merged <- df_merged[!is.na(df_merged$Count),]
df_merged$signficant <- df_merged$p.adjust <= 0.05
df_merged$Description <- tolower(gsub("GO_.._MM_","", df_merged$Description))
df_merged$Gene_Ratio <- sapply(df_merged$GeneRatio, function(x) round(eval(parse(text=x)), digits = 2))

# plot only significant ones
df_merged <- df_merged[df_merged$signficant,]

g <- ggplot(df_merged, aes(x = Gene_Ratio, y = forcats::fct_reorder(Description, Gene_Ratio, .desc = FALSE))) + 
  geom_point(aes(size = Count, col =  p.adjust )) +
  scale_colour_gradient("p.adjust", limits=c(0, 0.05), low="red", high = "blue") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("Gene Ratio")) +
  facet_wrap(~comp) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_GO_CHOSEN_ALL_Merged.pdf"), 
       height = 4, width = 8)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_GO_CHOSEN_ALL_Merged.svg"), 
       height = 4, width = 8)


# ####################################################################################################################################
# 11. Gene co-expression Network Analysis
# ####################################################################################################################################
n.cores <- 38                      # number of cores/threads to call for PCP
doPar <- TRUE                      # do we want to parallelize?
method <- "pearson"                # method for correlation. either pearson or spearman. 
FDR.cutoff <- 0.05                 # FDR threshold to define significant correlations upon shuffling samples. 
module.pval <- 0.05                # module significance p-value. Recommended is 0.05. 
hub.pval <- 0.05                   # connectivity significance p-value based random tetrahedral networks
cor.perm <- 10                     # number of permutations for calculating FDRs for all correlation pairs. 
hub.perm <- 100                    # number of permutations for calculating connectivity significance p-value. 
# annotation to be done on the downstream
annot.table <- NULL
id.col <- 1
symbol.col <- 2

# Load Expression Data
load("./Data/R_Data/Exp.RData", verbose = T)
# CPM & Log transform counts
data.Expr <- log2(edgeR::cpm(exp) + 1)
# https://github.com/songw01/MEGENA/issues/5
# You should run on all genes that have substantial variances across the samples 
# (e.g. coefficient of variation > 0.1 is a widely used threshold).
data.Expr <- data.Expr[matrixStats::rowVars(as.matrix(data.Expr)) > 0.2,]

# Calculate Pairwise Correlation
ijw <- MEGENA::calculate.correlation(data.Expr, 
                                     doPerm = cor.perm, 
                                     output.corTable = FALSE, 
                                     output.permFDR = FALSE,
                                     num.cores = n.cores)

# calculate PFN
el <- MEGENA::calculate.PFN(ijw[,1:3], 
                            doPar = doPar, 
                            num.cores = n.cores, 
                            keep.track = FALSE)

g <- graph.data.frame(el, directed = FALSE)

# perform MCA clustering.
MEGENA.output <- MEGENA::do.MEGENA(g, 
                                   mod.pval = module.pval,
                                   hub.pval = hub.pval,
                                   remove.unsig = TRUE,
                                   min.size = 10,
                                   max.size = vcount(g)/2,
                                   doPar = doPar,
                                   num.cores = n.cores,
                                   n.perm = hub.perm,
                                   save.output = FALSE)

summary.output <- MEGENA::MEGENA.ModuleSummary(MEGENA.output,
                                               mod.pvalue = module.pval,
                                               hub.pvalue = hub.pval,
                                               min.size = 10,
                                               max.size = vcount(g)/2,
                                               annot.table = annot.table,
                                               id.col = id.col,
                                               symbol.col = symbol.col,
                                               output.sig = TRUE)
# Examine Modules Obtained
print(summary.output$module.table)

pnet.obj <- MEGENA::plot_module(output.summary = summary.output,
                                PFN = g,
                                subset.module = "c1_2",
                                layout = "kamada.kawai",
                                label.hubs.only = TRUE,
                                gene.set = NULL,
                                color.code =  "grey",
                                output.plot = FALSE,
                                out.dir = "modulePlot",
                                col.names = c("magenta","green","cyan"),
                                label.scaleFactor = 20,
                                hubLabel.col = "black",
                                hubLabel.sizeProp = 1,
                                show.topn.hubs = Inf,
                                show.legend = TRUE)
# Generate Default Plots
pdf(paste0("./Plots/Co_Expression_Analysis_MEGENA/MEGENEA_Hub_Graph_Plot.pdf"), 
    height = 20, width = 20, useDingbats=FALSE)
print(pnet.obj[[1]])
dev.off()

module.table <- summary.output$module.table
# first column of module table must be labelled as "id"
colnames(module.table)[1] <- "id" 
hierarchy.obj <- MEGENA::plot_module_hierarchy(module.table = module.table,
                                               label.scaleFactor = 0.15,
                                               arrow.size = 0.03,
                                               node.label.color = "blue")
pdf(paste0("./Plots/Co_Expression_Analysis_MEGENA/MEGENEA_Hierarchy_Graph_Plot.pdf"), 
    height = 20, width = 20, useDingbats=FALSE)
print(hierarchy.obj[[1]])
dev.off()

MEGENA.result.table <-
  MEGENA::module_convert_to_table(MEGENA.output,
                                  mod.pval = 0.05,
                                  hub.pval = 0.05,
                                  min.size = 50,
                                  max.size = 2000)

# Save gmt object to use for downstream analyses
save(MEGENA.result.table, file = "./Data/R_Data/MEGENA_Transcriptomics_GMT.RData")

# ------------------------
# Transcriptomic Hub Genes
# ------------------------
load("./Data/R_Data/MEGENA_Transcriptomics_GMT.RData", verbose = T)

# Subset Hub genes from table and save to reports for downstream usage
Hub.genes <- MEGENA.result.table[!is.na(MEGENA.result.table$is.hub),]
Hub.genes$gene_name <- DESeq2.Results.Table[match(Hub.genes$id, DESeq2.Results.Table$gene_id), 'gene_name']

write.table(Hub.genes, 
            file = "./Reports/Co_Expression_Analysis_MEGENA/Transcriptomics_Hub_Genes.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# ----------------------
# Enrichment of Clusters
# ----------------------
load("./Data/R_Data/MEGENA_Transcriptomics_GMT.RData", verbose = T)
load("./Data/R_Data/Exp.RData", verbose = T)

genesets <- list()
genesets$GO <- tail(read.gmt("./Resources/GeneSets/MousePath_GO_gmt.gmt"), -1)
genesets$mGSKB <- read.gmt("./Resources/GeneSets/mGSKB_Ensembl.gmt")
genesets$Pathway <- tail(read.gmt("./Resources/GeneSets/MousePath_Pathway_gmt.gmt"), -1)
genesets$Metabolic <- tail(read.gmt("./Resources/GeneSets/MousePath_Metabolic_gmt.gmt"), -1)
genesets$TF <- tail(read.gmt("./Resources/GeneSets/MousePath_TF_gmt.gmt"), -1)

# Prepare gene groups
module.genes <- mapply(module = unique(MEGENA.result.table$module), function(module){
  MEGENA.result.table[MEGENA.result.table$module == module, "id"]
}, SIMPLIFY = FALSE)

gmt_GO <- tail(read.gmt("./Resources/GeneSets/MousePath_GO_gmt.gmt"), -1)
gmt_mGSKB <- read.gmt("./Resources/GeneSets/mGSKB_Ensembl.gmt")
gmt_Pathway <- tail(read.gmt("./Resources/GeneSets/MousePath_Pathway_gmt.gmt"), -1)
gmt_Metabolic <- tail(read.gmt("./Resources/GeneSets/MousePath_Metabolic_gmt.gmt"), -1)
gmt_TF <- tail(read.gmt("./Resources/GeneSets/MousePath_TF_gmt.gmt"), -1)

gmt_list <- list(gmt_GO, gmt_mGSKB, gmt_Pathway, gmt_Metabolic, gmt_TF)
gmt_list_names <- c("GO", "mGSKB", "Pathway", "Metabolic", "TF")

DE_list <- module.genes
DE_data_list <- names(module.genes)

En.list <-
  mapply(condition = DE_data_list, function(condition){
    # Specify Paths
    Plots_path = paste0('./Plots/Co_Expression_Analysis_MEGENA/Enrichment/', condition)
    Reports_path = paste0("./Reports/Co_Expression_Analysis_MEGENA/Enrichment/", condition)
    
    if(!dir.exists(Plots_path)){
      dir.create(Plots_path)
    }
    
    if(!dir.exists(Reports_path)){
      dir.create(Reports_path)
    }
    
    # Create Excel object
    wb <- openxlsx::createWorkbook()
    # Get gene ID
    gene.ID <- DE_list[[condition]]
    gene.Name <- DESeq2.Results.Table[match(gene.ID, DESeq2.Results.Table$gene_id), 'gene_name']
    
    en.tmp <-
      mapply(gmt = gmt_list_names, function(gmt){
        if(gmt == "mGSKB"){
          en <- enricher(gene = gene.ID,
                         TERM2GENE = gmt_list[[match(gmt, gmt_list_names)]])
          
          write.csv(data.frame(en@result), file = paste0(Reports_path, "/Enrichement_", condition, "_", gmt, ".csv"))
        } else {
          en <- enricher(gene = toupper(gene.Name),
                         TERM2GENE = gmt_list[[match(gmt, gmt_list_names)]])
          
          # export GE results
          openxlsx::addWorksheet(wb, gmt)
          openxlsx::writeData(wb, sheet = gmt, data.frame(en@result))
        }
        
        if(sum(en@result$p.adjust < 0.05) > 0){
          # Generate default Dot Plot
          en.dotplot <- dotplot(en, showCategory = 30)
          ggsave(plot = en.dotplot, 
                 paste0(Plots_path, "/", condition, "_GE_DotPlot_", gmt, ".pdf"), 
                 height = 20, width = 20)
          
          # Generate Upset Plot
          en.upsetplot <- enrichplot::upsetplot(en)
          ggsave(plot = en.upsetplot, 
                 file = paste0(Plots_path, "/", condition, "_GE_UpsetPlot_", gmt, ".pdf"),
                 width=20, height=8)
        }
        
        return(en@result[en@result$p.adjust < 0.05,])
      }, SIMPLIFY = FALSE)
    
    openxlsx::saveWorkbook(wb,
                           paste0(Reports_path, "/Enrichement_", 
                                  condition, "_",".xlsx"), 
                           overwrite = TRUE)
    return(en.tmp)
  }, SIMPLIFY = FALSE)

saveRDS(En.list, "./Data/R_Data/MEGENA_En.RDS")


# ####################################################################################################################################
# 12. Alternative Organoid Enrichment Plots
# ####################################################################################################################################
# Load analysis results
OrgT0_vs_InVivo_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/OrgT0_vs_InVivo_Enrichment_Overlap_Analysis.RDS")

OrgT0_vs_InVivo <- OrgT0_vs_InVivo_Enrichment_Overlap_Analysis$OrgT0.vs.InVivo_up$TF

df_merged <- OrgT0_vs_InVivo[1:30,]
df_merged <- df_merged[!is.na(df_merged$Count),]
df_merged$signficant <- df_merged$p.adjust <= 0.05
df_merged$Description <- tolower(gsub(".*_MM_","", df_merged$Description))
df_merged$Gene_Ratio <- sapply(df_merged$GeneRatio, function(x) round(eval(parse(text=x)), digits = 2))

# plot only significant ones
df_merged <- df_merged[df_merged$signficant,]

g <- ggplot(df_merged, aes(x = Gene_Ratio, y = forcats::fct_reorder(Description, Gene_Ratio, .desc = FALSE))) + 
  geom_point(aes(size = Count, col =  p.adjust )) +
  scale_colour_gradient("p.adjust", limits=c(0, 0.05), low="red", high = "blue") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("Gene Ratio")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )
ggsave(plot = g, paste0("./Plots/Alternative_Figures/OrgT0_vs_InVivo_Enrichment_UP_TF.pdf"), 
       height = 5, width = 4)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/OrgT0_vs_InVivo_Enrichment_UP_TF.svg"), 
       height = 5, width = 4)


# #####################################################################################################################################
# 13. GSEA KEGG Plot Top / Down 15 Changing Gene Sets
# #####################################################################################################################################
# Load GSEA results to plot
load("./Data/R_Data/Exp_DESeq2_GSEA.RData", verbose = T)
# Proper GSEA Plot
result <- DE.gsea$KEGG@result
result$core_enrichmentNgenes <- sapply(strsplit(unlist(result$core_enrichment), 
                                                split = "/", fixed = T), length)
result$core_enrichmentgenes <- unlist(sapply(strsplit(unlist(result$core_enrichment), 
                                                      split = "/", fixed = T), 
                                             function(x){paste(sort(x), collapse = " ")}))
result$signficant <- result$qvalues <= 0.01
result$label <- result$Description
result$gene_ratio <- result$core_enrichmentNgenes / result$setSize

# sort by abs value of NES and take top 15 values
tmp_data <- head(result[order(result$qvalues),], 15)

GSEA.plot <-
  ggplot(tmp_data, aes(x = NES, y = forcats::fct_reorder(label, NES, .desc = FALSE))) + 
  geom_point(aes(size = abs(gene_ratio), col =  qvalues )) +
  scale_colour_gradient("q-value", limits=c(0, 0.2), low="chartreuse3", high = "brown1") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("NES")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )

ggsave(plot = GSEA.plot, 
       paste0("Plots/Alternative_Figures/", "GSEA_Top_15_qValue_KEGG.pdf"), 
       height = 6, width = 7)


# #####################################################################################################################################
# 14. Joint Enrichment Plot for GO Enrichment HFD / DDC
# #####################################################################################################################################
# Load analysis results
DDC_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/DDC_Enrichment_Overlap_Analysis.RDS")

DDC.GO.Up <- DDC_Enrichment_Overlap_Analysis$DDC_up$GO
DESeq2.GO.Up <- DDC_Enrichment_Overlap_Analysis$DESeq2_up$GO

DDC.GO.Up$comp <- "DDC_Up"
DESeq2.GO.Up$comp <- "DESeq2_Up"

df_merged <- rbind(DDC.GO.Up[1:10,], DESeq2.GO.Up[1:10,])
df_merged$signficant <- df_merged$p.adjust <= 0.05
df_merged$Description <- tolower(gsub("GO_.._MM_","", df_merged$Description))
df_merged$Gene_Ratio <- sapply(df_merged$GeneRatio, function(x) round(eval(parse(text=x)), digits = 2))

# plot only significant ones
df_merged <- df_merged[df_merged$signficant,]

g <- ggplot(df_merged, aes(x = Gene_Ratio, y = forcats::fct_reorder(Description, Gene_Ratio, .desc = FALSE))) + 
  geom_point(aes(size = Count, col =  p.adjust )) +
  scale_colour_gradient("p.adjust", limits=c(0, 0.05), low="red", high = "blue") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("Gene Ratio")) +
  facet_wrap(~comp) +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_DDC_GO_DESeq2_Up_Merged_Top_10.pdf"), 
       height = 4, width = 8)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_DDC_GO_DESeq2_Up_Merged_Top_10.svg"), 
       height = 4, width = 8)

# DOWN
# Load analysis results
DDC_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/DDC_Enrichment_Overlap_Analysis.RDS")

DDC.GO.Down <- DDC_Enrichment_Overlap_Analysis$DDC_down$GO
DESeq2.GO.Down <- DDC_Enrichment_Overlap_Analysis$DESeq2_down$GO

DDC.GO.Down$comp <- "DDC_Down"
DESeq2.GO.Down$comp <- "DESeq2_Down"

df_merged <- rbind(DDC.GO.Down[1:10,], DESeq2.GO.Down[1:10,])
df_merged$signficant <- df_merged$p.adjust <= 0.05
df_merged$Description <- tolower(gsub("GO_.._MM_","", df_merged$Description))
df_merged$Gene_Ratio <- sapply(df_merged$GeneRatio, function(x) round(eval(parse(text=x)), digits = 2))

# plot only significant ones
df_merged <- df_merged[df_merged$signficant,]

g <- ggplot(df_merged, aes(x = Gene_Ratio, y = forcats::fct_reorder(Description, Gene_Ratio, .desc = FALSE))) + 
  geom_point(aes(size = Count, col =  p.adjust )) +
  scale_colour_gradient("p.adjust", limits=c(0, 0.05), low="red", high = "blue") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("Gene Ratio")) +
  facet_wrap(~comp) +
  
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_DDC_GO_DESeq2_Down_Merged_Top_10.pdf"), 
       height = 4, width = 8)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/Enrichment_Overlap_Analysis_DDC_GO_DESeq2_Down_Merged_Top_10.svg"), 
       height = 4, width = 8)


# #####################################################################################################################################
# 15. Boxplot of tet1 gene in CDvsHFD
# #####################################################################################################################################
# Load the expression data
load("./Data/R_Data/Exp.RData", verbose = T)

# CPM & Log transform counts
data.Expr <- log2(edgeR::cpm(exp) + 1)
# Get only expression of the gene of interest
data.Expr <- data.Expr['ENSMUSG00000047146',]
# Prepare dataframe for plotting
Expr.df.ncam1 <-
  data.frame(sample = names(data.Expr),
             `log(cpm + 1)` = data.Expr,
             condition = gsub('(.*)_.*', '\\1', names(data.Expr)))
g <- 
  ggplot(Expr.df.ncam1, aes(x=condition, y=`log.cpm...1.`, fill = condition)) + 
  geom_boxplot()+
  geom_dotplot(binaxis='y', stackdir='center')+
  labs(y="Log2(cpm + 1)", x = "Condition")+
  theme_classic() +
  ggpubr::stat_compare_means(method = "t.test")

ggsave(plot = g, paste0("./Plots/Alternative_Figures/Tet1_Boxplot.pdf"), 
       height = 4, width = 4)


# #####################################################################################################################################
# 16. canonical glycolysis      GO:0061621
# #####################################################################################################################################
# Load the expression data
OrgT0.vs.InVivo.DE <- readxl::read_xlsx("./Resources/Organoid_Paper/41556_2019_402_MOESM4_ESM.xlsx", 
                                        sheet = "8. T0 vs organoids", skip = 1)
colnames(OrgT0.vs.InVivo.DE) <- c("gene_name", "gene_id", colnames(OrgT0.vs.InVivo.DE)[3:36])
# Chosen Lipid genes to plot
Genes <- tolower(unique(readxl::read_xlsx("./Resources/GeneSets/GO_term_summary_20220519_070539.xlsx")$Symbol))

# Create dataframe with expression of genes
gene.list <- unique(Genes)
gene.list <- gene.list[gene.list %in% tolower(OrgT0.vs.InVivo.DE$gene_name)]
gene.list.exp <- OrgT0.vs.InVivo.DE[match(gene.list, tolower(OrgT0.vs.InVivo.DE$gene_name)),]
df <- 
  data.frame(gene = gene.list.exp$gene_name,
             Genes = gene.list.exp$b,
             p.value = gene.list.exp$qval,
             y = "Genes")
df$label <- ""
df$label[df$p.value < 0.05] <- "*"
df$Genes <- -1 * as.numeric(df$Genes)
# Order by decreasing value
df <- df[with(df, order(-Genes)),]
df$gene <- factor(df$gene, levels = df$gene)

g <- ggplot(df, aes(y, gene, fill= Genes)) + 
  geom_tile() +
  geom_text(aes(label= label),
            position = position_dodge(width=0),  size=5) +
  scale_fill_gradient2(low = "blue", mid = "white",high = "red") +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 12, angle = 90),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.spacing.x=unit(0, "lines"),
        panel.spacing.y=unit(0, "lines"),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 5),
        strip.text.x = element_text(size = 7, 
                                    face = "bold", color = "black"),
        legend.position="left") +
  coord_fixed(ratio = 1/4) + 
  xlab("") +
  scale_y_discrete(position = "right", expand = expand_scale(mult = c(0, 0))) +
  scale_x_discrete(expand = expand_scale(mult = c(0, 0)))

pdf("./Plots/Alternative_Figures/Selected_Lipid_Genes_Heatmap.pdf", 
    height = 3, width = 2.5, useDingbats=FALSE)
g
dev.off()


# #####################################################################################################################################
# 17. aerobic electron transport chain     GO:0019646
# #####################################################################################################################################
# Load the expression data
OrgT0.vs.InVivo.DE <- readxl::read_xlsx("./Resources/Organoid_Paper/41556_2019_402_MOESM4_ESM.xlsx", 
                                        sheet = "8. T0 vs organoids", skip = 1)
colnames(OrgT0.vs.InVivo.DE) <- c("gene_name", "gene_id", colnames(OrgT0.vs.InVivo.DE)[3:36])
# Chosen Lipid genes to plot
Genes <- tolower(unique(readxl::read_xlsx("./Resources/GeneSets/GO_term_summary_20220602_073529.xlsx")$Symbol))

# Create dataframe with expression of genes
gene.list <- unique(Genes)
gene.list <- gene.list[gene.list %in% tolower(OrgT0.vs.InVivo.DE$gene_name)]
gene.list.exp <- OrgT0.vs.InVivo.DE[match(gene.list, tolower(OrgT0.vs.InVivo.DE$gene_name)),]
df <- 
  data.frame(gene = gene.list.exp$gene_name,
             Genes = gene.list.exp$b,
             p.value = gene.list.exp$qval,
             y = "Genes")
df$label <- ""
df$label[df$p.value < 0.05] <- "*"
df$Genes <- -1 * as.numeric(df$Genes)
# Order by decreasing value
df <- df[with(df, order(-Genes)),]
df$gene <- factor(df$gene, levels = df$gene)

g <- ggplot(df, aes(y, gene, fill= Genes)) + 
  geom_tile() +
  geom_text(aes(label= label),
            position = position_dodge(width=0),  size=5) +
  scale_fill_gradient2(low = "blue", mid = "white",high = "red") +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 12, angle = 90),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.spacing.x=unit(0, "lines"),
        panel.spacing.y=unit(0, "lines"),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 5),
        strip.text.x = element_text(size = 7, 
                                    face = "bold", color = "black"),
        legend.position="left") +
  coord_fixed(ratio = 1/4) + 
  xlab("") +
  scale_y_discrete(position = "right", expand = expand_scale(mult = c(0, 0))) +
  scale_x_discrete(expand = expand_scale(mult = c(0, 0)))

pdf("./Plots/Alternative_Figures/GO_0019646_Heatmap.pdf", 
    height = 5, width = 3, useDingbats=FALSE)
g
dev.off()


# #####################################################################################################################################
# 18. regulation of oxidative phosphorylation     GO:0002082
# #####################################################################################################################################
# Load the expression data
OrgT0.vs.InVivo.DE <- readxl::read_xlsx("./Resources/Organoid_Paper/41556_2019_402_MOESM4_ESM.xlsx", 
                                        sheet = "8. T0 vs organoids", skip = 1)
colnames(OrgT0.vs.InVivo.DE) <- c("gene_name", "gene_id", colnames(OrgT0.vs.InVivo.DE)[3:36])
# Chosen Lipid genes to plot
Genes <- tolower(unique(readxl::read_xlsx("./Resources/GeneSets/GO_term_summary_20220602_073546.xlsx")$Symbol))

# Create dataframe with expression of genes
gene.list <- unique(Genes)
gene.list <- gene.list[gene.list %in% tolower(OrgT0.vs.InVivo.DE$gene_name)]
gene.list.exp <- OrgT0.vs.InVivo.DE[match(gene.list, tolower(OrgT0.vs.InVivo.DE$gene_name)),]
df <- 
  data.frame(gene = gene.list.exp$gene_name,
             Genes = gene.list.exp$b,
             p.value = gene.list.exp$qval,
             y = "Genes")
df$label <- ""
df$label[df$p.value < 0.05] <- "*"
df$Genes <- -1 * as.numeric(df$Genes)
# Order by decreasing value
df <- df[with(df, order(-Genes)),]
df$gene <- factor(df$gene, levels = df$gene)

g <- ggplot(df, aes(y, gene, fill= Genes)) + 
  geom_tile() +
  geom_text(aes(label= label),
            position = position_dodge(width=0),  size=5) +
  scale_fill_gradient2(low = "blue", mid = "white",high = "red") +
  theme(axis.text.x = element_text(face = "bold", color = "black", 
                                   size = 12, angle = 90),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.spacing.x=unit(0, "lines"),
        panel.spacing.y=unit(0, "lines"),
        axis.text.y = element_text(face = "bold", color = "black", 
                                   size = 5),
        strip.text.x = element_text(size = 7, 
                                    face = "bold", color = "black"),
        legend.position="left") +
  coord_fixed(ratio = 1/4) + 
  xlab("") +
  scale_y_discrete(position = "right", expand = expand_scale(mult = c(0, 0))) +
  scale_x_discrete(expand = expand_scale(mult = c(0, 0)))

pdf("./Plots/Alternative_Figures/GO_0002082_Heatmap.pdf", 
    height = 3, width = 2.5, useDingbats=FALSE)
g
dev.off()


# #####################################################################################################################################
# 19. HFD DE Volcano Plot
# #####################################################################################################################################
load("./Data/R_Data/Exp_DESeq2.RData", verbose = T)

allTopTables.df <- DESeq2.Results.Table
cont <- unique(allTopTables.df$Contrast)

saveDir <- paste0("./Plots/DESeq2/Volcano/")
if(!dir.exists(saveDir)){
  dir.create(saveDir, recursive = T)
}

dfPlot                  <- allTopTables.df[allTopTables.df$Contrast == cont, ]
dfPlot                  <- dfPlot[!is.na(dfPlot$log2FoldChange) & !is.na(dfPlot$padj), ]
dfPlot$color            <- "C1"
dfPlot$color[dfPlot$padj < 0.05] <- "C2"
dfPlot                  <- dfPlot[order(dfPlot$padj, decreasing = F), ]
dfPlot$labels           <- ""

sign <- abs(dfPlot$log2FoldChange) > 1 & dfPlot$padj < 0.05
dfPlot$labels[sign] <- dfPlot$gene_name[sign]
dfPlot <- dfPlot[with(dfPlot, order(padj)),]
dfPlot$labels[25:dim(dfPlot)[1]] <- ""

tmpPos1 <- which(dfPlot$log2FoldChange > 0 & dfPlot$labels != "")
tmpPos2 <- which(dfPlot$log2FoldChange < 0 & dfPlot$labels != "")

xLims_pos <- base::range(na.omit(dfPlot$log2FoldChange[dfPlot$log2FoldChange > 0]))
xLims_neg <- base::range(na.omit(dfPlot$log2FoldChange[dfPlot$log2FoldChange < 0]))

pl.volc <- ggplot(dfPlot, aes(x = log2FoldChange, y = -log10(padj), label = labels)) +
  geom_point(data = dfPlot[dfPlot$padj < 0.05 & dfPlot$log2FoldChange > 1, ], aes(id = padj), 
             size = 1.5, alpha = 1, color = '#ee442f') +
  geom_point(data = dfPlot[dfPlot$padj < 0.05 & dfPlot$log2FoldChange < -1, ], aes(id = padj), 
             size = 1.5, alpha = 1, color = '#63acbe') +
  geom_point(data = dfPlot[dfPlot$padj > 0.05, ], aes(id = padj), size = 1.5, alpha = 1, color = 'grey') +
  geom_point(data = dfPlot[dfPlot$padj < 0.05 & abs(dfPlot$log2FoldChange) < 1, ], aes(id = padj), 
             size = 1.5, alpha = 1, color = 'grey') +
  ggrepel::geom_text_repel(data               = dfPlot[tmpPos1, ],
                           nudge_x            = dfPlot[tmpPos1, ]$log2FoldChange + xLims_pos[2] * 0.1,
                           segment.size       = 0.2,
                           segment.color      = "grey50",
                           direction          = "y",
                           hjust              = 0,
                           size               = 1.5,
                           force              = 50,
                           min.segment.length = 0) +
  ggrepel::geom_text_repel(data               = dfPlot[tmpPos2, ],
                           nudge_x            = dfPlot[tmpPos2, ]$log2FoldChange + xLims_neg[1] * 0.1,
                           segment.size       = 0.2,
                           segment.color      = "grey50",
                           direction          = "y",
                           hjust              = 0,
                           size               = 1.5,
                           force              = 50,
                           min.segment.length = 0) +
  scale_color_manual(values = c("C1" = "#bebebe", "C2" = "#b0222b", "C3" = "#6799e4")) +
  scale_x_continuous(expand = expansion(mult = c(0.22, 0.22))) +
  scale_y_continuous(expand = expansion(mult = c(0.025, 0.05))) +
  xlab("log2 fold change") + 
  ylab(paste0("-log10(adj. pValue)")) +
  theme_classic() +
  theme(
    plot.title      = element_text(hjust = 0.5, size = 10)) +
  guides(color = FALSE)
ggsave(paste0(saveDir, "/volcanoPlot_contrast_withLabels_3x3.pdf"), 
       plot = pl.volc, width = 4, height = 4, useDingbats = F)


# #####################################################################################################################################
# 20. Custom PATHWAY / METABOLIC Plots
# #####################################################################################################################################
# Load analysis results
OrgT0_vs_InVivo_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/OrgT0_vs_InVivo_Enrichment_Overlap_Analysis.RDS")

OrgT0_vs_InVivo <- OrgT0_vs_InVivo_Enrichment_Overlap_Analysis$OrgT0.vs.InVivo_up$Pathway

df <- OrgT0_vs_InVivo[1:25,]
df <- df[!is.na(df$Count),]
df$signficant <- df$p.adjust <= 0.05
df$Description <- gsub("_", " ", tolower(gsub(".*_MM_","", df$Description)))
df$Gene_Ratio <- sapply(df$GeneRatio, function(x) round(eval(parse(text=x)), digits = 2))

# plot only significant ones
df <- df[df$signficant,]

g <- ggplot(df, aes(x = Gene_Ratio, y = forcats::fct_reorder(Description, Gene_Ratio, .desc = FALSE))) + 
  geom_point(aes(size = Count, col =  p.adjust )) +
  scale_colour_gradient("p.adjust", limits=c(0, 0.05), low="red", high = "blue") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("Gene Ratio")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )
ggsave(plot = g, paste0("./Plots/Alternative_Figures/OrgT0_vs_InVivo_Enrichment_UP_PATHWAY.pdf"), 
       height = 5, width = 7.5)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/OrgT0_vs_InVivo_Enrichment_UP_PATHWAY.svg"), 
       height = 5, width = 7.5)

# Load analysis results
OrgT0_vs_InVivo_Enrichment_Overlap_Analysis <- readRDS("./Data/R_Data/OrgT0_vs_InVivo_Enrichment_Overlap_Analysis.RDS")

OrgT0_vs_InVivo <- OrgT0_vs_InVivo_Enrichment_Overlap_Analysis$OrgT0.vs.InVivo_up$Metabolic
OrgT0_vs_InVivo <- OrgT0_vs_InVivo[grepl("^EHMN_MM|^KEGG_MM", OrgT0_vs_InVivo$ID),]

df <- OrgT0_vs_InVivo[1:25,]
df <- df[!is.na(df$Count),]
df$signficant <- df$p.adjust <= 0.05
df$Description <- gsub("_", " ", tolower(gsub(".*_MM_","", df$Description)))
df$Gene_Ratio <- sapply(df$GeneRatio, function(x) round(eval(parse(text=x)), digits = 2))

# plot only significant ones
df <- df[df$signficant,]

g <- ggplot(df, aes(x = Gene_Ratio, y = forcats::fct_reorder(Description, Gene_Ratio, .desc = FALSE))) + 
  geom_point(aes(size = Count, col =  p.adjust )) +
  scale_colour_gradient("p.adjust", limits=c(0, 0.05), low="red", high = "blue") +
  theme_bw(base_size = 14) +
  ylab(NULL) +
  xlab(paste0("Gene Ratio")) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
    axis.title = element_text(size = 19),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 13)
  )
ggsave(plot = g, paste0("./Plots/Alternative_Figures/OrgT0_vs_InVivo_Enrichment_UP_Metabolic.pdf"), 
       height = 5, width = 7.5)
ggsave(plot = g, paste0("./Plots/Alternative_Figures/OrgT0_vs_InVivo_Enrichment_UP_Metabolic.svg"), 
       height = 5, width = 7.5)


# #####################################################################################################################################
# 21. GEO Submission Metadata
# #####################################################################################################################################
load("./Data/R_Data/metadata.RData", verbose = T)

geo.meta <- 
  data.frame("Sample name" = metadata$Raw.file,
           "Title" = metadata$Raw.file,
           "source name" = metadata$Sample.type,
           "organism" = "Mus musculus",
           "Strain" = metadata$Mice.strain,
           "Age" = "8",
           "Treatment" = metadata$condition,
           "molecule" = "total RNA",
           "description" = metadata$SampleID,
           "processed" = paste0("refAlign.", metadata$Raw.file, ".ReadsPerGene.out.tab"),
           "raw1" = paste0(metadata$Raw.file, "_1.fq.gz"),
           "raw2" = paste0(metadata$Raw.file, "_2.fq.gz"))

# Generate md5sum for read count files
geo.counts.md5 <- data.frame("sample" = c(),
                             "md5sum" = c())
for(sample in paste0("refAlign.", metadata$Raw.file, ".ReadsPerGene.out.tab")){
  df_tmp <- data.frame(
    "sample" = sample,
    "type" = "raw read count",
    "md5sum" = gsub("(.*)\\s.*", "\\1", system(paste0("md5sum ", "./Data/Star_alignments/", sample), intern = TRUE)))
  geo.counts.md5 <- rbind(geo.counts.md5, df_tmp)
}

# Generate md5sum for raw files
geo.raw.md5 <- read.table("./Data/RAW/F20FTSEUHT0958_MOUipaE_12sample_lowyield_ece/Filter_SOAPnuke/md5sum_check.txt")
geo.raw.md5 <- data.frame(
  "sample" = gsub("Clean/.*/(.*.fq.gz$)", "\\1", geo.raw.md5$V2),
  "type" = "fq.gz",
  "md5" = geo.raw.md5$V1,
  "instr" = "BGISEQ-500 Transcriptome",
  "sing" = "paired-end")


# Create Excel object
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "geo_meta")
openxlsx::writeData(wb, sheet = "geo_meta", geo.meta)
openxlsx::addWorksheet(wb, "geo_counts_md5")
openxlsx::writeData(wb, sheet = "geo_counts_md5", geo.counts.md5)
openxlsx::addWorksheet(wb, "geo_raw_md5")
openxlsx::writeData(wb, sheet = "geo_raw_md5", geo.raw.md5)
openxlsx::saveWorkbook(wb,
                       paste0("./Reports/geo_metadata.xlsx"), 
                       overwrite = TRUE)

