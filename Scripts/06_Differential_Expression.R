# Load Metadata
load("./Data/R_Data/metadata.RData", verbose = T)
# Load Count Matrix
load("./Data/R_Data/count_mat.RData", verbose = T)
# Load GTF Object
load("./Data/R_Data/Ens_gtf.RData", verbose = T)

# Create DGEList Object
exp <- DGEList(counts = count_mat)
exp$samples <-cbind(exp$samples,
                    metadata[match(rownames(exp$samples),metadata$Raw.file),])
colnames(exp$counts) <- exp$samples$SampleID
exp$gtf <- Ens.gtf
# Calculate normalization factors to scale the library sizes
exp <- calcNormFactors(exp, method =c("TMM"))
# keep only genes with a minimal number of reads of > # of samples
exp$counts.Filtered <- exp$counts[rowSums(exp$counts) > length(metadata$SampleID),]
# specify the design matrix
Condition <- factor(exp$samples$condition)
Batch <- factor(exp$samples$Batch)
exp$design <- model.matrix(~ 0 + Condition + Batch)

# save Object
save(exp, file = "./Data/R_Data/Exp.RData")


# ----------------
# Voom DE Analysis
# ----------------
# load expression object 
load("./Data/R_Data/Exp.RData", verbose = T)

Voom <- voom(exp$counts.Filtered, design = exp$design, 
             lib.size = exp$samples$lib.size * exp$samples$norm.factors, 
             plot = T, normalize.method = "quantile")
# get the transpose to have the groups as the individuals (expected format for PCA function)
Voom$pca <- PCA(X = t(Voom$E), graph = F, scale.unit = FALSE)
# get the coordinates to plot the PCA customized using ggplot
Voom$pca$ind$coord <- cbind(Voom$pca$ind$coord, exp$samples)
Voom$pca$ind$coord <- Voom$pca$ind$coord[,unique(colnames(Voom$pca$ind$coord))]

# create the PCA plot using ggplot
g.PCA <-
ggplot(data = Voom$pca$ind$coord, aes(x= Dim.1, y = Dim.2, col = condition)) + 
  geom_point(shape = 1, size = 5, colour = "black")+
  geom_label_repel(aes(label = Voom$pca$ind$coord$SampleID),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50', 
                   size = 4, 
                   show.legend = FALSE) +
  theme_cowplot() +
  xlab(label = paste0("Dim.1 (", signif(Voom$pca$eig[1,2], 3), "%)")) +
  ylab(label = paste0("Dim.2 (", signif(Voom$pca$eig[2,2], 3), "%)")) +
  ggtitle('Transcriptomics PCA')+
  
  scale_color_brewer(palette ="Set1" )+
  theme(
    panel.background = element_blank(),
    axis.line = element_line(colour = "grey"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "plain")
  )+
  geom_rug('outside' = TRUE) +
  coord_cartesian(clip = "off")

ggsave(g.PCA, file = "./Plots/Voom/PCA.png", width = 6, height = 4)

# Construct the contrast matrix corresponding to specified contrasts
Voom$contrasts <- makeContrasts(Diet_Effect = ConditionHFD - ConditionCD,
                                levels = exp$design)
# fit a linear model
Voom$fit <- lmFit(Voom, design = exp$design)
Voom$fit.contrasts <- contrasts.fit(Voom$fit, contrasts = Voom$contrasts)
Voom$fit.eBayes <- eBayes(Voom$fit.contrasts, robust = TRUE)
Voom$fit.decideTests <- decideTests(Voom$fit.eBayes, lfc = log2(2), p.value = 0.05, adjust.method = "BH")
summary(Voom$fit.decideTests)

Voom$Results.Table <- topTable(Voom$fit.eBayes, n = Inf)
Voom$Results.Table$Contrast <- colnames(Voom$contrasts)[1]
Voom$Results.Table <-
cbind(Voom$Results.Table,
      exp$gtf[match(rownames(Voom$Results.Table), 
                    exp$gtf$gene_id), 
              c("gene_id", "gene_name", "gene_biotype", "gene_source", "start", "end")])

# save Results Table as excel sheet
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, colnames(Voom$contrasts)[1])
openxlsx::writeData(wb, sheet = colnames(Voom$contrasts)[1], Voom$Results.Table)
openxlsx::saveWorkbook(wb = wb, file = "./Reports/Voom_Results_Tables_DE.xlsx", overwrite = TRUE)

# save Voom Object
save(Voom, file = "./Data/R_Data/Exp_Voom.RData")


# -----------------
# EdgeR DE Analysis
# -----------------
# load expression object 
load("./Data/R_Data/Exp.RData", verbose = T)

EdgeR <- exp
EdgeR <- edgeR::estimateGLMCommonDisp(EdgeR,
                                      design = EdgeR$design,
                                      method = "CoxReid")
EdgeR <- edgeR::estimateGLMRobustDisp(EdgeR, 
                                      design = EdgeR$design,
                                      trend.method = 'auto')

EdgeR$glmfit <- edgeR::glmFit(EdgeR, design = EdgeR$design)
EdgeR$glmlrt <- edgeR::glmLRT(EdgeR$glmfit, 
                              contrast = makeContrasts(Diet_Effect = ConditionHFD - ConditionCD,
                                                                     levels = exp$design))
EdgeR$Results.Table <- EdgeR$glmlrt$table
EdgeR$Results.Table$adj.P.Val <- p.adjust(EdgeR$Results.Table$PValue, method = 'BH')
EdgeR$Results.Table$Contrast <- "Diet_Effect"
EdgeR$Results.Table <-
  cbind(EdgeR$Results.Table,
        exp$gtf[match(rownames(EdgeR$Results.Table), 
                      exp$gtf$gene_id), 
                c("gene_id", "gene_name", "gene_biotype", "gene_source", "start", "end")])


# save Results Table as excel sheet
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Diet_Effect")
openxlsx::writeData(wb, sheet = "Diet_Effect", EdgeR$Results.Table)
openxlsx::saveWorkbook(wb = wb, file = "./Reports/EdgeR_Results_Tables_DE.xlsx", overwrite = TRUE)

# save EdgeR Object
save(EdgeR, file = "./Data/R_Data/Exp_EdgeR.RData")


# ---------------
# DESeq2 Analysis
# ---------------
# load expression object 
load("./Data/R_Data/Exp.RData", verbose = T)

Info <- exp$samples
rownames(Info) <- Info$SampleID
DESeq2 <- DESeq2::DESeqDataSetFromMatrix(countData = exp$counts.Filtered,
                                         colData = Info, 
                                         design = exp$design)
DESeq2 <- DESeq2::DESeq(DESeq2)
DESeq2::resultsNames(DESeq2)
DESeq2.Results.Table <- data.frame(DESeq2::results(object = DESeq2, contrast = list(c("ConditionHFD"), c("ConditionCD"))))
DESeq2.Results.Table$Contrast <- "Diet_Effect"
DESeq2.Results.Table <-
  cbind(DESeq2.Results.Table,
        exp$gtf[match(rownames(DESeq2.Results.Table), 
                      exp$gtf$gene_id), 
                c("gene_id", "gene_name", "gene_biotype", "gene_source", "start", "end")])

# save Results Table as excel sheet
wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Diet_Effect")
openxlsx::writeData(wb, sheet = "Diet_Effect", DESeq2.Results.Table)
openxlsx::saveWorkbook(wb = wb, file = "./Reports/DESeq2_Results_Tables_DE.xlsx", overwrite = TRUE)

# save EdgeR Object
save(DESeq2, DESeq2.Results.Table, file = "./Data/R_Data/Exp_DESeq2.RData")


# ---------------------
# DE Results Comparison
# ---------------------
load("./Data/R_Data/Exp_Voom.RData", verbose = T)
load("./Data/R_Data/Exp_EdgeR.RData", verbose = T)
load("./Data/R_Data/Exp_DESeq2.RData", verbose = T)

DE.Voom <- Voom$Results.Table[abs(Voom$Results.Table$logFC) > 1 & Voom$Results.Table$adj.P.Val < 0.05,]$gene_id
DE.EdgeR <- EdgeR$Results.Table[abs(EdgeR$Results.Table$logFC) > 1 & EdgeR$Results.Table$adj.P.Val < 0.05,]$gene_id
DE.DESeq2 <- rownames(DESeq2.Results.Table[abs(DESeq2.Results.Table$log2FoldChange) > 1 & 
                                             !is.na(DESeq2.Results.Table$padj) & 
                                             DESeq2.Results.Table$padj < 0.05,])

DE.Venn <-
  VennDiagram::venn.diagram(
    x = list(DE.Voom, 
             DE.EdgeR,
             DE.DESeq2),
    category.names = c("Voom", "EdgeR", "DESeq2"),
    filename = NULL,
    
    # Output features
    imagetype="png" ,
    height = 480 , 
    width = 480 , 
    resolution = 300,
    compression = "lzw",
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("#abff7a", "#ff7a7a", "#ffe27a"),
    
    # Numbers
    cex = .6,
    fontface = "bold",
    fontfamily = "sans",
    
    # Set names
    cat.cex = 0.6,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.fontfamily = "sans",
  )

pdf("./Plots/DE_Comparison.pdf", 
    height = 6, width = 6, useDingbats=FALSE)
grid::grid.draw(DE.Venn, recording = F)
dev.off()

