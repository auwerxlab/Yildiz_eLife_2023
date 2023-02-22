# Retrieved: http://ge-lab.org/gskb/
gmt_GO <- tail(read.gmt("./Resources/GeneSets/MousePath_GO_gmt.gmt"), -1)
gmt_mGSKB <- read.gmt("./Resources/GeneSets/mGSKB_Ensembl.gmt")
gmt_Pathway <- tail(read.gmt("./Resources/GeneSets/MousePath_Pathway_gmt.gmt"), -1)
gmt_Metabolic <- tail(read.gmt("./Resources/GeneSets/MousePath_Metabolic_gmt.gmt"), -1)
gmt_TF <- tail(read.gmt("./Resources/GeneSets/MousePath_TF_gmt.gmt"), -1)

gmt_list <- list(gmt_GO, gmt_mGSKB, gmt_Pathway, gmt_Metabolic, gmt_TF)
gmt_list_names <- c("GO", "mGSKB", "Pathway", "Metabolic", "TF")

# Load Results from DE Analyses
load("./Data/R_Data/Exp_Voom.RData", verbose = T)
load("./Data/R_Data/Exp_EdgeR.RData", verbose = T)
load("./Data/R_Data/Exp_DESeq2.RData", verbose = T)

Voom_up <- Voom$Results.Table[Voom$Results.Table$logFC > 1 & Voom$Results.Table$adj.P.Val < 0.05,]
Voom_down <- Voom$Results.Table[Voom$Results.Table$logFC < -1 & Voom$Results.Table$adj.P.Val < 0.05,]

EdgeR_up <- EdgeR$Results.Table[EdgeR$Results.Table$logFC > 1 & EdgeR$Results.Table$adj.P.Val < 0.05,]
EdgeR_down <- EdgeR$Results.Table[EdgeR$Results.Table$logFC < -1 & EdgeR$Results.Table$adj.P.Val < 0.05,]

DESeq2_up <- DESeq2.Results.Table[DESeq2.Results.Table$log2FoldChange > 1 & 
                                    !is.na(DESeq2.Results.Table$padj) & 
                                    DESeq2.Results.Table$padj < 0.05,]
DESeq2_down <- DESeq2.Results.Table[DESeq2.Results.Table$log2FoldChange < -1 & 
                                      !is.na(DESeq2.Results.Table$padj) & 
                                      DESeq2.Results.Table$padj < 0.05,]

DE_list <- list(Voom_up, Voom_down, EdgeR_up, EdgeR_down, DESeq2_up, DESeq2_down)
DE_data_list <- c("Voom_up", "Voom_down", "EdgeR_up", "EdgeR_down", "DESeq2_up", "DESeq2_down")

for(DE_data in 1:6) {
  # Get condition name
  condition <- DE_data_list[DE_data]
  # Specify Paths
  Plots_path = paste0('./Plots/', gsub("(.*)_.*", "\\1", condition), "/Enrichement_Analysis/")
  Reports_path = paste0("./Reports/")
  # Create Excel object
  wb <- openxlsx::createWorkbook()
  # Get gene IDs
  gene.ID <- rownames(DE_list[[DE_data]])
  gene.Name <- DE_list[[DE_data]]$gene_name

  for(gmt in gmt_list_names){
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
             height = 6, width = 6)
      
      # Generate Upset Plot
      en.upsetplot <- enrichplot::upsetplot(en)
      ggsave(plot = en.upsetplot, 
             file = paste0(Plots_path, "/", condition, "_GE_UpsetPlot_", gmt, ".pdf"),
             width=20, height=8)
      
      # Generate Enrichmap Plot                
      en.enrichplot <- enrichplot::emapplot(enrichplot::pairwise_termsim(en),
                                            pie="count", cex_category=1.5, layout="kk")
      ggsave(en.enrichplot,
             file = paste0(Plots_path, "/", condition, "_GE_EnrichPlot_", gmt, ".pdf"),
             width=20, height=20)
    }
  }
  
  openxlsx::saveWorkbook(wb,
                         paste0("./Reports/Enrichement_", condition, "_", gmt,".xlsx"), 
                         overwrite = TRUE)
}

# ####################################################################################################################################
# 2. GSEA Analysis
# ####################################################################################################################################
# Retrieved: http://ge-lab.org/gskb/
gmt_GO <- tail(read.gmt("./Resources/GeneSets/MousePath_GO_gmt.gmt"), -1)
gmt_mGSKB <- read.gmt("./Resources/GeneSets/mGSKB_Ensembl.gmt")
gmt_Pathway <- tail(read.gmt("./Resources/GeneSets/MousePath_Pathway_gmt.gmt"), -1)
gmt_Metabolic <- tail(read.gmt("./Resources/GeneSets/MousePath_Metabolic_gmt.gmt"), -1)
gmt_TF <- tail(read.gmt("./Resources/GeneSets/MousePath_TF_gmt.gmt"), -1)
gmt_hallmark <- read.gmt("./Resources/GeneSets/h.all.v7.5.1.symbols.gmt")

gmt_list <- list(gmt_GO, gmt_mGSKB, gmt_Pathway, gmt_Metabolic, gmt_TF)
gmt_list_names <- c("GO", "mGSKB", "Pathway", "Metabolic", "TF", "gseGO", "KEGG")

# Load Results from DE Analyses
load("./Data/R_Data/Exp_Voom.RData", verbose = T)
load("./Data/R_Data/Exp_EdgeR.RData", verbose = T)
load("./Data/R_Data/Exp_DESeq2.RData", verbose = T)

geneList <- DESeq2.Results.Table$log2FoldChange
names(geneList) <- DESeq2.Results.Table$gene_name
# sort by decreasing order
geneList <- sort(geneList, decreasing = T)
# Get Entrez IDs
mart.ensembl <- biomaRt::useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
Ensembl.to.Entrez <- biomaRt::getBM(attributes = c("external_gene_name", 
                                     "entrezgene_id"), 
                      filters    = "external_gene_name", 
                      values     = names(geneList),
                      mart       = mart.ensembl,
                      uniqueRows = T,
                      useCache   = F)

Ensembl.to.Entrez$logFC <- geneList[match(Ensembl.to.Entrez$external_gene_name, names(geneList))]
# remove NA ones
Ensembl.to.Entrez <- Ensembl.to.Entrez[!is.na(Ensembl.to.Entrez$entrezgene_id),]
geneList.Entrez <- Ensembl.to.Entrez$logFC
names(geneList.Entrez) <- Ensembl.to.Entrez$entrezgene_id
# sort by decreasing order
geneList.Entrez <- sort(geneList.Entrez, decreasing = T)

# Prepare same list but with ensembl ids
geneList.Ensembl <- DESeq2.Results.Table$log2FoldChange
names(geneList.Ensembl) <- DESeq2.Results.Table$gene_id
# sort by decreasing order
geneList.Ensembl <- sort(geneList.Ensembl, decreasing = T)
names(geneList) <- toupper(names(geneList))

DE.gsea <-
mapply(geneset = 1:length(gmt_list_names), function(geneset){
  print(gmt_list_names[geneset])
  if(gmt_list_names[geneset] == "KEGG"){
    en.tmp <- clusterProfiler::gseKEGG(geneList = geneList.Entrez, organism = "mmu")
  } else if(gmt_list_names[geneset] == "gseGO"){
    en.tmp <- clusterProfiler::gseGO(geneList = geneList.Entrez, ont = "BP", OrgDb = org.Mm.eg.db)
  } else if (gmt_list_names[geneset] == "mGSKB") {
    en.tmp <- clusterProfiler::GSEA(geneList = geneList.Ensembl, 
                                    TERM2GENE = gmt_list[[geneset]])
  } else {
    en.tmp <- clusterProfiler::GSEA(geneList = geneList, 
                          TERM2GENE = gmt_list[[geneset]])
  }
  return(en.tmp)
}, SIMPLIFY = FALSE)
names(DE.gsea) <- gmt_list_names
    
save(DE.gsea, file = "./Data/R_Data/Exp_DESeq2_GSEA.RData")

# Print name of genesets related to fat / lipid
for(x in names(DE.gsea)){
  print(DE.gsea[[x]]@result$Description[grepl("METABOLISM", DE.gsea[[x]]@result$Description)])
}

# Create Excel object
wb <- openxlsx::createWorkbook()
for(geneset in names(DE.gsea)){
  comp.en.tmp <- DE.gsea[[geneset]]
    
  if(as.character(typeof(DE.gsea[[geneset]]) == "S4")){
    print(geneset)
    # Save results table to file
    openxlsx::addWorksheet(wb, geneset)
    df <- as.data.frame(comp.en.tmp@result)
    df <-
    df %>%
      mutate_all(stringi::stri_enc_toutf8)
    df <- df[validUTF8(df$ID),]
    openxlsx::writeData(wb, sheet = geneset, df)
      
    if(sum(comp.en.tmp@result$p.adjust < 0.05) > 0){
      # Generate default Dot Plot
      en.dotplot <- dotplot(comp.en.tmp, showCategory = 15)
      ggsave(plot = en.dotplot, 
              paste0("Plots/DESeq2/GSEA/", geneset, "_DotPlot.pdf"), 
              height = 12, width = 12)
        
      # Proper GSEA Plot
      result <- comp.en.tmp@result
      result$core_enrichmentNgenes <- sapply(strsplit(unlist(result$core_enrichment), 
                                                      split = "/", fixed = T), length)
      result$core_enrichmentgenes <- unlist(sapply(strsplit(unlist(result$core_enrichment), 
                                                            split = "/", fixed = T), 
                                                    function(x){paste(sort(x), collapse = " ")}))
      result$signficant <- result$qvalues <= 0.01
      result$label <- gsub("\\(.*\\)", "", gsub("_", " ", result$Description))
      result$gene_ratio <- result$core_enrichmentNgenes / result$setSize
        
      # sort by abs value of NES and take top 15 values
      tmp_data <- head(result[order(-abs(result$NES)),], 30)
        
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
              paste0("Plots/DESeq2/GSEA/", 
                    geneset, "_GSEA_DotPlot.pdf"), 
              height = 12, width = 12)
        
      # Generate Enrichmap Plot                
      en.enrichplot <- enrichplot::emapplot(enrichplot::pairwise_termsim(comp.en.tmp),
                                             pie="count", cex_category=1.5)
      ggsave(en.enrichplot,
              file = paste0("Plots/DESeq2/GSEA/", 
                            geneset, "_EnrichPlot.pdf"),
              width=20, height=20)
    }
  }
}
  
openxlsx::saveWorkbook(wb,
                        paste0("Reports/GSEA/", "DESeq2_GSEA", ".xlsx"), 
                        overwrite = TRUE)
