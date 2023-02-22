# Load DE data analysis
load("./Data/R_Data/Exp_Voom.RData", verbose = T)
load("./Data/R_Data/Exp_EdgeR.RData", verbose = T)
load("./Data/R_Data/Exp_DESeq2.RData", verbose = T)
load("./Data/R_Data/Exp.RData", verbose = T)

# Read data object of Sc_RNA-seq & DDC sample
ScRNA_Seq_DDC_vs_Control <- read.xlsx("./Resources/Sc_Data/1-s2.0-S1934590919301535-mmc4.xlsx", 
                                      sheet = "Sheet1", startRow = 2)[, 1:6]
ScRNA_Seq_DDC_vs_Control$gene_id = mapIds(org.Mm.eg.db,
                                            keys=ScRNA_Seq_DDC_vs_Control$Gene.Symbol, 
                                            column="ENSEMBL",
                                            keytype="SYMBOL",multiVals="first")
colnames(ScRNA_Seq_DDC_vs_Control) <- gsub("Gene.Symbol", "gene_name", colnames(ScRNA_Seq_DDC_vs_Control))

# Retrieved: http://ge-lab.org/gskb/
gmt_GO <- tail(read.gmt("./Resources/GeneSets/MousePath_GO_gmt.gmt"), -1)
gmt_mGSKB <- read.gmt("./Resources/GeneSets/mGSKB_Ensembl.gmt")
gmt_Pathway <- tail(read.gmt("./Resources/GeneSets/MousePath_Pathway_gmt.gmt"), -1)
gmt_Metabolic <- tail(read.gmt("./Resources/GeneSets/MousePath_Metabolic_gmt.gmt"), -1)
gmt_TF <- tail(read.gmt("./Resources/GeneSets/MousePath_TF_gmt.gmt"), -1)

gmt_list <- list(gmt_GO, gmt_mGSKB, gmt_Pathway, gmt_Metabolic, gmt_TF)
gmt_list_names <- c("GO", "mGSKB", "Pathway", "Metabolic", "TF")

# ---------------------------
# Enrichment Overlap Analysis
# ---------------------------
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

DDC_up <- ScRNA_Seq_DDC_vs_Control[ScRNA_Seq_DDC_vs_Control$`Fold.Change.(log2)` < -1 & 
                                       ScRNA_Seq_DDC_vs_Control$`p-value.*` < 0.05,]
DDC_down <- ScRNA_Seq_DDC_vs_Control[ScRNA_Seq_DDC_vs_Control$`Fold.Change.(log2)` > 1 & 
                                       ScRNA_Seq_DDC_vs_Control$`p-value.*` < 0.05,]

# Read data from "Epigenetic remodelling licences adult cholangiocytes for organoid formation and liver regeneration"
OrgT0.vs.InVivo.DE <- readxl::read_xlsx("./Resources/Organoid_Paper/41556_2019_402_MOESM4_ESM.xlsx", 
                                        sheet = "8. T0 vs organoids", skip = 1)
colnames(OrgT0.vs.InVivo.DE) <- c("gene_name", "gene_id", colnames(OrgT0.vs.InVivo.DE)[3:36])
OrgT0.vs.InVivo_up <- OrgT0.vs.InVivo.DE[OrgT0.vs.InVivo.DE$b < -1 & OrgT0.vs.InVivo.DE$qval < 0.05,]
OrgT0.vs.InVivo_down <- OrgT0.vs.InVivo.DE[OrgT0.vs.InVivo.DE$b > 1 & OrgT0.vs.InVivo.DE$qval < 0.05,]

DE_list <- list(Voom_up, Voom_down, EdgeR_up, EdgeR_down, DESeq2_up, DESeq2_down, DDC_up, DDC_down, OrgT0.vs.InVivo_up, OrgT0.vs.InVivo_down)
DE_data_list <- c("Voom_up", "Voom_down", "EdgeR_up", "EdgeR_down", "DESeq2_up", "DESeq2_down", "DDC_up", "DDC_down", "OrgT0.vs.InVivo_up", "OrgT0.vs.InVivo_down")

Use.Background <- FALSE
En.list <-
mapply(condition = DE_data_list, function(condition){
  x <- match(condition, DE_data_list)
  
  if(!Use.Background){
    # Specify Paths
    Plots_path = paste0('./Plots/DDC_Comparison/Enrichment_Overlap_Analysis_No_Background/', 
                        gsub("(.*)_.*", "\\1", condition))
    Reports_path = paste0("./Reports/DDC_Comparison/Enrichment_Overlap_Analysis_No_Background/", 
                          gsub("(.*)_.*", "\\1", condition))
  } else {
    # Specify Paths
    Plots_path = paste0('./Plots/DDC_Comparison/Enrichment_Overlap_Analysis/', 
                        gsub("(.*)_.*", "\\1", condition))
    Reports_path = paste0("./Reports/DDC_Comparison/Enrichment_Overlap_Analysis/", 
                          gsub("(.*)_.*", "\\1", condition))
  }
  
  if(!dir.exists(Plots_path)){
    dir.create(Plots_path, recursive = T)
  }
  
  if(!dir.exists(Reports_path)){
    dir.create(Reports_path, recursive = T)
  }
  
  # Create Excel object
  wb <- openxlsx::createWorkbook()
  # Get gene ID
  gene.ID <- DE_list[[x]]$gene_id
  gene.Name <- DE_list[[x]]$gene_name
  
  # Get background based on samples
  if(grepl("^DDC_", condition)){
    # Specify background from DDC file reported genes
    universe.ids <- ScRNA_Seq_DDC_vs_Control$gene_id
    universe.names <- ScRNA_Seq_DDC_vs_Control$gene_name
  } else if (grepl("^OrgT0.vs", condition)){
    # Specify background from DDC file reported genes
    universe.ids <- OrgT0.vs.InVivo.DE$gene_id
    universe.names <- OrgT0.vs.InVivo.DE$gene_name
  } else {
    # Specify background based on all measured genes
    universe.ids <- rownames(exp$counts)
    universe.names <- exp$gtf[match(universe.ids, exp$gtf$gene_id), "gene_name"] 
  }
  
  en.tmp <-
  mapply(gmt = gmt_list_names, function(gmt){
    
    if(Use.Background){
      if(gmt == "mGSKB"){
        en <- enricher(gene = gene.ID,
                       universe = universe.ids,
                       TERM2GENE = gmt_list[[match(gmt, gmt_list_names)]])
        
        write.csv(data.frame(en@result), 
                  file = paste0(Reports_path, "/Enrichement_", condition, "_", gmt, ".csv"))
      } else {
        en <- enricher(gene = toupper(gene.Name),
                       universe = toupper(unique(universe.names)),
                       TERM2GENE = gmt_list[[match(gmt, gmt_list_names)]])
        
        # export GE results
        openxlsx::addWorksheet(wb, gmt)
        openxlsx::writeData(wb, sheet = gmt, data.frame(en@result))
      } 
    } else {
      if(gmt == "mGSKB"){
        en <- enricher(gene = gene.ID,
                       TERM2GENE = gmt_list[[match(gmt, gmt_list_names)]])
        
        write.csv(data.frame(en@result), 
                  file = paste0(Reports_path, "/Enrichement_", condition, "_", gmt, ".csv"))
      } else {
        en <- enricher(gene = toupper(gene.Name),
                       TERM2GENE = gmt_list[[match(gmt, gmt_list_names)]])
        
        # export GE results
        openxlsx::addWorksheet(wb, gmt)
        openxlsx::writeData(wb, sheet = gmt, data.frame(en@result))
      }
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
      
      # Generate Enrichmap Plot                
      # en.enrichplot <- enrichplot::emapplot(enrichplot::pairwise_termsim(en),
      #                                       pie="count", cex_category=1.5, layout="kk")
      # ggsave(en.enrichplot,
      #        file = paste0(Plots_path, "/", condition, "_GE_EnrichPlot_", gmt, ".pdf"),
      #        width=20, height=20)
    }
    
    return(en@result[en@result$p.adjust < 0.05, "Description"])
  }, SIMPLIFY = FALSE)
  
  openxlsx::saveWorkbook(wb,
                         paste0(Reports_path, "/Enrichement_", 
                                condition, "_",".xlsx"), 
                         overwrite = TRUE)
  return(en.tmp)
}, SIMPLIFY = FALSE)

if(!Use.Background){
  saveRDS(En.list, "./Data/R_Data/DDC_Enrichment_Overlap_Analysis_No_Background.RDS")
} else {
  saveRDS(En.list, "./Data/R_Data/DDC_Enrichment_Overlap_Analysis.RDS")
}


# mapply(cond = c("up", "down"), function(cond){
#   mapply(geneset =  c("GO", "mGSKB", "Pathway", "Metabolic", "TF"), function(geneset){
#     DE.Venn <-
#       VennDiagram::venn.diagram(
#         x = list(En.list[[paste0("Voom_",cond)]][[geneset]], 
#                  En.list[[paste0("EdgeR_",cond)]][[geneset]],
#                  En.list[[paste0("DESeq2_",cond)]][[geneset]],
#                  En.list[[paste0("DDC_",cond)]][[geneset]]),
#         category.names = c("Voom", "EdgeR", "DESeq2", "DDC"),
#         filename = NULL,
#         
#         # Output features
#         imagetype="png" ,
#         height = 480 , 
#         width = 480 , 
#         resolution = 300,
#         compression = "lzw",
#         
#         # Circles
#         lwd = 2,
#         lty = 'blank',
#         fill = c("#abff7a", "#ff7a7a", "#ffe27a", "#7ac3ff"),
#         
#         # Numbers
#         cex = .6,
#         fontface = "bold",
#         fontfamily = "sans",
#         
#         # Set names
#         cat.cex = 0.6,
#         cat.fontface = "bold",
#         cat.default.pos = "outer",
#         cat.fontfamily = "sans",
#       )
#     
#     pdf(paste0("./Plots/DDC_Comparison/Enrichment_Overlap_Analysis/Venn_",cond, "_", geneset, ".pdf"), 
#         height = 6, width = 6, useDingbats=FALSE)
#     grid::grid.draw(DE.Venn, recording = F)
#     dev.off()
#     
#     write(En.list[[paste0("DDC_",cond)]][[geneset]][
#       En.list[[paste0("DDC_",cond)]][[geneset]] %in% c(En.list[[paste0("Voom_",cond)]][[geneset]],
#                                                        En.list[[paste0("EdgeR_",cond)]][[geneset]],
#                                                        En.list[[paste0("DESeq2_",cond)]][[geneset]])],
#       file = paste0("./Reports/DDC_Comparison/Enrichment_Overlap_Analysis/Venn_",cond, "_", geneset, ".txt"))
#   }, SIMPLIFY = FALSE)
# }, SIMPLIFY = FALSE)
