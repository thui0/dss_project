# perform ORA analysis for each DESeq2 results object
# uses all ENSEMBL IDs without filtering baseMean 0 genes
# DEG list cutoff of padj < 0.05 & | log2FoldChange | > 0.0

packages <- c("tidyverse", "magrittr", "org.Mm.eg.db", "clusterProfiler", "msigdbr", "enrichplot", "ggplot2", "ggupset", "cowplot")
lapply(packages, library, character.only = TRUE, quietly = TRUE)

# prepare term2gene (t2g) objects for gene set collections
t2g_hallmark <- msigdbr(species = "Mus musculus", category = "H") %>% 
  dplyr::select(gs_name, ensembl_gene)
t2g_oncogenic <- msigdbr(species = "Mus musculus", category = "C6") %>% 
  dplyr::select(gs_name, ensembl_gene)
t2g_C2_cgp <- msigdbr(species = "Mus musculus", category = "C2",subcategory = "CGP") %>% 
  dplyr::select(gs_name, ensembl_gene)
t2g_C2_kegg <- msigdbr(species = "Mus musculus", category = "C2",subcategory = "KEGG") %>% 
  dplyr::select(gs_name, ensembl_gene)
t2g_C2_reactome <- msigdbr(species = "Mus musculus", category = "C2",subcategory = "REACTOME") %>% 
  dplyr::select(gs_name, ensembl_gene)
t2g_C7_immuno <- msigdbr(species = "Mus musculus", category = "C7") %>% 
  dplyr::select(gs_name, ensembl_gene)

t2g_list <- list(
  'Hallmark' = t2g_hallmark,
  "C6_oncogenic" = t2g_oncogenic,
  "C2_CGP" = t2g_C2_cgp,
  "C2_KEGG" = t2g_C2_kegg,
  "C2_Reactome" = t2g_C2_reactome,
  "C7_immunologic" = t2g_C7_immuno
)

# log2FoldChange Threshold
lfc = 0.00

# prepare list of results objects and output directories
out = file.path("results", "ora")
dir.create(out, showWarnings = FALSE)

int_genotypeKO <- readRDS(file.path("results", "dge_analysis", "diffexp", "genotypeKO.treatmentDSS_interaction_term.rds")) %>% as.data.frame()
res_kodss_vs_wtdss <- readRDS(file.path("results", "dge_analysis", "diffexp", "KODSS_vs_WTDSS.rds")) %>% as.data.frame()
res_koc_vs_wtc <- readRDS(file.path("results", "dge_analysis", "diffexp", "KOnoDSS_vs_WTnoDSS.rds")) %>% as.data.frame()
res_ko <- readRDS(file.path("results", "dge_analysis", "diffexp", "KODSS_vs_KOnoDSS.rds")) %>% as.data.frame()
res_wt <- readRDS(file.path("results", "dge_analysis", "diffexp", "WTDSS_vs_WTnoDSS.rds")) %>% as.data.frame()

res_obj <- list(
  "res_int" = list(data = int_genotypeKO, dir = file.path(out,"res_int")),
  "KODSS_vs_WTDSS" = list(data = res_kodss_vs_wtdss, dir = file.path(out,"KODSS_vs_WTDSS")),
  "KOc_vs_WTc" = list(data = res_koc_vs_wtc, dir = file.path(out,"KOc_vs_WTc")),
  "KODSS_vs_KOc" = list(data = res_ko, dir = file.path(out,"KODSS_vs_KOc")),
  "WTDSS_vs_WTc" = list(data = res_wt, dir = file.path(out,"WTDSS_vs_WTc"))
)

res_obj %>% 
  map(~dir.create(.x$dir, showWarnings = FALSE))

# # keep only genes with nonzero counts
# res_obj <- res_obj %>%
#   map(~ list(data = .x$data %>% filter(baseMean > 0), dir = .x$dir))

# create function to perform all desired ORA analyses and output list of results objects
perform_ora <- function(df) {
  # select set of interesting genes and background genes for ORA
  pos_genes <- rownames(df[(!is.na(df$padj) & (df$padj < 0.05) & (df$log2FoldChange > lfc)), ])
  neg_genes <- rownames(df[(!is.na(df$padj) & (df$padj < 0.05) & (df$log2FoldChange < -lfc)), ])
  all_genes <- c(pos_genes,neg_genes)
  universe_genes <- rownames(df)
  
  # create wrapper functions to perform ORA with set parameters
  ## enrichGO wrapper function
  wrap_enrichGO <- function(gene, universe=universe_genes, ont) {
    tryCatch({
      enrichGO(
        gene = gene,
        universe = universe,
        OrgDb = 'org.Mm.eg.db',
        keyType = 'ENSEMBL',
        ont = ont,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500,
        readable = T
      )
    }, error = function(e) NULL)
  }
  ## enricher wrapper function
  wrap_enricher <- function(gene, universe=universe_genes, t2g) {
    tryCatch({
      enricher(
        gene = gene,
        universe = universe,
        pvalueCutoff = 0.05,
        pAdjustMethod = "BH",
        minGSSize = 10,
        maxGSSize = 500,
        qvalueCutoff = 0.2,
        TERM2GENE = t2g
      )
    }, error = function(e) NULL)
  }
  ora_res <- list(
    # GO enrichment analyses
    ## positive degs
    GOBP_pos = wrap_enrichGO(gene = pos_genes, ont = "BP"),
    GOCC_pos = wrap_enrichGO(gene = pos_genes, ont = "CC"),
    GOMF_pos = wrap_enrichGO(gene = pos_genes, ont = "MF"),
    ## negative degs
    GOBP_neg = wrap_enrichGO(gene = neg_genes, ont = "BP"),
    GOCC_neg = wrap_enrichGO(gene = neg_genes, ont = "CC"),
    GOMF_neg = wrap_enrichGO(gene = neg_genes, ont = "MF"),
    ## all degs
    GOBP_all = wrap_enrichGO(gene = all_genes, ont = "BP"),
    GOCC_all = wrap_enrichGO(gene = all_genes, ont = "CC"),
    GOMF_all = wrap_enrichGO(gene = all_genes, ont = "MF"),
    
    # MSigDB gene set collections
    ## positive degs
    Hallmark_pos = wrap_enricher(gene = pos_genes, t2g = t2g_hallmark),
    C6_oncogenic_pos = wrap_enricher(gene = pos_genes, t2g = t2g_oncogenic),
    C2_CGP_pos = wrap_enricher(gene = pos_genes, t2g = t2g_C2_cgp),
    C2_KEGG_pos = wrap_enricher(gene = pos_genes, t2g = t2g_C2_kegg),
    C2_Reactome_pos = wrap_enricher(gene = pos_genes, t2g = t2g_C2_reactome),
    C7_immunologic_pos = wrap_enricher(gene = pos_genes, t2g = t2g_C7_immuno),
    ## negative degs
    Hallmark_neg = wrap_enricher(gene = neg_genes, t2g = t2g_hallmark),
    C6_oncogenic_neg = wrap_enricher(gene = neg_genes, t2g = t2g_oncogenic),
    C2_CGP_neg = wrap_enricher(gene = neg_genes, t2g = t2g_C2_cgp),
    C2_KEGG_neg = wrap_enricher(gene = neg_genes, t2g = t2g_C2_kegg),
    C2_Reactome_neg = wrap_enricher(gene = neg_genes, t2g = t2g_C2_reactome),
    C7_immunologic_neg = wrap_enricher(gene = neg_genes, t2g = t2g_C7_immuno),
    ## all degs
    Hallmark_all = wrap_enricher(gene = all_genes, t2g = t2g_hallmark),
    C6_oncogenic_all = wrap_enricher(gene = all_genes, t2g = t2g_oncogenic),
    C2_CGP_all = wrap_enricher(gene = all_genes, t2g = t2g_C2_cgp),
    C2_KEGG_all = wrap_enricher(gene = all_genes, t2g = t2g_C2_kegg),
    C2_Reactome_all = wrap_enricher(gene = all_genes, t2g = t2g_C2_reactome),
    C7_immunologic_all = wrap_enricher(gene = all_genes, t2g = t2g_C7_immuno)
  )
  
  return(ora_res)
}

# create function to generate pdf report containing all plots for an enrichment result object
create_report <- function(enrich_res, analysis_name, path) {
  if (is.null(enrich_res)) {
    pdf(file = paste0(path, "/", analysis_name, "_plots_NULL_ERROR", ".pdf"))
    plot.new()
    text(
      x = 0.5,
      y = 0.5,
      paste0(
        "ERROR: Enrichment Result is NULL; no enrichment result found for ",
        analysis_name,
        ". \n No plots will be rendered for this analysis."
      )
    )
    dev.off()
    
    warning(paste0("Enrichment Result is NULL; no enrichment result found for ", analysis_name))
    return()
  }
  
  if (nrow(enrich_res) == 0) {
    pdf(file = paste0(path, "/", analysis_name, "_plots_ERROR", ".pdf"))
    plot.new()
    text(
      x = 0.5,
      y = 0.5,
      paste0(
        "ERROR: No significantly enriched terms found for ",
        analysis_name,
        ". \n No plots will be rendered for this analysis."
      )
    )
    dev.off()
    
    warning(paste0("No significantly enriched terms found for ", analysis_name))
    return()
  } else {
    # create safe wrappers for plotting
    safe_plot <- function(plot_call) {
      tryCatch({
        print(plot_call())
      }, error = function(e) {
        plot.new()
        text(x = 0.5, y = 0.5, labels = e$message, adj = c(0.5, 0.5), cex = 1.5)
        message("Error in generating plot: ", e$message)
      })
    }
    
    pdf(file = paste0(path, "/", analysis_name, "_plots", ".pdf"), width = 25, height = 30)
    
    # create goplot
    if ((enrich_res@ontology != "UNKNOWN") & (nrow(enrich_res) > 1)) {
      safe_plot(function() goplot(enrich_res, showCategory = 5, color = "p.adjust") + ggtitle(analysis_name))
    }
    
    # create barplot
    safe_plot(function() barplot(enrich_res, showCategory = 30, x = "Count", color = "p.adjust", title = analysis_name))
    
    # create dotplot
    safe_plot(function() dotplot(enrich_res, showCategory = 30, x = "Count", color = "p.adjust", title = analysis_name))
    
    # create cnetplots
    safe_plot(function() {
      cnet1 <- cnetplot(enrich_res, node_label = "category")
      cnet2 <- cnetplot(enrich_res, node_label = "gene")
      cnet3 <- cnetplot(enrich_res, node_label = "all")
      cowplot::plot_grid(cnet1, cnet2, cnet3, labels = c(paste0(analysis_name, "_category"), paste0(analysis_name, "_gene"), paste0(analysis_name, "_all")))
    })
    
    # create emapplot
    safe_plot(function() {
      emap2 <- emapplot(pairwise_termsim(enrich_res), layout = "nicely")
      aplot::plot_list(emap2, labels = paste0(analysis_name, "_nicely"))
    })
    
    if (nrow(enrich_res) >= 2) {
      # create treeplot
      safe_plot(function() {
        treeplot(pairwise_termsim(enrich_res), hclust_method = "average", nCluster = if (nrow(enrich_res) < 5) { nrow(enrich_res) %% 5 } else { 5 }) + ggtitle(analysis_name)
      })
      
      # create upsetplot
      safe_plot(function() upsetplot(enrich_res) + ggtitle(analysis_name))
      
    }
    dev.off()
  }
}

# create function to generate individual plots
create_plots <- function(enrich_res, analysis_name, path) {
  dir_path <- paste0(path, "/", analysis_name, "_plots")
  dir.create(dir_path)
  
  # check if enrich_res is NULL; enricher most likely could not find DEGs in the T2G list
  if (is.null(enrich_res)) {
    png(filename = paste0(path, "/", analysis_name, "_plots_NULL_ERROR", ".png"))
    plot.new()
    text(
      x = 0.5,
      y = 0.5,
      paste0(
        "ERROR: Enrichment Result is NULL; no enrichment result found for ",
        analysis_name,
        ". \n No plots will be rendered for this analysis."
      )
    )
    dev.off()
    
    warning(paste0("Enrichment Result is NULL; no enrichment result found for ", analysis_name))
    return()
  }
  
  
  # skip plot generation if no significant results are found
  if (nrow(enrich_res) == 0) {
    png(filename = paste0(dir_path, "/", "ERROR_no_sig_res", ".png"))
    plot.new()
    text(
      x = 0.5,
      y = 0.5,
      paste0(
        "ERROR: No significantly enriched terms found for ",
        analysis_name,
        ". \n No plots will be rendered for this analysis."
      )
    )
    dev.off()
    
    warning(paste0("No significantly enriched terms found for ", analysis_name))
    return()
  } else {
    # create safe wrappers for plotting
    safe_plot <- function(filename, width, height, plot_call) {
      tryCatch({
        png(filename,
            width = width,
            height = height,
            units = "in",
            res = 300)
        print(plot_call())
        dev.off()
      }, error = function(e) {
        error_filename <- sub("\\.png$", "_ERROR.png", filename)
        png(
          error_filename,
          width = 1000,
          height = 1200,
          res = 300
        )
        plot.new()
        text(
          x = 0.5,
          y = 0.5,
          labels = e$message,
          adj = c(0.5, 0.5),
          cex = 1.5
        )
        dev.off()
        message("Error in generating plot: ", e$message)
      })
    }
    
    # create goplot
    if (enrich_res@ontology != "UNKNOWN" & (nrow(enrich_res) > 1)) {
      safe_plot(
        filename = paste0(dir_path, "/", "goplot", ".png"),
        width = 25,
        height = 25,
        plot_call = function() {
          goplot(enrich_res,
                 showCategory = 5,
                 color = "p.adjust") + ggtitle(analysis_name)
        }
      )
    }
    
    # create barplot
    safe_plot(
      filename = paste0(dir_path, "/", "barplot", ".png"),
      width = 15,
      height = 20,
      plot_call = function() {
        barplot(
          enrich_res,
          showCategory = 30,
          x = "Count",
          color = "p.adjust",
          title = analysis_name
        )
      }
    )
    
    # create dotplot
    safe_plot(
      filename = paste0(dir_path, "/", "dotplot", ".png"),
      width = 15,
      height = 20,
      plot_call = function() {
        dotplot(
          enrich_res,
          showCategory = 30,
          x = "Count",
          color = "p.adjust",
          title = analysis_name
        )
      }
    )
    
    # create cnetplot
    safe_plot(
      filename = paste0(dir_path, "/", "cnetplot", ".png"),
      width = 25,
      height = 25,
      plot_call = function() {
        cnet1 <- cnetplot(enrich_res, node_label = "category")
        cnet2 <- cnetplot(enrich_res, node_label = "gene")
        cnet3 <- cnetplot(enrich_res, node_label = "all")
        cowplot::plot_grid(cnet1, cnet2, cnet3, labels = c(
          paste0(analysis_name, "_category"),
          paste0(analysis_name, "_gene"),
          paste0(analysis_name, "_all")
        ))
      }
    )
    
    # create emapplot
    safe_plot(
      filename = paste0(dir_path, "/", "emapplot", ".png"),
      width = 15,
      height = 15,
      plot_call = function() {
        emap2 <- emapplot(pairwise_termsim(enrich_res), layout = "nicely")
        aplot::plot_list(emap2, labels = paste0(analysis_name, "_enrichment_map"))
      }
    )
    
    if (nrow(enrich_res) >= 2) {
      # create treeplot
      safe_plot(
        filename = paste0(dir_path, "/", "treeplot", ".png"),
        width = 15,
        height = 15,
        plot_call = function() {
          treeplot(
            pairwise_termsim(enrich_res),
            hclust_method = "average",
            nCluster = if (nrow(enrich_res) < 5) {
              nrow(enrich_res) %% 5
            } else {
              5
            }
          ) + ggtitle(analysis_name)
        }
      )
      
      # create upsetplot
      safe_plot(
        filename = paste0(dir_path, "/", "upsetplot", ".png"),
        width = 15,
        height = 15,
        plot_call = function() {
          print(upsetplot(enrich_res) + ggtitle(analysis_name))
        }
      )
    }
  }
}

################################################################################

# perform ORA analyses

## loop through res_obj and perform_ora on each dataset
for (res in seq_along(res_obj)) {
  res_obj[[res]]$ora_res <- perform_ora(res_obj[[res]]$data)
}

## update readable parameter of each enrichResult Object
for (i in seq_along(res_obj)) {
  current <- res_obj[[i]]
  enrich_res_list <- current$ora_res
  
  # Iterate over each enrichResult in enrich_res_list to update
  ## only update if an results object is present; no sig ora res => return NULL
  for (j in seq_along(enrich_res_list)) {
    if (!is.null(enrich_res_list[[j]])) {
      enrich_res_list[[j]] <- setReadable(enrich_res_list[[j]], OrgDb = "org.Mm.eg.db", keyType = "ENSEMBL")
    }
  }
  
  # Update the enrich_res in res_obj with the modified list
  res_obj[[i]]$ora_res <- enrich_res_list
}

## create report pdf
for (i in seq_along(res_obj)) {
  res_name <- names(res_obj[i])
  current <- res_obj[[i]]
  enrich_res_list <- current$ora_res
  path <- current$dir
  
  # Iterate over each enrichResult in enrich_res_list
  for (j in seq_along(enrich_res_list)) {
    enrich_res <- enrich_res_list[[j]]
    analysis_name <- names(enrich_res_list)[j]
    
    create_report(enrich_res, analysis_name, path)
    
    print(paste("Processed create_report:", res_name, analysis_name))
  }
}

## create individual plots
for (i in seq_along(res_obj)) {
  res_name <- names(res_obj[i])
  current <- res_obj[[i]]
  enrich_res_list <- current$ora_res
  path <- current$dir
  
  # Iterate over each enrichResult in enrich_res_list
  for (j in seq_along(enrich_res_list)) {
    enrich_res <- enrich_res_list[[j]]
    analysis_name <- names(enrich_res_list)[j]
    
    create_plots(enrich_res, analysis_name, path)
    print(paste("Processed create_plot:", res_name, analysis_name))
  }
}

# save res_obj object as RDS file
saveRDS(res_obj, file.path(out,"res_obj.rds"))

# save each enrichResult to csv
for (i in seq_along(res_obj)) {
  res_name <- names(res_obj[i])
  current <- res_obj[[i]]
  enrich_res_list <- current$ora_res
  path <- current$dir
  
  # Iterate over each enrichResult in enrich_res_list
  for (j in seq_along(enrich_res_list)) {
    enrich_res <- enrich_res_list[[j]]
    analysis_name <- names(enrich_res_list)[j]
    filename = paste0(path,"/",analysis_name)
    
    # write results dataframe to csv and txt file
    ## Note: if enrich_res is NULL files will simply write "" (a pair of quotes) in the output file
    enrich_res_df <- as.data.frame(enrich_res)
    write.csv(enrich_res_df, file = paste0(filename,".csv"), quote = T)
    write.table(enrich_res_df, file = paste0(filename,".txt"), quote = F, sep = "\t")
    
    # write results summary to txt log file
    logfile = paste0(filename,"_summary.txt")
    sink(logfile, type = "output"); print(paste0(res_name,"_",analysis_name)); print(enrich_res); sink(type = "output")
  }
}

################################################################################

# Potential Errors to Catch For:

#####

# --> No gene can be mapped....
# --> Expected input gene ID: ENSMUSG00000048895,ENSMUSG00000089678,ENSMUSG00000098488,ENSMUSG00000028479,ENSMUSG00000029171,ENSMUSG00000033161
# --> return NULL...

#####
## The above error message occurs from clusterProfiler::enricher() if none of the genes in the `gene=` param is not found in TERM2GENE dataframe; e.g., test.genes %in% t2g_C2_kegg$ensembl_gene
