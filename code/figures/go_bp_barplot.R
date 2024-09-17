# re-create barplot of GO_BP ORA results for manuscript figure with appropriate fonts and details
# KODSS vs WTDSS ORA GO_BP_pos

library(tidyverse)
library(clusterProfiler)
library(msigdbr)
library(enrichplot)

# load rds results object
res_obj <- readRDS(file.path("results", "ora", "res_obj.rds"))

# extract relevant result from res_obj
ora_result <- res_obj$KODSS_vs_WTDSS$ora_res$GOBP_pos

p <- ora_result %>% barplot(showCategory = 20, x="Count", color = "p.adjust", title = "KODSS vs WTDSS GO:BP")
p

old <- p$theme

theme_set(old)

p2 <- p + theme_set(old) + theme(axis.text.y = element_text(size = 20))
# theme(axis.text = element_text(size = 30))
p2

out <- file.path("figures", "barplot")
dir.create(out, showWarnings = F, recursive = T)
ggsave(filename = file.path(out, "KODSS_vs_WTDSS_GOBP_pos_barplot.png"), plot = p2, width = 12, height = 16, units = "in", dpi = 300)
