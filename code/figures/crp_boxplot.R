# create boxplot of CRP assay results

library(tidyverse)
library(ggplot2)
library(ggprism)

data <- read.csv(file = file.path("data", "CRP_assay", "CRP kit Assay.csv"), header = TRUE, fill = TRUE)

out <- file.path("figures", "CRP_assay_boxplot")
dir.create(out, showWarnings = F, recursive = T)

long_data <- gather(data, key = "SampleID", value = "CRP_Levels", colnames(data))

long_data_filtered <- na.omit(long_data)
long_data_filtered$SampleID <- factor(long_data_filtered$SampleID, levels = c("WT", "WT_3days", "KO_3days", "WT_7days", "KO_7days"))

p <- ggplot(long_data_filtered, aes(x=SampleID, y=`CRP_Levels`)) + geom_boxplot() + theme_prism()

p

ggsave(filename = file.path(out, "CRP_boxplot.png"), plot = p, dpi = 300)



p2 <- long_data_filtered %>% 
  mutate(color_type = ifelse(SampleID=="KO_3days" | SampleID=="KO_7days", "KO", "WT")) %>% 
  ggplot(aes(x=SampleID, y=`CRP_Levels`, fill=color_type)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("grey", "white"), name="", labels=c(expression("Smad4"^"KO IEC"),"WT"))

old <- theme_set(theme_prism())

theme_update(legend.text.align = 0)

p2

ggsave(filename = file.path(out, "CRP_boxplot_grey.png"), plot = p2, dpi = 300)
