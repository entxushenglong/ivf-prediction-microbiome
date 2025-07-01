# Core plotting code for the study:
# "Machine Learning Improves IVF Success Prediction via Multi-site, Species-resolved Reproductive Microbiota"
# This script provides example code for generating main figures using microbiome data.


# Load required libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(microeco)
library(patchwork)

# Construct microtable object
dataset <- microtable$new(sample_table = sample_data,
                          otu_table = table, 
                          tax_table = taxonomy)
dataset$tidy_dataset()

# Format group labels
dataset$sample_table <- dataset$sample_table %>%
  mutate(Site = recode(Site, YD = "Vagina", GJ = "Cervix", YZG = "Uterus"),
         Implantation = ifelse(Implantation == "Yes", "Succeed", "Fail"))

# Define unified theme for all plots
mytheme <- theme_bw() + theme(
  panel.grid = element_blank(),
  legend.position = "top",
  axis.title = element_text(size = 12),
  axis.text = element_text(size = 12),
  axis.line = element_line(size = 0.5, color = "black")
)

# 1. Phylum-level barplot
p_phylum <- trans_abund$new(dataset, taxrank = "phylum", ntaxa = 8, groupmean = "Implantation")$plot_bar(
  others_color = "grey70", xtext_angle = 50, legend_text_italic = FALSE) +
  mytheme

# 2. Species-level barplot
p_species <- trans_abund$new(dataset, taxrank = "species", ntaxa = 8, groupmean = "Implantation")$plot_bar(
  others_color = "grey70", xtext_angle = 50, legend_text_italic = FALSE) +
  mytheme

# 3. Alpha diversity plots (Shannon and Inverse Simpson indices)
alpha_obj <- trans_alpha$new(dataset, group = "Implantation")
alpha_obj$cal_diff(method = "wilcox")

p_shannon <- alpha_obj$plot_alpha(measure = "Shannon", shape = "Implantation",
                                   color_values = c("#D2691E", "#4682B4")) + mytheme

p_invsimpson <- alpha_obj$plot_alpha(measure = "InvSimpson", shape = "Implantation",
                                     color_values = c("#D2691E", "#4682B4")) + mytheme

# 4. Beta diversity PCoA plot with PERMANOVA
dataset$cal_betadiv(unifrac = FALSE)
beta_obj <- trans_beta$new(dataset, group = "Implantation", measure = "bray")
beta_obj$cal_ordination(method = "PCoA")
beta_obj$cal_manova(manova_all = TRUE)

p_beta <- beta_obj$plot_ordination(plot_color = "Implantation", plot_shape = "Implantation",
                                   plot_type = c("point", "ellipse"),
                                   color_values = c("#D2691E", "#4682B4")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.5, vjust = 1.5,
           label = paste0("PERMANOVA p = ", beta_obj$res_manova$`Pr(>F)`[1]), size = 4.5) +
  mytheme

# 5. Differential abundance analysis (Wilcoxon rank-sum test)
diff_obj <- trans_diff$new(dataset, method = "wilcox", group = "Implantation", taxa_level = "species")
p_diff <- diff_obj$plot_diff_abund(use_number = 1:5, add_sig = TRUE, simplify_names = TRUE,
                                   color_values = c("#D2691E", "#4682B4")) + mytheme

# Example to save a plot
# ggsave("figure_alpha_shannon.pdf", p_shannon, width = 5, height = 4)
