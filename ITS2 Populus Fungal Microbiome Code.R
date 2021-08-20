setwd("~/Downloads")

# Download libraries
# CRAN Packages
library(cowplot)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(knitr)
library(readxl)
library(stringr)
library(vegan)
library(viridis)
library(tidyr)
library(doParallel)
library(tibble)
library(SpadeR)
library(nlme)
library(hillR)
library(yaml)
library(multcomp)
library(ape)
library(biomformat)
library(Biostrings)
library(DESeq2)
library(phyloseq)
library(microbiome)
library(ape)
library("dplyr") 
library("magrittr")

# Update BiocManager Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
BiocManager::install(c("DESeq2", "microbiome"))

# Install qiime2 and load qiime2 library 
if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")
library(qiime2R)

# Set theme for ggplot figures
theme_set(theme_cowplot())
theme_update(axis.title = element_text(size = 20), axis.text = element_text(size = 15), 
             legend.title = element_text(size = 20), legend.text = element_text(size = 15),
             legend.text.align = 0, 
             title = element_text(size = 20))


# Import Fungal Data
ASV_ITS <- read_qza("table_merge_ASV_ITS.qza")
library(phyloseq)

# Convert ASV_ITStable into .csv files 
## Will drop converted files into Downloads folder
ASV_ITStable <- ASV_ITS$data
write.csv(ASV_ITStable, file = "table_merge_ASV_ITS.csv")

# Convert Taxonomy_ITS table into .csv file 
Taxonomy_ITS <- read_qza("Taxonomy_ITS.qza")
write.csv(Taxonomy_ITS$table, file = "Taxonomy_ITS.csv", row.names = FALSE)

## Reformat taxonomy table
Taxnonomy_ITS <- Taxonomy_ITS$data %>% as.tibble() %>%
  mutate(Taxon=gsub("[a-z]__", "", Taxon)) %>%
  separate(Taxon, sep = ";",
           c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) %>%
  mutate(Phylum = replace_na(Phylum, "unidentified"))
write.csv(Taxnonomy_ITS, file = "Taxonomy_ITS.csv", row.names = FALSE)

# Concert Metadata_ITS into .csv file 
Metadata_ITS <- read.csv("Metadata_ITS.csv")

# Import data into phyloseq to create new phyloseq object
data_ITS <- read_csv2phyloseq(otu.file = "table_merge_ASV_ITS.csv",
                              taxonomy.file = "Taxonomy_ITS.csv",
                              metadata.file = "Metadata_ITS.csv")

# Pre-processing step: Trim data and make into proportion
data_ITS_counts <- data_ITS
data_ITS <- prune_taxa(taxa_sums(data_ITS_counts) > 0, data_ITS_counts)
data_ITS <- prune_samples(sample_sums(data_ITS)>=5000, data_ITS)
data_ITS <- transform_sample_counts(data_ITS, function(x) x/sum(x))

# Normalize data to 5000 reads per sample 
data_ITS_rae <- rarefy_even_depth(data_ITS_counts, sample.size = 5000)

# Extract metadata
OTU_ITS <- t(as(data_ITS_rae@otu_table, "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)
rich_ITS <- hill_taxa(OTU_ITS, q = 0)
rich_ITS <- merge(data_ITS_rae@sam_data, rich_ITS, by = 0)
shan_ITS <- hill_taxa(OTU_ITS, q = 1)
shan_ITS <- merge(data_ITS_rae@sam_data, shan_ITS, by = 0)
simp_ITS <- hill_taxa(OTU_ITS, q = 2)
simp_ITS <- merge(data_ITS_rae@sam_data, simp_ITS, by = 0)


## Alpha Diversity Figures 
# Bar Graph: Alpha Diversity (Richness)
png("Mycobiome_Richness.png", units="in", width=12, height=8, res=300)
rich <- rbind(rich_ITS %>% dplyr:: select(Species, Time_Sequential, y))
rich %>%
  filter(Time_Sequential > 0) %>%
  group_by(Species, Time_Sequential) %>%
  dplyr::summarise(mean = mean(y), se = sd(y)/sqrt(n())) %>%
  ggplot(aes(x = as.character(Time_Sequential), y = mean, fill = Species)) +
  theme(plot.title = element_text(size = 25, hjust = 0.5, vjust = 1.5)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), col = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, position = position_dodge(0.9)) +
  scale_x_discrete(labels=c("Winter '18", "Spring '18", "Summer '18", "Fall '18", "Winter '19",
                            "Spring '19", "Summer '19", "Fall '19")) +
  scale_fill_manual(values = c("#b66dff", "#b6dbff"),
                    labels = c(expression(italic("P. deltoides")),
                               c(expression(italic("P. trichocarpa"))))) +
  labs(title = "Mycobiome Richness Over Time", x = element_blank(), y = "Richness Score") 
dev.off()


# Bar Graph: Alpha Diversity (Shannon's)
png("Mycobiome_Shannon.png", units="in", width=12, height=8, res=300)
shan <- rbind(shan_ITS %>% dplyr:: select(Species, Time_Sequential, y))
shan %>% 
  filter(Time_Sequential > 0) %>%
  group_by(Species, Time_Sequential) %>%
  dplyr::summarise(mean = mean(y), se = sd(y)/sqrt(n())) %>%
  ggplot(aes(x = as.factor(Time_Sequential), y = mean, fill = Species)) +
  theme(plot.title = element_text(size = 25, hjust = 0.5, vjust = 1.5),
        legend.title = element_text(face = "bold"), legend.position = c(0.85,0.90)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), col = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, position = position_dodge(0.9)) +
  scale_x_discrete(labels=c("Winter '18", "Spring '18", "Summer '18", "Fall '18", "Winter '19",
                            "Spring '19", "Summer '19", "Fall '19")) +
  scale_fill_manual(values = c("#b66dff", "#b6dbff"),
                    labels = c(expression(italic("P. deltoides")),
                               c(expression(italic("P. trichocarpa"))))) +
  labs(title = "Mycobiome Diversity (H) Over Time", x = element_blank(), y = "Shannon Diversity Index (H)")
dev.off()


#Bar graph: Alpha Diversity (Simpson's)
png("Mycobiome_Simpson.png", units="in", width=12, height=8, res=300)
simp <- rbind(shan_ITS %>% dplyr:: select(Species, Time_Sequential, y))
simp %>%
  group_by(Species, Time_Sequential) %>%
  filter(Time_Sequential > 0) %>%
  dplyr::summarise(mean = mean(y), se = sd(y)/sqrt(n())) %>%
  ggplot(aes(x = as.factor(Time_Sequential), y = mean, fill = Species)) +
  theme(plot.title = element_text(size = 25, hjust = 0.5, vjust = 1.5)) +
  geom_bar(stat = "identity", position = position_dodge(0.9), col = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2, position = position_dodge(0.9)) +
  scale_x_discrete(labels=c("Winter '18", "Spring '18", "Summer '18", "Fall '18", "Winter '19",
                            "Spring '19", "Summer '19", "Fall '19")) +
  scale_fill_manual(values = c("#b66dff", "#b6dbff"),
                    labels = c(expression(italic("P. deltoides")),
                               c(expression(italic("P. trichocarpa"))))) +
  labs(title = "Mycobiome Diversity (D) Over Time", x = element_blank(), y = "Simpson's Diversity Index (D)")
dev.off()


# Interaction (Two-way ANOVA for Richness)
rich_anova <- aov(y~as.character(Time_Sequential)*Species,rich_ITS)
plot(rich_anova) #check that normality, constant variance, etc conditions are met
summary(rich_anova)

# Post hoc analysis (Richness)
rich_tukey <- aov(y~as.character(Time_Sequential)*Species,rich_ITS %>% filter(Time_Sequential >0)) 
plot(rich_tukey) #check that conditions are satisfied 
summary(rich_tukey)
TukeyHSD(rich_tukey, "as.character(Time_Sequential)")

# Interaction (Two-way ANOVA for Shannon's)
shan_anova <- aov(y~as.character(Time_Sequential)*Species, shan_ITS)
plot(shan_anova)
summary(shan_anova)

# Post hoc analysis (Shannon's)
shan_tukey <- aov(y~as.character(Time_Sequential)*Species,shan_ITS %>% filter(Time_Sequential >0)) 
plot(shan_tukey) #check that conditions are satisfied 
summary(shan_tukey)
TukeyHSD(shan_tukey)
TukeyHSD(shan_tukey, "as.character(Time_Sequential)")

# Interaction (Two-way ANOVA for Simpson's)
simp_anova <- aov(y~as.character(Time_Sequential)*Species, simp_ITS)
plot(simp_anova)
summary(simp_anova)

# Post hoc analysis (Simpson's)
simp_tukey <- aov(y~as.character(Time_Sequential)*Species,simp_ITS %>% filter(Time_Sequential >0)) 
plot(simp_tukey) #check that conditions are satisfied 
summary(simp_tukey)
TukeyHSD(simp_tukey)
TukeyHSD(simp_tukey, "as.character(Time_Sequential)")


## Beta Diversity Figures 
# PERMANOVA (Bray-Curtis)
# The vegan package contains all common ordination methods
library(vegan)
df_ITS <- data.frame(sample_data(data_ITS %>%
                                   subset_samples(Time_Sequential > 0)))
set.seed(12290)
OTU_ITS <- t(as(otu_table(data_ITS %>%
                            subset_samples(Time_Sequential > 0)), "matrix"))
OTU_ITS <- as.data.frame(OTU_ITS)
OTU_ITS <- OTU_ITS %>% 
  filter(`0c4e3bbd0210d8cf72949a1adab04172` >= 0) 
dist_ITS <- vegdist(OTU_ITS, "bray") 
df_ITS <- df_ITS[!(row.names(df_ITS)=="Besc103-r22-Aug2018-ITS"),]
set.seed(01221990)
mod <- adonis(dist_ITS ~ Species*as.character(Time_Sequential), data = df_ITS, permutations = 9999)


# Principal Coordinate Analysis (PCoA) Plot of Beta Diversity 
ordu <- ordinate(data_ITS %>%
                   subset_samples(Time_Sequential > 0),
                 "PCoA", "bray", weighted = T)

P <- plot_ordination(data_ITS %>%
                       subset_samples(Time_Sequential > 0),
                     ordu) 

png("Fung_PCOA.png", units="in", width=12, height=8, res=300)
PCOA_Plot <- ggplot(P$data, P$mapping) +
  ggtitle("Differences in Rhizosphere Fungal Community Composition") +
  geom_point(aes(color = as.character(Time_Sequential), shape = Species), size = 4, alpha = 0.5) +
  scale_x_continuous(name = paste0(paste0("PCoA 1 (", round(100*as.data.frame(ordu$values)[1,2], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("PCoA 2 (", round(100*as.data.frame(ordu$values)[2,2], 1)), "%)")) +
  theme(plot.title = element_text(size = 25, hjust = 0.5, vjust = 1.5),
        legend.title = element_text(face = "bold")) +
  guides(color = guide_legend(override.aes = list(pch = 16, alpha = 1), title = "Time\nSequential"),
         shape = guide_legend(title = "Species")) +
  scale_color_manual(labels = c("Winter '18", "Spring '18", "Summer '18", "Fall '18", "Winter '19",
                                "Spring '19", "Summer '19", "Fall '19"), 
                     values = c("#009292","#ff6db6","#ffff6d","#490092",
                                "#006ddb","#b66dff","#24ff24","#fd8d3c")) +
  scale_shape_manual(values = c(16,17),
                     labels = c(expression(italic("P. deltoides")), 
                                expression(italic("P. trichocarpa"))))
dev.off()

# Dispersion Analysis
df_ITS <- data.frame(sample_data(data_ITS))
df_ITS$Time_Sequential_char <- as.character(df_ITS$Time_Sequential)
dist_ITS <- phyloseq::distance(subset_samples(data_ITS), "bray")
set.seed(12290)
mod <- betadisper(d = dist_ITS, group = df_ITS$Time_Sequential_char) 
anova(mod)
TukeyHSD(mod)


# FUNGuild Analysis: Mycorrhizal Fungi
library(FSA)

#Ectomycorrhiza
ectmy <- read.csv("Ectomycorrhizal.csv")
ectmy_kruskal <- kruskal.test(Abundance ~as.character(Time_Sequential_number), ectmy %>% filter(Time_Sequential_number >0))
summary(ectmy_kruskal)

dunnTest(Abundance ~as.character(Time_Sequential_number), ectmy %>% filter(Time_Sequential_number >0))

# Bar graph: Relative Abundance Over Time
png("Ectmy_PCOA2.png", units="in", width=12, height=8, res=300)
ectmy %>%
  filter(Time_Sequential_number > 0) %>%
  group_by(as.character(Time_Sequential_number)) %>%
  dplyr::summarise(mean = mean(Abundance*100), se = sd(Abundance*100)/sqrt(n())) %>%
  
  ggplot(aes(x = `as.character(Time_Sequential_number)`, y = mean, 
             fill = `as.character(Time_Sequential_number)`)) +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 25)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  guides(color = guide_legend(override.aes = list(pch = 16, alpha = 1), title = element_blank()),
         shape = guide_legend(title = "Species")) +
  scale_x_discrete(labels = c("Winter '18", "Spring '18", "Summer '18", "Fall '18", "Winter '19",
                              "Spring '19", "Summer '19", "Fall '19")) +
  scale_fill_manual(values = c("#009292","#ff6db6","#ffff6d","#490092",
                               "#006ddb","#b66dff","#24ff24","#fd8d3c")) +
  labs(title = "Ectomycorrhizal Fungi (Beneficial Symbiont) \n Relative Abundance Over Time", 
       x = element_blank(), y = "Relative Abundance (%)")
dev.off()  

# Arbuscular Mycorrhiza 
arbmy <- data_ITS %>% 
  subset_taxa(Phylum == "Glomeromycota") %>% 
  tax_glom("Phylum") %>%
  psmelt()

arbmy_kruskal <- kruskal.test(Abundance ~as.character(Time_Sequential), arbmy %>% filter(Time_Sequential >0))
summary(arbmy_kruskal)

dunnTest(Abundance ~as.character(Time_Sequential), arbmy %>% filter(Time_Sequential >0))

## Bar graph: Relative Abundance Over Time
png("Arbmy_PCOA5.png", units="in", width=12, height=8, res=300)
arbmy %>%
  filter(Time_Sequential > 0) %>%
  group_by(as.character(Time_Sequential)) %>%
  dplyr::summarise(mean = mean(Abundance*100), se = sd(Abundance*100)/sqrt(n())) %>%
  ggplot(aes(x = `as.character(Time_Sequential)`, y = mean, 
             fill = `as.character(Time_Sequential)`))  +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5, vjust = 1.5, size = 25)) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  guides(color = guide_legend(override.aes = list(pch = 16, alpha = 1), title = element_blank()),
         shape = guide_legend(title = "Species")) +
  scale_x_discrete(labels = c("Winter '18", "Spring '18", "Summer '18", "Fall '18", "Winter '19",
                              "Spring '19", "Summer '19", "Fall '19")) +
  scale_fill_manual(values = c("#009292","#ff6db6","#ffff6d","#490092",
                               "#006ddb","#b66dff","#24ff24","#fd8d3c")) +
  labs(title = "Arbuscular Mycorrhizal Fungi (Beneficial Symbiont) n Relative Abundance Over Time", 
       x = element_blank(), y = "Relative Abundance (%)") 
dev.off()


# Fungal Taxonomic Assignment
abundant_taxa <- subset_taxa(data_ITS, Phylum == "Ascomycota" |
                               Phylum == "Basidiomycota"|
                               Phylum == "Mucoromycota"|
                               Phylum == "Mortierellomycota" |
                               Phylum == "Chytridiomycota" |
                               Phylum == "Rozellomycota" |
                               Phylum == "Entomophthoromycota" |
                               Phylum == "Glomeromycota")

colnames(tax_table(abundant_taxa))[1] <- "group"
tax_table(abundant_taxa)[,"group"] <- tax_table(abundant_taxa)[,"Phylum"]
abundant_taxa_glom <- tax_glom(abundant_taxa, taxrank = "group")
abundant_taxa_melt <- psmelt(abundant_taxa_glom)

sum2 <- abundant_taxa_melt %>% 
  group_by(Time_Sequential, sample_Species, group) %>%
  dplyr::summarise(means = mean(Abundance))

rare <- sum2 %>%
  group_by(Time_Sequential, sample_Species) %>%
  dplyr::summarise(means = 1- sum(means)) %>%
  mutate(group = "Other/Unidentified")
sum2 = rbind(sum2, rare)

sum2$group <- ordered(sum2$group, c("Ascomycota", "Basidiomycota", "Chytridiomycota", 
                                    "Entomophthoromycota", "Glomeromycota",
                                    "Mortierellomycota", "Mucoromycota", "Rozellomycota",
                                    "Other/Unidentified"))

Taxonomic.Assignment <- ggplot(sum2, aes(x = Time_Sequential, y = means, fill = group)) +
  geom_bar(position = "stack", stat = "identity", col = "black") +
  facet_wrap(~sample_Species, ncol = 3) +
  scale_fill_manual(values = c("#abd2ed", "#2ca02c", 
                               "#8c564b", "#d62728", "#9467bd", "yellow", "#e377c2", 
                               "#17becf", "#7f7f7f", "white")) +
  scale_y_continuous(name = NULL, breaks = NULL) +
  scale_x_discrete(name = "Time Sequential") +
  panel_border() +
  theme(strip.text = element_text(size = 10), axis.text.x = element_text(size = 10),
        axis.title.y = element_blank(), legend.position = "bottom", legend.title = element_blank(),
        legend.justification = "center") +
  guides(fill = guide_legend(ncol = 3)) 

png("Fungal Taxonomic Assignment.png", res = 300, height = 8, width = 12 , units = "in")
Taxonomic.Assignment
dev.off()
