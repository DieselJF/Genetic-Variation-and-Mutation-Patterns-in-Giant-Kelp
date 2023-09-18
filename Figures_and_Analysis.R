# R script to make all figures from paper. 
# Includes details on what input files to use
# These files can be easily made from Mutation Table

library(ggplot2)
library(dplyr)
library(tidyverse)
library(tidyr)
library(patchwork)
library(patchwork)
library(cowplot)

########
# Figure 1
########

mut <- read.table("ind_mut_for_R.txt", header = TRUE) # File containing the sum of all mutation in an individuals. Columns: 1- SampleID, 2- #Modifier mutations, 3- #Low impact mutations, 4- # Moderate mutations, 5- # High impact mutations

mut <- mut %>%
  pivot_longer(cols = c("MODIFIER", "LOW", "MODERATE", "HIGH"),
               names_to = "Mut_type",
               values_to = "Mut_number")

theme_set(theme_minimal(base_size = 14, base_family = "Arial"))
mut %>%
  filter(Mut_type != "MODIFIER") %>%
  ggplot(aes(Mut_number, fill = Mut_type)) +
  geom_histogram(position = "identity", alpha = 0.7) +
  theme(legend.position = c(0.5,0.552)) +
  ylab("Individuals") +
  xlab("# Mutations") +
  guides(fill=guide_legend(title="Mutation Impact"))

########
# Figure 2 - CI more high
########

# Input file = file containing mutations per gene, by individual. 
# Columns = 1- SampleID, 2- Gene, 3- #Modifier mutations, 4- #Low impact mutations, 5- #Moderate impact mutations, 6- #High impact mutations
mut_file <- read.table("IND_GENE_MODI_LOW_MODE_HIGH_ALL.txt", sep = "\t", header = TRUE)

# Phenotyping file provided by Osborn et al. 
mass <- read.table("tab_fixed_effects_and_phenotypes.tsv", header = T)
mass <- distinct(mass[,c("VCF_ID", "Population")])

mut_file <- merge(mut_file, mass, by = "VCF_ID")

theme_set(theme_minimal(base_size = 14, base_family = "Arial"))
mut_file %>%
  group_by(Population, VCF_ID) %>%
  summarize(sum_high = sum(ALL_HIGH) / (sum(ALL_LOW) + sum(ALL_MODERATE) + sum(ALL_HIGH))) %>%
  ggplot(aes(x = Population, y = sum_high)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5, color = "darkgray") +
  labs(x = "Population", y = "Ratio of High impact mutations") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "none"
  )

#############
# Figure 3 A and B - Size by annotation 
#############

# Input file = file containing mutations per gene, by individual. 
# Columns = 1- SampleID, 2- Gene, 3- #Modifier mutations, 4- #Low impact mutations, 5- #Moderate impact mutations, 6- #High impact mutations
# Same mut file from figure 2
mut_file <- read.table("IND_GENE_MODI_LOW_MODE_HIGH_ALL.txt", sep = "\t", header = TRUE)

kog_gene_file <- read.table("Function_gene_length.txt", sep = "\t") # Modified KOG gene pathway file. Columns = 1- KOG pathway annotated for a gene, 2- GeneID, 3- Gene length
colnames(kog_gene_file) <- c("kog", "gene", "Length")
kog_gene_file <- select(kog_gene_file, c("kog", "gene"))
kog_gene_file$kog <- "Annotated"
kog_gene_file <- distinct(kog_gene_file)


length_protID <- read.table("Gene_position.txt") # File for gene positions. Columns = 1- GeneID, 2- Scaffold, 3- Start nucleotide, 4- End nucleotide
colnames(length_protID) <- c("gene", "scaffold", "start", "end")
length_protID$gene_length <- length_protID$end - length_protID$start
length_protID <- select(length_protID, c("gene", "gene_length"))

rj_table <- right_join(kog_gene_file, length_protID)
rj_table[is.na(rj_table)] <- "Not annotated"


kog_length <- rj_table %>%
  group_by(kog) %>%
  summarize(total_func_size = sum(gene_length))

colnames(rj_table) <- c("kog", "GENE", "gene_length")

full_data <- merge(mut_file, rj_table, by = "GENE")
full_data <- merge(full_data, kog_length, by = "kog")

rj_table %>%
  ggplot(aes(x = kog, y = gene_length)) +
  geom_boxplot()

theme_set(theme_minimal(base_size = 14, base_family = "Arial"))

# Plot B
full_data %>%
  group_by(kog, VCF_ID) %>%
  summarize(ind_high_by_len = sum(ALL_HIGH)/total_func_size) %>% 
  distinct() %>%
  ggplot(aes(kog, ind_high_by_len)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  #geom_jitter(width = 0.2, alpha = 0.3, color = "darkgray")+
  ggtitle("B") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"
  ) +
  labs(x = "", y = "High impact mutation / pathway length")

# Plot A
theme_set(theme_minimal(base_size = 14, base_family = "Arial"))
rj_table %>%
  ggplot(aes(kog, gene_length)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  ggtitle("A") +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "none"
  ) +
  labs(x = "", y = "Gene length")

##############
# Figure 4 - Relative gene position
##############

# Input file containing gene of impact, relative gene position for a gene of 1000 nucleotides, and impact of mutation.
# Columns = 1- GeneID, 2- Effect, 3- Relative gene position, 4- Impact
data <- read.table("high_low_where.txt")
colnames(data) <- c("gene", "effect", "position", "impact")

# File containing strand information to correct for relative gene position 
# Columns = 1- GeneID, 2- Strand (+/-)
strand_data <- read.table("gene_strand.txt")
colnames(strand_data) <- c("gene", "strand")

data <- merge(data, strand_data, by = "gene")

data$relativePosition <- ifelse(data$strand == "+", data$position,  1000 - data$position)

kog_gene_file <- read.table("Function_gene_length.txt", sep = "\t")
colnames(kog_gene_file) <- c("kog", "gene", "Length")
kog_gene_file$annot <- "Annotated"
kog_gene_file <- select(kog_gene_file, c("annot", "gene"))
kog_gene_file <- distinct(kog_gene_file)


data <- right_join(kog_gene_file, data, by = "gene")
data[is.na(data)] <- "Not annotated"

data %>%
  ggplot(aes(relativePosition, fill = annot)) +
  geom_histogram(data = . %>% filter(annot == "Not annotated"), bins = 30, position = "identity") +
  geom_histogram(data = . %>% filter(annot == "Annotated"), bins = 30, position = "identity") +
  scale_fill_manual(values = c("gray77", "grey56")) +
  labs(x = "Relative gene position", y = "# Mutations") +
  labs(fill = "") +
  theme(legend.position = c(0.4,0.635))

#########
# Figure 5 A and B - Density of mutations & Ratios
#########

# Same file from Figure 2
mut_file <- read.table("IND_GENE_MODI_LOW_MODE_HIGH_ALL.txt", sep = "\t", header = TRUE)

# Same file from Figure 3
kog_gene_file <- read.table("Function_gene_length.txt", sep = "\t")
colnames(kog_gene_file) <- c("kog", "gene", "Length")
kog_gene_file <- select(kog_gene_file, c("kog", "gene"))

# File with gene positions
# Columns = 1- GenID, 2- Scaffold, 3- Gene start position, 4- Gene end position
length_protID <- read.table("Gene_position.txt")
colnames(length_protID) <- c("gene", "scaffold", "start", "end")
length_protID$gene_length <- length_protID$end - length_protID$start
length_protID <- select(length_protID, c("gene", "gene_length"))

gene_kog_length <- merge(kog_gene_file, length_protID, by = "gene")

length_protID$kog <- "All_genes"
length_protID <- length_protID[c("gene", "kog", "gene_length")]

gene_kog_length <- rbind(gene_kog_length, length_protID)

rm(length_protID, kog_gene_file)

kog_length <- gene_kog_length %>%
  group_by(kog) %>%
  summarize(total_func_size = sum(gene_length))

colnames(gene_kog_length) <- c("GENE", "kog", "gene_length")

full_data <- merge(mut_file, gene_kog_length, by = "GENE")

full_data <- merge(full_data, kog_length, by = "kog")

# #Mutations / total size
theme_set(theme_minimal(base_size = 14, base_family = "Arial"))
plot_A <- full_data %>%
  group_by(VCF_ID, kog) %>%
  summarise(mut_rate = (sum(ALL_HIGH) + sum(ALL_LOW) + sum(ALL_MODERATE))/total_func_size) %>%
  unique() %>%
  ggplot(aes(kog, mut_rate)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1)) +
  theme(axis.text=element_text(size=12)) +
  ylab("Mutations\nper nucleotide") +
  theme(axis.text.x = element_blank()) +
  xlab("")

# dN/dS copy
theme_set(theme_minimal(base_size = 14, base_family = "Arial"))
plot_B <- full_data %>%
  group_by(VCF_ID, kog) %>%
  summarise(dn_ds_exp = (sum(ALL_HIGH))/sum(ALL_LOW)) %>%
  unique() %>%
  ggplot(aes(kog, dn_ds_exp)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  ylim(c(0,0.3)) +
  theme(axis.text=element_text(size=12),
        plot.margin = margin(5, 5, 5, 50)) +
  ylab("dN/dS Proxy")+
  xlab("")

final_plot <- plot_A + plot_B + plot_layout(nrow=2)
ggsave("./final_plot_fig5.png", final_plot, width = 10, height = 10, bg = "white")

#############
# Figure 6 - TjD, Pi, Fst
#############

# Same from Figure 3
kog <- read.table("Function_gene_length.txt", sep = "\t")
colnames(kog) <- c("Function", "Gene", "Length")
kog <- select(kog, -Length)

# File generated after running "Statistics_pop_gen_per_gene"
pop_gen <- read.table("gene_pop_gen_stats.txt")
colnames(pop_gen) <- c("Gene", "scaffold", "start", "end", "length", "Fst", "Pi", "TjD")

final_table <- right_join(pop_gen, kog, by="Gene")
pop_gen$Function <- as.factor("Total")
final_table <- rbind(pop_gen, final_table)

theme_set(theme_minimal(base_size = 14, base_family = "Arial"))

# Pi plot
plot_A <- final_table %>%
  ggplot(aes(x = Function, y = Pi)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  labs(x = "Pi") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = NA)
  )+
  xlab("") +
  theme(axis.text.x = element_blank())

# TjD plot
plot_B <- final_table %>%
  ggplot(aes(x = Function, y = TjD)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  labs(y = "Tajima's D") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = NA)
  )+
  xlab("") +
  theme(axis.text.x = element_blank())

# Fst plot
plot_C <-final_table %>%
  ggplot(aes(x = Function, y = Fst)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  labs(x = "Fst") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1),  # Rotate x-axis labels
    legend.position = "none",
    panel.background = element_rect(fill = "white", colour = NA),
    plot.margin = margin(5, 5, 5, 50)
  ) +
  ylim(c(0,1))+
  xlab("")

final_plot <- plot_A + plot_B + plot_C + plot_layout(nrow=3)
ggsave("./final_plot_fig6.png", final_plot, width = 10, height = 10, bg = "white")

#############
# Figure 7 A through C - Biomass, Carbon, Nitrogen vs # high mutation
#############

# Same file from Figure 1
mut <- read.table("ind_mut_for_R.txt", header = T)

# Same file from figure 2, provided by Osborn et al.
mass <- read.table("tab_fixed_effects_and_phenotypes.tsv", header = T)
mass <- mass[,c("VCF_ID", "Population", "Biomass")]
colnames(mass) <- c("Sample", "Population", "Biomass")

mass <- na.omit(mass) %>%
  group_by(Sample, Population) %>%
  summarize(mean_B = mean(Biomass))

data <- merge(mut, mass, by = "Sample")

theme_set(theme_minimal(base_size = 14, base_family = "Arial"))
plot_A <- data %>%
  ggplot(aes(HIGH, mean_B)) +
  geom_point()+
  theme_minimal() +
  stat_smooth(method="lm", se = TRUE, color = "grey50") +
  xlab("# High impact mutations") +
  ylab("Mean Biomass") +
  xlab("")

# Carbon graph B
mut <- read.table("ind_mut_for_R.txt", header = T)

mass <- read.table("tab_fixed_effects_and_phenotypes.tsv", header = T)
mass <- mass[,c("VCF_ID", "Population", "Carbon")]
colnames(mass) <- c("Sample", "Population", "Carbon")

mass <- na.omit(mass) %>%
  group_by(Sample, Population) %>%
  summarize(mean_C = mean(Carbon))

data <- merge(mut, mass, by = "Sample")

theme_set(theme_minimal(base_size = 14, base_family = "Arial"))
plot_B <- data %>%
  ggplot(aes(HIGH, mean_C)) +
  geom_point()+
  theme_minimal() +
  stat_smooth(method="lm", se = TRUE, color = "grey50") +
  xlab("# High impact mutations") +
  ylab("Mean Carbon") +
  ylim(c(0,27)) +
  xlab("")

# Nitrogen graph C
mut <- read.table("ind_mut_for_R.txt", header = T)

mass <- read.table("tab_fixed_effects_and_phenotypes.tsv", header = T)
mass <- mass[,c("VCF_ID", "Population", "Nitrogen")]
colnames(mass) <- c("Sample", "Population", "Nitrogen")

mass <- na.omit(mass) %>%
  group_by(Sample, Population) %>%
  summarize(mean_N = mean(Nitrogen))

data <- merge(mut, mass, by = "Sample")

theme_set(theme_minimal(base_size = 14, base_family = "Arial"))
plot_C <- data %>%
  ggplot(aes(HIGH, mean_N)) +
  geom_point()+
  theme_minimal() +
  stat_smooth(method="lm", se = TRUE, color = "grey50") +
  xlab("# High impact mutations") +
  ylab("Mean Nitrogen")

final_plot <- plot_A + plot_B + plot_C + plot_layout(nrow=3)
ggsave("./final_plot_fig7.png", final_plot, width = 10, height = 10, bg = "white")

#############
# Fig 8 - AMD deaminase mutations
#############

# Same file from Figure 2
mut_file <- read.table("IND_GENE_MODI_LOW_MODE_HIGH_ALL.txt", sep = "\t", header = TRUE)

# Same file from Figure 2, provided by Osborne et al.
mass <- read.table("tab_fixed_effects_and_phenotypes.tsv", header = T)
mass <- mass[,c("VCF_ID", "Population", "Biomass", "Carbon")]
mass <- select(mass, c("VCF_ID", "Population", "Biomass", "Carbon"))
mass <- na.omit(mass) %>%
  group_by(VCF_ID) %>%
  summarize(mean_B = mean(Biomass), mean_C = mean(Carbon), Population) %>%
  distinct()

mut_file <- merge(mut_file, mass, by = "VCF_ID")

theme_set(theme_minimal(base_size = 14, base_family = "Arial"))
plot <- mut_file %>% 
  filter(GENE == "9620408") %>%
  mutate(HIGH = as.factor(ALL_HIGH)) %>%
  ggplot(aes(HIGH, mean_B)) + 
  geom_boxplot() +
  theme_minimal() +
  geom_point(aes(2,273), color = "red") + # ind 175
  geom_point(aes(2,282), color = "blue") + # ind 82
  ylab("Mean Biomass")  + 
  annotate("text", x = 2.15, y = 273, label = "Arg527His", color = "red")+ 
  annotate("text", x = 1.85, y = 282, label = "Gly11Glu", color = "blue") +
  xlab("") +
  scale_x_discrete(labels=c("Splice acceptor variant","wild type"))

ggsave("./final_plot_fig8.png", plot, width = 7, height = 6, bg = "white")
















