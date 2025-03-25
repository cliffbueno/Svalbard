# Network analysis for Svalbard Patches for Steve Schmidt
# By Cliff Bueno de Mesquita, Spring 2024.



#### 1. Overview ####
# Analysis Outline
## 1. This overview with background information
## 2. Setup (libraries, data import)
## 3. SPIEC-EASI (all data)
## 4. Spot Only analysis

# Data come from 2 different patches from the Midtre Lovenbreen forefield
# n = 12 samples from each patch (24 total)
# At each patch are 6 biological soil crust "spots" and 6 controls "bare"
# Control samples were always taken 2 cm away from spots
# 16S and 18S sequencing on Illumina MiSeq was performed in the Fierer Lab
# Data processed by Adam Solon with dada2 pipeline



#### 2. Setup ####
# Copied from other scripts; some but not all of these are needed
library(plyr) # Data manipulation
library(tidyverse) # Data manipulation
library(mctoolsr) # Microbial analyses
library(RColorBrewer) # Colors
library(vegan) # Multivariate analyses
library(indicspecies) # Indicator species
library(car) # Stats
library(FSA) # SE
library(magrittr) # Set names
library(PMCMRplus) # Stats
library(readxl) # Excel
library(writexl) # Excel
library(plotly) # Interactive plots
library(ggmap) # Maps
library(ggsn) # Maps
library(multcomp) # Tukey HSD and significance letters
library(emmeans) # Tukey HSD and significance letters
library(scales) # View colors
library(cowplot) # Multipanels
library(qvalue) # q values for indicator species
library(reshape2) # melt
library(gridExtra) # graphs
library(grid) # graphs
library(cowplot) # graphs
library(ggpubr) # graphs
library(ggExtra) # graphs
library(ggh4x) # graphs
library(dendextend) # graphs
library(corrplot) # correlation plots
library(pheatmap) # heatmaps
library(zCompositions) # CLR
library(compositions) # Aitchison
library(mobr) # rarefaction curves
library(plotly) # interactive graphs
library(pairwiseAdonis) # pairwise permanova
library(patchwork) # insets
library(ggbiplot) # ordinations
library(ggfortify) # ordinations
library(ggrepel) # repel text
library(usmap) # for map
library(SpiecEasi) # for networks
library(igraph) # for plotting networks
library(Matrix) # for networks
library(phyloseq) # if you want for networks
library(brainGraph) # for attack robustness (stability)
library(rnetcarto) # module identification
library(ggrepel) # labels
library(plotly) # Interactive graphs
library(networkD3) # Interactive networks
library(IgAScores) # For relative abundance

# Functions
find_hull <- function(df) df[chull(df$Axis01, df$Axis02),]
`%notin%` <- Negate(`%in%`)

# Repository path
setwd("~/Desktop/Svalbard")

# Shape functions
MyDiamond <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.frame.color <- params("vertex", "frame.color")
  if (length(vertex.frame.color) != 1 && !is.null(v)) {
    vertex.frame.color <- vertex.frame.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color, fg=vertex.frame.color,
          stars=cbind(vertex.size, vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
add_shape("diamond", clip=shapes("circle")$clip,
          plot=MyDiamond, parameters=list(vertex.frame.color="white",
                                          vertex.frame.width=1))

mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)

# Import and prepare data
tax_table_fp <- "~/Desktop/Svalbard/16S_seqtab_wTax_mctoolsr.txt"
map_fp <- "~/Desktop/Svalbard/Svalbard_Mapping.txt"
input_16S = load_taxa_table(tax_table_fp, map_fp) # 24 samples loaded

tax_table_fp <- "~/Desktop/Svalbard/18S_seqtab_wTax_mctoolsr.txt"
map_fp <- "~/Desktop/Svalbard/Svalbard_Mapping.txt"
input_18S = load_taxa_table(tax_table_fp, map_fp) # 24 samples loaded


# Filter data
input_16S_filt <- input_16S
input_16S_filt <- filter_taxa_from_input(input_16S_filt,
                                         taxa_to_remove = "Chloroplast") # 46 removed
input_16S_filt <- filter_taxa_from_input(input_16S_filt,
                                         taxa_to_remove = "Mitochondria") # 65 removed
input_16S_filt <- filter_taxa_from_input(input_16S_filt,
                                         taxa_to_remove = "Eukaryota") # 0 removed
input_16S_filt <- filter_taxa_from_input(input_16S_filt,
                                         taxa_to_remove = "NA",
                                         at_spec_level = 1) # 1 removed
singdoub_16S <- data.frame("count" = rowSums(input_16S_filt$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.))
input_16S_filt <- filter_taxa_from_input(input_16S_filt,
                                         taxa_IDs_to_remove = singdoub_16S$ASV) # 28 removed
nrow(input_16S_filt$data_loaded) # 2263

input_18S_filt <- input_18S
input_18S_filt <- filter_taxa_from_input(input_18S_filt,
                                         taxa_to_remove = "NA",
                                         at_spec_level = 1) # 1 removed
singdoub_18S <- data.frame("count" = rowSums(input_18S_filt$data_loaded)) %>%
  filter(count < 3) %>%
  mutate(ASV = rownames(.))
input_18S_filt <- filter_taxa_from_input(input_18S_filt,
                                         taxa_IDs_to_remove = singdoub_18S$ASV) # 3 removed
nrow(input_18S_filt$data_loaded) # 389

# Add sampleID
input_16S_filt$map_loaded$sampleID <- rownames(input_16S_filt$map_loaded)
input_18S_filt$map_loaded$sampleID <- rownames(input_18S_filt$map_loaded)

# Check sampleID match
sum(input_16S_filt$map_loaded$sampleID != input_18S_filt$map_loaded$sampleID) # 0, good.
sum(names(input_16S_filt$data_loaded) != names(input_18S$data_loaded)) # 0, good

# Add ASV ID to taxonomy table
input_16S_filt$taxonomy_loaded$taxonomy7 <- rownames(input_16S_filt$taxonomy_loaded)
input_18S_filt$taxonomy_loaded$taxonomy7 <- rownames(input_18S_filt$taxonomy_loaded)

# Need to create unique prok and euk ESV (ASV) IDs
# Unique ASV IDs are needed for rownames of $data_loaded and $taxonomy_loaded
input_16S_filt_toMerge <- input_16S_filt
num_16S <- seq(1:nrow(input_16S_filt_toMerge$data_loaded))
ASVlabel_16S <- rep("ASV_Prok_", nrow(input_16S_filt_toMerge$data_loaded))
IDs_16S <- paste(ASVlabel_16S, num_16S, sep = "")
rownames(input_16S_filt_toMerge$data_loaded) <- IDs_16S
rownames(input_16S_filt_toMerge$taxonomy_loaded) <- IDs_16S
input_16S_filt_toMerge$taxonomy_loaded$taxonomy7 <- IDs_16S

input_18S_filt_toMerge <- input_18S_filt
num_18S <- seq(1:nrow(input_18S_filt_toMerge$data_loaded))
ASVlabel_18S <- rep("ASV_Euk_", nrow(input_18S_filt_toMerge$data_loaded))
IDs_18S <- paste(ASVlabel_18S, num_18S, sep = "")
rownames(input_18S_filt_toMerge$data_loaded) <- IDs_18S
rownames(input_18S_filt_toMerge$taxonomy_loaded) <- IDs_18S
input_18S_filt_toMerge$taxonomy_loaded$taxonomy7 <- IDs_18S

input_filt_combined <- input_16S_filt_toMerge
input_filt_combined$data_loaded <- rbind(input_filt_combined$data_loaded,
                                         input_18S_filt_toMerge$data_loaded)
input_filt_combined$taxonomy_loaded <- rbind(input_filt_combined$taxonomy_loaded,
                                             input_18S_filt_toMerge$taxonomy_loaded)

nrow(input_filt_combined$data_loaded) # 2652 total taxa


### Prior to SPIEC-EASI, explore the data a little bit
# Not for publication, just for my own curiosity

### Venn Diagrams
# Check for ASV overlap in spots and controls for each patch
input_16S_filt$map_loaded <- input_16S_filt$map_loaded %>%
  mutate(Spot = substr(Replicate, start = 2, stop = 2)) %>%
  mutate(Spot = gsub("A", "Spot", Spot)) %>%
  mutate(Spot = gsub("B", "Bare", Spot)) %>%
  mutate(PatchName = gsub(1, "Patch1", Patch)) %>%
  mutate(PatchName = gsub(2, "Patch2", PatchName)) %>%
  mutate(PatchSpot = paste(PatchName, Spot, sep = "_"))
png("Venn16S.png", width = 7, height = 5, units = "in", res = 300)
plot_venn_diagram(input_16S_filt, "PatchSpot", 0.0000000000000000000000000001)
dev.off()

input_18S_filt$map_loaded <- input_18S_filt$map_loaded %>%
  mutate(Spot = substr(Replicate, start = 2, stop = 2)) %>%
  mutate(Spot = gsub("A", "Spot", Spot)) %>%
  mutate(Spot = gsub("B", "Bare", Spot)) %>%
  mutate(PatchName = gsub(1, "Patch1", Patch)) %>%
  mutate(PatchName = gsub(2, "Patch2", PatchName)) %>%
  mutate(PatchSpot = paste(PatchName, Spot, sep = "_"))
png("Venn18S.png", width = 7, height = 5, units = "in", res = 300)
plot_venn_diagram(input_18S_filt, "PatchSpot", 0.0000000000000000000000000001)
dev.off()

### Alpha
input_16S_filt$map_loaded$rich <- specnumber(input_16S_filt$data_loaded, 
                                             MARGIN = 2)
input_16S_filt$map_loaded$shannon <- vegan::diversity(input_16S_filt$data_loaded, 
                                                      index = "shannon", 
                                                      MARGIN = 2)
ggplot(input_16S_filt$map_loaded, aes(PatchSpot, rich, colour = PatchSpot)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.5) +
  labs(x = "Patch and Spot", y = "ASV Richness", colour = "PatchSpot") +
  scale_colour_manual(values = c("red1", "red3", "blue1", "blue4")) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14),
        legend.position = "none")

input_18S_filt$map_loaded$rich <- specnumber(input_18S_filt$data_loaded, 
                                             MARGIN = 2)
input_18S_filt$map_loaded$shannon <- vegan::diversity(input_18S_filt$data_loaded, 
                                                      index = "shannon", 
                                                      MARGIN = 2)
ggplot(input_18S_filt$map_loaded, aes(PatchSpot, rich, colour = PatchSpot)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(size = 3, alpha = 0.5) +
  labs(x = "Patch and Spot", y = "ASV Richness", colour = "PatchSpot") +
  scale_colour_manual(values = c("red1", "red3", "blue1", "blue4")) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16), 
        axis.text = element_text(size = 14),
        legend.position = "none")

### Beta
# First rarefy
sort(colSums(input_16S_filt$data_loaded))
set.seed(503)
input_16S_filt_rar = single_rarefy(input_16S_filt, 4129)

# Don't rarefy 18S, control samples (bare) are too low. Use relative abundance
sort(colSums(input_18S_filt$data_loaded))
input_18S_filt_rel <- input_18S_filt
input_18S_filt_rel$data_loaded <- relabund(input_18S_filt_rel$data_loaded)

# Quick ordination
dm_16S <- calc_dm(input_16S_filt_rar$data_loaded)
ord_16S <- calc_ordination(dm_16S, 'PCoA')
mctoolsr::plot_ordination(input_16S_filt_rar, ord_16S, 'PatchSpot', hulls = TRUE) +
  scale_colour_manual(values = c("red1", "red3", "blue1", "blue4")) +
  scale_fill_manual(values = c("red1", "red3", "blue1", "blue4"))

dm_18S <- calc_dm(input_18S_filt_rel$data_loaded)
ord_18S <- calc_ordination(dm_18S, 'PCoA')
mctoolsr::plot_ordination(input_18S_filt_rel, ord_18S, 'PatchSpot', hulls = TRUE) +
  scale_colour_manual(values = c("red1", "red3", "blue1", "blue4")) +
  scale_fill_manual(values = c("red1", "red3", "blue1", "blue4"))

# Barplots
tax_sum_16S_gen <- summarize_taxonomy(input_16S_filt_rar, level = 6, relative = T, report_higher_tax = F)
bars_gen_16S <- plot_taxa_bars(tax_sum_16S_gen,
                               input_16S_filt_rar$map_loaded,
                               "sampleID",
                               num_taxa = 12,
                               data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_16S_filt_rar$map_loaded, by = c("group_by" = "sampleID"))
top_gen_16S <- bars_gen_16S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_gen_16S <- bars_gen_16S %>%
  mutate(taxon = factor(taxon,
                        levels = c("NA", "Other", rev(top_gen_16S$taxon))))
ggplot(bars_gen_16S, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Sample", y = "Relative abundance (%)", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_nested_wrap(~ Spot + PatchName, scales = "free_x", ncol = 4) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2))

tax_sum_18S_gen <- summarize_taxonomy(input_18S_filt_rel, level = 6, report_higher_tax = F)
bars_gen_18S <- plot_taxa_bars(tax_sum_18S_gen,
                               input_18S_filt_rel$map_loaded,
                               "sampleID",
                               num_taxa = 12,
                               data_only = TRUE) %>%
  mutate(taxon = fct_rev(taxon)) %>%
  left_join(., input_18S_filt_rel$map_loaded, by = c("group_by" = "sampleID"))
top_gen_18S <- bars_gen_18S %>%
  group_by(taxon) %>%
  summarise(mean = mean(mean_value)) %>%
  filter(taxon != "Other") %>%
  filter(taxon != "NA") %>%
  arrange(-mean) %>%
  mutate(taxon = as.character(taxon))
bars_gen_18S <- bars_gen_18S %>%
  mutate(taxon = factor(taxon,
                        levels = c("NA", "Other", rev(top_gen_18S$taxon))))
ggplot(bars_gen_18S, aes(group_by, mean_value, fill = taxon)) +
  geom_bar(stat = "identity", colour = NA, linewidth = 0.25) +
  labs(x = "Sample", y = "Relative abundance (%)", fill = "Genus") +
  scale_fill_manual(values = c("grey75", "grey90", brewer.pal(12, "Paired")[12:1])) +
  scale_y_continuous(expand = c(0.01, 0.01)) + 
  facet_nested_wrap(~ Spot + PatchName, scales = "free_x", ncol = 4) +
  theme_classic() +
  theme(axis.title.y = element_text(face = "bold", size = 12),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        strip.text = element_text(size = 4),
        strip.background = element_rect(linewidth = 0.2))



#### 3. SPIEC-EASI ####
# Normalization is done internally, so use raw counts
# That will then be clr transformed in SPIEC-EASI
# Methods are glasso or mb (Meinshausen-Buhlmann's neighborhood selection)
# Follow the basic usage online
se.mb <- spiec.easi(as.matrix(input_filt_combined$data_loaded), 
                    method='mb', 
                    lambda.min.ratio=1e-2,
                    nlambda=20, 
                    pulsar.params=list(rep.num=50))
se.gl <- spiec.easi(as.matrix(input_filt_combined$data_loaded), 
                    method='glasso', 
                    lambda.min.ratio=1e-2,
                    nlambda=20, 
                    pulsar.params=list(rep.num=50))
sparcc <- sparcc(as.matrix(input_filt_combined$data_loaded))

## Define arbitrary threshold for SparCC correlation matrix for the graph
sparcc.graph <- abs(sparcc$Cor) >= 0.5
diag(sparcc.graph) <- 0
sparcc.graph <- Matrix(sparcc.graph, sparse=TRUE)

## Create igraph objects
ig.mb     <- adj2igraph(getRefit(se.mb))
ig.gl     <- adj2igraph(getRefit(se.gl))
ig.sparcc <- adj2igraph(sparcc.graph)

## set size of vertex proportional to clr-mean
vsize    <- rowMeans(clr(as.matrix(input_filt_combined$data_loaded), 1))+6
am.coord <- layout.fruchterman.reingold(ig.mb)

par(mfrow=c(1,3))
plot(ig.mb, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="MB")
plot(ig.gl, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="glasso")
plot(ig.sparcc, layout=am.coord, vertex.size=vsize, vertex.label=NA, main="sparcc")

# Edges
secor  <- cov2cor(getOptCov(se.gl))
sebeta <- symBeta(getOptBeta(se.mb), mode='maxabs')
elist.gl     <- summary(triu(secor*getRefit(se.gl), k=1))
elist.mb     <- summary(sebeta)
elist.sparcc <- summary(sparcc.graph*sparcc$Cor)

par(mfrow=c(1,1))
hist(elist.sparcc[,1], main='', xlab='edge weights')
hist(elist.mb[,3], col='forestgreen')
hist(elist.gl[,3], add=TRUE, col='red')

# Degree distributions
dd.gl     <- degree.distribution(ig.gl)
dd.mb     <- degree.distribution(ig.mb)
dd.sparcc <- degree.distribution(ig.sparcc)

plot(0:(length(dd.sparcc)-1), dd.sparcc, ylim=c(0,2), type='b',
     ylab="Frequency", xlab="Degree", main="Degree Distributions")
points(0:(length(dd.gl)-1), dd.gl, col="red" , type='b')
points(0:(length(dd.mb)-1), dd.mb, col="forestgreen", type='b')
legend("topright", c("MB", "glasso", "sparcc"),
       col=c("forestgreen", "red", "black"), pch=1, lty=1)

# glasso and mb were the same in this case
se.gl
se.gl$est
se.gl$select

# Good igraph tutorial here https://kateto.net/netscix2016.html
ig.gl
ig.gl[]
E(ig.gl) # 30 edges
V(ig.gl) # 39 vertices
edge_attr(ig.gl)
vertex_attr(ig.gl)
graph_attr(ig.gl)
edge_density(ig.gl, loops=F)
ecount(ig.gl)/(vcount(ig.gl)*(vcount(ig.gl)-1))
deg <- degree(ig.gl, mode="all")
plot(ig.gl, vertex.size=deg*3)
hist(deg, breaks=1:vcount(ig.gl)-1, main="Histogram of node degree")
deg.dist <- degree_distribution(ig.gl, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", 
      xlab="Degree", ylab="Cumulative Frequency")
degree(ig.gl, mode="in")
centr_degree(ig.gl, mode="in", normalized=T)
closeness(ig.gl, mode="all", weights=NA) 
centr_clo(ig.gl, mode="all", normalized=T) 
eigen_centrality(ig.gl, directed=T, weights=NA)
centr_eigen(ig.gl, directed=T, normalized=T) 
betweenness(ig.gl, directed=T, weights=NA)
edge_betweenness(ig.gl, directed=T, weights=NA)
centr_betw(ig.gl, directed=T, normalized=T)



#### _Abundant ####
# Filter to most abundant and prevalent taxa
# However, this needs to be done separately to get the top proks and top euks
nrow(input_filt_combined$data_loaded) # 2652
nrow(input_16S_filt$data_loaded) # 2263
nrow(input_18S_filt$data_loaded) # 389

# Let's try getting the top 100 of each, then merging.
top_16S <- as.data.frame(rowMeans(input_16S_filt$data_loaded)) %>%
  set_names("MeanCount") %>%
  arrange(desc(MeanCount)) %>%
  slice_head(n = 100)
input_16S_filt_abund <- filter_taxa_from_input(input_16S_filt,
                                               taxa_IDs_to_keep = rownames(top_16S),
                                               at_spec_level = 7)
nrow(input_16S_filt_abund$data_loaded) # 100

top_18S <- as.data.frame(rowMeans(input_18S_filt$data_loaded)) %>%
  set_names("MeanCount") %>%
  arrange(desc(MeanCount)) %>%
  slice_head(n = 100)
input_18S_filt_abund <- filter_taxa_from_input(input_18S_filt,
                                               taxa_IDs_to_keep = rownames(top_18S),
                                               at_spec_level = 7) # 46 removed
nrow(input_18S_filt_abund$data_loaded) # 100

# Need to create unique prok and euk ESV (ASV) IDs
# Unique ASV IDs are needed for rownames of $data_loaded and $taxonomy_loaded
input_16S_filt_abund_toMerge <- input_16S_filt_abund
num_16S <- seq(1:nrow(input_16S_filt_abund_toMerge$data_loaded))
ASVlabel_16S <- rep("ASV_Prok_", nrow(input_16S_filt_abund_toMerge$data_loaded))
IDs_16S <- paste(ASVlabel_16S, num_16S, sep = "")
rownames(input_16S_filt_abund_toMerge$data_loaded) <- IDs_16S
rownames(input_16S_filt_abund_toMerge$taxonomy_loaded) <- IDs_16S
input_16S_filt_abund_toMerge$taxonomy_loaded$taxonomy8 <- IDs_16S

input_18S_filt_abund_toMerge <- input_18S_filt_abund
num_18S <- seq(1:nrow(input_18S_filt_abund_toMerge$data_loaded))
ASVlabel_18S <- rep("ASV_Euk_", nrow(input_18S_filt_abund_toMerge$data_loaded))
IDs_18S <- paste(ASVlabel_18S, num_18S, sep = "")
rownames(input_18S_filt_abund_toMerge$data_loaded) <- IDs_18S
rownames(input_18S_filt_abund_toMerge$taxonomy_loaded) <- IDs_18S
input_18S_filt_abund_toMerge$taxonomy_loaded$taxonomy8 <- IDs_18S

input_filt_abund_combined <- input_16S_filt_abund_toMerge
input_filt_abund_combined$data_loaded <- rbind(input_filt_abund_combined$data_loaded,
                                         input_18S_filt_abund_toMerge$data_loaded)
input_filt_abund_combined$taxonomy_loaded <- rbind(input_filt_abund_combined$taxonomy_loaded,
                                             input_18S_filt_abund_toMerge$taxonomy_loaded)

nrow(input_filt_abund_combined$data_loaded) # 200 total taxa



#### _phyloseq ####
# Convert mctoolsr to phyloseq
names(input_filt_abund_combined$taxonomy_loaded) <- c("Domain", "Phylum", "Class", "Order",
                                            "Family", "Genus", "ESV_ID", "ESV_ID2")
otu <- otu_table(input_filt_abund_combined$data_loaded, taxa_are_rows = T)
tax <- tax_table(as.matrix(input_filt_abund_combined$taxonomy_loaded))
map <- sample_data(input_filt_abund_combined$map_loaded)
input.phy <- phyloseq(otu, tax, map)
se.mb2 <- spiec.easi(input.phy, method='mb', lambda.min.ratio=1e-2,
                     nlambda=20, pulsar.params=list(rep.num=50))
net <- adj2igraph(getRefit(se.mb2),  vertex.attr=list(name=taxa_names(input.phy)))

# Phyloseq plot
plot_network(net, 
             input.phy, 
             type='taxa', 
             color = "Phylum",
             shape = "Domain",
             point_size = 3,
             label = NULL)

plot_network(net, 
             input.phy, 
             type='taxa', 
             color = "Phylum",
             shape = "Domain",
             point_size = 3,
             label = NULL,
             layout.method = layout.circle)

pn <- plot_network(net, 
                   input.phy, 
                   type='taxa', 
                   color = "Phylum",
                   shape = "Domain",
                   point_size = 3,
                   label = NULL,
                   layout.method = layout.circle)
ggplotly(pn) # Doesn't show all edges! Glitch with ggplotly and networks?



#### _Igraph ####

# Stats
E(net) # 329
V(net) # 200
transitivity(net) # Average clustering coefficient. 0.06.
deg <- degree(net, mode="all")
mean(deg) # 3.29

# Add taxonomic information and color
# Add phylum, which is in the taxonomy_loaded table
# Check phylum numbers
table(input_filt_abund_combined$taxonomy_loaded$Phylum)

# Check match
sum(V(net)$name != rownames(input_filt_abund_combined$taxonomy_loaded)) # Check match

# Make factor
input_filt_abund_combined$taxonomy_loaded$Phylum <- as.factor(input_filt_abund_combined$taxonomy_loaded$Phylum)

# Confer to network
V(net)$phylum = input_filt_abund_combined$taxonomy_loaded$Phylum

# Check levels
levels(input_filt_abund_combined$taxonomy_loaded$Phylum) # There are 19 phyla

# Set n to number of levels
n <- length(levels(input_filt_abund_combined$taxonomy_loaded$Phylum))

# Save taxonomy and colors in tax
tax <- input_filt_abund_combined$taxonomy_loaded

# Get colors for n levels
#colrs <- hue_pal()(n) # ggplot hue palette
colrs <- colorRampPalette(brewer.pal(12, "Paired"))(n) # expanded Rcolorbrewer paired palette
tax$color <- tax$Phylum
levels(tax$Phylum)
tax$color <- recode_factor(tax$color,
                           "Acidobacteriota" = colrs[1],
                           "Actinobacteriota" = colrs[2],
                           "Amoebozoa" = colrs[3],
                           "Archaeplastida" = colrs[4],
                           "Armatimonadota" = colrs[5],
                           "Bacteroidota" = colrs[6],
                           "Bdellovibrionota" = colrs[7],
                           "Chloroflexi" = colrs[8],
                           "Cyanobacteria" = colrs[9],
                           "Deinococcota" = colrs[10],
                           "Excavata" = colrs[11],
                           "Incertae_Sedis" = colrs[12],
                           "Myxococcota" = colrs[13],
                           "NA" = colrs[14],
                           "Opisthokonta" = colrs[15],
                           "Planctomycetota" = colrs[16],
                           "Proteobacteria" = colrs[17],
                           "SAR" = colrs[18],
                           "Verrucomicrobiota" = colrs[19])
V(net)$color <- as.character(tax$color)

#### __Basic ####
# Plot. Layout in Circle (this is what Farrer et al. 2019 did)
par(mar = c(0,8,0,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color, 
     vertex.size = deg*2, 
     vertex.shape = "circle", 
     vertex.frame.color = "black",
     vertex.label = NA, 
     #edge.color = ifelse(cor.matrix$r > 0, "#619CFF","#F8766D"),
     edge.curved = 0.2,
     edge.width = 0.2,
     layout = layout_in_circle(net, order = order(V(net)$phylum)))
legend(x = -1.6, y = 0.5, levels(input_filt_abund_combined$taxonomy_loaded$Phylum), 
       pch = 21, col = "black", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1)

# Good. Basics are working.
# Now, refine taxonomy, shapes, color edges, other presentation things

# Check Modules
rng_adj <- igraph::get.adjacency(net, sparse = FALSE)
netcarto(rng_adj)
edge.betweenness.community(net)
fastgreedy.community(net)
walktrap.community(net)
spinglass.community(net)
leading.eigenvector.community(net)
label.propagation.community(net)
cluster_louvain(net)

# Color edge by weight
se.mb2$lambda
bm <- symBeta(getOptBeta(se.mb2), mode = "maxabs")
diag(bm) <- 0
weights <- Matrix::summary(t(bm))[,3]
E(net)$weight <- weights
E(net)[weight > 0]$color <- "red" # positive red
E(net)[weight < 0]$color <- "blue"  # negative blue
edge_attr(net)

# Vertex shape by role
adj <- as.matrix(as_adjacency_matrix(net))
role <- netcarto(adj)
role <- role[[1]]
v.names <- data.frame(name = V(net)$name)
v.names.ni <- v.names %>%
  filter(name %notin% role$name) %>%
  mutate(module = NA,
         connectivity = NA,
         participation = NA,
         role = NA)
role <- rbind(role, v.names.ni) %>%
  mutate(role = as.factor(role)) %>%
  droplevels() %>%
  mutate(shape = recode_factor(role,
                               "Peripheral Hub" = "square",
                               "Connector Hub" = "square",
                               "Kinless Hub" = "square",
                               "Connector" = "diamond",
                               "Kinless" = "triangle",
                               "Peripheral" = "circle",
                               "Ultra peripheral" = "circle")) %>%
  mutate(shape = as.character(shape)) %>%
  mutate(shape = replace_na(shape, "circle"))
vertex_attr(net)
V(net)$shape <- role$shape
hubcon <- role %>%
  filter(role == "Peripheral Hub" | role == "Connector Hub") %>%
  left_join(., input_filt_abund_combined$taxonomy_loaded, by = c("name" = "ESV_ID2"))

# Info for labels
length(V(net))
length(E(net))
round(mean(deg), 1)
round(transitivity(net), 3)

#### __+ -, Roles ####
# Plot Circle
par(mar = c(2,8,0,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color,
     vertex.size = deg*2,
     vertex.shape = V(net)$shape,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net, order = order(V(net)$phylum)))
legend(x = -1.7, y = 1.05, levels(input_filt_abund_combined$taxonomy_loaded$Phylum), 
       pch = 21, col = "black", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.adj = 0)
legend(x = -1.7, y = -0.3, c("Peripheral", "Connector", "Hub"), pch = c(21, 23, 22), 
       col = "black", pt.cex = 2, cex = 0.8, bty = "n", ncol = 1, 
       title = "Role", title.adj = 0)
legend(x = -1.7, y = -0.65, c("Positive", "Negative"), lwd = 2,
       col = c("red", "blue"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.adj = 0)
title(main = "Svalbard Patch Network (SPIEC-EASI)", adj = 0.5, line = -1.5, cex.main = 1.5)
mtext("Nodes = 200", side = 1, line = -2, cex = 0.8)
mtext("Edges = 329", side = 1, line = -1, cex = 0.8)
mtext("Mean Degrees = 3.3", side = 1, line = 0, cex = 0.8)
mtext("Clustering Coefficient = 0.062", side = 1, line = 1, cex = 0.8)

# Make Interactive? Can do ggplotly on the plot_network function (kind of) but not base graphic.

# Now just need to update taxonomy
# Need to label Bacteria vs. Eukaryota (B_, E_)
# Need to get Bacteria phylum and Eukaryota taxonomy4 ("Order")
View(input_filt_abund_combined$taxonomy_loaded)
input_filt_abund_combined$taxonomy_loaded <- input_filt_abund_combined$taxonomy_loaded %>%
  mutate(Phylum = as.character(Phylum)) %>%
  mutate(taxonomy8 = ifelse(Domain == "Bacteria",
                            Phylum,
                            Order)) %>%
  mutate(taxonomy = paste(Domain, taxonomy8, sep = "_")) %>%
  mutate(taxonomy = gsub("Bacteria", "B", taxonomy)) %>%
  mutate(taxonomy = gsub("Eukaryota", "E", taxonomy))

# Also make lichen functional group column. Will use that later.
View(input_filt_abund_combined$taxonomy_loaded)
input_filt_abund_combined$taxonomy_loaded <- input_filt_abund_combined$taxonomy_loaded %>%
  mutate(Group = ifelse(taxonomy == "B_Cyanobacteria",
                        "Cyanobacteria",
                        ifelse(taxonomy %in% c("E_Charophyta", "E_Chlorophyta"),
                                               "Algae",
                                               ifelse(taxonomy == "E_Fungi",
                                                      "Fungi",
                                                      NA))))

# Run as above but with taxonomy instead of Phylum
sum(V(net)$name != rownames(input_filt_abund_combined$taxonomy_loaded)) # Check match
input_filt_abund_combined$taxonomy_loaded$taxonomy <- as.factor(input_filt_abund_combined$taxonomy_loaded$taxonomy)
V(net)$taxonomy = input_filt_abund_combined$taxonomy_loaded$taxonomy
levels(input_filt_abund_combined$taxonomy_loaded$taxonomy) # There are 24 phyla
n <- length(levels(input_filt_abund_combined$taxonomy_loaded$taxonomy))
tax <- input_filt_abund_combined$taxonomy_loaded
colrs <- colorRampPalette(brewer.pal(12, "Paired"))(n) # expanded Rcolorbrewer paired palette
tax$color <- tax$taxonomy
levels(tax$taxonomy)
tax$color <- recode_factor(tax$color,
                           "B_Acidobacteriota" = colrs[1],
                           "B_Actinobacteriota" = colrs[2],
                           "B_Armatimonadota" = colrs[3],
                           "B_Bacteroidota" = colrs[4],
                           "B_Bdellovibrionota" = colrs[5],
                           "B_Chloroflexi" = colrs[6],
                           "B_Cyanobacteria" = colrs[7],
                           "B_Deinococcota" = colrs[8],
                           "B_Myxococcota" = colrs[9],
                           "B_Planctomycetota" = colrs[10],
                           "B_Proteobacteria" = colrs[11],
                           "B_Verrucomicrobiota" = colrs[12],
                           "E_Apusomonas" = colrs[13],
                           "E_Cercozoa" = colrs[14],
                           "E_Charophyta" = colrs[15],
                           "E_Chlorophyta" = colrs[16],
                           "E_Dictyamoeba" = colrs[17],
                           "E_Discicristata" = colrs[18],
                           "E_Fungi" = colrs[19],
                           "E_Labyrinthulomycetes" = colrs[20],
                           "E_MAST-12" = colrs[21],
                           "E_NA" = colrs[22],
                           "E_Ochrophyta" = colrs[23],
                           "E_Peronosporomycetes" = colrs[24])
V(net)$color <- as.character(tax$color)

# Also make a Group color variable
sum(V(net)$name != rownames(input_filt_abund_combined$taxonomy_loaded)) # Check match
input_filt_abund_combined$taxonomy_loaded$Group <- as.factor(input_filt_abund_combined$taxonomy_loaded$Group)
V(net)$Group = input_filt_abund_combined$taxonomy_loaded$Group
levels(input_filt_abund_combined$taxonomy_loaded$Group) # There are 3 levels
tax$groupcolor <- input_filt_abund_combined$taxonomy_loaded$Group
levels(tax$groupcolor)
tax$groupcolor <- recode_factor(tax$groupcolor,
                                "Algae" = "green",
                                "Cyanobacteria" = "green4",
                                "Fungi" = "brown")
V(net)$groupcolor <- as.character(tax$groupcolor)

# And add ESV ID as vertex attribute
V(net)$ESV_ID2 <- input_filt_abund_combined$taxonomy_loaded$ESV_ID2



#### __Good Taxonomy ####
# Replot, Circle, shape by Role
par(mar = c(3,8,1,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color,
     vertex.size = deg*2,
     vertex.shape = V(net)$shape,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net, order = order(V(net)$taxonomy)))
legend(x = -1.7, y = 1.2, levels(input_filt_abund_combined$taxonomy_loaded$taxonomy), 
       pch = 21, col = "black", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.adj = 0)
legend(x = -1.7, y = -0.5, c("Peripheral", "Connector", "Hub"), pch = c(21, 23, 22), 
       col = "black", pt.cex = 2, cex = 0.8, bty = "n", ncol = 1, 
       title = "Role", title.adj = 0)
legend(x = -1.7, y = -0.8, c("Positive", "Negative"), lwd = 2,
       col = c("red", "blue"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.adj = 0)
title(main = "Svalbard Patch Network (SPIEC-EASI)", adj = 0.5, line = -1, cex.main = 1.5)
mtext("Nodes = 200", side = 1, line = -1.5, cex = 0.8)
mtext("Edges = 329", side = 1, line = -0.5, cex = 0.8)
mtext("Mean Degrees = 3.3", side = 1, line = 0.5, cex = 0.8)
mtext("Clustering Coefficient = 0.062", side = 1, line = 1.5, cex = 0.8)

# Replot, Circle, shape by Domain
sum(V(net)$name != rownames(input_filt_abund_combined$taxonomy_loaded)) # Check match
V(net)$Domain = input_filt_abund_combined$taxonomy_loaded$Domain
V(net)$Domain <- gsub("Bacteria", "circle", V(net)$Domain)
V(net)$Domain <- gsub("Eukaryota", "triangle", V(net)$Domain)

png("Network.png", width = 7, height = 6, units = "in", res = 300)
par(mar = c(3,8,1,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color,
     vertex.size = deg*2,
     vertex.shape = V(net)$Domain,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net, order = order(V(net)$taxonomy)))
legend(x = -1.7, y = 1.25, levels(input_filt_abund_combined$taxonomy_loaded$taxonomy), 
       pch = c(rep(21, 12), rep(24, 12)), col = "black", pt.bg = colrs, pt.cex = 1.5, 
       cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.cex = 1.1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.7, y = -0.8, c("Bacteria", "Eukaryota"), pch = c(21, 24), 
       col = "black", pt.cex = 1.5, cex = 0.8, bty = "n", ncol = 1, 
       title = "Domain", title.cex = 1.1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.7, y = -1.1, c("Positive", "Negative"), lwd = 2,
       col = c("red", "blue"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.cex = 1.1, title.adj = 0)
title(main = "Svalbard Patch Network (SPIEC-EASI)", adj = 0.5, line = -1, cex.main = 1.5)
mtext("Nodes = 200", side = 1, line = -1.5, cex = 0.8)
mtext("Edges = 329", side = 1, line = -0.5, cex = 0.8)
mtext("Mean Degrees = 3.3", side = 1, line = 0.5, cex = 0.8)
mtext("Clustering Coefficient = 0.062", side = 1, line = 1.5, cex = 0.8)
dev.off()

# Can we highlight the potential lichen edges by making them super thick?
# Need to set the edge.width up to 3 on those instead of 0.5
edge_attr(net) # Have weight and color attribute, need to make size attribute
# Need to isolate the lichen edges
E(net_lich)
# 2 lichen edges. 
# ASV_Prok_64 --ASV_Euk_21 
# ASV_Euk_32 --ASV_Euk_27
E(net)
print(E(net))
elist <- get.edgelist(net) %>%
  as.data.frame() %>%
  set_names(c("V1", "V2")) %>%
  mutate(Lichen = ifelse(V1 == "ASV_Prok_64" & V2 == "ASV_Euk_21",
                         "Yes",
                         ifelse(V1 == "ASV_Euk_27" & V2 == "ASV_Euk_32",
                         "Yes",
                         "No")))
E(net)$Lichen <- elist$Lichen
E(net)$Size <- elist$Lichen
E(net)$Size <- gsub("Yes", 5, E(net)$Size)
E(net)$Size <- gsub("No", 0.25, E(net)$Size)
E(net)$Size <- as.numeric(E(net)$Size)

tax_lab <- input_filt_abund_combined$taxonomy_loaded %>%
  group_by(taxonomy) %>%
  slice_head(n = 1)

png("Network_LichenHighlight.png", width = 7, height = 6, units = "in", res = 300)
par(mar = c(2,8,0,0), xpd = TRUE) # bottom, left, top, right
plot(net,
     vertex.color = V(net)$color,
     vertex.size = deg*2,
     vertex.shape = V(net)$Domain,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0,
     edge.width = E(net)$Size,
     layout = layout_in_circle(net, order = order(V(net)$taxonomy)))
legend(x = -1.7, y = 1.2, tax_lab$taxonomy8, 
       pch = c(rep(21, 12), rep(24, 12)), col = "black", pt.bg = colrs, pt.cex = 1.5, 
       cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.cex = 1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.7, y = -0.76, c("Bacteria", "Eukaryota"), pch = c(21, 24), 
       col = "black", pt.cex = 1.5, cex = 0.8, bty = "n", ncol = 1, 
       title = "Domain", title.cex = 1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.7, y = -1.05, c("Positive", "Negative"), lwd = 2,
       col = c("red", "blue"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.cex = 1, title.adj = 0)
mtext("Nodes = 200, Edges = 329", side = 1, line = -1.5, cex = 0.8)
mtext("Mean Degrees = 3.3", side = 1, line = -0.75, cex = 0.8)
mtext("Clustering Coefficient = 0.062", side = 1, line = 0, cex = 0.8)
dev.off()



#### __Positive ####
net_pos <- delete.edges(net, E(net) [ weight < 0])

# Info for labels
length(V(net_pos))
length(E(net_pos))
deg_pos <- degree(net_pos, mode="all")
round(mean(deg_pos), 1)
round(transitivity(net_pos), 3)

# Plot Circle
png("Network_Pos.png", width = 7, height = 6, units = "in", res = 300)
par(mar = c(3,8,1,0), xpd = TRUE) # bottom, left, top, right
plot(net_pos,
     vertex.color = V(net_pos)$color,
     vertex.size = deg*2,
     vertex.shape = V(net_pos)$Domain,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net_pos, order = order(V(net_pos)$taxonomy)))
legend(x = -1.7, y = 1.25, levels(input_filt_abund_combined$taxonomy_loaded$taxonomy), 
       pch = c(rep(21, 12), rep(24, 12)), col = "black", pt.bg = colrs, pt.cex = 1.5, 
       cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.cex = 1.1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.7, y = -0.8, c("Bacteria", "Eukaryota"), pch = c(21, 24), 
       col = "black", pt.cex = 1.5, cex = 0.8, bty = "n", ncol = 1, 
       title = "Domain", title.cex = 1.1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.7, y = -1.1, c("Positive", "Negative"), lwd = 2,
       col = c("red", "blue"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.cex = 1.1, title.adj = 0)
title(main = "Svalbard Positive Connections (SPIEC-EASI)", adj = 0.5, line = -1, cex.main = 1.5)
mtext("Nodes = 200", side = 1, line = -1.5, cex = 0.8)
mtext("Edges = 226", side = 1, line = -0.5, cex = 0.8)
mtext("Mean Degrees = 2.3", side = 1, line = 0.5, cex = 0.8)
mtext("Clustering Coefficient = 0.079", side = 1, line = 1.5, cex = 0.8)
dev.off()



#### __Lichen ####
# Take the positive network and subset to just Fungi, Cyanobacteria, Algae
net_lich <- delete.vertices(net_pos, V(net_pos)[V(net_pos)$taxonomy %notin% c("B_Cyanobacteria",
                                                                              "E_Charophyta",
                                                                              "E_Chlorophyta",
                                                                              "E_Fungi")])

# Info for labels
length(V(net_lich)) # 44
length(E(net_lich)) # 7
deg_lich <- degree(net_lich, mode="all")
round(mean(deg_lich), 1) # 0.3
round(transitivity(net_lich), 3) # NA

# Plot Circle
png("Network_Lichen.png", width = 7, height = 6, units = "in", res = 300)
par(mar = c(4,8,2,0), xpd = TRUE) # bottom, left, top, right
plot(net_lich,
     vertex.color = V(net_lich)$groupcolor,
     vertex.size = (deg_lich+1)*5,
     vertex.shape = V(net_lich)$Domain,
     vertex.frame.color = "black",
     vertex.label = NA,
     #vertex.label = V(net_lich)$ESV_ID2, # Check IDs. Make note of lichens
     edge.curved = 0.5,
     edge.width = 2,
     layout = layout_in_circle(net_lich, order = order(V(net_lich)$Group)))
legend(x = -1.7, y = 0.8, levels(input_filt_abund_combined$taxonomy_loaded$Group), 
       pch = c(rep(21, 12), rep(24, 12)), col = "black", pt.bg = c("green", "green4", "brown"), 
       pt.cex = 2, cex = 0.8, bty = "n", ncol = 1,
       title = "Group", title.cex = 1.1, title.adj = 0, y.intersp = 1.2)
legend(x = -1.7, y = 0.3, c("Bacteria", "Eukaryota"), pch = c(21, 24), 
       col = "black", pt.cex = 1.5, cex = 0.8, bty = "n", ncol = 1, 
       title = "Domain", title.cex = 1.1, title.adj = 0, y.intersp = 1.2)
legend(x = -1.7, y = -0.2, c("Positive", "Negative"), lwd = 2,
       col = c("red", "blue"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.cex = 1.1, title.adj = 0)
title(main = "Cyano-Algae-Fungi Positive Connections", adj = 0.5, line = -0.5, cex.main = 1.5)
mtext("Nodes = 44", side = 1, line = -0.75, cex = 0.8)
mtext("Edges = 7; Potential lichens = 2", side = 1, line = 0.25, cex = 0.8)
mtext("Mean Degrees = 0.3", side = 1, line = 1.25, cex = 0.8)
mtext("Clustering Coefficient = 0", side = 1, line = 2.25, cex = 0.8)
dev.off()

# Lichens are:
View(input_filt_abund_combined$taxonomy_loaded)
# ASV_Prok_64 - ASV_Euk_21
# Cyanobacteria Cyanobacteriia Gloeobacterales Gloeobacteraceae Gloeobacter PCC-7421 ESV_45
# Opisthokonta Nucletmycea Fungi Herpotrichiellaceae NA ESV_7

# ASV_Euk_32 - ASV_Euk_27
# Eukaryota Archaeplastida Chloroplastida Chlorophyta Chlamydomonadales Chlamydomonas ESV_924
# Eukaryota Opisthokonta Nucletmycea Fungi NA NA ESV_497



#### _Betweenness ####
# Plot degree versus betweenness for each network (each state)
se <- se.mb2
net <- adj2igraph(getRefit(se),  vertex.attr=list(name=taxa_names(input.phy)))
bw <- data.frame("Degree" = degree(net, mode="all"),
                     "Betweenness" = betweenness(net)) %>%
  mutate("ASV" = rownames(.)) %>%
  left_join(., input_filt_abund_combined$taxonomy_loaded, by = c("ASV" = "ESV_ID2"))

nb.cols <- length(levels(bw$Phylum))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bw, aes(Degree, Betweenness, colour = Phylum)) +
  geom_jitter(shape = 16, size = 2, alpha = 1, width = 0.15) +
  labs(x = "Degree",
       y = "Betweenness",
       colour = "Phylum") +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(limits = c(0, 10),
                     breaks = c(0,1,2,3,4,5,6,7,8, 9, 10)) +
  theme_bw() +
  theme(legend.key.size = unit(0.1, "cm"))

# Check ASVs, add labels
ggplot(bw, aes(Degree, Betweenness, colour = Phylum)) +
  geom_jitter(shape = 16, size = 2, alpha = 1, width = 0.15) +
  geom_text(data = bw,
            aes(x = Degree, y = Betweenness, label = ASV), 
            size = 3, inherit.aes = F) +
  labs(x = "Degree",
       y = "Betweenness",
       colour = "Phylum") +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(limits = c(0, 10),
                     breaks = c(0,1,2,3,4,5,6,7,8, 9, 10)) +
  theme_bw() +
  theme(legend.key.size = unit(0.1, "cm"))

asvs_to_label <- c("ASV_Prok_37", "ASV_Euk_67", "ASV_Euk_91", "ASV_Prok_62")
bw_asvs <- bw %>%
  filter(ASV %in% asvs_to_label) %>%
  filter(Betweenness > 1500) %>%
  mutate("HighTax" = ifelse(Genus != "NA",
                            Genus,
                            ifelse(Family != "NA",
                                   Family,
                                   ifelse(Order != "NA",
                                          Order,
                                          ifelse(Class != "NA",
                                                 Class,
                                                 ifelse(Phylum != "NA",
                                                        Phylum,
                                                        Domain)))))) %>%
  mutate("HightTax_ESV" = paste(HighTax, ESV_ID, sep = "_"))

png("Degree_Betweenness.png", width = 7, height = 3, units = "in", res = 300)
ggplot(bw, aes(Degree, Betweenness, colour = Phylum)) +
  geom_jitter(shape = 16, size = 2, alpha = 1, 
              position = position_jitter(seed = 1, width = 0.15)) +
  geom_text_repel(data = bw_asvs,
                  #min.segment.length = 0,
                  aes(x = Degree, y = Betweenness, label = HightTax_ESV), 
                  size = 2, inherit.aes = F,
                  position = position_jitter(seed = 1, width = 0.15)) +
  labs(x = "Degree",
       y = "Betweenness",
       colour = "Phylum") +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(limits = c(0, 10),
                     breaks = c(0,1,2,3,4,5,6,7,8, 9, 10)) +
  theme_bw() +
  theme(legend.key.size = unit(0.1, "cm"))
dev.off()



#### _Participation ####
# Plot participation coefficient and within-module degree (Barnes et al.)
# Dashed lines at 0.61 and 2.2
# z-score (within module degree) is "connectivity"
# participation coefficient P is "participation
roles <- role %>%
  left_join(., input_filt_abund_combined$taxonomy_loaded, by = c("name" = "ESV_ID2")) %>%
  mutate(Phylum = as.factor(Phylum),
         role = as.factor(role)) %>%
  filter(is.na(module) == F) %>%
  droplevels()

nb.cols <- length(levels(roles$Phylum))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(roles, aes(participation, connectivity, colour = Phylum, shape = role)) +
  geom_vline(xintercept = 0.62, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_point(size = 2, alpha = 0.75) +
  labs(x = "Participation coefficient (P)",
       y = "Within module degree (z)",
       colour = "Phylum",
       shape = "Network Role") +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.key.size = unit(0.1, "cm"),
        panel.grid = element_blank())

# Check ASVs, add labels
ggplot(roles, aes(participation, connectivity, colour = Phylum, shape = role)) +
  geom_vline(xintercept = 0.62, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text(data = roles,
            aes(x = participation, y = connectivity, label = name), 
            size = 3, inherit.aes = F) +
  labs(x = "Participation coefficient (P)",
       y = "Within module degree (z)",
       colour = "Phylum",
       shape = "Network Role") +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.key.size = unit(0.1, "cm"),
        panel.grid = element_blank())

pz_asvs <- roles %>%
  filter(connectivity > 2.5 | participation > 0.62) %>%
  mutate("HighTax" = ifelse(Genus != "NA",
                            Genus,
                            ifelse(Family != "NA",
                                   Family,
                                   ifelse(Order != "NA",
                                          Order,
                                          Class))))

png("Participation_Connectivity.png", width = 7, height = 5, units = "in", res = 300)
ggplot(roles, aes(participation, connectivity, colour = Phylum, shape = role)) +
  geom_vline(xintercept = 0.62, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text_repel(data = pz_asvs,
                  min.segment.length = 0,
                  aes(x = participation, y = connectivity, label = HighTax), 
                  size = 2, inherit.aes = F) +
  labs(x = "Participation coefficient (P)",
       y = "Within module degree (z)",
       colour = "Phylum",
       shape = "Network Role") +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.key.size = unit(0.1, "cm"),
        legend.margin = margin(0,0,0,0),
        panel.grid = element_blank())
dev.off()

keystone <- roles %>%
  filter(participation > 0.62 | connectivity > 2.5)

#write_xlsx(keystone, format_headers = F, "keystone_taxa.xlsx")



#### 4. Spot Only ####
# Now repeat the above analysis but just for spots. Remove non-spot (control) samples
# n will now be 12 not 24.
spot_16S <- filter_data(input_16S_filt,
                        "sampleID",
                        keep_vals = c("P1.1A", "P1.2A", "P1.3A", "P1.4A", "P1.5A", "P1.6A",
                                      "P2.1A", "P2.2A", "P2.3A", "P2.4A", "P2.5A", "P2.6A"))
spot_18S <- filter_data(input_18S_filt,
                        "sampleID",
                        keep_vals = c("P1.1A", "P1.2A", "P1.3A", "P1.4A", "P1.5A", "P1.6A",
                                      "P2.1A", "P2.2A", "P2.3A", "P2.4A", "P2.5A", "P2.6A"))

#### _Abundant ####
nrow(spot_16S$data_loaded) # 1777
nrow(spot_18S$data_loaded) # 358

# Let's try getting the top 100 of each, then merging.
top_16S <- as.data.frame(rowMeans(spot_16S$data_loaded)) %>%
  set_names("MeanCount") %>%
  arrange(desc(MeanCount)) %>%
  slice_head(n = 100)
spot_16S_abund <- filter_taxa_from_input(spot_16S,
                                               taxa_IDs_to_keep = rownames(top_16S),
                                               at_spec_level = 7)
nrow(spot_16S_abund$data_loaded) # 100

top_18S <- as.data.frame(rowMeans(spot_18S$data_loaded)) %>%
  set_names("MeanCount") %>%
  arrange(desc(MeanCount)) %>%
  slice_head(n = 100)
spot_18S_abund <- filter_taxa_from_input(spot_18S,
                                               taxa_IDs_to_keep = rownames(top_18S),
                                               at_spec_level = 7) # 46 removed
nrow(spot_18S_abund$data_loaded) # 100

# Need to create unique prok and euk ESV (ASV) IDs
# Unique ASV IDs are needed for rownames of $data_loaded and $taxonomy_loaded
spot_16S_abund_toMerge <- spot_16S_abund
num_16S <- seq(1:nrow(spot_16S_abund_toMerge$data_loaded))
ASVlabel_16S <- rep("ASV_Prok_", nrow(spot_16S_abund_toMerge$data_loaded))
IDs_16S <- paste(ASVlabel_16S, num_16S, sep = "")
rownames(spot_16S_abund_toMerge$data_loaded) <- IDs_16S
rownames(spot_16S_abund_toMerge$taxonomy_loaded) <- IDs_16S
spot_16S_abund_toMerge$taxonomy_loaded$taxonomy8 <- IDs_16S

spot_18S_abund_toMerge <- spot_18S_abund
num_18S <- seq(1:nrow(spot_18S_abund_toMerge$data_loaded))
ASVlabel_18S <- rep("ASV_Euk_", nrow(spot_18S_abund_toMerge$data_loaded))
IDs_18S <- paste(ASVlabel_18S, num_18S, sep = "")
rownames(spot_18S_abund_toMerge$data_loaded) <- IDs_18S
rownames(spot_18S_abund_toMerge$taxonomy_loaded) <- IDs_18S
spot_18S_abund_toMerge$taxonomy_loaded$taxonomy8 <- IDs_18S

spot_abund_combined <- spot_16S_abund_toMerge
spot_abund_combined$data_loaded <- rbind(spot_abund_combined$data_loaded,
                                               spot_18S_abund_toMerge$data_loaded)
spot_abund_combined$taxonomy_loaded <- rbind(spot_abund_combined$taxonomy_loaded,
                                                   spot_18S_abund_toMerge$taxonomy_loaded)

nrow(spot_abund_combined$data_loaded) # 200 total taxa



#### _phyloseq ####
# Convert mctoolsr to phyloseq
names(spot_abund_combined$taxonomy_loaded) <- c("Domain", "Phylum", "Class", "Order",
                                                "Family", "Genus", "ESV_ID", "ESV_ID2")
otu <- otu_table(spot_abund_combined$data_loaded, taxa_are_rows = T)
tax <- tax_table(as.matrix(spot_abund_combined$taxonomy_loaded))
map <- sample_data(spot_abund_combined$map_loaded)
spot.phy <- phyloseq(otu, tax, map)
se.mb2_spot <- spiec.easi(spot.phy, method='mb', lambda.min.ratio=1e-2,
                     nlambda=20, pulsar.params=list(rep.num=50))
net_spot <- adj2igraph(getRefit(se.mb2_spot),  vertex.attr=list(name=taxa_names(spot.phy)))

# Phyloseq plot
plot_network(net_spot, 
             spot.phy, 
             type='taxa', 
             color = "Phylum",
             shape = "Domain",
             point_size = 3,
             label = NULL)

plot_network(net_spot, 
             spot.phy, 
             type='taxa', 
             color = "Phylum",
             shape = "Domain",
             point_size = 3,
             label = NULL,
             layout.method = layout.circle)

pn <- plot_network(net_spot, 
                   spot.phy, 
                   type='taxa', 
                   color = "Phylum",
                   shape = "Domain",
                   point_size = 3,
                   label = NULL,
                   layout.method = layout.circle)
ggplotly(pn) # Doesn't show all edges! Don't know why. Code is exactly the same as above. Too bad.



# Stats
E(net_spot) # 283
V(net_spot) # 200
transitivity(net_spot) # Average clustering coefficient. 0.093
deg_spot <- degree(net_spot, mode="all")
mean(deg_spot) # 2.83

# Add taxonomic information and color
# Add phylum, which is in the taxonomy_loaded table
# Check phylum numbers
table(spot_abund_combined$taxonomy_loaded$Phylum)

# Check match
sum(V(net_spot)$name != rownames(spot_abund_combined$taxonomy_loaded)) # Check match

# Make factor
spot_abund_combined$taxonomy_loaded$Phylum <- as.factor(spot_abund_combined$taxonomy_loaded$Phylum)

# Confer to network
V(net_spot)$phylum = spot_abund_combined$taxonomy_loaded$Phylum

# Check levels
levels(spot_abund_combined$taxonomy_loaded$Phylum) # There are 18 phyla

# Set n to number of levels
n <- length(levels(spot_abund_combined$taxonomy_loaded$Phylum))

# Save taxonomy and colors in tax
tax <- spot_abund_combined$taxonomy_loaded

# Get colors for n levels
#colrs <- hue_pal()(n) # ggplot hue palette
colrs <- colorRampPalette(brewer.pal(12, "Paired"))(n) # expanded Rcolorbrewer paired palette
tax$color <- tax$Phylum
levels(tax$Phylum)
tax$color <- recode_factor(tax$color,
                           "Acidobacteriota" = colrs[1],
                           "Actinobacteriota" = colrs[2],
                           "Amoebozoa" = colrs[3],
                           "Archaeplastida" = colrs[4],
                           "Armatimonadota" = colrs[5],
                           "Bacteroidota" = colrs[6],
                           #"Bdellovibrionota" = colrs[7],
                           "Chloroflexi" = colrs[7],
                           "Cyanobacteria" = colrs[8],
                           "Deinococcota" = colrs[9],
                           "Excavata" = colrs[10],
                           "Incertae_Sedis" = colrs[11],
                           "Myxococcota" = colrs[12],
                           "NA" = colrs[13],
                           "Opisthokonta" = colrs[14],
                           "Planctomycetota" = colrs[15],
                           "Proteobacteria" = colrs[16],
                           "SAR" = colrs[17],
                           "Verrucomicrobiota" = colrs[18])
V(net_spot)$color <- as.character(tax$color)



#### __Basic ####
# Plot. Layout in Circle (this is what Farrer et al. 2019 did)
par(mar = c(0,8,0,0), xpd = TRUE) # bottom, left, top, right
plot(net_spot,
     vertex.color = V(net_spot)$color, 
     vertex.size = deg_spot*2, 
     vertex.shape = "circle", 
     vertex.frame.color = "black",
     vertex.label = NA, 
     #edge.color = ifelse(cor.matrix$r > 0, "#619CFF","#F8766D"),
     edge.curved = 0.2,
     edge.width = 0.2,
     layout = layout_in_circle(net_spot, order = order(V(net_spot)$phylum)))
legend(x = -1.6, y = 0.5, levels(spot_abund_combined$taxonomy_loaded$Phylum), 
       pch = 21, col = "black", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1)

# Good. Basics are working.
# Now, refine taxonomy, shapes, color edges, other presentation things

# Check Modules
rng_adj_spot <- igraph::get.adjacency(net_spot, sparse = FALSE)
netcarto(rng_adj_spot)
edge.betweenness.community(net_spot)
fastgreedy.community(net_spot)
walktrap.community(net_spot)
spinglass.community(net_spot)
leading.eigenvector.community(net_spot)
label.propagation.community(net_spot)
cluster_louvain(net_spot)

# Color edge by weight
se.mb2_spot$lambda
bm_spot <- symBeta(getOptBeta(se.mb2_spot), mode = "maxabs")
diag(bm_spot) <- 0
weights_spot <- Matrix::summary(t(bm))[,3]
E(net_spot)$weight <- weights_spot
E(net_spot)[weight > 0]$color <-"red" # positive red
E(net_spot)[weight < 0]$color <-"blue"  # negative blue
edge_attr(net_spot)

# Vertex shape by role_spot
adj_spot <- as.matrix(as_adjacency_matrix(net_spot))
role_spot <- netcarto(adj_spot)
role_spot <- role_spot[[1]]
v.names_spot <- data.frame(name = V(net_spot)$name)
v.names.ni_spot <- v.names_spot %>%
  filter(name %notin% role_spot$name) %>%
  mutate(module = NA,
         connectivity = NA,
         participation = NA,
         role_spot = NA)
role_spot <- rbind(role_spot, v.names.ni_spot) %>%
  mutate(role = as.factor(role)) %>%
  droplevels() %>%
  mutate(shape = recode_factor(role,
                               "Peripheral Hub" = "square",
                               "Connector Hub" = "square",
                               "Kinless Hub" = "square",
                               "Connector" = "diamond",
                               "Kinless" = "triangle",
                               "Peripheral" = "circle",
                               "Ultra peripheral" = "circle")) %>%
  mutate(shape = as.character(shape)) %>%
  mutate(shape = replace_na(shape, "circle"))
vertex_attr(net_spot)
V(net_spot)$shape <- role_spot$shape
hubcon_spot <- role_spot %>%
  filter(role == "Peripheral Hub" | role == "Connector Hub") %>%
  left_join(., spot_abund_combined$taxonomy_loaded, by = c("name" = "ESV_ID2"))

# Info for labels
length(V(net_spot))
length(E(net_spot))
round(mean(deg_spot), 1)
round(transitivity(net_spot), 3)

#### __+ -, Roles ####
# Plot Circle
par(mar = c(2,8,0,0), xpd = TRUE) # bottom, left, top, right
plot(net_spot,
     vertex.color = V(net_spot)$color,
     vertex.size = deg_spot*2,
     vertex.shape = V(net_spot)$shape,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net_spot, order = order(V(net_spot)$phylum)))
legend(x = -1.7, y = 1.05, levels(spot_abund_combined$taxonomy_loaded$Phylum), 
       pch = 21, col = "black", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.adj = 0)
legend(x = -1.7, y = -0.3, c("Peripheral", "Connector", "Hub"), pch = c(21, 23, 22), 
       col = "black", pt.cex = 2, cex = 0.8, bty = "n", ncol = 1, 
       title = "Role", title.adj = 0)
legend(x = -1.7, y = -0.65, c("Positive", "Negative"), lwd = 2,
       col = c("red", "blue"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.adj = 0)
title(main = "Svalbard Spot Network (SPIEC-EASI)", adj = 0.5, line = -1.5, cex.main = 1.5)
mtext("Nodes = 200", side = 1, line = -2, cex = 0.8)
mtext("Edges = 283", side = 1, line = -1, cex = 0.8)
mtext("Mean Degrees = 2.8", side = 1, line = 0, cex = 0.8)
mtext("Clustering Coefficient = 0.093", side = 1, line = 1, cex = 0.8)

# Make Interactive? Can do ggplotly on the plot_network function but not base graphic.
# Can use networkD3 to plot an edgelist but can't format very well
net_spot_EL <- as.data.frame(as_edgelist(net_spot))
simpleNetwork(net_spot_EL,
              nodeColour = V(net_spot)$color)

# Now just need to update taxonomy
# Need to label Bacteria vs. Eukaryota (B_, E_)
# Need to get Bacteria phylum and Eukaryota taxonomy4 ("Order")
View(spot_abund_combined$taxonomy_loaded)
spot_abund_combined$taxonomy_loaded <- spot_abund_combined$taxonomy_loaded %>%
  mutate(Phylum = as.character(Phylum)) %>%
  mutate(taxonomy8 = ifelse(Domain == "Bacteria",
                            Phylum,
                            Order)) %>%
  mutate(taxonomy = paste(Domain, taxonomy8, sep = "_")) %>%
  mutate(taxonomy = gsub("Bacteria", "B", taxonomy)) %>%
  mutate(taxonomy = gsub("Eukaryota", "E", taxonomy))

# Also make lichen functional group column. Will use that later.
View(spot_abund_combined$taxonomy_loaded)
spot_abund_combined$taxonomy_loaded <- spot_abund_combined$taxonomy_loaded %>%
  mutate(Group = ifelse(taxonomy == "B_Cyanobacteria",
                        "Cyanobacteria",
                        ifelse(taxonomy %in% c("E_Charophyta", "E_Chlorophyta"),
                               "Algae",
                               ifelse(taxonomy == "E_Fungi",
                                      "Fungi",
                                      NA))))

# Run as above but with taxonomy instead of Phylum
sum(V(net_spot)$name != rownames(spot_abund_combined$taxonomy_loaded)) # Check match
spot_abund_combined$taxonomy_loaded$taxonomy <- as.factor(spot_abund_combined$taxonomy_loaded$taxonomy)
V(net_spot)$taxonomy = spot_abund_combined$taxonomy_loaded$taxonomy
levels(spot_abund_combined$taxonomy_loaded$taxonomy) # There are 24 phyla
n <- length(levels(spot_abund_combined$taxonomy_loaded$taxonomy))
tax <- spot_abund_combined$taxonomy_loaded
colrs <- colorRampPalette(brewer.pal(12, "Paired"))(n) # expanded Rcolorbrewer paired palette
tax$color <- tax$taxonomy
levels(tax$taxonomy)
tax$color <- recode_factor(tax$color,
                           "B_Acidobacteriota" = colrs[1],
                           "B_Actinobacteriota" = colrs[2],
                           "B_Armatimonadota" = colrs[3],
                           "B_Bacteroidota" = colrs[4],
                           #"B_Bdellovibrionota" = colrs[5],
                           "B_Chloroflexi" = colrs[5],
                           "B_Cyanobacteria" = colrs[6],
                           "B_Deinococcota" = colrs[7],
                           "B_Myxococcota" = colrs[8],
                           "B_Planctomycetota" = colrs[9],
                           "B_Proteobacteria" = colrs[10],
                           "B_Verrucomicrobiota" = colrs[11],
                           "E_Apusomonas" = colrs[12],
                           "E_Cercozoa" = colrs[13],
                           "E_Charophyta" = colrs[14],
                           "E_Chlorophyta" = colrs[15],
                           "E_Dictyamoeba" = colrs[16],
                           "E_Discicristata" = colrs[17],
                           "E_Fungi" = colrs[18],
                           "E_Labyrinthulomycetes" = colrs[19],
                           "E_MAST-12" = colrs[20],
                           "E_NA" = colrs[21],
                           "E_Ochrophyta" = colrs[22],
                           "E_Peronosporomycetes" = colrs[23])
V(net_spot)$color <- as.character(tax$color)

# Also make a Group color variable
sum(V(net_spot)$name != rownames(spot_abund_combined$taxonomy_loaded)) # Check match
spot_abund_combined$taxonomy_loaded$Group <- as.factor(spot_abund_combined$taxonomy_loaded$Group)
V(net_spot)$Group = spot_abund_combined$taxonomy_loaded$Group
levels(spot_abund_combined$taxonomy_loaded$Group) # There are 3 levels
tax$groupcolor <- spot_abund_combined$taxonomy_loaded$Group
levels(tax$groupcolor)
tax$groupcolor <- recode_factor(tax$groupcolor,
                                "Algae" = "green",
                                "Cyanobacteria" = "green4",
                                "Fungi" = "brown")
V(net_spot)$groupcolor <- as.character(tax$groupcolor)

# And add ESV ID as vertex attribute
V(net_spot)$ESV_ID2 <- spot_abund_combined$taxonomy_loaded$ESV_ID2



#### __Good Taxonomy ####
# Replot, Circle, shape by role_spot
par(mar = c(3,8,1,0), xpd = TRUE) # bottom, left, top, right
plot(net_spot,
     vertex.color = V(net_spot)$color,
     vertex.size = deg_spot*2,
     vertex.shape = V(net_spot)$shape,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net_spot, order = order(V(net_spot)$taxonomy)))
legend(x = -1.7, y = 1.2, levels(spot_abund_combined$taxonomy_loaded$taxonomy), 
       pch = 21, col = "black", pt.bg = colrs, pt.cex = 2, cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.adj = 0)
legend(x = -1.7, y = -0.5, c("Peripheral", "Connector", "Hub"), pch = c(21, 23, 22), 
       col = "black", pt.cex = 2, cex = 0.8, bty = "n", ncol = 1, 
       title = "Role", title.adj = 0)
legend(x = -1.7, y = -0.9, c("Positive", "Negative"), lwd = 2,
       col = c("red", "blue"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.adj = 0)
title(main = "Svalbard Spot Network (SPIEC-EASI)", adj = 0.5, line = -1, cex.main = 1.5)
mtext("Nodes = 200", side = 1, line = -1.5, cex = 0.8)
mtext("Edges = 283", side = 1, line = -0.5, cex = 0.8)
mtext("Mean Degrees = 2.8", side = 1, line = 0.5, cex = 0.8)
mtext("Clustering Coefficient = 0.093", side = 1, line = 1.5, cex = 0.8)

# Replot, Circle, shape by Domain
sum(V(net_spot)$name != rownames(spot_abund_combined$taxonomy_loaded)) # Check match
V(net_spot)$Domain = spot_abund_combined$taxonomy_loaded$Domain
V(net_spot)$Domain <- gsub("Bacteria", "circle", V(net_spot)$Domain)
V(net_spot)$Domain <- gsub("Eukaryota", "triangle", V(net_spot)$Domain)

png("Network_Spot.png", width = 7, height = 6, units = "in", res = 300)
par(mar = c(3,8,1,0), xpd = TRUE) # bottom, left, top, right
plot(net_spot,
     vertex.color = V(net_spot)$color,
     vertex.size = deg_spot*2,
     vertex.shape = V(net_spot)$Domain,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net_spot, order = order(V(net_spot)$taxonomy)))
legend(x = -1.7, y = 1.25, levels(spot_abund_combined$taxonomy_loaded$taxonomy), 
       pch = c(rep(21, 12), rep(24, 12)), col = "black", pt.bg = colrs, pt.cex = 1.5, 
       cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.cex = 1.1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.7, y = -0.8, c("Bacteria", "Eukaryota"), pch = c(21, 24), 
       col = "black", pt.cex = 1.5, cex = 0.8, bty = "n", ncol = 1, 
       title = "Domain", title.cex = 1.1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.7, y = -1.1, c("Positive", "Negative"), lwd = 2,
       col = c("red", "blue"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.cex = 1.1, title.adj = 0)
title(main = "Svalbard Spot Network (SPIEC-EASI)", adj = 0.5, line = -1, cex.main = 1.5)
mtext("Nodes = 200", side = 1, line = -1.5, cex = 0.8)
mtext("Edges = 283", side = 1, line = -0.5, cex = 0.8)
mtext("Mean Degrees = 2.8", side = 1, line = 0.5, cex = 0.8)
mtext("Clustering Coefficient = 0.093", side = 1, line = 1.5, cex = 0.8)
dev.off()



#### __Positive ####
net_spot_pos <- delete.edges(net_spot, E(net_spot) [ weight < 0])

# Info for labels
length(V(net_spot_pos))
length(E(net_spot_pos))
deg_pos_spot <- degree(net_spot_pos, mode="all")
round(mean(deg_pos_spot), 1)
round(transitivity(net_spot_pos), 3)

# Plot Circle
png("Network_Pos_Spot.png", width = 7, height = 6, units = "in", res = 300)
par(mar = c(3,8,1,0), xpd = TRUE) # bottom, left, top, right
plot(net_spot_pos,
     vertex.color = V(net_spot_pos)$color,
     vertex.size = deg_spot*2,
     vertex.shape = V(net_spot_pos)$Domain,
     vertex.frame.color = "black",
     vertex.label = NA, 
     edge.curved = 0.2,
     edge.width = 0.5,
     layout = layout_in_circle(net_spot_pos, order = order(V(net_spot_pos)$taxonomy)))
legend(x = -1.7, y = 1.25, levels(spot_abund_combined$taxonomy_loaded$taxonomy), 
       pch = c(rep(21, 12), rep(24, 12)), col = "black", pt.bg = colrs, pt.cex = 1.5, 
       cex = 0.8, bty = "n", ncol = 1,
       title = "Taxonomy", title.cex = 1.1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.7, y = -0.8, c("Bacteria", "Eukaryota"), pch = c(21, 24), 
       col = "black", pt.cex = 1.5, cex = 0.8, bty = "n", ncol = 1, 
       title = "Domain", title.cex = 1.1, title.adj = 0, y.intersp = 1.1)
legend(x = -1.7, y = -1.1, c("Positive", "Negative"), lwd = 2,
       col = c("red", "blue"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.cex = 1.1, title.adj = 0)
title(main = "Spot Positive Connections (SPIEC-EASI)", adj = 0.5, line = -1, cex.main = 1.5)
mtext("Nodes = 200", side = 1, line = -1.5, cex = 0.8)
mtext("Edges = 164", side = 1, line = -0.5, cex = 0.8)
mtext("Mean Degrees = 1.6", side = 1, line = 0.5, cex = 0.8)
mtext("Clustering Coefficient = 0.093", side = 1, line = 1.5, cex = 0.8)
dev.off()


#### __Lichen ####
# Take the positive network and subset to just Fungi, Cyanobacteria, Algae
net_lich_spot <- delete.vertices(net_spot_pos, 
                                 V(net_spot_pos)[V(net_spot_pos)$taxonomy %notin% c("B_Cyanobacteria",
                                                                                    "E_Charophyta",
                                                                                    "E_Chlorophyta",
                                                                                    "E_Fungi")])

# Info for labels
length(V(net_lich_spot)) # 48
length(E(net_lich_spot)) # 6
deg_lich_spot <- degree(net_lich_spot, mode="all")
round(mean(deg_lich_spot), 1) # 0.3
round(transitivity(net_lich_spot), 3) # NA

# Plot Circle
png("Network_Lichen_Spot.png", width = 7, height = 6, units = "in", res = 300)
par(mar = c(4,8,2,0), xpd = TRUE) # bottom, left, top, right
plot(net_lich_spot,
     vertex.color = V(net_lich_spot)$groupcolor,
     vertex.size = (deg_lich_spot+1)*5,
     vertex.shape = V(net_lich_spot)$Domain,
     vertex.frame.color = "black",
     vertex.label = NA,
     #vertex.label = V(net_lich_spot)$ESV_ID2, # Check IDs. Make note of lichens
     edge.curved = 0.5,
     edge.width = 2,
     layout = layout_in_circle(net_lich_spot, order = order(V(net_lich_spot)$Group)))
legend(x = -1.7, y = 0.8, levels(spot_abund_combined$taxonomy_loaded$Group), 
       pch = c(rep(21, 12), rep(24, 12)), col = "black", pt.bg = c("green", "green4", "brown"), 
       pt.cex = 2, cex = 0.8, bty = "n", ncol = 1,
       title = "Group", title.cex = 1.1, title.adj = 0, y.intersp = 1.2)
legend(x = -1.7, y = 0.3, c("Bacteria", "Eukaryota"), pch = c(21, 24), 
       col = "black", pt.cex = 1.5, cex = 0.8, bty = "n", ncol = 1, 
       title = "Domain", title.cex = 1.1, title.adj = 0, y.intersp = 1.2)
legend(x = -1.7, y = -0.2, c("Positive", "Negative"), lwd = 2,
       col = c("red", "blue"), cex = 0.8, bty = "n", ncol = 1, 
       title = "Association", title.cex = 1.1, title.adj = 0)
title(main = "Spot Cyano-Algae-Fungi Positive Connections", adj = 0.5, line = -0.5, cex.main = 1.25)
mtext("Nodes = 48", side = 1, line = -0.75, cex = 0.8)
mtext("Edges = 6; Potential lichens = 2", side = 1, line = 0.25, cex = 0.8)
mtext("Mean Degrees = 0.2", side = 1, line = 1.25, cex = 0.8)
mtext("Clustering Coefficient = 0", side = 1, line = 2.25, cex = 0.8)
dev.off()

# Lichens are:
View(spot_abund_combined$taxonomy_loaded)
# ASV_Prok_63 - ASV_Euk_20
# Cyanobacteria Cyanobacteriia Gloeobacterales Gloeobacteraceae Gloeobacter PCC-7421 ESV_45
# Opisthokonta Nucletmycea Fungi Herpotrichiellaceae NA ESV_7

# ASV_Euk_27 - ASV_Euk_21
# Eukaryota Archaeplastida Chloroplastida Chlorophyta NA NA ESV_524
# Eukaryota Opisthokonta Nucletmycea Fungi NA NA ESV_68

# Note, the first one is the same as the whole dataset. The second one is new.



#### _Betweenness ####
# Plot degree versus betweenness for each network (each state)
se_spot <- se.mb2_spot
net_spot <- adj2igraph(getRefit(se_spot),  vertex.attr=list(name=taxa_names(spot.phy)))
bw_spot <- data.frame("Degree" = degree(net_spot, mode="all"),
                      "Betweenness" = betweenness(net_spot)) %>%
  mutate("ASV" = rownames(.)) %>%
  left_join(., spot_abund_combined$taxonomy_loaded, by = c("ASV" = "ESV_ID2"))

nb.cols <- length(levels(as.factor(bw_spot$Phylum)))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(bw_spot, aes(Degree, Betweenness, colour = Phylum)) +
  geom_jitter(shape = 16, size = 2, alpha = 1, width = 0.15) +
  labs(x = "Degree",
       y = "Betweenness",
       colour = "Phylum") +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(limits = c(0, 10),
                     breaks = c(0,1,2,3,4,5,6,7,8, 9, 10)) +
  theme_bw() +
  theme(legend.key.size = unit(0.1, "cm"))

# Check ASVs, add labels
ggplot(bw_spot, aes(Degree, Betweenness, colour = Phylum)) +
  geom_jitter(shape = 16, size = 2, alpha = 1, width = 0.15) +
  geom_text(data = bw_spot,
            aes(x = Degree, y = Betweenness, label = ASV), 
            size = 3, inherit.aes = F) +
  labs(x = "Degree",
       y = "Betweenness",
       colour = "Phylum") +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(limits = c(0, 10),
                     breaks = c(0,1,2,3,4,5,6,7,8, 9, 10)) +
  theme_bw() +
  theme(legend.key.size = unit(0.1, "cm"))

asvs_to_label_spot <- c("ASV_Euk_65", "ASV_Euk_1", "ASV_Euk_2", "ASV_Prok_95", "ASV_Prok_26")
bw_asvs_spot <- bw_spot %>%
  filter(ASV %in% asvs_to_label_spot) %>%
  filter(Betweenness > 1500) %>%
  mutate("HighTax" = ifelse(Genus != "NA",
                            Genus,
                            ifelse(Family != "NA",
                                   Family,
                                   ifelse(Order != "NA",
                                          Order,
                                          ifelse(Class != "NA",
                                                 Class,
                                                 ifelse(Phylum != "NA",
                                                        Phylum,
                                                        Domain)))))) %>%
  mutate("HighTax_ESV" = paste(HighTax, ESV_ID, sep = "_"))

png("Degree_Betweenness_Spot.png", width = 7, height = 3, units = "in", res = 300)
ggplot(bw_spot, aes(Degree, Betweenness, colour = Phylum)) +
  geom_jitter(shape = 16, size = 2, alpha = 1, 
              position = position_jitter(seed = 1, width = 0.15)) +
  geom_text_repel(data = bw_asvs_spot,
                  #min.segment.length = 0,
                  aes(x = Degree, y = Betweenness, label = HighTax_ESV), 
                  size = 2, inherit.aes = F,
                  position = position_jitter(seed = 1, width = 0.15)) +
  labs(x = "Degree",
       y = "Betweenness",
       colour = "Phylum") +
  scale_color_manual(values = mycolors) +
  scale_x_continuous(limits = c(0, 10),
                     breaks = c(0,1,2,3,4,5,6,7,8, 9, 10)) +
  theme_bw() +
  theme(legend.key.size = unit(0.1, "cm"))
dev.off()



#### _Participation ####
# Plot participation coefficient and within-module degree (Barnes et al.)
# Dashed lines at 0.61 and 2.2
# z-score (within module degree) is "connectivity"
# participation coefficient P is "participation
roles_spot <- role_spot %>%
  left_join(., spot_abund_combined$taxonomy_loaded, by = c("name" = "ESV_ID2")) %>%
  mutate(Phylum = as.factor(Phylum),
         role = as.factor(role)) %>%
  filter(is.na(module) == F) %>%
  droplevels()

nb.cols <- length(levels(roles_spot$Phylum))
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
ggplot(roles_spot, aes(participation, connectivity, colour = Phylum, shape = role)) +
  geom_vline(xintercept = 0.62, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_point(size = 2, alpha = 0.75) +
  labs(x = "Participation coefficient (P)",
       y = "Within module degree (z)",
       colour = "Phylum",
       shape = "Network Role") +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.key.size = unit(0.1, "cm"),
        panel.grid = element_blank())

# Check ASVs, add labels
ggplot(roles_spot, aes(participation, connectivity, colour = Phylum, shape = role)) +
  geom_vline(xintercept = 0.62, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text(data = roles_spot,
            aes(x = participation, y = connectivity, label = name), 
            size = 3, inherit.aes = F) +
  labs(x = "Participation coefficient (P)",
       y = "Within module degree (z)",
       colour = "Phylum",
       shape = "Network Role") +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.key.size = unit(0.1, "cm"),
        panel.grid = element_blank())

pz_asvs_spot <- roles_spot %>%
  filter(connectivity > 2.5 | participation > 0.62) %>%
  mutate("HighTax" = ifelse(Genus != "NA",
                            Genus,
                            ifelse(Family != "NA",
                                   Family,
                                   ifelse(Order != "NA",
                                          Order,
                                          Class))))

png("Participation_Connectivity_Spot.png", width = 7, height = 5, units = "in", res = 300)
ggplot(roles_spot, aes(participation, connectivity, colour = Phylum, shape = role)) +
  geom_vline(xintercept = 0.62, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_hline(yintercept = 2.5, linetype = "dashed", colour = "gray", linewidth = 0.25) +
  geom_point(size = 2, alpha = 0.75) +
  geom_text_repel(data = pz_asvs_spot,
                  min.segment.length = 0,
                  aes(x = participation, y = connectivity, label = HighTax), 
                  size = 2, inherit.aes = F) +
  labs(x = "Participation coefficient (P)",
       y = "Within module degree (z)",
       colour = "Phylum",
       shape = "Network Role") +
  scale_color_manual(values = mycolors) +
  theme_bw() +
  theme(legend.key.size = unit(0.1, "cm"),
        legend.margin = margin(0,0,0,0),
        panel.grid = element_blank())
dev.off()

keystone_spot <- roles_spot %>%
  filter(participation > 0.62 | connectivity > 2.5)

#write_xlsx(keystone_spot, format_headers = F, "keystone_taxa_spot.xlsx")



#### End Script ####