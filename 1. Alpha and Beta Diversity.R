#This script analyzes alpha and beta diversity, and plots an ordination with environmental variables overlaid
  #Hannah Rodgers, Feb 2024

#packages
require(ggrepel)
require(phyloseq)
library(mgcv)
library(DHARMa)
require(vegan)
require(patchwork)
require(tidyverse)

#load in data
load("../Data/rhizo_phyloseq.RData")

#set which data you're working with as ps
ps <- ps_AMF_noNA_pruned

#### ALPHA DIVERSITY ####

#rarefy and check read depth
sort(colSums(otu_table(ps)))

ps_rarefy <- rarefy_even_depth(ps, rngseed = 1992)
sort(colSums(otu_table(ps_rarefy)))

# Calculate diversity metrics and then add in metadata
metadata <- data.frame(sample_data(ps))

Richness <- estimate_richness(ps_rarefy, split = TRUE, measures = NULL)
Richness$Breeding_Cycle <- metadata$Breeding_Cycle

# Plot all alpha diversity measures
plot_richness(ps_rarefy, "Breeding_Cycle") + geom_smooth() + geom_point()

#Run a GAM on observed diversity
mod <- gam(Simpson ~ s(Breeding_Cycle, bs = "cr", k=5), 
           family = gaussian (link = "identity"), 
           data = Richness)

# model checking with Dharma package
sim <- simulateResiduals(mod, plot = T)

# GET P VALUES
summary(mod)

# Predict values
mod_p <- tidygam::predict_gam(mod, length_out = 9)

#exponentiate to get back to normal scale #(only necessary if link = log or sqrt in GAM)
mod_p[c(2,4,5)] <- lapply(mod_p[c(2,4,5)], exp)
mod_p[c(2,4,5)] <- lapply(mod_p[c(2,4,5)], function(x) x^2)

#PLOT
(glom <- ggplot(Richness, aes(Breeding_Cycle, Observed)) +
    geom_point(aes(fill = as.factor(Breeding_Cycle)), 
               color = "black", pch=21, size = 4) +
  
   geom_smooth(data=mod_p, color = "black") +
  
   geom_ribbon(data = mod_p, aes(x= Breeding_Cycle, 
               ymax = upper_ci, ymin = lower_ci), alpha = 0.15) +
  
    theme_bw(base_size = 12) +
    theme(panel.background = element_rect(fill = '#F8F8F8')) +
    scale_x_continuous(breaks = seq(0, 9, by = 1)) +
    scale_fill_brewer(palette = "BrBG") +
    theme(legend.position = "none") +

    labs (x = "Breeding Cycle", y = "Number of OTUs", 
          title = "Glomeromyota Diversity")) #+
    #annotate("label", x = 8, y = 1200, label = "NS ", size = 4))

#patchwork
bac + fungi + glom + plot_annotation(tag_levels = "A")

#### Beta Diversity (Permanova and Ordination) ####
otu <- data.frame(otu_table(ps))
metadata <- data.frame(sample_data(ps))

# Calculate Bray Curtis distance matrix
adonis_dist <- vegdist(data.matrix(t(otu)), method="jaccard")

#calculate PERMANOVA
adonis2(adonis_dist ~ Breeding_Cycle, data = metadata, permutations = 50000)

#calculate and visualize dispersal
mod <- betadisper(adonis_dist, metadata$Breeding_Cycle, 
                  type = "centroid")
anova(mod)

### Plot Ordination with Environmental Variables Correlated ####

#perform NMDS. 
# if you don't set seed, randomization will yield slightly different results every time
set.seed(05261940)
nmds <- metaMDS(adonis_dist, distance = "bray") 
nmds # good stress value is below 0.2

# save NMDS scores
data.scores <- as.data.frame(scores(nmds$points))
  
#merge metadata with data.scores
data.scores$Sample_ID <- metadata$Sample_ID

  data.scores <- data.scores %>% 
  dplyr::left_join(metadata, by = 'Sample_ID')

#pull out NMDS scores
scores <- dplyr::select(data.scores, MDS1, MDS2) 

#select environmental variables to correlate to ordination. Can rename some
  # all variables to try: C.Cycling_Enzyme_Activity, PHOS, SUL, DOC, DON, PMN, Microbial_Biomass, AMF_Biomass
vars <- dplyr::select(data.scores, C.Cycling_Enzyme_Activity, SUL, DON, PHOS)

vars <- dplyr::rename(vars, "C Enzymes" = C.Cycling_Enzyme_Activity)

#use envfit to correlate environmental variables to ordination and check for significance
(en <- envfit(scores, vars, permutations = 10000))

# use results from envfit to create arrows (continuous variables)
arrows = as.data.frame(scores(en, "vectors"))

#rescale arrows to fit nicely within scatterplot of our data
arrows$var = rownames(arrows)
mult = max(abs(scores[, c('MDS1', 'MDS2')])) / max(abs(arrows[, 1:2])) / 2
arrows[, 1:2] = arrows[, 1:2] * mult

# PLOT ORDINATION
#arrows_amf <- arrows
#data.scores_bac <- data.scores

(ord_amf <- ggplot(data.scores_amf, aes(x = MDS1, y = MDS2)) + 
  
  geom_point(aes(fill = as.factor(Breeding_Cycle)), 
             color= "black", pch=21, size = 3.5) +
  geom_segment(data = arrows, 
               aes(x = 0, y = 0, xend = MDS1, yend = MDS2), 
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.02, "npc"))) +
  
  geom_label_repel(data = arrows, aes(x = MDS1, y = MDS2), 
                   colour = "black", label = row.names(arrows), 
                  nudge_x = 0.4, size = 3) +
  theme_bw(base_size = 12) +
  #theme(panel.background = element_rect(fill = 'gray95')) +
  scale_fill_brewer(palette = "BrBG") +
  
  guides(fill="none") +
  annotate("label", x = 0.8, y = -0.55, label = "p = 0.049", size = 3.5))

  #plot differences
      #AMF: x = 0.8, y = -0.6, label = "p = 0.049", fill = 'gray95'
      #bacteria: x = 0.8, y = 0.35, label = "p = 0.02"

## plot MDS1 above ordination ##

(scatter_amf <- ggplot(data.scores, aes(Breeding_Cycle, MDS1)) +
  
  geom_point(aes(fill = as.factor(Breeding_Cycle)), 
             color = "black", pch = 21, size = 3.5) +
    
  theme_bw(base_size = 12) +
 # theme(panel.background = element_rect(fill = '#F8F8F8')) +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 9, by = 1)) +
  scale_fill_brewer(palette = "BrBG") +
  facet_wrap(~ 'AM Fungi') +
  labs(x= "Breeding Cycle")) 


#patchwork ordination and scatterplot
png(file = "ordinations.png", width = 8, height = 6, units = "in", res = 1000)

(scatter_bac + scatter_amf)/(ord_bac + ord_amf) + plot_layout(heights = c(1,2.5)) + plot_annotation(tag_levels = "A")

dev.off()
