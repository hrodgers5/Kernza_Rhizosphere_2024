#This script uses the microeco package to analyze microbial functions and correlate environmental variables to taxa abundance and to  functions
  #Hannah Rodgers, Feb 2024

#packages
require(microeco)
require(phyloseq)
require(mgcv)
require(tidygam)
require(DHARMa)
require(patchwork)
require(vegan)
require(tidyverse)

#load in ps object
load("../Data/rhizo_phyloseq.RData")

#choose which dataset to use
ps <- ps_bac_noNA

#### CREATE MICROTABLE OBJECT ####
#pull out otu, taxonomy, and metadata tables
otu <- data.frame(otu_table(ps))
tax <- data.frame(tax_table(ps))
meta <- data.frame(sample_data(ps))
  
#Clean up tax table
tax$KINGDOM <- tax$DOMAIN
tax <- dplyr::rename(tax, "Kingdom" = KINGDOM, "Phylum" = PHYLUM, 
              "Class" = CLASS, "Order" = ORDER, "Family" = FAMILY, 
              "Genus" = GENUS, "Species" = SPECIES)


#create microtable object
micro <- microtable$new(otu, sample_table = meta, tax_table = tax)

#### ASSIGN FUNCTIONS ####
#create trans_func object
tf <- trans_func$new(micro)

#assign functions to OTUs using FAPROTAX or FUNGuilddatabase and view first few rows
  #for bacteria:
    tf$cal_spe_func(prok_database = "FAPROTAX")
  
    #for fungi:
    #tf$cal_spe_func(fungi_database = "FUNGuild")

#calculate percentages of each functional group for samples and view first few rows. 
    #abundance_weighted means use abundance of individuals, not taxa
tf$cal_spe_func_perc(abundance_weighted = TRUE)
tf$res_spe_func_perc[1:5, 1:2]

#save function percentages, and add in metadata
functions1 <- tf$res_spe_func_perc[]
meta$sampname <- rownames(meta)
functions1$sampname <- rownames(functions1)

meta2 <- select(meta, Breeding_Cycle, sampname)

functions1 <- functions1 %>% left_join(meta2, by = "sampname")

#run ordination/ permanova on all functions
pca <- functions1[,1:41] %>% 
  select(where(~ sum(.) != 0)) %>% 
  prcomp(center = TRUE, scale. = TRUE)

summary(pca)

#use PERMANOVA (adonis testing) to test for sig differences between groups
adonis2(functions1[,1:41] ~ Breeding_Cycle, data = functions1)

ord <- metaMDS(functions1[,1:41])

plot(ord)

require(ggbiplot)
pca <- prcomp(functions1[,1:41], scale. = TRUE)

ggbiplot(pca, obs.scale = 1, var.scale = 1,
         groups = as.factor(functions1$Breeding_Cycle), ellipse = TRUE, circle = TRUE) +
  scale_color_discrete(name = '') +
  theme(legend.direction = 'horizontal', legend.position = 'top')

#### PLOT FUNCTIONS WITH GGPLOT ####

#tidy that dataset, include only the functions we're interested in
functions2 <- functions1 %>% 
  dplyr::select("hydrocarbon_degradation", "methanotrophy", "nitrate_reduction", 
                         "nitrification", "nitrogen_fixation", "ureolysis", "Breeding_Cycle")
  
## run GAM analysis on functions 
mod <- gam(ureolysis ~ s(Breeding_Cycle, bs = "cr", k=5), 
           family = gaussian (link = "log"), 
           data = functions2)
#check mod
sim <- simulateResiduals(mod, plot = T)

#get p value
summary(mod)

# PREDICT 
mod_p <- tidygam::predict_gam(mod, length_out = 9)

#exponentiate to get back to normal scale (if link = log)
mod_p[c(2,4,5)] <- lapply(mod_p[c(2,4,5)], exp)

#PLOT
functions2$title <- "Ureolysis"
label <- "NS"

(ureolysis <- ggplot(functions2, aes(Breeding_Cycle, ureolysis)) +
    geom_point(aes(fill = as.factor(Breeding_Cycle)), 
               color = "black", pch = 21, size = 3.5) +
    
    geom_smooth(data=mod_p, color = "black") +
    
    geom_ribbon(data = mod_p, aes(x= Breeding_Cycle, 
                ymax = upper_ci, ymin = lower_ci), alpha = 0.15) +
    
    theme_bw(base_size = 12) +
    theme(panel.background = element_rect(fill = '#F8F8F8')) +
    scale_fill_brewer(palette = "BrBG") +
    scale_x_continuous(breaks = seq(0, 9, by = 1)) +
    facet_grid(. ~ title) +
    theme(legend.position = "none") +
    
    labs (x = "Breeding Cycle", y = "% of Total Bacteria") +
    theme(strip.background =element_rect(fill="darkolivegreen3")) +
    
    #theme(axis.title.x = element_blank()) +
    theme (axis.title.y = element_blank()) +
    annotate("label", x = Inf, y = Inf, label = label, 
             size = 3.5, vjust = 1.1, hjust = 1.1))

#colors: darkolivegreen3, tan3
#pathwork

png(file = "Figure7.png", width = 10, height = 5.5, units = "in", res = 1000)

hydrocarbon + methanotrophy + nitrogen + nitrification + nitrate + ureolysis + plot_annotation(tag_levels = "A")

dev.off()

#### CORRELATE FUNCTIONS TO ENVIRONMENTAL VARIABLES ####

#choose environmental variables
t3 <- trans_env$new(dataset = micro, env_cols = c(4,8,19,24,26,29:30))

#calculate correlation, choose which functions to include
t3$cal_cor(add_abund_table = tf$res_spe_func_perc[
    c("hydrocarbon_degradation", "methanotrophy", "nitrate_reduction", 
    "nitrification", "nitrogen_fixation", "ureolysis")], 
           cor_method = "spearman", p_adjust_method = "fdr", p_adjust_type = "Env")

#plot
t3$plot_cor(text_x_order = c("AMF_Biomass","Microbial_Biomass", 
                             "C.Cycling_Enzyme_Activity", "N.Cycling_Enzyme_Activity", 
                             "DOC", "DON", "PMN"))

#### CORRELATE TAXA TO ENVIRONMENTAL VARIABLES ####
#correlations between specific taxa and env variables
t1 <- trans_env$new(dataset = micro, add_data = meta[, c(4,8,19,24,26,29:30)])

micro$cal_abund(
  select_cols = NULL,
  rel = TRUE,
  merge_by = "|",
  split_group = FALSE,
  split_by = "&&",
  split_column = NULL
)


# 'p_adjust_type = "Env"' means p adjustment is performed for each environmental variable separately.

t1$cal_cor(use_data = "Family", p_adjust_method = "fdr", p_adjust_type = "Env")
t1$plot_cor()


#### DISSIMILARITY MATRIX / DISTANCE DECAY ####
#here, use the Mantel test to test for signficiant correlations between two dissimilarity matrices, instead of between two variables

#transform abundance data
ps_tr <- transform_sample_counts(ps, function(x) x / sum(x))
OTU<-t(data.frame(ps_tr@otu_table@.Data))

#pull out metadata
metadata <- data.frame(sample_data(ps_tr))

#select breeding cycle and environmental variables

BC <- select(metadata, Breeding_Cycle)
env <- select(metadata, AG, BG, BX, CBH, PHOS, SUL, N.Cycling_Enzyme_Activity, 
              DOC, DON, PMN, Protein, POXC, TotalN, SOC,
              Microbial_Biomass)

#convert each dataframe into distance matrices
#for abundance data, use Bray-Curtis similarity
otu.dist <- vegdist(OTU, method = "bray")

#for environmental data, use euclidean distance
bc.dist <- vegdist(BC, method = "euclidean")

env.dist <- vegdist(env, method = "euclidean")

#run the mantel test to look for sig correlations
mantel(otu.dist, bc.dist, method = "spearman", permutations = 9999, na.rm = TRUE)

mantel(otu.dist, env.dist, method = "spearman", permutations = 9999, na.rm = TRUE)

#new data frame with vectorized distance matrices
otu <- as.vector(otu.dist)
env <- as.vector(env.dist)
bc <- as.vector(bc.dist)

mat <- data.frame(otu, env, bc)

#GRAPHING

ggplot(mat, aes(y = otu, x = bc)) + 
  geom_point(size = 3, alpha = 0.9, shape = 21, 
             color = "black", aes(fill = bc)) +
  geom_smooth (method = "lm", color = "black") +
  labs(x = "Difference in Breeding Cycle", 
       y = "AMF Bray-Curtis Dissimilarity", 
       fill = "Difference in \nBreeding Cycle") +
  theme_bw() +
  scale_fill_distiller(palette = "PuBu")

