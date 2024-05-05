#This script creates PCA plots on PLFA data from Wyoming rhizosphere study and tests for differences between groups.

library(readxl)
library(vegan)
library(patchwork)
library(ggbiplot)
library(ggrepel)
library(ggcorrplot)
library(tidyverse)

#### KANSAS PLFA ORDINATION ####

#load in data
plfa <- read_excel("Rhizo_PLFA_Kansas.xlsx", sheet = "raw_peaks")
groups <- read_excel("Rhizo_PLFA_Kansas.xlsx", sheet = "named_peaks")

#calculate PCA on raw biomarkers (first remove any columns with all zeros)
pca <- plfa[,3:69] %>% 
  select(where(~ sum(.) != 0)) %>% 
  prcomp(center = TRUE, scale. = TRUE)

summary(pca)

#use PERMANOVA (adonis testing) to test for sig differences between groups
adonis2(plfa[,c(3:69)] ~ Breeding_Cycle, data = plfa)

#test for sig differences in dispersion (distance from centroid) between groups
disp2 <- betadisper(dist(plfa[,c(3:69)]),
                    #grouping variable
                    plfa$Breeding_Cycle, 
                    type = 'centroid')

anova(disp2)
boxplot(disp2)

#nicer dispersion plot
dispersion <- data.frame(disp2$distances, disp2$group)


(disp <- ggplot(dispersion, aes(as.numeric(as.character(disp2.group)), 
                                disp2.distances, 
                                fill = disp2.group)) +
  geom_point(color = "black", pch = 21, size = 3.5) +
    
    theme_light(base_size = 12) +
    theme(panel.background = element_rect(fill = '#F8F8F8')) +
    scale_fill_brewer(palette = "BrBG") +
   scale_x_continuous(breaks = seq(0, 9, by = 1)) +
    theme(legend.position = "none") +
    
    labs (x = "Breeding Cycle", y = "PLFA Biomarker Dispersion\n(Distance to Centroid)") +
   annotate("label", x = Inf, y = Inf, label = "p=0.02", size = 3.5, hjust = 1.1, vjust =1.1))

#ORDINATION PLOT

#put PC1 and PC2 into dataframe
plfa[, c('PC1', 'PC2')] <- pca$x[, 1:2]

#use envfit to correlate other PLFA groupings to ordination
groups$"AM Fungi" <- groups$AMF

vars <- select(groups, Fungi, Bacteria, "AM Fungi")

en <- vegan::envfit(pca, vars, permutations = 10000, na.rm = TRUE, strata = NULL)

# use results from envfit to create arrows (continuous variables)
arrows = as.data.frame(scores(en, "vectors"))

#rescale arrows to fit nicely within scatterplot of our data
arrows$var = rownames(arrows)

mult = max(abs(plfa[, c('PC1', 'PC2')])) / max(abs(arrows[, 1:2])) / 2
arrows[, 1:2] = arrows[, 1:2] * mult

#### PLOT IT! Use PLFAs for scatterplot, and arrows for arrows ####
(ord <-
  
  ggplot(plfa, aes(x = PC1, y = PC2, 
             color = as.factor(Breeding_Cycle))) +
  
  geom_point(aes(fill = as.factor(Breeding_Cycle)), 
             colour = "black", pch = 21, size = 3.5) +
  
  scale_fill_brewer(palette = "BrBG") +
  
  geom_segment(data = arrows, 
               aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.02, "npc"))) +
  
  geom_label_repel(data = arrows, 
                   aes(PC1 * 1.9, PC2 * 1.2, label = var), 
                   inherit.aes = FALSE, 
                   nudge_x = -0.5, 
                   box.padding = 0.1) +
  theme_bw(base_size = 12) +
   theme(panel.background = element_rect(fill = '#F8F8F8')) +
   theme(legend.position = "none") +
  
  ylim(-4,2.5) +
  
  labs(fill = "Breeding Cycle") +
  annotate("label", x = Inf, y = Inf, label = "p<0.001", size = 4, hjust = 1.1, vjust =1.1))

#save plot
png(file = "FigureS1.png", width = 7, height = 3.6, units = "in", res = 500)
ord + disp + plot_annotation(tag_levels = "A")
dev.off()





#Line Graph of GAM of just PC1
mod <- gam(PC1 ~ s(Breeding_Cycle, bs = "cr", k=5), 
           family = gaussian (link = "identity"), 
           data = plfa)

sim <- simulateResiduals(mod, plot = T)

summary(mod)

#plot
mod_p <- predict_gam(mod)

ggplot() +
  geom_point(data = plfa, 
             aes(y = PC1, x = Breeding_Cycle, 
                 color = as.factor(Breeding_Cycle)),
             size = 4) +
  geom_smooth(data=mod_p, aes(Breeding_Cycle, PC1), color = "black") +
  geom_ribbon(data = mod_p, aes(x= Breeding_Cycle, 
                                ymax = upper_ci, 
                                ymin = lower_ci), 
              alpha = 0.15) +
  theme_bw(base_size = 15) +
  theme(legend.position = "none") +
  scale_color_brewer(palette = "BrBG")
