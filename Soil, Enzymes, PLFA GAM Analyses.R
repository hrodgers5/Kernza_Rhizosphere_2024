#Analyzing Enzyme Activities and Organic Matter Pools as a function of Breeding Cycle using GAMs
#Builds correlogram

#load packages
library(mgcv)
library(readxl)
library(tidygam)
library(DHARMa)
#library(mgcViz)
library(patchwork)
#library(ggfortify)
library(tidyverse)

#read in data
enzymes <- read_excel("Rhizo_Enzymes_Kansas.xlsx")
soil <- read_excel("Rhizo_Soil_Kansas.xlsx")
plfa <- read_excel("Rhizo_PLFA_Kansas.xlsx")

#### DATA PROCESSING AND MEANS ####
#create scaled mean for C-cycling enzymes, and center it on 5
enzymes$mean_C <- enzymes %>% 
  select(AG, BG, BX, CBH) %>% 
  scale() %>% 
  rowMeans()

enzymes$mean_C <- enzymes$mean_C + 5

#calculate microbial efficiency
enzymes$Carbon <- enzymes$mean_C/ plfa$Total_Micro
enzymes$Nitrogen <- enzymes$LAP/ plfa$Total_Micro
enzymes$Phos <- enzymes$PHOS/ plfa$Total_Micro
enzymes$Sul <- enzymes$SUL/ plfa$Total_Micro

#calculate means
soil %>% filter(Breeding_Cycle == 9) %>% 
  select(PMN) %>% 
  colMeans()


#### BUILD GAMs and LINE GRAPHS ####
  #use k=5 since we have n~50, which is about what the mgcv package automatically chooses anyway
  #Family is Gamma because data is continuous and non-negative. Use gamma for data with neg values.
  #NOTE: Usually use link = log with Gamma, but instead using link = identity so that I don't have to backtransform. This can lead to negative predicted values, so watch out for that

#use link= "log" for gamma, link= identity for gaussian (LAP)

mod <- gam(Total_Micro ~ s(Breeding_Cycle, bs = "cr", k=5), 
           family = Gamma (link = "log"), 
           data = plfa)

# MODEL CHECKING with Dharma package
#tests residuals, checks for KS, dispersion, and outliers
sim <- simulateResiduals(mod, plot = T)

# GET P VALUES
summary(mod)

# PREDICT 
mod_p <- tidygam::predict_gam(mod, length_out = 9)

#exponentiate to get back to normal scale (if link = log)
    mod_p[c(2,4,5)] <- lapply(mod_p[c(2,4,5)], exp)

#PLOT
    
stats <- expression(atop("p<0.001", paste(R^2==0.28)))

(total <- ggplot(plfa, aes(Breeding_Cycle, Total_Micro)) +
  geom_point(aes(fill = as.factor(Breeding_Cycle)), 
             color = "black", pch = 21, size = 3.5) +
  
  geom_smooth(data=mod_p, color = "black") +
  
  geom_ribbon(data = mod_p, aes(x= Breeding_Cycle, 
              ymax = upper_ci, ymin = lower_ci), alpha = 0.15) +
  
  theme_light(base_size = 12) +
  theme(panel.background = element_rect(fill = '#F8F8F8')) +
  scale_fill_brewer(palette = "BrBG") +
  scale_x_continuous(breaks = seq(0, 9, by = 1)) +
  theme(legend.position = "none") +
  #facet_wrap (~ "Carbon Cycling") +
  theme(strip.text = element_text(colour = 'black')) +
  labs (x = "Breeding Cycle", y = "Total Microbial Biomass (nmol/g)") +
  annotate("label", x = Inf, y = Inf, label = stats, size = 3.5, vjust = 1.1, hjust = 1.1)
  ) 


#patchwork and save plots

png(file = "Figure2.png", width = 8, height = 3.6, units = "in", res = 1000)

total + amf + plot_annotation(tag_levels = "A")

dev.off()
