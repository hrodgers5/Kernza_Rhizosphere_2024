#This script tests if any taxa are differential abundant between groups using the LINDA package, then plots them
  #Hannah Rodgers, Feb 2024

#packages
require(phyloseq)
library(mia)
library(tidySummarizedExperiment)
library(MicrobiomeStat)
require(ggbreak)
require(patchwork)
require(tidyverse)

#load in data
load("../Data/rhizo_phyloseq.RData")

#set dataset as ps
ps <- ps_bac_noNA

#Optional: find % NAs and number of total groups at different taxonomic levels
  #PHYLUM, CLASS, ORDER, FAMILY, GENUS, SPECIES
taxa <- as.data.frame(tax_table(ps))
mean(is.na(taxa$GENUS)) #percent NAs
length(unique(taxa$FAMILY)) #number of total groups

#### LINDA TEST FOR DIFF REL ABUN ####
#Data preparation with TSE

## Transform PS object into a Mia object
tse <- makeTreeSummarizedExperimentFromPhyloseq(ps)

# Agglomerate by taxa of interest
tse_taxa <- subsetByPrevalentFeatures(tse, rank = "FAMILY")

##Test for Differential Relative Abundance with LinDA

res <- linda(
  feature.dat = as.data.frame(assay(tse_taxa)), 
  meta.dat = as.data.frame(colData(tse)), 
  formula = '~Breeding_Cycle', 
  alpha = 0.05, 
  prev.filter = 0, 
  mean.abund.filter = 0)

# save output
output2 <- as.data.frame(res$output)

#### GRAPH TAXA ####

#convert to relative abundance
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
colSums(otu_table(ps_rel))

#aggregate by taxa and convert to dataframe
ps_melt <- tax_glom(ps_rel, taxrank = "FAMILY", NArm = FALSE) %>% 
  psmelt()

#summarize (sum totals for each taxa/sample)
ps_melt2 <- ps_melt %>%
  group_by(FAMILY, Sample, Breeding_Cycle) %>%
  summarise(Abundance = sum(Abundance)) %>% 
  ungroup()

#FOR BAC PHYLA: select 10 most abundant taxa
  taxa_to_graph <- ps_melt2 %>%
    group_by(PHYLUM) %>% 
    summarise(Abundance = sum(Abundance)) %>% 
    arrange(dplyr::desc(Abundance)) %>% 
    dplyr::slice(1:10)
  
  taxa_to_graph <- taxa_to_graph$PHYLUM
  
  #rename all others "Other"
  ps_melt2 <- ps_melt2 %>% 
    mutate(PHYLUM = if_else(PHYLUM %in% taxa_to_graph, PHYLUM, "Other"))
  
#FOR FAMILIES SIG IMPACTED BY BREEDING: filter to just Glomeraceae, Kribbellaceae, or Pseudonocardiaceae
  ps_melt2 <- ps_melt2 %>% 
    filter(FAMILY %in% c("Glomeraceae", "Kribbellaceae", "Pseudonocardiaceae"))

   
#save bac and fungi data
   #phyla_bac <- ps_melt2
   #classes_amf <- ps_melt2
   families_bac <- ps_melt2
    families_amf <- ps_melt2
  
#merge families
  families_all <- bind_rows(families_bac, families_amf)

#GRAPH FAMILIES: line graph
  
label = parse(text = "p[adj]==0.047")

(kribb <- filter(families_all, FAMILY == "Kribbellaceae") %>% 
    
  ggplot(aes(x= Breeding_Cycle, y= Abundance)) +
  geom_point(aes(fill = as.factor(Breeding_Cycle)), 
             color = "black", pch = 21, size = 3.5) +
  stat_summary(fun.y=mean, colour="black", geom="line", size = 0.8) +
  
  scale_x_continuous(breaks = seq(0, 9, by = 1)) +
  theme_bw(base_size = 12) +
  theme(panel.background = element_rect(fill = '#F8F8F8'), 
        legend.position = "none") + 
  scale_fill_brewer(palette = "BrBG") +
    
  facet_wrap(~FAMILY) +
    
  annotate("label", x = Inf, y = Inf, label = label, size = 3.5, 
           vjust = 1.1, hjust = 1.1) +
  theme(axis.title.x=element_blank()) +
  theme(axis.title.y=element_blank()) +   
  labs(x= "Breeding Cycle", y = "Relative Abundance (%)"))

#patchwork and save
png(file = "Figure5.png", width = 11, height = 3.5, units = "in", res = 1000)

glom + pseudo + kribb + plot_annotation(tag_levels = "A")

dev.off()

#graph it: stream graph (first run streamgraph function below)
require(RColorBrewer)

#loess streamgraph

#graph for classes and phyla
classes_amf$CLASS <- fct_reorder(classes_amf$CLASS, classes_amf$Abundance, .desc = TRUE)

(class_amf <- 
  make_streamgraph_CLASS(classes_amf) +
  scale_x_continuous(breaks = seq(0, 9, by = 1), expand = c(0,0)) +
  scale_y_continuous(expand = c(0,0)) +
  theme_bw(base_size = 12) +
  scale_y_cut(breaks=c(0.01), which=c(2, 1), scales=c(1, 3)) +
  facet_wrap(~ "AM Fungal Classes") +
  theme(panel.background = element_rect(fill = '#F8F8F8'), axis.text.y=element_blank()) +  
  scale_fill_manual(values= rep(brewer.pal(11, "Set2"), times = 2), guide=guide_legend(reverse=T)) +
  labs(x= "Breeding Cycle", fill = "Class")
  )

#differences for AMF and bac: 

#patchwork

png(file = "Figure4.png", width = 10, height = 4, units = "in", res = 1000)

phyla_bac + class_amf + plot_annotation(tag_levels = "A") +  plot_layout(widths = c(1.15, 1))
  
dev.off()
  
#### FUNCTION TO MAKE STREAMGRAPH ####
#must change x and y variables in the code below

get_smooth_site <- function(site, data, span = 1.5, smooth.spline = FALSE){
  
  df <- data %>% filter(CLASS %in% site)
  date_min <- min(data$Breeding_Cycle)
  date_max <- max(data$Breeding_Cycle)
  if (!smooth.spline){
    pred_datenum <- seq(date_min, date_max)
    mod <- loess(Abundance ~ Breeding_Cycle, df, span = span)
    pred_y <- predict(mod, data.frame(Breeding_Cycle = pred_datenum, se = FALSE))
  }
  pred_date <- as.Date("1970-01-01") + ddays(pred_datenum)
  data.frame(CLASS = site,
             Breeding_Cycle = pred_datenum,
             Breeding_Cycle = pred_date,
             Abundance = pred_y)
}

if (FALSE){
  
  # Test for a single site
  test <- get_smooth_site("CS", data = testdata)
  head(test)
  ggplot(test, 
         aes(x = Breeding_Cycle, y = Abundance, color = fct_reorder(CLASS, Abundance))) + 
    geom_point()
  
  # For all sites: perform 'get_smooth_site' for each site (map_),
  #   then combine all results by stacking data frames (map_dfr)
  all_sites <- table(testdata$CLASS) %>% names() %>% rev()
  smooths_separate <- purrr::map_dfr(all_sites, 
                                     get_smooth_site, 
                                     data = testdata)
  ggplot(smooths_separate, aes(x = Breeding_Cycle, 
                               y = Abundance, 
                               fill = fct_reorder(CLASS, Abundance), 
                               colour = fct_reorder(CLASS, Abundance))) + 
    geom_col()
  
}


#
# get_smooths_cumul
#   Take smoothed values for all sites (result of 'purrr::map_dfr(sites, get_smooth_site)')
#   and make cumulative sum (in the order of 'groups)
# Output: data frame 
#

get_smooths_cumul <- function(smooths, sites){
  result <- smooths
  result$Abundance_lo <- 0   # for the first site
  for (i in 2:length(sites)){
    sel1 <- result$CLASS == sites[i-1]
    sel2 <- result$CLASS == sites[i]
    result$Abundance_lo[sel2] <- result$Abundance[sel1] 
    result$Abundance[sel2] <- result$Abundance[sel1] + result$Abundance[sel2] 
  }
  result
}

if (FALSE){
  
  # Test 
  smooths_cumul <- get_smooths_cumul(smooths_separate, all_sites)
  gg <- ggplot(smooths_cumul, aes(x = Breeding_Cycle, fill = fct_reorder(CLASS, Abundance))) + 
    geom_ribbon(aes(ymin = Abundance_lo, ymax = Abundance))
  gg
  gg + scale_y_log10()
  gg + geom_point(aes(y = Abundance))
  
  # The same function can be used to get the cumulative sums
  #   for the data as well which is useful for plotting:
  data_cumul <- get_smooths_cumul(testdata, all_sites)
  gg + geom_point(data = data_cumul, aes(y = Abundance))
  
  
}


get_data_cumul <- function(smooths, sites){
  result <- smooths
  result$Abundance_lo <- 0   # for the first site
  for (i in 2:length(sites)){
    sel1 <- result$CLASS == sites[i-1]
    sel2 <- result$CLASS == sites[i]
    result$Abundance_lo[sel2] <- result$Abundance[sel1] 
    result$Abundance[sel2] <- result$Abundance[sel1] + result$Abundance[sel2] 
  }
  result
}


#  
# Combine all steps
#

make_streamgraph_CLASS <- function(data,                      # must contain SITE_CODE, SAMPLE_DATE, WA_Avg
                             span = 0.8,                # lower number = more wiggly plot
                             logscale = FALSE,          # log10 y scale
                             brewer_palette = "RdYlGn", # color palette from colorbrewer2.org
                             brewer_direction = 1,      # color palette direction
                             add_datapoints = FALSE,    # adding points showing the data
                             add_databars = FALSE       # adding bars showing the data
){
  
  # Get sites (order in "opposite direction" )
  all_sites <- table(data$CLASS) %>% names() %>% rev()
  
  # Get smoothed curves for each site, combined in one data frame
  smooths_separate <- purrr::map_dfr(all_sites, 
                                     get_smooth_site, 
                                     data = data,
                                     span = span)
  
  # Get cumulative sums
  smooths_cumul <- get_smooths_cumul(smooths_separate, all_sites)
  
  # Make plot
  gg <- ggplot(smooths_cumul, aes(x = Breeding_Cycle, fill = fct_reorder(CLASS, Abundance))) + 
    geom_ribbon(aes(ymin = Abundance_lo, ymax = Abundance)) #+
  #scale_fill_brewer(palette = brewer_palette, direction = brewer_direction)
  
  if (logscale){
    gg + scale_y_log10()
  }
  
  if (add_datapoints){
    data_cumul <- get_smooths_cumul(data, all_sites)
    gg <- gg + geom_point(data = data_cumul, aes(y = Abundance))
  }
  
  if (add_databars){
    gg <- gg + 
      geom_col(data = data, aes(y = Abundance), color = "black", width = 2)
  }
  
  gg
  
}

if (FALSE){
  
  # Test
  make_streamgraph(testdata) 
  make_streamgraph(testdata, add_datapoints = TRUE)
  make_streamgraph(testdata, add_datapoints = TRUE, span = 1)
  
}

  
