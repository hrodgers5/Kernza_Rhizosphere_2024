# Load packages
require(tidyverse)
require(phyloseq)
require(Biostrings)
require(ape)
require(vegan)
require(mgcv)

#load environment
load("../Data/rhizo_phyloseq.RData")

#### MAKE PHYLOGENETIC TREE ####

#save data as ps
ps <- ps_AMF_noNA_pruned

#### save OTU table as df #
ESV<-t(as.data.frame(otu_table(ps)))

#add centroid = so that colnames match with fasta file
colnames(ESV) <- paste("centroid=", colnames(ESV), sep = "")

#### subset fasta file so that it has the same sequences as OTU table

#in local R, subset fasta file by the ESV dataframe

seqs<-Biostrings::readDNAStringSet("../Data/AMF_Data/zotus_nonchimeric.fa")
seqs@ranges@NAMES<-stringr::str_split(seqs@ranges@NAMES,pattern = ";", 3, simplify = T)[,1]

seqs_sub<-seqs[colnames(ESV)]
table(seqs_sub@ranges@NAMES %in% colnames(ESV))
table(colnames(ESV)%in%seqs_sub@ranges@NAMES)
rm(list=ls()[! ls() %in% c("ESV", "seqs_sub" )])

#save the pruned fasta file to run on Beartooth
writeXStringSet(seqs_sub, filepath = "ESVs_forBeartooth_AMF.fasta", format = "fasta")

#### BEARTOOTH STUFF ####
#FIRST, copy the fasta file over to the supercomputer using the code below on MobaXTerm
scp /drives/c/Users/hanna/Desktop/ESVs_forBeartooth.fasta hrodger3@beartooth.arcc.uwyo.edu:/project/microbiome/users/hrodger3/phylotree/
  
#NEXT, create a run the following bash script on Beartooth to use clustalo to align and then fasttree to make the tree
  
  align_maketree.sh
```{bash}
#!/bin/bash
#SBATCH --job-name phylotree
#SBATCH --mem=120GB
#SBATCH --time=6-00:00:00
#SBATCH --cpus-per-task=6
#SBATCH --account=microbiome
#SBATCH --output=phylotree_%A.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hrodger3@uwyo.edu

module load arcc/1.0 miniconda3/4.12.0
conda activate shotgun_env
cd /project/microbiome/users/hrodger3/phylotree

clustalo -i ESVs_Beartooth.fasta -o clustalo_out_16s.fa -v --threads=6

FastTree -nt clustalo_out_16s.fa > 16s_phylotree.nwk

```

#Finally, copy nwk file from beartooth to computer with following code:
#scp hrodger3@beartooth.arcc.uwyo.edu:/project/microbiome/users/hrodger3/phylotree/16s_phylotree.nwk /drives/c/Users/hanna/Desktop/
  
#### RUN UNIFRAC ####
#add in tree file (nwk) to ps object
tree <- ape::read.tree("../Data/AMF_Data/AMF_phylotree2.nwk")

#remove centroid =
taxa_names(tree) <- str_sub(taxa_names(tree), 10)
  
ps <- merge_phyloseq(ps, tree)

#run UniFrac analysis
require("doParallel")
require("foreach")
require("picante")

#plot a tree for a subset of the data
subset <- prune_taxa(taxa_names(ps)[1:1000], ps)
plot_tree(subset, color = "Breeding_Cycle")

#run Faith's phylogenetic diversity analysis. PD is faith's phylo div and SR is species richness
samp <- t(data.frame(otu_table(ps)))
tree <- phy_tree(ps)
meta <- data.frame(sample_data(ps))

pd <- pd(samp, tree, include.root=FALSE)

#merge with metadata
pd2 <- merge.data.frame(pd, meta, by = 'row.names')

#run gam
mod <- gam(PD ~ s(Breeding_Cycle, bs = "cr", k=5), 
           family = Gamma (link = "log"), 
           data = pd2)

summary(mod)

# PREDICT 
mod_p <- tidygam::predict_gam(mod, length_out = 9)

#exponentiate to get back to normal scale (if link = log)
mod_p[c(2,4,5)] <- lapply(mod_p[c(2,4,5)], exp)


ggplot(pd2, aes(x=Breeding_Cycle, y=PD)) +
  geom_point()







#run UniFrac analysis. Before doing unweighted, can rarefy
ps_rarefy <- rarefy_even_depth(ps, rngseed = 1992)

UF <- UniFrac(ps_rarefy, weighted=FALSE, normalized=TRUE, parallel=FALSE, fast=TRUE)

#run PERMANOVA
meta <- data.frame(sample_data(ps))

adonis2(UF ~ Breeding_Cycle, data = meta)

#use ordinate function to perform weighted UniFrac(which considers relative abundance of taxa, NOT just presence/absence) and then run PCA on that distance matrix
ordu = ordinate(ps, "PCoA", "unifrac", weighted=FALSE)

plot_ordination(ps, ordu, color= "Breeding_Group2") +
  geom_point(size = 5) +
  stat_ellipse()
