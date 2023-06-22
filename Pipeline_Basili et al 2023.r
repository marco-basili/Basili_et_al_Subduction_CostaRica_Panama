
## Dada2 read processing
library(dada2); packageVersion("dada2")

path <- "~/marco/merging/BMS_16S_fastq/16S_BMS_LLO/tot_database/"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_L00"), `[`, 1)
 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(245,200),
            maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, trimLeft=23,
            compress=TRUE, multithread=FALSE)
    head(out)


names(filtFs) <- sample.names
names(filtRs) <- sample.names
set.seed(100)
# Set the seed of R‘s random number generator, which is useful for creating simulations or random 
# objects that can be reproduced.

mergers_quar <- vector("list", length(sample.names))
names(mergers_quar) <- sample.names
for(sam in sample.names) {
  cat("Processing:", sam, "\n")
    errF_quar <- learnErrors(filtFs[[sam]], nbases=1e8, multithread=TRUE)
    derepF_quar <- derepFastq(filtFs[[sam]])
    ddF_quar <- dada(derepF_quar, err=errF_quar, multithread=TRUE)
    errR_quar <- learnErrors(filtRs[[sam]], nbases=1e8, multithread=TRUE)
    derepR_quar <- derepFastq(filtRs[[sam]])
    ddR_quar <- dada(derepR_quar, err=errR_quar, multithread=TRUE)
    merger_quar <- mergePairs(ddF_quar, derepF_quar, ddR_quar, derepR_quar)
    mergers_quar[[sam]] <- merger_quar
}

rm(derepF_quar); rm(derepR_quar); rm(ddR_quar); rm(ddF_quar); rm(errF_quar); rm(errR_quar) 

seqtab_quar <- makeSequenceTable(mergers_quar)
    
table(nchar(getSequences(seqtab_quar)))
seqtab_quar <- seqtab_quar[,nchar(colnames(seqtab_quar)) %in% seq(360, 380)]

seqtab_nochim_quar <- removeBimeraDenovo(seqtab_quar, method="consensus", multithread=TRUE)
seqtab_quar_t <- t(seqtab_nochim_quar)

getN <- function(x) sum(getUniques(x))
track_quar <- cbind(out, sapply(mergers_quar, getN), rowSums(seqtab_nochim_quar))
colnames(track_quar) <- c("input", "filtered", "merged", "nonchim")
rownames(track_quar) <- sample.names
head(track_quar)

taxa <- assignTaxonomy(seqtab_nochim_quar, "~/16S_BMS_LLO/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "~/16S_BMS_LLO/silva_species_assignment_v132.fa.gz")


# Load the packages
library(tidyverse)
library(sp)
library(sf)
library(rnaturalearth)
library(tidyr)
library(rgeos)
library(ggrepel)
library(ggplot2)
library("rnaturalearthdata")
library(ggspatial)
library(gghighlight)
library(ggsflabel)
library(pheatmap) 
library(rstatix)
#library(ggbiplot)
library("ggfortify")
library(phyloseq)
library(vegan)
library(geojsonio)
library(igraph)
library(repr)

library(gridExtra)
library(grid)
library(lattice)
library(cowplot) # for get_legend

library(ggthemes) # additional themes fro ggplot2
library(ggplotify)
library(factoextra)
library(microbiome)
library(ggplot2)
library(ggdendro)
library(gtable) # for align the plot
library(dplyr)

#visualizzazione barplot
library(wesanderson)
library(ggplot2)

#  creare barplot con ggplot
library("dplyr")
library("magrittr")
library("plotly")
library("ggplot2")
library("tidyr")
library("cowplot")
library(ggplotify)
library(pals)

#  devtools::install_github("nefosl/ggtern")
#  Viridis color library

# other structure and plot
library(RColorBrewer) # nice color options
library(dplyr) # data handling
library(network) # networks
library(intergraph)  # networks
library(GGally)   # network plotting with ggplot
library(igraph)  # networks
library(gridExtra) # gridding plots
library(grid)
library(ape) # importing and handling phylogenetic trees

library(magrittr) #
library(rioja) # plotting packages for tabular bubbleplots (inkspot)
library(ggpubr)
# library(ggtern) # ternary plots for geochemistry
library(plyr)
library(coda.base)
library(tidyverse)
library(propr)
library(missForest) # Imputing missing values in dataframes using Random Forests
library(VSURF)
library(ggnet)   # network plotting with ggplot
library(igraph) # networks
library(graphics)
library(microbiomeutilities)
library(ggthemes)

# set size of the plot
options(jupyter.plot_scale=1.5)
library(viridis)
# plot size
options(repr.plot.width=10, repr.plot.height=20)

# laod the sample table
coord <- read.csv("~/Desktop/Merging_2017_2018/paper_17_18/File_usati_per_R/coordinate.csv")
row.names(coord)
# Create data frame of only longitude and latitude values
coords <- select(coord, longitude, latitude)

# Create sf object with geo_data data frame and CRS
points_sf <- st_as_sf(coord, coords = c("longitude", "latitude"), crs = 4326, agr = "Altitude")
points_sf

# Create sf object with coast and country data
coast_sf <- ne_coastline(scale = "medium", returnclass = "sf")
countries_sf <- ne_countries(scale = "medium", returnclass = "sf")


#########    Figure  1 
#########    Suppl. Figure 1

 #  Map with sample point colored in accordance to variables 
ggplot() + 
  geom_sf(data = countries_sf) + 
  geom_sf(data = points_sf, 
          aes(fill =rc_ra, size = DSub ,stroke = 0.8), shape = 21, 
          alpha = 0.9,
          show.legend = "point") +
    scale_fill_viridis() +
  #  geom_text_repel(data = coord, size = 2.5, 
  #               aes(x = longitude, y = latitude, label = code)) +
  coord_sf(xlim = c(-86, -79), ylim = c(7, 12),expand = FALSE,
           datum = NA) + # removes graticules
  # scale_size(range = c(4, 15))+
  #geom_text(data = states, aes(X, Y, label = ID), size = 5) +
  labs(title = "Subduction zone",
       x = NULL,
       y = NULL) +
  annotation_scale(location = "bl", width_hint = 0.4) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                         style = north_arrow_fancy_orienteering) +
  theme(panel.grid.major = element_line(color = gray(0.5), linetype = "dashed", 
                                        size = 0.1), panel.background = element_rect(fill = "aliceblue"))+
  guides(color = guide_legend(override.aes = list(size = 6))) 
#  ggsave('map_tas.png', width = 6, height= 4)

###################################

# Split variables in categorical groups 

coord_geophisic <- coord[,c(113:132)]
coord_geochemic_major <- coord[,c(70:79)]
coord_geochemic_trace<- coord[,c(80:112)]
coord_environment<-coord[,c(31:60)]
coord[c(7:9, 31:60,  70:132)]
coord_paper<-coord[,c(31:60,113:132)]
  
colnames(coord_geophisic)
colnames(coord_geochemic_major)
colnames(coord_geochemic_trace)
colnames(coord_environment)
colnames(coord_paper)
colnames(coord)

cor_env_matrix2 <- cor(coord_geophisic, use="pairwise.complete.obs")
pheatmap(cor_env_matrix2,
         border_color = "grey60",
         cellwidth = 10, cellheight = 5,
         clustering_distance_rows = "euclidean",  
         clustering_distance_cols = "euclidean", 
         clustering_method = "average", 
         fontsize_row = 5,
         fontsize_col = 5)
         
         
################################################

###########       Figure  1 -  PCA

pc <- prcomp(na.omit(coord[,c(113:132)]),
             center = TRUE,
            scale. = TRUE)

# remove unwanted samples (47) because na.omit removes the entire row with NA

groups <- as.factor(coord[-47,13])
shape <-  as.factor(coord[-47,20])

fviz_pca_biplot(pc, obs.scale = 1, var.scale = 1, alpha=0.8,
                label = "var",
                col.var = "contrib") + 
  geom_point(aes(shape = factor(shape), fill = factor(groups), size = 4)) +
   scale_shape_manual(values=c(21,22,23, 24))+
   scale_fill_viridis(discrete=TRUE)+
  guides(shape = guide_legend(title = "Sub"),
         fill = guide_legend(title = "Area"))
         
##########################################

# correlation between variable and Groups and statistical analisys (wilcox)
library(ggpubr)
library(rstatix)
library(ggprism)
library(patchwork)
library(magrittr)

#########  Suppl. figure 2

p_cor <- ggplot(coord, aes(location3,DSub))+
	geom_boxplot(aes(fill=location3))+
	scale_fill_viridis(discrete=TRUE) +
	theme(axis.text.x = element_text(angle = 45, hjust=1))+
	theme(legend.position = "none",
	panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
	panel.background=element_rect(fill="white"),
	panel.border = element_rect(fill = NA, colour = rgb(100, 100, 100, maxColorValue = 255)))

# ggsave('box1.svg', width = 6, height= 4)
df_p_val <- rstatix::wilcox_test(coord,DSub ~ location3) %>% 
	rstatix::add_xy_position()
p_cor <-p_cor+ stat_pvalue_manual(df_p_val, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE)
p_cor
     
     
#######   Figure 1         
#######   Suppl. Figure 1
# Scatterplot 

ggplot(coord, aes(DSub, rc_ra, fill=location3))+
  #geom_text_repel(data = coord, aes(label = code), size = 3, max.overlaps = 40) +
  geom_point(shape = 21, size = 4, alpha = 0.9)+
  scale_fill_viridis(discrete=T)+
  scale_size_continuous(range = c(0.1, 15))+
  #geom_smooth(method="lm", formula= (y~log(-x)), se=FALSE, color=2)+ 
  theme(legend.position = "right") +
  theme_bw()+
  scale_shape_manual(values=c(21,22,23,24,25)) 

cor.test(coord$rc_ra, coord$DSub, method = "pearson")




#######################################################
#####    TAS Diagram - Figure 2
  
p <- ggplot(data=coord, mapping=aes(x=sio2_d,y=na2o_d+k2o_d,color=tectonics4)) +
    geom_blank() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_y_continuous(limits=c(0,15), expand = c(0, 0)) + 
    scale_x_continuous(limits=c(40,80), expand = c(0, 0)) +
    labs(y=expression(Na[2]*O + K[2]*O*~ wt~'%'), x=expression(SiO[2]*~ wt~'%'))+
    annotate("segment", x=45, xend=45, y=1, yend=5)+
    annotate("segment", x=45, xend=52, y=5, yend=5)+
    annotate("segment", x=52, xend=69, y=5, yend=8)+
    annotate("segment", x=76.5, xend=69, y=1, yend=8)+
    annotate("segment", x=69, xend=69, y=8, yend=13)+
    annotate("segment", x=45, xend=61.32, y=5, yend=13.7)+
    annotate("segment", x=52, xend=52, y=1, yend=5)+
    annotate("segment", x=57, xend=57, y=1, yend=5.9)+
    annotate("segment", x=63, xend=63, y=1, yend=6.9)+
    annotate("segment", x=52, xend=49.4, y=5, yend=7.3)+
    annotate("segment", x=57, xend=53.05, y=5.9, yend=9.25)+
    annotate("segment", x=63, xend=57.6, y=6.9, yend=11.7)+
    annotate("segment", x=41, xend=45, y=3, yend=3)+
    annotate("segment", x=41, xend=41, y=1, yend=3)+
    annotate("segment", x=41, xend=41, y=3, yend=7, linetype="dashed")+
    annotate("segment", x=41, xend=45, y=7, yend=9.4, linetype="dashed")+
    annotate("segment", x=45, xend=52.5, y=9.4, yend=14)+
    annotate("segment", x=49.4, xend=45, y=7.3, yend=9.4)+
    annotate("segment", x=53, xend=48.4, y=9.3, yend=11.5)+
    annotate("segment", x=57.6, xend=50.3, y=11.7, yend=15)
  
tas <- p + annotate("text", label = "Basalt", x = 48.5, y = 1.5, size=3)+
    annotate("text", label = "Basaltic\n andesite", x = 54.8, y = 1.5, size=3)+
    annotate("text", label = "Andesite", x = 60, y = 1.5, size=3)+
    annotate("text", label = "Dacite", x = 67.5, y = 1.5, size=3)+
    annotate("text", label = "Rhyolite", x = 72, y = 8, size=3)+
    annotate("text", label = "Trachy- \n basalt", x = 47.8, y = 5.7, size=3)+
    annotate("text", label = "Basaltic \n trachy- \n andesite", x = 52.5, y = 7, size=3)+
    annotate("text", label = "Trachy- \n andesite", x = 57.8, y = 8.2, size=3)+
    annotate("text", label = "Trachydacite", x = 65, y = 9, size=3)+
    annotate("text", label = "Trachyte", x = 62.5, y = 11.5, size=3)+
    annotate("text", label = "Picro- \n basalt", x = 43, y = 1.5, size=3)+
    annotate("text", label = "Basanite ", x = 44, y = 6, size=3)+
    annotate("text", label = "Tephrite ", x = 43.5, y = 7, size=3)+
    annotate("text", label = "Phono- \n tephrite", x = 48.5, y = 9.5, size=3)+
    annotate("text", label = "Tephri- \n phonolite", x = 52.5, y = 11.5, size=3)+
    annotate("text", label = "Phonolite", x = 57, y = 14, size=3)+
    annotate("text", label = "Foidite", x = 42, y = 12, size=3)+
    geom_point(data = coord, aes(x=sio2_d,y=na2o_d+k2o_d, fill = location3) , shape = 21, size=7, stroke=0.8, alpha=0.8, color ="black")+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=0))+
    # geom_label_repel(mapping = aes(label = code), size = 2, color = "black", 
      #               fill= "white", force = 1, seed = 10, max.overlaps = 30)+ 
    geom_text_repel(aes(label = code), size = 3, color = "black", force = 5, max.overlaps = 30)+
    scale_fill_viridis(discrete=TRUE)
  tas
  
  

###################################################

###   Microbial analysis

# load otu table from DADA2
otu_table <- read.table("~/ASV_tab_merg.csv", sep = ",", header = TRUE)
row.names(otu_table) <- as.character(unlist(otu_table[,1]))
otu_table = otu_table[,-1 ]
colnames(otu_table)

# load taxa table
tax_table <- read.table("~/taxa_merg.csv", sep = ",", header = TRUE)
row.names(tax_table) <- as.character(unlist(row.names(otu_table)))
tax_table= tax_table[,-1]
tax_table = as.matrix(tax_table)

# load metadata table
database <- read.csv("~/data_sample.csv", sep = ",")
row.names(database)= database[,5]
row.names(database)

# creat phyloseq object
otu_table <- as.matrix(otu_table)
tax_table <- as.matrix(tax_table)
taxa_names(otu_table)
ps <- phyloseq(otu_table(otu_table, taxa_are_rows=TRUE), tax_table(tax_table), sample_data(database))
sum(sampleSums(ps))

# create the sequences object and name ASV
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
sample_sums(ps)
head(refseq(ps))

# remove unwanted sample
ps_filt <- subset_samples(ps, sample_names(ps) != "LB180410.f_S44262")
ps_filt <- subset_samples(ps_filt, sample_names(ps_filt) != "PGF.PG170224_S38589")
ps_filt <- subset_samples(ps_filt, sample_names(ps_filt) != "PGS.PG170224_S38581")
ps_filt  <- subset_samples(ps_filt, sample_names(ps_filt) != "Extraction.blank_S44240")
ps_filt  <- subset_samples(ps_filt, sample_names(ps_filt) != "S33.S.33_S38573")
ps_filt  <- subset_samples(ps_filt, sample_names(ps_filt) != "S14.S.14_S38590")
ps_filt <- subset_samples(ps_filt, sample_names(ps_filt) != "PC140111.f_S44238")
ps_filt  <- subset_samples(ps_filt, sample_names(ps_filt) != "DS170420.f_S44276")
ps_filt  <- subset_samples(ps_filt, sample_names(ps_filt) != "AG170420.f_S44247")
ps_filt


## Saving the phyloseq objects that can be shared and deposited as QC_data for downstearm applications
saveRDS(ps_filt, "prok_data_row.rds") # Basic phyloseq object for future applications

# removed of Eukaryota, Chloroplast, Mitochondria
ps_noECM <- subset_taxa(ps_filt,  (Kingdom != "Eukaryota") | is.na(Kingdom))
ps_noECM <- subset_taxa(ps_noECM, (Order!="Chloroplast") | is.na(Order))
ps_noECM <- subset_taxa(ps_noECM, (Family!="Mitochondria") | is.na(Family))
ps_noECM

# Removing the known DNA Extraction contaminants from Sheik et al., 2018
# Assuming you are starting from a phyloseq object called prok_data_raw, for the contaminant removal (cr) and ends with an object called prok_data_cr. This can be further processed to
# remove the pathogens, starting from prok_data_cr and ending with
# a prok_data obejct. Apply this before normalization.

prok_data_cr <- subset_taxa(ps_noECM,  (Genus != "Afipia") | is.na(Genus))
##prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Aquabacterium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Asticcacaulis") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Aurantimonas") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Beijerinckia") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Bosea") | is.na(Genus))
##prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Bradyrhizobium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Brevundimonas") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Caulobacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Craurococcus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Devosia") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Hoefleae") | is.na(Genus))
##prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Mesorhizobium") | is.na(Genus))
##prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Methylobacterium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Novosphingobium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Ochrobactrum") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Paracoccus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Pedomicrobium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Phyllobacterium") | is.na(Genus))
##prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Rhizobium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Sphingobium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Sphingomonas") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Sphingopyxis") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Acidovorax") | is.na(Genus))
##prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Azoarcus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Azospira") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Burkholderia") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Comamonas") | is.na(Genus))
##prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Cupriavidus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Curvibacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Delftiae") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Duganella") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Herbaspirillum") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Janthinobacterium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Kingella") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Leptothrix") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Limnobacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Massilia") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Methylophilus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Methyloversatilis") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Oxalobacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Pelomonas") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Polaromonas") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Neisseria") | is.na(Genus))
##prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Ralstonia") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Schlegelella") | is.na(Genus))
##prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Sulfuritalea") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Undibacterium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Variovorax") | is.na(Genus))
##prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Acinetobacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Enhydrobacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Enterobacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Escherichia") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Nevskia") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Pasteurella") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Pseudoxanthomonas") | is.na(Genus))
#prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Psychrobacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Stenotrophomonas") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Xanthomonas") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Aeromicrobium") | is.na(Genus))
#prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Actinomyces") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Arthrobacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Beutenbergia") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Brevibacterium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Corynebacterium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Curtobacterium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Dietzia") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Janibacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Kocuria") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Microbacterium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Micrococcus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Microlunatus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Patulibacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Propionibacterium") | is.na(Genus))
#prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Rhodococcus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Tsukamurella") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Chryseobacterium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Dyadobacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Flavobacterium") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Hydrotalea") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Niastella") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Parabacteroides") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Pedobacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Prevotella") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Wautersiella") | is.na(Genus))
#prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Deinococcus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Abiotrophia") | is.na(Genus))
#prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Bacillus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Brevibacillus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Brochothrix") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Facklamia") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Olivibacter") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Lactobacillus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Paenibacillus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Ruminococcus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Staphylococcus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Streptococcus") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Veillonella") | is.na(Genus))
prok_data_cr <- subset_taxa(prok_data_cr,  (Genus != "Fusobacterium") | is.na(Genus))
prok_data_cr

readcount(otu_table(prok_data_cr))/readcount(otu_table(ps_noECM))

# Removing the potential human pathogens and contaminants
# This step needs to be evaluated with attention since many of these genera
# might be relevant in many environmental settings. Usually it is better to compare
# before-after removal to see what and how much you are removing. Feel free to
# experiment with the different groups and evaluate the results.

prok_data <- subset_taxa(prok_data_cr, (Genus != "Abiotrophia") |  is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Achromobacter") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Acinetobacter") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Actinobacillus") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Arcanobacterium") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Babesia") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Bifidobacterium") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Bartonella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Bordetella") |  is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Borrelia") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Brodetella") |  is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Brucella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Capnocytophaga") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Chlamydia") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Citrobacter") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Comamonas") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Corynebacterium_1") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Corynebacterium") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Coxiella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Cronobacter") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Cutibacterium") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Dermatophilus") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Ehrlichia") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Order != "Enterobacteriales") | is.na(Order))
prok_data <- subset_taxa(prok_data, (Genus != "Enterococcus") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Erysipelothrix") |  is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Escherichia") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Escherichia/Shigella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Francisella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Gardnerella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Granulicatella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Haemophilus") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Hafnia") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Helicobacter") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Klebsiella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Kocuria") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Lactococcus") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Lactobacillus") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Lawsonia") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Legionella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Leptospira") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Listeria") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Merkel_cell") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Micrococcus") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Morganella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Mycoplasma") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Neisseria") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Nocardia") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Pasteurella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Plesiomonas") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Propionibacterium") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Proteus") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Providencia") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Pseudomonas") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Rhodococcus") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Rickettsiae") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Roseomonas") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Rothia") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Salmonella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Serratia") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Shewanella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Shigella") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Sphaerophorus") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Staphylococcus") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Stenotrophomonas") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Streptococcus") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Treponema") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Vibrio") | is.na(Genus))
prok_data <- subset_taxa(prok_data, (Genus != "Yersinia") | is.na(Genus))
prok_data

# remove low abundance ASV
prok_data_prune <- prune_taxa(taxa_sums(prok_data) > 5, prok_data)
prok_data_prune

# After these filtering step remove the ASV and samples that are left with all zeros
prok_data_prune = filter_taxa(prok_data_prune, function(x) sum(x) > 0, TRUE)
sum(sample_sums(prok_data_prune))

# Plot Richness on not normalized data
p_richness <- plot_richness(prok_data, x = "location3" , measures =c("Chao1")) +
   geom_boxplot(aes(fill=location3)) +
    scale_fill_viridis(discrete=TRUE)+
#xlab("Hot spring group") +
ylab("Chao1") +
   theme_classic()+
   facet_wrap(~sample_type) +
   theme_classic()
   
########  Suppl. Figure 4
p_richness$data$tectonics4 <- factor(p_r$data$location3, levels=order)
p_richness


# Normalize the counts across the different samples by converting the abundance to
# relative abundance and multiply by the median library size
ps_filt_norm <- transform_sample_counts(prok_data_prune, function(x) ((x / sum(x))*median(readcount(prok_data_prune))))
readcount(ps_filt_norm)

##################################################

##############   Figure 4
##############   Supll- Figure 6
# NMDS

# Create the weighted and unweighted distance matrices
ps_j_un <- distance(ps_filt_norm, method = "jaccard", binary = TRUE)
ps_j_w <- distance(ps_filt_norm, method = "jaccard")

## nMDS with Jaccard weighted and unweighted
prok_nmds_w <- ordinate(ps_filt_norm, ps_j_w, method = "NMDS", trymax=100)
prok_nmds_uw <- ordinate(ps_filt_norm, ps_j_un, method = "NMDS", trymax=100)

#prok_mds_w <- ordinate(ps_filt_norm, ps_j_w, method = "MDS", trymax=100)
#  prok_mds_uw <- ordinate(ps_filt_norm, ps_j_un, method = "MDS", trymax=100)

# set order for the visualization of the plot
order_2 <- c("Backarc", "Guanacaste_region","Cordillera central","Outer Forearc", "Panama")


NMDS<-plot_ordination(ps_filt_norm, prok_nmds_w, 
                      color = "Ca_Set", shape = "sample_type") 
#scale_color_manual(values = c("darkslateblue", "turquoise1", "dodgerblue2", "green2", "yellow"))
# get_legend
leg_NMDS<-get_legend(NMDS)
leg<-as.ggplot(leg_NMDS)

#############################################


plot_ordination(ps_filt_norm, prok_nmds_w, type = "sample") +
  geom_point(aes(fill=location3, shape=sample_type), size=5) +
  scale_shape_manual(values = c(21,24))+
  stat_ellipse(aes(color=location3), type="t", linetype=2 , show.legend=FALSE, level=0.95) +
  stat_stars(aes(color=location3), geom="segment", show.legend=FALSE, alpha=0.2)+
  #geom_encircle(aes(fill=location3), alpha=0.1, s_shape=1, expand=0, show.legend=FALSE) +
  #scale_fill_discrete(labels = c("Baia Terranova", "Ross Sea")) +
  #scale_color_discrete(labels = c("Baia Terranova", "Ross Sea")) +
  labs(fill = "Sampling area",
       color = "95% confidence interval")+
  # caption = 'Giovannelli Lab / @d_giovannelli') + 
  theme_glab(base_size = 20)   

# correlation NMDS dimension and varibles
scores(prok_nmds_w)
tab <- as.data.frame(scores(prok_nmds_w))
tab$Temp <- sample_data(ps_filt_norm)$Temp
tab$pH <- sample_data(ps_filt_norm)$pH
tab$DSub <- sample_data(ps_filt_norm)$DSub
tab$rc_ra <- sample_data(ps_filt_norm)$rc_ra

library(rsq)
# linear model
ggplot(tab,aes(x=rc_ra,y=NMDS2)) +
  geom_point( shape = 21, size=4, stroke=1, alpha=1) +
  geom_smooth(method=lm)


#######################################

#    ENVFIT

## Create two functions for vector fitting in Vegan
# convert the sample_data() within a phyloseq object to a vegan compatible data object
pssd2veg <- function(physeq) {
  sd <- sample_data(physeq)
  return(as(sd,"data.frame"))
}
# Create an environmental data object only
environmental <- pssd2veg(ps_filt_norm)
colnames(environmental)

env_geophisic <- environmental[,c(119:132,134:139)]
env_environment <- environmental[,c(38:41,43:45,47:51,56:62,8:10,32:36)]
env_geochemic_major<-environmental[,c(71:80)]
env_geochemic_trace<-environmental[,c(82:88,90:94,96:98,100:105,107:111,114,116:118)]

# use ordination NMDS Jacc
env_fitting <- envfit(prok_nmds_w, env_environment, perm=999, na.rm =TRUE)
env_fitting
p.adjust(env_fitting$vectors$pvals, method="bonferroni")

env_fitting_crust <- envfit(prok_nmds_w, env_geophisic, perm=999, na.rm =TRUE)
env_fitting_crust
p.adjust(env_fitting_crust$vectors$pvals, method="holm")

env_fitting_maj <- envfit(prok_nmds_w, env_geochemic_major, perm=999, na.rm =TRUE)
env_fitting_maj
p.adjust(env_fitting_maj$vectors$pvals, method="holm")

env_fitting_min <- envfit(prok_nmds_w, env_geochemic_trace, perm=999, na.rm =TRUE)
env_fitting_min
p.adjust(env_fitting_min$vectors$pvals, method="holm")

# plot nmds 
plot(prok_nmds_w, display = "sites", main = "nMDS weighted Jac") +
    geom_point(size=8, fill="black")
   
# plot envfit
plot(env_fitting, p.max = 0.05, col = "red")
plot(env_fitting_maj, p.max =0.05, col = "gold")
plot(env_fitting_min, p.max = 0.05, col = "darkblue")
plot(env_fitting_crust, p.max = 0.05, col = "darkgreen")
plot(env_fitting_tot, p.max = 0.05, col = "grey")
plot(env_fitting_pochi, p.max = 0.05, col = "purple")

########################################################

######################################################

######    Figure 3

#  hierarchical cluster plot
 
dist.bray<-vegdist(t(otu_table(ps_filt_norm)),method="jaccard")
hc_bray <- hclust(dist.bray, method = "complete") # heirarchal clustering
# the agglomeration method to be used. 
# This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#convert cluster object to use with ggplot
sample_norm <- sample_data(ps_filt_norm)
dendr    <- dendro_data(hc_bray, type="rectangle") # convert for ggplot

## !! ## 
#  add the grouping factor
clust.df <- data.frame(label = sample_norm$code, cluster =factor(sample_norm$location3))
# clust    <- cutree(hc_bray,k=3)                    # find 2 clusters
# dendr[["labels"]] has the labels, merge with clust.df based on label column
dendr[["labels"]] <- merge(dendr[["labels"]],clust.df, by="label")

# plot the dendrogram; note use of color=cluster in geom_text(...)
hclust_filt <- ggplot() + 
  geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend)) + 
  # geom_text(data=label(dendr), aes(x, y, label=label, hjust=0, color=cluster),   size=2.5) +
  geom_text(data=dendr$label, aes(x, y, label=label, hjust=-0.5, fontface = "bold"),  size=2) +
  geom_point(data=dendr$label, aes(x, y,  fill=cluster, shape = sample_norm$tectonics4),  size=2) +
  scale_fill_manual(values = c("darkslateblue", "turquoise1", "dodgerblue2", "green2", "yellow"))+
  coord_flip() + scale_y_reverse(expand=c(0.3, 0)) + 
  theme(axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.title.x = element_blank(),
        axis.text = element_text(size = (8)),
        panel.background=element_rect(fill="white"),
        plot.margin = unit(c(0, 0.5,0, 0), "cm"),
        panel.grid=element_blank())+
  theme(plot.margin = unit(c(0.5,-10,0.5,2), "pt"))+
  scale_shape_manual(values=c(21,22,23,24,25))+
  theme(legend.position = "none")
  
hclust_filt

#    Barplot
#    script for plotbar with phylum > 1% of abundance
#    be careful of the otu-table orientation

ps_norm_phylum_prune <- prune_taxa(taxa_sums(ps_phylum_rel) > sum(readcount(ps_phylum_rel))/100, ps_phylum_rel)
phy_norm <- data.frame(otu_table(ps_norm_phylum_prune))
Others_norm <- c(colSums(subset(otu_table(ps_phylum_rel), taxa_sums(ps_phylum_rel) < sum(readcount(ps_phylum_rel)/100))))
phy_norm<- rbind(phy_norm, Others_norm)
list <- c(data.frame(tax_table(ps_norm_phylum_prune))$Phylum, "Others")
list[is.na(list)] <- "NA"

row.names(phy_norm) <- list
# process Proteobacteria abundances, in order to visualise at Class level
phy_norm_prune_proteoB<- phy_norm[!(row.names(phy_norm) %in% c("Proteobacteria")), ]
ps_proteo <- subset_taxa(ps_class_rel, Phylum == "Proteobacteria")
ps_proteo_class = tax_glom(ps_proteo, taxrank="Class", NArm=FALSE)
phy_proteo <- data.frame(otu_table(ps_proteo_class))
colSums(phy_proteo)
colSums(phy_norm_prune_proteoB)
Proteo_Class <- c(data.frame(tax_table(ps_proteo_class))$Class)
# if is present an "NA" group in Proteobacteria
Proteo_Class[is.na(Proteo_Class)] = "NA Proteobacteria"
row.names(phy_proteo) <- Proteo_Class
phy_tab <- rbind(phy_proteo, phy_norm_prune_proteoB)

# Gathering data  !!!
# for the ggplot barplot

phy_tab$Taxa <- row.names(phy_tab)
aggdata_phy <- gather(phy_tab, key = "sample", value = "Abundance", -Taxa)
# Converte to factor for preserving sequence in our visualisation
# same for the sample
aggdata_phy$Taxa<- factor(aggdata_phy$Taxa,  levels = unique(aggdata_phy$Taxa))
aggdata_phy$sample<- factor(aggdata_phy$sample,  levels = unique(aggdata_phy$sample))
aggdata_phy

# use the hclust order for the barplot 
dendr_for_barplot <- as.dendrogram(hc_bray)
labels(dendr_for_barplot)
aggdata_phy$sample <- factor(aggdata_phy$sample, levels = labels(dendr_for_barplot))

mainplot <- ggplot(aggdata_phy, aes(width = 0.8, fill=Taxa, y=Abundance, x=sample)) +
	geom_bar(position="stack", stat="identity") +
	#scale_x_discrete(guide = guide_axis(angle = 0)) +
	#scale_x_discrete(limits = rev)+
	labs(y = "Relative abundance (%)")+
	scale_fill_manual(values=mycols4)+ coord_flip() +  # per girarlo orizzontalmente
	theme(plot.margin = unit(c(18,2,5,0), "pt"))+
	theme(axis.title.y = element_blank(), 
	panel.background=element_rect(fill="white"),
	plot.margin = unit(c(0.3, 0.5,0, 0), "cm"),
	axis.text = element_text(size = (8)),
	panel.grid=element_blank(),
	axis.text.y=element_blank(),
	axis.ticks.y=element_blank())+
	theme(legend.position = "none")

# visualization of several variable using the same order of the cluster
Temp<- sample_data(ps_filt_norm)$Temp
df_temp <- data.frame(code, Temp )  
df_temp$Temp_code <- "Temp"
df_temp$code <- factor(df_temp$code, levels = labels(dendr_for_barplot))

col_temp<-ggplot(df_temp, aes(Temp_code, code, fill= Temp)) + 
	scale_fill_gradient2(low = "skyblue3", mid = "white", midpoint = 47, high = "orangered3" ) +
	geom_tile() +
 	theme(legend.position = "none",axis.text.y = element_blank(), axis.text.x = element_blank(),
        axis.ticks.y =  element_blank(),axis.title = element_blank(),
        axis.ticks.x =  element_blank(),
        panel.grid.major = element_blank(),  panel.grid.minor = element_blank(),
        panel.background=element_rect(fill="white"),
        plot.margin = unit(c(0.3, 0,0.8, 0), "cm"),
        panel.border = element_rect(fill = NA, colour = rgb(100, 100, 100, maxColorValue = 255)))+ 
	theme(legend.text = element_text(size=3),
        legend.key = element_blank())
col_temp
col_temp + scale_fill_viridis(option = "C")

#  add legend

leg_clust<-as.ggplot(get_legend(ggplot(coord, aes(Temp, location3))+
	geom_point(aes(color=location3), size = 5)+
	theme(title = element_blank())+
	scale_color_manual(values = c("darkslateblue", "turquoise1", "dodgerblue2", "green2", "yellow"))+
	theme(legend.key = element_blank())))
leg_clust

leg_bar<-as.ggplot(get_legend(ggplot(aggdata_phy, aes(sample, Abundance))+
	geom_point(aes(color=Taxa), size = 5)+
	theme(legend.key = element_blank())+
	scale_color_manual(values=mycols4)+
	guides(color = guide_legend(title = "Phylum",label.theme = element_text(size = 7,angle = 0), ncol = 4)) + 
	theme(legend.key.height = unit(1, "mm"))))
leg_bar

ggarrange(hclust_filt , mainplot, 
	ggarrange(col_temp, widths = c(0.1,0.1,0.1,0.1,0.1), ncol = 5, nrow = 1),
	leg_clust, leg_bar, 
	heights = c(2, 0.6), widths = c(1,2,1), ncol = 3, nrow = 2)

