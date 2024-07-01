library(tidyverse)
library(phyloseq)
source(file = "Scripts/0_Util_fonctions.R")

# Colors

color_code = c(
"Carnivora"="#f87575",
"Rodentia"="#d68c45",
"Cetartiodactyla"="#a7c957",
"Primates"="#13315c",
  
"Ape"="#13315c",
"Lemur"="#b8b8ff",
"NWM"="#e574bc",
"OWM"="#81667a",

"Europe"="#00b450",
"Australia"="#8da9c4",
"North_America"="#00b489",
"Asia"="#a8dadc",
"Africa"="#13315c",
"South_America"="#00b4d8",

"Industrial"="#8da9c4",
"Non-industrial"="#a7c957",
"Paleofeces"= "darkred"
)

# Load Amato et al. 
rawOTU_table = read.table("Data/Primates/qs850_qs10052_qs1448_qs11212_qs10317_qs11358_nobloom_nodoubletons_min1k.biom_tax_gut_adult_reduced_even10010_reduced_nochlor_nomit_even9870_fixed.txt",
                          check.names=F, header=T)

PT = rawOTU_table %>% 
  select(taxonomy,OTU_ID) %>% 
  separate(taxonomy, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep=";") %>% 
  as.matrix()

ASV_S <- rawOTU_table %>% 
  select(-OTU_ID,-taxonomy)
rownames(ASV_S) <- rawOTU_table$OTU_ID

# rectify sample names for american gut project 
AGP = read_tsv("Data/Primates/filereport_read_run_PRJEB11419_tsv.txt") %>%
  select(secondary_sample_accession,sample_title) %>% 
  distinct()

MT = read.table("Data/Primates/human_NHP_mapping.txt", header = T) %>% 
  left_join(AGP,by=c("Accession_number"="secondary_sample_accession")) %>% 
  mutate(SampleID_correction = ifelse(is.na(sample_title),SampleID,sample_title))

# Make phyloseq object
rownames(MT)=MT$SampleID_correction
rownames(PT)=PT[,'OTU_ID']

primate_PS=phyloseq(tax_table(PT),
                    otu_table(ASV_S,taxa_are_rows=T),
                    sample_data(MT)) %>%
  rarefy_even_depth(2000)

primate_PS


#################
# Human panel 16S #
#################

human_PS <- primate_PS %>% 
  subset_samples(Species=="Homo_sapiens")
human_PS

betaM="jaccard" # , binary = TRUE

# Ordination 
beta_values = human_PS %>% 
  distance(betaM, binary = TRUE)
ordination = beta_values %>% vegan::metaMDS() 

Ordination_scores_humans <-  ordination$points %>% 
  as_tibble(rownames = "SampleID_correction") %>% 
  left_join(MT) 

Panel_humans = Ordination_scores_humans  %>% 
  ggplot(aes(y=MDS1,x=MDS2,col=geographic_region,shape=life_style_2)) + #ylim(-.3,.3)+
  geom_point(alpha=.7,size=1)  +
  MyTheme +
  scale_color_manual(values=color_code) 

Panel_humans 

unique(Ordination_scores_humans$geographic_region)

#################
# Human panel metaG + feces #
#################

coprolites <- readxl::read_xlsx("Data/Human_Wibowo/41586_2021_3532_MOESM6_ESM.xlsx", 
                                sheet = 3, col_names = F)

#metadata
copro_MT <- as.data.frame(t(coprolites[1:2,-1]))
colnames(copro_MT) <- c("Type","SampleID")
row.names(copro_MT) <- copro_MT$SampleID

#reformat sp table
copro_sp <- coprolites[3:1201,]; sp_list <- copro_sp[,1]
copro_sp <- matrix(as.numeric(as.matrix(copro_sp[,-1])), ncol = 799) 
rownames(copro_sp) <- sp_list$...1; colnames(copro_sp) <- copro_MT$SampleID

#phyloseq
wibowo_PS=phyloseq(otu_table(copro_sp, taxa_are_rows=T),
                    sample_data(copro_MT)) %>% 
  subset_samples(!Type=="Soil")

wibowo_PS
sample_sums(wibowo_PS)

# plot
betaM="jaccard" # , binary = TRUE

# Ordination 
beta_values_wibo = wibowo_PS %>% 
  distance(betaM, binary = TRUE)
ordination_wibo = beta_values_wibo %>% vegan::metaMDS() 

Ordination_scores_humans_wibo <-  ordination_wibo$points %>% 
  as_tibble(rownames = "SampleID") %>% 
  left_join(copro_MT) %>% 
  mutate(sizeP = ifelse(Type=="Paleofeces",1.1,1))

Panel_humans_wibo = Ordination_scores_humans_wibo  %>% 
  mutate(Type =  factor(Type, levels = c("Industrial","Non-industrial","Paleofeces"))) %>%
  arrange(Type)   %>% 
  ggplot(aes(x=MDS1,y=MDS2,col=Type)) + #ylim(-.3,.3)+
  geom_point(aes(size=sizeP,alpha=sizeP))  +
  MyTheme +
  scale_color_manual(values=color_code)  + scale_size(range = c(1, 2)) +
  scale_alpha(range = c(0.7, 1))
  
Panel_humans_wibo


#################
# Primate panel #
#################


betaM="jaccard" # , binary = TRUE

# Ordination 
beta_values = primate_PS %>% 
  distance(betaM, binary = TRUE)
ordination = beta_values %>% vegan::metaMDS() 

Ordination_scores_primates <-  ordination$points %>% 
  as_tibble(rownames = "SampleID_correction") %>% 
  left_join(MT) 

Panel_primates = Ordination_scores_primates  %>% 
  ggplot(aes(y=MDS1,x=MDS2,col=phyl_group2)) + #ylim(-.3,.3)+
  geom_point(alpha=1,size=1)  +
  MyTheme  +
  scale_color_manual(values=color_code) 

Panel_primates



#################
# Mammals panel #
#################

# Load Song et al. 
ASV_S_mammals = read.table("Data/Mammals/Song_ASVs_counts_filt.txt")
PT_mammals = as.matrix(read.table("Data/Mammals/Song_ASVs_taxonomy_filt.txt"))
MT_mammals = read.table("Data/Mammals/Song_Metadata.txt")

# Make phyloseq object
rownames(MT_mammals)=MT_mammals$SampleID
rownames(PT_mammals)=PT[,'ASV']
mammals_PS=phyloseq(tax_table(PT_mammals),otu_table(ASV_S_mammals,taxa_are_rows=T),sample_data(MT_mammals))
mammals_PS

# Filter primate samples 
set.seed(30052024)

main_orders <- mammals_PS %>% 
  subset_samples(Order=c("Carnivora","Rodentia","Cetartiodactyla","Priamtes")) %>%
  subset_samples(SampleID %in% subset_sp) %>%
  rarefy_even_depth(2000) 

main_orders <- mammals_PS %>% 
  subset_samples(Order%in%c("Carnivora","Rodentia","Cetartiodactyla","Primates")) %>%
  rarefy_even_depth(2000) 
main_orders

#only wild specimens
subset_sp <- sample_data(main_orders) %>%
  subset(captive_wild=="wild") %>%
  pull(SampleID)

main_orders <- main_orders %>% 
  subset_samples(SampleID %in% subset_sp)
main_orders

betaM="jaccard" # , binary = TRUE

beta_values = main_orders  %>% 
  distance(betaM)#, binary = TRUE)
ordination = beta_values %>% vegan::metaMDS() 

Ordination_scores_mammals <-  ordination$points %>% 
  as_tibble(rownames = "SampleID") %>% 
  left_join(MT) %>% 
  mutate(Species=Time_Tree_curated)

Panel_mammals = Ordination_scores_mammals %>% 
  ggplot(aes(y=MDS1,x=MDS2,col=Order)) + #ylim(-.3,.3)+
  geom_point(alpha=1,size=1)  +
  MyTheme  +
  scale_color_manual(values=color_code) 

Panel_mammals


#Save panels 
ggsave("Figs/Prep/Panel_humans.pdf",plot = Panel_humans,
       width = 100, height = 70, device = 'pdf',units="mm")

ggsave("Figs/Prep/Panel_humans_wibowo.pdf",plot = Panel_humans_wibo,
       width = 100, height = 70, device = 'pdf',units="mm")

ggsave("Figs/Prep/Panel_mammals.pdf",plot = Panel_mammals,
       width = 100, height = 70, device = 'pdf',units="mm")

ggsave("Figs/Prep/Panel_primates.pdf",plot = Panel_primates,
       width = 100, height = 70, device = 'pdf',units="mm")

Panel_humans
Panel_mammals
Panel_primates






