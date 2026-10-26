
######## Experiment I field surveys ----------------------------------------------------------------------------------------------------------------
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf") 
######## normality and homogeneity of variance  --------------------------------
library(ggpubr) 
library(car)
data = read.csv("field_data_155_35sp.csv",header=TRUE, row.names=1)

ggqqplot(data$P.concentrations, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$Root.diameter, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$Root.length.colonization, color = "blue", main="Normal Q-Q Plot")


data = read.csv("field_data_391otu_raw_155.csv",header=TRUE, row.names=1)

ggqqplot(data$P.concentrations, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$Root.diameter, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$Root.length.colonization, color = "blue", main="Normal Q-Q Plot")

leveneTest( data$P.concentrations  ~ Plant.species , data = data) 
leveneTest( data$Root.diameter  ~ Plant.species , data = data) 
leveneTest( data$Root.length.colonization  ~ Plant.species , data = data) 



########  v.phylomaker tree  35sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf") 
df <- read.csv("35sp.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "1_35sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL
#write.tree(sp_tr, "1-35sp-phylomakers1-mega.nwk")

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  

####
library(dplyr)
tree <- read.tree("1_35sp_phylomakers1.nwk")
name_map <- c(
  "Elymus_nutans" = "En",
  "Poa_pratensis" = "Pp",
  "Poa_poophagorum" = "Po",
  "Koeleria_macrantha" = "Km",
  "Carex_capillifolia" = "Cp",
  "Trichophorum_pumilum" = "Tp",
  "Carex_kokanica" = "Ck",
  "Aster_souliei" = "Au",
  "Leontopodium_leontopodioides" = "Ll",
  "Swertia_tetraptera" = "St",
  "Ligularia_virgaurea" = "Lv",
  "Anemone_obtusiloba" = "Ao",
  "Trollius_farreri" = "Tf",
  "Ranunculus_tanguticus" = "Rt",
  "Halenia_elliptica" = "He",
  "Elsholtzia_densa" = "Ed",
  "Deschampsia_cespitosa" = "Dc",
  "Delphinium_kamaonense" = "Dk",
  "Potentilla_fragarioides" = "Pf",
  "Argentina_anserina" = "Aa",
  "Saussurea_nigrescens" = "Sn",
  "Anaphalis_lactea" = "Al",
  "Festuca_ovina" = "Fo",
  "Thalictrum_alpinum" = "Tl",
  "Artemisia_gmelinii" = "Ag",
  "Salvia_roborowskii" = "Sr",
  "Allium_sikkimense" = "As",
  "Anemone_trullifolia" = "Ar",
  "Veronica_vandellioides" = "Vv",
  "Lagotis_brachystachya" = "Lb",
  "Cerastium_arvense" = "Ca",
  "Ajania_tenuifolia" = "At",
  "Gymnaconitum_gymnandrum" = "Gg",
  "Artemisia_princeps" = "Ap",
  "Agrostis_hugoniana" = "Ah"
)


setdiff(tree$tip.label, names(name_map))  

tree$tip.label <- name_map[tree$tip.label]
any(is.na(tree$tip.label))

write.tree( tree, file = "1_35sp_phylomakers1_abbr.nwk")




########  v.phylomaker tree  31sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf") 
df <- read.csv("31sp.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "1_31sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL
#write.tree(sp_tr, "1-35sp-phylomakers1-mega.nwk")

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  

####
library(dplyr)
tree <- read.tree("1_31sp_phylomakers1.nwk")
name_map <- c(
  "Elymus_nutans" = "En",
  "Poa_pratensis" = "Pp",
  "Poa_poophagorum" = "Po",
  "Koeleria_macrantha" = "Km",
  "Carex_capillifolia" = "Cp",
  "Trichophorum_pumilum" = "Tp",
  "Carex_kokanica" = "Ck",
  "Aster_souliei" = "Au",
  "Leontopodium_leontopodioides" = "Ll",
  "Swertia_tetraptera" = "St",
  "Ligularia_virgaurea" = "Lv",
  "Anemone_obtusiloba" = "Ao",
  "Trollius_farreri" = "Tf",
  "Ranunculus_tanguticus" = "Rt",
  "Halenia_elliptica" = "He",
  "Elsholtzia_densa" = "Ed",
  "Deschampsia_cespitosa" = "Dc",
  "Delphinium_kamaonense" = "Dk",
  "Potentilla_fragarioides" = "Pf",
  "Argentina_anserina" = "Aa",
  "Saussurea_nigrescens" = "Sn",
  "Anaphalis_lactea" = "Al",
  "Festuca_ovina" = "Fo",
  "Thalictrum_alpinum" = "Tl",
  "Artemisia_gmelinii" = "Ag",
  "Salvia_roborowskii" = "Sr",
  "Allium_sikkimense" = "As",
  "Anemone_trullifolia" = "Ar",
  "Veronica_vandellioides" = "Vv",
  "Lagotis_brachystachya" = "Lb",
  "Cerastium_arvense" = "Ca",
  "Ajania_tenuifolia" = "At",
  "Gymnaconitum_gymnandrum" = "Gg",
  "Artemisia_princeps" = "Ap",
  "Agrostis_hugoniana" = "Ah"
)


setdiff(tree$tip.label, names(name_map))  

tree$tip.label <- name_map[tree$tip.label]
any(is.na(tree$tip.label))

write.tree( tree, file = "1_31sp_phylomakers1_abbr.nwk")




######## -- Fig S1 -- 155 35sp mean calculation & Plant phylogenetic signals of plant traits, RLC ----------------
library(plyr) 
library(ggtree)
library(ape) 
library(picante) 
library(phylosignal) 
library(phylobase) 
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf") 
data = read.csv("field_data_391otu_raw_155.csv",header=TRUE, row.names=1)

data_mean <- ddply(data, c("Plant.species","Plant.species0",  "Plant.family"),summarise,   
                   P.concentrations = mean(P.concentrations,na.rm=TRUE),
                   Root.diameter = mean(Root.diameter,na.rm=TRUE),
                   Root.length.colonization = mean(Root.length.colonization,na.rm=TRUE)
)

data_mean

data.table::fwrite(data_mean,"field_data_155_35sp.csv",  row.names = TRUE, bom = TRUE)

data_mean1 = data_mean
rownames(data_mean1) <- data_mean[, 1]
data_mean1 <- data_mean1[ , -c(1,2,3)]

##
rooted_tree <- read.tree("1_35sp_phylomakers1_abbr.nwk") 
plot(rooted_tree)

tree_p <- as(rooted_tree, "phylo4")
class(tree_p)
p4d <- phylo4d(tree_p, data_mean1)
phy.sig <- phyloSignal(p4d, reps = 9999)
phy.sig ## 35sp 155
#write.csv(phy.sig, '1_field_data_physig.csv')  

##
library(ggplot2) 
library(ggpubr)
library(ggtreeExtra)
library(ggstar) 
library(ggtree) 
library(treeio)
library(ggnewscale) 
library(aplot)
library (reshape2)
library(vegan)
library(ape)
library(plyr)
library( RColorBrewer)
library(gghalves)
library(doBy)

se <- function(x) sd(x)/sqrt(length(x))
conf_95 <- function(x) t.test(x, conf.level = 0.95)$conf.int[2]

sumdata <- summaryBy( P.concentrations + Root.diameter + Root.length.colonization ~ Plant.species, data, FUN = c(mean, sd, se, conf_95))
sumdata
sumdata$Family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
sumdata$P.concentrations <- sumdata$P.concentrations.mean
sumdata$Root.diameter <- sumdata$Root.diameter.mean
sumdata$Root.length.colonization <- sumdata$Root.length.colonization.mean
colnames(sumdata)

group_file <- sumdata[,c(1, 14)]
row.names(group_file) <- group_file$Plant.species
groupInfo_list <- list()
groupInfo_list[["Family"]] <- split(row.names(group_file), group_file$Family)
for (group_name in names(groupInfo_list)) {
  rooted_tree <- groupOTU(rooted_tree, groupInfo_list[[group_name]], group_name)
}
rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)
rooted_tree.data = fortify(rooted_tree)
head(rooted_tree.data)

cols <- c("Amaryllidaceae"= "#5F97C6", "Cyperaceae" = "#A3B150", "Poaceae" = "#8FB9D1", 
          "Ranunculaceae" = "#B85A68", 
          # 
          "Rosaceae" =  "#FFA64D", "Fabaceae" = "#5880ba" , 
          "Brassicaceae"= "#9e5c93", "Onagraceae" = "#7B8F5E", 
          # 
          "Caryophyllaceae" = "#A83270" ,  "Polygonaceae"= "#9BBF6E",  
          "Gentianaceae"= "#E0A16F", "Lamiaceae"  = "#B07AA8", "Plantaginaceae" = "#497d75", 
          "Apiaceae" = "#BFA87F", "Asteraceae" = "#E9746A")

p00 <- ggtree( rooted_tree , layout="rectangular", ladderize = T,  branch.length = "none", size=0.7 ) +  
  geom_tippoint(aes(color= Family ),  size = 2 , alpha = 1, position = position_nudge(x = 0.01) ) + 
  geom_tiplab(aes(color= Family), size=4, fontface = 3, hjust = 0, alpha = 1, offset=0.2, show.legend=F ) + 
  scale_color_manual(values =cols) +  
  theme(legend.position="none", legend.title=element_text(size=11), legend.text =element_text(size=9),
        axis.title = element_text(size=10)) + labs(color = "Plant family") + geom_rootedge(rootedge = 0) +
  guides(color = guide_legend(override.aes = list(size = 1.5, shape = 16), order = 1, nrow = 12))   
p00


p1 <- ggplot(sumdata, aes(x = Plant.species, y= P.concentrations, fill = Family)) +  
  geom_col(width = 0.6, color = "black", alpha = 1) +
  geom_errorbar(aes(ymin= P.concentrations - P.concentrations.se, ymax= P.concentrations + P.concentrations.se), width=0 ,linewidth = 0.6)+ 
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.text.y= element_blank(),
        axis.text=element_text(size=11), legend.position = "none",  axis.title=element_text(size=12))+
  scale_fill_manual(values =cols) +
  labs(x="",y="Shoot P concentration (mg/g)")+
  coord_flip()
p1


p2 <- ggplot(sumdata, aes(x = Plant.species, y= Root.diameter , fill = Family)) +  
  geom_col(width = 0.6, color = "black", alpha = 1) +
  geom_errorbar(aes(ymin= Root.diameter - Root.diameter.se , ymax= Root.diameter + Root.diameter.se ), width=0 ,linewidth = 0.6 )+ 
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), axis.text.y= element_blank(),
        axis.text=element_text(size=11),legend.position = "none",  axis.title=element_text(size=12))+ 
  scale_fill_manual(values = cols) +
  labs(x="",y="Root diameter (mm)")+ 
  coord_flip()
p2


p3 <- ggplot(sumdata, aes(x = Plant.species, y= Root.length.colonization, fill = Family)) +  
  geom_col(width = 0.6, color = "black", alpha = 1) +
  geom_errorbar(aes(ymin=Root.length.colonization - Root.length.colonization.se, ymax= Root.length.colonization + Root.length.colonization.se), width=0 ,linewidth = 0.6)+ 
  theme_classic()+
  scale_y_continuous(expand = c(0, 0))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(), axis.text.y= element_blank(),
        axis.text=element_text(size=11),legend.position = "none",  axis.title=element_text(size=12))+
  scale_fill_manual(values = cols) +
  labs(x="",y="Root length colonization (%)")+ 
  coord_flip()
p3

#2_fig.s1_10*9
figs1 <- p1 %>% insert_left(p00, width = 0.8) %>% insert_right(p2,width = 1) %>% 
  insert_right(p3,width = 1)

figs1

#eoffice::topptx( figs1 , width = 10 ,  height = 9 , "1_155_35sp_figs1_10_9.pptx" ) 


######## -- Table S3  -- 31sp Sensitivity analysis between the all-species dataset and the typical-AM-species dataset----------------
library(car)
library(ggplot2)
library(ggpmisc)
library(scales)
library(ggpubr)
library(cowplot)
library(ppcor)
library(doBy)
library(caper)
rm(list = ls()) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf") 
data = read.csv("field_data_155_31sp.csv",header=TRUE, row.names=1) 

##
rooted_tree <- read.tree("1_31sp_phylomakers1.nwk")
plot(rooted_tree)

data$Plant.species0 <- gsub(" ", "_", data$Plant.species0)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data , names.col = "Plant.species0", vcv = TRUE, vcv.dim = 2)
nrow(comp_data$data)
colnames(data)


###
lm_model <- lm( Root.length.colonization ~ Root.diameter , data = data) 
summary(lm_model)
pgls_model <- pgls(  Root.length.colonization ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model3 <- predict(pgls_model, newdata = data)  


##
lm_model <- lm( P.concentrations ~ Root.length.colonization, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model2 <- predict(pgls_model, newdata = data )  



########################################################################################################################################
######## Experiment II greenhouse experiments --------------------------------------------------------------------------------------------------------
######## normality and homogeneity of variance  --------------------------------
library(ggpubr) 
library(car)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf")
data = read.csv("greenhouse_data_mean.csv",header=TRUE, row.names=1)
colnames(data)

ggqqplot(data$P.concentrations, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$Total.biomass, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$Root.diameter, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$Root.length.colonization, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$MPR, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$MGR, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$MRR, color = "blue", main="Normal Q-Q Plot")


data = read.csv("greenhouse_data_2.csv",header=TRUE, row.names=1)
colnames(data)

ggqqplot(data$P.concentrations, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$Root.diameter, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$Root.length.colonization, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$MPR, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$MGR, color = "blue", main="Normal Q-Q Plot")
ggqqplot(data$MRR, color = "blue", main="Normal Q-Q Plot")


leveneTest( data$P.concentrations  ~ Plant.species , data = data) 
leveneTest( data$Root.diameter  ~ Plant.species , data = data) 
leveneTest( data$Root.length.colonization  ~ Plant.species , data = data) 

leveneTest( data$MPR  ~ Plant.species , data = data) 
leveneTest( data$MGR  ~ Plant.species , data = data) 
leveneTest( data$MRR  ~ Plant.species , data = data) 

##  v.phylomaker tree    45sp --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf") 
df <- read.csv("45sp.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "1_45sp_phylomakers1.nwk")


is.rooted(sp_tr)       
is.binary.tree(sp_tr) 


####
library(ape)
library(dplyr)
tree <- read.tree("1-45sp-phylomakers1.nwk")
name_map <- c(
  "Agrostis_hugoniana" = "Ah",
  "Anaphalis_lactea" = "Al",
  "Angelica_nitida" = "An",
  "Aster_souliei" = "Au",
  "Astragalus_polycladus" = "Ay",
  "Cerastium_arvense" = "Ca",
  "Capsella_bursa-pastoris" = "Cb",
  "Carum_carvi" = "Cc",
  "Calamagrostis_flavens" = "Cf",
  "Carduus_nutans" = "Cn",
  "Carex_capillifolia" = "Cp",
  "Chamaesium_thalictrifolium" = "Ct",
  "Deschampsia_cespitosa" = "Dc",
  "Delphinium_kamaonense" = "Dk",
  "Erigeron_breviscapus" = "Eb",
  "Elymus_nutans" = "En",
  "Epilobium_palustre" = "Ep",
  "Festuca_ovina" = "Fo",
  "Festuca_sinensis" = "Fs",
  "Gymnaconitum_gymnandrum" = "Gg",
  "Koeleria_macrantha" = "Km",
  "Lepidium_apetalum" = "La",
  "Leontopodium_leontopodioides" = "Ll",
  "Ligularia_virgaurea" = "Lv",
  "Oxytropis_kansuensis" = "Ok",
  "Plantago_depressa" = "Pd",
  "Potentilla_fragarioides" = "Pf",
  "Poa_poophagorum" = "Po",
  "Poa_pratensis" = "Pp",
  "Rumex_patientia" = "Rp",
  "Ranunculus_tanguticus" = "Rt",
  "Saussurea_alpina" = "Sa",
  "Silene_gonosperma" = "Sg",
  "Stipa_aliena" = "Sl",
  "Stellaria_media" = "Sm",
  "Saussurea_nigrescens" = "Sn",
  "Silene_aprica" = "Sp",
  "Saussurea_superba" = "Ss",
  "Stipa_purpurea" = "Su",
  "Thlaspi_arvense" = "Ta",
  "Trollius_farreri" = "Tf",
  "Tibetia_himalaica" = "Th",
  "Thalictrum_tenue" = "Tt",
  "Taraxacum_lugubre" = "Tu",
  "Veronica_eriogyne" = "Ve"
)


setdiff(tree$tip.label, names(name_map))  

tree$tip.label <- name_map[tree$tip.label]
any(is.na(tree$tip.label))

write.tree( tree, file = "1_45sp_phylomakers1_abbr.nwk")

##  v.phylomaker tree    37sp --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf") 
df <- read.csv("37sp.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "1_37sp_phylomakers1.nwk")

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  

####
library(ape)
library(dplyr)
tree <- read.tree("1_37sp_phylomakers1.nwk")
name_map <- c(
  "Agrostis_hugoniana" = "Ah",
  "Anaphalis_lactea" = "Al",
  "Angelica_nitida" = "An",
  "Aster_souliei" = "Au",
  "Astragalus_polycladus" = "Ay",
  "Cerastium_arvense" = "Ca",
  "Capsella_bursa-pastoris" = "Cb",
  "Carum_carvi" = "Cc",
  "Calamagrostis_flavens" = "Cf",
  "Carduus_nutans" = "Cn",
  "Carex_capillifolia" = "Cp",
  "Chamaesium_thalictrifolium" = "Ct",
  "Deschampsia_cespitosa" = "Dc",
  "Delphinium_kamaonense" = "Dk",
  "Erigeron_breviscapus" = "Eb",
  "Elymus_nutans" = "En",
  "Epilobium_palustre" = "Ep",
  "Festuca_ovina" = "Fo",
  "Festuca_sinensis" = "Fs",
  "Gymnaconitum_gymnandrum" = "Gg",
  "Koeleria_macrantha" = "Km",
  "Lepidium_apetalum" = "La",
  "Leontopodium_leontopodioides" = "Ll",
  "Ligularia_virgaurea" = "Lv",
  "Oxytropis_kansuensis" = "Ok",
  "Plantago_depressa" = "Pd",
  "Potentilla_fragarioides" = "Pf",
  "Poa_poophagorum" = "Po",
  "Poa_pratensis" = "Pp",
  "Rumex_patientia" = "Rp",
  "Ranunculus_tanguticus" = "Rt",
  "Saussurea_alpina" = "Sa",
  "Silene_gonosperma" = "Sg",
  "Stipa_aliena" = "Sl",
  "Stellaria_media" = "Sm",
  "Saussurea_nigrescens" = "Sn",
  "Silene_aprica" = "Sp",
  "Saussurea_superba" = "Ss",
  "Stipa_purpurea" = "Su",
  "Thlaspi_arvense" = "Ta",
  "Trollius_farreri" = "Tf",
  "Tibetia_himalaica" = "Th",
  "Thalictrum_tenue" = "Tt",
  "Taraxacum_lugubre" = "Tu",
  "Veronica_eriogyne" = "Ve"
)

setdiff(tree$tip.label, names(name_map))  

tree$tip.label <- name_map[tree$tip.label]
any(is.na(tree$tip.label))

write.tree( tree, file = "1_37sp_phylomakers1_abbr.nwk")

######## -- Table S2  -- 225 45sp  mean calculation & Plant phylogenetic signals of plant traits and AM fungal variables ----------------
library(plyr)
library(ggtree)
library(ape)
library(picante)
library(phylosignal)
library(phylobase)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf") 
data = read.csv("greenhouse_data_2.csv",header=TRUE, row.names=1)
colnames(data)


data_mean <- ddply(data, c("Plant.species0"),summarise,   
                   Root.diameter = mean(Root.diameter,na.rm=TRUE),
                   Root.diameter_NM = mean(Root.diameter_NM,na.rm=TRUE),
                   P.concentrations = mean(P.concentrations,na.rm=TRUE),
                   P.concentrations_NM = mean(P.concentrations_NM,na.rm=TRUE),
                   Total.biomass = mean(Total.biomass,na.rm=TRUE),
                   Total.biomass_NM= mean(Total.biomass_NM,na.rm=TRUE),
                   
                   Root.length.colonization = mean(Root.length.colonization,na.rm=TRUE),
                   MGR = mean(MGR, na.rm=TRUE), 
                   MPR = mean(MPR,na.rm=TRUE),
                   MRR = mean(MRR,na.rm=TRUE) )

data_mean

#write.csv(data_mean,"greenhouse_data_mean.csv")

data_mean1=data_mean
rownames(data_mean1) <- data_mean[, 1]
data_mean1 <- data_mean1[, -1]

##
rooted_tree <- read.tree("1_45sp_phylomakers1.nwk") 
plot(rooted_tree)

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

tree_p <- as(rooted_tree, "phylo4")
class(tree_p)
p4d <- phylo4d(tree_p, data_mean1)
phy.sig <- phyloSignal(p4d, reps = 9999)
#v.phylomaker tree 
phy.sig

write.csv(phy.sig, '1_greenhouse_data_v1physig.csv')  

######## -- Table S4  -- 37sp Regression relationships of MPR and MGR with root length colonization, root diameter, and MRR ----------------
library(car) 
library(ggplot2)
library(ggpmisc)
library(scales)
library(ggpubr)
library(cowplot)
library(ppcor)
library(doBy)
library(caper)
rm(list = ls()) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf")
data = read.csv("greenhouse_data_mean_37sp.csv",header=TRUE, row.names=1) 

##
rooted_tree <- read.tree("1_37sp_phylomakers1.nwk")
plot(rooted_tree)

data$Plant.species0 <- gsub(" ", "_", data$Plant.species0)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data , names.col = "Plant.species0", vcv = TRUE, vcv.dim = 2)
nrow(comp_data$data)
colnames(data)

###
lm_model <- lm( Root.length.colonization ~ Root.diameter , data = data) 
summary(lm_model)
pgls_model <- pgls(  Root.length.colonization ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model3 <- predict(pgls_model, newdata = data)  


###
lm_model <- lm( P.concentrations ~ Root.length.colonization , data = data) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model1 <- predict(pgls_model, newdata = data)  

###
lm_model <- lm( Total.biomass ~ Root.length.colonization , data = data) 
summary(lm_model)
pgls_model <- pgls(  Total.biomass  ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model1 <- predict(pgls_model, newdata = data)  

###
lm_model <- lm( MPR ~ Root.length.colonization , data = data) 
summary(lm_model)
pgls_model <- pgls(  MPR ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model1 <- predict(pgls_model, newdata = data)  


###
lm_model <- lm( MGR ~ Root.length.colonization , data = data) 
summary(lm_model)
pgls_model <- pgls(  MGR ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model2 <- predict(pgls_model, newdata = data)  



###
lm_model <- lm( MPR ~ Root.diameter , data = data) 
summary(lm_model)
pgls_model <- pgls(  MPR ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model3 <- predict(pgls_model, newdata = data)  


###
lm_model <- lm( MGR ~ Root.diameter , data = data) 
summary(lm_model)
pgls_model <- pgls(  MGR ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model4 <- predict(pgls_model, newdata = data)  


###
lm_model <- lm( MPR ~ MRR , data = data) 
summary(lm_model)
pgls_model <- pgls(  MPR ~ MRR , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model5 <- predict(pgls_model, newdata = data)  


###
lm_model <- lm( MGR ~ MRR , data = data) 
summary(lm_model)
pgls_model <- pgls(  MGR ~ MRR , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model6 <- predict(pgls_model, newdata = data)  


######## -- Fig 4  -- 225 45sp Plant’s responses to mycorrhizal inoculation ----------------
library(ggplot2)
library(ggpubr)
library(eoffice)
library(ggtreeExtra) 
library(ggstar) 
library(ggtree) 
library(treeio)
library(ggnewscale) 
library(aplot)
library(reshape2)
library(vegan)
library(ape)
library(plyr)
library(dplyr)
library(doBy)
library(RColorBrewer)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf")
data = read.csv("greenhouse_data_2.csv",header=TRUE, row.names=1)

##
rooted_tree <- read.tree("1_45sp_phylomakers1_abbr.nwk")
plot(rooted_tree)

is.rooted(rooted_tree)    
is.binary.tree(rooted_tree)  

se <- function(x) sd(x)/sqrt(length(x))
conf_95 <- function(x) t.test(x, conf.level = 0.95)$conf.int[2]
sumdata <- summaryBy( MPR ~ Plant.species, data, FUN = c(mean, sd, se, conf_95))
sumdata$Plant.family <- data$Plant.family[match(sumdata$Plant.species, data$Plant.species)]
group_file <- sumdata[,c(1, 6)]
row.names(group_file) <- group_file$Plant.species
groupInfo_list <- list()
groupInfo_list[["Plant.family"]] <- split(row.names(group_file), group_file$Plant.family)

for (group_name in names(groupInfo_list)) {
  tree <- groupOTU(rooted_tree, groupInfo_list[[group_name]], group_name)
}

tree.data = fortify(tree)
head(tree.data)

cols <- c("Amaryllidaceae"= "#5F97C6", "Cyperaceae" = "#A3B150", "Poaceae" = "#8FB9D1", 
          "Ranunculaceae" = "#B85A68", 
          # 
          "Rosaceae" =  "#FFA64D", "Fabaceae" = "#5880ba" , 
          "Brassicaceae"= "#9e5c93", "Onagraceae" = "#7B8F5E", 
          # 
          "Caryophyllaceae" = "#A83270" ,  "Polygonaceae"= "#9BBF6E",  
          "Gentianaceae"= "#E0A16F", "Lamiaceae"  = "#B07AA8", "Plantaginaceae" = "#497d75", 
          "Apiaceae" = "#BFA87F", "Asteraceae" = "#E9746A") 


p00 <- ggtree( tree, layout="rectangular", ladderize = T, branch.length = "none",size=0.7) +  
  geom_tippoint(aes(color= Plant.family ),  size = 2 , alpha = 1, position = position_nudge(x = 0.01)) + 
  geom_tiplab(aes(color= Plant.family), size=4, fontface = 3, hjust = 0, alpha = 1, offset=0.2, show.legend=F) + 
  scale_color_manual(values =cols) +  
  theme(legend.position="none", legend.title=element_text(size=11), legend.text =element_text(size=9),
        axis.title =element_text(size=10)) + labs( color ="Plant family")+ geom_rootedge(rootedge = 0) +
  guides(color = guide_legend(override.aes = list(size = 1.5, shape = 16), order = 1, nrow = 12))  

p00

results1 <- data %>% 
  dplyr::group_by(Plant.species, Plant.family) %>%
  dplyr::summarise(
    mean_MPR = mean(MPR, na.rm = TRUE),
    se_MPR = sd(MPR, na.rm = TRUE) / sqrt(n()),
    t_test = list(t.test(MPR, mu = 0)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(mean_MPR) %>%
  mutate(Plant.species = factor(Plant.species, levels = Plant.species))


p1 <- ggplot(results1, aes(x = Plant.species, y = mean_MPR, fill = Plant.family)) +
  geom_col(width = 0.6, color = "black", alpha = 1) +
  geom_errorbar(aes(ymin = mean_MPR - se_MPR, ymax = mean_MPR + se_MPR),
                width = 0 , linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6, linetype = "solid") +
  geom_text(aes(y = mean_MPR * 1.4, label = significance),
            size = 4.5, color = "black", vjust = 0.8) +
  scale_fill_manual(values = cols) +
  labs(x = "", y = "MPR") +
  coord_flip() +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        legend.position = "none") 

p1


results2 <- data %>%
  dplyr::group_by(Plant.species, Plant.family) %>%
  dplyr::summarise(
    mean_MGR = mean(MGR, na.rm = TRUE),
    se_MGR = sd(MGR, na.rm = TRUE) / sqrt(n()),
    t_test = list(t.test(MGR, mu = 0)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(mean_MGR) %>%
  mutate(Plant.species = factor(Plant.species, levels = Plant.species))

p2 <- ggplot(results2, aes(x = Plant.species, y = mean_MGR, fill = Plant.family)) +
  geom_col(width = 0.6, color = "black", alpha = 1) +
  geom_errorbar(aes(ymin = mean_MGR - se_MGR, ymax = mean_MGR + se_MGR),
                width = 0 , linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6, linetype = "solid") +
  geom_text(aes(y = mean_MGR * 1.5, label = significance),
            size = 4.5, color = "black", vjust = 0.8) +
  scale_fill_manual(values = cols) +
  labs(x = "", y = "MGR") +
  coord_flip() +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        legend.position = "none") 

p2


results3 <- data %>%
  dplyr::group_by(Plant.species, Plant.family) %>%
  dplyr::summarise(
    mean_MRR = mean(MRR, na.rm = TRUE),
    se_MRR = sd(MRR, na.rm = TRUE) / sqrt(n()),
    t_test = list(t.test(MRR, mu = 0)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = sapply(t_test, function(x) x$p.value),
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01 ~ "**",
      p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(mean_MRR) %>%
  mutate(Plant.species = factor(Plant.species, levels = Plant.species))

p3 <- ggplot(results3, aes(x = Plant.species, y = mean_MRR, fill = Plant.family)) +
  geom_col(width = 0.6, color = "black", alpha = 1) +
  geom_errorbar(aes(ymin = mean_MRR - se_MRR, ymax = mean_MRR + se_MRR),
                width = 0, linewidth = 0.6) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6, linetype = "solid") +
  geom_text(aes(y = mean_MRR * 1.3, label = significance),
            size = 4.5, color = "black", vjust = 0.8) +
  scale_fill_manual(values = cols) +
  labs(x = "", y = "MRR") +
  coord_flip() +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 12),
        legend.position = "none") 
p3


blank <- ggplot() + theme_void()



fig4 <- p1 %>%
  insert_left(p00, width = 0.8) %>%          
  insert_right(blank, width = 0.03) %>%      
  insert_right(p2,   width = 1) %>%         
  insert_right(blank, width = 0.03) %>%     
  insert_right(p3,   width = 1)

fig4

#eoffice::topptx( fig4 , width = 10 ,  height = 9 , "1_225_45sp_fig4_10_9.pptx" ) 


######## -- Regression relationships of MPR and MGR with root diameter under nonmycorrhizal conditions ----------------
library(car) 
library(ggplot2)
library(ggpmisc)
library(scales)
library(ggpubr)
library(cowplot)
library(ppcor)
library(doBy)
library(caper)
rm(list = ls())  
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf")
data = read.csv("greenhouse_data_mean.csv",header=TRUE, row.names=1) 

##
rooted_tree <- read.tree("1_45sp_phylomakers1.nwk")
plot(rooted_tree)

data$Plant.species0 <- gsub(" ", "_", data$Plant.species0)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data , names.col = "Plant.species0", vcv = TRUE, vcv.dim = 2)
nrow(comp_data$data)
colnames(data)

###
lm_model <- lm( MPR ~ Root.diameter_NM , data = data) 
summary(lm_model)
pgls_model <- pgls(  MPR ~ Root.diameter_NM , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model1 <- predict(pgls_model, newdata = data)  


###
lm_model <- lm( MGR ~ Root.diameter_NM , data = data) 
summary(lm_model)
pgls_model <- pgls(  MGR ~ Root.diameter_NM , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model2 <- predict(pgls_model, newdata = data)  


##
lm_model <- lm( P.concentrations ~ Root.length.colonization, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model2 <- predict(pgls_model, newdata = data )  


########################################################################################################################################
######## Experiment III global synthesis ----------------------------
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf/global synthesis")
######## Median calculation ----------------
library(plyr)
library(dplyr)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf/global synthesis/median calculation") 
data = read.csv("rd_4101records.csv",header=TRUE, row.names=1)
colnames(data) 

check_family <- data %>%
  group_by(PlantSpecies) %>%
  summarise(
    n_family = n_distinct(PlantFamily),
    families = paste(unique(PlantFamily), collapse = "; ")
  ) %>%
  filter(n_family > 1)
check_family

data_median <- ddply(
  data,
  "PlantSpecies",
  summarise,
  PlantFamily = unique(PlantFamily),
  RootDiameter_median = median(RootDiameter, na.rm = TRUE),
  FirstQuartile = quantile(RootDiameter, 0.25, na.rm = TRUE),
  ThirdQuartile = quantile(RootDiameter, 0.75, na.rm = TRUE),
  IQR = IQR(RootDiameter, na.rm = TRUE),
  SampleSize = sum(!is.na(RootDiameter)) )
data_median

write.csv(data_median,"0_rd_898sp.csv")


###
data = read.csv("rlc_5683records.csv",header=TRUE, row.names=1)
colnames(data) 

check_family <- data %>%
  group_by(PlantSpecies) %>%
  summarise(
    n_family = n_distinct(PlantFamily),
    families = paste(unique(PlantFamily), collapse = "; ")
  ) %>%
  filter(n_family > 1)
check_family


data_median <- ddply(
  data,
  "PlantSpecies",
  summarise,
  PlantFamily = unique(PlantFamily),
  RootLengthColonization_median = median(RootLengthColonization, na.rm = TRUE),
  FirstQuartile = quantile(RootLengthColonization, 0.25, na.rm = TRUE),
  ThirdQuartile = quantile(RootLengthColonization, 0.75, na.rm = TRUE),
  IQR = IQR(RootLengthColonization, na.rm = TRUE),
  SampleSize = sum(!is.na(RootLengthColonization)))
data_median

write.csv(data_median,"0_rlc_1794sp.csv")


###
data = read.csv("mgr_2104records.csv",header=TRUE, row.names=1)
colnames(data) 

check_family <- data %>%
  group_by(PlantSpecies) %>%
  summarise(
    n_family = n_distinct(PlantFamily),
    families = paste(unique(PlantFamily), collapse = "; ")
  ) %>%
  filter(n_family > 1)
check_family

data_median <- ddply(
  data,
  "PlantSpecies",
  summarise,
  PlantFamily = unique(PlantFamily),
  MycorrhizalGrowthResponse_median = median(MycorrhizalGrowthResponse, na.rm = TRUE),
  FirstQuartile = quantile(MycorrhizalGrowthResponse, 0.25, na.rm = TRUE),
  ThirdQuartile = quantile(MycorrhizalGrowthResponse, 0.75, na.rm = TRUE),
  IQR = IQR(MycorrhizalGrowthResponse, na.rm = TRUE),
  SampleSize = sum(!is.na(MycorrhizalGrowthResponse)) )
data_median

write.csv(data_median,"0_mgr_348sp.csv")

######## Data pairing ------------------------------
library(dplyr)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf/global synthesis/median calculation") 
rd <- read.csv("0_rd_898sp.csv")
rlc <- read.csv("0_rlc_1794sp.csv")

# 提取两个数据集共有物种，并合并数据
rd_rlc_common <- inner_join(
  rd %>% dplyr::select(PlantSpecies, PlantFamily, RootDiameter_median,
                       FirstQuartile, ThirdQuartile, IQR, SampleSize),
  rlc %>% dplyr::select(PlantSpecies, PlantFamily, RootLengthColonization_median,
                        FirstQuartile, ThirdQuartile, IQR, SampleSize),
  by = "PlantSpecies",
  suffix = c("_RD", "_RLC")
)

# 查看结果
head(rd_rlc_common)
nrow(rd_rlc_common)

write.csv( rd_rlc_common,  "1_rd_rlc_395sp.csv",  row.names = FALSE)


###
rd <- read.csv("0_rd_898sp.csv")
mgr <- read.csv("0_mgr_348sp.csv")

# 提取两个数据集共有物种，并合并数据
rd_rlc_common <- inner_join(
  rd %>% dplyr::select(PlantSpecies, PlantFamily, RootDiameter_median,
                       FirstQuartile, ThirdQuartile, IQR, SampleSize),
  mgr %>% dplyr::select(PlantSpecies, PlantFamily, MycorrhizalGrowthResponse_median,
                        FirstQuartile, ThirdQuartile, IQR, SampleSize),
  by = "PlantSpecies",
  suffix = c("_RD", "_MGR")
)

# 查看结果
head(rd_rlc_common)
nrow(rd_rlc_common)

write.csv( rd_rlc_common,  "1_rd_mgr_162sp.csv",  row.names = FALSE)


###
rlc <- read.csv("0_rlc_1794sp.csv")
mgr <- read.csv("0_mgr_348sp.csv")

# 提取两个数据集共有物种，并合并数据
rd_rlc_common <- inner_join(
  rlc %>% dplyr::select(PlantSpecies, PlantFamily, RootLengthColonization_median,
                        FirstQuartile, ThirdQuartile, IQR, SampleSize),
  mgr %>% dplyr::select(PlantSpecies, PlantFamily, MycorrhizalGrowthResponse_median,
                        FirstQuartile, ThirdQuartile, IQR, SampleSize),
  by = "PlantSpecies",
  suffix = c("_RLC", "_MGR")
)

# 查看结果
head(rd_rlc_common)
nrow(rd_rlc_common)

write.csv( rd_rlc_common,  "1_rlc_mgr_211sp.csv",  row.names = FALSE)

##  v.phylomaker tree rlc 1794sp rd 898sp mgr 348sp --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf/global synthesis") 
df <- read.csv("1794sp.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "0_rlc_1794sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  


####
rm(list = ls())
df <- read.csv("898sp.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "0_rd_898sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL


is.rooted(sp_tr)       
is.binary.tree(sp_tr)  


####
rm(list = ls())

df <- read.csv("348sp.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "0_mgr_348sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL


is.rooted(sp_tr)       
is.binary.tree(sp_tr)  


##  v.phylomaker tree   395sp  211sp  162sp --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf/global synthesis") 
df <- read.csv("395sp.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "1_rd_rlc_395sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  


####
rm(list = ls())

df <- read.csv("211sp.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "1_rlc_mgr_211sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  


####
rm(list = ls())
df <- read.csv("162sp.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "1_rd_mgr_162sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  


######## -- Table S2 -- Plant phylogenetic signals ----------------
library(ggtree)
library(ape)
library(picante)
library(phylosignal)
library(phylobase)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf/global synthesis")  
data <- read.csv("0_mgr_348sp.csv", header=TRUE, row.names=1)
rooted_tree <- read.tree("0_mgr_348sp_phylomakers1.nwk")

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

#all(sort(rooted_tree$label) == sort(rownames(data)))
tree_p <- as(rooted_tree, "phylo4")
class(tree_p)

p4d <- phylo4d(tree_p, data)
phy.sig <- phyloSignal(p4d, reps = 9999)
phy.sig


####
data <- read.csv("0_rd_898sp.csv", header=TRUE, row.names=1)
rooted_tree <- read.tree("0_rd_898sp_phylomakers1.nwk")

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

#all(sort(rooted_tree$label) == sort(rownames(data)))

tree_p <- as(rooted_tree, "phylo4")
class(tree_p)

p4d <- phylo4d(tree_p, data)
phy.sig <- phyloSignal(p4d, reps = 9999)
phy.sig


####
data <- read.csv("0_rlc_1794sp.csv", header=TRUE, row.names=1)
rooted_tree <- read.tree("0_rlc_1794sp_phylomakers1.nwk")

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

#all(sort(rooted_tree$label) == sort(rownames(data)))

tree_p <- as(rooted_tree, "phylo4")
class(tree_p)

p4d <- phylo4d(tree_p, data)
phy.sig <- phyloSignal(p4d, reps = 9999)
phy.sig


######## global synthesis Sensitivity analysis between the all-species dataset and the typical-AM-species dataset ----------------------------
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf/global synthesis/global synthesis_d3") 
##  v.phylomaker tree    326sp 203sp  155sp --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf/global synthesis/global synthesis_d3") 
df <- read.csv("326_rd_rlc.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "1_rd_rlc_326sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  


####
rm(list = ls())
df <- read.csv("203_rlc_mgr.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "1_rlc_mgr_203sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  



####
rm(list = ls())
df <- read.csv("155_rd_mgr.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1")
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "1_rd_mgr_155sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  


######## -- Table S3 -- Regression relationships among root diameter, root length colonization, and mycorrhizal growth response ----------------
library(readxl) 
library(pheatmap) 
library(psych)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(car)
library(caper)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf/global synthesis/global synthesis_d3") 
### 395 sp  
data1 = read.csv("1_rd_rlc_326sp.csv", header=TRUE)
##
rooted_tree <- read.tree("1_rd_rlc_326sp_phylomakers1.nwk")
plot(rooted_tree)

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data1 , names.col = "PlantSpecies", vcv = TRUE, vcv.dim = 2)
nrow(comp_data$data)

colnames(data1)

### log(RootLengthColonization+1) ~ log(RootDiameter)
# degree = 2 表示二次多项式
#lm_model <- lm( log(RootLengthColonization+1) ~ poly(log(RootDiameter), 2, raw = TRUE), data = data1 )  
#pgls_model <- pgls( log(RootLengthColonization+1) ~ poly(log(RootDiameter), 2, raw = TRUE),   data = comp_data, lambda = "ML")
#lm_model <- lm( RootLengthColonization ~ RootDiameter , data = data1 ) 
#pgls_model <- pgls( RootLengthColonization ~ RootDiameter , data = comp_data,  lambda = "ML") 

lm_model <- lm( log(RootLengthColonization+1) ~ log(RootDiameter) , data = data1 ) 
summary(lm_model)
pgls_model <- pgls( log(RootLengthColonization+1) ~ log(RootDiameter) , data = comp_data,  lambda = "ML") 
summary(pgls_model)
pgls_model1 <- predict(pgls_model, newdata = data1 )  


### 211 sp 
data2 = read.csv("1_rlc_mgr_203sp.csv", header=TRUE)
##
rooted_tree <- read.tree("1_rlc_mgr_203sp_phylomakers1.nwk")
plot(rooted_tree)

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data2 <- comparative.data( phy = rooted_tree, data = data2 , names.col = "PlantSpecies", vcv = TRUE, vcv.dim = 2)
nrow(comp_data2$data)
colnames(data2)

###
lm_model <- lm( MycorrhizalGrowthResponse ~ log(RootLengthColonization+1)  , data = data2 ) 
summary(lm_model)
pgls_model <- pgls( MycorrhizalGrowthResponse ~ log(RootLengthColonization+1)  , data = comp_data2,  lambda = "ML") 
summary(pgls_model)
pgls_model2 <- predict(pgls_model, newdata = data2 )  


### 162 sp 
data3 = read.csv("1_rd_mgr_155sp.csv", header=TRUE)
##
rooted_tree <- read.tree("1_rd_mgr_155sp_phylomakers1.nwk")
plot(rooted_tree)

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data3 <- comparative.data( phy = rooted_tree, data = data3 , names.col = "PlantSpecies", vcv = TRUE, vcv.dim = 2)
nrow(comp_data3$data)
colnames(data3)

###
lm_model <- lm( MycorrhizalGrowthResponse ~ log(RootDiameter) , data = data3 ) 
summary(lm_model)
pgls_model <- pgls( MycorrhizalGrowthResponse ~ log(RootDiameter) , data = comp_data3,  lambda = "ML") 
summary(pgls_model)
pgls_model3 <- predict(pgls_model, newdata = data3 )    



########################################################################################################################################
######## -- Fig 1  -- Relationship between plant traits and AM fungal colonization  ---------------
library(ggplot2) 
library(gridExtra)
library(reshape)
library(ggpmisc)
library(car)
library(caper)
library(doBy)
library(dplyr) 
library(ape)
rm(list = ls()) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf") 
data1 = read.csv("field_data_155_35sp.csv",header=TRUE, row.names=1) 

rooted_tree <- read.tree("1_35sp_phylomakers1_abbr.nwk") 
plot(rooted_tree) 

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data1 , names.col = "Plant.species", vcv = TRUE, vcv.dim = 2)
nrow(comp_data$data)
colnames(data1)

##
lm_model <- lm( Root.length.colonization ~ Root.diameter , data = data1 ) 
summary(lm_model)
pgls_model <- pgls(  Root.length.colonization ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model1 <- predict(pgls_model, newdata = data1 )  

phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Root.diameter), max(dat0$Root.diameter), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( Root.length.colonization ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

pgls_ci11 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci1)

fig1a <- ggplot( data1 , aes(y = Root.length.colonization , x = Root.diameter  ))+   labs(title = "Field: 35 species") +   
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "Root length colonization (%)" ) +
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22",  fill = "grey88",  formula = y ~ poly(x, 1, raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  #geom_ribbon(data = pgls_ci11,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model1 ), linewidth = 0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0), axis.text = element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) + 
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.0, vjust = - 2, color = "grey22",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.132, 0.032 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.0, vjust = - 0.85, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.073, 0.116 ), size = 4, inherit.aes = FALSE )
fig1a


###
lm_model <- lm( P.concentrations ~ Root.length.colonization , data = data1 ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model1 <- predict(pgls_model, newdata = data1 )  



######## -- Fig 1 45sp -- Relationship between plant traits and AM fungal colonization 
data = read.csv("greenhouse_data_mean.csv", header=TRUE, row.names=1)

##
rooted_tree <- read.tree("1_45sp_phylomakers1.nwk")
plot(rooted_tree)

data$Plant.species0 <- gsub(" ", "_", data$Plant.species0)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data , names.col = "Plant.species0", vcv = TRUE, vcv.dim = 2)
nrow(comp_data$data)
colnames(data)

###
lm_model <- lm( Root.length.colonization ~ Root.diameter , data = data) 
summary(lm_model)
pgls_model <- pgls(  Root.length.colonization ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model1 <- predict(pgls_model, newdata = data)  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species0 <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Root.diameter), max(dat0$Root.diameter), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species0", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( Root.length.colonization ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci1 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig1b <- ggplot(data, aes(y = Root.length.colonization , x = Root.diameter  ))+  labs(title = "Greenhouse: 45 species") +     
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "Root length colonization (%)" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  #geom_ribbon(data = pgls_ci1,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model1 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 1.02, colour = "grey22", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.098, 0.036 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 2.1, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.005, 0.629 ), size = 4, inherit.aes = FALSE )

fig1b


####
### 395 sp 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf/global synthesis") 
data1 = read.csv("1_rd_rlc_395sp.csv",header=TRUE)
##
rooted_tree <- read.tree("1_rd_rlc_395sp_phylomakers1.nwk")
plot(rooted_tree)

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data1 , names.col = "PlantSpecies", vcv = TRUE, vcv.dim = 2)
nrow(comp_data$data)
colnames(data1)

###
lm_model <- lm( log(RootLengthColonization+1) ~ log(RootDiameter) , data = data1 ) 
summary(lm_model)
pgls_model <- pgls( log(RootLengthColonization+1) ~ log(RootDiameter) , data = comp_data,  lambda = "ML") 
summary(pgls_model)
pgls_model1 <- predict(pgls_model, newdata = data1 )  


##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$PlantSpecies <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$RootDiameter, na.rm = TRUE), max(dat0$RootDiameter, na.rm = TRUE), length.out = 100 )

boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))

for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "PlantSpecies", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try(  pgls(log(RootLengthColonization+1) ~ log(RootDiameter),   data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( RootDiameter = pred_grid) ) } } }

# CI
pgls_ci1 <- data.frame( RootDiameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig1c <- ggplot( data1, aes( y = log(RootLengthColonization+1) , x = log(RootDiameter)  ))+    labs(title = "Global synthesis: 395 species") +    
  scale_y_continuous(name ="ln(Root length colonization [%] +1)") + scale_x_continuous(name ="ln(Root diameter [mm])") + 
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 1 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci1 ,  aes(x = log(RootDiameter) , ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model1 ), linewidth=1, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33", face = "plain", hjust = 0), axis.text=element_text(size = 10.5 ), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size = 10.5 ), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1 ))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 1.0, colour = "grey22",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.077, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 2.05, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'<'~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.001, 0.729 ), size = 4, inherit.aes = FALSE )

fig1c


#
fig1 <- cowplot::plot_grid( fig1a,  fig1b, fig1c,  
                            ncol = 3, nrow = 1, labels = c('(a)', "(b)", "(c) ",'(d)', "(e)"), 
                            label_fontface = "bold", label_x = -0.03, label_y = 1.01 , align = c("hv"))
fig1


setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf")
eoffice::topptx( fig1 , width = 9 ,  height= 3.5 , "1_fig1_9_3.5.pptx" ) 




######## -- Fig 5 -- Regression relationships of MPR and MGR with root length colonization, root diameter, and MRR ----------------
library(car) 
library(ggplot2)
library(ggpmisc)
library(scales)
library(ggpubr)
library(cowplot)
library(ppcor)
library(doBy)
library(caper)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf") 
data = read.csv("greenhouse_data_mean.csv",header=TRUE, row.names=1) 

##
rooted_tree <- read.tree("1_45sp_phylomakers1.nwk")
plot(rooted_tree)

data$Plant.species0 <- gsub(" ", "_", data$Plant.species0)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data , names.col = "Plant.species0", vcv = TRUE, vcv.dim = 2)
nrow(comp_data$data)
colnames(data)

###
lm_model <- lm( MPR ~ Root.length.colonization , data = data) 
summary(lm_model)
pgls_model <- pgls(  MPR ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model1 <- predict(pgls_model, newdata = data)  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species0 <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Root.length.colonization), max(dat0$Root.length.colonization), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species0", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( MPR ~ Root.length.colonization ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.length.colonization = pred_grid) ) } } }

# CI
pgls_ci1 <- data.frame( Root.length.colonization = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig5a <- ggplot(data, aes(y = MPR , x = Root.length.colonization  )) +  labs(title = "Greenhouse: 45 species") +    
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = "MPR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  geom_ribbon(data = pgls_ci1,  aes(x = Root.length.colonization, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model1 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33",  face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 1.01 ,  colour = "grey22", 
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.299, 0.001 ), size = 4 , inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 2.05 , color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.211, 0.001 ), size = 4 , inherit.aes = FALSE )

fig5a


###
lm_model <- lm( MGR ~ Root.length.colonization , data = data) 
summary(lm_model)
pgls_model <- pgls(  MGR ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model2 <- predict(pgls_model, newdata = data)  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species0 <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Root.length.colonization), max(dat0$Root.length.colonization), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species0", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( MGR ~ Root.length.colonization ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.length.colonization = pred_grid) ) } } }

# CI
pgls_ci2 <- data.frame( Root.length.colonization = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig5b <- ggplot(data, aes(y = MGR , x = Root.length.colonization  ))+     labs(title = "Greenhouse: 45 species") + 
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = "MGR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  #geom_ribbon(data = pgls_ci2,  aes(x = Root.length.colonization, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model2 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33",  face = "plain", hjust = 0),axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text( size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 1.01 , color = "grey22",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.258, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 2.05 , color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.019, 0.362 ), size = 4, inherit.aes = FALSE )

fig5b


###
lm_model <- lm( MPR ~ Root.diameter , data = data) 
summary(lm_model)
pgls_model <- pgls(  MPR ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model3 <- predict(pgls_model, newdata = data)  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species0 <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Root.diameter), max(dat0$Root.diameter), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species0", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( MPR ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci3 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig5c <- ggplot(data, aes(y = MPR , x = Root.diameter  )) + labs(title = "Greenhouse: 45 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "MPR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  geom_ribbon(data = pgls_ci3,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model3 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33",  face = "plain", hjust = 0), axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 1.01, color = "grey22",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.197, 0.002 ), size = 4 , inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 2.05, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.120, 0.020 ), size = 4 , inherit.aes = FALSE )

fig5c


###
lm_model <- lm( MGR ~ Root.diameter , data = data) 
summary(lm_model)
pgls_model <- pgls(  MGR ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model4 <- predict(pgls_model, newdata = data)  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species0 <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Root.diameter), max(dat0$Root.diameter), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species0", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( MGR ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci4 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig5d <- ggplot(data, aes(y = MGR , x = Root.diameter  ))+      labs(title = "Greenhouse: 45 species") +
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "MGR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  #geom_ribbon(data = pgls_ci4,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model4 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33",  face = "plain", hjust = 0), axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.01, vjust = - 1.4, color = "grey22",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.258, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.01, vjust = - 0.25, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.061, 0.101 ), size = 4, inherit.aes = FALSE )

fig5d



###
lm_model <- lm( MPR ~ MRR , data = data) 
summary(lm_model)
pgls_model <- pgls(  MPR ~ MRR , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model5 <- predict(pgls_model, newdata = data)  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species0 <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$MRR), max(dat0$MRR), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE) 
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species0", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( MPR ~ MRR ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( MRR = pred_grid) ) } } }

# CI
pgls_ci5 <- data.frame( MRR = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig5e <- ggplot(data, aes(y = MPR , x = MRR  )) +      labs(title = "Greenhouse: 45 species") +
  scale_x_continuous(name ="MRR") + scale_y_continuous(name = "MPR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  #geom_ribbon(data = pgls_ci5,  aes(x = MRR, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model5 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33",  face = "plain", hjust = 0), axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 1.0, color = "grey22",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.193, 0.003 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 2.05, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.072, 0.074 ), size = 4, inherit.aes = FALSE )

fig5e


###
lm_model <- lm( MGR ~ MRR , data = data) 
summary(lm_model)
pgls_model <- pgls(  MGR ~ MRR , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model6 <- predict(pgls_model, newdata = data)  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species0 <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$MRR), max(dat0$MRR), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species0", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( MGR ~ MRR ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( MRR = pred_grid) ) } } }

# CI
pgls_ci6 <- data.frame( MRR = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig5f <- ggplot(data, aes(y = MGR , x = MRR  ))+      labs(title = "Greenhouse: 45 species") +
  scale_x_continuous(name ="MRR") + scale_y_continuous(name = "MGR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  geom_ribbon(data = pgls_ci6,  aes(x = MRR, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model6 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33",  face = "plain", hjust = 0), axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.15, size = 13),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", nrow= 12, label.hjust = 0))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.01, vjust = - 1.4, color = "grey22",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.617, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.01, vjust = - 0.25, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.489, 0.001 ), size = 4, inherit.aes = FALSE )

fig5f


######## -- Fig 5 -- Regression relationships among root diameter, root length colonization,  and mycorrhizal growth response 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf/global synthesis")

### 211 sp 
data2 = read.csv("1_rlc_mgr_211sp.csv",header=TRUE)
#data2.1 = read.csv("1_rlc_mgr_211sp_2.csv",header=TRUE)
##
rooted_tree <- read.tree("1_rlc_mgr_211sp_phylomakers1.nwk")
plot(rooted_tree)

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data2 <- comparative.data( phy = rooted_tree, data = data2 , names.col = "PlantSpecies", vcv = TRUE, vcv.dim = 2)
nrow(comp_data2$data)
colnames(data2)

###
lm_model <- lm( MycorrhizalGrowthResponse ~ log(RootLengthColonization+1)  , data = data2 ) 
summary(lm_model)
pgls_model <- pgls( MycorrhizalGrowthResponse ~ log(RootLengthColonization+1)  , data = comp_data2,  lambda = "ML") 
summary(pgls_model)
pgls_model7 <- predict(pgls_model, newdata = data2 )  

##
phy  <- comp_data2$phy
dat0 <- comp_data2$data
dat0$PlantSpecies <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$RootLengthColonization, na.rm = TRUE), max(dat0$RootLengthColonization, na.rm = TRUE), length.out = 100 )

boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))

for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "PlantSpecies", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try(  pgls( MycorrhizalGrowthResponse ~ log(RootLengthColonization+1) ,   data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( RootLengthColonization = pred_grid) ) } } }

# CI
pgls_ci7 <- data.frame( RootLengthColonization = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )

fig5g <- ggplot(data2, aes(y = MycorrhizalGrowthResponse  , x = log(RootLengthColonization+1)  ))+    labs(title = "Global synthesis: 211 species") +   
  scale_y_continuous(name ="MGR") + scale_x_continuous(name ="ln(Root length colonization [%] +1)") + 
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, colour = "grey22", fill = "grey88", data = data2, formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci7 ,  aes(x = log(RootLengthColonization+1) , ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model7 ), linewidth= 0.8 , color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33",  face = "plain", hjust = 0), axis.text=element_text( size=10.5 ), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text( size=10.5 ), legend.title = element_text(hjust=0.8),  axis.title=element_text( size = 12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 0.95,  colour = "grey22",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.071, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.02, vjust = 2.0, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.052, 0.001 ), size = 4, inherit.aes = FALSE )

fig5g


### 162 sp 
data3 = read.csv("1_rd_mgr_162sp.csv",header=TRUE)
##
rooted_tree <- read.tree("1_rd_mgr_162sp_phylomakers1.nwk")
plot(rooted_tree)

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data3 <- comparative.data( phy = rooted_tree, data = data3 , names.col = "PlantSpecies", vcv = TRUE, vcv.dim = 2)
nrow(comp_data3$data)
colnames(data3)

###
lm_model <- lm( MycorrhizalGrowthResponse ~ log(RootDiameter) , data = data3 ) 
summary(lm_model)
pgls_model <- pgls( MycorrhizalGrowthResponse ~ log(RootDiameter) , data = comp_data3,  lambda = "ML") 
summary(pgls_model)
pgls_model8 <- predict(pgls_model, newdata = data3 )    

##
phy  <- comp_data3$phy
dat0 <- comp_data3$data
dat0$PlantSpecies <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$RootDiameter, na.rm = TRUE), max(dat0$RootDiameter, na.rm = TRUE), length.out = 100 )

boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))

for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "PlantSpecies", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try(  pgls( MycorrhizalGrowthResponse ~ log(RootDiameter),   data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict( m ,  newdata = data.frame( RootDiameter = pred_grid) ) } } }

# CI
pgls_ci8 <- data.frame( RootDiameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig5h <- ggplot( data3, aes(y = MycorrhizalGrowthResponse , x = log(RootDiameter)  ))+    labs(title = "Global synthesis: 162 species") +      
  scale_y_continuous(name ="MGR") + scale_x_continuous(name ="ln(Root diameter [mm])") + 
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, colour = "grey22", fill = "grey88", data = data3, formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci8 ,  aes(x = log(RootDiameter), ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes( y = pgls_model8 ), linewidth = 0.8, color = "#CB4641", linetype = "solid" ) + 
  theme_classic() + theme(plot.title = element_text( size = 12.5, color= "grey33",  face = "plain", hjust = 0), axis.text = element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) + 
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = -1.25, color = "grey22",
             label = sprintf("'OLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.092, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = -0.2, color = "#CB4641",
             label = sprintf("'PGLS'*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.040, 0.011 ), size = 4, inherit.aes = FALSE )

fig5h



fig5 <- cowplot::plot_grid( fig5a,  fig5b, fig5c, 
                            fig5d, fig5e, fig5f,
                            fig5g, fig5h,
                            ncol = 3, nrow = 3, labels = c( "(a)", "(b)",'(c)',
                                                            "(d)", "(e)","(f)",  
                                                            "(g)", "(h)", "(i)",
                                                            "(j)", "(k)", "(l)" ), 
                            label_fontface = "bold", label_x = -0.0, label_y = 1.01 , align = c("hv"))
fig5


setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf") 
eoffice::topptx( fig5 , width = 9 ,  height = 9.6 , "1_fig5_9_9.6.pptx" ) 


########################################################################################################################################
###### Calculate the coefficient of variation of root diameter --------------
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf")
data1 = read.csv("field_data_155_35sp.csv",header=TRUE, row.names=1) 
colnames(data1)

data1$Root.diameter

mean_RD <- mean(data1$Root.diameter, na.rm = TRUE)
sd_RD <- sd(data1$Root.diameter, na.rm = TRUE)
CV_RD <- sd_RD / mean_RD * 100
CV_RD # [1] 62.72251


data2 = read.csv("greenhouse_data_mean.csv",header=TRUE, row.names=1) 
colnames(data2)

mean_RD <- mean(data2$Root.diameter, na.rm = TRUE)
sd_RD <- sd(data2$Root.diameter, na.rm = TRUE)
CV_RD <- sd_RD / mean_RD * 100
CV_RD # [1] 25.36745


####
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/1_amf/global synthesis") 
### 395 sp 
data3 = read.csv("1_rd_rlc_395sp.csv",header=TRUE)
colnames(data3)

mean_RD <- mean(data3$RootDiameter, na.rm = TRUE)
sd_RD <- sd(data3$RootDiameter, na.rm = TRUE)
CV_RD <- sd_RD / mean_RD * 100
CV_RD #    [1] 65.41678
