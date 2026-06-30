######## Determination of the threshold for sample exclusion: 20th percentile -------------------------
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
otu_all = read.csv("field_data_391otu_raw_155.csv", header=TRUE, row.names=1)
colnames(otu_all)

otu = otu_all[ , -c( 1:12 )]


read_counts <- rowSums(otu)
summary(read_counts)
hist(read_counts, breaks = 100, main = "Sequencing depth distribution")


five_quantiles <- quantile(read_counts, probs = c(0, 0.2, 0.4, 0.6, 0.8, 1))
print(five_quantiles)

# 也可以在直方图上标记这些分位数 
hist(read_counts, breaks = 20 , main = "Sequencing depth distribution")
abline(v = five_quantiles, col = "red", lty = 2, lwd =  0.8 )
legend("topright", legend = names(five_quantiles), col = "red", lty = 2, lwd = 1.5, cex = 0.7)


# 1. 计算每个样品的总序列数
sample_reads <- rowSums(otu)
sample_reads <- sort(sample_reads, decreasing = TRUE)


library(ggplot2)
read_df <- data.frame(Sample = names(sample_reads), Reads = sample_reads)

ggplot(read_df, aes(x = reorder(Sample, -Reads), y = Reads)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_hline(yintercept = 9360, color = "red", linetype = "dashed", linewidth = 0.5) +
  labs(x = "Sample", y = "Total reads", 
       title = "Sequencing depth across samples") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



########  v.phylomaker tree  35sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
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


########  v.phylomaker tree  28sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
df <- read.csv("28sp.csv")
head(df)

sp_df <- data.frame(
  Species = df$Species,
  Genus = df$Genus,
  Family = df$Family )

tree <- phylo.maker( sp = sp_df, tree = GBOTB.extended.TPL, scenarios = "S1", output.sp.list = TRUE)
sp_tr <- tree$scenario.1
class(sp_tr)
sp_tr
plot(sp_tr)

write.tree(sp_tr, "1_28sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL
#write.tree(sp_tr, "1-35sp-phylomakers1-mega.nwk")

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  

####
library(dplyr)
tree <- read.tree("1_28sp_phylomakers1.nwk")
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

write.tree( tree, file = "1_28sp_phylomakers1_abbr.nwk")


setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys")  
data = read.csv("field_data_118.csv", header=TRUE, row.names=1)

data_mean <-  plyr::ddply(data, c("Plant.species", "Plant.species0", "Plant.family"),summarise,   
                          P.concentrations = mean(P.concentrations,na.rm=TRUE),
                          Root.diameter = mean(Root.diameter,na.rm=TRUE),
                          Root.length.colonization = mean(Root.length.colonization,na.rm=TRUE),
                          
                          Family.richness = mean(Family.richness,na.rm=TRUE),
                          Genus.richness = mean(Genus.richness,na.rm=TRUE),
                          OTU.richness = mean(OTU.richness,na.rm=TRUE),

                          Glomeraceae = mean( Glomeraceae ,na.rm=TRUE),
                          Acaulosporaceae = mean( Acaulosporaceae ,na.rm=TRUE),
                          Archaeosporaceae = mean( Archaeosporaceae, na.rm=TRUE),
                          Entrophosporaceae = mean( Entrophosporaceae, na.rm=TRUE),
                          Diversisporaceae = mean( Diversisporaceae, na.rm=TRUE),
                          
                          Dominikia = mean( Dominikia, na.rm=TRUE),
                          Rhizophagus = mean( Rhizophagus , na.rm=TRUE),
                          Acaulospora = mean( Acaulospora, na.rm=TRUE),
                          Polonospora = mean( Polonospora , na.rm=TRUE),
                          Funneliformis = mean( Funneliformis, na.rm=TRUE),
                          Entrophospora = mean( Entrophospora , na.rm=TRUE),
                          Kamienskia = mean( Kamienskia , na.rm=TRUE),
                          Microdominikia = mean( Microdominikia, na.rm=TRUE),
                          Others = mean( Others , na.rm=TRUE),
                          
                          
                          Module.1 = mean( Module.1, na.rm=TRUE),
                          Module.2 = mean( Module.2, na.rm=TRUE),
                          Module.3 = mean( Module.3, na.rm=TRUE),
                          Module.4 = mean( Module.4, na.rm=TRUE),
                          Module.5 = mean( Module.5, na.rm=TRUE),
                          
                          OTU1 = mean( OTU1, na.rm=TRUE),
                          OTU2 = mean( OTU2, na.rm=TRUE),
                          OTU3 = mean( OTU3, na.rm=TRUE),
                          OTU4 = mean( OTU4, na.rm=TRUE),
                          OTU5 = mean( OTU5, na.rm=TRUE)
                          
)

data_mean

data.table::fwrite(data_mean,"field_data_118_28sp.csv",  row.names = TRUE, bom = TRUE)


######## RA & Alpha_diversity_index calculation ---------------
library(amplicon)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys")
otu_all <- read.csv('field_data_386otu_rare_118.csv', header=T,row.names = 1)
colnames(otu_all)
rownames(otu_all) <- otu_all[, 12]
otu_t <- t(otu_all[ , -c(1:12)])

otu_ra <- sweep( otu_t, 2, colSums(otu_t),  "/" )

data.table::fwrite( otu_ra, "field_data_386otu_ra_118.csv",  row.names = TRUE, bom = TRUE)

##
Tax <- read.csv('field_data_386otu_taxon.csv', header=T, row.names = 1)

####
####
otu_all <- read.csv('field_data_386otu_rare_118.csv', header=T,row.names = 1)
colnames(otu_all)
rownames(otu_all) <- otu_all[, 12]
otu_t <- t(otu_all[ , -c(1:12)])

Tax <- read.csv('field_data_386otu_taxon.csv', header=T, row.names = 1)

format2stamp(otu_t, Tax, thre = 0.0, prefix = "RA_118") 

####
####
otu_all <- read.csv('field_data_386otu_rare_118_28sp.csv', header=T,row.names = 1)
colnames(otu_all)
rownames(otu_all) <- otu_all[, 2]
otu_t <- t(otu_all[ , -c(1:3)])

Tax <- read.csv('field_data_386otu_taxon.csv', header=T, row.names = 1)

format2stamp(otu_t, Tax, thre = 0.0, prefix = "RA_118_28sp") 


###
library(vegan)
library(picante)
otu_rare<- read.csv("field_data_386otu_rare_118.csv",header=TRUE, row.names=1)
colnames(otu_rare)
otu_rare_s <- otu_rare[,-c(1:13)] 

Alpha_diversity_index <- function(x, tree = NULL, base = exp(1)) {
  x <- as.matrix(x)
  
  # Richness: 丰度 > 0 的物种数（适用于小数）
  Obs <- rowSums(x > 0)
  
  # Shannon 和 Simpson (vegan 支持小数)
  Shannon <- vegan::diversity(x, index = 'shannon', base = base)
  Simpson <- vegan::diversity(x, index = 'simpson')
  
  # Pielou 均匀度
  Pielou <- Shannon / log(Obs, base)
  
  # Good's coverage: 对于小数数据不适用，改为计算相对丰度 > 0.001% 的物种占比
  # 或者直接跳过这个指标
  low_abundance <- rowSums(x > 0 & x < 0.0001)  # 丰度极低的物种数
  goods_coverage <- 1 - low_abundance / Obs
  
  # 组装结果
  result <- rbind(Obs, Shannon, Simpson, Pielou, goods_coverage)
  result <- as.data.frame(t(result))
  colnames(result) <- c("Richness", "Shannon", "Simpson", "Pielou", "Goods_coverage")
  
  if (!is.null(tree)) {
    # PD 需要整数，可以将小数转为 presence/absence
    x_pa <- ifelse(x > 0, 1, 0)
    Pd <- pd(x_pa, tree, include.root = FALSE)
    result$PD <- Pd$PD
  }
  
  return(result)
}

alpha_diversity <- Alpha_diversity_index(otu_rare_s)
dat <- cbind(otu_rare[,c(1:13)], alpha_diversity)

data.table::fwrite( dat,"1_field_data_386otu_rare_118_alpha.diversity.csv",  row.names = TRUE, bom = TRUE)


######## Co-occurrence network analysis -- 118samples ra -------------
library(psych)
library(reshape2)
library(igraph)
library(Hmisc)
rm(list = ls()) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
otu_a = read.csv("field_data_386otu_ra_118.csv",header=TRUE, row.names=1) 
otu = otu_a[ , -c(1:12)]
#sapply(otu, class)

set.seed(123)
cor <- rcorr(as.matrix(otu), type = "spearman")

r_raw <- cor$r
p_matrix <- cor$P
p_vector <- as.vector(p_matrix)
p <- p.adjust(p_vector, method = 'BH')

r <- r_raw
r[p > 0.05 | abs(r_raw) < 0.6] <- 0

diag(r) <- 0 

r <- (r + t(r)) / 2
r[is.na(r)] <- 0

write.csv(data.frame(as.data.frame(r), check.names = FALSE), '1_network_386otu_118_corrmatrix.csv')     

r[is.na(r)] <- 0
g = graph_from_adjacency_matrix(r, mode ="undirected",weighted=TRUE,diag = FALSE)
g
g = delete_vertices(g, names(igraph::degree(g)[igraph::degree(g) == 0]))

E(g)$corr = E(g)$weight       
E(g)$weight = abs(E(g)$weight)   

tax = read.csv('field_data_386otu_taxon.csv', row.names=1, header=T)    

tax = tax[as.character(V(g)$name), ]          
V(g)$Kingdom = tax$Kingdom                          
V(g)$Phylum = tax$Phylum                           
V(g)$Class = tax$Class                              
V(g)$Order = tax$Order                           
V(g)$Family = tax$Family                           
V(g)$Genus = tax$Genus                               
V(g)$Species = tax$Species 
V(g)$ra = tax$ra

node_list = data.frame(
  label = names(V(g)),
  kingdom = V(g)$Kingdom,
  phylum = V(g)$Phylum,
  class = V(g)$Class,
  order = V(g)$Order,
  family = V(g)$Family,
  genus=V(g)$Genus,
  species = V(g)$Species,
  ra = V(g)$ra  )   

head(node_list)

write.csv(node_list, '1_118_network_node.csv')   

edge = data.frame(as_edgelist(g)) 
edge_list = data.frame(
  source = edge[[1]],
  target = edge[[2]],
  weight = E(g)$weight,
  correlation = E(g)$corr
)

head(edge_list)

write.csv(edge_list, '1_118_network_edge.csv')  

write_graph(g, '1_network_386otu_118.graphml', format = 'graphml')                     

####
####
otu_a = read.csv("field_data_386otu_ra_118.csv", header=TRUE, row.names=1) 
rownames(otu_a) <- otu_a[ , 12 ]
otu = otu_a[ , -c(1:12)]
otu_t = t(otu)
head(otu_t)
module0 = read.csv("1_0624_network_218_358_module.csv", header=TRUE, row.names=1)

# Extract the OTU names from module0
otu_list <- module0$v_name

# Extract the corresponding OTU row from otu_t
otu_module <- otu_t[rownames(otu_t) %in% otu_list, ]

write.csv( otu_module , "1_0624_network_218_358_module0.csv", row.names = T)


####
module = read.csv("1_0624_network_218_358_module1.csv", header=TRUE, row.names=1)

module_by_species <- aggregate(. ~ modularity_class, data = module , FUN = sum , na.rm = TRUE)

head(module_by_species[, 1:5])

write.csv(module_by_species, "1_0624_network_218_358_module2.csv", row.names = FALSE)


######## -- Fig S1 -- rrarefy & rarecurve plot —— 118samples 28species  ---------------
library(vegan)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys")
otu_all = read.csv("field_data_386otu_raw_118.csv", header=TRUE, row.names=1)
otu1 = otu_all[ , -c( 1:12 )]

par(mar = c(5, 5, 1, 2))
rarecurve( otu1 ,  step = 10,  cex.lab = 1.2,  main = "122samples" ,
           label = F ,  lines = TRUE, col = "steelblue", lwd = 0.8,
           xlab = "Number of reads per sample",  ylab = "Number of OTU" )

###
set.seed(12345)
otu_rrarefy = rrarefy(otu1, min(rowSums(otu1)))  

otu_all_rrarefy = cbind(otu_all[1:12], otu_rrarefy)

data.table::fwrite(otu_all_rrarefy, "field_data_386otu_rare_118.csv",  row.names = TRUE, bom = TRUE)

##2_rarecurve_118samples_6_6
par(mar = c(5, 5, 1, 2))
rarecurve(otu_rrarefy,  step = 10,  cex.lab = 1.2,  main = NULL,
          label = FALSE,  lines = TRUE, col = "steelblue", lwd = 0.8,
          xlab = "Number of reads per sample",  ylab = "Number of OTU" )


# 按物种合并OTU数据（计算平均值）
otu2 = otu_all_rrarefy[ , -c(1:4, 7:12 )]

otu_by_species <- aggregate(. ~ Plant.species, 
                            data = otu2[, !names(otu2) %in% "Plant.family"], 
                            FUN = mean, na.rm = TRUE)


family_info <- unique(otu2[, c("Plant.species", "Plant.family")])
otu_by_species <- merge(otu_by_species, family_info, by = "Plant.species", all.x = TRUE)


otu_by_species <- otu_by_species[, c("Plant.species", "Plant.family", 
                                     setdiff(names(otu_by_species), 
                                             c("Plant.species", "Plant.family")))]
head(otu_by_species[, 1:5])

data.table::fwrite( otu_by_species , "field_data_386otu_rare_118_28sp.csv",  row.names = TRUE, bom = TRUE)


otu = otu_by_species
rownames(otu) <- otu[, 1]
otu <- otu[ , -c(1,2)]

otu_int <- round(otu)
##2_rarecurve_118_28species_6_6
rarecurve(otu_int,  step = 10,  cex.lab = 1.2,  main = NULL,
          label = FALSE,  lines = TRUE, col = "tomato", lwd = 0.8,
          xlab = "Number of reads per plant species",  ylab = "Number of OTU" )

######## -- Fig S2 -- Plant phylogenetic signals of plant traits, rlc  155 35species ----------------
library(plyr) 
library(ggtree)
library(ape) 
library(picante) 
library(phylosignal) 
library(phylobase) 
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys")  
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
          # 超蔷薇类
          "Rosaceae" =  "#FFA64D", "Fabaceae" = "#5880ba" , 
          "Brassicaceae"= "#9e5c93", "Onagraceae" = "steelblue", 
          # 超菊类
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
  labs(x="",y="Tissue P concentration (mg/g)")+
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


######## -- Table S2 -- field Relationship between plant traits and AM fungal colonization  ---------------
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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF2/field surveys_d3") 
data1 = read.csv("field_data_155_31sp.csv",header=TRUE, row.names=1) 

##
rooted_tree <- read.tree("1_31sp_phylomakers1_abbr.nwk") 
plot(rooted_tree) 

#rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)
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
pgls_model1 <- predict( pgls_model, newdata = data1 )  

##
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

# CI
pgls_ci11 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci)

fig1a <- ggplot( data1 , aes( y = Root.length.colonization , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "Root length colonization (%)" ) +
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "grey22",  fill = "grey88",  formula = y ~ poly(x, 1, raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  #geom_ribbon(data = pgls_ci11,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model1 ), linewidth = 0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text = element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.05, vjust = - 2.8, color = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.082, 0.117 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.05, vjust = - 1.3, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.072, 0.146 ), size = 4, inherit.aes = FALSE )
fig1a


###
lm_model <- lm( P.concentrations ~ Root.length.colonization, data = data1 ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model2 <- predict(pgls_model, newdata = data1 )  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Root.length.colonization), max(dat0$Root.length.colonization), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ Root.length.colonization ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.length.colonization = pred_grid) ) } } }

# CI
pgls_ci22 <- data.frame( Root.length.colonization = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig1b <- ggplot( data1, aes( y = P.concentrations , x = Root.length.colonization  )) +      
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  #geom_ribbon(data = pgls_ci22,  aes(x = Root.length.colonization, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model2 ), linewidth= 0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 4, color = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.069, 0.155 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 4.9, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.036, 0.304 ), size = 4, inherit.aes = FALSE )


fig1b


fig1.1 <- cowplot::plot_grid( fig1a,  fig1b, 
                              ncol = 2, nrow = 1, labels = c('(a)', "(b)", "(c)",'(d)', "(e)"), 
                              label_fontface = "bold", label_x = -0.02, label_y = 1.02 , align = c("hv"))
fig1.1


###
data = read.csv("field_data_109_26sp.csv",header=TRUE, row.names=1) 
colnames(data)

##
rooted_tree <- read.tree("1_26sp_phylomakers1_abbr.nwk") 
plot(rooted_tree)

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data, names.col = "Plant.species", vcv = T, vcv.dim = 2 )
nrow(comp_data$data)
nrow(data)
setdiff(data$Plant.species, rownames(comp_data$data))

##
lm_model <- lm( Rhizophagus ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Rhizophagus ~ Root.diameter , data = comp_data, lambda = "ML" )  #, kappa = "ML", delta = "ML"
summary(pgls_model)
pgls_model11 <- predict(pgls_model, newdata = data )  

##
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( Rhizophagus ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci1 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci1)

fig1d <- ggplot( data , aes(y = Rhizophagus , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = expression("RA of "*italic(Rhizophagus)*" (%)")) +
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci1,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model11 ), linewidth = 0.8 , color = "#CB4641", linetype = "solid" ) + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "grey22", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 1.3, color = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.329, 0.002 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.03, vjust = 2.5, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.329, 0.002 ), size = 4, inherit.aes = FALSE )
fig1d


###
lm_model <- lm( P.concentrations ~ Rhizophagus, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Rhizophagus , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model22 <- predict(pgls_model, newdata = data )  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Rhizophagus), max(dat0$Rhizophagus), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ Rhizophagus ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Rhizophagus = pred_grid) ) } } }

# CI
pgls_ci2 <- data.frame( Rhizophagus = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig1e <- ggplot( data, aes(y = P.concentrations , x = Rhizophagus  ))+      
  scale_x_continuous(name = expression("RA of "*italic(Rhizophagus)*" (%)")) + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "dashed" )+ 
  #geom_ribbon(data = pgls_ci2,  aes(x = Rhizophagus, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  # geom_line(aes(y = pgls_model22 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  #guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.1, colour = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.128, 0.073 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.3, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.099, 0.118 ), size = 4, inherit.aes = FALSE )

fig1e


##
lm_model <- lm( mpd.uw ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  mpd.uw ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model3 <- predict(pgls_model, newdata = data )  

##
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( mpd.uw ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci3 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci)

fig1f <- ggplot( data , aes(y = mpd.uw , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "MPD" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci3 ,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model3 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.13, vjust = - 5.2, colour = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.013, 0.585 ), size = 4 , inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.12, vjust = - 3.4, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.011, 0.614 ), size = 4 , inherit.aes = FALSE )
fig1f


###
lm_model <- lm( P.concentrations ~ mpd.uw, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ mpd.uw , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model4 <- predict(pgls_model, newdata = data )  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$mpd.uw), max(dat0$mpd.uw), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ mpd.uw ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( mpd.uw = pred_grid) ) } } }

# CI
pgls_ci4 <- data.frame( mpd.uw = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )

fig1g <- ggplot( data, aes(y = P.concentrations , x = mpd.uw  ))+      
  scale_x_continuous(name ="MPD") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci4,  aes(x = mpd.uw, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model4 ), linewidth = 0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 4.1, colour = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.177, 0.032 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 5.2, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.231, 0.013 ), size = 4, inherit.aes = FALSE )

fig1g


##
lm_model <- lm( PCoA1 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  PCoA1 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model5 <- predict(pgls_model, newdata = data )  

##
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( PCoA1 ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci5 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci5)

fig1i <- ggplot( data , aes(y = PCoA1 , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "PCoA axis1" ) +
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci5,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model5 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust = - 6.8, colour = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.004, 0.754 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.09, vjust = - 5, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.004, 0.754 ), size = 4, inherit.aes = FALSE )
fig1i


###
lm_model <- lm( P.concentrations ~ PCoA1, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ PCoA1 , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model6 <- predict(pgls_model, newdata = data )  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$PCoA1), max(dat0$PCoA1), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ PCoA1 ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( PCoA1 = pred_grid) ) } } }

# CI
pgls_ci6 <- data.frame( PCoA1 = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig1j <- ggplot( data, aes(y = P.concentrations ,  x = PCoA1  )) +      
  scale_x_continuous(name ="PCoA axis1") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci6,  aes(x = PCoA1, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model6 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust =  2.2, colour = "grey22", 
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.166, 0.039 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.09, vjust =  3.4, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.290, 0.005 ), size = 4, inherit.aes = FALSE )

fig1j


##
lm_model <- lm( PCoA2 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  PCoA2 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model7 <- predict(pgls_model, newdata = data )  

##
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( PCoA2 ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci7 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci7)

fig1k <- ggplot( data , aes(y = PCoA2 , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "PCoA axis2" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "dashed"  )+ 
  #geom_ribbon(data = pgls_ci7,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model7 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.2, vjust =  -3 , colour = "grey22", 
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.144, 0.056 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.18, vjust = -1.5, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.070, 0.192 ), size = 4, inherit.aes = FALSE )
fig1k


###
lm_model <- lm( P.concentrations ~ PCoA2, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ PCoA2 , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model8 <- predict(pgls_model, newdata = data )  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$PCoA2), max(dat0$PCoA2), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ PCoA2 ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( PCoA2 = pred_grid) ) } } }

# CI
pgls_ci8 <- data.frame( PCoA2 = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig1l <- ggplot( data, aes(y = P.concentrations , x = PCoA2  )) +      
  scale_x_continuous(name ="PCoA axis2") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci8,  aes(x = PCoA2, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model8 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text = element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1 , vjust =  - 1.8 , colour = "grey22", 
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.274, 0.006 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.09 , vjust =  - 0.4 , color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.274, 0.006 ), size = 4, inherit.aes = FALSE )


fig1l


##
lm_model <- lm( Module.2 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Module.2 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model111 <- predict(pgls_model, newdata = data )  

##
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( Module.2 ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci9 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci1)

fig2m <- ggplot( data , aes( y = Module.2 , x = Root.diameter )) +      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "RA of Module #2 (%)" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci9,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model111 ), linewidth = 0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 4.6, colour = "grey22", 
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.248, 0.010 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.03, vjust = 5.65, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.248, 0.010 ), size = 4, inherit.aes = FALSE )
fig2m


###
lm_model <- lm( P.concentrations ~ Module.2, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Module.2 , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model222 <- predict(pgls_model, newdata = data )  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Module.2), max(dat0$Module.2), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ Module.2 ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Module.2 = pred_grid) ) } } }

# CI
pgls_ci10 <- data.frame( Module.2 = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig2n <- ggplot( data, aes(y = P.concentrations , x = Module.2  ))+      
  scale_x_continuous(name ="RA of Module #2 (%)") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci10,  aes(x = Module.2, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model222 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1, colour = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.300, 0.004 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.3, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.280, 0.005 ), size = 4, inherit.aes = FALSE )

fig2n


fig1.1 <- cowplot::plot_grid( fig1a, fig1b, fig1d, fig1e,  fig1f, fig1g,  fig1i,  fig1j, fig1k,  fig1l, fig2m, fig2n,
                              ncol = 3, nrow = 4, labels = c( "(a)", "(b)",'(c)', "(d)", "(e)","(f)",  "(g)",  "(h)", 
                                                              "(i)", "(j)","(k)", "(l)", "(m)", "(n)", "(o)", "(p)",  
                                                              "(q)","(r)","(s)", "(t)"), 
                              label_fontface = "bold", label_x = -0.02, label_y = 1.0 , align = c("hv"))
fig1.1


eoffice::topptx( fig1.1 , width = 10 ,  height = 14 , "1_0617_109_26sp_fig1_10*14.pptx" ) 


fig1.3 <- cowplot::plot_grid( fig1h,  fig1k,  fig1l,
                              ncol = 3, nrow = 1, labels = c( "", "(j)", "(k)" ), 
                              rel_widths  = c(1.015, 0.5, 0.5), 
                              label_fontface = "bold", label_x = -0.02, label_y = 1.0 , align = c(""))
fig1.3



#1_0529_fig.1_14*11
fig1 <- cowplot::plot_grid(fig1.1,  fig1.3, 
                           ncol = 1, nrow = 2, 
                           rel_heights = c(2, 1),  # 上面是下面的一半
                           label_fontface = "bold", 
                           label_x = -0.02, 
                           label_y = 1.0, 
                           align = c("v"))
fig1


eoffice::topptx( fig1 , width = 11.52 ,  height = 12.4 , "1_118_28sp_fig1_11.52*12.4.pptx" ) 


######## -- Table S3 -- Plant phylogenetic signals of AM fungal variables 28species ----------------
library(plyr) 
library(ggtree)
library(ape) 
library(picante) 
library(phylosignal) 
library(phylobase) 
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys")  
data_mean = read.csv("field_data_118_28sp.csv",header=TRUE, row.names=1) 

data_mean1 = data_mean
rownames(data_mean1) <- data_mean[, 1]
data_mean1 <- data_mean1[ , -c(1,2,3)]

##
rooted_tree <- read.tree("1_28sp_phylomakers1_abbr.nwk") 
plot(rooted_tree)

tree_p <- as(rooted_tree, "phylo4")
class(tree_p)
p4d <- phylo4d(tree_p, data_mean1)
phy.sig <- phyloSignal(p4d, reps = 9999)
phy.sig ## 28sp 118

write.csv(phy.sig, '1_0624_118_28sp_physig.csv')  



######## -- Fig S4 a -- AM fungal genus ---------------
library(ggalluvial)
library(ggplot2)
library(tidyr)
library(dplyr)
rm(list = ls()) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys")
data = read.csv("field_data_118_28sp.csv",header=TRUE, row.names=1) 
colnames(data)

genus = data[, c(15:23)]
head(genus)
rownames(genus) <- data[, 1]
genus$Sample <- rownames(genus)

genus_long <- genus  %>%  pivot_longer(cols = -Sample, names_to = "Genus", values_to = "Abundance")

genus_order <- c("Dominikia", "Rhizophagus", "Acaulospora", "Polonospora", 
                 "Funneliformis", "Entrophospora",   "Kamienskia", "Microdominikia", "Others")

genus_long$Genus <- factor(genus_long$Genus, levels = genus_order)

genus_labels <- setNames(
  c(expression(italic("Dominikia")), 
    expression(italic("Rhizophagus")), 
    expression(italic("Acaulospora")), 
    expression(italic("Polonospora")), 
    expression(italic("Funneliformis")), 
    expression(italic("Entrophospora")), 
    expression(italic("Microdominikia")), 
    expression(italic("Kamienskia")), 
    "Others"),    c("Dominikia", "Rhizophagus", "Acaulospora", "Polonospora", 
                    "Funneliformis", "Entrophospora",   "Kamienskia", "Microdominikia", "Others") )


p_sankey1 <- ggplot(genus_long, aes(x = Sample, y = Abundance, alluvium = Genus,  stratum = Genus, fill = Genus) ) +
  geom_alluvium(aes(fill = Genus),   alpha = 0.3,  width = 0.5, show.legend = F ) +
  geom_stratum(width = 0.5, alpha = 0.7, color = "black",show.legend = T ) +
  scale_fill_manual(values = c("Dominikia" = "#5880ba",  "Rhizophagus" = "#E9746A", "Acaulospora" = "#7B8F5E",
                               "Polonospora" = "#B07AA8","Funneliformis" = "#497d75", "Entrophospora" = "#FFA64D",
                               "Kamienskia" = "#A65928",   "Microdominikia" = "#8FB9D1","Others" = "#999999"), 
                    breaks = genus_order, labels = genus_labels) +  
  labs(x = "", y = "Relative abundance (%)",  title = "AM Fungal Genus Composition Across species", fill = "Genus") +
  theme_minimal() + theme(axis.text.x = element_text(size=11, face = "italic", color = "black"), axis.text.y = element_text(size=10, angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0, size=12),  
                          axis.title.x = element_blank(), axis.title.y = element_text(size=12.5), plot.title = element_blank(),
                          plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.55) )

p_sankey1


figs5a <- cowplot::plot_grid( p_sankey1 ,
                              ncol = 1, nrow = 1, labels = c( "(a)",'(d)', "(e)"), 
                              label_fontface = "bold", label_x = -0.005, label_y = 1.01 , align = c("hv"))
figs5a

#eoffice::topptx( figs5.2 , width = 14 ,  height= 8.5 , "1_figs4_14*8.5.pptx" )


######## -- Fig S4 b -- pcoa ------------------------------------
library(vegan)
library(ggplot2)
library(ggrepel)
library(ape)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys")
data = read.csv("field_data_118_28sp.csv",header=TRUE, row.names=1) 
otu_a <- read.csv("field_data_386otu_ra_118_28sp.csv",header=TRUE, row.names=1)
#otu = otu_a[, -c(1,2,3)] 

vare.dist <- vegdist(otu_a, method="bray")

pcoa_res <- pcoa(vare.dist)

pcoa_axis <- as.data.frame(pcoa_res$vectors[, 1:2])


colnames(pcoa_axis) <- c("PCoA1", "PCoA2")

eig_percent <- round(pcoa_res$values$Relative_eig[1:2] * 100, 2)
eig_percent

#write.csv( pcoa_axis , '1_rare_noh_118_28sp_pcoa12.csv')  

set.seed(123)
# envfit
fit <- envfit( pcoa_res$vectors , data$Root.diameter, permutations = 9999)
fit
fit <- envfit( pcoa_res$vectors , data$P.concentrations, permutations = 9999)
fit

adonis2(vare.dist ~ data$Root.diameter, permutations = 9999) 
adonis2(vare.dist ~ data$P.concentrations, permutations = 9999) 

fit <- envfit( pcoa_res$vectors, data[, c("Root.diameter", "P.concentrations")], permutations = 9999 ) 

# 
fit_df <- data.frame( Axis1 = fit$vectors$arrows[, 1], Axis2 = fit$vectors$arrows[, 2],
                      r2 = fit$vectors$r,
                      P = fit$vectors$pvals )

fit_df$Variable <- rownames(fit$vectors$arrows)

# 调整箭头长度（可根据图形再微调）
arrow_mul <- 0.6
fit_df$PCoA1 <- fit_df$Axis.1 * arrow_mul
fit_df$PCoA2 <- fit_df$Axis.2 * arrow_mul
fit_df$Variable <- c("Root diameter", "Tissue P concentration")

# 准备数据
pcoa_axis$Sample <- rownames(otu_a)
eig_percent_formatted <- sprintf("%.2f", eig_percent)

# 基础 PCoA 图
pcoa <- ggplot(pcoa_axis, aes(x = PCoA1, y = PCoA2, )) + 
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_text_repel(aes(label = Sample), size = 4.5, max.overlaps = 20) +
  #geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  #geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = paste0("PCoA axis1 (", eig_percent_formatted[1], "%)"), y = paste0("PCoA axis2 (", eig_percent_formatted[2], "%)"), title = "AM fungal communities" ) +
  labs(x = paste0("PCoA axis1 (", eig_percent_formatted[1], "%)"), y = paste0("PCoA axis2 (", eig_percent_formatted[2], "%)"), title = "AM fungal communities" ) +
  theme_bw() + # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  theme(panel.grid = element_blank(), axis.text.y = element_text(angle = 90, hjust = 0.5), 
        plot.title = element_blank(), #element_text(hjust = 0, size = 16), 
        legend.text = element_text(size=10), legend.title = element_text(hjust=0.5, size = 14), legend.position = c(0.8,0.8),
        axis.title = element_text(size = 12.5), axis.text = element_text(size = 10.5)) 
# geom_segment( data = fit_df, aes(x = 0, y = 0, xend = PCoA1, yend = PCoA2), arrow = arrow(length = unit(0.25, "cm")), linewidth = 0.8, color = "black", inherit.aes = FALSE ) +
#  geom_text( data = fit_df,aes(x = PCoA1 * 0.9, y = PCoA2 *1.2, label = Variable ),
#             size = 3.1, color = "black", nudge_x = 0.02,  nudge_y = 0.02,inherit.aes = FALSE ) +
#  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), 
#             size = 3, inherit.aes = FALSE,  parse = TRUE, hjust = -1.85, vjust = 23.3, label = sprintf("italic(R^2)~'='~%.3f*', '*italic(P)~'='~%.3f", 0.231, 0.036) ) +
#  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y),
#             size = 3, inherit.aes = FALSE , parse = TRUE, hjust = -1.3, vjust = 25.8, label = sprintf("italic(R^2)~'='~%.3f*', '*italic(P)~'<'~%.3f", 0.453, 0.001) )

#  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y),   parse = TRUE, hjust = -2.6, vjust = 2.3, label = "'Root diameter:'",  size = 4, inherit.aes = FALSE ) + 
# geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), size = 4, inherit.aes = FALSE,  parse = TRUE, hjust = -2.2, vjust = 2.9, label = sprintf("italic(R^2)~'='~%.3f*', '*italic(P)~'='~%.3f", 0.231, 0.036) ) +
#geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y),  parse = TRUE, hjust = -1.64, vjust = 6.7,  label = "'Tissue P concentration:'",  size = 4, inherit.aes = FALSE ) + 
# geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), size = 4, inherit.aes = FALSE , parse = TRUE, hjust = -2.2, vjust = 5.8, label = sprintf("italic(R^2)~'='~%.3f*', '*italic(P)~'<'~%.3f", 0.453, 0.001) )

pcoa


######## -- Fig S4 c --  genus ----------------------------------
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
data1 = read.csv("1_0624_network_module5.csv",header=TRUE, row.names=1) 
head(data1)

Module = c("Module #1", "Module #2", "Module #3", "Module #4", "Module #5")
data1$Module <- Module 
Module <- rownames(data1) 


data_long <- data1 %>%
  pivot_longer(cols = -Module, names_to = "Genus", values_to = "Abundance")

# 
genus_order <- c("Dominikia", "Rhizophagus", "Acaulospora", "Polonospora", 
                 "Funneliformis", "Entrophospora", "Kamienskia",  "Microdominikia", "Others")

# 
data_long$Genus <- factor(data_long$Genus, levels = genus_order)

genus_labels <- setNames(
  c(expression(italic("Dominikia")), 
    expression(italic("Rhizophagus")), 
    expression(italic("Acaulospora")), 
    expression(italic("Polonospora")), 
    expression(italic("Funneliformis")), 
    expression(italic("Entrophospora")), 
    expression(italic("Kamienskia")), 
    expression(italic("Microdominikia")), 
    "Others"),  # Others 不加 italic
  c("Dominikia", "Rhizophagus", "Acaulospora", "Polonospora", 
    "Funneliformis", "Entrophospora", "Kamienskia", "Microdominikia", "Others")
)


p_sankey2 <- ggplot(data_long, aes(x = Module, y = Abundance, alluvium = Genus,  stratum = Genus, fill = Genus) ) +
  geom_alluvium(aes(fill = Genus),   alpha = 0.3,  width = 0.3, show.legend = F ) +
  geom_stratum(width = 0.3, alpha = 0.7, show.legend = T ) +
  scale_fill_manual(values = c("Dominikia" = "#5880ba",  "Rhizophagus" = "#E9746A", "Acaulospora" = "#7B8F5E",
                               "Polonospora" = "#B07AA8","Funneliformis" = "#497d75", "Entrophospora" = "#FFA64D",
                               "Kamienskia" = "#A65928",   "Microdominikia" = "#8FB9D1","Others" = "#999999"), 
                    breaks = genus_order,labels = genus_labels)  +  scale_x_discrete(expand = c(0.08, 0.08)) + 
  labs(x = "", y = "Relative abundance (%)",  title = "AM Fungal Genus Composition Across module", fill = "Genus") +
  theme_minimal() + theme(axis.text.x = element_text(size=11, face = "plain", vjust = -1,  color = "black"), axis.text.y = element_text(size=10, angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0, size=12),  legend.position = "none",
                          axis.title.x = element_blank(), axis.title.y = element_text(size=12.5), plot.title = element_blank(),
                          plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.55) )

p_sankey2


###
figs5.1 <- cowplot::plot_grid( pcoa, p_sankey2, 
                               ncol = 2, nrow = 1, labels = c( "(b)",'(c)', "(d)"), 
                               label_fontface = "bold", label_x = -0.02, label_y = 1.03 , align = c("hv"))
figs5.1

#1_0529_figs5_14_8.5
figs5.2 <- cowplot::plot_grid( figs5a , figs5.1 ,
                               ncol = 1, nrow = 2,  
                               label_fontface = "bold", label_x = -0.01, label_y = 1.01 , align = c("hv"))
figs5.2


eoffice::topptx( figs5.2 , width = 14 ,  height= 8.5 , "1_0624_figs4_14*8.5.pptx" )


######## -- Fig S5 -- Distribution and sharing patterns of AM fungal OTUs ---------------
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
data_all = read.csv("field_data_386otu_rare_118_28sp.csv", header=TRUE, row.names=1) 
colnames(data_all)
rownames(data_all) <- data_all[, 2]
data_all <- data_all[ , -c(1,2,3)]
data = t( data_all )

data <- data[rowSums(data) > 0, ]

otu_counts <- colSums(data > 0) 

df_plant <- data.frame( 
  Plant = names(otu_counts), 
  OTU_Count = as.numeric(otu_counts) ) 

gg_plant_otu <- ggplot(df_plant, aes(x = reorder(Plant, -OTU_Count), y = OTU_Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = OTU_Count), vjust = -0.5, size = 3.5) +
  labs(x = "Plant species", y = "Number of OTUs", title = "Number of OTUs per plant species") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 11),
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "italic"),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    axis.title = element_text(size = 12),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 11),
    legend.title = element_text(hjust = 0.8)
  )

gg_plant_otu

plants_per_otu <- rowSums(data > 0)

otu_dist_table <- table(plants_per_otu)

df_otu_dist <- data.frame(
  Num_Plants = as.numeric(names(otu_dist_table)),
  Num_OTUs = as.numeric(otu_dist_table) )

gg_otu_dist <- ggplot(df_otu_dist, aes(x = Num_Plants, y = Num_OTUs)) +
  geom_bar(stat = "identity", fill = "tomato") +
  geom_text(aes(label = Num_OTUs), vjust = -0.5, size = 3.5) +
  scale_x_continuous(breaks = seq(min(df_otu_dist$Num_Plants), max(df_otu_dist$Num_Plants), by = 3)) +
  labs(x = "Number of plant species", y = "Number of shared OTUs",
       title = "Frequency distribution of OTUs shared across plant species") +
  theme_classic() +
  theme(
    axis.text = element_text(size = 11),
    axis.text.y = element_text(angle = 90, hjust = 0.5),
    axis.title = element_text(size = 12),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.text = element_text(size = 11),
    legend.title = element_text(hjust = 0.8)
  )
gg_otu_dist


# Extract OTUs that appear in only one plant
otu_single <- names(plants_per_otu[plants_per_otu == 1])

# Obtain the presence matrix of these OTUs from the raw data
single_otu_data <- data[otu_single, , drop = FALSE]  

# Identify the single plant name (column name) corresponding to each OTU
plant_for_otu <- apply(single_otu_data, 1, function(row) {
  colnames(single_otu_data)[which(row > 0)]
})


result_df <- data.frame( OTU = otu_single, Plant = plant_for_otu, stringsAsFactors = FALSE )

# Add genus annotations
library(dplyr)
library(ggalluvial)
tax = read.csv("field_data_386otu_taxon.csv", header= T, row.names=1 ) 
tax$OTU_ID = rownames(tax) 

tax_sub <- tax %>% filter(OTU_ID %in% otu_single) %>%
  select(OTU_ID, Genus)

result_df <- merge(result_df, tax_sub, by.x = "OTU", by.y = "OTU_ID", all.x = TRUE)

head(result_df)

genus_plant_table <- as.data.frame(table(result_df$Genus, result_df$Plant))
colnames(genus_plant_table) <- c("Genus", "Plant", "Count")

genus_plant_counts <- genus_plant_table[genus_plant_table$Count > 0, ]

genus_plant_counts <- genus_plant_counts[order(-genus_plant_counts$Count), ]

sankey_data <- genus_plant_counts
sankey_data$Flow <- 1:nrow(sankey_data)  

sankey_data$Flow <- seq_len(nrow(sankey_data))

head(sankey_data)
length(unique(sankey_data$Genus)) # 17 genus
unique(sankey_data$Genus)
length(unique(sankey_data$Plant)) # 26species

library(ggplot2)
library(ggalluvial)
library(dplyr)

genus_order <- c("Dominikia", "Rhizophagus", "Acaulospora", 
                 "Funneliformis", "Entrophospora","Microdominikia", 
                 "Kamienskia", "Newgenus1",  "Diversispora", 
                 "Sacculospora",  "Newgenus3","Glomus",    
                 "Paraglomus",  "Archaeospora", "Scutellospora",  
                 "Newgenus2", "Pacispora"  )

# 
sankey_data$Genus <- factor(sankey_data$Genus, levels = genus_order)

# 
genus_labels <- setNames(
  c(expression(italic("Dominikia")), 
    expression(italic("Rhizophagus")), 
    expression(italic("Acaulospora")), 
    expression(italic("Funneliformis")), 
    expression(italic("Entrophospora")), 
    expression(italic("Microdominikia")), 
    
    expression(italic("Kamienskia")), 
    expression(italic("Newgenus1")), 
    expression(italic("Diversispora")), 
    expression(italic("Sacculospora")), 
    expression(italic("Newgenus3")), 
    expression(italic("Glomus")), 
    
    expression(italic("Paraglomus")), 
    expression(italic("Archaeospora")), 
    expression(italic("Scutellospora")), 
    expression(italic("Newgenus2")), 
    expression(italic("Pacispora")) ),
  
  c("Dominikia", "Rhizophagus", "Acaulospora", 
    "Funneliformis", "Entrophospora","Microdominikia", 
    "Kamienskia", "Newgenus1",  "Diversispora", 
    "Sacculospora",  "Newgenus3","Glomus",    
    "Paraglomus",  "Archaeospora", "Scutellospora",  
    "Newgenus2", "Pacispora"  )  )

# 桑基图
p_sankey <- ggplot( sankey_data, aes(axis1 = Genus,  axis2 = Plant,   y = Count) ) +
  geom_alluvium( aes(fill = Genus), width = 0.25, alpha = 0.5 ) +
  geom_stratum( aes(fill = after_stat(stratum)), width = 0.25, color = "grey99", alpha = 0.8 ) +
  geom_text( stat = "stratum",aes(label = after_stat(stratum)), size = 2) +
  scale_x_discrete( limits = c("AM fungi", "Plant"), expand = c(.08, .08) ) +
  scale_fill_manual(values = c( "Dominikia" = "#5880ba", "Rhizophagus" = "#E9746A","Acaulospora" = "steelblue",
                                "Funneliformis" = "#FFA64D", "Entrophospora" = "#497d75", "Microdominikia" = "#8FB9D1", "Kamienskia" = "#A65928",
                                "Newgenus1" = "#4A6B7C",  "Diversispora" = "#B07AA8",  "Sacculospora" = "#999999",  "Newgenus3" = "#AF4D2A", "Glomus" = "#BF8F98", "Paraglomus" = "#8ecad0",  
                                "Archaeospora" = "#AA67B3", "Scutellospora" = "#7FC1B5",  "Newgenus2" = "#11b7Ca", "Pacispora" = "#bf887f" ),
                    breaks = genus_order,  labels = genus_labels,  na.value = NA ) +
  labs( x = NULL, y = "Interaction frequency", fill = "Genus", title = "Distribution of 93 plant-specific AM fungal OTUs" ) +
  theme_minimal() + theme( axis.text.x = element_text(size = 10, face = "bold"), axis.text.y = element_blank(), 
                           axis.title = element_text(size = 12),  axis.ticks = element_blank(), legend.position = "none",
                           legend.text = element_text(size = 10), legend.title = element_text(size = 11),
                           panel.grid = element_blank(), panel.border = element_rect( fill = NA,  color = "grey80", linewidth = 0.6 ),
                           plot.background = element_blank() )  + geom_text( stat = "alluvium", data = subset(sankey_data, Count >= 0),
                                                                             aes(label = paste0( Count)), size = 1 , nudge_x = 0.1)

p_sankey


#2_figs4_8*14
figs4 <- cowplot::plot_grid( gg_plant_otu,  gg_otu_dist , p_sankey, ncol = 1, nrow = 3, rel_heights = c(1, 1, 2.5),
                             labels = c('(a)', "(b)", "(c)",'(d)', "(e)", "(f)", "(g)", "(h)", 
                                        "(i)", "(j)","(k)", "(l)", "(m)","(n)","(o)", "(p)", 
                                        "(q)","(r)","(s)", "(t)"), 
                             label_fontface = "bold", label_x = 0 , label_y = 1.01 , align = c("h"))

figs4


eoffice::topptx( figs4 , width = 12 ,  height = 16 , "1_118_28sp_figs4_8*20.pptx" ) 


######## Nestedness and modularity of plant-AM fungal bipartite network ---------------
library(bipartite) 
library(vegan)
library(dplyr) 
rm(list = ls()) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
# nestedness of plant-AM fungal bipartite network
otu_a = read.csv("field_data_386otu_ra_118_28sp.csv",header=TRUE, row.names=1) 

binary_matrix <- ifelse(otu_a > 0, 1, 0)

set.seed(1234)
null_model <- nullmodel( binary_matrix , N = 1000, method = "r2dtable")

nested_NODF <- function(x) nested(x, method = "NODF")

result <- oecosimu( binary_matrix, nested_NODF, method = "r2dtable", nsimul = 1000 )
result

result$oecosimu$statistic
#      statistic  SES   mean   2.5%   50%   97.5%   Pr(sim.)    
# NODF    31.659 20.254 23.972 23.227 23.988 24.684 0.000999 ***
result$oecosimu$means
# 23.97194
sd(result$oecosimu$simulated)
# 0.3795319


# modularity of plant-AM fungal bipartite network
set.seed(1234)
otu_a = read.csv("field_data_386otu_ra_118_28sp.csv",header=TRUE, row.names=1) 

binary_matrix <- ifelse(otu_a > 0, 1, 0)

mod_result <- computeModules(binary_matrix)
plotModuleWeb(mod_result)

mod_obs <- mod_result
mod_obs@likelihood
Q_obs <- mod_obs@likelihood
Q_obs
# [1] 0.2147847

n_sim <- 1000
Q_null <- numeric(n_sim)


for (i in 1:n_sim) {
  null_mat <- nullmodel(binary_matrix, N = 1, method = "r2dtable")[[1]]
  mod_null <- try(computeModules(null_mat), silent = TRUE)
  if (!inherits(mod_null, "try-error")) {
    Q_null[i] <- mod_null@likelihood
  } else {
    Q_null[i] <- NA
  }
}

Q_null <- Q_null[!is.na(Q_null)]

p_value <- (sum(Q_null >= Q_obs) + 1) / (length(Q_null) + 1)
p_value

mean(Q_null)
# [1] 0.237242
sd(Q_null)
# [1] 0.003127572

SES <- (Q_obs - mean(Q_null))/sd(Q_null)
SES
# [1] -7.180437
p_two <- (sum(abs(Q_null - mean(Q_null)) >= abs(Q_obs - mean(Q_null))) + 1)/(n_sim + 1)
p_two
# [1] 0.000999001
p_upper <- (sum(Q_null >= Q_obs) + 1 )/(n_sim+1)
p_upper
# [1] 1



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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
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
head(pgls_ci)

fig1a <- ggplot( data1 , aes(y = Root.length.colonization , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "Root length colonization (%)" ) +
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22",  fill = "grey88",  formula = y ~ poly(x, 1, raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  #geom_ribbon(data = pgls_ci11,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model1 ), linewidth = 0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text = element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.05, vjust = - 3, color = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.132, 0.032 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.05, vjust = - 1.6, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.073, 0.116 ), size = 4, inherit.aes = FALSE )
fig1a

##
lm_model <- lm( P.concentrations ~ Root.length.colonization, data = data1 ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model2 <- predict(pgls_model, newdata = data1 )  

phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Root.length.colonization), max(dat0$Root.length.colonization), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ Root.length.colonization ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.length.colonization = pred_grid) ) } } }

# CI
pgls_ci22 <- data.frame( Root.length.colonization = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )

fig1b <- ggplot( data1, aes(y = P.concentrations , x = Root.length.colonization  ))+      
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  #geom_ribbon(data = pgls_ci22,  aes(x = Root.length.colonization, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model2 ), linewidth= 0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 3.6, color = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.119, 0.042), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 4.7, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.057, 0.167), size = 4, inherit.aes = FALSE )


fig1b


######## -- Fig 1 -- Relationship between plant traits and AM fungal variables 
library(ggplot2)
library(gridExtra)
library(reshape)
library(ggpmisc)
library(car)
library(caper)
library(doBy)
library(dplyr) 
library(ape)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
data = read.csv("field_data_118_28sp.csv",header=TRUE, row.names=1) 
colnames(data)

b<- data[, c(4:35)]
str(b)

newrow <- list()
for (i in 1:length(b)) {  m <- cor.test(b[,i],b$Root.diameter , method="pearson")
newrow[[i]] <- cbind(m[4],m[3]) }
result <- do.call(rbind,newrow)
colnames(result) <- c("pearson r","p.value")
rownames(result) <- colnames(b)
sig.result <- result[result[,2]<=0.1,]
sig.result

##
rooted_tree <- read.tree("1_28sp_phylomakers1_abbr.nwk") 
plot(rooted_tree)

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data, names.col = "Plant.species", vcv = T, vcv.dim = 2 )
nrow(comp_data$data)
nrow(data)
setdiff(data$Plant.species, rownames(comp_data$data))

##
lm_model <- lm( Rhizophagus ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Rhizophagus ~ Root.diameter , data = comp_data, lambda = "ML" )  
summary(pgls_model)
pgls_model11 <- predict(pgls_model, newdata = data )  

coef(lm_model)
coef(pgls_model)
all.equal( coef(lm_model), coef(pgls_model) )

##
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( Rhizophagus ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci1 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci1)

fig1d <- ggplot( data , aes(y = Rhizophagus , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = expression("RA of "*italic(Rhizophagus)*" (%)")) +
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci1,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model11 ), linewidth = 0.8 , color = "#CB4641", linetype = "solid" ) + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "grey22", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 1.3, color = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.375, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.03, vjust = 2.5, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.375, 0.001 ), size = 4, inherit.aes = FALSE )
fig1d



###
lm_model <- lm( P.concentrations ~ Rhizophagus, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Rhizophagus , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model22 <- predict(pgls_model, newdata = data )  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Rhizophagus), max(dat0$Rhizophagus), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ Rhizophagus ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Rhizophagus = pred_grid) ) } } }

# CI
pgls_ci2 <- data.frame( Rhizophagus = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig1e <- ggplot( data, aes(y = P.concentrations , x = Rhizophagus  ))+      
  scale_x_continuous(name = expression("RA of "*italic(Rhizophagus)*" (%)")) + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci2,  aes(x = Rhizophagus, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model22 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  #guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.1, colour = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.169, 0.030), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.3, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.117, 0.074), size = 4, inherit.aes = FALSE )

fig1e


##
lm_model <- lm( OTU.richness ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  OTU.richness ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model3 <- predict(pgls_model, newdata = data )  

##
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( OTU.richness ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci3 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci)

fig1f <- ggplot( data , aes(y = OTU.richness , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "OTU richness" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci3 ,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model3 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = - 4.5, colour = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.013, 0.560 ), size = 4 , inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = - 3, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.013, 0.560 ), size = 4 , inherit.aes = FALSE )
fig1f


###
lm_model <- lm( P.concentrations ~ OTU.richness, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ OTU.richness , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model4 <- predict(pgls_model, newdata = data )  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$OTU.richness), max(dat0$OTU.richness), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ OTU.richness ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( OTU.richness = pred_grid) ) } } }

# CI
pgls_ci4 <- data.frame( OTU.richness = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )

fig1g <- ggplot( data, aes(y = P.concentrations , x = OTU.richness  ))+      
  scale_x_continuous(name ="OTU richness") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid" )+ 
  #geom_ribbon(data = pgls_ci4,  aes(x = OTU.richness, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model4 ), linewidth = 0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 4, colour = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'<'~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.001, 0.972 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 5, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.002, 0.843 ), size = 4, inherit.aes = FALSE )

fig1g


##
lm_model <- lm( Family.richness ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Family.richness ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model5 <- predict(pgls_model, newdata = data )  

##
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( Family.richness ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci5 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci5)

fig1i <- ggplot( data , aes(y = Family.richness , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "Family richness" ) +
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci5,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model5 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
 # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = 1.1 , colour = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.005, 0.713 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = 2.2 , color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.005, 0.713 ), size = 4, inherit.aes = FALSE )
fig1i


###
lm_model <- lm( P.concentrations ~ Family.richness, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Family.richness , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model6 <- predict(pgls_model, newdata = data )  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Family.richness), max(dat0$Family.richness), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ Family.richness ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Family.richness = pred_grid) ) } } }

# CI
pgls_ci6 <- data.frame( Family.richness = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig1j <- ggplot( data, aes(y = P.concentrations ,  x = Family.richness  )) +      
  scale_x_continuous(name ="Family richness") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci6,  aes(x = Family.richness, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model6 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust =  3.7, colour = "grey22", 
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.161, 0.034 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.03, vjust =  4.7, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.154, 0.039 ), size = 4, inherit.aes = FALSE )

fig1j


##
lm_model <- lm( PCoA2 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  PCoA2 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model7 <- predict(pgls_model, newdata = data )  

##
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( PCoA2 ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci7 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci7)

fig1k <- ggplot( data , aes( y = PCoA2 , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "PCoA axis2 of AM fungal community" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci7,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model7 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust =  2.2 , colour = "grey22", 
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.224, 0.011 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust =  3.3, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.130, 0.059 ), size = 4, inherit.aes = FALSE )
fig1k


###
lm_model <- lm( P.concentrations ~ PCoA2, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ PCoA2 , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model8 <- predict(pgls_model, newdata = data )  

coef(lm_model)
coef(pgls_model)
all.equal( coef(lm_model), coef(pgls_model) )

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$PCoA2), max(dat0$PCoA2), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ PCoA2 ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( PCoA2 = pred_grid) ) } } }

# CI
pgls_ci8 <- data.frame( PCoA2 = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig1l <- ggplot( data, aes(y = P.concentrations , x = PCoA2  )) +      
  scale_x_continuous(name ="PCoA axis2 of AM fungal community") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22",  fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci8,  aes(x = PCoA2, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model8 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text = element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  3.8, colour = "grey22", 
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.315, 0.002 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust =  4.8, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.315, 0.002 ), size = 4, inherit.aes = FALSE )


fig1l


##
lm_model <- lm( Module.2 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Module.2 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model111 <- predict(pgls_model, newdata = data )  

coef(lm_model)
coef(pgls_model)
all.equal( coef(lm_model), coef(pgls_model) )


##
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( Module.2 ~ Root.diameter ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.diameter = pred_grid) ) } } }

# CI
pgls_ci9 <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )
head(pgls_ci1)

fig2m <- ggplot( data , aes( y = Module.2 , x = Root.diameter ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "RA of Module #2 (%)" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth = 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci9,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model111 ), linewidth = 0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.06, vjust = 0.8, colour = "grey22", 
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.285, 0.003 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 2, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.285, 0.003 ), size = 4, inherit.aes = FALSE )
fig2m


###
lm_model <- lm( P.concentrations ~ Module.2, data = data ) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Module.2 , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model222 <- predict(pgls_model, newdata = data )  

##
phy  <- comp_data$phy
dat0 <- comp_data$data
dat0$Plant.species <- rownames(dat0)
colnames(dat0)
set.seed(123)
nboot <- 1000
pred_grid <- seq(min(dat0$Module.2), max(dat0$Module.2), length.out = 100)
boot_mat <- matrix(NA, nrow = nboot, ncol = length(pred_grid))
for (i in 1:nboot) { idx <- sample(1:nrow(dat0), size = floor(0.8 * nrow(dat0)),  replace = FALSE)
boot_data <- dat0[idx, ]
comp_boot <- try( comparative.data( phy = phy,  data = boot_data, names.col = "Plant.species", vcv = TRUE, warn.dropped = FALSE  ), silent = TRUE)
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ Module.2 ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Module.2 = pred_grid) ) } } }

# CI
pgls_ci10 <- data.frame( Module.2 = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig2n <- ggplot( data, aes(y = P.concentrations , x = Module.2  )) +      
  scale_x_continuous(name ="RA of Module #2 (%)") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  geom_ribbon(data = pgls_ci10,  aes(x = Module.2, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model222 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1, colour = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.335, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.1, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.300, 0.003 ), size = 4, inherit.aes = FALSE )

fig2n


#1_0529_fig.1_14*11
fig1.1 <- cowplot::plot_grid( fig1a,  fig1f, fig1k, 
                              fig1b, fig1g, fig1l,
                              fig1d, fig1i, fig2m,
                              fig1e, fig1j,fig2n,
                              
                              ncol = 3, nrow = 4, labels = c( "(a)", "(e)",'(i)',
                                                              "(b)", "(f)","(j)",  
                                                              "(c)", "(g)", "(k)",
                                                              "(d)", "(h)", "(l)" ), 
                              label_fontface = "bold", label_x = -0.02, label_y = 1.0 , align = c("hv"))
fig1.1

eoffice::topptx( fig1.1 , width = 10 ,  height = 14.4 , "1_0624_118_28sp_fig1_10*14.4.pptx" ) 

######## -- r p ols  pgls   ---------------
library(ggplot2)
library(gridExtra)
library(reshape)
library(ggpmisc)
library(car)
library(caper)
library(doBy)
library(dplyr) 
library(ape)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/field surveys") 
data = read.csv("field_data_118_28sp.csv",header=TRUE, row.names=1) 
colnames(data)

b<- data[, c(3:99)]
str(b)

newrow <- list()
for (i in 1:length(b)) {  m <- cor.test(b[,i],b$P.concentrations , method="pearson")
newrow[[i]] <- cbind(m[4],m[3]) }
result <- do.call(rbind,newrow)
colnames(result) <- c("pearson r","p.value")
rownames(result) <- colnames(b)
sig.result <- result[result[,2]<=0.1,]
sig.result

##
rooted_tree <- read.tree("1_28sp_phylomakers1_abbr.nwk")
plot(rooted_tree)

#rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data , names.col = "Plant.species", vcv = TRUE, vcv.dim = 2)
nrow(comp_data$data)

##
lm_model <- lm( PCoA1 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  PCoA1 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model11 <- predict(pgls_model, newdata = data )  


###
lm_model <- lm( P.concentrations ~ PCoA1 , data = data )
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ PCoA1 , data = comp_data, lambda = "ML")
summary(pgls_model)
pgls_model22 <- predict(pgls_model, newdata = data )


##
lm_model <- lm( Polonospora ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Polonospora ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model11 <- predict(pgls_model, newdata = data )  


###
lm_model <- lm( P.concentrations ~ Polonospora , data = data )
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Polonospora , data = comp_data, lambda = "ML")
summary(pgls_model)
pgls_model22 <- predict(pgls_model, newdata = data )


##
lm_model <- lm( Module.4 ~ Root.diameter , data = data ) 
summary(lm_model)
pgls_model <- pgls(  Module.4 ~ Root.diameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model11 <- predict(pgls_model, newdata = data )  


###
lm_model <- lm( P.concentrations ~ Module.4 , data = data )
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Module.4 , data = comp_data, lambda = "ML")
summary(pgls_model)
pgls_model22 <- predict(pgls_model, newdata = data )


