
######## Median calculation ----------------
library(plyr)
library(dplyr)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis/median calculation") 
data = read.csv("rd_981records.csv",header=TRUE, row.names=1)
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

#write.csv(data_median,"rd_981records_median_900sp.csv")


###
data = read.csv("rlc_2396records.csv",header=TRUE, row.names=1)
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

#write.csv(data_median,"rlc_2396records_median_1796sp.csv")


###
data = read.csv("mgr_2108records.csv",header=TRUE, row.names=1)
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
  SampleSize = sum(!is.na(MycorrhizalGrowthResponse))
)
data_median

#write.csv(data_median,"mgr_2108records_median_349sp.csv")

######## Plant phylogenetic signals ----------------
########  root diameter ----------------
library(ggtree)
library(ape)
library(picante)
library(phylosignal)
library(phylobase)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis") 
data <- read.csv("0_rd_781sp.csv", header=TRUE, row.names=1)
tree <- read.tree("0_rd_781sp.nwk")

rooted_tree1 <- root(tree, outgroup = "Equisetum", node = NULL, resolve.root = T, interactive = FALSE, edgelabel = FALSE)
plot(rooted_tree1)
rooted_tree  <- drop.tip(rooted_tree1, "Equisetum")
plot(rooted_tree)
sum(rooted_tree$edge.length == 0)
summary(rooted_tree$edge.length)
rooted_tree$edge.length[rooted_tree$edge.length == 0] <- 1e-15

all(sort(rooted_tree$label) == sort(rownames(data)))

tree_p <- as(rooted_tree, "phylo4")
class(tree_p)

p4d <- phylo4d(tree_p, data)
phy.sig <- phyloSignal(p4d, reps = 9999)
phy.sig


########  root length colonization ----------------
library(ggtree)
library(ape)
library(picante)
library(phylosignal)
library(phylobase)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis") 
data <- read.csv("0_rlc_1214sp.csv", header=TRUE, row.names=1)
tree <- read.tree("0_rlc_1214sp.nwk")

rooted_tree1 <- root(tree, outgroup = "Equisetum", node = NULL, resolve.root = T, interactive = FALSE, edgelabel = FALSE)
plot(rooted_tree1)
rooted_tree  <- drop.tip(rooted_tree1, "Equisetum")
plot(rooted_tree)
sum(rooted_tree$edge.length == 0)
summary(rooted_tree$edge.length)

all(sort(rooted_tree$label) == sort(rownames(data)))

tree_p <- as(rooted_tree, "phylo4")
class(tree_p)

p4d <- phylo4d(tree_p, data)
phy.sig <- phyloSignal(p4d, reps = 9999)
phy.sig


########  mycorrhizal growth response ----------------
library(ggtree)
library(ape)
library(picante)
library(phylosignal)
library(phylobase)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis") 
data <- read.csv("0_mgr_317sp.csv", header=TRUE, row.names=1)
tree <- read.tree("0_mgr_317sp.nwk")

rooted_tree1 <- root(tree, outgroup = "Equisetum", node = NULL, resolve.root = T, interactive = FALSE, edgelabel = FALSE)
plot(rooted_tree1)
rooted_tree  <- drop.tip(rooted_tree1, "Equisetum")
plot(rooted_tree)
sum(rooted_tree$edge.length == 0)
summary(rooted_tree$edge.length)

all(sort(rooted_tree$label) == sort(rownames(data)))

tree_p <- as(rooted_tree, "phylo4")
class(tree_p)

p4d <- phylo4d(tree_p, data)
phy.sig <- phyloSignal(p4d, reps = 9999)
phy.sig


########  Regression relationships among root diameter, root length colonization, and mycorrhizal growth response ----------------
library(readxl) 
library(pheatmap) 
library(psych)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggpmisc)
library(car)
library(caper)
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis") 

cols <- c(
  "Asteraceae" = "#E9746A",  "Campanulaceae" = "#E9446A",
  "Apiaceae" = "#B85A78", "Caprifoliaceae" = "#C45B7E", 
  "Plantaginaceae" = "#A83248",  "Orobanchaceae" = "#BC356C", 
  "Scrophulariaceae" = "#AF4D2A", "Lamiaceae" = "#D94F30", 
  "Verbenaceae" = "#C94F70", "Mazaceae" = "#bc4235", "Acanthaceae" = "#cf6155",  "Solanaceae" = "#AA67B3", 
  "Rubiaceae" = "#9e5c93", "Gentianaceae" = "#9e8c93",  "Apocynaceae" = "#9e6c83", 
  "Balsaminaceae" = "#DE4C7A",  "Primulaceae" ="#DE6C7A",
  "Caryophyllaceae" = "#7FC1B5", "Polygonaceae" = "#6BA89F",  "Plumbaginaceae" = "#497d75",
  "Amaranthaceae" = "#3A9E6E", "Fabaceae" = "#FFA64D", "Rosaceae" = "#E0A16F", "Urticaceae" = "#E0A92F",
  "Onagraceae" = "#9A8A8A",  "Brassicaceae" = "#B09AA8",
  "Lythraceae" = "#15969D", "Hypericaceae" = "#55948D",  "Violaceae"= "#9BAD9E", 
  "Euphorbiaceae" = "#A9B267", "Linaceae" = "#A3B150", "Oxalidaceae" = "#9BBF6E",  
  "Geraniaceae" = "#D96C4F", "Ranunculaceae" = "#BFA87F", "Berberidaceae"="#bf887f", "Papaveraceae" = "#BCB87F", 
  "Boraginaceae" = "#BFA93F", "Biebersteiniaceae" ="#BF8F98", 
  "Poaceae" = "#8FB9D1",  "Cyperaceae" = "#5F97C6", "Juncaceae" = "#8ecad0", 
  "Amaryllidaceae" = "#6C8FA6",  "Asparagaceae" = "#8FB6D4",   "Hypoxidaceae" ="#5F77C6",
  "Juncaginaceae" = "#4A6B7C", "Araceae" = "#4A6B9C",  "Tofieldiaceae"= "#4A6B9E",
  "Commelinaceae" ="#58b1ba", "Liliaceae"="#11b7Ca",
  "Musaceae"= "#5880ba", "Equisetaceae" = "#7B8F5E", "Dennstaedtiaceae" = "#5F8E4C", "Onocleaceae" = "#567F5D")


###
data1 = read.csv("1_rd_rlc_395sp.csv",header=TRUE)

fig5a <- ggplot(data1, aes(y = log(RootLengthColonization+1), x = log(RootDiameter)   )) +     
  geom_point(size =  2.5, alpha = 1, shape = 16,  aes(color = PlantFamily), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")),  formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.02,label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  scale_y_continuous(name ="ln(Root length colonization [%] +1)") + scale_x_continuous(name ="ln(Root diameter [mm])") + 
  labs(title ="395 plant species") +  scale_color_manual(values = cols )+  labs(color="Plant family")+
  theme_bw()+ theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                    legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.15),  axis.title=element_text(size=12)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", label.hjust = 0))

fig5a


####
data2 = read.csv("1_rlc_mgr_212sp.csv",header=TRUE)

fig5b <- ggplot(data2, aes( y = MycorrhizalGrowthResponse, x = log(RootLengthColonization +1) )) +     
  geom_point(size =  2.5, alpha = 1, shape = 16,  aes(color = PlantFamily), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")),  formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.02,label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  scale_y_continuous(name ="MGR")+ scale_x_continuous(name ="ln(Root length colonization [%] +1)") + 
  labs(title ="212 plant species") + scale_color_manual(values = cols )+  labs(color="Plant family")+
  theme_bw()+ theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                    legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.15),  axis.title=element_text(size=12)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", label.hjust = 0))

fig5b


###
data3 = read.csv("1_rd_mgr_162sp.csv", header=TRUE)

fig5c <- ggplot(data3, aes( y = MycorrhizalGrowthResponse, x = log(RootDiameter)  )) +    
  geom_point(size =  2.5, alpha = 1, shape = 16,  aes(color =PlantFamily), show.legend = F ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")),  formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.02,label.y = 0.98, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  scale_y_continuous(name ="MGR")+ scale_x_continuous(name ="ln(Root diameter [mm])") + 
  labs(title ="162 plant species") +  scale_color_manual(values = cols )+  labs(color="Plant family")+
  theme_bw()+ theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                    legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0.15),  axis.title=element_text(size=12)) +
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", label.hjust = 0))

fig5c

#1_fig5_10*4
fig5 <- cowplot::plot_grid( fig5a, fig5b, fig5c,
                  ncol = 3, nrow = 1, labels = c('(a)', "(b)", "(c)"), 
                  label_fontface = "bold", label_x = -0.0, label_y = 1.01 , align = c("hv"))
fig5

###
data_all <- bind_rows(data1, data2, data3)

legend <- ggplot(data_all, aes( y = MycorrhizalGrowthResponse , x = log(RootLengthColonization+1) , color = PlantFamily  )) +      
  geom_point(size =  2.5, alpha = 1, shape = 16,  aes(color = PlantFamily ), show.legend = T ) + 
  geom_smooth(method = lm, colour = "black", fill = "grey77",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" )+ 
  stat_poly_eq( aes(label = paste("atop(",after_stat(rr.label),",",after_stat(p.value.label),")", sep = "")),  formula = y ~ poly(x, 1,raw = TRUE), parse = TRUE,label.x = 0.98,label.y = 0.02, colour = "black",size=4.2, rr.digits = 3, p.digits = 3 )+
  scale_y_continuous(name ="")+ scale_x_continuous(name ="") + 
  labs(title ="") +  scale_color_manual(values = cols, name = "Plant family" ,  
                                         breaks = c( "Asteraceae", "Campanulaceae", "Apiaceae", "Caprifoliaceae","Plantaginaceae", "Orobanchaceae", "Scrophulariaceae", "Lamiaceae",
                                         "Verbenaceae", "Mazaceae", "Acanthaceae", "Solanaceae", "Rubiaceae", "Gentianaceae", "Apocynaceae", "Balsaminaceae", 
                                         "Primulaceae", "Caryophyllaceae", "Polygonaceae","Plumbaginaceae", "Amaranthaceae", "Fabaceae",  "Rosaceae" , "Urticaceae", 
                                         "Onagraceae", "Brassicaceae", "Lythraceae", "Hypericaceae" , "Violaceae", "Euphorbiaceae", "Linaceae", "Oxalidaceae",
                                         "Geraniaceae", "Ranunculaceae", "Berberidaceae", "Papaveraceae", "Boraginaceae", "Biebersteiniaceae", "Poaceae", "Cyperaceae", 
                                         "Juncaceae", "Amaryllidaceae",  "Asparagaceae", "Hypoxidaceae", "Juncaginaceae", "Araceae", "Tofieldiaceae", "Commelinaceae", 
                                         "Musaceae", "Equisetaceae", "Dennstaedtiaceae", "Onocleaceae" ) ) + 
  theme_bw()+ theme(axis.text=element_text(size=11), axis.text.y = element_text(angle = 90, hjust = 0.5), legend.position = "bottom",
                    legend.text=element_text(size=11), legend.title = element_text(size=12, hjust=0),  axis.title=element_text(size=12)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "left", label.position = "right", ncol = 6  ))
legend




######## Phylogenetic generalized least squares regressions ----------------
library(car)
library(caper) 
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis") 

###
data = read.csv("2_rd_rlc_376sp.csv",header=TRUE, row.names=1)
tree <- read.tree("2_rd_rlc_376sp.nwk")

tree$tip.label <- gsub("_", " ", tree$tip.label)
rooted_tree1 <- root(tree, outgroup = "Equisetum", node = NULL, resolve.root = T, interactive = FALSE, edgelabel = FALSE)
plot(rooted_tree1)
rooted_tree  <- drop.tip(rooted_tree1, "Equisetum")
plot(rooted_tree)

all(sort(rooted_tree$label) == sort(rownames(data)))

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data , names.col = "PlantSpecies", vcv = TRUE, vcv.dim = 2)

outlier_test = outlierTest(lm( RootLengthColonization ~  RootDiameter , data =  data)) 
outlier_test
pgls_model <- pgls( RootLengthColonization  ~ RootDiameter , data = comp_data, lambda = "ML") 
summary(pgls_model)


###
data = read.csv("2_rlc_mgr_205sp.csv",header=TRUE, row.names=1)
tree <- read.tree("2_rlc_mgr_205sp.nwk")

tree$tip.label <- gsub("_", " ", tree$tip.label)
rooted_tree1 <- root(tree, outgroup = "Equisetum", node = NULL, resolve.root = T, interactive = FALSE, edgelabel = FALSE)
plot(rooted_tree1)
rooted_tree  <- drop.tip(rooted_tree1, "Equisetum")
plot(rooted_tree)

all(sort(rooted_tree$label) == sort(rownames(data)))

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data , names.col = "PlantSpecies", vcv = TRUE, vcv.dim = 2)

outlier_test = outlierTest(lm( MycorrhizalGrowthResponse ~  RootLengthColonization , data =  data)) 
outlier_test
pgls_model <- pgls( MycorrhizalGrowthResponse ~  RootLengthColonization , data = comp_data, lambda = "ML") 
summary(pgls_model)


###
data = read.csv("2_rd_mgr_159sp.csv",header=TRUE, row.names=1)
tree <- read.tree("2_rd_mgr_159sp.nwk")

tree$tip.label <- gsub("_", " ", tree$tip.label)
rooted_tree1 <- root(tree, outgroup = "Equisetum", node = NULL, resolve.root = T, interactive = FALSE, edgelabel = FALSE)
plot(rooted_tree1)
rooted_tree  <- drop.tip(rooted_tree1, "Equisetum")
plot(rooted_tree)

all(sort(rooted_tree$label) == sort(rownames(data)))

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data , names.col = "PlantSpecies", vcv = TRUE, vcv.dim = 2)

outlier_test = outlierTest(lm( MycorrhizalGrowthResponse ~  RootDiameter , data =  data)) 
outlier_test
pgls_model <- pgls( MycorrhizalGrowthResponse ~  RootDiameter , data = comp_data, lambda = "ML") 
summary(pgls_model)
