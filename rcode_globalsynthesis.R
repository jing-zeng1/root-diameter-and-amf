
######## Median calculation ----------------
library(plyr)
library(dplyr)
rm(list = ls())
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

#write.csv(data_median,"0_rd_900sp.csv")


###
data = read.csv("rlc_2392records.csv",header=TRUE, row.names=1)
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

#write.csv(data_median,"0_rlc_1795sp.csv")


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

#write.csv(data_median,"0_mgr_348sp.csv")

##  v.phylomaker tree  rlc  1795sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis") 
df <- read.csv("1795sp.csv")
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

#write.tree(sp_tr, "0_rlc_1795sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL

is.rooted(sp_tr)       
is.binary.tree(sp_tr)  


##  v.phylomaker tree   rd  900sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis") 
df <- read.csv("900sp.csv")
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

#write.tree(sp_tr, "0_rd_900sp_phylomakers1.nwk")

#sp_tr$node.label <- NULL


is.rooted(sp_tr)       
is.binary.tree(sp_tr)  


##  v.phylomaker tree   mgr  348sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis") 
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


##  v.phylomaker tree    395sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis") 
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


##  v.phylomaker tree    211sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis")
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



##  v.phylomaker tree    162sp  --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis") 
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


######## -- Table S2 -- Regression relationships among root diameter, root length colonization,  and mycorrhizal growth response ----------------
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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF2/global synthesis_d3")

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


fig5a <- ggplot( data1, aes( y = log(RootLengthColonization+1) , x = log(RootDiameter)  ))+      
  scale_y_continuous(name ="ln(Root length colonization [%] +1)") + scale_x_continuous(name ="ln(Root diameter [mm])") + 
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  #geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 1 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci1 ,  aes(x = log(RootDiameter) , ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model1 ), linewidth=1, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size = 10.5 ), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size = 10.5 ), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1 ))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.1, vjust = -3.1, colour = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.003, 0.324 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.09, vjust = -1.6, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'<'~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.001, 0.815 ), size = 4, inherit.aes = FALSE )

fig5a


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
pgls_ci2 <- data.frame( RootLengthColonization = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )

fig5b <- ggplot(data2, aes(y = MycorrhizalGrowthResponse  , x = log(RootLengthColonization+1)  ))+      
  scale_y_continuous(name ="MGR") + scale_x_continuous(name ="ln(Root length colonization [%] +1)") + 
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, colour = "grey22", fill = "grey88", data = data2, formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci2 ,  aes(x = log(RootLengthColonization+1) , ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model2 ), linewidth= 0.8 , color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text( size=10.5 ), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text( size=10.5 ), legend.title = element_text(hjust=0.8),  axis.title=element_text( size = 12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 1.5,  colour = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.061, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.03, vjust = 2.7, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.054, 0.001 ), size = 4, inherit.aes = FALSE )

fig5b


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
pgls_ci3 <- data.frame( RootDiameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig5c <- ggplot( data3, aes(y = MycorrhizalGrowthResponse , x = log(RootDiameter)  ))+      
  scale_y_continuous(name ="MGR") + scale_x_continuous(name ="ln(Root diameter [mm])") + 
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, colour = "grey22", fill = "grey88", data = data3, formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci3 ,  aes(x = log(RootDiameter), ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes( y = pgls_model3 ), linewidth = 0.8, color = "#CB4641", linetype = "solid" ) + 
  theme_classic() + theme(axis.text = element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) + 
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.1, vjust = 2.4, color = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.113, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.09, vjust = 3.6, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.069, 0.001 ), size = 4, inherit.aes = FALSE )

fig5c


#1_fig5_10.8*4
fig5 <- cowplot::plot_grid( fig5a, fig5b, fig5c, 
                            ncol = 3, nrow = 1, labels = c('(a)', "(b)", "(c)"), 
                            label_fontface = "bold", label_x = -0.02, label_y = 1.02 , align = c("hv"))
fig5

eoffice::topptx( fig5 , width = 10.8 ,  height= 4 , "1_fig5_d3_10.8*4.pptx" ) 



######## -- Table S3 -- Plant phylogenetic signals   ----------------
library(ggtree)
library(ape)
library(picante)
library(phylosignal)
library(phylobase)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis") 
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
data <- read.csv("0_rd_900sp.csv", header=TRUE, row.names=1)
rooted_tree <- read.tree("0_rd_900sp_phylomakers1.nwk")

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

#all(sort(rooted_tree$label) == sort(rownames(data)))

tree_p <- as(rooted_tree, "phylo4")
class(tree_p)

p4d <- phylo4d(tree_p, data)
phy.sig <- phyloSignal(p4d, reps = 9999)
phy.sig


####
data <- read.csv("0_rlc_1795sp.csv", header=TRUE, row.names=1)
rooted_tree <- read.tree("0_rlc_1795sp_phylomakers1.nwk")

rooted_tree$tip.label <- gsub("_", " ", rooted_tree$tip.label)

#all(sort(rooted_tree$label) == sort(rownames(data)))

tree_p <- as(rooted_tree, "phylo4")
class(tree_p)

p4d <- phylo4d(tree_p, data)
phy.sig <- phyloSignal(p4d, reps = 9999)
phy.sig


######## -- Fig 4 -- Regression relationships among root diameter, root length colonization,  and mycorrhizal growth response ----------------
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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/global synthesis") 

### 395 sp 
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


fig5a <- ggplot( data1, aes( y = log(RootLengthColonization+1) , x = log(RootDiameter)  ))+      
  scale_y_continuous(name ="ln(Root length colonization [%] +1)") + scale_x_continuous(name ="ln(Root diameter [mm])") + 
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 1 , linetype = "solid"  )+ 
  #geom_ribbon(data = pgls_ci1 ,  aes(x = log(RootDiameter) , ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model1 ), linewidth=1, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size = 10.5 ), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size = 10.5 ), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1 ))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 1.4, colour = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.077, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.03, vjust = 2.6, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'<'~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.001, 0.981 ), size = 4, inherit.aes = FALSE )

fig5a


### 211 sp 
data2 = read.csv("1_rlc_mgr_211sp.csv",header=TRUE)
data2.1 = read.csv("1_rlc_mgr_211sp_2.csv",header=TRUE)
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
pgls_model2 <- predict(pgls_model, newdata = data2 )  

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
pgls_ci2 <- data.frame( RootLengthColonization = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )

fig5b <- ggplot(data2, aes(y = MycorrhizalGrowthResponse  , x = log(RootLengthColonization+1)  ))+      
  scale_y_continuous(name ="MGR") + scale_x_continuous(name ="ln(Root length colonization [%] +1)") + 
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, colour = "grey22", fill = "grey88", data = data2, formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci2 ,  aes(x = log(RootLengthColonization+1) , ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model2 ), linewidth= 0.8 , color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text( size=10.5 ), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text( size=10.5 ), legend.title = element_text(hjust=0.8),  axis.title=element_text( size = 12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 1.5,  colour = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.081, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.03, vjust = 2.7, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.060, 0.001 ), size = 4, inherit.aes = FALSE )

fig5b


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
pgls_model3 <- predict(pgls_model, newdata = data3 )    

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
pgls_ci3 <- data.frame( RootDiameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig5c <- ggplot( data3, aes(y = MycorrhizalGrowthResponse , x = log(RootDiameter)  ))+      
  scale_y_continuous(name ="MGR") + scale_x_continuous(name ="ln(Root diameter [mm])") + 
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, colour = "grey22", fill = "grey88", data = data3, formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth= 0.8 , linetype = "solid"  )+ 
  geom_ribbon(data = pgls_ci3 ,  aes(x = log(RootDiameter), ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes( y = pgls_model3 ), linewidth = 0.8, color = "#CB4641", linetype = "solid" ) + 
  theme_classic() + theme(axis.text = element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
 # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) + 
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = -2, color = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.125, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1, vjust = -0.6, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.060, 0.002 ), size = 4, inherit.aes = FALSE )

fig5c


#1_0529_fig5_10.8*3.8
fig5 <- cowplot::plot_grid( fig5a, fig5b, fig5c, 
                            ncol = 3, nrow = 1, labels = c('(a)', "(b)", "(c)"), 
                            label_fontface = "bold", label_x = -0.02, label_y = 1.01 , align = c("hv"))
fig5

eoffice::topptx( fig5 , width = 9.2 ,  height= 3.6 , "1_fig4_9.2*3.6.pptx" ) 



