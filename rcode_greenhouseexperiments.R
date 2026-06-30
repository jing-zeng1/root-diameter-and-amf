
##  v.phylomaker tree    45sp --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments") 
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

#write.tree(sp_tr, "1_45sp_phylomakers1.nwk")

sp_tr$node.label <- NULL
#write.tree(sp_tr, "1_45sp_phylomakers1_mega.nwk")

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

#write.tree( tree, file = "1_45sp_phylomakers1_abbr.nwk")

##  v.phylomaker tree    37sp --------------------------
#devtools::install_github("jinyizju/V.PhyloMaker2")
library(picante)
library(V.PhyloMaker2)
library(ggtree)
library(ape)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments") 
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

sp_tr$node.label <- NULL
#write.tree(sp_tr, "1_45sp_phylomakers1_mega.nwk")

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

######## -- Table S2 37sp -- Regression relationships of MPR and MGR with root length colonization, root diameter, and MRR ----------------
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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments")  
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
pgls_ci <- data.frame( Root.length.colonization = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig4a <- ggplot(data, aes(y = MPR , x = Root.length.colonization  )) +      
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = "MPR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  geom_ribbon(data = pgls_ci,  aes(x = Root.length.colonization, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model1 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.5,  colour = "grey22", 
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.124, 0.033 ), size = 4 , inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.6, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.084, 0.081 ), size = 4 , inherit.aes = FALSE )

fig4a


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
pgls_ci <- data.frame( Root.length.colonization = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig4b <- ggplot(data, aes(y = MGR , x = Root.length.colonization  ))+      
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = "MGR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  #geom_ribbon(data = pgls_ci,  aes(x = Root.length.colonization, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model2 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text( size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.4 , color = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.141, 0.022 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 2.6 , color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.002, 0.805 ), size = 4, inherit.aes = FALSE )

fig4b


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


fig4c <- ggplot(data, aes(y = MPR , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "MPR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  geom_ribbon(data = pgls_ci3 ,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model3 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.2, color = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.165, 0.013 ), size = 4 , inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.4, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.119, 0.037 ), size = 4 , inherit.aes = FALSE )

fig4c


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
pgls_ci <- data.frame( Root.diameter = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


fig4d <- ggplot(data, aes(y = MGR , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "MGR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  #geom_ribbon(data = pgls_ci,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model4 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.03, vjust =  1.2, color = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.240, 0.002 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = - 0.02, vjust =  2.5, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.043, 0.220 ), size = 4, inherit.aes = FALSE )

fig4d


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


fig4e <- ggplot(data, aes(y = MPR , x = MRR  )) +      
  scale_x_continuous(name ="MRR") + scale_y_continuous(name = "MPR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  geom_ribbon(data = pgls_ci5,  aes(x = MRR, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model5 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.4, color = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.139, 0.023 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.6, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.093, 0.066 ), size = 4, inherit.aes = FALSE )

fig4e


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


fig4f <- ggplot(data, aes(y = MGR , x = MRR  ))+      
  scale_x_continuous(name ="MRR") + scale_y_continuous(name = "MGR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  geom_ribbon(data = pgls_ci6,  aes(x = MRR, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model6 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                          legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.15, size = 13),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", nrow= 12, label.hjust = 0))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.05, vjust = - 1.8 , color = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.603, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.05, vjust = - 0.4 , color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.495, 0.001 ), size = 4, inherit.aes = FALSE )

fig4f


# 1_0529_fig4_10.8*7.6
fig4 <- cowplot::plot_grid( fig4a, fig4c, fig4e, fig4b, fig4d, fig4f, 
                            ncol = 3, nrow = 2, labels = c('(a)', "(c)", "(e)",'(b)', "(d)", "(f)"), 
                            label_fontface = "bold", label_x = -0.02, label_y = 1.02 , align = c("hv"))
fig4


#1_fig4_10.8*7.2
#fig4 <- cowplot::plot_grid( pp1, legend, align = "", axis = "tpl", rel_widths = c( 1.7, 0.3 )) 
#fig4

#eoffice::topptx( fig4 , width = 8.64 ,  height = 6.5 , "1_37sp_fig4_8.64*6.5.pptx" ) 








######## -- Table S3 45sp -- Plant phylogenetic signals of plant traits and AM fungal variables ----------------
library(plyr) 
library(ggtree)
library(ape)
library(picante)
library(phylosignal)
library(phylobase)
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments")  
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

write.csv(data_mean,"greenhouse_data_mean.csv")

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

######## -- Fig S6 45sp -- Relationship between plant traits and AM fungal colonization ---------------
library(ggplot2)
library(gridExtra)
library(reshape)
library(ggpmisc)
library(caper)
library(car) 
library(doBy)
library(dplyr) 
rm(list = ls())
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments")  
data = read.csv("greenhouse_data_mean.csv",header=TRUE, row.names=1)

##
rooted_tree <- read.tree("1_45sp_phylomakers1.nwk")
plot(rooted_tree)

data$Plant.species0 <- gsub(" ", "_", data$Plant.species0)

rooted_tree$node.label
rooted_tree$node.label <- NULL
comp_data <- comparative.data( phy = rooted_tree, data = data , names.col = "Plant.species0", vcv = TRUE, vcv.dim = 2)

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


figs5a <- ggplot(data, aes(y = Root.length.colonization , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "Root length colonization (%)" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  #geom_ribbon(data = pgls_ci1,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model1 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
 # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.05, vjust = - 3.9, colour = "grey22", 
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.098, 0.036 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.05, vjust = - 2.3, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.005, 0.629 ), size = 4, inherit.aes = FALSE )

figs5a


###
lm_model <- lm( P.concentrations ~ Root.length.colonization , data = data) 
summary(lm_model)
pgls_model <- pgls(  P.concentrations ~ Root.length.colonization , data = comp_data, lambda = "ML") 
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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( P.concentrations ~ Root.length.colonization ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.length.colonization = pred_grid) ) } } }

# CI
pgls_ci2 <- data.frame( Root.length.colonization = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


figs5b <- ggplot(data, aes(y = P.concentrations , x = Root.length.colonization  ))+      
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = expression("Tissue P concentration (mg g"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "dashed" , data = data )+ 
  #geom_ribbon(data = pgls_ci2,  aes(x = Root.length.colonization, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model2 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  #theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 2.2, color = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.071, 0.076 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 3.5, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.019, 0.371 ), size = 4, inherit.aes = FALSE )
  
figs5b


###
lm_model <- lm( Total.biomass ~ Root.length.colonization , data = data) 
summary(lm_model)
pgls_model <- pgls(  Total.biomass ~ Root.length.colonization , data = comp_data, lambda = "ML") 
summary(pgls_model)
pgls_model3 <- predict(pgls_model, newdata = data)  

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
if (!inherits(comp_boot, "try-error")) {  m <- try( pgls( Total.biomass ~ Root.length.colonization ,  data = comp_boot, lambda = "ML"), silent = TRUE )
if (!inherits(m, "try-error")) {  boot_mat[i, ] <- predict(  m,  newdata = data.frame( Root.length.colonization = pred_grid) ) } } }

# CI
pgls_ci3 <- data.frame( Root.length.colonization = pred_grid, fit = apply(boot_mat, 2, mean, na.rm = TRUE), lwr = apply(boot_mat, 2, quantile, 0.025, na.rm = TRUE), upr = apply(boot_mat, 2, quantile, 0.975, na.rm = TRUE) )


figs5c <- ggplot(data, aes(y = Total.biomass , x = Root.length.colonization  ))+ 
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = expression("Total biomass (g plant"^{-1}*")")  )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  #geom_smooth(method = lm, colour = "black", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "dashed" , data = data )+ 
  #geom_ribbon(data = pgls_ci3,  aes(x = Root.length.colonization, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model3 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.15, size =13 ),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
 # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", nrow= 12, label.hjust = 0))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 2.3,color = "grey22",
             label = sprintf("italic('lm')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.018, 0.380 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 3.5, color = "#CB4641", 
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.005, 0.638 ), size = 4, inherit.aes = FALSE )

figs5c


#1_0529_figs5_10.8*3.6
figs5 <- cowplot::plot_grid( figs5a,  figs5b, figs5c,  
                             ncol = 3, nrow = 1, labels = c('(a)', "(b)", "(c) ",'(d)', "(e)"), 
                             label_fontface = "bold", label_x = -0.03, label_y = 1.01 , align = c("hv"))
figs5

eoffice::topptx( figs5 , width = 8.64 ,  height= 3.1 , "1_45sp_figs5_8.64*3.1.pptx" ) 


######## -- Fig 2 45sp -- Plant’s responses to mycorrhizal inoculation ----------------
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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments")  
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

#1_45sp_fig3_10*9
fig3 <- p1 %>%
  insert_left(p00, width = 0.8) %>%          
  insert_right(blank, width = 0.03) %>%      
  insert_right(p2,   width = 1) %>%         
  insert_right(blank, width = 0.03) %>%     
  insert_right(p3,   width = 1)

fig3

eoffice::topptx( fig3 , width = 10 ,  height = 9 , "1_118_28sp_fig2_10*9.pptx" ) 


######## -- Fig 3 45sp -- Regression relationships of MPR and MGR with root length colonization, root diameter, and MRR ----------------
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
setwd("/Users/zengjing/Desktop/Rootdiameter and AMF/greenhouse experiments")  
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


fig4a <- ggplot(data, aes(y = MPR , x = Root.length.colonization  )) +      
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = "MPR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, colour = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  geom_ribbon(data = pgls_ci1,  aes(x = Root.length.colonization, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model1 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
 # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.4,  colour = "grey22", 
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.299, 0.001 ), size = 4 , inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.6, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.211, 0.001 ), size = 4 , inherit.aes = FALSE )

fig4a


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


fig4b <- ggplot(data, aes(y = MGR , x = Root.length.colonization  ))+      
  scale_x_continuous(name ="Root length colonization (%)") + scale_y_continuous(name = "MGR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  #geom_ribbon(data = pgls_ci2,  aes(x = Root.length.colonization, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model2 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text( size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
#  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.2 , color = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.258, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 2.4 , color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.019, 0.362 ), size = 4, inherit.aes = FALSE )

fig4b


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


fig4c <- ggplot(data, aes(y = MPR , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "MPR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  geom_ribbon(data = pgls_ci3,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model3 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5 )) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
 # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.1, color = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.197, 0.002 ), size = 4 , inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.3, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.120, 0.020 ), size = 4 , inherit.aes = FALSE )

fig4c


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


fig4d <- ggplot(data, aes(y = MGR , x = Root.diameter  ))+      
  scale_x_continuous(name ="Root diameter (mm)") + scale_y_continuous(name = "MGR" )+
    geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  #geom_ribbon(data = pgls_ci4,  aes(x = Root.diameter, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  #geom_line(aes(y = pgls_model4 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
#  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.05, vjust = - 1.7, color = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.258, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.05, vjust = - 0.4, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.061, 0.101 ), size = 4, inherit.aes = FALSE )

fig4d



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


fig4e <- ggplot(data, aes(y = MPR , x = MRR  )) +      
  scale_x_continuous(name ="MRR") + scale_y_continuous(name = "MPR" )+
    geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) + 
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  geom_ribbon(data = pgls_ci5,  aes(x = MRR, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model5 ), linewidth=0.8, color = "#CB4641", linetype = "dashed") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.8),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
 # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "left", label.hjust = 1))+
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.05, vjust = 1.4, color = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.193, 0.003 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = -Inf, y = Inf),  aes(x = x, y = y), parse = TRUE, hjust = -0.04, vjust = 2.6, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'='~'%.3f'", 0.072, 0.074 ), size = 4, inherit.aes = FALSE )

fig4e


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


fig4f <- ggplot(data, aes(y = MGR , x = MRR  ))+      
  scale_x_continuous(name ="MRR") + scale_y_continuous(name = "MGR" )+
  geom_point( size = 2.5,  shape = 21,  fill = "#d3deea" , color = "grey22", show.legend = F ) + scale_size(range = c(1, 6)) +
  geom_smooth(method = lm, color = "grey22", fill = "grey88",  formula = y ~ poly(x, 1,raw = TRUE), se = T, linewidth=0.8 , linetype = "solid" , data = data )+ 
  geom_ribbon(data = pgls_ci6,  aes(x = MRR, ymin = lwr, ymax = upr), fill = "#CB4641",  alpha = 0.1, inherit.aes = FALSE) +
  geom_line(aes(y = pgls_model6 ), linewidth=0.8, color = "#CB4641", linetype = "solid") + 
  theme_classic() + theme(axis.text=element_text(size=10.5), axis.text.y = element_text(angle = 90, hjust = 0.5), 
                     legend.text=element_text(size=10.5), legend.title = element_text(hjust=0.15, size = 13),  axis.title=element_text(size=12.5)) + 
  theme(plot.background = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
 # theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55)) +
  guides(color = guide_legend(override.aes = list(size = 2.5), title.position = "top", label.position = "right", nrow= 12, label.hjust = 0))+
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.05, vjust = - 1.7, color = "grey22",
             label = sprintf("italic('ols')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.617, 0.001 ), size = 4, inherit.aes = FALSE ) + 
  geom_text( data = data.frame(x = Inf, y = -Inf),  aes(x = x, y = y), parse = TRUE, hjust = 1.05, vjust = - 0.4, color = "#CB4641",
             label = sprintf("italic('pgls')*': '*italic(R^2)~'='~'%.3f'*', '*italic(P)~'<'~'%.3f'", 0.489, 0.001 ), size = 4, inherit.aes = FALSE )

fig4f


# 1_0529_fig4_10.8*7.6
fig4 <- cowplot::plot_grid( fig4a, fig4c, fig4e,  fig4b,  fig4d,  fig4f, 
                   ncol = 3, nrow = 2, labels = c('(a)', "(c)", "(e)",'(b)', "(d)", "(f)"), 
                   label_fontface = "bold", label_x = -0.02, label_y = 1.02 , align = c("hv"))
fig4


#1_fig4_10.8*7.2
#fig4 <- cowplot::plot_grid( pp1, legend, align = "", axis = "tpl", rel_widths = c( 1.7, 0.3 )) 
#fig4

eoffice::topptx( fig4 , width = 8.64 ,  height = 6.5 , "1_45sp_fig4_8.64*6.5.pptx" ) 

 


### Regression relationships of MPR and MGR with root diameter under nonmycorrhizal conditions
colnames(data)
outlier_test = outlierTest(lm( MPR ~ Root.diameter_NM , data =  data )) 
outlier_test
cor_test_result <- cor.test( data$MPR, data$Root.diameter_NM, method = "pearson")
cor_test_result$estimate
cor_test_result$p.value 

outlier_test = outlierTest(lm( MGR ~ Root.diameter_NM , data =  data )) 
outlier_test
cor_test_result <- cor.test( data$MGR, data$Root.diameter_NM, method = "pearson")
cor_test_result$estimate
cor_test_result$p.value 



